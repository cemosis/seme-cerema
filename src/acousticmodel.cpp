
#include <acousticmodel.hpp>

namespace Feel
{

AcousticModel::AcousticModel()
{
    M_verbose = boption(_name="acoustic-model.verbose");

    M_vSoundVelocity = doption(_name="sound-velocity");
    M_refPressure = doption(_name="reference-pressure");
    M_airDensity = doption(_name="air-density");
    M_M = doption(_name="coeff.M");
    M_alpha = doption(_name="coeff.alpha");
    M_d_prob = doption(_name="coeff.d-prob");//0.5;

    M_useOneMatVecByDirection = boption(_name="use-one-matvec-by-direction");
    M_useOneBackendByDirection = boption(_name="use-one-backend-by-direction");
    M_hasBuildCstPart = false;

    this->createThetaPhi();

}

void
AcousticModel::createThetaPhi()
{
    // in parallel start by master rank (write msh file), and load after others process
    std::string fileNameMeshDescThetaPhi = Environment::expand( soption(_name="gmsh.filename-thetaPhi") );
    if ( Environment::isMasterRank() )
        M_meshThetaPhi = Feel::loadMesh(_mesh=new mesh_thetaphi_type, _filename=fileNameMeshDescThetaPhi,
                                        _worldcomm=Environment::worldCommSeq());
    Environment::worldComm().barrier();
    if ( !Environment::isMasterRank() )
    {
        std::string fileNameMshThetaPhi = fs::path( fileNameMeshDescThetaPhi ).stem().string()+".msh";
        //std::cout << "filename to load " << fileNameMshThetaPhi << "\n";
        M_meshThetaPhi = Feel::loadMesh(_mesh=new mesh_thetaphi_type, _filename=fileNameMshThetaPhi,
                                        _worldcomm=Environment::worldCommSeq());
    }

    // build function space for theta phi
    M_XhThetaPhi = space_thetaphi_type::New(_mesh=M_meshThetaPhi,
                                            _worldscomm=std::vector<WorldComm>(1,Environment::worldCommSeq()) );



    // find doublon in thetaphi functionspace and compute dof's periodicity
    //std::vector<bool> dofIsUsed_thetaPhi(M_XhThetaPhi->nDof(),true);
    M_dofIsUsed_thetaPhi.resize(M_XhThetaPhi->nDof(),true);
    //std::map<size_type,size_type> mapPeriodicDofThetaPhi;

    bool findThetaPhiLeft=false,findThetaPhiRight=false;
    size_type idThetaPhiLeft = invalid_size_type_value, idThetaPhiRight = invalid_size_type_value;
    for( size_type dof_thetaPhi=0 ; dof_thetaPhi<M_XhThetaPhi->nDof() ; ++dof_thetaPhi )
    {
        auto const& mythetaPhiPt = M_XhThetaPhi->dof()->dofPoint(dof_thetaPhi).get<0>();

        if ( mythetaPhiPt[1] < (2*M_PI-0.00001) )
        {

            if ( std::abs( mythetaPhiPt[0] ) < 1e-9 )
            {
                if ( !findThetaPhiLeft )
                {
                    M_dofToUseInThetaPhi.insert(dof_thetaPhi);
                    findThetaPhiLeft=true;
                    idThetaPhiLeft = dof_thetaPhi;
                }
                else
                {
                    M_dofIsUsed_thetaPhi[dof_thetaPhi] = false;
                    M_mapPeriodicDofThetaPhi[ dof_thetaPhi ] = idThetaPhiLeft;
                }
            }
            else if ( std::abs( mythetaPhiPt[0] - M_PI ) < 1e-9 )
            {
                if ( !findThetaPhiRight )
                {
                    M_dofToUseInThetaPhi.insert(dof_thetaPhi);
                    findThetaPhiRight=true;
                    idThetaPhiRight = dof_thetaPhi;
                }
                else
                {
                    M_dofIsUsed_thetaPhi[dof_thetaPhi] = false;
                    M_mapPeriodicDofThetaPhi[ dof_thetaPhi ] = idThetaPhiRight;
                }
            }
            else
            {
                M_dofToUseInThetaPhi.insert(dof_thetaPhi);
            }
        }
        else
        {
            M_dofIsUsed_thetaPhi[dof_thetaPhi] = false;
            // search corresponding dof
            bool hasFindPeriodicDof=false;
            for( size_type dof_thetaPhi2=0 ; dof_thetaPhi2<M_XhThetaPhi->nDof() && !hasFindPeriodicDof ; ++dof_thetaPhi2 )
            {
                auto const& mythetaPhiPt2 = M_XhThetaPhi->dof()->dofPoint(dof_thetaPhi2).get<0>();
                if ( std::abs( mythetaPhiPt[0]-mythetaPhiPt2[0] ) < 1e-9 && std::abs( mythetaPhiPt2[1] ) < 1e-9 )
                {
                    M_mapPeriodicDofThetaPhi[ dof_thetaPhi ] = dof_thetaPhi2;
                    hasFindPeriodicDof=true;
                }
            }
        }

    }

    if ( Environment::isMasterRank() )
        std::cout << "XhThetaPhi->nDof : = " << M_XhThetaPhi->nDof() << " (size with no doublon : " << M_dofToUseInThetaPhi.size() << ")\n";
}

void
AcousticModel::loadMesh( mesh_ptrtype mesh)
{
    M_mesh=mesh;

    this->createTransportFE();
    this->createPrecomputeBC();
}

void
AcousticModel::createTransportFE()
{
    // build space
    M_Xh = space_type::New(_mesh=M_mesh);
    //ULoc = Xh->elementPtr(); // ASUP
    //ULocSumAllDirection = Xh->elementPtr();
    M_ULocIntegrateAllDirection = M_Xh->elementPtr();
    M_projBCDirichlet = M_Xh->elementPtr();
    M_projSourceTerm = M_Xh->elementPtr();
    M_soundPressureLevel = M_Xh->elementPtr();

    M_ULoc_thetaPhi.resize(M_XhThetaPhi->nDof() );
    for( size_type dofThetaPhi : M_dofToUseInThetaPhi )
    {
        M_ULoc_thetaPhi[dofThetaPhi] = M_Xh->elementPtr();
        //*ULoc_thetaPhi[dofThetaPhi] = *ULoc;
    }


    M_bdfDensity.resize( M_XhThetaPhi->nDof() );
    for( size_type dofThetaPhi : M_dofToUseInThetaPhi )
        M_bdfDensity[dofThetaPhi] = bdf( _vm=Environment::vm(), _space=M_Xh,
                                  _name=(boost::format("energy%1%-%2%_%3%")%dofThetaPhi %Environment::worldComm().rank() %Environment::worldComm().size() ).str());

    if ( Environment::isMasterRank() )
        std::cout << "Xh->nDof : = " << M_Xh->nDof() << std::endl;

    M_exporter = exporter( _mesh=M_mesh,_name="myexporter" );


}


void
AcousticModel::createPrecomputeBC()
{
    auto pcf =  M_mesh->gm()->preComputeOnFaces( M_mesh->gm(), M_mesh->gm()->referenceConvex().barycenterFaces() );
    CHECK( M_mesh->beginFace()!=M_mesh->endFace() ) << "no face on mesh for this proc : must be take into account!\n";
    auto firstFacePtr = M_mesh->beginFace();
    auto ctx = M_mesh->gm()->context<vm::POINT|vm::NORMAL|vm::KB|vm::JACOBIAN>( firstFacePtr->element0(),
                                                                                pcf,
                                                                                firstFacePtr->pos_first() );
    //auto wBCintegrate = Xh_thetaPhi->elementPtr();
    //vitesse de propagation
    std::vector<double> vDirection(3,0.);

    M_vxProj=M_XhThetaPhi->elementPtr();
    M_vyProj=M_XhThetaPhi->elementPtr();
    M_vzProj=M_XhThetaPhi->elementPtr();
    M_vxProj->on( _range=elements(M_meshThetaPhi), _expr=M_vSoundVelocity*sin(Px())*cos(Py()) );
    M_vyProj->on( _range=elements(M_meshThetaPhi), _expr=M_vSoundVelocity*sin(Px())*sin(Py()) );
    M_vzProj->on( _range=elements(M_meshThetaPhi), _expr=M_vSoundVelocity*cos(Px()) );


    std::map<size_type,std::map<size_type,std::pair<node_type,node_type> > > mapThetaPhiAndHatThetaPhiNEW;
    //std::map<size_type,context_ptrtype> M_ctxEvalHatThetaPhi;//( new context_type( Xh_thetaPhi ) );
    //std::map<size_type,std::map<size_type,size_type> > M_mapEvalFromCtxThetaPhi;
    //std::map<size_type,std::map<size_type,bool> > M_doComputeHatThetaPhi;
    //std::map<size_type, node_type> mapUnitNormalByFace;
    //std::map<size_type,vector_ptrtype > M_mapOperatorDiffuseBC;
    std::set< boost::tuple<double,double,double,size_type> > setOfNormalVectorStored;

    size_type nBoundaryFaces = nelements(boundaryfaces(M_mesh),false);
    size_type cptBoundaryFaces = 0;
    bool useStorageNormalBoundaryFaces= boption(_name="use-storage-from-normal-boundary-faces"); // )true;
    for ( auto const& face : boundaryfaces(M_mesh) )
    {
        ++cptBoundaryFaces;
        std::cout << "\r";
        if ( Environment::isMasterRank() )
            std::cout << "cptBoundaryFaces = " << std::setw(10) << cptBoundaryFaces << "/" << nBoundaryFaces << std::flush;// << "\n";

        size_type faceId = face.id();
        ctx->update( face.element0(), face.pos_first() );
        auto unitNormal = ctx->unitNormal();

        if ( useStorageNormalBoundaryFaces )
        {
            bool findNormal = false;
            size_type faceIdReference = invalid_size_type_value;
            for ( auto const& normalVecStored : setOfNormalVectorStored )
            {
                if ( std::abs( normalVecStored.get<0>() - unitNormal[0] ) < 1e-8 &&
                     std::abs( normalVecStored.get<1>() - unitNormal[1] ) < 1e-8 &&
                     std::abs( normalVecStored.get<2>() - unitNormal[2] ) < 1e-8 )
                {
                    findNormal = true;
                    faceIdReference = normalVecStored.get<3>();
                    break;
                }
            }

            if ( findNormal )
            {
                M_ctxEvalHatThetaPhi[faceId] = M_ctxEvalHatThetaPhi.find(faceIdReference)->second;
                M_mapUnitNormalByFace[ faceId ] = unitNormal;
                M_doComputeHatThetaPhi[faceId] = M_doComputeHatThetaPhi.find(faceIdReference)->second;
                mapThetaPhiAndHatThetaPhiNEW[ faceId ] = mapThetaPhiAndHatThetaPhiNEW.find(faceIdReference)->second;
                M_mapEvalFromCtxThetaPhi[faceId] = M_mapEvalFromCtxThetaPhi.find(faceIdReference)->second;
                M_mapOperatorDiffuseBC[faceId] = M_mapOperatorDiffuseBC.find(faceIdReference)->second;
                // next face
                continue;
            }

            // say that data for this normal vector will be computed
            setOfNormalVectorStored.insert( boost::make_tuple( unitNormal[0],unitNormal[1],unitNormal[2],faceId ) );
        }


        M_ctxEvalHatThetaPhi[faceId].reset( new context_type( M_XhThetaPhi ) );
        M_mapUnitNormalByFace[ faceId ] = unitNormal;

        size_type cptNodeCtx = 0;
        for( size_type dof_thetaPhi : M_dofToUseInThetaPhi )
        {
            vDirection[0] = M_vxProj->operator()(dof_thetaPhi);
            vDirection[1] = M_vyProj->operator()(dof_thetaPhi);
            vDirection[2] = M_vzProj->operator()(dof_thetaPhi);

            double valVdotN = vDirection[0]*unitNormal[0] + vDirection[1]*unitNormal[1] + vDirection[2]*unitNormal[2];

            // on prend les vecteur speculaire partant de la paroi (on calcul les vecteurs arrivants sur la paroi )
            if ( valVdotN > -1e-6 )
            {
                M_doComputeHatThetaPhi[faceId][dof_thetaPhi]=false;
                continue;
            }

            M_doComputeHatThetaPhi[faceId][dof_thetaPhi]=true;

            // compute specular vector in cartesian coord (thank to tangent plane projector)
            node_type hatThetaPhiCartesian(3);
            hatThetaPhiCartesian[0] = vDirection[0] - 2*valVdotN*unitNormal[0];
            hatThetaPhiCartesian[1] = vDirection[1] - 2*valVdotN*unitNormal[1];
            hatThetaPhiCartesian[2] = vDirection[2] - 2*valVdotN*unitNormal[2];
            double normHatThetaPhiCartesian = math::sqrt(math::pow(hatThetaPhiCartesian[0],2)+math::pow(hatThetaPhiCartesian[1],2)+math::pow(hatThetaPhiCartesian[2],2));
            for ( uint16_type k=0;k<3;++k)
                hatThetaPhiCartesian[k] /= normHatThetaPhiCartesian;
            normHatThetaPhiCartesian = math::sqrt(math::pow(hatThetaPhiCartesian[0],2)+math::pow(hatThetaPhiCartesian[1],2)+math::pow(hatThetaPhiCartesian[2],2));
            CHECK( math::abs(normHatThetaPhiCartesian-1) < 1e-9 ) << "normHatThetaPhiCartesian = " << normHatThetaPhiCartesian;

            auto thetaPhiPt = M_XhThetaPhi->dof()->dofPoint(dof_thetaPhi).get<0>();
            node_type hatThetaPhi(2);
            hatThetaPhi[0] = math::acos(hatThetaPhiCartesian[2]/1);
            CHECK( hatThetaPhi[0] >= -1e-8 && hatThetaPhi[0] <= (M_PI+1e-8) ) << "hatTheta is not in [0,pi] : "  << hatThetaPhi[0];
            // see http://en.wikipedia.org/wiki/Atan2
            hatThetaPhi[1] = std::atan2(hatThetaPhiCartesian[1],hatThetaPhiCartesian[0]);
            // translate phi if not in [0,2*pi]
            if ( hatThetaPhi[1] < 0 )
                hatThetaPhi[1] += 2*M_PI;
            CHECK( hatThetaPhi[1] >= -1e-8 && hatThetaPhi[1] <= (2*M_PI+1e-8) ) << " hatPhi is not in [0,2pi] : "  << hatThetaPhi[1];

            node_type hatThetaPhiCheck(3);
            hatThetaPhiCheck[0] = math::sin(hatThetaPhi[0])*math::cos(hatThetaPhi[1]);
            hatThetaPhiCheck[1] = math::sin(hatThetaPhi[0])*math::sin(hatThetaPhi[1]);
            hatThetaPhiCheck[2] = math::cos(hatThetaPhi[0]);
            CHECK( ( std::abs( hatThetaPhiCartesian[0]-hatThetaPhiCheck[0] ) < 1e-9 ) &&
                   ( std::abs( hatThetaPhiCartesian[1]-hatThetaPhiCheck[1] ) < 1e-9 ) &&
                   ( std::abs( hatThetaPhiCartesian[2]-hatThetaPhiCheck[2] ) < 1e-9 ) ) << "aie : \n"
                                                                                        << "x =  " << hatThetaPhiCartesian[0] << " vs " << hatThetaPhiCheck[0] << "\n"
                                                                                        << "y =  " << hatThetaPhiCartesian[1] << " vs " << hatThetaPhiCheck[1] << "\n"
                                                                                        << "z =  " << hatThetaPhiCartesian[2] << " vs " << hatThetaPhiCheck[2] << "\n";

            //std::cout << " thetaPhiPt vs hatThetaPhi " << thetaPhiPt << " vs " << hatThetaPhi << "\n";
            mapThetaPhiAndHatThetaPhiNEW[ faceId ][ dof_thetaPhi ] = std::make_pair( thetaPhiPt,hatThetaPhi );

            M_ctxEvalHatThetaPhi[faceId]->add( hatThetaPhi );
            M_mapEvalFromCtxThetaPhi[faceId][dof_thetaPhi]=cptNodeCtx;
            ++cptNodeCtx;

        } // for (dof_thetaPhi)

        if ( true )
        {
            M_mapOperatorDiffuseBC[faceId] = backend()->newVector(M_XhThetaPhi);

            auto thetaPrime = Px();
            auto phiPrime = Py();
            auto vPrime = M_vSoundVelocity*vec( sin(thetaPrime)*cos(phiPrime), sin(thetaPrime)*sin(phiPrime), cos(thetaPrime) );
            auto vPrimeDotN = vPrime(0,0)*unitNormal[0] + vPrime(1,0)*unitNormal[1] + vPrime(2,0)*unitNormal[2];
#if 0
            auto chiVPrimeDotN = chi( vPrimeDotN > cst(-1e-6) );// not compile with linearform(maybe bilinear also
#else
            auto chiVPrimeDotN = chi( ( sin(thetaPrime)*cos(phiPrime)*unitNormal[0] + sin(thetaPrime)*sin(phiPrime)*unitNormal[1] + cos(thetaPrime)*unitNormal[2] )  > cst(-1e-6) );
#endif
            auto uThetaPhi = M_XhThetaPhi->element();
            auto bcIntegrateExpr = (M_alpha*M_d_prob/(M_PI*M_vSoundVelocity))*vPrimeDotN*id(uThetaPhi)*sin(thetaPrime);
            form1(_test=M_XhThetaPhi,_vector=M_mapOperatorDiffuseBC[faceId] )
                = integrate(_range=elements(M_meshThetaPhi),
                            _expr=bcIntegrateExpr*chiVPrimeDotN
                            );
            M_mapOperatorDiffuseBC[faceId]->close();
        }


    } // faces

    if ( Environment::isMasterRank() )
        std::cout << "\nwaiting other process (can be take a long time)" << std::flush;
    Environment::worldComm().barrier();
    if ( Environment::isMasterRank() )
        std::cout << "\rfinish precompute (size of normalVectorStored.size : " << setOfNormalVectorStored.size() << ")\n";


}


void
AcousticModel::init()
{
    if ( M_bdfDensity[*M_dofToUseInThetaPhi.begin()]->isRestart() )
    {
        for( size_type dofThetaPhi : M_dofToUseInThetaPhi )
        {
            M_bdfDensity[dofThetaPhi]->restart();
            // load a previous solution as current solution
            *M_ULoc_thetaPhi[dofThetaPhi] = M_bdfDensity[dofThetaPhi]->unknown(0);
        }
        // restart exporter
        if ( M_exporter->doExport() )
            M_exporter->restart( M_bdfDensity[*M_dofToUseInThetaPhi.begin()]->timeInitial() );

        if ( Environment::isMasterRank() )
            std::cout << "restart time step and export done\n";
    }
    else
    {
        for( size_type dofThetaPhi : M_dofToUseInThetaPhi )
            M_bdfDensity[dofThetaPhi]->start( *M_ULoc_thetaPhi[dofThetaPhi] /* *M_projSourceTerm*/);
    }


    //M_backend = backend();
    //M_mat = backend()->newMatrix(_test=M_Xh,_trial=M_Xh);
    //M_rhs = backend()->newVector(M_Xh);


    M_mat.resize( M_XhThetaPhi->nDof() );
    M_rhs.resize( M_XhThetaPhi->nDof() );
    if ( M_useOneMatVecByDirection )
    {
        for( size_type dofThetaPhi : M_dofToUseInThetaPhi )
        {
            M_mat[dofThetaPhi] = backend()->newMatrix(_test=M_Xh,_trial=M_Xh);
            M_rhs[dofThetaPhi] = backend()->newVector(M_Xh);
        }

    }
    else
    {
        auto mymat = backend()->newMatrix(_test=M_Xh,_trial=M_Xh);
        auto myrhs = backend()->newVector(M_Xh);
        for( size_type dofThetaPhi : M_dofToUseInThetaPhi )
        {
            M_mat[dofThetaPhi] = mymat;
            M_rhs[dofThetaPhi] = myrhs;
        }
    }

    M_backends.resize( M_XhThetaPhi->nDof() );
    M_preconditioners.resize( M_XhThetaPhi->nDof() );
    if ( M_useOneBackendByDirection )
    {
        for( size_type dofThetaPhi : M_dofToUseInThetaPhi )
        {
            M_backends[dofThetaPhi] = backend_type::build( Environment::vm(), "", Environment::worldComm() );
            M_preconditioners[dofThetaPhi] = preconditioner( _prefix=M_backends[dofThetaPhi]->prefix(),_matrix=M_mat[dofThetaPhi],_pc=M_backends[dofThetaPhi]->pcEnumType(),
                                                             _pcfactormatsolverpackage=M_backends[dofThetaPhi]->matSolverPackageEnumType(), _backend=M_backends[dofThetaPhi] );
        }
    }
    else
    {
        auto mybackend = backend();
        auto myprec = preconditioner( _prefix=mybackend->prefix(),_matrix=M_mat[*M_dofToUseInThetaPhi.begin()],_pc=mybackend->pcEnumType(),
                                      _pcfactormatsolverpackage=mybackend->matSolverPackageEnumType(), _backend=mybackend );
        for( size_type dofThetaPhi : M_dofToUseInThetaPhi )
        {
            M_backends[dofThetaPhi] = mybackend;
            M_preconditioners[dofThetaPhi] = myprec;
        }
    }

    //M_mat = M_backends[*M_dofToUseInThetaPhi.begin()]->newMatrix(_test=M_Xh,_trial=M_Xh);
    //M_rhs = M_backends[*M_dofToUseInThetaPhi.begin()]->newVector(M_Xh);
}

bool
AcousticModel::hasFinishedTimeStep()
{
    return M_bdfDensity[*M_dofToUseInThetaPhi.begin()]->isFinished();
}
void
AcousticModel::nextTimeStep()
{
    for( size_type dofThetaPhi : M_dofToUseInThetaPhi )
        M_bdfDensity[dofThetaPhi]->next( *M_ULoc_thetaPhi[dofThetaPhi] );
}
double
AcousticModel::time()
{
    return M_bdfDensity[*M_dofToUseInThetaPhi.begin()]->time();
}


void
AcousticModel::exportResults( double time )
{
    //M_exporter->step( time )->add( "sound-intensity-sum", *ULocSumAllDirection );
    M_exporter->step( time )->add( "sound-density", *M_ULocIntegrateAllDirection );
    bool doExportSolForEachVecDir = false;
    if ( doExportSolForEachVecDir )// (dof_thetaPhi % 5) == 0 )
        for( size_type dofThetaPhi : M_dofToUseInThetaPhi )
            M_exporter->step(time)->add( (boost::format("sound-density%1%")%dofThetaPhi).str(), *(M_ULoc_thetaPhi[dofThetaPhi]) );

    /*double vSoundVelocity=0.;
    double refPressure=0.;
    double airDensity=0.;*/
    //auto chiPositiveExpr = chi( idv(ULocIntegrateAllDirection) > 1e-8 );
    //auto ensurePositiveExpr = 1e-5*(1-chiPositiveExpr);
    auto soundPressureLevelExpr = 10*vf::log( idv(M_ULocIntegrateAllDirection)*M_airDensity*M_vSoundVelocity/cst(std::pow(M_refPressure,2) ) );
    M_soundPressureLevel->on(_range=elements(M_Xh->mesh()),_expr=soundPressureLevelExpr );
    M_exporter->step( time )->add( "sound-pressure-level", *M_soundPressureLevel );

    M_exporter->save();

}



void
AcousticModel::solve()
{
    size_type cpt_dof_thetaPhi = 0;
    for( size_type dofThetaPhi : M_dofToUseInThetaPhi )
    {
        if ( Environment::isMasterRank() )
        {
            if ( M_verbose )
                std::cout << "dof_thetaPhi : " << ++cpt_dof_thetaPhi << "/" << M_dofToUseInThetaPhi.size() << "\n";
            else
            {
                std::cout << "\r";
                std::cout << "dof_thetaPhi : " << ++cpt_dof_thetaPhi << "/" << M_dofToUseInThetaPhi.size() << std::flush;
            }
        }

        //M_mat[dofThetaPhi]->zero();
        //M_rhs[dofThetaPhi]->zero();

        this->updateAssemblyTransportModel( M_mat[dofThetaPhi],M_rhs[dofThetaPhi],dofThetaPhi, !M_hasBuildCstPart );
        this->updateAssemblyBC( M_mat[dofThetaPhi],M_rhs[dofThetaPhi],dofThetaPhi, !M_hasBuildCstPart );

        boost::mpi::timer mytimer;

        if ( M_useOneMatVecByDirection && !M_useOneBackendByDirection )
        {
            // a particular case
#if 0
            auto myprec = preconditioner( _prefix=M_backends[dofThetaPhi]->prefix(),_matrix=M_mat[dofThetaPhi],_pc=M_backends[dofThetaPhi]->pcEnumType(),
                                          _pcfactormatsolverpackage=M_backends[dofThetaPhi]->matSolverPackageEnumType(), _backend=M_backends[dofThetaPhi],_rebuild=true );
            M_backends[dofThetaPhi]->solve(_matrix=M_mat[dofThetaPhi], _rhs=M_rhs[dofThetaPhi], _solution=*M_ULoc_thetaPhi[dofThetaPhi],_prec=myprec);
#endif
            backend(_rebuild=true)->solve(_matrix=M_mat[dofThetaPhi], _rhs=M_rhs[dofThetaPhi], _solution=*M_ULoc_thetaPhi[dofThetaPhi]);
        }
        else
            M_backends[dofThetaPhi]->solve(_matrix=M_mat[dofThetaPhi], _rhs=M_rhs[dofThetaPhi], _solution=*M_ULoc_thetaPhi[dofThetaPhi],_prec=M_preconditioners[dofThetaPhi] );

        double tElapsed = mytimer.elapsed();
        if ( M_verbose && Environment::isMasterRank() )
            std::cout << " solve done in " << tElapsed << "s\n";

    } // for( size_type dof_thetaPhi : M_dofToUseInThetaPhi )

    // optimisation in this case only
    if ( M_useOneMatVecByDirection )
        M_hasBuildCstPart=true;

    if ( Environment::isMasterRank() && !M_verbose )
        std::cout << "\n";


    // compute field of interest
    M_ULocIntegrateAllDirection->zero();
    auto energyProjOnThetaPhi = M_XhThetaPhi->elementPtr();
    for ( size_type dof=0;dof<M_Xh->nLocalDof();++dof)
    {
        // nodal projection on ThetaPhi space and integrate
        for( size_type dof_thetaPhi=0 ; dof_thetaPhi<M_XhThetaPhi->nDof() ; ++dof_thetaPhi )
        {
            if ( M_dofIsUsed_thetaPhi[dof_thetaPhi] )
                energyProjOnThetaPhi->set(dof_thetaPhi,M_ULoc_thetaPhi[dof_thetaPhi]->operator()( dof ) );
            else
                energyProjOnThetaPhi->set(dof_thetaPhi,M_ULoc_thetaPhi[M_mapPeriodicDofThetaPhi.find(dof_thetaPhi)->second]->operator()( dof ) );
        }
        auto theta = Px();
        double nodalVal = integrate(_range=elements(M_meshThetaPhi),
                                    _expr=idv(energyProjOnThetaPhi)*sin(theta) ).evaluate(false)(0,0);
        M_ULocIntegrateAllDirection->set( dof, nodalVal );
    }

}


void
AcousticModel::updateAssemblyTransportModel( sparse_matrix_ptrtype & mat,vector_ptrtype & rhs, size_type dofThetaPhi, bool buildCstPart )
{
    if ( buildCstPart )
        mat->zero();
    rhs->zero();

    // vector transport direction
    std::vector<double> vDirection(3,0.);
    vDirection[0] = M_vxProj->operator()(dofThetaPhi);
    vDirection[1] = M_vyProj->operator()(dofThetaPhi);
    vDirection[2] = M_vzProj->operator()(dofThetaPhi);
    //if ( Environment::isMasterRank() )
    //    std::cout << "vDirection[0] " << vDirection[0] << " vDirection[1] " << vDirection[1] << " vDirection[2] " << vDirection[2] << "\n";
    auto u = M_ULoc_thetaPhi[dofThetaPhi];
    //-----------------------------------------------------------------//
    boost::mpi::timer mytimer;
   //mytimer.restart();
    // time discretisation and source term
    if ( buildCstPart )
        form2(_test=M_Xh,_trial=M_Xh,_matrix=mat) +=
            integrate(_range=elements(M_mesh),
                      _expr=M_bdfDensity[dofThetaPhi]->polyDerivCoefficient(0)*idt(u)*id(u) );
    auto polyDerivInTime = M_bdfDensity[dofThetaPhi]->polyDeriv();
    auto rhsExpr = idv(M_projSourceTerm)+idv(polyDerivInTime);

    form1(_test=M_Xh,_vector=rhs) +=
        integrate(_range=elements(M_mesh),
                  _expr=rhsExpr*id(u) );

    // propagation term
    auto vDirExpr = vec( cst(vDirection[0]),cst(vDirection[1]),cst(vDirection[2]));
    if ( buildCstPart )
        form2(_test=M_Xh,_trial=M_Xh,_matrix=mat) +=
            integrate(_range=elements(M_mesh),
                      _expr=gradt(u)*vDirExpr*id(u) );

    double tElapsed = mytimer.elapsed();
    if ( M_verbose && Environment::isMasterRank() )
        std::cout << " assemblage base done in " << tElapsed << "s\n";

    // stabilisation
    bool doStab = true;
    if ( doStab )
    {
        mytimer.restart();
        // SUPG stab
        auto coeff = vf::h()/( 2*norm2( vDirExpr ) );
        auto residualStabForm2 = M_bdfDensity[dofThetaPhi]->polyDerivCoefficient(0)*idt(u) + gradt(u)*vDirExpr;
        auto residualStabForm1 = rhsExpr;//idv(polyDerivInTime);
        if ( buildCstPart )
            form2(_test=M_Xh,_trial=M_Xh,_matrix=mat) +=
                integrate(_range=elements(M_mesh),
                          _expr= coeff*(grad(u)*vDirExpr)*residualStabForm2 );
        form1(_test=M_Xh,_vector=rhs) +=
            integrate(_range=elements(M_mesh),
                      _expr=coeff*(grad(u)*vDirExpr)*residualStabForm1 );
        tElapsed = mytimer.elapsed();
        if ( M_verbose && Environment::isMasterRank() )
            std::cout << " assemblage stab done in " << tElapsed << "s\n";
    }


}
void
AcousticModel::updateAssemblyBC( sparse_matrix_ptrtype & mat,vector_ptrtype & rhs, size_type dofThetaPhi, bool buildCstPart )
{
    boost::mpi::timer mytimer;

    M_projBCDirichlet->zero();
    typedef std::vector<boost::reference_wrapper<typename MeshTraits<mesh_type>::face_type const> > cont_range_type;
    boost::shared_ptr<cont_range_type> myelts( new cont_range_type );
    myelts->clear();//resize(0);

    std::vector<double> vDirection(3,0.);
    vDirection[0] = M_vxProj->operator()(dofThetaPhi);
    vDirection[1] = M_vyProj->operator()(dofThetaPhi);
    vDirection[2] = M_vzProj->operator()(dofThetaPhi);

    auto wBCintegrate = M_XhThetaPhi->elementPtr();

    // close rhs in this case
    if ( !buildCstPart )
        rhs->close();

    bool extrapUseBdf=false;

    std::vector<bool> dofdone(M_Xh->nLocalDof(),false);
    for ( auto const& face : boundaryfaces(M_mesh) )
    {
        size_type faceId = face.id();

        auto unitNormal = M_mapUnitNormalByFace[ faceId ];
        double valVdotN = vDirection[0]*unitNormal[0] + vDirection[1]*unitNormal[1] + vDirection[2]*unitNormal[2];
        if ( std::abs(valVdotN) < 1e-6 ) continue;
        if ( !M_doComputeHatThetaPhi[faceId][dofThetaPhi] ) continue;


        if ( M_doComputeHatThetaPhi[faceId][dofThetaPhi] ) //valVdotN < 0 )
            myelts->push_back(boost::cref(face));

        for ( uint16_type fDof = 0 ; fDof<M_Xh->dof()->nLocalDofOnFace(true) ; ++fDof)
        {
            const size_type thedof = M_Xh->dof()->faceLocalToGlobal(face.id(),fDof,0).get<0>();
            if ( dofdone[thedof] ) continue;

            // nodal projection of w(x,t,theta,phi) on x,t -> w(theta,phi;x,t)
            for( size_type dof_thetaPhi2=0 ; dof_thetaPhi2<M_XhThetaPhi->nDof() ; ++dof_thetaPhi2 )
            {
                if ( !extrapUseBdf )
                {
                    //use solution at last time step
                    if ( M_dofIsUsed_thetaPhi[dof_thetaPhi2] )
                        wBCintegrate->set(dof_thetaPhi2,M_ULoc_thetaPhi[dof_thetaPhi2]->operator()( thedof ) );
                    else
                        wBCintegrate->set(dof_thetaPhi2,M_ULoc_thetaPhi[M_mapPeriodicDofThetaPhi.find(dof_thetaPhi2)->second]->operator()( thedof ) );
                }
                else
                {
                    // extrapolation using bdf
                    if ( M_dofIsUsed_thetaPhi[dof_thetaPhi2] )
                        wBCintegrate->set(dof_thetaPhi2,M_bdfDensity[dof_thetaPhi2]->poly()( thedof ) );
                    else
                        wBCintegrate->set(dof_thetaPhi2,M_bdfDensity[M_mapPeriodicDofThetaPhi.find(dof_thetaPhi2)->second]->poly()( thedof ) );
                }
            }


            bool useFirstTermInBC = true;
            bool useSecondTermInBC = true;
            if ( M_doComputeHatThetaPhi[faceId][dofThetaPhi] )//valVdotN < 0 )
            {
                //boundary condition with d=0 : w(r,theta,phi,t)=-alpha*w(r,\hat{theta},\hat{phi},t)
                if ( useFirstTermInBC )
                {
                    double wHatThetaPhi = 0;
                    auto evalAtHatThetaPhi = evaluateFromContext( _context=*M_ctxEvalHatThetaPhi[faceId],
                                                                  _expr=idv(wBCintegrate) );
                    wHatThetaPhi = evalAtHatThetaPhi( M_mapEvalFromCtxThetaPhi[faceId][dofThetaPhi],0 );

                    // update projBCDirichlet at this dof
                    // warning : not minus!
                    M_projBCDirichlet->add( thedof, M_alpha*(1-M_d_prob)*wHatThetaPhi );
                }

                if ( useSecondTermInBC )
                {
                    double bcIntegrate = inner_product(*M_mapOperatorDiffuseBC.find(faceId)->second,wBCintegrate->container());
#if !defined(NDEBUG)
                    auto thetaPrime = Px();
                    auto phiPrime = Py();
                    auto vPrime = vSoundVelocity*vec( sin(thetaPrime)*cos(phiPrime), sin(thetaPrime)*sin(phiPrime), cos(thetaPrime) );
                    auto vPrimeDotN = vPrime(0,0)*unitNormal[0] + vPrime(1,0)*unitNormal[1] + vPrime(2,0)*unitNormal[2];
                    auto chiVPrimeDotN = chi( vPrimeDotN > cst(-1e-6) );
                    auto bcIntegrateExpr = (alpha*d_prob/(M_PI*vSoundVelocity))*vPrimeDotN*idv(wBCintegrate)*sin(thetaPrime);

                    double bcIntegrate2 = integrate(_range=elements(M_meshThetaPhi),
                                                    _expr=bcIntegrateExpr*chiVPrimeDotN ).evaluate(false)(0,0);
                    CHECK( std::abs( bcIntegrate - bcIntegrate2 ) < 1e-8 ) << " error between opt " << bcIntegrate <<" and eval : " << bcIntegrate2 << "\n";
#endif
                    M_projBCDirichlet->add( thedof, bcIntegrate );
                }


                if ( !buildCstPart )
                    rhs->set( thedof,M_projBCDirichlet->operator()(thedof) );

            } // if M_doComputeHatThetaPhi


            dofdone[thedof] = true;

        } // for ( uint16_type fDof = 0 ; fDof<Xh->dof()->nLocalDofOnFace(true) ; ++fDof)

    } // for ( auto const& face : boundaryfaces(mesh))

    if ( buildCstPart )
    {
        auto myMarkedFacesBCRange = boost::make_tuple( mpl::size_t<MESH_FACES>(),
                                                       myelts->begin(),
                                                       myelts->end(),
                                                       myelts );

        // up mat and rhs for bc Dirichlet condition
        form2(_test=M_Xh,_trial=M_Xh,_matrix=mat) +=
            on(_range=myMarkedFacesBCRange,
               _rhs=rhs, _element=*M_ULoc_thetaPhi[dofThetaPhi], _expr=idv(M_projBCDirichlet) );

    }

    double tElapsed = mytimer.elapsed();
    if ( M_verbose && Environment::isMasterRank() )
        std::cout << " assemblage bc done in " << tElapsed << "s\n";


}














Feel::po::options_description
acousticModel_options()
{
	po::options_description accousticoptions( "Accoustic options" );
	accousticoptions.add_options()

		( "coeff.M", po::value<double>()->default_value( 0.01 ), "atmospheric attenuation coefficient" )
        ( "coeff.alpha", po::value<double>()->default_value( 0.5 ), "reflexion coefficient" )
        ( "coeff.d-prob", po::value<double>()->default_value( 0.5 ), "probability that reflexion is non specular" )
        ( "sound-velocity", po::value<double>()->default_value( 343 ), " sound velocity [m/s]" )
        ( "air-density", po::value<double>()->default_value( 1.225 ), "air-density [kg/m^3]" )
        ( "reference-pressure", po::value<double>()->default_value( 2e-5 ), "reference-pressure [Pa]" )

		( "stab", po::value<bool>()->default_value( true ), "use SUPG stab with fem solver" )
		//( "stab-rho", po::value<double>()->default_value( 0.25 ), "coeff" )
        ( "gmsh.filename-thetaPhi", po::value<std::string>(), "name for theta phi mesh" )

		( "use-first-term-bc", po::value<bool>()->default_value( true ), "use specular bc" )
		( "use-second-term-bc", po::value<bool>()->default_value( true ), "use non specular bc" )

		( "use-source-as-initial-solution", po::value<bool>()->default_value( true ), "use-source-as-initial-solution" )

		( "scaling-coeff", po::value<double>()->default_value( 5. ), "coeff" )
        ( "center_x", po::value<double>()->default_value( 1 ), "x-position of source" )
        ( "center_y", po::value<double>()->default_value( 0.5 ), "y-position of source" )
        ( "center_z", po::value<double>()->default_value( 0.5 ), "z-position of source" )
        ( "radius", po::value<double>()->default_value( 0.1 ), "radius of source" )

		( "do-export-foreach-sol", po::value<bool>()->default_value( false ), "coeff" )

        ( "extrapolation.use-bdf", po::value<bool>()->default_value( false ), "coeff" )


        ( "acoustic-model.verbose", po::value<bool>()->default_value( false ), "print some additional info" )

        ( "acoustic-model.export-propagation-vectors-in-geofile", po::value<bool>()->default_value( false ), "export propagation vectors in geofile" )
        ( "acoustic-model.export-specular-vectors-in-geofile", po::value<bool>()->default_value( false ), "export speculars vectors in geofile" )
        ( "specular-meshexport-nFaces", po::value<int>()->default_value( 1 ), "coeff" )

        ( "use-storage-from-normal-boundary-faces", po::value<bool>()->default_value( true ), "use-storage-from-normal-boundary-faces" )

        ( "use-one-backend-by-direction", po::value<bool>()->default_value( false ), "use-one-backend-by-direction" )
        ( "use-one-matvec-by-direction", po::value<bool>()->default_value( false ), "use-one-matvec-by-direction" )

		;

    return accousticoptions;
}

} // namespace Feel
