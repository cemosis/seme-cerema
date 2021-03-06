/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:set fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <feel/feelcore/environment.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/projector.hpp>

using namespace Feel;



int main(int argc, char**argv )
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
        //( "nThetaPhi", po::value<int>()->default_value( 10000 ), "max number of direction" )

		( "use-first-term-bc", po::value<bool>()->default_value( true ), "use specular bc" )
		( "use-second-term-bc", po::value<bool>()->default_value( true ), "use non specular bc" )

		( "use-source-as-initial-solution", po::value<bool>()->default_value( true ), "use-source-as-initial-solution" )

		( "scaling-coeff", po::value<double>()->default_value( 5. ), "coeff" )
        ( "center_x", po::value<double>()->default_value( 1 ), "x-position of source" )
        ( "center_y", po::value<double>()->default_value( 0.5 ), "y-position of source" )
        ( "center_z", po::value<double>()->default_value( 0.5 ), "z-position of source" )
        ( "radius", po::value<double>()->default_value( 0.1 ), "radius of source" )

		( "do-export-foreach-sol", po::value<bool>()->default_value( false ), "coeff" )
        //( "coupling-direction", po::value<bool>()->default_value( false ), "coeff" )

        ( "extrapolation.use-bdf", po::value<bool>()->default_value( false ), "coeff" )

        //( "dof_thetaPhi" , po::value<int>()->default_value( 0 ), "coeff" )

        ( "acoustic-model.verbose", po::value<bool>()->default_value( false ), "print some additional info" )

        ( "acoustic-model.export-propagation-vectors-in-geofile", po::value<bool>()->default_value( false ), "export propagation vectors in geofile" )
        ( "acoustic-model.export-specular-vectors-in-geofile", po::value<bool>()->default_value( false ), "export speculars vectors in geofile" )
        ( "specular-meshexport-nFaces", po::value<int>()->default_value( 1 ), "coeff" )

        ( "use-storage-from-normal-boundary-faces", po::value<bool>()->default_value( true ), "use-storage-from-normal-boundary-faces" )
		;

	Environment env( _argc=argc, _argv=argv,
                     _desc=accousticoptions,
                     _desc_lib=feel_options().add( backend_options( "l2proj") ),
                     _about=about(_name="seme-acoustic",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    static const uint16_type nDim = 3;

    typedef Mesh<Hypercube<1,1,nDim>> mesh_1d_type;
    typedef boost::shared_ptr<mesh_1d_type> mesh_1d_ptrtype;
    typedef FunctionSpace<mesh_1d_type, bases<Lagrange<1, Scalar> > > space_1d_type;
    typedef boost::shared_ptr<space_1d_type> space_1d_ptrtype;
    typedef space_1d_type::element_type element_1d_type;
    typedef boost::shared_ptr<element_1d_type> element_1d_ptrtype;

    //typedef Mesh<Hypercube<nDim>> mesh_type;
    typedef Mesh<Simplex<nDim>> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef FunctionSpace<mesh_type, bases<Lagrange<1, Scalar> > > space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef space_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;

    typedef Bdf<space_type> bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrtype;

    typedef Mesh<Hypercube<2>> mesh_2d_type;
    typedef boost::shared_ptr<mesh_2d_type> mesh_2d_ptrtype;
    typedef FunctionSpace<mesh_2d_type, bases<Lagrange<1, Scalar> > > space_2d_type;
    typedef boost::shared_ptr<space_2d_type> space_2d_ptrtype;
    typedef space_2d_type::element_type element_2d_type;
    typedef boost::shared_ptr<element_2d_type> element_2d_ptrtype;

    typedef Backend<double> backend_type;
    //typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef backend_type::vector_ptrtype vector_ptrtype;


    double doStab = boption(_name="stab");
    //double rho = doption(_name="stab-rho");
    double M = doption(_name="coeff.M");
    double alpha = doption(_name="coeff.alpha");
    double d_prob = doption(_name="coeff.d-prob");//0.5;
    double scalCoef = doption(_name="scaling-coeff");
    double vSoundVelocity = doption(_name="sound-velocity");//1;
    bool useFirstTermInBC = boption(_name="use-first-term-bc");
    bool useSecondTermInBC = boption(_name="use-second-term-bc");
    bool doExportSolForEachVecDir = boption(_name="do-export-foreach-sol");
    //bool couplingDirection = boption(_name="coupling-direction");//ASUP
    bool extrapUseBdf = boption(_name="extrapolation.use-bdf");

    double airDensity = doption(_name="air-density");
    double refPressure = doption(_name="reference-pressure");

    bool exportPropVecInGeoFile = boption(_prefix="acoustic-model",_name="export-propagation-vectors-in-geofile");
    bool exportSpecVecInGeoFile = boption(_prefix="acoustic-model",_name="export-specular-vectors-in-geofile");

    //-----------------------------------------------------------------//
    // create FE mesh/space and 1d meshes/spaces for bases dim
    mesh_ptrtype mesh;
    space_ptrtype Xh;
    element_ptrtype ULoc,ULocSumAllDirection,ULocIntegrateAllDirection,projBCDirichlet,projSourceTerm;
    element_ptrtype soundPressureLevel;


    //-----------------------------------------------------------------//
    // build mesh for theta phi
    mesh_2d_ptrtype mesh_thetaPhi;
    // in // start by master rank (write msh file), and load after others process
    std::string fileNameMeshDescThetaPhi = Environment::expand( soption(_name="gmsh.filename-thetaPhi") );
    if ( Environment::isMasterRank() )
        mesh_thetaPhi = loadMesh(_mesh=new mesh_2d_type, _filename=fileNameMeshDescThetaPhi,
                                 _worldcomm=Environment::worldCommSeq());
    Environment::worldComm().barrier();
    if ( !Environment::isMasterRank() )
    {
        std::string fileNameMshThetaPhi = fs::path( fileNameMeshDescThetaPhi ).stem().string()+".msh";
        //std::cout << "filename to load " << fileNameMshThetaPhi << "\n";
        mesh_thetaPhi = loadMesh(_mesh=new mesh_2d_type, _filename=fileNameMshThetaPhi,
                                 _worldcomm=Environment::worldCommSeq());
    }

    // build function space for theta phi
    space_2d_ptrtype Xh_thetaPhi;
    Xh_thetaPhi = space_2d_type::New(_mesh=mesh_thetaPhi,
                                     _worldscomm=std::vector<WorldComm>(1,Environment::worldCommSeq()) );

    // find doublon in thetaphi functionspace and compute dof's periodicity
    std::set<size_type> dofToUse_thetaPhi;
    std::vector<bool> dofIsUsed_thetaPhi(Xh_thetaPhi->nDof(),true);
    std::map<size_type,size_type> mapPeriodicDof_thetaPhi;

    bool findThetaPhiLeft=false,findThetaPhiRight=false;
    size_type idThetaPhiLeft = invalid_size_type_value, idThetaPhiRight = invalid_size_type_value;
    for( size_type dof_thetaPhi=0 ; dof_thetaPhi<Xh_thetaPhi->nDof() ; ++dof_thetaPhi )
    {
        auto const& mythetaPhiPt = Xh_thetaPhi->dof()->dofPoint(dof_thetaPhi).get<0>();

        if ( mythetaPhiPt[1] < (2*M_PI-0.00001) )
        {

            if ( std::abs( mythetaPhiPt[0] ) < 1e-9 )
            {
                if ( !findThetaPhiLeft )
                {
                    dofToUse_thetaPhi.insert(dof_thetaPhi);
                    findThetaPhiLeft=true;
                    idThetaPhiLeft = dof_thetaPhi;
                }
                else
                {
                    dofIsUsed_thetaPhi[dof_thetaPhi] = false;
                    mapPeriodicDof_thetaPhi[ dof_thetaPhi ] = idThetaPhiLeft;
                }
            }
            else if ( std::abs( mythetaPhiPt[0] - M_PI ) < 1e-9 )
            {
                if ( !findThetaPhiRight )
                {
                    dofToUse_thetaPhi.insert(dof_thetaPhi);
                    findThetaPhiRight=true;
                    idThetaPhiRight = dof_thetaPhi;
                }
                else
                {
                    dofIsUsed_thetaPhi[dof_thetaPhi] = false;
                    mapPeriodicDof_thetaPhi[ dof_thetaPhi ] = idThetaPhiRight;
                }
            }
            else
            {
                dofToUse_thetaPhi.insert(dof_thetaPhi);
            }
        }
        else
        {
            dofIsUsed_thetaPhi[dof_thetaPhi] = false;
            // search corresponding dof
            bool hasFindPeriodicDof=false;
            for( size_type dof_thetaPhi2=0 ; dof_thetaPhi2<Xh_thetaPhi->nDof() && !hasFindPeriodicDof ; ++dof_thetaPhi2 )
            {
                auto const& mythetaPhiPt2 = Xh_thetaPhi->dof()->dofPoint(dof_thetaPhi2).get<0>();
                if ( std::abs( mythetaPhiPt[0]-mythetaPhiPt2[0] ) < 1e-9 && std::abs( mythetaPhiPt2[1] ) < 1e-9 )
                {
                    mapPeriodicDof_thetaPhi[ dof_thetaPhi ] = dof_thetaPhi2;
                    hasFindPeriodicDof=true;
                }
            }
        }

    }

    if ( Environment::isMasterRank() )
        std::cout << "Xh_thetaPhi->nDof : = " << Xh_thetaPhi->nDof() << std::endl;

    //-----------------------------------------------------------------//
    // build mesh
    mesh = loadMesh(_mesh=new mesh_type);

    //-----------------------------------------------------------------//
    // build space
    Xh = space_type::New(_mesh=mesh);
    ULoc = Xh->elementPtr(); // ASUP
    ULocSumAllDirection = Xh->elementPtr();
    ULocIntegrateAllDirection = Xh->elementPtr();
    projBCDirichlet = Xh->elementPtr();
    projSourceTerm = Xh->elementPtr();
    soundPressureLevel = Xh->elementPtr();

    std::vector<bdf_ptrtype> myBdf( Xh_thetaPhi->nDof() );
    for( size_type dof_thetaPhi : dofToUse_thetaPhi )
        myBdf[dof_thetaPhi] = bdf( _vm=Environment::vm(), _space=Xh,
                                   _name=(boost::format("energy%1%-%2%_%3%")%dof_thetaPhi %Environment::worldComm().rank() %Environment::worldComm().size() ).str());

    if ( Environment::isMasterRank() )
        std::cout << "Xh->nDof : = " << Xh->nDof() << std::endl;

    //-----------------------------------------------------------------//
    // source expression
    double center_x=doption(_name="center_x"),center_y=doption(_name="center_y"),center_z=doption(_name="center_z");
    double radius=doption(_name="radius");
    bool useSourceAsInitialSolution = boption("use-source-as-initial-solution");


    auto chiInitSol1d = chi( pow(Px()-center_x,2) <= cst(std::pow(radius,2)) );
    auto chiInitSol2d = chi( pow(Px()-center_x,2) + pow(Py()-center_y,2) <= cst(std::pow(radius,2)) );
    auto chiInitSol3d = chi( pow(Px()-center_x,2) + pow(Py()-center_y,2) + pow(Pz()-center_z,2) <= cst(std::pow(radius,2)) );
    auto sourceTerm1dExpr = -scalCoef*
        (Px()-(center_x-radius))*(Px()-(center_x+radius))*chiInitSol1d;
    auto sourceTerm2dExpr = -scalCoef*
        (Px()-(center_x-radius))*(Px()-(center_x+radius))*
        (Py()-(center_y-radius))*(Py()-(center_y+radius))*chiInitSol2d;

    auto sourceTerm3dExpr = scalCoef*chiInitSol3d;
    /*auto sourceTerm3dExpr = -scalCoef*
        (Px()-(center_x-radius))*(Px()-(center_x+radius))*
        (Py()-(center_y-radius))*(Py()-(center_y+radius))*
        (Pz()-(center_z-radius))*(Pz()-(center_z+radius))*chiInitSol3d;*/
#if 0
    if ( nDim == 1 )
        ULoc->on( _range=elements(mesh), _expr=-scalCoef*Px()*(Px()-cst(1.))*chi(Px()<=1) );
    else if ( nDim == 2 )
        ULoc->on( _range=elements(mesh), _expr=-scalCoef*Px()*(Px()-cst(1.))*chi(Px()<=1)*Py()*(1.-Py()) );
    else if ( nDim == 3 )
        ULoc->on( _range=elements(mesh), _expr=sourceTerm3dExpr );
#endif

#if 0
    // L2 projection of source expression
    auto l2proj = opProjection(_domainSpace=Xh,_imageSpace=Xh,
                               _backend=backend_type::build( Environment::vm(), "l2proj", Xh->worldComm() ),
                               _type=::Feel::L2 );
    if ( nDim == 1 )
        *projSourceTerm = l2proj->operator()( sourceTerm1dExpr );
    else if ( nDim == 2 )
        *projSourceTerm = l2proj->operator()( sourceTerm2dExpr );
    else if ( nDim == 3 )
        *projSourceTerm = l2proj->operator()( sourceTerm3dExpr );
#else
    *projSourceTerm = vf::project( _space=Xh,_expr=sourceTerm3dExpr );
#endif
    // take initial solution as source term
    if ( useSourceAsInitialSolution )
    {
        *ULoc=*projSourceTerm;
        // put zero in order to have no source in simulation
        projSourceTerm->zero();
    }

    std::vector<element_ptrtype> ULoc_thetaPhi(Xh_thetaPhi->nDof() );
    for( size_type dof_thetaPhi : dofToUse_thetaPhi )
    {
        ULoc_thetaPhi[dof_thetaPhi] = Xh->elementPtr();
        *ULoc_thetaPhi[dof_thetaPhi] = *ULoc;
    }
    //for ( size_type dof=0;dof<Xh->nLocalDof();++dof)
    //ULoc_thetaPhi[dof_thetaPhi]->set(dof, ULoc->operator()(dof) );

    //-----------------------------------------------------------------//
    //int nThetaPhiUsed = ioption(_name="nThetaPhi");
    //size_type nDof_thetaPhi = std::min( size_type(nThetaPhiUsed),dofToUse_thetaPhi.size() );//ASUP

    // build exporter
    auto e = exporter( _mesh=mesh,_name="myexporter" );

    //-----------------------------------------------------------------//

    if ( !myBdf[0]->isRestart() )
    {
        for( size_type dof_thetaPhi : dofToUse_thetaPhi )
            myBdf[dof_thetaPhi]->start(*ULoc);

        // export initial solution
        e->step(0)->add( "sound-intensity-sum", *ULoc );
        e->step(0)->add( "sound-intensity-integrate", *ULoc );

        if ( doExportSolForEachVecDir )
            for( size_type dof_thetaPhi : dofToUse_thetaPhi )
                e->step(0)->add( (boost::format("NEWULoc%1%")%dof_thetaPhi).str(), *ULoc );

        e->step(0)->add( "sound-pressure-level", *soundPressureLevel );
        e->save();

        if ( Environment::isMasterRank() )
            std::cout << "init time step and export initial solution done\n";
    }
    else
    {
        for( size_type dof_thetaPhi : dofToUse_thetaPhi )
        {
            myBdf[dof_thetaPhi]->restart();
            // load a previous solution as current solution
            *ULoc_thetaPhi[dof_thetaPhi] = myBdf[dof_thetaPhi]->unknown(0);
        }
        // restart exporter
        if ( e->doExport() )
            e->restart( myBdf[0]->timeInitial() );

        if ( Environment::isMasterRank() )
            std::cout << "restart time step and export done\n";
    }


    //-----------------------------------------------------------------//


    //-----------------------------------------------------------------//
    // create algebraic structure
    auto rhs = backend()->newVector( Xh );
    auto USol = backend()->newVector( Xh );
    auto mat = backend()->newMatrix(_test=Xh,_trial=Xh);

    //std::vector<boost::shared_ptr<Backend<double> > > mybackends(nDof_thetaPhi,backend(_rebuild=true));
    auto mybackend = backend();
    auto myprec = preconditioner( _prefix=mybackend->prefix(),_matrix=mat,_pc=mybackend->pcEnumType(),
                                  _pcfactormatsolverpackage=mybackend->matSolverPackageEnumType(), _backend=mybackend );

    // copy initial solution to vector solution
    *USol = *ULoc;
    //-----------------------------------------------------------------//
    //vitesse de propagation
    std::vector<double> vDirection(nDim,0.);

    element_2d_ptrtype vx_proj, vy_proj, vz_proj;
    vx_proj=Xh_thetaPhi->elementPtr();
    vy_proj=Xh_thetaPhi->elementPtr();
    vz_proj=Xh_thetaPhi->elementPtr();
    vx_proj->on( _range=elements(mesh_thetaPhi), _expr=vSoundVelocity*sin(Px())*cos(Py()) );
    vy_proj->on( _range=elements(mesh_thetaPhi), _expr=vSoundVelocity*sin(Px())*sin(Py()) );
    vz_proj->on( _range=elements(mesh_thetaPhi), _expr=vSoundVelocity*cos(Px()) );

    // export each vector direction in geo file for Gmsh
    if ( exportPropVecInGeoFile && !myBdf[0]->isRestart() && Environment::isMasterRank() )
    {
        std::ofstream fileAllDir("allvectordirection.geo", std::ios::out);
        fileAllDir << "myh=0.5;\n";
        fileAllDir << "Point(1) = {0, 0, 0, myh};\n";
        for( size_type dof_thetaPhi : dofToUse_thetaPhi )
            fileAllDir << "Point(" << dof_thetaPhi+2 <<") = {"
                       << vx_proj->operator()(dof_thetaPhi) << ","
                       << vy_proj->operator()(dof_thetaPhi) << ","
                       << vz_proj->operator()(dof_thetaPhi) << ", myh};\n";

        for( size_type dof_thetaPhi : dofToUse_thetaPhi )
            fileAllDir << "Line(" << dof_thetaPhi+1 <<") = {1,"<< dof_thetaPhi+2 << "};\n";
        fileAllDir.close();
    }

    //-----------------------------------------------------------------//

    auto pcf =  mesh->gm()->preComputeOnFaces( mesh->gm(), mesh->gm()->referenceConvex().barycenterFaces() );
    CHECK( mesh->beginFace()!=mesh->endFace() ) << "no face on mesh for this proc : must be take into account!\n";
    auto firstFacePtr = mesh->beginFace();
    auto ctx = mesh->gm()->context<vm::POINT|vm::NORMAL|vm::KB|vm::JACOBIAN>( firstFacePtr->element0(),
                                                                              pcf,
                                                                              firstFacePtr->pos_first() );
    auto wBCintegrate = Xh_thetaPhi->elementPtr();

    //-----------------------------------------------------------------//
    //-----------------------------------------------------------------//
    //-----------------------------------------------------------------//
    //-----------------------------------------------------------------//
    // localisation \hat{theta},\hat{phi}
    typedef space_2d_type::Context context_type;
    typedef boost::shared_ptr<context_type> context_ptrtype;

    auto locMeshThetaPhi = mesh_thetaPhi->tool_localization();
    locMeshThetaPhi->updateForUse();
    std::vector<size_type> localisationHatThetaPhi( Xh_thetaPhi->nDof(), invalid_size_type_value );

    //-----------------------------------------------------------------//
    //-----------------------------------------------------------------//
    // NEW
    std::map<size_type,std::map<size_type,std::pair<node_type,node_type> > > mapThetaPhiAndHatThetaPhiNEW;
    std::map<size_type,context_ptrtype> ctxEvalHatThetaPhiNEW;//( new context_type( Xh_thetaPhi ) );
    std::map<size_type,std::map<size_type,size_type> > mapEvalFromCtxThetaPhiNEW;
    std::map<size_type,std::map<size_type,bool> > doComputeHatThetaPhi;
    std::map<size_type, node_type> mapUnitNormalByFace;
    std::map<size_type,vector_ptrtype > mapOperatorDiffuseBC;
    std::set< boost::tuple<double,double,double,size_type> > setOfNormalVectorStored;

    size_type nBoundaryFaces = nelements(boundaryfaces(mesh),false);
    size_type cptBoundaryFaces = 0;
    bool useStorageNormalBoundaryFaces= boption(_name="use-storage-from-normal-boundary-faces"); // )true;
    for ( auto const& face : boundaryfaces(mesh) )
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
                ctxEvalHatThetaPhiNEW[faceId] = ctxEvalHatThetaPhiNEW.find(faceIdReference)->second;
                mapUnitNormalByFace[ faceId ] = unitNormal;
                doComputeHatThetaPhi[faceId] = doComputeHatThetaPhi.find(faceIdReference)->second;
                mapThetaPhiAndHatThetaPhiNEW[ faceId ] = mapThetaPhiAndHatThetaPhiNEW.find(faceIdReference)->second;
                mapEvalFromCtxThetaPhiNEW[faceId] = mapEvalFromCtxThetaPhiNEW.find(faceIdReference)->second;
                mapOperatorDiffuseBC[faceId] = mapOperatorDiffuseBC.find(faceIdReference)->second;
                // next face
                continue;
            }

            // say that data for this normal vector will be computed
            setOfNormalVectorStored.insert( boost::make_tuple( unitNormal[0],unitNormal[1],unitNormal[2],faceId ) );
        }


        ctxEvalHatThetaPhiNEW[faceId].reset( new context_type( Xh_thetaPhi ) );
        mapUnitNormalByFace[ faceId ] = unitNormal;

        size_type cptNodeCtx = 0;
        for( size_type dof_thetaPhi : dofToUse_thetaPhi )
        {
            vDirection[0] = vx_proj->operator()(dof_thetaPhi);
            if( nDim >= 2 )
                vDirection[1] = vy_proj->operator()(dof_thetaPhi);
            if( nDim == 3 )
                vDirection[2] = vz_proj->operator()(dof_thetaPhi);

            double valVdotN = vDirection[0]*unitNormal[0] + vDirection[1]*unitNormal[1] + vDirection[2]*unitNormal[2];

            // on prend les vecteur speculaire partant de la paroi (on calcul les vecteurs arrivants sur la paroi )
            if ( valVdotN > -1e-6 )
            {
                doComputeHatThetaPhi[faceId][dof_thetaPhi]=false;
                continue;
            }

            doComputeHatThetaPhi[faceId][dof_thetaPhi]=true;

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

            auto thetaPhiPt = Xh_thetaPhi->dof()->dofPoint(dof_thetaPhi).get<0>();
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

            ctxEvalHatThetaPhiNEW[faceId]->add( hatThetaPhi );
            mapEvalFromCtxThetaPhiNEW[faceId][dof_thetaPhi]=cptNodeCtx;
            ++cptNodeCtx;

        } // for (dof_thetaPhi)

        if ( true )
        {
            mapOperatorDiffuseBC[faceId] = backend()->newVector(Xh_thetaPhi);

            auto thetaPrime = Px();
            auto phiPrime = Py();
            auto vPrime = vSoundVelocity*vec( sin(thetaPrime)*cos(phiPrime), sin(thetaPrime)*sin(phiPrime), cos(thetaPrime) );
            auto vPrimeDotN = vPrime(0,0)*unitNormal[0] + vPrime(1,0)*unitNormal[1] + vPrime(2,0)*unitNormal[2];
#if 0
            auto chiVPrimeDotN = chi( vPrimeDotN > cst(-1e-6) );// not compile with linearform(maybe bilinear also
#else
            auto chiVPrimeDotN = chi( ( sin(thetaPrime)*cos(phiPrime)*unitNormal[0] + sin(thetaPrime)*sin(phiPrime)*unitNormal[1] + cos(thetaPrime)*unitNormal[2] )  > cst(-1e-6) );
#endif
            auto uThetaPhi = Xh_thetaPhi->element();
            auto bcIntegrateExpr = (alpha*d_prob/(M_PI*vSoundVelocity))*vPrimeDotN*id(uThetaPhi)*sin(thetaPrime);
            form1(_test=Xh_thetaPhi,_vector=mapOperatorDiffuseBC[faceId] )
                = integrate(_range=elements(mesh_thetaPhi),
                            _expr=bcIntegrateExpr*chiVPrimeDotN
                            );
            mapOperatorDiffuseBC[faceId]->close();
        }


    } // faces

    if ( Environment::isMasterRank() )
        std::cout << "\nwaiting other process (can be take a long time)" << std::flush;
    Environment::worldComm().barrier();
    if ( Environment::isMasterRank() )
        std::cout << "\rfinish precompute (size of normalVectorStored.size : " << setOfNormalVectorStored.size() << ")\n";





    //-----------------------------------------------------------------//
    //-----------------------------------------------------------------//
    //-----------------------------------------------------------------//
    //-----------------------------------------------------------------//
    //-----------------------------------------------------------------//

    // export each vector speculaire direction in geo file for Gmsh
    if ( !myBdf[0]->isRestart() && Environment::isMasterRank() )
    {
        std::ofstream fileAllDir("speculaire.geo", std::ios::out);

        int cptGmshPt=1,cptGmshLine=1;

        std::vector<int> listMarkerLinesAllThetaPhi,listMarkerLinesAllHatThetaPhi,listMarkerBoundaryFaces,listMarkerNormal;
        std::map<size_type,/*std::list<int>*/std::vector<int> > listMarkerLinesThetaPhiFixed, listMarkerLinesHatThetaPhiFixed;
        fileAllDir << "myh=0.5;\n";

        int nFaceExported = ioption(_name="specular-meshexport-nFaces");
        size_type cptFaces = 0;
        for ( auto const& face : boundaryfaces(mesh) )
        {
            ++cptFaces;
            if ( nFaceExported > 0 && cptFaces > nFaceExported )
                break;

            size_type faceId = face.id();
            uint16_type nPt = face.nPoints();
            auto bary = ublas::column( face.G(), 0 );

            ctx->update( face.element0(), face.pos_first() );
            auto unitNormal = ctx->unitNormal();

            int cptPtStart=cptGmshPt;
            for ( uint16_type k = 0;k< nPt; ++k,++cptGmshPt)
            {
                auto pt = ublas::column( face.G(), k );
                fileAllDir << "Point(" << cptGmshPt <<") = {"
                           << pt[0] << "," << pt[1] << "," << pt[2] << ", myh};\n";
            }
            //for ( uint16_type k = 0;k< nPt; ++k,++cptGmshPt)
            //    fileAllDir << "Line(" << dof_thetaPhi+1 <<") = {1,"<< dof_thetaPhi+2 << "};\n";

            for ( uint16_type k = 0;k< nPt-1; ++k,++cptGmshPt,++cptGmshLine)
            {
                fileAllDir << "Line(" << cptGmshLine <<") = {"<< cptPtStart+k << ","<< cptPtStart+k+1 << "};\n";
                listMarkerBoundaryFaces.push_back( cptGmshLine );
            }
            //++cptGmshLine;
            fileAllDir << "Line(" << ++cptGmshLine <<") = {"<< cptPtStart+nPt-1 << ","<< cptPtStart << "};\n";
            listMarkerBoundaryFaces.push_back( cptGmshLine );


            fileAllDir << "Point(" << ++cptGmshPt <<") = {"
                       << bary[0] << "," << bary[1] << "," << bary[2] << ", myh};\n";
            int saveIdBary = cptGmshPt;

            fileAllDir << "Point(" << ++cptGmshPt <<") = {"
                       << bary[0]+unitNormal[0] << ","
                       << bary[1]+unitNormal[1] << ","
                       << bary[2]+unitNormal[2] << ", myh};\n";
            fileAllDir << "Line(" << ++cptGmshLine <<") = {"<< saveIdBary << ","<< cptGmshPt << "};\n";
            listMarkerNormal.push_back( cptGmshLine );



            for( size_type dof_thetaPhi : dofToUse_thetaPhi )
            {
                if ( !doComputeHatThetaPhi[faceId][dof_thetaPhi] ) continue;
                //if ( dof_thetaPhi != ioption(_name="dof_thetaPhi") ) continue;
                //CHECK( doComputeHatThetaPhi[faceId][dof_thetaPhi] ) << "HatThetaPhi not computed";
#if 0
                double thetheta = mapThetaPhiAndHatThetaPhi[dof_thetaPhi].first[0];
                double thephi = mapThetaPhiAndHatThetaPhi[dof_thetaPhi].first[1];
                double thehattheta = mapThetaPhiAndHatThetaPhi[dof_thetaPhi].second[0];
                double thehatphi = mapThetaPhiAndHatThetaPhi[dof_thetaPhi].second[1];
#else
                double thetheta = mapThetaPhiAndHatThetaPhiNEW[ faceId ][dof_thetaPhi].first[0];
                double thephi = mapThetaPhiAndHatThetaPhiNEW[ faceId ][dof_thetaPhi].first[1];
                double thehattheta = mapThetaPhiAndHatThetaPhiNEW[ faceId ][dof_thetaPhi].second[0];
                double thehatphi = mapThetaPhiAndHatThetaPhiNEW[ faceId ][dof_thetaPhi].second[1];
#endif

                double theVx = math::sin(thetheta)*math::cos(thephi)+bary[0];
                double theVy = math::sin(thetheta)*math::sin(thephi)+bary[1];
                double theVz = math::cos(thetheta)+bary[2];
                double theHatVx = math::sin(thehattheta)*math::cos(thehatphi)+bary[0];
                double theHatVy = math::sin(thehattheta)*math::sin(thehatphi)+bary[1];
                double theHatVz = math::cos(thehattheta)+bary[2];
                //fileAllDir << "Point(" << ++cptGmshPt <<") = {"
                //           << bary[0] << "," << bary[1] << "," << bary[2] << ", myh};\n";
                fileAllDir << "Point(" << ++cptGmshPt <<") = {"
                           << theVx << "," << theVy << "," << theVz << ", myh};\n";
                fileAllDir << "Point(" << ++cptGmshPt <<") = {"
                           << theHatVx << "," << theHatVy << "," << theHatVz << ", myh};\n";

                fileAllDir << "Line(" << ++cptGmshLine <<") = {"<< saveIdBary << ","<< cptGmshPt-1 << "};\n";
                listMarkerLinesAllThetaPhi.push_back(cptGmshLine);
                listMarkerLinesThetaPhiFixed[dof_thetaPhi].push_back(cptGmshLine);

                fileAllDir << "Line(" << ++cptGmshLine <<") = {"<< saveIdBary << ","<< cptGmshPt << "};\n";
                listMarkerLinesAllHatThetaPhi.push_back(cptGmshLine);
                listMarkerLinesHatThetaPhiFixed[dof_thetaPhi].push_back(cptGmshLine);

                ++cptGmshPt, ++cptGmshLine;
            }

        } // faces

        int nMarker = listMarkerLinesAllThetaPhi.size();
        fileAllDir << "Physical Line(\"ThetaPhi\") = {";
        if (nMarker > 1 )
            for ( int k = 0;k<(nMarker-1);++k )
                fileAllDir << listMarkerLinesAllThetaPhi[k] << ",";
        fileAllDir << listMarkerLinesAllThetaPhi[nMarker-1] << "};\n";
        fileAllDir << "Physical Line(\"HatThetaPhi\") = {";
        if (nMarker > 1 )
            for ( int k = 0;k<(nMarker-1);++k )
                fileAllDir << listMarkerLinesAllHatThetaPhi[k] << ",";
        fileAllDir << listMarkerLinesAllHatThetaPhi[nMarker-1] << "};\n";

        int nMarkerBoundaryFaces = listMarkerBoundaryFaces.size();
        fileAllDir << "Physical Line(\"BoundaryFaces\") = {";
        if (nMarker > 1 )
            for ( int k = 0;k<(nMarkerBoundaryFaces-1);++k )
                fileAllDir << listMarkerBoundaryFaces[k] << ",";
        fileAllDir << listMarkerBoundaryFaces[nMarkerBoundaryFaces-1] << "};\n";

        int nMarkerNormal = listMarkerNormal.size();
        fileAllDir << "Physical Line(\"Normal\") = {";
        if (nMarker > 1 )
            for ( int k = 0;k<(nMarkerNormal-1);++k )
                fileAllDir << listMarkerNormal[k] << ",";
        fileAllDir << listMarkerNormal[nMarkerNormal-1] << "};\n";


        for ( auto const& myMark : listMarkerLinesThetaPhiFixed )
        {
            size_type theDofThetaPhi = myMark.first;
            int nMarker2 = myMark.second.size();

            fileAllDir << "Physical Line(\"ThetaPhi" << theDofThetaPhi <<"\") = {";
            if (nMarker2 > 1 )
                for ( int k = 0;k<(nMarker2-1);++k )
                    fileAllDir << myMark.second[k] << ",";
                    fileAllDir << myMark.second[nMarker2-1] << "};\n";

            fileAllDir << "Physical Line(\"HatThetaPhi" << theDofThetaPhi <<"\") = {";
            if (nMarker2 > 1 )
                for ( int k = 0;k<(nMarker2-1);++k )
                    fileAllDir << listMarkerLinesHatThetaPhiFixed[theDofThetaPhi][k] << ",";
                    fileAllDir << listMarkerLinesHatThetaPhiFixed[theDofThetaPhi][nMarker2-1] << "};\n";
        }

        fileAllDir.close();
    }





#if 1 // commment up to end







    typedef std::vector<boost::reference_wrapper<typename MeshTraits<mesh_type>::face_type const> > cont_range_type;
    boost::shared_ptr<cont_range_type> myelts( new cont_range_type );


    if ( false && Environment::isMasterRank() )
        for( size_type dof_thetaPhi : dofToUse_thetaPhi )
        {
            // vecteur direction
            vDirection[0] = vx_proj->operator()(dof_thetaPhi);
            if( nDim >= 2 )
                vDirection[1] = vy_proj->operator()(dof_thetaPhi);
            if( nDim == 3 )
                vDirection[2] = vz_proj->operator()(dof_thetaPhi);
                std::cout << "vDirection[0] " << vDirection[0] << " vDirection[1] " << vDirection[1] << " vDirection[2] " << vDirection[2] << "\n";

        }

    //-----------------------------------------------------------------//
    //-----------------------------------------------------------------//
    //-----------------------------------------------------------------//
    boost::mpi::timer mytimer;
    bool verbose = boption(_name="acoustic-model.verbose");
    // time loop
    while ( !myBdf[0]->isFinished() )
    {
        double time = myBdf[0]->time();
        if ( Environment::isMasterRank() )
            std::cout << "\n-----------------------------------\n"
                      << "time " << time << "\n";

        size_type cpt_dof_thetaPhi = 0;
        for( size_type dof_thetaPhi : dofToUse_thetaPhi )
        {
            if ( Environment::isMasterRank() )
            {
                if ( verbose )
                    std::cout << "dof_thetaPhi : " << ++cpt_dof_thetaPhi << "/" << dofToUse_thetaPhi.size() << "\n";
                else
                {
                    std::cout << "\r";
                    std::cout << "dof_thetaPhi : " << ++cpt_dof_thetaPhi << "/" << dofToUse_thetaPhi.size() << std::flush;
                }
            }

            mat->zero();
            rhs->zero();
            //USol->zero();
            *USol = *ULoc_thetaPhi[dof_thetaPhi];

            // vector transport direction
            vDirection[0] = vx_proj->operator()(dof_thetaPhi);
            if( nDim >= 2 )
                vDirection[1] = vy_proj->operator()(dof_thetaPhi);
            if( nDim == 3 )
                vDirection[2] = vz_proj->operator()(dof_thetaPhi);
            //if ( Environment::isMasterRank() )
            //    std::cout << "vDirection[0] " << vDirection[0] << " vDirection[1] " << vDirection[1] << " vDirection[2] " << vDirection[2] << "\n";

            //-----------------------------------------------------------------//
            mytimer.restart();
            // time discretisation and source term
            form2(_test=Xh,_trial=Xh,_matrix=mat) +=
                integrate(_range=elements(mesh),
                          _expr=myBdf[dof_thetaPhi]->polyDerivCoefficient(0)*idt(ULoc)*id(ULoc) );

            //auto polyDerivInTime = myBdf[dof_thetaPhi]->polyDeriv();
            //auto polyDerivInTime = (couplingDirection)? bdfEnergyDensity->polyDeriv() : myBdf[dof_thetaPhi]->polyDeriv();
            auto polyDerivInTime = myBdf[dof_thetaPhi]->polyDeriv();
            auto rhsExpr = idv(projSourceTerm)+idv(polyDerivInTime);
            //auto rhsExpr = idv(polyDerivInTime);

            form1(_test=Xh,_vector=rhs) +=
                integrate(_range=elements(mesh),
                          _expr=rhsExpr*id(ULoc) );

            // propagation term
            auto vDirExpr = vec( cst(vDirection[0]),cst(vDirection[1]),cst(vDirection[2]));
            form2(_test=Xh,_trial=Xh,_matrix=mat) +=
                integrate(_range=elements(mesh),
                          _expr=gradt(ULoc)*vDirExpr*id(ULoc) );

            double tElapsed = mytimer.elapsed();
            if ( verbose && Environment::isMasterRank() )
                std::cout << " assemblage base done in " << tElapsed << "s\n";

            //-----------------------------------------------------------------//
            // stabilisation
            if ( doStab )
            {
                mytimer.restart();
                // SUPG stab
                auto coeff = vf::h()/( 2*norm2( vDirExpr ) );
                auto residualStabForm2 = myBdf[dof_thetaPhi]->polyDerivCoefficient(0)*idt(ULoc) + gradt(ULoc)*vDirExpr;
                auto residualStabForm1 = rhsExpr;//idv(polyDerivInTime);
                form2(_test=Xh,_trial=Xh,_matrix=mat) +=
                    integrate(_range=elements(mesh),
                              _expr= coeff*(grad(ULoc)*vDirExpr)*residualStabForm2 );
                form1(_test=Xh,_vector=rhs) +=
                    integrate(_range=elements(mesh),
                              _expr=coeff*(grad(ULoc)*vDirExpr)*residualStabForm1 );
                tElapsed = mytimer.elapsed();
                if ( verbose && Environment::isMasterRank() )
                    std::cout << " assemblage stab done in " << tElapsed << "s\n";
            }

            //-----------------------------------------------------------------//
            // boundary conditions
            mytimer.restart();

            projBCDirichlet->zero();
            myelts->clear();//resize(0);
            std::vector<bool> dofdone(Xh->nLocalDof(),false);
            for ( auto const& face : boundaryfaces(mesh) )
            {
                size_type faceId = face.id();
#if 0
                ctx->update( face.element0(), face.pos_first() );
                auto unitNormal = ctx->unitNormal();
#endif
                auto unitNormal = mapUnitNormalByFace[ faceId ];
                double valVdotN = vDirection[0]*unitNormal[0] + vDirection[1]*unitNormal[1] + vDirection[2]*unitNormal[2];
                if ( std::abs(valVdotN) < 1e-6 ) continue;
                //continue;

#if !defined(NDEBUG)
                if ( doComputeHatThetaPhi[faceId][dof_thetaPhi] )
                {
                    double thetheta = mapThetaPhiAndHatThetaPhiNEW[faceId][dof_thetaPhi].first[0];
                    double thephi = mapThetaPhiAndHatThetaPhiNEW[faceId][dof_thetaPhi].first[1];
                    double thehattheta = mapThetaPhiAndHatThetaPhiNEW[faceId][dof_thetaPhi].second[0];
                    double thehatphi = mapThetaPhiAndHatThetaPhiNEW[faceId][dof_thetaPhi].second[1];

                    double theVx = math::sin(thetheta)*math::cos(thephi);
                    double theVy = math::sin(thetheta)*math::sin(thephi);
                    double theVz = math::cos(thetheta);
                    double theHatVx = math::sin(thehattheta)*math::cos(thehatphi);
                    double theHatVy = math::sin(thehattheta)*math::sin(thehatphi);
                    double theHatVz = math::cos(thehattheta);

                    double res = theVx*unitNormal[0] + theVy*unitNormal[1] + theVz*unitNormal[2] + theHatVx*unitNormal[0] + theHatVy*unitNormal[1] + theHatVz*unitNormal[2];
                    CHECK ( std::abs(res)<1e-12 ) << "invalid vn != hatv n : " << res << "\n";
                }

                //std::cout << "face id " << face.id() << "normal " << face.element0().normal(face.pos_first())
                //          <<  " ctx->normal(0) " << ctx->unitNormal()//ctx->normal()
                //          << "\n";
#endif

                if ( doComputeHatThetaPhi[faceId][dof_thetaPhi] ) //valVdotN < 0 )
                    myelts->push_back(boost::cref(face));
                //continue;

                for ( uint16_type fDof = 0 ; fDof<Xh->dof()->nLocalDofOnFace(true) ; ++fDof)
                {
                    const size_type thedof = Xh->dof()->faceLocalToGlobal(face.id(),fDof,0).get<0>();
                    if ( dofdone[thedof] ) continue;

                    //if ( doComputeHatThetaPhi[faceId][dof_thetaPhi] )//valVdotN < 0 )
                    //{
                    //for( size_type dof_thetaPhi2 : dofToUse_thetaPhi )
                    //    wBCintegrate->set( dof_thetaPhi2,ULoc_thetaPhi[ dof_thetaPhi2 ]->operator()( thedof ) );

                    for( size_type dof_thetaPhi2=0 ; dof_thetaPhi2<Xh_thetaPhi->nDof() ; ++dof_thetaPhi2 )
                    {
                        if ( !extrapUseBdf )
                        {
                            //use solution at last time step
                            if ( dofIsUsed_thetaPhi[dof_thetaPhi2] )
                            {
                                //if ( dof_thetaPhi == 1 && dof_thetaPhi2 == 1 && std::abs(ublas::column(face.G(),0)[1] -6. ) < 1.1  && std::abs(ublas::column(face.G(),0)[2] - 12. ) < 1e-5  )
                                //std::cout << "val1 :" << ULoc_thetaPhi[dof_thetaPhi2]->operator()( thedof ) << "\n";
                                wBCintegrate->set(dof_thetaPhi2,ULoc_thetaPhi[dof_thetaPhi2]->operator()( thedof ) );
                            }
                            else
                            {
                                //if ( dof_thetaPhi == 1 && dof_thetaPhi2 == 1 && std::abs(ublas::column(face.G(),0)[1] -6. ) < 1.1  && std::abs(ublas::column(face.G(),0)[2] - 12. ) < 1e-5  )
                                //std::cout << "val2 :" << ULoc_thetaPhi[mapPeriodicDof_thetaPhi.find(dof_thetaPhi2)->second]->operator()( thedof ) << "\n";
                                wBCintegrate->set(dof_thetaPhi2,ULoc_thetaPhi[mapPeriodicDof_thetaPhi.find(dof_thetaPhi2)->second]->operator()( thedof ) );
                            }
                        }
                        else
                        {
                            // extrapolation using bdf
                            if ( dofIsUsed_thetaPhi[dof_thetaPhi2] )
                                wBCintegrate->set(dof_thetaPhi2,myBdf[dof_thetaPhi2]->poly()( thedof ) );
                            else
                                wBCintegrate->set(dof_thetaPhi2,myBdf[mapPeriodicDof_thetaPhi.find(dof_thetaPhi2)->second]->poly()( thedof ) );
                        }
                    }

                    if ( doComputeHatThetaPhi[faceId][dof_thetaPhi] )//valVdotN < 0 )
                    {
                        //boundary condition with d=0 : w(r,theta,phi,t)=-alpha*w(r,\hat{theta},\hat{phi},t)
                        if ( useFirstTermInBC )
                        {
                            double wHatThetaPhi = 0;
                            if ( false )
                            {
                                size_type idElt = localisationHatThetaPhi[dof_thetaPhi];
                                //std::cout << " idElt "<< idElt << "\n";
                                for ( uint16_type lDof = 0 ; lDof<Xh_thetaPhi->dof()->nLocalDof(true) ; ++lDof)
                                {
                                    const size_type thedof2 = Xh_thetaPhi->dof()->localToGlobal(idElt,lDof,0).index();
                                    if ( dofIsUsed_thetaPhi[thedof2] )
                                        wHatThetaPhi+= ULoc_thetaPhi[ thedof2 ]->operator()( thedof );
                                    else
                                        wHatThetaPhi+= ULoc_thetaPhi[ mapPeriodicDof_thetaPhi.find(thedof2)->second ]->operator()( thedof );
                                }
                                wHatThetaPhi/=Xh_thetaPhi->dof()->nLocalDof(true);
                            }
                            else
                            {
                                auto evalAtHatThetaPhi = evaluateFromContext( _context=*ctxEvalHatThetaPhiNEW[faceId],
                                                                              _expr=idv(wBCintegrate) );
                                wHatThetaPhi = evalAtHatThetaPhi( mapEvalFromCtxThetaPhiNEW[faceId][dof_thetaPhi],0 );
                                /*if ( dof_thetaPhi == 1 && std::abs(ublas::column(face.G(),0)[1] -6. ) < 1.1  && std::abs(ublas::column(face.G(),0)[2] - 0. ) < 1e-5 &&
                                     ublas::column(face.G(),0)[0] > 6 && ublas::column(face.G(),0)[0] < 10  )
                                     std::cout << " wHatThetaPhi " << wHatThetaPhi << "\n";*/
                            }

                            // update projBCDirichlet at this dof
                            // warning : not minus!
                            projBCDirichlet->add( thedof, alpha*(1-d_prob)*wHatThetaPhi );
                        }



                        if ( useSecondTermInBC )
                        {
                            double bcIntegrate = inner_product(*mapOperatorDiffuseBC.find(faceId)->second,wBCintegrate->container());
#if !defined(NDEBUG)
                            auto thetaPrime = Px();
                            auto phiPrime = Py();
                            auto vPrime = vSoundVelocity*vec( sin(thetaPrime)*cos(phiPrime), sin(thetaPrime)*sin(phiPrime), cos(thetaPrime) );
                            auto vPrimeDotN = vPrime(0,0)*unitNormal[0] + vPrime(1,0)*unitNormal[1] + vPrime(2,0)*unitNormal[2];
                            auto chiVPrimeDotN = chi( vPrimeDotN > cst(-1e-6) );
                            auto bcIntegrateExpr = (alpha*d_prob/(M_PI*vSoundVelocity))*vPrimeDotN*idv(wBCintegrate)*sin(thetaPrime);

                            double bcIntegrate2 = integrate(_range=elements(mesh_thetaPhi),
                                                            _expr=bcIntegrateExpr*chiVPrimeDotN ).evaluate(false)(0,0);
                            CHECK( std::abs( bcIntegrate - bcIntegrate2 ) < 1e-8 ) << " error between opt " << bcIntegrate <<" and eval : " << bcIntegrate2 << "\n";
#endif
                            projBCDirichlet->add( thedof, bcIntegrate );
                        }


                    } // if ( valVdotN < 0 )

                    dofdone[thedof] = true;
                } // for ( fDof ... )
            } // for ( auto const& face : boundaryfaces(mesh) )

            if ( true )
            {
                auto myMarkedFacesBCRange = boost::make_tuple( mpl::size_t<MESH_FACES>(),
                                                               myelts->begin(),
                                                               myelts->end(),
                                                               myelts );

                // up mat and rhs for bc Dirichlet condition
                form2(_test=Xh,_trial=Xh,_matrix=mat) +=
                    on(_range=myMarkedFacesBCRange,//boundaryfaces(mesh),
                       _rhs=rhs, _element=*ULoc, _expr=idv(projBCDirichlet) );
            }

            tElapsed = mytimer.elapsed();
            if ( verbose && Environment::isMasterRank() )
                std::cout << " assemblage bc done in " << tElapsed << "s\n";

            //-----------------------------------------------------------------//
            mytimer.restart();
            mybackend->solve(_matrix=mat, _rhs=rhs, _solution=USol,_prec=myprec );
            //if ( Environment::isMasterRank() )
            //    std::cout << " USol->l2Norm() " << USol->l2Norm() << "\n";

            // copy algebraic vector into finite element function
            *ULoc_thetaPhi[dof_thetaPhi] = *USol;

            tElapsed = mytimer.elapsed();
            if ( verbose && Environment::isMasterRank() )
                std::cout << " solve done in " << tElapsed << "s\n";


        } // dof_thetaPhi

        //-----------------------------------------------------------------//
        // compute field of interest
        ULocSumAllDirection->zero();
        ULocIntegrateAllDirection->zero();
        auto energyProjOnThetaPhi = Xh_thetaPhi->elementPtr();
        for ( size_type dof=0;dof<Xh->nLocalDof();++dof)
        {
            // type 1 : sum all direction method
            for( size_type dof_thetaPhi : dofToUse_thetaPhi )
                ULocSumAllDirection->add( dof, ULoc_thetaPhi[dof_thetaPhi]->operator()( dof ) );
            //ULocSumAllDirection->set( dof, ULocSumAllDirection->operator()(dof)/dofToUse_thetaPhi.size()  );

            // type 2 : nodal projection on ThetaPhi space and integrate
            for( size_type dof_thetaPhi=0 ; dof_thetaPhi<Xh_thetaPhi->nDof() ; ++dof_thetaPhi )
            {
                if ( dofIsUsed_thetaPhi[dof_thetaPhi] )
                    energyProjOnThetaPhi->set(dof_thetaPhi,ULoc_thetaPhi[dof_thetaPhi]->operator()( dof ) );
                else
                    energyProjOnThetaPhi->set(dof_thetaPhi,ULoc_thetaPhi[mapPeriodicDof_thetaPhi.find(dof_thetaPhi)->second]->operator()( dof ) );
            }
            auto theta = Px();
            double nodalVal = integrate(_range=elements(mesh_thetaPhi),
                                        _expr=idv(energyProjOnThetaPhi)*sin(theta) ).evaluate(false)(0,0);
            ULocIntegrateAllDirection->set( dof, nodalVal );
        }

        //-----------------------------------------------------------------//
        // export solutions
        e->step( time )->add( "sound-intensity-sum", *ULocSumAllDirection );
        e->step( time )->add( "sound-intensity-integrate", *ULocIntegrateAllDirection );
        if ( doExportSolForEachVecDir )// (dof_thetaPhi % 5) == 0 )
            for( size_type dof_thetaPhi : dofToUse_thetaPhi )
                e->step(time)->add( (boost::format("NEWULoc%1%")%dof_thetaPhi).str(), *(ULoc_thetaPhi[dof_thetaPhi]) );

        auto chiPositiveExpr = chi( idv(ULocIntegrateAllDirection) > 1e-8 );
        auto ensurePositiveExpr = 1e-5*(1-chiPositiveExpr);
        auto soundPressureLevelExpr = 10*vf::log( idv(ULocIntegrateAllDirection)*airDensity*vSoundVelocity/cst(std::pow(refPressure,2) ) );
        //auto soundPressureLevelExpr = 10*vf::log( ( idv(ULocIntegrateAllDirection)*airDensity*vSoundVelocity/cst(std::pow(refPressure,2) ) )*chiPositiveExpr + ensurePositiveExpr );
        //auto soundPressureLevelExpr = 10*vf::log( (idv(ULocIntegrateAllDirection)+cst(1.))*airDensity*vSoundVelocity/cst(std::pow(refPressure,2) ) );
        //*soundPressureLevel = vf::project( _space=Xh,_expr=soundPressureLevelExpr );
        soundPressureLevel->on(_range=elements(Xh->mesh()),_expr=soundPressureLevelExpr );
        e->step( time )->add( "sound-pressure-level", *soundPressureLevel );

        e->save();

        //-----------------------------------------------------------------//
        // update time step
        for( size_type dof_thetaPhi : dofToUse_thetaPhi )
            myBdf[dof_thetaPhi]->next( *ULoc_thetaPhi[dof_thetaPhi] );

        //ULocIntegrateAllDirection->scale(1./(4*M_PI)/*dofToUse_thetaPhi.size()*/);
        //bdfEnergyDensity->next( *ULocIntegrateAllDirection );
        //ULocIntegrateAllDirection->scale(4*M_PI/*dofToUse_thetaPhi.size()*/);
        //bdfEnergyDensity->next( *ULocSumAllDirection );


    } // time loop

#endif

    return 0;

}
