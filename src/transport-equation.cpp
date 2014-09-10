#include <feel/feel.hpp>


namespace Feel
{

template <uint16_type DIM>
class TransportEquation
{
public :
    static const uint16_type nDim = DIM;
    typedef Mesh<Simplex<nDim> > mesh_type;
    typedef Lagrange<1, Scalar,Continuous> basis_type;
    typedef FunctionSpace<mesh_type, bases<basis_type> > space_type;
    typedef typename space_type::element_type element_type;

    typedef Bdf<space_type> bdf_type;
    typedef Exporter<mesh_type,1> export_type;

    typedef Backend<double> backend_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    TransportEquation()
        :
        M_vDirection( nDim ),
        M_solverType( soption(_name="solver-type") ),
        M_hasInitAlgebraicDataStructuresSolverFEM( false ),
        M_hasInitAlgebraicDataStructuresSolverFEMFCT( false )
    {
        M_vDirection[0]=1;
    }

    void loadMesh( boost::shared_ptr<mesh_type> mesh )
    {
        M_mesh=mesh;
        M_Xh = space_type::New( mesh );// Pch<1>( M_mesh );
        M_u = M_Xh->elementPtr();
        M_bdf = bdf( _space=M_Xh, _name="mybdf" );
        M_exporter = exporter( _mesh=mesh );
    }

    template < typename ExprT >
    void updateInitialSolution( vf::Expr<ExprT> const& expr )
    {
        M_u->on( _range=elements(M_mesh),_expr=expr );
    }

    void startTimeStep()
    {
        M_bdf->start(*M_u);
    }

    bool hasFinishedTimeStep()
    {
        return M_bdf->isFinished();
    }
    void nextTimeStep()
    {
        M_bdf->next(*M_u);
    }
    double time()
    {
        return M_bdf->time();
    }
    void exportResults( double time )
    {
        M_exporter->step(time)->add( "u", *M_u );
        M_exporter->save();
    }

    std::string solverType() const { return M_solverType; }

    void solve()
    {
        if ( this->solverType() == "fem" )
            this->solveFEM();
        else if ( this->solverType() == "fem-fct" )
            this->solveFEMFCT();
    }

private :
    void initAlgebraicDataStructuresSolverFEM()
    {
        M_mat = backend()->newMatrix(_test=M_Xh,_trial=M_Xh);
        M_rhs = backend()->newVector(M_Xh);
        M_hasInitAlgebraicDataStructuresSolverFEM = true;
    }
    void initAlgebraicDataStructuresSolverFEMFCT()
    {
        M_mat = backend()->newMatrix(_test=M_Xh,_trial=M_Xh);
        M_rhs = backend()->newVector(M_Xh);


        M_matMc = backend()->newMatrix(_test=M_Xh,_trial=M_Xh);
        M_matMl = backend()->newMatrix(_test=M_Xh,_trial=M_Xh);
        M_matK = backend()->newMatrix(_test=M_Xh,_trial=M_Xh);
        M_matD = backend()->newMatrix(_test=M_Xh,_trial=M_Xh);


        form2(_test=M_Xh,_trial=M_Xh,_matrix=M_matMc) +=
            integrate(_range=elements(M_mesh),
                      _expr= M_bdf->polyDerivCoefficient(0)*idt(M_u)*id(M_u) );
        M_matMc->close();

        auto vDirExpr = vec( cst(M_vDirection[0]) );//,cst(vDirection[1]),cst(vDirection[2]));
        form2(_test=M_Xh,_trial=M_Xh,_matrix=M_matK) +=
            integrate(_range=elements(M_mesh),
                      _expr=-gradt(M_u)*vDirExpr*id(M_u) );
        M_matK->close();

        //for ( size_type i=0;i<M_Xh->nLocalDof();++i)
        //    M_matMl->set(i,i, M_matMc->operator()(i,i) );
        for ( auto const& graphData : *M_matMl->graph() )
        {
            size_type rowId = graphData.first;
            double sumCol=0;
            for ( size_type colId : graphData.second.get<2>() )
                sumCol += M_matMc->operator()( rowId, colId );
            M_matMl->set( rowId, rowId, sumCol );
        }
        M_matMl->close();

        for ( auto const& graphData : *M_matK->graph() )
        {
            size_type rowId = graphData.first;
            for (  size_type colId : graphData.second.get<2>() )
            {
                if ( rowId < colId )
                {
                    double kij = std::max(0., std::max( -M_matK->operator()( rowId,colId ), -M_matK->operator()( colId,rowId ) ) );
                    M_matD->set( rowId, colId, kij );
                    M_matD->set( colId, rowId, kij );
                }
            }
        }
        for ( auto const& graphData : *M_matK->graph() )
        {
            size_type rowId = graphData.first;
            double sumCol=0;
            for ( size_type colId : graphData.second.get<2>() )
                if ( rowId != colId )
                    sumCol += M_matD->operator()( rowId, colId );
            M_matD->set( rowId, rowId, -sumCol );
        }
        M_matD->close();

        M_hasInitAlgebraicDataStructuresSolverFEMFCT = true;
    }

    void solveFEM()
    {
        if ( !M_hasInitAlgebraicDataStructuresSolverFEM )
            this->initAlgebraicDataStructuresSolverFEM();

        this->updateAssemblySolverFEM();

        backend()->solve(_matrix=M_mat,_rhs=M_rhs,_solution=*M_u);
    }

    void solveFEMFCT()
    {
        if ( !M_hasInitAlgebraicDataStructuresSolverFEMFCT )
            this->initAlgebraicDataStructuresSolverFEMFCT();

        for ( int k=0;k<15;++k )
        {
            this->updateAssemblySolverFEMFCT();

            backend()->solve(_matrix=M_mat,_rhs=M_rhs,_solution=*M_u);
        }
    }

    void updateAssemblySolverFEM()
    {
        bool doStab = boption(_name="do-stab");
        double theta=doption(_name="theta");
        M_mat->zero();
        M_rhs->zero();
        auto bf = form2(_test=M_Xh,_trial=M_Xh,_matrix=M_mat);
        auto lf = form1(_test=M_Xh,_vector=M_rhs);
        auto u = M_u;
        auto v = M_u;
        auto vDirExpr = vec( cst(M_vDirection[0]) );//,cst(vDirection[1]),cst(vDirection[2]));

        // time discretisation and source term
        bf += integrate(_range=elements(M_mesh),
                        _expr=M_bdf->polyDerivCoefficient(0)*idt(u)*id(v) );

        auto polyDerivInTime = M_bdf->polyDeriv();
        auto rhsExpr = /*idv(projSourceTerm)+*/idv(polyDerivInTime);

        lf += integrate(_range=elements(M_mesh),
                        _expr=rhsExpr*id(v) );

        // advection term
        bf += integrate(_range=elements(M_mesh),
                        _expr=theta*gradt(u)*vDirExpr*id(v) );

        if ( std::abs(theta-1)>1e-8 )
        {
            lf += integrate(_range=elements(M_mesh),
                            _expr=-(1.-theta)*gradv(u)*vDirExpr*id(v) );
        }

        // stabilisation
        if ( doStab )
        {
            // SUPG stab
            auto coeff = vf::h()/( 2*norm2( vDirExpr ) );
            auto residualStabForm2 = M_bdf->polyDerivCoefficient(0)*idt(u) + gradt(u)*vDirExpr;
            //auto residualStabForm2 = M_bdf->polyDerivCoefficient(0)*idt(u) + theta*gradt(u)*vDirExpr;
            auto residualStabForm1 = rhsExpr;
            //auto residualStabForm1 = rhsExpr -(1.-theta)*gradv(u)*vDirExpr;
            bf += integrate(_range=elements(M_mesh),
                            _expr= coeff*(/*theta**/grad(v)*vDirExpr)*residualStabForm2 );
            lf += integrate(_range=elements(M_mesh),
                            _expr=coeff*(/*theta**/grad(v)*vDirExpr)*residualStabForm1 );
        }

    }

    void updateAssemblySolverFEMFCT()
    {
        double theta=doption(_name="theta");
        auto uold = M_bdf->unknown(0);

        auto UoldVec = backend()->newVector(M_Xh);
        auto UoldVec2 = backend()->newVector(M_Xh);

        *UoldVec = uold;// *M_u;
        M_mat->zero();
        M_rhs->zero();

        // time discretisation and source term
        M_mat->addMatrix( 1., M_matMl );
        M_rhs->addVector( UoldVec, M_matMl );
        //mat->addMatrix( 1., matMc );
        //rhs->addVector( Uold, matMc );

        M_mat->addMatrix( -theta, M_matK );
        M_mat->addMatrix( -theta, M_matD );

        if ( std::abs(theta-1)>1e-8 )
        {
            *UoldVec2 = uold;//*M_u;
            UoldVec2->scale(1.-theta );
            M_rhs->addVector( UoldVec2, M_matK );
            //Uold2->scale(-1.);
            M_rhs->addVector( UoldVec2, M_matD );
        }

        auto M_matF = backend()->newMatrix(_test=M_Xh,_trial=M_Xh);
        auto PiPlusVec = backend()->newVector(_test=M_Xh), PiMinusVec= backend()->newVector(_test=M_Xh);
        auto QiPlusVec = backend()->newVector(_test=M_Xh), QiMinusVec= backend()->newVector(_test=M_Xh);
        auto RiPlusVec = backend()->newVector(_test=M_Xh), RiMinusVec= backend()->newVector(_test=M_Xh);
        auto FiVec = backend()->newVector(_test=M_Xh);

        for ( auto const& graphData : *M_matMl->graph() )
        {
            size_type rowId = graphData.first;
            for ( size_type colId : graphData.second.get<2>() )
            {
                if ( colId == rowId ) continue;
                double ui = M_u->operator()(rowId);
                double uj = M_u->operator()(colId);
                if ( std::abs(ui-uj) < 1e-6 ) continue;

                double uoldi = uold(rowId);
                double uoldj = uold(colId);
                double Fmij = M_matMc->operator()(rowId,colId)*(ui-uj-uoldi+uoldj);
                double Dij = M_matD->operator()(rowId,colId);
                double Fdij = theta*Dij*(ui-uj) + (1.-theta)*Dij*(uoldi-uoldj);

                double Fij = Fmij + Fdij;
                if ( Fij*(uj-ui) > 0 ) Fij=0; // prelimiting
                M_matF->set( rowId,colId, Fij);
            }
        }

        for ( auto const& graphData : *M_matMl->graph() )
        {
            size_type rowId = graphData.first;
            double PiPlus=0,PiMinus=0;
            double QiPlus=0,QiMinus=0;
            //double uimin=INT_MAX,uimax=INT_MIN;

            for ( size_type colId : graphData.second.get<2>() )
            {
                if ( colId == rowId ) continue;
                double Fij = M_matF->operator()( rowId, colId );
                PiPlus+= std::max(0., Fij);
                PiMinus+= std::min(0., Fij);

                double ui = M_u->operator()(rowId);
                double uj = M_u->operator()(colId);

                QiPlus = std::max( QiPlus, M_matMl->operator()(rowId,rowId)*(uj-ui) );
                QiMinus = std::min( QiMinus, M_matMl->operator()(rowId,rowId)*(uj-ui) );
            }
            PiPlusVec->add(rowId,PiPlus);
            PiMinusVec->add(rowId,PiMinus);
            QiPlusVec->set(rowId,QiPlus);
            QiMinusVec->set(rowId,QiMinus);

        }
        for ( auto const& graphData : *M_matMl->graph() )
        {
            size_type rowId = graphData.first;
            if ( std::abs(PiPlusVec->operator()(rowId) ) < 1e-8 || std::abs(PiMinusVec->operator()(rowId) ) < 1e-8 ) continue;
            double RiPlus = std::min(1., QiPlusVec->operator()(rowId)/PiPlusVec->operator()(rowId) );
            double RiMinus = std::min(1., QiMinusVec->operator()(rowId)/PiMinusVec->operator()(rowId) );
            RiPlusVec->set(rowId,RiPlus);
            RiMinusVec->set(rowId,RiMinus);
        }

        for ( auto const& graphData : *M_matMl->graph() )
        {
            size_type rowId = graphData.first;
            double FiSum=0;
            for ( size_type colId : graphData.second.get<2>() )
            {
                if ( colId == rowId ) continue;

                double Fij = M_matF->operator()( rowId,colId );
                if ( std::abs(Fij) < 1e-6 ) continue;

                double alphaij=0;
                if ( Fij > 0 )
                    alphaij = std::min( RiPlusVec->operator()(rowId), RiMinusVec->operator()(colId) );
                else
                    alphaij = std::min( RiMinusVec->operator()(rowId), RiPlusVec->operator()(colId) );
                FiSum += alphaij*Fij;
            }
            FiVec->set( rowId, FiSum );
        }
        M_rhs->add(1.,FiVec);

#if 0
        auto uold = M_bdf->unknown(0);
        for ( auto const& graphData : *M_matMl->graph() )
        {
            size_type rowId = graphData.first;
            double sumCol=0;
            for ( size_type colId : graphData.second.get<2>() )
            {
                double ui = M_u->operator()(rowId);
                double uj = M_u->operator()(colId);
                if ( std::abs(ui-uj) < 1e-6 ) continue;
                double uoldi = uold(rowId);
                double uoldj = uold(colId);
                double Fmij = M_matMc->operator()(rowId,colId)*(ui-uj-uoldi+uoldj);
                double Fdij = M_matD->operator()(rowId,colId)*(ui-uj);
                double pij = (Fdij+Fmij)/(uj-ui);
                double Fij = std::min(0,pij)*(uj-ui);
                double Fprimeij = std::min(-pij, M_matK->operator()(rowId,colId)+M_matD->operator()(rowId,colId) )*(ui-uj);
                double DeltaFij = Fij-Fprimeij;

                //M_matDeltaF->set( rowId, colId, Fij - Fprimeij );
            }
        }
        //auto M_matFstar
#endif
    }

private :

    boost::shared_ptr<mesh_type> M_mesh;
    boost::shared_ptr<space_type> M_Xh;
    boost::shared_ptr<element_type> M_u;
    boost::shared_ptr<bdf_type> M_bdf;
    boost::shared_ptr<export_type> M_exporter;

    node_type M_vDirection;

    std::string M_solverType;
    bool M_hasInitAlgebraicDataStructuresSolverFEM;
    bool M_hasInitAlgebraicDataStructuresSolverFEMFCT;

    sparse_matrix_ptrtype M_mat;
    vector_ptrtype M_rhs;

    sparse_matrix_ptrtype M_matMc,M_matMl,M_matK,M_matD;
    sparse_matrix_ptrtype M_matFstar,M_matDeltaFstar;

};


} // namespace Feel


namespace Feel
{

boost::shared_ptr<Mesh<Simplex<1>>>
createMesh( mpl::int_<1> /**/ )
{
    typedef Mesh<Simplex<1>> mesh_type;
    double meshSize1d = doption(_name="gmsh.hsize");
    GeoTool::Line geoX( meshSize1d,"structH",
                        GeoTool::Node(0),GeoTool::Node(1));
    geoX.setMarker(_type="point",_name="BoundaryLeft",_marker1=true);
    geoX.setMarker(_type="point",_name="BoundaryRight",_marker2=true);
    geoX.setMarker(_type="line",_name="OmegaX",_markerAll=true);
    auto mesh = geoX.createMesh(_mesh=new mesh_type,
                                _name = "meshX" );
    return mesh;
}

template <uint16_type DIM>
boost::shared_ptr<Mesh<Simplex<DIM>>>
createMesh()
{
    return createMesh( mpl::int_<DIM>() );
}

} // namespace Feel


int main(int argc, char**argv )
{
    //# marker1 #
    using namespace Feel;
	po::options_description transportequation( "transport-equation options" );
	transportequation.add_options()
        ( "mu", po::value<double>()->default_value( 1.0 ), "coeff" )
        ( "use-mesh-relation", po::value<bool>()->default_value( true ), "coeff" )
        ( "do-stab", po::value<bool>()->default_value( true ), "coeff" )
        ( "theta", po::value<double>()->default_value( 0.5 ), "coeff" )
        ( "solver-type", po::value<std::string>()->default_value( "fem" ), "coeff" )
		;

	Environment env( _argc=argc, _argv=argv,
                   _desc=transportequation,
                   _about=about(_name="transport-equation",
                                _author="Feel++ Consortium",
                                _email="feelpp-devel@feelpp.org"));


    TransportEquation<1> transportEq;
    // load mesh
    auto mesh = createMesh<1>();
    transportEq.loadMesh( mesh );
    // load initial solution from an expression
    auto solInitialExpr = sqrt( 1- pow((Px()-0.2)/0.15 ,2)*chi(abs(Px()-0.2)<=0.15) )*chi(abs(Px()-0.2)<=0.15);
    transportEq.updateInitialSolution( solInitialExpr );
    // export initial solution
    transportEq.exportResults(0);

    // start time loop
    for ( transportEq.startTimeStep() ; !transportEq.hasFinishedTimeStep() ; transportEq.nextTimeStep() )
    {
        if ( Environment::isMasterRank() )
        {
            std::cout << "------------------------------------------------------------\n";
            std::cout << "Time " << transportEq.time() << "s\n";
        }
        // solve current time step
        transportEq.solve();
        // export visualisation results (with paraview)
        transportEq.exportResults( transportEq.time() );
    }

    return 0;
}
