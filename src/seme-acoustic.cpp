/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:set fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <feel/feelcore/environment.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>

using namespace Feel;


std::vector< boost::shared_ptr< Mesh<Hypercube<1,1,1> > > >
createFDSubMesh1D( boost::shared_ptr< Mesh<Hypercube<1,1,1> > > mesh )
{
    std::vector< boost::shared_ptr< Mesh<Hypercube<1,1,1> > > > res(1);
    res[0] = mesh;
    return res;
}

std::vector< boost::shared_ptr< Mesh<Hypercube<1,1,2> > > >
createFDSubMesh1D( boost::shared_ptr< Mesh<Hypercube<2,1,2> > > mesh )
{
    std::vector< boost::shared_ptr< Mesh<Hypercube<1,1,2> > > > res(2);
    res[0] = createSubmesh(mesh,markedfaces(mesh,"bordX") );
    res[1] = createSubmesh(mesh,markedfaces(mesh,"bordY") );
    return res;
}

std::vector< boost::shared_ptr< Mesh<Hypercube<1,1,3> > > >
createFDSubMesh1D( boost::shared_ptr< Mesh<Hypercube<3,1,3> > > mesh )
{
    std::vector< boost::shared_ptr< Mesh<Hypercube<1,1,3> > > > res(3);
    res[0] = createSubmesh(mesh,markededges(mesh,"bordX") );
    res[1] = createSubmesh(mesh,markededges(mesh,"bordY") );
    res[2] = createSubmesh(mesh,markededges(mesh,"bordZ") );
    return res;
}

int main(int argc, char**argv )
{
	po::options_description accousticoptions( "Accoustic options" );
	accousticoptions.add_options()
		( "meshsizeX", po::value<double>()->default_value( 0.1 ), "coeff" )
		( "meshsizeY", po::value<double>()->default_value( 0.1 ), "coeff" )
		( "meshsizeZ", po::value<double>()->default_value( 0.1 ), "coeff" )

		( "time-step", po::value<double>()->default_value( 0.1 ), "coeff" )
		( "time-final", po::value<double>()->default_value( 0.1 ), "coeff" )
		( "scaling-coeff", po::value<double>()->default_value( 5. ), "coeff" )
		( "M", po::value<double>()->default_value( 0.01 ), "coeff" )
		( "stab", po::value<bool>()->default_value( true ), "coeff" )
		( "stab-rho", po::value<double>()->default_value( 0.25 ), "coeff" )
        ( "gmsh.filename-thetaPhi", po::value<std::string>(), "name for theta phi mesh" )
        ( "alpha", po::value<double>()->default_value( 0.5 ), "alpha coeff" )
        ( "d-prob", po::value<double>()->default_value( 0.5 ), "d coeff" )
        ( "sound-velocity", po::value<double>()->default_value( 343 ), "alpha coeff" )
		( "use-second-term-bc", po::value<bool>()->default_value( true ), "coeff" )


        ( "center_x", po::value<double>()->default_value( 1 ), "alpha coeff" )
        ( "center_y", po::value<double>()->default_value( 0.5 ), "alpha coeff" )
        ( "center_z", po::value<double>()->default_value( 0.5 ), "alpha coeff" )
        ( "radius", po::value<double>()->default_value( 0.1 ), "alpha coeff" )
        ( "vel0", po::value<double>()->default_value( 5 ), "vel0" )
        ( "vel1", po::value<double>()->default_value( 0 ), "vel1" )
        ( "nvel", po::value<int>()->default_value( 1 ), "nvel" )
		;

	Environment env( _argc=argc, _argv=argv,
                     _desc=accousticoptions.add(feel_options()).add(backend_options("prec0")).add(backend_options("prec1")),
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

    typedef Mesh<Hypercube<nDim>> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef FunctionSpace<mesh_type, bases<Lagrange<1, Scalar> > > space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef space_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;
    typedef FunctionSpace<mesh_type, bases<Lagrange<1, Vectorial> > > space_vec_type;
    typedef boost::shared_ptr<space_vec_type> space_vec_ptrtype;
    typedef space_vec_type::element_type element_vec_type;
    typedef boost::shared_ptr<element_vec_type> element_vec_ptrtype;

    typedef Mesh<Hypercube<2>> mesh_2d_type;
    typedef boost::shared_ptr<mesh_2d_type> mesh_2d_ptrtype;
    typedef FunctionSpace<mesh_2d_type, bases<Lagrange<1, Scalar> > > space_2d_type;
    typedef boost::shared_ptr<space_2d_type> space_2d_ptrtype;
    typedef space_2d_type::element_type element_2d_type;
    typedef boost::shared_ptr<element_2d_type> element_2d_ptrtype;


    double meshSizeX = doption(_name="meshsizeX");
    double meshSizeY = doption(_name="meshsizeY");
    double meshSizeZ = doption(_name="meshsizeZ");
    double timeStep = doption(_name="time-step");
    double timeFinal = doption(_name="time-final");
    double doStab = boption(_name="stab");
    double rho = doption(_name="stab-rho");
    double M = doption(_name="M");
    double scalCoef = doption(_name="scaling-coeff");
    double alpha = doption(_name="alpha");
    double d_prob = doption(_name="d-prob");//0.5;
    bool useSecondTermInBC = boption(_name="use-second-term-bc");
    //-----------------------------------------------------------------//
    // create FE mesh/space and 1d meshes/spaces for bases dim
    mesh_ptrtype mesh;
    space_ptrtype Xh;
    element_ptrtype ULoc,ULocSumAllDirection;
    mesh_1d_ptrtype meshX,meshY,meshZ;
    mesh_2d_ptrtype mesh_thetaPhi;
    space_1d_ptrtype XhX, XhY,XhZ;
    space_2d_ptrtype Xh_thetaPhi;
    space_vec_ptrtype XhVec;
    element_vec_ptrtype ULocVec;

    mesh_thetaPhi = loadMesh(_mesh=new mesh_2d_type, _filename=soption(_name="gmsh.filename-thetaPhi"));
    Xh_thetaPhi = space_2d_type::New(_mesh=mesh_thetaPhi);

    if (nDim == 1 )
    {
        GeoTool::Node x1(0);
        GeoTool::Node x2(5);
        GeoTool::Line geoX( meshSizeX,"structH",x1,x2);
        geoX.setMarker(_type="point",_name="Boundary",_markerAll=true);
        geoX.setMarker(_type="line",_name="OmegaX",_markerAll=true);
        //meshX = geoX.createMesh(_mesh=new mesh_1d_type,
        //                        _name = "meshX" );
        mesh = geoX.createMesh(_mesh=new mesh_type,
                               _name = "meshX" );
        Xh = space_type::New(_mesh=mesh);
        ULoc = Xh->elementPtr();

        auto submeshes1d = createFDSubMesh1D(mesh);
        auto meshX = submeshes1d[0];
        XhX = space_1d_type::New(_mesh=meshX);
        //UX = XhX->elementPtr();
    }
    else if (nDim >= 2 )
    {
        mesh = loadMesh(_mesh=new mesh_type);
        Xh = space_type::New(_mesh=mesh);
        ULoc = Xh->elementPtr();
        ULocSumAllDirection = Xh->elementPtr();

        auto submeshes1d = createFDSubMesh1D(mesh);
        auto meshX = submeshes1d[0];//createSubmesh(mesh,markedfaces(mesh,"bordX") );
        XhX = space_1d_type::New(_mesh=meshX);
        auto meshY = submeshes1d[1];//createSubmesh(mesh,markedfaces(mesh,"bordY") );
        XhY = space_1d_type::New(_mesh=meshY);
        if (nDim == 3 )
        {
            auto meshZ = submeshes1d[2];
            XhZ = space_1d_type::New(_mesh=meshZ);
            //XhVec = space_vec_type::New(_mesh=mesh);
            //ULocVec = XhVec->elementPtr();
        }
    }

    //-----------------------------------------------------------------//

    size_type nDofX = 5./meshSizeX+1;
    CHECK( XhX->nDof() == nDofX ) << "invalid XhX->nDof() "<< XhX->nDof() << " != nDofX " << nDofX << "\n";
    size_type nDofY = 1./meshSizeY+1;
    if (nDim >= 2 )
        CHECK( XhY->nDof() == nDofY ) << "invalid XhY->nDof() "<< XhY->nDof() << " != nDofY " << nDofY << "\n";
    size_type nDofZ = 1./meshSizeZ+1;
    if (nDim >= 3 )
        CHECK( XhZ->nDof() == nDofZ ) << "invalid XhZ->nDof() "<< XhZ->nDof() << " != nDofZ " << nDofZ << "\n";

    std::cout << "nDofX : " << nDofX << "\n";
    if (nDim >= 2 )
        std::cout << "nDofY : " << nDofY << "\n";
    if (nDim >= 3 )
        std::cout << "nDofZ : " << nDofZ << "\n";

    //size_type startDofX = 0, startDofY = nDofX;
    if (nDim >= 2 )
    {
        size_type theFDndof = 0;
        /**/ if (nDim == 1 ) theFDndof = nDofX;
        else if (nDim == 2 ) theFDndof = nDofX*nDofY;
        else if (nDim == 3 ) theFDndof = nDofX*nDofY*nDofZ;

        CHECK( theFDndof == mesh->numPoints() ) << "invalid theFDndof "
                                                  <<  theFDndof << " != " << mesh->numPoints()
                                                  << " with nDofX " << nDofX
                                                  << " with nDofY " << nDofY
                                                  << " with nDofZ " << nDofZ
                                                  << " \n";
        CHECK( theFDndof == Xh->nDof() ) << "invalid nDofX*nDofY " << theFDndof << " Xh->nDof() "  << Xh->nDof() << "\n";
    }

    std::map<size_type,size_type> mapFEtoFD,mapFDtoFE;
    for ( size_type dof=0;dof<Xh->nDof();++dof)
    {
        auto ptDof = Xh->dof()->dofPoint(dof).get<0>();
        size_type idX = std::floor( ptDof[0]/meshSizeX+0.01*meshSizeX );
        size_type idY = (nDim >= 2 ) ? std::floor( ptDof[1]/meshSizeY+0.01*meshSizeY ) : invalid_size_type_value;
        size_type idZ = (nDim == 3 ) ? std::floor( ptDof[2]/meshSizeZ+0.01*meshSizeZ ) : invalid_size_type_value;

        size_type dofDF = invalid_size_type_value;
        if (nDim == 1 )
            dofDF = idX;
        else if ( nDim == 2 )
            dofDF = idX+idY*nDofX;
        else if ( nDim == 3 )
            dofDF = idX+idY*nDofX+idZ*nDofX*nDofY;

        CHECK( dofDF < Xh->nDof() ) << "invalid dofDF " << dofDF << " and Xh->nDof() " << Xh->nDof()
                                    << " idX " << idX << " idY " << idY
                                    << " ptDof " << ptDof
                                    << "\n";

        CHECK( mapFEtoFD.find(dof) == mapFEtoFD.end() ) << "aie1\n";
        CHECK( mapFDtoFE.find(dofDF) == mapFDtoFE.end() ) << "aie2\n";
        mapFEtoFD[dof]=dofDF;
        mapFDtoFE[dofDF]=dof;
    }

    //-----------------------------------------------------------------//
    // initial solution
    //double center_x=1,center_y=0.5,center_z=0.5;
    double center_x=doption(_name="center_x"),center_y=doption(_name="center_y"),center_z=doption(_name="center_z");
    //double center_x=10*meshSizeX,center_y=5*meshSizeX,center_z=5*meshSizeX;
    //double radius=2*meshSizeX;
    double radius=doption(_name="radius");//2*meshSizeX;
#if 0
    if ( nDim == 1 )
        ULoc->on( _range=elements(mesh), _expr=-scalCoef*Px()*(Px()-cst(1.))*chi(Px()<=1) );
    else if ( nDim == 2 )
        ULoc->on( _range=elements(mesh), _expr=-scalCoef*Px()*(Px()-cst(1.))*chi(Px()<=1)*Py()*(1.-Py()) );
    else if ( nDim == 3 )
        ULoc->on( _range=elements(mesh), _expr=-scalCoef*Px()*(Px()-cst(1.))*chi(Px()<=1)*Py()*(1.-Py())*Pz()*(1.-Pz()) );
#else
    auto chiInitSol3d = chi( pow(Px()-center_x,2) + pow(Py()-center_y,2) + pow(Pz()-center_z,2) <= cst(std::pow(radius,2)) );
    if ( nDim == 1 )
        ULoc->on( _range=elements(mesh), _expr=-scalCoef*Px()*(Px()-cst(1.))*chi(Px()<=1) );
    else if ( nDim == 2 )
        ULoc->on( _range=elements(mesh), _expr=-scalCoef*Px()*(Px()-cst(1.))*chi(Px()<=1)*Py()*(1.-Py()) );
    else if ( nDim == 3 )
        ULoc->on( _range=elements(mesh), _expr=-scalCoef*(Px()-(center_x-radius))*(Px()-(center_x+radius))*(Py()-(center_y-radius))*(Py()-(center_y+radius))*
                  (Pz()-(center_z-radius))*(Pz()-(center_z+radius))*chiInitSol3d );
#endif
    //-----------------------------------------------------------------//
    // export initial solution

    auto e2 = exporter( _mesh=mesh,_name="myexporter2" );
    auto e3 = exporter( _mesh=mesh,_name="myexporter" );
#if 0
    //e->step(0)->add( "ULoc", *ULoc );
    for( size_type dof_thetaPhi=0 ; dof_thetaPhi<Xh_thetaPhi->nDof() ; ++dof_thetaPhi )
        e->step(0)->add( (boost::format("ULoc%1%")%dof_thetaPhi).str(), *ULoc );

    e->step(0)->add( "ULocSumAllDirection", *ULocSumAllDirection );

    if (nDim == 3 && false )
        e->step(0)->add( "ULocVec", *ULocVec );
    e->save();
#endif

    std::cout << "export done" << std::endl;

    //-----------------------------------------------------------------//
    // create algebraic structure
    auto rhs = backend()->newVector( Xh->nDof(), Xh->nDof() );/*nDofX*nDofY*/
    auto USol = backend()->newVector( Xh->nDof(), Xh->nDof() );
    auto mat = backend()->newMatrix(Xh,Xh);
    // copy initial solution to vector solution
    for ( size_type k=0;k<Xh->nDof();++k)
        USol->set( k, ULoc->operator()( mapFDtoFE[k] ) );

    //-----------------------------------------------------------------//
    //vitesse de propagation
    std::vector<std::vector<double>> vDirection(Xh_thetaPhi->nDof(),std::vector<double>(nDim,0.));
    // ???
    double vStrange = doption(_name="sound-velocity");//1;

    element_2d_ptrtype vx_proj, vy_proj, vz_proj;
    vx_proj=Xh_thetaPhi->elementPtr();
    vy_proj=Xh_thetaPhi->elementPtr();
    vz_proj=Xh_thetaPhi->elementPtr();
    //vx_proj->on( _range=elements(mesh_thetaPhi), _expr=cos(Px()) );
    //vy_proj->on( _range=elements(mesh_thetaPhi), _expr=sqrt(pow(sin(Px()),2))*cos(Py()) );
    //vz_proj->on( _range=elements(mesh_thetaPhi), _expr=sqrt(pow(sin(Px()),2))*sin(Py()) );
    //vx_proj->on( _range=elements(mesh_thetaPhi), _expr=cos(Px()) );
    //vy_proj->on( _range=elements(mesh_thetaPhi), _expr=sin(Px())*cos(Py()) );
    //vz_proj->on( _range=elements(mesh_thetaPhi), _expr=sin(Px())*sin(Py()) );
    vx_proj->on( _range=elements(mesh_thetaPhi), _expr=vStrange*sin(Px())*cos(Py()) );
    vy_proj->on( _range=elements(mesh_thetaPhi), _expr=vStrange*sin(Px())*sin(Py()) );
    vz_proj->on( _range=elements(mesh_thetaPhi), _expr=vStrange*cos(Px()) );


    std::ofstream fileAllDir("allvectordirection.geo", std::ios::out);
    fileAllDir << "myh=0.5;\n";
    fileAllDir << "Point(1) = {0, 0, 0, myh};\n";
    for( size_type dof_thetaPhi=0 ; dof_thetaPhi<Xh_thetaPhi->nDof() ; ++dof_thetaPhi )
        fileAllDir << "Point(" << dof_thetaPhi+2 <<") = {"
                   << vx_proj->operator()(dof_thetaPhi) << ","
                   << vy_proj->operator()(dof_thetaPhi) << ","
                   << vz_proj->operator()(dof_thetaPhi) << ", myh};\n";

    for( size_type dof_thetaPhi=0 ; dof_thetaPhi<Xh_thetaPhi->nDof() ; ++dof_thetaPhi )
        fileAllDir << "Line(" << dof_thetaPhi+1 <<") = {1,"<< dof_thetaPhi+2 << "};\n";
    fileAllDir.close();



    std::cout << "Xh_thetaPhi->nDof() = " << Xh_thetaPhi->nDof() << std::endl;
    std::vector<element_ptrtype> ULoc_thetaPhi(Xh_thetaPhi->nDof(), Xh->elementPtr() );
    for( size_type dof_thetaPhi=0 ; dof_thetaPhi<Xh_thetaPhi->nDof() ; ++dof_thetaPhi )
        for ( size_type dof=0;dof<Xh->nDof();++dof)
            ULoc_thetaPhi[dof_thetaPhi]->set(dof, ULoc->operator()(dof) );

    auto pcf =  mesh->gm()->preComputeOnFaces( mesh->gm(), mesh->gm()->referenceConvex().barycenterFaces() );
    auto firstFacePtr = mesh->beginFace();
    auto ctx = mesh->gm()->context<vm::POINT|vm::NORMAL|vm::KB|vm::JACOBIAN>( firstFacePtr->element0(),
                                                                              pcf,
                                                                              firstFacePtr->pos_first() );
    auto wBCintegrate = Xh_thetaPhi->elementPtr();

    auto locMeshThetaPhi = mesh_thetaPhi->tool_localization();
    //locMeshThetaPhi->init();
    std::vector<size_type> localisationHatThetaPhi( Xh_thetaPhi->nDof(), invalid_size_type_value );
    for( size_type dof_thetaPhi=0; dof_thetaPhi<Xh_thetaPhi->nDof(); ++dof_thetaPhi)
    {
        auto thetaPhiPt = Xh_thetaPhi->dof()->dofPoint(dof_thetaPhi).get<0>();
        node_type hatThetaPhi(2);
        hatThetaPhi[0] = M_PI - thetaPhiPt[0]; // [0,pi]
        hatThetaPhi[1] = ( thetaPhiPt[1]<=M_PI ) ? M_PI -thetaPhiPt[1] : M_PI -thetaPhiPt[1] + 2*M_PI; // [0,2*pi]
        CHECK( hatThetaPhi[0] >= 0 &&  hatThetaPhi[0] <= M_PI ) << " hhihih " << hatThetaPhi[0] << "\n";
        CHECK( hatThetaPhi[1] >= 0 &&  hatThetaPhi[1] <= 2*M_PI ) << " hhihih" << hatThetaPhi[1] << "\n";
#if 0
        auto resLocalisation = locMeshThetaPhi->searchElement(hatThetaPhi);
        CHECK( resLocalisation.get<0>() ) << " aieee pas trouver!!\n";
        localisationHatThetaPhi[dof_thetaPhi] = resLocalisation.get<1>();
        std::cout << " the pt " << hatThetaPhi << " is localised in element " << mesh_thetaPhi->element( resLocalisation.get<1>() ).G() << "\n";
#else
        bool hasFind=false;
        for ( auto const& elt : elements(mesh_thetaPhi) )
        {
            std::vector<double> ptX,ptY;
            for (int k=0;k<4;++k)
            {
                auto const& pt = ublas::column(elt.G(),k);
                ptX.push_back( pt[0] );
                ptY.push_back( pt[1] );
            }
            double xmin = *std::min_element( ptX.begin(),ptX.end() );
            double xmax = *std::max_element( ptX.begin(),ptX.end() );
            double ymin = *std::min_element( ptY.begin(),ptY.end() );
            double ymax = *std::max_element( ptY.begin(),ptY.end() );

            double eps=1e-3;
            if ( hatThetaPhi[0] >= (xmin-eps) && hatThetaPhi[0] <= (xmax+eps) &&
                 hatThetaPhi[1] >= (ymin-eps) && hatThetaPhi[1] <= (ymax+eps) )
            {
                localisationHatThetaPhi[dof_thetaPhi] = elt.id();
                hasFind=true;
                break;
            }
        }
        //if( hasFind ) std::cout <<  "has Find the pt " << hatThetaPhi << "\n";
        CHECK( hasFind ) <<  "hasNotFind the pt " << hatThetaPhi << "\n";

#endif
    }


    std::vector<boost::shared_ptr<Backend<double> > > mybackends(Xh_thetaPhi->nDof(),backend(_rebuild=true));


    //for( size_type dof_thetaPhi=0; dof_thetaPhi<ioption("nvel")/*Xh_thetaPhi->nDof()*/; ++dof_thetaPhi)
    ULocSumAllDirection->zero();
    int N = (ioption("nvel" )==0)?Xh_thetaPhi->nDof():ioption("nvel" );
    //for( size_type dof_thetaPhi=0; dof_thetaPhi<Xh_thetaPhi->nDof(); ++dof_thetaPhi)
    for( size_type dof_thetaPhi=0; dof_thetaPhi<N; ++dof_thetaPhi)
        {
            auto rhs = backend()->newVector( Xh->nDof(), Xh->nDof() );/*nDofX*nDofY*/
            auto USol = backend()->newVector( Xh->nDof(), Xh->nDof() );
            auto mat = backend()->newMatrix(Xh,Xh);
            //mat->zero();
            //rhs->zero();
           //USol->zero();
            //rhs = backend()->newVector( Xh->nDof(), Xh->nDof() );
            //auto mat2 = backend()->newMatrix(Xh,Xh);
            //auto rhs2 = backend()->newVector( Xh->nDof(), Xh->nDof() );
            //auto USol2 = backend()->newVector( Xh->nDof(), Xh->nDof() );
            auto mat2 = mat;
            auto rhs2 = rhs;
            auto USol2 = USol;
            auto mstr = (boost::format("mat%1%.m")%dof_thetaPhi).str();
            auto str = (boost::format("myexporter%1%")%dof_thetaPhi).str();
            auto e1 = exporter( _mesh=mesh,_name=str );

        for ( size_type dof=0;dof<Xh->nDof();++dof)
            ULoc_thetaPhi[dof_thetaPhi]->set(dof, ULoc->operator()(dof) );
    // time loop
    for ( double time=timeStep ; time<timeFinal ; time+=timeStep )
    {
        std::cout << "time " << time << "\n";
        mat->zero();
        rhs->zero();
        USol->zero();

        if ( ioption("nvel" ) == 0 )
        {
            vDirection[dof_thetaPhi][0] = vx_proj->operator()(dof_thetaPhi);
            if( nDim >= 2 )
                vDirection[dof_thetaPhi][1] = vy_proj->operator()(dof_thetaPhi);
            if( nDim == 3 )
                vDirection[dof_thetaPhi][2] = vz_proj->operator()(dof_thetaPhi);
        }
        else
        {
            //vDirection[0] = ( dof_thetaPhi==0 )?8:0;
            //if ( (dof_thetaPhi == 0) && (time < 2*timeStep+1e-6))
            if ( (dof_thetaPhi ==  0 ))
                vDirection[dof_thetaPhi][0] = doption("vel0");
            else
                vDirection[dof_thetaPhi][0] =  doption("vel1");
            vDirection[dof_thetaPhi][1] = 0;
            vDirection[dof_thetaPhi][2] = 0;
        }
            std::cout << "vDirection[0] " << vDirection[dof_thetaPhi][0] << " vDirection[1] " << vDirection[dof_thetaPhi][1] << " vDirection[2] " << vDirection[dof_thetaPhi][2] << "\n";

            std::vector<bool> dofdone(Xh->nDof(),false);
            for ( auto const& face : boundaryfaces(mesh) )
            {
                ctx->update( face.element0(), face.pos_first() );
                auto unitNormal = ctx->unitNormal();
                double valVdotN = vDirection[dof_thetaPhi][0]*unitNormal[0]+ vDirection[dof_thetaPhi][1]*unitNormal[1] + vDirection[dof_thetaPhi][2]*unitNormal[2];
#if 0
                std::cout << "face id " << face.id() << "normal " << face.element0().normal(face.pos_first())
                          <<  " ctx->normal(0) " << ctx->unitNormal()//ctx->normal()
                          << "\n";
#endif
                for ( uint16_type fDof = 0 ; fDof<Xh->dof()->nLocalDofOnFace(true) ; ++fDof)
                {
                    const size_type thedof = Xh->dof()->faceLocalToGlobal(face.id(),fDof,0).get<0>();
                    if ( dofdone[thedof] ) continue;

                    for(size_type dof_thetaPhi2=0; dof_thetaPhi2<Xh_thetaPhi->nDof(); ++dof_thetaPhi2)
                        wBCintegrate->set( dof_thetaPhi2,ULoc_thetaPhi[ dof_thetaPhi2 ]->operator()( thedof ) );


                    size_type dofDF = mapFEtoFD[thedof];
                    // impose dirichlet condition in matrix
                    mat2->add( dofDF, dofDF, 1. );
                    //boundary condition with d=0 : w(r,theta,phi,t)=-alpha*w(r,\hat{theta},\hat{phi},t)
                    if ( valVdotN < 0 )
                    {
                        size_type idElt = localisationHatThetaPhi[dof_thetaPhi];
                        //std::cout << " idElt "<< idElt << "\n";
                        double wHatThetaPhi = 0;
                        for ( uint16_type lDof = 0 ; lDof<Xh_thetaPhi->dof()->nLocalDof(true) ; ++lDof)
                        {
                            const size_type thedof2 = Xh_thetaPhi->dof()->localToGlobal(idElt,lDof,0).index();
                            //const size_type thedof2 = Xh_thetaPhi->dof()->localToGlobal( *mesh_thetaPhi->elementIterator( idElt ),lDof,0).index();
                            //std::cout << " thedof2 "<< thedof2 << "\n";

                            wHatThetaPhi+= ULoc_thetaPhi[ thedof2 ]->operator()( thedof );
                        }
                        wHatThetaPhi/=Xh_thetaPhi->dof()->nLocalDof(true);
                        //double d_prob = 0.5;
                        rhs2->add( dofDF,-alpha*(1-d_prob)*wHatThetaPhi );

                        if ( useSecondTermInBC )
                        {
                            auto thetaPrime = Px();
                            auto phiPrime = Py();
                            auto vPrime = vec(cos(thetaPrime), sin(thetaPrime)*cos(phiPrime), sin(thetaPrime)*sin(phiPrime) );
                            //auto chiVdotN = cst(1.);//chi( inner(vPrime,N() ) < cst(0.) );
                            //auto vDootN = vPrime(0,0)*Nx() + vPrime(1,0)*Nx() + vPrime(2,0)*Nz();
                            auto vDootN = vPrime(0,0)*unitNormal[0] + vPrime(1,0)*unitNormal[1] + vPrime(2,0)*unitNormal[2];
                            auto chiVdotN = chi( vDootN < cst(0.) );
                            auto bcIntegrateExpr = alpha*d_prob/(M_PI*vStrange)*vDootN*idv(wBCintegrate)*sin(thetaPrime);
                            double bcIntegrate = integrate(_range=elements(mesh_thetaPhi),
                                                           _expr=bcIntegrateExpr*chiVdotN ).evaluate()(0,0);
                            rhs2->add( dofDF, bcIntegrate );
                        }
                    }
                    else
                    {
                        rhs2->add( dofDF, 0. );
                    }

                    dofdone[thedof] = true;
                }
            } // for ( auto const& face : boundaryfaces(mesh) )


            for ( size_type i=1;i< nDofX-1;++i)
            {
                size_type nDofInternalY = (nDim >= 2 )?nDofY-1 : 2;
                for ( size_type j=1;j< nDofInternalY;++j)
                {
                    size_type nDofInternalZ = (nDim == 3 )?nDofZ-1 : 2;
                    for ( size_type k=1;k< nDofInternalZ;++k)
                    {
                        size_type thedof = invalid_size_type_value;
                        if ( nDim == 1 )
                            thedof = i;
                        else if ( nDim == 2 )
                            thedof = i+j*nDofX;
                        else if ( nDim == 3 )
                            thedof = i+j*nDofX+k*nDofX*nDofY;
                        //std::cout << " thedof " << thedof << "with i " << i << " j " << j <<" nDofX " << nDofX << "\n";

                        // transient + reaction terms
                        mat2->add( thedof, thedof, 1./timeStep + M*vStrange );
                        //rhs2->add( thedof, 1./timeStep*(USol->operator()(thedof)) );
                        rhs2->add( thedof, (1./timeStep)*ULoc_thetaPhi[dof_thetaPhi]->operator()( mapFDtoFE[thedof] ) );


                        // grad term
#if 1
                        double valGradX = vDirection[dof_thetaPhi][0]/(2*meshSizeX);
                        mat2->add( thedof, thedof-1, -valGradX );
                        mat2->add( thedof, thedof+1, valGradX );
                        if ( nDim>=2 )
                        {
                            double valGradY = vDirection[dof_thetaPhi][1]/(2*meshSizeY);
                            mat2->add( thedof, thedof-nDofX, -valGradY );
                            mat2->add( thedof, thedof+nDofX, valGradY );
                        }
                        if ( nDim==3 )
                        {
                            double valGradZ = vDirection[dof_thetaPhi][2]/(2*meshSizeZ);
                            mat2->add( thedof, thedof-nDofX*nDofY, -valGradZ );
                            mat2->add( thedof, thedof+nDofX*nDofY, valGradZ );
                        }

#else
                        double valGrad = vDirection[dof_thetaPhi][0]/(meshSizeX);
                        //mat2->add( i, i-1, -valGrad );
                        mat2->add( i, i, -valGrad );
                        mat2->add( i, i+1, valGrad );
#endif
                        if ( doStab )
                        {
                            // artificial diffusion term
                            //double valDiff = rho*meshSizeX/std::pow(meshSizeX,2);
                            double valDiff = rho*std::min(meshSizeX,meshSizeY);
                            for(int i_dim=0; i_dim<nDim; ++i_dim)
                            {
                                if( std::abs(vDirection[dof_thetaPhi][i_dim]) > 1e-5 ) //if( v[i_dim] != 0)
                                    valDiff *= vDirection[dof_thetaPhi][i_dim];
                            }
                            double valDiffX = valDiff/std::pow(meshSizeX,2);
                            double valDiffY = valDiff/std::pow(meshSizeY,2);
                            double valDiffZ = valDiff/std::pow(meshSizeZ,2);
                            mat2->add( thedof, thedof-1, -valDiffX );
                            mat2->add( thedof, thedof, 2*valDiffX );
                            mat2->add( thedof, thedof+1, -valDiffX );

                            if ( nDim>=2 )
                            {
                                mat2->add( thedof, thedof-nDofX, -valDiffY );
                                mat2->add( thedof, thedof, 2*valDiffY );
                                mat2->add( thedof, thedof+nDofX, -valDiffY );
                            }
                            if ( nDim==3 )
                            {
                                mat2->add( thedof, thedof-nDofX*nDofY, -valDiffZ );
                                mat2->add( thedof, thedof, 2*valDiffZ );
                                mat2->add( thedof, thedof+nDofX*nDofY, -valDiffZ );
                            }


                        } // do stab

                    } // k
                } // j
            } // i
            if ( 0 )
            {
                mat->printMatlab(mstr);


                auto pstr = (boost::format("prec%1%")%dof_thetaPhi).str();
                auto b = backend(_rebuild=true,_name=pstr);
                auto p=preconditioner( _prefix=pstr,_matrix=mat,_pc=b->pcEnumType()/*LU_PRECOND*/,
                                       _pcfactormatsolverpackage=b->matSolverPackageEnumType(), _backend=b, _rebuild = true );
                b->solve(_matrix=mat, _rhs=rhs, _solution=USol, _prec = p );
                p.reset();
                b.reset();
            }
            else
            {
                auto b = backend(_rebuild=true );
                b->solve(_matrix=mat, _rhs=rhs, _solution=USol);
            }
            //mybackends[dof_thetaPhi]->solve(_matrix=mat2, _rhs=rhs2, _solution=USol2 );
            std::cout << " rhs2->l2Norm() " << rhs2->l2Norm() << "\n";
            std::cout << " USol->l2Norm() " << USol->l2Norm() << "\n";

            //for ( size_type k=0;k<Xh->nDof();++k)
            //    ULoc->set( mapFDtoFE[k], USol->operator()( k ) );//mapFEtoFD,mapFDtoFE

            for ( size_type k=0;k<Xh->nDof();++k)
                ULoc_thetaPhi[dof_thetaPhi]->set( mapFDtoFE[k], USol2->operator()( k ) );//mapFEtoFD,mapFDtoFE

            //if (dof_thetaPhi == 0 )
            //{
                std::cout << " do export\n";
                //*ULoc_thetaPhi[dof_thetaPhi] = *ULoc;
                //e->step(time)->add( "ULoc", *ULoc );
#if 0
                if ( dof_thetaPhi ==  0 )
                {
                    e1->step(time)->add( (boost::format("ULoc%1%")%dof_thetaPhi).str(), *(ULoc_thetaPhi[dof_thetaPhi])/* *ULoc*/ );
                    e1->save();
                }
                if ( dof_thetaPhi ==  1 )
                {
                    e2->step(time)->add( (boost::format("ULoc%1%")%dof_thetaPhi).str(), *(ULoc_thetaPhi[dof_thetaPhi])/* *ULoc*/ );
                    e2->save();
                }
#else
                e1->step(time)->add( (boost::format("ULoc%1%")%dof_thetaPhi).str(), *(ULoc_thetaPhi[dof_thetaPhi])/* *ULoc*/ );
                e1->save();
#endif
                //}
            //break;
        } // dof_thetaPhi

#if 0
        // compute field of interest
        ULocSumAllDirection->zero();
        for ( size_type dof=0;dof<Xh->nDof();++dof)
        {
            //size_type cptVecProp=0;
            for( size_type dof_thetaPhi=0 ; dof_thetaPhi<Xh_thetaPhi->nDof() ; ++dof_thetaPhi )
            {
                //auto mythetaPhiPt = Xh_thetaPhi->dof()->dofPoint(dof_thetaPhi).get<0>();
                //if ( mythetaPhiPt[1] >= (2*M_PI-0.00001) ) continue;
                //++cptVecProp;
                ULocSumAllDirection->add( dof, ULoc_thetaPhi[dof_thetaPhi]->operator()( dof ) );
            }
            ULocSumAllDirection->set( dof, ULocSumAllDirection->operator()(dof)/Xh_thetaPhi->nDof()  );
        }

        e3->step( time )->add( "ULocSumAllDirection", *ULocSumAllDirection );
#else
        for ( size_type dof=0;dof<Xh->nDof();++dof)
        {
            //size_type cptVecProp=0;
            for( size_type dof_thetaPhi=0 ; dof_thetaPhi<Xh_thetaPhi->nDof() ; ++dof_thetaPhi )
            {
                //auto mythetaPhiPt = Xh_thetaPhi->dof()->dofPoint(dof_thetaPhi).get<0>();
                //if ( mythetaPhiPt[1] >= (2*M_PI-0.00001) ) continue;
                //++cptVecProp;
                ULocSumAllDirection->add( dof, ULoc_thetaPhi[dof_thetaPhi]->operator()( dof ) );
            }
        }

#endif
        //e->save();

    } // time loop
    for ( size_type dof=0;dof<Xh->nDof();++dof)
    {
        ULocSumAllDirection->set( dof, ULocSumAllDirection->operator()(dof)/Xh_thetaPhi->nDof()  );
    }
    e3->add( "ULocSumAllDirection", *ULocSumAllDirection );
    e3->save();
    return 0;

}
