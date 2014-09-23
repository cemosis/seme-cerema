
#include <acousticmodel.hpp>


using namespace Feel;

Feel::po::options_description
acousticSimulation_options()
{
	po::options_description accousticoptions( "Accoustic-Simulation options" );
	accousticoptions.add_options()

        ( "use-initial-solution", po::value<bool>()->default_value( true ), "use-initial-solution" )
        ( "initial-solution.scaling-coeff", po::value<double>()->default_value( 5. ), "coeff" )
        ( "initial-solution.center_x", po::value<double>()->default_value( 1 ), "x-position of source" )
        ( "initial-solution.center_y", po::value<double>()->default_value( 0.5 ), "y-position of source" )
        ( "initial-solution.center_z", po::value<double>()->default_value( 0.5 ), "z-position of source" )
        ( "initial-solution.radius", po::value<double>()->default_value( 0.1 ), "radius of source" )

        ( "use-sound-source", po::value<bool>()->default_value( false ), "use-sound-source" )
        ( "sound-source.scaling-coeff", po::value<double>()->default_value( 5. ), "coeff" )
        ( "sound-source.center_x", po::value<double>()->default_value( 1 ), "x-position of source" )
        ( "sound-source.center_y", po::value<double>()->default_value( 0.5 ), "y-position of source" )
        ( "sound-source.center_z", po::value<double>()->default_value( 0.5 ), "z-position of source" )
        ( "sound-source.radius", po::value<double>()->default_value( 0.1 ), "radius of source" )

        ;
    return accousticoptions;

}
int main(int argc, char**argv )
{
	Environment env( _argc=argc, _argv=argv,
                     _desc=acousticModel_options().add( acousticSimulation_options() ),
                     _desc_lib=feel_options(),
                     _about=about(_name="acoustic_simulation",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    AcousticModel myAcousticModel;

    auto mesh = loadMesh(_mesh=new AcousticModel::mesh_type);
    myAcousticModel.loadMesh(mesh);

    // source expression
    double initSolScalCoef = doption(_name="initial-solution.scaling-coeff");
    double initSolCenterX=doption(_name="initial-solution.center_x");
    double initSolCenterY=doption(_name="initial-solution.center_y");
    double initSolCenterZ=doption(_name="initial-solution.center_z");
    double initSolRadius=doption(_name="initial-solution.radius");
    bool useSourceAsInitialSolution = boption("use-initial-solution");
    auto initSolChi = chi( pow(Px()-initSolCenterX,2) + pow(Py()-initSolCenterY,2) + pow(Pz()-initSolCenterZ,2) <= cst(std::pow(initSolRadius,2)) );
    auto initSolExpr = initSolScalCoef*initSolChi;
    /*auto sourceTerm3dExpr = -scalCoef*
        (Px()-(center_x-radius))*(Px()-(center_x+radius))*
        (Py()-(center_y-radius))*(Py()-(center_y+radius))*
        (Pz()-(center_z-radius))*(Pz()-(center_z+radius))*chiInitSol3d;*/
    if ( useSourceAsInitialSolution )
        myAcousticModel.updateInitialSolution( initSolExpr );

    double SourceScalCoef = doption(_name="sound-source.scaling-coeff");
    double SourceCenterX=doption(_name="sound-source.center_x");
    double SourceCenterY=doption(_name="sound-source.center_y");
    double SourceCenterZ=doption(_name="sound-source.center_z");
    double SourceRadius=doption(_name="sound-source.radius");
    bool useSoundSource = boption("use-sound-source");
    auto soundSourceChi = chi( pow(Px()-SourceCenterX,2) + pow(Py()-SourceCenterY,2) + pow(Pz()-SourceCenterZ,2) <= cst(std::pow(SourceRadius,2)) );
    auto soundSourceExpr = SourceScalCoef*soundSourceChi*(pow(Px()-SourceCenterX,2) + pow(Py()-SourceCenterY,2) + pow(Pz()-SourceCenterZ,2) - cst(std::pow(SourceRadius,2)) );
    if ( useSoundSource )
        myAcousticModel.updateSoundSource( soundSourceExpr );



    myAcousticModel.init();

    if ( !myAcousticModel.isRestart() )
        myAcousticModel.exportResults(0.);



    // start time loop
    for ( ; !myAcousticModel.hasFinishedTimeStep() ; myAcousticModel.nextTimeStep() )
    {
        if ( Environment::isMasterRank() )
        {
            std::cout << "------------------------------------------------------------\n";
            std::cout << "Time " << myAcousticModel.time() << "s\n";
        }
        // solve current time step
        myAcousticModel.solve();
        // export visualisation results (with paraview)
        myAcousticModel.exportResults( myAcousticModel.time() );
    }


    return 0;
}
