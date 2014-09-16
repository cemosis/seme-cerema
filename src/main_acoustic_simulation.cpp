
#include <acousticmodel.hpp>


using namespace Feel;

int main(int argc, char**argv )
{
	Environment env( _argc=argc, _argv=argv,
                     _desc=acousticModel_options(),
                     _desc_lib=feel_options(),
                     _about=about(_name="acoustic_simulation",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    AcousticModel myAcousticModel;

    auto mesh = loadMesh(_mesh=new AcousticModel::mesh_type);
    myAcousticModel.loadMesh(mesh);

    // source expression
    double scalCoef = doption(_name="scaling-coeff");
    double center_x=doption(_name="center_x"),center_y=doption(_name="center_y"),center_z=doption(_name="center_z");
    double radius=doption(_name="radius");
    bool useSourceAsInitialSolution = boption("use-source-as-initial-solution");

    auto chiInitSol3d = chi( pow(Px()-center_x,2) + pow(Py()-center_y,2) + pow(Pz()-center_z,2) <= cst(std::pow(radius,2)) );
    auto sourceTerm3dExpr = scalCoef*chiInitSol3d;
    /*auto sourceTerm3dExpr = -scalCoef*
        (Px()-(center_x-radius))*(Px()-(center_x+radius))*
        (Py()-(center_y-radius))*(Py()-(center_y+radius))*
        (Pz()-(center_z-radius))*(Pz()-(center_z+radius))*chiInitSol3d;*/

    if ( useSourceAsInitialSolution )
        myAcousticModel.updateInitialSolution( sourceTerm3dExpr );

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
