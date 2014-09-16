#ifndef __ACOUSTICMODEL_H
#define __ACOUSTICMODEL_H 1

#include <feel/feelalg/backend.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>


namespace Feel
{


class AcousticModel
{
public :

    //typedef Mesh<Hypercube<nDim>> mesh_type;
    typedef Mesh<Simplex<3>> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef FunctionSpace<mesh_type, bases<Lagrange<1, Scalar> > > space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef space_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;

    typedef Bdf<space_type> bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrtype;
    // exporter
    typedef Exporter<mesh_type,1> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    typedef Mesh<Hypercube<2>> mesh_thetaphi_type;
    typedef boost::shared_ptr<mesh_thetaphi_type> mesh_thetaphi_ptrtype;
    typedef FunctionSpace<mesh_thetaphi_type, bases<Lagrange<1, Scalar> > > space_thetaphi_type;
    typedef boost::shared_ptr<space_thetaphi_type> space_thetaphi_ptrtype;
    typedef space_thetaphi_type::element_type element_thetaphi_type;
    typedef boost::shared_ptr<element_thetaphi_type> element_thetaphi_ptrtype;
    // evaluation context
    typedef space_thetaphi_type::Context context_type;
    typedef boost::shared_ptr<context_type> context_ptrtype;

    typedef Backend<double> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef backend_type::vector_ptrtype vector_ptrtype;
    typedef Preconditioner<double> preconditioner_type;
    typedef boost::shared_ptr<preconditioner_type> preconditioner_ptrtype;

    AcousticModel();
    void loadMesh( mesh_ptrtype mesh);
    void init();
    void exportResults( double time );

    bool hasFinishedTimeStep();
    void nextTimeStep();
    double time();
    bool isRestart() { return M_bdfDensity[*M_dofToUseInThetaPhi.begin()]->isRestart(); }

    template < typename ExprT >
    void updateInitialSolution( vf::Expr<ExprT> const& expr )
    {
        // only for visu
        M_ULocIntegrateAllDirection->on( _range=elements(M_mesh),_expr=expr );

        for( size_type dofThetaPhi : M_dofToUseInThetaPhi )
            *M_ULoc_thetaPhi[dofThetaPhi] = *M_ULocIntegrateAllDirection;
    }

    void solve();
private :
    void createThetaPhi();
    void createTransportFE();
    void createPrecomputeBC();

    void updateAssemblyTransportModel( sparse_matrix_ptrtype & mat,vector_ptrtype &rhs,size_type dofThetaPhi );
    void updateAssemblyBC( sparse_matrix_ptrtype & mat,vector_ptrtype & rhs, size_type dofThetaPhi );
private :

    bool M_verbose;

    mesh_thetaphi_ptrtype M_meshThetaPhi;
    space_thetaphi_ptrtype M_XhThetaPhi;
    std::set<size_type> M_dofToUseInThetaPhi;
    element_thetaphi_ptrtype M_vxProj, M_vyProj, M_vzProj;
    std::map<size_type,size_type> M_mapPeriodicDofThetaPhi;
    std::vector<bool> M_dofIsUsed_thetaPhi;

    mesh_ptrtype M_mesh;
    space_ptrtype M_Xh;
    element_ptrtype /*ULoc,ULocSumAllDirection,*/M_ULocIntegrateAllDirection,M_projBCDirichlet,M_projSourceTerm;
    element_ptrtype M_soundPressureLevel;
    std::vector<element_ptrtype> M_ULoc_thetaPhi;

    std::vector<bdf_ptrtype> M_bdfDensity;
    export_ptrtype M_exporter;

    double M_vSoundVelocity;
    double M_refPressure;
    double M_airDensity;
    double M_M;
    double M_alpha;
    double M_d_prob;


    std::map<size_type, node_type> M_mapUnitNormalByFace;
    std::map<size_type,std::map<size_type,bool> > M_doComputeHatThetaPhi;
    std::map<size_type,context_ptrtype> M_ctxEvalHatThetaPhi;
    std::map<size_type,std::map<size_type,size_type> > M_mapEvalFromCtxThetaPhi;
    std::map<size_type,vector_ptrtype > M_mapOperatorDiffuseBC;

    std::vector<backend_ptrtype> M_backends;
    std::vector<preconditioner_ptrtype> M_preconditioners;
    std::vector<sparse_matrix_ptrtype> M_mat;
    std::vector<vector_ptrtype> M_rhs;
    bool M_useOneMatVecByDirection, M_useOneBackendByDirection;

}; // class AcousticModel


Feel::po::options_description acousticModel_options();


} // namespace Feel

#endif //__ACOUSTICMODEL_H
