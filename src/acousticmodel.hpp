#ifndef __ACOUSTICMODEL_H
#define __ACOUSTICMODEL_H 1

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/pch.hpp>

#include <feel/feelvf/vf.hpp>


namespace Feel
{


class AcousticModel
{
public :

    typedef Mesh<Hypercube<2>> mesh_2d_type;
    typedef boost::shared_ptr<mesh_2d_type> mesh_2d_ptrtype;
    typedef FunctionSpace<mesh_2d_type, bases<Lagrange<1, Scalar> > > space_2d_type;
    typedef boost::shared_ptr<space_2d_type> space_2d_ptrtype;
    typedef space_2d_type::element_type element_2d_type;
    typedef boost::shared_ptr<element_2d_type> element_2d_ptrtype;

    AcousticModel() {}

    void init();

    void initThetaPhi();

private :

    mesh_2d_ptrtype M_meshThetaPhi;


}; // class AcousticModel

} // namespace Feel

#endif //__ACOUSTICMODEL_H
