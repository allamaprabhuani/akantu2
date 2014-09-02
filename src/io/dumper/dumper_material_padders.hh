
#ifndef __AKANTU_DUMPER_MATERIAL_PADDERS_HH__
#define __AKANTU_DUMPER_MATERIAL_PADDERS_HH__
/* -------------------------------------------------------------------------- */
#include "dumper_padding_helper.hh"
/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__
__BEGIN_AKANTU_DUMPER__
/* -------------------------------------------------------------------------- */


template<class T, class R>
class MaterialPadder : public PadderGeneric<Vector<T>, R > {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:
  MaterialPadder(const SolidMechanicsModel & model) : 
    model(model),
    element_index_by_material(model.getElementIndexByMaterial()) { }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

  /// return the material from the global element index
  const Material & getMaterialFromGlobalIndex(Element global_index){

    UInt index = global_index.getIndex();
    UInt material_id = element_index_by_material(global_index.getType())(index);
    const Material & material = model.getMaterial(material_id);
    return material;
  }

  /// return the type of the element from global index
  const ElementType getElementTypeFromGlobalIndex(Element global_index){
    return global_index.getType();
  }

protected:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

  /// all material padders probably need access to solid mechanics model
  const SolidMechanicsModel & model;
  /// they also need an access to the map from global ids to material id and local ids
  const ElementTypeMapArray<UInt>  & element_index_by_material;
  /// the number of data per element
  const ElementTypeMapArray<UInt>  nb_data_per_element;

};

/* -------------------------------------------------------------------------- */

template <UInt spatial_dimension>
class StressPadder : 
  public MaterialPadder<Real,Matrix<Real> > {

public:
  StressPadder(const SolidMechanicsModel & model) :
    MaterialPadder<Real, Matrix<Real> >(model){
    this->setPadding(3,3);
  }

  inline Matrix<Real> func(const Vector<Real> & in, Element global_element_id){

    UInt nrows = spatial_dimension;    
    UInt ncols = in.size() / nrows;
    UInt nb_data = in.size() / ncols;

    Matrix<Real> stress = this->pad(in, nrows,ncols, nb_data);
    const Material & material = this->getMaterialFromGlobalIndex(global_element_id);
    bool plane_strain = !material.getParam<bool>("Plane_Stress");
    if(plane_strain) {
      Real nu = material.getParam<Real>("nu");
      for (UInt d = 0; d < nb_data; ++d) {
	stress(2, 2 + 3*d) = nu * (stress(0, 0 + 3*d) + stress(1, 1 + 3*d));
      }
    }
    return stress;
  }

  UInt getDim(){return 9;};

  UInt getNbComponent(UInt old_nb_comp){
    return this->getDim();
  };

};

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
class StrainPadder : public MaterialPadder<Real, Matrix<Real> > {
public:
  StrainPadder(const SolidMechanicsModel & model) : 
    MaterialPadder<Real, Matrix<Real> >(model) {
    this->setPadding(3,3);
  }

  inline Matrix<Real> func(const Vector<Real> & in, Element global_element_id){

    UInt nrows = spatial_dimension;    
    UInt ncols = in.size() / nrows;
    UInt nb_data = in.size() / ncols;

    Matrix<Real> strain = this->pad(in, nrows,ncols, nb_data);
    const Material & material = this->getMaterialFromGlobalIndex(global_element_id);
    bool plane_stress = material.getParam<bool>("Plane_Stress");
    if(plane_stress) {
      Real nu = material.getParam<Real>("nu");
      for (UInt d = 0; d < nb_data; ++d) {
	strain(2, 2 + 3*d) = nu / (nu - 1) * (strain(0, 0 + 3*d) + strain(1, 1 + 3*d));
      }
    }
    return strain;
  }

  UInt getDim(){return 9;};

  UInt getNbComponent(UInt old_nb_comp){
    return this->getDim();
  };

};

/* -------------------------------------------------------------------------- */

__END_AKANTU_DUMPER__
__END_AKANTU__


#endif /* __AKANTU_DUMPER_MATERIAL_PADDERS_HH__ */
