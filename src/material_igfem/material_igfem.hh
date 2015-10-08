/**
 * @file   element_class_igfem.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 *
 * @brief  Parent material for IGFEM
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material.hh"
#include "igfem_internal_field.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_IGFEM_HH__
#define __AKANTU_MATERIAL_IGFEM_HH__

/* -------------------------------------------------------------------------- */
namespace akantu {
  class SolidMechanicsModelIGFEM;
}

__BEGIN_AKANTU__


class MaterialIGFEM : public virtual Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  typedef FEEngineTemplate<IntegratorGauss,
		      ShapeLagrange, _ek_igfem> MyFEEngineIGFEMType;

  MaterialIGFEM(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialIGFEM();

protected:
 
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  inline void interpolateInternal(const Element new_el,
				  const Element old_el,
				  Vector<Real> & interpolated,
				  const Vector<Real> & internal,
				  const UInt nb_quads_new,
				  const UInt nb_quads_old);

  virtual void computeAllStresses(GhostType ghost_type = _not_ghost);


  /* ------------------------------------------------------------------------ */
  /* MeshEventHandler inherited members                                       */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  virtual void computeQuadraturePointsCoordinates(ElementTypeMapArray<Real> & quadrature_points_coordinates,
						  const GhostType & ghost_type) const;
  // virtual void onElementsAdded(const Array<Element> & element_list,
  //                              const NewElementsEvent & event) {};

  // virtual void onElementsRemoved(const Array<Element> & element_list,
  //                                const ElementTypeMapArray<UInt> & new_numbering,
  //                                const RemovedElementsEvent & event) {};
protected:

  /// constitutive law
  virtual void computeStress(__attribute__((unused)) ElementType el_type,
                             __attribute__((unused)) GhostType ghost_type = _not_ghost)  {
   
  }
  void initialize();

  
  virtual ElementTypeMap<UInt> getInternalDataPerElem(const ID & id,
						      const ElementKind & element_kind,
						      const ID & fe_engine_id) const;

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
public:
  virtual inline UInt getNbDataForElements(const Array<Element> & elements,
                                           SynchronizationTag tag) const;

  virtual inline void packElementData(CommunicationBuffer & buffer,
                                      const Array<Element> & elements,
                                      SynchronizationTag tag) const;

  virtual inline void unpackElementData(CommunicationBuffer & buffer,
                                        const Array<Element> & elements,
                                        SynchronizationTag tag);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

protected:

  const UInt nb_sub_materials = 2;

  /// pointer to the solid mechanics model for igfem elements
  SolidMechanicsModelIGFEM * model;

  ///internal field of bool to know to which sub-material a quad point belongs
  IGFEMInternalField<UInt> sub_material; 

  /// material name of first sub-material
  std::string name_sub_mat_1;

  /// material name of first sub-material
  std::string name_sub_mat_2;

};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_igfem_inline_impl.cc"



__END_AKANTU__

#include "igfem_internal_field_tmpl.hh"

#endif /* __AKANTU_MATERIAL_IGFEM_HH__ */
