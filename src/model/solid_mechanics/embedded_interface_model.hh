/**
 * @file   embedded_interface_model.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Mon Mar 9 2015
 * @date last modification: Mon Mar 9 2015
 *
 * @brief  Model of Solid Mechanics with embedded interfaces
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_EMBEDDED_INTERFACE_MODEL_HH__
#define __AKANTU_EMBEDDED_INTERFACE_MODEL_HH__

#include "aka_common.hh"

#include "solid_mechanics_model.hh"
#include "mesh.hh"

#include "embedded_interface_intersector.hh"

#include <CGAL/Cartesian.h>

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

typedef CGAL::Cartesian<Real> K;

/**
 * @brief Solid mechanics model using the embedded model.
 *
 * This SolidMechanicsModel subclass implements the embedded model,
 * a method used to represent 1D elements in a finite elements model
 * (eg. reinforcements in concrete).
 *
 * In addition to the SolidMechanicsModel properties, this model has
 * a mesh of the 1D elements embedded in the model, and an instance of the
 * MeshGeomContainer class for the computation of the intersections of the
 * 1D elements with the background (bulk) mesh.
 */
class EmbeddedInterfaceModel : public SolidMechanicsModel {

  typedef FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_regular> MyFEEngineType;


  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  /// Constructor
  EmbeddedInterfaceModel(Mesh & mesh,
                         Mesh & primitive_mesh,
                         UInt spatial_dimension = _all_dimensions,
                         const ID & id = "embedded_interface_model",
                         const MemoryID & memory_id = 0);

  /// Destructor
  virtual ~EmbeddedInterfaceModel();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// Initialise the model
  virtual void initFull(const ModelOptions & options = default_solid_mechanics_model_options);

  /// Initialise the materials
  virtual void initMaterials();

  /// Dump
  virtual dumper::Field * createElementalField(const std::string & field_name,
					       const std::string & group_name,
					       bool padding_flag,
					       const ElementKind & kind,
                 const std::string & fe_engine_id = "");

#ifndef SWIG
  /// give the amount of data per element
  ElementTypeMap<UInt> getInternalDataPerElem(const std::string & field_name,
					     const ElementKind & kind,
               const std::string & fe_engine_id = "");
#endif


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get interface mesh
  AKANTU_GET_MACRO(InterfaceMesh, *interface_mesh, Mesh &);

  /// get associated elements
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(InterfaceAssociatedElements,
                                   interface_mesh->getData<Element>("associated_element"),
                                   Element);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// Interface mesh (weak reference)
  Mesh * interface_mesh;

  /// Mesh used to create the CGAL primitives for intersections
  Mesh & primitive_mesh;

  /// Material selector for interface
  MaterialSelector * interface_material_selector;
  
  /// Intersector object to build the interface mesh
  EmbeddedInterfaceIntersector intersector;

};

template<typename T>
class InterfaceMeshDataMaterialSelector : public MaterialSelector {};

template<>
class InterfaceMeshDataMaterialSelector<std::string> : public MaterialSelector {
public:
  InterfaceMeshDataMaterialSelector(const std::string & name, const EmbeddedInterfaceModel & model) :
    names(model.getInterfaceMesh().getData<std::string>(name)), model(model) {}

  UInt operator() (const Element & element) {
    try {
      DebugLevel dbl = debug::getDebugLevel();
      debug::setDebugLevel(dblError);
      std::string material_name = names(element.type, element.ghost_type)(element.element);
      debug::setDebugLevel(dbl);

      return model.getMaterialIndex(material_name);
    } catch (...) {
      return MaterialSelector::operator()(element);
    }
  }

private:
  const ElementTypeMapArray<std::string> & names;
protected:
  const EmbeddedInterfaceModel & model;
};


__END_AKANTU__

#endif // __AKANTU_EMBEDDED_INTERFACE_MODEL_HH__
