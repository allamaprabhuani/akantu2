/**
 * @file   solid_mechanics_model_cohesive.hh
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Thu Apr 19 10:34:00 2012
 *
 * @brief  Solid mechanics model for cohesive elements
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "solid_mechanics_model.hh"

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SOLID_MECHANICS_MODEL_COHESIVE_HH__
#define __AKANTU_SOLID_MECHANICS_MODEL_COHESIVE_HH__

__BEGIN_AKANTU__

class SolidMechanicsModelCohesive : public SolidMechanicsModel {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  typedef FEMTemplate< IntegratorCohesive<IntegratorGauss>, ShapeCohesive<ShapeLagrange> > MyFEMCohesiveType;

  SolidMechanicsModelCohesive(Mesh & mesh,
			      UInt spatial_dimension = 0,
			      const ID & id = "solid_mechanics_model_cohesive",
			      const MemoryID & memory_id = 0);

  //  virtual ~SolidMechanicsModelCohesive();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /// function to perform a stress check on each facet and insert
  /// cohesive elements if needed
  void checkCohesiveStress();

  /// function to insert cohesive elements on the selected facets
  void insertCohesiveElements(const Vector<UInt> & facet_insertion);

  /// initialize completely the model
  void initFull(std::string material_file,
                AnalysisMethod method = _explicit_dynamic);

  /// initialize completely the model for extrinsic elements
  void initExtrinsic();

  /// initialize the model
  void initModel();

  /// initialize cohesive material
  void initCohesiveMaterial();

  /// compute reversible energy
  Real getReversibleEnergy();

  /// compute dissipated energy
  Real getDissipatedEnergy();

private:

  /// function to double a given facet and update the list of doubled
  /// nodes
  void doubleFacet(Element & facet, Vector<UInt> & doubled_nodes,
		   Vector<UInt> & doubled_facets);

  /// function to double a subfacet given start and end index for
  /// local facet_to_subfacet vector, and update the list of doubled
  /// nodes
  void doubleSubfacet(const Element & subfacet, UInt start, UInt end,
		      Vector<UInt> & doubled_nodes);

  /// double middle nodes if facets are _segment_3
  void doubleMiddleNode(Vector<UInt> & doubled_nodes,
			const Vector<UInt> & doubled_facets);

  /// function to update nodal parameters for doubled nodes
  void updateDoubledNodes(const Vector<UInt> & doubled_nodes);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /// get the cohesive element type
  AKANTU_GET_MACRO(CohesiveElementType, type_cohesive, ElementType);

  /// get the cohesive index
  AKANTU_GET_MACRO(CohesiveIndex, cohesive_index, UInt);

  /// get the facet type
  AKANTU_GET_MACRO(FacetType, type_facet, ElementType);

  /// get the facet mesh
  AKANTU_GET_MACRO(MeshFacets, mesh_facets, const Mesh &);

  /// get the sigma limit vector for automatic insertion
  AKANTU_GET_MACRO(SigmaLimit, sigma_lim, const Vector<Real> &);
  AKANTU_GET_MACRO_NOT_CONST(SigmaLimit, sigma_lim, Vector<Real> &);

  /// get the facets check vector
  AKANTU_GET_MACRO_NOT_CONST(FacetsCheck, facets_check, Vector<bool> &);

  /// THIS HAS TO BE CHANGED
  AKANTU_GET_MACRO(Tangents, tangents, const Vector<Real> &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// cohesive element type
  ElementType type_cohesive;

  /// facet type
  ElementType type_facet;

  /// cohesive material index in materials vector
  UInt cohesive_index;

  /// mesh containing facets and their data structures
  Mesh mesh_facets;

  /// vector containing a sigma limit for automatic insertion
  Vector<Real> sigma_lim;

  /// vector containing facets in which cohesive elements can be automatically inserted
  Vector<bool> facets_check;

  /// @todo store tangents when normals are computed:
  Vector<Real> tangents;

  /// quadrature points coordinates by element
  ByElementTypeReal quad_elements;

  /// list of facet quadrature points positions by element
  ByElementTypeReal elements_quad_facets;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "solid_mechanics_model_cohesive_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const SolidMechanicsModelCohesive & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__


#endif /* __AKANTU_SOLID_MECHANICS_MODEL_COHESIVE_HH__ */
