/**
 * @file   test_data_accessor.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Tue Apr  2 14:25:20 2013
 *
 * @brief  Data Accessor class for testing
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

#include "aka_common.hh"
#include "mesh.hh"
#include "data_accessor.hh"


/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class TestAccessor : public DataAccessor {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  inline TestAccessor(const Mesh & mesh, const ByElementTypeReal & barycenters);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Barycenter, barycenters, Real);

  /* ------------------------------------------------------------------------ */
  /* Ghost Synchronizer inherited members                                     */
  /* ------------------------------------------------------------------------ */
protected:
  inline UInt getNbDataForElements(const Array<Element> & elements,
				   SynchronizationTag tag) const;
  inline void packElementData(CommunicationBuffer & buffer,
			      const Array<Element> & elements,
			      SynchronizationTag tag) const;
  inline void unpackElementData(CommunicationBuffer & buffer,
				const Array<Element> & elements,
				SynchronizationTag tag);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  const ByElementTypeReal & barycenters;
  const Mesh & mesh;
};


/* -------------------------------------------------------------------------- */
/* TestSynchronizer implementation                                            */
/* -------------------------------------------------------------------------- */
inline TestAccessor::TestAccessor(const Mesh & mesh,
				  const ByElementTypeReal & barycenters)
  : barycenters(barycenters), mesh(mesh) { }

inline UInt TestAccessor::getNbDataForElements(const Array<Element> & elements,
					       __attribute__ ((unused)) SynchronizationTag tag) const {
  if(elements.getSize())
    // return Mesh::getSpatialDimension(elements(0).type) * sizeof(Real) * elements.getSize();
    return mesh.getSpatialDimension() * sizeof(Real) * elements.getSize();
  else
    return 0;
}

inline void TestAccessor::packElementData(CommunicationBuffer & buffer,
					  const Array<Element> & elements,
					  __attribute__ ((unused)) SynchronizationTag tag) const {
  UInt spatial_dimension = mesh.getSpatialDimension();
  Array<Element>::const_iterator<Element> bit  = elements.begin();
  Array<Element>::const_iterator<Element> bend = elements.end();
  for (; bit != bend; ++bit) {
    const Element & element = *bit;

    Vector<Real> bary(this->barycenters(element.type, element.ghost_type).storage()
		      + element.element * spatial_dimension,
		      spatial_dimension);
    buffer << bary;
  }
}

inline void TestAccessor::unpackElementData(CommunicationBuffer & buffer,
					    const Array<Element> & elements,
					    __attribute__ ((unused)) SynchronizationTag tag) {
  UInt spatial_dimension = mesh.getSpatialDimension();
  Array<Element>::const_iterator<Element> bit  = elements.begin();
  Array<Element>::const_iterator<Element> bend = elements.end();
  for (; bit != bend; ++bit) {
    const Element & element = *bit;

    Vector<Real> barycenter_loc(this->barycenters(element.type, element.ghost_type).storage()
				+ element.element * spatial_dimension,
				spatial_dimension);

    Vector<Real> bary(spatial_dimension);
    buffer >> bary;
    std::cout << element << barycenter_loc << std::endl;
    Real tolerance = 1e-15;
    for (UInt i = 0; i < spatial_dimension; ++i) {
      if(!(std::abs(bary(i) - barycenter_loc(i)) <= tolerance))
        AKANTU_DEBUG_ERROR("Unpacking an unknown value for the element: "
                           << element
                           << "(barycenter[" << i << "] = " << barycenter_loc(i)
                           << " and buffer[" << i << "] = " << bary(i) << ") - tag: " << tag);
    }
  }
}


__END_AKANTU__
