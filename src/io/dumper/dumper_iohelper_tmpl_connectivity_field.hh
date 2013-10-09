/**
 * @file   dumper_iohelper_tmpl_connectivity_field.hh
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Mon May  6 10:40:29 2013
 *
 * @brief  Connectivity field dumper
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


/* -------------------------------------------------------------------------- */
class DumperIOHelper::cohesive_connectivity_field_iterator : public element_iterator<UInt, UInt, Vector, cohesive_connectivity_field_iterator, false> {
public:
  typedef element_iterator<UInt, UInt,
			   Vector, cohesive_connectivity_field_iterator, false> parent;

  typedef parent::it_type     it_type;
  typedef parent::data_type   data_type;
  typedef parent::return_type return_type;
  typedef parent::field_type  field_type;
  typedef parent::internal_iterator internal_iterator;
public:
  cohesive_connectivity_field_iterator(const field_type & field,
				       __attribute__((unused)) UInt n,
				       const field_type::type_iterator & t_it,
				       const field_type::type_iterator & t_it_end,
				       const internal_iterator & it,
				       ElementType element_type,
				       const GhostType ghost_type = _not_ghost,
				       const ByElementTypeArray<UInt> * filter = NULL,
				       UInt * fit = NULL) :
    parent(field, 0, t_it, t_it_end, it, element_type, ghost_type, filter, fit) {

    write_order[_cohesive_3d_12].push_back(0);
    write_order[_cohesive_3d_12].push_back(1);
    write_order[_cohesive_3d_12].push_back(2);
    write_order[_cohesive_3d_12].push_back(6);
    write_order[_cohesive_3d_12].push_back(7);
    write_order[_cohesive_3d_12].push_back(8);
    write_order[_cohesive_3d_12].push_back(3);
    write_order[_cohesive_3d_12].push_back(4);
    write_order[_cohesive_3d_12].push_back(5);
    write_order[_cohesive_3d_12].push_back(9);
    write_order[_cohesive_3d_12].push_back(10);
    write_order[_cohesive_3d_12].push_back(11);

    write_order[_cohesive_3d_6].push_back(0);
    write_order[_cohesive_3d_6].push_back(1);
    write_order[_cohesive_3d_6].push_back(2);
    write_order[_cohesive_3d_6].push_back(3);
    write_order[_cohesive_3d_6].push_back(4);
    write_order[_cohesive_3d_6].push_back(5);

    write_order[_cohesive_2d_6].push_back(0);
    write_order[_cohesive_2d_6].push_back(2);
    write_order[_cohesive_2d_6].push_back(1);
    write_order[_cohesive_2d_6].push_back(4);
    write_order[_cohesive_2d_6].push_back(5);
    write_order[_cohesive_2d_6].push_back(3);

    write_order[_cohesive_2d_4].push_back(0);
    write_order[_cohesive_2d_4].push_back(1);
    write_order[_cohesive_2d_4].push_back(3);
    write_order[_cohesive_2d_4].push_back(2);

  }

  return_type operator*() {
    ElementType type = *tit;
    const Vector<UInt> & conn = *vit;
    Vector<UInt> new_conn(conn.size());

    for (UInt n = 0; n < conn.size(); ++n)
      new_conn(n) = conn(write_order[type][n]);

    return new_conn;
  }

protected:

  std::map<ElementType, std::vector<UInt> > write_order;

};


/* -------------------------------------------------------------------------- */
class DumperIOHelper::CohesiveConnectivityField : public GenericElementalField<UInt,
									       cohesive_connectivity_field_iterator,
									       Vector,
									       false> {
public:
  typedef cohesive_connectivity_field_iterator iterator;
private:
  typedef GenericElementalField<UInt, iterator, Vector, false> parent;
public:
  /* ------------------------------------------------------------------------ */
  CohesiveConnectivityField(const Mesh & mesh,
			    UInt spatial_dimension = 0,
			    GhostType ghost_type = _not_ghost) :
    parent(mesh.getConnectivities(), 0, spatial_dimension, ghost_type, _ek_cohesive) { }

  virtual void registerToDumper(__attribute__((unused)) const std::string & id,
			   iohelper::Dumper & dumper) {
    dumper.addElemDataField("connectivities", *this);
  }

};
