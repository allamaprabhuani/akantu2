/**
 * @file   contact_2d_explicit.cc
 *
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 *
 * @date   Fri Nov 19 14:23:18 2010
 *
 * @brief  Contact class for explicit contact in 2d based on DCR algorithm
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

#include "contact_2d_explicit.hh"
#include "contact_search.hh"
#include "aka_array.hh"

#define PROJ_TOL 1.E-8

__BEGIN_AKANTU__

Contact2dExplicit::Contact2dExplicit(const SolidMechanicsModel & model,
				     const ContactType & type,
				     const ContactID & id,
				     const MemoryID & memory_id) :
  Contact(model, type, id, memory_id), friction_coefficient(0.),
  coefficient_of_restitution(0.) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model.getFEM().getMesh().getSpatialDimension();
  if(spatial_dimension != 2)
    AKANTU_DEBUG_ERROR("Wrong ContactType for contact in 2d!");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Contact2dExplicit::~Contact2dExplicit()
{
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}


// void Contact2dExplicit::defineSurfaces() {
//   AKANTU_DEBUG_IN();

//   Mesh & mesh = model.fem->getMesh();

//   // const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
//   // UInt nb_types = type_list.size();

//   /// Declare vectors
//   // Array<UInt> node_offset;
//   // Array<UInt> node_to_elem
//     ;
//   UInt * facet_conn_val = NULL;
//   UInt nb_facet_elements;

//   /// Consider only linear facet elements
//   // for(it = type_list.begin(); it != type_list.end(); ++it) {
//   ElementType facet_type;
//     // if(type == _line_1) {
//   nb_facet_elements = mesh.getNbElement(facet_type);
//   facet_conn_val = mesh.getConnectivity(facet_type).values;
//       // buildNode2ElementsByElementType(mesh, facet_type, node_offset, node_to_elem);

//   /// Ricominciamo usando le funzioni di più alto livello in mesh

//   const UInt nb_surfaces = getNbSurfaces(facet_type);
//   nb_facet_elements = mesh.getNbElement(facet_type);
//   const Array<UInt> * elem_to_surf_val = getSurfaceValues(facet_type).values;

//   for (UInt i = 0; i < nb_surfaces; ++i)
//     addMasterSurface(i);

//   /// OK ora riordinare le superfici con i rispettivi nodi ed elementi (in linea)
//   surf_to_node_offset.resize(nb_surfaces+1);

//   for (UInt i = 0; i < nb_facet_elements; ++i)
//     surf_to_elem_offset.val[elem_to_surf_val[i]]++;

//   ///  rearrange surf_to_elem_offset
//   for (UInt i = 1; i < nb_surfaces; ++i) surf_to_elem_offset.values[i] += surf_to_elem_offset.values[i-1];
//   for (UInt i = nb_surfaces; i > 0; --i) surf_to_elem_offset.values[i]  = surf_to_elem_offset[i-1];
//   surf_to_elem_offset.values[0] = 0;

//   for (UInt i = 0; i < nb_facet_elements; ++i)
//     surf_to_elem[elem_to_surf_val[]]


//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */
void Contact2dExplicit::solveContact() {
  AKANTU_DEBUG_IN();

  /// Loop over master surfaces
  std::vector<Surface>::iterator it;
  std::vector<Surface> master_surfaces = getMasterSurfaces();
  for (it = master_surfaces.begin(); it != master_surfaces.end(); ++it) {

    /// Get penetration list (find node to segment penetration)
    std::stringstream sstr; sstr << id << ":penetration_list";
    PenetrationList pen_list(sstr.str());
    contact_search->findPenetration(*it, pen_list);

    if(pen_list.penetrating_nodes.getSize() > 0) {
      Array<Real> vel_norm(0, 2);
      Array<Real> vel_fric(0, 2);
      Array<UInt> nodes_index(0, 2);

      /// Remove node to segment penetrations
      projectNodesOnSegments(pen_list, nodes_index);

      /// Compute normal velocities of impacting nodes according to the normal linear momenta
      computeNormalVelocities(pen_list, nodes_index, vel_norm);

      /// Compute friction velocities
      computeFrictionVelocities(pen_list, nodes_index, vel_norm, vel_fric);

      /// Update post-impact velocities
      updatePostImpactVelocities(pen_list, nodes_index, vel_norm, vel_fric);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Contact2dExplicit::projectNodesOnSegments(PenetrationList & pen_list, Array<UInt> & nodes_index) {
  AKANTU_DEBUG_IN();

  UInt dim = 2;
  ElementType el_type = _segment_2;

  const SolidMechanicsModel & model = getModel();
  Real * inc_val = model.getIncrement().values;
  Real * disp_val = model.getDisplacement().values;
  Real * pos_val = model.getCurrentPosition().values;

  Array<UInt> & pen_nodes = pen_list.penetrating_nodes;
  Array<UInt> & pen_facets_off = pen_list.penetrated_facets_offset(el_type, _not_ghost);
  Array<UInt> & pen_facets = pen_list.penetrated_facets(el_type, _not_ghost);
  Array<Real> & projected_positions = pen_list.projected_positions(el_type, _not_ghost);
  Array<Real> & gaps = pen_list.gaps(el_type, _not_ghost);
  Array<Real> & facets_normals = pen_list.facets_normals(el_type, _not_ghost);

  Real * delta = new Real[dim*pen_nodes.getSize()];

  UInt * conn_val = model.getFEM().getMesh().getConnectivity(el_type, _not_ghost).values;
  UInt elem_nodes = Mesh::getNbNodesPerElement(el_type);

  nodes_index.resize(3*pen_nodes.getSize());

  UInt * index = nodes_index.values;

  for (UInt n = 0; n < pen_nodes.getSize(); ++n) {
    const UInt i_node = pen_nodes(n);
    for (UInt el = pen_facets_off(n); el < pen_facets_off(n+1); ++el) {
      /// Project node on segment according to its projection
      Real proj = projected_positions(n);

      index[n*3]   = conn_val[elem_nodes*pen_facets(n)];
      index[n*3+1] = conn_val[elem_nodes*pen_facets(n) + 1];
      index[n*3+2] = i_node;

      if(proj > PROJ_TOL && proj < 1.-PROJ_TOL) {
	delta[2*n]   = -gaps(n)*facets_normals(el, 0);
	delta[2*n+1] = -gaps(n)*facets_normals(el, 1);
      }

      else if(proj <= PROJ_TOL) {
	delta[2*n]   = pos_val[dim*index[n*3]]   - pos_val[dim*i_node];
	delta[2*n+1] = pos_val[dim*index[n*3]+1] - pos_val[dim*i_node+1];
	if(proj < -PROJ_TOL) {
	  Real modulus = sqrt(delta[2*n]*delta[2*n]+delta[2*n+1]*delta[2*n+1]);
	  facets_normals(el, 0) = delta[2*n]/modulus;
	  facets_normals(el, 1) = delta[2*n+1]/modulus;
	}
	proj = 0.;
      }

      else { /* proj >= 1.- PROJ_TOL */
	delta[2*n]   = pos_val[dim*index[n*3+1]]   - pos_val[dim*i_node];
	delta[2*n+1] = pos_val[dim*index[n*3+1]+1] - pos_val[dim*i_node+1];
	if(proj > 1.+PROJ_TOL) {
	  Real modulus = sqrt(delta[2*n]*delta[2*n]+delta[2*n+1]*delta[2*n+1]);
	  facets_normals(el, 0) = delta[2*n]/modulus;
	  facets_normals(el, 0) = delta[2*n+1]/modulus;
	  // Note from Nicolas  : I put the same  index during the modification,
	  // but shouldn't it be facets_normals(el, 0) on the second line ?
	}
	proj = 1.;
      }
    }
  }

  /// Update displacement and increment of the projected nodes
  for (UInt n = 0; n < pen_nodes.getSize(); ++n) {
    //    UInt i_node = pen_nodes.values[n];

    pos_val[dim*index[n*3+2]]    += delta[2*n];
    pos_val[dim*index[n*3+2]+1]  += delta[2*n+1];
    disp_val[dim*index[n*3+2]]   += delta[2*n];
    disp_val[dim*index[n*3+2]+1] += delta[2*n+1];
    inc_val[dim*index[n*3+2]]    += delta[2*n];
    inc_val[dim*index[n*3+2]+1]  += delta[2*n+1];
  }

  delete [] delta;

  AKANTU_DEBUG_OUT();
}

void Contact2dExplicit::computeNormalVelocities(PenetrationList & pen_list,
						Array<UInt> & nodes_index,
						Array<Real> & vel_norm) {
  AKANTU_DEBUG_IN();

  const UInt dim = 2;
  const SolidMechanicsModel & model = getModel();
  const ElementType el_type = _segment_2;
  //  UInt * conn_val = model.getFEM().getMesh().getConnectivity(el_type, _not_ghost).values;
  //  const UInt elem_nodes = Mesh::getNbNodesPerElement(el_type);

  Array<Real> & projected_positions = pen_list.projected_positions(el_type, _not_ghost);
  Array<Real> & facets_normals = pen_list.facets_normals(el_type, _not_ghost);

  Real * vel_val = model.getVelocity().values;
  Real * mass_val = model.getMass().values;
  Real * pos_val = model.getCurrentPosition().values;

  vel_norm.resize(6*pen_list.penetrating_nodes.getSize());
  Real * v_n = vel_norm.values;
  UInt * index = nodes_index.values;
  Real dg[6];

  /// Loop over projected nodes
  for (UInt n = 0; n < pen_list.penetrating_nodes.getSize(); ++n) {

    // for (UInt el = pen_list.penetrated_facets_offset[el_type]->values[n]; el < pen_list.penetrated_facets_offset[el_type]->values[n+1]; ++el) {

    Real proj = projected_positions(n);

      /// Fill gap derivative
      if (proj < PROJ_TOL) {
	dg[0] =  facets_normals(n, 0);
	dg[1] =  facets_normals(n, 1);
	dg[2] = 0.;
	dg[3] = 0.;
	dg[4] = -facets_normals(n, 0);
	dg[5] = -facets_normals(n, 1);
      }
      else if (proj > 1.-PROJ_TOL) {
	dg[0] = 0.;
	dg[1] = 0.;
	dg[2] =  facets_normals(n, 0);
	dg[3] =  facets_normals(n, 1);
	dg[4] = -facets_normals(n, 0);
	dg[5] = -facets_normals(n, 1);
      }
      else {
	dg[0] =  pos_val[index[n*3+2]*dim+1] - pos_val[index[n*3+1]*dim+1];
	dg[1] =  pos_val[index[n*3+1]*dim]   - pos_val[index[n*3+2]*dim];
	dg[2] =  pos_val[index[n*3]*dim+1]   - pos_val[index[n*3+2]*dim+1];
	dg[3] =  pos_val[index[n*3+2]*dim]   - pos_val[index[n*3]*dim];
	dg[4] =  pos_val[index[n*3+1]*dim+1] - pos_val[index[n*3]*dim+1];
	dg[5] =  pos_val[index[n*3]*dim]     - pos_val[index[n*3+1]*dim];
      }

      /// Compute normal velocities exchanged during impact
      Real temp1 = 0., temp2 = 0.;
      for (UInt i=0; i<3; ++i)
	for (UInt j = 0; j < dim; ++j) {
	  temp1 += dg[i*dim+j]*vel_val[dim*index[3*n+i]+j];
	  temp2 += dg[i*dim+j]*dg[i*dim+j]/mass_val[index[3*n+i]];
	}
      temp1 /= temp2;

      for (UInt i=0; i<3; ++i)
	for (UInt j = 0; j < dim; ++j)
	  v_n[n*6+i*dim+j] = temp1*dg[i*dim+j]/mass_val[index[3*n+i]];

    // }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
 void Contact2dExplicit::computeFrictionVelocities(PenetrationList & pen_list,
						  Array<UInt> & nodes_index,
						  Array<Real> & vel_norm,
						  Array<Real> & vel_fric) {

  AKANTU_DEBUG_IN();

  vel_fric.resize(6*pen_list.penetrating_nodes.getSize());
  if (getFrictionCoefficient() == 0) {
    memset(vel_fric.values, 0, (6*pen_list.penetrating_nodes.getSize())*sizeof(Real));
    return;
  }

  const ElementType el_type = _segment_2;
  const SolidMechanicsModel & model = getModel();
  //  UInt * conn_val = model.getFEM().getMesh().getConnectivity(el_type, _not_ghost).values;
  //  const UInt elem_nodes = Mesh::getNbNodesPerElement(el_type);
  const UInt dim = 2;

  Array<Real> & projected_positions = pen_list.projected_positions(el_type, _not_ghost);

  Real * vel_val = model.getVelocity().values;
  Real * mass_val = model.getMass().values;

  Real * v_n = vel_norm.values;
  Real * v_f = vel_fric.values;
  UInt * index = nodes_index.values;

  /// Loop over projected nodes
  for (UInt n = 0; n < pen_list.penetrating_nodes.getSize(); ++n) {

    // for (UInt el = pen_list.penetrated_facets_offset[el_type]->values[n]; el < pen_list.penetrated_facets_offset[el_type]->values[n+1]; ++el) {


    const Real h[3] = {1.-projected_positions(n),
		       projected_positions(n), -1.};
    Real temp[3] = {0., 0., 0.};

    for (UInt i=0; i<3; ++i) {
      temp[0] += h[i]*vel_val[dim*index[i+3*n]];
      temp[1] += h[i]*vel_val[dim*index[i+3*n]+1];
      temp[2] += h[i]*h[i]/mass_val[index[i+3*n]];
    }

    /// Compute non-fixed components of velocities
    for (UInt i=0; i<2; ++i)
      for (UInt j=0; j<3; ++j)
	v_f[n*6+j*dim+i] = temp[i]*h[j]/(temp[2]*mass_val[index[3*n+j]]);

    /// Compute slide components of velocities
    for (UInt i=0; i<6; ++i)
      v_f[n*6+i] = v_f[n*6+i] - v_n[n*6+i];

    /// Compute final friction velocities */
    memset(temp, 0, 3*sizeof(Real));
    for (UInt i=0; i<3; ++i)
      for (UInt j = 0; j < dim; ++j) {
	temp[0] += v_n[n*6+i*dim+j]*v_n[n*6+i*dim+j]/mass_val[index[3*n+i]];
	temp[1] += v_f[n*6+i*dim+j]*v_f[n*6+i*dim+j]/mass_val[index[3*n+i]];
      }
    Real mu = sqrt(temp[0]/temp[1])*getFrictionCoefficient();

    if (mu < 1.)
      for (UInt i=0; i<6; ++i)
	v_f[i] *= mu;
    // }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Contact2dExplicit::updatePostImpactVelocities(PenetrationList & pen_list,
						   Array<UInt> & nodes_index,
						   Array<Real> & vel_norm,
						   Array<Real> & vel_fric) {

  AKANTU_DEBUG_IN();

  const UInt dim = 2;

  Real * v = model.getVelocity().values;
  Real * v_n = vel_norm.values;
  Real * v_f = vel_fric.values;
  UInt * index = nodes_index.values;

  /// Loop over projected nodes
  for (UInt n = 0; n < pen_list.penetrating_nodes.getSize(); ++n)
    for (UInt i = 0; i < 3; ++i)
      for (UInt j = 0; j < dim; ++j)
	v[dim*index[3*n+i]+j] -= ((1.+coefficient_of_restitution)*v_n[n*6+dim*i+j] + v_f[n*6+dim*i+j]);

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
