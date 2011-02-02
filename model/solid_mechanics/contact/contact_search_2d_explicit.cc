/**
 * @file   contact_search_2d_explicit.cc
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @date   Wed Nov  3 15:06:52 2010
 *
 * @brief  Contact 
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique fédérale de Lausanne)
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

#include "contact_search_2d_explicit.hh"
#include "contact.hh"
#include "contact_neighbor_structure.hh"
#include "aka_memory.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
ContactSearch2dExplicit::ContactSearch2dExplicit(Contact & contact,
						 const ContactNeighborStructureType & neighbors_structure_type,
						 const ContactSearchType & type,
						 const ContactSearchID & id) :
  ContactSearch(contact, neighbors_structure_type, type, id) {

  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
ContactSearch2dExplicit::~ContactSearch2dExplicit() {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactSearch2dExplicit::findPenetration(const Surface & master_surface, PenetrationList & pen_list) {
  AKANTU_DEBUG_IN();

  //const ContactNeighborStructureType & neighbor_type = getContactNeighborStructureType();
  const ContactNeighborStructure & structure = getContactNeighborStructure(master_surface);
  const NeighborList & neigh_list = const_cast<ContactNeighborStructure&>(structure).getNeighborList();

  // Real * inc_val = contact.getModel().getIncrement().values;
  // Real * disp_val = contact.getModel().getDisplacement().values;
  // Real * vel_val = contact.getModel().getVelocity().values;
  const Contact & contact = getContact();
  Real * pos_val = contact.getModel().getCurrentPosition().values;

  //contact.getModel().updateCurrentPosition();
  /// compute current_position = initial_position + displacement temp
  // Real * coord = contact.getModel().getFEM().getMesh().getNodes().values;
  // Real * disp = contact.getModel().getDisplacement().values;
  // Real * pos_val = new Real[contact.getModel().getFEM().getMesh().getNodes().getSize()];
  // for (UInt n = 0; n < contact.getModel().getFEM().getMesh().getNodes().getSize(); ++n) 
  //   pos_val[n] = coord[n] + disp[n];

  Mesh & mesh = contact.getModel().getFEM().getMesh();

  ElementType el_type = _segment_2; /* Only linear element at the moment */
  const ContactSearchID & id = getID();

   /// Alloc space for Penetration class members
  std::stringstream sstr_pfo; sstr_pfo << id << ":penetrated_facets_offset:" << el_type;
  pen_list.penetrated_facets_offset[el_type] = new Vector<UInt>(0, 1, sstr_pfo.str());
  std::stringstream sstr_pf; sstr_pf << id << ":penetrated_facet:" << el_type;
  pen_list.penetrated_facets[el_type] = new Vector<UInt>(0 ,1, sstr_pf.str());
  std::stringstream sstr_fn; sstr_fn << id << ":facet_normals:" << el_type;
  pen_list.facets_normals[el_type] = new Vector<Real>(0, 2, sstr_fn.str());
  std::stringstream sstr_g; sstr_g << id << ":gaps:" << el_type;
  pen_list.gaps[el_type] = new Vector<Real>(0, 1, sstr_g.str());
  std::stringstream sstr_pp; sstr_pp << id << ":projected_positions:" << el_type;
  pen_list.projected_positions[el_type] = new Vector<Real>(0, 1, sstr_pp.str());

  //pen_list.nb_nodes = 0;

  UInt * facets_off_val = neigh_list.facets_offset[el_type]->values;
  UInt * facets_val = neigh_list.facets[el_type]->values;

  UInt * node_to_el_off_val = contact.getNodeToElementsOffset(el_type).values;
  UInt * node_to_el_val = contact.getNodeToElements(el_type).values;

  UInt * conn_val = mesh.getConnectivity(el_type).values;
  UInt elem_nodes = Mesh::getNbNodesPerElement(el_type);

  /// Loop over the impactor nodes to detect possible penetrations
  for (UInt i = 0; i < neigh_list.impactor_nodes.getSize(); ++i) {

    UInt i_node = neigh_list.impactor_nodes.values[i];
    Real *x3 = &pos_val[i_node*2];

    /// Loop over elements nearby the impactor node    
    for (UInt i_el = facets_off_val[i]; i_el < facets_off_val[i+1]; ++i_el) {

      UInt facet = facets_val[i_el];
      Real * x1 = &pos_val[2*conn_val[facet*elem_nodes]];
      Real * x2 = &pos_val[2*conn_val[facet*elem_nodes+1]];

      Real vec_surf[2];
      Real vec_normal[2];
      Real vec_dist[2];
      Real length = Math::distance_2d(x1, x2);
      Math::vector_2d(x1, x2, vec_surf);
      Math::vector_2d(x1, x3, vec_dist);
      Math::normal2(vec_surf, vec_normal);
      
      /* Projection normalized over length*/
      Real projection = Math::vectorDot2(vec_surf, vec_dist)/(length*length); 
      Real gap = Math::vectorDot2(vec_dist, vec_normal);

      bool find_proj = false;

      /// Penetration has occurred
      if(gap < -PEN_TOL) {

	InterType test_pen = Detect_Intersection(conn_val[facet*elem_nodes], conn_val[facet*elem_nodes+1],
			    i_node, vec_surf, vec_dist, gap);

	/// Node has intersected segment
	if(test_pen != _no) {

	  Real proj = Math::vectorDot2(vec_surf, vec_dist)/(length*length);

	  UInt c_facet = facet;

	  /// Projection on neighbor facet which shares to node1
	  if(test_pen == _node_1 || (test_pen==_yes && proj < 0.-PROJ_TOL)) {

	    // Find index of neighbor facet
	    for (UInt i = node_to_el_off_val[conn_val[facet*elem_nodes]]; i < node_to_el_off_val[conn_val[facet*elem_nodes]+1]; ++i)
	      if(node_to_el_val[i] != facet /* &&  on same surface? */) {
		c_facet = node_to_el_val[i];
		break;
	      }

	    if(c_facet != facet) {
	      find_proj = checkProjectionAdjacentFacet(pen_list, facet, c_facet, i_node, proj, el_type);
	      if(find_proj == true)
		break;
	    }
	    if (proj < 0.-PROJ_TOL && test_pen == _node_1)
	      break;
	  }

	  /// Projection on neighbor facet which shares to node2
	  if(test_pen == _node_2 || (test_pen==_yes && proj > 1.+PROJ_TOL)) {

	    // Find index of neighbor facet
	    for (UInt i = node_to_el_off_val[conn_val[facet*elem_nodes]]; i < node_to_el_off_val[conn_val[facet*elem_nodes]+1]; ++i)
	      if(node_to_el_val[i] != facet /* &&  on same surface? */) {
		c_facet = node_to_el_val[i];
		break;
	      }

	    if(c_facet != facet) {
	      find_proj = checkProjectionAdjacentFacet(pen_list, facet, c_facet, i_node, proj, el_type);
	      if(find_proj == true)
		break;
	    }
	    if (proj > 0.+PROJ_TOL && test_pen == _node_2)
	      break;
	  }

	  /// Save data for projection on this facet -> Check if projection outside in solve?
	  pen_list.penetrated_facets_offset[el_type]->push_back(pen_list.penetrating_nodes.getSize());
	  pen_list.penetrating_nodes.push_back(i_node);	  
	  pen_list.penetrated_facets[el_type]->push_back(facet);
	  pen_list.facets_normals[el_type]->push_back(vec_normal); /* Correct ? */
	  pen_list.gaps[el_type]->push_back(gap);
	  pen_list.projected_positions[el_type]->push_back(proj);
	  //pen_list.penetrating_nodes.getSize()++;
	}
      }
    }
  }
  pen_list.penetrated_facets_offset[el_type]->push_back(pen_list.penetrating_nodes.getSize());

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
InterType ContactSearch2dExplicit::Detect_Intersection(UInt node1, UInt node2, UInt node3,
					       Real *vec_surf, Real *vec_dist, Real gap)
{
AKANTU_DEBUG_IN();

  const Real eps = 1.e-12, /* Tolerance to set discriminant equal to zero */
    eps2 = 1e6, /* */
    eps3 = 1e-6; /* Tolerance for space intersection "l" outside segment (tolerance for projection) */
    // eps4=1e-2;
  /* Upper tolerance for time: penetration may occurred in previous time_step but had been not detected */ 
  const Real eps4 = 1.0001*PEN_TOL/(-gap-PEN_TOL);

  Real a[2], b[2], c[2], d[2], den, delta, t[2]={-1.,-1.}, l[2]={-1.,-1.}, k, j;

    Real * inc_val = getContact().getModel().getIncrement().values;
  
  /* Initialize vectors of Equation to solve : 0 = a + b*l +c*t + d*t*l */
  for(UInt i=0; i<2; i++) {
    a[i] = -vec_dist[i];
    b[i] = vec_surf[i];
    c[i] = inc_val[2*node3+i]-inc_val[2*node1+i];
    d[i] = inc_val[2*node1+i]-inc_val[2*node2+i];
  }

  /* Compute therm of quadratic equation: t²+k*t+j=0 */

  k = (a[1]*d[0]+c[1]*b[0]-c[0]*b[1]-d[1]*a[0])/(c[1]*d[0]-d[1]*c[0]);

  j = (a[1]*b[0]-b[1]*a[0])/(c[1]*d[0]-d[1]*c[0]);


  if(isnan(k)!=true && isnan(j)!=true) {

    if(k < 1.e12 && j < 1.e12) { /* Equation quadratic */ /* If quadratic therm smaller than eps*other therm exclude them */

    /* Compute discriminant */
    delta = k*k-4.*j;

    /* Compute solution of quadratic equation */
    do {

      /* Ceck if disc<0 */
      if(delta < 0.) {
	if(delta < -eps)
	  return _no;
	else  // delta close to the machine precision -> only one solution
	  /* delta = 0. */
	  t[0] = -k/2.;
	  break;
      }

      /* Discriminant bigger-equal than zero */
      t[0] = (-k-getSign(k)*sqrt(delta))/2.;
      t[1] = (-k+getSign(k)*sqrt(delta))/2.;
      if(abs(t[1]) < 2.*eps)
	t[1] = j/t[0];

    } while(0);

    /* Ceck solution */
    for(UInt i=0; i<2; i++)
      if(t[i]>= 0.-eps3 && t[i]<=1.+eps4) { // Within time step
	l[i] = -(a[0]+c[0]*t[i])/(b[0]+t[i]*d[0]);
	if(l[i]>= 0.-eps3 && l[i]<= 1.+eps3) {

	  if(/*t[i]<1.-eps &&*/ l[i]>0.+1.e-3 && l[i]<1.-1.e-3)
	    return _yes;

	  else if(l[i]<0.5) { /* Node 3 close at the beginning to node 1*/
	    return _node_1;
	  }
	  else { /* Node 3 close at the beginning to node 2 */
	    return _node_2;
	  }
	}
      }

    return _no;
    }
    else { /* New Equation: t²+k*t+j=0 -> k*t+j=0*/

      t[0]=-(a[1]*b[0]-b[1]*a[0])/(a[1]*d[0]+c[1]*b[0]-c[0]*b[1]-d[1]*a[0]);

      if(t[0]>= 0.-eps3 && t[0]<=1.+eps4) { // Within time step
	l[0] = -(a[0]+c[0]*t[0])/(b[0]+t[0]*d[0]);
	if(l[0]>= 0.-eps3 && l[0]<= 1.+eps3) {
	  if(/*t[i]<1.-eps &&*/ l[0]>0.+1.e-3 && l[0]<1.-1.e-3)
	    return _yes;
	  
	  else if(l[0]<0.5)  /* Node 3 close at the beginning to node 1*/
	    return _node_1;
	    
	  else  /* Node 3 close at the beginning to node 2 */
	    return _node_2;
	}
      }
      return _no;
    }
  }

  else if(abs(c[1]/d[0])>eps2 && abs(c[0]/d[1])>eps2) { /* neglect array d */

    den = b[0]*c[1]-b[1]*c[0];
    t[0] = (a[0]*b[1]-b[0]*a[1])/den;
    if(isnan(t[0])) /* Motion parallel */
      return _no; /* ? */

    /* Check solution */
    if(t[0]>= 0.-eps3 && t[0]<=1.+eps4) { // Possible Intersection Within time step
      l[0] = -(c[1]*a[0]-c[0]*a[1])/den;
      if(l[0]>= 0.-eps3 && l[0]<= 1.+eps3)
	return _yes;
    }
    return _no;
  }

  else if(abs(d[0]/c[1])>eps2 && abs(d[1]/c[0])>eps2) { /* neglect array c */

   t[0] = -(-a[1]*b[0]+b[1]*a[0])/(d[1]*a[0]-d[0]*a[1]);
   if(isnan(t[0]))
     return _no; /* ? */

   /* Ceck solution */
   if(t[0]>= 0.-eps3 && t[0]<=1.+eps4) { // Possible Intersection Within time step
     l[0] = -a[0]/(b[0]+t[0]*d[0]);
     if(l[0]>= 0.-eps3 && l[0]<= 1.+eps3)
       return _yes;
   }
   return _no;
  }

  /* Neglect array c and d*/
  return _no; /* undefined (either no motion or parallel motion)*/

AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
bool ContactSearch2dExplicit::checkProjectionAdjacentFacet(PenetrationList & pen_list, UInt facet, UInt c_facet, UInt i_node, Real old_proj, ElementType el_type)
 {
   AKANTU_DEBUG_IN();
   
   UInt * conn_val = contact.getModel().getFEM().getMesh().getConnectivity(el_type).values;
   UInt elem_nodes = Mesh::getNbNodesPerElement(el_type);

   UInt node1 = conn_val[c_facet*elem_nodes];
   UInt node2 = conn_val[c_facet*elem_nodes+1];

   Real * pos_val = contact.getModel().getCurrentPosition().values;
   Real * x1 = &pos_val[2*node1];
   Real * x2 = &pos_val[2*node2];
   Real * x3 = &pos_val[2*i_node]; 

   Real vec_surf[2];
   Real vec_normal[2];
   Real vec_dist[2];
   Real length = Math::distance_2d(x1, x2);
   Math::vector_2d(x1, x2, vec_surf);
   Math::vector_2d(x1, x3, vec_dist);
   Math::normal2(vec_surf, vec_normal);
 
   Real gap = Math::vectorDot2(vec_dist, vec_normal);

   if(gap < 0.-PEN_TOL) {
     Real proj = Math::vectorDot2(vec_surf, vec_dist)/(length*length);
     if(proj >= 0. && proj <= 1.) { /* Project on c_facet or node */
 
       /// Save data on penetration list
       if(old_proj < 0. || old_proj > 1.) { /* (project on adjacent segment) */
	 
	 pen_list.penetrated_facets_offset[el_type]->push_back(pen_list.penetrating_nodes.getSize());
	 pen_list.penetrating_nodes.push_back(i_node);
	 pen_list.penetrated_facets[el_type]->push_back(c_facet);
	 pen_list.facets_normals[el_type]->push_back(vec_normal); /* Correct ? */
	 pen_list.gaps[el_type]->push_back(gap);
	 pen_list.projected_positions[el_type]->push_back(proj);
	 //pen_list.penetrating_nodes.getSize()++;
	 return true;
       }

       InterType test_pen = Detect_Intersection(node1, node2, i_node, vec_surf, vec_dist, gap);
       if(test_pen == _no)
	 return false;
 
       /* project on node (compute new normal) */
       Real new_normal[2] = {0.,0.};
       if(old_proj < 0.5) {
	 Real * x4 = &pos_val[2*conn_val[elem_nodes*facet+1]];
	 Math::vector_2d(x1, x4, vec_surf);
	 Math::normal2(vec_surf, vec_normal);
	 // new_normal[0] = -x3[0]+x2[0];
	 // new_normal[1] = -x3[1]+x2[1];
	 pen_list.projected_positions[el_type]->push_back(0.);
	 gap = Math::distance_2d(x3, x2);
       }
       else {
	 Real * x4 = &pos_val[2*conn_val[elem_nodes*facet]];
	 Math::vector_2d(x4, x2, vec_surf);
	 Math::normal2(vec_surf, vec_normal);
	 // new_normal[0] = -x3[0]+x1[0];
	 // new_normal[1] = -x3[1]+x1[1];
	 pen_list.projected_positions[el_type]->push_back(1.);
	 gap = Math::distance_2d(x3, x1);
       }
       // Real new_gap = sqrt(new_normal[0]*new_normal[0]+new_normal[1]*new_normal[1]);
       pen_list.penetrated_facets_offset[el_type]->push_back(pen_list.penetrating_nodes.getSize());
       pen_list.penetrating_nodes.push_back(i_node);
       pen_list.penetrated_facets[el_type]->push_back(facet);
       pen_list.facets_normals[el_type]->push_back(vec_normal);
       pen_list.gaps[el_type]->push_back(gap);
       //pen_list.penetrating_nodes.getSize()++;
       return true;
     }
   }
   return false;
   AKANTU_DEBUG_OUT();
 }


__END_AKANTU__
