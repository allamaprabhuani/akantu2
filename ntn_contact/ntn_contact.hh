/**
 * @file   ntn_contact.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon Feb 20 15:13:23 2012
 *
 * @brief  contact for node to node discretization
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

#ifndef __AST_NTN_CONTACT_HH__
#define __AST_NTN_CONTACT_HH__

/* -------------------------------------------------------------------------- */
// simtools
#include "ntn_base_contact.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
using namespace akantu;

/* -------------------------------------------------------------------------- */
class NTNContact : public NTNBaseContact {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  NTNContact(SolidMechanicsModel & model,
	     const ContactID & id = "contact",
	     const MemoryID & memory_id = 0);
  virtual ~NTNContact() {};
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// add surface pair and pair nodes according to the surface normal
  void addSurfacePair(const Surface & slave, 
		      const Surface & master, 
		      UInt surface_normal_dir);

  /// fills the pairs vector with interface node pairs (*,0)=slaves, (*,1)=masters
  static void pairInterfaceNodes(const SubBoundary & slave_boundary, 
				 const SubBoundary & master_boundary,
				 UInt surface_normal_dir,
				 const Mesh & mesh,
				 Array<UInt> & pairs);

  // add node pairs from a list with pairs(*,0)=slaves and pairs(*,1)=masters
  void addNodePairs(const Array<UInt> & pairs);

  /// add node pair
  virtual void addSplitNode(UInt slave, UInt master);

  /// update (compute the normals on the master nodes)
  virtual void updateNormals();

  /// update the lumped boundary B matrix
  virtual void updateLumpedBoundary();

  /// update the impedance matrix
  virtual void updateImpedance();

  /// impose the normal contact force
  virtual void applyContactPressure();

  /// dump restart file
  virtual void dumpRestart(const std::string & file_name) const;

  /// read restart file
  virtual void readRestart(const std::string & file_name);

  /// compute the normal gap
  virtual void computeNormalGap(Array<Real> & gap) const {
    this->computeRelativeNormalField(this->model.getCurrentPosition(),
				     gap);
  };

  /// compute relative normal field (only value that has to be multiplied with the normal)
  /// relative to master nodes
  virtual void computeRelativeNormalField(const Array<Real> & field,
					  Array<Real> & rel_normal_field) const;

  /// compute relative tangential field (complet array)
  /// relative to master nodes
  virtual void computeRelativeTangentialField(const Array<Real> & field,
					      Array<Real> & rel_tang_field) const;
  
  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

protected:
  /// synchronize arrays
  virtual void syncArrays(SyncChoice sync_choice);

  /* ------------------------------------------------------------------------ */
  /* Dumpable                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  virtual void addDumpFieldToDumper(const std::string & dumper_name,
				    const std::string & field_id);
  //  virtual void addDumpFieldVector(const std::string & field_id);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Masters,                               masters, const SynchronizedArray<UInt> &)
  AKANTU_GET_MACRO(LumpedBoundaryMasters, lumped_boundary_masters, const SynchronizedArray<Real> &)

  /// get interface node pairs (*,0) are slaves, (*,1) are masters
  void getNodePairs(Array<UInt> & pairs) const;

  /// get index of node in either slaves or masters array
  /// if node is in neither of them, return -1
  Int getNodeIndex(UInt node) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// array of master nodes
  SynchronizedArray<UInt> masters;
  /// lumped boundary of master nodes
  SynchronizedArray<Real> lumped_boundary_masters;

  // element list for dump and lumped_boundary
  ByElementTypeArray<UInt> master_elements;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "ntn_contact_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const NTNContact & _this)
{
  _this.printself(stream);
  return stream;
}

__END_SIMTOOLS__

#endif /* __AST_NTN_CONTACT_HH__ */
