/**
 * @file   mesh_io_diana.hh
 * @author Alodie Schneuwly <alodie.schneuwly@epfl.ch>
 * @date   Thu Mar 10 11:24:27 2011
 *
 * @brief
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

#ifndef __AKANTU_MESH_IO_DIANA_HH__
#define __AKANTU_MESH_IO_DIANA_HH__

/* -------------------------------------------------------------------------- */
#include "mesh_io.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class MeshIODiana : public MeshIO {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MeshIODiana();
  virtual ~MeshIODiana();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// read a mesh from the file
  virtual void read(const std::string & filename, Mesh&mesh);

  /// write a mesh to a file
  virtual void write(const std::string & filename, const Mesh & mesh);

private:
  std::string readCoordinates(std::ifstream & infile,
			      Mesh & mesh,
			      UInt & first_node_number);

  std::string readElements(std::ifstream & infile,
			   Mesh & mesh,
			   std::vector<Element> & global_to_local_index,
			   UInt first_node_number);

  std::string readGroups(std::ifstream & infile, Mesh & mesh,
			 std::vector<Element> & global_to_local_index,
			 UInt first_node_number);

  std::string readConnectivity(std::ifstream & infile,
			       Mesh & mesh,
			       std::vector<Element> & global_to_local_index,
			       UInt first_node_number);

  std::string readMaterialElement(std::ifstream & infile,
				  Mesh & mesh,
				  std::vector<Element> & global_to_local_index);

  std::string readMaterial(std::ifstream & infile,
			   Mesh & mesh,
			   const std::string & filename);

  UInt readInterval(std::stringstream & line,
		    std::set<UInt> & interval);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  const Vector<UInt> & getNodeGroup(const std::string & group_name) const {
    std::map<std::string, Vector<UInt> *>::const_iterator it = node_groups.find(group_name);
    AKANTU_DEBUG_ASSERT(it != node_groups.end(), "There is no nodes group named : " << group_name);
    return *it->second;
  }

  const std::vector<Element> & getElementGroup(const std::string & group_name) const {
    std::map<std::string, std::vector<Element> *>::const_iterator it = element_groups.find(group_name);
    AKANTU_DEBUG_ASSERT(it != element_groups.end(), "There is no elements group named : " << group_name);
    return *it->second;
  }

  std::vector<std::string> getNodeGroupsNames() const {
    std::vector<std::string> names;
    std::map<std::string, Vector<UInt> *>::const_iterator it;
    for(it = node_groups.begin(); it != node_groups.end(); ++it)
      names.push_back(it->first);

    return names;
  }

  std::vector<std::string> getElementGroupsNames() const {
    std::vector<std::string> names;
    std::map<std::string, std::vector<Element> *>::const_iterator it;
    for(it = element_groups.begin(); it != element_groups.end(); ++it)
      names.push_back(it->first);

    return names;
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  std::map<std::string, ElementType> _diana_to_akantu_element_types;

  std::map<std::string, std::string> _diana_to_akantu_mat_prop;


  std::map<std::string, Vector<UInt> *> node_groups;
  std::map<std::string, std::vector<Element> *> element_groups;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "mesh_io_diana_inline_impl.cc"
// /// standard output stream operator
// inline std::ostream & operator <<(std::ostream & stream, const MeshIODiana & _this)
// {
//   _this.printself(stream);
//   return stream;
// }


__END_AKANTU__

#endif /* __AKANTU_MESH_IO_DIANA_HH__ */
