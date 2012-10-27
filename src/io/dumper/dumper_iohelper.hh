/**
 * @file   dumper_iohelper.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Oct 24 13:20:16 2012
 *
 * @brief  Define the akantu dumper interface for IOhelper dumpers
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
#ifndef __AKANTU_DUMPER_IOHELPER_HH__
#define __AKANTU_DUMPER_IOHELPER_HH__

#include "aka_common.hh"
#include "aka_types.hh"
#include "mesh.hh"


#include <io_helper.hh>

__BEGIN_AKANTU__

class DumperIOHelper {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  DumperIOHelper();
  virtual ~DumperIOHelper();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  class Field;

  void registerMesh(const Mesh & mesh, UInt spatial_dimension = 0,
		    const GhostType & ghost_type = _not_ghost,
		    const ElementKind & element_kind = _ek_not_defined);
  void registerField(const std::string & field_id, Field * field);
  void unRegisterField(const std::string & field_id);

  void dump();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Dumper, *dumper, iohelper::Dumper &)

  static iohelper::ElemType getIOHelperType(ElementType type);
protected:
  template <ElementType type> static iohelper::ElemType getIOHelperType();

public:
  /* ------------------------------------------------------------------------ */
  /* Field descriptors                                                        */
  /* ------------------------------------------------------------------------ */
  /// Field interface
  class Field {
  public:
    Field() : padding_n(0), padding_m(0) {}
    virtual ~Field() {};
    virtual void registerToDumper(const std::string & id, iohelper::Dumper & dupmer) = 0;
    virtual void setPadding(UInt n, UInt m = 1) {
      padding_n = n;
      padding_m = m;
    }
  protected:
    UInt padding_n, padding_m;
  };


  /* ------------------------------------------------------------------------ */
  template<class T, template<class> class R>
  class iterator_helper;

  /* ------------------------------------------------------------------------ */
  /* Nodal field wrapper */
  template<typename T, template<typename> class return_type = types::Vector>
  class NodalField;

  /* ------------------------------------------------------------------------ */
  /* Generic class used as interface for the others */
  template<typename i_type, typename d_type,
	   template<typename> class ret_type, class daughter>
  class element_iterator;

  template<typename i_type, typename d_type,
	   template<typename> class ret_type, class daughter>
  class quadrature_point_iterator;

  template<typename T, class iterator_type,
	   template<typename> class ret_type>
  class GenericElementalField;

  template<typename T, class iterator_type,
	   template<typename> class ret_type>
  class GenericQuadraturePointsField;

  /* ------------------------------------------------------------------------ */
  /* Elemental Fields wrapper */
  class element_type_field_iterator;

  template<typename T, template<class> class ret_type>
  class elemental_field_iterator;

  class ElementTypeField;

  template<typename T, template<typename> class ret_type = types::Vector>
  class ElementalField;

  /* ------------------------------------------------------------------------ */
  /* Material Field wrapper */
  template<typename T, template<class> class ret_type>
  class internal_material_field_iterator;

  template<typename T, template<class> class ret_type = types::Vector>
  class InternalMaterialField;

  /* ------------------------------------------------------------------------ */
  /* Field homogenizing wrapper */
  template<typename T, class Container, template<class> class sub_type>
  class PaddingHomogenizingFunctor;

  template<typename T, class Container, template<class> class sub_type>
  class AvgHomogenizingFunctor;

  template<typename T, template< typename, template<class> class> class Container,
	   template<typename, class, template<class> class> class Funct = AvgHomogenizingFunctor,
	   template<typename> class ret_type = types::Vector>
  class HomogenizedField;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// internal iohelper::Dumper
  iohelper::Dumper * dumper;

  typedef std::map<std::string, Field *> Fields;
  /// list of registered fields to dump
  Fields fields;

  /// dump counter
  UInt count;

  /// filename
  std::string filename;
};

#include "dumper_iohelper_tmpl.hh"

__END_AKANTU__

#endif /* __AKANTU_DUMPER_IOHELPER_HH__ */
