/**
 * Copyright (©) 2013-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include <type_traits>

/* -------------------------------------------------------------------------- */
#ifndef AKANTU_DUMPER_IOHELPER_TMPL_VARIABLE_HH_
#define AKANTU_DUMPER_IOHELPER_TMPL_VARIABLE_HH_
/* -------------------------------------------------------------------------- */

namespace akantu {
namespace dumpers {
  /* --------------------------------------------------------------------------
   */

  /// Variable interface
  class VariableBase {
  public:
    VariableBase() = default;
    virtual ~VariableBase() = default;
    virtual void registerToDumper(const std::string & id,
                                  iohelper::Dumper & dumper) = 0;
  };

  /* --------------------------------------------------------------------------
   */

  template <typename T, bool is_scal = std::is_arithmetic<T>::value>
  class Variable : public VariableBase {
  public:
    Variable(const T & t) : vari(t) {}

    void registerToDumper(const std::string & id,
                          iohelper::Dumper & dumper) override {
      dumper.addVariable(id, *this);
    }

    const T & operator[](Idx i) const { return vari[i]; }

    Int getDim() { return vari.size(); }
    iohelper::DataType getDataType() { return iohelper::getDataType<T>(); }

  protected:
    const T & vari;
  };

  /* --------------------------------------------------------------------------
   */
  template <typename T> class Variable<Vector<T>, false> : public VariableBase {
  public:
    Variable(const Vector<T> & t) : vari(t) {}

    void registerToDumper(const std::string & id,
                          iohelper::Dumper & dumper) override {
      dumper.addVariable(id, *this);
    }

    const T & operator[](Idx i) const { return vari[i]; }

    Int getDim() { return vari.size(); }
    iohelper::DataType getDataType() { return iohelper::getDataType<T>(); }

  protected:
    const Vector<T> & vari;
  };

  /* --------------------------------------------------------------------------
   */

  template <typename T> class Variable<T, true> : public VariableBase {
  public:
    Variable(const T & t) : vari(t) {}

    void registerToDumper(const std::string & id,
                          iohelper::Dumper & dumper) override {
      dumper.addVariable(id, *this);
    }

    const T & operator[](Idx /*i*/) const { return vari; }

    Int getDim() { return 1; }
    iohelper::DataType getDataType() { return iohelper::getDataType<T>(); }

  protected:
    const T & vari;
  };

} // namespace dumpers
} // namespace akantu

#endif /* AKANTU_DUMPER_IOHELPER_TMPL_VARIABLE_HH_ */
