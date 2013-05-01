/**
 * @file   weight_function.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Apr 13 16:49:32 2012
 *
 * @brief  Weight functions for non local materials
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
#include "aka_types.hh"
#include "solid_mechanics_model.hh"
#include <cmath>
#if defined(AKANTU_DEBUG_TOOLS)
#include "aka_debug_tools.hh"
#include <string>
#endif

/* -------------------------------------------------------------------------- */
#include <vector>


#ifndef __AKANTU_WEIGHT_FUNCTION_HH__
#define __AKANTU_WEIGHT_FUNCTION_HH__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/*  Normal weight function                                                    */
/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
class BaseWeightFunction {
public:
  BaseWeightFunction(Material & material) : material(material) {}

  virtual ~BaseWeightFunction() {}

  virtual void init() {};

  virtual void updateInternals(__attribute__((unused)) const ByElementTypeReal & quadrature_points_coordinates) {};

  /* ------------------------------------------------------------------------ */
  inline void setRadius(Real radius) { R = radius; R2 = R * R; }

  /* ------------------------------------------------------------------------ */
  inline void selectType(__attribute__((unused)) ElementType type1,
                         __attribute__((unused)) GhostType   ghost_type1,
                         __attribute__((unused)) ElementType type2,
                         __attribute__((unused)) GhostType   ghost_type2) {
  }

  /* ------------------------------------------------------------------------ */
  inline Real operator()(Real r,
                         __attribute__((unused)) QuadraturePoint & q1,
                         __attribute__((unused)) QuadraturePoint & q2) {
    Real w = 0;
    if(r <= R) {
      Real alpha = (1. - r*r / R2);
      w = alpha * alpha;
      // *weight = 1 - sqrt(r / radius);
    }
    return w;
  }

  /* ------------------------------------------------------------------------ */
  bool parseParam(__attribute__((unused)) const std::string & key,
                  __attribute__((unused)) const std::string & value) {
    return false;
  }

  /* ------------------------------------------------------------------------ */
  virtual void printself(std::ostream & stream) const {
    stream << "BaseWeightFunction";
  }

public:
  virtual UInt getNbDataForElements(__attribute__((unused)) const Array<Element> & elements,
                                    __attribute__((unused)) SynchronizationTag tag) const {
    return 0;
  }

  virtual inline void packElementData(__attribute__((unused)) CommunicationBuffer & buffer,
                                      __attribute__((unused)) const Array<Element> & elements,
                                      __attribute__((unused)) SynchronizationTag tag) const {}

  virtual inline void unpackElementData(__attribute__((unused)) CommunicationBuffer & buffer,
                                        __attribute__((unused)) const Array<Element> & elements,
                                        __attribute__((unused)) SynchronizationTag tag) {}
protected:
  Material & material;

  Real R;
  Real R2;
};

/* -------------------------------------------------------------------------- */
/*  Damage weight function                                                    */
/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
class DamagedWeightFunction : public BaseWeightFunction<spatial_dimension> {
public:
  DamagedWeightFunction(Material & material) : BaseWeightFunction<spatial_dimension>(material) {}

  inline void selectType(__attribute__((unused)) ElementType type1,
                         __attribute__((unused)) GhostType ghost_type1,
                         ElementType type2,
                         GhostType ghost_type2) {
    selected_damage = &this->material.getArray("damage", type2, ghost_type2);
  }

  inline Real operator()(Real r, __attribute__((unused)) QuadraturePoint & q1, QuadraturePoint & q2) {
    UInt quad = q2.global_num;
    Real D = (*selected_damage)(quad);
    Real Radius = (1.-D)*(1.-D) * this->R2;
    if(Radius < Math::getTolerance()) {
      Radius = 0.01 * 0.01 * this->R2;
    }
    Real alpha = std::max(0., 1. - r*r / Radius);
    Real w = alpha * alpha;
    return w;
  }

  /* ------------------------------------------------------------------------ */
  bool parseParam(__attribute__((unused)) const std::string & key,
                  __attribute__((unused)) const std::string & value) {
    return false;
  }

  /* ------------------------------------------------------------------------ */
  virtual void printself(std::ostream & stream) const {
    stream << "DamagedWeightFunction";
  }

private:
  const Array<Real> * selected_damage;
};

/* -------------------------------------------------------------------------- */
/*  Remove damaged weight function                                            */
/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
class RemoveDamagedWeightFunction : public BaseWeightFunction<spatial_dimension> {
public:
  RemoveDamagedWeightFunction(Material & material) : BaseWeightFunction<spatial_dimension>(material) {}

  inline void selectType(__attribute__((unused)) ElementType type1,
                         __attribute__((unused)) GhostType ghost_type1,
                         ElementType type2,
                         GhostType ghost_type2) {
    selected_damage = &this->material.getArray("damage", type2, ghost_type2);
  }

  inline Real operator()(Real r, __attribute__((unused)) QuadraturePoint & q1, QuadraturePoint & q2) {
    UInt quad = q2.global_num;

    if(q1 == q2) return 1.;

    Real D = (*selected_damage)(quad);
    Real w = 0.;
    if(D < damage_limit) {
      Real alpha = std::max(0., 1. - r*r / this->R2);
      w = alpha * alpha;
    }
    return w;
  }

  /* ------------------------------------------------------------------------ */
  bool parseParam(const std::string & key,
                  const std::string & value) {
    std::stringstream sstr(value);
    if(key == "damage_limit") { sstr >> damage_limit; }
    else return BaseWeightFunction<spatial_dimension>::parseParam(key, value);
    return true;
  }

  /* ------------------------------------------------------------------------ */
  virtual void printself(std::ostream & stream) const {
    stream << "RemoveDamagedWeightFunction [damage_limit: " << damage_limit << "]";
  }

  virtual UInt getNbDataForElements(const Array<Element> & elements,
                                    SynchronizationTag tag) const {
    if(tag == _gst_mnl_weight)
      return this->material.getModel().getNbQuadraturePoints(elements) * sizeof(Real);

    return 0;
  }

  virtual inline void packElementData(CommunicationBuffer & buffer,
                                      const Array<Element> & elements,
                                      SynchronizationTag tag) const {
    if(tag == _gst_mnl_weight) {
      MaterialDamage<spatial_dimension> & mat = dynamic_cast<MaterialDamage<spatial_dimension> &>(this->material);
      this->material.packElementDataHelper(mat.getDamage(),
                                           buffer,
                                           elements);
#if defined(AKANTU_DEBUG_TOOLS)
#if defined(AKANTU_CORE_CXX11)
      debug::element_manager.print(debug::_dm_material,
                                   [&elements, &mat](const Element & el)->std::string {
                                     std::stringstream out;
                                     UInt pos = elements.find(el);
                                     if(pos != UInt(-1)) {
                                       Real d = mat.getArray("damage", el.type, el.ghost_type)(el.element);
                                       if(d > 0.3)
                                         out << " damage sent: " << d;
                                     }
                                     return out.str();
                                   });
#else
      debug::element_manager.printData(debug::_dm_material, "RemoveDamagedWeightFunction: packElementData",
                                       mat.getDamage(), mat.getElementFilter());
#endif
#endif
    }
  }

  virtual inline void unpackElementData(CommunicationBuffer & buffer,
                                        const Array<Element> & elements,
                                        SynchronizationTag tag) {
    if(tag == _gst_mnl_weight) {
      MaterialDamage<spatial_dimension> & mat = dynamic_cast<MaterialDamage<spatial_dimension> &>(this->material);
      this->material.unpackElementDataHelper(mat.getDamage(),
                                             buffer,
                                             elements);
#if defined(AKANTU_DEBUG_TOOLS)
#if defined(AKANTU_CORE_CXX11)
      debug::element_manager.print(debug::_dm_material,
                                   [&elements, &mat](const Element & el)->std::string {
                                     std::stringstream out;
                                     UInt pos = elements.find(el);
                                     if(pos != UInt(-1)) {
                                       Real d = mat.getArray("damage", el.type, el.ghost_type)(el.element);
                                       if(d > 0.3)
                                         out << " damage recv: " << d;
                                     }
                                     return out.str();
                                   });
#else
      debug::element_manager.printData(debug::_dm_material, "RemoveDamagedWeightFunction: unpackElementData",
                                       mat.getDamage(), mat.getElementFilter());
#endif
#endif

    }
  }


private:
  /// limit at which a point is considered as complitely broken
  Real damage_limit;

  /// internal pointer to the current damage vector
  const Array<Real> * selected_damage;
};
/* -------------------------------------------------------------------------- */
/* Remove damaged with damage rate weight function                                             */
/* -------------------------------------------------------------------------- */

template<UInt spatial_dimension>
class RemoveDamagedWithDamageRateWeightFunction : public BaseWeightFunction<spatial_dimension> {
public:
  RemoveDamagedWithDamageRateWeightFunction(Material & material) : BaseWeightFunction<spatial_dimension>(material) {}

  inline void selectType(__attribute__((unused)) ElementType type1,
                         __attribute__((unused)) GhostType ghost_type1,
                         ElementType type2,
                         GhostType ghost_type2) {
    selected_damage_with_damage_rate = &(this->material.getArray("damage",type2, ghost_type2));
    selected_damage_rate_with_damage_rate = &(this->material.getArray("damage-rate",type2, ghost_type2));
  }

  inline Real operator()(Real r, __attribute__((unused)) QuadraturePoint & q1, QuadraturePoint & q2) {
    UInt quad = q2.global_num;

    if(q1.global_num == quad) return 1.;

    Real D = (*selected_damage_with_damage_rate)(quad);
    Real w = 0.;
    Real alphaexp = 1.;
    Real betaexp = 2.;
    if(D < damage_limit_with_damage_rate) {
      Real alpha = std::max(0., 1. - pow((r*r / this->R2),alphaexp));
      w = pow(alpha, betaexp);
    }

    return w;
  }
  /* ------------------------------------------------------------------------ */
  bool setParam(const std::string & key,
                const std::string & value) {
    std::stringstream sstr(value);
    if(key == "damage_limit") { sstr >> damage_limit_with_damage_rate; }
    else return BaseWeightFunction<spatial_dimension>::setParam(key, value);
    return true;
  }

  /* ------------------------------------------------------------------------ */
  virtual void printself(std::ostream & stream) const {
    stream << "RemoveDamagedWithDamageRateWeightFunction [damage_limit: " << damage_limit_with_damage_rate << "]";
  }

private:
  /// limit at which a point is considered as complitely broken
  Real damage_limit_with_damage_rate;

  /// internal pointer to the current damage vector
  const Array<Real> * selected_damage_with_damage_rate;

  /// internal pointer to the current damage rate vector
  const Array<Real> * selected_damage_rate_with_damage_rate;
};

/* -------------------------------------------------------------------------- */
/* Stress Based Weight                                                        */
/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
class StressBasedWeightFunction : public BaseWeightFunction<spatial_dimension> {
public:
  StressBasedWeightFunction(Material & material);

  void init();

  virtual void updateInternals(__attribute__((unused)) const ByElementTypeReal & quadrature_points_coordinates) {
    updatePrincipalStress(_not_ghost);
    updatePrincipalStress(_ghost);
  };

  void updatePrincipalStress(GhostType ghost_type);

  inline void updateQuadraturePointsCoordinates(ByElementTypeReal & quadrature_points_coordinates);

  inline void selectType(ElementType type1, GhostType ghost_type1,
                         ElementType type2, GhostType ghost_type2);

  inline Real operator()(Real r, QuadraturePoint & q1, QuadraturePoint & q2);

  inline Real computeRhoSquare(Real r,
                               Vector<Real> & eigs,
                               Matrix<Real> & eigenvects,
                               Vector<Real> & x_s);

  /* ------------------------------------------------------------------------ */
  bool parseParam(const std::string & key, const std::string & value) {
    std::stringstream sstr(value);
    if(key == "ft") { sstr >> ft; }
    else return BaseWeightFunction<spatial_dimension>::parseParam(key, value);
    return true;
  }

  /* ------------------------------------------------------------------------ */
  virtual void printself(std::ostream & stream) const {
    stream << "StressBasedWeightFunction [ft: " << ft << "]";
  }

private:
  Real ft;

  ByElementTypeReal stress_diag;
  Array<Real> * selected_stress_diag;
  ByElementTypeReal stress_base;
  Array<Real> * selected_stress_base;


  ByElementTypeReal quadrature_points_coordinates;
  Array<Real> * selected_position_1;
  Array<Real> * selected_position_2;

  ByElementTypeReal characteristic_size;
  Array<Real> * selected_characteristic_size;
};

template<UInt spatial_dimension>
inline std::ostream & operator <<(std::ostream & stream,
                                  const BaseWeightFunction<spatial_dimension> & _this)
{
  _this.printself(stream);
  return stream;
}


template<UInt spatial_dimension>
inline std::ostream & operator <<(std::ostream & stream,
                                  const BaseWeightFunction<spatial_dimension> * _this)
{
  _this->printself(stream);
  return stream;
}


template<UInt d>
inline std::ostream & operator >>(std::ostream & stream,
                                  __attribute__((unused)) const BaseWeightFunction<d> * _this)
{  AKANTU_DEBUG_TO_IMPLEMENT(); return stream; }
template<UInt d>
inline std::ostream & operator >>(std::ostream & stream,
                                  __attribute__((unused)) const RemoveDamagedWeightFunction<d> * _this)
{  AKANTU_DEBUG_TO_IMPLEMENT(); return stream; }
template<UInt d>
inline std::ostream & operator >>(std::ostream & stream,
                                  __attribute__((unused)) const DamagedWeightFunction<d> * _this)
{  AKANTU_DEBUG_TO_IMPLEMENT(); return stream; }
template<UInt d>
inline std::ostream & operator >>(std::ostream & stream,
                                  __attribute__((unused)) const StressBasedWeightFunction<d> * _this)
{  AKANTU_DEBUG_TO_IMPLEMENT(); return stream; }


#include "weight_function_tmpl.hh"

__END_AKANTU__

#endif /* __AKANTU_WEIGHT_FUNCTION_HH__ */
