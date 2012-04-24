/**
 * @file   weight_function.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Mar  8 16:17:00 2012
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

/* -------------------------------------------------------------------------- */
#include <vector>


#ifndef __AKANTU_WEIGHT_FUNCTION_HH__
#define __AKANTU_WEIGHT_FUNCTION_HH__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
class BaseWeightFunction {
public:
  BaseWeightFunction(const Material & material) : material(material) {}

  virtual ~BaseWeightFunction() {}

  virtual void init() {};

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
  bool setParam(__attribute__((unused)) const std::string & key,
		__attribute__((unused)) const std::string & value) {
    return false;
  }

  /* ------------------------------------------------------------------------ */
  virtual void printself(std::ostream & stream) const {
    stream << "BaseWeightFunction";
  }

protected:
  const Material & material;

  Real R;
  Real R2;
};

/* -------------------------------------------------------------------------- */
class DamagedWeightFunction : public BaseWeightFunction {
public:
  DamagedWeightFunction(const Material & material) : BaseWeightFunction(material) {}

  inline void selectType(__attribute__((unused)) ElementType type1,__attribute__((unused)) GhostType ghost_type1,
			 ElementType type2, GhostType ghost_type2) {
    selected_damage = &(dynamic_cast<const MaterialDamage &>(material).getDamage(type2, ghost_type2));
  }

  inline Real operator()(Real r, __attribute__((unused)) QuadraturePoint & q1, QuadraturePoint & q2) {
    UInt quad = q2.global_num;
    Real D = (*selected_damage)(quad);
    Real Radius = (1.-D)*(1.-D) * R2;
    if(Radius < Math::getTolerance()) {
      Radius = 0.01 * 0.01 * R2;
    }
    Real alpha = std::max(0., 1. - r*r / Radius);
    Real w = alpha * alpha;
    return w;
  }

  /* ------------------------------------------------------------------------ */
  bool setParam(__attribute__((unused)) const std::string & key,
		__attribute__((unused)) const std::string & value) {
    return false;
  }

  /* ------------------------------------------------------------------------ */
  virtual void printself(std::ostream & stream) const {
    stream << "DamageWeightFunction";
  }

private:
  const Vector<Real> * selected_damage;
};


/* -------------------------------------------------------------------------- */
/* Stress Based Weight                                                        */
/* -------------------------------------------------------------------------- */
class StressBasedWeightFunction : public BaseWeightFunction {
public:
  StressBasedWeightFunction(const Material & material);

  void init();

  void updatePrincipalStress(GhostType ghost_type);

  inline void updateQuadraturePointsCoordinates(ByElementTypeReal & quadrature_points_coordinates);

  inline void selectType(ElementType type1, GhostType ghost_type1,
			 ElementType type2, GhostType ghost_type2);

  inline Real operator()(Real r, QuadraturePoint & q1, QuadraturePoint & q2);

  /* ------------------------------------------------------------------------ */
  bool setParam(const std::string & key, const std::string & value) {
    std::stringstream sstr(value);
    if(key == "ft") { sstr >> ft; }
    else return false;
    return true;
  }

  /* ------------------------------------------------------------------------ */
  virtual void printself(std::ostream & stream) const {
    stream << "StressBasedWeightFunction [ft : " << ft << "]";
  }

private:
  Real ft;
  UInt spatial_dimension;

  ByElementTypeReal stress_diag;
  Vector<Real> * selected_stress_diag;
  ByElementTypeReal stress_base;
  Vector<Real> * selected_stress_base;


  ByElementTypeReal quadrature_points_coordinates;
  Vector<Real> * selected_position_1;
  Vector<Real> * selected_position_2;

  ByElementTypeReal characteristic_size;
  Vector<Real> * selected_characteristic_size;
};


inline std::ostream & operator <<(std::ostream & stream, const BaseWeightFunction & _this)
{
  _this.printself(stream);
  return stream;
}

#include "weight_function_tmpl.hh"

__END_AKANTU__

#endif /* __AKANTU_WEIGHT_FUNCTION_HH__ */
