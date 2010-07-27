/**
 * @file   material_elastic.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 11:53:52 2010
 *
 * @brief  Specialization of the material class for the elastic material
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */

template<> class Material<_elastic> : public MaterialBase {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Material(SolidMechanicsModel & model, const MaterialID & id = "")  :
    MaterialBase(model, id) {
    AKANTU_DEBUG_IN();

    rho = 1;
    young_modulus = 1;
    nu = 1;

    AKANTU_DEBUG_OUT();
  }

  virtual ~Material() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {
    std::string space;
    for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

    stream << space << "Material<_elastic> [" << std::endl;
    stream << space << " + id              : " << id << std::endl;
    stream << space << " + density         : " << rho << std::endl;
    stream << space << " + Young modulus   : " << young_modulus << std::endl;
    stream << space << " + Poisson's ratio : " << nu << std::endl;
    stream << space << "]" << std::endl;
  }


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// the young modulus
  UInt young_modulus;

  /// Poisson coefficient
  UInt nu;
};
