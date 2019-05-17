/**
 * @file   resolution.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon Jan 7 2019
 * @date last modification: Mon Jan 7 2019
 *
 * @brief  Mother class for all contact resolutions
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_factory.hh"
#include "aka_memory.hh"
#include "parsable.hh"
#include "parser.hh"
#include "fe_engine.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_RESOLUTION_HH__
#define __AKANTU_RESOLUTION_HH__

/* -------------------------------------------------------------------------- */
namespace akantu {
  class Model;
  class ContactMechanicsModel;
} // namespace akantu


namespace akantu {

/**
 * Interface of all contact resolutions
 * Prerequisites for a new resolution
 * - inherit from this class
 * - implement the following methods:
 * \code
 *
 *  virtual void computeNormalForce();
 *  virtual void computeFrictionForce();
 *
 * \endcode
 *
 */
class Resolution : public Memory,
		   public Parsable {

  /* ------------------------------------------------------------------------ */
  /* Constructor/Destructor                                                   */
  /* ------------------------------------------------------------------------ */
public:
  /// instantiate contact resolution with defaults
  Resolution(ContactMechanicsModel & model, const ID & id = "");

  /// Destructor
  ~Resolution() override;

protected:
  void initialize();

  /// computes coordinates of a given element
  void computeCoordinates(const Element & , Matrix<Real> &);

  /// computes tangents
  void computeTangents(Matrix<Real> & ,
		       Matrix<Real> & ,
		       Matrix<Real> & );

  /// computes metric tensor (@f$m_{\alpha\beta}@f$) where @f$\alpha,
  /// \beta@f$ are surface directions
  void computeMetricTensor(Matrix<Real> &,
			   Matrix<Real> &);

  /// computes N array
  void computeN(Vector<Real> &,
		Vector<Real> &,
		Vector<Real> &);

  /// computes (@f$N_{\alpha}@f$) where \alpha is surface dimension
  void computeNalpha(Array<Real>  &,
		     Matrix<Real> &,
		     Vector<Real> &);

  /// computes (@f$T_{\alpha}@f$) where @f$\alpha@f$ is surface dimension
  void computeTalpha(Array<Real> & ,
		     Vector<Real> & ,
		     Matrix<Real> & );

  /// computes (@f$D_{\alpha}@f$) where @f$\alpha@f$ is surface dimension
  void computeDalpha(Array<Real> & ,
		     Array<Real> & ,
		     Array<Real> & ,
		     Matrix<Real> &,
		     Real & ); 

  /* ------------------------------------------------------------------------ */
  /* Functions that resolutions can/should reimplement                        */
  /* ------------------------------------------------------------------------ */
protected:
  /// computes the normal force
  virtual void computeNormalForce(__attribute__((unused)) Vector<Real> & /* force */,
				  __attribute__((unused)) Vector<Real>  & /* n */,
				  __attribute__((unused)) Real & /* gap */) {
    AKANTU_TO_IMPLEMENT();
  }

  /// computes the friction force
  virtual void computeFrictionForce(__attribute__((unused)) Vector<Real> & /* force */,
				    __attribute__((unused)) Array<Real>  & /* d_alpha  */,
				    __attribute__((unused)) Matrix<Real> &,
				    __attribute__((unused)) Vector<Real> &,
				    __attribute__((unused)) Real & /* gap */) {
    AKANTU_TO_IMPLEMENT();
  }

  /// compute the tangent moduli
  virtual void computeTangentModuli(__attribute__((unused)) Matrix<Real> &,
				    __attribute__((unused)) Vector<Real> & /* n */,
				    __attribute__((unused)) Array<Real> & /* n_alpha */,
				    __attribute__((unused)) Array<Real> & /* d_alpha */,
				    __attribute__((unused)) Matrix<Real> & /* A */,
				    __attribute__((unused)) Real & /* gap */) {
    AKANTU_TO_IMPLEMENT();
  }

  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */ 
public:
  /// assemble the residual for this resolution
  void assembleInternalForces(GhostType ghost_type);

  /// assemble the stiffness matrix for this resolution
  void assembleStiffnessMatrix(GhostType ghost_type);

public:
  /// function to print the contain of the class
  void printself(std::ostream & stream, int indent = 0) const override;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:  
  /// Link to the fem object in the model
  FEEngine & fem;
  
  /// resolution name
  std::string name;
  
  /// model to which the resolution belong
  ContactMechanicsModel & model;

  /// friciton coefficient : mu
  Real mu;
  
  /// spatial dimension
  UInt spatial_dimension;
				   
};

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const Resolution & _this) {
  _this.printself(stream);
  return stream;
}

  
} // namespace akantu

/* -------------------------------------------------------------------------- */
namespace akantu {
using ResolutionFactory =
    Factory<Resolution, ID, UInt, const ID &, ContactMechanicsModel &, const ID &>;

/// macaulay bracket to convert  positive gap to zero  
template <typename T>
T macaulay(T var) {return var < 0 ? 0 : var; }

template <typename T>
T heaviside(T var) {return var < 0 ? 0 : 1;  }
} // namespace akantu

#define INSTANTIATE_RESOLUTION_ONLY(res_name)                                    \
  class res_name                                                  
 
#define RESOLUTION_DEFAULT_PER_DIM_ALLOCATOR(id, res_name)                       \
  [](UInt dim, const ID &, ContactMechanicsModel & model,                        \
     const ID & id) -> std::unique_ptr<Resolution> {                             \
    switch (dim) {							\
    case 1:								\
      return std::make_unique<res_name>(model, id);			\
    case 2:								\
      return std::make_unique<res_name>(model, id);			\
    case 3:								\
      return std::make_unique<res_name>(model, id);			\
    default:								\
      AKANTU_EXCEPTION("The dimension "					\
                       << dim << "is not a valid dimension for the contact resolution "	\
                       << #id);						\
    }									\
  }


#define INSTANTIATE_RESOLUTION(id, res_name)                                     \
  INSTANTIATE_RESOLUTION_ONLY(res_name);                                         \
  static bool resolution_is_alocated_##id[[gnu::unused]] =                       \
      ResolutionFactory::getInstance().registerAllocator(                        \
          #id, RESOLUTION_DEFAULT_PER_DIM_ALLOCATOR(id, res_name))

#endif /* __AKANTU_RESOLUTION_HH__  */


