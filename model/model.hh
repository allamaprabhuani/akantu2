/**
 * @file   model.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 22 11:02:42 2010
 *
 * @brief  Interface of a model
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MODEL_HH__
#define __AKANTU_MODEL_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_memory.hh"
#include "mesh.hh"
#include "fem.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class Model : public Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Model(UInt spatial_dimension,
	const ModelID & id = "model",
	const MemoryID & memory_id = 0) :
    Memory(memory_id), id(id), spatial_dimension(spatial_dimension) {
    AKANTU_DEBUG_IN();
    std::stringstream sstr; sstr << id << ":fem";
    this->fem = new FEM(spatial_dimension, sstr.str(), memory_id);
    AKANTU_DEBUG_OUT();
  };

  Model(Mesh & mesh,
	UInt spatial_dimension = 0,
	const ModelID & id = "model",
	const MemoryID & memory_id = 0) :
    Memory(memory_id), id(id) {
    AKANTU_DEBUG_IN();
    this->spatial_dimension = (spatial_dimension == 0) ? mesh.getSpatialDimension() : spatial_dimension;
    std::stringstream sstr; sstr << id << ":fem";
    this->fem = new FEM(mesh, spatial_dimension, sstr.str(), memory_id);
    AKANTU_DEBUG_OUT();
  };

  virtual ~Model() {
    AKANTU_DEBUG_IN();
    delete fem;
    AKANTU_DEBUG_OUT();
  };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  virtual void initModel() = 0;

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const = 0;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                 */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);

  AKANTU_GET_MACRO(FEM, *fem, const FEM &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// id
  ModelID id;

  /// the spatial dimension
  UInt spatial_dimension;

  /// the main fem object present in all  models
  FEM * fem;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

///#include "model_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const Model & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__

#endif /* __AKANTU_MODEL_HH__ */
