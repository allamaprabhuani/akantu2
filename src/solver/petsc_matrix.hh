/**
 * @file   petsc_matrix.hh
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Mon Jul 21 14:49:49 2014
 *
 * @brief  Interface for PETSc matrices
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

#ifndef __AKANTU_PETSC_MATRIX_HH__
#define __AKANTU_PETSC_MATRIX_HH__

/* -------------------------------------------------------------------------- */
#include "sparse_matrix.hh"

struct Mat;
struct AO;

/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__



class PetscMatrix : SparseMatrix {
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */  

  class PetscMatrixValue {

  public:

    PetscMatrixValue(){}

    inline PetscMatrixValue operator = (Real v){
      MatSetValue(*mat,i,j,v,INSERT_VALUES);
      return *this;
    };


    inline PetscMatrixValue operator += (Real v){
      //     std::cout << "set " << std::abs(Real(i)-Real(j)) << std::endl;
      MatSetValue(*mat,i,j,v,ADD_VALUES);
      return *this;
    };

    inline PetscMatrixValue operator -= (Real v){
      MatSetValue(*mat,i,j,-1.*v,ADD_VALUES);
      return *this;
    };

    friend PetscMatrix;

  private:

    UInt i,j;
    Mat * mat;
  };

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  PetscMatrix(UInt size,
	       const SparseMatrixType & sparse_matrix_type,
	       UInt nb_degree_of_freedom,
	       const ID & id = "petsc_matrix",
	       const MemoryID & memory_id = 0);

  PetscMatrix(const PetscMatrix & matrix,
	       const ID & id = "petsc_matrix",
	       const MemoryID & memory_id = 0);

  virtual ~PetscMatrix();

  typedef std::pair<UInt, UInt> KeyCOO;
  typedef unordered_map<KeyCOO, UInt>::type coordinate_list_map;


  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
  

private:

  Mat * mat;
  AO  * ao;


};

__END_AKANTU__

#endif /* __AKANTU_PETSCI_MATRIX_HH__ */
