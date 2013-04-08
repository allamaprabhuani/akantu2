/**
 * @file   material_finite_deformation.cc
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 * @date   Wed Feb 22 16:31:20 2012
 *
 * @brief  Specialization of the material class for finite deformation
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
#include "material_finite_deformation.hh"
#include "solid_mechanics_model.hh"
#include "sparse_matrix.hh"
#include "dof_synchronizer.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MaterialFiniteDeformation::MaterialFiniteDeformation(SolidMechanicsModel & model, const ID & id) :
Material(model, id),
delta_stress("delta_stress", id),
delta_strain("delta_strain", id),
stress_at_t("stress_at_t", id) {
    AKANTU_DEBUG_IN();
    
    finite_deformation=true;

    /// allocate strain stress for local elements
    initInternalArray(delta_strain, spatial_dimension * spatial_dimension);
    initInternalArray(delta_stress, spatial_dimension * spatial_dimension);
    initInternalArray(stress_at_t, spatial_dimension * spatial_dimension);

    AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
MaterialFiniteDeformation::~MaterialFiniteDeformation() {
    AKANTU_DEBUG_IN();

    AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialFiniteDeformation::initMaterial() {
    AKANTU_DEBUG_IN();

    Material::initMaterial();
    resizeInternalArray(delta_stress);
    resizeInternalArray(delta_strain);
    resizeInternalArray(stress_at_t);

    AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

/**
 * Compute  the  residual  by  assembling  @f$\int_{e}  t_e  N_e dS @f$
 *
 * @param[in] displacements nodes displacements
 * @param[in] ghost_type compute the residual for _ghost or _not_ghost element
 */

void MaterialFiniteDeformation::UpdateStressesAtT(GhostType ghost_type) {
    AKANTU_DEBUG_IN();

    //resizeInternalArray(stress_at_t);

    UInt spatial_dimension = model->getSpatialDimension();

    Mesh::type_iterator it = model->getFEM().getMesh().firstType(spatial_dimension, ghost_type);
    Mesh::type_iterator last_type = model->getFEM().getMesh().lastType(spatial_dimension, ghost_type);

    for (; it != last_type; ++it) {
        UInt dim= spatial_dimension;
        
        Array<Real>::iterator< Matrix<Real> > stress_it = stress(*it, ghost_type).begin(dim, dim);
        Array<Real>::iterator< Matrix<Real> > stress_end = stress(*it, ghost_type).end(dim, dim);
        Array<Real>::iterator< Matrix<Real> > stress_at_t_it = stress_at_t(*it, ghost_type).begin(dim, dim);
        Array<Real>::iterator< Matrix<Real> > strain_it = strain(*it, ghost_type).begin(dim, dim);

        
        Matrix<Real> F(spatial_dimension, spatial_dimension);
        Matrix<Real> FtS(spatial_dimension, spatial_dimension);


        for (; stress_it != stress_end ; ++stress_it, ++stress_at_t_it, ++strain_it){
            Matrix<Real> & grad_u = *strain_it;
            Matrix<Real> & s_t = *stress_it;
            Matrix<Real> & cauchy_t = *stress_at_t_it;
            
            Real J;
            switch (spatial_dimension){
                case 3:
                    gradUToF<3 > (grad_u, F);
                    J=Math::det3(F.storage());
                    break;
                case 2:
                    gradUToF<2 > (grad_u, F);
                    J=Math::det2(F.storage());
                    break;
            }
            
            
            FtS.mul<true, false>(F,s_t);
            cauchy_t.mul<false, false>(FtS,F);
            cauchy_t *= 1.0/J;
            for(UInt i=0;i<spatial_dimension;i++)
                for(UInt j=0;j<spatial_dimension;j++)
                    cauchy_t(i,j)=s_t(i,j);
            
            
        }
            

    }

    AKANTU_DEBUG_OUT();
}

/*
 * Incremental computation of the stress and the strain using Updated Lagrangian Formulation
 */
void MaterialFiniteDeformation::computeAllStresses(GhostType ghost_type) {
    AKANTU_DEBUG_IN();

    resizeInternalArray(delta_stress);
    resizeInternalArray(delta_strain);

    UInt spatial_dimension = model->getSpatialDimension();

    Mesh::type_iterator it = model->getFEM().getMesh().firstType(spatial_dimension, ghost_type);
    Mesh::type_iterator last_type = model->getFEM().getMesh().lastType(spatial_dimension, ghost_type);

    for (; it != last_type; ++it) {
        Array<UInt> & elem_filter = element_filter(*it, ghost_type);
        Array<Real> & delta_strain_vect = delta_strain(*it, ghost_type);
        Array<Real> & strain_vect = strain(*it, ghost_type);

        /// compute @f$\nabla u@f$
        model->getFEM().gradientOnQuadraturePoints(model->getDisplacement(), strain_vect,
                spatial_dimension,
                *it, ghost_type, &elem_filter);


        /// compute @f$\nabla u@f$
        model->getFEM().gradientOnQuadraturePoints(model->getIncrement(), delta_strain_vect,
                spatial_dimension,
                *it, ghost_type, &elem_filter);

        /// compute @f$\mathbf{\sigma}_q@f$ from @f$\nabla u@f$
        computeStress(*it, ghost_type);
    }

    AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialFiniteDeformation::assembleResidual(GhostType ghost_type) {
    AKANTU_DEBUG_IN();

    UInt spatial_dimension = model->getSpatialDimension();

    Array<Real> & residual = const_cast<Array<Real> &> (model->getResidual());

    Mesh & mesh = model->getFEM().getMesh();
    Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type);
    Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type);
    for (; it != last_type; ++it) {
        const Array<Real> & shapes_derivatives = model->getFEM().getShapesDerivatives(*it, ghost_type);

        Array<UInt> & elem_filter = element_filter(*it, ghost_type);
        Array<Real> & strain_vect = strain(*it, ghost_type);


        model->getFEM().gradientOnQuadraturePoints(model->getDisplacement(), strain_vect,
                spatial_dimension, *it, ghost_type, &elem_filter);


        UInt nb_element = elem_filter.getSize();
        UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);
        UInt nb_quadrature_points = model->getFEM().getNbQuadraturePoints(*it, ghost_type);


        Array<Real> * shapes_derivatives_filtered = new Array<Real > (nb_element * nb_quadrature_points,
                spatial_dimension * nb_nodes_per_element,
                "shapes derivatives filtered");


        Array<Real>::const_iterator< Matrix<Real> > shapes_derivatives_it = shapes_derivatives.begin(spatial_dimension,
                nb_nodes_per_element);

        Array<Real>::iterator< Matrix<Real> > shapes_derivatives_filtered_it = shapes_derivatives_filtered->begin(spatial_dimension,
                nb_nodes_per_element);
        UInt * elem_filter_val = elem_filter.storage();
        for (UInt e = 0; e < nb_element; ++e, ++elem_filter_val)
            for (UInt q = 0; q < nb_quadrature_points; ++q, ++shapes_derivatives_filtered_it)
                *shapes_derivatives_filtered_it = shapes_derivatives_it[*elem_filter_val * nb_quadrature_points + q];


        //Set stress vectors

        UInt stress_size = getCauchyStressArraySize(spatial_dimension);

        Array<Real> * stress_voigt_vector =
                new Array<Real > (nb_element*nb_quadrature_points, stress_size,
                "tangent_stiffness_matrix");

        stress_voigt_vector->clear();

        Array<Real>::iterator< Matrix<Real> > S_it =
                this->stress(*it, ghost_type).begin(spatial_dimension, spatial_dimension);

        Array<Real>::iterator< Matrix<Real> > S_end =
                this->stress(*it, ghost_type).end(spatial_dimension, spatial_dimension);

        Array<Real>::iterator< Matrix<Real> > stress_vect_it = stress_voigt_vector->begin(stress_size, 1);


        for (; S_it != S_end; ++S_it, ++stress_vect_it) {
            Matrix<Real> & S = *S_it;
            Matrix<Real> & S_vect = *stress_vect_it;
            switch (spatial_dimension) {
                case 1:
                    SetCauchyStressArray < 1 > (S, S_vect);
                    break;
                case 2:
                    SetCauchyStressArray < 2 > (S, S_vect);
                    break;
                case 3:
                    SetCauchyStressArray < 3 > (S, S_vect);
                    break;
            }
        }

        //Set matrices B and BNL*


        /// compute @f$\mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
        UInt bt_s_size = spatial_dimension * nb_nodes_per_element;

        Array<Real> * bt_s = new Array<Real > (nb_element * nb_quadrature_points,
                bt_s_size,
                "B^t*S");

        Array<Real> * B = new Array<Real > (1, stress_size * bt_s_size, "B");
        Array<Real> * B2 = new Array<Real > (1, stress_size * bt_s_size, "B2");

        Array<Real>::iterator< Matrix<Real> > B_it = B->begin(stress_size, bt_s_size);
        Array<Real>::iterator< Matrix<Real> > B2_it = B2->begin(stress_size, bt_s_size);

        stress_vect_it = stress_voigt_vector->begin(stress_size, 1);
        Array<Real>::iterator< Matrix<Real> > stress_vect_end = stress_voigt_vector->end(stress_size, 1);
        Array<Real>::iterator< Matrix<Real> > grad_u_it = this->strain(*it, ghost_type).begin(spatial_dimension, spatial_dimension);

        shapes_derivatives_filtered_it = shapes_derivatives_filtered->begin(spatial_dimension, nb_nodes_per_element);


        Array<Real>::iterator< Matrix<Real> > bt_s_it = bt_s->begin(bt_s_size, 1);

        Matrix<Real> & B_tensor = *B_it;
        Matrix<Real> & B2_tensor = *B2_it;

        for (; stress_vect_it != stress_vect_end; ++stress_vect_it, ++shapes_derivatives_filtered_it, ++bt_s_it) {
            Matrix<Real> & grad_u = *grad_u_it;
            Matrix<Real> & r_it = *bt_s_it;
            Matrix<Real> & S_it = *stress_vect_it;

            switch (spatial_dimension) {
                case 1:
                {
                    transferBMatrixToBL1 < 1 > (*shapes_derivatives_filtered_it, B_tensor, nb_nodes_per_element);
                    transferBMatrixToBL2 < 1 > (*shapes_derivatives_filtered_it, grad_u, B2_tensor, nb_nodes_per_element);
                    break;
                }
                case 2:
                {
                    transferBMatrixToBL1 < 2 > (*shapes_derivatives_filtered_it, B_tensor, nb_nodes_per_element);
                    transferBMatrixToBL2 < 2 > (*shapes_derivatives_filtered_it, grad_u, B2_tensor, nb_nodes_per_element);
                    break;
                }
                case 3:
                {
                    transferBMatrixToBL1 < 3 > (*shapes_derivatives_filtered_it, B_tensor, nb_nodes_per_element);
                    transferBMatrixToBL2 < 3 > (*shapes_derivatives_filtered_it, grad_u, B2_tensor, nb_nodes_per_element);
                    break;
                }
            }

            B_tensor += B2_tensor;

            r_it.mul < true, false > (B_tensor, S_it);
        }

        delete B;
        delete B2;

        /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
        Array<Real> * r_e = new Array<Real > (nb_element,
                bt_s_size, "r_e");

        model->getFEM().integrate(*bt_s, *r_e,
                bt_s_size,
                *it, ghost_type,
                &elem_filter);

        delete bt_s;

        model->getFEM().assembleArray(*r_e, residual,
                model->getDOFSynchronizer().getLocalDOFEquationNumbers(),
                residual.getNbComponent(),
                *it, ghost_type, &elem_filter, -1);

        delete r_e;
        
    }
    AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialFiniteDeformation::assembleStiffnessMatrix(GhostType ghost_type) {

    AKANTU_DEBUG_IN();

    //Material::assembleStiffnessMatrix(ghost_type);
    UInt spatial_dimension = model->getSpatialDimension();

    Mesh & mesh = model->getFEM().getMesh();
    Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type);
    Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type);
    for (; it != last_type; ++it) {
        switch (spatial_dimension) {
            case 1:
            {
                assembleStiffnessMatrixNL < 1 > (*it, ghost_type);
                assembleStiffnessMatrixL2 < 1 > (*it, ghost_type);
                break;
            }
            case 2:
            {
                assembleStiffnessMatrixNL < 2 > (*it, ghost_type);
                assembleStiffnessMatrixL2 < 2 > (*it, ghost_type);
                break;
            }
            case 3:
            {
                assembleStiffnessMatrixNL < 3 > (*it, ghost_type);
                assembleStiffnessMatrixL2 < 3 > (*it, ghost_type);
                break;
            }
        }
    }
    
    AKANTU_DEBUG_OUT();
}

template<UInt dim>
void MaterialFiniteDeformation::assembleStiffnessMatrixNL(const ElementType & type,
        GhostType ghost_type) {
    AKANTU_DEBUG_IN();

    SparseMatrix & K = const_cast<SparseMatrix &> (model->getStiffnessMatrix());

    const Array<Real> & shapes_derivatives = model->getFEM().getShapesDerivatives(type, ghost_type);

    Array<UInt> & elem_filter = element_filter(type, ghost_type);
    //Array<Real> & strain_vect = delta_strain(type, ghost_type);

    UInt nb_element = elem_filter.getSize();
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    UInt nb_quadrature_points = model->getFEM().getNbQuadraturePoints(type, ghost_type);

    //strain_vect.resize(nb_quadrature_points * nb_element);

    // model->getFEM().gradientOnQuadraturePoints(model->getIncrement(), strain_vect,
    //        dim, type, ghost_type, &elem_filter);

    UInt cauchy_matrix_size = getCauchyStressMatrixSize(dim);

    Array<Real> * cauchy_stress_matrix =
            new Array<Real > (nb_element*nb_quadrature_points, cauchy_matrix_size * cauchy_matrix_size,
            "cauchy_stress_matrix");

    cauchy_stress_matrix->clear();

    SetCauchyStressMatrix(type, *cauchy_stress_matrix, ghost_type); //TODO: replace by setcauchystressmatrix()


    Array<Real> * shapes_derivatives_filtered = new Array<Real > (nb_element * nb_quadrature_points,
            dim * nb_nodes_per_element,
            "shapes derivatives filtered");


    Array<Real>::const_iterator< Matrix<Real> > shapes_derivatives_it = shapes_derivatives.begin(spatial_dimension,
            nb_nodes_per_element);

    Array<Real>::iterator< Matrix<Real> > shapes_derivatives_filtered_it = shapes_derivatives_filtered->begin(spatial_dimension,
            nb_nodes_per_element);
    UInt * elem_filter_val = elem_filter.storage();
    for (UInt e = 0; e < nb_element; ++e, ++elem_filter_val)
        for (UInt q = 0; q < nb_quadrature_points; ++q, ++shapes_derivatives_filtered_it)
            *shapes_derivatives_filtered_it = shapes_derivatives_it[*elem_filter_val * nb_quadrature_points + q];

    /// compute @f$\mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
    UInt bt_d_b_size = dim * nb_nodes_per_element;

    Array<Real> * bt_d_b = new Array<Real > (nb_element * nb_quadrature_points,
            bt_d_b_size * bt_d_b_size,
            "B^t*D*B");

    Matrix<Real> B(cauchy_matrix_size, bt_d_b_size);
    Matrix<Real> Bt_D(bt_d_b_size, cauchy_matrix_size);

    shapes_derivatives_filtered_it = shapes_derivatives_filtered->begin(spatial_dimension, nb_nodes_per_element);

    Array<Real>::iterator< Matrix<Real> > Bt_D_B_it = bt_d_b->begin(bt_d_b_size,
            bt_d_b_size);

    Array<Real>::iterator< Matrix<Real> > D_it = cauchy_stress_matrix->begin(cauchy_matrix_size,
            cauchy_matrix_size);
    Array<Real>::iterator< Matrix<Real> > D_end = cauchy_stress_matrix->end(cauchy_matrix_size,
            cauchy_matrix_size);


    for (; D_it != D_end; ++D_it, ++Bt_D_B_it, ++shapes_derivatives_filtered_it) {
        Matrix<Real> & D = *D_it;
        Matrix<Real> & Bt_D_B = *Bt_D_B_it;

        transferBMatrixToBNL< dim > (*shapes_derivatives_filtered_it, B, nb_nodes_per_element);
        Bt_D.mul < true, false > (B, D);
        Bt_D_B.mul < false, false > (Bt_D, B);
    }

    delete cauchy_stress_matrix;
    delete shapes_derivatives_filtered;

    /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
    Array<Real> * K_e = new Array<Real > (nb_element,
            bt_d_b_size * bt_d_b_size,
            "K_e");

    model->getFEM().integrate(*bt_d_b, *K_e,
            bt_d_b_size * bt_d_b_size,
            type, ghost_type,
            &elem_filter);

    delete bt_d_b;

    model->getFEM().assembleMatrix(*K_e, K, spatial_dimension, type, ghost_type, &elem_filter);
    delete K_e;

    AKANTU_DEBUG_OUT();
}

template<UInt dim>
void MaterialFiniteDeformation::assembleStiffnessMatrixL2(const ElementType & type,
        GhostType ghost_type) {
    AKANTU_DEBUG_IN();

    SparseMatrix & K = const_cast<SparseMatrix &> (model->getStiffnessMatrix());

    const Array<Real> & shapes_derivatives = model->getFEM().getShapesDerivatives(type, ghost_type);

    Array<UInt> & elem_filter = element_filter(type, ghost_type);
    Array<Real> & strain_vect = strain(type, ghost_type);

    UInt nb_element = elem_filter.getSize();
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    UInt nb_quadrature_points = model->getFEM().getNbQuadraturePoints(type, ghost_type);

    strain_vect.resize(nb_quadrature_points * nb_element);

    model->getFEM().gradientOnQuadraturePoints(model->getDisplacement(), strain_vect,
            dim, type, ghost_type, &elem_filter);

    UInt tangent_size = getTangentStiffnessVoigtSize(dim);

    Array<Real> * tangent_stiffness_matrix =
            new Array<Real > (nb_element*nb_quadrature_points, tangent_size * tangent_size,
            "tangent_stiffness_matrix");

    tangent_stiffness_matrix->clear();

    computeTangentModuli(type, *tangent_stiffness_matrix, ghost_type);


    Array<Real> * shapes_derivatives_filtered = new Array<Real > (nb_element * nb_quadrature_points,
            dim * nb_nodes_per_element,
            "shapes derivatives filtered");


    Array<Real>::const_iterator< Matrix<Real> > shapes_derivatives_it = shapes_derivatives.begin(spatial_dimension,
            nb_nodes_per_element);

    Array<Real>::iterator< Matrix<Real> > shapes_derivatives_filtered_it = shapes_derivatives_filtered->begin(spatial_dimension,
            nb_nodes_per_element);
    UInt * elem_filter_val = elem_filter.storage();
    for (UInt e = 0; e < nb_element; ++e, ++elem_filter_val)
        for (UInt q = 0; q < nb_quadrature_points; ++q, ++shapes_derivatives_filtered_it)
            *shapes_derivatives_filtered_it = shapes_derivatives_it[*elem_filter_val * nb_quadrature_points + q];

    /// compute @f$\mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
    UInt bt_d_b_size = dim * nb_nodes_per_element;

    Array<Real> * bt_d_b = new Array<Real > (nb_element * nb_quadrature_points,
            bt_d_b_size * bt_d_b_size,
            "B^t*D*B");

    Matrix<Real> B(tangent_size, dim * nb_nodes_per_element);
    Matrix<Real> B2(tangent_size, dim * nb_nodes_per_element);
    Matrix<Real> Bt_D(dim * nb_nodes_per_element, tangent_size);

    shapes_derivatives_filtered_it = shapes_derivatives_filtered->begin(spatial_dimension, nb_nodes_per_element);

    Array<Real>::iterator< Matrix<Real> > Bt_D_B_it = bt_d_b->begin(dim*nb_nodes_per_element,
            dim * nb_nodes_per_element);

    Array<Real>::iterator< Matrix<Real> > grad_u_it = strain_vect.begin(dim, dim);

    Array<Real>::iterator< Matrix<Real> > D_it = tangent_stiffness_matrix->begin(tangent_size,
            tangent_size);
    Array<Real>::iterator< Matrix<Real> > D_end = tangent_stiffness_matrix->end(tangent_size,
            tangent_size);


    for (; D_it != D_end; ++D_it, ++Bt_D_B_it, ++shapes_derivatives_filtered_it, ++grad_u_it) {
        Matrix<Real> & grad_u = *grad_u_it;
        Matrix<Real> & D = *D_it;
        Matrix<Real> & Bt_D_B = *Bt_D_B_it;

        transferBMatrixToBL1<dim > (*shapes_derivatives_filtered_it, B, nb_nodes_per_element);
        transferBMatrixToBL2< dim > (*shapes_derivatives_filtered_it, grad_u, B2, nb_nodes_per_element);
        B += B2;
        Bt_D.mul < true, false > (B, D);
        Bt_D_B.mul < false, false > (Bt_D, B);
    }

    delete tangent_stiffness_matrix;
    delete shapes_derivatives_filtered;

    /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
    Array<Real> * K_e = new Array<Real > (nb_element,
            bt_d_b_size * bt_d_b_size,
            "K_e");

    model->getFEM().integrate(*bt_d_b, *K_e,
            bt_d_b_size * bt_d_b_size,
            type, ghost_type,
            &elem_filter);

    delete bt_d_b;

    model->getFEM().assembleMatrix(*K_e, K, spatial_dimension, type, ghost_type, &elem_filter);
    delete K_e;

    AKANTU_DEBUG_OUT();
}

void MaterialFiniteDeformation::SetCauchyStressMatrix(const ElementType & type, Array<Real> & stress_matrix, GhostType ghost_type) {
    AKANTU_DEBUG_IN();

    Array<Real> & stress_vect = stress(type, ghost_type);
    //Array<Real> & stress_vect = stress_at_t(type, ghost_type);

    UInt cauchy_matrix_size = getCauchyStressMatrixSize(spatial_dimension);

    UInt dim_stress_tensor = spatial_dimension;

    Array<Real>::iterator< Matrix<Real> > S_it = stress_vect.begin(dim_stress_tensor,
            dim_stress_tensor);
    Array<Real>::iterator< Matrix<Real> > S_end = stress_vect.end(dim_stress_tensor,
            dim_stress_tensor);

    Array<Real>::iterator< Matrix<Real> > Stress_matrix_it = stress_matrix.begin(cauchy_matrix_size,
            cauchy_matrix_size);


    for (; S_it != S_end; ++S_it, ++Stress_matrix_it) {
        Matrix<Real> & S = *S_it;
        Matrix<Real> & S_Matrix = *Stress_matrix_it;
        switch (spatial_dimension) {
            case 1:
                SetCauchyStressMatrix < 1 > (S, S_Matrix);
                break;
            case 2:
                SetCauchyStressMatrix < 2 > (S, S_Matrix);
                break;
            case 3:
                SetCauchyStressMatrix < 3 > (S, S_Matrix);
                break;
        }
    }

    AKANTU_DEBUG_OUT();
}

void MaterialFiniteDeformation::SetCauchyStressArray(const ElementType & type, Array<Real> & stress_vect, GhostType ghost_type) {
    AKANTU_DEBUG_IN();

    Array<Real> & stress_matrix = stress(type, ghost_type);

    UInt cauchy_vector_size = getCauchyStressArraySize(spatial_dimension);

    UInt dim_stress_tensor = 3;

    Array<Real>::iterator< Matrix<Real> > S_it = stress_matrix.begin(dim_stress_tensor,
            dim_stress_tensor);
    Array<Real>::iterator< Matrix<Real> > S_end = stress_matrix.end(dim_stress_tensor,
            dim_stress_tensor);

    Array<Real>::iterator< Matrix<Real> > stress_vect_it = stress_vect.begin(cauchy_vector_size, 1);


    for (; S_it != S_end; ++S_it, ++stress_vect_it) {
        Matrix<Real> & S = *S_it;
        Matrix<Real> & S_vect = *stress_vect_it;
        switch (spatial_dimension) {
            case 1:
                SetCauchyStressMatrix < 1 > (S, S_vect);
                break;
            case 2:
                SetCauchyStressMatrix < 2 > (S, S_vect);
                break;
            case 3:
                SetCauchyStressMatrix < 3 > (S, S_vect);
                break;
        }
    }

    AKANTU_DEBUG_OUT();
}

__END_AKANTU__
