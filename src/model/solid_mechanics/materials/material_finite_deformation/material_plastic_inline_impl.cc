/**
 * @file   material_plastic_inline_impl.cc
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 * @date   Wed Aug 04 10:58:42 2010
 *
 * @brief  Implementation of the inline functions of the material elastic
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

#include <cmath>
#include "material_plastic.hh"


/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void MaterialPlastic<dim>::computeDeltaStressOnQuad(const Matrix<Real> & grad_u, const Matrix<Real> & grad_delta_u,
        Matrix<Real> & delta_S){

  /*Real J = 1.0;
    switch (dim){
        case 2:
            J = Math::det2(F.storage());
            break;
        case 3:
            J = Math::det3(F.storage());
            break;
    }
    if (std::abs(J) > Math::getTolerance()) {*/

        /*
         * Updated Lagrangian approach
         *
         * ^{t+\Delta t}S_{ij}:         2nd Piola Kirchhoff tensor at t+\Delta t
         * ^{t}\sigma_{ij}:             Cauchy stress tensor at t
         * \Delta S_{ij}:               Increment of the
         *
         * ^{t+\Delta t}S_{ij}=^{t}\sigma_{ij}+\Delta S_{ij}2nd Piola Kirchhoff tensor
         *
         */
        //cauchy += S;

        /*
         * ^{t+\Delta t}\sigma_{ij}:    Cauchy stress tensor at t+\Delta t
         * F:                           Deformation Gradient
         *
         * ^{t+\Delta t}\sigma_{ij}=1/det(F) F_ir ^{t+\Delta t}S_{rs} F_{sj}
         *
         */

        /* Matrix<Real> Aux(3, 3);
        Matrix<Real> Aux2(3, 3);

        Aux.mul < false, false > (F, cauchy);
        Aux2.mul < false, true > (Aux, F, 1.0 / J);

        //cauchy.mul < false, false > (F, cauchy);

        //cauchy.mul < false, true > (cauchy, F, 1.0/J );

        for (UInt i = 0; i < dim; i++)
            for (UInt j = 0; j < dim; j++)
                cauchy(i, j) = Aux2(i, j);

        //std::cout << cauchy << std::endl;

    } else
    cauchy.clear();*/

}

template<UInt dim>
inline void MaterialPlastic<dim>::computeStressOnQuad(Matrix<Real> & grad_u,
                                                           Matrix<Real> & sigma) {
    //Neo hookean book
    Matrix<Real> F(dim, dim);
    Matrix<Real> C(dim, dim);//Right green
    Matrix<Real> Cminus(dim, dim);//Right green


    this->template gradUToF<dim > (grad_u, F);
    this->rightCauchy(F, C);
    Real J = 1.0;
    switch (dim) {
        case 2:
            Math::inv2(C.storage(), Cminus.storage());
            J = Math::det2(F.storage());
            break;
        case 3:
            Math::inv3(C.storage(), Cminus.storage());
            J = Math::det3(F.storage());
            break;
    }


    for (UInt i = 0; i < dim; ++i)
        for (UInt j = 0; j < dim; ++j)
            sigma(i, j) = (i == j) *  mu  + (lambda * log(J) - mu) * Cminus(i, j);
    //Neo hookean book
    /*Matrix<Real> F(dim, dim);
    Matrix<Real> C(dim, dim);//Right green


    this->template gradUToF<dim > (grad_u, F);
    this->rightCauchy(F, C);
    Real J = 1.0;
    switch (dim) {
        case 2:
            J = Math::det2(F.storage());
            break;
        case 3:
            J = Math::det3(F.storage());
            break;
    }

    Real trace_green = 0.5 * ( C.trace() - dim);

    for (UInt i = 0; i < dim; ++i)
        for (UInt j = 0; j < dim; ++j)
          sigma(i, j) = (i == j) * trace_green * lambda/J  + ( mu - lambda * log(J) ) / J
          * (0.5 * ( C(i, j) - (i==j) ) + 0.5 * ( C(j, i) - (j==i) ) );*/

}

template<UInt dim>
inline void MaterialPlastic<dim>::computePiolaKirchhoffOnQuad(const Matrix<Real> & E,
        Matrix<Real> & S) {

    Real trace = E.trace(); /// trace = (\nabla u)_{kk}

    /// \sigma_{ij} = \lambda * (\nabla u)_{kk} * \delta_{ij} + \mu * (\nabla u_{ij} + \nabla u_{ji})
    for (UInt i = 0; i < dim; ++i)
        for (UInt j = 0; j < dim; ++j)
            S(i, j) = (i == j) * lambda * trace + 2.0 * mu * E(i, j);

}


/**************************************************************************************/
/*  Computation of the potential energy for a this neo hookean material */
template<UInt dim>
inline void MaterialPlastic<dim>::computePotentialEnergyOnQuad(Matrix<Real> & grad_u,
                                                               Matrix<Real> & sigma,
                                                               Real & epot){
    Matrix<Real> F(dim, dim);
    Matrix<Real> C(dim, dim);//Right green

    this->template gradUToF<dim > (grad_u, F);
    this->rightCauchy(F, C);
    Real J = 1.0;
    switch (dim) {
        case 2:
            J = Math::det2(F.storage());
            break;
        case 3:
            J = Math::det3(F.storage());
            break;
    }

    epot=0.5*lambda*pow(log(J),2.)+ mu * (-log(J)+0.5*(C.trace()-dim));

}

/*template<UInt spatial_dimension>
inline void MaterialPlastic<spatial_dimension>::updateStressOnQuad(const Matrix<Real> & sigma,
        Matrix<Real> & cauchy_sigma) {

    for (UInt i = 0; i < spatial_dimension; ++i)
        for (UInt j = 0; j < spatial_dimension; ++j)
            cauchy_sigma(i, j) += sigma(i, j);

}*/

/* -------------------------------------------------------------------------- */
/*template<>
inline void MaterialPlastic < 1 > ::computeStressOnQuad(const Matrix<Real> & F, const Matrix<Real> & S,
        Matrix<Real> & cauchy) {
    cauchy(0, 0) = E * F(0, 0);
}*/

/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void MaterialPlastic<dim>::computeTangentModuliOnQuad(Matrix<Real> & tangent, Matrix<Real> & grad_u) {

    //Neo hookean book
    UInt cols = tangent.cols();
    UInt rows = tangent.rows();
    Matrix<Real> F(dim, dim);
    Matrix<Real> C(dim, dim);
    Matrix<Real> Cminus(dim, dim);
    this->template gradUToF<dim > (grad_u, F);
    this->rightCauchy(F, C);
    Real J = 1.0;
    switch (dim){
        case 2:
            Math::inv2(C.storage(), Cminus.storage());
            J = Math::det2(F.storage());
            break;
        case 3:
            Math::inv3(C.storage(), Cminus.storage());
            J = Math::det3(F.storage());
            break;
    }

    for (UInt m = 0; m < rows; m++) {
        UInt i, j;
        if (m < dim) {
            i = m;
            j = m;
        } else {
            if (dim == 3) {
                if (m == 3) {
                    i = 1;
                    j = 2;
                } else if (m == 4) {
                    i = 0;
                    j = 2;
                } else if (m == 5) {
                    i = 0;
                    j = 1;
                }
            } else if (dim == 2) {
                if (m == 2) {
                    i = 0;
                    j = 1;
                }
            }
        }

        for (UInt n = 0; n < cols; n++) {
            UInt k, l;
            if (n < dim) {
                k = n;
                l = n;
            } else {
                if (dim == 3) {
                    if (n == 3) {
                        k = 1;
                        l = 2;
                    } else if (n == 4) {
                        k = 0;
                        l = 2;
                    } else if (n == 5) {
                        k = 0;
                        l = 1;
                    }
                } else if (dim == 2) {
                    if (n == 2) {
                        k = 0;
                        l = 1;
                    }
                }
            }

            //book Bathe
            /*            tangent(m, n) = J*( lambda * Cminus(i, j) * Cminus(k, l) +
                          mu * (Cminus(i, k) * Cminus(j, l) + Cminus(i, l) * Cminus(k, j)));*/

            //Linear elastic
            /*            tangent(m, n) = lambda * (i==j) *  (k==l) +
                          mu * ((i==k) * (j==l) + (i==l) *(k==j));*/

            //book belytchko
            tangent(m, n) = lambda * Cminus(i, j) * Cminus(k, l) +
              (mu - lambda * log(J)) * (Cminus(i, k) * Cminus(j, l) + Cminus(i, l) * Cminus(k, j));

        }
    }

    /*
    UInt cols = tangent.cols();
    UInt rows = tangent.rows();
    Matrix<Real> F(dim, dim);
    Matrix<Real> C(dim, dim);
    this->template gradUToF<dim > (grad_u, F);
    this->rightCauchy(F, C);
    Real J = 1.0;
    switch (dim){
        case 2:
            J = Math::det2(F.storage());
            break;
        case 3:
            J = Math::det3(F.storage());
            break;
    }

    Real mu_NH = (mu - lambda * log(J))/J;
    Real lambda_NH = lambda/J;

    for (UInt m = 0; m < rows; m++) {
        UInt i, j;
        if (m < dim) {
            i = m;
            j = m;
        } else {
            if (dim == 3) {
                if (m == 3) {
                    i = 0;
                    j = 1;
                } else if (m == 4) {
                    i = 1;
                    j = 2;
                } else if (m == 5) {
                    i = 2;
                    j = 0;
                }
            } else if (dim == 2) {
                if (m == 2) {
                    i = 0;
                    j = 1;
                }
            }
        }

        for (UInt n = 0; n < cols; n++) {
            UInt k, l;
            if (n < dim) {
                k = n;
                l = n;
            } else {
                if (dim == 3) {
                    if (n == 3) {
                        k = 0;
                        l = 1;
                    } else if (n == 4) {
                        k = 1;
                        l = 2;
                    } else if (n == 5) {
                        k = 2;
                        l = 0;
                    }
                } else if (dim == 2) {
                    if (n == 2) {
                        k = 0;
                        l = 1;
                    }
                }
            }

            tangent(m, n) = lambda_NH * (i==j) * (k==l) + mu_NH * ((i==k) * (j==l) + (i==l) * (k==j));

        }
        }*/

}

/* -------------------------------------------------------------------------- */
/*template<>
inline void MaterialPlastic < 1 > ::computeTangentModuliOnQuad(Matrix<Real> & tangent) {
    tangent(0, 0) = E;
}*/

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline Real MaterialPlastic<spatial_dimension>::getStableTimeStep(Real h,
        __attribute__((unused)) const Element & element) {
    return (h / getPushWaveSpeed());
}
