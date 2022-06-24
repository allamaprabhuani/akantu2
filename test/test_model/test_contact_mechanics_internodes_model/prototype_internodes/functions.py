#!/usr/bin/env python
# coding: utf-8

import akantu as aka
import numpy as np
import scipy as sp
import scipy.sparse.linalg as spla
from scipy.sparse.linalg import spsolve

def nodes_to_dofs(nodes):
    """Returns dof numbers for nodes in 2D. """
    return np.column_stack([2*nodes, 2*nodes+1]).reshape([-1, 1])

def assemble_interface_masses(model, dofs1i, dofs2i):
    """Assembles mass matrices for interfaces.

    :param model: akantu solid mechanics model
    :param dofs1i: dofs of body 1 interface (master)
    :param dofs2i: dofs of body 2 interface (slave)
    :returns: sparse matrices of interface mass
    """
    # assemble and store mass
    model.assembleMass()
    M_aka = model.getDOFManager().getMatrix('M')
    M = sp.sparse.lil_matrix(aka.AkantuSparseMatrix(M_aka))

    # select only the dofs of the interface
    M1i = M[np.ix_(dofs1i, dofs1i)]
    M2i = M[np.ix_(dofs2i, dofs2i)]

    M1i = sp.sparse.csc_matrix(M1i)
    M2i = sp.sparse.csc_matrix(M2i)

    return M1i, M2i

def assemble_Bs(M1i, M2i, R12, R21, indx1i, indx2i,
        nb_free_dofs, nb_constraint_dofs):
    """Assembles B and B_tilde for internodes matrix.

    :param M1i: sparse matrix interface mass 1 
    :param M2i: sparse matrix interface mass 2 
    :param R12: interpolation matrix master to slave
    :param R21: interpolation matrix slave to master
    :param indx1i: global indices of interface 1 (master)
    :param indx2i: global indices of interface 2 (slave)
    :param nb_free_dofs: number of non blocked dofs 
    :param nb_constraint_dofs: number of master dofs which act as constraint 
    :returns: sparse blocke matrices for internodes formulation 
    """
    B = sp.sparse.lil_matrix(np.zeros((nb_free_dofs, nb_constraint_dofs)))
    B_tilde = sp.sparse.lil_matrix(np.zeros((nb_constraint_dofs, nb_free_dofs)))
    C = sp.sparse.lil_matrix(np.zeros((nb_constraint_dofs, nb_constraint_dofs)))

    B[indx1i, :] = - M1i
    B[indx2i, :] = M2i * R21

    B_tilde[:, indx1i] = sp.sparse.eye(nb_constraint_dofs)
    B_tilde[:, indx2i] = - R12

    return B, B_tilde, C

def assemble_A_explicit(K_free, B, B_tilde, C):
    """Explicitly assembles internodes matrix.

    :param K_free: stiffness matrix with non blocked dofs
    :param B: B matrix of internodes formulation
    :param B_tilde: B_tilde matrix of internodes formulation
    :param C: all zero matrix of internodes formulation
    """
    A1 = sp.sparse.hstack([K_free, B])
    A2 = sp.sparse.hstack([B_tilde, C])

    A = sp.sparse.vstack([A1, A2])
    A = sp.sparse.csr_matrix(A)

    return A


def assemble_b(f_free, R12_normal, K, displacements, free_dofs, blocked_dofs, 
        positions1i, positions2i, nb_constraint_dofs, rescaling):
    """Assemble right hand side of Ax=b.

    :param f_free: force vector of free dofs 
    :param R12_normal: nodal interpolation matrix master to slave
    :param K: stiffness matrix
    :param displacements: displacements vector
    :param free_dofs: dofs numbers of non blocked dofs
    :param blocked_dofs: dofs numbers of blocked dofs
    :param positions1i: initial positions of interface 1 nodes
    :param positions2i: initial positions of interface 2 nodes
    :param nb_constraint_dofs: number of master dofs which act as constraint
    :param rescaling: rescaling stiffness matrix
    """
    nb_free_dofs = len(free_dofs)
    b = np.zeros(nb_free_dofs + nb_constraint_dofs)

    # Dirichlet displacements
    K_blocked = K[free_dofs, :]
    K_blocked = K_blocked[:, blocked_dofs]
    dirichlet = K_blocked.dot(displacements[blocked_dofs])

    b[0:nb_free_dofs] = 1/rescaling * (f_free - dirichlet) 
    b[nb_free_dofs:] = (R12_normal.dot(positions2i) - positions1i).ravel()

    return b

def solve_direct(A, b, positions, displacements, free_dofs,
        nb_constraint_dofs, nb_dofs, rescaling):
    """Direct solve of internodes equation Ax=b.

    :param A: sparse matrix internodes
    :param b: right hand side vector
    :param positions: initial positions of all nodes
    :param displacements: displacements vector
    :param free_dofs: dofs numbers of non blocked dofs
    :param nb_constraint_dofs: number of master dofs which act as constraint
    :param nb_dofs: total number of dofs
    :param rescaling: rescaling stiffness matrix
    :returns: new positions, displacements and lambdas resulting from constraints
    """
    nb_free_dofs = len(free_dofs)

    x = spsolve(A, b)
    displacements[free_dofs] = x[:nb_free_dofs]

    positions_new = positions + displacements.reshape([-1, 2])
    lambdas = rescaling * x[nb_free_dofs:nb_free_dofs+nb_constraint_dofs].reshape([-1, 2])

    return positions_new, displacements, lambdas

# all functions below are needed for iterative solve 
def solve_iterative(A_op_spla, b, M_inv_op_spla, positions, displacements,
       free_dofs, ordering, nb_constraint_dofs, nb_dofs, rescaling):
    """Iterative solve of internodes equation Ax=b.

    :param A_op_spla: scipy linear operator for A
    :param b: right hand side vector
    :param M_inv_op_spla: scipy linear operator for M
    :param positions: initial positions of all nodes
    :param displacements: displacements vector
    :param free_dofs: dofs numbers of non blocked dofs
    :param ordering: reordered stiffness matrix due to cholesky factorization
    :param nb_constraint_dofs: number of master dofs which act as constraint
    :param nb_dofs: total number of dofs
    :param rescaling: rescaling stiffness matrix
    :returns: new positions and lambdas resulting from constraints
    """
    nb_free_dofs = len(free_dofs)

    x, exitCode = spla.gmres(A_op_spla, b, tol=1e-07, restart=5, maxiter=3, M=M_inv_op_spla)

    disp_order = x[:nb_free_dofs]
    inv_ordering = np.argsort(ordering)
    displacements[free_dofs] = disp_order[inv_ordering]

    positions_new = positions + displacements.reshape([-1, 2])
    lambdas = rescaling * x[nb_free_dofs:nb_free_dofs+nb_constraint_dofs].reshape([-1, 2])

    return positions_new, displacements, lambdas

def linearoperator_A(x, K_free_order, B, B_tilde, nb_free_dofs):
    """Creates linear operator for internodes matrix.

    :param K_free_order: ordered stiffness matrix with non blocked dofs
    :param B: B matrix of internodes formulation
    :param B_tilde: B_tilde matrix of internodes formulation
    :param nb_free_dofs: number of non blocked dofs 
    :param nb_constraint_dofs: number of master dofs which act as constraint 
    """

    # A @ [x1; x2] = [y1; y2]
    x1 = x[:nb_free_dofs]
    x2 = x[nb_free_dofs:]

    y1 = K_free_order.dot(x1) + B.dot(x2)
    y2 = B_tilde.dot(x1)

    y = np.concatenate((y1, y2))

    return y 

def prepare_precond(M1i, M2i, R12, R21, L22, L22L22_T, indx1i, indx2i,
        nb_boundary_dofs):
    """Prepares matrix V needed for construction of preconditioner.

    :param M1i: sparse matrix interface mass 1 
    :param M2i: sparse matrix interface mass 2 
    :param R12: interpolation matrix master to slave
    :param R21: interpolation matrix slave to master
    :param L22: block 22 of cholesky factorization
    :param L22L22_T: precomputed L_22 * L_22^T
    :param indx1i: indices of interface 1 (master)
    :param indx2i: indices of interface 2 (slave)
    :param nb_boundary_dofs: number of initial boundary dofs
    :return: matrix V needed for preconditioner
    """
    indxi = np.append(indx1i, indx2i)
 
    # assemble Y1*Y2^T
    Q1 = -M2i.dot((R21.dot(spla.inv(M1i))))
    Q2 = -R12.T

    Y1 = sp.sparse.vstack([sp.sparse.identity(len(indx1i)), Q1])
    Y2 = sp.sparse.vstack([sp.sparse.identity(len(indx1i)), Q2])
    Y1Y2_T = Y1.dot(Y2.T)
    
    Z = sp.sparse.lil_matrix((nb_boundary_dofs, nb_boundary_dofs))
    Z[np.ix_(indxi, indxi)] = Y1Y2_T

    # assemble V = L22*L22^T - alpha*Y1*Y2^T
    alpha = 7.922 # fixed for preconditioner
    V = L22L22_T - alpha*Z

    return V

def linearoperator_precond(y, L22, V, B, M1i, chol_factor,
        nb_free_dofs, nb_boundary_dofs):
    """Creates a linear operator for inverse of preconditioner M^-1
    , where MAx=Mb and Ax=y.

    :param y: solution to solve for
    :param L22: submatrice of chelosky for perturbed stiffness matrix
    :param V: prepared matrix for preconditioner
    :param B: B matrix of internodes formulation
    :param M1i: sparse matrix interface mass 1
    :param chol_factor: sksparse factor to do operations with L
    :param nb_free_dofs: number of non blocked dofs 
    :param nb_boundary_dofs: indices of initial boundary interface dofs
    """
    boundary_start = nb_free_dofs - nb_boundary_dofs 

    # M @ [x1; x2] = [y1; y2]
    y1 = y[:nb_free_dofs]
    y2 = y[nb_free_dofs:]

    alpha = 7.922 # fixed for preconditioner
    x2 = spla.spsolve(M1i, -alpha*y2)

    y1_tilde = y1 - 2*B.dot(x2)
    x1 = chol_factor.solve_L(y1_tilde, use_LDLt_decomposition=False)

    x1a = x1[:boundary_start]

    x1b = L22.dot(x1[boundary_start:])

    # solving for V could be done
    # with GMRES as well + preconditiong
    x1b = spla.spsolve(V, x1b)
    x1b = L22.T.dot(x1b)

    x1 = np.concatenate((x1a, x1b))
    x1 = chol_factor.solve_Lt(x1, use_LDLt_decomposition=False)

    x = np.concatenate((x1, x2))

    return x

def solve_iterative(A_op_spla, b, M_inv_op_spla, positions, displacements,
       free_dofs, ordering, nb_constraint_dofs, nb_dofs, rescaling):
    """Iterative solve of internodes equation Ax=b.

    :param A_op_spla: scipy linear operator for A
    :param b: right hand side vector
    :param M_inv_op_spla: scipy linear operator for B
    :param positions: initial positions of all nodes
    :param displacements: displacements vector
    :param free_dofs: dofs numbers of non blocked dofs
    :param ordering: reordered stiffness matrix due to cholesky factorization
    :param nb_constraint_dofs: number of master dofs which act as constraint
    :param nb_dofs: total number of dofs
    :param rescaling: rescaling stiffness matrix
    :returns: new positions and lambdas resulting from constraints
    """
    nb_free_dofs = len(free_dofs)

    x, exitCode = spla.gmres(A_op_spla, b, tol=1e-07, restart=5, maxiter=3, M=M_inv_op_spla)

    disp_order = x[:nb_free_dofs]
    inv_ordering = np.argsort(ordering)
    displacements[free_dofs] = disp_order[inv_ordering]

    positions_new = positions + displacements.reshape([-1, 2])
    lambdas = rescaling * x[nb_free_dofs:nb_free_dofs+nb_constraint_dofs].reshape([-1, 2])

    return positions_new, displacements, lambdas
