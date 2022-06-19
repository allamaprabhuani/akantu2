#!/usr/bin/env python
# coding: utf-8

import akantu as aka
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from scipy.sparse.linalg import LinearOperator
from sksparse.cholmod import cholesky

import sys
sys.path.append("..")
from prototype_internodes.functions import * 
from prototype_internodes.functions_contact_probl import * 
from prototype_internodes.init_model import init_model


# example
def main():
    mesh_file = 'contact.msh'
    material_file = 'material.dat'

    aka.parseInput(material_file)
    spatial_dimension = 2

    mesh = aka.Mesh(spatial_dimension)
    mesh.read(mesh_file)

    # initialize model
    model = aka.SolidMechanicsModel(mesh)
    model.initFull(_analysis_method=aka._implicit_dynamic)

    # boundary conditions
    displacements = np.zeros(mesh.getNbNodes()*spatial_dimension)

    model.applyBC(aka.FixedValue(0., aka._x), 'lower_bottom')
    model.applyBC(aka.FixedValue(0., aka._y), 'lower_bottom')

    # Dirichlet
    # nodes_top = mesh.getElementGroup('upper_top').getNodeGroup().getNodes().ravel()
    # displacements = displacements.reshape([-1, 2])
    # displacements[nodes_top, 1] = -0.1
    # displacements = displacements.ravel()

    # Neumann (K is not invertible)
    traction = np.zeros(spatial_dimension)
    traction[1] = -1e9
    model.applyBC(aka.FromTraction(traction), "upper_top")

    # init and solve
    model, data = init_model(model, mesh, mesh_file, material_file, spatial_dimension,
            displacements, "upper_bottom", "lower_top")

    solve_step_iterative(model, data, nb_max_iter=10, plot=True)

def solve_step_iterative(model, data, nb_max_iter=10, plot=False):
    #---------- assemble stiffness matrices: K, K_free ----------
    model.assembleStiffnessMatrix()

    K_aka = model.dof_manager.getMatrix("K")
    K = sp.sparse.csc_matrix(aka.AkantuSparseMatrix(K_aka))

    # K for all non blocked dofs
    rescaling = 30e9 # with E
    K_free = K[np.ix_(data['free_dofs'], data['free_dofs'])] / rescaling

    # calculate perturbed K free
    epsilon = 1e-8
    K_tilde = K_free + epsilon*sp.sparse.identity(data['nb_free_dofs'])

    #---------- ordering and indexing of K and K_tilde ----------
    sorter_free_dofs = np.argsort(data['free_dofs'])
    indx1b = sorter_free_dofs[np.searchsorted(data['free_dofs'], data['dofs1b'], sorter=sorter_free_dofs)]
    indx2b = sorter_free_dofs[np.searchsorted(data['free_dofs'], data['dofs2b'], sorter=sorter_free_dofs)]
    indxb = np.append(indx1b, indx2b)

    # ordering for cholesky of K_tilde with scikit-sparse
    factor = cholesky(K_tilde, beta=0)

    permut = factor.P()
    permut_inv = np.argsort(permut)

    indx_permut = np.setdiff1d(permut, indxb, assume_unique=True)
    ordering = np.append(indx_permut, indxb)
    ordering_inv = np.argsort(ordering)

    K_tilde_order = K_tilde[np.ix_(ordering, ordering)]
    K_free_order = K_free[np.ix_(ordering, ordering)]

    #---------- cholesky decomposition of K_tilde with scikit-sparse ----------
    # this is not very efficient in sksparse
    chol_factor = cholesky(K_tilde_order, beta=0, ordering_method="natural")
    L = chol_factor.L()
    L = sp.sparse.csr_matrix(L)

    # explicit calculations of blocks
    boundary_start = data['nb_free_dofs'] - len(indxb)
    L22 = L[boundary_start:, boundary_start:] 
    L22L22_T = L22.dot(L22.T)

    #---------- assemble external forces: f, f_ordering ----------
    f = model.getExternalForce().ravel()

    # f for all non blocked dofs
    f_free = f[data['free_dofs']]
    f_free_order = f_free[ordering]

    #---------- initialize interface ----------
    # 1i denotes interface 1, 2i denotes inteface 2

    # take all boundary initialy as interface nodes
    # possible selection with radius calculation
    nodes1i = data['nodes1b']
    positions1i = data['positions1b']
    nodes2i = data['nodes2b'] 
    positions2i = data['positions2b']

    if plot:
        plt.figure()
        plt.triplot(data['positions'][:, 0], data['positions'][:, 1], data['connectivity'])
        plt.title('mesh')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.axis('scaled')
        plt.title('Initial')
        # plt.ylim([0.9, 1.1])
        plt.show()

    #---------- internodes iterations ----------
    # run until converged or nb_max_iter is attained
    for i in range(nb_max_iter):
        print('--- iteration', i+1, '---')

        # select nodes belonging to interface
        nodes1i, nodes2i, positions1i, positions2i, radiuses1, radiuses2 = find_contact_nodes(nodes1i, nodes2i,
                positions1i, positions2i)

        dofs1i = nodes_to_dofs(nodes1i).ravel()
        dofs2i = nodes_to_dofs(nodes2i).ravel()
        nb_constraint_dofs = len(dofs1i)

        # local index of the interface among the boundary dofs
        sorter_boundary_dofs = np.argsort(data['boundary_dofs'])
        indx1i_loc = sorter_boundary_dofs[np.searchsorted(data['boundary_dofs'], dofs1i, sorter=sorter_boundary_dofs)]
        indx2i_loc = sorter_boundary_dofs[np.searchsorted(data['boundary_dofs'], dofs2i, sorter=sorter_boundary_dofs)]

        # global index of the interface among all dofs
        boundary_start = data['nb_free_dofs'] - len(indxb)
        indx1i = indx1i_loc + boundary_start
        indx2i = indx2i_loc + boundary_start

        # subassemble matrices
        R12_normal, R21_normal, R12, R21 = assemble_Rijs(positions1i, positions2i, radiuses1, radiuses2)    
        M1i, M2i = assemble_interface_masses(model, dofs1i, dofs2i)
        B, B_tilde, C = assemble_Bs(M1i, M2i, R12, R21, indx1i, indx2i, data['nb_free_dofs'], nb_constraint_dofs)

        b = assemble_b(f_free_order, R12_normal, K, data['displacements'], data['free_dofs'],
                data['blocked_dofs'], positions1i, positions2i, nb_constraint_dofs, rescaling)

        A_size = data['nb_free_dofs'] + nb_constraint_dofs
        A_op = lambda x: linearoperator_A(x, K_tilde_order, B, B_tilde, data['nb_free_dofs'])
        A_op_spla = LinearOperator((A_size, A_size), matvec=A_op)

        # prepare preconditioner
        V = prepare_precond(M1i, M2i, R12, R21, L22, L22L22_T,
                indx1i_loc, indx2i_loc, data['nb_boundary_dofs'])
        M_inv_op = lambda y: linearoperator_precond(y, L22, V, B, M1i, chol_factor,
                data['nb_free_dofs'], data['nb_boundary_dofs'])
        M_inv_op_spla = LinearOperator((A_size, A_size), matvec=M_inv_op)

        # solve
        positions_new, displacements, lambdas = solve_iterative(A_op_spla, b, M_inv_op_spla,
                data['positions'], data['displacements'], data['free_dofs'], ordering,
                nb_constraint_dofs, data['nb_dofs'], rescaling)

        if plot:
            plt.figure()
            plt.triplot(positions_new[:, 0], positions_new[:, 1], data['connectivity'])
            plt.title('mesh')
            plt.xlabel('x')
            plt.ylabel('y')
            plt.axis('scaled')
            plt.title('Iteration ' + str(i+1))
            # plt.ylim([0.9, 1.1])
            plt.show()

        # add or remove nodes
        nodes1i, nodes2i, diff_nb_nodes1i, diff_nb_nodes2i = remove_traction(positions_new,
                data['connectivity1b'], data['connectivity2b'], data['connectivity1b_body'], data['connectivity2b_body'],
                nodes1i, nodes2i, data['nodes1b'], data['nodes2b'], lambdas, R12, R21)

        positions1i = data['positions'][nodes1i, :]
        positions2i = data['positions'][nodes2i, :]

        if np.abs(diff_nb_nodes1i)+np.abs(diff_nb_nodes2i) == 0:
            print()
            print('successfully converged in', i+1, 'iterations')
            break

    return positions_new, displacements


if __name__ == '__main__':
    main()
