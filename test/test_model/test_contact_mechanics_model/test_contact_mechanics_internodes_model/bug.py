#!/usr/bin/env python
# coding: utf-8

import akantu as aka
import numpy as np
import scipy
import matplotlib.pyplot as plt

from prototype_internodes.functions_contact_probl import * 
from prototype_internodes.functions_direct import * 
from prototype_internodes.example_direct import \
        init_direct, solve_step_direct


def reference_setup():
    mesh_file = 'bug.msh'
    material_file = 'prototype_internodes/material.dat'

    aka.parseInput(material_file)
    spatial_dimension = 2

    mesh = aka.Mesh(spatial_dimension)
    mesh.read(mesh_file)

    # initialize model
    model = aka.SolidMechanicsModel(mesh)
    model.initFull(_analysis_method=aka._implicit_dynamic)

    # boundary conditions
    model.applyBC(aka.FixedValue(0., aka._x), 'lower_bottom')
    model.applyBC(aka.FixedValue(0., aka._y), 'lower_bottom')

    # Dirichlet
    # model.applyBC(aka.FixedValue(-0.1, aka._y), 'upper_top')

    # Neumann (K is not invertible)
    traction = np.zeros(spatial_dimension)
    traction[1] = -1e9
    model.applyBC(aka.FromTraction(traction), "upper_top")
    model.applyBC(aka.FixedValue(0., aka._x), 'upper_top')

    traction = np.zeros(spatial_dimension)
    traction[1] = -1e9
    model.applyBC(aka.FromTraction(traction), "upper_top")

    # init
    model, data = init_direct(model, mesh, mesh_file, material_file, spatial_dimension,
            "upper_bottom", "lower_top")

    return model, mesh, data

def test_assembleInternodesMatrix():
    # reference
    ref_model, mesh, ref_data = reference_setup()



    nodes1i, nodes2i, positions1i, positions2i, radiuses1, radiuses2 = \
        find_contact_nodes(ref_data['nodes1b'], ref_data['nodes2b'],
            ref_data['positions1b'], ref_data['positions2b'])

    dofs1i = nodes_to_dofs(nodes1i).ravel()
    dofs2i = nodes_to_dofs(nodes2i).ravel()

    spatial_dimension = 2
    nb_constraint_dofs = len(dofs1i)
    nb_free_dofs = spatial_dimension*len(ref_data['nodes'])

    _, _, R12, R21 = assemble_Rijs(positions1i, positions2i, radiuses1, radiuses2)
    M1i, M2i = assemble_interface_masses(ref_model, dofs1i, dofs2i)
    B, B_tilde, C = assemble_Bs(M1i, M2i, R12, R21, dofs1i, dofs2i, nb_free_dofs, nb_constraint_dofs) 


    ref_model.assembleStiffnessMatrix()
    K_aka = ref_model.dof_manager.getMatrix("K")
    K = sp.sparse.csc_matrix(aka.AkantuSparseMatrix(K_aka))

    A_ref = assemble_A_explicit(K, B, B_tilde, C).todense()

    # akantu
    contact_model = aka.ContactMechanicsInternodesModel(mesh)
    detector = contact_model.getContactDetectorInternodes();

    contact_model.initFull(_analysis_method=aka._static)
    contact_model.assembleInternodesMatrix()

    A_aka = contact_model.dof_manager.getMatrix("K")
    A = sp.sparse.csc_matrix(aka.AkantuSparseMatrix(A_aka)).todense()

    print(A)
    print(A_ref)

    B_master = A[:nb_free_dofs, nb_free_dofs:]
    B_master_ref = A_ref[:nb_free_dofs, nb_free_dofs:]

    B_slave = A[:nb_free_dofs, nb_free_dofs:]
    B_slave_ref = A_ref[:nb_free_dofs, nb_free_dofs:]

    B_tilde_master = A[nb_free_dofs:, :nb_free_dofs]
    B_tilde_master_ref = A_ref[nb_free_dofs:, :nb_free_dofs]

    B_tilde_slave = A[nb_free_dofs:, :nb_free_dofs]
    B_tilde_slave_ref = A_ref[nb_free_dofs:, :nb_free_dofs]

    np.testing.assert_allclose(B_master, B_master_ref)
    np.testing.assert_allclose(B_slave, B_slave_ref)
    np.testing.assert_allclose(B_tilde_master, B_tilde_master_ref)
    np.testing.assert_allclose(B_tilde_slave, B_tilde_slave_ref)

    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if A[i, j] != A_ref[i, j]:
                rel_error = np.abs((A[i, j]-A_ref[i, j]) / A_ref[i, j])
                if (rel_error > 1e-7):
                    print("(%d, %d): %0.20f !=and %0.20f" % (i, j, A[i, j], A_ref[i, j]))

    np.testing.assert_allclose(A[:nb_free_dofs, :nb_free_dofs],
            A_ref[:nb_free_dofs, :nb_free_dofs])

if __name__ == '__main__':
    test_assembleInternodesMatrix()
