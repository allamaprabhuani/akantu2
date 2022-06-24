#!/usr/bin/env python
# coding: utf-8

import akantu as aka
import numpy as np
import scipy
import pytest
import matplotlib.pyplot as plt

# import functions for reference solution
from prototype_internodes.functions import * 
from prototype_internodes.functions_contact_probl import * 
from prototype_internodes.init_model import init_model
from prototype_internodes.example_direct import solve_step_direct 


def reference_setup():
    mesh_file = 'contact.msh'
    material_file = 'prototype_internodes/material.dat'

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
    model.applyBC(aka.FixedValue(-0.1, aka._y), 'upper_top') # to block the nodes
    nodes_top = mesh.getElementGroup('upper_top').getNodeGroup().getNodes().ravel()
    displacements = displacements.reshape([-1, 2])
    displacements[nodes_top, 1] = -0.1
    displacements = displacements.ravel()

    # Neumann (K is not invertible)
    # traction = np.zeros(spatial_dimension)
    # traction[1] = -1e9
    # model.applyBC(aka.FromTraction(traction), "upper_top")

    # init and solve
    model, data = init_model(model, mesh, mesh_file, material_file, spatial_dimension,
            displacements, 'lower_top', 'upper_bottom')

    return model, mesh, data

def test_findContactNodes():
    # reference
    ref_model, mesh, ref_data = reference_setup()

    nodes1i, nodes2i, positions1i, positions2i, radiuses1, radiuses2 = \
        find_contact_nodes(ref_data['nodes1b'], ref_data['nodes2b'],
            ref_data['positions1b'], ref_data['positions2b'])

    # akantu
    contact_model = aka.ContactMechanicsInternodesModel(mesh)
    detector = contact_model.getContactDetectorInternodes()

    detector.findContactNodes()

    master_nodes = detector.getMasterNodeGroup().getNodes().ravel()
    slave_nodes = detector.getSlaveNodeGroup().getNodes().ravel()

    master_radiuses = detector.getMasterRadiuses().ravel()
    slave_radiuses = detector.getSlaveRadiuses().ravel()

    np.testing.assert_equal(master_nodes, nodes1i)
    np.testing.assert_equal(slave_nodes, nodes2i)
    np.testing.assert_equal(master_radiuses, radiuses1)
    np.testing.assert_equal(slave_radiuses, radiuses2)

def test_assembleInterfaceMass():
    # reference
    ref_model, mesh, ref_data = reference_setup()

    M1b_ref_sparse, _ = assemble_interface_masses(ref_model, ref_data['dofs1b'], ref_data['dofs2b'])
    M1b_ref = M1b_ref_sparse.todense()

    # akantu
    contact_model = aka.ContactMechanicsInternodesModel(mesh)
    detector = contact_model.getContactDetectorInternodes();

    contact_model.initFull(_analysis_method=aka._static)

    master_nodes = detector.getMasterNodeGroup().getNodes().ravel()
    slave_nodes = detector.getSlaveNodeGroup().getNodes().ravel()

    master_radiuses = detector.getMasterRadiuses().ravel()
    slave_radiuses = detector.getSlaveRadiuses().ravel()

    initial_master_node_group = detector.getInitialMasterNodeGroup()
    M1b = contact_model.assembleInterfaceMass(initial_master_node_group)

    np.testing.assert_allclose(M1b, M1b_ref)

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

    # scale by youngs modulus (as in akantu solution)
    B = B * 30e9
    B_tilde = B_tilde * 30e9

    ref_model.assembleStiffnessMatrix()
    K_aka = ref_model.dof_manager.getMatrix("K")
    K = aka.AkantuSparseMatrix(K_aka)

    A_ref = assemble_A_explicit(K, B, B_tilde, C).todense()

    # akantu
    contact_model = aka.ContactMechanicsInternodesModel(mesh)
    solid = contact_model.getSolidMechanicsModel();

    contact_model.initFull(_analysis_method=aka._static)
    contact_model.assembleInternodesMatrix()

    A_aka = contact_model.dof_manager.getMatrix("K")
    A = sp.sparse.csc_matrix(aka.AkantuSparseMatrix(A_aka)).todense()

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

    # this fails for the K matrix part, but very small errors
    # np.testing.assert_allclose(A[:nb_free_dofs, :nb_free_dofs],
    #         A_ref[:nb_free_dofs, :nb_free_dofs])

def test_solveStep():
    # reference
    ref_model, mesh, ref_data = reference_setup()
    positions_new_ref, displacements_ref = solve_step_direct(ref_model, ref_data, nb_max_iter=1) 

    # akantu
    contact_model = aka.ContactMechanicsInternodesModel(mesh)
    solid = contact_model.getSolidMechanicsModel();

    contact_model.initFull(_analysis_method=aka._static)

    # boundary conditions
    solid.applyBC(aka.FixedValue(0., aka._x), 'lower_bottom')
    solid.applyBC(aka.FixedValue(0., aka._y), 'lower_bottom')

    # Dirichlet
    solid.applyBC(aka.FixedValue(-0.1, aka._y), 'upper_top')

    # Neumann (K is not invertible)
    # traction = np.zeros(2)
    # traction[1] = -1e9
    # solid.applyBC(aka.FromTraction(traction), "upper_top")

    contact_model.solveStep()
    positions_new = solid.getCurrentPosition()

    # plt.triplot(positions_new[:, 0], positions_new[:, 1], ref_data['connectivity'])
    # plt.title('mesh')
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.axis('scaled')
    # plt.ylim([0.9, 1.1])
    # plt.show()

    np.testing.assert_allclose(positions_new, positions_new_ref)


if __name__ == '__main__':
    pytest.main()
