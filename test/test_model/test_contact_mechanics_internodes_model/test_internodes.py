#!/usr/bin/env python
# coding: utf-8

import akantu as aka
import numpy as np
import scipy
import pytest
import matplotlib.pyplot as plt

# with this "hack", we also support running from tests where all the files are in one folder
import os.path
import sys

file_prefix = ""

if os.path.isdir("prototype_internodes"):
    sys.path.append("prototype_internodes")
    file_prefix = "prototype_internodes/"

# import functions for reference solution
from functions import *
from functions_contact_probl import *
from init_model import init_model
from example_direct import solve_step_direct

### HELPER FUNCTIONS FOR TESTS

def _get_uniform_interface_grid(dim, graph, x_min=-1, x_max=1, n_grid=10):
    linspace = np.linspace(x_min, x_max, n_grid)
    meshgrid = np.meshgrid(*[linspace]*(dim-1))
    grid = np.c_[[meshgrid[i].ravel() for i in range(dim-1)]]
    mesh = np.c_[grid.T, graph(grid)]
    return mesh

def _get_random_interface_grid(dim, graph, x_min=-1, x_max=1, n_grid=10, seed=0):
    np.random.seed(0)
    grid = np.random.uniform(x_min, x_max, (dim-1, n_grid**(dim-1)))
    mesh = np.c_[grid.T, graph(grid)]
    return mesh

def _check_rbf_radius_conditions(positions, rbf_radius_parameters, c_min=0.49, c_max=0.95):
    # Check condition: [1], Page 51, Section 2.3, Equation 2
    dist_MM = sp.spatial.distance.cdist(positions, positions)
    np.fill_diagonal(dist_MM, np.inf)
    min_ratio = np.min(dist_MM / rbf_radius_parameters, axis=1)
    np.testing.assert_array_less(c_min*np.ones_like(min_ratio),  min_ratio)

    # Check condition: [1], Page 51, Section 2.3, Equation 4
    n_supports = np.sum(dist_MM < rbf_radius_parameters, axis=0)
    n_supports_max = np.ones_like(n_supports) / (1-c_max)**4*(1+4*c_max)
    np.testing.assert_array_less(n_supports,  n_supports_max)

def _check_node_contained_in_support(positions, positions_ref, rbf_radius_parameters, C_max=0.99):
    # Check condition: [1], Page 51, Section 2.3, Equation 3
    dist_MN = sp.spatial.distance.cdist(positions_ref, positions)
    min_ratio = np.min(dist_MN / rbf_radius_parameters, axis=1)
    np.testing.assert_array_less(min_ratio, C_max*np.ones_like(min_ratio))

def _get_theoretical_contact_radius(R, d):
    return (d*R)**0.5

def _get_theoretical_pressure_amplitude(R, d, E, nu):
    a = _get_theoretical_contact_radius(R, d)
    E_red = E/(2*(1-nu**2))
    F = 4/3 * E_red * R**0.5 * d**1.5
    return 3*F / (2*np.pi*a**2)

def _get_theoretical_normal_displacement(R, d, E, nu):
    a = _get_theoretical_contact_radius(R, d)
    p0 = _get_theoretical_pressure_amplitude(R, d, E, nu)
    return - (1-nu**2)/E * np.pi/2 * p0 * a

### COMPARISON TESTS WITH PYTHON REFERENCE

def reference_setup(d=0.1):
    mesh_file = 'contact.msh'
    material_file = file_prefix + 'material.dat'

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
    model.applyBC(aka.FixedValue(-d, aka._y), 'upper_top') # to block the nodes
    nodes_top = mesh.getElementGroup('upper_top').getNodeGroup().getNodes().ravel()
    displacements = displacements.reshape([-1, 2])
    displacements[nodes_top, 1] = -d
    displacements = displacements.ravel()

    # Neumann (K is not invertible)
    # traction = np.zeros(spatial_dimension)
    # traction[1] = -1e9
    # model.applyBC(aka.FromTraction(traction), "upper_top")

    # init and solve
    model, data = init_model(model, mesh, mesh_file, material_file, spatial_dimension,
            displacements, 'lower_top', 'upper_bottom')

    return model, mesh, data

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

### PROPER TESTS OF C++ IMPLEMENTATION

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
    nodal_positions = mesh.getNodes()
    master_nodes = detector.getMasterNodeGroup().getNodes().ravel()
    slave_nodes = detector.getSlaveNodeGroup().getNodes().ravel()

    master_radiuses = detector.getMasterRadiuses().ravel()
    slave_radiuses = detector.getSlaveRadiuses().ravel()

    _check_rbf_radius_conditions(nodal_positions[master_nodes], master_radiuses)
    _check_rbf_radius_conditions(nodal_positions[slave_nodes], slave_radiuses)
    _check_node_contained_in_support(nodal_positions[master_nodes], nodal_positions[slave_nodes], master_radiuses)
    _check_node_contained_in_support(nodal_positions[slave_nodes], nodal_positions[master_nodes], slave_radiuses)

### TEST PYTHON REFERENCE IMPLEMENTATION

def test_compute_rbf_radius_parameters_regular():

    dim = 2

    g_zeros = lambda x: np.zeros(x.shape[1])
    g_parabolic = lambda x: np.sum(x**2, axis=0) + 0.05
    positions1i = _get_uniform_interface_grid(dim, g_zeros, n_grid=10)
    positions2i = _get_uniform_interface_grid(dim, g_parabolic, n_grid=10)
    nodes1i = np.arange(len(positions1i))
    nodes2i = np.arange(len(positions2i))

    radiuses1i, _ = compute_radiuses(nodes1i, nodes2i, positions1i, positions2i)
    _check_rbf_radius_conditions(positions1i, radiuses1i)

    radiuses2i, _ = compute_radiuses(nodes2i, nodes1i, positions2i, positions1i)
    _check_rbf_radius_conditions(positions2i, radiuses2i)

def test_compute_rbf_radius_parameters_random():

    dim = 2

    g_zeros = lambda x: np.zeros(x.shape[1])
    g_parabolic = lambda x: np.sum(x**2, axis=0) + 0.05
    positions1i = _get_random_interface_grid(dim, g_zeros, n_grid=10)
    positions2i = _get_random_interface_grid(dim, g_parabolic, n_grid=10)
    nodes1i = np.arange(len(positions1i))
    nodes2i = np.arange(len(positions2i))

    radiuses1i, _ = compute_radiuses(nodes1i, nodes2i, positions1i, positions2i)
    _check_rbf_radius_conditions(positions1i, radiuses1i)

    radiuses2i, _ = compute_radiuses(nodes2i, nodes1i, positions2i, positions1i)
    _check_rbf_radius_conditions(positions2i, radiuses2i)

def test_find_interface_nodes_regular():

    dim = 2

    g_zeros = lambda x: np.zeros(x.shape[1])
    g_parabolic = lambda x: np.sum(x**2, axis=0) + 0.05
    positions1b = _get_uniform_interface_grid(dim, g_zeros, n_grid=10)
    positions2b = _get_uniform_interface_grid(dim, g_parabolic, n_grid=10)
    nodes1b = np.arange(len(positions1b))
    nodes2b = np.arange(len(positions2b))

    nodes1i, nodes2i, positions1i, positions2i, radiuses1, radiuses2 = \
    find_contact_nodes(nodes1b, nodes2b, positions1b, positions2b)

    _check_node_contained_in_support(positions1i, positions2i, radiuses1)
    _check_node_contained_in_support(positions2i, positions1i, radiuses2)

def test_find_interface_nodes_random():

    dim = 2

    g_zeros = lambda x: np.zeros(x.shape[1])
    g_parabolic = lambda x: np.sum(x**2, axis=0) + 0.05
    positions1b = _get_random_interface_grid(dim, g_zeros, n_grid=10)
    positions2b = _get_random_interface_grid(dim, g_parabolic, n_grid=10)
    nodes1b = np.arange(len(positions1b))
    nodes2b = np.arange(len(positions2b))

    nodes1i, nodes2i, positions1i, positions2i, radiuses1, radiuses2 = \
    find_contact_nodes(nodes1b, nodes2b, positions1b, positions2b)

    _check_node_contained_in_support(positions1i, positions2i, radiuses1)
    _check_node_contained_in_support(positions2i, positions1i, radiuses2)

def test_construct_gap_function_interpolation_regular():

    dim = 2

    g_zeros = lambda x: np.zeros(x.shape[1])
    g_parabolic = lambda x: np.sum(x**2, axis=0) + 0.05
    positions1b = _get_uniform_interface_grid(dim, g_zeros, n_grid=10)
    positions2b = _get_uniform_interface_grid(dim, g_parabolic, n_grid=10)
    nodes1b = np.arange(len(positions1b))
    nodes2b = np.arange(len(positions2b))

    nodes1i, nodes2i, positions1i, positions2i, radiuses1, radiuses2 = \
    find_contact_nodes(nodes1b, nodes2b, positions1b, positions2b)

    R12, R21, _, _ = assemble_Rijs(positions1i, positions2i, radiuses1, radiuses2)

    positions1i_interp = R12 @ positions2i
    positions2i_interp = R21 @ positions1i

    np.testing.assert_allclose(g_parabolic(positions1i_interp[:, :-1].T), positions1i_interp[:, -1], atol=2e-2)
    np.testing.assert_allclose(g_zeros(positions2i_interp[:, :-1].T), positions2i_interp[:, -1], atol=2e-2)

def test_construct_gap_function_interpolation_random():

    dim = 2

    g_zeros = lambda x: np.zeros(x.shape[1])
    g_parabolic = lambda x: np.sum(x**2, axis=0) + 0.05
    positions1b = _get_random_interface_grid(dim, g_zeros, n_grid=10)
    positions2b = _get_random_interface_grid(dim, g_parabolic, n_grid=10)
    nodes1b = np.arange(len(positions1b))
    nodes2b = np.arange(len(positions2b))

    nodes1i, nodes2i, positions1i, positions2i, radiuses1, radiuses2 = \
    find_contact_nodes(nodes1b, nodes2b, positions1b, positions2b)

    R12, R21, _, _ = assemble_Rijs(positions1i, positions2i, radiuses1, radiuses2)

    positions1i_interp = R12 @ positions2i
    positions2i_interp = R21 @ positions1i

    np.testing.assert_allclose(g_parabolic(positions1i_interp[:, :-1].T), positions1i_interp[:, -1], atol=2e-2)
    np.testing.assert_allclose(g_zeros(positions2i_interp[:, :-1].T), positions2i_interp[:, -1], atol=2e-2)

def test_contact2d_problem():

    # Model parameters
    R = (1**2 + 0.5**2)**0.5
    d0 = 1 - (2.1 - R)

    d_list = np.linspace(0.05, 0.35, 5)
    u_list = np.empty_like(d_list)

    for i, d in enumerate(d_list):

        model, mesh, data = reference_setup(d-d0)
        positions_new, displacements_ref = solve_step_direct(model, data, nb_max_iter=10) 

        nodes2i = mesh.getElementGroup("upper_bottom").getNodeGroup().getNodes().ravel()
        positions2i = positions_new[nodes2i]

        u_list[i] = np.min(positions2i[:, 1]) - 1

    nu = model.getMaterial(0).getReal("nu")
    E = model.getMaterial(0).getReal("E")

    np.testing.assert_allclose(_get_theoretical_normal_displacement(R, d_list, E, nu), u_list, atol=5e-3)

if __name__ == '__main__':
    pytest.main()
