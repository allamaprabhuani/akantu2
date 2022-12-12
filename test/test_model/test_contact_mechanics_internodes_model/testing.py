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

def reference_setup():
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

def _get_uniform_interface_grid(dim, graph, x_min=-1, x_max=1, n_grid=10):
    linspace = np.linspace(x_min, x_max, n_grid)
    meshgrid = np.meshgrid(*[linspace]*(dim-1))
    grid = np.c_[[meshgrid[i].ravel() for i in range(dim-1)]]
    mesh = np.c_[grid.T, graph(grid)]
    return mesh

def _get_random_interface_grid(dim, graph, x_min=-1, x_max=1, n_grid=10, seed=1):
    np.random.seed(0)
    grid = np.random.uniform(x_min, x_max, (dim-1, n_grid**(dim-1)))
    mesh = np.c_[grid.T, graph(grid)]
    return mesh

def _check_node_contained_in_support(positions, positions_ref, rbf_radius_parameters, C_max=0.99):
    # Check condition: [1], Page 51, Section 2.3, Equation 3
    dist_MN = scipy.spatial.distance.cdist(positions_ref, positions)
    min_ratio = np.min(dist_MN / rbf_radius_parameters, axis=1)
    np.testing.assert_array_less(min_ratio, C_max*np.ones_like(min_ratio))

dim = 2

g_zeros = lambda x: np.zeros(x.shape[1])
g_parabolic = lambda x: np.sum(x**2, axis=0) + 0.05
positions_primary = _get_random_interface_grid(dim, g_zeros, n_grid=10)
positions_secondary = _get_random_interface_grid(dim, g_parabolic, n_grid=10)
nodes_primary = np.arange(len(positions_primary))
nodes_secondary = np.arange(len(positions_secondary))

ref_model, mesh, ref_data = reference_setup()

print(ref_data['nodes1b'].shape, nodes_primary.shape)
print(ref_data['positions1b'].shape, positions_primary.shape)

nodes1i, nodes2i, positions1i, positions2i, radiuses1, radiuses2 = \
    find_contact_nodes(ref_data['nodes1b'], ref_data['nodes2b'],
        ref_data['positions1b'], ref_data['positions2b'])

print(nodes1i.shape)
print(positions1i.shape)
print(radiuses1.shape)

nodes1i, nodes2i, positions1i, positions2i, radiuses1, radiuses2 = \
find_contact_nodes(nodes_primary, nodes_secondary,
    positions_primary, positions_secondary)

print(positions1i, positions2i,  radiuses1, radiuses2)