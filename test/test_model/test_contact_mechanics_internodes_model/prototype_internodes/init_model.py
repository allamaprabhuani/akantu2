#!/usr/bin/env python
# coding: utf-8

import akantu as aka
import numpy as np

import sys
sys.path.append("..")
from prototype_internodes.functions import nodes_to_dofs


def init_model(model, mesh, mesh_file, material_file, spatial_dimension,
        displacements, gmsh_interface1, gmsh_interface2):
    # dictonary to store mesh data
    data = {}

    #---------- initialize mesh data ----------
    # displacements
    data['displacements'] = displacements

    # data['positions] of nodes
    data['nb_nodes'] = mesh.getNbNodes()
    data['nodes'] = np.arange(data['nb_nodes'])
    data['positions'] = mesh.getNodes()
    data['dofs'] = nodes_to_dofs(data['nodes'])
    data['nb_dofs'] = data['nb_nodes'] * spatial_dimension

    # connectivity
    data['connectivity'] = mesh.getConnectivity(aka._triangle_3)
    data['connectivity_boundary'] = mesh.getConnectivity(aka._segment_2)

    # coordinates, nodes and dofs of body 1
    data['nodes1'] = mesh.getElementGroup("body_lower").getNodeGroup().getNodes().ravel()
    data['positions1'] = data['positions'][data['nodes1']]
    data['dofs1'] = nodes_to_dofs(data['nodes1']).ravel()

    # coordinates, nodes and dofs of body 2 
    data['nodes2'] = mesh.getElementGroup("body_upper").getNodeGroup().getNodes().ravel()
    data['positions2'] = data['positions'][data['nodes2']]
    data['dofs2'] = nodes_to_dofs(data['nodes2']).ravel()


    #---------- initialize boundary ----------
    # 1b denotes boundary 1, 2b denotes boundary 2
    # boundary in this case means the potential interface nodes!

    # get nodes from selected surfaces
    data['nodes1b'] = mesh.getElementGroup(gmsh_interface1).getNodeGroup().getNodes().ravel()
    data['nodes2b'] = mesh.getElementGroup(gmsh_interface2).getNodeGroup().getNodes().ravel()

    # or get nodes via akantu boundary algorithm
    # _ = mesh.createBoundaryGroupFromGeometry()
    # data['nodes1b'] = mesh.getElementGroup("boundary_0").getNodeGroup().getNodes().ravel()
    # data['nodes2b'] = mesh.getElementGroup("boundary_1").getNodeGroup().getNodes().ravel()

    data['positions1b'] = data['positions'][data['nodes1b']]
    data['dofs1b'] = nodes_to_dofs(data['nodes1b']).ravel()

    data['positions2b'] = data['positions'][data['nodes2b']]
    data['dofs2b'] = nodes_to_dofs(data['nodes2b']).ravel()

    data['boundary_dofs'] = np.append(data['dofs1b'], data['dofs2b'])
    data['nb_boundary_dofs'] = len(data['boundary_dofs']) 

    # connectivity of segements on the boundary
    data['connectivity1b'] = data['connectivity_boundary'][np.in1d(data['connectivity_boundary'],
            data['nodes1b']).reshape(data['connectivity_boundary'].shape).any(axis=1)]
    data['connectivity2b'] = data['connectivity_boundary'][np.in1d(data['connectivity_boundary'],
            data['nodes2b']).reshape(data['connectivity_boundary'].shape).any(axis=1)]

    # connectivity of elements belonging to boundary
    data['connectivity1b_body'] = data['connectivity'][np.in1d(data['connectivity'],
            data['nodes1b']).reshape(data['connectivity'].shape).any(axis=1)]
    data['connectivity2b_body'] = data['connectivity'][np.in1d(data['connectivity'],
            data['nodes2b']).reshape(data['connectivity'].shape).any(axis=1)]


    #---------- boundary conditions ----------
    # get free dofs 
    blocked_dofs_mask = model.getBlockedDOFs().ravel()
    data['blocked_dofs'] = data['dofs'][blocked_dofs_mask].ravel()
    data['free_dofs'] = data['dofs'][~blocked_dofs_mask].ravel()
    data['nb_free_dofs'] = len(data['free_dofs'])

    # remove nodes wich are blocked (if akantu boundary algorithm)
    # data['dofs1b'] = data['dofs1b'][np.in1d(data['dofs1b'], data['free_dofs'])]
    # data['dofs2b'] = data['dofs2b'][np.in1d(data['dofs2b'], data['free_dofs'])]

    return model, data
