#!/usr/bin/env python
# coding: utf-8

import numpy as np
import scipy as sp
from scipy.linalg import lstsq 
from scipy import spatial

from functions import nodes_to_dofs


def find_contact_nodes(nodes1i, nodes2i, positions1i, positions2i):
    """Contact nodes algorithm.

    :param nodes1i: nodes of body 1 interface (master)
    :param nodes2i: nodes of body 2 interface (slave)
    :param positions1i: initial positions of body 1 interface (master)
    :param positions2i: initial positions of body 2 interface (slave)
    :returns: selected nodes and radiuses
    """
    radiuses1, nnzR21 = compute_radiuses(nodes1i, nodes2i, positions1i, positions2i)
    radiuses2, nnzR12 = compute_radiuses(nodes2i, nodes1i, positions2i, positions1i)

    nodes1i_mask = nnzR12 > 0
    nodes2i_mask = nnzR21 > 0

    while np.any(nodes1i_mask == False) or np.any(nodes2i_mask == False):
        nodes1i = nodes1i[nodes1i_mask]
        nodes2i = nodes2i[nodes2i_mask]

        positions1i = positions1i[nodes1i_mask]
        positions2i = positions2i[nodes2i_mask]

        dofs1i = nodes_to_dofs(nodes1i).ravel()
        dofs2i = nodes_to_dofs(nodes2i).ravel()

        radiuses1, nnzR21 = compute_radiuses(nodes1i, nodes2i, positions1i, positions2i)
        radiuses2, nnzR12 = compute_radiuses(nodes2i, nodes1i, positions2i, positions1i)

        nodes1i_mask = nnzR12 > 0
        nodes2i_mask = nnzR21 > 0

    return nodes1i, nodes2i, positions1i, positions2i, radiuses1, radiuses2

def compute_radiuses(nodes1i, nodes2i, positions1i, positions2i):
    """Computes the radiuses of attack.

    :param nodes1i: nodes of body 1 interface (master)
    :param nodes2i: nodes of body 2 interface (slave)
    :param positions1i: initial positions of body 1 interface (master)
    :param positions2i: initial positions of body 2 interface (slave)
    :returns: radiuses, and resulting nonzero components of R matrix
    """
    # fixed parameters
    c = 0.5 # conditition (2)
    C = 0.95 # condition (3)
    n = 1 # consider 1 nearest neighboors
    d = 0.05 # tolerance, for radius of "attack" estimation

    M = len(positions1i)
    N = len(positions2i)

    radiuses = np.zeros(M)
    nnzRMM = np.zeros(M)
    nnzRNM = np.zeros(N)
    nnzCMM = np.zeros(M)
    nnzCNM = np.zeros(M)

    maxS = np.inf
    f = 0
    niter = 0
    maxiter = 10

    while maxS > f and niter < maxiter-1:
        f = np.floor(1/(np.power(1-c, 4)*(1+4*c))) # maximum number of supports

        for k in range(M):
            point = positions1i[k, :].reshape(1, -1)
            neighbors = positions1i.copy()
            neighbors[k, :] = np.inf
            distMM = spatial.distance.cdist(neighbors, point).ravel()
            distMN = spatial.distance.cdist(positions2i, point).ravel()

            rMM = np.min(distMM)
            rNM = np.sqrt(d*d + 0.25*np.power(rMM, 2))
            radius = np.maximum(rMM, rNM)

            if radius > rMM/c:
                radius = rMM/c

            s1 = distMM < radius
            s2 = distMN < C*radius

            nnzRMM[s1] = nnzRMM[s1] + 1
            nnzRNM[s2] = nnzRNM[s2] + 1
            nnzCMM[k] = np.sum(s1)
            nnzCNM[k] = np.sum(s2)

            radiuses[k] = radius

        maxS = np.max(nnzRMM)

        if maxS > f:
            c = 0.5 * (1 + c);
            nnzRMM = np.zeros(M)
            nnzRNM = np.zeros(N)
            nnzCMM = np.zeros(M)
            nnzCNM = np.zeros(M)
            niter = niter+1

    return radiuses, nnzRNM

def wendland(dists, radiuses):
    """ Compute the Beckert & Wendland RBF

    :param dists: distances 
    :param radiuses: radius for each distance (same length)
    """
    result = np.zeros(len(dists))

    mask = dists <= radiuses
    result[mask] = np.power(1-dists[mask]/radiuses[mask], 4) * (1+4*dists[mask]/radiuses[mask])
    return result

def phi_constructor(positions_i, positions_j, radiuses_j, rad_func):
    """Construct Phi for RBF interpolation.

    :param positions_i: positions of evaluation points
    :param positions_j: positions of reference points 
    :param radiuses_j: radiuses of reference points
    """
    N = len(positions_i)
    M = len(positions_j)

    dists = spatial.distance.cdist(positions_i, positions_j)
    radiuses_j = np.tile(radiuses_j, N)
    phi = rad_func(dists.ravel(), radiuses_j.ravel())

    return phi.reshape([N, M])

def Rij_constructor(positions_i, positions_j, radiuses_j):
    """Construct a RBF interpolation.

    :param positions_i: positions of evaluation points
    :param positions_j: positions of reference points 
    :param radiuses_j: radiuses of reference points
    """
    phiMM = phi_constructor(positions_j, positions_j, radiuses_j, wendland)
    phiNM = phi_constructor(positions_i, positions_j, radiuses_j, wendland)

    Rij = phiNM.dot(np.linalg.inv(phiMM))
    g = Rij.dot(np.ones((Rij.shape[1], 1)))
    Rij_norm = Rij * (1/g)
    return Rij_norm

def assemble_Rijs(positions1i, positions2i, radiuses1, radiuses2):
    """ Assembles the Radial Basis function (RBF) interpolation matrices.

    :param positions_i: positions of evaluation points
    :param positions_j: positions of reference points 
    :param radiuses_i: radiuses of evaluation points
    :param radiuses_j: radiuses of reference points
    """
    R12_normal = Rij_constructor(positions1i, positions2i, radiuses2)
    R21_normal = Rij_constructor(positions2i, positions1i, radiuses1)

    R21 = sp.sparse.csr_matrix(extend_to_2D(R21_normal))
    R12 = sp.sparse.csr_matrix(extend_to_2D(R12_normal))

    return R12_normal, R21_normal, R12, R21

def extend_to_2D(R):
    """ Extend the RBF matrices to 2D. """
    R_extended = np.repeat(np.repeat(R,2,axis=1), 2, axis=0)
    R_extended[1::2,::2] = 0
    R_extended[::2,1::2] = 0
    return R_extended

def remove_traction(positions_new, connectivity1b, connectivity2b, connectivity1b_body, connectivity2b_body,
        nodes1i, nodes2i, nodes1b, nodes2b, lambda1, R12, R21):

    """ Check if there is still traction between the bodies.

    :param positions_new: new positions after internodes step
    :param connectivity1b: connectivity of boundary surface 1 (master) 
    :param connectivity2b: connectivity of boundary surface 2 (slave) 
    :param connectivity1b_body: connectivity of boundary elements 1 (master) 
    :param connectivity2b_body: connectivity of boundary elements 2 (slave) 
    :param nodes1i: nodes of interface 1
    :param nodes2i: nodes of interface 2
    :param nodes1b: nodes of boundary 1 (potential interface nodes)
    :param nodes2b: nodes of boundary 2 (potential interface nodes)
    :param lambda1: lambdas of interface 1 (solution)
    :param R12: RBF interpolation master to slave
    :param R21: RBF interpolation slave to master
    """
    normals1b = compute_normals(positions_new, nodes1b, connectivity1b, connectivity1b_body)
    normals2b = compute_normals(positions_new, nodes2b, connectivity2b, connectivity2b_body)

    normals1i = normals1b[np.in1d(nodes1b, nodes1i)]
    normals2i = normals2b[np.in1d(nodes2b, nodes2i)]

    lambda2 = -R21.dot(lambda1.reshape([-1, 1])).reshape([-1, 2])

    scalar1 = np.sum(lambda1*normals1i, axis=1)
    scalar2 = np.sum(lambda2*normals2i, axis=1)

    nodes1i_dump = nodes1i[scalar1>0]
    nodes2i_dump = nodes2i[scalar2>0]

    if len(nodes1i_dump) == 0 and len(nodes2i_dump) == 0:
        # gap verification
        nodes1i_add, nodes2i_add = detect_gaps(positions_new, nodes1i, nodes2i, normals1i, normals2i)

        nodes1i, diff_nb_nodes1i = update_interface(nodes1i_add, nodes1i, 'add')
        nodes2i, diff_nb_nodes2i = update_interface(nodes2i_add, nodes2i, 'add')

        print(diff_nb_nodes1i, ' nodes added to interface 1')
        print(diff_nb_nodes2i, ' nodes added to interface 2')
    else:
        nodes1i, diff_nb_nodes1i = update_interface(nodes1i_dump, nodes1i, 'dump')
        nodes2i, diff_nb_nodes2i = update_interface(nodes2i_dump, nodes2i, 'dump')

        print(diff_nb_nodes1i, ' nodes removed from interface 1')
        print(diff_nb_nodes2i, ' nodes removed from interface 2')

    return nodes1i, nodes2i, diff_nb_nodes1i, diff_nb_nodes2i

def update_interface(new_nodes, nodesi, case):
    """ Remove of add interfaces after traction check. """
    if case == 'dump':
        nodesi_new = nodesi[~np.in1d(nodesi, new_nodes)]
    if case == 'add':
        nodesi_new = np.union1d(nodesi, new_nodes)

    diff_nb_nodes = len(nodesi) -len(nodesi_new)
    return nodesi_new, diff_nb_nodes

def compute_normals(positions_new, nodesb, connectivityb, connectivityb_body):
    """Comput normals on interface surface.

    :param positions_new: new positions after internodes step
    :param nodesb: nodes of boundary(potential interface nodes)
    :param connectivityb: connectivity of boundary surface
    """
    n = len(nodesb)
    m = len(connectivityb)

    connectivityi_body = connectivityb_body[np.in1d(connectivityb_body, nodesb).reshape(connectivityb_body.shape).any(axis=1)]
    nodesb_body = np.unique(connectivityi_body[~np.isin(connectivityb_body, nodesb)])

    tangents = positions_new[connectivityb[:, 1]] - positions_new[connectivityb[:, 0]]
    lengths = np.linalg.norm(tangents, axis=1).reshape([-1,1])
    tangents = tangents/lengths

    normals = np.zeros((m, 2))
    normals[:, 0] = -tangents[:, 1]
    normals[:, 1] = tangents[:, 0]

    normals_avg = np.zeros((n, 2))
    gamma = 1e-3 # step size
    for j in range(n):
        node = nodesb[j]
        coord = positions_new[node, :]
        id = np.in1d(connectivityb, node).reshape(connectivityb.shape).any(axis=1)
        length = lengths[id]
        normal_avg = 1/np.sum(length)*np.sum(normals[id, :]*length, axis=0)

        tang_plus = (coord + gamma*normal_avg).reshape([-1, 2])
        tang_minus = (coord - gamma*normal_avg).reshape([-1, 2])

        min_plus = np.min(spatial.distance.cdist(positions_new[nodesb_body, :], tang_plus).ravel())
        min_minus = np.min(spatial.distance.cdist(positions_new[nodesb_body, :], tang_minus).ravel())

        if min_plus > min_minus:
            normals_avg[j, :] = normal_avg
        else:
            normals_avg[j, :] = -normal_avg

    norms = np.linalg.norm(normals_avg, axis=1).reshape([-1, 1])
    normals_avg = (normals_avg/norms).reshape([-1, 2])
    return normals_avg

def detect_gaps(positions_new, nodes1i, nodes2i, normals1i, normals2i):
    """Detect gaps between interfaces.

    :param positions_new: new positions after internodes step
    :param nodes1i: nodes of interface 1
    :param nodes2i: nodes of interface 2
    :param normals1i: normal vectors of nodes on interface 1
    :param normals2i: normal vectors of nodes on interface 2
    """
    tol = 0.9 # tolerance for gap detection (could be changed as input)
    h = 0.05 # mesh size (shouldn't be fixed!)
    positions1i = positions_new[nodes1i, :]
    positions2i = positions_new[nodes2i, :]

    nodes1i, nodes2i, positions1i, positions2i, radiuses1, radiuses2 = find_contact_nodes(nodes1i, nodes2i, positions1i, positions2i)

    R21_normal = Rij_constructor(positions2i, positions1i, radiuses1)
    R21 = sp.sparse.csr_matrix(extend_to_2D(R21_normal))

    R12_normal = Rij_constructor(positions1i, positions2i, radiuses2)
    R12 = sp.sparse.csr_matrix(extend_to_2D(R12_normal))

    diffs1 = R12.dot(positions2i.reshape([-1, 1])) - positions1i.reshape([-1, 1])
    diffs1 = diffs1.reshape([-1, 2])
    diffs2 = R21.dot(positions1i.reshape([-1, 1])) - positions2i.reshape([-1, 1])
    diffs2 = diffs2.reshape([-1, 2])

    scalar1 = np.sum(diffs1*normals1i, axis=1)
    scalar2 = np.sum(diffs2*normals2i, axis=1)

    threshold = -tol*h

    nodes1i_add = nodes1i[scalar1<threshold]
    nodes2i_add = nodes2i[scalar2<threshold]

    return nodes1i_add, nodes2i_add
