"""
References
----------

[1] Y. Voet et. al.: The INTERNODES method for applications in contact mechanics
    and dedicated preconditioning techniques.
    Computers & Mathematics with Applications, vol. 127, 2022, pp. 48-64
[2] Y. Voet: On the preconditioning of the INTERNODES matrix for applications
    in contact mechanics.
    Master's thesis EPFL, 2021.
"""

import numpy as np
import scipy as sp
import akantu as aka

def nodes_to_dofs(nodes, dim):
    """Obtain the DOFs that correspond to the indices in nodes"""
    return dim*np.repeat(nodes, dim) + np.tile(range(dim), nodes.shape[0])

def expand_to_dim(matrix, dim):
    """Expand matrix entries to (dim x dim)-blocks corresponding to the DOFs"""
    return sp.sparse.kron(matrix, np.eye(dim))

def remove_rows_without_elements(matrix, elements):
    """Remove rows from matrix that are not contained in elements"""
    return matrix[np.isin(matrix, elements).any(axis=1)]

def wendland(delta):
    """Evaluate Wendland C2 function (see [1], Table 1)"""
    return (1 - delta)**4 * (1 + 4*delta) * (delta <= 1)

def wendland_rbf(distances, radii):
    """Evaluate Wendland radial basis function (see [1], p. 49)"""
    return wendland(distances / radii)

def construct_rbf_matrix(positions, positions_ref, radii_ref, rbf):
    """Construct radial basis matrix $\Phi_{NM}$ (see [1], p. 49)"""
    distance_matrix = sp.spatial.distance.cdist(positions, positions_ref)
    return rbf(distance_matrix, radii_ref)

def construct_interpolation_matrix(positions, positions_ref, radii_ref, rbf):
    """Construct interpolation matrix $\R_{NM}$ (see [1], p. 49-50)"""
    rbf_matrix_self = construct_rbf_matrix(positions_ref, positions_ref, radii_ref, rbf)
    rbf_matrix = construct_rbf_matrix(positions, positions_ref, radii_ref, rbf)
    
    interpolation_matrix = np.linalg.solve(rbf_matrix_self.T, rbf_matrix.T).T
    normalization = np.sum(interpolation_matrix, axis=1)[:, np.newaxis]
    return sp.sparse.csr_matrix(interpolation_matrix / normalization)

def compute_radii(positions, positions_ref, lower_precision_bound=0.5, upper_precision_bound=0.95, attack_radius_tolerance=0.05, max_iter=10):
    """Get radius parameters for radial basis functions (see [1], p. 49, pt. 1)

    lower_precision_bound 
        [1], eq. (2)
    upper_precision_bound 
        [1], eq. (3)
    attack_radius_tolerance
        TODO: makes no sense...
    """
    radii = np.zeros(positions.shape[0])

    # Criterion for strict diagonal dominance by rows ([1], eq. (4))
    max_n_supports = np.floor(1 / wendland(lower_precision_bound))

    for _ in range(max_iter):
        nnzRMM = np.zeros(positions.shape[0])
        nnzRNM = np.zeros(positions_ref.shape[0])

        for k in range(positions.shape[0]):
            distMM = np.linalg.norm(positions - positions[k, None, :], axis=-1)
            distMM[k] = np.inf
            distMN = np.linalg.norm(positions_ref - positions[k, None, :], axis=-1)

            rMM = np.min(distMM)
            # TODO: Magic number 0.25
            rNM = (attack_radius_tolerance**2 + 0.25*rMM**2)**0.5
            radii[k] = min(max(rMM, rNM), rMM / lower_precision_bound)

            nnzRMM[distMM < radii[k]] += 1
            nnzRNM[distMN < upper_precision_bound*radii[k]] += 1

        if np.max(nnzRMM) <= max_n_supports:
            break

        lower_precision_bound = (1 + lower_precision_bound) / 2

    return radii, nnzRNM

class ContactMechanicsInternodes(object):

    def __init__(self, dim, mesh, model, gmsh_interface1, gmsh_interface2, rbf=wendland_rbf):
        self.dim = dim
        self.mesh = mesh
        self.model = model
        self.nodes_interface_master = mesh.getElementGroup(gmsh_interface1).getNodeGroup().getNodes().ravel()
        self.nodes_interface_slave = mesh.getElementGroup(gmsh_interface2).getNodeGroup().getNodes().ravel()
        self.nodes_boundary_master = mesh.getElementGroup(gmsh_interface1).getNodeGroup().getNodes().ravel()
        self.nodes_boundary_slave = mesh.getElementGroup(gmsh_interface2).getNodeGroup().getNodes().ravel()
        
        self.positions_master = mesh.getNodes()[self.nodes_interface_master]
        self.positions_slave = mesh.getNodes()[self.nodes_interface_slave]

        self.dofs_master = nodes_to_dofs(self.nodes_interface_master, dim=self.dim)
        self.dofs_slave = nodes_to_dofs(self.nodes_interface_slave, dim=self.dim)
        self.dofs = np.arange(mesh.getNbNodes()*self.dim)

        blocked_dofs_mask = model.getBlockedDOFs().ravel()
        self.dofs_blocked = self.dofs[blocked_dofs_mask]
        self.dofs_free = self.dofs[~blocked_dofs_mask]

        self.n_dofs_free = len(self.dofs_free)
        self.n_constraint_dofs = len(self.dofs_master)

        connectivity_boundary = mesh.getConnectivity(aka._segment_2)
        self.connectivity_boundary_master = remove_rows_without_elements(connectivity_boundary, self.nodes_interface_master)
        self.connectivity_boundary_slave = remove_rows_without_elements(connectivity_boundary, self.nodes_interface_slave)

        connectivity_body = mesh.getConnectivity(aka._triangle_3)
        self.connectivity_boundary_body_master = remove_rows_without_elements(connectivity_body, self.nodes_interface_master)
        self.connectivity_boundary_body_slave = remove_rows_without_elements(connectivity_body, self.nodes_interface_slave)

        self.rbf = rbf
        self.radii_master = None
        self.radii_slave = None

        self.normals_interface_master = None
        self.normals_interface_slave = None
        
        self.rescaling_factor = 30e+9  # TODO: Tunable, Young's modulus
        self.K = None  # Stiffness matrix

        self.R12 = None  # Interpolation matrix from slave to master
        self.R21 = None  # Interpolation matrix from master to slave

        self.M1 = None  # Interface mass matrix of mater
        self.M2 = None  # Interface mass matrix of slave

        self.B = None  # Block component of INTERNODES matrix
        self.B_tilde = None  # Block component of INTERNODES matrix

        self.internodes_matrix = None  # INTERNODES matrix (see [1], eq. 13)
        self.force_term = None  # Force term b (see [1], eq. 13)

    def find_contact_nodes(self):
        """Identify contact nodes (see [2], p. 22, algorithm 1)"""
        while True:
            self.radii_master, nnzR21 = compute_radii(self.positions_master, self.positions_slave)
            self.radii_slave, nnzR12 = compute_radii(self.positions_slave, self.positions_master)

            nodes_mask_master = nnzR12 > 0
            nodes_mask_slave = nnzR21 > 0

            self.nodes_interface_master = self.nodes_interface_master[nodes_mask_master]
            self.nodes_interface_slave = self.nodes_interface_slave[nodes_mask_slave]
            self.positions_master = self.positions_master[nodes_mask_master]
            self.positions_slave = self.positions_slave[nodes_mask_slave]

            if np.all(nodes_mask_master) and np.all(nodes_mask_slave):
                break

        self.dofs_master = nodes_to_dofs(self.nodes_interface_master, dim=self.dim).ravel()
        self.dofs_slave = nodes_to_dofs(self.nodes_interface_slave, dim=self.dim).ravel()
        self.n_constraint_dofs = len(self.dofs_master)

    def assemble_interpolation_matrices(self):
        """Assemble the interpolation matrices $\R_{NM}$ (see [1], p. 53)"""
        self.R12 = construct_interpolation_matrix(self.positions_master, self.positions_slave, self.radii_slave, self.rbf)
        self.R21 = construct_interpolation_matrix(self.positions_slave, self.positions_master, self.radii_master, self.rbf)

    def assemble_interface_mass_matrices(self):
        """Assemble the interface mass matrices $M_1, M_2$ (see [1], eq. 12)"""
        self.model.assembleMass()
        mass_matrix = self.model.getDOFManager().getMatrix('M')
        mass_matrix = aka.AkantuSparseMatrix(mass_matrix).toarray()
        mass_matrix = sp.sparse.csr_matrix(mass_matrix)

        self.M1 = mass_matrix[np.ix_(self.dofs_master, self.dofs_master)]
        self.M2 = mass_matrix[np.ix_(self.dofs_slave, self.dofs_slave)]

    def assemble_stiffness_matrix(self):
        """Assemble the global stiffness matrix $M_1, M_2$ (see [1], p. 54)"""
        self.model.assembleStiffnessMatrix()
        self.K = self.model.getDOFManager().getMatrix('K')
        self.K = aka.AkantuSparseMatrix(self.K).toarray()
        self.K = sp.sparse.csr_matrix(self.K)

    def assemble_B_matrices(self):
        """Assemble block components of INTERNODES matrix (see [1], eq. 13)"""
        sorter_free_dofs = np.argsort(self.dofs_free)
        indx1i = sorter_free_dofs[np.searchsorted(self.dofs_free, self.dofs_master, sorter=sorter_free_dofs)]
        indx2i = sorter_free_dofs[np.searchsorted(self.dofs_free, self.dofs_slave, sorter=sorter_free_dofs)]

        self.B = sp.sparse.csr_matrix((self.n_dofs_free, self.n_constraint_dofs), dtype=np.float64)
        self.B_tilde = sp.sparse.csr_matrix((self.n_constraint_dofs, self.n_dofs_free), dtype=np.float64)

        self.B[indx1i, :] = - self.M1
        self.B[indx2i, :] = self.M2 * expand_to_dim(self.R21, dim=self.dim)

        self.B_tilde[:, indx1i] = sp.sparse.eye(self.n_constraint_dofs, dtype=np.float64, format="csr")
        self.B_tilde[:, indx2i] = - expand_to_dim(self.R12, dim=self.dim)

    def assemble_internodes_matrix(self):
        """Assemble the INTERNODES matrix (see [1], eq. 13)"""
        K_free = self.K[np.ix_(self.dofs_free, self.dofs_free)] / self.rescaling_factor
        self.internodes_matrix = sp.sparse.vstack([
            sp.sparse.hstack([K_free, self.B]),
            sp.sparse.hstack([self.B_tilde, sp.sparse.csr_matrix(
                (self.B_tilde.shape[0], self.B.shape[1]), dtype=np.float64
            )])
        ])

    def assemble_force_term(self, f_free, displacements):
        """Assemble the force term (see [1], eq. 13)
        
        Parameters
        ----------
        f_free : np.ndarray
            Force applied to free DOFs
        displacements : np.ndarray
            Displacements of the free DOFs
        """
        dirichlet_displacements = self.K[np.ix_(self.dofs_free, self.dofs_blocked)] * displacements[self.dofs_blocked]
        self.force_term = np.concatenate([(f_free - dirichlet_displacements) / self.rescaling_factor,
                                          (self.R12 * self.positions_slave - self.positions_master).ravel()])

    def solve_direct(self, displacements):
        """Solve the INTERNODES system of equations, Ax = b (see [1], eq. 13)
        
        Parameters
        ----------
        displacements : np.ndarray
            Displacements of the free DOFs before the solve step 

        Returns
        -------
        positions_new : np.ndarray
            Positions after the solve step 
        displacements : np.ndarray
            Displacements of the free DOFs after the solve step 
        lambdas : np.ndarray
            Lagrange multipliers after the solve step
        """
        x = sp.sparse.linalg.spsolve(self.internodes_matrix, self.force_term)
        displacements[self.dofs_free] = x[:self.n_dofs_free]
        positions_new = self.mesh.getNodes() + displacements.reshape((-1, self.dim))
        lambdas = self.rescaling_factor * x[self.n_dofs_free:].reshape((-1, self.dim))

        return positions_new, displacements, lambdas

    def remove_traction(self, positions_new, lambda1):
        # TODO: Document
        normals1b = self.compute_normals(positions_new, self.nodes_boundary_master, self.connectivity_boundary_master, self.connectivity_boundary_body_master)
        normals2b = self.compute_normals(positions_new, self.nodes_boundary_slave, self.connectivity_boundary_slave, self.connectivity_boundary_body_slave)

        self.normals_interface_master = normals1b[np.in1d(self.nodes_boundary_master, self.nodes_interface_master)]
        self.normals_interface_slave = normals2b[np.in1d(self.nodes_boundary_slave, self.nodes_interface_slave)]

        lambda2 = - self.R21 * lambda1

        scalar1 = np.sum(lambda1 * self.normals_interface_master, axis=1)
        scalar2 = np.sum(lambda2 * self.normals_interface_slave, axis=1)

        nodes1i_dump = self.nodes_interface_master[scalar1 > 0]
        nodes2i_dump = self.nodes_interface_slave[scalar2 > 0]

        if len(nodes1i_dump) == len(nodes2i_dump) == 0:
            nodes1i_add, nodes2i_add = self.detect_gaps(positions_new)

            self.nodes_interface_master, diff_nb_nodes1i = self.update_interface(nodes1i_add, self.nodes_interface_master, 'add')
            self.nodes_interface_slave, diff_nb_nodes2i = self.update_interface(nodes2i_add, self.nodes_interface_slave, 'add')

            print(diff_nb_nodes1i, ' nodes added to interface 1')
            print(diff_nb_nodes2i, ' nodes added to interface 2')
        else:
            self.nodes_interface_master, diff_nb_nodes1i = self.update_interface(nodes1i_dump, self.nodes_interface_master, 'dump')
            self.nodes_interface_slave, diff_nb_nodes2i = self.update_interface(nodes2i_dump, self.nodes_interface_slave, 'dump')

            print(diff_nb_nodes1i, ' nodes removed from interface 1')
            print(diff_nb_nodes2i, ' nodes removed from interface 2')

        self.positions_master = self.mesh.getNodes()[self.nodes_interface_master]
        self.positions_slave = self.mesh.getNodes()[self.nodes_interface_slave]

        converged = (diff_nb_nodes1i == diff_nb_nodes2i == 0)
        return converged

    def update_interface(self, new_nodes, nodesi, case):
        # TODO: Document and maybe merge with remove_traction
        if case == 'dump':
            nodesi_new = nodesi[~np.in1d(nodesi, new_nodes)]
        elif case == 'add':
            nodesi_new = np.union1d(nodesi, new_nodes)

        diff_nb_nodes = len(nodesi) - len(nodesi_new)
        return nodesi_new, diff_nb_nodes

    def compute_normals(self, positions_new, nodesb, connectivity_boundary, connectivity_boundary_body, step_size=1e-3):
        # TODO: Document and refactor
        # TODO: This only works in 2D -> Rewrite for arbitrary dimension
        connectivity_interface_body = connectivity_boundary_body[np.in1d(connectivity_boundary_body, nodesb).reshape(connectivity_boundary_body.shape).any(axis=1)]
        nodesb_body = np.unique(connectivity_interface_body[~np.isin(connectivity_boundary_body, nodesb)])

        tangents = positions_new[connectivity_boundary[:, 1]] - positions_new[connectivity_boundary[:, 0]]
        lengths = np.linalg.norm(tangents, axis=1)[:, np.newaxis]
        tangents /= lengths

        normals = np.zeros_like(tangents)
        normals[:, 0] = -tangents[:, 1]
        normals[:, 1] = tangents[:, 0]

        n = len(nodesb)
        normals_avg = np.zeros((n, self.dim))
        for j in range(n):
            node = nodesb[j]
            coord = positions_new[node, :]
            id = np.in1d(connectivity_boundary, node).reshape(connectivity_boundary.shape).any(axis=1)
            length = lengths[id]
            normal_avg = 1 / np.sum(length) * np.sum(normals[id, :]*length, axis=0)

            tang_plus = (coord + step_size*normal_avg).reshape((-1, self.dim))
            tang_minus = (coord - step_size*normal_avg).reshape((-1, self.dim))

            min_plus = sp.spatial.distance.cdist(positions_new[nodesb_body, :], tang_plus).min()
            min_minus = sp.spatial.distance.cdist(positions_new[nodesb_body, :], tang_minus).min()

            normals_avg[j, :] = normal_avg * (1 if min_plus > min_minus else -1)

        normals_avg /= np.linalg.norm(normals_avg, axis=1)[:, np.newaxis]
        return normals_avg

    def detect_gaps(self, positions_new, tolerance=0.9, mesh_size=0.05):
        # TODO: Document
        self.positions_master = positions_new[self.nodes_interface_master, :]
        self.positions_slave = positions_new[self.nodes_interface_slave, :]

        self.find_contact_nodes()

        diffs1 = self.R12 * self.positions_slave - self.positions_master
        diffs2 = self.R21 * self.positions_master - self.positions_slave

        scalar1 = np.sum(diffs1 * self.normals_interface_master, axis=1)
        scalar2 = np.sum(diffs2 * self.normals_interface_slave, axis=1)

        threshold = - tolerance * mesh_size

        nodes1i_add = self.nodes_interface_master[scalar1 < threshold]
        nodes2i_add = self.nodes_interface_slave[scalar2 < threshold]

        return nodes1i_add, nodes2i_add
