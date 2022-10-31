"""
References
----------

[1] Y. Voet et. al.: The INTERNODES method for applications in contact mechanics
    and dedicated preconditioning techniques.
    Computers & Mathematics with Applications, vol. 127, 2022, pp. 48-64
    https://doi.org/10.1016/j.camwa.2022.09.019
[2] Y. Voet: On the preconditioning of the INTERNODES matrix for applications
    in contact mechanics.
    Master's thesis EPFL, 2021.
"""

import numpy as np
import scipy as sp
import akantu as aka

def nodes_to_dofs(nodes, dim):
    """Obtain the DOFs that correspond to the node indices
    
    Parameters
    ----------
    nodes : 1D list or np.ndarray
        List of node indices.
    dim : int
        Spatial dimension.
    
    Returns
    -------
    DOFs : np.ndarray
        Array with the degrees of freedom corresponding to the nodes.

    Example
    -------
    >>> nodes =  np.array([1, 2, 5])
    >>> nodes_to_dofs(nodes, dim=2)
    >>>   = array([2, 3, 4, 5, 10, 11])
    """
    return dim*np.repeat(nodes, dim) + np.tile(range(dim), len(nodes))

def expand_to_dim(matrix, dim):
    """Expand matrix entries to (dim x dim)-blocks corresponding to the DOFs.

    Parameters
    ----------
    matrix : 2D np.ndarray, sparse or dense scipy matrix
        A matrix.
    dim : int
        The spatial dimension of the problem.

    Returns
    -------
    matrix_expanded : sparse scipy matrix
        The expanded matrix.

    Example
    -------
    >>> matrix = np.array([[1., 2.],
    >>>                    [3., 4.]])
    >>> expand_to_dim(matrix, dim=2).todense()
    >>>   = matrix([[1., 0., 2., 0.],
    >>>             [0., 1., 0., 2.],
    >>>             [3., 0., 4., 0.],
    >>>             [0., 3., 0., 4.]])
    """
    return sp.sparse.kron(matrix, np.eye(dim), format='csr')

def remove_rows_without_items(matrix, items):
    """Remove rows from matrix that are not contained in items.

    Parameters
    ----------
    matrix : 2D np.ndarray
        A matrix.
    items : 1D np.ndarray or list
        A list with items.

    Returns
    -------
    matrix_removed : 2D np.ndarray
        The matrix where all rows without an item have been eliminated.

    Example
    -------
    >>> matrix = np.array([[1, 2],
    >>>                    [3, 4],
    >>>                    [5, 6],
    >>>                    [7, 8]])
    >>> items = [3, 8]
    >>> remove_rows_without_elements(matrix, items)
    >>>   = array([[3, 4],
                   [7, 8]])
    """
    return matrix[np.isin(matrix, items).any(axis=1)]

def wendland(delta):
    """Evaluate Wendland C2 function.
    
    Reference
    ---------
    [1], Page 49, Section 2.1, Table 1
    """
    return (1 - delta)**4 * (1 + 4*delta) * (delta <= 1)

def wendland_rbf(distances, radiuses):
    """Evaluate Wendland radial basis function.
    
    Reference
    ---------
    [1], Page 49, Section 2.1, Bottom left
    """
    return wendland(distances / radiuses)

def construct_rbf_matrix(positions, positions_ref, radiuses_ref, rbf):
    """Construct radial basis matrix $\Phi_{NM}$.
    
    Reference
    ---------
    [1], Page 49, Section 2.1, Point 1
    """
    return rbf(sp.spatial.distance.cdist(positions, positions_ref), radiuses_ref)

def construct_interpolation_matrix(positions, positions_ref, radiuses_ref, rbf):
    """Construct interpolation matrix $\R_{NM}$.

    Reference
    ---------
    [1], Page 49, Section 2.1, Bottom right
    """
    rbf_matrix_MM = construct_rbf_matrix(positions_ref, positions_ref, radiuses_ref, rbf)
    rbf_matrix_NM = construct_rbf_matrix(positions, positions_ref, radiuses_ref, rbf)

    # Compute raw interpolation matrix without rescaling
    interpolation_matrix = np.linalg.solve(rbf_matrix_MM.T, rbf_matrix_NM.T).T

    # Compute the diagonal rescaling factors (diagonal entries of $D_{NN}$)
    rescaling_factors = np.sum(interpolation_matrix, axis=1)[:, np.newaxis]
    return sp.sparse.csr_matrix(interpolation_matrix / rescaling_factors)

def compute_rbf_radius_parameters(positions, positions_ref, c=0.5, C=0.95, attack_radius_tolerance=0.05, max_iter=10):
    """Compute radius parameters for radial basis functions.

    Parameters
    ----------
    lower_distance_bound : float in (0, 1), default is 0.5
        [1], Page 51, Section 2.3, Equation 2
        (The empirical default value is given at the right on same page)
    upper_distance_bound : float in (c, 1), default is 0.5
        [1], Page 51, Section 2.3, Equation 3
        (The empirical default value is given at the bottom right on same page)
    attack_radius_tolerance : float, default is 0.05
        TODO: makes no sense...
    max_iter : int, default is 10
        Maximum number of iterations for finding suitable radius parameters.

    Returns
    -------
    rbf_radius_parameters : np.ndarray
        The radius parameters to use in the radial basis function.
    n_nnz_elements_PHI_NM
        TODO: No idea what this is.

    Reference
    ---------
    TODO: Where is this algorithm exactly coming from?! Heuristic?
    """

    # Iteratively increase parameter 'c' until good radius parameters are found
    for _ in range(max_iter):

        # Criterion for strict diagonal dominance by rows
        # [1], Page 51, Section 2.3, Equation 4
        max_n_supports = np.floor(1 / wendland(c))

        # Compute distance matrices between positions and reference positions
        distance_matrix_MM = sp.spatial.distance.cdist(positions, positions)
        distance_matrix_NM = sp.spatial.distance.cdist(positions_ref, positions)

        # Set self-distance (distance between a position and itself) to infinite
        np.fill_diagonal(distance_matrix_MM, np.inf)

        # Minimum distance between two distinct positions
        min_distance_MM = np.min(distance_matrix_MM, axis=1)
    
        # Minimum distance between a position and a reference position?!?!
        # TODO: Find out what this magic does. Is 0.25 = c^2? And why not np.min?!
        min_distance_NM = (attack_radius_tolerance**2 + 0.25*min_distance_MM**2)**0.5

        # Set radius parameters to the largest value satisfying the conditions
        # Candidate radiuses are maximum between closest interpolation point
        # and TODO: what is min_distance_NM?!?!
        # [1], Page 51, Section 2.3, right center
        candidate_radiuses = np.maximum(min_distance_MM, min_distance_NM)
        rbf_radius_parameters = np.minimum(candidate_radiuses, min_distance_MM / c)

        # Number of non-zero off-diagonal elements in rows of $\Phi_{MM}$ 
        n_nnz_elements_PHI_MM = np.sum(distance_matrix_MM < rbf_radius_parameters, axis=1)

        # TODO: What exactly is this?
        n_nnz_elements_PHI_NM = np.sum(distance_matrix_NM < C*rbf_radius_parameters, axis=1)

        # Check if criterion for strict diagonal dominance by rows is satisfied
        if np.max(n_nnz_elements_PHI_MM) <= max_n_supports:
            break

        # Increase c
        c = (c + 1) / 2

    return rbf_radius_parameters, n_nnz_elements_PHI_NM

class ContactMechanicsInternodes(object):

    def __init__(self, dim, mesh, model, name_candidate_interface_master, name_candidate_interface_slave, rbf=wendland_rbf):
        self.dim = dim
        self.mesh = mesh
        self.model = model
    
        # Candidate nodes for contact interface
        self.nodes_candidate_master = mesh.getElementGroup(name_candidate_interface_master).getNodeGroup().getNodes().ravel()
        self.nodes_candidate_slave = mesh.getElementGroup(name_candidate_interface_slave).getNodeGroup().getNodes().ravel()
        
        # Nodes, positions, and corresponding dofs of master/slave interface
        self.nodes_interface_master = mesh.getElementGroup(name_candidate_interface_master).getNodeGroup().getNodes().ravel()
        self.nodes_interface_slave = mesh.getElementGroup(name_candidate_interface_slave).getNodeGroup().getNodes().ravel()
        self.positions_interface_master = mesh.getNodes()[self.nodes_interface_master]
        self.positions_interface_slave = mesh.getNodes()[self.nodes_interface_slave]
        self.dofs_interface_master = nodes_to_dofs(self.nodes_interface_master, dim=self.dim)
        self.dofs_interface_slave = nodes_to_dofs(self.nodes_interface_slave, dim=self.dim)
        self.n_dofs_interface_master = len(self.dofs_interface_master)

        # All dofs, blocked dofs, and free dobs
        self.dofs = np.arange(mesh.getNbNodes()*self.dim)
        self.dofs_blocked = self.dofs[model.getBlockedDOFs().ravel()]
        self.dofs_free = self.dofs[~model.getBlockedDOFs().ravel()]
        self.n_dofs_free = len(self.dofs_free)

        # Connectivity of the model (line segments and triangular elements)
        self.connectivity_segments = mesh.getConnectivity(aka._segment_2)
        self.connectivity_triangles = mesh.getConnectivity(aka._triangle_3)

        # Radial basis function and the corresponding radius parameters
        self.rbf = rbf
        self.rbf_radius_parameters_master = None
        self.rbf_radius_parameters_slave = None
        
        # Objects for use in the solution process
        self.rescaling_factor = 30e+9  # TODO: Tunable, Young's modulus
        self.K = None  # Stiffness matrix

        self.R12 = None  # Interpolation matrix from slave to master
        self.R21 = None  # Interpolation matrix from master to slave

        self.M1 = None  # Interface mass matrix of master
        self.M2 = None  # Interface mass matrix of slave

        self.B = None  # Block component of INTERNODES matrix
        self.B_tilde = None  # Block component of INTERNODES matrix

        self.internodes_matrix = None  # INTERNODES matrix
        self.force_term = None  # Force term b

    def find_contact_nodes(self):
        """Find contact/interface nodes while trying to satisfy the constraints.
       
        Reference
        ---------
        [1], Page 51, Section 2.3, Equation 2 and 3
        """
        while True:
            # Determine the radial basis function radius parameters
            self.rbf_radius_parameters_master, nnzR21 = compute_rbf_radius_parameters(self.positions_interface_master, self.positions_interface_slave)
            self.rbf_radius_parameters_slave, nnzR12 = compute_rbf_radius_parameters(self.positions_interface_slave, self.positions_interface_master)

            # Update interface positions
            self.positions_interface_master = self.positions_interface_master[nnzR12 > 0]
            self.positions_interface_slave = self.positions_interface_slave[nnzR21 > 0]
            
            # Update interface nodes
            self.nodes_interface_master = self.nodes_interface_master[nnzR12 > 0]
            self.nodes_interface_slave = self.nodes_interface_slave[nnzR21 > 0]
            
            # Stop algorithm if for all slave and master interface nodes
            # [1], Page 51, Section 2.3, Equation 3 is satisfied
            if np.all(nnzR12 > 0) and np.all(nnzR21 > 0):
                break

        # Update the corresponding degrees of freedom
        self.dofs_interface_master = nodes_to_dofs(self.nodes_interface_master, dim=self.dim).ravel()
        self.dofs_interface_slave = nodes_to_dofs(self.nodes_interface_slave, dim=self.dim).ravel()
        self.n_dofs_interface_master = len(self.dofs_interface_master)

    def assemble_interface_mass_matrices(self):
        """Assemble the interface mass matrices $M_1, M_2$.

        Reference
        ---------
        [1], Page 53, Section 3.4, Equation 12
        """
        self.model.assembleMass()
        mass_matrix = self.model.getDOFManager().getMatrix('M')
        mass_matrix = aka.AkantuSparseMatrix(mass_matrix).toarray()
        mass_matrix = sp.sparse.csr_matrix(mass_matrix)

        # Slice global mass matrix into master and slave matrices by dofs
        self.M1 = mass_matrix[np.ix_(self.dofs_interface_master, self.dofs_interface_master)]
        self.M2 = mass_matrix[np.ix_(self.dofs_interface_slave, self.dofs_interface_slave)]

    def assemble_interpolation_matrices(self):
        """Assemble the interpolation matrices $R_{12}$ and $R_{21}$.
        
        Reference
        ---------
        [1], Page 53, Section 3.4, Bottom right
        """
        self.R12 = construct_interpolation_matrix(self.positions_interface_master, self.positions_interface_slave, self.rbf_radius_parameters_slave, self.rbf)
        self.R21 = construct_interpolation_matrix(self.positions_interface_slave, self.positions_interface_master, self.rbf_radius_parameters_master, self.rbf)

    def assemble_B_matrices(self):
        """Assemble block components of INTERNODES matrix.

        Reference
        ---------
        [1], Page 54, Section 4, Bottom left
        """
        # Find indices in 'dofs_free' corresponding to interface of master and slave
        idx_interface_master = np.argwhere(np.in1d(self.dofs_free, self.dofs_interface_master)).squeeze()
        idx_interface_slave = np.argwhere(np.in1d(self.dofs_free, self.dofs_interface_slave)).squeeze()

        # Initialize the matrices
        self.B = sp.sparse.csr_matrix((self.n_dofs_free, self.n_dofs_interface_master), dtype=np.float64)
        self.B_tilde = sp.sparse.csr_matrix((self.n_dofs_interface_master, self.n_dofs_free), dtype=np.float64)

        # Fill in the blocks as described in the reference [1]
        self.B[idx_interface_master, :] = - self.M1
        self.B[idx_interface_slave, :] = self.M2 * expand_to_dim(self.R21, dim=self.dim)

        self.B_tilde[:, idx_interface_master] = sp.sparse.eye(self.n_dofs_interface_master, dtype=np.float64, format="csr")
        self.B_tilde[:, idx_interface_slave] = - expand_to_dim(self.R12, dim=self.dim)

    def assemble_stiffness_matrix(self):
        """Assemble the global stiffness matrix $K$."""
        self.model.assembleStiffnessMatrix()
        self.K = self.model.getDOFManager().getMatrix('K')
        self.K = aka.AkantuSparseMatrix(self.K).toarray()
        self.K = sp.sparse.csr_matrix(self.K)

    def assemble_internodes_matrix(self):
        """Assemble the INTERNODES matrix.

        Reference
        ---------
        [1], Page 54, Section 4, Equation 13
        """
        # Rescale stiffness matrix restricted to free dofs by Young's modulus
        K_free = self.K[np.ix_(self.dofs_free, self.dofs_free)] / self.rescaling_factor
        self.internodes_matrix = sp.sparse.vstack([
            sp.sparse.hstack([K_free, self.B]),
            sp.sparse.hstack([self.B_tilde, sp.sparse.csr_matrix(
                (self.B_tilde.shape[0], self.B.shape[1]), dtype=np.float64
            )])
        ])

    def assemble_force_term(self, f_free, displacements):
        """Assemble the force term.
        
        Parameters
        ----------
        f_free : np.ndarray
            Force applied to free DOFs.
        displacements : np.ndarray
            Displacements of the free DOFs.

        Reference
        ---------
        [1], Page 54, Section 4, Equation 13
        """
        # Compute displacements for Dirichlet boundary condition offset
        dirichlet_offset = self.K[np.ix_(self.dofs_free, self.dofs_blocked)] * displacements[self.dofs_blocked]
        
        # First component $f$ of force term is a adjusted by offset and rescaled
        virtual_force = (f_free - dirichlet_offset) / self.rescaling_factor

        # Nodal gaps $d$ between interpolated and true positions of master nodes
        nodal_gaps = self.R12 * self.positions_interface_slave - self.positions_interface_master
    
        self.force_term = np.concatenate([virtual_force, nodal_gaps.ravel()])

    def solve_direct(self, displacements):
        """Solve the INTERNODES system of equations.
        
        Parameters
        ----------
        displacements : np.ndarray
            Displacements of the free DOFs before the solve step 

        Returns
        -------
        positions_new : np.ndarray
            Positions after the solve step 
        displacements : np.ndarray
            Displacements $u$ of the free DOFs after the solve step 
        lambdas : np.ndarray
            Lagrange multipliers $\lambda$ after the solve step
        
        Reference
        ---------
        [1], Page 54, Section 4, Equation 13
        """
        x = sp.sparse.linalg.spsolve(self.internodes_matrix, self.force_term)

        # FIll in the computed displacements at the free dofs
        displacements[self.dofs_free] = x[:self.n_dofs_free]
        positions_new = self.mesh.getNodes() + displacements.reshape((-1, self.dim))

        # Reshape and revert the rescaling of the Lagrange multipliers $\lambda$
        lambdas = x[self.n_dofs_free:].reshape((-1, self.dim)) * self.rescaling_factor

        return positions_new, displacements, lambdas

    def update_interface(self, positions_new, lambdas):
        """Update the interfaces according to penetrations and tension
        
        Parameters
        ----------
        positions_new : np.ndarray
            New positions of nodes after solving the INTERNODES system.
        lambdas : np.ndarray
            Lagrange multipliers obtained from solving the INTERNODES system.

        Returns
        -------
        converged : bool
            If all Lagrange multipliers negative and no penetration is detected

        Reference
        ---------
        [2], Page 23, Algorithm 2, Lines 5-16
        """
        # Get the connectivities of master and slave
        connectivity_segments_master = remove_rows_without_items(self.connectivity_segments, self.nodes_interface_master)
        connectivity_segments_slave = remove_rows_without_items(self.connectivity_segments, self.nodes_interface_slave)
        connectivity_triangles_master = remove_rows_without_items(self.connectivity_triangles, self.nodes_interface_master)
        connectivity_triangles_slave = remove_rows_without_items(self.connectivity_triangles, self.nodes_interface_slave)

        # Compute the normals
        normals_candidate_master = self.compute_normals(positions_new, self.nodes_candidate_master, connectivity_segments_master, connectivity_triangles_master)
        normals_candidate_slave = self.compute_normals(positions_new, self.nodes_candidate_slave, connectivity_segments_slave, connectivity_triangles_slave)
        normals_interface_master = normals_candidate_master[np.in1d(self.nodes_candidate_master, self.nodes_interface_master)]
        normals_interface_slave = normals_candidate_slave[np.in1d(self.nodes_candidate_slave, self.nodes_interface_slave)]

        # Interpolate the Lagrange multipliers of the slave
        lambdas_slave = - self.R21 * lambdas

        # Mark nodes with positive projected Lagrange multipliers for dumping
        positive_lambda_proj_master = np.sum(lambdas * normals_interface_master, axis=1) > 0
        positive_lambda_proj_slave = np.sum(lambdas_slave * normals_interface_slave, axis=1) > 0
        nodes_to_dump_master = self.nodes_interface_master[positive_lambda_proj_master]
        nodes_to_dump_slave = self.nodes_interface_slave[positive_lambda_proj_slave]

        # If projected Lagrange multipliers are all negative
        if len(nodes_to_dump_master) == len(nodes_to_dump_slave) == 0:
            nodes_to_add_master, nodes_to_add_slave = self.detect_penetration_nodes(positions_new, normals_interface_master, normals_interface_slave)

            # If no penetration is detected, then we have convergence
            if (np.all(np.in1d(nodes_to_add_master, self.nodes_interface_master))
             or np.all(np.in1d(nodes_to_add_slave, self.nodes_interface_slave))):
                return True

            # Add new nodes to master and slave
            self.nodes_interface_master = np.union1d(self.nodes_interface_master, nodes_to_add_master)
            self.nodes_interface_slave = np.union1d(self.nodes_interface_slave, nodes_to_add_slave)

        else:
            # Dump nodes from master and slave
            self.nodes_interface_master = np.setdiff1d(self.nodes_interface_master, nodes_to_dump_master)
            self.nodes_interface_slave = np.setdiff1d(self.nodes_interface_slave, nodes_to_dump_slave)

        # Update interface positions of master and slave
        self.positions_interface_master = self.mesh.getNodes()[self.nodes_interface_master]
        self.positions_interface_slave = self.mesh.getNodes()[self.nodes_interface_slave]

        return False

    def compute_normals(self, positions_new, nodes, connectivity_segments, connectivity_triangles, step_size=1e-3):
        """Compute the normals.
        TODO: Require reference and merge 2d/3d cases.
        
        Parameters
        ----------
        positions_new : np.ndarray
            New positions of nodes after solving the INTERNODES system.
        nodes : np.ndarray 
            Nodes for which the normals are computed for.
        connectivity_segments : np.ndarray
            Connectivity of the line segments in the mesh.
        connectivity_triangles : np.ndarray
            Connectivity of the triangular elements in the mesh.
    
        Returns
        -------
        normals_avg : np.ndarray
        """
        # Determine boundary nodes that aren't interface elements
        nodesb_body = np.setdiff1d(connectivity_triangles, nodes)

        if self.dim == 2:
            tangent1 = positions_new[connectivity_segments[:, 1]] - positions_new[connectivity_segments[:, 0]]
            tangent2 = [0, 0, -1]
        elif self.dim == 3:  # Do this for only interface edges
            tangent1 = positions_new[connectivity_triangles[:, 1]] - positions_new[connectivity_triangles[:, 0]]
            tangent2 = positions_new[connectivity_triangles[:, 2]] - positions_new[connectivity_triangles[:, 0]]

        # Compute normal vectors
        normals = - np.cross(tangent1, tangent2)[:, :self.dim]

        normals_avg = np.zeros((len(nodes), self.dim))
        for j, node in enumerate(nodes):
            if self.dim == 2:
                id = np.isin(connectivity_segments, node).any(axis=1)
            elif self.dim == 3:
                id = np.isin(connectivity_triangles, node).any(axis=1)

            # Compute average normal on interface boundary
            normal_avg = np.sum(normals[id], axis=0) / np.sum(np.linalg.norm(normals[id], axis=1))

            # Extend positions by +/- a step size in the normal direction
            positions_plus = (positions_new[node] + step_size*normal_avg)
            positions_minus = (positions_new[node] - step_size*normal_avg)

            # Find minimum distance of new positions of non-interface nodes to WHAT?!?!
            min_plus = sp.spatial.distance.cdist(positions_new[nodesb_body], positions_plus.reshape((-1, self.dim))).min()
            min_minus = sp.spatial.distance.cdist(positions_new[nodesb_body], positions_minus.reshape((-1, self.dim))).min()

            # Change sign of the average normal vector
            normals_avg[j] = normal_avg * (1 if min_plus > min_minus else -1)

        normals_avg /= np.linalg.norm(normals_avg, axis=1)[:, np.newaxis]
        return normals_avg

    def detect_penetration_nodes(self, positions_new, normals_interface_master, normals_interface_slave, tolerance=0.9, mesh_size=0.05):
        """Detect the nodes on slave and master interface that penetrate.
        TODO: Require reference and make mesh size configurable!!

        Parameters
        ----------
        positions_new : np.ndarray
            New positions of nodes after solving the INTERNODES system.
        normals_interface_master : np.ndarray
            Normal vectors of master interface.
        normals_interface_slave : np.ndarray
            Normal vectors of slave interface.
        tolerance : float in (0, 1), default is 0.9
            Tolerance for what counts as penetration or not.
        mesh_size : float, default is 0.05 (TODO: No default but adaptive)
            Representative size of mesh to determine when penetration happens.
        
        Returns
        -------
        nodes_penetration_master : np.ndarray
            Nodes of master interface where penetration is observed.
        nodes_penetration_slave : np.ndarray
            Nodes of slave interface where penetration is observed.
        """

        # Update positions of master and slave nodes along interface
        self.positions_interface_master = positions_new[self.nodes_interface_master]
        self.positions_interface_slave = positions_new[self.nodes_interface_slave]

        # Find contact nodes with contact detection algorithm
        self.find_contact_nodes()

        # Determine the size of the nodal gaps on interface
        nodal_gaps_master = self.R12 * self.positions_interface_slave - self.positions_interface_master
        nodal_gaps_slave = self.R21 * self.positions_interface_master - self.positions_interface_slave

        # Penetration if projected nodal gap onto normal is sufficiently small
        threshold = - tolerance * mesh_size
        penetrates_slave = np.sum(nodal_gaps_master * normals_interface_master, axis=1) < threshold
        penetrates_master = np.sum(nodal_gaps_slave * normals_interface_slave, axis=1) < threshold

        # Add the nodes where penetration is observed to the interface
        nodes_penetration_master = self.nodes_interface_master[penetrates_slave]
        nodes_penetration_slave = self.nodes_interface_slave[penetrates_master]

        return nodes_penetration_master, nodes_penetration_slave
