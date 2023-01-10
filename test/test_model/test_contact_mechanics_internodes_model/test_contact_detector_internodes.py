import pytest
import numpy as np
import akantu as aka
import scipy as sp

def _wendland_rbf(distances, radiuses):
    """
    Evaluate the wendland radial basis function.

        phi(d, r) = max{(1 - d/r)^4 * (1 + 4*d/r), 0}

    Parameters
    ----------
    distances : np.ndarray
        Distances between nodes.
    radiuses : np.ndarray
        Radius parameters for radial basis function.
    """
    delta = distances / radiuses
    return (1 - delta) ** 4 * (1 + 4 * delta) * (delta <= 1)

def _check_rbf_radius_conditions(positions, rbf_radius_parameters, c_min=0.05, c_max=0.95, n_samples=100):
    """
    Checks the conditions of the radial basis function are satisfied.

    Parameters
    ----------
    positions : np.ndarray
        Positions of interface nodes.
    rbf_radius_parameters : np.ndarray
        Radius parameters corresponding to the interface nodes.
    c_min : float 
        Minimum value to start search for 'c' satisfying the conditions.
        'c' defined as used in [1], Page 51, Section 2.3, Equation 2.
    c_max : float
        Minimum value to stoop search for 'c' satisfying the conditions.
    n_samples : int
        Number of uniform steps from c_min to c_max in search for 'c'.

    Reference
    ---------
    [1] Y. Voet et. al.: The INTERNODES method for applications in contact mechanics
        and dedicated preconditioning techniques.
        Computers & Mathematics with Applications, vol. 127, 2022, pp. 48-64
        https://doi.org/10.1016/j.camwa.2022.09.019
    """
    dist_MM = sp.spatial.distance.cdist(positions, positions)
    np.fill_diagonal(dist_MM, np.inf)
    min_ratio = np.min(dist_MM / rbf_radius_parameters, axis=1)
    n_supports = np.sum(dist_MM < rbf_radius_parameters, axis=0)

    for c in np.linspace(c_min, c_max, n_samples):

        # Check condition: [1], Page 51, Section 2.3, Equation 2
        condition_1 = np.all(c*np.ones_like(min_ratio) <  min_ratio)

        # Check condition: [1], Page 51, Section 2.3, Equation 4
        n_supports_max = np.ones_like(n_supports) / _wendland_rbf(c, 1)
        condition_2 = np.all(n_supports < n_supports_max)

        if condition_1 and condition_2:
            return

    assert False, "No 'c' was found which satisfies radius parameter conditions"

def _check_node_contained_in_support(positions, positions_ref, rbf_radius_parameters, C_max=0.99):
    """
    Checks the conditions of the radial basis function are satisfied.

    Parameters
    ----------
    positions : np.ndarray
        Positions of interface nodes.
    positions_ref : np.ndarray
        Position where the interpolant will be evaluated (opposing interface).
    rbf_radius_parameters : np.ndarray
        Radius parameters corresponding to the interface nodes.
    C_max : float 
        Maximum value constant 'C' in [1], Page 51, Section 2.3, Equation 3

    Reference
    ---------
    [1] Y. Voet et. al.: The INTERNODES method for applications in contact mechanics
        and dedicated preconditioning techniques.
        Computers & Mathematics with Applications, vol. 127, 2022, pp. 48-64
        https://doi.org/10.1016/j.camwa.2022.09.019
    """
    # Check condition: [1], Page 51, Section 2.3, Equation 3
    dist_MN = sp.spatial.distance.cdist(positions_ref, positions)
    min_ratio = np.min(dist_MN / rbf_radius_parameters, axis=1)
    np.testing.assert_array_less(min_ratio, C_max*np.ones_like(min_ratio))

def test_findContactNodes_2d():
    """
    Test if the node search finds radius parameters and interface nodes that
    satisfy the conditions in [1] for a 2d problem of a plane interfacing with
    a circle.

    Reference
    ---------
    [1] Y. Voet et. al.: The INTERNODES method for applications in contact mechanics
        and dedicated preconditioning techniques.
        Computers & Mathematics with Applications, vol. 127, 2022, pp. 48-64
        https://doi.org/10.1016/j.camwa.2022.09.019
    """
    # Set up contact problem
    mesh_file = 'contact2d_plane_circle.msh'
    material_file = 'material.dat'
    spatial_dimension = 2
    aka.parseInput(material_file)

    # Read mesh
    mesh = aka.Mesh(spatial_dimension)
    mesh.read(mesh_file)

    # Find contact nodes using detector
    detector = aka.ContactDetectorInternodes(mesh)
    detector.findContactNodes(detector.getMasterNodeGroup(), detector.getSlaveNodeGroup())

    # Get interface nodes
    nodes_interface_master = detector.getMasterNodeGroup().getNodes().ravel()
    nodes_interface_slave = detector.getSlaveNodeGroup().getNodes().ravel()

    # Get radius parameters of radial basis functions
    radiuses_interface_master = detector.getMasterRadiuses().ravel()
    radiuses_interface_slave = detector.getSlaveRadiuses().ravel()

    # Get nodal positions
    positions = mesh.getNodes()

    # Check if radius parameters satisfy the conditions from the paper
    _check_rbf_radius_conditions(positions[nodes_interface_master], radiuses_interface_master)
    _check_rbf_radius_conditions(positions[nodes_interface_slave], radiuses_interface_slave)

    # Check if all reference nodes are supported by an interpolation node
    _check_node_contained_in_support(positions[nodes_interface_master], positions[nodes_interface_slave], radiuses_interface_master)
    _check_node_contained_in_support(positions[nodes_interface_slave], positions[nodes_interface_master], radiuses_interface_slave)

def test_findContactNodes_3d():
    """
    Test if the node search finds radius parameters and interface nodes that
    satisfy the conditions in [1] for a 3d problem of a plane interfacing with
    a sphere.

    Reference
    ---------
    [1] Y. Voet et. al.: The INTERNODES method for applications in contact mechanics
        and dedicated preconditioning techniques.
        Computers & Mathematics with Applications, vol. 127, 2022, pp. 48-64
        https://doi.org/10.1016/j.camwa.2022.09.019
    """
    # Set up contact problem
    mesh_file = 'contact3d_plane_sphere.msh'
    material_file = 'material.dat'
    spatial_dimension = 3
    aka.parseInput(material_file)

    # Read mesh
    mesh = aka.Mesh(spatial_dimension)
    mesh.read(mesh_file)

    # Find contact nodes using detector
    detector = aka.ContactDetectorInternodes(mesh)
    detector.findContactNodes(detector.getMasterNodeGroup(), detector.getSlaveNodeGroup())

    # Get interface nodes
    nodes_interface_master = detector.getMasterNodeGroup().getNodes().ravel()
    nodes_interface_slave = detector.getSlaveNodeGroup().getNodes().ravel()

    # Get radius parameters of radial basis functions
    radiuses_interface_master = detector.getMasterRadiuses().ravel()
    radiuses_interface_slave = detector.getSlaveRadiuses().ravel()

    # Get nodal positions
    positions = mesh.getNodes()

    # Check if radius parameters satisfy the conditions from the paper
    _check_rbf_radius_conditions(positions[nodes_interface_master], radiuses_interface_master)
    _check_rbf_radius_conditions(positions[nodes_interface_slave], radiuses_interface_slave)

    # Check if all reference nodes are supported by an interpolation node
    _check_node_contained_in_support(positions[nodes_interface_master], positions[nodes_interface_slave], radiuses_interface_master)
    _check_node_contained_in_support(positions[nodes_interface_slave], positions[nodes_interface_master], radiuses_interface_slave)

if __name__ == '__main__':
    pytest.main()