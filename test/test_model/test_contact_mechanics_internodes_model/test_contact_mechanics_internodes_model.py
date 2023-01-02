import pytest
import numpy as np
import akantu as aka
import scipy as sp

def _get_theoretical_contact_radius(R, d):
    return (d * R) ** 0.5

def _get_theoretical_pressure_amplitude(R, d, E, nu):
    a = _get_theoretical_contact_radius(R, d)
    E_red = E / (2 * (1 - nu**2))
    F = 4 / 3 * E_red * R**0.5 * d**1.5
    return 3 * F / (2 * np.pi * a**2)

def _get_theoretical_normal_displacement(R, d, E, nu):
    a = _get_theoretical_contact_radius(R, d)
    p0 = _get_theoretical_pressure_amplitude(R, d, E, nu)
    K = np.pi / 2
    return -(1 - nu**2) / E * K * p0 * a

def test_contact2d_plane_circle():

    # Set up contact problem
    mesh_file = 'contact2d_plane_circle.msh'
    material_file = 'contact2d_material.dat'
    spatial_dimension = 2
    aka.parseInput(material_file)

    d0 = 0.05  # Initial overlap between primary and secondary mesh
    R = 0.5  # Radius of circle
    d = 0.1  # Displacement enforced on top of circle

    # Read mesh
    mesh = aka.Mesh(spatial_dimension)
    mesh.read(mesh_file)

    model = aka.ContactMechanicsInternodesModel(mesh)
    detector = model.getContactDetectorInternodes()
    solid = model.getSolidMechanicsModel()

    model.initFull(_analysis_method=aka._static)

    # Apply boundary conditions
    model.applyBC(aka.FixedValue(0., aka._x), 'primary_fixed')
    model.applyBC(aka.FixedValue(0., aka._y), 'primary_fixed')
    model.applyBC(aka.FixedValue(0., aka._x), 'secondary_fixed')
    model.applyBC(aka.FixedValue(-d+d0, aka._y), 'secondary_fixed')

    model.solveStep()

    positions = model.getSolidMechanicsModel().getCurrentPosition()
    positions_interface_secondary = positions[detector.getSlaveNodeGroup().getNodes().ravel()]

    a = np.max(sp.spatial.distance.cdist(positions_interface_secondary, positions_interface_secondary))/2
    u = np.min(positions_interface_secondary[:, 1])

    E = solid.getMaterial(0).getReal("E")
    nu = solid.getMaterial(0).getReal("nu")

    # Check if normal displacement corresponds to theoretical expectation
    np.testing.assert_allclose(_get_theoretical_normal_displacement(R, d, E, nu), u, atol=5e-3)

    # Check if contact radius corresponds to theoretical expectation
    np.testing.assert_allclose(_get_theoretical_contact_radius(R, d), a, atol=5e-2)

def test_contact2d_circle_circle():

    # Set up contact problem
    mesh_file = 'contact2d_circle_circle.msh'
    material_file = 'contact2d_material.dat'
    spatial_dimension = 2
    aka.parseInput(material_file)

    d0 = 0.05  # Initial overlap between primary and secondary mesh
    d = 0.1  # Displacement enforced on top of circle

    mesh = aka.Mesh(spatial_dimension)
    mesh.read(mesh_file)

    model = aka.ContactMechanicsInternodesModel(mesh)
    detector = model.getContactDetectorInternodes()

    model.initFull(_analysis_method=aka._static)

    # Apply boundary conditions
    model.applyBC(aka.FixedValue(0., aka._x), 'primary_fixed')
    model.applyBC(aka.FixedValue((d-d0)/2., aka._y), 'primary_fixed')
    model.applyBC(aka.FixedValue(0., aka._x), 'secondary_fixed')
    model.applyBC(aka.FixedValue((-d+d0)/2, aka._y), 'secondary_fixed')

    model.solveStep()

    positions = model.getSolidMechanicsModel().getCurrentPosition()
    positions_interface_primary = positions[detector.getMasterNodeGroup().getNodes().ravel()]
    positions_interface_secondary = positions[detector.getSlaveNodeGroup().getNodes().ravel()]

    # Check if primary interface nodes are close to y=0
    np.testing.assert_allclose(positions_interface_primary[:, 1], np.zeros(len(positions_interface_primary)), atol=1e-2)

    # Check if secondary interface nodes are close to y=0
    np.testing.assert_allclose(positions_interface_secondary[:, 1], np.zeros(len(positions_interface_secondary)), atol=1e-2)

if __name__ == '__main__':
    pytest.main()