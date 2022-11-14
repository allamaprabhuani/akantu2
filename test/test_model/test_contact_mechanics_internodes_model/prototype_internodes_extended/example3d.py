import akantu as aka
import numpy as np
import matplotlib.pyplot as plt
from contact_mechanics_internodes import ContactMechanicsInternodes

def plot_mesh(positions, triangle_indices):
    plt.figure()
    plt.triplot(positions[:, 0], positions[:, 1], triangle_indices)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('scaled')
    plt.show()

mesh_file = 'contact3d.msh'
material_file = 'material.dat'
spatial_dimension = 3
aka.parseInput(material_file)

mesh = aka.Mesh(spatial_dimension)
mesh.read(mesh_file)

model = aka.SolidMechanicsModel(mesh)
model.initFull(_analysis_method=aka._implicit_dynamic)

model.applyBC(aka.FixedValue(0., aka._x), 'lower_bottom')
model.applyBC(aka.FixedValue(0., aka._y), 'lower_bottom')
model.applyBC(aka.FixedValue(0., aka._z), 'lower_bottom')
model.applyBC(aka.FixedValue(0., aka._x), 'upper_top')
model.applyBC(aka.FixedValue(-0.01, aka._y), 'upper_top')
model.applyBC(aka.FixedValue(0., aka._z), 'upper_top')

nodes_top = mesh.getElementGroup('upper_top').getNodeGroup().getNodes().ravel()
displacements = np.zeros((mesh.getNbNodes(), spatial_dimension))
displacements[nodes_top, 1] = -0.01
displacements = displacements.ravel()

internodes_model = ContactMechanicsInternodes(spatial_dimension, mesh, model, 'lower_top', 'upper_bottom', blocked_nodes_name='blocked_nodes')

f_free = model.getExternalForce().ravel()
f_free = f_free[internodes_model.dofs_free]

# Plot initial configuration
plot_mesh(internodes_model.mesh.getNodes(), internodes_model.mesh.getConnectivity(aka._triangle_3))

for i in range(10):
    print("----> Starting iteration", i+1, "<----")

    # Find the interface nodes
    internodes_model.find_interface_nodes()
    print("Interface nodes primary: ", internodes_model.nodes_interface_primary)
    print("Interface nodes secondary: ", internodes_model.nodes_interface_secondary)

    # Assemble model
    internodes_model.assemble_interpolation_matrices()
    internodes_model.assemble_stiffness_matrix()
    internodes_model.assemble_interface_mass_matrices()
    internodes_model.assemble_B_matrices()
    internodes_model.assemble_internodes_matrix()
    internodes_model.assemble_force_term(f_free, displacements)

    # Solve model
    positions_new, displacements, lambdas = internodes_model.solve_direct(displacements)

    # Plot the obtained solution
    plot_mesh(positions_new, internodes_model.mesh.getConnectivity(aka._triangle_3))

    # Update the interface nodes and check if it converged
    converged = internodes_model.update_interface(positions_new, lambdas)

    if converged:
        print('\nsuccessfully converged in', i+1, 'iterations')
        break