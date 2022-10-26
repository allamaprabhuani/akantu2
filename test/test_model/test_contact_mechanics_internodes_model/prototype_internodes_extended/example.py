import akantu as aka
import numpy as np
import matplotlib.pyplot as plt
from contact_mechanics_internodes import ContactMechanicsInternodes

mesh_file = 'contact.msh'
material_file = 'material.dat'

aka.parseInput(material_file)
spatial_dimension = 2

mesh = aka.Mesh(spatial_dimension)
mesh.read(mesh_file)

model = aka.SolidMechanicsModel(mesh)
model.initFull(_analysis_method=aka._implicit_dynamic)

displacements = np.zeros(mesh.getNbNodes()*spatial_dimension)

model.applyBC(aka.FixedValue(0., aka._x), 'lower_bottom')
model.applyBC(aka.FixedValue(0., aka._y), 'lower_bottom')

model.applyBC(aka.FixedValue(-0.1, aka._y), 'upper_top')
nodes_top = mesh.getElementGroup('upper_top').getNodeGroup().getNodes().ravel()
displacements = displacements.reshape([-1, 2])
displacements[nodes_top, 1] = -0.1
displacements = displacements.ravel()

internodes_model = ContactMechanicsInternodes(2, mesh, model, 'lower_top', 'upper_bottom')

f_free = model.getExternalForce().ravel()
f_free = f_free[internodes_model.dofs_free]

plt.figure()
plt.triplot(internodes_model.mesh.getNodes()[:, 0], internodes_model.mesh.getNodes()[:, 1], internodes_model.mesh.getConnectivity(aka._triangle_3))
plt.title('mesh')
plt.xlabel('x')
plt.ylabel('y')
plt.axis('scaled')
plt.title('Initial')
plt.show()

#---------- internodes iterations ----------
# run until converged or nb_max_iter is attained
for i in range(10):
    print('--- iteration', i+1, '---')

    # select nodes belonging to interface
    internodes_model.find_contact_nodes()

    # assemble model matrices
    internodes_model.assemble_stiffness_matrix()
    internodes_model.assemble_interface_mass_matrices()
    internodes_model.assemble_interpolation_matrices()
    internodes_model.assemble_B_matrices()
    internodes_model.assemble_internodes_matrix()
    internodes_model.assemble_force_term(f_free, displacements)

    # solve
    positions_new, displacements, lambdas = internodes_model.solve_direct(displacements)

    plt.figure()
    plt.triplot(positions_new[:, 0], positions_new[:, 1], internodes_model.mesh.getConnectivity(aka._triangle_3))
    plt.title('mesh')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('scaled')
    plt.title('Iteration ' + str(i+1))
    plt.show()

    # add or remove nodes
    converged = internodes_model.remove_traction(positions_new, lambdas)

    if converged:
        print('\nsuccessfully converged in', i+1, 'iterations')
        break