import sys
sys.path.append("/home/bp/akantu/build/python/")
sys.path.append("/home/bp/akantu/test")
sys.path.append("/home/bp/akantu/test/test_fe_engine")
sys.path.append("/home/bp/akantu/build/src")

import numpy as np
import akantu as aka
import pytest
import matplotlib.pyplot as plt

# with this "hack", we also support running from tests where all the files are in one folder
import os.path
import sys

file_prefix = ""

if os.path.isdir("prototype_internodes"):
    sys.path.append("prototype_internodes")
    file_prefix = "prototype_internodes/"

def test_contact2d_circle():
    mesh_file = 'contact2d_circle.msh'
    material_file = file_prefix + 'material.dat'

    aka.parseInput(material_file)
    spatial_dimension = 2

    mesh = aka.Mesh(spatial_dimension)
    mesh.read(mesh_file)
    triangles = mesh.getConnectivity(aka._triangle_3)

    contact = aka.ContactMechanicsInternodesModel(mesh)
    detector = contact.getContactDetectorInternodes()
    solid = contact.getSolidMechanicsModel()

    def print_node_groups():
        print("Current master node group: ", np.array(detector.getMasterNodeGroup().getNodes()).flatten())
        print("Current slave node group: ", np.array(detector.getSlaveNodeGroup().getNodes()).flatten())

    print("Initial master node group: ", np.array(detector.getInitialMasterNodeGroup().getNodes()).flatten())
    print("Initial slave node group: ", np.array(detector.getInitialSlaveNodeGroup().getNodes()).flatten())
    #print_node_groups()
    contact.initFull(_analysis_method=aka._static)
    #print_node_groups()

    # boundary conditions
    contact.applyBC(aka.FixedValue(0., aka._x), 'primary_fixed')
    contact.applyBC(aka.FixedValue(0., aka._y), 'primary_fixed')
    contact.applyBC(aka.FixedValue(0., aka._x), 'secondary_fixed')
    contact.applyBC(aka.FixedValue(-0.1, aka._y), 'secondary_fixed')

    # setup debugging
    def plot_mesh():
        print_node_groups()

        plot_mesh.iter += 1

        current_positions = contact.getSolidMechanicsModel().getCurrentPosition()

        plt.triplot(current_positions[:, 0], current_positions[:, 1], triangles, label="mesh")
        interface_nodes = np.zeros(len(current_positions), dtype=bool)
        interface_nodes[np.array(detector.getMasterNodeGroup().getNodes()).flatten()] = 1
        interface_nodes[np.array(detector.getSlaveNodeGroup().getNodes()).flatten()] = 1
        plt.scatter(current_positions[interface_nodes, 0], current_positions[interface_nodes, 1], color="blue", label="interface")
        plt.scatter(current_positions[8, 0], current_positions[8, 1], color='red', label="8")
        #plt.scatter(current_positions[66, 0], current_positions[66, 1], color='green', label="66")
        plt.legend()
        plt.title("After iteration %d" % plot_mesh.iter)
        ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
        plt.show()
    plot_mesh.iter = 0

    contact.setStartIterationCallback(plot_mesh)

    # Actually solve
    contact.solveStep()

    # Plot final mesh
    plot_mesh()

if __name__ == '__main__':
    #pytest.main()
    test_contact2d_circle()