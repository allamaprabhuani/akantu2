#!/usr/bin/env python3

import akantu as aka
import time


spatial_dimension = 2
aka.parseInput('./materials/detection-explicit.dat')

mesh = aka.Mesh(spatial_dimension)
mesh.read('./mesh/detection-explicit.msh')

model = aka.ContactMechanicsModel(mesh)

model.initFull(_analysis_method=aka._explicit_lumped_mass)
surface_selector = aka.PhysicalSurfaceSelector(mesh)
model.getContactDetector().setSurfaceSelector(surface_selector)


model.setBaseName("detection-explicit")
model.addDumpFieldVector("normals")
model.addDumpField("gaps")
model.addDumpField("areas")


start_time = time.time()
model.search()
finish_time = time.time()
print('Search time = %s seconds', finish_time - start_time)

model.dump()

# by default the contact model creates a group named contact_surface
contact_surface = mesh.getElementGroup("contact_surface")
