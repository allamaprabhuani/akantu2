import akantu as akantu
import time

spatial_dimension = 2
akantu.parseInput('./materials/detection-explicit.dat')

mesh = akantu.Mesh(spatial_dimension)
mesh.read('./mesh/detection-explicit.msh')

model = akantu.ContactMechanicsModel(mesh)

model.initFull(akantu.ContactMechanicsModelOptions(akantu._explicit_lumped_mass))
surface_selector = akantu.PhysicalSurfaceSelector(mesh)
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
