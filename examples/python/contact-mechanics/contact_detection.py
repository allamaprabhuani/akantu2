import py11_akantu as akantu

spatial_dimension = 2
akantu.parseInput('./material_press.dat')

mesh = akantu.Mesh(spatial_dimension)
mesh.read('./press_2d.msh')

detector = akantu.ContactMechanicsModel(mesh)

detector.initFull(akantu.ContactMechanicsModelOptions(akantu._explicit_dynamic_contact))
surface_selector = akantu.PhysicalSurfaceSelector(mesh)
