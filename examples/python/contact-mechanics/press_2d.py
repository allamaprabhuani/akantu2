import py11_akantu as akantu

max_steps = 20000
imposing_steps = 10000
damping_interval = 10
max_displacement = 0.01
damping_ratio = 0.9

spatial_dimension = 2
akantu.parseInput('./material_press.dat')

mesh = akantu.Mesh(spatial_dimension)
mesh.read('./press_2d.msh')

coupler = akantu.CouplerSolidContact(mesh)
solid = coupler.getSolidMechanicsModel()
contact = coupler.getContactMechanicsModel()

selector = akantu.MeshDataMaterialSelectorString("physical_names", solid)

solid.initFull(akantu._explicit_lumped_mass)
contact.initFull(akantu._explicit_dynamic_contact)

solid.applyBC(akantu.FixedValue(0.0, akantu._x), "bottom")
solid.applyBC(akantu.FixedValue(0.0, akantu._y), "bottom")

coupler.initFull(akantu._explicit_dynamic_contact)

time_step = solid.getStableTimeStep()
time_step *= 0.1
coupler.setTimeStep(time_step)

coupler.setBaseName("python-press")
coupler.addDumpFieldVector("displacement")
coupler.addDumpFieldVector("normals")
coupler.addDumpFieldVector("tangents")
coupler.addDumpFieldVector("contact_force")
coupler.addDumpFieldVector("external_force")
coupler.addDumpFieldVector("internal_force")
coupler.addDumpField("gaps")
coupler.addDumpField("areas")
coupler.addDumpField("blocked_dofs")
coupler.addDumpField("grad_u")
coupler.addDumpField("stress")

coupler.dump()

velocity = solid.getVelocity()

for s in range(0, max_steps):

    if s > 0 and s <= imposing_steps/2:
        increment = max_displacement / (imposing_steps/2.0)
        solid.applyBC(akantu.IncrementValue(-increment, akantu._y), "top")

    if s >= imposing_steps/2 and s <= max_steps/2:
        increment = max_displacement / (max_steps/2.0 - imposing_steps/2.0)
        solid.applyBC(akantu.IncrementValue(-increment, akantu._x), "top")

    coupler.solveStep()

    if s % damping_interval == 0:
        velocity *= damping_ratio

    if s % 100 == 0:
        coupler.dump()
