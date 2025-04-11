-- Mesh
meshgen = mesh.MeshGenerator.Create({
  inputs = {
    mesh.FromFileMeshGenerator.Create({
      filename = "opensn.msh",
    }),
  },
})
grid = meshgen:Execute()

-- Material
Ng = 1
xs_source = xs.CreateSimpleOneGroup(0.1, 0.1)
xs_block = xs.CreateSimpleOneGroup(1.0e-6, 0.001)

-- Source
src_strength = {}
for g = 1, Ng do
  src_strength[g] = 0.0
end
src_strength[1] = 100.0
src = lbs.VolumetricSource.Create({ block_ids = {0}, group_strength = src_strength })

-- Quadrature
Npolar = 8
Nazimuthal = 16
pquad = aquad.CreateGLCProductQuadrature3DXYZ(Npolar, Nazimuthal)

-- Set up solver
lbs_block = {
  mesh = grid,
  num_groups = Ng,
  groupsets = {
    {
      groups_from_to = { 0, Ng -1 },
      angular_quadrature = pquad,
      angle_aggregation_type = "single",
      inner_linear_method = "petsc_gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 500,
    },
  },
  xs_map = {
    { block_ids = {0}, xs = xs_source },
    { block_ids = {1}, xs = xs_block },
  },
}

lbs_options = {
  scattering_order = 0,
  volumetric_sources = { src },
}

phys = lbs.DiscreteOrdinatesProblem.Create(lbs_block)
phys:SetOptions(lbs_options)
ss_solver = lbs.SteadyStateSolver.Create({ lbs_problem = phys })

-- Solve
ss_solver:Initialize()
ss_solver:Execute()

fflist,count = lbs.GetScalarFieldFunctionList(phys)
fieldfunc.ExportToVTKMulti(fflist,'flux')
