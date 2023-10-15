-- 2D Transport test with Vacuum BC and Isotropic external source.

--############################################### Setup mesh
meshgen = chi_mesh.MeshGenerator.Create
({
  inputs =
  {
    chi_mesh.FromFileMeshGenerator.Create
    ({
      filename="./meshes/tri_simple_square.obj"
    }),
  },
  partitioner = chi.KBAGraphPartitioner.Create
  ({
    nx = 2, ny=2, nz=1,
    xcuts = {0.5}, ycuts = {0.5},
  })
})

chi_mesh.MeshGenerator.Execute(meshgen)
chiMeshHandlerExportMeshToVTK("ZObjMesh")
--############################################### Set Material IDs

--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("mat1");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
  CHI_XSFILE,"xs_1g_MatA.cxs")

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
num_groups = 1
src={}
for g=1,num_groups do
  src[g] = 1.0
end
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--############################################### Setup Physics
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 1)
chiOptimizeAngularQuadratureForPolarSymmetry(pquad, 4.0*math.pi)

lbs_block =
{
  num_groups = num_groups,
  groupsets =
  {
    {
      groups_from_to = {0, 0},
      angular_quadrature_handle = pquad,
      angle_aggregation_num_subsets = 1,
      inner_linear_method = "gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 300,
      gmres_restart_interval = 30,
    }
  }
}

lbs_options =
{
  scattering_order = 0,
}

phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
lbs.SetOptions(phys, lbs_options)

--############################################### Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({lbs_solver_handle = phys})

chiSolverInitialize(ss_solver)
chiSolverExecute(ss_solver)

--############################################### Get field functions
fflist,count = chiLBSGetScalarFieldFunctionList(phys)
chiExportFieldFunctionToVTK(fflist[1],"T2D1g_read")
