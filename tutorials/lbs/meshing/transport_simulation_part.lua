--[[ @doc
\page readmesh_rest The Rest of the Transport Simulation Input
\tableofcontents

\section mat_rest Materials and Sources
We create two materials and add two properties to it:
+ TRANSPORT_XSECTIONS for the transport cross sections, and
+ ISOTROPIC_MG_SOURCE for the isotropic volumetric source

We assign cross sections per material by loading the file containing the cross sections.

We create lua tables containing the volumetric multigroup sources and assign them to their
 respective material by passing an array.

-- @end ]]
--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("material-1");
materials[2] = chiPhysicsAddMaterial("material-2");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
  CHI_XSFILE,"xs_1g_MatA.cxs")
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)
chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
  CHI_XSFILE,"xs_1g_MatB.cxs")

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)
num_groups = 1
src1={}
src2={}
for g=1,num_groups do
  src1[g] = 0.0
  src2[g] = 300.0
end
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src1)
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src2)

--[[ @doc
\section ang_rest Angular Quadrature

We call a product Gauss-Legendre-Chebyshev quadrature and pass the number of **positive** polar cosines (here ```npolar = 2```) and the number of azimuthal subdivisions in **one quadrant** (```nazimu = 1```). This creates a 3D angular quadrature.

We finish by optimizing the quadrature to only use the positive hemisphere for 2D simulations.
-- @end ]]
--############################################### Setup the Angular Quadrature
nazimu = 4
npolar = 2
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,nazimu,npolar)
chiOptimizeAngularQuadratureForPolarSymmetry(pquad, 4.0*math.pi)

--############################################### Setup LBS parameters
--[[ @doc
\section lbs_solver_rest Linear Boltzmann Solver
\subsection options_solver_rest Options for the Linear Boltzmann Solver (LBS)
In the LBS block, we provide
+ the number of energy groups,
+ the groupsets (with 0-indexing), the handle for the angular quadrature, the angle aggregation, the solver type, tolerances, and other solver options.
-- @end ]]
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
--[[ @doc
\subsection more_options_rest Further Options for the Linear Boltzmann Solver
In the LBS options, we pass the maximum scattering order to be employed (should be less than the one supplied the cross section file)
-- @end ]]
lbs_options =
{
  scattering_order = 0,
}
--[[ @doc
\subsection putting_together_rest Putting the Linear Boltzmann Solver Together
We create the physics solver, initialize it, and execute it,
-- @end ]]
phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
lbs.SetOptions(phys, lbs_options)

--############################################### Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({lbs_solver_handle = phys})

chiSolverInitialize(ss_solver)
chiSolverExecute(ss_solver)

--[[ @doc
\section postprocessing_rest Post-Processing via Field Functions
We extract the scalar flux (i.e., the first entry in the field function list; recall that lua indexing starts at 1) and export it to a VTK file whose name is supplied by the user.

- @end ]]
--############################################### Get field functions
fflist,count = chiLBSGetScalarFieldFunctionList(phys)
vtk_basename = "read_obj"
chiExportFieldFunctionToVTK(fflist[1],vtk_basename)

