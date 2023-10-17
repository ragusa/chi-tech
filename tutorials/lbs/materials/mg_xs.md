# Multigroup Cross Sections
___
## Background
Chi-Tech is not provided with cross-section libraries. Users are expected to supply their own multigroup cross-section data.
One may use open-source software to generate this data (e.g., NJOY, Dragon, OpenMC).

Jump to [the Chi-Tech cross-section format description](./mg_xs.md#XS)

## Mesh
A simple orthogonal 2D mesh.

We create a right parallelepiped logical volume that contains the entire mesh and we assign a 0 for material ID to all cells found inside the logical volume. Logical volumes are quite powerful, see subsequent tutorials on meshing.

```
--############################################### Setup mesh
nodes={}
n_cells=4
length=2.
xmin = -length/2.
dx = length/n_cells
for i=1,(n_cells+1) do
    k=i-1
    nodes[i] = xmin + k*dx
end

meshgen = chi_mesh.OrthogonalMeshGenerator.Create
({
    node_sets = {nodes,nodes},
})

chi_mesh.MeshGenerator.Execute(meshgen)

--############################################### Set Material IDs
vol0 = chi_mesh.RPPLogicalVolume.Create({infx=true, infy=true, infz=true})
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

```
## Materials
We create a material and add two properties to it:
+ TRANSPORT_XSECTIONS for the transport cross sections, and
+ ISOTROPIC_MG_SOURCE for the isotropic volumetric source
```
--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Material_A");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)

```
## <a name="XS"></a> Cross Sections
We assign the cross sections to the material by loading the file containing the cross sections.

The ```CHI_XSFILE``` format is as follows:
```
# Add your comments
NUM_GROUPS ng
NUM_MOMENTS nmom

SIGMA_T_BEGIN
0 value
.
.
ng-1 value
SIGMA_T_END

SIGMA_A_BEGIN
0 value
.
.
ng-1 value
SIGMA_A_END

TRANSFER_MOMENTS_BEGIN
# Add your comments (in the repo, you will typically see before each moment: Zeroth moment (l=moment) )
M_GPRIME_G_VAL 0 0 0 value
.
M_GPRIME_G_VAL moment gprime g value
.
M_GPRIME_G_VAL nmom-1 ng-1 ng-1 value
TRANSFER_MOMENTS_END

```
```
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        CHI_XSFILE,"xs_1g_MatA.cxs")

```
## Volumetric Source
We create a lua table containing the volumetric multigroup source and assign it to the material by passing that array.
```
num_groups = 1
src={}
for g=1,num_groups do
    src[g] = 1.0
end
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

```
## Angular Quadrature
We call a product Gauss-Legendre-Chebyshev quadrature and pass the number of **positive** polar cosines (here ```npolar = 2```) and the number of azimuthal subdivisions in **one quadrant** (```nazimu = 1```). This creates a 3D angular quadrature.

We finish by optimizing the quadrature to only use the positive hemisphere for 2D simulations.
```
--############################################### Setup the Angular Quadrature
nazimu = 1
npolar = 2
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,nazimu,npolar)
chiOptimizeAngularQuadratureForPolarSymmetry(pquad, 4.0*math.pi)

```
## Linear Boltzmann Solver
In the LBS block, we provide
+ the number of energy groups,
+ the groupsets (with 0-indexing), the handle for the angular quadrature, the angle aggregation, the solver type, tolerances, and other solver options.

In the LBS options, we pass the maximum scattering order to be employed (should be less than the one supplied the cross section file)

We then create the physics solver, initialize it, and execute it.

```
--############################################### Setup LBS parameters
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

```
## Post-Processing via Field Functions
We extract the scalar flux (i.e., the first entry in the field function list; recall that lua indexing starts at 1) and export it to a VTK file whose name is supplied by the user.

The resulting scalar flux is shown below:
![Scalar_flux](images/first_example_scalar_flux.png)
```
--############################################### Get field functions
fflist,count = chiLBSGetScalarFieldFunctionList(phys)
vtk_basename = "first_example"
chiExportFieldFunctionToVTK(fflist[1],vtk_basename)

```
## Possible Extensions:
1. Change the number of MPI processes (you may want to delete the safeguard at the top of the input file to run with any number of MPI ranks);
2. Change the spatial resolution by increasing or decreasing the number of cells;
3. Change the angular resolution by increasing or decreasing the number of polar and azimuthal subdivisions.
___
## The complete input is below:
You can copy/paste it or look in the file named ```./materials/mg_xs.lua```:
```
--############################################### Setup mesh
nodes={}
n_cells=4
length=2.
xmin = -length/2.
dx = length/n_cells
for i=1,(n_cells+1) do
    k=i-1
    nodes[i] = xmin + k*dx
end

meshgen = chi_mesh.OrthogonalMeshGenerator.Create
({
    node_sets = {nodes,nodes},
})

chi_mesh.MeshGenerator.Execute(meshgen)

--############################################### Set Material IDs
vol0 = chi_mesh.RPPLogicalVolume.Create({infx=true, infy=true, infz=true})
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Material_A");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)

chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        CHI_XSFILE,"xs_1g_MatA.cxs")

num_groups = 1
src={}
for g=1,num_groups do
    src[g] = 1.0
end
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--############################################### Setup the Angular Quadrature
nazimu = 1
npolar = 2
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,nazimu,npolar)
chiOptimizeAngularQuadratureForPolarSymmetry(pquad, 4.0*math.pi)

--############################################### Setup LBS parameters
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
vtk_basename = "first_example"
chiExportFieldFunctionToVTK(fflist[1],vtk_basename)

```
___
Back to [**Tutorial Home**](../tutorials_transport.md#first_example)