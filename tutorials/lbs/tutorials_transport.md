# Tutorials for the Linear Boltzmann Solver (LBS)
```lbs``` is the main root for the LBS tutorials

--------------------------------------
## <a name="first_example"></a> **A first example for an LBS input**
- [**First Example with Orthogonal Grids and KBA Partitioning**](./first/first_example.md)

--------------------------------------
## <a name="Meshing"></a> **Meshing**
*Note*:  Examples, denoted by ```#```, use KBA partitioning.
- [**A Simple 1D-slab Grid: Reed's Problem**](./meshing/reed.md)
- 2D Unstructured Grids:
  - [```#```**Creating and Reading a 2D .gmsh file with material IDs and Boundary IDs**](./meshing/read_2D_gmsh.md)
  - [```#```**Creating and Reading a 2D .obj file with material IDs and Boundary IDs**](./meshing/read_2D_obj_mesh.md)
  - [```#```**Creating and Reading a 2D spiderweb grid file with material IDs and Boundary IDs**](./meshing/2D_spider.md)
- 3D Extrusion of 2D Unstructured Grids:
  - [```#```**3D Extrusion of 2D Unstructured Grid with material IDs and Boundary IDs**](./meshing/3D_extrusion.md)
- 3D Unstructured Grids:
  - [```#```**Reading a 3D Grid with material IDs and Boundary IDs**](./meshing/3D_mesh.md)
- Logical Volumes:
  - [```#```**Logical Volumes in 2D Orthogonal Grids**](./meshing/logical_volumes_2Dortho.md)
  - [```#```**Logical Volumes in 3D Orthogonal Grids**](./meshing/logical_volumes_3Dortho.md)
  - [```#```**Logical Volumes in 2D Extruded Unstructured Grids**](./meshing/logical_volumes_2D_extruded.md)
- Using SplitMesh:
  - [```#```**Using SplitMesh**](./meshing/splitmesh.md)

--------------------------------------
## <a name="Angular"></a> **Angular Quadrature**
- [**Gauss-Legendre-Chebyshev Product Quadrature**](./angular_quad/productGLC.md)
- [**Simplified LDFE Quadrature**](./angular_quad/ldfe_quadrature.md)
- [**Quadrature for Cylindrical Geometry**](./angular_quad/cyl_quad.md)

--------------------------------------
## <a name="Materials"></a> **Materials**
i.e., **Cross Sections and Volumetric Sources**
- [**Multigroup Cross-Section Sets and Sources**](./materials/mg_xs.md)
- [**Multi-materials**](./materials/multi_materials.md)
- [**Anisotropic Scattering**](./materials/anisotropic.md)

--------------------------------------
## <a name="Boundary"></a> **Boundary Conditions and Boundary Sources**
- [**Reflective BC**](./boundary/reflective_bc.md)
- [**Boundary Source**](./boundary/bd_src.md)

--------------------------------------
## <a name="Steady_state"></a> **Steady-state Source-driven Runs in XYZ Geometries**
- [**1g Run with Source Iteration**](./steady_state/1g_SI.md)
- [**1g Run with Unpreconditioned GMRes**](./steady_state/1g_GMRes.md)
- [**1g Run with Source Iteration + DSA**](./steady_state/1_g_SI_DSA.md)
- [**1g Run with Preconditioned GMRes**](./steady_state/1g_Precond_GMRes.md)

- [**Multigroup Run with GroupSets**](./steady_state/mg_groupsets.md)
- [**Two-grid acceleration in Graphite**](./steady_state/mg_two_grid.md)

--------------------------------------
## <a name="Keigenvalue"></a> **K-eigenvalue Runs in XYZ Geometries**

--------------------------------------
## <a name="Adjoint"></a> **Adjoint Runs and Post-processors**

--------------------------------------
## <a name="Steady_state_cyl"></a> **Steady-state Source-driven Runs in RZ Cylindrical Geometries**

--------------------------------------
## <a name="Transient"></a> **Transient Runs**
