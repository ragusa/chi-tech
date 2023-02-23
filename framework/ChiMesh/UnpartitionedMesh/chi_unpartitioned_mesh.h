#ifndef CHI_MESH_UNPARTITIONED_MESH_H
#define CHI_MESH_UNPARTITIONED_MESH_H

#include "ChiMesh/chi_mesh.h"
#include "ChiMesh/Cell/cell.h"

class vtkCell;
class vtkUnstructuredGrid;
template<class T>
class vtkSmartPointer;

//###################################################################
/**This object is intented for unpartitioned meshes that still require
 * partitioning.*/
class chi_mesh::UnpartitionedMesh
{
public:
  struct LightWeightFace
  {
    std::vector<uint64_t> vertex_ids;
    bool has_neighbor = false;
    uint64_t neighbor=0;

    LightWeightFace() = default;
    explicit
    LightWeightFace(std::vector<uint64_t> in_vertex_ids) :
      vertex_ids(std::move(in_vertex_ids)) {}
  };
  struct LightWeightCell
  {
    const chi_mesh::CellType type;
    const chi_mesh::CellType sub_type;
    chi_mesh::Vertex centroid;
    int material_id=-1;
    std::vector<uint64_t> vertex_ids;
    std::vector<LightWeightFace> faces;

    explicit
    LightWeightCell(chi_mesh::CellType in_type,
                    chi_mesh::CellType in_sub_type) :
                    type(in_type),
                    sub_type(in_sub_type) {}
  };

public:
  std::vector<chi_mesh::Vertex>    vertices;
  std::vector<LightWeightCell*>    raw_cells;
  std::vector<LightWeightCell*>    raw_boundary_cells;
  std::vector<std::set<uint64_t>>  vertex_cell_subscriptions;

  MeshAttributes attributes = NONE;

public:
  struct Options
  {
    std::string file_name;
    std::string material_id_fieldname;
    std::string boundary_id_fieldname;
    double scale=1.0;
    size_t ortho_Nx = 0;
    size_t ortho_Ny = 0;
    size_t ortho_Nz = 0;
  }mesh_options;

  struct BoundBox
  {
    double xmin=0.0, xmax=0.0,
           ymin=0.0, ymax=0.0,
           zmin=0.0, zmax=0.0;
  } bound_box;

protected:
  static LightWeightCell* CreateCellFromVTKPolyhedron(vtkCell* vtk_cell);
  static LightWeightCell* CreateCellFromVTKPolygon(vtkCell* vtk_cell);
  static LightWeightCell* CreateCellFromVTKLine(vtkCell* vtk_cell);
  static LightWeightCell* CreateCellFromVTKVertex(vtkCell* vtk_cell);

  typedef vtkSmartPointer<vtkUnstructuredGrid> vtkUGridPtr;
  static int FindHighestDimension(std::vector<vtkUGridPtr>& ugrid_blocks);
  static vtkUGridPtr
    ConsolidateAndCleanBlocks(std::vector<vtkUGridPtr>& ugrid_blocks,
                              int desired_dimension);
  void CopyUGridCellsAndPoints(vtkUnstructuredGrid& ugrid,
                               double scale);


  static std::vector<uint64_t>
    BuildBlockCellExtents(std::vector<vtkUGridPtr>& ugrid_blocks,
                          int desired_dimension);
  void SetMaterialIDsFromBlocks(const std::vector<uint64_t>& block_mat_ids);

  std::vector<int>
    BuildCellMaterialIDsFromField(vtkUGridPtr &ugrid,
                                  const std::string& field_name,
                                  const std::string& file_name) const;
  void SetMaterialIDsFromList(const std::vector<int>& material_ids);


public:
  void BuildMeshConnectivity();
  void ComputeCentroidsAndCheckQuality();

  void ReadFromVTU(const Options& options);
  void ReadFromEnsightGold(const Options& options);
  void ReadFromWavefrontOBJ(const Options& options);

  void ReadFromMsh(const Options& options);

  void ReadFromExodus(const Options& options);

  void PushProxyCell(const std::string& type_str,
                     const std::string& sub_type_str,
                     int cell_num_faces,
                     int cell_material_id,
                     const std::vector<std::vector<uint64_t>>& proxy_faces);

  ~UnpartitionedMesh()
  {
    for (auto& cell : raw_cells)          delete cell;
    for (auto& cell : raw_boundary_cells) delete cell;
  }
  void CleanUp()
  {
    for (auto& cell : raw_cells)          delete cell;
    for (auto& cell : raw_boundary_cells) delete cell;
    vertices.clear();
    raw_cells.clear();
    raw_boundary_cells.clear();
    vertex_cell_subscriptions.clear();
  }
};


#endif //CHI_MESH_UNPARTITIONED_MESH_H
