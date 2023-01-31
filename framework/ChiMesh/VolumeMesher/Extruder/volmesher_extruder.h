#ifndef VOLUME_MESHER_EXTRUDER_H
#define VOLUME_MESHER_EXTRUDER_H

#include <utility>

#include "ChiMesh/VolumeMesher/chi_volumemesher.h"
#include "ChiMesh/Cell/cell.h"



//###################################################################
/**An extruder mesher taking a flat surface and extruding it.*/
class chi_mesh::VolumeMesherExtruder : public chi_mesh::VolumeMesher
{
public:
  enum class TemplateType : int
  {
    SURFACE_MESH       = 1,
    UNPARTITIONED_MESH = 2
  };
  struct MeshLayer
  {
    std::string name;
    double height;
    int    sub_divisions;
  };
private:
  const TemplateType template_type;
  std::shared_ptr<SurfaceMesh> template_surface_mesh = nullptr;
  std::shared_ptr<UnpartitionedMesh> template_unpartitioned_mesh = nullptr;
public:
  std::vector<MeshLayer> input_layers;
  std::vector<double> vertex_layers;
  size_t node_z_index_incr=0;

public:
  explicit
  VolumeMesherExtruder(std::shared_ptr<chi_mesh::SurfaceMesh> in_surface_mesh) :
    VolumeMesher(VolumeMesherType::EXTRUDER),
    template_type(TemplateType::SURFACE_MESH),
    template_surface_mesh(std::move(in_surface_mesh))
  {}
  explicit
  VolumeMesherExtruder(std::shared_ptr<chi_mesh::UnpartitionedMesh> in_unpartitioned_mesh) :
    VolumeMesher(VolumeMesherType::EXTRUDER),
    template_type(TemplateType::UNPARTITIONED_MESH),
    template_unpartitioned_mesh(std::move(in_unpartitioned_mesh))
  {}

  void Execute() override;

private:
  //utils
  void SetupLayers(int default_layer_count=1);

  chi_mesh::Vector3 ProjectCentroidToLevel(const chi_mesh::Vector3& centroid,
                                           size_t level);
  int GetCellKBAPartitionIDFromCentroid(chi_mesh::Vector3& centroid);


  bool HasLocalScope(
    const chi_mesh::Cell& template_cell,
    const chi_mesh::MeshContinuum& template_continuum,
    size_t z_level);

  std::unique_ptr<chi_mesh::Cell>
  MakeExtrudedCell(const chi_mesh::Cell& template_cell,
                   const chi_mesh::MeshContinuum& grid,
                   size_t z_level,
                   uint64_t cell_global_id,
                   int partition_id,
                   size_t num_template_cells);

  void CreateLocalNodes(chi_mesh::MeshContinuum& template_grid,
                        chi_mesh::MeshContinuum& grid);

  void ExtrudeCells(chi_mesh::MeshContinuum& template_grid,
                    chi_mesh::MeshContinuum& grid);

};



#endif //VOLUME_MESHER_EXTRUDER_H