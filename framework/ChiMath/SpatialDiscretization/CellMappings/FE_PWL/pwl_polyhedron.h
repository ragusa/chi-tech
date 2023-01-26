#ifndef PWL_POLYHEDRON_VALUES_H
#define PWL_POLYHEDRON_VALUES_H

#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_cellbase.h"
#include "ChiMath/Quadratures/Tetrahedron/quadrature_tetrahedron.h"
#include "ChiMath/Quadratures/Triangle/quadrature_triangle.h"
#include "ChiMesh/Cell/cell.h"
#include "pwl_cellbase.h"



//###################################################################
namespace chi_math
{
  /**Object for handling piecewise linear
     * shape functions on polyhedron shaped 3D cells.
     *
     * This object has a whitepaper associated with it
     * (<a target="_blank"
     * href="../../whitepages/FEM/PWLPolyhedron/PWLPolyhedron.pdf">
     * here</a>)
     *
     * Some notes on indexing:\n
     *  - IntV_shapeI_gradshapeJ, given i and j results in a vector.
     *  - IntS_shapeI_shapeJ, requires f, then cell i, then cell j
     *  - face_dof_mappings, is as follows face_dof_mappings[f][fi] and
     *    returns cell i.
     * */
  class PolyhedronMappingFE_PWL : public chi_math::CellMappingFE_PWL
  {
  private:
    /**Stores the data for each side's tetrahedron. */
    struct FEside_data3d
    {
      double                    detJ = 0.0;
      double                    detJ_surf = 0.0;
      std::vector<uint64_t>     v_index;
      chi_mesh::Vector3         v0;
      chi_mesh::Matrix3x3       J;
      chi_mesh::Matrix3x3       Jinv;
      chi_mesh::Matrix3x3       JTinv;
    };
    //Goes into
    /**Stores data for each face.*/
    struct FEface_data
    {
      std::vector<FEside_data3d> sides;
      chi_mesh::Vector3 normal;
    };
    //Goes int face_data


    /**Lowest level of mapping dof i.*/
    struct FEnodeSideMap
    {
      int index = -1;
      bool part_of_face = false;
    };
    //Goes into
    /**Intermediate level of mapping.*/
    struct FEnodeFaceMap
    {
      std::vector<FEnodeSideMap> side_map;
    };
    //Goes into
    /**Node map per face.*/
    struct FEnodeMap
    {
      std::vector<FEnodeFaceMap> face_map;
    };
    //Goes into node_maps
    // node n
    // face f
    // side s
    // node_maps[n]->face_map[f]->side_map[s]

  private:
    std::vector<double>            face_betaf;     ///< Face Beta-factor.
    double                         alphac;         ///< Cell alpha-factor.

    std::vector<FEface_data>       face_data;      ///< Holds determinants and data tet-by-tet.
  private:
    std::vector<FEnodeMap>         node_side_maps; ///< Maps nodes to side tets.

  private:
    const QuadratureTetrahedron& volume_quadrature;
    const QuadratureTriangle&    surface_quadrature;

  public:
    //00_constrdestr.cc
    PolyhedronMappingFE_PWL(const chi_mesh::Cell& polyh_cell,
                            const chi_mesh::MeshContinuumConstPtr& ref_grid,
                            const QuadratureTetrahedron& volume_quadrature,
                            const QuadratureTriangle&    surface_quadrature);

    void InitializeVolumeQuadraturePointData(
      finite_element::InternalQuadraturePointData& internal_data) const override;

    void InitializeFaceQuadraturePointData(
      unsigned int face,
      finite_element::FaceQuadraturePointData& faces_qp_data) const override;

    //################################################## Define standard
    //                                                   tetrahedron linear shape
    //                                                   functions
    //01a_reftet.cc
  private:
    static
    double TetShape(unsigned int index,
                    const chi_mesh::Vector3& qpoint,
                    bool on_surface = false);
    static double TetGradShape_x(unsigned int index);
    static double TetGradShape_y(unsigned int index);
    static double TetGradShape_z(unsigned int index);

    //################################################## Shape functions per face-side
    //01b_sidevalues.cc
  private:
    double FaceSideShape(unsigned int face_index,
                         unsigned int side_index,
                         unsigned int i,
                         const chi_mesh::Vector3& qpoint,
                         bool on_surface = false) const;
    double FaceSideGradShape_x(unsigned int face_index,
                               unsigned int side_index,
                               unsigned int i) const;
    double FaceSideGradShape_y(unsigned int face_index,
                               unsigned int side_index,
                               unsigned int i) const;
    double FaceSideGradShape_z(unsigned int face_index,
                               unsigned int side_index,
                               unsigned int i) const;

    //############################################### Actual shape functions
    //                                                as function of cartesian
    //                                                coordinates
  public:
    double ShapeValue(int i, const chi_mesh::Vector3& xyz) const override;
    chi_mesh::Vector3 GradShapeValue(
      int i,
      const chi_mesh::Vector3& xyz) const override;


    void ShapeValues(const chi_mesh::Vector3& xyz,
                     std::vector<double>& shape_values) const override;

    void GradShapeValues(
      const chi_mesh::Vector3& xyz,
      std::vector<chi_mesh::Vector3>& gradshape_values) const override;

  };
}

#endif
