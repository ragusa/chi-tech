#ifndef SPATIAL_DISCRETIZATION_PWL_BASE_H
#define SPATIAL_DISCRETIZATION_PWL_BASE_H

#include "ChiMath/SpatialDiscretization/FiniteElement/spatial_discretization_FE.h"

#include "ChiMath/Quadratures/Line/quadrature_line.h"
#include "ChiMath/Quadratures/Triangle/quadrature_triangle.h"
#include "ChiMath/Quadratures/Quadrilateral/quadrature_quadrilateral.h"
#include "ChiMath/Quadratures/Tetrahedron/quadrature_tetrahedron.h"

namespace chi_math
{

  class SpatialDiscretization_PWLBase : public chi_math::SpatialDiscretization_FE
  {
  protected:
    QuadratureLine          line_quad_order_arbitrary;
    QuadratureTriangle      tri_quad_order_arbitrary;
    QuadratureQuadrilateral quad_quad_order_arbitrary;
    QuadratureTetrahedron   tet_quad_order_arbitrary;

  protected:
    explicit
    SpatialDiscretization_PWLBase(chi_mesh::MeshContinuumPtr& in_grid,
                                  finite_element::SetupFlags in_setup_flags,
                                  QuadratureOrder qorder,
                                  SDMType in_type,
                                  CoordinateSystemType in_cs_type
                                  );
    //01
  protected:
    void PreComputeCellSDValues();
    void PreComputeNeighborCellSDValues();

    void CreateCellMappings();

    //02
    //Child specialized:
    //OrderNodes

    //03
    //Child specialized:
    //BuildSparsityPattern

    //04
    //Child specialized:
    //MapDOF
    //MapDOFLocal

    //05
  public:
    size_t GetNumLocalDOFs(const UnknownManager& unknown_manager) const override;
    size_t GetNumGlobalDOFs(const UnknownManager& unknown_manager) const override;

    //Child Specialized:
    //GetNumGhostDOFs
    //GetGhostDOFIndices

    size_t GetCellNumNodes(const chi_mesh::Cell& cell) const override;

    std::vector<chi_mesh::Vector3>
    GetCellNodeLocations(const chi_mesh::Cell& cell) const override;

    void LocalizePETScVector(Vec petsc_vector,
                             std::vector<double>& local_vector,
                             const UnknownManager& unknown_manager)
    const override;

    //FE-utils
    const finite_element::UnitIntegralData&
    GetUnitIntegrals(const chi_mesh::Cell& cell) override;

    const finite_element::InternalQuadraturePointData&
    GetQPData_Volumetric(const chi_mesh::Cell& cell) override;

    const finite_element::FaceQuadraturePointData&
    GetQPData_Surface(const chi_mesh::Cell& cell,
                      unsigned int face_index) override;
  };

}//namespace chi_math

#endif //SPATIAL_DISCRETIZATION_PWL_BASE_H
