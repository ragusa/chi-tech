#include "chi_lua.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"
#include "ChiMath/Quadratures/angular_quadrature_base.h"
#include "ChiMath/Quadratures/angular_product_quadrature.h"
#include "ChiMath/chi_math_range.h"

#include "ChiPhysics/FieldFunction/fieldfunction.h"
#include "ChiPhysics/PhysicsMaterial/transportxsections/material_property_transportxsections.h"

#include "ChiDataTypes/ndarray.h"

namespace chi_unit_sim_tests
{

/**PWLD Sweep. */
int chiSimTest91_PWLD(lua_State* L)
{
  const std::string fname = "chiSimTest91_PWLD";
  const int num_args = lua_gettop(L);
  chi::log.Log() << "chiSimTest91_PWLD num_args = " << num_args;

  if (chi::mpi.process_count != 1)
    throw std::logic_error(fname + ": Is serial only.");

  //============================================= Get grid
  auto grid_ptr = chi_mesh::GetCurrentHandler().GetGrid();
  const auto& grid = *grid_ptr;

  chi::log.Log() << "Global num cells: " << grid.GetGlobalNumberOfCells();

  //============================================= Make Orthogonal mapping
  const auto ijk_info = grid.GetIJKInfo();
  const auto& ijk_mapping = grid.MakeIJKToGlobalIDMapping();


  const auto Nx = static_cast<int64_t>(ijk_info[0]);
  const auto Ny = static_cast<int64_t>(ijk_info[1]);
  const auto Nz = static_cast<int64_t>(ijk_info[2]);

  const auto Dim1 = chi_mesh::DIMENSION_1;
  const auto Dim2 = chi_mesh::DIMENSION_2;
  const auto Dim3 = chi_mesh::DIMENSION_3;

  int    dimension = 0;
  if (grid.Attributes() & Dim1) dimension = 1;
  if (grid.Attributes() & Dim2) dimension = 2;
  if (grid.Attributes() & Dim3) dimension = 3;

  //============================================= Make SDM
  typedef std::shared_ptr<chi_math::SpatialDiscretization> SDMPtr;
  SDMPtr sdm_ptr = chi_math::SpatialDiscretization_PWLD::New(grid_ptr);
  const auto& sdm = *sdm_ptr;

  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;

  const size_t num_local_nodes = sdm.GetNumLocalDOFs(OneDofPerNode);
  const size_t num_globl_nodes = sdm.GetNumGlobalDOFs(OneDofPerNode);

  chi::log.Log() << "Num local nodes: " << num_local_nodes;
  chi::log.Log() << "Num globl nodes: " << num_globl_nodes;

  //============================================= Make an angular quadrature
  std::shared_ptr<chi_math::AngularQuadrature> quadrature;
  if (dimension == 1)
    quadrature = std::make_shared<chi_math::AngularQuadratureProdGL>(8);
  else if (dimension == 2)
  {
    quadrature = std::make_shared<chi_math::AngularQuadratureProdGLC>(8,8);
//    quadrature->OptimizeForPolarSymmetry(4.0*M_PI);
  }
  else if (dimension == 3)
    quadrature = std::make_shared<chi_math::AngularQuadratureProdGLC>(8,8);
  else
    throw std::logic_error(fname + "Error with the dimensionality "
                                   "of the mesh.");
  chi::log.Log() << "Quadrature created." << std::endl;

  //============================================= Set/Get params
  const size_t scat_order = 1;
  const size_t num_groups = 20;

  quadrature->BuildMomentToDiscreteOperator(scat_order,dimension);
  quadrature->BuildDiscreteToMomentOperator(scat_order,dimension);

  const auto& m2d = quadrature->GetMomentToDiscreteOperator();
  const auto& d2m = quadrature->GetDiscreteToMomentOperator();
  const auto& m_ell_em_map = quadrature->GetMomentToHarmonicsIndexMap();

  const size_t num_moments = m_ell_em_map.size();
  const size_t num_dirs = quadrature->omegas.size();

  chi::log.Log() << "End Set/Get params." << std::endl;
  chi::log.Log() << "Num Moments: " << num_moments << std::endl;

  //============================================= Make Unknown Managers
  const auto VecN = chi_math::UnknownType::VECTOR_N;
  using Unknown = chi_math::Unknown;

  std::vector<Unknown> phi_uks(num_moments, Unknown(VecN, num_groups));
  std::vector<Unknown> psi_uks(num_dirs,    Unknown(VecN, num_groups));

  const chi_math::UnknownManager phi_uk_man(phi_uks);
  const chi_math::UnknownManager psi_uk_man(psi_uks);

  const size_t num_local_phi_dofs = sdm.GetNumLocalDOFs(phi_uk_man);
  const size_t num_local_psi_dofs = sdm.GetNumLocalDOFs(psi_uk_man);

  chi::log.Log() << "End ukmanagers." << std::endl;

  //============================================= Make XSs
  chi_physics::TransportCrossSections xs;
  xs.MakeFromChiXSFile("tests/xs_graphite_pure.cxs");

  //============================================= Initializes vectors
  std::vector<double> phi_old(num_local_phi_dofs,0.0);
  std::vector<double> psi_old(num_local_psi_dofs, 0.0);
  auto source_moments = phi_old;
  auto phi_new        = phi_old;
  auto q_source       = phi_old;

  chi::log.Log() << "End vectors." << std::endl;

  //============================================= Make material source term
  for (const auto& cell : grid.local_cells)
  {
    const auto& cc = cell.centroid;
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    if (cc.x < 0.5 and cc.y < 0.5 and cc.z < 0.5 and
        cc.x >-0.5 and cc.y >-0.5 and cc.z >-0.5)
    {
      for (size_t i=0; i<num_nodes; ++i)
      {
        const int64_t dof_map = sdm.MapDOFLocal(cell,i,phi_uk_man,0,0);

        q_source[dof_map] = 1.0;
      }
    }
  }

  //============================================= Precompute cell matrices
  typedef chi_mesh::Vector3 Vec3;
  typedef std::vector<Vec3> VecVec3;
  typedef std::vector<VecVec3> MatVec3;
  typedef std::vector<double> VecDbl;
  typedef std::vector<VecDbl> MatDbl;
  typedef std::vector<MatDbl> VecMatDbl;

  std::vector<MatVec3>   cell_Gmatrices;
  std::vector<MatDbl>    cell_Mmatrices;
  std::vector<VecMatDbl> cell_faceMmatrices;

  for (const auto& cell : grid.local_cells)
  {
    const auto&  cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes    = cell_mapping.NumNodes();
    const auto   vol_qp_data  = cell_mapping.MakeVolumeQuadraturePointData();

    MatVec3 IntV_shapeI_gradshapeJ(num_nodes, VecVec3(num_nodes,Vec3(0,0,0)));
    MatDbl  IntV_shapeI_shapeJ(num_nodes, VecDbl(num_nodes,0.0));

    for (unsigned int i = 0; i < num_nodes; ++i)
      for (unsigned int j = 0; j < num_nodes; ++j)
        for (const auto& qp : vol_qp_data.QuadraturePointIndices())
        {
          IntV_shapeI_gradshapeJ[i][j]
            += vol_qp_data.ShapeValue(i, qp) *
               vol_qp_data.ShapeGrad(j, qp) *
               vol_qp_data.JxW(qp);

          IntV_shapeI_shapeJ[i][j]
            += vol_qp_data.ShapeValue(i, qp) *
               vol_qp_data.ShapeValue(j, qp) *
               vol_qp_data.JxW(qp);
        }// for qp

    cell_Gmatrices.push_back(std::move(IntV_shapeI_gradshapeJ));
    cell_Mmatrices.push_back(std::move(IntV_shapeI_shapeJ));

    const size_t num_faces = cell.faces.size();
    VecMatDbl faces_Mmatrices;
    for (size_t f=0; f<num_faces; ++f)
    {
      const auto face_qp_data = cell_mapping.MakeFaceQuadraturePointData(f);
      MatDbl IntS_shapeI_shapeJ(num_nodes,VecDbl(num_nodes,0.0));
      for (unsigned int i = 0; i < num_nodes; ++i)
        for (unsigned int j = 0; j < num_nodes; ++j)
          for (const auto& qp : face_qp_data.QuadraturePointIndices())
            IntS_shapeI_shapeJ[i][j]
              += face_qp_data.ShapeValue(i, qp) *
                 face_qp_data.ShapeValue(j, qp) *
                 face_qp_data.JxW(qp);

      faces_Mmatrices.push_back(std::move(IntS_shapeI_shapeJ));
    }//for face f

    cell_faceMmatrices.push_back(std::move(faces_Mmatrices));

  }//for cell

  chi::log.Log() << "End cell matrices." << std::endl;

  //============================================= Make Grid internal face mapping
  const auto cell_adj_mapping = sdm.MakeInternalFaceNodeMappings();

  //============================================= Define sweep chunk
  auto SweepChunk = [&ijk_mapping, &grid, &sdm,
                     &num_moments,
                     &phi_uk_man, &psi_uk_man,
                     &m2d,&d2m,
                     &phi_new, &source_moments, &psi_old,
                     &cell_Gmatrices, &cell_Mmatrices, &cell_faceMmatrices,
                     &cell_adj_mapping]
    (const std::array<int64_t,3>& ijk,
     const Vec3& omega,
     const size_t d,
     const chi_physics::TransportCrossSections& cell_xs)
  {
    const auto cell_global_id = ijk_mapping.MapNDtoLin(ijk[1],ijk[0],ijk[2]);
    const auto& cell = grid.cells[cell_global_id];
    const auto cell_local_id = cell.local_id;
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const size_t num_faces = cell.faces.size();

    const std::vector<double> zero_vector(num_groups,0.0);

    const auto& G = cell_Gmatrices[cell_local_id];
    const auto& M = cell_Mmatrices[cell_local_id];

    MatDbl A(num_nodes, VecDbl(num_nodes, 0.0));
    MatDbl b(num_groups, VecDbl(num_nodes, 0.0));

    //================================= Gradient matrix
    for (size_t i = 0; i < num_nodes; ++i)
      for (size_t j = 0; j < num_nodes; ++j)
        A[i][j] = omega.Dot(G[i][j]);

    //================================= Surface integrals
    for (size_t f=0; f<num_faces; ++f)
    {
      const auto& face = cell.faces[f];
      const double mu = omega.Dot(face.normal);

      if (mu < 0.0)
      {
        const auto& M_surf = cell_faceMmatrices[cell_local_id][f];

        const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);
        for (size_t fi=0; fi<num_face_nodes; ++fi)
        {
          const int i = cell_mapping.MapFaceNode(f,fi);
          for (size_t fj=0; fj<num_face_nodes; ++fj)
          {
            const int j = cell_mapping.MapFaceNode(f,fj);

            const double* upwind_psi = zero_vector.data();
            if (face.has_neighbor)
            {
              const auto& adj_cell = grid.cells[face.neighbor_id];
              const int aj = cell_adj_mapping[cell.local_id][f][fj];
              const int64_t ajmap = sdm.MapDOFLocal(adj_cell,aj,psi_uk_man,d,0);
              upwind_psi = &psi_old[ajmap];
            }

            const double mu_Nij = -mu * M_surf[i][j];
            A[i][j] += mu_Nij;
            for (int g=0; g<num_groups; ++g)
              b[g][i] += upwind_psi[g]*mu_Nij;
          }//for fj
        }//for fi
      }//if internal incident face
    }//for face

    for (size_t g=0; g<num_groups; ++g)
    {
      auto Atemp = A;
      VecDbl source(num_nodes, 0.0);
      //Nodal source moments
      for (size_t i=0; i<num_nodes; ++i)
      {
        double temp_src = 0.0;
        for (size_t m=0; m<num_moments; ++m)
        {
          const int64_t dof_map = sdm.MapDOFLocal(cell,i,phi_uk_man,m,g);
          temp_src += m2d[m][d]*source_moments[dof_map];
        }//for m
        source[i] = temp_src;
      }//for i


      //Mass Matrix and Source
      const double sigma_tg = cell_xs.sigma_t[g];

      for (int i = 0; i < num_nodes; ++i)
      {
        double temp = 0.0;
        for (int j = 0; j < num_nodes; ++j)
        {
          const double Mij = M[i][j];
          Atemp[i][j] = A[i][j] + Mij*sigma_tg;
          temp += Mij*source[j];
        }//for j
        b[g][i] += temp;
      }//for i

      // ============================= Solve system
      chi_math::GaussElimination(Atemp, b[g], static_cast<int>(num_nodes));
    }//for g

    //Accumulate flux-moments
    for (size_t m=0; m<num_moments; ++m)
    {
      const double wn_d2m = d2m[m][d];
      for (size_t i=0; i<num_nodes; ++i)
      {
        const int64_t dof_map = sdm.MapDOFLocal(cell,i,phi_uk_man,m,0);
        for (size_t g=0; g<num_groups; ++g)
          phi_new[dof_map + g] += wn_d2m * b[g][i];
      }
    }

    //Save angular fluxes
    for (size_t i=0; i<num_nodes; ++i)
    {
      const int64_t dof_map = sdm.MapDOFLocal(cell,i,psi_uk_man,d,0);
      for (size_t g=0; g<num_groups; ++g)
        psi_old[dof_map + g] = b[g][i];
    }
  };


  //============================================= Define sweep for all dirs
  auto Sweep = [&num_dirs,&quadrature,Nx,Ny,Nz,&SweepChunk,&xs]()
  {
    for (size_t d=0; d<num_dirs; ++d)
    {
      const auto &omega = quadrature->omegas[d];
      const auto &weight = quadrature->weights[d];

      std::vector<int64_t> iorder, jorder, korder;
      if (omega.x > 0.0) iorder = chi_math::Range<int64_t>(0, Nx);
      else               iorder = chi_math::Range<int64_t>(Nx - 1, -1, -1);
      if (omega.y > 0.0) jorder = chi_math::Range<int64_t>(0, Ny);
      else               jorder = chi_math::Range<int64_t>(Ny - 1, -1, -1);
      if (omega.z > 0.0) korder = chi_math::Range<int64_t>(0, Nz);
      else               korder = chi_math::Range<int64_t>(Nz - 1, -1, -1);

      for (auto i: iorder)
        for (auto j: jorder)
          for (auto k: korder)
            SweepChunk({i,j,k}, omega, d, xs);
    }//for d
  };

  //============================================= Define SetSource routine
  auto SetSource = [&source_moments, &phi_old, &q_source, &grid, &sdm,
                    &m_ell_em_map, &xs, num_moments,
                    &phi_uk_man]()
  {
    for (const auto& cell : grid.local_cells)
    {
      const auto& cell_mapping = sdm.GetCellMapping(cell);
      const size_t num_nodes = cell_mapping.NumNodes();
      const auto& S = xs.transfer_matrices;

      for (size_t i=0; i<num_nodes; ++i)
      {
        for (size_t m=0; m<num_moments; ++m)
        {
          const int64_t dof_map = sdm.MapDOFLocal(cell,i,phi_uk_man,m,0);
          const auto ell = m_ell_em_map[m].ell;

          for (size_t g=0; g<num_groups; ++g)
          {
            //Fixed source
            source_moments[dof_map + g] = q_source[dof_map + g];

            //Inscattering
            if (ell < S.size())
            {
              double inscat_g = 0.0;
              for (const auto& [row_g, gprime, sigma_sm] : S[ell].Row(g))
                inscat_g += sigma_sm * phi_old[dof_map + gprime];

              source_moments[dof_map + g] += inscat_g;
            }
          }//for g
        }//for m
      }//for node i
    }//for cell
  };

  //============================================= Define L-infinite-norm
  auto ComputeRelativePWChange = [&grid,&sdm,
                                  &num_moments,
                                  &phi_uk_man]
    (const std::vector<double>& in_phi_new,
     const std::vector<double>& in_phi_old)
  {
    double pw_change = 0.0;

    for (const auto& cell : grid.local_cells)
    {
      const auto& cell_mapping = sdm.GetCellMapping(cell);
      const size_t num_nodes = cell_mapping.NumNodes();

      for (size_t i=0; i<num_nodes; ++i)
      {
        //Get scalar moments
        const int64_t m0_map = sdm.MapDOFLocal(cell,i,phi_uk_man,0,0);

        const double* phi_new_m0 = &in_phi_new[m0_map];
        const double* phi_old_m0 = &in_phi_old[m0_map];
        for (size_t m=0; m<num_moments; ++m)
        {
          const int64_t m_map = sdm.MapDOFLocal(cell,i,phi_uk_man,m,0);

          const double* phi_new_m = &in_phi_new[m_map];
          const double* phi_old_m = &in_phi_old[m_map];

          for (size_t g=0; g<num_groups; ++g)
          {
            const double abs_phi_new_g_m0 = std::fabs(phi_new_m0[g]);
            const double abs_phi_old_g_m0 = std::fabs(phi_old_m0[g]);

            const double max_denominator = std::max(abs_phi_new_g_m0,
                                                    abs_phi_old_g_m0);

            const double delta_phi = std::fabs(phi_new_m[g] - phi_old_m[g]);

            if (max_denominator >= std::numeric_limits<double>::min())
              pw_change = std::max(delta_phi/max_denominator,pw_change);
            else
              pw_change = std::max(delta_phi,pw_change);
          }//for g
        }//for m
      }//for i
    }//for cell

    return pw_change;
  };

  //============================================= Classic Richardson iteration
  chi::log.Log() << "Starting iterations" << std::endl;
  for (size_t iter=0; iter<200; ++iter)
  {
    phi_new.assign(phi_new.size(), 0.0);
    //Build rhs
    SetSource();
    Sweep();

    const double rel_change = ComputeRelativePWChange(phi_new, phi_old);

    std::stringstream outstr;
    outstr << "Iteration " << std::setw(5) << iter << " ";
    {
      char buffer[100];
      snprintf(buffer,100, "%11.3e\n", rel_change);
      outstr << buffer;
    }

    chi::log.Log() << outstr.str();

    phi_old = phi_new;

    if (rel_change < 1.0e-6 and iter > 0)
      break;
  }//for iteration

  //============================================= Create Field Functions
  std::vector<std::shared_ptr<chi_physics::FieldFunction>> ff_list;

  ff_list.push_back(std::make_shared<chi_physics::FieldFunction>(
    "Phi",                                           //Text name
    sdm_ptr,                                         //Spatial Discr.
    chi_math::Unknown(chi_math::UnknownType::VECTOR_N,num_groups) //Unknown
  ));

  const std::vector<std::string> dim_strings = {"x","y","z"};
  for (const std::string& dim : dim_strings)
    ff_list.push_back(std::make_shared<chi_physics::FieldFunction>(
      "J-"+dim,                                        //Text name
      sdm_ptr,                                         //Spatial Discr.
      chi_math::Unknown(chi_math::UnknownType::VECTOR_N,num_groups) //Unknown
    ));

  //============================================= Localize zeroth moment
  //This routine extracts a single moment vector
  //from the vector that contains multiple moments
  const chi_math::UnknownManager m0_uk_man(
    {chi_math::Unknown(chi_math::UnknownType::VECTOR_N,num_groups)});
  const size_t num_m0_dofs = sdm.GetNumLocalDOFs(m0_uk_man);

  std::vector<double> m0_phi(num_m0_dofs, 0.0);
  std::vector<double> mx_phi(num_m0_dofs, 0.0); //Y(1,1)  - X-component
  std::vector<double> my_phi(num_m0_dofs, 0.0); //Y(1,-1) - Y-component
  std::vector<double> mz_phi(num_m0_dofs, 0.0); //Y(1,0)  - Z-component

  sdm.CopyVectorWithUnknownScope(phi_old,     //from vector
                                 m0_phi,      //to vector
                                 phi_uk_man,  //from dof-structure
                                 0,           //from unknown-id
                                 m0_uk_man,   //to dof-structure
                                 0);          //to unknown-id

  ff_list[0]->UpdateFieldVector(m0_phi);

  std::array<unsigned int,3> j_map = {0,0,0};
  if (dimension == 1 and num_moments >= 2) j_map = {0,0,1};
  if (dimension == 2 and num_moments >= 3) j_map = {2,1,0};
  if (dimension == 3 and num_moments >= 4) j_map = {3,1,2};

  sdm.CopyVectorWithUnknownScope(
    phi_old, mx_phi, phi_uk_man, j_map[0], m0_uk_man, 0);
  sdm.CopyVectorWithUnknownScope(
    phi_old, my_phi, phi_uk_man, j_map[1], m0_uk_man, 0);
  sdm.CopyVectorWithUnknownScope(
    phi_old, mz_phi, phi_uk_man, j_map[2], m0_uk_man, 0);

  ff_list[1]->UpdateFieldVector(mx_phi);
  ff_list[2]->UpdateFieldVector(my_phi);
  ff_list[3]->UpdateFieldVector(mz_phi);


  //============================================= Update field function
  chi_physics::FieldFunction::FFList const_ff_list;
  for (const auto& ff_ptr : ff_list)
    const_ff_list.push_back(ff_ptr);
  chi_physics::FieldFunction::ExportMultipleToVTK("SimTest_91a_PWLD",
                                                  const_ff_list);

  return 0;
}

}//namespace chi_unit_tests