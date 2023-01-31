#ifndef CHI_PHYSICS_TRANSPORT_CROSS_SECTIONS_H
#define CHI_PHYSICS_TRANSPORT_CROSS_SECTIONS_H

#include "ChiPhysics/PhysicsMaterial/material_property_base.h"
#include "ChiMath/SparseMatrix/chi_math_sparse_matrix.h"

typedef std::vector<std::pair<double,double>> Tvecdbl_vecdbl;

/**\defgroup LuaTransportXSs Transport Cross Sections
 * \ingroup LuaPhysics*/

namespace chi_physics
{

//###################################################################
/** Class for handling Transport-Theory related cross sections.*/
class TransportCrossSections : public chi_physics::MaterialProperty
{
protected:
  using TransferMatrix = chi_math::SparseMatrix;
  using EmissionSpectra = std::vector<std::vector<double>>;

public:
  /// A struct containing data for a delayed neutron precursor.
  struct Precursor
  {
    double decay_constant = 0.0;
    double fractional_yield = 0.0;
    std::vector<double> emission_spectrum;
  };

public:
  unsigned int num_groups = 0;       ///< Total number of groups
  unsigned int scattering_order = 0; ///< Legendre scattering order
  unsigned int num_precursors = 0;   ///< Number of precursors

  bool is_fissionable = false;
  bool is_fission_scaled = false;

  /// Energy bin boundaries in MeV
  std::vector<std::vector<double>> e_bounds;

  std::vector<double> sigma_t;  ///< Total cross section
  std::vector<double> sigma_a;  ///< Absorption cross section
  std::vector<double> sigma_f;  ///< Fission cross section

  std::vector<double> nu_sigma_f;
  std::vector<double> nu_prompt_sigma_f;
  std::vector<double> nu_delayed_sigma_f;

  std::vector<double> inv_velocity;

  std::vector<TransferMatrix> transfer_matrices;
  std::vector<std::vector<double>> production_matrix;

  std::vector<Precursor> precursors;

  //Diffusion quantities
public:
  bool diffusion_initialized = false;
public:
  std::vector<double> diffusion_coeff; ///< Transport corrected diffusion coeff
  std::vector<double> sigma_removal;   ///< Removal cross section
  std::vector<double> sigma_s_gtog;    ///< Within-group scattering xs

  //Monte-Carlo quantities
public:
  bool scattering_initialized = false;
private:
  std::vector<std::vector<double>>         cdf_gprime_g;
  std::vector<std::vector<Tvecdbl_vecdbl>> scat_angles_gprime_g;

public:
  //00
  TransportCrossSections();

private:
  void Reset();

public:
  void MakeSimple0(int n_grps, double sigma);
  void MakeSimple1(int n_grps, double sigma, double c);
  void MakeCombined(std::vector<std::pair<int,double>>& combinations);

private:
  void ComputeAbsorption();

public:
  void ScaleFissionData(double k);

public:
  //01
  void MakeFromChiXSFile(const std::string &file_name);

  //02
  void ComputeDiffusionParameters();

  //05
  void PushLuaTable(lua_State* L) override;

  //06
  void ExportToChiXSFile(const std::string& file_name);


};

}//namespace chi_physics


#endif