#include "angular_quadrature_base.h"

#include "ChiMath/Quadratures/LegendrePoly/legendrepoly.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <iomanip>
#include <numeric>

//###################################################################
/**Optimizes the angular quadrature for polar symmetry by removing
 * all the direction with downward pointing polar angles.
 *
 * \param normalization float. (Optional) The default is a negative number
 *                             which does not apply any normalization. If a
 *                             positive number is provided, the weights will be
 *                             normalized to sum to this number.*/
//void chi_math::AngularQuadrature::
//  OptimizeForPolarSymmetry(const double normalization)
//{
//  std::vector<chi_mesh::Vector3>                 new_omegas;
//  std::vector<chi_math::QuadraturePointPhiTheta> new_abscissae;
//  std::vector<double>                            new_weights;
//
//  const size_t num_dirs = omegas.size();
//  double weight_sum = 0.0;
//  for (size_t d=0; d<num_dirs; ++d)
//    if (omegas[d].z > 0.0)
//    {
//      new_abscissae.emplace_back(abscissae[d]);
//      new_weights.emplace_back(weights[d]);
//      new_omegas.emplace_back(omegas[d]);
//      weight_sum += weights[d];
//    }
//
//  if (normalization > 0.0)
//    for (double& w : new_weights)
//      w *= normalization/weight_sum;
//
//  abscissae = std::move(new_abscissae);
//  weights   = std::move(new_weights);
//  omegas    = std::move(new_omegas);
//}

//###################################################################
/**Populates a map of moment m to the angular_Spherical Harmonic indices
 * required.*/
void chi_math::AngularQuadrature::
  MakeHarmonicIndices(unsigned int scattering_order, int dimension)
{
  m_to_ell_em_map.clear();

  if (dimension == 1)
    for (int ell=0; ell<=scattering_order; ell++)
      m_to_ell_em_map.emplace_back(ell,0);
  else if (dimension == 2)
    for (int ell=0; ell<=scattering_order; ell++)
      for (int m=-ell; m<=ell; m+=2)
        m_to_ell_em_map.emplace_back(ell,m);
  else if (dimension == 3)
    for (int ell=0; ell<=scattering_order; ell++)
      for (int m=-ell; m<=ell; m++)
        m_to_ell_em_map.emplace_back(ell,m);
}

//###################################################################
/**Computes the discrete to moment operator.*/
void chi_math::AngularQuadrature::
  BuildDiscreteToMomentOperator(unsigned int scattering_order, int dimension)
{
  if (d2m_op_built) return;

  d2m_op.clear();
  MakeHarmonicIndices(scattering_order,dimension);

  const size_t num_angles = abscissae.size();
  const size_t num_moms = m_to_ell_em_map.size();

  for (const auto& ell_em : m_to_ell_em_map)
  {
    std::vector<double> cur_mom;
    cur_mom.reserve(num_angles);

    for (int n=0; n<num_angles; n++)
    {
      const auto& cur_angle = abscissae[n];
      double value = chi_math::Ylm(ell_em.ell,ell_em.m,
                                   cur_angle.phi,
                                   cur_angle.theta);
      double w = weights[n];
      cur_mom.push_back(value*w);
    }

    d2m_op.push_back(cur_mom);
  }
  d2m_op_built = true;

  //=================================== Verbose printout
  std::stringstream outs;
  outs
    << "\nQuadrature d2m operator:\n";
  for (int n=0; n<num_angles; n++)
  {
    outs << std::setw(5) << n;
    for (int m=0; m<num_moms; m++)
    {
      outs
        << std::setw(15) << std::left << std::fixed
        << std::setprecision(10) << d2m_op[m][n] << " ";
    }
    outs << "\n";
  }

  chi::log.Log0Verbose1() << outs.str();
}

//###################################################################
/**Computes the moment to discrete operator.*/
void chi_math::AngularQuadrature::
  BuildMomentToDiscreteOperator(unsigned int scattering_order, int dimension)
{
  if (m2d_op_built) return;

  m2d_op.clear();
  MakeHarmonicIndices(scattering_order,dimension);

  const size_t num_angles = abscissae.size();
  const size_t num_moms = m_to_ell_em_map.size();

  const auto normalization =
    std::accumulate(weights.begin(), weights.end(), 0.0);

  for (const auto& ell_em : m_to_ell_em_map)
  {
    std::vector<double> cur_mom;
    cur_mom.reserve(num_angles);

    for (int n=0; n<num_angles; n++)
    {
      const auto& cur_angle = abscissae[n];
      double value = ((2.0*ell_em.ell+1.0)/normalization)*
                     chi_math::Ylm(ell_em.ell,ell_em.m,
                                   cur_angle.phi,
                                   cur_angle.theta);
      cur_mom.push_back(value);
    }

    m2d_op.push_back(cur_mom);
  }//for m
  m2d_op_built = true;

  //=================================== Verbose printout
  std::stringstream outs;

  outs
    << "\nQuadrature m2d operator:\n";
  for (int n=0; n<num_angles; n++)
  {
    outs << std::setw(5) << n;
    for (int m=0; m<num_moms; m++)
    {
      outs
        << std::setw(15) << std::left << std::fixed
        << std::setprecision(10) << m2d_op[m][n] << " ";
    }
    outs << "\n";
  }

  chi::log.Log0Verbose1() << outs.str();
}

//###################################################################
/**Returns a reference to the precomputed d2m operator. This will
 * throw a std::logic_error if the operator has not been built yet.
 * The operator is accessed as [m][d], where m is the moment index
 * and d is the direction index.*/
std::vector<std::vector<double>> const&
  chi_math::AngularQuadrature::GetDiscreteToMomentOperator() const
{
  const std::string fname = __FUNCTION__;
  if (not d2m_op_built)
    throw std::logic_error(fname + ": Called but D2M operator not yet built. "
           "Make a call to BuildDiscreteToMomentOperator before using this.");
  return d2m_op;
}

//###################################################################
/**Returns a reference to the precomputed m2d operator. This will
 * throw a std::logic_error if the operator has not been built yet.
 * The operator is accessed as [m][d], where m is the moment index
 * and d is the direction index.*/
std::vector<std::vector<double>> const&
  chi_math::AngularQuadrature::GetMomentToDiscreteOperator() const
{
  const std::string fname = __FUNCTION__;
  if (not m2d_op_built)
    throw std::logic_error(fname + ": Called but M2D operator not yet built. "
           "Make a call to BuildMomentToDiscreteOperator before using this.");
  return m2d_op;
}

//###################################################################
/**Returns a reference to the precomputed harmonic index map. This will
 * throw a std::logic_error if the map has not been built yet.*/
const std::vector<chi_math::AngularQuadrature::HarmonicIndices>&
  chi_math::AngularQuadrature::GetMomentToHarmonicsIndexMap() const
{
  const std::string fname = __FUNCTION__;
  if (not (d2m_op_built or m2d_op_built))
    throw std::logic_error(fname + ": Called but map not yet built. "
           "Make a call to either BuildDiscreteToMomentOperator or"
           "BuildMomentToDiscreteOperator before using this.");
  return m_to_ell_em_map;
}

//###################################################################
/**Constructor using custom directions.*/
chi_math::AngularQuadratureCustom::
  AngularQuadratureCustom(std::vector<double> &azimuthal,
                          std::vector<double> &polar,
                          std::vector<double> &in_weights, bool verbose)
{
  size_t Na = azimuthal.size();
  size_t Np = polar.size();
  size_t Nw = in_weights.size();

  if ((Na-Np != 0) or (Na-Nw != 0))
  {
    chi::log.LogAllError()
      << "chi_math::AngularQuadrature::InitializeWithCustom: supplied"
         " vectors need to be of equal length.";
    chi::Exit(EXIT_FAILURE);
  }

  //================================================== Create angle pairs
  std::stringstream ostr;
  double weight_sum = 0.0;

  for (unsigned int i=0; i<Na; i++)
  {
    const auto abscissa =
      chi_math::QuadraturePointPhiTheta(azimuthal[i], polar[i]);

    abscissae.push_back(abscissa);

    const double weight = in_weights[i];
    weights.push_back(weight);
    weight_sum += weight;

    if (verbose)
    {
      char buf[200];
      snprintf(buf,200,"Varphi=%.2f Theta=%.2f Weight=%.3e\n",
              abscissa.phi*180.0/M_PI,
              abscissa.theta*180.0/M_PI,
              weight);
      ostr << buf;
    }
  }

  //================================================== Create omega list
  for (const auto& qpoint : abscissae)
  {
    chi_mesh::Vector3 new_omega;
    new_omega.x = sin(qpoint.theta)*cos(qpoint.phi);
    new_omega.y = sin(qpoint.theta)*sin(qpoint.phi);
    new_omega.z = cos(qpoint.theta);

    omegas.push_back(new_omega);
  }

  if (verbose)
  {
    chi::log.Log()
      << ostr.str() << "\n"
      << "Weight sum=" << weight_sum;
  }
}