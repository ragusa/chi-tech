#include "angular_quadrature_base2.h"

#include "ChiMath/Quadratures/LegendrePoly/legendrepoly.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <iomanip>
#include <numeric>

/**Populates a map of moment m to the angular_Spherical Harmonic indices
 * required.*/
void chi_math::AngularQuadrature::
  MakeHarmonicIndices(unsigned int scattering_order)
{
  // short cut
  unsigned int L = scattering_order;

  m_to_ell_em_map.clear();
  if (not galerkin) {
    if (dimension == 1)
    {
      num_moms = L + 1;
      for (int ell = 0; ell <= scattering_order; ++ell)
        m_to_ell_em_map.emplace_back(ell, 0);
    }
    else if (dimension == 2)
    {
      num_moms = (L + 1) * (L + 2) / 2 ;
      for (int ell = 0; ell <= scattering_order; ++ell)
        for (int m = -ell; m <= ell; m += 2)
          m_to_ell_em_map.emplace_back(ell, m);
    }
    else if (dimension == 3)
    {
      num_moms = (L + 1) * (L + 1) ;
      for (int ell = 0; ell <= scattering_order; ++ell)
        for (int m = -ell; m <= ell; ++m)
          m_to_ell_em_map.emplace_back(ell, m);
    }
  }
  else
  {
    chi::log.Log0Error() << "Galerkin not coded";
    chi::Exit(123);
  }
}

//###################################################################
/**Computes the discrete to moment operator.*/
void chi_math::AngularQuadrature::
  BuildDiscreteToMomentUsual(unsigned int scattering_order)
{
  d2m.resize(num_moms, std::vector<double>(num_dirs));

  unsigned int mom = 0;
  for (const auto& ell_em : m_to_ell_em_map)
  {
    std::vector<double> cur_mom;
    cur_mom.resize(num_dirs, 0.);

    for (int d=0; d<num_dirs; ++d)
    {
      const auto& cur_angle = abscissae[d];
      double sph_value = chi_math::Ylm(ell_em.ell,ell_em.m,
                                   cur_angle.phi,
                                   cur_angle.theta);
      cur_mom[d] = sph_value * weights[d];
    }//for d
    d2m[mom] = cur_mom;
    mom += 1;
  }

  //=================================== Verbose printout
  std::stringstream outs;
  outs
    << "\nQuadrature d2m operator:\n";
  for (int d=0; d<num_dirs; ++d)
  {
    outs << std::setw(5) << d;
    for (int m=0; m<num_moms; ++m)
    {
      outs
        << std::setw(15) << std::left << std::fixed
        << std::setprecision(10) << d2m[m][d] << " ";
    }
    outs << "\n";
  }

  chi::log.Log0Verbose1() << outs.str();
}

//###################################################################
/**Computes the moment to discrete operator.*/
void chi_math::AngularQuadrature::
  BuildMomentToDiscreteUsual(unsigned int scattering_order)
{
  m2d.resize(num_dirs, std::vector<double>(num_moms));

//  const auto normalization =
//    std::accumulate(weights.begin(), weights.end(), 0.0);

  for (int d=0; d<num_dirs; ++d)
  {
    const auto& cur_angle = abscissae[d];

    std::vector<double> cur_dir;
    cur_dir.reserve(num_moms);

    unsigned int mom = 0;
    for (const auto& ell_em : m_to_ell_em_map)
    {
      double value = ((2.0*ell_em.ell+1.0)/normalization)*
                     chi_math::Ylm(ell_em.ell,ell_em.m,
                                   cur_angle.phi,
                                   cur_angle.theta);
      cur_dir[mom] = value;
    }
    m2d[d] = cur_dir;
  }//for d

  //=================================== Verbose printout
  std::stringstream outs;

  outs
    << "\nQuadrature m2d operator:\n";
  for (int d=0; d<num_dirs; ++d)
  {
    outs << std::setw(5) << d;
    for (int m=0; m<num_moms; ++m)
    {
      outs
        << std::setw(15) << std::left << std::fixed
        << std::setprecision(10) << m2d[d][m] << " ";
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
  if (not d2m_built)
    throw std::logic_error(fname + ": Called but D2M operator not yet built. "
           "Make a call to BuildDiscreteToMomentOperator before using this.");
  return d2m;
}

//###################################################################
/**Returns a reference to the precomputed m2d operator. This will
 * throw a std::logic_error if the operator has not been built yet.
 * The operator is accessed as [d][m], where m is the moment index
 * and d is the direction index.*/
std::vector<std::vector<double>> const&
  chi_math::AngularQuadrature::GetMomentToDiscreteOperator() const
{
  const std::string fname = __FUNCTION__;
  if (not m2d_built)
    throw std::logic_error(fname + ": Called but M2D operator not yet built. "
           "Make a call to BuildMomentToDiscreteOperator before using this.");
  return m2d;
}

//###################################################################
/**Returns a reference to the precomputed harmonic index map. This will
 * throw a std::logic_error if the map has not been built yet.*/
const std::vector<chi_math::AngularQuadrature::HarmonicIndices>&
  chi_math::AngularQuadrature::GetMomentToHarmonicsIndexMap() const
{
  const std::string fname = __FUNCTION__;
  if (not (d2m_built or m2d_built))
    throw std::logic_error(fname + ": Called but map not yet built. "
           "Make a call to either BuildDiscreteToMomentOperator or"
           "BuildMomentToDiscreteOperator before using this.");
  return m_to_ell_em_map;
}

//###################################################################
///**Constructor using custom directions.*/
//chi_math::AngularQuadratureCustom::
//  AngularQuadratureCustom(std::vector<double> &azimuthal,
//                          std::vector<double> &polar,
//                          std::vector<double> &in_weights, bool verbose)
//{
//  size_t Na = azimuthal.size();
//  size_t Np = polar.size();
//  size_t Nw = in_weights.size();
//
//  if ((Na-Np != 0) or (Na-Nw != 0))
//  {
//    chi::log.LogAllError()
//      << "chi_math::AngularQuadrature::InitializeWithCustom: supplied"
//         " vectors need to be of equal length.";
//    chi::Exit(EXIT_FAILURE);
//  }
//
//  //================================================== Create angle pairs
//  std::stringstream ostr;
//  double weight_sum = 0.0;
//
//  for (unsigned int i=0; i<Na; i++)
//  {
//    const auto abscissa =
//      chi_math::QuadraturePointPhiTheta(azimuthal[i], polar[i]);
//
//    abscissae.push_back(abscissa);
//
//    const double weight = in_weights[i];
//    weights.push_back(weight);
//    weight_sum += weight;
//
//    if (verbose)
//    {
//      char buf[200];
//      snprintf(buf,200,"Varphi=%.2f Theta=%.2f Weight=%.3e\n",
//              abscissa.phi*180.0/M_PI,
//              abscissa.theta*180.0/M_PI,
//              weight);
//      ostr << buf;
//    }
//  }
//
//  //================================================== Create omega list
//  for (const auto& qpoint : abscissae)
//  {
//    chi_mesh::Vector3 new_omega;
//    new_omega.x = sin(qpoint.theta)*cos(qpoint.phi);
//    new_omega.y = sin(qpoint.theta)*sin(qpoint.phi);
//    new_omega.z = cos(qpoint.theta);
//
//    omegas.push_back(new_omega);
//  }
//
//  if (verbose)
//  {
//    chi::log.Log()
//      << ostr.str() << "\n"
//      << "Weight sum=" << weight_sum;
//  }
//}