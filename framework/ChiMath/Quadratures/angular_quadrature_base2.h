#ifndef _angular_quadrature_base_h
#define _angular_quadrature_base_h

#include <vector>

#include "ChiMesh/chi_mesh.h"

namespace chi_math
{
  struct QuadraturePointPhiTheta;

  enum class AngularQuadratureType
  {
    GaussLegendre        = 1,
    GaussLobatto         = 2, //jcr
    ProductQuadrature    = 3,
    TriangularQuadrature = 4,
    SLDFESQ              = 5
  };
  class AngularQuadrature;
}

/**Simple structure to add angle values to directions.*/
struct chi_math::QuadraturePointPhiTheta
{
  double phi=0.0;
  double theta=0.0;
  QuadraturePointPhiTheta(const double phi, const double theta)
  : phi(phi), theta(theta)
  {
  }
};

//################################################################### Class def
/**Base class for angular quadratures.*/
class chi_math::AngularQuadrature
{
public:
  const chi_math::AngularQuadratureType type;
public:
  // direction cosines
  std::vector<chi_mesh::Vector3>                 omegas;
  // weights
  std::vector<double>                            weights;
  // angle values phi and theta for each omega (used in aggregation)
  std::vector<chi_math::QuadraturePointPhiTheta> abscissae;
  // dimension of the problem
  int dimension;
  // number of directions
  int num_dirs;
  // number of moments
  int num_moms;
  // is galerkin_quad ?
  bool galerkin = false;
  // normalization
  double normalization;

  // mapping from single index to {ell,m} pairs
  struct HarmonicIndices
  {
    unsigned int ell = 0;
    int          m   = 0;

    HarmonicIndices()=default;
    HarmonicIndices(unsigned int in_ell, int in_m) : ell(in_ell),m(in_m)
    {
    }

    bool operator==(const HarmonicIndices& other) const
    {
      return (ell == other.ell and m == other.m);
    }
  };

protected:
  std::vector<std::vector<double>> d2m;
  std::vector<std::vector<double>> m2d;
  std::vector<HarmonicIndices>     m_to_ell_em_map;
  bool                             d2m_built = false;
  bool                             m2d_built = false;

public:
  // constructor
  explicit
  AngularQuadrature(chi_math::AngularQuadratureType in_type) :
    type(in_type)
  {}

  // destructor
  virtual ~AngularQuadrature() = default;


  // d2m and m2d operators
  virtual void BuildDiscreteToMomentUsual(unsigned int scattering_order);
  virtual void BuildMomentToDiscreteUsual(unsigned int scattering_order);

  // getters
  std::vector<std::vector<double>> const&
  GetDiscreteToMomentOperator() const;

  std::vector<std::vector<double>> const&
  GetMomentToDiscreteOperator() const;

  const std::vector<HarmonicIndices>&
  GetMomentToHarmonicsIndexMap() const;

protected:
  virtual void MakeHarmonicIndices(unsigned int scattering_order);

};

#endif
