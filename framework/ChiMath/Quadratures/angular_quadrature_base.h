#ifndef _angular_quadrature_base_h
#define _angular_quadrature_base_h

#include <vector>

#include "ChiMesh/chi_mesh.h"

namespace chi_math
{
  struct QuadraturePointPhiTheta;

  enum class AngularQuadratureType
  {
    Arbitrary         = 1,
    ProductQuadrature = 2,
    SLDFESQ           = 3
  };
  class AngularQuadrature;
  class AngularQuadratureCustom;
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
  // angle values phi and theta for each omega (used in aggregation)
  std::vector<chi_math::QuadraturePointPhiTheta> abscissae;
  // weights. when d2m is created, make sure weights are overwritten
  // we keep weights for compatibility.
  std::vector<double>                            weights;
  // dimension of the problem
  int dim;
  // number of directions
  int num_dirs;

  // mapping from single index to {ell,m} pairs
  struct HarmonicIndices
  {
    unsigned int ell=0;
    int          m=0;

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
  std::vector<std::vector<double>> d2m_op;
  std::vector<std::vector<double>> m2d_op;
  std::vector<HarmonicIndices>     m_to_ell_em_map;
  bool                             d2m_op_built = false;
  bool                             m2d_op_built = false;

public:
  // constructors
  AngularQuadrature() :
  type(chi_math::AngularQuadratureType::Arbitrary)
  {}

  explicit
  AngularQuadrature(chi_math::AngularQuadratureType in_type) :
    type(in_type)
  {}

  // destructor
  virtual ~AngularQuadrature() = default;

  // should never be used
//  virtual void OptimizeForPolarSymmetry(double normalization);

  // d2m and m2d operators
  virtual void BuildDiscreteToMomentOperator(unsigned int scattering_order,
                                             int dimension);
  virtual void BuildMomentToDiscreteOperator(unsigned int scattering_order,
                                             int dimension);

  // getters
  std::vector<std::vector<double>> const&
  GetDiscreteToMomentOperator() const;

  std::vector<std::vector<double>> const&
  GetMomentToDiscreteOperator() const;

  const std::vector<HarmonicIndices>&
  GetMomentToHarmonicsIndexMap() const;

protected:
  virtual void MakeHarmonicIndices(unsigned int scattering_order, int dimension);

};

//----------------------- the custom angular quadrature class
class chi_math::AngularQuadratureCustom : public chi_math::AngularQuadrature
{
public:
  AngularQuadratureCustom(std::vector<double>& azimuthal,
                          std::vector<double>& polar,
                          std::vector<double>& in_weights,
                          bool verbose);
};

#endif
