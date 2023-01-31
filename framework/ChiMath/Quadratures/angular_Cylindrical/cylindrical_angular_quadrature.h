#ifndef CYLINDRICAL_ANGULAR_QUADRATURE_H
#define CYLINDRICAL_ANGULAR_QUADRATURE_H

#include "ChiMath/Quadratures/angular_curvilinear_quadrature.h"

#include "ChiMath/Quadratures/quadrature.h"

namespace chi_math
{
  class CylindricalAngularQuadrature;
}

/** angular_Cylindrical product angular quadrature. */
class chi_math::CylindricalAngularQuadrature : public chi_math::CurvilinearAngularQuadrature
{
//  Methods
public:
  /** Effective constructor. Initialize with one-dimensional quadratures:
   *  the azimuthal quadrature is applied at each polar level.
   *  If not already present in the azimuthal quadrature, the method inserts
   *  the azimuthal starting directions. */
  CylindricalAngularQuadrature(
    const chi_math::Quadrature& quad_polar,
    const chi_math::Quadrature& quad_azimu,
    const bool verbose=false);
  /** Effective constructor. Initialize with one-dimensional quadratures:
   *  a possibly diverse azimuthal quadrature is applied at each polar level.
   *  If not already present in the azimuthal quadrature, the method inserts
   *  the azimuthal starting directions. */
  CylindricalAngularQuadrature(
    const chi_math::Quadrature& quad_polar,
    const std::vector<chi_math::Quadrature>& quad_azimu_vec,
    const bool verbose=false);
  /** Default destructor. */
  virtual ~CylindricalAngularQuadrature() = default;

  void MakeHarmonicIndices(unsigned int scattering_order, int dimension) override;
private:
  /** Initialize with one-dimensional quadratures: a polar quadrature and
   *  a possibly unique azimuthal quadrature for each polar level. */
  void Initialize(const chi_math::Quadrature& quad_polar,
                  const std::vector<chi_math::Quadrature>& quad_azimu_vec,
                  const bool verbose=false);
  /** Initialize parametrizing factors of the cylindrical angular quadrature,
   *  starting from a fully initialized underlying product quadrature. */
  void InitializeParameters();
};

#endif // CYLINDRICAL_ANGULAR_QUADRATURE_H
