#ifndef CURVILINEAR_ANGULAR_QUADRATURE_H
#define CURVILINEAR_ANGULAR_QUADRATURE_H

#include "ChiMath/Quadratures/angular_product_quadrature.h"

namespace chi_math
{
  class CurvilinearAngularQuadrature;
}

/** Base class for curvilinear angular quadratures (product angular
 *  quadratures with additional direction-dependent parameters).
 */
class chi_math::CurvilinearAngularQuadrature : public chi_math::ProductQuadrature
{
//  Attributes
protected:
  /** Factor to account for angular diamond differencing. */
  std::vector<double> fac_diamond_difference;
  /** Factor to account for discretisation of the component of the streaming
   *  operator that contains the angular derivative. */
  std::vector<double> fac_streaming_operator;

//  Methods
public:
  /** Return constant reference to fac_diamond_difference. */
  const std::vector<double>& GetDiamondDifferenceFactor() const
  { return fac_diamond_difference; }
  /** Return constant reference to fac_streaming_operator. */
  const std::vector<double>& GetStreamingOperatorFactor() const
  { return fac_streaming_operator; }
protected:
  /** Default constructor. */
  CurvilinearAngularQuadrature() = default;
  /** Default destructor. */
  virtual ~CurvilinearAngularQuadrature() = default;
};

#endif // CURVILINEAR_ANGULAR_QUADRATURE_H
