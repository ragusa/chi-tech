#include "quadrature_gausslegendre.h"
#include "ChiMath/Quadratures/LegendrePoly/legendrepoly.h"
#include <cmath>

#include "chi_runtime.h"
#include "chi_log.h"

#include <algorithm>

//###################################################################
/**Populates the abscissae and weights for a Gauss-Legendre
 * quadrature given the degree \f$ p \f$ of the mononomial such that
 * the quadrature rule integrates exactly the weighted integrand
 * \f$ \rho(x) x^{p} \f$, with \f$ \rho(x) := 1 \f$,
 * on the interval \f$ [-1;+1] \f$.
 * The number of points generated will be ceil((O+1)/2).*/
chi_math::QuadratureGaussLegendre::
  QuadratureGaussLegendre(QuadratureOrder in_order,
                          bool verbose,
                          unsigned int max_iters,
                          double tol)
  : chi_math::Quadrature(in_order)
{
  const unsigned int N = std::ceil( ((int)order + 1)/2.0 );
  Initialize(N, verbose, max_iters, tol);
}

//###################################################################
/**Populates the abscissae and weights for a Gauss-Legendre
 * quadrature given the number of desired quadrature points. The
 * order of the quadrature will be 2N-1.*/
chi_math::QuadratureGaussLegendre::
  QuadratureGaussLegendre(unsigned int N,
                          bool verbose,
                          unsigned int max_iters,
                          double tol)
  : chi_math::Quadrature((QuadratureOrder)(2*N-1))
{
  Initialize(N, verbose, max_iters, tol);
}

//###################################################################
/**Populates the abscissae and weights for a Gauss-Legendre
 * quadrature given the number of desired quadrature points.*/
void
chi_math::QuadratureGaussLegendre::Initialize(unsigned int N, bool verbose,
                                              unsigned int max_iters, double tol)
{
  switch (order)
  {
    default:
    {
      if (verbose)
        chi::log.Log() << "Initializing Gauss-Legendre Quadrature "
                         "with " << N << " q-points";

      //========================= Compute the roots
      auto roots = FindRoots(N, max_iters, tol);
      for (auto v : roots) qpoints.emplace_back(v);

      //========================= Compute the weights
      weights.resize(N,1.0);
      for (size_t k=0; k < qpoints.size(); k++)
      {
        weights[k] =
          2.0 * (1.0 - qpoints[k][0] * qpoints[k][0]) /
          ((N + 1) * (N + 1) *
           Legendre(N+1, qpoints[k][0]) * Legendre(N + 1, qpoints[k][0]) );

        if (verbose)
          chi::log.Log()
            << "root[" << k << "]=" << qpoints[k][0]
            << ", weight=" << weights[k];
      }//for abscissae

      break;
    }
  }//switch order

  range = {-1, +1};
}


//###################################################################
/** Finds the roots of the Legendre polynomial.
 *
 * The algorithm is that depicted in:
 *
 * [1] Barrera-Figueroa, et al., "Multiple root finder algorithm for Legendre
 *     and Chebyshev polynomials via Newton's method", Annales Mathematicae et
 *     Informaticae, 33 (2006) pp. 3-13.
 *
 * \param N Is the order of the polynomial.
 * \param roots Is a reference to the roots.
 * \param max_iters Maximum newton iterations to perform for each root.
 *        Default: 1000.
 * \param tol Tolerance at which the newton iteration will be terminated.
 *        Default: 1.0e-12.
 *
 * \author Jan*/
std::vector<double> chi_math::QuadratureGaussLegendre::FindRoots(
  unsigned int N, unsigned int max_iters, double tol)
{
  //======================================== Populate init guess
  //This initial guess proved to be quite important
  //at higher N since the roots start to get
  //squeezed to -1 and 1.
  int num_search_intvls = 1000;
  if (N>64)
    num_search_intvls *= 10;
  if (N>256)
    num_search_intvls *= 10;
  if (N>768)
    num_search_intvls *= 10;

  if (N>2056)
  {
    num_search_intvls *= 10;
    chi::log.Log0Warning()
      << "chi_math::QuadratureGaussLegendre::FindRoots: "
      << "The order of the polynomial for which to find the roots is "
      << "greater than 2056. Accuracy of the root finder will be diminished "
      << "along with a reduction in stability.";
  }

  // For this code we simply check to see where the
  // polynomial changes sign.
  double delta = 2.0/num_search_intvls;
  std::vector<double> xk(N, 0.0);
  int counter = -1;
  for(size_t i=0; i<num_search_intvls; i++)
  {
    double x_i = -1.0 + i*delta;
    double x_ip1 = x_i + delta;

    if (Legendre(N,x_i)*Legendre(N,x_ip1) < 0.0)
      xk[++counter] = (x_ip1 + x_i) / 2.0;
  }

  //======================================== Apply algorithm
  // Refer to equation 4.3 in [1]. Sum 1 (S1) is used in the
  // computation of B at x_k. Sum 2 (S2) is used in equation 4.3.
  // Equation 4.3 is broken up into pieces as follows:
  //  - a = block bracket containing the second derivative
  //  - b = denominator
  //  - c = everything but xold
  for (int k=0; k<N; k++)
  {
    for (size_t iteration=0; iteration<max_iters; iteration++)
    {
      double xold = xk[k];
      double f   = Legendre(N,xold);      //Function evaluation
      double fp  = dLegendredx(N,xold);   //First derivative
      double fpp = d2Legendredx2(N,xold); //Second derivative

      //===================== Compute sum 1
      double S1 = 0.0;
      for (int i=0; i<=(k-1); i++)
        S1 += 1.0/(xk[k] - xk[i]);

      //===================== Compute B at x_k
      double B_xk = fp - f*S1;

      //===================== Compute sum 2
      double S2 = 0.0;
      for (int i=0; i<=(k-1); i++)
        S2 += 1.0 / (xk[k] - xk[i]) / (xk[k] - xk[i]);

      //===================== Compute final formula
      double a    = fpp + f*S2;
      double b    = B_xk*B_xk + fp*fp - f*a;
      double c    = 2.0*f*B_xk/b;

      xk[k] = xold - c;

      if (std::fabs(xk[k] - xold) < tol)
        break;
    }//for iteration
  }//for k

  std::stable_sort(xk.begin(), xk.end());

  return xk;
}