#include "ChiLua/chi_lua.h"

#include "chi_runtime.h"

#include "ChiMath/Quadratures/angular_product_quadrature.h"

#include "chi_log.h"

#include <memory>

//########################################################## Create empty system
/** Creates a Product-quadrature.
 *
\param QuadratureType int Quadrature identifier.
\param values varying Varying options based on the quadrature type.

##_

###QuadratureType:
GAUSS_LEGENDRE\n
 Gauss-Legendre quadrature for the polar angles and no quadrature rule
 for the azimuthal angle. Suitable only for 1D simulations. Expects
 to be followed by the number of angles Np. Optionally a verbosity flag
 can be added.\n\n

GAUSS_LEGENDRE_LEGENDRE\n
 Gauss-Legendre quadrature for both the polar and azimuthal dimension.
 Expects to be followed by number of Azimuthal and Polar angles.
 Optionally a verbosity flag can be added.\n\n

GAUSS_LEGENDRE_CHEBYSHEV\n
 Gauss-Legendre quadrature for the polar angle but Gauss-Chebyshev
 for the azimuthal angle.
 Expects to be followed by number of Azimuthal and Polar angles.
 Optionally a verbosity flag can be added.\n\n

CUSTOM_QUADRATURE\n
 Expects to be followed by three lua tables. The first table is an array,
 of length Na, of the azimuthal angles (radians). The second table is an array,
 of length Np, of the polar angles (radians). The third table is an array, of
 length Na*Np, and contains the weight associated with each angle pair.
 Optionally a verbosity flag can be added.\n\n


\return Returns a unique handle to the created product quadrature rule

\ingroup LuaQuadrature
\author Jan*/
int chiCreateProductQuadrature(lua_State *L)
{
  int num_args = lua_gettop(L);
  //============================================= Parse argument
  int ident = lua_tonumber(L,1);
  bool verbose = false;



  if (ident == (int)chi_math::ProductQuadratureType::GAUSS_LEGENDRE)
  {
    if (num_args<2)
      LuaPostArgAmountError("chiCreateProductQuadrature",2,num_args);

    int Np = lua_tonumber(L,2);
    if (num_args == 3)
      verbose = lua_toboolean(L,3);

    chi::log.Log() << "Creating Gauss-Legendre Quadrature\n";

    auto new_quad =
      std::make_shared<chi_math::AngularQuadratureProdGL>(Np, verbose);

    chi::angular_quadrature_stack.push_back(new_quad);
    const size_t index = chi::angular_quadrature_stack.size() - 1;
    lua_pushnumber(L,static_cast<lua_Number>(index));

    if (verbose)
    {
      chi::log.Log()
        << "Created Gauss-Legendre Quadrature with "
        << new_quad->azimu_ang.size()
        << " azimuthal angles and "
        << new_quad->polar_ang.size()
        << " polar angles.";
    }

    return 1;
  }
//  else if (ident == (int)chi_math::ProductQuadratureType::GAUSS_LEGENDRE_LEGENDRE)
//  {
//    if (num_args<3)
//      LuaPostArgAmountError("chiCreateProductQuadrature",3,num_args);
//
//    int Np = lua_tonumber(L,2);
//    int Na = lua_tonumber(L,3);
//    if (num_args == 4)
//      verbose = lua_toboolean(L,4);
//
//    chi::log.Log() << "Creating Gauss-Legendre-Legendre Quadrature\n";
//
//    auto new_quad =
//      std::make_shared<chi_math::AngularQuadratureProdGLL>(Np,Na,verbose);
//
//    chi::angular_quadrature_stack.push_back(new_quad);
//    const size_t index = chi::angular_quadrature_stack.size() - 1;
//    lua_pushnumber(L,static_cast<lua_Number>(index));
//
//    if (verbose)
//    {
//      chi::log.Log()
//        << "Created Gauss-Legendre-Legendre Quadrature with "
//        << new_quad->azimu_ang.size()
//        << " azimuthal angles and "
//        << new_quad->polar_ang.size()
//        << " polar angles.";
//    }
//
//    return 1;
//  }
  else if (ident == (int)chi_math::ProductQuadratureType::GAUSS_LEGENDRE_CHEBYSHEV)
  {
    if (num_args<3)
      LuaPostArgAmountError("chiCreateProductQuadrature",3,num_args);

    int Np = lua_tonumber(L,2);
    int Na = lua_tonumber(L,3);
    if (num_args == 4)
      verbose = lua_toboolean(L,4);

    chi::log.Log() << "Creating Gauss-Legendre-ChebyShev Quadrature\n";

    auto new_quad =
      std::make_shared<chi_math::AngularQuadratureProdGLC>(Np,Na,verbose);

    chi::angular_quadrature_stack.push_back(new_quad);
    const size_t index = chi::angular_quadrature_stack.size() - 1;
    lua_pushnumber(L,static_cast<lua_Number>(index));

    if (verbose)
    {
      chi::log.Log()
      << "Created Gauss-Legendre-Chebyshev Quadrature with "
      << new_quad->azimu_ang.size()
      << " azimuthal angles and "
      << new_quad->polar_ang.size()
      << " polar angles.";
    }

    return 1;
  }
  else if (ident == (int)chi_math::ProductQuadratureType::CUSTOM_QUADRATURE)
  {
    if (num_args<4)
      LuaPostArgAmountError("chiCreateProductQuadrature:CUSTOM_QUADRATURE",3,num_args);

    if (not lua_istable(L,2))
    {
      chi::log.LogAllError()
        << "chiCreateProductQuadrature:CUSTOM_QUADRATURE, second argument must "
        << "be a lua table.";
     chi::Exit(EXIT_FAILURE);
    }
    if (not lua_istable(L,3))
    {
      chi::log.LogAllError()
        << "chiCreateProductQuadrature:CUSTOM_QUADRATURE, third argument must "
        << "be a lua table.";
     chi::Exit(EXIT_FAILURE);
    }
    if (not lua_istable(L,4))
    {
      chi::log.LogAllError()
        << "chiCreateProductQuadrature:CUSTOM_QUADRATURE, fourth argument must "
        << "be a lua table.";
     chi::Exit(EXIT_FAILURE);
    }
    if (num_args == 5)
      verbose = lua_toboolean(L,4);

    size_t Na = lua_rawlen(L,2);
    size_t Np = lua_rawlen(L,3);
    size_t Nw = lua_rawlen(L,4);

    std::vector<double> azimuthal(Na,0.0);
    std::vector<double> polar(Np,0.0);
    std::vector<double> weights(Nw,0.0);

    for (int n=1; n<=Na; ++n)
    {
      lua_pushnumber(L,n);
      lua_gettable(L,2);
      azimuthal[n-1] = lua_tonumber(L,-1);
      lua_pop(L,1);
    }
    for (int n=1; n<=Np; ++n)
    {
      lua_pushnumber(L,n);
      lua_gettable(L,3);
      polar[n-1] = lua_tonumber(L,-1);
      lua_pop(L,1);
    }
    for (int n=1; n<=Nw; ++n)
    {
      lua_pushnumber(L,n);
      lua_gettable(L,4);
      weights[n-1] = lua_tonumber(L,-1);
      lua_pop(L,1);
    }

    chi::log.Log() << "Creating custom product quadrature Quadrature\n";

    chi::log.Log() << Na << " " << Np << " " << Nw;

    auto new_quad = std::make_shared<chi_math::AngularQuadratureProdCustom>(
      azimuthal, polar, weights, verbose);

    chi::angular_quadrature_stack.push_back(new_quad);
    const size_t index = chi::angular_quadrature_stack.size() - 1;
    lua_pushnumber(L,static_cast<lua_Number>(index));

    if (verbose)
    {
      chi::log.Log()
        << "Created Custom Quadrature with "
        << new_quad->azimu_ang.size()
        << " azimuthal angles and "
        << new_quad->polar_ang.size()
        << " polar angles.";
    }

    return 1;
  }
  else
  {
    chi::log.LogAllError()
      << "In call to chiCreateProductQuadrature. Unsupported quadrature type"
         " supplied. Given: " << ident;
   chi::Exit(EXIT_FAILURE);
  }
  return 0;
}