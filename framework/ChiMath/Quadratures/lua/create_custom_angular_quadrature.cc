#include "ChiLua/chi_lua.h"

#include "chi_runtime.h"

#include "ChiMath/Quadratures/angular_quadrature_base.h"

#include "chi_log.h"

//########################################################## Create empty system
/** Creates a custom angular quadrature.
 *
\param azimuthal_angles array A lua table with N entries each being an azimuthal
                              angle.
\param polar_angles array A lua table with N entries each being a polar angle
                              angle.
\param weight array A lua table with N entries each being a quadrature weight.


\return Returns a unique handle to the created angular quadrature

\ingroup LuaQuadrature
\author Jan*/
int chiCreateCustomAngularQuadrature(lua_State *L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  if (num_args != 3)
    LuaPostArgAmountError(fname,3,num_args);

  LuaCheckNilValue(fname,L,1);
  LuaCheckNilValue(fname,L,2);
  LuaCheckNilValue(fname,L,3);

  LuaCheckTableValue(fname,L,1);
  LuaCheckTableValue(fname,L,2);
  LuaCheckTableValue(fname,L,3);

  size_t Na = lua_rawlen(L,1);
  size_t Np = lua_rawlen(L,2);
  size_t Nw = lua_rawlen(L,3);

  if ((Na-Np != 0) or (Na-Nw !=0))
  {
    chi::log.LogAllError()
      << fname +": Tables lengths supplied "
         "are not of equal lengths.";
   chi::Exit(EXIT_FAILURE);
  }

  std::vector<double> azi_angles(Na,0.0);
  std::vector<double> pol_angles(Na,0.0);
  std::vector<double> weights(Na,0.0);

  for (int n=1; n<=Na; ++n)
  {
    lua_pushnumber(L,n);
    lua_gettable(L,1);
    azi_angles[n-1] = lua_tonumber(L,-1);
    lua_pop(L,1);
  }

  for (int n=1; n<=Na; ++n)
  {
    lua_pushnumber(L,n);
    lua_gettable(L,2);
    pol_angles[n-1] = lua_tonumber(L,-1);
    lua_pop(L,1);
  }

  for (int n=1; n<=Na; ++n)
  {
    lua_pushnumber(L,n);
    lua_gettable(L,3);
    weights[n-1] = lua_tonumber(L,-1);
    lua_pop(L,1);
  }

  chi::log.Log() << "Creating Custom Angular Quadrature\n";

  auto new_quad = std::make_shared<chi_math::AngularQuadratureCustom>(
    azi_angles, pol_angles, weights, false);

  chi::angular_quadrature_stack.push_back(new_quad);
  size_t index = chi::angular_quadrature_stack.size() - 1;
  lua_pushnumber(L,static_cast<lua_Number>(index));

  return 1;
}