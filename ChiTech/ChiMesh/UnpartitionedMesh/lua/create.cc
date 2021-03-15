#include "../../../ChiLua/chi_lua.h"

#include "../chi_unpartitioned_mesh.h"
#include "ChiMesh/MeshHandler/chi_meshhandler.h"

/** \defgroup LuaUnpartitionedMesh Unpartitioned Mesh-Reader
 * \ingroup LuaMesh
 */

//###################################################################
/**Creates an unpartitioned mesh from VTK Unstructured mesh files.
 *
 * \image html "InProgressImage.png" width=200px
 *
 * \param file_name char Filename of the .vtu file.
 *
 * \ingroup LuaUnpartitionedMesh
 *
 * \return A handle to the newly created UnpartitionedMesh*/
int chiUnpartitionedMeshFromVTU(lua_State* L)
{
  const char func_name[] = "chiUnpartitionedMeshFromVTU";
  int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError(func_name,1,num_args);

  const char* temp = lua_tostring(L,1);
  auto new_object = new chi_mesh::UnpartitionedMesh;

  chi_mesh::UnpartitionedMesh::Options options;
  options.file_name = std::string(temp);

  new_object->ReadFromVTU(options);

  auto handler = chi_mesh::GetCurrentHandler();
  handler->unpartitionedmesh_stack.push_back(new_object);

  lua_pushnumber(L,handler->unpartitionedmesh_stack.size()-1);

  return 1;
}

//###################################################################
/**Creates an unpartitioned mesh from starccm+ exported
 * Ensight Gold mesh files.
 *
 * \param file_name char Filename of the .case file.
 * \param scale float Scale to apply to the mesh
 *
 * \ingroup LuaUnpartitionedMesh
 *
 * \return A handle to the newly created UnpartitionedMesh*/
int chiUnpartitionedMeshFromEnsightGold(lua_State* L)
{
  const char func_name[] = "chiUnpartitionedMeshFromEnsightGold";
  int num_args = lua_gettop(L);
  if (num_args <1)
    LuaPostArgAmountError(func_name,1,num_args);

  const char* temp = lua_tostring(L,1);
  double scale = 1.0;
  if (num_args >= 2) scale = lua_tonumber(L,2);
  auto new_object = new chi_mesh::UnpartitionedMesh;

  chi_mesh::UnpartitionedMesh::Options options;
  options.file_name = std::string(temp);
  options.scale = scale;

  new_object->ReadFromEnsightGold(options);

  auto handler = chi_mesh::GetCurrentHandler();
  handler->unpartitionedmesh_stack.push_back(new_object);

  lua_pushnumber(L,handler->unpartitionedmesh_stack.size()-1);

  return 1;
}

//###################################################################
/**Creates an unpartitioned mesh from a wavefront .obj file.
 *
 * \param file_name char Filename of the .case file.
 *
 * \ingroup LuaUnpartitionedMesh
 *
 * \return A handle to the newly created UnpartitionedMesh*/
int chiUnpartitionedMeshFromWavefrontOBJ(lua_State* L)
{
  const char func_name[] = "chiUnpartitionedMeshFromWavefrontOBJ";
  int num_args = lua_gettop(L);
  if (num_args <1)
    LuaPostArgAmountError(func_name,1,num_args);

  const char* temp = lua_tostring(L,1);

  auto new_object = new chi_mesh::UnpartitionedMesh;

  chi_mesh::UnpartitionedMesh::Options options;
  options.file_name = std::string(temp);

  new_object->ReadFromWavefrontOBJ(options);

  auto handler = chi_mesh::GetCurrentHandler();
  handler->unpartitionedmesh_stack.push_back(new_object);

  lua_pushnumber(L,handler->unpartitionedmesh_stack.size()-1);

  return 1;
}

//###################################################################
/**Creates an unpartitioned mesh from a .msh file.
 *
 * \param file_name char Filename of the .msh file.
 *
 * \ingroup LuaUnpartitionedMesh
 *
 * \return A handle to the newly created UnpartitionedMesh*/
int chiUnpartitionedMeshFromMshFormat(lua_State* L)
{
  const char func_name[] = "chiUnpartitionedMeshFromMshFormat";
  int num_args = lua_gettop(L);
  if (num_args <1)
    LuaPostArgAmountError(func_name,1,num_args);

  const char* temp = lua_tostring(L,1);

  auto new_object = new chi_mesh::UnpartitionedMesh;

  chi_mesh::UnpartitionedMesh::Options options;
  options.file_name = std::string(temp);

  new_object->ReadFromMsh(options);

  auto handler = chi_mesh::GetCurrentHandler();
  handler->unpartitionedmesh_stack.push_back(new_object);

  lua_pushnumber(L,handler->unpartitionedmesh_stack.size()-1);

  return 1;
}
