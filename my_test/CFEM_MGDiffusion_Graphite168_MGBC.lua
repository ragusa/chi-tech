function printTable(t, f)
    local function printTableHelper(obj, cnt)
        local cnt = cnt or 0
        if type(obj) == "table" then
            io.write("\n", string.rep("\t", cnt), "{\n")
            cnt = cnt + 1
            for k,v in pairs(obj) do
                if type(k) == "string" then
                    io.write(string.rep("\t",cnt), '["'..k..'"]', ' = ')
                end
                if type(k) == "number" then
                    io.write(string.rep("\t",cnt), "["..k.."]", " = ")
                end
                printTableHelper(v, cnt)
                io.write(",\n")
            end
            cnt = cnt-1
            io.write(string.rep("\t", cnt), "}")
        elseif type(obj) == "string" then
            io.write(string.format("%q", obj))
        else
            io.write(tostring(obj))
        end
    end

    if f == nil then
        printTableHelper(t)
    else
        io.output(f)
        io.write("return")
        printTableHelper(t)
        io.output(io.stdout)
    end
end

--############################################### Setup mesh
chiMeshHandlerCreate()

mesh={}
N=2
L=20
xmin = -L/2
dx = L/N
for i=1,(N+1) do
    k=i-1
    mesh[i] = xmin + k*dx
end

chiMeshCreateUnpartitioned2DOrthoMesh(mesh,mesh)
chiVolumeMesherExecute();
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)

--############################################### Set Material IDs
 chiVolumeMesherSetMatIDToAll(0)


-- Setboundary IDs
-- xmin,xmax,ymin,ymax,zmin,zmax
e_vol = chiLogicalVolumeCreate(RPP,0.99999,1000,-1000,1000,-1000,1000)
w_vol = chiLogicalVolumeCreate(RPP,-1000,-0.9999,-1000,1000,-1000,1000)
n_vol = chiLogicalVolumeCreate(RPP,-1000,1000,0.99999,1000,-1000,1000)
s_vol = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,-0.99999,-1000,1000)

e_bndry = 0
w_bndry = 1
n_bndry = 2
s_bndry = 3

chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,e_vol,e_bndry)
chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,w_vol,w_bndry)
chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,n_vol,n_bndry)
chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,s_vol,s_bndry)

--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Mat_outer");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)

num_groups = 2
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        CHI_XSFILE,"../my_test/xs_2g_mat1up.cxs")
num_groups = 168
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        CHI_XSFILE,"../ChiTest/xs_graphite_pure.cxs")

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)

--############################################### Add external src
src={}
for g=1,num_groups do
    src[g] = 0.0
end
src[1] = 1.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--############################################### Add material properties
--#### CFEM stuff
phys1 = chiCFEMMGDiffusionSolverCreate()

chiSolverSetBasicOption(phys1, "residual_tolerance", 1E-8)
chiSolverSetBasicOption(phys1, "thermal_flux_error", 1E-7)
chiSolverSetBasicOption(phys1, "max_thermal_iters", 1)
chiSolverSetBasicOption(phys1, "verbose_level", 1)

atab = {}
btab = {}
ftab = {}
for g=1,num_groups do
    atab[g] = 0.25
    btab[g] = 0.5
    ftab[g] = 0.0
end

--chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",e_bndry,"neumann",atab)
--chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",w_bndry,"neumann",atab)
--chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",n_bndry,"neumann",atab)
--chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",s_bndry,"neumann",atab)

--chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",e_bndry,"vacuum")
--chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",w_bndry,"vacuum")
--chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",n_bndry,"vacuum")
--chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",s_bndry,"vacuum")

chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",e_bndry,"reflecting")
chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",w_bndry,"reflecting")
chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",n_bndry,"reflecting")
chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",s_bndry,"reflecting")

--chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",e_bndry,"robin",atab,btab,ftab)
--chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",w_bndry,"robin",atab,btab,ftab)
--chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",n_bndry,"robin",atab,btab,ftab)
--ftab[1] = 1.0
--chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",s_bndry,"robin",atab,btab,ftab)

chiSolverInitialize(phys1)
chiSolverExecute(phys1)

----############################################### Visualize the field function
fflist,count = chiGetFieldFunctionList(phys1)
chiExportMultiFieldFunctionToVTK(fflist,"square_dirG")
--for g=1,num_groups do
--    g_string=string.format("%03d",g)
--    chiExportFieldFunctionToVTK(fflist[g],"square_up_flx"..g_string,"Flux_Diff"..g_string)
--end

ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1
chiFFInterpolationSetProperty(curffi,OPERATION,OP_AVG)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[1])
--chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[2])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
print(chiFFInterpolationGetValue(curffi))

--printTable(fflist)