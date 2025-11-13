help([[
This module loads libraries for MPAS-Model
]])

whatis([===[Loads libraries for MPAS-Model ]===])
prepend_path("MODULEPATH", "/ncrc/proj/epic/spack-stack/c6/spack-stack-1.9.3/envs/ue-oneapi-2024.2.1/install/modulefiles/Core")

load("stack-oneapi/2024.2.1")
load("stack-cray-mpich/8.1.32")
load("cmake/3.27.9")
load("parallel-netcdf/1.12.3")
load("parallelio/2.6.2")

if mode() == "load" then
  setenv("PNETCDF", os.getenv("parallel_netcdf_ROOT"))
end
if mode() == "unload" then
  unsetenv("PNETCDF")
end

setenv("CMAKE_C_COMPILER", "mpicc")
setenv("CMAKE_CXX_COMPILER", "mpic++")
setenv("CMAKE_Fortran_COMPILER", "mpifort")
