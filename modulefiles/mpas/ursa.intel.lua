help([[
This module loads libraries for MPAS-Model
]])

prepend_path("MODULEPATH",
  "/contrib/spack-stack/spack-stack-1.9.3/envs/ue-oneapi-2024.2.1/install/modulefiles/Core")

load("stack-oneapi/2024.2.1")
load("stack-intel-oneapi-mpi/2021.13")
load("intel-oneapi-mkl/2024.2.1")
load("cmake/3.27.9")
load("parallel-netcdf/1.12.3")
load("parallelio/2.6.2")

setenv("CMAKE_C_COMPILER", "mpiicc")
setenv("CMAKE_CXX_COMPILER", "mpiicpc")
setenv("CMAKE_Fortran_COMPILER", "mpiifort")

if mode() == "load" then

  local pnetcdf_root = os.getenv("parallel_netcdf_ROOT")
  local netcdf_c_root = os.getenv("netcdf_c_ROOT")
  local netcdf_fortran_root = os.getenv("netcdf_fortran_ROOT")

  if pnetcdf_root then
    setenv("PNETCDF", pnetcdf_root)
  end

  if netcdf_c_root then
    setenv("NetCDF_C_ROOT", netcdf_c_root)
  end

  if netcdf_fortran_root then
    setenv("NetCDF_FORTRAN_ROOT", netcdf_fortran_root)
  end

end

if mode() == "unload" then
  unsetenv("PNETCDF")
  unsetenv("NetCDF_C_ROOT")
  unsetenv("NetCDF_FORTRAN_ROOT")
end
