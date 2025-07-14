!>\file  mpas_smoke_config.F90
!! This file contains the configuration for MPAS-Smoke.
!
! Haiqin.Li@noaa.gov  
! 10/2024
!
module mpas_smoke_config

  use mpas_kind_types
  
  implicit none

  !-- constant paramters
  real(RKIND), parameter :: epsilc     = 1.e-12
  real(RKIND), parameter :: pi         = 3.1415926
  !-- aerosol module configurations
  integer :: seas_opt = 1  ! GOCART scheme by default (2 = NGAC) 
  integer :: addsmoke_flag = 1
  integer :: n_dbg_lines = 3
  integer :: nfire_types = 5
  logical :: dbg_opt       = .false.
  logical :: aero_ind_fdb  = .false.
  ! --
  integer, parameter :: CHEM_OPT_GOCART= 1
  integer, parameter :: num_moist=2
  integer, parameter :: num_emis_seas=5, num_emis_dust=5   
  integer, parameter :: p_dust_1=10, p_dust_2=11, p_dust_3=12, p_dust_4=13, p_dust_5=14
  integer, parameter :: p_edust1=1,p_edust2=2,p_edust3=3,p_edust4=4,p_edust5=5
  integer, parameter :: ndvel = 20
  integer :: numgas = 0

end module mpas_smoke_config
