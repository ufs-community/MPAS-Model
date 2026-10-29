module mpas_gfs_base_pressure_option
 implicit none
 private
 public :: gfs_fixed_base_mode
 logical, save :: initialized=.false., fixed=.false.
contains
 logical function gfs_fixed_base_mode()
   character(len=32) :: value
   integer :: status
   !$omp critical(gfs_pressure_mode_init)
   if (.not.initialized) then
     call get_environment_variable('MPAS_GFS_PRESSURE_MODE',value,status=status)
     if (status==0) then
       select case(trim(adjustl(value)))
       case('0','original'); fixed=.false.
       case('1','fixed_base'); fixed=.true.
       case default
         error stop 'MPAS_GFS_PRESSURE_MODE must be 0 or 1'
       end select
     endif
     initialized=.true.
   endif
   gfs_fixed_base_mode=fixed
   !$omp end critical(gfs_pressure_mode_init)
 end function
end module
