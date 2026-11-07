module mpas_gfs_base_pressure_option
 implicit none
 private
 public :: gfs_pressure_mode, gfs_fixed_base_mode
 integer, save :: mode = 0
 logical, save :: initialized = .false.
contains
 integer function gfs_pressure_mode()
   character(len=32) :: value
   integer :: status
   !$omp critical(gfs_pressure_mode_init)
   if (.not.initialized) then
     value = ''
     call get_environment_variable('MPAS_GFS_PRESSURE_MODE', value, status=status)
     if (status == 0) then
       select case (trim(adjustl(value)))
       case ('0','original')
         mode = 0
       case ('1','fixed_base')
         mode = 1
       case ('2','mpas_full')
         mode = 2
       case default
         error stop 'MPAS_GFS_PRESSURE_MODE must be 0, 1, or 2'
       end select
     else if (status /= 1) then
       error stop 'MPAS_GFS_PRESSURE_MODE environment read failed'
     endif
     initialized = .true.
   endif
   gfs_pressure_mode = mode
   !$omp end critical(gfs_pressure_mode_init)
 end function gfs_pressure_mode

 logical function gfs_fixed_base_mode()
   gfs_fixed_base_mode = (gfs_pressure_mode() == 1)
 end function gfs_fixed_base_mode
end module mpas_gfs_base_pressure_option
