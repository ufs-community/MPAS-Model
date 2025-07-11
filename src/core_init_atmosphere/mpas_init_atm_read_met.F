! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
#define _METGRID

module init_atm_read_met

   use iso_fortran_env, only: real32

   integer, parameter :: MAX_FILENAME_LEN = 1024

   real (kind=real32), parameter :: EARTH_RADIUS_M = 6370000.   ! same as MM5 system

   ! Projection codes for proj_info structure:
   INTEGER, PRIVATE, PARAMETER  :: PROJ_LATLON = 0
   INTEGER, PRIVATE, PARAMETER  :: PROJ_LC = 1
   INTEGER, PRIVATE, PARAMETER  :: PROJ_PS = 2
   INTEGER, PRIVATE, PARAMETER  :: PROJ_MERC = 3
   INTEGER, PRIVATE, PARAMETER  :: PROJ_GAUSS = 4
   INTEGER, PRIVATE, PARAMETER  :: PROJ_CASSINI = 6

   ! Derived types
   type met_data
      integer                       :: version, nx, ny, iproj
      real (kind=real32)                 :: xfcst, xlvl, startlat, startlon, starti, startj, &
                                       deltalat, deltalon, dx, dy, xlonc, &
                                       truelat1, truelat2, earth_radius
      real (kind=real32)                 :: centerlat, centerlon, pole_lat, pole_lon
      real (kind=real32), pointer, dimension(:,:) :: slab
      logical                       :: is_wind_grid_rel
      character (len=9)             :: field
      character (len=24)            :: hdate
      character (len=25)            :: units
      character (len=32)            :: map_source
      character (len=46)            :: desc
   end type met_data


   ! State variables?
   integer :: input_unit
   character (len=MAX_FILENAME_LEN) :: filename
 
   contains
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: read_met_init
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine read_met_init(fg_source, source_is_constant, datestr, istatus)

      use mpas_derived_types, only : MPAS_LOG_ERR
      use mpas_log, only : mpas_log_write
      use mpas_io_units, only : mpas_new_unit
 
      implicit none
  
      ! Arguments
      integer, intent(out) :: istatus
      logical, intent(in) :: source_is_constant
      character (len=*), intent(in) :: fg_source
      character (len=*), intent(in) :: datestr
  
      ! Local variables
      integer :: io_status
      logical :: is_used

      istatus = 0
    
      !  1) BUILD FILENAME BASED ON TIME 
      filename = ' '
      if (.not. source_is_constant) then 
         write(filename, '(a)') trim(fg_source)//':'//trim(datestr)
      else
         write(filename, '(a)') trim(fg_source)
      end if
  
      !  2) OPEN FILE
      call mpas_new_unit(input_unit, unformatted = .true.)
      if (input_unit < 0) call mpas_log_write('In read_met_init(), couldn''t find an available Fortran unit.', messageType=MPAS_LOG_ERR)
      open(unit=input_unit, file=trim(filename), status='old', form='unformatted', iostat=io_status)

      if (io_status > 0) istatus = 1

      return
  
 
   end subroutine read_met_init
 
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: read_next_met_field
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine read_next_met_field(fg_data, istatus)

      use mpas_derived_types, only : MPAS_LOG_ERR
      use mpas_log, only : mpas_log_write
 
      implicit none
  
      ! Arguments
      type (met_data), intent(inout) :: fg_data
      integer, intent(out) :: istatus
  
      ! Local variables
      character (len=8) :: startloc
  
      istatus = 1
  
      !  1) READ FORMAT VERSION
      read(unit=input_unit,err=1001,end=1001) fg_data % version
  
      ! PREGRID
      if (fg_data % version == 3) then

         read(unit=input_unit) fg_data % hdate, &
                               fg_data % xfcst, &
                               fg_data % field, &
                               fg_data % units, &
                               fg_data % desc,  &
                               fg_data % xlvl,  &
                               fg_data % nx,    &
                               fg_data % ny,    &
                               fg_data % iproj

         fg_data % map_source = ' '

         if (fg_data % field == 'HGT      ') fg_data % field = 'GHT      '

         fg_data % starti = 1.0
         fg_data % startj = 1.0
     
         ! Cylindrical equidistant
         if (fg_data % iproj == 0) then
            fg_data % iproj = PROJ_LATLON
            read(unit=input_unit,err=1001,end=1001) fg_data % startlat, &
                                                    fg_data % startlon, &
                                                    fg_data % deltalat, &
                                                    fg_data % deltalon
     
         ! Mercator
         else if (fg_data % iproj == 1) then
            fg_data % iproj = PROJ_MERC
            read(unit=input_unit,err=1001,end=1001) fg_data % startlat, &
                                                    fg_data % startlon, &
                                                    fg_data % dx,       &
                                                    fg_data % dy,       &
                                                    fg_data % truelat1
     
         ! Lambert conformal
         else if (fg_data % iproj == 3) then
            fg_data % iproj = PROJ_LC
            read(unit=input_unit,err=1001,end=1001) fg_data % startlat, &
                                                    fg_data % startlon, &
                                                    fg_data % dx,       &
                                                    fg_data % dy,       &
                                                    fg_data % xlonc,    &
                                                    fg_data % truelat1, &
                                                    fg_data % truelat2
     
         ! Polar stereographic
         else if (fg_data % iproj == 5) then
            fg_data % iproj = PROJ_PS
            read(unit=input_unit,err=1001,end=1001) fg_data % startlat, &
                                                    fg_data % startlon, &
                                                    fg_data % dx,       &
                                                    fg_data % dy,       &
                                                    fg_data % xlonc,    &
                                                    fg_data % truelat1

            call mpas_log_write('Unrecognized projection code $i when reading from '//trim(filename), intArgs=(/fg_data % iproj/), messageType=MPAS_LOG_ERR)
     
         end if
     
         fg_data % earth_radius = EARTH_RADIUS_M / 1000.

#if (defined _GEOGRID) || (defined _METGRID)
         fg_data % dx = fg_data % dx * 1000.
         fg_data % dy = fg_data % dy * 1000.

         if (fg_data % xlonc    > 180.) fg_data % xlonc    = fg_data%xlonc    - 360.

         if (fg_data % startlon > 180.) fg_data % startlon = fg_data%startlon - 360.
  
         if (fg_data % startlat < -90.) fg_data % startlat = -90.
         if (fg_data % startlat >  90.) fg_data % startlat = 90.
#endif
     
         fg_data % is_wind_grid_rel = .true.
     
         allocate(fg_data % slab(fg_data % nx, fg_data % ny))
         read(unit=input_unit,err=1001,end=1001) fg_data % slab
     
         istatus = 0 
    
      ! GRIB_PREP
      else if (fg_data % version == 4) then
  
         read(unit=input_unit) fg_data % hdate,      &
                               fg_data % xfcst,      &
                               fg_data % map_source, &
                               fg_data % field,      &
                               fg_data % units,      &
                               fg_data % desc,       &
                               fg_data % xlvl,       &
                               fg_data % nx,         &
                               fg_data % ny,         &
                               fg_data % iproj
  
         if (fg_data % field == 'HGT      ') fg_data % field = 'GHT      '
  
         ! Cylindrical equidistant
         if (fg_data % iproj == 0) then
            fg_data % iproj = PROJ_LATLON
            read(unit=input_unit,err=1001,end=1001) startloc, &
                                                    fg_data % startlat, &
                                                    fg_data % startlon, &
                                                    fg_data % deltalat, &
                                                    fg_data % deltalon

         ! Mercator
         else if (fg_data % iproj == 1) then
            fg_data % iproj = PROJ_MERC
            read(unit=input_unit,err=1001,end=1001) startloc, &
                                                    fg_data % startlat, &
                                                    fg_data % startlon, &
                                                    fg_data % dx,       &
                                                    fg_data % dy,       &
                                                    fg_data % truelat1

         ! Lambert conformal
         else if (fg_data % iproj == 3) then
            fg_data % iproj = PROJ_LC
            read(unit=input_unit,err=1001,end=1001) startloc, &
                                                    fg_data % startlat, &
                                                    fg_data % startlon, &
                                                    fg_data % dx,       &
                                                    fg_data % dy,       &
                                                    fg_data % xlonc,    &
                                                    fg_data % truelat1, &
                                                    fg_data % truelat2

         ! Polar stereographic
         else if (fg_data % iproj == 5) then
            fg_data % iproj = PROJ_PS
            read(unit=input_unit,err=1001,end=1001) startloc, &
                                                    fg_data % startlat, &
                                                    fg_data % startlon, &
                                                    fg_data % dx,       &
                                                    fg_data % dy,       &
                                                    fg_data % xlonc,    &
                                                    fg_data % truelat1
     
         ! ?????????
         else
            call mpas_log_write('Unrecognized projection code $i when reading from '//trim(filename), intArgs=(/fg_data % iproj/), messageType=MPAS_LOG_ERR)
     
         end if
  
         if (startloc == 'CENTER  ') then
            fg_data % starti = real(fg_data % nx)/2.
            fg_data % startj = real(fg_data % ny)/2.
         else if (startloc == 'SWCORNER') then
            fg_data % starti = 1.0
            fg_data % startj = 1.0
         end if

         fg_data % earth_radius = EARTH_RADIUS_M / 1000.

#if (defined _GEOGRID) || (defined _METGRID)
         fg_data % dx = fg_data % dx * 1000.
         fg_data % dy = fg_data % dy * 1000.

         if (fg_data % xlonc    > 180.) fg_data % xlonc    = fg_data % xlonc    - 360.

         if (fg_data % startlon > 180.) fg_data % startlon = fg_data % startlon - 360.
  
         if (fg_data % startlat < -90.) fg_data % startlat = -90.
         if (fg_data % startlat >  90.) fg_data % startlat = 90.
#endif
         
         fg_data % is_wind_grid_rel = .true.
      
         allocate(fg_data % slab(fg_data % nx, fg_data % ny))
         read(unit=input_unit,err=1001,end=1001) fg_data % slab
      
         istatus = 0

      ! WPS
      else if (fg_data % version == 5) then
  
         read(unit=input_unit) fg_data % hdate,      &
                               fg_data % xfcst,      &
                               fg_data % map_source, &
                               fg_data % field,      &
                               fg_data % units,      &
                               fg_data % desc,       &
                               fg_data % xlvl,       &
                               fg_data % nx,         &
                               fg_data % ny,         &
                               fg_data % iproj
  
         if (fg_data % field == 'HGT      ') fg_data % field = 'GHT      '
  
         ! Cylindrical equidistant
         if (fg_data % iproj == 0) then
            fg_data % iproj = PROJ_LATLON
            read(unit=input_unit,err=1001,end=1001) startloc, &
                                                    fg_data % startlat, &
                                                    fg_data % startlon, &
                                                    fg_data % deltalat, &
                                                    fg_data % deltalon, &
                                                    fg_data % earth_radius
            fg_data % dx = 0.0
            fg_data % dy = 0.0

         ! Mercator
         else if (fg_data % iproj == 1) then
            fg_data % iproj = PROJ_MERC
            read(unit=input_unit,err=1001,end=1001) startloc, &
                                                    fg_data % startlat, &
                                                    fg_data % startlon, &
                                                    fg_data % dx,       &
                                                    fg_data % dy,       &
                                                    fg_data % truelat1, &
                                                    fg_data % earth_radius

         ! Lambert conformal
         else if (fg_data % iproj == 3) then
            fg_data % iproj = PROJ_LC
            read(unit=input_unit,err=1001,end=1001) startloc, &
                                                    fg_data % startlat, &
                                                    fg_data % startlon, &
                                                    fg_data % dx,       &
                                                    fg_data % dy,       &
                                                    fg_data % xlonc,    &
                                                    fg_data % truelat1, &
                                                    fg_data % truelat2, &
                                                    fg_data % earth_radius

         ! Gaussian
         else if (fg_data % iproj == 4) then
            fg_data % iproj = PROJ_GAUSS
            read(unit=input_unit,err=1001,end=1001) startloc, &
                                                    fg_data % startlat, &
                                                    fg_data % startlon, &
                                                    fg_data % deltalat, &
                                                    fg_data % deltalon, &
                                                    fg_data % earth_radius
            fg_data % dx = 0.0
            fg_data % dy = 0.0

         ! Polar stereographic
         else if (fg_data % iproj == 5) then
            fg_data % iproj = PROJ_PS
            read(unit=input_unit,err=1001,end=1001) startloc, &
                                                    fg_data % startlat, &
                                                    fg_data % startlon, &
                                                    fg_data % dx,       &
                                                    fg_data % dy,       &
                                                    fg_data % xlonc,    &
                                                    fg_data % truelat1, &
                                                    fg_data % earth_radius
         ! CASSINI
         else if (fg_data % iproj == 6) then
            fg_data % iproj = PROJ_CASSINI
            read(unit=input_unit,err=1001,end=1001) startloc, &
                                                    fg_data % startlat, &
                                                    fg_data % startlon, &
                                                    fg_data % dx,       &
                                                    fg_data % dy,       &
                                                    fg_data % centerlat, &
                                                    fg_data % centerlon, &
                                                    fg_data % earth_radius
            if ( fg_data % centerlat > 0. ) then
                fg_data % pole_lat = 90. - fg_data % centerlat
                fg_data % pole_lon = 180.
                fg_data % xlonc = -fg_data % centerlon
            else
                fg_data % pole_lat = 90. + fg_data % centerlat
                fg_data % pole_lon = 0.
                fg_data % xlonc = 180. - fg_data % centerlon
            end if
            fg_data % deltalon = fg_data % dx
            fg_data % deltalat = fg_data % dy
     
         ! ?????????
         else
            call mpas_log_write('Unrecognized projection code $i when reading from '//trim(filename), intArgs=(/fg_data % iproj/), messageType=MPAS_LOG_ERR)
     
         end if
  
         if (startloc == 'CENTER  ') then
            fg_data % starti = real(fg_data % nx)/2.
            fg_data % startj = real(fg_data % ny)/2.
         else if (startloc == 'SWCORNER') then
            fg_data % starti = 1.0
            fg_data % startj = 1.0
         end if

#if (defined _GEOGRID) || (defined _METGRID)
         fg_data % dx = fg_data % dx * 1000.
         fg_data % dy = fg_data % dy * 1000.

         if (fg_data % xlonc    > 180.) fg_data % xlonc    = fg_data % xlonc    - 360.

         if (fg_data % startlon > 180.) fg_data % startlon = fg_data % startlon - 360.
         
         if (fg_data % startlat < -90.) fg_data % startlat = -90.
         if (fg_data % startlat >  90.) fg_data % startlat =  90.
#endif
 
         read(unit=input_unit,err=1001,end=1001) fg_data % is_wind_grid_rel
      
         allocate(fg_data % slab(fg_data % nx, fg_data % ny))
         read(unit=input_unit,err=1001,end=1001) fg_data % slab
      
         istatus = 0

      else
         call mpas_log_write('Didn''t recognize format version of data in '//trim(filename)//'.', messageType=MPAS_LOG_ERR)
         call mpas_log_write('Found version $i but expected either 3, 4, or 5. This could be an endian problem.', intArgs=(/fg_data % version/), messageType=MPAS_LOG_ERR)
      end if
  
      return
 
   1001 return
 
   end subroutine read_next_met_field
 
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: read_met_close
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine read_met_close()
 
      use mpas_io_units, only : mpas_release_unit

      implicit none
  
      close(unit=input_unit)
      call mpas_release_unit(input_unit)
      filename = 'UNINITIALIZED_FILENAME'
  
   end subroutine read_met_close

end module init_atm_read_met
