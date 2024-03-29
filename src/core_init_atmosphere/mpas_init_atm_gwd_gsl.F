! Module mpas_init_atm_gwd_gsl
!
! This program calls subroutines which calculate the parameters
! required for the GSL subgrid-scale orographic gravity-wave drag (GWDO)
! suite on the MPAS mesh.  These parameters are for the small-scale
! GWD (Tsiringakis et al., 2017) and turbulent orographic form drag (TOFD)
! (Beljaars et al., 2004) schemes of the GSL drag suite.
! The output fields are:
! - stddev      standard deviation of subgrid-scale topograpy
! - convexity   convexity (kurtosis) of subgrid-scale topography
! - ol{1,2,3,4} orographic effective lengths of subgrid-scale topography
!   for 4 orientations: 1-westerly, 2-southerly, 3-southwesterly, 4-northwesterly
! - oa{1,2,3,4} orographic asymmetries of subgrid-scale topography
!   for 4 orientations: 1-westerly, 2-southerly, 3-southwesterly, 4-northwesterly
!
! Based on code by Michael Duda provided by NCAR/MMM
!
! Brief description of program:  Creates orographic (oro_data) files
! needed by the GSL drag suite physics parameterization 
!
! Author: Michael Toy, NOAA/GSL
!
module mpas_init_atm_gwd_gsl

! Temporary comment
! use gsl_oro_data_sm_scale, only: calc_gsl_oro_data_sm_scale
! use gsl_oro_data_lg_scale, only: calc_gsl_oro_data_lg_scale

   use iso_c_binding, only : c_char, c_int, c_float, c_ptr, c_loc

   use mpas_derived_types, only : MPAS_LOG_ERR
   use mpas_framework
   use mpas_timekeeping
   use mpas_log, only : mpas_log_write
   use mpas_c_interfacing, only : mpas_f_to_c_string

   public :: calc_gsl_oro_data

   private

   integer, parameter :: I1KIND = selected_int_kind(2)


   contains


   subroutine calc_gsl_oro_data(domain,iErr)

      use mpas_derived_types
      use mpas_kind_types
      use mpas_timer
      use mpas_stream_manager

      implicit none

      type (domain_type), intent(inout) :: domain
      integer, intent(inout) :: iErr

      type (mpas_pool_type), pointer :: mesh, state
      integer :: iCell
      real (kind=RKIND) :: cellArea
      real (kind=RKIND), pointer :: config_gwd_cell_scaling
      integer, pointer :: nCells
      integer, pointer :: nEdges
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell
      logical :: onUnitSphere
      real (kind=RKIND), pointer :: sphere_radius
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell, dcEdge
      character(len=StrKIND), pointer :: config_geog_data_path
      character(len=StrKIND), pointer :: config_topo_data
      character(len=StrKIND) :: geog_sub_path
      character(len=StrKIND+1) :: geog_data_path      ! same as config_geog_data_path, but guaranteed to have a trailing slash

      integer (kind=I1KIND), dimension(:), pointer :: hlanduse ! Dominant land mask (0 or 1)

      logical :: duplicate_oro_data_file   ! flag for whether oro_data_ls file is a duplicate
                                           ! of oro_data_ss due to minimum grid size being less than 7.5km

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', state)

      call mpas_log_write('In sub calc_gsl_oro_data #1')

      call mpas_pool_get_config(mesh, 'sphere_radius', sphere_radius)
      call mpas_log_write('In sub calc_gsl_oro_data #2')
      call mpas_pool_get_config(domain % configs, 'config_geog_data_path', config_geog_data_path)
      call mpas_log_write('In sub calc_gsl_oro_data #3')
      call mpas_pool_get_config(domain % configs, 'config_topo_data', config_topo_data)
      call mpas_log_write('In sub calc_gsl_oro_data #4')
      call mpas_pool_get_config(domain % configs, 'config_gwd_cell_scaling', config_gwd_cell_scaling)

      call mpas_log_write('In sub calc_gsl_oro_data #5')

      !
      ! Retrieve pointers to arrays holding the latitudes and longitudes of cells
      !
      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_array(mesh, 'latCell', latCell)
      call mpas_pool_get_array(mesh, 'lonCell', lonCell)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)

      allocate(hlanduse(nCells+1))    ! +1, since we access hlanduse(cellsOnCell(i,iCell)) later on for iCell=1,nCells

      call mpas_log_write('In sub calc_gsl_oro_data #6')

      !
      ! It is possible that this code is called before the mesh fields have been scaled
      ! up to "Earth-sized". Because we need "Earth" distances to cut out bounding
      ! boxes from topography, we try here to detect whether we are on an unscaled
      ! unit sphere or not: if the maximum dcEdge value is less than 1.0, assume this
      ! is the case.
      !
      if (maxval(dcEdge(1:nEdges)) < 1.0_RKIND) then
         call mpas_log_write('Computing GWD statistics on a unit sphere')
         onUnitSphere = .true.
      else
         onUnitSphere = .false.
      end if

      if (config_gwd_cell_scaling /= 1.0) then
         call mpas_log_write('Using effective cell diameters scaled by a factor of $r', realArgs=(/config_gwd_cell_scaling/))
         call mpas_log_write('in the computation of GWD static fields.')
      end if

      call mpas_log_write('In sub calc_gsl_oro_data #7')


      !
      ! Loop to compute approximate area of each MPAS horizontal grid cell
      !
      do iCell=1,nCells

         print *, "iCell =", iCell

      end do


      ! Temporary comment
      ! call calc_gsl_oro_data_sm_scale(,duplicate_oro_data_file)


      print *, "duplicate_oro_data_file =", duplicate_oro_data_file
      print *

      ! Temporary comment
      ! if ( .not.duplicate_oro_data_file ) then
      !    call calc_gsl_oro_data_lg_scale(tile_num,res_indx,halo)
      ! end if


      call mpas_log_write('')
      call mpas_log_write('End module gsl_oro_data')
      call mpas_log_write('')

      iErr = 0

      return
   end subroutine calc_gsl_oro_data

end module mpas_init_atm_gwd_gsl
