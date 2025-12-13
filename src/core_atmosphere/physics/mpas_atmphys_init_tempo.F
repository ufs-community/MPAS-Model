! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
#define DM_BCAST_MACRO(A) call mpas_dmpar_bcast_real4s(dminfo,size(A),A)

!=================================================================================================================
 module mpas_atmphys_init_tempo
 use mpas_dmpar
 use mpas_kind_types
 use mpas_log
 use mpas_pool_routines

 use mpas_atmphys_utilities
 use module_mp_tempo_params, only: nc_ml_input, nc_ml_nodes, nc_ml_output, &
      nr_ml_input, nr_ml_nodes, nr_ml_output, &
      nc_ml_trans_mean, nc_ml_trans_var, &
      nc_ml_w00, nc_ml_w01, nc_ml_b00, nc_ml_b01
 use module_mp_tempo_ml, only: MLdata, tempo_save_or_read_ml_data

 implicit none
 private
 public:: init_tempo_aerosols_forMPAS, tempo_ml_init

!MPAS main initialization of the TEMPO parameterization of cloud microphysics with nucleation of cloud
!droplets based on distributions of CCNs and INs (aerosol-aware parameterization).

 contains

!=================================================================================================================
 subroutine tempo_ml_init(configs)
! Called once to initial data for tempo_ml
!=================================================================================================================

!input arguments:
type(mpas_pool_type),intent(in):: configs

!local variables and pointers:
 type(MLdata), dimension(2) :: tempo_ml_data
 logical,pointer:: config_tempo_ml_nc_pbl,config_tempo_ml_nc,config_tempo_ml_nr

 call mpas_pool_get_config(configs, 'config_tempo_ml_nc_pbl', config_tempo_ml_nc_pbl)
 call mpas_pool_get_config(configs, 'config_tempo_ml_nc', config_tempo_ml_nc)
 call mpas_pool_get_config(configs, 'config_tempo_ml_nr', config_tempo_ml_nr)

 if(.not. (config_tempo_ml_nc_pbl .or. config_tempo_ml_nc .or. config_tempo_ml_nr)) then
    call mpas_log_write('--- All configuration flags for TEMPO ML are false... TEMPO ML will not be used')
    return
 endif

 if(config_tempo_ml_nc .or. config_tempo_ml_nr) then
    call mpas_log_write('--- TEMPO ML for nc and nr prediction not yet working... These flags will be ignored')
 endif

 if(config_tempo_ml_nc_pbl) then
    call mpas_log_write('--- Using TEMPO ML prediction to give life to the PBL clouds')
 endif

 ! Cloud water
 tempo_ml_data(1)%input_size = nc_ml_input
 tempo_ml_data(1)%node_size = nc_ml_nodes
 tempo_ml_data(1)%output_size = nc_ml_output

 if (.not.allocated(tempo_ml_data(1)%transform_mean)) allocate(tempo_ml_data(1)%transform_mean(nc_ml_input))
 if (.not.allocated(tempo_ml_data(1)%transform_var)) allocate(tempo_ml_data(1)%transform_var(nc_ml_input))

 tempo_ml_data(1)%transform_mean = nc_ml_trans_mean
 tempo_ml_data(1)%transform_var = nc_ml_trans_var

 if (.not.allocated(tempo_ml_data(1)%weights00)) allocate(tempo_ml_data(1)%weights00(nc_ml_nodes,nc_ml_input))
 if (.not.allocated(tempo_ml_data(1)%weights01)) allocate(tempo_ml_data(1)%weights01(nc_ml_output,nc_ml_nodes))
 if (.not.allocated(tempo_ml_data(1)%bias00)) allocate(tempo_ml_data(1)%bias00(nc_ml_nodes))
 if (.not.allocated(tempo_ml_data(1)%bias01)) allocate(tempo_ml_data(1)%bias01(nc_ml_output))

 tempo_ml_data(1)%weights00 = reshape(nc_ml_w00, (/nc_ml_nodes, nc_ml_input/))
 tempo_ml_data(1)%weights01 = reshape(nc_ml_w01, (/nc_ml_output, nc_ml_nodes/))
 tempo_ml_data(1)%bias00 = nc_ml_b00
 tempo_ml_data(1)%bias01 = nc_ml_b01

 ! Rain water
 tempo_ml_data(2)%input_size = nc_ml_input
 tempo_ml_data(2)%node_size = nc_ml_nodes
 tempo_ml_data(2)%output_size = nc_ml_output

 if (.not.allocated(tempo_ml_data(2)%transform_mean)) allocate(tempo_ml_data(2)%transform_mean(nc_ml_input))
 if (.not.allocated(tempo_ml_data(2)%transform_var)) allocate(tempo_ml_data(2)%transform_var(nc_ml_input))

 tempo_ml_data(2)%transform_mean = nc_ml_trans_mean
 tempo_ml_data(2)%transform_var = nc_ml_trans_var

 if (.not.allocated(tempo_ml_data(2)%weights00)) allocate(tempo_ml_data(2)%weights00(nc_ml_nodes,nc_ml_input))
 if (.not.allocated(tempo_ml_data(2)%weights01)) allocate(tempo_ml_data(2)%weights01(nc_ml_output,nc_ml_nodes))
 if (.not.allocated(tempo_ml_data(2)%bias00)) allocate(tempo_ml_data(2)%bias00(nc_ml_nodes))
 if (.not.allocated(tempo_ml_data(2)%bias01)) allocate(tempo_ml_data(2)%bias01(nc_ml_output))

 tempo_ml_data(2)%weights00 = reshape(nc_ml_w00, (/nc_ml_nodes, nc_ml_input/))
 tempo_ml_data(2)%weights01 = reshape(nc_ml_w01, (/nc_ml_output, nc_ml_nodes/))
 tempo_ml_data(2)%bias00 = nc_ml_b00
 tempo_ml_data(2)%bias01 = nc_ml_b01

 ! Save neural network
 call tempo_save_or_read_ml_data(ml_data_in=tempo_ml_data)

 end subroutine tempo_ml_init

!=================================================================================================================
 subroutine init_tempo_aerosols_forMPAS(do_restart,dminfo,mesh,state,time_lev,diag_physics)
!=================================================================================================================

!input variables:
 type(dm_info),intent(in):: dminfo
 type(mpas_pool_type),intent(in):: mesh
 logical,intent(in):: do_restart
 integer,intent(in):: time_lev

!inout variables:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: state

!local variables and pointers:
 integer,pointer:: nCellsSolve,nVertLevels
 integer,pointer:: index_nifa,index_nwfa

 real(kind=RKIND),dimension(:),pointer    :: nifa2d,nwfa2d
 real(kind=RKIND),dimension(:,:),pointer  :: nifa,nwfa
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

 character(len=StrKIND):: mess

 integer:: iCell, k

 real(kind=RKIND):: max_test
 real(kind=RKIND):: nifa_max,nifa_min,global_nifa_max,global_nifa_min
 real(kind=RKIND):: nwfa_max,nwfa_min,global_nwfa_max,global_nwfa_min

!-----------------------------------------------------------------------------------------------------------------
!... if do_restart is true, then we do not need to check the initialization of nwfa, nifa, and nwfa2d. If false,
!    then, we proceed with the initialization:
 if(do_restart) return

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_array(diag_physics,'nifa2d',nifa2d)
 call mpas_pool_get_array(diag_physics,'nwfa2d',nwfa2d)

 call mpas_pool_get_dimension(state,'index_nifa'  ,index_nifa  )
 call mpas_pool_get_dimension(state,'index_nwfa'  ,index_nwfa  )

 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 nifa   => scalars(index_nifa,:,:)
 nwfa   => scalars(index_nwfa,:,:)

!... initialize the distribution of hygroscopic ("water friendly") aerosols if not already initialized using
!    GOCART data:
 global_nwfa_min = 0._RKIND
 global_nwfa_max = 0._RKIND
 nwfa_min = minval(nwfa(:,1:nCellsSolve))
 nwfa_max = maxval(nwfa(:,1:nCellsSolve))
 call mpas_dmpar_min_real(dminfo,nwfa_min,global_nwfa_min)
 call mpas_dmpar_max_real(dminfo,nwfa_max,global_nwfa_max)
 call mpas_log_write('--- global_nwfa_min = $r',realArgs=(/global_nwfa_min/))
 call mpas_log_write('--- global_nwfa_max = $r',realArgs=(/global_nwfa_max/))

 if(global_nwfa_min /= 0._RKIND .and. global_nwfa_max /= 0._RKIND) then
   call mpas_log_write('--- tempo() info: initialized nwfa using the climatological GOCART data.')
 endif

!... initialize the distribution of nonhygroscopic ("ice friendly") aerosols if not already initialized using
!    GOCART data:
 global_nifa_min = 0._RKIND
 global_nifa_max = 0._RKIND
 nifa_min = minval(nifa(:,1:nCellsSolve))
 nifa_max = maxval(nifa(:,1:nCellsSolve))
 call mpas_dmpar_min_real(dminfo,nifa_min,global_nifa_min)
 call mpas_dmpar_max_real(dminfo,nifa_max,global_nifa_max)
 call mpas_log_write('--- global_nifa_min = $r',realArgs=(/global_nifa_min/))
 call mpas_log_write('--- global_nifa_max = $r',realArgs=(/global_nifa_max/))

 if(global_nifa_min == 0._RKIND .and. global_nifa_max == 0._RKIND) then
   call mpas_log_write('--- initialize nifa using the climatological GOCART data.')
 endif

 end subroutine init_tempo_aerosols_forMPAS

!=================================================================================================================
 end module mpas_atmphys_init_tempo
!=================================================================================================================
