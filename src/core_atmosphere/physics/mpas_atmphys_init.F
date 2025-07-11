! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_init
 use mpas_kind_types
 use mpas_pool_routines
 use mpas_timekeeping

 use mpas_atmphys_driver_convection,only: init_convection
 use mpas_atmphys_driver_lsm,only: init_lsm
 use mpas_atmphys_driver_microphysics,only: init_microphysics
 use mpas_atmphys_driver_pbl,only: init_pbl
 use mpas_atmphys_driver_radiation_lw,only: init_radiation_lw
 use mpas_atmphys_driver_radiation_sw,only: init_radiation_sw
 use mpas_atmphys_driver_sfclayer,only: init_sfclayer
 use mpas_atmphys_vars,only: f_qc,f_qr,f_qi,f_qs,f_qg,f_qoz,f_nc,f_ni,f_nifa,f_nwfa,f_nbca

 use mpas_atmphys_landuse
 use mpas_atmphys_o3climatology
 use mpas_atmphys_lsm_noahmpinit,only: init_lsm_noahmp

 use bl_ugwpv1_ngw, only: ugwpv1_ngw_init

 implicit none
 private
 public:: physics_init

!MPAS main initialization subroutine for all physics parameterizations.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines in mpas_atmphys_init:
! ---------------------------------
! physics_init    : call initialization of individual physics parameterizations.
! init_dir_forphys: needed for initialization of "reconstruct" subroutines.
! r3_normalize    : needed for initialization of "reconstruct" subroutines.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * added structure diag in calls to subroutine init_radiation_lw and init_radiation_sw.
!   Laura D. Fowler (laura@ucar.edu) / 2013-07-01.
! * added call to subroutine init_o3climatology. reads monthly-mean climatological ozone data and interpolates
!   ozone data to the MPAS grid.
!   Laura D. Fowler (laura@ucar.edu) / 2013-07-03.
! * added the calculation of the mean distance between cell centers.
!   Laura D. Fowler (laura@ucar.edu) / 2013-08-22.
! * added initialization of variable xicem.
!   Laura D. Fowler (laura@ucar.edu) / 2013-08-24.
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * added initialization of the accumulated surface pressure. Added initialization of the tendency and the
!   accumulated tendency of the surface pressure.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * renamed config_conv_deep_scheme to config_convection_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2014-09-18.
! * changed the argument list in the call to subroutine microphysics_init, needed to include the Thompson
!   parameterization of cloud microphysics.
!   Laura D. Fowler (laura@ucar.edu) / 2015-03-28.
! * modified the initialization of i_rainc and i_rainnc, now that the convection and cloud microphysics
!   parameterizations are in "packages."
!   Laura D. Fowler (laura@ucar.edu) / 2106-04-13.
! * removed the calculation of the variable dcEdge_m which is no longer needed in the different physics
!   parameterizations.
!   Laura D. Fowler (laura@ucar.edu) / 2016-10-18.
! * added the subroutine init_physics_flags to initialize f_qc,f_qr,f_qi,f_qs,f_qg,f_nc,and f_ni.
!   Laura D. Fowler (laura@ucar.edu) / 2024-02-14.
! * added call to subroutine init_lsm_noahmp to initialize the Noah-MP land surface scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2024-03-11.


 contains


!=================================================================================================================
 subroutine physics_init(dminfo,stream_manager,clock,configs,mesh,diag,tend,state,time_lev,diag_physics,     &
                         diag_physics_noahmp,ngw_input,atm_input,sfc_input,output_noahmp)
!=================================================================================================================

use mpas_stream_manager

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(MPAS_streamManager_type),intent(inout):: stream_manager
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs
 type(MPAS_Clock_type),intent(in):: clock

 integer,intent(in):: time_lev

!inout arguments:
 type(mpas_pool_type),intent(inout):: state
 type(mpas_pool_type),intent(inout):: diag
 type(mpas_pool_type),intent(inout):: tend
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: diag_physics_noahmp
 type(mpas_pool_type),intent(inout):: ngw_input
 type(mpas_pool_type),intent(inout):: atm_input
 type(mpas_pool_type),intent(inout):: sfc_input
 type(mpas_pool_type),intent(inout):: output_noahmp

!local pointers:
 logical,pointer:: config_do_restart,    &
                   config_do_DAcycling,  &
                   config_o3climatology, &
                   config_oml1d

 character(len=StrKIND),pointer::            &
                   config_convection_scheme, &
                   config_lsm_scheme,        &
                   config_microp_scheme,     &
                   config_pbl_scheme,        &
                   config_sfclayer_scheme,   &
                   config_radt_lw_scheme,    &
                   config_radt_sw_scheme

 integer,pointer:: nCellsSolve,nLags
 integer,dimension(:),pointer:: i_rainc,i_rainnc
 integer,dimension(:),pointer:: i_acswdnb,i_acswdnbc,i_acswdnt,i_acswdntc, &
                                i_acswupb,i_acswupbc,i_acswupt,i_acswuptc, &
                                i_aclwdnb,i_aclwdnbc,i_aclwdnt,i_aclwdntc, &
                                i_aclwupb,i_aclwupbc,i_aclwupt,i_aclwuptc

 real(kind=RKIND),dimension(:),pointer:: acswdnb,acswdnbc,acswdnt,acswdntc, &
                                         acswupb,acswupbc,acswupt,acswuptc, &
                                         aclwdnb,aclwdnbc,aclwdnt,aclwdntc, &
                                         aclwupb,aclwupbc,aclwupt,aclwuptc
 real(kind=RKIND),dimension(:),pointer:: nsteps_accum,ndays_accum,tday_accum, &
                                         tyear_accum,tyear_mean
 real(kind=RKIND),dimension(:),pointer:: sst,sstsk,tmn,xice,xicem
 real(kind=RKIND),dimension(:,:),pointer:: tlag

 real(kind=RKIND),pointer:: config_oml_hml0
 real(kind=RKIND),dimension(:),pointer:: t_oml,t_oml_initial,t_oml_200m_initial
 real(kind=RKIND),dimension(:),pointer:: h_oml,h_oml_initial,hu_oml,hv_oml
 real(kind=RKIND),dimension(:),pointer:: rdzw,dzu
 real(kind=RKIND),pointer:: config_dt
 integer,pointer:: nCells,nVertLevels
 character(len=StrKIND),pointer:: gwdo_scheme
 logical,pointer:: ngw_scheme
 integer,pointer:: ntau_d1y
 real(kind=RKIND),pointer:: knob_ugwp_tauamp
 real(kind=RKIND),dimension(:),pointer:: ugwp_taulat
 integer,dimension(:),pointer:: jindx1_tau, jindx2_tau
 real(kind=RKIND),dimension(:),pointer:: ddy_j1tau, ddy_j2tau
 real(kind=RKIND),dimension(:),pointer:: latCell

!local variables and arrays:
 type(MPAS_Time_Type):: currTime

 logical:: init_done
 integer:: ierr,julday 
 integer:: iCell,iLag

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('')
!call mpas_log_write('--- enter subroutine physics_init:')

 call mpas_pool_get_config(configs,'config_do_restart'       ,config_do_restart       )
 call mpas_pool_get_config(configs,'config_do_DAcycling'     ,config_do_DAcycling     )
 call mpas_pool_get_config(configs,'config_o3climatology'    ,config_o3climatology    )
 call mpas_pool_get_config(configs,'config_convection_scheme',config_convection_scheme)
 call mpas_pool_get_config(configs,'config_lsm_scheme'       ,config_lsm_scheme       )
 call mpas_pool_get_config(configs,'config_microp_scheme'    ,config_microp_scheme    )
 call mpas_pool_get_config(configs,'config_pbl_scheme'       ,config_pbl_scheme       )
 call mpas_pool_get_config(configs,'config_sfclayer_scheme'  ,config_sfclayer_scheme  )
 call mpas_pool_get_config(configs,'config_radt_lw_scheme'   ,config_radt_lw_scheme   )
 call mpas_pool_get_config(configs,'config_radt_sw_scheme'   ,config_radt_sw_scheme   )

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nLags'      ,nLags      )

 call mpas_pool_get_array(diag_physics,'i_acswdnb'   ,i_acswdnb   )
 call mpas_pool_get_array(diag_physics,'i_acswdnbc'  ,i_acswdnbc  )
 call mpas_pool_get_array(diag_physics,'i_acswdnt'   ,i_acswdnt   )
 call mpas_pool_get_array(diag_physics,'i_acswdntc'  ,i_acswdntc  )
 call mpas_pool_get_array(diag_physics,'i_acswupb'   ,i_acswupb   )
 call mpas_pool_get_array(diag_physics,'i_acswupbc'  ,i_acswupbc  )
 call mpas_pool_get_array(diag_physics,'i_acswupt'   ,i_acswupt   )
 call mpas_pool_get_array(diag_physics,'i_acswuptc'  ,i_acswuptc  )
 call mpas_pool_get_array(diag_physics,'i_aclwdnb'   ,i_aclwdnb   )
 call mpas_pool_get_array(diag_physics,'i_aclwdnbc'  ,i_aclwdnbc  )
 call mpas_pool_get_array(diag_physics,'i_aclwdnt'   ,i_aclwdnt   )
 call mpas_pool_get_array(diag_physics,'i_aclwdntc'  ,i_aclwdntc  )
 call mpas_pool_get_array(diag_physics,'i_aclwupb'   ,i_aclwupb   )
 call mpas_pool_get_array(diag_physics,'i_aclwupbc'  ,i_aclwupbc  )
 call mpas_pool_get_array(diag_physics,'i_aclwupt'   ,i_aclwupt   )
 call mpas_pool_get_array(diag_physics,'i_aclwuptc'  ,i_aclwuptc  )

 call mpas_pool_get_array(diag_physics,'acswdnb'     ,acswdnb     )
 call mpas_pool_get_array(diag_physics,'acswdnbc'    ,acswdnbc    )
 call mpas_pool_get_array(diag_physics,'acswdnt'     ,acswdnt     )
 call mpas_pool_get_array(diag_physics,'acswdntc'    ,acswdntc    )
 call mpas_pool_get_array(diag_physics,'acswupb'     ,acswupb     )
 call mpas_pool_get_array(diag_physics,'acswupbc'    ,acswupbc    )
 call mpas_pool_get_array(diag_physics,'acswupt'     ,acswupt     )
 call mpas_pool_get_array(diag_physics,'acswuptc'    ,acswuptc    )
 call mpas_pool_get_array(diag_physics,'aclwdnb'     ,aclwdnb     )
 call mpas_pool_get_array(diag_physics,'aclwdnbc'    ,aclwdnbc    )
 call mpas_pool_get_array(diag_physics,'aclwdnt'     ,aclwdnt     )
 call mpas_pool_get_array(diag_physics,'aclwdntc'    ,aclwdntc    )
 call mpas_pool_get_array(diag_physics,'aclwupb'     ,aclwupb     )
 call mpas_pool_get_array(diag_physics,'aclwupbc'    ,aclwupbc    )
 call mpas_pool_get_array(diag_physics,'aclwupt'     ,aclwupt     )
 call mpas_pool_get_array(diag_physics,'aclwuptc'    ,aclwuptc    )

 call mpas_pool_get_array(diag_physics,'nsteps_accum',nsteps_accum)
 call mpas_pool_get_array(diag_physics,'ndays_accum' ,ndays_accum )
 call mpas_pool_get_array(diag_physics,'tday_accum'  ,tday_accum  )
 call mpas_pool_get_array(diag_physics,'tyear_accum' ,tyear_accum )
 call mpas_pool_get_array(diag_physics,'tyear_mean'  ,tyear_mean  )
 call mpas_pool_get_array(diag_physics,'tlag'        ,tlag        )
 call mpas_pool_get_array(diag_physics,'sstsk'       ,sstsk       )
 call mpas_pool_get_array(diag_physics,'xicem'       ,xicem       )

 call mpas_pool_get_array(sfc_input,'sst' ,sst )
 call mpas_pool_get_array(sfc_input,'tmn' ,tmn )
 call mpas_pool_get_array(sfc_input,'xice',xice)

 call mpas_pool_get_array(diag_physics,'t_oml'             ,t_oml)
 call mpas_pool_get_array(diag_physics,'t_oml_initial'     ,t_oml_initial)
 call mpas_pool_get_array(diag_physics,'t_oml_200m_initial',t_oml_200m_initial)
 call mpas_pool_get_array(diag_physics,'h_oml'             ,h_oml)
 call mpas_pool_get_array(diag_physics,'h_oml_initial'     ,h_oml_initial)
 call mpas_pool_get_array(diag_physics,'hu_oml'            ,hu_oml)
 call mpas_pool_get_array(diag_physics,'hv_oml'            ,hv_oml)
 call mpas_pool_get_config(configs,'config_oml1d'          ,config_oml1d  )
 call mpas_pool_get_config(configs,'config_oml_hml0'       ,config_oml_hml0  )
 call mpas_pool_get_config(configs,'config_gwdo_scheme'    ,gwdo_scheme )
 call mpas_pool_get_config(configs,'config_ngw_scheme'     ,ngw_scheme  )
 call mpas_pool_get_config(configs,'config_dt'             ,config_dt  )
 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels  )
 call mpas_pool_get_array(mesh,'rdzw'   ,rdzw   )
 call mpas_pool_get_array(mesh,'dzu'    ,dzu    )

 currTime = mpas_get_clock_time(clock,MPAS_NOW,ierr)
 call mpas_get_time(curr_time=currTime,DoY=julday,ierr=ierr)

!initialization of east-north directions to convert u-tendencies from cell centers to cell
!edges:
 call init_dirs_forphys(mesh)

!initialization of logical flags for cloud mixing ratios and number concentrations, and aerosols
!number concentrations from the Thompson cloud microphysics:
 call init_physics_flags(state,f_qc,f_qr,f_qi,f_qs,f_qg,f_qoz,f_nc,f_ni,f_nifa,f_nwfa,f_nbca)

!initialization of counters i_rainc and i_rainnc. i_rainc and i_rainnc track the number of
!times the accumulated convective (rainc) and grid-scale (rainnc) rain exceed the prescribed
!threshold value:
 if(.not. config_do_restart .and. config_convection_scheme.ne.'off') then
    call mpas_pool_get_array(diag_physics,'i_rainc',i_rainc)
    do iCell = 1, nCellsSolve
       i_rainc(iCell)  = 0
    enddo
 endif
 if(.not. config_do_restart .and. config_microp_scheme.ne.'off') then
    call mpas_pool_get_array(diag_physics,'i_rainnc',i_rainnc)
    do iCell = 1, nCellsSolve
       i_rainnc(iCell) = 0
    enddo
 endif

!initialization of counters i_acsw* and i_aclw*. i_acsw* and i_aclw* track the number of times
!the accumulated long and short-wave radiation fluxes exceed their prescribed theshold values.
 if(.not. config_do_restart) then
    do iCell = 1, nCellsSolve
       i_acswdnb(iCell)  = 0
       i_acswdnbc(iCell) = 0
       i_acswdnt(iCell)  = 0
       i_acswdntc(iCell) = 0
       i_acswupb(iCell)  = 0
       i_acswupbc(iCell) = 0
       i_acswupt(iCell)  = 0
       i_acswuptc(iCell) = 0

       i_aclwdnb(iCell)  = 0
       i_aclwdnbc(iCell) = 0
       i_aclwdnt(iCell)  = 0
       i_aclwdntc(iCell) = 0
       i_aclwupb(iCell)  = 0
       i_aclwupbc(iCell) = 0
       i_aclwupt(iCell)  = 0
       i_aclwuptc(iCell) = 0

       acswdnb(iCell)  = 0._RKIND
       acswdnbc(iCell) = 0._RKIND
       acswdnt(iCell)  = 0._RKIND
       acswdntc(iCell) = 0._RKIND
       acswupb(iCell)  = 0._RKIND
       acswupbc(iCell) = 0._RKIND
       acswupt(iCell)  = 0._RKIND
       acswuptc(iCell) = 0._RKIND

       aclwdnb(iCell)  = 0._RKIND
       aclwdnbc(iCell) = 0._RKIND
       aclwdnt(iCell)  = 0._RKIND
       aclwdntc(iCell) = 0._RKIND
       aclwupb(iCell)  = 0._RKIND
       aclwupbc(iCell) = 0._RKIND
       aclwupt(iCell)  = 0._RKIND
       aclwuptc(iCell) = 0._RKIND
    enddo
 endif

!initialization of xicem:
 if(.not.config_do_restart .and. .not. config_do_DAcycling) then
!   call mpas_log_write('--- initialization of xicem:')
    do iCell = 1, nCellsSolve
       xicem(iCell) = xice(iCell)
    enddo
 endif

!initialization of the local sea-surface temperature when a diurnal cycle of the
!sea-surface temperature is applied. This avoids having the array sstsk equal to
!zero over land:
 if(.not. config_do_restart .and. .not. config_do_DAcycling) then
!   call mpas_log_write('--- initialization of sstsk:')
    do iCell = 1, nCellsSolve
       sstsk(iCell) = sst(iCell)
    enddo
 endif

!initialized the 1D ocean mixed-layer model (code from wrf module_sf_oml):
 if(config_oml1d) then
    if(.not. config_do_restart) then
       call mpas_log_write('--- initialization of 1D ocean mixed layer model ')
       do iCell = 1, nCellsSolve
          t_oml(iCell) = sst(iCell)
          t_oml_initial(iCell) = sst(iCell)
       enddo
       if(config_oml_hml0 .gt. 0) then
          do iCell = 1, nCellsSolve
             h_oml(iCell) = config_oml_hml0
             h_oml_initial(iCell) = config_oml_hml0
             hu_oml(iCell) = 0.
             hv_oml(iCell) = 0.
             t_oml_200m_initial(iCell) = sst(iCell) - 5.
          enddo
       elseif(config_oml_hml0 .eq. 0) then
! initializing with climatological mixed layer depth only:
          do iCell = 1, nCellsSolve
             h_oml(iCell) = h_oml_initial(iCell)
             hu_oml(iCell) = 0.
             hv_oml(iCell) = 0.
             t_oml_200m_initial(iCell) = sst(iCell) - 5.
          enddo
       else
          do iCell = 1, nCellsSolve
             h_oml(iCell) = h_oml_initial(iCell)
             ! WRF COMMENT:
             ! fill in near coast area with SST: 200 K was set as missing value in ocean pre-processing code
             if( (t_oml_200m_initial(iCell) > 200.) .and. (t_oml_200m_initial(iCell) <= 200.) )  &
                  t_oml_200m_initial(iCell) = sst(iCell)
          enddo
       endif
    endif
 endif

!initialization of temperatures needed for updating the deep soil temperature:
 if(.not. config_do_restart) then
    do iCell = 1, nCellsSolve
       nsteps_accum(iCell) = 0._RKIND
       ndays_accum(iCell)  = 0._RKIND
       tday_accum(iCell)   = 0._RKIND
       tyear_accum(iCell)  = 0._RKIND
       tyear_mean(iCell)   = tmn(iCell)
       do iLag = 1, nLags
          tlag(iLag,iCell) = tmn(iCell)
       enddo
    enddo
 endif

!read the input files that contain the monthly-mean ozone climatology on fixed pressure levels:
 if(config_o3climatology) call init_o3climatology(mesh,atm_input)

!initialization of global surface properties. set here for now, but may be moved when time
!manager is implemented:
 call landuse_init_forMPAS(dminfo,julday,mesh,configs,diag_physics,sfc_input)

!initialization of parameterized convective processes:
 if(config_convection_scheme .ne. 'off') &
    call init_convection(dminfo,mesh,configs,diag_physics)

!initialization of cloud microphysics processes:
 if(config_microp_scheme .ne. 'off') &
    call init_microphysics(dminfo,configs,mesh,state,time_lev,sfc_input,diag_physics)

!initialization of PBL processes:
 if(config_pbl_scheme .ne. 'off') call init_pbl(configs)

!initialization of surface layer processes:
 if(config_sfclayer_scheme .ne. 'off') call init_sfclayer(configs,diag_physics,sfc_input)

!initialization of land-surface model:
 if(config_lsm_scheme .ne. 'off') then
    if(config_lsm_scheme .eq. 'sf_noah' .or. config_lsm_scheme .eq. 'sf_ruc') then
       call init_lsm(dminfo,mesh,configs,diag_physics,sfc_input)
    elseif(config_lsm_scheme .eq. 'sf_noahmp') then
       call init_lsm_noahmp(configs,mesh,clock,diag_physics,diag_physics_noahmp,output_noahmp,sfc_input)
    endif
 endif

!initialization of shortwave radiation processes:
 init_done = .false.
 if(config_radt_sw_scheme.ne.'off') then
    if(trim(config_radt_sw_scheme) .eq. 'cam_sw') then
       call init_radiation_sw(dminfo,configs,mesh,atm_input,diag,diag_physics,state,time_lev)
       init_done = .true.
    else
       call init_radiation_sw(dminfo,configs)
    endif
 endif

!initialization of longwave radiation processes: if we run the CAM radiation codes, the initia
!lization of the longwave and shortwave parameterizations is the same, and needs to be called
!only once:
 if(config_radt_lw_scheme.ne.'off') then    
    if(trim(config_radt_lw_scheme) .eq. 'cam_lw') then
       if(.not. init_done) then
          call init_radiation_lw(dminfo,configs,mesh,atm_input,diag,diag_physics,state,time_lev)
       else
!          call mpas_log_write('')
!          call mpas_log_write('--- camrad lw initialization done above')
       endif
    else
       call init_radiation_lw(dminfo,configs)
    endif
 endif

!initialization of non-stationary gravity wave drag
 if(trim(gwdo_scheme).eq.'bl_ugwp_gwdo') then
    ! Read in ugwp_oro_data
    call MPAS_stream_mgr_read(stream_manager, streamID='ugwp_oro_data_in', whence=MPAS_STREAM_NEAREST, ierr=ierr)
    call MPAS_stream_mgr_reset_alarms(stream_manager, streamID='ugwp_oro_data_in', direction=MPAS_STREAM_INPUT, ierr=ierr)
    if (ngw_scheme) then
       call mpas_log_write('Initializing non-stationary GWD scheme',masterOnly=.true.)
       call MPAS_stream_mgr_read(stream_manager, streamID='ugwp_ngw_in', whence=MPAS_STREAM_NEAREST, ierr=ierr)
       call MPAS_stream_mgr_reset_alarms(stream_manager, streamID='ugwp_ngw_in', direction=MPAS_STREAM_INPUT, ierr=ierr)
       call mpas_pool_get_config(configs,'config_knob_ugwp_tauamp', knob_ugwp_tauamp)
       call mpas_pool_get_dimension(mesh,'lat',ntau_d1y)
       call mpas_pool_get_array(ngw_input,'LATS',ugwp_taulat)
       call mpas_pool_get_array(ngw_input,'jindx1_tau'  ,jindx1_tau  )
       call mpas_pool_get_array(ngw_input,'jindx2_tau'  ,jindx2_tau  )
       call mpas_pool_get_array(ngw_input,'ddy_j1tau'   ,ddy_j1tau   )
       call mpas_pool_get_array(ngw_input,'ddy_j2tau'   ,ddy_j2tau   )
       call mpas_pool_get_array(mesh,'latCell',latCell)
       call mpas_log_write('--- Initializing UGWP non-stationary GWD parameters ---',masterOnly=.true.)
       call ugwpv1_ngw_init(latCell,nVertLevels,config_dt,rdzw,dzu,ntau_d1y,          &
               knob_ugwp_tauamp,ugwp_taulat,jindx1_tau,jindx2_tau,ddy_j1tau,ddy_j2tau)
       call mpas_log_write('--- UGWP non-stationary GWD parameters initialized ---',masterOnly=.true.)
    endif
 endif

! call mpas_log_write('')
! call mpas_log_write('--- end subroutine physics_init')
! call mpas_log_write('')

 end subroutine physics_init

!=================================================================================================================
 subroutine init_physics_flags(state,f_qc,f_qr,f_qi,f_qs,f_qg,f_qoz,f_nc,f_ni,f_nifa,f_nwfa,f_nbca)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: state

!output arguments:
 logical,intent(out):: f_qc,f_qr,f_qi,f_qs,f_qg,f_qoz
 logical,intent(out):: f_nc,f_ni,f_nifa,f_nwfa,f_nbca

!local pointers:
 integer,pointer:: index_qc,index_qr,index_qi,index_qs,index_qg
 integer,pointer:: index_nc,index_ni,index_nifa,index_nwfa

!-----------------------------------------------------------------------------------------------------------------

!initializes the logicals assigned to mixing ratios:
 f_qc  = .false.
 f_qr  = .false.
 f_qi  = .false.
 f_qs  = .false.
 f_qg  = .false.
 f_qoz = .false. !qoz is not defined in Registry.xml and f_qoz is initialized to false.
 call mpas_pool_get_dimension(state,'index_qc',index_qc)
 call mpas_pool_get_dimension(state,'index_qr',index_qr)
 call mpas_pool_get_dimension(state,'index_qi',index_qi)
 call mpas_pool_get_dimension(state,'index_qs',index_qs)
 call mpas_pool_get_dimension(state,'index_qg',index_qg)

 if(index_qc .gt. -1) f_qc = .true.
 if(index_qr .gt. -1) f_qr = .true.
 if(index_qi .gt. -1) f_qi = .true.
 if(index_qs .gt. -1) f_qs = .true.
 if(index_qg .gt. -1) f_qg = .true.

!initializes the logical assigned to number concentrations:
 f_nc   = .false.
 f_ni   = .false.
 f_nifa = .false.
 f_nwfa = .false.
 f_nbca = .false. !nbca is not defined in Registry.xml - therefore f_nc is initialized to false.
 call mpas_pool_get_dimension(state,'index_nc'  ,index_nc  )
 call mpas_pool_get_dimension(state,'index_ni'  ,index_ni  )
 call mpas_pool_get_dimension(state,'index_nifa',index_nifa)
 call mpas_pool_get_dimension(state,'index_nwfa',index_nwfa)

 if(index_nc   .gt. -1) f_nc   = .true.
 if(index_ni   .gt. -1) f_ni   = .true.
 if(index_nifa .gt. -1) f_nifa = .true.
 if(index_nwfa .gt. -1) f_nwfa = .true.

 end subroutine init_physics_flags

!=================================================================================================================
 subroutine init_dirs_forphys(mesh)
!=================================================================================================================

!inout arguments:
!----------------
 type(mpas_pool_type),intent(in):: mesh

!local pointers:
 integer,pointer:: nCells
 real(kind=RKIND),dimension(:),pointer:: latCell,lonCell
 real(kind=RKIND),dimension(:,:),pointer:: east,north

!local variables:
 integer:: iCell

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_dimension(mesh,'nCells',nCells)

 call mpas_pool_get_array(mesh,'latCell',latCell)
 call mpas_pool_get_array(mesh,'lonCell',lonCell)
 call mpas_pool_get_array(mesh,'east'   ,east   )
 call mpas_pool_get_array(mesh,'north'  ,north  )

!Compute unit vectors in east and north directions for each cell:
 do iCell = 1, nCells

    east(1,iCell) = -sin(lonCell(iCell))
    east(2,iCell) =  cos(lonCell(iCell))
    east(3,iCell) =  0.0
    call r3_normalize(east(1,iCell), east(2,iCell), east(3,iCell))

    north(1,iCell) = -cos(lonCell(iCell))*sin(latCell(iCell))
    north(2,iCell) = -sin(lonCell(iCell))*sin(latCell(iCell))
    north(3,iCell) =  cos(latCell(iCell))
    call r3_normalize(north(1,iCell), north(2,iCell), north(3,iCell))

 end do

 end subroutine init_dirs_forphys

!=================================================================================================================
 subroutine r3_normalize(ax, ay, az)
!=================================================================================================================
!normalizes the vector (ax, ay, az)

 real (kind=RKIND), intent(inout) :: ax, ay, az
 real (kind=RKIND) :: mi

!-----------------------------------------------------------------------------------------------------------------

 mi = 1.0 / sqrt(ax**2 + ay**2 + az**2)
 ax = ax * mi
 ay = ay * mi
 az = az * mi

 end subroutine r3_normalize

!=================================================================================================================
 end module mpas_atmphys_init
!=================================================================================================================
