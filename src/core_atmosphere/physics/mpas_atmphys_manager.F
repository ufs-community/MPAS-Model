! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_manager
 use mpas_kind_types
 use mpas_pool_routines
 use mpas_timekeeping
 use mpas_stream_manager

 use mpas_atmphys_constants
 use mpas_atmphys_vars
 use mpas_atmphys_o3climatology
 use mpas_atmphys_update_surface
 use mpas_atmphys_utilities

 implicit none
 private
 public:: physics_timetracker,physics_run_init

 integer, public:: year                 !Current year.
 integer, public:: month                !Current month.
 integer, public:: day                  !Current day of the month.
 integer, public:: julday               !Initial Julian day.
 real(kind=RKIND), public:: curr_julday !Current Julian day (= 0.0 at 0Z on January 1st).
 real(kind=RKIND), public:: gmt         !Greenwich mean time hour of model start (hr)

!defines alarms for calling the long- and short-wave radiation codes, for calling the convec-
!tion scheme:
 character(len=*), parameter:: radtlwAlarmID = 'radtlw'
 character(len=*), parameter:: radtswAlarmID = 'radtsw'
 character(len=*), parameter:: convAlarmID   = 'conv'
 character(len=*), parameter:: pblAlarmID    = 'pbl'

!defines alarm to update the surface boundary conditions:
 character(len=*), parameter:: sfcbdyAlarmID = 'sfcbdy'

!defines alarm to update the background surface albedo and the greeness fraction:
 character(len=*), parameter:: greenAlarmID  = 'green'

!defines alarm to update the ozone path length,the trace gas path length,the total emissivity,
!and the total absorptivity in the "CAM" long-wave radiation codes. The default time interval
!between updates is 6 hours and is set with config_camrad_abs_update (00:30:00).
 character(len=*), parameter:: camAlarmID    = 'cam'

!defines alarm to save the CAM arrays absnst, absnxt, and emstot to restart files. When the
!alarm rings, the local arrays absnt_p, absnxt_p, and emstot_p are copied to the MPAS arrays
!for writing to restart files at the bottom of the time-step:
 character(len=*), parameter:: camlwAlarmID  = 'camlw'
 type(MPAS_TimeInterval_Type):: camlwTimeStep

!defines alarm to check if the accumulated rain due to cloud microphysics and convection is
!greater than its maximum allowed value:
 character(len=*), parameter:: acrainAlarmID   = 'acrain'
 type(MPAS_TimeInterval_Type):: acrainTimeStep

!defines alarm to check if the accumulated radiation diagnostics due to long- and short-wave
!radiation is greater than its maximum allowed value:
 character(len=*), parameter:: acradtAlarmID   = 'acradt'
 type(MPAS_TimeInterval_Type):: acradtTimeStep

!defines alarm to compute some physics diagnostics, such as radar reflectivity:
 character(len=*), parameter:: diagAlarmID     = 'diag'

 integer :: h, m, s, s_n, s_d, DoY, yr
 real(kind=RKIND) :: utc_h


!Main time manager for physics parameterizations.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines in mpas_atmphys_manager:
! ------------------------------------
! physics_timetracker: check alarms and update boundary conditions if needed.
! physics_run_init   : create and initialize alarms used for physics parameterizations.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * added initialization of variable sf_surface_physics in subroutine physics_run_init. see
!   definition of sf_surface_physics in mpas_atmphys_vars.F
!   Laura D. Fowler (laura@ucar.edu) / 2013-03-11.
! * removed the namelist option config_conv_shallow_scheme and associated sourcecode.
!   Laura D. Fowler (laura@ucar.edu) / 2013-05-29.
! * added call to subroutine o3climatology_from_MPAS to interpolate the climatological ozone
!   mixing ratios to the current julian day.
!   Laura D. Fowler (laura@ucar.edu) / 2013-07-03.
! * added domain%dminfo in call to subroutine physics_update_sst to print local and global
!   min and max values of the updated sea-surface temperatures and sea-ice fractions.
!   Laura D. Fowler (laura@ucar.edu) / 2013-07-24.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * renamed config_conv_deep_scheme to config_convection_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2014-09-18.
! * added initialization of variables has_reqc,has_reqi,and has_reqs needed in the calls to radiation codes
!   rrtmg_lwrad and rrmtg_swrad.
!   Laura D. Fowler (laura@ucar.edu) / 2016-07-07.
! * in subroutine physics_run_init, removed the initialization of the local variable gwdo_scheme. gwdo_scheme
!   is no longer needed and can be replaced with config_gwdo_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2017-02-16.
! * in subroutine physics_run_init, removed the initialization of the local variable lsm_scheme. lsm_scheme
!   is no longer needed and can be replaced with config_lsm_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2017-02-16.
! * in subroutine physics_run_init, removed the initialization of the local variable sfclayer_scheme.
!   sfclayer_scheme is no longer needed and can be replaced with config_sfclayer_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2017-02-16.
! * in subroutine physics_run_init, removed the initialization of the local variable pbl_scheme. pbl_scheme
!   is no longer needed and can be replaced with config_pbl_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2017-02-16.
! * in subroutine physics_run_init, removed the initialization of the local variable radt_cld_scheme.
!   radt_cld_scheme is no longer needed and can be replaced with config_radt_cld_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2017-02-16.
! * in subroutine physics_run_init, removed the initialization of the local variable radt_lw_scheme.
!   radt_lw_scheme is no longer needed and can be replaced with config_radt_lw_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2017-02-16.
! * in subroutine physics_run_init, removed the initialization of the local variable radt_sw_scheme.
!   radt_sw_scheme is no longer needed and can be replaced with config_radt_sw_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2017-02-16.
! * in subroutine physics_run_init, removed the initialization of the local variable convection_scheme.
!   convection_scheme is no longer needed and can be replaced with config_convection_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2017-02-16.
! * in subroutine physics_run_init, removed the initialization of the local variable microp_scheme.
!   microp_scheme is no longer needed and can be replaced with config_microp_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2017-02-16.
! * replaced the option "noah" with "sf_noah" to run the NOAH land surface scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2022-02-18.


 contains


!=================================================================================================================
 subroutine physics_timetracker(domain,dt,clock,itimestep,xtime_s)
!=================================================================================================================

 integer,intent(in):: itimestep
 real(kind=RKIND),intent(in):: dt
 real(kind=RKIND),intent(in) :: xtime_s

!inout arguments:
 type(MPAS_Clock_type),intent(inout):: clock
 type(domain_type),intent(inout)    :: domain

!local pointers:
 logical,pointer:: config_frac_seaice,        &
                   config_o3climatology,      &
                   config_sfc_albedo,         &
                   config_sst_update,         &
                   config_gvf_update,         &
                   config_sstdiurn_update,    &
                   config_deepsoiltemp_update

 character(len=StrKIND),pointer:: config_convection_scheme, &
                                  config_radt_lw_scheme,    &
                                  config_radt_sw_scheme

 character(len=StrKIND),pointer:: config_conv_interval,    &
                                  config_radtlw_interval,  &
                                  config_radtsw_interval

 type(block_type),pointer    :: block
 type(mpas_pool_type),pointer:: mesh
 type(mpas_pool_type),pointer:: configs
 type(mpas_pool_type),pointer:: diag_physics
 type(mpas_pool_type),pointer:: sfc_input
 type(mpas_pool_type),pointer:: atm_input

!local variables:

 type(MPAS_Time_Type):: currTime
 type (MPAS_TimeInterval_type) :: dtInterval

 logical:: LeapYear
 character(len=StrKIND)   :: timeStamp
 integer:: ierr
 real(kind=RKIND):: utc_s
 real(kind=RKIND):: xtime_m

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('')
!call mpas_log_write('--- enter subroutine physics_timetracker: itimestep = $i', intArgs=(/itimestep/))

 call mpas_pool_get_config(domain%blocklist%configs,'config_convection_scheme',config_convection_scheme)
 call mpas_pool_get_config(domain%blocklist%configs,'config_radt_lw_scheme'   ,config_radt_lw_scheme   )
 call mpas_pool_get_config(domain%blocklist%configs,'config_radt_sw_scheme'   ,config_radt_sw_scheme   )

 call mpas_pool_get_config(domain%blocklist%configs,'config_conv_interval'  ,config_conv_interval  )
 call mpas_pool_get_config(domain%blocklist%configs,'config_radtlw_interval',config_radtlw_interval)
 call mpas_pool_get_config(domain%blocklist%configs,'config_radtsw_interval',config_radtsw_interval)

 call mpas_pool_get_config(domain%blocklist%configs,'config_frac_seaice'        ,config_frac_seaice        )
 call mpas_pool_get_config(domain%blocklist%configs,'config_o3climatology'      ,config_o3climatology      )
 call mpas_pool_get_config(domain%blocklist%configs,'config_sfc_albedo'         ,config_sfc_albedo         )
 call mpas_pool_get_config(domain%blocklist%configs,'config_sst_update'         ,config_sst_update         )
 call mpas_pool_get_config(domain%blocklist%configs,'config_gvf_update'         ,config_gvf_update         )
 call mpas_pool_get_config(domain%blocklist%configs,'config_sstdiurn_update'    ,config_sstdiurn_update    )
 call mpas_pool_get_config(domain%blocklist%configs,'config_deepsoiltemp_update',config_deepsoiltemp_update)

!update the current julian day and current year:

 currTime = mpas_get_clock_time(clock,MPAS_NOW,ierr)
 call mpas_get_time(curr_time=currTime,dateTimeString=timeStamp,YYYY=yr,MM=month,DD=day,H=h,M=m, &
                    S=s,S_n=s_n,S_d=s_d,DoY=DoY,ierr=ierr)

 utc_h = real(h) + real(m) / 60.0 + real(s + s_n / s_d) / 3600.0
 utc_s = real(s + s_n + s_d) 
 year   = yr
 julday = DoY
 curr_julday = real(julday-1) + utc_h / 24.0
 LeapYear = isLeapYear(year)
!call mpas_log_write('     YEAR        = $i', intArgs=(/year/))
!call mpas_log_write('     JULDAY      = $i', intArgs=(/julday/))
!call mpas_log_write('     GMT         = $r', realArgs=(/gmt/))
!call mpas_log_write('     UTC_H       = $r', realArgs=(/utc_h/))
!call mpas_log_write('     CURR_JULDAY = $r', realArgs=(/curr_julday/))
!call mpas_log_write('     LEAP_YEAR   = $l', logicArgs=(/LeapYear/))
!call mpas_log_write('     TIME STAMP  = '//trim(timeStamp))

 block => domain % blocklist
 do while(associated(block))

    call mpas_pool_get_subpool(block%structs,'mesh'        ,mesh        )
    call mpas_pool_get_subpool(block%structs,'sfc_input'   ,sfc_input   )
    call mpas_pool_get_subpool(block%structs,'atm_input'   ,atm_input   )
    call mpas_pool_get_subpool(block%structs,'diag_physics',diag_physics)


    !update the background surface albedo and greeness of vegetation: interpolation of input
    !monthly values to current day:
    if(mpas_is_alarm_ringing(clock,greenAlarmID,ierr=ierr)) then
       call mpas_reset_clock_alarm(clock,greenAlarmID,ierr=ierr)
       call mpas_log_write('--- time to update background surface albedo, greeness fraction.')
       call physics_update_surface(timeStamp,config_sfc_albedo,config_gvf_update,mesh,sfc_input,diag_physics)
    endif

    !update surface boundary conditions with input sea-surface temperatures and fractional
    !sea-ice coverage:
    if(mpas_is_alarm_ringing(clock,sfcbdyAlarmID,ierr=ierr)) then
       call mpas_reset_clock_alarm(clock,sfcbdyAlarmID,ierr=ierr)
       if(config_sst_update) &
          call physics_update_sst(domain%dminfo,config_frac_seaice,mesh,sfc_input,diag_physics)
    endif   

    !apply a diurnal cycle to the sea-surface temperature:
    if(config_sstdiurn_update) &
       call physics_update_sstskin(dt_dyn,mesh,diag_physics,sfc_input)

    !update the deep soil temperature:
    if(config_deepsoiltemp_update) &
       call physics_update_deepsoiltemp(LeapYear,dt_dyn,curr_julday,mesh, &
                    sfc_input,diag_physics)

    block => block % next
 end do 

!check to see if it is time to run the longwave and shortwave radiation codes:
 if(trim(config_radt_lw_scheme) /= "off") then
    l_radtlw = .false.

    if(config_radtlw_interval /= "none") then
       if(mpas_is_alarm_ringing(clock,radtlwAlarmID,ierr=ierr)) then
          call mpas_reset_clock_alarm(clock,radtlwAlarmID,ierr=ierr)
          l_radtlw = .true.
       endif
    elseif(config_radtlw_interval == "none") then
       l_radtlw = .true.
    endif
    call mpas_log_write('--- time to run the LW radiation scheme L_RADLW = $l',logicArgs=(/l_radtlw/))
 endif

 if(trim(config_radt_sw_scheme) /= "off") then
    l_radtsw = .false.

    if(config_radtsw_interval /= "none") then
       if(mpas_is_alarm_ringing(clock,radtswAlarmID,ierr=ierr)) then
          call mpas_reset_clock_alarm(clock,radtswAlarmID,ierr=ierr)
          l_radtsw = .true.
       endif
    elseif(config_radtsw_interval == "none") then
       l_radtsw = .true.
    endif
    call mpas_log_write('--- time to run the SW radiation scheme L_RADSW = $l',logicArgs=(/l_radtsw/))
 endif

!check to see if it is time to run the parameterization of convection:
 if(trim(config_convection_scheme) /= "off") then
    l_conv = .false.

    if(config_conv_interval /= "none") then
       if(mpas_is_alarm_ringing(clock,convAlarmID,ierr=ierr)) then
          call mpas_reset_clock_alarm(clock,convAlarmID,ierr=ierr)
          l_conv = .true.
       endif
    elseif(config_conv_interval == "none") then
       l_conv = .true.
    endif
    call mpas_log_write('--- time to run the convection scheme L_CONV    = $l',logicArgs=(/l_conv/))
 endif

!check to see if it is time to update ozone to the current julian day in the RRTMG radiation codes:
 if(config_o3climatology) then
    block => domain % blocklist
    do while(associated(block))

       if((l_radtlw .and. trim(config_radt_lw_scheme) .eq. "rrtmg_lw") .or. &
           l_radtsw .and. trim(config_radt_sw_scheme) .eq. "rrtmg_sw" ) then
          call o3climatology_from_MPAS(curr_julday,mesh,atm_input,diag_physics)
          call mpas_log_write('--- time to update the ozone climatology for RRTMG radiation codes')
       endif
    
       block => block % next
    end do 
 endif

!check to see if it is time to update the ozone trace gas path lengths,the total emissivity,
!and the total absorptivity in the "CAM" long-wave radiation codes.
 if(trim(config_radt_lw_scheme) .eq. "cam_lw" .or. &
    trim(config_radt_sw_scheme) .eq. "cam_sw" ) then

    doabsems = .false.
    if(mpas_is_alarm_ringing(clock,camAlarmID,ierr=ierr)) then
       call mpas_reset_clock_alarm(clock,camAlarmID,ierr=ierr)
       doabsems = .true.
    endif
    call mpas_log_write('--- time to update CAM absorptivity and emissivity arrays DOABSEMS =$l',logicArgs=(/doabsems/))
 endif

!check to see if it is time to save the local CAM arrays absnst_p, absnxt_p, and emstot_p to
!the MPAS arrays:
 if(trim(config_radt_lw_scheme) .eq. "cam_lw") then
    l_camlw = .false.
    if(mpas_is_alarm_ringing(clock,camlwAlarmID,camlwTimeStep,ierr=ierr)) then
       call mpas_reset_clock_alarm(clock,camlwAlarmID,camlwTimeStep,ierr=ierr)
       l_camlw = .true.
    endif
    call mpas_log_write('--- time to write local CAM arrays to MPAS arrays L_CAMLW          = $l',logicArgs=(/l_camlw/))
 endif

!check to see if it is time to apply limit to the accumulated rain due to cloud microphysics
!and convection:
 if(trim(config_convection_scheme) /= "off") then
    l_acrain = .false.
    if(mpas_is_alarm_ringing(clock,acrainAlarmID,acrainTimeStep,ierr=ierr)) then
       call mpas_reset_clock_alarm(clock,acrainAlarmID,acrainTimeStep,ierr=ierr)
       l_acrain = .true.
    endif
    call mpas_log_write('--- time to apply limit to accumulated rainc and rainnc L_ACRAIN   = $l',logicArgs=(/l_acrain/))
 endif

!check to see if it is time to apply limit to the accumulated radiation diagnostics due to
!long- and short-wave radiation:
 if(trim(config_radt_lw_scheme) /= "off" .or. trim(config_radt_sw_scheme) /= "off") then
    l_acradt = .false.
    if(mpas_is_alarm_ringing(clock,acradtAlarmID,acradtTimeStep,ierr=ierr)) then
       call mpas_reset_clock_alarm(clock,acradtAlarmID,acradtTimeStep,ierr=ierr)
       l_acradt = .true.
    endif
    call mpas_log_write('--- time to apply limit to accumulated radiation diags. L_ACRADT   = $l',logicArgs=(/l_acradt/))
 endif

!check to see if it is time to calculate additional physics diagnostics:
 l_diags = .false.
 if (mpas_is_alarm_ringing(clock,diagAlarmID,ierr=ierr)) then
    call mpas_reset_clock_alarm(clock,diagAlarmID,ierr=ierr)
 end if
 call mpas_set_timeInterval(dtInterval,dt=dt,ierr=ierr)
 if (mpas_is_alarm_ringing(clock,diagAlarmID,interval=dtInterval,ierr=ierr)) then
    l_diags = .true.
 end if
 call mpas_log_write('--- time to calculate additional physics_diagnostics               = $l',logicArgs=(/l_diags/))

 end subroutine physics_timetracker

!=================================================================================================================
 subroutine physics_run_init(configs,mesh,state,clock,stream_manager)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: state
 type(MPAS_Clock_type):: clock
 type (MPAS_streamManager_type), intent(inout) :: stream_manager

!local pointers:
 character(len=StrKIND),pointer:: config_convection_scheme, &
                                  config_lsm_scheme,        &
                                  config_microp_scheme,     &
                                  config_radt_lw_scheme,    &
                                  config_radt_sw_scheme

 character(len=StrKIND),pointer:: config_conv_interval,     &
                                  config_pbl_interval,      &
                                  config_radtlw_interval,   &
                                  config_radtsw_interval,   &
                                  config_bucket_update,     &
                                  config_camrad_abs_update, &
                                  config_greeness_update

 logical,pointer:: config_sst_update
 logical,pointer:: config_frac_seaice
 logical,pointer:: config_microp_re


 integer,pointer:: cam_dim1
 integer,pointer:: nMonths
 integer,pointer:: nAerosols,nAerLevels,nOznLevels
 integer,pointer:: nCellsSolve,nSoilLevels,nVertLevels,nlcat,nscat

 real(kind=RKIND):: dt

!local variables:
 type(MPAS_Time_Type):: startTime,alarmStartTime
 type(MPAS_TimeInterval_Type):: alarmTimeStep, alarmTimeStepHi
 character(len=StrKIND) :: stream_interval
 integer:: ierr

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('')
!call mpas_log_write('--- enter subroutine physics_run_init:')

 call mpas_pool_get_config(configs,'config_convection_scheme',config_convection_scheme)
 call mpas_pool_get_config(configs,'config_lsm_scheme'       ,config_lsm_scheme       )
 call mpas_pool_get_config(configs,'config_microp_scheme'    ,config_microp_scheme    )
 call mpas_pool_get_config(configs,'config_radt_lw_scheme'   ,config_radt_lw_scheme   )
 call mpas_pool_get_config(configs,'config_radt_sw_scheme'   ,config_radt_sw_scheme   )

 call mpas_pool_get_config(configs,'config_conv_interval'    ,config_conv_interval    )
 call mpas_pool_get_config(configs,'config_pbl_interval'     ,config_pbl_interval     )
 call mpas_pool_get_config(configs,'config_radtlw_interval'  ,config_radtlw_interval  )
 call mpas_pool_get_config(configs,'config_radtsw_interval'  ,config_radtsw_interval  )
 call mpas_pool_get_config(configs,'config_bucket_update'    ,config_bucket_update    )
 call mpas_pool_get_config(configs,'config_camrad_abs_update',config_camrad_abs_update)
 call mpas_pool_get_config(configs,'config_greeness_update'  ,config_greeness_update  )
 call mpas_pool_get_config(configs,'config_sst_update'       ,config_sst_update       )
 call mpas_pool_get_config(configs,'config_frac_seaice'      ,config_frac_seaice      )
 call mpas_pool_get_config(configs,'config_microp_re'        ,config_microp_re        )

 call mpas_get_timeInterval(mpas_get_clock_timestep(clock, ierr), dt=dt)


 call mpas_pool_get_dimension(mesh,'cam_dim1'   ,cam_dim1   )
 call mpas_pool_get_dimension(mesh,'nMonths'    ,nMonths    )
 call mpas_pool_get_dimension(mesh,'nAerLevels' ,nAerLevels )
 call mpas_pool_get_dimension(mesh,'nOznLevels' ,nOznLevels )
 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nSoilLevels',nSoilLevels)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nlcat'      ,nlcat      )
 call mpas_pool_get_dimension(mesh,'nscat'      , nscat     )
 call mpas_pool_get_dimension(state,'num_aerosols',nAerosols)

!initialization of gmt, julian day, and alarms:

 startTime = mpas_get_clock_time(clock, MPAS_START_TIME, ierr)
 call mpas_get_time(curr_time=startTime,YYYY=yr,H=h,M=m,S=s,S_n=s_n,S_d=s_d,DoY=DoY,ierr=ierr)
 utc_h  = real(h) + real(m) / 60.0 + real(s + s_n / s_d) / 3600.0
 year   = yr
 gmt    = utc_h
 julday = DoY
 curr_julday = real(julday-1)
! call mpas_log_write('     YEAR        =$i', intArgs=(/year/))
! call mpas_log_write('     JULDAY      =$i', intArgs=(/julday/))
! call mpas_log_write('     GMT         =$r', realArgs=(/gmt/))
! call mpas_log_write('     UTC_H       =$r', realArgs=(/utc_h/))
! call mpas_log_write('     CURR_JULDAY =$r', realArgs=(/curr_julday/))

!set alarms for calling the longwave and shortwave radiation schemes, the convection schemes,
!and the PBL schemes at intervals different (greater) than the dynamical time-step:

 if(trim(config_radtlw_interval) /= "none") then
    call mpas_set_timeInterval(alarmTimeStep,timeString=config_radtlw_interval,ierr=ierr)
    alarmStartTime = startTime
    call mpas_add_clock_alarm(clock,radtlwAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_run_init: error creating radtlwAlarmID')

    call mpas_get_timeInterval(interval=alarmTimeStep,dt=dt_radtlw,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_run_init: error defining dt_radtlw')

 elseif(trim(config_radtlw_interval) == "none") then
    dt_radtlw = dt

 else
    call physics_error_fatal('subroutine physics_run_init: dt_radtlw is not defined')

 endif

 if(trim(config_radtsw_interval) /= "none") then
    call mpas_set_timeInterval(alarmTimeStep,timeString=config_radtsw_interval,ierr=ierr)
    alarmStartTime = startTime
    call mpas_add_clock_alarm(clock,radtswAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_run_init: error creating alarm radtsw')

    call mpas_get_timeInterval(interval=alarmTimeStep,dt=dt_radtsw,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_run_init: error defining radtswAlarmID')

 elseif(trim(config_radtsw_interval) == "none") then
    dt_radtsw = dt

 else
    call physics_error_fatal('subroutine physics_run_init: dt_radtsw is not defined')

 endif

 if(trim(config_conv_interval) /= "none") then
    call mpas_set_timeInterval(alarmTimeStep,timeString=config_conv_interval,ierr=ierr)
    alarmStartTime = startTime
    call mpas_add_clock_alarm(clock,convAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_init: error creating convAlarmID')

    call mpas_get_timeInterval(interval=alarmTimeStep,dt=dt_cu,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_run_init: error defining dt_cu')

 elseif(trim(config_conv_interval) == "none") then
    dt_cu = dt

 else
    call physics_error_fatal('subroutine physics_run_init: dt_cu is not defined')

 endif

 if(trim(config_pbl_interval) /= "none") then
    call mpas_set_timeInterval(alarmTimeStep,timeString=config_pbl_interval,ierr=ierr)
    alarmStartTime = startTime
    call mpas_add_clock_alarm(clock,pblAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_init: error creating pblAlarmID')

    call mpas_get_timeInterval(interval=alarmTimeStep,dt=dt_pbl,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_run_init: error defining dt_pbl')

 elseif(trim(config_pbl_interval) == "none") then
    dt_pbl = dt

 else
    call physics_error_fatal('subroutine physics_run_init: dt_pbl is not defined')

 endif

!set alarm for updating the background surface albedo and the greeness fraction:
 call mpas_set_timeInterval(alarmTimeStep,timeString=config_greeness_update,ierr=ierr)
 alarmStartTime = startTime
 call mpas_add_clock_alarm(clock,greenAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_init: error creating alarm greeness')

!set alarm for updating the surface boundary conditions:
 if (config_sst_update) then
    call MPAS_stream_mgr_get_property(stream_manager, 'surface', MPAS_STREAM_PROPERTY_RECORD_INTV, stream_interval, &
                                      direction=MPAS_STREAM_INPUT, ierr=ierr)
    call mpas_set_timeInterval(alarmTimeStep,timeString=stream_interval,ierr=ierr)
    alarmStartTime = startTime
    call mpas_add_clock_alarm(clock,sfcbdyAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_init: error creating alarm sfcbdy')
 endif

!set alarm to update the ozone path length, the trace gas path length, the total emissivity,
!and the total absorptivity in the "CAM" long-wave radiation codes.
 if(trim(config_radt_lw_scheme) .eq. "cam_lw" .or. &
    trim(config_radt_sw_scheme) .eq. "cam_sw" ) then
    call mpas_set_timeInterval(alarmTimeStep,timeString=config_camrad_abs_update,ierr=ierr)
    alarmStartTime = startTime
    call mpas_add_clock_alarm(clock,camAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_init: error creating alarm CAM')
 endif

!set alarm to write the "CAM" local arrays absnst_p, absnxt_p, and emstot_p to the MPAS arrays
!for writing to the restart file at the bottom of the time-step:
 if(trim(config_radt_lw_scheme) .eq. "cam_lw" ) then
    call mpas_set_timeInterval(camlwTimeStep,dt=dt,ierr=ierr)
    call MPAS_stream_mgr_get_property(stream_manager, 'restart', MPAS_STREAM_PROPERTY_RECORD_INTV, stream_interval, &
                                      direction=MPAS_STREAM_OUTPUT, ierr=ierr)
    if(trim(stream_interval) /= 'none') then
       call mpas_set_timeInterval(alarmTimeStep,timeString=stream_interval,ierr=ierr)
       alarmStartTime = startTime + alarmTimeStep
       call mpas_add_clock_alarm(clock,camlwAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
       if(ierr /= 0) &
          call physics_error_fatal('subroutine physics_init: error creating alarm CAMLW')
    endif
 endif

!set alarm to check if the accumulated rain due to cloud microphysics and convection is
!greater than its maximum allowed value:
 if(config_bucket_update /= "none") then
    call mpas_set_timeInterval(acrainTimeStep,dt=dt,ierr=ierr)
    call mpas_set_timeInterval(alarmTimeStep,timeString=config_bucket_update,ierr=ierr)
    alarmStartTime = startTime + alarmTimeStep
    call mpas_add_clock_alarm(clock,acrainAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
       if(ierr /= 0) &
          call physics_error_fatal('subroutine physics_init: error creating alarm rain limit')
 endif

!set alarm to check if the accumulated radiation diagnostics due to long- and short-wave radiation
!is greater than its maximum allowed value:
 if(config_bucket_update /= "none") then
    call mpas_set_timeInterval(acradtTimeStep,dt=dt,ierr=ierr)
    call mpas_set_timeInterval(alarmTimeStep,timeString=config_bucket_update,ierr=ierr)
    alarmStartTime = startTime + alarmTimeStep
    call mpas_add_clock_alarm(clock,acradtAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
       if(ierr /= 0) &
          call physics_error_fatal('subroutine physics_init: error creating alarm radiation limit')
 endif

!set alarm to calculate physics diagnostics on IO outpt only:
 call MPAS_stream_mgr_get_property(stream_manager, 'output', MPAS_STREAM_PROPERTY_RECORD_INTV, stream_interval, &
                                   direction=MPAS_STREAM_OUTPUT, ierr=ierr)
 if(trim(stream_interval) /= 'none') then
    call mpas_set_timeInterval(alarmTimeStep,timeString=stream_interval,ierr=ierr)
    
    ! If the diagnostic interval is higher in resolution, use it instead
    call MPAS_stream_mgr_get_property(stream_manager, 'diagnostics', MPAS_STREAM_PROPERTY_RECORD_INTV, stream_interval, &
                                      direction=MPAS_STREAM_OUTPUT, ierr=ierr)
    if(trim(stream_interval) /= 'none') then
       call mpas_set_timeInterval(alarmTimeStepHi,timeString=stream_interval,ierr=ierr)
       if (alarmTimeStepHi < alarmTimeStep) then
          alarmTimeStep = alarmTimeStepHi
       end if
    end if
    alarmStartTime = startTime
    call mpas_add_clock_alarm(clock,diagAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_init: error creating alarm diag')
 else
    call MPAS_stream_mgr_get_property(stream_manager, 'diagnostics', MPAS_STREAM_PROPERTY_RECORD_INTV, stream_interval, &
                                      direction=MPAS_STREAM_OUTPUT, ierr=ierr)
    if(trim(stream_interval) /= 'none') then
       call mpas_set_timeInterval(alarmTimeStep,timeString=stream_interval,ierr=ierr)
       alarmStartTime = startTime
       call mpas_add_clock_alarm(clock,diagAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
       if(ierr /= 0) &
          call physics_error_fatal('subroutine physics_init: error creating alarm diag')
    end if
 endif

 call mpas_log_write(' ')
 call mpas_log_write('DT_RADTLW = $r',realArgs=(/dt_radtlw/))
 call mpas_log_write('DT_RADTSW = $r',realArgs=(/dt_radtsw/))
 call mpas_log_write('DT_CU     = $r',realArgs=(/dt_cu/))
 call mpas_log_write('DT_PBL    = $r',realArgs=(/dt_pbl/))

!initialization of physics dimensions to mimic a rectangular grid:
 ims=1   ; ime = nCellsSolve
 jms=1   ; jme=1
 kms=1   ; kme = nVertLevels+1

 ids=ims ; ide=ime + 1
 jds=jms ; jde=jme + 1
 kds=kms ; kde=kme 

 its=ims ; ite = ime 
 jts=jms ; jte = jme
 kts=kms ; kte = kme-1

 call mpas_log_write(' ')
 call mpas_log_write('IMS = $i   IME = $i',intArgs=(/ims,ime/))
 call mpas_log_write('JMS = $i   JME = $i',intArgs=(/jms,jme/))
 call mpas_log_write('KMS = $i   KME = $i',intArgs=(/kms,kme/))
 call mpas_log_write('IDS = $i   IDE = $i',intArgs=(/ids,ide/))
 call mpas_log_write('JDS = $i   JDE = $i',intArgs=(/jds,jde/))
 call mpas_log_write('KDS = $i   KDE = $i',intArgs=(/kds,kde/))
 call mpas_log_write('ITS = $i   ITE = $i',intArgs=(/its,ite/))
 call mpas_log_write('JTS = $i   JTE = $i',intArgs=(/jts,jte/))
 call mpas_log_write('KTS = $i   KTE = $i',intArgs=(/kts,kte/))

!initialization local physics variables:
 num_months = nMonths
 num_soils  = nSoilLevels
 num_landcat = nlcat
 num_soilcat = nscat

 if(trim(config_lsm_scheme) .eq. "sf_noah") sf_surface_physics = 2

!initialization of local physics time-steps:
!... dynamics:
 dt_dyn    = dt
!... cloud microphysics:
 dt_microp = dt_dyn
 n_microp  = 1
 if(trim(config_microp_scheme)=='mp_thompson' .or. &
    trim(config_microp_scheme)=='mp_tempo'    .or. &
    trim(config_microp_scheme)=='mp_thompson_aerosols') then
    dt_microp = 90._RKIND
    n_microp  = max(nint(dt_dyn/dt_microp),1)
    dt_microp = dt_dyn / n_microp
    if(dt_dyn <= dt_microp) dt_microp = dt_dyn
 endif
if( trim(config_microp_scheme)=='mp_nssl2m' ) then
    IF ( dt_dyn > 1.25001*60._RKIND) THEN ! max 75s dt
      n_microp  = max(nint(dt_dyn/60._RKIND),2)
      dt_microp = dt_dyn / n_microp
    ENDIF
 endif
 call mpas_log_write(' ')
 call mpas_log_write('--- specifics on cloud microphysics option microp_scheme = '//trim(config_microp_scheme))
 call mpas_log_write('--- dt_microp = $r', realArgs=(/dt_microp/))
 call mpas_log_write('--- n_microp  = $i', intArgs=(/n_microp/))
!... convection:
 l_conv = .false.
 n_cu   = nint(dt_cu/dt_dyn)
 n_cu   = max(n_cu,1)
!... radiation:
 l_radtlw = .false.
 l_radtsw = .false.
!... others:
 l_diags  = .false.
 l_camlw  = .false.
 l_acrain = .false.
 l_acradt = .false.

!initialization for CAM radiation schemes only:
 if(trim(config_radt_lw_scheme) .eq. "cam_lw" .or. &
    trim(config_radt_sw_scheme) .eq. "cam_sw" ) then

    !this logical initiates the calculation of "pre-dfined arrays in the longwave
    !radiation code. these arrays are calculated once if it is the beginning of a
    !new run or if it is a restart run:
    doabsems     = .false.
    cam_abs_dim1 = cam_dim1
    cam_abs_dim2 = nVertLevels + 1

    !initializes the number of aerosols, and the prescribed vertical dimensions for
    !aerosols and ozone mixing ratios:
    call mpas_log_write('')
    call mpas_log_write('--- doabsems     = $l', logicArgs=(/doabsems/))
    call mpas_log_write('--- cam_abs_dim1 = $i', intArgs=(/cam_abs_dim1/))
    call mpas_log_write('--- cam_abs_dim2 = $i', intArgs=(/cam_abs_dim2/))
    call mpas_log_write('--- nAerosols    = $i', intArgs=(/nAerosols/))
    call mpas_log_write('--- naerLevels   = $i', intArgs=(/naerLevels/))
    num_aerosols  = nAerosols
    num_aerLevels = naerLevels

 endif 

!initialization of number of ozone levels:
 num_oznlevels = noznLevels

!initialization of sea-ice threshold:
 if(.not. config_frac_seaice) then
    xice_threshold = 0.5
 elseif(config_frac_seaice) then
    xice_threshold = 0.02
 endif

!initialization for the calculation of the cloud effective radii of cloud water, cloud ice, and snow:
 has_reqc = 0
 has_reqi = 0
 has_reqs = 0
 if(config_microp_re) then
    if(trim(config_microp_scheme)=='mp_thompson'          .or. &
       trim(config_microp_scheme)=='mp_tempo'             .or. &
       trim(config_microp_scheme)=='mp_thompson_aerosols' .or. &
       trim(config_microp_scheme)=='mp_wsm6' .or. trim(config_microp_scheme)=='mp_nssl2m') then
       if(trim(config_radt_lw_scheme)=='rrtmg_lw' .and. trim(config_radt_sw_scheme)=='rrtmg_sw') then
          has_reqc = 1
          has_reqi = 1
          has_reqs = 1
       endif
    endif
 endif
 call mpas_log_write('--- has_reqc  = $i', intArgs=(/has_reqc/))
 call mpas_log_write('--- has_reqi  = $i', intArgs=(/has_reqi/))
 call mpas_log_write('--- has_reqs  = $i', intArgs=(/has_reqs/))
 call mpas_log_write(' ')

 end subroutine physics_run_init

!=================================================================================================================
 end module mpas_atmphys_manager
!=================================================================================================================
 
