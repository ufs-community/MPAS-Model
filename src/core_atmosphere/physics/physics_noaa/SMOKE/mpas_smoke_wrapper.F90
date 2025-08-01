!>\file mpas_smoke_wrapper.F90
!! This file is MPAS Smoke wrapper
!! Haiqin.Li@noaa.gov 09/2024

module mpas_smoke_wrapper

   use mpas_kind_types
   use mpas_pool_routines
   use mpas_constants,        only: cp
   use mpas_smoke_config
   use mpas_smoke_init 
   use module_plumerise,      only : ebu_driver
   use module_add_emiss_burn, only : add_emis_burn
   use dep_dry_simple_mod,    only : dry_dep_driver_simple
   use dep_dry_mod_emerson,   only : dry_dep_driver_emerson
   use module_wetdep_ls,      only : wetdep_ls

   implicit none

   private

   public :: mpas_smoke_driver

contains

    subroutine mpas_smoke_driver(                                                            &
           num_chem              , chemistry_start             , chem           ,            &    
           kemit                 , kbio, kfire, kvol, index_smoke_fine,                      &
           index_e_bb_in_smoke_fine, index_e_bb_out_smoke_fine,                              &
           frp_in                , frp_out,    fre_in, fre_out, hwp,                         &
           totprcp_prev24        , hwp_prev24     , frp_prev24,    fre_prev24,               &
           hfx_bb                , qfx_bb         ,  frac_grid_burned    ,                   &
           min_bb_plume          , max_bb_plume,                                             &
           coef_bb_dc            , e_bb_in, e_bb_out,    num_e_bb_in, num_e_bb_out,          &
           ddvel                 , wetdep_resolved       , tend_chem_settle      ,           & 
           do_mpas_smoke         ,                                                           &
           hwp_method            , hwp_alpha             , wetdep_ls_opt         ,           &
           wetdep_ls_alpha       , plumerise_opt         , plume_wind_eff       ,            &
           plume_alpha           , bb_emis_scale_factor, ebb_dcycle             ,            &
           drydep_opt            , pm_settling           , add_fire_heat_flux   ,            &
           add_fire_moist_flux   , plumerisefire_frq     ,                                   &
           ktau                  , dt                    , dxcell               ,            &
           area                  ,                                                           & 
           xland                 , u10                   , v10                  ,            &
           ust                   , xlat                  , xlong                ,            &
           tskin                 , pblh                  , t2m                  ,            &
           p8w                   , dz8w                  , z_at_w               ,            &
           p_phy                 , t_phy                 , u_phy                ,            &
           v_phy                 , qv                    , vvel                 ,            &
           pi_phy                , rho_phy               , kpbl                 ,            &
           nsoil                 , smois                 , tslb                 ,            &
           ivgtyp                , isltyp                , nlcat                ,            &
           swdown                , z0                    , snowh                ,            &
           julian                , rmol                  , raincv               ,            &
           rainncv               , dpt2m                 , znt                  ,            &
           mavail                , g                     , vegfra               ,            &
           landusef              , cldfrac               , ktop_deep            ,            &
           cp                    , rd                    , gmt                  ,            &
           ids       , ide       , jds       , jde       , kds       , kde      ,            &
           ims       , ime       , jms       , jme       , kms       , kme      ,            &
           its       , ite       , jts       , jte       , kts       , kte                   &
                                                                                             )
        
    implicit none

! intent arguments:
! array indexes
    integer,intent(in):: ids,ide,jds,jde,kds,kde,        &
                         ims,ime,jms,jme,kms,kme,        &
                         its,ite,jts,jte,kts,kte
! Timestep, day, constants
    real(RKIND),intent(in) :: dt, julian, g, cp, rd, gmt
! Time step #
    integer,intent(in):: ktau
! Dimensions and indexes
    integer,intent(in):: nsoil, nlcat, num_chem, chemistry_start
    integer,intent(in):: kemit, kbio, kfire, kvol
    integer,intent(in):: num_e_bb_in, num_e_bb_out
! 2D mesh arguments
    real(RKIND),intent(in), dimension(ims:ime, jms:jme)             :: xlat, xlong, dxcell, area, xland   ! grid
! 2D Met input
    integer,intent(in), dimension(ims:ime, jms:jme)                 :: isltyp, ivgtyp ! domainant soil, vegetation type
    integer,intent(in), dimension(ims:ime, jms:jme)                 :: kpbl          ! k-index of PBLH
    integer,intent(in), dimension(ims:ime, jms:jme),optional        :: ktop_deep
    real(RKIND),intent(in), dimension(ims:ime, jms:jme)             :: u10, v10      ! 10-m winds
    real(RKIND),intent(in), dimension(ims:ime, jms:jme)             :: tskin, t2m, dpt2m            ! temperature
    real(RKIND),intent(in), dimension(ims:ime, jms:jme)             :: pblh              ! PBL height [m]
    real(RKIND),intent(in), dimension(ims:ime, jms:jme)             :: vegfra
    real(RKIND),intent(in), dimension(ims:ime, jms:jme)             :: swdown, z0, snowh, znt
    real(RKIND),intent(in), dimension(ims:ime, jms:jme)             :: raincv, rainncv, mavail                    
    real(RKIND),intent(inout), dimension(ims:ime, jms:jme)          :: rmol, ust
! BB forecast input (previous 24 hours)
    real(RKIND),intent(in), dimension(ims:ims, jms:jme, 24),        & 
                                                   optional         :: totprcp_prev24, hwp_prev24, frp_prev24, fre_prev24
! 3D Met input 
    real(RKIND),intent(in), dimension(ims:ime, kms:kme, jms:jme)    :: p8w,    dz8w,    z_at_w, cldfrac,   &
                                                                       p_phy,  t_phy,   u_phy,  v_phy,     &
                                                                       pi_phy, rho_phy, vvel   
! 3D emission input
    real(RKIND),intent(in), dimension(ims:ime,1:kfire,jms:jme,1:num_e_bb_in),optional  :: e_bb_in
! JLS - TODO, if we update QV via moist flux, we will need to update the scalar in the driver
    real(RKIND),intent(inout), dimension(ims:ime, kms:kme, jms:jme)           :: qv
    real(RKIND),intent(in), dimension(ims:ime,1:nsoil, jms:jme)               :: smois, tslb
    real(RKIND),intent(in), dimension(ims:ime,1:nlcat, jms:jme)               :: landusef
! Chemistry indexes into MPAS scalar array
    integer, intent(in) :: index_smoke_fine
    integer, intent(in) :: index_e_bb_in_smoke_fine
    integer, intent(in) :: index_e_bb_out_smoke_fine
! 2D chemistry input (only) arrays 
    real(RKIND),intent(in),dimension(ims:ime, jms:jme),optional  :: frp_in, fre_in      ! Fire input
! 2D input/output arrays
    real(RKIND),intent(inout),dimension(ims:ime, jms:jme),optional :: frp_out, fre_out
    real(RKIND),intent(inout),dimension(ims:ime, jms:jme),optional :: hwp, coef_bb_dc
    real(RKIND),intent(inout),dimension(ims:ime, jms:jme),optional :: hfx_bb, qfx_bb, frac_grid_burned
    integer,intent(inout),dimension(ims:ime,jms:jme),optional      :: min_bb_plume, max_bb_plume
! 2D + chem output arrays
    real(RKIND),intent(inout), dimension(ims:ime, jms:jme, 1:num_chem),optional :: wetdep_resolved
    real(RKIND),intent(inout), dimension(ims:ime, jms:jme, 1:num_chem),optional :: ddvel
! 3D output arrays
    real(RKIND),intent(inout), dimension(ims:ime, kms:kme, jms:jme,1:num_e_bb_out),optional    :: e_bb_out
! 3D + chem output arrays
    real(RKIND),intent(inout), dimension(ims:ime, kms:kme, jms:jme, 1:num_chem),optional :: tend_chem_settle
    real(RKIND),intent(inout), dimension(ims:ime, kms:kme, jms:jme, 1:num_chem) :: chem
!----------------------------------
!>-- Local Variables
!>-- 3D met
    real(RKIND), dimension(ims:ime, kms:kme, jms:jme) :: rri,     &
                     wind_phy,theta_phy,zmid,t8w,relhum
!>-- indexes, time
    integer :: julday
!>- dust & chemistry variables
    real(RKIND), dimension(ims:ime, 1:nlcat, jms:jme) :: vegfrac
    ! JLS, temporary, need to read in like SMOKE_RRFS/MPAS
    real(RKIND), dimension(ims:ime, jms:jme) :: total_flashrate
!>-- Namelist options
     logical,intent(in)                :: do_mpas_smoke
     integer,intent(in)                :: hwp_method
     real(RKIND),intent(in)            :: hwp_alpha
     integer,intent(in)                :: wetdep_ls_opt
     real(kind=RKIND),intent(in)       :: wetdep_ls_alpha
     integer,intent(in)                :: plumerise_opt
     integer,intent(in)                :: plume_wind_eff
     real(kind=RKIND),intent(in)       :: plume_alpha
     real(kind=RKIND),intent(in)       :: bb_emis_scale_factor
     integer,intent(in)                :: ebb_dcycle
     integer,intent(in)                :: drydep_opt
     integer,intent(in)                :: pm_settling
     logical,intent(in)                :: add_fire_heat_flux
     logical,intent(in)                :: add_fire_moist_flux
     integer,intent(in)                :: plumerisefire_frq

!>- plume variables
    ! -- buffers
    real(RKIND), dimension(ims:ime, kms:kme, jms:jme) :: ebu!,ebu_coarse
    real(RKIND), dimension(ims:ime, jms:jme)          :: flam_frac,                               &
                                                         fire_hist, peak_hr,                      &
                                                         fire_end_hr, hwp_day_avg,                &
                                                         uspdavg2d, hpbl2d 
    real(RKIND), dimension(ims:ime, jms:jme)          :: lu_nofire, lu_qfire, lu_sfire
    integer,     dimension(ims:ime, jms:jme)          :: fire_type
    integer,     dimension(ims:ime, jms:jme)          :: kpbl_thetav
    logical                                           :: call_plume
    real(RKIND), parameter                            :: conv_frpi   = 1.e-06_RKIND  ! FRP conversion factor, MW to W
    real(RKIND), parameter                            :: conv_frei   = 1.e-06_RKIND  ! FRE conversion factor, MW-s to W-s
!>- Dry deposition - temporary - move to output 
    real(RKIND), dimension(ims:ime, jms:jme, 1:num_chem)          ::  drydep_flux_local
    real(RKIND), dimension(ims:ime, kms:kme, jms:jme, 1:num_chem) :: vgrav     ! gravitational settling velocit
    real(RKIND), dimension(ims:ime, kms:kme, jms:jme) :: thetav
!> -- other
    real(RKIND)    :: theta
    real(RKIND)    :: curr_secs
    integer        :: nbegin, nv
    integer        :: i, j, k, kp, n
    character(100) :: errmsg
    integer        :: errflg
    logical        :: do_plumerise
! FRP/plumerise related thresholds
    real(RKIND), parameter :: frp_min        = 1.e+7     ! Minimum FRP (Watts) to distribute smoke in PBL, 10MW
    real(RKIND), parameter :: frp_max        = 2.e+10    ! Maximum FRP over 3km Pixel, 20,000 MW
    real(RKIND), parameter :: fre_min        = -999.     ! Minimum FRE (Watt seconds) to distrubute smoke in PBL, TODO
    real(RKIND), parameter :: fre_max        = -999.     ! Maximumm FRE (Watt seconds) "                       ", TODO     
    real(RKIND), parameter :: zpbl_threshold = 2.e+3     ! Minimum PBL depth to have plume rise 
    real(RKIND), parameter :: uspd_threshold = 5.        ! Wind speed averaged across PBL depth to control smoke release levels 
    real(RKIND), parameter :: frp_wthreshold = 1.e+9     ! Minimum FRP (Watts) to have plume rise in windy conditions
    real(RKIND), parameter :: fre_wthreshold = 1.e+9     ! Minimum FRE (Watt-seconds) to have plume rise in windy conditions
    real(RKIND), parameter :: ebb_min        = 1.e-3     ! Minimum smoke emissions (ug/m2/s)

    errmsg = ''
    errflg = 0
 
  ! If not simulating smoke or pollen, get outta here...
    if ( .not. do_mpas_smoke) return

    uspdavg2d   = 0._RKIND
    hpbl2d      = 0._RKIND
    peak_hr     = 0._RKIND
    flam_frac   = 0._RKIND
    fire_type   = 0

    curr_secs = ktau * dt
    julday = int(julian)

    do_plumerise = .false.
    if (plumerise_opt .gt. 0 ) do_plumerise = .true. 
    ! plumerise frequency in minutes set up by the namelist input
    call_plume       = (do_plumerise .and. (plumerisefire_frq > 0))
    if (call_plume) call_plume = (mod(int(curr_secs), max(1, 60*plumerisefire_frq)) == 0) .or. (ktau == 2)
! 
!   Reorder chemistry indices -- TODO if ktau = 1?
!
    call set_scalar_indices(chemistry_start,                             &
                    index_smoke_fine)
!
!
!
    call mpas_log_write( ' Calling smoke prep')
    !>- get ready for chemistry run
    call mpas_smoke_prep(                                                   &
        do_mpas_smoke,                                                      &
        ktau, nlcat,cp,ebb_dcycle,ebb_min,                                  &
        xland,xlat,xlong,ivgtyp,isltyp,landusef,                            & ! JLS TODO LANDUSEF /= VEGFRA?
        snowh,u10,v10,t2m,dpt2m,mavail,hwp,hwp_day_avg,                     &
        index_e_bb_in_smoke_fine,num_e_bb_in,kfire,e_bb_in,                 &
        t_phy,u_phy,v_phy,p_phy,pi_phy,z_at_w,                              &     
        rho_phy,qv,relhum,rri,                                              &
        total_flashrate,                                                    &
        wind_phy,theta_phy,zmid,kpbl_thetav,                                &
        peak_hr,coef_bb_dc,fire_hist,                                       &
        lu_nofire, lu_qfire, lu_sfire, fire_type,                           &
        ids,ide, jds,jde, kds,kde,                                          &
        ims,ime, jms,jme, kms,kme,                                          &
        its,ite, jts,jte, kts,kte)
        
   if ( do_mpas_smoke ) then

! HERE IS WHERE WE'LL CALCULATE THE EMISSIONS OR READ IT IN
!     if ( calc_bb_emis_online ) then
!     if ( ebb_dcycle .eq. 1 ) then
!


    if (ktau==1) then
      do j=jts,jte
      do i=its,ite
        ebu(i,kts,j)= e_bb_in(i,kts,j,index_e_bb_in_smoke_fine)
        !ebu_coarse(i,kts,j) = e_bb_in(i,kts,j,index_e_bb_in_smoke_coarse)
        do k=kts+1,kte
         ebu(i,k,j)= 0._RKIND
         !ebu_coarse(i,k,j)= 0._RKIND
        enddo
      enddo
      enddo
    else
      do j=jts,jte
      do k=kts,kte
      do i=its,ite
      ! ebu is divided by coef_bb_dc since it is applied in the output
        ebu(i,k,j) = e_bb_out(i,k,j,index_e_bb_out_smoke_fine) / &
                     MAX(1.E-4_RKIND,coef_bb_dc(i,j))
        !ebu_coarse(i,k,j) = e_bb_out(i,k,j,index_e_bb_out_smoke_coarse) / &
        !             MAX(1.E-4_RKIND,coef_bb_dc(i,j))
      enddo
      enddo
      enddo
    endif
  ! Compute the heat/moisture fluxes
    if ( add_fire_heat_flux ) then
     do j = jts,jte
     do i = its,ite
       if ( coef_bb_dc(i,j)*frp_in(i,j) .ge. 1.E7_RKIND ) then
          hfx_bb(i,j)           = min(max(0._RKIND,0.88_RKIND * coef_bb_dc(i,j)*frp_in(i,j) / &
                                  0.55_RKIND / area(i,j)) ,5000._RKIND) ! W m-2 [0 - 10,000]
          frac_grid_burned(i,j) = min(max(0._RKIND, 1.3_RKIND * 0.0006_RKIND * &
                                  coef_bb_dc(i,j)*frp_in(i,j)/area(i,j) ), &
                                  1._RKIND)
       else
          hfx_bb(i,j)           = 0._RKIND
          frac_grid_burned(i,j) = 0._RKIND
       endif
     enddo
     enddo
    endif
   ! JLS, input emissions or scale?
    if (add_fire_moist_flux) then
      do j = jts,jte
      do i = its,ite
        if ( coef_bb_dc(i,j)*frp_in(i,j) .ge. 1.E7_RKIND ) then
           qfx_bb(i,j)           = 0._RKIND
        else
           qfx_bb(i,j)           = 0._RKIND
        endif
      enddo
      enddo
    endif

    ! compute wild-fire plumes
    if (call_plume) then
      ! Apply the diurnal cycle coefficient to frp_out ()
 ! JLS -- Should this be moved outside the "call_plume" IF block?
      do j=jts,jte
      do i=its,ite
        if ( fire_type(i,j) .eq. 4 ) then ! only apply scaling factor to wildfires
           frp_out(i,j) = min(bb_emis_scale_factor*frp_in(i,j)*coef_bb_dc(i,j),frp_max)
        else
           frp_out(i,j) = min(frp_in(i,j)*coef_bb_dc(i,j),frp_max)
        endif
      enddo
      enddo

      call mpas_log_write( ' Calling ebu_driver')
      call ebu_driver (                                               &
                 flam_frac,kfire,                                     &
                 e_bb_in(:,:,:,index_e_bb_in_smoke_fine),             &
                 ebu,                                                 &
                 !e_bb_in(:,:,:,index_e_bb_in_smoke_coarse),           &
                 !ebu_coarse,                                          &
                 theta_phy,qv,                                        &
                 rho_phy,vvel,u_phy,v_phy,pi_phy,wind_phy,            &
                 z_at_w,zmid,g,cp,rd,                                 &
                 frp_out, min_bb_plume, max_bb_plume,                 &
                 plume_wind_eff,                                      &
                 do_plumerise,                                        &
                 kpbl_thetav,kpbl,curr_secs,                          &
                 xlat, xlong, uspdavg2d, hpbl2d, plume_alpha,         &
                 frp_min, frp_wthreshold,                             &
                 zpbl_threshold, uspd_threshold,                      &
                 ids,ide, jds,jde, kds,kde,                           &
                 ims,ime, jms,jme, kms,kme,                           &
                 its,ite, jts,jte, kts,kte, errmsg, errflg            )
      if(errflg/=0) return
    end if

    ! -- add biomass burning emissions at every timestep
    if (addsmoke_flag == 1) then
     call mpas_log_write( ' Calling add_emis_burn')
     call add_emis_burn(dt,dz8w,rho_phy,pi,ebb_min,                   &
                        chem(:,:,:,p_smoke_fine),                     &
                        !chem(:,:,:,p_smoke_coarse),                   &
                        julday,gmt,xlat,xlong,                        &
                        fire_end_hr, peak_hr,curr_secs,               &
                        coef_bb_dc,fire_hist,hwp,hwp_day_avg,         &
                        !swdown,ebb_dcycle,ebu,ebu_coarse,fire_type,   &
                        swdown,ebb_dcycle,ebu,fire_type,   &
                        qv, add_fire_moist_flux,                      &
                        bb_emis_scale_factor,                         &   
                        ids,ide, jds,jde, kds,kde,                    &
                        ims,ime, jms,jme, kms,kme,                    &
                        its,ite, jts,jte, kts,kte                     )
    endif
  endif ! if do_mpas_smoke

    !>-- compute dry deposition, based on Emerson et al., (2020)
    if (drydep_opt == 1) then
     call mpas_log_write( ' Calling dry_dep_driver_emerson')
     call dry_dep_driver_emerson(rmol,ust,znt,num_chem,ddvel,         &
        vgrav,chem,dz8w,snowh,t_phy,p_phy,rho_phy,ivgtyp,g,dt,        &
        drydep_flux_local,tend_chem_settle,dbg_opt,                   &
        pm_settling,                                                  &
        ids,ide, jds,jde, kds,kde,                                    &
        ims,ime, jms,jme, kms,kme,                                    &
        its,ite, jts,jte, kts,kte, curr_secs, xlat, xlong             )
    !>-- compute dry deposition based on simple parameterization (HRRR-Smoke)
    elseif (drydep_opt == 2) then
     call dry_dep_driver_simple(rmol,ust,ndvel,ddvel,                 &
        ids,ide, jds,jde, kds,kde,                                    &
        ims,ime, jms,jme, kms,kme,                                    &
        its,ite, jts,jte, kts,kte                                     )
    else
        ddvel=0._RKIND
    endif

 !>- large-scale wet deposition
    if (wetdep_ls_opt == 1) then
       call mpas_log_write( ' Calling wetdep_ls')
       call  wetdep_ls(dt,chem,rainncv,qv,                            &
                     rho_phy,num_chem,dz8w,vvel,p_phy,                &
                     wetdep_ls_alpha,                                 &
                     wetdep_resolved,                                 &
                     ids,ide, jds,jde, kds,kde,                       &
                     ims,ime, jms,jme, kms,kme,                       &
                     its,ite, jts,jte, kts,kte                        )
    endif
    
    !>-- output of MPAS-Smoke
    do j=jts,jte
    do k=kts,kte
    do i=its,ite
       if (do_mpas_smoke) then
          e_bb_out(i,k,j,index_e_bb_out_smoke_fine)=ebu(i,k,j) * coef_bb_dc(i,j)
          !e_bb_out(i,k,j,index_e_bb_out_smoke_coarse)=ebu_coarse(i,k,j) * coef_bb_dc(i,j)
       endif
       ! TODO, floating point exception in debug
       !do nv = 1,num_chem
       !   chem(i,k,j,nv) = max(1.e-12_RKIND,min(5000._RKIND,chem(i,k,j,nv)))
       !enddo
    enddo
    enddo
    enddo
    
 end subroutine mpas_smoke_driver

 subroutine mpas_smoke_prep(                                                &
        do_mpas_smoke, &
        ktau, nlcat,cp,ebb_dcycle,ebb_min,                                  &
        xland,xlat,xlong,ivgtyp,isltyp,vegfrac,                             &
        snowh,u10,v10,t2m,dpt2m,wetness,hwp,hwp_day_avg,                    &
        index_e_bb_in_smoke_fine,num_e_bb_in,kfire,e_bb_in,                 &
        t_phy,u_phy,v_phy,p_phy,pi_phy,z_at_w,                              &
        rho_phy,qv,relhum,rri,                                              &
        total_flashrate,                                                    &
        wind_phy,theta_phy,zmid,kpbl_thetav,                                &
        peak_hr,coef_bb_dc,fire_hist,                                       &
        lu_nofire, lu_qfire, lu_sfire, fire_type,                           &
        ids,ide, jds,jde, kds,kde,                                          &
        ims,ime, jms,jme, kms,kme,                                          &
        its,ite, jts,jte, kts,kte)

    !intent arguments:
    integer,intent(in):: ids,ide,jds,jde,kds,kde,                           &
                         ims,ime,jms,jme,kms,kme,                           &
                         its,ite,jts,jte,kts,kte

    integer,intent(in)     :: ktau, nlcat,ebb_dcycle,                       &
                              kfire, num_e_bb_in, index_e_bb_in_smoke_fine
    logical,intent(in)     :: do_mpas_smoke
    real(RKIND),intent(in) :: cp, ebb_min
    integer,intent(in), dimension(ims:ime, jms:jme) :: isltyp, ivgtyp
    real(RKIND),intent(in),   dimension(ims:ime, nlcat, jms:jme) :: vegfrac

    real(RKIND),intent(in),   dimension(ims:ime, jms:jme) :: xland, xlat, xlong,                   &
                                        snowh, u10, v10, t2m, dpt2m, wetness
    real(RKIND),intent(in),   dimension(ims:ime, kms:kme, jms:jme) :: qv, z_at_w,                  &
                                        p_phy, t_phy, u_phy, v_phy, pi_phy, rho_phy
    real(RKIND),intent(out),  dimension(ims:ime, kms:kme, jms:jme) :: zmid, wind_phy,theta_phy,   &
                                       relhum, rri
    integer    ,intent(out),  dimension(ims:ime, jms:jme) :: fire_type, kpbl_thetav
    real(RKIND),intent(out),  dimension(ims:ime, jms:jme) :: lu_nofire, lu_qfire, lu_sfire
    real(RKIND),intent(out),  dimension(ims:ime, jms:jme) :: fire_hist, peak_hr, hwp_day_avg
    real(RKIND),intent(inout),dimension(ims:ime,jms:jme),optional :: hwp, coef_bb_dc
    real(RKIND),intent(out),  dimension(ims:ime, jms:jme) :: total_flashrate
    real(RKIND),intent(in),   dimension(ims:ime,1:kfire,jms:jme,1:num_e_bb_in),optional :: e_bb_in

    !local variables
    real(RKIND), parameter :: delta_theta4gust = 0.5
    real(RKIND) :: theta, wdgust, snoweq, term2, term3
!    real(RKIND), dimension(ims:ime, jms:jme) :: 
    real(RKIND), dimension(ims:ime, kms:kme, jms:jme) :: thetav

    integer :: i, j, k, k1, nv

    if ( do_mpas_smoke ) then
       if (ktau==1) then
         do j=jts,jte
         do i=its,ite
           fire_hist   (i,j) = 1._RKIND
           coef_bb_dc  (i,j) = 1._RKIND
           if (xlong(i,j)<230.) then
               peak_hr(i,j)= 0.0* 3600.     ! peak at 24 UTC, fires in Alaska
           elseif(xlong(i,j)<245.) then
               peak_hr(i,j)= 23.0* 3600.
           elseif (xlong(i,j)<260.) then
               peak_hr(i,j)= 22.0* 3600.    ! peak at 22 UTC, fires in the western US
           elseif (xlong(i,j)<275.) then
               peak_hr(i,j)= 21.0* 3600.
           elseif (xlong(i,j)<290.) then    ! peak at 20 UTC, fires in the eastern US
               peak_hr(i,j)= 20.0* 3600.
           else
               peak_hr(i,j)= 19.0* 3600.
           endif
         enddo
         enddo
       endif
    endif

    do j=jts,jte
    do k=kts,kte
    do i=its,ite
      zmid(i,k,j)=z_at_w(i,k,j)
    enddo
    enddo
    enddo

    do j=jts,jte
    do k=kts,kte
    do i=its,ite
      theta_phy(i,k,j) = t_phy(i,k,j)/pi_phy(i,k,j)*cp
      wind_phy(i,k,j) = sqrt(u_phy(i,k,j)**2. + v_phy(i,k,j)**2.)
      rri(i,k,j) = 1._RKIND/rho_phy(i,k,j)
    enddo
    enddo
    enddo

    ! JLS, TODO, lightning flashrate
    do j=jts,jte
    do i=its,ite
       total_flashrate(i,j) = 0._RKIND
    enddo
    enddo
    
    ! Calculate relative humidity [0.1 -- 0.95]
    do j=jts,jte
    do k=kts,kte
    do i=its,ite
         relhum(i,k,j) = max(.1_RKIND,MIN( .95_RKIND, qv(i,k,j) / &
               (3.80_RKIND*exp(17.27_RKIND*(t_phy(i,k,j)-273._RKIND)/ &
               (t_phy(i,k,j)-36._RKIND))/(.01_RKIND*p_phy(i,k,j)))))
    enddo
    enddo
    enddo

    !---- Calculate PBLH and K-PBL based on virtual potential temperature profile
    !---- First calculate THETAV
    do j = jts,jte
    do k = kts,kte
    do i = its,ite
       theta = t_phy(i,k,j) * (1.E5_RKIND/p_phy(i,k,j))**0.286_RKIND
       thetav(i,k,j) = theta * (1._RKIND + 0.61_RKIND * qv(i,k,j))
    enddo
    enddo
    enddo
    !---- Now use the UPP code to deterimine the height and level
    do i = its, ite
    do j = jts, jte
       if ( thetav(i,kts+1,j) .lt. ( thetav(i,kts,j) + delta_theta4gust) ) then
          do k = kts+1, kte
             k1 = k
!--- give theta-v at the sfc a 0.5K boost in the PBLH definition
             if ( thetav(i,kts+k-1,j) .gt. ( thetav(i,kts,j) + delta_theta4gust) ) then
                exit
             endif
          enddo
          kpbl_thetav(i,j) = k1
       else
          kpbl_thetav(i,j) = kts + 1
       endif
   enddo
   enddo

    if (do_mpas_smoke) then
      !RAR: change this to the fractional LU type; fire_type: 0- no fires, 1- Ag
      ! or urban fires, 2- prescribed fires in wooded area, 3- wildfires
       if (ebb_dcycle==2) then
         do j=jts,jte
         do i=its,ite
           if (e_bb_in(i,1,j,index_e_bb_in_smoke_fine)<ebb_min) then
              fire_type(i,j) = 0
              lu_nofire(i,j) = 1.0
           else
             ! Permanent wetlands, snow/ice, water, barren tundra:
             lu_nofire(i,j)= vegfrac(i,11,j) + vegfrac(i,15,j) + vegfrac(i,17,j) + vegfrac(i,20,j)
             ! cropland, urban, cropland/natural mosaic, barren and sparsely
             ! vegetated and non-vegetation areas:
             lu_qfire(i,j) = lu_nofire(i,j) + vegfrac(i,12,j) + vegfrac(i,13,j) + vegfrac(i,14,j) + vegfrac(i,16,j)
             ! Savannas and grassland fires, these fires last longer than the Ag fires:
             lu_sfire(i,j) = lu_qfire(i,j) + vegfrac(i,8,j) + vegfrac(i,9,j) + vegfrac(i,10,j)
             if (lu_nofire(i,j)>0.95) then ! no fires
               fire_type(i,j) = 0
             else if (lu_qfire(i,j)>0.9) then   ! Ag. and urban fires
               fire_type(i,j) = 1
             else if (xlong(i,j)>260. .AND. xlat(i,j)>25. .AND. xlat(i,j)<41.) then
               fire_type(i,j) = 2    ! slash burn and wildfires in the east, eastern temperate forest ecosystem
             else if (lu_sfire(i,j)>0.8) then
               fire_type(i,j) = 3    ! savanna and grassland fires
             else
               fire_type(i,j) = 4    ! potential wildfires
             end if
           end if
         end do
         end do
       endif ! ebb_dycycle == 2
   
       !>-- HWP: Pre-release of RRFSv1 method - using wind gust calculated via UPP Method
       do i=its, ite
       do j=jts, jte
         wdgust  =max(sqrt(u10(i,j)**2.+v10(i,j)**2.),3._RKIND)
         snoweq  =max((25._RKIND - snowh(i,j))/25._RKIND,0._RKIND)
!         term2   = max(t2m(i,j)-dpt2m(i,j),15._RKIND)**1.03 ! TODO, floating point exception
         term2   = 15._RKIND
         term3   = (1._RKIND-wetness(i,j))**0.4
         hwp(i,j)= 0.177_RKIND * wdgust**0.97 * term2 * ((1._RKIND-wetness(i,j))**0.4) * snoweq 
         hwp_day_avg(i,j)=hwp(i,j)
       enddo
       enddo
    endif ! do_mpas_smoke
   

  end subroutine mpas_smoke_prep
  

!> @}
  end module mpas_smoke_wrapper
