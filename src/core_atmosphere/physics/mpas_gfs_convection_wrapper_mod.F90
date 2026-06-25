!> MPAS -> GFS SAS convection bridge using the old-Tiedtke-style vertical
!> packing pattern.
!>
!> Scientific conventions used here:
!>   * The wrapper is stateless with respect to vertical level identity.
!>     Every physics call rebuilds a current-time MPAS<->SAS k mapping from
!>     the current pressure column, analogous to the old Tiedtke bridge.
!>   * MPAS and the GFS SAS/shalcnv interface are treated as BOTTOM-to-TOP:
!>     k=1 is the lowest model layer and k=km is the model top.
!>     The wrapper still rebuilds kmap every call, but kmap(i,k)=k.
!>     Returned tendencies are unpacked with the same current-call mapping.
!>   * psp, prslp, delp are in Pa, following sascnvn.meta/shalcnv.meta.
!>   * phil is geopotential, m2 s-2.
!>   * dot is pressure vertical velocity, Pa s-1, computed from MPAS w as
!>       omega ~= -rho*g*w .
!>   * shalcnv heat/evap are kinematic fluxes:
!>       heat = HFX/(rho_sfc*cp)  [K m s-1]
!>       evap = QFX/rho_sfc       [kg kg-1 m s-1]
!>   * SAS updates local T/q/u/v/qc/qi; MPAS receives tendencies only.
!>   * MPAS input t is physical temperature in K, not potential temperature.
!>     SAS receives t directly.
!>   * SAS dT/dt is converted to MPAS potential-temperature tendency for rthcuten.
!>   * Least-artificial dt=300 test: no artificial heat/moisture/momentum caps.
!>     Keep only NaN/bad-column protection and GFS deep/shallow parameter split.
!>   * No arbitrary tendency caps are applied. Columns with nonphysical
!>     temperature inputs are skipped rather than aborting the run.

! TEST CONFIG: latest SAMF deep+shallow, driver-staged T/qv,
! TKE-EDMF coupling through tkeh/maxMF/do_mynnedmf, kbot/ktop guard.
module mpas_gfs_convection_wrapper_mod

   use mpas_kind_types, only: RKIND
   use mpas_atmphys_constants, only: gravity, cp, R_d, R_v, xlv

   implicit none
   private
   public :: mpas_call_gfs_convection

   real(kind=RKIND), parameter :: t0 = 273.15_RKIND
   real(kind=RKIND), parameter :: qmin = 1.0e-12_RKIND

contains

   subroutine mpas_call_gfs_convection(im, km, dt,                              &
        pres_mid, pres_int, z_int, dx, u, v, t, qv, qc, qi, w,                      &
        hpbl, hfx, qfx, xland,                                                  &
        tkeh_in, maxmf_in, do_mynnedmf_in,                                      &
        rthcuten, rqvcuten, rqccuten, rqicuten, rucuten, rvcuten,               &
        raincv, pratec, cubot, cutop, ierr, errmsg)

      use samfdeepcnv, only: samfdeepcnv_run
      use samfshalcnv, only: samfshalcnv_run

      implicit none

      integer, intent(in) :: im, km
      real(kind=RKIND), intent(in) :: dt

      real(kind=RKIND), intent(in) :: pres_mid(im,km)      ! Pa, MPAS layer pressure
      real(kind=RKIND), intent(in) :: pres_int(im,km+1)    ! Pa, MPAS interface pressure
      real(kind=RKIND), intent(in) :: z_int(im,km+1)       ! m, interface height
      real(kind=RKIND), intent(in) :: dx(im)              ! cell length scale, m
      real(kind=RKIND), intent(in) :: u(im,km), v(im,km)   ! m s-1
      real(kind=RKIND), intent(in) :: t(im,km)             ! physical temperature, K; staged by MPAS driver
      real(kind=RKIND), intent(in) :: qv(im,km)            ! kg kg-1; staged by MPAS driver
      real(kind=RKIND), intent(in) :: qc(im,km), qi(im,km) ! kg kg-1
      real(kind=RKIND), intent(in) :: w(im,km)             ! m s-1
      real(kind=RKIND), intent(in) :: hpbl(im)             ! m
      real(kind=RKIND), intent(in) :: hfx(im), qfx(im)     ! W m-2, kg m-2 s-1
      real(kind=RKIND), intent(in) :: xland(im)            ! MPAS/WRF convention: 1 land, 2 water
      real(kind=RKIND), intent(in) :: tkeh_in(im,km+1)     ! interface TKE from TKE-EDMF, m2 s-2
      real(kind=RKIND), intent(in) :: maxmf_in(im)         ! max EDMF mass flux from PBL; 0 if unavailable
      logical, intent(in) :: do_mynnedmf_in                ! true when coupled to TKE-EDMF/MYNN-EDMF PBL

      real(kind=RKIND), intent(inout) :: rthcuten(im,km)   ! theta tendency, K s-1
      real(kind=RKIND), intent(inout) :: rqvcuten(im,km)   ! kg kg-1 s-1
      real(kind=RKIND), intent(inout) :: rqccuten(im,km)   ! kg kg-1 s-1
      real(kind=RKIND), intent(inout) :: rqicuten(im,km)   ! kg kg-1 s-1
      real(kind=RKIND), intent(inout) :: rucuten(im,km)    ! m s-2
      real(kind=RKIND), intent(inout) :: rvcuten(im,km)    ! m s-2
      real(kind=RKIND), intent(inout) :: raincv(im)        ! mm per physics step
      real(kind=RKIND), intent(inout) :: pratec(im)        ! mm s-1
      real(kind=RKIND), intent(inout) :: cubot(im)         ! diagnostic index
      real(kind=RKIND), intent(inout) :: cutop(im)         ! diagnostic index

      integer, intent(out) :: ierr
      character(len=*), intent(out) :: errmsg

      integer :: i, k, kk, kkdbg
      integer :: kmap(im,km)
      integer :: jcap, ncloud, mp_phys, mp_phys_mg
      integer :: kbot(im), ktop(im), kcnv(im), islimsk(im)
      integer :: errflg_deep, errflg_shal
      integer :: surf_k(im)
      logical :: bad_col(im)
      character(len=256) :: emsg_deep, emsg_shal

      real(kind=RKIND) :: hvap, cvap, cliq, fv, eps, epsm1, kappa
      real(kind=RKIND) :: clam_deep, c0_deep, c1_deep, betal_deep, betas_deep, evfact_deep, evfactl_deep, pgcon_deep
      real(kind=RKIND) :: clam_shal, c0_shal, c1_shal, pgcon_shal, evef_shal
      real(kind=RKIND) :: tv, rho, theta_fac
      real(kind=RKIND) :: raw_t
      real(kind=RKIND) :: dtemp_sas, rth_test, dtdt_test
      real(kind=RKIND) :: rqv_raw, rqc_raw, rqi_raw, ru_raw, rv_raw
      real(kind=RKIND) :: pif_surf, pif_top, pbot, ptop

      real(kind=RKIND) :: psp(im), delp(im,km), prslp(im,km), phil(im,km)
      real(kind=RKIND) :: dot(im,km)
      real(kind=RKIND) :: qlc(im,km), qli(im,km)
      real(kind=RKIND) :: q1(im,km), t1(im,km), u1(im,km), v1(im,km)
      real(kind=RKIND) :: q1_old(im,km), t1_old(im,km), u1_old(im,km), v1_old(im,km)
      real(kind=RKIND) :: qlc_old(im,km), qli_old(im,km)
      real(kind=RKIND) :: cldwrk(im), rn_deep(im), rn_shal(im)
      real(kind=RKIND) :: ud_mf(im,km), dd_mf(im,km), dt_mf(im,km)
      real(kind=RKIND) :: ud_mf_sh(im,km), dt_mf_sh(im,km)
      real(kind=RKIND) :: cnvw(im,km), cnvc(im,km)

      real(kind=RKIND) :: qlcn(im,km), qicn(im,km), w_upi(im,km), cf_upi(im,km)
      real(kind=RKIND) :: cnv_mfd(im,km), cnv_dqldt(im,km), clcn(im,km)
      real(kind=RKIND) :: cnv_fice(im,km), cnv_ndrop(im,km), cnv_nice(im,km)

      real(kind=RKIND) :: heat(im), evap(im)
      ! Latest CCPP SAMF deep/shallow interface work arrays.
      integer :: nn_samf, itc_samf, ntc_samf, ntk_samf, ntr_samf
      logical :: first_time_step_samf, restart_samf
      logical :: hwrf_samfdeep, hwrf_samfshal
      logical :: progsigma_samf, progomega_samf
      logical :: do_ca, ca_closure, ca_entr, ca_trigger
      logical :: do_mynnedmf, sigmab_coldstart
      real(kind=RKIND) :: asolfac, cscale, nthresh
      real(kind=RKIND) :: betadcu, betamcu, betascu
      real(kind=RKIND) :: cat_adj_deep, cat_adj_shal
      real(kind=RKIND) :: garea(im), fscav(1), maxMF(im), ca_deep(im), rainevap(im)
      real(kind=RKIND) :: tmf(im,km,1), qmicro(im,km), prevsq(im,km)
      real(kind=RKIND) :: sigmain(im,km), sigmaout(im,km)
      real(kind=RKIND) :: omegain(im,km), omegaout(im,km)
      real(kind=RKIND) :: qtr(im,km,2), qtr_sh(im,km,2)
      real(kind=RKIND) :: dqtr_deep(im,km,2), dqtr_shal(im,km,2)
      real(kind=RKIND) :: ten_t_deep(im,km), ten_u_deep(im,km), ten_v_deep(im,km)
      real(kind=RKIND) :: ten_t_shal(im,km), ten_u_shal(im,km), ten_v_shal(im,km)
      real(kind=RKIND) :: ten_q_deep(im,km,1), ten_q_shal(im,km,1)
      real(kind=RKIND) :: t1_sh(im,km), q1_sh(im,km), u1_sh(im,km), v1_sh(im,km)
      real(kind=RKIND) :: tkeh(im,km+1)

      ierr = 0
      errmsg = ' '
      emsg_deep = ' '
      emsg_shal = ' '
      errflg_deep = 0
      errflg_shal = 0

      hvap  = xlv
      cvap  = 1870._RKIND
      cliq  = 4190._RKIND
      fv    = R_v/R_d - 1._RKIND
      eps   = R_d/R_v
      epsm1 = eps - 1._RKIND
      kappa = R_d / cp

      ! TEST VERSION: all caps plus nonzero SAS momentum limiter +/-5.0e-3 m s-2.
      ! GFS SAS control parameters.
      ! Deep and shallow SAS intentionally use different GFS constants.
      ! Deep mass-flux convection:
      !   clam_deep    = 0.1
      !   c0s_deep     = 0.002
      !   c1_deep      = 0.002
      !   betal_deep   = 0.05
      !   betas_deep   = 0.05
      !   evfact_deep  = 0.3
      !   evfactl_deep = 0.3
      !   pgcon_deep   = 0.55
      ! Shallow mass-flux convection:
      !   clam_shal    = 0.3
      !   c0s_shal     = 0.002
      !   c1_shal      = 5.e-4
      !   pgcon_shal   = 0.55
      !   evef_shal    = 0.09 in CCPP/GFS, but this standalone shalcnv_run
      !                  interface does not expose evef; it is kept here as
      !                  documentation for a later shalcnv interface update.
      jcap = 0
      ncloud = 1
      mp_phys    = 0
      mp_phys_mg = -999

      clam_deep    = 0.10_RKIND
      c0_deep      = 0.002_RKIND
      c1_deep      = 0.002_RKIND
      betal_deep   = 0.05_RKIND
      betas_deep   = 0.05_RKIND
      evfact_deep  = 0.30_RKIND
      evfactl_deep = 0.30_RKIND
      pgcon_deep   = 0.55_RKIND

      clam_shal    = 0.30_RKIND
      c0_shal      = 0.002_RKIND
      c1_shal      = 5.0e-4_RKIND
      pgcon_shal   = 0.55_RKIND
      evef_shal    = 0.09_RKIND

      ! Latest SAMF interface controls.  Disable optional aerosol/prognostic
      ! sigma/omega/TKE pathways for this MPAS bridge test; provide valid
      ! arrays anyway because the CCPP routines use assumed-shape dummies.
      nn_samf = 2
      itc_samf = 0
      ntc_samf = 0
      ntk_samf = 0
      ntr_samf = 0
      first_time_step_samf = .false.
      restart_samf = .false.
      hwrf_samfdeep = .false.
      hwrf_samfshal = .false.
      progsigma_samf = .false.
      progomega_samf = .false.
      do_ca = .false.
      ca_closure = .false.
      ca_entr = .false.
      ca_trigger = .false.
      do_mynnedmf = do_mynnedmf_in
      sigmab_coldstart = .false.
      asolfac = 0.958_RKIND
      cscale  = 1.0_RKIND
      nthresh = 0.0_RKIND
      betadcu = 2.0_RKIND
      betamcu = 1.0_RKIND
      betascu = 8.0_RKIND
      cat_adj_deep = 1.0_RKIND
      cat_adj_shal = 1.0_RKIND
      fscav(1) = 0.0_RKIND
      do i = 1, im
         maxMF(i) = max(0.0_RKIND, maxmf_in(i))
      enddo
      ca_deep(:) = 0.0_RKIND
      rainevap(:) = 0.0_RKIND
      tmf(:,:,:) = 0.0_RKIND
      qmicro(:,:) = 0.0_RKIND
      prevsq(:,:) = 0.0_RKIND
      sigmain(:,:) = 0.0_RKIND
      sigmaout(:,:) = 0.0_RKIND
      omegain(:,:) = 0.0_RKIND
      omegaout(:,:) = 0.0_RKIND
      qtr(:,:,:) = 0.0_RKIND
      qtr_sh(:,:,:) = 0.0_RKIND
      dqtr_deep(:,:,:) = 0.0_RKIND
      dqtr_shal(:,:,:) = 0.0_RKIND
      ten_t_deep(:,:) = 0.0_RKIND
      ten_u_deep(:,:) = 0.0_RKIND
      ten_v_deep(:,:) = 0.0_RKIND
      ten_q_deep(:,:,:) = 0.0_RKIND
      ten_t_shal(:,:) = 0.0_RKIND
      ten_u_shal(:,:) = 0.0_RKIND
      ten_v_shal(:,:) = 0.0_RKIND
      ten_q_shal(:,:,:) = 0.0_RKIND
      t1_sh(:,:) = 0.0_RKIND
      q1_sh(:,:) = 0.0_RKIND
      u1_sh(:,:) = 0.0_RKIND
      v1_sh(:,:) = 0.0_RKIND
      do k = 1, km+1
      do i = 1, im
         if (tkeh_in(i,k) == tkeh_in(i,k)) then
            tkeh(i,k) = max(0.0_RKIND, tkeh_in(i,k))
         else
            tkeh(i,k) = 0.0_RKIND
         endif
      enddo
      enddo

      ! Build current-call vertical mapping every call.
      ! Do not carry any k identity from previous calls.
      ! MPAS and GFS SAS/shalcnv are both treated as bottom-up:
      ! k=1 is the lowest model layer, k=km is the model top.
      do i = 1, im
         surf_k(i) = 1
         do k = 1, km
            kmap(i,k) = k
         enddo

         psp(i) = pres_int(i,surf_k(i))

         bad_col(i) = .false.

         if (psp(i) /= psp(i) .or. psp(i) <= 1000._RKIND .or. psp(i) > 120000._RKIND) then
            ierr = 11
            write(errmsg,'(a,i8,1x,es14.6)') 'bad SAS surface pressure psp at i=', i, psp(i)
            return
         endif

         kbot(i) = 0
         ktop(i) = 0
         kcnv(i) = 0
         cldwrk(i) = 0._RKIND
         rn_deep(i) = 0._RKIND
         rn_shal(i) = 0._RKIND

         if (xland(i) < 1.5_RKIND) then
            islimsk(i) = 1       ! land
         else
            islimsk(i) = 0       ! ocean; add sea-ice category later when available
         endif
         garea(i) = max(1.0_RKIND, dx(i)*dx(i))

         ! Surface density for shallow-convection flux conversion.
         ! The driver must pass physical temperature t_phy_p, not potential temperature.
         raw_t = t(i,surf_k(i))
         if (raw_t /= raw_t .or. raw_t < 150._RKIND .or. raw_t > 360._RKIND) then
            bad_col(i) = .true.
            raw_t = min(360._RKIND, max(150._RKIND, raw_t))
         endif
         tv = raw_t * (1._RKIND + fv * max(qmin, qv(i,surf_k(i))))

         rho = pres_mid(i,surf_k(i)) / (R_d * tv)
         if (rho /= rho .or. rho <= 0._RKIND) then
            ierr = 13
            write(errmsg,'(a,i8,3(1x,es14.6))') 'bad SAS surface rho at i,p,tv,rho=', &
                 i, pres_mid(i,surf_k(i)), tv, rho
            return
         endif
         heat(i) = hfx(i) / (rho * cp)
         evap(i) = qfx(i) / rho
      enddo

      do k = 1, km
      do i = 1, im
         ! SAS/shalcnv index k maps to the current MPAS index kk.
         ! kmap is rebuilt every call and preserves bottom-up ordering.
         kk = kmap(i,k)

         ! Robust layer pressure construction.  Some MPAS interface arrays are not
         ! strictly monotonic locally after remapping/physics adjustment, so never
         ! infer bottom/top from the raw index alone.  SAS only needs positive mass
         ! thickness in Pa and a layer-center pressure in Pa.
         pbot = max(pres_int(i,kk), pres_int(i,kk+1))
         ptop = min(pres_int(i,kk), pres_int(i,kk+1))
         delp(i,k)  = pbot - ptop
         prslp(i,k) = 0.5_RKIND * (pbot + ptop)
         phil(i,k)  = gravity * 0.5_RKIND * (z_int(i,kk) + z_int(i,kk+1))

         if (delp(i,k) /= delp(i,k) .or. delp(i,k) <= 0._RKIND) then
            ierr = 13
            write(errmsg,'(a,2i8,1x,l1,3(1x,es14.6))') &
                 'bad SAS delp after pbot/ptop at i,k,pbot,ptop,delp=', &
                 i, k, pbot, ptop, delp(i,k)
            return
         endif

         if (prslp(i,k) /= prslp(i,k) .or. prslp(i,k) <= 0._RKIND .or. prslp(i,k) > 120000._RKIND) then
            ierr = 14
            write(errmsg,'(a,2i8,1x,es14.6)') 'bad SAS prslp at i,k=', i, k, prslp(i,k)
            return
         endif

         ! MPAS driver passes physical temperature.  GFS SAS expects physical T.
         raw_t = t(i,kk)
         if (raw_t /= raw_t .or. raw_t < 150._RKIND .or. raw_t > 360._RKIND) then
            bad_col(i) = .true.
            t1(i,k) = min(360._RKIND, max(150._RKIND, raw_t))
         else
            t1(i,k) = raw_t
         endif

         tv  = t1(i,k) * (1._RKIND + fv * max(qmin, qv(i,kk)))
         rho = prslp(i,k) / (R_d * tv)
         if (rho /= rho .or. rho <= 0._RKIND) then
            ierr = 16
            write(errmsg,'(a,2i8,3(1x,es14.6))') 'bad SAS rho at i,k,p,tv,rho=', &
                 i, k, prslp(i,k), tv, rho
            return
         endif

         dot(i,k) = -rho * gravity * w(i,kk)

         qlc(i,k) = max(0._RKIND, qc(i,kk))
         qli(i,k) = max(0._RKIND, qi(i,kk))
         qlc_old(i,k) = qlc(i,k)
         qli_old(i,k) = qli(i,k)
         qtr(i,k,1) = qli(i,k)
         qtr(i,k,2) = qlc(i,k)
         qtr_sh(i,k,1) = qtr(i,k,1)
         qtr_sh(i,k,2) = qtr(i,k,2)

         q1(i,k) = max(qmin, qv(i,kk))
         u1(i,k) = u(i,kk)
         v1(i,k) = v(i,kk)

         q1_old(i,k) = q1(i,k)
         t1_old(i,k) = t1(i,k)
         u1_old(i,k) = u1(i,k)
         v1_old(i,k) = v1(i,k)

         ud_mf(i,k) = 0._RKIND
         dd_mf(i,k) = 0._RKIND
         dt_mf(i,k) = 0._RKIND
         ud_mf_sh(i,k) = 0._RKIND
         dt_mf_sh(i,k) = 0._RKIND
         cnvw(i,k) = 0._RKIND
         cnvc(i,k) = 0._RKIND

         qlcn(i,k) = 0._RKIND
         qicn(i,k) = 0._RKIND
         w_upi(i,k) = 0._RKIND
         cf_upi(i,k) = 0._RKIND
         cnv_mfd(i,k) = 0._RKIND
         cnv_dqldt(i,k) = 0._RKIND
         clcn(i,k) = 0._RKIND
         cnv_fice(i,k) = 0._RKIND
         cnv_ndrop(i,k) = 0._RKIND
         cnv_nice(i,k) = 0._RKIND
      enddo
      enddo

      !-----------------------------------------------------------------
      ! Latest CCPP SAMF deep convection interface.
      ! The latest scheme returns tendencies (ten_t/ten_q/ten_u/ten_v)
      ! rather than updating t1/q1/u1/v1 in place.
      !-----------------------------------------------------------------
      call samfdeepcnv_run(                                                   &
           im, km, nn_samf, first_time_step_samf, restart_samf,               &
           tmf, qmicro, itc_samf, ntc_samf, cliq, cp, cvap,                  &
           eps, epsm1, fv, gravity, hvap, R_d, R_v,                           &
           t0, dt, ntk_samf, ntr_samf, delp,                                  &
           ten_t_deep, ten_u_deep, ten_v_deep, ten_q_deep,                    &
           prslp, psp, phil, tkeh, qtr, dqtr_deep, prevsq,                    &
           q1, q1, t1, u1, v1, fscav,                                         &
           hwrf_samfdeep, progsigma_samf, progomega_samf,                    &
           cldwrk, rn_deep, kbot, ktop, kcnv,                                 &
           islimsk, garea, dot, ncloud, hpbl, ud_mf, dd_mf, dt_mf,            &
           cnvw, cnvc, qlcn, qicn, w_upi, cf_upi, cnv_mfd,                   &
           cnv_dqldt, clcn, cnv_fice, cnv_ndrop, cnv_nice,                   &
           mp_phys, mp_phys_mg,                                               &
           clam_deep, c0_deep, c1_deep, betal_deep, betas_deep,               &
           evfact_deep, pgcon_deep, asolfac, cscale,                         &
           do_ca, ca_closure, ca_entr, ca_trigger, nthresh, ca_deep,          &
           rainevap, sigmain, sigmaout, omegain, omegaout,                    &
           betadcu, betamcu, betascu, maxMF,                                  &
           do_mynnedmf, sigmab_coldstart, cat_adj_deep,                      &
           emsg_deep, errflg_deep)

      if (errflg_deep /= 0) then
         ierr = errflg_deep
         errmsg = 'samfdeepcnv_run failed: '//trim(emsg_deep)
         return
      endif

      ! Stage the input to shallow convection with deep-convection tendencies,
      ! matching the CCPP sequence behavior where shallow sees the updated state.
      do k = 1, km
      do i = 1, im
         t1_sh(i,k) = t1(i,k) + dt * ten_t_deep(i,k)
         q1_sh(i,k) = max(qmin, q1(i,k) + dt * ten_q_deep(i,k,1))
         u1_sh(i,k) = u1(i,k) + dt * ten_u_deep(i,k)
         v1_sh(i,k) = v1(i,k) + dt * ten_v_deep(i,k)
         qtr_sh(i,k,1) = qtr(i,k,1) + dt * dqtr_deep(i,k,1)
         qtr_sh(i,k,2) = qtr(i,k,2) + dt * dqtr_deep(i,k,2)
      enddo
      enddo

      call samfshalcnv_run(                                                  &
           im, km, nn_samf, itc_samf, ntc_samf, cliq, cp, cvap,             &
           eps, epsm1, fv, gravity, hvap, R_d, R_v,                          &
           t0, dt, ntk_samf, ntr_samf, delp,                                 &
           first_time_step_samf, restart_samf,                               &
           tmf, qmicro, progsigma_samf, progomega_samf,                     &
           prslp, psp, phil, tkeh, qtr_sh, dqtr_shal, prevsq,                &
           q1_sh, q1_sh, t1_sh, u1_sh, v1_sh, fscav,                         &
           rn_shal, kbot, ktop, kcnv, islimsk, garea, cscale,                &
           ten_t_shal, ten_u_shal, ten_v_shal, ten_q_shal,                   &
           dot, ncloud, hpbl, ud_mf_sh, dt_mf_sh, cnvw, cnvc,                &
           clam_shal, c0_shal, c1_shal, evef_shal, pgcon_shal,               &
           asolfac, hwrf_samfshal,                                           &
           sigmain, sigmaout, omegain, omegaout,                             &
           betadcu, betamcu, betascu, cat_adj_shal,                          &
           emsg_shal, errflg_shal)

      if (errflg_shal /= 0) then
         ierr = errflg_shal
         errmsg = 'samfshalcnv_run failed: '//trim(emsg_shal)
         return
      endif

      do k = 1, km
      do i = 1, im
         kk = kmap(i,k)

         if (bad_col(i)) then
            rthcuten(i,kk) = 0._RKIND
            rqvcuten(i,kk) = 0._RKIND
            rqccuten(i,kk) = 0._RKIND
            rqicuten(i,kk) = 0._RKIND
            rucuten(i,kk)  = 0._RKIND
            rvcuten(i,kk)  = 0._RKIND
         else
            theta_fac = (100000._RKIND / prslp(i,k)) ** kappa

            rth_test = (ten_t_deep(i,k) + ten_t_shal(i,k)) * theta_fac
            rqv_raw  =  ten_q_deep(i,k,1) + ten_q_shal(i,k,1)
            rqc_raw  =  dqtr_deep(i,k,2) + dqtr_shal(i,k,2)
            rqi_raw  =  dqtr_deep(i,k,1) + dqtr_shal(i,k,1)
            ru_raw   =  ten_u_deep(i,k) + ten_u_shal(i,k)
            rv_raw   =  ten_v_deep(i,k) + ten_v_shal(i,k)

            if (rth_test /= rth_test .or. rqv_raw /= rqv_raw .or.             &
                rqc_raw /= rqc_raw .or. rqi_raw /= rqi_raw .or.              &
                ru_raw /= ru_raw .or. rv_raw /= rv_raw) then
               write(0,*) 'NaN tendency from latest SAMF GFS convection'
               write(0,*) 'i,k,kk=', i,k,kk
               call flush(0)
               ierr = 99
               errmsg = 'NaN tendency from latest SAMF GFS convection'
               return
            endif

            if (abs(rth_test) > 5.0e-2_RKIND .or. abs(rqv_raw) > 2.0e-5_RKIND .or. &
                abs(rqc_raw) > 2.0e-5_RKIND .or. abs(rqi_raw) > 2.0e-5_RKIND .or. &
                abs(ru_raw)  > 5.0e-3_RKIND .or. abs(rv_raw)  > 5.0e-3_RKIND) then
               write(0,*) 'LARGE LATEST SAMF GFS SAS TENDENCY -- REPORT ONLY'
               write(0,*) 'i,k,kk             = ', i, k, kk
               write(0,*) 'rth,rqv,rqc,rqi    = ', rth_test, rqv_raw, rqc_raw, rqi_raw
               write(0,*) 'ru,rv              = ', ru_raw, rv_raw
               write(0,*) 'deep rn,shal rn    = ', rn_deep(i), rn_shal(i)
               write(0,*) 'kbot,ktop,kcnv     = ', kbot(i), ktop(i), kcnv(i)
               call flush(0)
            endif

            rthcuten(i,kk) = rth_test
            rqvcuten(i,kk) = rqv_raw
            rqccuten(i,kk) = rqc_raw
            rqicuten(i,kk) = rqi_raw
            rucuten(i,kk)  = ru_raw
            rvcuten(i,kk)  = rv_raw
         endif

         if (rthcuten(i,kk) /= rthcuten(i,kk)) rthcuten(i,kk) = 0._RKIND
         if (rqvcuten(i,kk) /= rqvcuten(i,kk)) rqvcuten(i,kk) = 0._RKIND
         if (rqccuten(i,kk) /= rqccuten(i,kk)) rqccuten(i,kk) = 0._RKIND
         if (rqicuten(i,kk) /= rqicuten(i,kk)) rqicuten(i,kk) = 0._RKIND
         if (rucuten(i,kk)  /= rucuten(i,kk) ) rucuten(i,kk)  = 0._RKIND
         if (rvcuten(i,kk)  /= rvcuten(i,kk) ) rvcuten(i,kk)  = 0._RKIND
      enddo
      enddo

      do i = 1, im
         if (rn_deep(i) /= rn_deep(i)) rn_deep(i) = 0._RKIND
         if (rn_shal(i) /= rn_shal(i)) rn_shal(i) = 0._RKIND
         rn_deep(i) = max(0._RKIND, rn_deep(i))
         rn_shal(i) = max(0._RKIND, rn_shal(i))

         if (bad_col(i)) then
            rn_deep(i) = 0._RKIND
            rn_shal(i) = 0._RKIND
            raincv(i) = 0._RKIND
            pratec(i) = 0._RKIND
            kbot(i) = 0
            ktop(i) = 0
         else
            raincv(i) = raincv(i) + 1000._RKIND * (rn_deep(i) + rn_shal(i))
            pratec(i) = 1000._RKIND * (rn_deep(i) + rn_shal(i)) / dt
         endif

         if (kbot(i) >= 1 .and. kbot(i) <= km) then
            cubot(i) = real(kmap(i,kbot(i)), RKIND)
         else
            cubot(i) = real(km+1, RKIND)
         endif

         if (ktop(i) >= 1 .and. ktop(i) <= km) then
            cutop(i) = real(kmap(i,ktop(i)), RKIND)
         else
            cutop(i) = 1._RKIND
         endif
      enddo

   end subroutine mpas_call_gfs_convection

end module mpas_gfs_convection_wrapper_mod
