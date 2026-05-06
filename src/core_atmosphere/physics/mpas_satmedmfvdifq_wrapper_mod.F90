module mpas_satmedmfvdifq_wrapper_mod

  use mpas_kind_types, only: RKIND
  use satmedmfvdifq, only: satmedmfvdifq_run
  implicit none

  type mpas_satmedmfvdifq_config_type
    logical :: sa3dtke      = .false.
    logical :: tte_edmf     = .false.
    logical :: dspheat      = .false.
    logical :: use_oceanuv  = .false.
    logical :: do_canopy    = .false.
    logical :: cplaqm       = .false.
    logical :: gen_tend     = .false.
    logical :: ldiag3d      = .false.

    integer :: sfc_rlm = 0
    integer :: tc_pbl  = 0
    integer :: use_lpt = 0

    real(kind=RKIND) :: xkzm_m = 1.0_RKIND
    real(kind=RKIND) :: xkzm_h = 1.0_RKIND
    real(kind=RKIND) :: xkzm_s = 0.7_RKIND
    real(kind=RKIND) :: dspfac = 1.0_RKIND
    real(kind=RKIND) :: bl_upfr = 1.0_RKIND
    real(kind=RKIND) :: bl_dnfr = 1.0_RKIND
    real(kind=RKIND) :: rlmx = 300.0_RKIND
    real(kind=RKIND) :: elmx = 300.0_RKIND
  end type mpas_satmedmfvdifq_config_type

contains

  subroutine mpas_call_satmedmfvdifq(nCells, nVertLevels, ntrac, dt,      &
                                     z_mid, z_int, areaCell,             &
                                     u_mpas, v_mpas, t_mpas,             &
                                     qv_mpas, qc_mpas, qi_mpas, tke_mpas,&
                                     p_mid, p_int, exner_mid,            &
                                     sw_heat, lw_heat, coszen,           &
                                     skin_temp, shflx, lhflx, stress_in, &
                                     z0_mpas, u10, v10, vegfra_in, rb_in, fm_in, fh_in,    &
                                     hpbl_out, kpbl_out,                 &
                                     ten_t_out, ten_u_out, ten_v_out,    &
                                     cfg, errmsg, errflg)

    integer, intent(in) :: nCells, nVertLevels, ntrac
    real(kind=RKIND), intent(in) :: dt
    real(kind=RKIND), intent(in) :: z_mid(nCells,nVertLevels)
    real(kind=RKIND), intent(in) :: z_int(nCells,nVertLevels+1)
    real(kind=RKIND), intent(in) :: areaCell(nCells)
    real(kind=RKIND), intent(in) :: u_mpas(nCells,nVertLevels)
    real(kind=RKIND), intent(in) :: v_mpas(nCells,nVertLevels)
    real(kind=RKIND), intent(in) :: t_mpas(nCells,nVertLevels)
    real(kind=RKIND), intent(in) :: qv_mpas(nCells,nVertLevels)
    real(kind=RKIND), intent(in) :: qc_mpas(nCells,nVertLevels)
    real(kind=RKIND), intent(in) :: qi_mpas(nCells,nVertLevels)
    real(kind=RKIND), intent(inout) :: tke_mpas(nCells,nVertLevels)
    real(kind=RKIND), intent(in) :: p_mid(nCells,nVertLevels)
    real(kind=RKIND), intent(in) :: p_int(nCells,nVertLevels+1)
    real(kind=RKIND), intent(in) :: exner_mid(nCells,nVertLevels)
    real(kind=RKIND), intent(in) :: sw_heat(nCells,nVertLevels)
    real(kind=RKIND), intent(in) :: lw_heat(nCells,nVertLevels)
    real(kind=RKIND), intent(in) :: coszen(nCells)
    real(kind=RKIND), intent(in) :: skin_temp(nCells)
    real(kind=RKIND), intent(in) :: shflx(nCells)
    real(kind=RKIND), intent(in) :: lhflx(nCells)
    real(kind=RKIND), intent(in) :: stress_in(nCells)
    real(kind=RKIND), intent(in) :: z0_mpas(nCells)
    real(kind=RKIND), intent(in) :: u10(nCells), v10(nCells)
    real(kind=RKIND), intent(in) :: vegfra_in(nCells),rb_in(nCells)
    real(kind=RKIND), intent(in) :: fm_in(nCells), fh_in(nCells)
    real(kind=RKIND), intent(out) :: hpbl_out(nCells)
    integer, intent(out) :: kpbl_out(nCells)
    real(kind=RKIND), intent(out) :: ten_t_out(nCells,nVertLevels)
    real(kind=RKIND), intent(out) :: ten_u_out(nCells,nVertLevels)
    real(kind=RKIND), intent(out) :: ten_v_out(nCells,nVertLevels)
    type(mpas_satmedmfvdifq_config_type), intent(in) :: cfg
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    integer :: im, km
    integer :: i, k
    integer :: ntqv, ntcw, ntiw, ntrw, ntke
    integer :: index_of_temperature, index_of_x_wind
    integer :: index_of_y_wind, index_of_process_pbl

    real(kind=RKIND), parameter :: grav  = 9.80665_RKIND
    real(kind=RKIND), parameter :: pi    = 3.14159265358979323846_RKIND
    real(kind=RKIND), parameter :: rd    = 287.0_RKIND
    real(kind=RKIND), parameter :: cp    = 1004.0_RKIND
    real(kind=RKIND), parameter :: rv    = 461.5_RKIND
    real(kind=RKIND), parameter :: hvap  = 2.5e6_RKIND
    real(kind=RKIND), parameter :: hfus  = 3.3358e5_RKIND
    real(kind=RKIND), parameter :: fv    = rv/rd - 1.0_RKIND
    real(kind=RKIND), parameter :: eps   = rd/rv
    real(kind=RKIND), parameter :: epsm1 = eps - 1.0_RKIND
    real(kind=RKIND), parameter :: z0lo = 0.1_RKIND
    real(kind=RKIND), parameter :: z0up = 1.0_RKIND

    real(kind=RKIND), allocatable :: rtg(:,:,:), q1(:,:,:)
    real(kind=RKIND), allocatable :: u1(:,:), v1(:,:), t1(:,:)
    real(kind=RKIND), allocatable :: swh(:,:), hlw(:,:), xmu(:)
    real(kind=RKIND), allocatable :: garea(:), zvfun(:), sigmaf(:)
    real(kind=RKIND), allocatable :: psk(:), rbsoil(:), zorl(:)
    real(kind=RKIND), allocatable :: u10m(:), v10m(:), fm(:), fh(:)
    real(kind=RKIND), allocatable :: tsea(:), heat(:), evap(:)
    real(kind=RKIND), allocatable :: stress(:), spd1(:)
    real(kind=RKIND), allocatable :: prsi(:,:), del(:,:), prsl(:,:), prslk(:,:)
    real(kind=RKIND), allocatable :: phii(:,:), phil(:,:)
    real(kind=RKIND), allocatable :: dusfc(:), dvsfc(:), dtsfc(:), dqsfc(:)
    real(kind=RKIND), allocatable :: hpbl(:), dkt(:,:), dku(:,:), tkeh(:,:)
    real(kind=RKIND), allocatable :: def_1(:,:), def_2(:,:), def_3(:,:)
    real(kind=RKIND), allocatable :: dku3d_h(:,:), dku3d_e(:,:)
    real(kind=RKIND), allocatable :: ten_t(:,:), ten_u(:,:), ten_v(:,:)
    real(kind=RKIND), allocatable :: dv(:,:), du(:,:), tdt(:,:)
    real(kind=RKIND), allocatable :: dtend(:,:,:)
    real(kind=RKIND), allocatable :: claie(:), cfch(:), cfrt(:), cclu(:), cpopu(:)
    real(kind=RKIND) :: tem1, tem2
    integer, allocatable :: kpbl(:), kinver(:), dtidx(:,:)

    im = nCells
    km = nVertLevels

    errmsg = ''
    errflg = 0

    ! Tracer layout for the UFS routine.
    ntqv = 1
    ntcw = 2
    ntiw = 3
    ntrw = 0
    ntke = ntrac

    if (ntrac < 4) then
      errmsg = 'mpas_call_satmedmfvdifq: ntrac must be at least 4: qv,qc,qi,tke'
      errflg = 1
      return
    endif

    allocate(rtg(im,km,ntrac), q1(im,km,ntrac))
    allocate(u1(im,km), v1(im,km), t1(im,km))
    allocate(swh(im,km), hlw(im,km), xmu(im))
    allocate(garea(im), zvfun(im), sigmaf(im))
    allocate(psk(im), rbsoil(im), zorl(im))
    allocate(u10m(im), v10m(im), fm(im), fh(im))
    allocate(tsea(im), heat(im), evap(im), stress(im), spd1(im))
    allocate(prsi(im,km+1), del(im,km), prsl(im,km), prslk(im,km))
    allocate(phii(im,km+1), phil(im,km))
    allocate(dusfc(im), dvsfc(im), dtsfc(im), dqsfc(im))
    allocate(hpbl(im), dkt(im,km), dku(im,km), tkeh(im,km))
    allocate(def_1(im,km), def_2(im,km), def_3(im,km))
    allocate(dku3d_h(im,km), dku3d_e(im,km))
    allocate(ten_t(im,km), ten_u(im,km), ten_v(im,km))
    allocate(dv(im,km), du(im,km), tdt(im,km))
    allocate(kpbl(im), kinver(im))
    allocate(dtidx(3,1), dtend(im,km,3))
    allocate(claie(im), cfch(im), cfrt(im), cclu(im), cpopu(im))

    rtg = 0.0_RKIND
    q1 = 0.0_RKIND
    dkt = 0.0_RKIND
    dku = 0.0_RKIND
    tkeh = 0.0_RKIND
    def_1 = 0.0_RKIND
    def_2 = 0.0_RKIND
    def_3 = 0.0_RKIND
    dku3d_h = 0.0_RKIND
    dku3d_e = 0.0_RKIND
    ten_t = 0.0_RKIND
    ten_u = 0.0_RKIND
    ten_v = 0.0_RKIND
    dv = 0.0_RKIND
    du = 0.0_RKIND
    tdt = 0.0_RKIND
    dtend = 0.0_RKIND
    dtidx = 0
    claie = 0.0_RKIND
    cfch  = 0.0_RKIND
    cfrt  = 0.0_RKIND
    cclu  = 0.0_RKIND
    cpopu = 0.0_RKIND

    do k = 1, km
      do i = 1, im
        u1(i,k) = u_mpas(i,k)
        v1(i,k) = v_mpas(i,k)
        t1(i,k) = t_mpas(i,k)

        q1(i,k,ntqv) = max(qv_mpas(i,k), 1.0e-12_RKIND)
        q1(i,k,ntcw) = max(qc_mpas(i,k), 0.0_RKIND)
        q1(i,k,ntiw) = max(qi_mpas(i,k), 0.0_RKIND)
        q1(i,k,ntke) = max(tke_mpas(i,k), 1.0e-9_RKIND)

        prsl(i,k)  = p_mid(i,k)
        prslk(i,k) = exner_mid(i,k)

        ! UFS routine expects geopotential, not geometric height.
        phil(i,k) = grav * z_mid(i,k)

        swh(i,k) = sw_heat(i,k)
        hlw(i,k) = lw_heat(i,k)
      enddo
    enddo

    do k = 1, km+1
      do i = 1, im
        prsi(i,k) = p_int(i,k)

        ! UFS routine internally does zi=phii/grav, so pass geopotential.
        phii(i,k) = grav * z_int(i,k)
      enddo
    enddo

    do k = 1, km
      do i = 1, im
        del(i,k) = abs(prsi(i,k) - prsi(i,k+1))
      enddo
    enddo

    do i = 1, im
      garea(i) = areaCell(i)

      xmu(i) = max(coszen(i), 0.0_RKIND)

      ! UFS code uses z0 = 0.01*zorl, so zorl is in cm.
      zorl(i) = max(z0_mpas(i), 1.0e-6_RKIND) * 100.0_RKIND

      tsea(i) = skin_temp(i)
      heat(i) = shflx(i)

      ! If MPAS gives latent heat flux W m-2, convert to kg m-2 s-1.
      evap(i) = lhflx(i) / hvap

      stress(i) = max(stress_in(i), 0.0_RKIND)
      spd1(i) = max(sqrt(u1(i,1)**2 + v1(i,1)**2), 0.1_RKIND)

      u10m(i) = u10(i)
      v10m(i) = v10(i)
      rbsoil(i) = rb_in(i)
      fm(i) = max(fm_in(i), 1.0e-6_RKIND)
      fh(i) = max(fh_in(i), 1.0e-6_RKIND)

      ! If MPAS does not have rbsoil, start neutral.
!     rbsoil(i) = 0.0_RKIND

      ! psk is surface Exner. Use lowest layer as fallback.
      psk(i) = prslk(i,1)

      ! Need real mappings later from vegetation/roughness data.
      sigmaf(i) = vegfra_in(i)
!     zvfun(i) = 1.0_RKIND
      !-------------------------------------------------------
! Compute zvfun: function of surface roughness and vegetation
!-------------------------------------------------------

! z0 from MPAS is in meters; UFS uses zorl in cm
! Here we stay consistent with MPAS units (meters)
      tem1 = (z0_mpas(i) - z0lo) / (z0up - z0lo)
! limit between 0 and 1
      tem1 = min(max(tem1, 0.0_RKIND), 1.0_RKIND)
! ensure minimum vegetation fraction
      tem2 = max(sigmaf(i), 0.1_RKIND)
! final function
      zvfun(i) = sqrt(tem1 * tem2)

      ! No inversion limiter initially.
      kinver(i) = km

      kpbl(i) = 1
      hpbl(i) = 0.0_RKIND
    enddo

    index_of_temperature = 1
    index_of_x_wind = 2
    index_of_y_wind = 3
    index_of_process_pbl = 1

    call satmedmfvdifq_run(im, km, ntrac, ntcw, ntrw, ntiw, ntke,       &
         grav, pi, rd, cp, rv, hvap, hfus, fv, eps, epsm1,             &
         def_1, def_2, def_3, cfg%sa3dtke, dku3d_h, dku3d_e,           &
         dv, du, tdt, rtg, u1, v1, t1, q1,                             &
         swh, hlw, xmu, garea, zvfun, sigmaf,                          &
         psk, rbsoil, zorl, u10m, v10m, fm, fh,                        &
         tsea, heat, evap, stress, spd1, kpbl,                         &
         prsi, del, prsl, prslk, phii, phil, dt, cfg%tte_edmf,         &
         cfg%dspheat, dusfc, dvsfc, dtsfc, dqsfc, hpbl, dkt, dku, tkeh,&
         kinver, cfg%xkzm_m, cfg%xkzm_h, cfg%xkzm_s, cfg%dspfac,       &
         cfg%bl_upfr, cfg%bl_dnfr, cfg%rlmx, cfg%elmx,                 &
         cfg%sfc_rlm, cfg%tc_pbl, cfg%use_lpt,                         &
         cfg%do_canopy, cfg%cplaqm, claie, cfch, cfrt, cclu, cpopu,    &
         ntqv, dtend, dtidx, index_of_temperature,                     &
         index_of_x_wind, index_of_y_wind, index_of_process_pbl,       &
         cfg%gen_tend, cfg%ldiag3d, errmsg, errflg)

    if (errflg /= 0) return

    do k = 1, km
      do i = 1, im
        ten_t_out(i,k) = tdt(i,k)
        ten_u_out(i,k) = du(i,k)
        ten_v_out(i,k) = dv(i,k)

        tke_mpas(i,k) = max(q1(i,k,ntke), 0.0_RKIND)
      enddo
    enddo

    do i = 1, im
      hpbl_out(i) = hpbl(i)
      kpbl_out(i) = kpbl(i)
    enddo

  end subroutine mpas_call_satmedmfvdifq

end module mpas_satmedmfvdifq_wrapper_mod
