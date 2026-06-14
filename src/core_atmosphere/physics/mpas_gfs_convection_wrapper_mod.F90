!> MPAS wrapper for GFS SAS deep convection (sascnvn) + GFS shallow convection
!> (shalcnv).  This wrapper is intended to be called from
!> mpas_atmphys_driver_convection.F after convection_from_MPAS has filled the
!> *_p work arrays, and before convection_to_MPAS copies tendencies back.
!
! Integration notes:
!  1. Add these uses in mpas_atmphys_driver_convection.F:
!        use mpas_gfs_convection_wrapper_mod, only: mpas_call_gfs_convection
!     and add the GFS modules if they are module-wrapped in your tree, e.g.:
!        use sascnvn, only: sascnvn_run
!        use shalcnv, only: shalcnv_run
!  2. Add a convection case, e.g. case("cu_gfs_sas"), and call this wrapper.
!  3. Add cu_gfs_sas to Registry.xml package/namelist checks and to
!     physics_namelist_check legal values.
!
module mpas_gfs_convection_wrapper_mod

   use mpas_kind_types,       only: RKIND
   use mpas_atmphys_constants, only: gravity, cp, R_d, R_v, xlv

   implicit none
   private
   public :: mpas_call_gfs_convection

   ! MPAS constants module does not export t0 in this code base.
   ! GFS SAS expects the freezing-point reference temperature.
   real(kind=RKIND), parameter :: t0 = 273.15_RKIND

contains

   subroutine mpas_call_gfs_convection(im, km, dt,                              &
        pres_mid, pres_int, z_int, u, v, t, qv, qc, qi, w,                      &
        hpbl, hfx, qfx, xland,                                                  &
        rthcuten, rqvcuten, rqccuten, rqicuten, rucuten, rvcuten,               &
        raincv, pratec, cubot, cutop, ierr, errmsg)

      ! GFS routines.  If your GFS physics source files are not module files,
      ! put explicit interfaces here or module-wrap sascnvn.F and shalcnv.F.
      use sascnvn, only: sascnvn_run
      use shalcnv, only: shalcnv_run

      implicit none

      integer, intent(in) :: im, km
      real(kind=RKIND), intent(in) :: dt

      ! MPAS/WRF-style arrays: horizontal x vertical.
      real(kind=RKIND), intent(in)    :: pres_mid(im,km)     ! Pa, layer pressure
      real(kind=RKIND), intent(in)    :: pres_int(im,km+1)   ! Pa, interface pressure
      real(kind=RKIND), intent(in)    :: z_int(im,km+1)      ! m, interface height
      real(kind=RKIND), intent(inout) :: u(im,km), v(im,km)  ! m s-1
      real(kind=RKIND), intent(inout) :: t(im,km)            ! K
      real(kind=RKIND), intent(inout) :: qv(im,km)           ! kg kg-1
      real(kind=RKIND), intent(inout) :: qc(im,km), qi(im,km)! kg kg-1
      real(kind=RKIND), intent(in)    :: w(im,km)            ! m s-1, approximate
      real(kind=RKIND), intent(in)    :: hpbl(im)            ! m
      real(kind=RKIND), intent(in)    :: hfx(im), qfx(im)    ! W m-2, kg m-2 s-1
      real(kind=RKIND), intent(in)    :: xland(im)           ! MPAS/WRF: 1 land, 2 water

      ! Tendencies returned to MPAS, units per second.
      real(kind=RKIND), intent(inout) :: rthcuten(im,km)
      real(kind=RKIND), intent(inout) :: rqvcuten(im,km)
      real(kind=RKIND), intent(inout) :: rqccuten(im,km)
      real(kind=RKIND), intent(inout) :: rqicuten(im,km)
      real(kind=RKIND), intent(inout) :: rucuten(im,km)
      real(kind=RKIND), intent(inout) :: rvcuten(im,km)
      real(kind=RKIND), intent(inout) :: raincv(im)          ! mm over dt
      real(kind=RKIND), intent(inout) :: pratec(im)          ! mm s-1
      real(kind=RKIND), intent(inout) :: cubot(im), cutop(im)

      integer, intent(out) :: ierr
      character(len=*), intent(out) :: errmsg

      integer :: i, k
      integer :: jcap, ncloud, mp_phys, mp_phys_mg
      integer :: kbot(im), ktop(im), kcnv(im), islimsk(im)
      integer :: errflg
      real(kind=RKIND) :: hvap, cvap, cliq, fv, eps, epsm1
      real(kind=RKIND) :: clam, c0, c1, betal, betas, evfact, evfactl, pgcon
      real(kind=RKIND) :: psp(im), delp(im,km), prslp(im,km), phil(im,km)
      real(kind=RKIND) :: dot(im,km), q1_old(im,km), t1_old(im,km)
      real(kind=RKIND) :: u1_old(im,km), v1_old(im,km)
      real(kind=RKIND) :: qlc(im,km), qli(im,km), cnvw(im,km), cnvc(im,km)
      real(kind=RKIND) :: cldwrk(im), rn_deep(im), rn_shal(im)
      real(kind=RKIND) :: ud_mf(im,km), dd_mf(im,km), dt_mf(im,km)
      real(kind=RKIND) :: ud_mf_sh(im,km), dt_mf_sh(im,km)
      real(kind=RKIND) :: heat(im), evap(im)
      character(len=256) :: emsg

      ierr = 0
      errmsg = ' '

      hvap  = xlv
      cvap  = 1870._RKIND
      cliq  = 4190._RKIND
      fv    = R_v/R_d - 1._RKIND
      eps   = R_d/R_v
      epsm1 = eps - 1._RKIND

      ! Conservative default SAS tunables from GFS-style interface.  Move these
      ! to namelist after the first clean build/science test.
      jcap    = 0
      ncloud  = 1
      mp_phys    = 0
      mp_phys_mg = -999
      clam    = 0.10_RKIND
      c0      = 0.002_RKIND
      c1      = 0.002_RKIND
      betal   = 0.05_RKIND
      betas   = 0.05_RKIND
      evfact  = 0.30_RKIND
      evfactl = 0.30_RKIND
      pgcon   = 0.55_RKIND

      do i = 1, im
         psp(i) = pres_int(i,1)
         kbot(i) = 0
         ktop(i) = 0
         kcnv(i) = 0
         islimsk(i) = merge(1, 0, xland(i) < 1.5_RKIND)  ! GFS: 1 land, 0 sea
         heat(i) = hfx(i)
         evap(i) = qfx(i)
         rn_deep(i) = 0._RKIND
         rn_shal(i) = 0._RKIND
         cldwrk(i) = 0._RKIND
      end do

      do k = 1, km
      do i = 1, im
         delp(i,k)  = max(1._RKIND, pres_int(i,k) - pres_int(i,k+1))
         prslp(i,k) = pres_mid(i,k)
         phil(i,k)  = gravity * 0.5_RKIND * (z_int(i,k) + z_int(i,k+1))
         dot(i,k)   = -w(i,k) * delp(i,k) / max(1._RKIND, z_int(i,k+1)-z_int(i,k))

         qlc(i,k) = max(0._RKIND, qc(i,k))
         qli(i,k) = max(0._RKIND, qi(i,k))
         cnvw(i,k) = 0._RKIND
         cnvc(i,k) = 0._RKIND

         q1_old(i,k) = qv(i,k)
         t1_old(i,k) = t(i,k)
         u1_old(i,k) = u(i,k)
         v1_old(i,k) = v(i,k)

         ud_mf(i,k) = 0._RKIND
         dd_mf(i,k) = 0._RKIND
         dt_mf(i,k) = 0._RKIND
         ud_mf_sh(i,k) = 0._RKIND
         dt_mf_sh(i,k) = 0._RKIND
      end do
      end do
      call sascnvn_run(gravity, cp, hvap, R_v, fv, t0, R_d, cvap, cliq, eps, epsm1, &
           im, km, jcap, dt, delp, prslp, psp, phil, qlc, qli,                    &
           qv, t, u, v, cldwrk, rn_deep, kbot, ktop, kcnv, islimsk,                &
           dot, ncloud, ud_mf, dd_mf, dt_mf, cnvw, cnvc,                           &
           mp_phys=mp_phys, mp_phys_mg=mp_phys_mg,                                &
           clam=clam, c0=c0, c1=c1, betal=betal, betas=betas,                      &
           evfact=evfact, evfactl=evfactl, pgcon=pgcon,                            &
           errmsg=emsg, errflg=errflg)
      if (errflg /= 0) then
         ierr = errflg
         errmsg = 'sascnvn_run failed: '//trim(emsg)
         return
      end if

      call shalcnv_run(gravity, cp, hvap, R_v, fv, t0, R_d, cvap, cliq, eps, epsm1, &
           im, km, jcap, dt, delp, prslp, psp, phil, qlc, qli,                    &
           qv, t, u, v, rn_shal, kbot, ktop, kcnv, islimsk,                       &
           dot, ncloud, hpbl, heat, evap, ud_mf_sh, dt_mf_sh, cnvw, cnvc,          &
           clam, c0, c1, pgcon, emsg, errflg)
      if (errflg /= 0) then
         ierr = errflg
         errmsg = 'shalcnv_run failed: '//trim(emsg)
         return
      end if

      do k = 1, km
      do i = 1, im
         rthcuten(i,k) = rthcuten(i,k) + (t(i,k)  - t1_old(i,k)) / dt
         rqvcuten(i,k) = rqvcuten(i,k) + (qv(i,k) - q1_old(i,k)) / dt
         rucuten(i,k)  = rucuten(i,k)  + (u(i,k)  - u1_old(i,k)) / dt
         rvcuten(i,k)  = rvcuten(i,k)  + (v(i,k)  - v1_old(i,k)) / dt
         rqccuten(i,k) = rqccuten(i,k) + (qlc(i,k) - qc(i,k)) / dt
         rqicuten(i,k) = rqicuten(i,k) + (qli(i,k) - qi(i,k)) / dt
      end do
      end do

      do i = 1, im
         raincv(i) = raincv(i) + rn_deep(i) + rn_shal(i)
         pratec(i) = (rn_deep(i) + rn_shal(i)) / dt
         cubot(i) = real(kbot(i), RKIND)
         cutop(i) = real(ktop(i), RKIND)
      end do

   end subroutine mpas_call_gfs_convection

end module mpas_gfs_convection_wrapper_mod
