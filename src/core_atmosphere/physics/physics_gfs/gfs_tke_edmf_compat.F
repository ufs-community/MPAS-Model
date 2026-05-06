!> \file gfs_tke_edmf_compat.F
!! Minimal compatibility layer for GFS/UFS TKE-EDMF routines in MPAS.
!!
!! Fixed-form safe version.  Keep lines short to avoid 72-column truncation.

      module gfs_tke_edmf_machine
      use mpas_kind_types, only: RKIND
      implicit none
      integer, parameter :: kind_phys = RKIND
      end module gfs_tke_edmf_machine

      module physcons
      use gfs_tke_edmf_machine, only: kind_phys
      implicit none
      real(kind=kind_phys), parameter :: con_g = 9.80665_kind_phys
      real(kind=kind_phys), parameter :: con_cp = 1004.0_kind_phys
      real(kind=kind_phys), parameter :: con_rd = 287.0_kind_phys
      real(kind=kind_phys), parameter :: con_rv = 461.5_kind_phys
      real(kind=kind_phys), parameter :: con_hvap = 2.5e6_kind_phys
      real(kind=kind_phys), parameter :: con_hfus =
     &     3.3358e5_kind_phys
      real(kind=kind_phys), parameter :: con_eps = con_rd / con_rv
      real(kind=kind_phys), parameter :: con_epsm1 =
     &     con_eps - 1.0_kind_phys
      real(kind=kind_phys), parameter :: con_fvirt =
     &     con_rv / con_rd - 1.0_kind_phys
      end module physcons

      module funcphys
      use gfs_tke_edmf_machine, only: kind_phys
      implicit none

      contains

      real(kind=kind_phys) function fpvs(t)
      real(kind=kind_phys), intent(in) :: t
      real(kind=kind_phys) :: tc
      real(kind=kind_phys) :: a1
      real(kind=kind_phys) :: a2
      real(kind=kind_phys) :: a3

      tc = t - 273.15_kind_phys
      a1 = 611.2_kind_phys
      a2 = 17.67_kind_phys * tc
      a3 = tc + 243.5_kind_phys
      fpvs = a1 * exp(a2 / a3)

      end function fpvs

      end module funcphys
