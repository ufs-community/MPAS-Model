!> \file gfs_tke_edmf_compat.F
!! Minimal compatibility layer for using selected GFS/UFS TKE-EDMF routines
!! inside the MPAS physics build.
!!
!! This file supplies the small GFS-style modules referenced by the imported
!! routines: machine, physcons, and funcphys.  It should be used only if these
!! modules are not already provided by the MPAS build.

      module machine
        use mpas_kind_types, only: RKIND
        implicit none
        integer, parameter :: kind_phys = RKIND
      end module machine

      module physcons
        use machine, only: kind_phys
        implicit none
        real(kind=kind_phys), parameter :: con_g     = 9.80665_kind_phys
        real(kind=kind_phys), parameter :: con_cp    = 1004.0_kind_phys
        real(kind=kind_phys), parameter :: con_rd    = 287.0_kind_phys
        real(kind=kind_phys), parameter :: con_rv    = 461.5_kind_phys
        real(kind=kind_phys), parameter :: con_hvap  = 2.5e6_kind_phys
        real(kind=kind_phys), parameter :: con_hfus  = 3.3358e5_kind_phys
        real(kind=kind_phys), parameter :: con_eps   = con_rd / con_rv
        real(kind=kind_phys), parameter :: con_epsm1 = con_eps - 1.0_kind_phys
      end module physcons

      module funcphys
        use machine, only: kind_phys
        implicit none
      contains

        real(kind=kind_phys) function fpvs(t)
          real(kind=kind_phys), intent(in) :: t
          real(kind=kind_phys) :: tc

          ! Simple saturation vapor pressure over liquid water, returned in Pa.
          ! This is intended as a compatibility fallback.  If your MPAS build
          ! already provides the GFS funcphys module, use that one instead.
          tc = t - 273.15_kind_phys
          fpvs = 611.2_kind_phys * exp(17.67_kind_phys * tc / &
                 (tc + 243.5_kind_phys))
        end function fpvs

      end module funcphys
