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
      real(kind=kind_phys), parameter :: con_hfus = 3.3358e5_kind_phys
      real(kind=kind_phys), parameter :: con_eps = con_rd/con_rv
      real(kind=kind_phys), parameter :: con_epsm1 = con_eps-1
      real(kind=kind_phys), parameter :: con_fvirt = con_rv/con_rd-1
!> minimum aerosol concentration
      real(kind=kind_phys),parameter:: qamin = 1.e-16_kind_phys
      end module physcons

      module funcphys
      use gfs_tke_edmf_machine, only: kind_phys
      use physcons, only: con_rv, con_hvap, con_hfus
      implicit none
      private
      public :: fpvs
!     GFS mixed-phase formula from module_gfs_funcphys.F fpvsx.
!     GFS physical constants; keep in sync with source physcons.
      real(kind=kind_phys), parameter :: ttp=273.16_kind_phys
      real(kind=kind_phys), parameter :: psat=610.78_kind_phys
      real(kind=kind_phys), parameter :: cvap=1846._kind_phys
      real(kind=kind_phys), parameter :: cliq=4185.5_kind_phys
      real(kind=kind_phys), parameter :: csol=2106._kind_phys
      contains
      real(kind=kind_phys) function fpvs(t)
      real(kind=kind_phys), intent(in) :: t
      real(kind=kind_phys), parameter :: tice=ttp-20._kind_phys
      real(kind=kind_phys), parameter :: dldtl=cvap-cliq
      real(kind=kind_phys), parameter :: dldti=cvap-csol
      real(kind=kind_phys), parameter :: xponal=-dldtl/con_rv
      real(kind=kind_phys), parameter :: xponbl=xponal+
     &                         con_hvap/(con_rv*ttp)
      real(kind=kind_phys), parameter :: xponai=-dldti/con_rv
      real(kind=kind_phys), parameter :: xponbi=xponai+
     &                   (con_hvap+con_hfus)/(con_rv*ttp)
      real(kind=kind_phys) :: tr,w,pvl,pvi
      tr=ttp/t
      if (t.ge.ttp) then
         fpvs=psat*(tr**xponal)*exp(xponbl*(1._kind_phys-tr))
      else if (t.lt.tice) then
         fpvs=psat*(tr**xponai)*exp(xponbi*(1._kind_phys-tr))
      else
         w=(t-tice)/(ttp-tice)
         pvl=psat*(tr**xponal)*exp(xponbl*(1._kind_phys-tr))
         pvi=psat*(tr**xponai)*exp(xponbi*(1._kind_phys-tr))
         fpvs=w*pvl+(1._kind_phys-w)*pvi
      endif
      end function fpvs
      end module funcphys
