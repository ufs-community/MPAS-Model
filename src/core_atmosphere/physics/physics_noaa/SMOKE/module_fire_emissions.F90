!>\file  module_wildfire_smoke_emissions.F90
!! This file contains the MPAS-Aerosols/RRFS wildfire emission module

module module_wildfire_smoke_emissions
!
!  This module developed by Johana Romero-Alvarez and Jordan Schnell (NOAA GSL)
!  For serious questions contact johana.romero-alvarez@noaa.gov
!
  use mpas_kind_types
  use mpas_smoke_init, only : p_smoke_fine, p_smoke_coarse
  use mpas_smoke_config

  implicit none

  private

  public :: caclulate_smoke_emissions

contains


  subroutine calculate_smoke_emissions(dt,gmt,julday,nlcat,           &
                                       xlat,xlong,xland,ivgtyp,               &
                                       ef_map, &
                                       dz8w,rho_phy,hwp,               &
                                       frp_in,fre_in,fire_type,fire_end_hr,peak_hr                                                  &
                           ids,ide, jds,jde, kds,kde,               &
                           ims,ime, jms,jme, kms,kme,               &
                           its,ite, jts,jte, kts,kte                )

   IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: julday,nlcat,                     &
                                  ids,ide, jds,jde, kds,kde,        &
                                  ims,ime, jms,jme, kms,kme,        &
                                  its,ite, jts,jte, kts,kte

   REAL(RKIND), INTENT(IN    ) :: dt,gmt

   REAL(RKIND),DIMENSION(ims:ime,jms:jme),INTENT(IN) :: xlat,xlong,xland,ivgtyp
   REAL(RKIND),DIMENSION(ims:ims,kms:kme,jms:jme), &
                                          INTENT(IN) :: dz8w,rho_phy
   REAL(RKIND),DIMENSION(ims:ime,jms:jme),INTENT(IN) :: hwp, frp_in, fre_in
   REAL(RKIND),DIMENSION(ims:ime,jms:jme),INTENT(IN) :: fire_type,fire_end_hr,peak_hr
   REAL(RKIND),DIMENSION(ims:ime,1:nlcat,jms:jme),INTENT(IN) :: ef_map


   REAL(RKIND),DIMENSION(ims:ime,kms:kme,jms:jme) :: e_bb_out_smoke_fine, &
                                                     e_bb_out_smoke_coarse
                   
  ! local
   INTEGER :: i,j,k,n
   REAL(RKIND) :: conv
   REAL(RKIND), PARAMETER :: BETA !! TODO, NAMELIST


   do j = jts,jte
   do i = its,ite



   enddo
   enddo


  end subroutine calculate_smoke_emisisons

end module module_wildfire_smoke_emissions
