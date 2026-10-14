! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the
! LICENSE file
! distributed with this code, or at
! http://mpas-dev.github.com/license.html
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODULE INTERP_MODULE
!
! This module provides routines for vertical interpolation.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module init_atm_vinterp
   use mpas_kind_types

   contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: interp_array_from_string
   !
   ! Purpose:
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real (kind=RKIND) function vertical_interp(target_z, nz, zf, order, extrap, surface_val, sealev_val, ierr)

      use mpas_log, only : mpas_log_write
      use mpas_derived_types, only : MPAS_LOG_ERR

      implicit none

      real (kind=RKIND), intent(in) :: target_z
      integer, intent(in) :: nz
      real (kind=RKIND), dimension(2,nz), intent(in) :: zf      ! zf(1,:) is column of vertical coordinate values, zf(2,:) is column of field values
      integer, intent(in), optional :: order
      integer, intent(in), optional :: extrap                   ! can take values 0 = constant, 1 = linear (default), 2 = lapse-rate
      real (kind=RKIND), intent(in), optional :: surface_val
      real (kind=RKIND), intent(in), optional :: sealev_val
      integer, intent(out), optional :: ierr
      
      integer :: k, lm, lp
      real (kind=RKIND) :: wm, wp
      real (kind=RKIND) :: slope,z_temp, T_temp, dz_temp

      integer :: interp_order, extrap_type
      real (kind=RKIND) :: surface, sealevel

      if (present(ierr)) ierr = 0
      
      if (present(order)) then
         interp_order = order
      else
         interp_order = 2
      end if

      if (present(extrap)) then
         extrap_type = extrap
      else
         extrap_type = 1
      end if

      if (present(surface_val)) then
         surface = surface_val
      else
         surface = 200100.0
      end if

      if (present(sealev_val)) then
         sealevel = sealev_val
      else
         sealevel = 201300.0
      end if

      !
      ! Extrapolation required
      !
      if (target_z < zf(1,1)) then
         if (extrap_type == 0) then
            vertical_interp = zf(2,1)
         else if (extrap_type == 1) then
            slope = (zf(2,2) - zf(2,1)) / (zf(1,2) - zf(1,1))
            vertical_interp = zf(2,1) + slope * (target_z - zf(1,1))
         else if (extrap_type == 2) then
            vertical_interp = zf(2,1) - (target_z - zf(1,1))*0.0065
         else if (extrap_type == 3) then
            vertical_interp = zf(2,1) - (target_z - zf(1,1))*0.0065
         end if
         return
      end if
      if (target_z >= zf(1,nz)) then
         if (extrap_type == 0) then
            vertical_interp = zf(2,nz)
         else if (extrap_type == 1) then
            slope = (zf(2,nz) - zf(2,nz-1)) / (zf(1,nz) - zf(1,nz-1))
            vertical_interp = zf(2,nz) + slope * (target_z - zf(1,nz))
         else if (extrap_type == 2) then
             z_temp = zf(1,nz)
             T_temp = zf(2,nz)

             ! Layer 1: Troposphere (up to 11,000 m)
             ! Gamma = -6.5 K/km (Cooling)
             if (z_temp < 11000.0 .and. target_z > z_temp) then
                dz_temp = min(target_z, 11000.0) - z_temp
                T_temp = T_temp - 0.0065 * dz_temp
                z_temp = z_temp + dz_temp
             end if

             ! Layer 2: Lower Stratosphere / Tropopause (11km to 20km)
             ! Gamma = 0.0 K/km (Isothermal)
             if (z_temp >= 11000.0 .and. z_temp < 20000.0 .and. target_z > z_temp) then
                dz_temp = min(target_z, 20000.0) - z_temp
                T_temp = T_temp - 0.0 * dz_temp
                z_temp = z_temp + dz_temp
             end if

             ! Layer 3: Middle Stratosphere (20km to 32km)
             ! Gamma = +1.0 K/km (Warming)
             if (z_temp >= 20000.0 .and. z_temp < 32000.0 .and. target_z > z_temp) then
                dz_temp = min(target_z, 32000.0) - z_temp
                T_temp = T_temp - (-0.001) * dz_temp
                z_temp = z_temp + dz_temp
             end if

             ! Layer 4: Upper Stratosphere (32km to 47km)
             ! Gamma = +2.8 K/km (Warming)
             if (z_temp >= 32000.0 .and. z_temp < 47000.0 .and. target_z > z_temp) then
                dz_temp = min(target_z, 47000.0) - z_temp
                T_temp = T_temp - (-0.0028) * dz_temp
                z_temp = z_temp + dz_temp
             end if

             ! Layer 5: Stratopause (47km to 51km)
             ! Gamma = 0.0 K/km (Isothermal)
             if (z_temp >= 47000.0 .and. z_temp < 51000.0 .and. target_z > z_temp) then
                dz_temp = min(target_z, 51000.0) - z_temp
                T_temp = T_temp - 0.0 * dz_temp
                z_temp = z_temp + dz_temp
             end if

             ! Layer 6: Lower Mesosphere (51km to 71km)
             ! Gamma = -2.8 K/km (Cooling)
             if (z_temp >= 51000.0 .and. z_temp < 71000.0 .and. target_z > z_temp) then
                dz_temp = min(target_z, 71000.0) - z_temp
                T_temp = T_temp - 0.0028 * dz_temp
                z_temp = z_temp + dz_temp
             end if

             ! Layer 7: Upper Mesosphere (71km to 85km)
             ! Gamma = -2.0 K/km (Cooling)
             if (z_temp >= 71000.0 .and. z_temp < 85000.0 .and. target_z > z_temp) then
                dz_temp = min(target_z, 85000.0) - z_temp
                T_temp = T_temp - 0.002 * dz_temp
                z_temp = z_temp + dz_temp
             end if

             if (z_temp >= 85000.0 .and. target_z > z_temp) then
                call mpas_log_write('target_z >= 85km. Stop!', messageType=MPAS_LOG_ERR)
                if (present(ierr)) ierr = 1
                return
             end if

             ! Layer 8: Lower Thermosphere (85km to 110km)
             ! Gamma = +2.0 K/km (Warming)
             !if (z_temp >= 85000.0 .and. z_temp < 110000.0 .and. target_z > z_temp) then
             !   dz_temp = min(target_z, 110000.0) - z_temp
             !   T_temp = T_temp - (-0.002) * dz_temp
             !   z_temp = z_temp + dz_temp
             !end if

         else if (extrap_type == 3) then
             vertical_interp = zf(2,nz)
         end if
         return
      end if


      !
      ! No extrapolation required
      !
      do k=1,nz-1
         if (target_z >= zf(1,k) .and. target_z < zf(1,k+1)) then
            lm = k
            lp = k+1
            wm = (zf(1,k+1) - target_z) / (zf(1,k+1) - zf(1,k))
            wp = (target_z - zf(1,k)) / (zf(1,k+1) - zf(1,k))
            exit
         end if
      end do

      vertical_interp = wm*zf(2,lm) + wp*zf(2,lp)

      return

   end function vertical_interp

end module init_atm_vinterp
