! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_driver
 use mpas_kind_types
 use mpas_pool_routines

 use mpas_atmphys_driver_cloudiness
 use mpas_atmphys_driver_convection
 use mpas_atmphys_driver_gwdo
 use mpas_atmphys_driver_lsm
 use mpas_atmphys_driver_lsm_noahmp
 use mpas_atmphys_driver_pbl
 use mpas_atmphys_driver_radiation_lw
 use mpas_atmphys_driver_radiation_sw 
 use mpas_atmphys_driver_seaice,only: allocate_seaice,deallocate_seaice,driver_seaice
 use mpas_atmphys_driver_sfclayer
 use mpas_atmphys_driver_oml
 use mpas_atmphys_driver_smoke
 use mpas_atmphys_constants
 use mpas_atmphys_interface
 use mpas_gfs_base_pressure_option, only: gfs_fixed_base_mode
 use mpas_atmphys_sfc_diagnostics,only: atmphys_sfc_diagnostics
 use mpas_atmphys_update
 use mpas_atmphys_vars, only: l_camlw,l_conv,l_radtlw,l_radtsw,do_chemistry, &
      pres_base_mid_gfs_p, pres_base_int_gfs_p, &
      ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte
 use mpas_timer

 implicit none
 private
 public:: physics_driver

   ! Provides definition of halo_exchange_routine
#include "mpas_halo_interface.inc"

!MPAS top physics driver.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutine physics_driver is the top physics driver from which separate drivers for all physics
! parameterizations, except cloud microphysics parameterizations are called.
!
! subroutines called in mpas_atmphys_driver:
! ------------------------------------------
! allocate_forall_physics     : allocate local arrays defining atmospheric soundings (pressure,..)
! allocate_cloudiness         : allocate all local arrays used in driver_cloudiness.
! allocate_convection         : allocate all local arrays used in driver_convection.
! allocate_gwdo               : allocate all local arrays used in driver_gwdo.
! allocate_lsm                : allocate all local arrays used in driver_lsm.
! allocate_pbl                : allocate all local arrays used in driver_pbl.
! allocate_radiation_lw       : allocate all local arrays used in driver_radiation_lw.
! allocate_radiation_sw       : allocate all local arrays used in driver_radiation_sw.
! allocate_sfclayer           : allocate all local arrays used in driver_sfclayer.
!
! deallocate_forall_physics   : deallocate local arrays defining atmospheric soundings.
! deallocate_cloudiness       : dedeallocate all local arrays used in driver_cloudiness.
! deallocate_convection       : deallocate all local arrays used in driver_convection.
! deallocate_gwdo             : deallocate all local arrays used in driver_gwdo.
! deallocate_lsm              : deallocate all local arrays used in driver_lsm.
! deallocate_pbl              : deallocate all local arrays used in driver_pbl.
! deallocate_radiation_lw     : deallocate all local arrays used in driver_radiation_lw.
! deallocate_radiation_sw     : deallocate all local arrays used in driver_radiation_sw.
! deallocate_sfclayer         : deallocate all local arrays used in driver_sfclayer.
!
! MPAS_to_physics             :
! driver_cloudiness           : driver for parameterization of fractional cloudiness.
! driver_convection           : driver for parameterization of convection.
! driver_gwdo                 : driver for parameterization of gravity wave drag over orography.
! driver_lsm                  : driver for land-surface scheme.
! driver_pbl                  : driver for planetary boundary layer scheme.
! driver_radiation_sw         : driver for short wave radiation schemes.
! driver_radiation_lw         : driver for long wave radiation schemes.
! driver_sfclayer             : driver for surface layer scheme.
! update_convection_step1     : updates lifetime of deep convective clouds in Kain-Fritsch scheme.
! update_convection_step2     : updates accumulated precipitation output from convection schemes.
! update_radiation_diagnostics: updates accumualted radiation diagnostics from radiation schemes.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * removed call to calculate atmospheric soundings for the hydrostatic dynamical core.
!   Laura D. Fowler (2013-05-06).
! * removed the namelist option config_eddy_scheme and associated sourcecode.
! * removed the namelist option config_conv_shallow_scheme and associated sourcecode.
!   Laura D. Fowler (birch.ucar.edu) / 2013-05-29.
! * added block%atm_input in calls to subroutines driver_radiation_lw amd driver_radiation_lw.
!   Laura D. Fowler (laura@ucar.edu) / 2013-07-03.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * renamed config_conv_deep_scheme to config_convection_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2014-09-18.
! * in the call to driver_convection, added block%configs needed for the implementation of the
!   Grell-Freitas convection scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-30.
! * modified the call to the subroutines driver_sfclayer and driver_pbl for the implementation
!   of the MYNN surface layer scheme and PBL schemes. itimestep and block%configs are added to
!   the argument list.
!   Laura D. Fowler (laura@ucar.edu) / 2015-01-06.
! * now only call subroutine update_convection_step2 when config_convection_scheme is not off.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-13.
! * modified call to driver_cloudiness to accomodate the calculation of the cloud fraction with the Thompson
!   cloud microphysics scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2016-06-04.
! * added call to the Noah-MP land surface scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2024-03-11.


 contains


!=================================================================================================================
 subroutine physics_driver(domain,itimestep,xtime_s,exchange_halo_group)
!=================================================================================================================

!input arguments:
 integer,intent(in):: itimestep
 real(kind=RKIND),intent(in):: xtime_s
 procedure (halo_exchange_routine) :: exchange_halo_group

!inout arguments:
 type(domain_type),intent(inout):: domain

!local pointers:
 type(mpas_pool_type),pointer::  configs,             &
                                 mesh,                &
                                 state,               &
                                 diag,                &
                                 diag_physics,        &
                                 diag_physics_noahmp, &
                                 output_noahmp,       &
                                 tend_physics,        &
                                 atm_input,           &
                                 ngw_input,           &
                                 sfc_input,           &
                                 emission_input

 logical,pointer:: config_frac_seaice

 character(len=StrKIND),pointer:: config_bucket_update,     &
                                  config_convection_scheme, &
                                  config_gwdo_scheme,       &
                                  config_lsm_scheme,        &
                                  config_pbl_scheme,        &
                                  config_radt_lw_scheme,    &
                                  config_radt_sw_scheme,    &
                                  config_sfclayer_scheme

 logical, pointer:: config_oml1d
 real(kind=RKIND),pointer:: config_bucket_radt

!local variables:
 type(block_type),pointer:: block

 integer:: time_lev
 integer:: thread

 integer,pointer:: nThreads
 integer,dimension(:),pointer:: cellSolveThreadStart, cellSolveThreadEnd

!=================================================================================================================
!call mpas_log_write('')
!call mpas_log_write('--- enter subroutine mpas_atmphys_driver:')

 call mpas_timer_start('physics driver')

 call mpas_pool_get_config(domain%configs,'config_convection_scheme',config_convection_scheme)
 call mpas_pool_get_config(domain%configs,'config_gwdo_scheme'      ,config_gwdo_scheme      )
 call mpas_pool_get_config(domain%configs,'config_lsm_scheme'       ,config_lsm_scheme       )
 call mpas_pool_get_config(domain%configs,'config_pbl_scheme'       ,config_pbl_scheme       )
 call mpas_pool_get_config(domain%configs,'config_radt_lw_scheme'   ,config_radt_lw_scheme   )
 call mpas_pool_get_config(domain%configs,'config_radt_sw_scheme'   ,config_radt_sw_scheme   )
 call mpas_pool_get_config(domain%configs,'config_sfclayer_scheme'  ,config_sfclayer_scheme  )
 call mpas_pool_get_config(domain%configs,'config_bucket_radt'      ,config_bucket_radt      )
 call mpas_pool_get_config(domain%configs,'config_bucket_update'    ,config_bucket_update    )
 call mpas_pool_get_config(domain%configs,'config_frac_seaice'      ,config_frac_seaice      ) 
 call mpas_pool_get_config(domain%configs,'config_oml1d'            ,config_oml1d            )

 if(config_convection_scheme .ne. 'off' .or. &
    config_lsm_scheme        .ne. 'off' .or. &
    config_pbl_scheme        .ne. 'off' .or. & 
    config_radt_lw_scheme    .ne. 'off' .or. &
    config_radt_sw_scheme    .ne. 'off' .or. &
    config_sfclayer_scheme   .ne. 'off') then

 block => domain % blocklist
 do while(associated(block))

    call mpas_pool_get_subpool(block%structs,'mesh'               ,mesh               )
    call mpas_pool_get_subpool(block%structs,'state'              ,state              )
    call mpas_pool_get_subpool(block%structs,'diag'               ,diag               )
    call mpas_pool_get_subpool(block%structs,'diag_physics'       ,diag_physics       )
    call mpas_pool_get_subpool(block%structs,'diag_physics_noahmp',diag_physics_noahmp)
    call mpas_pool_get_subpool(block%structs,'output_noahmp'      ,output_noahmp      )
    call mpas_pool_get_subpool(block%structs,'atm_input'          ,atm_input          )
    call mpas_pool_get_subpool(block%structs,'sfc_input'          ,sfc_input          )
    call mpas_pool_get_subpool(block%structs,'ngw_input'          ,ngw_input          )
    call mpas_pool_get_subpool(block%structs,'tend_physics'       ,tend_physics       )
    if ( do_chemistry ) then
       call mpas_pool_get_subpool(block%structs,'emission_input',emission_input)
    endif

    call mpas_pool_get_dimension(block%dimensions,'nThreads',nThreads)

    call mpas_pool_get_dimension(block%dimensions,'cellSolveThreadStart',cellSolveThreadStart)
    call mpas_pool_get_dimension(block%dimensions,'cellSolveThreadEnd',cellSolveThreadEnd)

    !allocate arrays shared by all physics parameterizations:
    call allocate_forall_physics(block%configs)

    !physics prep step:
    time_lev = 1

!$OMP PARALLEL DO
    do thread=1,nThreads
       call MPAS_to_physics(block%configs,mesh,state,time_lev,diag,diag_physics, &
                            cellSolveThreadStart(thread),cellSolveThreadEnd(thread))
    end do
!$OMP END PARALLEL DO

    ! GFS-specific base-pressure preparation (mode 1 only).
    ! Keep all reconstruction outside the common MPAS_to_physics interface.
    ! Modes 0 and 2, and non-GFS physics, do not use these arrays.
    if ((trim(config_pbl_scheme) == 'bl_tke_edmf' .or. &
         trim(config_convection_scheme) == 'cu_gfs_sas') .and. &
         gfs_fixed_base_mode()) then

       if (.not.allocated(pres_base_mid_gfs_p)) then
          allocate(pres_base_mid_gfs_p(ims:ime,kms:kme,jms:jme))
       endif
       if (.not.allocated(pres_base_int_gfs_p)) then
          allocate(pres_base_int_gfs_p(ims:ime,kms:kme,jms:jme))
       endif

!$OMP PARALLEL DO
       do thread=1,nThreads
          call prepare_gfs_base_pressure(mesh,diag, &
               cellSolveThreadStart(thread),cellSolveThreadEnd(thread))
       enddo
!$OMP END PARALLEL DO
    endif

    !call to cloud scheme:
    if(l_radtlw .or. l_radtsw) then
       call allocate_cloudiness(block%configs)
!$OMP PARALLEL DO
       do thread=1,nThreads
          call driver_cloudiness(itimestep,block%configs,mesh,diag_physics,sfc_input, &
                                 cellSolveThreadStart(thread),cellSolveThreadEnd(thread))
       end do
!$OMP END PARALLEL DO
    endif

    !call to short wave radiation scheme:
    if(l_radtsw) then
       time_lev = 1
       call allocate_radiation_sw(block%configs,xtime_s)
!$OMP PARALLEL DO
       do thread=1,nThreads
          call driver_radiation_sw(itimestep,block%configs,mesh,state,time_lev,diag_physics, &
                                   atm_input,sfc_input,tend_physics,xtime_s, &
                                   cellSolveThreadStart(thread),cellSolveThreadEnd(thread))
       end do
!$OMP END PARALLEL DO
    endif

    !call to long wave radiation scheme:
    if(l_radtlw) then
       time_lev = 1
       call allocate_radiation_lw(block%configs,xtime_s)
!$OMP PARALLEL DO
       do thread=1,nThreads
          call driver_radiation_lw(xtime_s,block%configs,mesh,state,time_lev,diag_physics, &
                                   atm_input,sfc_input,tend_physics, &
                                   cellSolveThreadStart(thread),cellSolveThreadEnd(thread))
       end do
!$OMP END PARALLEL DO
    endif

    !call to accumulate long- and short-wave diagnostics if needed:
     if(config_bucket_update /= 'none' .and. config_bucket_radt .gt. 0._RKIND) then
!$OMP PARALLEL DO
       do thread=1,nThreads
          call update_radiation_diagnostics(block%configs,mesh,diag_physics, &
                                            cellSolveThreadStart(thread),cellSolveThreadEnd(thread))
       end do
!$OMP END PARALLEL DO
     endif

    !deallocate all radiation arrays:
    if(config_radt_sw_scheme.ne.'off' .or. config_radt_lw_scheme.ne.'off') &
       call deallocate_cloudiness(block%configs)
    if(config_radt_sw_scheme.ne.'off') call deallocate_radiation_sw(block%configs)
    if(config_radt_lw_scheme.ne.'off') call deallocate_radiation_lw(block%configs)

    !call to surface-layer scheme:
    if(config_sfclayer_scheme .ne. 'off') then
       call allocate_sfclayer(block%configs)
!$OMP PARALLEL DO
       do thread=1,nThreads
          call driver_sfclayer(itimestep,block%configs,mesh,diag_physics,sfc_input, &
                               cellSolveThreadStart(thread),cellSolveThreadEnd(thread))
       end do
!$OMP END PARALLEL DO
       call deallocate_sfclayer(block%configs)
    endif

    !call to 1d ocean mixed-layer model
    if(config_oml1d) call driver_oml1d(block%configs,mesh,diag,diag_physics,sfc_input)

    !call to land-surface scheme:
    if(config_lsm_scheme .ne. 'off') then
       if(config_lsm_scheme .eq. 'sf_noah' .or. config_lsm_scheme .eq. 'sf_ruc') then
          call allocate_lsm(block%configs)
!$OMP PARALLEL DO
          do thread=1,nThreads
             call driver_lsm(itimestep,block%configs,mesh,diag_physics,sfc_input, &
                             cellSolveThreadStart(thread),cellSolveThreadEnd(thread))
          end do
!$OMP END PARALLEL DO
       call deallocate_lsm(block%configs)

       elseif(config_lsm_scheme == 'sf_noahmp') then
          do thread=1,nThreads
             call driver_lsm_noahmp(block%configs,mesh,state,time_lev,diag,diag_physics,   &
                                    diag_physics_noahmp,output_noahmp,sfc_input,itimestep, &
                                    cellSolveThreadStart(thread),cellSolveThreadEnd(thread))
          enddo
       endif

       call allocate_seaice(block%configs)
!$OMP PARALLEL DO
       do thread=1,nThreads
          call driver_seaice(itimestep,block%configs,mesh,diag_physics,sfc_input, &
                             cellSolveThreadStart(thread),cellSolveThreadEnd(thread))
       enddo
!$OMP END PARALLEL DO
       call deallocate_seaice(block%configs)

!$OMP PARALLEL DO
       do thread=1,nThreads
          call atmphys_sfc_diagnostics(block%configs,mesh,diag,diag_physics,sfc_input,output_noahmp, &
                           cellSolveThreadStart(thread),cellSolveThreadEnd(thread))
       enddo
!$OMP END PARALLEL DO
    endif

    !call to mpas smoke model:
    if ( do_chemistry ) then
       time_lev = 1
       call allocate_smoke(block%configs)
!$OMP PARALLEL DO
       do thread=1,nThreads
          call driver_smoke(itimestep,time_lev,emission_input,state,block%configs,mesh,      &
                            sfc_input,diag_physics,tend_physics,                         &
                            cellSolveThreadStart(thread),cellSolveThreadEnd(thread)      )
       end do
!$OMP END PARALLEL DO
       call deallocate_smoke(block%configs)
    endif

    !call to pbl schemes:
    if(config_pbl_scheme .ne. 'off' .and. config_sfclayer_scheme .ne. 'off') then
       call allocate_pbl(block%configs)
!$OMP PARALLEL DO
       do thread=1,nThreads
          call driver_pbl(itimestep,block%configs,state,mesh,sfc_input,diag_physics,tend_physics, &
                          cellSolveThreadStart(thread),cellSolveThreadEnd(thread))
       end do
!$OMP END PARALLEL DO
       call deallocate_pbl(block%configs)

    endif

    !call to gravity wave drag over orography scheme:
    if(config_gwdo_scheme .ne. 'off') then
       call allocate_gwdo(block%configs)
!$OMP PARALLEL DO
       do thread=1,nThreads
          call driver_gwdo(itimestep,block%configs,mesh,sfc_input,ngw_input,diag_physics,       &
                           tend_physics,cellSolveThreadStart(thread),cellSolveThreadEnd(thread))
       end do
!$OMP END PARALLEL DO
       call deallocate_gwdo(block%configs)
    endif

    !call to convection scheme:
!$OMP PARALLEL DO
    do thread=1,nThreads
       call update_convection_step1(block%configs,diag_physics,tend_physics, &
                            cellSolveThreadStart(thread),cellSolveThreadEnd(thread))
    end do
!$OMP END PARALLEL DO
    if(l_conv) then
       call allocate_convection(block%configs)
!$OMP PARALLEL DO
       do thread=1,nThreads
          call driver_convection(itimestep,block%configs,state,mesh,sfc_input,diag_physics,tend_physics,diag,domain%dminfo, &
                                 cellSolveThreadStart(thread),cellSolveThreadEnd(thread),exchange_halo_group,block)
       end do
!$OMP END PARALLEL DO
       call deallocate_convection(block%configs)
    endif
    !update diagnostics:
    if(config_convection_scheme .ne. 'off') then
!$OMP PARALLEL DO
       do thread=1,nThreads
          call update_convection_step2(block%configs,diag_physics, &
                               cellSolveThreadStart(thread),cellSolveThreadEnd(thread))
       end do
!$OMP END PARALLEL DO
    end if

    ! Release mode-1 GFS pressure arrays after all physics drivers finish.
    ! No mode check is needed: these tests are safe for modes 0 and 2.
    if (allocated(pres_base_mid_gfs_p)) deallocate(pres_base_mid_gfs_p)
    if (allocated(pres_base_int_gfs_p)) deallocate(pres_base_int_gfs_p)

    !deallocate arrays shared by all physics parameterizations:
    call deallocate_forall_physics(block%configs)

    block => block % next
 end do 

 endif

 call mpas_timer_stop('physics driver')

!call mpas_log_write('--- enter subroutine mpas_atmphys_driver:')
!call mpas_log_write('')

 end subroutine physics_driver

!=================================================================================================================
!=================================================================================================================
!=================================================================================================================
! GFS-only base-pressure preparation. Called by physics driver, not MPAS_to_physics.
! This routine is tile-safe when each OpenMP thread owns disjoint [its:ite].
 subroutine prepare_gfs_base_pressure(mesh,diag,its,ite)
   type(mpas_pool_type), intent(in) :: mesh,diag
   integer, intent(in) :: its,ite
   real(kind=RKIND), dimension(:,:), pointer :: pressure_b,zgrid
   integer :: i,j

   call mpas_pool_get_array(diag,'pressure_base',pressure_b)
   call mpas_pool_get_array(mesh,'zgrid',zgrid)

   if (kte+1 > kme) error stop 'GFS_BASE_PRESSURE: kme must include kte+1'
   if (.not.allocated(pres_base_mid_gfs_p) .or. &
       .not.allocated(pres_base_int_gfs_p)) &
       error stop 'GFS_BASE_PRESSURE: reference arrays not allocated'
   do j=jts,jte
      do i=its,ite
         call gfs_base_pressure_column( &
              pressure_b(kts:kte,i), zgrid(kts:kte+1,i), &
              pres_base_mid_gfs_p(i,kts:kte,j), &
              pres_base_int_gfs_p(i,kts:kte+1,j))
      enddo
   enddo
 end subroutine prepare_gfs_base_pressure

! Construct strictly monotone reference interfaces using geometric-height
! interpolation/extrapolation of ln(pressure_base). MPAS vertical index increases
! upward. The routine rejects invalid reference states rather than repairing them.
 subroutine gfs_base_pressure_column(pc, zi, pm, pi)
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 real(kind=RKIND), intent(in) :: pc(:),zi(:)
 real(kind=RKIND), intent(out) :: pm(:),pi(:)
 real(kind=RKIND) :: zc(size(pc)), slope(size(pc)-1), sl, zr
 integer :: n,k
 n=size(pc)
 if(n<2 .or. size(zi)/=n+1) error stop 'GFS_BASE_PRESSURE: bad column dimensions'
 do k=1,n
   if (.not.ieee_is_finite(pc(k)) .or. pc(k)<=0._RKIND) error stop 'GFS_BASE_PRESSURE: invalid base pressure'
   if (.not.ieee_is_finite(zi(k)) .or. .not.ieee_is_finite(zi(k+1))) &
      error stop 'GFS_BASE_PRESSURE: nonfinite zgrid'
   if (zi(k+1)<=zi(k)) error stop 'GFS_BASE_PRESSURE: nonpositive dz'
   zc(k)=0.5_RKIND*(zi(k)+zi(k+1))
 enddo
 do k=1,n-1
   if (pc(k+1)>=pc(k)) error stop 'GFS_BASE_PRESSURE: nonmonotone pressure_base'
   slope(k)=(log(pc(k+1))-log(pc(k)))/(zc(k+1)-zc(k))
 enddo
 pi(1)=exp(log(pc(1))+slope(1)*(zi(1)-zc(1)))
 do k=2,n
   sl=(log(pc(k))-log(pc(k-1)))/(zc(k)-zc(k-1))
   pi(k)=exp(log(pc(k-1))+sl*(zi(k)-zc(k-1)))
 enddo
 pi(n+1)=exp(log(pc(n))+slope(n-1)*(zi(n+1)-zc(n)))
 do k=1,n
   if (pi(k)<=pi(k+1)) error stop 'GFS_BASE_PRESSURE: nonpositive delp'
   pm(k)=pc(k)
 enddo
 end subroutine gfs_base_pressure_column

 end module mpas_atmphys_driver
!=================================================================================================================
