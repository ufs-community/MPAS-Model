! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_driver_cloudiness
 use mpas_kind_types
 use mpas_pool_routines
 use mpas_timer, only : mpas_timer_start, mpas_timer_stop

 use mpas_atmphys_constants, only: ep_2
 use mpas_atmphys_vars
 use module_mp_thompson_cldfra3
 
 implicit none
 private
 public:: allocate_cloudiness,   &
          deallocate_cloudiness, &
          driver_cloudiness

!MPAS driver for parameterization of the diagnostic cloud fraction.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines in mpas_atmphys_driver_cloudiness:
! ----------------------------------------------
! allocate_cloudiness  : allocate local arrays for parameterization of diagnostic cloudiness.
! deallocate_cloudiness: deallocate local arrays for parameterization of diagnostic cloudiness.
! cloudiness_from_MPAS : initialize local arrays.
! cloudiness_to_MPAS   : copy local arrays to MPAS arrays.
! driver_cloudiness    : main driver (called from subroutine physics_driver).
! calc_cldincidence    : calculates the cloud fraction as 0 or 1, depending on cloud condensates.
! calc_cldfraction     : calculates the cloud fraction as a function of the relative humidity.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * added subroutine cal_cldfra3 which contains the calculation of the radiative cloud fraction in conjunction
!   with the Thompson cloud microphysics.
!   Laura D. Fowler (laura@ucar.edu) / 2016-06-04.
! * initialize the local "radiation" water vapor, cloud water, cloud ice, and snow mixing ratios. When using the
!   option "cld_fraction_thompson", the "radiative" cloud water and ice paths are updated in conjunction with
!   cloud formation, but changes to the cloud water and cloud ice mixing ratios only affect the long wave and
!   short wave radiation codes.
!   Laura D. Fowler (laura@ucar.edu) / 2016-07-05.
! * since we removed the local variable radt_cld_scheme from mpas_atmphys_vars.F, now defines radt_cld_scheme
!   as a pointer to config_radt_cld_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2017-02-16.
! * this is a bug fix. dx_p is converted from meters to kilometers prior to calling the thompson parameterization
!   of the cloud fraction.
!   Laura D. Fowler (laura@ucar.edu) / 2024-03-23.


 contains


!=================================================================================================================
 subroutine allocate_cloudiness(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in)   :: configs
!local pointer:
 character(len=StrKIND),pointer:: radt_cld_scheme
 character(len=StrKIND),pointer:: convection_scheme

 call mpas_pool_get_config(configs,'config_radt_cld_scheme',radt_cld_scheme)
 call mpas_pool_get_config(configs,'config_convection_scheme',convection_scheme)

 if(.not.allocated(dx_p)     ) allocate(dx_p(ims:ime,jms:jme)             ) 
 if(.not.allocated(xland_p)  ) allocate(xland_p(ims:ime,jms:jme)          )
 if(.not.allocated(cldfrac_p)) allocate(cldfrac_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(qvrad_p)  ) allocate(qvrad_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(qcrad_p)  ) allocate(qcrad_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(qirad_p)  ) allocate(qirad_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(qsrad_p)  ) allocate(qsrad_p(ims:ime,kms:kme,jms:jme)  )

 cld_fraction_select: select case (trim(radt_cld_scheme))

    case("cld_fraction_mynn")
       if(.not.allocated(qibl_p)) allocate(qibl_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(qcbl_p)) allocate(qcbl_p(ims:ime,kms:kme,jms:jme))
    case default

 end select cld_fraction_select

 if (convection_scheme=="cu_grell_freitas_li") then
    if(.not.allocated(qicu_p)) allocate(qicu_p(ims:ime,kms:kme,jms:jme))
    if(.not.allocated(qccu_p)) allocate(qccu_p(ims:ime,kms:kme,jms:jme))
 endif

 end subroutine allocate_cloudiness

!=================================================================================================================
 subroutine deallocate_cloudiness(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in)   :: configs
!local pointer:
 character(len=StrKIND),pointer:: radt_cld_scheme
 character(len=StrKIND),pointer:: convection_scheme

 call mpas_pool_get_config(configs,'config_radt_cld_scheme',radt_cld_scheme)
 call mpas_pool_get_config(configs,'config_convection_scheme',convection_scheme)

 if(allocated(dx_p)     ) deallocate(dx_p     )
 if(allocated(xland_p)  ) deallocate(xland_p  )
 if(allocated(cldfrac_p)) deallocate(cldfrac_p)
 if(allocated(qvrad_p)  ) deallocate(qvrad_p  )
 if(allocated(qcrad_p)  ) deallocate(qcrad_p  )
 if(allocated(qirad_p)  ) deallocate(qirad_p  )
 if(allocated(qsrad_p)  ) deallocate(qsrad_p  )

 cld_fraction_select: select case (trim(radt_cld_scheme))

    case("cld_fraction_mynn")

       if(allocated(qibl_p)) deallocate(qibl_p)
       if(allocated(qcbl_p)) deallocate(qcbl_p)
    case default

 end select cld_fraction_select

 if (convection_scheme=="cu_grell_freitas_li") then
    if(allocated(qicu_p)) deallocate(qicu_p)
    if(allocated(qccu_p)) deallocate(qccu_p)
 endif


 end subroutine deallocate_cloudiness

!=================================================================================================================
 subroutine cloudiness_from_MPAS(itimestep,configs,mesh,diag_physics,sfc_input,its,ite)
!=================================================================================================================

!input arguments:
 integer,intent(in):: its,ite
 integer,intent(in):: itimestep

!input and inout arguments:
 type(mpas_pool_type),intent(in)   :: configs
 type(mpas_pool_type),intent(in)   :: mesh
 type(mpas_pool_type),intent(in)   :: sfc_input

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local variables and pointers:
 integer:: i,j,k

 real(kind=RKIND),pointer:: len_disp
 real(kind=RKIND),dimension(:),pointer:: areaCell,meshDensity
 real(kind=RKIND),dimension(:),pointer:: xland

 real(kind=RKIND),dimension(:,:),pointer:: cldfrac_bl,qi_bl,qc_bl,cldfrac_cu,qc_cu,qi_cu

 character(len=StrKIND),pointer:: radt_cld_scheme
 character(len=StrKIND),pointer:: convection_scheme

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_len_disp',len_disp)
 call mpas_pool_get_config(configs,'config_convection_scheme',convection_scheme)

 call mpas_pool_get_array(mesh,'areaCell'   ,areaCell   )
 call mpas_pool_get_array(mesh,'meshDensity',meshDensity)

 call mpas_pool_get_array(sfc_input,'xland',xland)
 call mpas_pool_get_config(configs,'config_radt_cld_scheme',radt_cld_scheme)

 do j = jts,jte
    do i = its,ite
       dx_p(i,j)    = len_disp / meshDensity(i)**0.25
       !conversion of dx_p from meters to kilometers.
       dx_p(i,j)    = dx_p(i,j)*0.001
       xland_p(i,j) = xland(i)
    enddo

    !--- initialize the radiative cloud fraction and water vapor, cloud water, cloud ice,
    !    and snow mixing ratios:
    do k = kts,kte
    do i = its,ite
       qvrad_p(i,k,j)   = qv_p(i,k,j)
       qcrad_p(i,k,j)   = qc_p(i,k,j)
       qirad_p(i,k,j)   = qi_p(i,k,j)
       qsrad_p(i,k,j)   = qs_p(i,k,j)
       cldfrac_p(i,k,j) = 0._RKIND
    enddo
    enddo
 enddo

 cld_fraction_select: select case (trim(radt_cld_scheme))
 
    case("cld_fraction_mynn")
       call mpas_pool_get_array(diag_physics,'cldfrac_bl',cldfrac_bl)
       call mpas_pool_get_array(diag_physics,'qi_bl',qi_bl)
       call mpas_pool_get_array(diag_physics,'qc_bl',qc_bl)
       if(itimestep > 0) then
          do j = jts,jte
          do k = kts,kte
          do i = its,ite
             if (cldfrac_bl(k,i) .gt. 0.001) then
                cldfrac_p(i,k,j) = cldfrac_bl(k,i)
             endif
          enddo
          enddo
          enddo
       else
          do j = jts,jte
          do k = kts,kte
          do i = its,ite
             if	(qcrad_p(i,k,j) .gt. 1e-6 .or. &
                 qirad_p(i,k,j) .gt. 1e-9 .or. &
                 qsrad_p(i,k,j) .gt. 1e-8 ) then
                cldfrac_p(i,k,j) = 1._RKIND
      	     endif
          enddo
          enddo
          enddo
       endif

       do j = jts,jte
       do k = kts,kte
       do i = its,ite
          qcbl_p(i,k,j) = qc_bl(k,i)
          qibl_p(i,k,j) = qi_bl(k,i)
       enddo
       enddo
       enddo

    case default

 end select cld_fraction_select  
 
 if (convection_scheme=="cu_grell_freitas_li") then
    call mpas_pool_get_array(diag_physics,'cldfrac_cu',cldfrac_cu)
    call mpas_pool_get_array(diag_physics,'qi_cu',qi_cu)
    call mpas_pool_get_array(diag_physics,'qc_cu',qc_cu)
    do j = jts,jte
    do k = kts,kte
    do i = its,ite
       if (cldfrac_cu(k,i) .gt. 0.001) then
          cldfrac_p(i,k,j) = min(1._RKIND,cldfrac_p(i,k,j) + cldfrac_cu(k,i))
       endif
       qccu_p(i,k,j) = qc_cu(k,i)
       qicu_p(i,k,j) = qi_cu(k,i)
    enddo
    enddo
    enddo
 endif

 end subroutine cloudiness_from_MPAS

!=================================================================================================================
 subroutine cloudiness_to_MPAS(diag_physics,its,ite)
!=================================================================================================================

!input arguments:
 integer,intent(in):: its,ite

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local variables and pointers:
 integer:: i,j,k

 real(kind=RKIND),dimension(:,:),pointer:: cldfrac

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_array(diag_physics,'cldfrac',cldfrac)

 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    cldfrac(k,i) = cldfrac_p(i,k,j)
 enddo
 enddo
 enddo
 
 end subroutine cloudiness_to_MPAS

!=================================================================================================================
 subroutine driver_cloudiness(itimestep,configs,mesh,diag_physics,sfc_input,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in)   :: configs
 type(mpas_pool_type),intent(in)   :: mesh
 type(mpas_pool_type),intent(in)   :: sfc_input

 integer,intent(in):: its,ite
 integer,intent(in):: itimestep

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local variables and pointers:
 character(len=StrKIND),pointer:: radt_cld_scheme
 character(len=StrKIND),pointer:: convection_scheme

 integer:: i,j,k

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('')
!call mpas_log_write('--- enter subroutine driver_cloudiness:')

 call mpas_pool_get_config(configs,'config_radt_cld_scheme',radt_cld_scheme)
 call mpas_pool_get_config(configs,'config_convection_scheme',convection_scheme)

!copy MPAS arrays to local arrays:
 call cloudiness_from_MPAS(itimestep,configs,mesh,diag_physics,sfc_input,its,ite)

 cld_fraction_select: select case (trim(radt_cld_scheme))
    case("cld_incidence")
      call mpas_timer_start('calc_cldincidence')
      call calc_cldincidence(cldfrac_p,qcrad_p,qirad_p,f_qc,f_qi,its,ite)
      call mpas_timer_stop('calc_cldincidence')

    case("cld_fraction")
      call mpas_timer_start('calc_cldfraction')
      !print*, minval(cldfrac_p), minval(t_p), minval(pres_p), minval(qvrad_p), minval(qcrad_p), minval(qirad_p)
      call calc_cldfraction(cldfrac_p,t_p,pres_p,qvrad_p,qcrad_p,qirad_p,qs_p,its,ite)
      call mpas_timer_stop('calc_cldfraction')

    case("cld_fraction_thompson")
      call mpas_timer_start('cal_cldfra3')
      call cal_cldfra3( &
           cldfra = cldfrac_p , qv     = qvrad_p    , qc = qcrad_p , qi  = qirad_p , &
           qs     = qsrad_p   , p      = pres_hyd_p , t  = t_p     , rho = rho_p   , &
           xland  = xland_p   , gridkm = dx_p       ,                                &
           ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,   &
           ims = ims , ime = ime , jms = jms , jme = jme , kms = kds , kme = kme ,   &
           its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte     &
                      )
      call mpas_timer_stop('cal_cldfra3')

    case("cld_fraction_mynn")
      call mpas_timer_start('cal_mynncld')
      call cal_mynncld(cldfrac_p , qcrad_p   , qirad_p    , &
                       qcbl_p    , qibl_p    ,              &
                       its,  ite , jts , jte , kts , kte  )
      call mpas_timer_stop('cal_mynncld')

    case default

 end select cld_fraction_select

 if (convection_scheme=="cu_grell_freitas_li") then
    call mpas_timer_start('cal_gflcld')
    call cal_gflcld(cldfrac_p , qcrad_p   , qirad_p    , &
                   qccu_p    , qicu_p    ,              &
                   its,  ite , jts , jte , kts , kte  )
    call mpas_timer_stop('cal_gflcld')
 endif

!copy local arrays to MPAS grid:
 call cloudiness_to_MPAS(diag_physics,its,ite)

!call mpas_log_write('--- end subroutine driver_cloudiness.')

 end subroutine driver_cloudiness

!=================================================================================================================
 subroutine calc_cldincidence(cldfrac,qc,qi,f_qc,f_qi,its,ite)
!=================================================================================================================

!input arguments:
 logical,intent(in):: f_qc,f_qi
 integer,intent(in):: its,ite
 real(kind=RKIND),intent(in),dimension(ims:ime,kms:kme,jms:jme):: qc,qi

!output arguments:
 real(kind=RKIND),intent(out),dimension(ims:ime,kms:kme,jms:jme):: cldfrac

!local variables:
 integer:: i,j,k

 real(kind=RKIND),parameter:: thresh = 1.e-06

!-----------------------------------------------------------------------------------------------------------------

 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    cldfrac(i,k,j) = 0._RKIND
 enddo
 enddo
 enddo

 if(f_qc .and. f_qi) then
    do j = jts,jte
    do k = kts,kte
    do i = its,ite
       if(qc(i,k,j)+qi(i,k,j) .gt. thresh) cldfrac(i,k,j) = 1.0_RKIND
    enddo
    enddo
    enddo
 elseif(f_qc) then
    do j = jts,jte
    do k = kts,kte
    do i = its,ite
       if(qc(i,k,j) .gt. thresh) cldfrac(i,k,j) = 1.0_RKIND
    enddo
    enddo
    enddo
 endif
    
 end subroutine calc_cldincidence

!=================================================================================================================
 subroutine calc_cldfraction(cldfrac,t_p,pres_p,qv,qc,qi,qs,its,ite)
!=================================================================================================================

!input arguments:
 integer,intent(in):: its,ite
 real(kind=RKIND),intent(in),dimension(ims:ime,kms:kme,jms:jme):: qv,qc,qi,qs
 real(kind=RKIND),intent(in),dimension(ims:ime,kms:kme,jms:jme):: t_p,pres_p

!output arguments:
 real(kind=RKIND),intent(out),dimension(ims:ime,kms:kme,jms:jme):: cldfrac

!local variables:
 integer:: i,j,k

 real(kind=RKIND),parameter:: alpha0  = 100.
 real(kind=RKIND),parameter:: gamma   = 0.49
 real(kind=RKIND),parameter:: qcldmin = 1.e-12
 real(kind=RKIND),parameter:: pexp    = 0.25
 real(kind=RKIND),parameter:: rhgrid  = 1.0

 real(kind=RKIND),parameter:: svp1  = 0.61078
 real(kind=RKIND),parameter:: svp2  = 17.2693882
 real(kind=RKIND),parameter:: svpi2 = 21.8745584
 real(kind=RKIND),parameter:: svp3  = 35.86
 real(kind=RKIND),parameter:: svpi3 = 7.66
 real(kind=RKIND),parameter:: svpt0 = 273.15

 real(kind=RKIND):: esi,esw,qvsi,qvsw
 real(kind=RKIND):: arg,denom,qcld,qvs,rhum,subsat,weight

!-----------------------------------------------------------------------------------------------------------------

 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    cldfrac(i,k,j) = 0._RKIND
 enddo
 enddo
 enddo
 
 do j = jts,jte
 do k = kts,kte
 do i = its,ite

!... calculation of the saturation mixing ratios over water and over ice (Murray, 1966):
    esw = 1000. * svp1 * exp(svp2 * (t_p(i,k,j) - svpt0) / (t_p(i,k,j) - svp3))
    esi = 1000. * svp1 * exp(svpi2 * (t_p(i,k,j) - svpt0) / (t_p(i,k,j) - svpi3))

    qvsw = ep_2 * esw / (pres_p(i,k,j) - esw)
    qvsi = ep_2 * esi / (pres_p(i,k,j) - esi)

    qcld = qc(i,k,j) + qi(i,k,j) + qs(i,k,j)
    if(qcld .lt. qcldmin) then
       weight = 0.
    else
       weight = (qi(i,k,j) + qs(i,k,j)) / qcld
    endif

    qvs = (1-weight) * qvsw + weight * qvsi
    rhum = qv(i,k,j) / qvs

    if(qcld .lt. qcldmin) then

       !assume that the cloud fraction is equal to 0. when the cloudy mixing ratio equals 0.
       cldfrac(i,k,j) = 0.

    elseif(rhum .ge. rhgrid) then
       !assume that the cloud fraction is equal to 1. when the relative humidity equal 100%.
       cldfrac(i,k,j) = 1.

    else
       !computation of the cloud fraction:
       subsat = max(1.e-10,rhgrid*qvs-qv(i,k,j))
       denom  = subsat**gamma
       arg    = max(-6.9,-alpha0*qcld/denom) ! exp(-6.9) = 0.001

       rhum = max(1.e-10,rhum)
       cldfrac(i,k,j) = (rhum/rhgrid)**pexp*(1.-exp(arg))
       if(cldfrac(i,k,j) .lt. 0.01) cldfrac(i,k,j) = 0.

    endif

 enddo
 enddo
 enddo

 end subroutine calc_cldfraction

!=================================================================================================================
 subroutine cal_mynncld(cldfrac,qc,qi,qc_bl,qi_bl,its,ite,jts,jte,kts,kte)
!=================================================================================================================

!input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte
 real(kind=RKIND),intent(inout),dimension(ims:ime,kms:kme,jms:jme):: qc,qi
 real(kind=RKIND),intent(in),   dimension(ims:ime,kms:kme,jms:jme):: qc_bl,qi_bl
 real(kind=RKIND),intent(in),   dimension(ims:ime,kms:kme,jms:jme):: cldfrac

!local variables:
 integer:: i,j,k

 real(kind=RKIND),parameter:: qc_thresh  = 1.e-06
 real(kind=RKIND),parameter:: qi_thresh  = 1.e-09
 real(kind=RKIND),parameter:: cld_thresh = 1.e-03

 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    
    if (qc(i,k,j) < qc_thresh .AND. cldfrac(i,k,j) > cld_thresh) THEN
       !call mpas_log_write("Adding in cldfra_bl to qc at $i,$i,$i", intArgs=(/i,j,k/))
       qc(i,k,j)=qc(i,k,j) + qc_bl(i,k,j)
    endif
    if (qi(i,k,j) < qi_thresh .AND. cldfrac(i,k,j) > cld_thresh) THEN
       !call mpas_log_write("Adding in cldfra_bl to qi at $i,$i,$i", intArgs=(/i,j,k/))
       qi(i,k,j)=qi(i,k,j) + qi_bl(i,k,j)
    endif
 enddo
 enddo
 enddo

 end subroutine cal_mynncld
!=================================================================================================================
!=================================================================================================================
 subroutine cal_gflcld(cldfrac,qc,qi,qc_cu,qi_cu,its,ite,jts,jte,kts,kte)
!=================================================================================================================

!input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte
 real(kind=RKIND),intent(inout),dimension(ims:ime,kms:kme,jms:jme):: qc,qi
 real(kind=RKIND),intent(in),   dimension(ims:ime,kms:kme,jms:jme):: qc_cu,qi_cu
 real(kind=RKIND),intent(in),   dimension(ims:ime,kms:kme,jms:jme):: cldfrac

!local variables:
 integer:: i,j,k

 real(kind=RKIND),parameter:: qc_thresh  = 1.e-06
 real(kind=RKIND),parameter:: qi_thresh  = 1.e-09
 real(kind=RKIND),parameter:: cld_thresh = 1.e-03

 do j = jts,jte
 do k = kts,kte
 do i = its,ite

    if (qc(i,k,j) < qc_thresh .AND. cldfrac(i,k,j) > cld_thresh) THEN
       !call mpas_log_write("Adding in cldfra_bl to qc at $i,$i,$i", intArgs=(/i,j,k/))
       qc(i,k,j)=qc(i,k,j) + qc_cu(i,k,j)
    endif
    if (qi(i,k,j) < qi_thresh .AND. cldfrac(i,k,j) > cld_thresh) THEN
       !call mpas_log_write("Adding in cldfra_bl to qi at $i,$i,$i", intArgs=(/i,j,k/))
       qi(i,k,j)=qi(i,k,j) + qi_cu(i,k,j)
    endif
 enddo
 enddo
 enddo

 end subroutine cal_gflcld
!=================================================================================================================
 end module mpas_atmphys_driver_cloudiness
!=================================================================================================================
