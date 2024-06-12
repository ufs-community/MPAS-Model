! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!=================================================================================================================
 module mpas_atmphys_init_microphysics
 use mpas_dmpar
 use mpas_kind_types
 use mpas_pool_routines

 use mpas_atmphys_vars
 use mpas_atmphys_constants
 use mpas_atmphys_utilities
 use module_mp_thompson_params, only: Nt_c_l, Nt_c_o, aero_max, nwfa_default
 use module_mp_thompson_utils, only: make_hydrometeor_number_concentrations
!use module_mp_thompson, only: is_aerosol_aware,naCCN0,naCCN1,naIN0,naIN1,ntb_arc,ntb_arw,ntb_art,ntb_arr, &
!                              ntb_ark,tnccn_act

 implicit none
 private
 public:: init_thompson_clouddroplets_forMPAS

!MPAS main initialization of the Thompson parameterization of cloud microphysics with nucleation of cloud
!droplets based on distributions of CCNs and INs (aerosol-aware parameterization).
!Laura D. Fowler (send comments to laura@ucar.edu).
!2016-03-28.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * added "use mpas_dmpar" at the top of the module.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-04.


 contains


!=================================================================================================================
 subroutine init_thompson_clouddroplets_forMPAS(mesh,sfc_input,diag,diag_physics,state,configs)
!=================================================================================================================

!input variables:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: sfc_input

!inout variables:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: state
 type(mpas_pool_type),intent(in):: diag
 type(mpas_pool_type),intent(in):: configs

!local variables and pointers:
 integer,pointer:: index_qc, index_nc, index_qr, index_nr
 integer,pointer:: index_qi, index_ni, index_qv, index_nwfa
 integer,pointer:: nCellsSolve
 integer,pointer:: nVertLevels
 integer,dimension(:),pointer:: landmask

 logical,pointer:: thompson_aerosol_aware

 real(kind=RKIND),dimension(:),pointer:: nt_c,mu_c,nwfa_emissions
 real(kind=RKIND),dimension(:,:),pointer:: pressure_b, pressure_p, theta_m, exner, zgrid
 real(kind=RKIND),dimension(:,:),pointer  :: qv,qc,qr,qi,nc,nr,ni,nwfa
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

 ! local variables:
 integer:: i,k,j,iCell
 logical:: do_init_number_concentrations
 integer,parameter:: k_sfc = 1
 real(kind=RKIND):: deltaz

 do_init_number_concentrations = .false.
!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('')
!call mpas_log_write('--- enter subroutine init_thompson_clouddroplets_forMPAS:')

 call mpas_pool_get_config(configs,'config_thompson_aerosol_aware',thompson_aerosol_aware)

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)

 call mpas_pool_get_array(mesh,'zgrid',zgrid)
 call mpas_pool_get_array(state,'scalars',scalars)
 call mpas_pool_get_dimension(state,'index_nwfa',index_nwfa)
 call mpas_pool_get_array(diag_physics,'nwfa_emissions',nwfa_emissions)
 
 do iCell = 1, nCellsSolve
    deltaz = zgrid(k_sfc+1,iCell)-zgrid(k_sfc,iCell)
    if (deltaz < 9.0) then
       deltaz = 9.0
    endif
    nwfa_emissions(iCell) = max(nwfa_default, min(aero_max, scalars(index_nwfa,k_sfc,iCell))) * &
         0.000196 * 50.0 / deltaz
 enddo
 call physics_message('Setting surface wfa emissions for tempo microphysics')

 if ((thompson_aerosol_aware) .and. (do_init_number_concentrations)) then
    call mpas_pool_get_array(state,'scalars',scalars)
    call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

    call mpas_pool_get_array(diag,'exner',exner)
    call mpas_pool_get_array(diag,'pressure_p',pressure_p)
    call mpas_pool_get_array(diag,'pressure_base',pressure_b)
    call mpas_pool_get_array(state,'theta_m',theta_m,1)

    call mpas_pool_get_dimension(state,'index_qc',index_qc)
    call mpas_pool_get_dimension(state,'index_nc',index_nc)
    call mpas_pool_get_dimension(state,'index_qr',index_qr)
    call mpas_pool_get_dimension(state,'index_nr',index_nr)
    call mpas_pool_get_dimension(state,'index_qi',index_qi)
    call mpas_pool_get_dimension(state,'index_ni',index_ni)
    call mpas_pool_get_dimension(state,'index_qv',index_qv)
    call mpas_pool_get_dimension(state,'index_nwfa',index_nwfa)

    if(.not.allocated(pres_p) ) allocate(pres_p(ims:ime,kms:kme,jms:jme) )
    if(.not.allocated(qv_p) ) allocate(qv_p(ims:ime,kms:kme,jms:jme) )
    if(.not.allocated(th_p) ) allocate(th_p(ims:ime,kms:kme,jms:jme) )
    if(.not.allocated(t_p) ) allocate(t_p(ims:ime,kms:kme,jms:jme) )
    if(.not.allocated(rho_p) ) allocate(rho_p(ims:ime,kms:kme,jms:jme) )

    if(.not.allocated(qc_p) ) allocate(qc_p(ims:ime,kms:kme,jms:jme) )
    if(.not.allocated(nc_p) ) allocate(nc_p(ims:ime,kms:kme,jms:jme) )
    if(.not.allocated(qr_p) ) allocate(qr_p(ims:ime,kms:kme,jms:jme) )
    if(.not.allocated(nr_p) ) allocate(nr_p(ims:ime,kms:kme,jms:jme) )
    if(.not.allocated(qi_p) ) allocate(qi_p(ims:ime,kms:kme,jms:jme) )
    if(.not.allocated(ni_p) ) allocate(ni_p(ims:ime,kms:kme,jms:jme) )
    if(.not.allocated(qv_p) ) allocate(qv_p(ims:ime,kms:kme,jms:jme) )
    if(.not.allocated(nwfa_p) ) allocate(nwfa_p(ims:ime,kms:kme,jms:jme) )

    qc => scalars(index_qc,:,:)
    nc => scalars(index_nc,:,:)
    qr => scalars(index_qr,:,:)
    nr => scalars(index_nr,:,:)
    qi => scalars(index_qi,:,:)
    ni => scalars(index_ni,:,:)
    qv => scalars(index_qv,:,:)
    nwfa => scalars(index_nwfa,:,:)

    do j = jts, jte
       do k = kts, kte
          do i = its, ite
             pres_p(i,k,j) = pressure_p(k,i) + pressure_b(k,i)
             qv_p(i,k,j) = qv(k,i)
             th_p(i,k,j) = theta_m(k,i) / (1._RKIND + R_v/R_d * qv_p(i,k,j))
             t_p(i,k,j) = th_p(i,k,j)*exner(k,i)
             rho_p(i,k,i) = ep_2*pres_p(i,k,j) / (R_d*t_p(i,k,j)*(qv_p(i,k,j)+ep_2))

             qc_p(i,k,j) = qc(k,i)
             nc_p(i,k,j) = nc(k,i)
             qr_p(i,k,j) = qr(k,i)
             nr_p(i,k,j) = nr(k,i)
             qi_p(i,k,j) = qi(k,i)
             ni_p(i,k,j) = ni(k,i)
             qv_p(i,k,j) = qv(k,i)
             nwfa_p(i,k,j) = nwfa(k,i)

          enddo
       enddo
    enddo

    call make_hydrometeor_number_concentrations(qc_p, qr_p, qi_p, nwfa_p, t_p, rho_p, nc_p, nr_p, ni_p)

     do j = jts, jte
       do k = kts, kte
          do i = its, ite
             nc(k,i) = nc_p(i,k,j)
             nr(k,i) = nr_p(i,k,j)
             ni(k,i) = ni_p(i,k,j)
          enddo
       enddo
    enddo

    if(allocated(pres_p)  ) deallocate(pres_p  )
    if(allocated(qv_p)  ) deallocate(qv_p  )
    if(allocated(th_p)  ) deallocate(th_p  )
    if(allocated(t_p)  ) deallocate(t_p  )
    if(allocated(rho_p)  ) deallocate(rho_p  )

    if(allocated(qc_p)  ) deallocate(qc_p  )
    if(allocated(nc_p)  ) deallocate(nc_p  )
    if(allocated(qr_p)  ) deallocate(qr_p  )
    if(allocated(nr_p)  ) deallocate(nr_p  )
    if(allocated(qi_p)  ) deallocate(qi_p  )
    if(allocated(nr_p)  ) deallocate(ni_p  )
    if(allocated(qv_p)  ) deallocate(qv_p  )
    if(allocated(nwfa_p)  ) deallocate(nwfa_p  )

 else ! non aerosol-aware initial number concentrations

    call mpas_pool_get_array(sfc_input,'landmask',landmask)
    call mpas_pool_get_array(diag_physics,'nt_c',nt_c)
    call mpas_pool_get_array(diag_physics,'mu_c',mu_c)

!... initialize the prescribed number of cloud droplets, and mu_c (parameter in the exponential of the generalized
!gamma distribution) as a function of the land-cean mask. as set in the thompson cloud microphysics scheme, nt_c
!is set to 100 per cc (100.E6 m^-3) for maritime cases and 300 per cc (300.E6 m^-3) for continental cases.
    do iCell = 1, nCellsSolve
       if(landmask(iCell) .eq. 1) then
          nt_c(iCell) = Nt_c_l
       elseif(landmask(iCell) .eq. 0) then
          nt_c(iCell) = Nt_c_o
       endif
       mu_c(iCell) = MIN(15., (1000.e6/nt_c(iCell) + 2.))
    enddo

    call physics_message('calling init_thompson_clouddroplets_forMPAS() without aerosols')
 endif

! AAJ - Set initial value of nc to 10% of the available WFA.
! This is a quick fix (01 April 2024) since the values of nc
! that are set above are to high to use as intial conditions
! with the aerosol-aware scheme.
! do k = 1, nVertLevels
!    do iCell = 1, nCellsSolve
!       scalars(index_nc,k,iCell) = max(10.e6, min(500.e6, (0.1 * scalars(index_nwfa,k,iCell))))
!    enddo
! enddo

!call mpas_log_write('--- end subroutine init_thompson_clouddroplets_forMPAS.')

 end subroutine init_thompson_clouddroplets_forMPAS 

!=================================================================================================================
 end module mpas_atmphys_init_microphysics
!=================================================================================================================

 
 
