! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_driver_convection
 use mpas_kind_types
 use mpas_pool_routines
 use mpas_timer, only : mpas_timer_start, mpas_timer_stop

 use mpas_atmphys_constants
 use mpas_atmphys_utilities
 use mpas_atmphys_vars

!wrf physics:
 use module_cu_gf
 use module_cu_kfeta
 use module_cu_tiedtke
 use module_cu_ntiedtke,only: cu_ntiedtke_driver
 use mpas_atm_dimensions, only : nVertLevels
 use mpas_dmpar, only : IO_NODE

 implicit none
 private
 public:: allocate_convection,     &
          deallocate_convection,   &
          init_convection,         &
          driver_convection,       &
          update_convection_step1, &
          update_convection_step2
 ! Abstract interface for routine used to communicate halos of fields
 ! in a named group
 !
 abstract interface
      subroutine halo_exchange_routine(domain, halo_group, ierr)

         use mpas_derived_types, only : domain_type

         type (domain_type), intent(inout) :: domain
         character(len=*), intent(in) :: halo_group
         integer, intent(out), optional :: ierr

      end subroutine halo_exchange_routine
 end interface


!MPAS driver for parameterization of convection.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines in mpas_atmphys_driver_convection:
! ----------------------------------------------
! allocate_convection     : allocate local arrays for parameterization of convection.
! deallocate_convection   : deallocate local arrays for parameterization of convection.
! init_convection         : initialization of individual convection scheme.
! driver_convection       : main driver (called from subroutine physics_driver).
! convection_from_MPAS    : initialize local arrays.
! convection_to_MPAS      : copy local arrays to MPAS arrays.
! update_convection_step1 : updates lifetime of deep convective clouds in Kain-Fritsch scheme.
! update_convection_step2 : updates accumulated precipitation output from convection schemes.
!
! WRF physics called from driver_convection:
! ------------------------------------------
! * module_cu_kfeta   : Kain-Fritsch convection scheme.
! * module_cu_tiedtke : Tiedtke convection scheme.

! add-ons and modifications to sourcecode:
! ----------------------------------------
! * removed the pre-processor option "do_hydrostatic_pressure" before call to the subroutines kf_eta_cps
!   and tiedtke. 
! * removed call to the updated Kain-Fritsch convection scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2013-05-29.
! * added the mean distance between cell centers in the call to the Kain-Fritsch parameterization of convection.
!   Laura D. Fowler (laura@ucar.edu) / 2013-08-22. 
! * in call to subroutine kf_eta_cps, replaced the variable g (that originally pointed to gravity) with gravity,
!   for simplicity.
!   Laura D. Fowler (laura@ucar.edu) / 2014-03-21.
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * Modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * in ./physics_wrf, updated the Tiedtke convection scheme to that in WRF version 3.6.1. The call to subroutine 
!   cu_tiedtke has been updated accordingly to include the sensible heat flux.
!   Laura D. Fowler (laura@ucar.edu) / 2014-09-11.
! * renamed "kain_fritsch" with "cu_kain_fritsch" and "tiedtke" with "cu_tiedtke".
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-22.
! * added the implementation of the Grell-Freitas convection scheme (option cu_grell_freitas).
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-30.
!   code cleaner, now that we are using physics packages.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-01.
! * in the call to subroutine cu_grell_freitas, added the array pratec for the calculation of rainc.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-14.
! * in the call to subroutine cu_grell_freitas, removed the argument ktop_deep.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-20.
! * updated the call to the Tiedtke parameterization of convection to WRF version 3.8.1. This option is
!   triggered using the option cu_tiedtke.
!   Laura D. Fowler (laura@ucar.edu) / 2016-08-18.
! * added the call to the "new" Tiedtke parameterization of convection from WRF version 3.8.1. This option is
!   triggered using the option cu_ntiedtke.
!   Laura D. Fowler (laura@ucar.edu) / 2016-09-20.
! * for the kain_fritsch parameterization of convection, change the definition of dx_p to match that used in the
!   Grell-Freitas and "new Tiedtke" parameterization.
!   Laura D. Fowler (laura@ucar.edu) / 2016-10-18.
! * since we removed the local variable convection_scheme from mpas_atmphys_vars.F, now defines convection_scheme
!   as a pointer to config_convection_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2017-02-16.
! * removed f_qv,f_qr,and f_qs in the calls to cu_tiedtke and cu_ntiedtke. removed f_qv and f_qc in the call to
!   kf_eta_cps.
!   Laura D. Fowler (laura@ucar.edu) / 2024-02-13.


 contains


!=================================================================================================================
 subroutine allocate_convection(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

!local variables and pointers:
 character(len=StrKIND),pointer:: convection_scheme

 integer,pointer:: config_gf_sub3d
 integer:: i,k,j

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_convection_scheme',convection_scheme)

 if(.not.allocated(cu_act_flag)) allocate(cu_act_flag(ims:ime,jms:jme)       )
 if(.not.allocated(rthcuten_p) ) allocate(rthcuten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rqvcuten_p) ) allocate(rqvcuten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rqccuten_p) ) allocate(rqccuten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rqicuten_p) ) allocate(rqicuten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(pratec_p)   ) allocate(pratec_p(ims:ime,jms:jme)          )
 if(.not.allocated(raincv_p)   ) allocate(raincv_p(ims:ime,jms:jme)          )
 if(.not.allocated(rainncv_p)  ) allocate(rainncv_p(ims:ime,jms:jme)         )

 do i = its,ite
 do j = jts,jte
    pratec_p(i,j) = 0._RKIND
    raincv_p(i,j) = 0._RKIND
 enddo
 enddo

 do i = its,ite
 do k = kts,kte
 do j = jts,jte
    rthcuten_p(i,k,j) = 0._RKIND
    rqvcuten_p(i,k,j) = 0._RKIND
    rqccuten_p(i,k,j) = 0._RKIND
    rqicuten_p(i,k,j) = 0._RKIND
 enddo
 enddo
 enddo

 convection_select: select case(convection_scheme)

    case ("cu_grell_freitas")
       if(.not.allocated(ht_p)           ) allocate(ht_p(ims:ime,jms:jme)               )
       if(.not.allocated(cubot_p)        ) allocate(cubot_p(ims:ime,jms:jme)            )
       if(.not.allocated(cutop_p)        ) allocate(cutop_p(ims:ime,jms:jme)            )
       if(.not.allocated(kpbl_p)         ) allocate(kpbl_p(ims:ime,jms:jme)             )
       if(.not.allocated(k22_shallow_p)  ) allocate(k22_shallow_p(ims:ime,jms:jme)      )
       if(.not.allocated(kbcon_shallow_p)) allocate(kbcon_shallow_p(ims:ime,jms:jme)    )
       if(.not.allocated(ktop_shallow_p) ) allocate(ktop_shallow_p(ims:ime,jms:jme)     )
       if(.not.allocated(kbot_shallow_p) ) allocate(kbot_shallow_p(ims:ime,jms:jme)     )
       if(.not.allocated(ktop_deep_p)    ) allocate(ktop_deep_p(ims:ime,jms:jme)        )

       if(.not.allocated(bilbc_p)        ) allocate(bilbc_p(ims:ime,jms:jme)            )
       if(.not.allocated(dx_p)           ) allocate(dx_p(ims:ime,jms:jme)               )
       if(.not.allocated(area_p)         ) allocate(area_p(ims:ime,jms:jme)             )
       if(.not.allocated(gsw_p)          ) allocate(gsw_p(ims:ime,jms:jme)              )
       if(.not.allocated(hfx_p)          ) allocate(hfx_p(ims:ime,jms:jme)              )
       if(.not.allocated(qfx_p)          ) allocate(qfx_p(ims:ime,jms:jme)              )
       if(.not.allocated(maxmfbl_p)      ) allocate(maxmfbl_p(ims:ime,jms:jme)          )
       if(.not.allocated(xland_p)        ) allocate(xland_p(ims:ime,jms:jme)            )

       if(.not.allocated(sigma_deep_p)   ) allocate(sigma_deep_p(ims:ime,jms:jme)       )
       if(.not.allocated(sigma_deep_f_p) ) allocate(sigma_deep_f_p(ims:ime,jms:jme)       )

       if(.not.allocated(rthblten_p)     ) allocate(rthblten_p(ims:ime,kms:kme,jms:jme) )
       if(.not.allocated(rthdynten_p)    ) allocate(rthdynten_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(rthraten_p)     ) allocate(rthraten_p(ims:ime,kms:kme,jms:jme) )
       if(.not.allocated(rqvblten_p)     ) allocate(rqvblten_p(ims:ime,kms:kme,jms:jme) )
       if(.not.allocated(rqvdynten_p)    ) allocate(rqvdynten_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(rucuten_p)      ) allocate(rucuten_p(ims:ime,kms:kme,jms:jme)  )
       if(.not.allocated(rvcuten_p)      ) allocate(rvcuten_p(ims:ime,kms:kme,jms:jme)  )

       if(.not.allocated(xmb_total_p)    ) allocate(xmb_total_p(ims:ime,jms:jme)        )
       if(.not.allocated(xmb_shallow_p)  ) allocate(xmb_shallow_p(ims:ime,jms:jme)      )

       if(.not.allocated(qccu_p)         ) allocate(qccu_p(ims:ime,kms:kme,jms:jme)     )
       if(.not.allocated(qicu_p)         ) allocate(qicu_p(ims:ime,kms:kme,jms:jme)     )

       if(.not.allocated(vcpool_p)       ) allocate(vcpool_p(ims:ime,jms:jme)           )
       if(.not.allocated(cldfrac_cu_p)   ) allocate(cldfrac_cu_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(refl10cm_cu_p)  ) allocate(refl10cm_cu_p(ims:ime,kms:kme,jms:jme))

       do i = its,ite
       do j = jts,jte
          sigma_deep_p(i,j) = 0._RKIND
          sigma_deep_f_p(i,j) = 0._RKIND
       enddo
       enddo
       call mpas_pool_get_config(configs,'config_gf_sub3d',config_gf_sub3d)
            if(.not.allocated(sub3d_rthcuten_p)) allocate(sub3d_rthcuten_p(ims:ime,kms:kme,jms:jme)  )
            if(.not.allocated(sub3d_rqvcuten_p)) allocate(sub3d_rqvcuten_p(ims:ime,kms:kme,jms:jme)  )
            if(.not.allocated(sub3d_rucuten_p))  allocate(sub3d_rucuten_p (ims:ime,kms:kme,jms:jme)  )
            if(.not.allocated(sub3d_rvcuten_p))  allocate(sub3d_rvcuten_p (ims:ime,kms:kme,jms:jme)  )
            do i = its,ite
            do k = kts,kte
            do j = jts,jte
               sub3d_rthcuten_p(i,k,j)  = 0._RKIND
               sub3d_rqvcuten_p(i,k,j)  = 0._RKIND
               sub3d_rucuten_p(i,k,j)   = 0._RKIND
               sub3d_rvcuten_p(i,k,j)   = 0._RKIND
            enddo
            enddo
            enddo
    case ("cu_kain_fritsch")
       if(.not.allocated(dx_p)       ) allocate(dx_p(ims:ime,jms:jme)              )
       if(.not.allocated(area_p)     ) allocate(area_p(ims:ime,jms:jme)            )
       if(.not.allocated(nca_p)      ) allocate(nca_p(ims:ime,jms:jme)             )
       if(.not.allocated(cubot_p)    ) allocate(cubot_p(ims:ime,jms:jme)           )
       if(.not.allocated(cutop_p)    ) allocate(cutop_p(ims:ime,jms:jme)           ) 
       if(.not.allocated(w0avg_p)    ) allocate(w0avg_p(ims:ime,kms:kme,jms:jme)   )
       if(.not.allocated(rqrcuten_p) ) allocate(rqrcuten_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(rqscuten_p) ) allocate(rqscuten_p(ims:ime,kms:kme,jms:jme))

       do i = its,ite
       do j = jts,jte
          cubot_p(i,j) = DBLE(kte+1)
          cutop_p(i,j) = DBLE(kts)
       enddo
       enddo

       do i = its,ite
       do k = kts,kte
       do j = jts,jte
          rqrcuten_p(i,k,j) = 0._RKIND
          rqscuten_p(i,k,j) = 0._RKIND
       enddo
       enddo
       enddo

    case ("cu_tiedtke","cu_ntiedtke")
       if(.not.allocated(hfx_p)        ) allocate(hfx_p(ims:ime,jms:jme)                )
       if(.not.allocated(qfx_p)        ) allocate(qfx_p(ims:ime,jms:jme)                )
       if(.not.allocated(xland_p)      ) allocate(xland_p(ims:ime,jms:jme)              )
       if(.not.allocated(rucuten_p)    ) allocate(rucuten_p(ims:ime,kms:kme,jms:jme)    )
       if(.not.allocated(rvcuten_p)    ) allocate(rvcuten_p(ims:ime,kms:kme,jms:jme)    )

       do i = its,ite
       do j = jts,jte
          hfx_p(i,j)   = 0._RKIND
          qfx_p(i,j)   = 0._RKIND
          xland_p(i,j) = 0._RKIND
       enddo
       enddo

       do i = its,ite
       do k = kts,kte
       do j = jts,jte
          rucuten_p(i,k,j)     = 0._RKIND
          rvcuten_p(i,k,j)     = 0._RKIND
       enddo
       enddo
       enddo

       cu_tiedtke_select: select case(convection_scheme)

          case ("cu_tiedtke")
             if(.not.allocated(rqvdynten_p)  ) allocate(rqvdynten_p(ims:ime,kms:kme,jms:jme)  )
             if(.not.allocated(rqvdynblten_p)) allocate(rqvdynblten_p(ims:ime,kms:kme,jms:jme))

             do i = its,ite
             do k = kts,kte
             do j = jts,jte
                rqvdynten_p(i,k,j)   = 0._RKIND
                rqvdynblten_p(i,k,j) = 0._RKIND
             enddo
             enddo
             enddo

          case("cu_ntiedtke")
             if(.not.allocated(dx_p)     ) allocate(dx_p(ims:ime,jms:jme)             )
             if(.not.allocated(rqvften_p)) allocate(rqvften_p(ims:ime,kms:kme,jms:jme))
             if(.not.allocated(rthften_p)) allocate(rthften_p(ims:ime,kms:kme,jms:jme))

             do i = its,ite
             do k = kts,kte
             do j = jts,jte
                rqvften_p(i,k,j) = 0._RKIND
                rthften_p(i,k,j) = 0._RKIND
             enddo
             enddo
             enddo

           case default

       end select cu_tiedtke_select

    case default

 end select convection_select

 end subroutine allocate_convection

!=================================================================================================================
 subroutine deallocate_convection(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

!local pointers:
 character(len=StrKIND),pointer:: convection_scheme

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_convection_scheme',convection_scheme)

 if(allocated(cu_act_flag)) deallocate(cu_act_flag)
 if(allocated(rthcuten_p) ) deallocate(rthcuten_p )
 if(allocated(rqvcuten_p) ) deallocate(rqvcuten_p )
 if(allocated(rqccuten_p) ) deallocate(rqccuten_p )
 if(allocated(rqicuten_p) ) deallocate(rqicuten_p )
 if(allocated(pratec_p)   ) deallocate(pratec_p   )
 if(allocated(raincv_p)   ) deallocate(raincv_p   )
 if(allocated(rainncv_p)  ) deallocate(rainncv_p  )

 convection_select: select case(convection_scheme)

    case ("cu_grell_freitas")
       if(allocated(ht_p)           ) deallocate(ht_p           )
       if(allocated(cubot_p)        ) deallocate(cubot_p        )
       if(allocated(cutop_p)        ) deallocate(cutop_p        )
       if(allocated(kpbl_p)         ) deallocate(kpbl_p         )
       if(allocated(k22_shallow_p)  ) deallocate(k22_shallow_p  )
       if(allocated(kbcon_shallow_p)) deallocate(kbcon_shallow_p)
       if(allocated(ktop_shallow_p) ) deallocate(ktop_shallow_p )
       if(allocated(kbot_shallow_p) ) deallocate(kbot_shallow_p )
       if(allocated(ktop_deep_p)    ) deallocate(ktop_deep_p    )

       if(allocated(bilbc_p)        ) deallocate(bilbc_p        )
       if(allocated(dx_p)           ) deallocate(dx_p           )
       if(allocated(area_p)         ) deallocate(area_p         )
       if(allocated(gsw_p)          ) deallocate(gsw_p          )
       if(allocated(hfx_p)          ) deallocate(hfx_p          )
       if(allocated(qfx_p)          ) deallocate(qfx_p          )
       if(allocated(maxmfbl_p)      ) deallocate(maxmfbl_p      )
       if(allocated(xland_p)        ) deallocate(xland_p        )
       if(allocated(rthblten_p)     ) deallocate(rthblten_p     )
       if(allocated(rthdynten_p)    ) deallocate(rthdynten_p    )
       if(allocated(rthraten_p)     ) deallocate(rthraten_p     )
       if(allocated(rqvblten_p)     ) deallocate(rqvblten_p     )
       if(allocated(rqvdynten_p)    ) deallocate(rqvdynten_p    )
       if(allocated(rucuten_p)      ) deallocate(rucuten_p      )
       if(allocated(rvcuten_p)      ) deallocate(rvcuten_p      )

       if(allocated(xmb_total_p)    ) deallocate(xmb_total_p    )
       if(allocated(xmb_shallow_p)  ) deallocate(xmb_shallow_p  )

       if(allocated(qccu_p)         ) deallocate(qccu_p         )
       if(allocated(qicu_p)         ) deallocate(qicu_p         )

       if(allocated(vcpool_p)       ) deallocate(vcpool_p       )
       if(allocated(cldfrac_cu_p)   ) deallocate(cldfrac_cu_p   )
       if(allocated(refl10cm_cu_p)  ) deallocate(refl10cm_cu_p  )
       if(allocated(sigma_deep_p)   ) deallocate(sigma_deep_p   )
       if(allocated(sigma_deep_f_p) ) deallocate(sigma_deep_f_p   )
       if(allocated(sub3d_rthcuten_p)) deallocate(sub3d_rthcuten_p)
       if(allocated(sub3d_rqvcuten_p)) deallocate(sub3d_rqvcuten_p)
       if(allocated(sub3d_rucuten_p) ) deallocate(sub3d_rucuten_p )
       if(allocated(sub3d_rvcuten_p) ) deallocate(sub3d_rvcuten_p )

    case ("cu_kain_fritsch")
       if(allocated(dx_p)         ) deallocate(dx_p         )
       if(allocated(area_p)       ) deallocate(area_p       )
       if(allocated(nca_p)        ) deallocate(nca_p        )
       if(allocated(cubot_p)      ) deallocate(cubot_p      )
       if(allocated(cutop_p)      ) deallocate(cutop_p      ) 
       if(allocated(w0avg_p)      ) deallocate(w0avg_p      )
       if(allocated(rqrcuten_p)   ) deallocate(rqrcuten_p   )
       if(allocated(rqscuten_p)   ) deallocate(rqscuten_p   )

    case ("cu_tiedtke","cu_ntiedtke")
       if(allocated(hfx_p)        ) deallocate(hfx_p        )
       if(allocated(qfx_p)        ) deallocate(qfx_p        )
       if(allocated(xland_p)      ) deallocate(xland_p      )
       if(allocated(rucuten_p)    ) deallocate(rucuten_p    )
       if(allocated(rvcuten_p)    ) deallocate(rvcuten_p    )

       cu_tiedtke_select: select case(convection_scheme)

          case ("cu_tiedtke")
             if(allocated(rqvdynten_p)  ) deallocate(rqvdynten_p  )
             if(allocated(rqvdynblten_p)) deallocate(rqvdynblten_p)

          case ("cu_ntiedtke")
             if(allocated(dx_p)     ) deallocate(dx_p)
             if(allocated(rqvften_p)) deallocate(rqvften_p)
             if(allocated(rthften_p)) deallocate(rthften_p)

          case default

       end select cu_tiedtke_select

    case default

 end select convection_select

 end subroutine deallocate_convection

!=================================================================================================================
 subroutine init_convection(dminfo,mesh,configs,diag_physics)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs
 type(dm_info),intent(in):: dminfo

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local pointers:
 logical,pointer:: config_do_restart
 character(len=StrKIND),pointer:: convection_scheme

 integer,pointer:: nCells
 real(kind=RKIND),dimension(:),pointer:: nca

!local variables:
 integer:: iCell

 integer, pointer :: config_gf_sub3d
 integer :: mynum, init_stat


!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_pool_get_config(configs,'config_do_restart'       ,config_do_restart)
 call mpas_pool_get_config(configs,'config_convection_scheme',convection_scheme)

 convection_select: select case(convection_scheme)

    case ("cu_kain_fritsch")
       call mpas_pool_get_array(diag_physics,'nca',nca)
       if(.not. config_do_restart) then
          do iCell = 1, nCells
             nca(iCell) = -100._RKIND
          enddo
       endif
       call kf_lutab(svp1,svp2,svp3,svpt0)

    case default

 end select convection_select

 end subroutine init_convection

!=================================================================================================================
 subroutine driver_convection(itimestep,configs,mesh,sfc_input,diag_physics,tend_physics,diag,dminfo,its,ite &
                             ,exchange_halo_group,block)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: sfc_input
 type(mpas_pool_type),intent(in):: diag
 type(dm_info),intent(in):: dminfo
 type(block_type), intent(in) :: block
 procedure (halo_exchange_routine) :: exchange_halo_group

 integer,intent(in):: its,ite
 integer,intent(in):: itimestep

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

!local variables and arrays:
 logical:: log_convection

 integer:: i,j
 integer:: icount,initflag,sub3d

 real(kind=RKIND):: dx

!local pointers:
 logical,pointer:: config_do_restart
 integer,pointer:: gfconv_closure_deep,gfconv_closure_shallow,config_gf_sub3d
 character(len=StrKIND),pointer:: convection_scheme
 character(len=StrKIND),pointer:: pbl_scheme
 real(kind=RKIND),pointer:: len_disp

!variables specific to Kain_Fritsch parameterization:
 logical:: warm_rain,adapt_step_flag
 integer:: ktau
 real(kind=RKIND):: curr_secs
 real(kind=RKIND):: cudt
 real(kind=RKIND):: cudtacttime

!CCPP-compliant flags:
 character(len=StrKIND):: errmsg
 integer:: errflg

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('')
!call mpas_log_write('--- enter subroutine driver_convection:')

!initialization of CCPP-compliant flags:
 errmsg = ' '
 errflg = 0

 call mpas_pool_get_config(configs,'config_gfconv_closure_deep',gfconv_closure_deep)
 call mpas_pool_get_config(configs,'config_gf_sub3d',config_gf_sub3d)
 sub3d=0
 if(config_gf_sub3d > 0) sub3d=config_gf_sub3d  !sub3d=1
 call mpas_pool_get_config(configs,'config_gfconv_closure_shallow',gfconv_closure_shallow)
 call mpas_pool_get_config(configs,'config_len_disp'         ,len_disp         )
 call mpas_pool_get_config(configs,'config_do_restart'       ,config_do_restart)
 call mpas_pool_get_config(configs,'config_convection_scheme',convection_scheme)
 call mpas_pool_get_config(configs,'config_pbl_scheme',pbl_scheme)
!call mpas_pool_get_array(mesh, 'bdyMaskCell', bdyMaskCell)
!call mpas_pool_get_array(mesh, 'specZoneMaskCell',specZoneMaskCell )

!initialize instantaneous precipitation, and copy convective tendencies from the dynamics to
!the physics grid:
 call convection_from_MPAS(dt_dyn,configs,mesh,sfc_input,diag_physics,tend_physics,diag,its,ite)

!... convert the convection time-step to minutes:
 cudt = dt_cu/60.

!... call to convection schemes:
 curr_secs   = -1
 cudtacttime = -1
 adapt_step_flag = .false.
 do j = jts, jte
 do i = its, ite
    cu_act_flag(i,j) = .false.
 enddo
 enddo

!... initialization of initflag needed in the calls cu_tiedtke and cu_ntiedtke:
 initflag = 1
 if(config_do_restart .or. itimestep > 1) initflag = 0

 convection_select: select case(convection_scheme)

    case ("cu_grell_freitas")
       call mpas_timer_start('Grell-Freitas')
       call cu_grell_freitas( &
             itimestep       = itimestep              , dt            = dt_dyn          , &
             dxCell          = dx_p                   , areaCell      = area_p          , &
             u               = u_p                    , v             = v_p             , &
             w               = w_p                    , t             = t_p             , &
             q               = qv_p                   , rho           = rho_p           , &
             p               = pres_hyd_p             , pi            = pi_p            , &
             p8w             = pres2_hyd_p            , dz8w          = dz_p            , &
             ht              = ht_p                   , xland         = xland_p         , &
             gsw             = gsw_p                  , xlv           = xlv             , &
             cp              = cp                     , g             = gravity         , &
             r_v             = R_v                    , hfx           = hfx_p           , &
             qfx             = qfx_p                  , rthblten      = rthblten_p      , &
             rqvblten        = rqvblten_p             , rthften       = rthdynten_p     , &
             rqvften         = rqvdynten_p            , rthraten      = rthraten_p      , &
             kpbl            = kpbl_p                 , raincv        = raincv_p        , &
             pratec          = pratec_p               , htop          = cutop_p         , &
             hbot            = cubot_p                , k22_shallow   = k22_shallow_p   , & 
             kbcon_shallow   = kbcon_shallow_p        , ktop_shallow  = ktop_shallow_p  , & 
             xmb_total       = xmb_total_p            , xmb_shallow   = xmb_shallow_p   , &
             gdc             = qccu_p                 , gdc2          = qicu_p          , &
             rthcuten        = rthcuten_p             , rqvcuten      = rqvcuten_p      , & 
             rqccuten        = rqccuten_p             , rqicuten      = rqicuten_p      , &
             rucuten         = rucuten_p              , rvcuten       = rvcuten_p       , &
             sub3d_rthcuten  = sub3d_rthcuten_p       , sub3d_rqvcuten= sub3d_rqvcuten_p, &
             sub3d_rucuten   = sub3d_rucuten_p        , sub3d_rvcuten = sub3d_rvcuten_p , &
             ichoice_deep    = gfconv_closure_deep    , sig_deep      = sigma_deep_p    , &
             sig_deep_far    = sigma_deep_f_p                                           , &
             ichoice_shallow = gfconv_closure_shallow , bilbc         = bilbc_p         , &
             ishallow_g3     = ishallow               , sub3d         = sub3d           , &
             pbl_scheme      = pbl_scheme             , maxmf         = maxmfbl_p       , &
             qc3d            = qc_p                   , qi3d          = qi_p            , &
             qc_cu           = qccu_p                 , qi_cu         = qicu_p          , &
             cldfrac_cu      = cldfrac_cu_p           , refl10cm_cu   = refl10cm_cu_p   , &
             rainncv         = rainncv_p              , vcpool        = vcpool_p        , &
             ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde      , &
             ims = ims , ime = ime , jms = jms , jme = jme , kms = kds , kme = kme      , &
             its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte        &
                            )
       call mpas_timer_stop('Grell-Freitas')

    case ("cu_kain_fritsch")
       if(itimestep == 1) then
          ktau = itimestep
       else
          ktau = itimestep + 1
       endif

       call mpas_timer_start('Kain-Fritsch')
       call kf_eta_cps ( &
             pcps            = pres_hyd_p      , t         = t_p        ,            &
             dt              = dt_dyn          , ktau      = ktau       ,            &
             dxCell          = dx_p            , areaCell  = area_p     ,            &
             cudt            = cudt            , curr_secs = curr_secs  ,            & 
             adapt_step_flag = adapt_step_flag , rho       = rho_p      ,            &
             raincv          = raincv_p        , pratec    = pratec_p   ,            &
             nca             = nca_p           , u         = u_p        ,            &
             v               = v_p             , th        = th_p       ,            &
             pi              = pi_p            , w         = w_p        ,            &
             dz8w            = dz_p            , w0avg     = w0avg_p    ,            &
             xlv0            = xlv0            , xlv1      = xlv1       ,            &
             xls0            = xls0            , xls1      = xls1       ,            &
             cp              = cp              , r         = r_d        ,            &
             g               = gravity         , ep1       = ep_1       ,            &
             ep2             = ep_2            , svp1      = svp1       ,            &
             svp2            = svp2            , svp3      = svp3       ,            &
             svpt0           = svpt0           , stepcu    = n_cu       ,            &
             cu_act_flag     = cu_act_flag     , warm_rain = warm_rain  ,            &
             cutop           = cutop_p         , cubot     = cubot_p    ,            &
             qv              = qv_p            , f_qr      = f_qr       ,            &
             f_qi            = f_qi            , f_qs      = f_qs       ,            &
             rthcuten        = rthcuten_p      , rqvcuten  = rqvcuten_p ,            &
             rqccuten        = rqccuten_p      , rqrcuten  = rqrcuten_p ,            &
             rqicuten        = rqicuten_p      , rqscuten  = rqscuten_p ,            &
             ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde , &
             ims = ims , ime = ime , jms = jms , jme = jme , kms = kds , kme = kme , &
             its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte   &
                        )
       call mpas_timer_stop('Kain-Fritsch')

    case("cu_tiedtke")
       call mpas_timer_start('Tiedtke')
       call cu_tiedtke( &
             pcps        = pres_hyd_p  , p8w       = pres2_hyd_p   ,                 &
             znu         = znu_hyd_p   , t3d       = t_p           ,                 &
             dt          = dt_dyn      , itimestep = initflag      ,                 &
             stepcu      = n_cu        , raincv    = raincv_p      ,                 &
             pratec      = pratec_p    , qfx       = qfx_p         ,                 &
             u3d         = u_p         , v3d       = v_p           ,                 &
             w           = w_p         , qv3d      = qv_p          ,                 &
             qc3d        = qc_p        , qi3d      = qi_p          ,                 &
             pi3d        = pi_p        , rho3d     = rho_p         ,                 &
             qvften      = rqvdynten_p , qvpblten  = rqvdynblten_p ,                 &
             dz8w        = dz_p        , xland     = xland_p       ,                 &
             cu_act_flag = cu_act_flag , f_qc      = f_qc          ,                 &
             f_qi        = f_qi        , rthcuten  = rthcuten_p    ,                 &
             rqvcuten    = rqvcuten_p  , rqccuten  = rqccuten_p    ,                 &
             rqicuten    = rqicuten_p  , rucuten   = rucuten_p     ,                 &
             rvcuten     = rvcuten_p   ,                                             &
             ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde , &
             ims = ims , ime = ime , jms = jms , jme = jme , kms = kds , kme = kme , &
             its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte   &
                      )
       call mpas_timer_stop('Tiedtke')

    case("cu_ntiedtke")
       call mpas_timer_start('cu_ntiedtke')
       call cu_ntiedtke_driver( &
             pcps     = pres_hyd_p , p8w         = pres2_hyd_p ,                     &
             t3d      = t_p        , dz8w        = dz_p        ,                     &
             dt       = dt_dyn     , itimestep   = initflag    ,                     &
             stepcu   = n_cu       , raincv      = raincv_p    ,                     &
             pratec   = pratec_p   , qfx         = qfx_p       ,                     &
             hfx      = hfx_p      , xland       = xland_p     ,                     &
             dx       = dx_p       , u3d         = u_p         ,                     &
             v3d      = v_p        , w           = w_p         ,                     &
             qv3d     = qv_p       , qc3d        = qc_p        ,                     &
             qi3d     = qi_p       , pi3d        = pi_p        ,                     &
             rho3d    = rho_p      , qvften      = rqvften_p   ,                     &
             thften   = rthften_p  , cu_act_flag = cu_act_flag ,                     &
             f_qc     = f_qc       , f_qi        = f_qi        ,                     &
             rthcuten = rthcuten_p , rqvcuten    = rqvcuten_p  ,                     &
             rqccuten = rqccuten_p , rqicuten    = rqicuten_p  ,                     &
             rucuten  = rucuten_p  , rvcuten     = rvcuten_p   ,                     &
             grav     = gravity    , xlf         = xlf         ,                     &
             xls      = xls        , xlv         = xlv         ,                     &
             rd       = R_d        , rv          = R_v         ,                     &
             cp       = cp         ,                                                 &
             errmsg   = errmsg     , errflg      = errflg      ,                     &
             ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde , &
             ims = ims , ime = ime , jms = jms , jme = jme , kms = kds , kme = kme , &
             its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte   &
                              )
       call mpas_timer_stop('cu_ntiedtke')

    case default

 end select convection_select

!copy instantaneous and accumulated precipitation, convective tendencies, and "other" arrays
!specific to convection parameterization back to the dynamics grid:
 call convection_to_MPAS(configs,mesh,diag_physics,tend_physics,its,ite,exchange_halo_group,block)

!call mpas_log_write('--- end subroutine driver_convection.')
!call mpas_log_write('')

 end subroutine driver_convection

!=================================================================================================================
 subroutine convection_from_MPAS(dt_dyn,configs,mesh,sfc_input,diag_physics,tend_physics,diag,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: sfc_input
 type(mpas_pool_type),intent(in):: diag_physics
 type(mpas_pool_type),intent(in):: tend_physics
 type(mpas_pool_type),intent(in):: diag

 integer,intent(in):: its,ite
 real(kind=RKIND),intent(in):: dt_dyn

!local pointers:
 character(len=StrKIND),pointer:: convection_scheme
 integer,dimension(:),pointer:: kpbl,k22_shallow,kbcon_shallow,ktop_shallow,ktop_deep
 integer,dimension(:),pointer:: bdyMaskCell
 real(kind=RKIND),dimension(:),pointer  :: areaCell,meshDensity
 real(kind=RKIND),dimension(:),pointer  :: nca,cubot,cutop,cuprec,rainncv,raincv,vcpool
 real(kind=RKIND),dimension(:),pointer  :: gsw,hfx,qfx,xland
 real(kind=RKIND),dimension(:),pointer  :: xmb_total,xmb_shallow,sigma_deep,sigma_deep_far
 real(kind=RKIND),dimension(:,:),pointer:: zgrid
 real(kind=RKIND),dimension(:,:),pointer:: qc_cu,qi_cu
 real(kind=RKIND),dimension(:,:),pointer:: w0avg
 real(kind=RKIND),dimension(:,:),pointer:: rthcuten,rqvcuten,rqccuten,rqicuten,rqrcuten,rqscuten 
 real(kind=RKIND),dimension(:,:),pointer:: rthblten,rthdynten,rthratenlw,rthratensw
 real(kind=RKIND),dimension(:,:),pointer:: rqvblten,rqvdynten,rucuten,rvcuten
 real(kind=RKIND),dimension(:,:),pointer:: cldfrac_cu,refl10cm_cu

!local variables:
 integer:: i,j,k
 integer:: iEdge
 real(kind=RKIND),pointer:: len_disp
 integer,pointer:: config_gf_sub3d
 logical result

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_convection_scheme',convection_scheme)

 call mpas_pool_get_array(diag_physics,'cuprec',cuprec)
 call mpas_pool_get_array(diag_physics,'raincv',raincv)
 call mpas_pool_get_array(diag_physics,'rainncv',rainncv)

 call mpas_pool_get_array(tend_physics,'rthcuten',rthcuten)
 call mpas_pool_get_array(tend_physics,'rqvcuten',rqvcuten)
 call mpas_pool_get_array(tend_physics,'rqccuten',rqccuten)
 call mpas_pool_get_array(tend_physics,'rqicuten',rqicuten)

 do j = jts,jte
 do i = its,ite
    raincv_p(i,j) = raincv(i)
    rainncv_p(i,j)= rainncv(i)
    pratec_p(i,j) = cuprec(i)
    do k = kts,kte
       rthcuten_p(i,k,j) = rthcuten(k,i)
       rqvcuten_p(i,k,j) = rqvcuten(k,i)
       rqccuten_p(i,k,j) = rqccuten(k,i)
       rqicuten_p(i,k,j) = rqicuten(k,i)
    enddo
 enddo
 enddo

 convection_select: select case(convection_scheme)

    case ("cu_grell_freitas")
       call mpas_pool_get_config(configs,'config_gf_sub3d',config_gf_sub3d)
       call mpas_pool_get_config(configs,'config_len_disp',len_disp)

       call mpas_pool_get_array(mesh,'areaCell'   ,areaCell   )
       call mpas_pool_get_array(mesh,'meshDensity',meshDensity)
       call mpas_pool_get_array(mesh,'zgrid'      ,zgrid      )
       call mpas_pool_get_array(mesh,'bdyMaskCell',bdyMaskCell) 
       call mpas_pool_get_array(sfc_input,'xland',xland)

       call mpas_pool_get_array(diag_physics,'gsw'          ,gsw          )
       call mpas_pool_get_array(diag_physics,'hfx'          ,hfx          )
       call mpas_pool_get_array(diag_physics,'qfx'          ,qfx          )
       call mpas_pool_get_array(diag_physics,'kpbl'         ,kpbl         )
       call mpas_pool_get_array(diag_physics,'cubot'        ,cubot        )
       call mpas_pool_get_array(diag_physics,'cutop'        ,cutop        )
       call mpas_pool_get_array(diag_physics,'xmb_total'    ,xmb_total    )
       call mpas_pool_get_array(diag_physics,'xmb_shallow'  ,xmb_shallow  )
       call mpas_pool_get_array(diag_physics,'k22_shallow'  ,k22_shallow  )
       call mpas_pool_get_array(diag_physics,'kbcon_shallow',kbcon_shallow)
       call mpas_pool_get_array(diag_physics,'ktop_shallow' ,ktop_shallow )
       call mpas_pool_get_array(diag_physics,'ktop_deep'    ,ktop_deep    )
       call mpas_pool_get_array(diag_physics,'qc_cu'        ,qc_cu        )
       call mpas_pool_get_array(diag_physics,'qi_cu'        ,qi_cu        )
       call mpas_pool_get_array(diag_physics,'vcpool'       ,vcpool       )
       call mpas_pool_get_array(diag_physics,'cldfrac_cu'   ,cldfrac_cu   )
       call mpas_pool_get_array(diag_physics,'refl10cm_cu'  ,refl10cm_cu  )
       call mpas_pool_get_array(diag_physics,'sigma_deep'   ,sigma_deep   )
       call mpas_pool_get_array(diag_physics,'sigma_deep_far'   ,sigma_deep_far )

       call mpas_pool_get_array(tend_physics,'rqvblten'  ,rqvblten  )
       call mpas_pool_get_array(tend_physics,'rqvdynten' ,rqvdynten )
       call mpas_pool_get_array(tend_physics,'rthblten'  ,rthblten  )
       call mpas_pool_get_array(tend_physics,'rthdynten' ,rthdynten )
       call mpas_pool_get_array(tend_physics,'rthratenlw',rthratenlw)
       call mpas_pool_get_array(tend_physics,'rthratensw',rthratensw)
       call mpas_pool_get_array(tend_physics,'rucuten'   ,rucuten   )
       call mpas_pool_get_array(tend_physics,'rvcuten'   ,rvcuten   )

       do j = jts,jte
       do i = its,ite
          dx_p(i,j) = len_disp / meshDensity(i)**0.25
          area_p(i,j)          = areaCell(i)
          ht_p(i,j)            = zgrid(1,i)
          xland_p(i,j)         = xland(i)
          gsw_p(i,j)           = gsw(i)
          hfx_p(i,j)           = hfx(i)
          qfx_p(i,j)           = qfx(i)
          kpbl_p(i,j)          = kpbl(i)
          cubot_p(i,j)         = cubot(i)
          cutop_p(i,j)         = cutop(i)
          bilbc_p(i,j)         = bdyMaskCell(i)
          xmb_total_p(i,j)     = xmb_total(i)
          xmb_shallow_p(i,j)   = xmb_shallow(i)

          k22_shallow_p(i,j)   = k22_shallow(i)
          kbcon_shallow_p(i,j) = kbcon_shallow(i)
          ktop_shallow_p(i,j)  = ktop_shallow(i)
          ktop_deep_p(i,j)     = ktop_deep(i)
          vcpool_p(i,j)        = vcpool(i)
          sigma_deep_p(i,j)    = sigma_deep(i)
          sigma_deep_f_p(i,j)  = sigma_deep_far(i)

          do k = kts,kte
             qccu_p(i,k,j)      = qc_cu(k,i)
             qicu_p(i,k,j)      = qi_cu(k,i)

             rthblten_p(i,k,j)  = rthblten(k,i)
             rthdynten_p(i,k,j) = rthdynten(k,i)
             rthraten_p(i,k,j)  = rthratenlw(k,i) + rthratensw(k,i)
             rqvblten_p(i,k,j)  = rqvblten(k,i)
             rqvdynten_p(i,k,j) = rqvdynten(k,i)
             rucuten_p(i,k,j)   = rucuten(k,i)
             rvcuten_p(i,k,j)   = rvcuten(k,i)

             cldfrac_cu_p(i,k,j)= cldfrac_cu(k,i)
             refl10cm_cu_p(i,k,j)= refl10cm_cu(k,i)
          enddo
       enddo
       enddo

    case ("cu_kain_fritsch")
       call mpas_pool_get_config(configs,'config_len_disp',len_disp)

       call mpas_pool_get_array(mesh,'areaCell',areaCell)
       call mpas_pool_get_array(mesh,'meshDensity',meshDensity)

       call mpas_pool_get_array(diag_physics,'cubot'   ,cubot   )
       call mpas_pool_get_array(diag_physics,'cutop'   ,cutop   )
       call mpas_pool_get_array(diag_physics,'nca'     ,nca     )
       call mpas_pool_get_array(diag_physics,'w0avg'   ,w0avg   )

       call mpas_pool_get_array(tend_physics,'rqrcuten',rqrcuten)
       call mpas_pool_get_array(tend_physics,'rqscuten',rqscuten)

       do j = jts,jte
       do i = its,ite
          !area of grid-cell:       
          dx_p(i,j)    = len_disp / meshDensity(i)**0.25
          area_p(i,j)  = areaCell(i)
          cubot_p(i,j) = cubot(i)
          cutop_p(i,j) = cutop(i)

          do k = kts,kte
             rqrcuten_p(i,k,j) = rqrcuten(k,i)
             rqscuten_p(i,k,j) = rqscuten(k,i)
          enddo

          !decreases the characteristic time period that convection remains active. When nca_p
          !becomes less than the convective timestep, convective tendencies and precipitation
          !are reset to zero (note that this is also done in subroutine kf_eta_cps).
          nca_p(i,j) = nca(i)

          do k = kts,kte
             w0avg_p(i,k,j) = w0avg(k,i)
          enddo
       enddo
       enddo

    case ("cu_tiedtke","cu_ntiedtke")
       call mpas_pool_get_array(sfc_input,'xland',xland)
       call mpas_pool_get_array(diag_physics,'hfx',hfx)
       call mpas_pool_get_array(diag_physics,'qfx',qfx)

       call mpas_pool_get_array(tend_physics,'rqvblten' ,rqvblten )
       call mpas_pool_get_array(tend_physics,'rqvdynten',rqvdynten)
       call mpas_pool_get_array(tend_physics,'rucuten'  ,rucuten  )
       call mpas_pool_get_array(tend_physics,'rvcuten'  ,rvcuten  )

       do j = jts,jte
       do i = its,ite
          xland_p(i,j) = xland(i)
          hfx_p(i,j)   = hfx(i)
          qfx_p(i,j)   = qfx(i)
       enddo

       do k = kts,kte
       do i = its,ite
          rucuten_p(i,k,j) = rucuten(k,i)
          rvcuten_p(i,k,j) = rvcuten(k,i)
       enddo
       enddo
       enddo

       cu_tiedtke_select: select case(convection_scheme)

          case ("cu_tiedtke")
             do j = jts,jte
             do k = kts,kte
             do i = its,ite
                rqvdynblten_p(i,k,j) = rqvblten(k,i) 
                rqvdynten_p(i,k,j)   = rqvdynten(k,i)
             enddo
             enddo
             enddo
 
          case("cu_ntiedtke")
             call mpas_pool_get_config(configs,'config_len_disp',len_disp)
             call mpas_pool_get_array(mesh,'meshDensity',meshDensity)

             call mpas_pool_get_array(tend_physics,'rqvdynten' ,rqvdynten )
             call mpas_pool_get_array(tend_physics,'rqvblten'  ,rqvblten  )
             call mpas_pool_get_array(tend_physics,'rthdynten' ,rthdynten )
             call mpas_pool_get_array(tend_physics,'rthblten'  ,rthblten  )
             call mpas_pool_get_array(tend_physics,'rthratenlw',rthratenlw)
             call mpas_pool_get_array(tend_physics,'rthratensw',rthratensw)
             
             do j = jts,jte
             do i = its,ite
                dx_p(i,j) = len_disp / meshDensity(i)**0.25
             enddo

             do k = kts,kte
             do i = its,ite
                rqvften_p(i,k,j) = rqvdynten(k,i) + rqvblten(k,i)
                rthften_p(i,k,j) = (rthdynten(k,i) + rthblten(k,i) &
                                 + rthratenlw(k,i) + rthratensw(k,i)) * pi_p(i,k,j) 
             enddo
             enddo
             enddo

          case default

       end select cu_tiedtke_select

    case default

 end select convection_select

 end subroutine convection_from_MPAS

!=================================================================================================================
 subroutine convection_to_MPAS(configs,mesh,diag_physics,tend_physics,its,ite,exchange_halo_group,block)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(block_type), intent(in) :: block
 procedure (halo_exchange_routine) :: exchange_halo_group
 integer,intent(in):: its,ite

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

!local variables:
 integer:: i,k,j
 integer:: iCell,ineighbor,ineighbor_far

!local pointers:
 character(len=StrKIND),pointer:: convection_scheme
 integer,dimension(:),pointer:: k22_shallow,kbcon_shallow,ktop_shallow,ktop_deep
 real(kind=RKIND),dimension(:),pointer  :: nca,cubot,cutop,cuprec,rainncv,raincv,vcpool
 real(kind=RKIND),dimension(:),pointer  :: sigma_deep,sigma_deep_far
 real(kind=RKIND),dimension(:),pointer  :: xmb_total,xmb_shallow
 real(kind=RKIND),dimension(:,:),pointer:: qc_cu,qi_cu
 real(kind=RKIND),dimension(:,:),pointer:: w0avg
 real(kind=RKIND),dimension(:,:),pointer:: rthcuten,rqvcuten,rqccuten,rqicuten,rqrcuten,rqscuten
 real(kind=RKIND),dimension(:,:),pointer:: rucuten,rvcuten
 real(kind=RKIND),dimension(:,:),pointer:: cldfrac_cu,refl10cm_cu
 real(kind=RKIND),dimension(:,:),pointer:: sub3d_rucuten,sub3d_rvcuten,sub3d_rthcuten,sub3d_rqvcuten
 integer,pointer:: config_gf_sub3d,nCellsSolve
 integer, dimension(:), pointer:: nEdgesOnCell,bdyMaskCell,ncellsOnCellsOnCell
 integer, dimension(:,:), pointer:: cellsOnCell,cellsOnCellsOnCell
 logical result

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_convection_scheme',convection_scheme)

 call mpas_pool_get_array(diag_physics,'cuprec',cuprec)
 call mpas_pool_get_array(diag_physics,'raincv',raincv)
 call mpas_pool_get_array(diag_physics,'rainncv',rainncv)

 call mpas_pool_get_array(tend_physics,'rthcuten',rthcuten)
 call mpas_pool_get_array(tend_physics,'rqvcuten',rqvcuten)
 call mpas_pool_get_array(tend_physics,'rqccuten',rqccuten)
 call mpas_pool_get_array(tend_physics,'rqicuten',rqicuten)

 do j = jts,jte
 do i = its,ite
    raincv(i) = raincv_p(i,j)
    cuprec(i) = pratec_p(i,j)
    do k = kts, kte
       rthcuten(k,i) = rthcuten_p(i,k,j)
       rqvcuten(k,i) = rqvcuten_p(i,k,j)
       rqccuten(k,i) = rqccuten_p(i,k,j)
       rqicuten(k,i) = rqicuten_p(i,k,j)
    enddo
 enddo
 enddo

 convection_select: select case(convection_scheme)

    case ("cu_grell_freitas")
       call mpas_pool_get_array(diag_physics,'cubot'        ,cubot        )
       call mpas_pool_get_array(diag_physics,'cutop'        ,cutop        )
       call mpas_pool_get_array(diag_physics,'xmb_total'    ,xmb_total    )
       call mpas_pool_get_array(diag_physics,'xmb_shallow'  ,xmb_shallow  )
       call mpas_pool_get_array(diag_physics,'k22_shallow'  ,k22_shallow  )
       call mpas_pool_get_array(diag_physics,'kbcon_shallow',kbcon_shallow)
       call mpas_pool_get_array(diag_physics,'ktop_shallow' ,ktop_shallow )
       call mpas_pool_get_array(diag_physics,'ktop_deep'    ,ktop_deep    )
       call mpas_pool_get_array(diag_physics,'qc_cu'        ,qc_cu        )
       call mpas_pool_get_array(diag_physics,'qi_cu'        ,qi_cu        )
       call mpas_pool_get_array(diag_physics,'cldfrac_cu'   ,cldfrac_cu   )
       call mpas_pool_get_array(diag_physics,'refl10cm_cu'  ,refl10cm_cu  )

       call mpas_pool_get_array(diag_physics,'vcpool'       ,vcpool       )
       call mpas_pool_get_array(diag_physics,'sigma_deep'  ,sigma_deep    )
       call mpas_pool_get_array(diag_physics,'sigma_deep_far',sigma_deep_far  )
       call mpas_pool_get_array(tend_physics,'rucuten',rucuten)
       call mpas_pool_get_array(tend_physics,'rvcuten',rvcuten)
       call mpas_pool_get_array(mesh,'bdyMaskCell',bdyMaskCell)

       do j = jts,jte
       do i = its,ite
          cubot(i)         = cubot_p(i,j)
          cutop(i)         = cutop_p(i,j)
          xmb_total(i)     = xmb_total_p(i,j)
          xmb_shallow(i)   = xmb_shallow_p(i,j)
          k22_shallow(i)   = k22_shallow_p(i,j)
          kbcon_shallow(i) = kbcon_shallow_p(i,j)
          ktop_shallow(i)  = ktop_shallow_p(i,j)
          ktop_deep(i)     = ktop_deep_p(i,j)
          vcpool(i)        = vcpool_p(i,j)
          sigma_deep(i)    = sigma_deep_p(i,j)
          sigma_deep_far(i)= sigma_deep_f_p(i,j)

          do k = kts,kte
             qc_cu(k,i)    = qccu_p(i,k,j)
             qi_cu(k,i)    = qicu_p(i,k,j)
             rucuten(k,i)  = rucuten_p(i,k,j)
             rvcuten(k,i)  = rvcuten_p(i,k,j)

             cldfrac_cu (k,i)  = cldfrac_cu_p(i,k,j)
             refl10cm_cu(k,i)  = refl10cm_cu_p(i,k,j)
          enddo
       enddo
       enddo

       !-- including option of the env subsidence lateral spreading
       call mpas_pool_get_config(configs,'config_gf_sub3d',config_gf_sub3d)
       call mpas_pool_get_array(tend_physics,'sub3d_rucuten'  ,sub3d_rucuten )
       call mpas_pool_get_array(tend_physics,'sub3d_rvcuten'  ,sub3d_rvcuten )
       call mpas_pool_get_array(tend_physics,'sub3d_rqvcuten' ,sub3d_rqvcuten)
       call mpas_pool_get_array(tend_physics,'sub3d_rthcuten' ,sub3d_rthcuten)

            do j = jts,jte
            do i = its,ite
               do k = kts,kte
                  sub3d_rthcuten(k,i) = sub3d_rthcuten_p(i,k,j)
                  sub3d_rqvcuten(k,i) = sub3d_rqvcuten_p(i,k,j)
                  sub3d_rucuten(k,i)  = sub3d_rucuten_p(i,k,j)
                  sub3d_rvcuten(k,i)  = sub3d_rvcuten_p(i,k,j)
               enddo
            enddo
            enddo
            !
            !
       call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
       call mpas_pool_get_array(mesh, 'cellsOnCellsOnCell', cellsOnCellsonCell)
       call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
       call mpas_pool_get_array(mesh, 'nCellsOnCellsOnCell',ncellsOnCellsonCell)
       call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
              if(config_gf_sub3d > 0) then
            !
!
            !--- part 2: broadcast the subsidence tendencies

            call exchange_halo_group(block % domain, 'physics:sub3d_cuten')

            !-- part 3: spread the subsidence tendencies to the 1st neighbors
            do iCell = 1, nCellsSolve
            ! do not do this at the boundary of the regional domain!
            !if(bdyMaskCell(icell) == 7 )cycle
            if(bdyMaskCell(iCell) > 4 )cycle
              do k = 1, nVertLevels
                 rthcuten(k,iCell) = rthcuten(k,iCell) + (1._RKIND - sigma_deep(iCell)) * sub3d_rthcuten(k,iCell)
                 rqvcuten(k,iCell) = rqvcuten(k,iCell) + (1._RKIND - sigma_deep(iCell)) * sub3d_rqvcuten(k,iCell)
                 rucuten(k,iCell)  = rucuten(k,iCell)  + (1._RKIND - sigma_deep(iCell)) * sub3d_rucuten(k,iCell)
                 rvcuten(k,iCell)  = rvcuten(k,iCell)  + (1._RKIND - sigma_deep(iCell)) * sub3d_rvcuten(k,iCell)
                 do i = 1, nEdgesOnCell(iCell)
                   ineighbor = cellsOnCell(i,iCell)
                   ! do not do this if there is no convection in that neighbor cell!
                   ! ineighbor here is the i-th cell from the first row, which should have CellsOnCells elements
                   if(sigma_deep(ineighbor) == 0.)cycle
! do this for first row cells
! 1-sigma_deep-sigma_deep_far what comes in from the first row neighbor,
! 1-sigma_deep_far            what xomes in from second row neighbor

                   rthcuten(k,iCell) = rthcuten(k,iCell) + &
                                       (sigma_deep(ineighbor)- sigma_deep_far(ineighbor)) &
                                       *sub3d_rthcuten(k,ineighbor)/ &
                                       real(nEdgesOnCell(icell),RKIND)
                   rqvcuten(k,iCell) = rqvcuten(k,iCell) + &
                                       (sigma_deep(ineighbor)- sigma_deep_far(ineighbor)) &
                                       *sub3d_rqvcuten(k,ineighbor)/ &
                                       real(nEdgesOnCell(icell),RKIND)
                   rucuten(k,iCell) = rucuten(k,iCell) + &
                                       (sigma_deep(ineighbor) - sigma_deep_far(ineighbor)) &
                                      *sub3d_rucuten(k,ineighbor)/ &
                                       real(nEdgesOnCell(icell),RKIND)
                   rvcuten(k,iCell) = rvcuten(k,iCell) + &
                                       (sigma_deep(ineighbor)- sigma_deep_far(ineighbor))  &
                                       *sub3d_rvcuten(k,ineighbor)/ &
                                       real(nEdgesOnCell(icell),RKIND)
                   enddo ! i, first row
! now for second row for this cell 
!
                   do j=1,nCellsOnCellsOnCell(icell)
                   ! ineighbor_far here is the j-th cell from the second row, which should have nCellsOnCellsOnCell elements
                   ineighbor_far = cellsOnCellsOnCell(j,iCell)
                   rthcuten(k,iCell) = rthcuten(k,iCell) + &
                                       sigma_deep_far(ineighbor_far)*sub3d_rthcuten(k,ineighbor_far)/ &
                                       real(nCellsOnCellsOnCell(icell),RKIND)                               
                   rqvcuten(k,iCell) = rqvcuten(k,iCell) + &
                                       sigma_deep_far(ineighbor_far)*sub3d_rqvcuten(k,ineighbor_far)/ &
                                       real(nCellsOnCellsOnCell(icell),RKIND)
                   rucuten(k,iCell) = rucuten(k,iCell) + &
                                       sigma_deep_far(ineighbor_far)*sub3d_rucuten(k,ineighbor_far)/ &
                                       real(nCellsOnCellsOnCell(icell),RKIND)
                   rvcuten(k,iCell) = rvcuten(k,iCell) + &
                                       sigma_deep_far(ineighbor_far)*sub3d_rvcuten(k,ineighbor_far)/ &
                                       real(nCellsOnCellsOnCell(icell),RKIND)
                   enddo  ! j - second row
              enddo ! k
            enddo ! icell
       else !if (bdyMaskCell(icell) > 4  .or. config_gf_sub3d == 0 ) then
            do i = 1, nCellsSolve
            do k = 1, nVertLevels
                 rucuten(k,i)  = rucuten(k,i)  + sub3d_rucuten(k,i)
                 rvcuten(k,i)  = rvcuten(k,i)  + sub3d_rvcuten(k,i)
                 rthcuten(k,i) = rthcuten(k,i) + sub3d_rthcuten(k,i)
                 rqvcuten(k,i) = rqvcuten(k,i) + sub3d_rqvcuten(k,i)
            enddo
            enddo

       endif

    case ("cu_kain_fritsch")
       call mpas_pool_get_array(diag_physics,'cubot',cubot)
       call mpas_pool_get_array(diag_physics,'cutop',cutop)
       call mpas_pool_get_array(diag_physics,'nca'  ,nca  )
       call mpas_pool_get_array(diag_physics,'w0avg',w0avg)

       call mpas_pool_get_array(tend_physics,'rqrcuten',rqrcuten)
       call mpas_pool_get_array(tend_physics,'rqscuten',rqscuten)

       do j = jts,jte
       do i = its,ite
          cubot(i) = cubot_p(i,j)
          cutop(i) = cutop_p(i,j)
          nca(i)   = nca_p(i,j)
          do k = kts, kte
             w0avg(k,i)    = w0avg_p(i,k,j)
             rqrcuten(k,i) = rqrcuten_p(i,k,j)
             rqscuten(k,i) = rqscuten_p(i,k,j)
          enddo                          
       enddo
       enddo

    case ("cu_tiedtke","cu_ntiedtke")
       call mpas_pool_get_array(tend_physics,'rucuten',rucuten)
       call mpas_pool_get_array(tend_physics,'rvcuten',rvcuten)

       do j = jts,jte
       do k = kts,kte
       do i = its,ite
          rucuten(k,i) = rucuten_p(i,k,j)
          rvcuten(k,i) = rvcuten_p(i,k,j)
       enddo
       enddo
       enddo

    case default

 end select convection_select

 end subroutine convection_to_MPAS

!=================================================================================================================
 subroutine update_convection_step1(configs,diag_physics,tend_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 integer,intent(in):: its,ite

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

!local pointers:
 character(len=StrKIND),pointer:: convection_scheme
 real(kind=RKIND),dimension(:),pointer  :: nca,cubot,cutop,cuprec,raincv
 real(kind=RKIND),dimension(:,:),pointer:: rthcuten,rqvcuten,rqccuten,rqicuten,rqrcuten,rqscuten

!local variables and arrays:
 integer:: i,k

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_convection_scheme',convection_scheme)

 convection_select: select case(convection_scheme)

    case ("cu_kain_fritsch")
      call mpas_pool_get_array(diag_physics,'nca'   ,nca   )
      call mpas_pool_get_array(diag_physics,'cubot' ,cubot )
      call mpas_pool_get_array(diag_physics,'cutop' ,cutop )
      call mpas_pool_get_array(diag_physics,'cuprec',cuprec)
      call mpas_pool_get_array(diag_physics,'raincv',raincv)

      call mpas_pool_get_array(tend_physics,'rthcuten',rthcuten)
      call mpas_pool_get_array(tend_physics,'rqvcuten',rqvcuten)
      call mpas_pool_get_array(tend_physics,'rqccuten',rqccuten)
      call mpas_pool_get_array(tend_physics,'rqicuten',rqicuten)
      call mpas_pool_get_array(tend_physics,'rqrcuten',rqrcuten)
      call mpas_pool_get_array(tend_physics,'rqscuten',rqscuten)

       do i = its, ite
          !decreases the characteristic time period that convection remains active. When nca_p
          !becomes less than the convective timestep, convective tendencies and precipitation
          !are reset to zero (note that this is also done in subroutine kf_eta_cps).
          if(nca(i) .gt. 0.) then
             nca(i) = nca(i) - dt_dyn
             
             if(nca(i) .lt. 0.5*dt_dyn) then
                do k = kts,kte
                   rthcuten(k,i) = 0._RKIND
                   rqvcuten(k,i) = 0._RKIND
                   rqccuten(k,i) = 0._RKIND
                   rqrcuten(k,i) = 0._RKIND
                   rqicuten(k,i) = 0._RKIND
                   rqscuten(k,i) = 0._RKIND
                enddo
                raincv(i) = 0._RKIND
                cuprec(i) = 0._RKIND
                cubot(i)  = kte+1
                cutop(i)  = kts
             endif
          endif
       enddo

    case default

 end select convection_select

 end subroutine update_convection_step1

!=================================================================================================================
 subroutine update_convection_step2(configs,diag_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 integer,intent(in):: its,ite

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local pointers:
 integer,dimension(:),pointer:: i_rainc

 real(kind=RKIND),pointer:: bucket_rainc
 real(kind=RKIND),dimension(:),pointer:: cuprec,rainc

!local variables and arrays:
 integer:: i

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_bucket_rainc',bucket_rainc)

 call mpas_pool_get_array(diag_physics,'i_rainc',i_rainc)
 call mpas_pool_get_array(diag_physics,'cuprec' ,cuprec )
 call mpas_pool_get_array(diag_physics,'rainc'  ,rainc  )

!update the accumulated precipitation at the end of each dynamic time-step:
 do i = its, ite
    rainc(i) = rainc(i) + cuprec(i) * dt_dyn

    if(l_acrain .and. bucket_rainc.gt.0._RKIND .and. &
       rainc(i).gt.bucket_rainc) then
       i_rainc(i) = i_rainc(i) + 1
       rainc(i)   = rainc(i) - bucket_rainc
    endif
 enddo

 end subroutine update_convection_step2

!=================================================================================================================
 end module mpas_atmphys_driver_convection
!=================================================================================================================
