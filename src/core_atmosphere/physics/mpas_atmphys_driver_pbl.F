! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_driver_pbl
 use mpas_kind_types
 use mpas_pool_routines
 use mpas_timer, only : mpas_timer_start, mpas_timer_stop

 use mpas_atmphys_constants
 use mpas_atmphys_vars

!wrf physics:
 use module_bl_mynn_wrapper
 use module_bl_ysu

 implicit none
 private
 public:: allocate_pbl,   &
          deallocate_pbl, &
          driver_pbl

!MPAS driver for parameterization of Planetary Boundary Layer (PBL) processes.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines in mpas_atmphys_driver_pbl:
! ---------------------------------------
! allocate_pbl  : allocate local arrays for parameterization of PBL processes.
! deallocate_pbl: deallocate local arrays for parameterization of PBL processes.
! driver_pbl    : main driver (called from subroutine physics_driver).
! pbl_from_MPAS : initialize local arrays.
! pbl_to_MPAS   : copy local arrays to MPAS arrays.
!
! WRF physics called from driver_pbl:
! -----------------------------------
! * module_bl_ysu : YSU PBL scheme.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * removed the pre-processor option "do_hydrostatic_pressure" before call to subroutine ysu.
!   Laura D. Fowler (birch.ucar.edu) / 2013-05-29.
! * in call to subroutine ysu, replaced the variable g (that originally pointed to gravity)
!   with gravity, for simplicity.
!   Laura D. Fowler (laura@ucar.edu) / 2014-03-21.
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * renamed "ysu" with "bl_ysu".
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-25.
! * added the implementation of the MYNN PBL scheme from WRF 3.6.1.
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-30.
! * corrected the initialization of sh3d for the mynn parameterization.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-13.
! * for the mynn parameterization, change the definition of dx_p to match that used in other physics
!   parameterizations.
!   Laura D. Fowler (laura@ucar.edu) / 2016-10-18.
! * updated the call to subroutine ysu in conjunction with updating module_bl_ysu.F from WRF version 3.6.1 to
!   WRF version 3.8.1
!   Laura D. Fowler (laura@ucar.edu) / 2016-10-27.
! * since we removed the local variable pbl_scheme from mpas_atmphys_vars.F, now defines pbl_scheme as a pointer
!   to config_pbl_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2917-02-16.
! * after updating module_bl_ysu.F to WRF version 4.0.3, corrected call to subroutine ysu to output diagnostics of
!   exchange coefficients exch_h and exch_m.
!   Laura D. Fowler (laura@ucar.edu) / 2019-03-12.
! * updated the call to subroutine ysu after updating the YSU PBL scheme to that in WRF 4.4.1. added the flags
!   errmsg and errflg in the call to subroutine ysu for compliance with the CCPP framework. also removed local
!   variable regime_p which is no longer needed in the call to subroutine ysu.
!   Laura D. Fowler (laura@ucar.edu) / 2023-05-15.


 contains


!=================================================================================================================
 subroutine allocate_pbl(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

!local pointers:
 character(len=StrKIND),pointer:: pbl_scheme
 integer,pointer :: config_mynn_output
 integer,pointer :: config_tkebudget

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_pbl_scheme',pbl_scheme)

 if(.not.allocated(hfx_p)  ) allocate(hfx_p(ims:ime,jms:jme)  )
 if(.not.allocated(qfx_p)  ) allocate(qfx_p(ims:ime,jms:jme)  )
 if(.not.allocated(ust_p)  ) allocate(ust_p(ims:ime,jms:jme)  )
 if(.not.allocated(wspd_p) ) allocate(wspd_p(ims:ime,jms:jme) )
 if(.not.allocated(xland_p)) allocate(xland_p(ims:ime,jms:jme))
 if(.not.allocated(hpbl_p) ) allocate(hpbl_p(ims:ime,jms:jme) )
 if(.not.allocated(kpbl_p) ) allocate(kpbl_p(ims:ime,jms:jme) )
 if(.not.allocated(znt_p)  ) allocate(znt_p(ims:ime,jms:jme)  )
 if(.not.allocated(delta_p)) allocate(delta_p(ims:ime,jms:jme))
 if(.not.allocated(wstar_p)) allocate(wstar_p(ims:ime,jms:jme))
 if(.not.allocated(uoce_p) ) allocate(uoce_p(ims:ime,jms:jme) )
 if(.not.allocated(voce_p) ) allocate(voce_p(ims:ime,jms:jme) )


 !tendencies:
 if(.not.allocated(rublten_p) ) allocate(rublten_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(rvblten_p) ) allocate(rvblten_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(rthblten_p)) allocate(rthblten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rqvblten_p)) allocate(rqvblten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rqcblten_p)) allocate(rqcblten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rqiblten_p)) allocate(rqiblten_p(ims:ime,kms:kme,jms:jme))

 !exchange coefficients:
 if(.not.allocated(kzh_p)) allocate(kzh_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(kzm_p)) allocate(kzm_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(kzq_p)) allocate(kzq_p(ims:ime,kms:kme,jms:jme))

 pbl_select: select case (trim(pbl_scheme))

    case("bl_ysu")
       !from surface-layer model:
       if(.not.allocated(br_p)    ) allocate(br_p(ims:ime,jms:jme)          )
       if(.not.allocated(ctopo_p) ) allocate(ctopo_p(ims:ime,jms:jme)       )
       if(.not.allocated(ctopo2_p)) allocate(ctopo2_p(ims:ime,jms:jme)      )
       if(.not.allocated(psih_p)  ) allocate(psih_p(ims:ime,jms:jme)        )
       if(.not.allocated(psim_p)  ) allocate(psim_p(ims:ime,jms:jme)        )
       if(.not.allocated(u10_p)   ) allocate(u10_p(ims:ime,jms:jme)         )
       if(.not.allocated(v10_p)   ) allocate(v10_p(ims:ime,jms:jme)         )
       if(.not.allocated(exch_p)  ) allocate(exch_p(ims:ime,kms:kme,jms:jme))
       !from radiation schemes:
       if(.not.allocated(rthraten_p)) allocate(rthraten_p(ims:ime,kms:kme,jms:jme))

    case("bl_mynn")
       call mpas_pool_get_config(configs,'config_mynn_output',config_mynn_output)
       call mpas_pool_get_config(configs,'config_tkebudget',config_tkebudget)
         
       if(.not.allocated(dx_p)        ) allocate(dx_p(ims:ime,jms:jme)   )
       if(.not.allocated(ch_p)        ) allocate(ch_p(ims:ime,jms:jme)   )
       if(.not.allocated(qcg_p)       ) allocate(qcg_p(ims:ime,jms:jme)  )
       if(.not.allocated(qsfc_p)      ) allocate(qsfc_p(ims:ime,jms:jme) )
       if(.not.allocated(rmol_p)      ) allocate(rmol_p(ims:ime,jms:jme) )
       if(.not.allocated(tsk_p)       ) allocate(tsk_p(ims:ime,jms:jme)  )
       if(.not.allocated(vdfg_p)      ) allocate(vdfg_p(ims:ime,jms:jme) )
       if(.not.allocated(maxmf_p)     ) allocate(maxmf_p(ims:ime,jms:jme))
       if(.not.allocated(maxwidth_p)  ) allocate(maxwidth_p(ims:ime,jms:jme))
       if(.not.allocated(ztop_plume_p)) allocate(ztop_plume_p(ims:ime,jms:jme))

       if(.not.allocated(cov_p)   ) allocate(cov_p(ims:ime,kms:kme,jms:jme)   )
       if(.not.allocated(qke_p)   ) allocate(qke_p(ims:ime,kms:kme,jms:jme)   )
       if(.not.allocated(qsq_p)   ) allocate(qsq_p(ims:ime,kms:kme,jms:jme)   )
       if(.not.allocated(tsq_p)   ) allocate(tsq_p(ims:ime,kms:kme,jms:jme)   )
       if(.not.allocated(qkeadv_p)) allocate(qkeadv_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(elpbl_p) ) allocate(elpbl_p(ims:ime,kms:kme,jms:jme) )
       if(.not.allocated(sh3d_p)  ) allocate(sh3d_p(ims:ime,kms:kme,jms:jme)  )
       if(.not.allocated(sm3d_p)  ) allocate(sm3d_p(ims:ime,kms:kme,jms:jme)  )

       if (config_tkebudget .gt. 0) then
          if(.not.allocated(dqke_p)  ) allocate(dqke_p(ims:ime,kms:kme,jms:jme)  )
          if(.not.allocated(qbuoy_p) ) allocate(qbuoy_p(ims:ime,kms:kme,jms:jme) )
          if(.not.allocated(qdiss_p) ) allocate(qdiss_p(ims:ime,kms:kme,jms:jme) )
          if(.not.allocated(qshear_p)) allocate(qshear_p(ims:ime,kms:kme,jms:jme))
          if(.not.allocated(qwt_p)   ) allocate(qwt_p(ims:ime,kms:kme,jms:jme)   )
       else
          if(.not.allocated(dqke_p)  ) allocate(dqke_p(1,1,1)  )
          if(.not.allocated(qbuoy_p) ) allocate(qbuoy_p(1,1,1) )
          if(.not.allocated(qdiss_p) ) allocate(qdiss_p(1,1,1) )
          if(.not.allocated(qshear_p)) allocate(qshear_p(1,1,1))
          if(.not.allocated(qwt_p)   ) allocate(qwt_p(1,1,1)   )
       endif

       if (config_mynn_output .gt. 0) then
          if(.not.allocated(edmf_a_p)  ) allocate(edmf_a_p(ims:ime,kms:kme,jms:jme)  )
          if(.not.allocated(edmf_w_p)  ) allocate(edmf_w_p(ims:ime,kms:kme,jms:jme)  )
          if(.not.allocated(edmf_thl_p)) allocate(edmf_thl_p(ims:ime,kms:kme,jms:jme))
          if(.not.allocated(edmf_qt_p) ) allocate(edmf_qt_p(ims:ime,kms:kme,jms:jme) )
          if(.not.allocated(edmf_qc_p) ) allocate(edmf_qc_p(ims:ime,kms:kme,jms:jme) )
          if(.not.allocated(edmf_ent_p)) allocate(edmf_ent_p(ims:ime,kms:kme,jms:jme))
          if(.not.allocated(sub_thl_p) ) allocate(edmf_w_p(ims:ime,kms:kme,jms:jme)  )
          if(.not.allocated(sub_sqv_p) ) allocate(edmf_thl_p(ims:ime,kms:kme,jms:jme))
          if(.not.allocated(det_thl_p) ) allocate(edmf_qt_p(ims:ime,kms:kme,jms:jme) )
          if(.not.allocated(det_sqv_p) ) allocate(edmf_qc_p(ims:ime,kms:kme,jms:jme) )
       else
          if(.not.allocated(edmf_a_p)  ) allocate(edmf_a_p(1,1,1)  )
          if(.not.allocated(edmf_w_p)  ) allocate(edmf_w_p(1,1,1)  )
          if(.not.allocated(edmf_thl_p)) allocate(edmf_thl_p(1,1,1))
          if(.not.allocated(edmf_qt_p) ) allocate(edmf_qt_p(1,1,1) )
          if(.not.allocated(edmf_qc_p) ) allocate(edmf_qc_p(1,1,1) )
          if(.not.allocated(edmf_ent_p)) allocate(edmf_ent_p(1,1,1))
          if(.not.allocated(sub_thl_p) ) allocate(sub_thl_p(1,1,1) )
          if(.not.allocated(sub_sqv_p) ) allocate(sub_sqv_p(1,1,1) )
          if(.not.allocated(det_thl_p) ) allocate(det_thl_p(1,1,1) )
          if(.not.allocated(det_sqv_p) ) allocate(det_sqv_p(1,1,1) )
       endif
       
       if(.not.allocated(rqsblten_p)  ) allocate(rqsblten_p(ims:ime,kms:kme,jms:jme)  )
       if(.not.allocated(rniblten_p)  ) allocate(rniblten_p(ims:ime,kms:kme,jms:jme)  )
       if(.not.allocated(rncblten_p)  ) allocate(rncblten_p(ims:ime,kms:kme,jms:jme)  )
       if(.not.allocated(rnwfablten_p)) allocate(rnwfablten_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(rnifablten_p)) allocate(rnifablten_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(rnbcablten_p)) allocate(rnbcablten_p(ims:ime,kms:kme,jms:jme)) 
       if(.not.allocated(nbca_p)      ) allocate(nbca_p(ims:ime,kms:kme,jms:jme)      )
       
       if(.not.allocated(qc_bl_p)    ) allocate(qc_bl_p(ims:ime,kms:kme,jms:jme)    )
       if(.not.allocated(qi_bl_p)    ) allocate(qi_bl_p(ims:ime,kms:kme,jms:jme)    )
       if(.not.allocated(cldfrac_bl_p))allocate(cldfrac_bl_p(ims:ime,kms:kme,jms:jme))

       !from radiation schemes:
       if(.not.allocated(rthraten_p)) allocate(rthraten_p(ims:ime,kms:kme,jms:jme))

       !spp:
       if(.not.allocated(pattern_spp_pbl_p)) allocate(pattern_spp_pbl_p(ims:ime,kms:kme,jms:jme))
    case default

 end select pbl_select

 end subroutine allocate_pbl

!=================================================================================================================
 subroutine deallocate_pbl(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

!local pointers:
 character(len=StrKIND),pointer:: pbl_scheme

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_pbl_scheme',pbl_scheme)

 if(allocated(hfx_p)  ) deallocate(hfx_p  )
 if(allocated(qfx_p)  ) deallocate(qfx_p  )
 if(allocated(ust_p)  ) deallocate(ust_p  )
 if(allocated(wspd_p) ) deallocate(wspd_p )
 if(allocated(xland_p)) deallocate(xland_p)
 if(allocated(hpbl_p) ) deallocate(hpbl_p )
 if(allocated(kpbl_p) ) deallocate(kpbl_p )
 if(allocated(znt_p)  ) deallocate(znt_p  )
 if(allocated(delta_p)) deallocate(delta_p)
 if(allocated(wstar_p)) deallocate(wstar_p)
 if(allocated(uoce_p) ) deallocate(uoce_p )
 if(allocated(voce_p) ) deallocate(voce_p )

 !tendencies:
 if(allocated(rublten_p) ) deallocate(rublten_p )
 if(allocated(rvblten_p) ) deallocate(rvblten_p )
 if(allocated(rthblten_p)) deallocate(rthblten_p)
 if(allocated(rqvblten_p)) deallocate(rqvblten_p)
 if(allocated(rqcblten_p)) deallocate(rqcblten_p)
 if(allocated(rqiblten_p)) deallocate(rqiblten_p)

 !exchange coefficients:
 if(allocated(kzh_p)) deallocate(kzh_p)
 if(allocated(kzm_p)) deallocate(kzm_p)
 if(allocated(kzq_p)) deallocate(kzq_p)

 pbl_select: select case (trim(pbl_scheme))

    case("bl_ysu")
       !from surface-layer model:
       if(allocated(br_p)    ) deallocate(br_p    )
       if(allocated(ctopo_p) ) deallocate(ctopo_p )
       if(allocated(ctopo2_p)) deallocate(ctopo2_p)
       if(allocated(psih_p)  ) deallocate(psih_p  )
       if(allocated(psim_p)  ) deallocate(psim_p  )
       if(allocated(u10_p)   ) deallocate(u10_p   )
       if(allocated(v10_p)   ) deallocate(v10_p   )
       if(allocated(exch_p)  ) deallocate(exch_p  )
       !from radiation schemes:
       if(allocated(rthraten_p)) deallocate(rthraten_p)

    case("bl_mynn")
       if(allocated(dx_p)        ) deallocate(dx_p   )
       if(allocated(ch_p)        ) deallocate(ch_p   )
       if(allocated(qcg_p)       ) deallocate(qcg_p  )
       if(allocated(qsfc_p)      ) deallocate(qsfc_p )
       if(allocated(rmol_p)      ) deallocate(rmol_p )
       if(allocated(tsk_p)       ) deallocate(tsk_p  )
       if(allocated(vdfg_p)      ) deallocate(vdfg_p )
       if(allocated(maxmf_p)     ) deallocate(maxmf_p)
       if(allocated(maxwidth_p)  ) deallocate(maxwidth_p)
       if(allocated(ztop_plume_p)) deallocate(ztop_plume_p)

       if(allocated(cov_p)   ) deallocate(cov_p   )
       if(allocated(qke_p)   ) deallocate(qke_p   )
       if(allocated(qsq_p)   ) deallocate(qsq_p   )
       if(allocated(tsq_p)   ) deallocate(tsq_p   )
       if(allocated(qkeadv_p)) deallocate(qkeadv_p)
       if(allocated(elpbl_p) ) deallocate(elpbl_p )
       if(allocated(sh3d_p)  ) deallocate(sh3d_p  )
       if(allocated(sm3d_p)  ) deallocate(sm3d_p  )
       if(allocated(dqke_p)  ) deallocate(dqke_p  )
       if(allocated(qbuoy_p) ) deallocate(qbuoy_p )
       if(allocated(qdiss_p) ) deallocate(qdiss_p )
       if(allocated(qshear_p)) deallocate(qshear_p)
       if(allocated(qwt_p)   ) deallocate(qwt_p   )

       if(allocated(rqsblten_p)  ) deallocate(rqsblten_p  )
       if(allocated(rniblten_p)  ) deallocate(rniblten_p  )
       if(allocated(rncblten_p)  ) deallocate(rncblten_p  )
       if(allocated(rnwfablten_p)) deallocate(rnwfablten_p)
       if(allocated(rnifablten_p)) deallocate(rnifablten_p)
       if(allocated(rnbcablten_p)) deallocate(rnbcablten_p)
       if(allocated(nbca_p)      ) deallocate(nbca_p      )
       
       if(allocated(qc_bl_p)    ) deallocate(qc_bl_p    )
       if(allocated(qi_bl_p)    ) deallocate(qi_bl_p    )
       if(allocated(cldfrac_bl_p))deallocate(cldfrac_bl_p)
       if(allocated(rthraten_p) ) deallocate(rthraten_p )

       if(allocated(pattern_spp_pbl_p)) deallocate(pattern_spp_pbl_p) 
    case default

 end select pbl_select

 end subroutine deallocate_pbl

!=================================================================================================================
 subroutine pbl_from_MPAS(configs,mesh,sfc_input,diag_physics,tend_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: diag_physics
 type(mpas_pool_type),intent(in):: sfc_input
 type(mpas_pool_type),intent(in):: tend_physics

 integer,intent(in):: its,ite

!local variables:
 integer:: i,k,j

!local pointers:
 character(len=StrKIND),pointer:: pbl_scheme
 integer,pointer:: config_tkebudget
 integer,pointer:: config_mynn_mixlength
 integer,pointer:: config_mynn_cloudpdf
 integer,pointer:: config_mynn_edmf
 integer,pointer:: config_mynn_edmf_mom
 integer,pointer:: config_mynn_edmf_tke
 integer,pointer:: config_mynn_mixscalars
 integer,pointer:: config_mynn_mixqt
 integer,pointer:: config_mynn_output
 integer,pointer:: config_icloud_bl
 real,   pointer:: config_mynn_closure
 logical,pointer:: config_mynn_tkeadvect
      
 real(kind=RKIND),dimension(:),pointer:: hfx,hpbl,qfx,ust,wspd,xland,znt
 real(kind=RKIND),dimension(:),pointer:: delta,wstar

!local pointers for YSU scheme:
 logical,pointer:: config_ysu_pblmix
 real(kind=RKIND),dimension(:),pointer:: br,fh,fm,u10,v10
 real(kind=RKIND),dimension(:,:),pointer:: rthratenlw,rthratensw

!local pointers for MYNN scheme:
 real(kind=RKIND),pointer:: len_disp
 real(kind=RKIND),dimension(:),pointer  :: meshDensity
 real(kind=RKIND),dimension(:),pointer  :: ch,qsfc,qcg,rmol,skintemp
 real(kind=RKIND),dimension(:,:),pointer:: cov,qsq,tsq,sh3d,sm3d
 real(kind=RKIND),dimension(:,:),pointer:: qke,qke_adv,el_pbl
 real(kind=RKIND),dimension(:,:),pointer:: cldfrac_bl,qi_bl,qc_bl
 real(kind=RKIND),dimension(:,:),pointer:: dqke,qdiss,qshear,qbuoy,qwt
 real(kind=RKIND),dimension(:,:),pointer:: edmf_a,edmf_w,edmf_thl
 real(kind=RKIND),dimension(:,:),pointer:: edmf_qt,edmf_qc,edmf_ent
 real(kind=RKIND),dimension(:,:),pointer:: sub_thl,sub_sqv,det_thl,det_sqv
      
!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_pbl_scheme',pbl_scheme)

 call mpas_pool_get_array(diag_physics,'hfx'   ,hfx   )
 call mpas_pool_get_array(diag_physics,'hpbl'  ,hpbl  )
 call mpas_pool_get_array(diag_physics,'qfx'   ,qfx   )
 call mpas_pool_get_array(diag_physics,'ust'   ,ust   )
 call mpas_pool_get_array(diag_physics,'wspd'  ,wspd  )
 call mpas_pool_get_array(diag_physics,'znt'   ,znt   )
 call mpas_pool_get_array(diag_physics,'wstar' ,wstar )
 call mpas_pool_get_array(sfc_input,'skintemp' ,skintemp)
 call mpas_pool_get_array(sfc_input   ,'xland' ,xland )

 do j = jts,jte
 do i = its,ite
    !from surface-layer model:
    hfx_p(i,j)    = hfx(i)
    hpbl_p(i,j)   = hpbl(i)
    qfx_p(i,j)    = qfx(i)
    ust_p(i,j)    = ust(i)
    wspd_p(i,j)   = wspd(i)
    xland_p(i,j)  = xland(i)
    tsk_p(i,j)    = skintemp(i)
    kpbl_p(i,j)   = 1
    znt_p(i,j)    = znt(i)
    wstar_p(i,j)  = wstar(i)
    !... ocean currents are set to zero:
    uoce_p(i,j)  = 0._RKIND
    voce_p(i,j)  = 0._RKIND
 enddo
 enddo

 pbl_select: select case (trim(pbl_scheme))

    case("bl_ysu")
       call mpas_pool_get_config(configs,'config_ysu_pblmix',config_ysu_pblmix)

       call mpas_pool_get_array(diag_physics,'br'    ,br    )
       call mpas_pool_get_array(diag_physics,'fm'    ,fm    )
       call mpas_pool_get_array(diag_physics,'fh'    ,fh    )
       call mpas_pool_get_array(diag_physics,'u10'   ,u10   )
       call mpas_pool_get_array(diag_physics,'v10'   ,v10   )
       call mpas_pool_get_array(diag_physics,'delta' ,delta )
       
       call mpas_pool_get_array(tend_physics,'rthratenlw',rthratenlw)
       call mpas_pool_get_array(tend_physics,'rthratensw',rthratensw)

       ysu_pblmix = 0
       if(config_ysu_pblmix) ysu_pblmix = 1

       do j = jts,jte
       do i = its,ite
          !from surface-layer model:
          br_p(i,j)     = br(i)
          psim_p(i,j)   = fm(i)
          psih_p(i,j)   = fh(i)
          u10_p(i,j)    = u10(i)
          v10_p(i,j)    = v10(i)
          delta_p(i,j)  = delta(i)
          !initialization for YSU PBL scheme:
          ctopo_p(i,j)  = 1._RKIND
          ctopo2_p(i,j) = 1._RKIND
       enddo
       enddo

       do j = jts,jte
       do k = kts,kte
       do i = its,ite
          exch_p(i,k,j)     = 0._RKIND
          rthraten_p(i,k,j) = rthratenlw(k,i) + rthratensw(k,i) 
       enddo
       enddo
       enddo

    case("bl_mynn")
       call mpas_pool_get_config(configs,'config_len_disp'       ,len_disp)
       call mpas_pool_get_config(configs,'config_mynn_output'    ,config_mynn_output)
       call mpas_pool_get_config(configs,'config_tkebudget'      ,config_tkebudget)
       call mpas_pool_get_config(configs,'config_mynn_cloudpdf'  ,config_mynn_cloudpdf)
       call mpas_pool_get_config(configs,'config_mynn_mixlength' ,config_mynn_mixlength)
       call mpas_pool_get_config(configs,'config_mynn_edmf'      ,config_mynn_edmf)
       call mpas_pool_get_config(configs,'config_mynn_edmf_mom'  ,config_mynn_edmf_mom)
       call mpas_pool_get_config(configs,'config_mynn_edmf_tke'  ,config_mynn_edmf_tke)
       call mpas_pool_get_config(configs,'config_mynn_mixscalars',config_mynn_mixscalars)
       call mpas_pool_get_config(configs,'config_mynn_mixqt'     ,config_mynn_mixqt)
       call mpas_pool_get_config(configs,'config_mynn_closure'   ,config_mynn_closure)
       call mpas_pool_get_config(configs,'config_mynn_tkeadvect' ,config_mynn_tkeadvect)
       call mpas_pool_get_config(configs,'config_icloud_bl'      ,config_icloud_bl) 
       
       call mpas_pool_get_array(mesh,'meshDensity',meshDensity)

       call mpas_pool_get_array(sfc_input   ,'skintemp',skintemp)
       call mpas_pool_get_array(diag_physics,'ch'      ,ch      )
       call mpas_pool_get_array(diag_physics,'qcg'     ,qcg     )
       call mpas_pool_get_array(diag_physics,'qsfc'    ,qsfc    )
       call mpas_pool_get_array(diag_physics,'rmol'    ,rmol    )

       call mpas_pool_get_array(diag_physics,'el_pbl' ,el_pbl   )
       call mpas_pool_get_array(diag_physics,'cov'    ,cov      )
       call mpas_pool_get_array(diag_physics,'qke'    ,qke      )
       call mpas_pool_get_array(diag_physics,'qke_adv',qke_adv  )
       call mpas_pool_get_array(diag_physics,'qsq'    ,qsq      )
       call mpas_pool_get_array(diag_physics,'tsq'    ,tsq      )
       call mpas_pool_get_array(diag_physics,'sh3d'   ,sh3d     )
       call mpas_pool_get_array(diag_physics,'sm3d'   ,sm3d     )
       call mpas_pool_get_array(diag_physics,'cldfrac_bl',cldfrac_bl)
       call mpas_pool_get_array(diag_physics,'qc_bl'  ,qc_bl    )
       call mpas_pool_get_array(diag_physics,'qi_bl'  ,qi_bl    )
       
       if (config_tkebudget .gt. 0) then
          call mpas_pool_get_array(diag_physics,'dqke'  ,dqke   )
          call mpas_pool_get_array(diag_physics,'qdiss' ,qdiss  )
          call mpas_pool_get_array(diag_physics,'qwt'   ,qwt    )
          call mpas_pool_get_array(diag_physics,'qbuoy' ,qbuoy  )
          call mpas_pool_get_array(diag_physics,'qshear',qshear )
       endif

       if (config_mynn_output .gt. 0) then
          call mpas_pool_get_array(diag_physics,'edmf_a'  ,edmf_a  )
          call mpas_pool_get_array(diag_physics,'edmf_w'  ,edmf_w  )
          call mpas_pool_get_array(diag_physics,'edmf_thl',edmf_thl)
          call mpas_pool_get_array(diag_physics,'edmf_qt' ,edmf_qt )
          call mpas_pool_get_array(diag_physics,'edmf_qc' ,edmf_qc )
          call mpas_pool_get_array(diag_physics,'edmf_ent',edmf_ent)
          call mpas_pool_get_array(diag_physics,'sub_thl' ,sub_thl )
          call mpas_pool_get_array(diag_physics,'sub_sqv' ,sub_sqv )
          call mpas_pool_get_array(diag_physics,'det_thl' ,det_thl )
          call mpas_pool_get_array(diag_physics,'det_sqv' ,det_sqv )
       endif

       do j = jts,jte
       do i = its,ite
          dx_p(i,j)    = len_disp / meshDensity(i)**0.25
          ch_p(i,j)    = ch(i)
          qcg_p(i,j)   = qcg(i)
          qsfc_p(i,j)  = qsfc(i)
          rmol_p(i,j)  = rmol(i)
          tsk_p(i,j)   = skintemp(i)
          !... no gravitational settling of fog/cloud droplets (grav_settling = 0):
          vdfg_p(i,j)  = 0._RKIND
       enddo
       enddo

       do j = jts,jte
       do k = kts,kte
       do i = its,ite
          elpbl_p(i,k,j)      = el_pbl(k,i)
          cov_p(i,k,j)        = cov(k,i)
          qke_p(i,k,j)        = qke(k,i)
          qsq_p(i,k,j)        = qsq(k,i)
          tsq_p(i,k,j)        = tsq(k,i)
          qkeadv_p(i,k,j)     = qke_adv(k,i)
          sh3d_p(i,k,j)       = sh3d(k,i)
          cldfrac_bl_p(i,k,j) = cldfrac_bl(k,i)
          qi_bl_p(i,k,j)      = qi_bl(k,i)
          qc_bl_p(i,k,j)      = qc_bl(k,i)

          rqsblten_p(i,k,j)   = 0._RKIND
          rniblten_p(i,k,j)   = 0._RKIND
          rncblten_p(i,k,j)   = 0.0_RKIND
          rnwfablten_p(i,k,j) = 0.0_RKIND
          rnifablten_p(i,k,j) = 0.0_RKIND
          rnbcablten_p(i,k,j) = 0.0_RKIND
          nbca_p(i,k,j)       = 0.0_RKIND
          
          pattern_spp_pbl_p(i,k,j)=0.0_RKIND
       enddo
       enddo
       enddo

       !... outputs:
       dqke_p        = 0._RKIND
       qbuoy_p       = 0._RKIND
       qdiss_p       = 0._RKIND
       qshear_p      = 0._RKIND
       qwt_p         = 0._RKIND

       edmf_a_p      = 0._RKIND
       edmf_w_p      = 0._RKIND
       edmf_thl_p    = 0._RKIND
       edmf_qt_p     = 0._RKIND
       edmf_qc_p     = 0._RKIND
       edmf_ent_p    = 0._RKIND
       sub_thl_p     = 0._RKIND
       sub_sqv_p     = 0._RKIND
       det_thl_p     = 0._RKIND
       det_sqv_p     = 0._RKIND
       
    case default

 end select pbl_select

 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    rublten_p(i,k,j)  = 0._RKIND
    rvblten_p(i,k,j)  = 0._RKIND
    rthblten_p(i,k,j) = 0._RKIND
    rqvblten_p(i,k,j) = 0._RKIND
    rqcblten_p(i,k,j) = 0._RKIND
    rqiblten_p(i,k,j) = 0._RKIND

    kzh_p(i,k,j)    = 0._RKIND
    kzm_p(i,k,j)    = 0._RKIND
    kzq_p(i,k,j)    = 0._RKIND
 enddo
 enddo
 enddo

 end subroutine pbl_from_MPAS

!=================================================================================================================
 subroutine pbl_to_MPAS(configs,diag_physics,tend_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

 integer,intent(in):: its,ite

!local variables:
 integer:: i,k,j

!local pointers:
 character(len=StrKIND),pointer:: pbl_scheme
 integer,pointer::config_tkebudget
 integer,pointer::config_mynn_output

 integer,dimension(:),pointer:: kpbl

 real(kind=RKIND),dimension(:),pointer  :: hpbl
 real(kind=RKIND),dimension(:,:),pointer:: kzh,kzm,kzq
 real(kind=RKIND),dimension(:,:),pointer:: rublten,rvblten,rthblten
 real(kind=RKIND),dimension(:,:),pointer:: rqvblten,rqcblten,rqiblten
 real(kind=RKIND),dimension(:,:),pointer:: rniblten,rncblten
 real(kind=RKIND),dimension(:,:),pointer:: rnwfablten,rnifablten

!local pointers for YSU scheme:
 real(kind=RKIND),dimension(:,:),pointer:: exch_h

!local pointers for MYNN scheme:
 real(kind=RKIND),dimension(:),pointer  :: delta,wstar
 real(kind=RKIND),dimension(:,:),pointer:: cov,qke,qsq,tsq,sh3d,sm3d
 real(kind=RKIND),dimension(:,:),pointer:: qke_adv,el_pbl
 real(kind=RKIND),dimension(:,:),pointer:: cldfrac_bl,qi_bl,qc_bl
 real(kind=RKIND),dimension(:,:),pointer:: dqke,qdiss,qshear,qbuoy,qwt
 real(kind=RKIND),dimension(:,:),pointer:: edmf_a,edmf_w,edmf_thl
 real(kind=RKIND),dimension(:,:),pointer:: edmf_qt,edmf_qc,edmf_ent
 real(kind=RKIND),dimension(:,:),pointer:: sub_thl,sub_sqv,det_thl,det_sqv
!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_pbl_scheme'  ,pbl_scheme)
 call mpas_pool_get_config(configs,'config_tkebudget'   ,config_tkebudget)
 call mpas_pool_get_config(configs,'config_mynn_output' ,config_mynn_output)

 call mpas_pool_get_array(diag_physics,'kpbl' ,kpbl )
 call mpas_pool_get_array(diag_physics,'hpbl' ,hpbl )
 call mpas_pool_get_array(diag_physics,'kzh'  ,kzh  )
 call mpas_pool_get_array(diag_physics,'kzm'  ,kzm  )
 call mpas_pool_get_array(diag_physics,'kzq'  ,kzq  )
 call mpas_pool_get_array(diag_physics,'wstar',wstar)

 call mpas_pool_get_array(tend_physics,'rublten' ,rublten )
 call mpas_pool_get_array(tend_physics,'rvblten' ,rvblten )
 call mpas_pool_get_array(tend_physics,'rthblten',rthblten)
 call mpas_pool_get_array(tend_physics,'rqvblten',rqvblten)
 call mpas_pool_get_array(tend_physics,'rqcblten',rqcblten)
 call mpas_pool_get_array(tend_physics,'rqiblten',rqiblten)

 do j = jts,jte
 do i = its,ite
    hpbl(i)  = hpbl_p(i,j)
    kpbl(i)  = kpbl_p(i,j)
    wstar(i) = wstar_p(i,j)
 enddo
 enddo

 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    rublten(k,i)  = rublten_p(i,k,j)
    rvblten(k,i)  = rvblten_p(i,k,j)
    rthblten(k,i) = rthblten_p(i,k,j)
    rqvblten(k,i) = rqvblten_p(i,k,j)
    rqcblten(k,i) = rqcblten_p(i,k,j)
    rqiblten(k,i) = rqiblten_p(i,k,j)

    kzh(k,i) = kzh_p(i,k,j)
    kzm(k,i) = kzm_p(i,k,j)
    kzq(k,i) = kzh_p(i,k,j)
 enddo
 enddo
 enddo

 pbl_select: select case (trim(pbl_scheme))

    case("bl_ysu")
       call mpas_pool_get_array(diag_physics,'exch_h',exch_h)
       call mpas_pool_get_array(diag_physics,'delta',delta)
       
       do j = jts,jte
       do k = kts,kte
       do i = its,ite
          exch_h(k,i) = exch_p(i,k,j)
       enddo
       enddo
       enddo

       do j = jts,jte
       do i = its,ite
          delta(i) = delta_p(i,j)
       enddo
       enddo
       
    case("bl_mynn")
       call mpas_pool_get_array(diag_physics,'el_pbl'  , el_pbl   )
       call mpas_pool_get_array(diag_physics,'cov'     , cov      )
       call mpas_pool_get_array(diag_physics,'qke'     , qke      )
       call mpas_pool_get_array(diag_physics,'qke_adv' , qke_adv  )
       call mpas_pool_get_array(diag_physics,'qsq'     , qsq      )
       call mpas_pool_get_array(diag_physics,'tsq'     , tsq      )
       call mpas_pool_get_array(diag_physics,'sh3d'    , sh3d     )
       call mpas_pool_get_array(diag_physics,'sm3d'    , sm3d     )
       if (config_tkebudget .gt. 0) then
          call mpas_pool_get_array(diag_physics,'dqke'    , dqke     )
          call mpas_pool_get_array(diag_physics,'qbuoy'   , qbuoy    )
          call mpas_pool_get_array(diag_physics,'qdiss'   , qdiss    )
          call mpas_pool_get_array(diag_physics,'qshear'  , qshear   )
          call mpas_pool_get_array(diag_physics,'qwt'     , qwt      )
       endif
       if (config_mynn_output .gt. 0) then
          call mpas_pool_get_array(diag_physics,'edmf_a'   , edmf_a   )
          call mpas_pool_get_array(diag_physics,'edmf_w'   , edmf_w   )
          call mpas_pool_get_array(diag_physics,'edmf_thl' , edmf_thl )
          call mpas_pool_get_array(diag_physics,'edmf_qt'  , edmf_qt  )
          call mpas_pool_get_array(diag_physics,'edmf_qc'  , edmf_qc  )
          call mpas_pool_get_array(diag_physics,'edmf_ent' , edmf_ent )
          call mpas_pool_get_array(diag_physics,'sub_thl'  , sub_thl  )
          call mpas_pool_get_array(diag_physics,'sub_sqv'  , sub_sqv  )
          call mpas_pool_get_array(diag_physics,'det_thl'  , det_thl  )
          call mpas_pool_get_array(diag_physics,'det_sqv'  , det_sqv  )
       endif
       call mpas_pool_get_array(diag_physics,'cldfrac_bl',cldfrac_bl)
       call mpas_pool_get_array(diag_physics,'qc_bl'    ,qc_bl    )
       call mpas_pool_get_array(diag_physics,'qi_bl'    ,qi_bl    )
       call mpas_pool_get_array(tend_physics,'rniblten', rniblten )
       call mpas_pool_get_array(tend_physics,'rncblten', rncblten )
       call mpas_pool_get_array(tend_physics,'rnwfablten', rnwfablten )
       call mpas_pool_get_array(tend_physics,'rnifablten', rnifablten )
       
       do j = jts,jte
       do k = kts,kte
       do i = its,ite
          el_pbl(k,i)   = elpbl_p(i,k,j)
          cov(k,i)      = cov_p(i,k,j)
          qke(k,i)      = qke_p(i,k,j)
          qsq(k,i)      = qsq_p(i,k,j)
          tsq(k,i)      = tsq_p(i,k,j)
          sh3d(k,i)     = sh3d_p(i,k,j)
          sm3d(k,i)     = sm3d_p(i,k,j)
          qke_adv(k,i)  = qkeadv_p(i,k,j)

          cldfrac_bl(k,i)= cldfrac_bl_p(i,k,j)
          qc_bl(k,i)    = qc_bl_p(i,k,j)
          qi_bl(k,i)    = qi_bl_p(i,k,j)

          rncblten(k,i) = rncblten_p(i,k,j)
          rniblten(k,i) = rniblten_p(i,k,j)
          rnwfablten(k,i) = rnwfablten_p(i,k,j)
          rnifablten(k,i) = rnifablten_p(i,k,j)
       enddo
       enddo
       enddo

       !optional output
       if (config_tkebudget .gt. 0) then
          do j = jts,jte
          do k = kts,kte
          do i = its,ite
             dqke(k,i)     = dqke_p(i,k,j)
             qbuoy(k,i)    = qbuoy_p(i,k,j)
             qdiss(k,i)    = qdiss_p(i,k,j)
             qshear(k,i)   = qshear_p(i,k,j)
             qwt(k,i)      = qwt_p(i,k,j)
          enddo
          enddo
          enddo
       endif

       if (config_mynn_output .gt. 0) then
          do j = jts,jte
          do k = kts,kte
          do i = its,ite
             edmf_a(k,i)     = edmf_a_p(i,k,j)
             edmf_w(k,i)     = edmf_w_p(i,k,j)
             edmf_thl(k,i)   = edmf_thl_p(i,k,j)
             edmf_qt(k,i)    = edmf_qt_p(i,k,j)
             edmf_qc(k,i)    = edmf_qc_p(i,k,j)
             edmf_ent(k,i)   = edmf_ent_p(i,k,j)
             det_thl(k,i)    = sub_thl_p(i,k,j)
             det_sqv(k,i)    = sub_sqv_p(i,k,j)
             det_thl(k,i)    = det_thl_p(i,k,j)
             det_sqv(k,i)    = det_sqv_p(i,k,j)
          enddo
          enddo
          enddo
       endif
      
    case default

 end select pbl_select

 end subroutine pbl_to_MPAS
 
!=================================================================================================================
 subroutine driver_pbl(itimestep,configs,mesh,sfc_input,diag_physics,tend_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh

 integer,intent(in):: its,ite
 integer,intent(in):: itimestep

!inout arguments:
 type(mpas_pool_type),intent(inout):: sfc_input
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

!local pointers:
 logical,pointer:: config_do_restart
 character(len=StrKIND),pointer:: pbl_scheme
 integer,pointer:: config_icloud_bl
 integer,pointer:: config_tkebudget
 integer,pointer:: config_mynn_mixlength
 integer,pointer:: config_mynn_cloudpdf
 integer,pointer:: config_mynn_edmf
 integer,pointer:: config_mynn_edmf_mom
 integer,pointer:: config_mynn_edmf_tke
 integer,pointer:: config_mynn_mixscalars
 integer,pointer:: config_mynn_mixqt
 integer,pointer:: config_mynn_cloudmix
 integer,pointer:: config_mynn_output
 integer,pointer:: config_spp_pbl
 real   ,pointer:: config_mynn_closure
 logical,pointer:: config_mynn_tkeadvect
 logical,pointer:: config_do_DAcycling
      
!local variables:
 integer:: initflag
 integer:: i,k,j

!CCPP-compliant flags:
 character(len=StrKIND):: errmsg
 integer:: errflg
 integer, pointer :: indexToCellID(:)
!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('')
!call mpas_log_write('--- enter subroutine driver_pbl:')

!initialization of CCPP-compliant flags:
 errmsg = ' '
 errflg = 0

 call mpas_pool_get_config(configs,'config_do_restart'     ,config_do_restart     )
 call mpas_pool_get_config(configs,'config_do_DAcycling'   ,config_do_DAcycling   )
 call mpas_pool_get_config(configs,'config_pbl_scheme'     ,pbl_scheme            )
 call mpas_pool_get_config(configs,'config_icloud_bl'      ,config_icloud_bl      )
 call mpas_pool_get_config(configs,'config_mynn_output'    ,config_mynn_output    )
 call mpas_pool_get_config(configs,'config_mynn_cloudpdf'  ,config_mynn_cloudpdf  )
 call mpas_pool_get_config(configs,'config_mynn_mixlength' ,config_mynn_mixlength )
 call mpas_pool_get_config(configs,'config_mynn_edmf'      ,config_mynn_edmf      )
 call mpas_pool_get_config(configs,'config_mynn_edmf_mom'  ,config_mynn_edmf_mom  )
 call mpas_pool_get_config(configs,'config_mynn_edmf_tke'  ,config_mynn_edmf_tke  )
 call mpas_pool_get_config(configs,'config_mynn_mixscalars',config_mynn_mixscalars)
 call mpas_pool_get_config(configs,'config_mynn_mixqt'     ,config_mynn_mixqt     )
 call mpas_pool_get_config(configs,'config_mynn_cloudmix'  ,config_mynn_cloudmix  ) 
 call mpas_pool_get_config(configs,'config_mynn_closure'   ,config_mynn_closure   )
 call mpas_pool_get_config(configs,'config_mynn_tkeadvect' ,config_mynn_tkeadvect )
 call mpas_pool_get_config(configs,'config_tkebudget'      ,config_tkebudget )
 call mpas_pool_get_config(configs,'config_spp_pbl'        ,config_spp_pbl )
 call mpas_pool_get_array(mesh,'indexToCellID',indexToCellID)
!copy MPAS arrays to local arrays:
 call pbl_from_MPAS(configs,mesh,sfc_input,diag_physics,tend_physics,its,ite)

 initflag = 1
 if(config_do_restart .or. itimestep > 1) initflag = 0

 pbl_select: select case (trim(pbl_scheme))

    case("bl_ysu")
       call mpas_timer_start('bl_ysu')
       call ysu ( &
                 p3d      = pres_hyd_p , p3di     = pres2_hyd_p , psfc     = psfc_p     , &
                 t3d      = t_p        , dz8w     = dz_p        , pi3d     = pi_p       , &
                 u3d      = u_p        , v3d      = v_p         , qv3d     = qv_p       , &
                 qc3d     = qc_p       , qi3d     = qi_p        , rublten  = rublten_p  , &
                 rvblten  = rvblten_p  , rthblten = rthblten_p  , rqvblten = rqvblten_p , &
                 rqcblten = rqcblten_p , rqiblten = rqiblten_p  , flag_qc  = f_qc       , &
                 flag_qi  = f_qi       , cp       = cp          , g        = gravity    , &
                 rovcp    = rcp        , rd       = R_d         , rovg     = rdg        , & 
                 ep1      = ep_1       , ep2      = ep_2        , karman   = karman     , &
                 xlv      = xlv        , rv       = R_v         , znt      = znt_p      , &
                 ust      = ust_p      , hpbl     = hpbl_p      , psim     = psim_p     , &
                 psih     = psih_p     , xland    = xland_p     , hfx      = hfx_p      , &
                 qfx      = qfx_p      , wspd     = wspd_p      , br       = br_p       , &
                 dt       = dt_pbl     , kpbl2d   = kpbl_p      , exch_h   = kzh_p      , &
                 exch_m   = kzm_p      , wstar    = wstar_p     , delta    = delta_p    , &
                 uoce     = uoce_p     , voce     = voce_p      , rthraten = rthraten_p , &
                 u10      = u10_p      , v10      = v10_p       , ctopo    = ctopo_p    , &
                 ctopo2   = ctopo2_p   , flag_bep = flag_bep    , idiff    = idiff      , &
                 ysu_topdown_pblmix = ysu_pblmix ,                                        &
                 errmsg   = errmsg     , errflg   = errflg      ,                         &
                 ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde  , &
                 ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme  , &
                 its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte    &
                )
       call mpas_timer_stop('bl_ysu')

    case("bl_mynn")
       call mpas_timer_start('bl_mynnedmf')
       call mynnedmf_wrapper_run ( &
                 !3d input/output
                 p        = pres_hyd_p    , exner    = pi_p          , t3d      = t_p          , &
                 th       = th_p          , dz       = dz_p          , rho      = rho_p        , &
                 u        = u_p           , v        = v_p           , w        = w_p          , &
                 qv       = qv_p          , qc       = qc_p          , qi       = qi_p         , &
                 qs       = qs_p          , qni      = ni_p          , qnc      = nc_p         , &
                 qnwfa    = nwfa_p        , qnifa    = nifa_p        , qnbca    = nbca_p       , &
                 qc_bl    = qc_bl_p       , qi_bl    = qi_bl_p       , cldfra_bl= cldfrac_bl_p , &
                 pattern_spp_pbl= pattern_spp_pbl_p,                                             &
                 !tendency arrays
                 rublten  = rublten_p     , rvblten  = rvblten_p     , rthblten = rthblten_p   , &
                 rqvblten = rqvblten_p    , rqcblten = rqcblten_p    , rqiblten = rqiblten_p   , &
                 rqsblten = rqsblten_p    , rqncblten = rncblten_p   , rqniblten = rniblten_p  , &
                 rqnwfablten=rnwfablten_p , rqnifablten=rnifablten_p , rqnbcablten=rnbcablten_p, &
                 rthraten = rthraten_p    ,                                                      &
                 !moist flags
                 flag_qc  = f_qc          , flag_qnc = f_qnc         , flag_qs  = f_qs         , &
                 flag_qi  = f_qi          , flag_qni = f_qni         , flag_qnwfa=f_nwfa       , &
                 flag_qnifa=f_nifa        , flag_qnbca=.false.       ,                           &
                 !2d arrays
                 kpbl     = kpbl_p        , pblh     = hpbl_p        , xland    = xland_p      , &
                 hfx      = hfx_p         , qfx      = qfx_p         , ch       = ch_p         , &
                 qsfc     = qsfc_p        , ust      = ust_p         , rmol     = rmol_p       , &
                 wspd     = wspd_p        , dxc      = dx_p          , delt     = dt_pbl       , &
                 uoce     = uoce_p        , voce     = voce_p        , ts       = tsk_p        , &
                 ps       = psfc_p        , znt      = znt_p         , maxmf    = maxmf_p      , &
                 maxwidth = maxwidth_p    , ztop_plume=ztop_plume_p  ,                           &
                 !turbulence arrays
                 qke      = qke_p         , qke_adv  = qkeadv_p      , el_pbl   = elpbl_p      , &
                 sh3d     = sh3d_p        , sm3d     = sm3d_p        , cov      = cov_p        , &
                 tsq      = tsq_p         , qsq      = qsq_p         , kh       = kzh_p        , &
                 km       = kzm_p         ,                                                      &
                 !tke budget arrays
                 dqke     = dqke_p        , qwt      = qwt_p         ,                           &
                 qshear   = qshear_p      , qbuoy    = qbuoy_p       , qdiss    = qdiss_p      , &
                 !optional output
      	      	 edmf_a   = edmf_a_p      , edmf_w   = edmf_w_p      , edmf_thl = edmf_thl_p   , &
                 edmf_qc  = edmf_qc_p     , edmf_qt  = edmf_qt_p     , edmf_ent = edmf_ent_p   , &
                 sub_thl  = sub_thl_p     , sub_sqv  = sub_sqv_p     , det_thl  = det_thl_p    , &
                 det_sqv  = det_sqv_p     ,                                                      &
#if (WRF_CHEM == 1)
                 !smoke/chem---not yet added
                 mix_chem = mix_chem      , chem3d   = chem_p        , vd3d     = vd_p         , &
                 nchem    = nchem         , kdvel    = kdvel         , ndvel    = ndvel        , &
                 num_vert_mix=numvertmix  , frp_mean = frp_mean_p    , emis_ant_no=emis_ant_no_p,&
                 enh_mix  = enh_mix       ,                                                      &
#endif
                 !configuration parameters
                 initflag                 = initflag                 ,                           &
                 restart                  = config_do_restart        ,                           &
                 cycling                  = config_do_DAcycling      ,                           &
                 bl_mynn_cloudpdf         = config_mynn_cloudpdf     ,                           &
                 bl_mynn_tkeadvect        = config_mynn_tkeadvect    ,                           &
                 tke_budget               = config_tkebudget         ,                           &
                 bl_mynn_output           = config_mynn_output       ,                           &
                 bl_mynn_mixlength        = config_mynn_mixlength    ,                           &
                 bl_mynn_edmf             = config_mynn_edmf         ,                           &
                 bl_mynn_edmf_mom         = config_mynn_edmf_mom     ,                           &
                 bl_mynn_edmf_tke         = config_mynn_edmf_tke     ,                           &
                 bl_mynn_mixscalars       = config_mynn_mixscalars   ,                           &
                 bl_mynn_mixqt            = config_mynn_mixqt        ,                           &
                 bl_mynn_cloudmix         = config_mynn_cloudmix     ,                           &
                 bl_mynn_closure          = config_mynn_closure      ,                           &
                 icloud_bl                = config_icloud_bl         ,                           &
                 spp_pbl                  = config_spp_pbl           ,                           &
                 ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde         , &
                 ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme         , &
                 its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte         ) 
       call mpas_timer_stop('bl_mynnedmf')

    case default

 end select pbl_select

!copy local arrays to MPAS grid:
 call pbl_to_MPAS(configs,diag_physics,tend_physics,its,ite)

!call mpas_log_write('--- end subroutine driver_pbl.')

 end subroutine driver_pbl

!=================================================================================================================
 end module mpas_atmphys_driver_pbl
!=================================================================================================================
