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
 use mpas_timer,only: mpas_timer_start,mpas_timer_stop
 use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

 use mpas_atmphys_constants
 use mpas_atmphys_vars
 use mpas_satmedmfvdifq_wrapper_mod, only : mpas_satmedmfvdifq_config_type, &
                                            mpas_call_satmedmfvdifq

 use bl_mynn,only: bl_mynn_init
 use module_bl_mynn,only: mynn_bl_driver
 use module_bl_mynnedmf_driver,only: mynnedmf_driver
 use module_bl_ysu
 use module_bl_myjpbl

 implicit none
 private
 public:: allocate_pbl,   &
          deallocate_pbl, &
          init_pbl,       &
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
! * in the call to subroutine mynn_bl_driver,renamed f_qnc to f_nc, and f_qni to f_ni.
!   Laura D. Fowler (laura@ucar.edu) / 2024-02-14.
! * updated the MYNN PBL scheme to the sourcecode from WRF version 4.6.
!   Laura D. Fowler (laura@ucar.edu) / 2024-02.15.


 contains


!=================================================================================================================
 subroutine allocate_pbl(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

!local pointers:
 character(len=StrKIND),pointer:: pbl_scheme

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
 if(.not.allocated(uoce_p) ) allocate(uoce_p(ims:ime,jms:jme) )
 if(.not.allocated(voce_p) ) allocate(voce_p(ims:ime,jms:jme) )

 !tendencies:
 if(.not.allocated(rublten_p) ) allocate(rublten_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(rvblten_p) ) allocate(rvblten_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(rthblten_p)) allocate(rthblten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rqvblten_p)) allocate(rqvblten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rqcblten_p)) allocate(rqcblten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rqiblten_p)) allocate(rqiblten_p(ims:ime,kms:kme,jms:jme))

 if(.not.allocated(rthraten_p)) allocate(rthraten_p(ims:ime,kms:kme,jms:jme))

 !exchange coefficients:
 if(.not.allocated(kzh_p)) allocate(kzh_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(kzm_p)) allocate(kzm_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(kzq_p)) allocate(kzq_p(ims:ime,kms:kme,jms:jme))

 !allocation of chemistry arrays
 if (do_chemistry) then
    if(.not.allocated(ddvel_p)) allocate(ddvel_p(ims:ime,jms:jme,1:num_chem))
    if(.not.allocated(tend_chem_settle_p)) allocate(tend_chem_settle_p(ims:ime,kms:kme,jms:jme,1:num_chem))
    if(.not.allocated(chem_p)) allocate(chem_p(ims:ime,kms:kme,jms:jme,1:num_chem))
    if(.not.allocated(frp_out_p)) allocate(frp_out_p(ims:ime,jms:jme))
    if(.not.allocated(aero_emis_for_enhmix_p)) allocate(aero_emis_for_enhmix_p(ims:ime,jms:jme))
 endif

 pbl_select: select case (trim(pbl_scheme))

    case("bl_ysu")
       !from surface-layer model:
       if(.not.allocated(br_p)    ) allocate(br_p(ims:ime,jms:jme)          )
       if(.not.allocated(ctopo_p) ) allocate(ctopo_p(ims:ime,jms:jme)       )
       if(.not.allocated(ctopo2_p)) allocate(ctopo2_p(ims:ime,jms:jme)      )
       if(.not.allocated(delta_p) ) allocate(delta_p(ims:ime,jms:jme)       )
       if(.not.allocated(psih_p)  ) allocate(psih_p(ims:ime,jms:jme)        )
       if(.not.allocated(psim_p)  ) allocate(psim_p(ims:ime,jms:jme)        )
       if(.not.allocated(u10_p)   ) allocate(u10_p(ims:ime,jms:jme)         )
       if(.not.allocated(v10_p)   ) allocate(v10_p(ims:ime,jms:jme)         )
       if(.not.allocated(exch_p)  ) allocate(exch_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(wstar_p) ) allocate(wstar_p(ims:ime,jms:jme)       )


    case("bl_tke_edmf")
       ! Minimal MPAS application additions for UFS TKE-EDMF adapter.
       ! Reuse existing local arrays to stage MPAS fields for the wrapper:
       !   dx_p       : areaCell
       !   tsk_p      : skin temperature
       !   psim_p     : fm
       !   psih_p     : fh
       !   rthraten_p : SW radiative heating
       !   exch_p     : LW radiative heating
       if(.not.allocated(vegfra_p)  ) allocate(vegfra_p(ims:ime,jms:jme)      )
       if(.not.allocated(br_p)      ) allocate(br_p(ims:ime,jms:jme)          )
       if(.not.allocated(dx_p)      ) allocate(dx_p(ims:ime,jms:jme)          )
       if(.not.allocated(tsk_p)     ) allocate(tsk_p(ims:ime,jms:jme)         )
       if(.not.allocated(psim_p)    ) allocate(psim_p(ims:ime,jms:jme)        )
       if(.not.allocated(psih_p)    ) allocate(psih_p(ims:ime,jms:jme)        )
       if(.not.allocated(u10_p)     ) allocate(u10_p(ims:ime,jms:jme)         )
       if(.not.allocated(v10_p)     ) allocate(v10_p(ims:ime,jms:jme)         )
       if(.not.allocated(exch_p)    ) allocate(exch_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(kzh_p)     ) allocate(kzh_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(kzm_p)     ) allocate(kzm_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(kzq_p)     ) allocate(kzq_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(rthraten_p)) allocate(rthraten_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(tkepbl_p)  ) allocate(tkepbl_p(ims:ime,kms:kme,jms:jme)    )

    case("bl_mynn")
       if(.not.allocated(kbl_plume_p) ) allocate(kbl_plume_p(ims:ime,jms:jme)         )
       if(.not.allocated(dx_p)        ) allocate(dx_p(ims:ime,jms:jme)                )
       if(.not.allocated(ch_p)        ) allocate(ch_p(ims:ime,jms:jme)                )
       if(.not.allocated(qsfc_p)      ) allocate(qsfc_p(ims:ime,jms:jme)              )
       if(.not.allocated(rmol_p)      ) allocate(rmol_p(ims:ime,jms:jme)              )
       if(.not.allocated(tsk_p)       ) allocate(tsk_p(ims:ime,jms:jme)               )
       if(.not.allocated(maxwidthbl_p)) allocate(maxwidthbl_p(ims:ime,jms:jme)        )
       if(.not.allocated(maxmfbl_p)   ) allocate(maxmfbl_p(ims:ime,jms:jme)           )
       if(.not.allocated(zbl_plume_p) ) allocate(zbl_plume_p(ims:ime,jms:jme)         )
       if(.not.allocated(cov_p)       ) allocate(cov_p(ims:ime,kms:kme,jms:jme)       )
       if(.not.allocated(qke_p)       ) allocate(qke_p(ims:ime,kms:kme,jms:jme)       )
       if(.not.allocated(qsq_p)       ) allocate(qsq_p(ims:ime,kms:kme,jms:jme)       )
       if(.not.allocated(tsq_p)       ) allocate(tsq_p(ims:ime,kms:kme,jms:jme)       )
       if(.not.allocated(qkeadv_p)    ) allocate(qkeadv_p(ims:ime,kms:kme,jms:jme)    )
       if(.not.allocated(elpbl_p)     ) allocate(elpbl_p(ims:ime,kms:kme,jms:jme)     )
       if(.not.allocated(tkepbl_p)    ) allocate(tkepbl_p(ims:ime,kms:kme,jms:jme)    )
       if(.not.allocated(sh3d_p)      ) allocate(sh3d_p(ims:ime,kms:kme,jms:jme)      )
       if(.not.allocated(sm3d_p)      ) allocate(sm3d_p(ims:ime,kms:kme,jms:jme)      )
       if(.not.allocated(dqke_p)      ) allocate(dqke_p(ims:ime,kms:kme,jms:jme)      )
       if(.not.allocated(qbuoy_p)     ) allocate(qbuoy_p(ims:ime,kms:kme,jms:jme)     )
       if(.not.allocated(qdiss_p)     ) allocate(qdiss_p(ims:ime,kms:kme,jms:jme)     )
       if(.not.allocated(qshear_p)    ) allocate(qshear_p(ims:ime,kms:kme,jms:jme)    )
       if(.not.allocated(qwt_p)       ) allocate(qwt_p(ims:ime,kms:kme,jms:jme)       )
       if(.not.allocated(qcbl_p)      ) allocate(qcbl_p(ims:ime,kms:kme,jms:jme)      )
       if(.not.allocated(qibl_p)      ) allocate(qibl_p(ims:ime,kms:kme,jms:jme)      )
       if(.not.allocated(cldfrabl_p)  ) allocate(cldfrabl_p(ims:ime,kms:kme,jms:jme)  )
       if(.not.allocated(edmfa_p)     ) allocate(edmfa_p(ims:ime,kms:kme,jms:jme)     )
       if(.not.allocated(edmfw_p)     ) allocate(edmfw_p(ims:ime,kms:kme,jms:jme)     )
       if(.not.allocated(edmfqt_p)    ) allocate(edmfqt_p(ims:ime,kms:kme,jms:jme)    )
       if(.not.allocated(edmfthl_p)   ) allocate(edmfthl_p(ims:ime,kms:kme,jms:jme)   )
       if(.not.allocated(edmfent_p)   ) allocate(edmfent_p(ims:ime,kms:kme,jms:jme)   )
       if(.not.allocated(edmfqc_p)    ) allocate(edmfqc_p(ims:ime,kms:kme,jms:jme)    )
       if(.not.allocated(subthl_p)    ) allocate(subthl_p(ims:ime,kms:kme,jms:jme)    )
       if(.not.allocated(subqv_p)     ) allocate(subqv_p(ims:ime,kms:kme,jms:jme)     )
       if(.not.allocated(detthl_p)    ) allocate(detthl_p(ims:ime,kms:kme,jms:jme)    )
       if(.not.allocated(detqv_p)     ) allocate(detqv_p(ims:ime,kms:kme,jms:jme)     )

       !additional tendencies:
       if(.not.allocated(rqsblten_p)  ) allocate(rqsblten_p(ims:ime,kms:kme,jms:jme)  )
       if(.not.allocated(rncblten_p)  ) allocate(rncblten_p(ims:ime,kms:kme,jms:jme)  )
       if(.not.allocated(rniblten_p)  ) allocate(rniblten_p(ims:ime,kms:kme,jms:jme)  )
       if(.not.allocated(rnifablten_p)) allocate(rnifablten_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(rnwfablten_p)) allocate(rnwfablten_p(ims:ime,kms:kme,jms:jme))

       !allocation of additional arrays:
       if(.not.allocated(pattern_spp_pbl)) allocate(pattern_spp_pbl(ims:ime,kms:kme,jms:jme))
       
    case("bl_mynnedmf") !GSL version
       if(.not.allocated(dx_p)          ) allocate(dx_p(ims:ime,jms:jme)                )
       if(.not.allocated(ch_p)          ) allocate(ch_p(ims:ime,jms:jme)                )
       if(.not.allocated(qsfc_p)        ) allocate(qsfc_p(ims:ime,jms:jme)              )
       if(.not.allocated(tsk_p)         ) allocate(tsk_p(ims:ime,jms:jme)               )
       if(.not.allocated(maxwidthbl_p)  ) allocate(maxwidthbl_p(ims:ime,jms:jme)        )
       if(.not.allocated(maxmfbl_p)     ) allocate(maxmfbl_p(ims:ime,jms:jme)           )
       if(.not.allocated(zbl_plume_p)   ) allocate(zbl_plume_p(ims:ime,jms:jme)         )
       if(.not.allocated(excessh_p)     ) allocate(excessh_p(ims:ime,jms:jme)           )
       if(.not.allocated(excessq_p)     ) allocate(excessq_p(ims:ime,jms:jme)           ) 
       
       if(.not.allocated(cov_p)         ) allocate(cov_p(ims:ime,kms:kme,jms:jme)       )
       if(.not.allocated(qke_p)         ) allocate(qke_p(ims:ime,kms:kme,jms:jme)       )
       if(.not.allocated(qsq_p)         ) allocate(qsq_p(ims:ime,kms:kme,jms:jme)       )
       if(.not.allocated(tsq_p)         ) allocate(tsq_p(ims:ime,kms:kme,jms:jme)       )
       if(.not.allocated(qkeadv_p)      ) allocate(qkeadv_p(ims:ime,kms:kme,jms:jme)    )
       if(.not.allocated(elpbl_p)       ) allocate(elpbl_p(ims:ime,kms:kme,jms:jme)     )
       if(.not.allocated(sh3d_p)        ) allocate(sh3d_p(ims:ime,kms:kme,jms:jme)      )
       if(.not.allocated(sm3d_p)        ) allocate(sm3d_p(ims:ime,kms:kme,jms:jme)      )
       if(.not.allocated(dqke_p)        ) allocate(dqke_p(ims:ime,kms:kme,jms:jme)      )
       if(.not.allocated(qbuoy_p)       ) allocate(qbuoy_p(ims:ime,kms:kme,jms:jme)     )
       if(.not.allocated(qdiss_p)       ) allocate(qdiss_p(ims:ime,kms:kme,jms:jme)     )
       if(.not.allocated(qshear_p)      ) allocate(qshear_p(ims:ime,kms:kme,jms:jme)    )
       if(.not.allocated(qwt_p)         ) allocate(qwt_p(ims:ime,kms:kme,jms:jme)       )
       if(.not.allocated(qcbl_p)        ) allocate(qcbl_p(ims:ime,kms:kme,jms:jme)      )
       if(.not.allocated(qibl_p)        ) allocate(qibl_p(ims:ime,kms:kme,jms:jme)      )
       if(.not.allocated(cldfrabl_p)    ) allocate(cldfrabl_p(ims:ime,kms:kme,jms:jme)  )
       if(.not.allocated(edmfa_p)       ) allocate(edmfa_p(ims:ime,kms:kme,jms:jme)     )
       if(.not.allocated(edmfw_p)       ) allocate(edmfw_p(ims:ime,kms:kme,jms:jme)     )
       if(.not.allocated(edmfqt_p)      ) allocate(edmfqt_p(ims:ime,kms:kme,jms:jme)    )
       if(.not.allocated(edmfthl_p)     ) allocate(edmfthl_p(ims:ime,kms:kme,jms:jme)   )
       if(.not.allocated(edmfent_p)     ) allocate(edmfent_p(ims:ime,kms:kme,jms:jme)   )
       if(.not.allocated(edmfqc_p)      ) allocate(edmfqc_p(ims:ime,kms:kme,jms:jme)    )
       if(.not.allocated(subthl_p)      ) allocate(subthl_p(ims:ime,kms:kme,jms:jme)    )
       if(.not.allocated(subqv_p)       ) allocate(subqv_p(ims:ime,kms:kme,jms:jme)     )
       if(.not.allocated(detthl_p)      ) allocate(detthl_p(ims:ime,kms:kme,jms:jme)    )
       if(.not.allocated(detqv_p)       ) allocate(detqv_p(ims:ime,kms:kme,jms:jme)     )

       !additional tendencies:
       if(.not.allocated(rqsblten_p)    ) allocate(rqsblten_p(ims:ime,kms:kme,jms:jme)  )
       if(.not.allocated(rncblten_p)    ) allocate(rncblten_p(ims:ime,kms:kme,jms:jme)  )
       if(.not.allocated(rniblten_p)    ) allocate(rniblten_p(ims:ime,kms:kme,jms:jme)  )
       if(.not.allocated(rnifablten_p)  ) allocate(rnifablten_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(rnwfablten_p)  ) allocate(rnwfablten_p(ims:ime,kms:kme,jms:jme))

       !allocation of additional arrays:
       if(.not.allocated(pattern_spp_pbl)) allocate(pattern_spp_pbl(ims:ime,kms:kme,jms:jme))

    case("bl_myj")
       if(.not.allocated(ht_p)         ) allocate(ht_p(ims:ime,jms:jme)      )
       if(.not.allocated(chlowq_p)     ) allocate(chlowq_p(ims:ime,jms:jme)  )
       if(.not.allocated(thz0_p)       ) allocate(thz0_p(ims:ime,jms:jme)    )
       if(.not.allocated(uz0_p)        ) allocate(uz0_p(ims:ime,jms:jme)     )
       if(.not.allocated(vz0_p)        ) allocate(vz0_p(ims:ime,jms:jme)     )
       if(.not.allocated(qz0_p)        ) allocate(qz0_p(ims:ime,jms:jme)     )
       if(.not.allocated(tsk_p)        ) allocate(tsk_p(ims:ime,jms:jme)     )
       if(.not.allocated(qsfc_p)       ) allocate(qsfc_p(ims:ime,jms:jme)    )
       if(.not.allocated(xice_p)       ) allocate(xice_p(ims:ime,jms:jme)    )
       if(.not.allocated(akhs_p)       ) allocate(akhs_p(ims:ime,jms:jme)    )
       if(.not.allocated(akms_p)       ) allocate(akms_p(ims:ime,jms:jme)    )
       if(.not.allocated(ct_p)         ) allocate(ct_p(ims:ime,jms:jme)      )
       if(.not.allocated(lowlyr_p)     ) allocate(lowlyr_p(ims:ime,jms:jme)  )
       if(.not.allocated(mixht_p)      ) allocate(mixht_p(ims:ime,jms:jme)   )
       if(.not.allocated(lh_p)         ) allocate(lh_p(ims:ime,jms:jme)      ) 
       if(.not.allocated(snow_p)       ) allocate(snow_p(ims:ime,jms:jme)   )
       if(.not.allocated(tkepbl_p)     ) allocate(tkepbl_p(ims:ime,kms:kme,jms:jme)  )
       if(.not.allocated(elpbl_p)      ) allocate(elpbl_p(ims:ime,kms:kme,jms:jme)   )
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
 if(allocated(uoce_p) ) deallocate(uoce_p )
 if(allocated(voce_p) ) deallocate(voce_p )

 !tendencies:
 if(allocated(rublten_p) ) deallocate(rublten_p )
 if(allocated(rvblten_p) ) deallocate(rvblten_p )
 if(allocated(rthblten_p)) deallocate(rthblten_p)
 if(allocated(rqvblten_p)) deallocate(rqvblten_p)
 if(allocated(rqcblten_p)) deallocate(rqcblten_p)
 if(allocated(rqiblten_p)) deallocate(rqiblten_p)

 if(allocated(rthraten_p)) deallocate(rthraten_p)

 !exchange coefficients:
 if(allocated(kzh_p)) deallocate(kzh_p)
 if(allocated(kzm_p)) deallocate(kzm_p)
 if(allocated(kzq_p)) deallocate(kzq_p)

 !deallocation of chemistry arrays
 if (do_chemistry) then
    if(allocated(ddvel_p)            ) deallocate(ddvel_p             )
    if(allocated(chem_p)             ) deallocate(chem_p              )
    if(allocated(tend_chem_settle_p) ) deallocate(tend_chem_settle_p  )
    if(allocated(frp_out_p)          ) deallocate(frp_out_p           )
    if(allocated(aero_emis_for_enhmix_p)) deallocate(aero_emis_for_enhmix_p )
 endif

 pbl_select: select case (trim(pbl_scheme))

    case("bl_ysu")
       !from surface-layer model:
       if(allocated(br_p)    ) deallocate(br_p    )
       if(allocated(ctopo_p) ) deallocate(ctopo_p )
       if(allocated(ctopo2_p)) deallocate(ctopo2_p)
       if(allocated(delta_p) ) deallocate(delta_p )
       if(allocated(psih_p)  ) deallocate(psih_p  )
       if(allocated(psim_p)  ) deallocate(psim_p  )
       if(allocated(u10_p)   ) deallocate(u10_p   )
       if(allocated(v10_p)   ) deallocate(v10_p   )
       if(allocated(exch_p)  ) deallocate(exch_p  )
       if(allocated(wstar_p) ) deallocate(wstar_p )


    case("bl_tke_edmf")
       if(allocated(vegfra_p)  ) deallocate(vegfra_p  )
       if(allocated(br_p)      ) deallocate(br_p      )
       if(allocated(dx_p)      ) deallocate(dx_p      )
       if(allocated(tsk_p)     ) deallocate(tsk_p     )
       if(allocated(psim_p)    ) deallocate(psim_p    )
       if(allocated(psih_p)    ) deallocate(psih_p    )
       if(allocated(u10_p)     ) deallocate(u10_p     )
       if(allocated(v10_p)     ) deallocate(v10_p     )
       if(allocated(exch_p)    ) deallocate(exch_p    )
       if(allocated(kzh_p)     ) deallocate(kzh_p     )
       if(allocated(kzm_p)     ) deallocate(kzm_p     )
       if(allocated(kzq_p)     ) deallocate(kzq_p     )
       if(allocated(rthraten_p)) deallocate(rthraten_p)
       if(allocated(tkepbl_p)  ) deallocate(tkepbl_p    )

    case("bl_mynn")
       if(allocated(kbl_plume_p) ) deallocate(kbl_plume_p )
       if(allocated(dx_p)        ) deallocate(dx_p        )
       if(allocated(ch_p)        ) deallocate(ch_p        )
       if(allocated(qsfc_p)      ) deallocate(qsfc_p      )
       if(allocated(rmol_p)      ) deallocate(rmol_p      )
       if(allocated(tsk_p)       ) deallocate(tsk_p       )
       if(allocated(maxwidthbl_p)) deallocate(maxwidthbl_p)
       if(allocated(maxmfbl_p)   ) deallocate(maxmfbl_p   )
       if(allocated(zbl_plume_p) ) deallocate(zbl_plume_p )

       if(allocated(cov_p)       ) deallocate(cov_p       )
       if(allocated(qke_p)       ) deallocate(qke_p       )
       if(allocated(qsq_p)       ) deallocate(qsq_p       )
       if(allocated(tsq_p)       ) deallocate(tsq_p       )
       if(allocated(qkeadv_p)    ) deallocate(qkeadv_p    )
       if(allocated(elpbl_p)     ) deallocate(elpbl_p     )
       if(allocated(tkepbl_p)    ) deallocate(tkepbl_p    )
       if(allocated(sh3d_p)      ) deallocate(sh3d_p      )
       if(allocated(sm3d_p)      ) deallocate(sm3d_p      )
       if(allocated(dqke_p)      ) deallocate(dqke_p      )
       if(allocated(qbuoy_p)     ) deallocate(qbuoy_p     )
       if(allocated(qdiss_p)     ) deallocate(qdiss_p     )
       if(allocated(qshear_p)    ) deallocate(qshear_p    )
       if(allocated(qwt_p)       ) deallocate(qwt_p       )
       if(allocated(qcbl_p)      ) deallocate(qcbl_p      )
       if(allocated(qibl_p)      ) deallocate(qibl_p      )
       if(allocated(cldfrabl_p)  ) deallocate(cldfrabl_p  )
       if(allocated(edmfa_p)     ) deallocate(edmfa_p     )
       if(allocated(edmfw_p)     ) deallocate(edmfw_p     )
       if(allocated(edmfqt_p)    ) deallocate(edmfqt_p    )
       if(allocated(edmfthl_p)   ) deallocate(edmfthl_p   )
       if(allocated(edmfent_p)   ) deallocate(edmfent_p   )
       if(allocated(edmfqc_p)    ) deallocate(edmfqc_p    )
       if(allocated(subthl_p)    ) deallocate(subthl_p    )
       if(allocated(subqv_p)     ) deallocate(subqv_p     )
       if(allocated(detthl_p)    ) deallocate(detthl_p    )
       if(allocated(detqv_p)     ) deallocate(detqv_p     )

       !additional tendencies:
       if(allocated(rqsblten_p)  ) deallocate(rqsblten_p  )
       if(allocated(rncblten_p)  ) deallocate(rncblten_p  )
       if(allocated(rniblten_p)  ) deallocate(rniblten_p  )
       if(allocated(rnifablten_p)) deallocate(rnifablten_p)
       if(allocated(rnwfablten_p)) deallocate(rnwfablten_p)

       !deallocation of additional arrays:
       if(allocated(pattern_spp_pbl)) deallocate(pattern_spp_pbl)
       
    case("bl_mynnedmf")
       if(allocated(dx_p)          ) deallocate(dx_p          )
       if(allocated(ch_p)          ) deallocate(ch_p          )
       if(allocated(qsfc_p)        ) deallocate(qsfc_p        )
       if(allocated(tsk_p)         ) deallocate(tsk_p         )
       if(allocated(maxwidthbl_p)  ) deallocate(maxwidthbl_p  )
       if(allocated(maxmfbl_p)     ) deallocate(maxmfbl_p     )
       if(allocated(zbl_plume_p)   ) deallocate(zbl_plume_p   )
       if(allocated(excessh_p)     ) deallocate(excessh_p     )                                                                                                         
       if(allocated(excessq_p)     ) deallocate(excessq_p     )
       
       if(allocated(cov_p)       ) deallocate(cov_p       )
       if(allocated(qke_p)       ) deallocate(qke_p       )
       if(allocated(qsq_p)       ) deallocate(qsq_p       )
       if(allocated(tsq_p)       ) deallocate(tsq_p       )
       if(allocated(qkeadv_p)    ) deallocate(qkeadv_p    )
       if(allocated(elpbl_p)     ) deallocate(elpbl_p     )
       if(allocated(sh3d_p)      ) deallocate(sh3d_p      )
       if(allocated(sm3d_p)      ) deallocate(sm3d_p      )
       if(allocated(dqke_p)      ) deallocate(dqke_p      )
       if(allocated(qbuoy_p)     ) deallocate(qbuoy_p     )
       if(allocated(qdiss_p)     ) deallocate(qdiss_p     )
       if(allocated(qshear_p)    ) deallocate(qshear_p    )
       if(allocated(qwt_p)       ) deallocate(qwt_p       )
       if(allocated(qcbl_p)      ) deallocate(qcbl_p      )
       if(allocated(qibl_p)      ) deallocate(qibl_p      )
       if(allocated(cldfrabl_p)  ) deallocate(cldfrabl_p  )
       if(allocated(edmfa_p)     ) deallocate(edmfa_p     )
       if(allocated(edmfw_p)     ) deallocate(edmfw_p     )
       if(allocated(edmfqt_p)    ) deallocate(edmfqt_p    )
       if(allocated(edmfthl_p)   ) deallocate(edmfthl_p   )
       if(allocated(edmfent_p)   ) deallocate(edmfent_p   )
       if(allocated(edmfqc_p)    ) deallocate(edmfqc_p    )
       if(allocated(subthl_p)    ) deallocate(subthl_p    )
       if(allocated(subqv_p)     ) deallocate(subqv_p     )
       if(allocated(detthl_p)    ) deallocate(detthl_p    )
       if(allocated(detqv_p)     ) deallocate(detqv_p     )

       !additional tendencies:
       if(allocated(rqsblten_p)  ) deallocate(rqsblten_p  )
       if(allocated(rncblten_p)  ) deallocate(rncblten_p  )
       if(allocated(rniblten_p)  ) deallocate(rniblten_p  )
       if(allocated(rnifablten_p)) deallocate(rnifablten_p)
       if(allocated(rnwfablten_p)) deallocate(rnwfablten_p)

       !deallocation of additional arrays:
       if(allocated(pattern_spp_pbl)) deallocate(pattern_spp_pbl)

   case("bl_myj")
       if(allocated(ht_p)         ) deallocate(ht_p      )
       if(allocated(chlowq_p)     ) deallocate(chlowq_p  )
       if(allocated(thz0_p)       ) deallocate(thz0_p    )
       if(allocated(uz0_p)        ) deallocate(uz0_p     )
       if(allocated(vz0_p)        ) deallocate(vz0_p     )
       if(allocated(qz0_p)        ) deallocate(qz0_p     )
       if(allocated(tsk_p)        ) deallocate(tsk_p     )
       if(allocated(xice_p)       ) deallocate(xice_p    )
       if(allocated(snow_p)       ) deallocate(snow_p    )
       if(allocated(tkepbl_p)     ) deallocate(tkepbl_p  )
       if(allocated(elpbl_p)      ) deallocate(elpbl_p   )
       if(allocated(akhs_p)       ) deallocate(akhs_p    )
       if(allocated(akms_p)       ) deallocate(akms_p    )
       if(allocated(ct_p)         ) deallocate(ct_p      )
       if(allocated(lowlyr_p)     ) deallocate(lowlyr_p  )
       if(allocated(mixht_p)      ) deallocate(mixht_p   )
       if(allocated(lh_p)         ) deallocate(lh_p      )
       if(allocated(qsfc_p)       ) deallocate(qsfc_p    )

    case default

 end select pbl_select

 end subroutine deallocate_pbl

!=================================================================================================================
 subroutine pbl_from_MPAS(configs,state,mesh,sfc_input,diag_physics,tend_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: diag_physics
 type(mpas_pool_type),intent(in):: sfc_input
 type(mpas_pool_type),intent(in):: tend_physics
 type(mpas_pool_type),intent(in):: state

 integer,intent(in):: its,ite

!local variables:
 integer:: i,k,j,n

!local pointers:
 character(len=StrKIND),pointer:: pbl_scheme

 real(kind=RKIND),dimension(:),pointer:: hfx,hpbl,qfx,ust,wspd,xland,znt
 real(kind=RKIND),dimension(:),pointer:: delta,wstar

!local pointers for YSU scheme:
 logical,pointer:: config_ysu_pblmix
 real(kind=RKIND),dimension(:),pointer:: br,fh,fm,u10,v10
 real(kind=RKIND),dimension(:,:),pointer:: rthratenlw,rthratensw

!local pointers for bl_tke_edmf adapter:
 real(kind=RKIND),dimension(:),pointer:: areaCell,skintemp,vegfra

!local pointers for MYNN scheme:
 real(kind=RKIND),pointer:: len_disp
 real(kind=RKIND),dimension(:),pointer  :: meshDensity
 real(kind=RKIND),dimension(:),pointer  :: ch,qsfc,rmol
 real(kind=RKIND),dimension(:,:),pointer:: cov,qke,qsq,tsq,sh3d,sm3d,tke_pbl,qke_adv,el_pbl
 real(kind=RKIND),dimension(:,:),pointer:: cldfrac_bl,qc_bl,qi_bl
 real(kind=RKIND),dimension(:,:),pointer:: edmf_a,edmf_ent,edmf_qc,edmf_qt,edmf_thl,edmf_w
 real(kind=RKIND),dimension(:,:),pointer:: sub_thl,sub_qv,det_thl,det_qv

!local pointers for tke_edmf
! tke is stored as scalars(5,k,i) for bl_tke_edmf_in.

!local pointers for chemistry
 integer,pointer                           :: chemistry_start,chemistry_end
 real(kind=RKIND),dimension(:,:,:),pointer :: scalars
 real(kind=RKIND),dimension(:,:,:),pointer :: tend_chem_settle
 real(kind=RKIND),dimension(:,:  ),pointer :: ddvel
 real(kind=RKIND),dimension(:,:,:),pointer :: chem
 real(kind=RKIND),dimension(:),pointer     :: frp_out
 real(kind=RKIND),dimension(:),pointer     :: aero_emis_for_enhmix
 integer,pointer:: drydep_opt
 logical,pointer:: enh_mix
 character(len=StrKIND),pointer :: config_smoke_scheme, &
                                   config_anthro_scheme, &
                                   config_rwc_scheme
      
!local pointers for MYJ scheme:
 real(kind=RKIND),dimension(:),pointer   :: chlowq,thz0,qz0,uz0,vz0,ct,akhs,akms,lh,mixht
 real(kind=RKIND),dimension(:,:),pointer :: zgrid
 integer,dimension(:),pointer            :: lowlyr
!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_pbl_scheme',pbl_scheme)

 call mpas_pool_get_array(diag_physics,'hfx' ,hfx )
 call mpas_pool_get_array(diag_physics,'hpbl',hpbl)
 call mpas_pool_get_array(diag_physics,'qfx' ,qfx )
 call mpas_pool_get_array(diag_physics,'ust' ,ust )
 call mpas_pool_get_array(diag_physics,'wspd',wspd)
 call mpas_pool_get_array(diag_physics,'znt' ,znt )

 call mpas_pool_get_array(tend_physics,'rthratenlw',rthratenlw)
 call mpas_pool_get_array(tend_physics,'rthratensw',rthratensw)

 call mpas_pool_get_array(sfc_input,'xland',xland)

 do j = jts,jte
 do i = its,ite
    !from surface-layer model:
    hfx_p(i,j)   = hfx(i)
    hpbl_p(i,j)  = hpbl(i)
    qfx_p(i,j)   = qfx(i)
    ust_p(i,j)   = ust(i)
    wspd_p(i,j)  = wspd(i)
    xland_p(i,j) = xland(i)
    kpbl_p(i,j)  = 1
    znt_p(i,j)   = znt(i)
    !... ocean currents are set to zero:
    uoce_p(i,j)  = 0._RKIND
    voce_p(i,j)  = 0._RKIND
 enddo
 do k = kts,kte
 do i = its,ite
    rthraten_p(i,k,j) = rthratenlw(k,i) + rthratensw(k,i)
 enddo
 enddo
 enddo

 if (do_chemistry) then
    call mpas_pool_get_array(state,'scalars',scalars,1)
    call mpas_pool_get_config(configs,'drydep_opt',drydep_opt)
    if ( drydep_opt .ne. 0 ) then
       call mpas_pool_get_array(tend_physics,'tend_chem_settle',tend_chem_settle)
       call mpas_pool_get_array(diag_physics,'ddvel',ddvel)
    endif
    call mpas_pool_get_dimension(state, 'chemistry_start', chemistry_start)
    call mpas_pool_get_dimension(state, 'chemistry_end', chemistry_end)
    chem => scalars(chemistry_start:chemistry_end,:,:)
    call mpas_pool_get_config(configs,'config_mynn_enh_mix',enh_mix)

    if (enh_mix) then
       call mpas_pool_get_config(configs,'config_smoke_scheme',config_smoke_scheme)
       if ( trim(config_smoke_scheme) .eq. 'on' ) then
          call mpas_pool_get_array(diag_physics,'frp_out',frp_out)
          do j = jts,jte
          do i = its,ite
                frp_out_p(i,j)     = frp_out(i)
          enddo
          enddo
       else
          do j = jts,jte
          do i = its,ite
                frp_out_p(i,j)     = 0._RKIND
          enddo
          enddo
       endif

       call mpas_pool_get_config(configs,'config_anthro_scheme',config_anthro_scheme)
       call mpas_pool_get_config(configs,'config_rwc_scheme',config_rwc_scheme)
       if ( config_anthro_scheme .ne. 'off' .or. config_rwc_scheme .ne. 'off') then
          call mpas_pool_get_array(diag_physics,'aero_emis_for_enhmix',aero_emis_for_enhmix)
          do j = jts,jte
          do i = its,ite
             aero_emis_for_enhmix_p(i,j) = aero_emis_for_enhmix(i)
          enddo
          enddo
       else
          do j = jts,jte
          do i = its,ite
             aero_emis_for_enhmix_p(i,j) = 0._RKIND
          enddo
          enddo
       endif
    endif ! If enh


    do n = 1,num_chem
    do j = jts,jte
    do i = its,ite
       if ( drydep_opt .ne. 0 ) then
          ddvel_p(i,j,n) = ddvel(n,i)
       else
          ddvel_p(i,j,:) = 0._RKIND
       endif
    enddo
    enddo
    enddo

    do n = 1,num_chem
    do j = jts,jte
    do k = kts,kte
    do i = its,ite
       chem_p(i,k,j,n) = chem(n,k,i)
       if ( drydep_opt .ne. 0 ) then
          tend_chem_settle_p(i,k,j,n) = tend_chem_settle(n,k,i)
       else
          tend_chem_settle_p(i,k,j,n) = 0._RKIND
       endif
    enddo
    enddo
    enddo
    enddo
 endif ! chem

 pbl_select: select case (trim(pbl_scheme))

    case("bl_ysu")
       call mpas_pool_get_config(configs,'config_ysu_pblmix',config_ysu_pblmix)

       call mpas_pool_get_array(diag_physics,'br'   ,br   )
       call mpas_pool_get_array(diag_physics,'delta',delta)
       call mpas_pool_get_array(diag_physics,'fm'   ,fm   )
       call mpas_pool_get_array(diag_physics,'fh'   ,fh   )
       call mpas_pool_get_array(diag_physics,'u10'  ,u10  )
       call mpas_pool_get_array(diag_physics,'v10'  ,v10  )
       call mpas_pool_get_array(diag_physics,'wstar',wstar)

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
          wstar_p(i,j)  = wstar(i)
          !initialization for YSU PBL scheme:
          ctopo_p(i,j)  = 1._RKIND
          ctopo2_p(i,j) = 1._RKIND
       enddo
       enddo

       do j = jts,jte
       do k = kts,kte
       do i = its,ite
          exch_p(i,k,j) = 0._RKIND
       enddo
       enddo
       enddo


    case("bl_tke_edmf")
       ! Stage MPAS fields needed by the UFS TKE-EDMF adapter.
       call mpas_pool_get_array(sfc_input,   'vegfra'    ,vegfra    )
       call mpas_pool_get_array(mesh        ,'areaCell'  ,areaCell  )
       call mpas_pool_get_array(sfc_input   ,'skintemp'  ,skintemp  )
       call mpas_pool_get_array(diag_physics,'br'        ,br        )
       call mpas_pool_get_array(diag_physics,'fm'        ,fm        )
       call mpas_pool_get_array(diag_physics,'fh'        ,fh        )
       call mpas_pool_get_array(diag_physics,'u10'       ,u10       )
       call mpas_pool_get_array(diag_physics,'v10'       ,v10       )
       call mpas_pool_get_array(tend_physics,'rthratenlw',rthratenlw)
       call mpas_pool_get_array(tend_physics,'rthratensw',rthratensw)
! ------------------------------------------------------------------
! Prognostic TKE for GFS TKE-EDMF
!
! In the Registry, tke is a member of var_array scalars:
!   qv, qc, qr, qi, tke, ...
! Therefore it cannot be retrieved with mpas_pool_get_array(scalars,'tke',...).
! Get the parent state scalar array and use scalar index 5 for tke.
! ------------------------------------------------------------------
       call mpas_pool_get_array(state,'scalars',scalars,1)

       do j = jts,jte
       do i = its,ite
          vegfra_p(i,j) = vegfra(i)
          dx_p(i,j)   = areaCell(i)
          tsk_p(i,j)  = skintemp(i)
          br_p(i,j)   = br(i)
          psim_p(i,j) = fm(i)
          psih_p(i,j) = fh(i)
          u10_p(i,j)  = u10(i)
          v10_p(i,j)  = v10(i)
       enddo
       do k = kts,kte
       do i = its,ite
          rthraten_p(i,k,j) = rthratensw(k,i)
          exch_p(i,k,j)     = rthratenlw(k,i)
! Do not read tke_pbl from MPAS here.  In some MPAS builds this
! diagnostic is not allocated before the TKE-EDMF call, and reading
! tke_pbl(k,i) causes the SIGSEGV seen near this block.  The wrapper
! will update tkepbl_p and pbl_to_MPAS will write it back afterward.
!         tkepbl_p(i,k,j)   = 0.0_RKIND
          tkepbl_p(i,k,j)   = max(1.0e-9_RKIND, scalars(5,k,i))
       enddo
       enddo
       enddo

    case("bl_mynn")
       call mpas_pool_get_config(configs,'config_len_disp',len_disp)
       call mpas_pool_get_array(mesh,'meshDensity',meshDensity)

       call mpas_pool_get_array(sfc_input,'skintemp',skintemp)
       call mpas_pool_get_array(diag_physics,'ch'        ,ch        )
       call mpas_pool_get_array(diag_physics,'qsfc'      ,qsfc      )
       call mpas_pool_get_array(diag_physics,'rmol'      ,rmol      )
       call mpas_pool_get_array(diag_physics,'el_pbl'    ,el_pbl    )
       call mpas_pool_get_array(diag_physics,'cov'       ,cov       )
       call mpas_pool_get_array(diag_physics,'qke'       ,qke       )
       call mpas_pool_get_array(diag_physics,'qke_adv'   ,qke_adv   )
       call mpas_pool_get_array(diag_physics,'qsq'       ,qsq       )
       call mpas_pool_get_array(diag_physics,'tsq'       ,tsq       )
       call mpas_pool_get_array(diag_physics,'tke_pbl'   ,tke_pbl   )
       call mpas_pool_get_array(diag_physics,'sh3d'      ,sh3d      )
       call mpas_pool_get_array(diag_physics,'sm3d'      ,sm3d      )
       call mpas_pool_get_array(diag_physics,'cldfrac_bl',cldfrac_bl)
       call mpas_pool_get_array(diag_physics,'qc_bl'     ,qc_bl     )
       call mpas_pool_get_array(diag_physics,'qi_bl'     ,qi_bl     )
       call mpas_pool_get_array(diag_physics,'edmf_a'    ,edmf_a    )
       call mpas_pool_get_array(diag_physics,'edmf_ent'  ,edmf_ent  )
       call mpas_pool_get_array(diag_physics,'edmf_qc'   ,edmf_qc   )
       call mpas_pool_get_array(diag_physics,'edmf_qt'   ,edmf_qt   )
       call mpas_pool_get_array(diag_physics,'edmf_thl'  ,edmf_thl  )
       call mpas_pool_get_array(diag_physics,'edmf_w'    ,edmf_w    )
       call mpas_pool_get_array(diag_physics,'sub_thl'   ,sub_thl   )
       call mpas_pool_get_array(diag_physics,'sub_qv'    ,sub_qv    )
       call mpas_pool_get_array(diag_physics,'det_thl'   ,det_thl   )
       call mpas_pool_get_array(diag_physics,'det_qv'    ,det_qv    )

       do j = jts,jte
       do i = its,ite
          dx_p(i,j)   = len_disp / meshDensity(i)**0.25
          ch_p(i,j)   = ch(i)
          qsfc_p(i,j) = qsfc(i)
          rmol_p(i,j) = rmol(i)
          tsk_p(i,j)  = skintemp(i)
       enddo
       enddo

       do j = jts,jte
       do k = kts,kte
       do i = its,ite
          elpbl_p(i,k,j)    = el_pbl(k,i)
          cov_p(i,k,j)      = cov(k,i)
          qke_p(i,k,j)      = qke(k,i)
          qsq_p(i,k,j)      = qsq(k,i)
          tsq_p(i,k,j)      = tsq(k,i)
          tkepbl_p(i,k,j)   = tke_pbl(k,i)
          qkeadv_p(i,k,j)   = qke_adv(k,i)
          sh3d_p(i,k,j)     = sh3d(k,i)
          sm3d_p(i,k,j)     = sm3d(k,i)
          cldfrabl_p(i,k,j) = cldfrac_bl(k,i)
          qcbl_p(i,k,j)     = qc_bl(k,i)
          qibl_p(i,k,j)     = qi_bl(k,i)
          edmfa_p(i,k,j)    = edmf_a(k,i)
          edmfent_p(i,k,j)  = edmf_ent(k,i)
          edmfqc_p(i,k,j)   = edmf_qc(k,i)
          edmfqt_p(i,k,j)   = edmf_qt(k,i)
          edmfthl_p(i,k,j)  = edmf_thl(k,i)
          edmfw_p(i,k,j)    = edmf_w(k,i)
          subthl_p(i,k,j)   = sub_thl(k,i)
          subqv_p(i,k,j)    = sub_qv(k,i)
          detthl_p(i,k,j)   = det_thl(k,i)
          detqv_p(i,k,j)    = det_qv(k,i)
          dqke_p(i,k,j)     = 0._RKIND
          qbuoy_p(i,k,j)    = 0._RKIND
          qdiss_p(i,k,j)    = 0._RKIND
          qshear_p(i,k,j)   = 0._RKIND
          qwt_p(i,k,j)      = 0._RKIND

          rqsblten_p(i,k,j)   = 0._RKIND
          rncblten_p(i,k,j)   = 0._RKIND
          rniblten_p(i,k,j)   = 0._RKIND
          rnifablten_p(i,k,j) = 0._RKIND
          rnwfablten_p(i,k,j) = 0._RKIND

          pattern_spp_pbl(i,k,j) = 0._RKIND
       enddo
       enddo
       do i = its,ite
          kbl_plume_p(i,j)  = 0
          maxwidthbl_p(i,j) = 0._RKIND
          maxmfbl_p(i,j)    = 0._RKIND
          zbl_plume_p(i,j)  = 0
       enddo
       enddo
       
    case("bl_mynnedmf")
       call mpas_pool_get_config(configs,'config_len_disp',len_disp)
       call mpas_pool_get_array(mesh,'meshDensity',meshDensity)

       call mpas_pool_get_array(sfc_input,'skintemp',skintemp)
       call mpas_pool_get_array(diag_physics,'ch'        ,ch        )
       call mpas_pool_get_array(diag_physics,'qsfc'      ,qsfc      )
       call mpas_pool_get_array(diag_physics,'el_pbl'    ,el_pbl    )
       call mpas_pool_get_array(diag_physics,'cov'       ,cov       )
       call mpas_pool_get_array(diag_physics,'qke'       ,qke       )
       call mpas_pool_get_array(diag_physics,'qke_adv'   ,qke_adv   )
       call mpas_pool_get_array(diag_physics,'qsq'       ,qsq       )
       call mpas_pool_get_array(diag_physics,'tsq'       ,tsq       )
       call mpas_pool_get_array(diag_physics,'sh3d'      ,sh3d      )
       call mpas_pool_get_array(diag_physics,'sm3d'      ,sm3d      )
       call mpas_pool_get_array(diag_physics,'cldfrac_bl',cldfrac_bl)
       call mpas_pool_get_array(diag_physics,'qc_bl'     ,qc_bl     )
       call mpas_pool_get_array(diag_physics,'qi_bl'     ,qi_bl     )
       call mpas_pool_get_array(diag_physics,'edmf_a'    ,edmf_a    )
       call mpas_pool_get_array(diag_physics,'edmf_ent'  ,edmf_ent  )
       call mpas_pool_get_array(diag_physics,'edmf_qc'   ,edmf_qc   )
       call mpas_pool_get_array(diag_physics,'edmf_qt'   ,edmf_qt   )
       call mpas_pool_get_array(diag_physics,'edmf_thl'  ,edmf_thl  )
       call mpas_pool_get_array(diag_physics,'edmf_w'    ,edmf_w    )
       call mpas_pool_get_array(diag_physics,'sub_thl'   ,sub_thl   )
       call mpas_pool_get_array(diag_physics,'sub_qv'    ,sub_qv    )
       call mpas_pool_get_array(diag_physics,'det_thl'   ,det_thl   )
       call mpas_pool_get_array(diag_physics,'det_qv'    ,det_qv    )

       do j = jts,jte
       do i = its,ite
          dx_p(i,j)   = len_disp / meshDensity(i)**0.25
          ch_p(i,j)   = ch(i)
          qsfc_p(i,j) = qsfc(i)
          tsk_p(i,j)  = skintemp(i)
       enddo
       enddo

       do j = jts,jte
       do k = kts,kte
       do i = its,ite
          elpbl_p(i,k,j)    = el_pbl(k,i)
          cov_p(i,k,j)      = cov(k,i)
          qke_p(i,k,j)      = qke(k,i)
          qsq_p(i,k,j)      = qsq(k,i)
          tsq_p(i,k,j)      = tsq(k,i)
          qkeadv_p(i,k,j)   = qke_adv(k,i)
          sh3d_p(i,k,j)     = sh3d(k,i)
          sm3d_p(i,k,j)     = sm3d(k,i)
          cldfrabl_p(i,k,j) = cldfrac_bl(k,i)
          qcbl_p(i,k,j)     = qc_bl(k,i)
          qibl_p(i,k,j)     = qi_bl(k,i)
          edmfa_p(i,k,j)    = edmf_a(k,i)
          edmfent_p(i,k,j)  = edmf_ent(k,i)
          edmfqc_p(i,k,j)   = edmf_qc(k,i)
          edmfqt_p(i,k,j)   = edmf_qt(k,i)
          edmfthl_p(i,k,j)  = edmf_thl(k,i)
          edmfw_p(i,k,j)    = edmf_w(k,i)
          subthl_p(i,k,j)   = sub_thl(k,i)
          subqv_p(i,k,j)    = sub_qv(k,i)
          detthl_p(i,k,j)   = det_thl(k,i)
          detqv_p(i,k,j)    = det_qv(k,i)
          dqke_p(i,k,j)     = 0._RKIND
          qbuoy_p(i,k,j)    = 0._RKIND
          qdiss_p(i,k,j)    = 0._RKIND
          qshear_p(i,k,j)   = 0._RKIND
          qwt_p(i,k,j)      = 0._RKIND

          rqsblten_p(i,k,j)   = 0._RKIND
          rncblten_p(i,k,j)   = 0._RKIND
          rniblten_p(i,k,j)   = 0._RKIND
          rnifablten_p(i,k,j) = 0._RKIND
          rnwfablten_p(i,k,j) = 0._RKIND

          pattern_spp_pbl(i,k,j) = 0._RKIND

       enddo
       enddo
       enddo

       do j = jts,jte
       do i = its,ite
          maxwidthbl_p(i,j) = 0._RKIND
          maxmfbl_p(i,j)    = 0._RKIND
          zbl_plume_p(i,j)  = 0._RKIND
          excessh_p(i,j)    = 0._RKIND
          excessq_p(i,j)    = 0._RKIND
       enddo
       enddo

    case("bl_myj")
       !input variables:
       call mpas_pool_get_array(mesh        ,'zgrid' ,zgrid)
       call mpas_pool_get_array(diag_physics, 'lowlyr', lowlyr)
       call mpas_pool_get_array(diag_physics, 'chlowq', chlowq)
       call mpas_pool_get_array(diag_physics, 'lh', lh)
       call mpas_pool_get_array(sfc_input   ,'skintemp',skintemp)
       !inout
       call mpas_pool_get_array(diag_physics, 'thz0', thz0)
       call mpas_pool_get_array(diag_physics, 'uz0', uz0)
       call mpas_pool_get_array(diag_physics, 'vz0', vz0)
       call mpas_pool_get_array(diag_physics, 'qz0', qz0)
       call mpas_pool_get_array(diag_physics, 'akhs', akhs)
       call mpas_pool_get_array(diag_physics, 'akms', akms)
       call mpas_pool_get_array(diag_physics,'tke_pbl',tke_pbl  )
       call mpas_pool_get_array(diag_physics,'qsfc',qsfc)
       !out
       call mpas_pool_get_array(diag_physics, 'ct', ct)
       call mpas_pool_get_array(diag_physics, 'mixht', mixht)

       do j = jts, jte
       do k = kts, kte
       do i = its, ite
          if(k == kts) then
             !input variables:
             ht_p(i,j) = zgrid(1,i)
             lowlyr_p(i,j) = lowlyr(i)
             chlowq_p(i,j) = chlowq(i)
             lh_p(i,j) = lh(i)
             tsk_p(i,j) = skintemp(i)
             !inout variables
             thz0_p(i,j) = thz0(i)
             qz0_p(i,j) = qz0(i)
             uz0_p(i,j) = uz0(i)
             vz0_p(i,j) = vz0(i)
             akms_p(i,j) = akms(i)
             akhs_p(i,j) = akhs(i)
             qsfc_p(i,j) = qsfc(i)
             !out variables
             ct_p(i,j) = 0._RKIND
             mixht_p(i,j) = 0._RKIND
          endif
          tkepbl_p(i,k,j) = tke_pbl(k,i)
       enddo
       enddo
       enddo
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

    kzh_p(i,k,j) = 0._RKIND
    kzm_p(i,k,j) = 0._RKIND
    kzq_p(i,k,j) = 0._RKIND
 enddo
 enddo
 enddo

 end subroutine pbl_from_MPAS

!=================================================================================================================
 subroutine pbl_to_MPAS(configs,state,diag_physics,tend_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics
 type(mpas_pool_type),intent(inout):: state

 integer,intent(in):: its,ite

!local variables:
 integer:: i,k,j,n

!local pointers:
 character(len=StrKIND),pointer:: pbl_scheme

 integer,dimension(:),pointer:: kpbl

 real(kind=RKIND),dimension(:),pointer  :: hpbl
 real(kind=RKIND),dimension(:,:),pointer:: kzh,kzm,kzq
 real(kind=RKIND),dimension(:,:),pointer:: rublten,rvblten,rthblten,rqvblten,rqcblten,rqiblten,rqsblten
 real(kind=RKIND),dimension(:,:),pointer:: rncblten,rniblten,rnifablten,rnwfablten

!local pointers for YSU scheme:
 real(kind=RKIND),dimension(:,:),pointer:: exch_h
 real(kind=RKIND),dimension(:),pointer  :: delta,wstar
      
!local pointers for MYNN scheme:
 real(kind=RKIND),dimension(:),pointer  :: maxmf,maxwidth,ztop_plume,excess_h,excess_q
 real(kind=RKIND),dimension(:,:),pointer:: cov,qke,qsq,tsq,sh3d,sm3d,tke_pbl,qke_adv,el_pbl,dqke,qbuoy, &
                                           qdiss,qshear,qwt
 real(kind=RKIND),dimension(:,:),pointer:: cldfrac_bl,qc_bl,qi_bl
 real(kind=RKIND),dimension(:,:),pointer:: edmf_a,edmf_ent,edmf_qc,edmf_qt,edmf_thl,edmf_w
 real(kind=RKIND),dimension(:,:),pointer:: sub_thl,sub_qv,det_thl,det_qv

!local pointers for tke_edmf
! tke is stored as scalars(5,k,i) for bl_tke_edmf_in.

!local pointers for chemistry
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars
 real(kind=RKIND),dimension(:,:,:),pointer:: chem
 real(kind=RKIND),dimension(:,:,:),pointer:: tend_chem_settle
 real(kind=RKIND),dimension(:,:),pointer  :: ddvel
 integer,pointer :: chemistry_start,chemistry_end

!local pointers for MYJ scheme:
 real(kind=RKIND),dimension(:),pointer   :: thz0,qz0,uz0,vz0,akhs,akms,ct,mixht,qsfc,ust,znt

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_pbl_scheme',pbl_scheme)

 call mpas_pool_get_array(diag_physics,'kpbl' ,kpbl )
 call mpas_pool_get_array(diag_physics,'hpbl' ,hpbl )
 call mpas_pool_get_array(diag_physics,'kzh'  ,kzh  )
 call mpas_pool_get_array(diag_physics,'kzm'  ,kzm  )
 call mpas_pool_get_array(diag_physics,'kzq'  ,kzq  )

 call mpas_pool_get_array(tend_physics,'rublten' ,rublten )
 call mpas_pool_get_array(tend_physics,'rvblten' ,rvblten )
 call mpas_pool_get_array(tend_physics,'rthblten',rthblten)
 call mpas_pool_get_array(tend_physics,'rqvblten',rqvblten)
 call mpas_pool_get_array(tend_physics,'rqcblten',rqcblten)
 call mpas_pool_get_array(tend_physics,'rqiblten',rqiblten)

 do j = jts,jte
 do i = its,ite
    hpbl(i) = hpbl_p(i,j)
    kpbl(i) = kpbl_p(i,j)
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

!   kzh(k,i) = kzh_p(i,k,j)
!   kzm(k,i) = kzm_p(i,k,j)
!   kzq(k,i) = kzh_p(i,k,j)
     ! For bl_tke_edmf, do not touch MPAS K-diagnostic fields here.
! The TKE-EDMF wrapper used in this branch does not return a native
! MPAS-safe exch_h/kzh/kzm/kzq diagnostic.  Writing these fields caused
! a later crash in mpas_halo.  Keep only the physical tendencies.
if (trim(pbl_scheme) /= "bl_tke_edmf") then
   kzh(k,i) = kzh_p(i,k,j)
   kzm(k,i) = kzm_p(i,k,j)
   kzq(k,i) = kzh_p(i,k,j)
endif
 enddo
 enddo
 enddo

 if (do_chemistry ) then
    call mpas_pool_get_array(state,'scalars',scalars,1)
    call mpas_pool_get_dimension(state, 'chemistry_start', chemistry_start)
    call mpas_pool_get_dimension(state, 'chemistry_end', chemistry_end)
    chem => scalars(chemistry_start:chemistry_end,:,:)
    do n = 1,num_chem
    do j = jts,jte
    do k = kts,kte
    do i = its,ite
       chem(n,k,i)    = chem_p(i,k,j,n)
    enddo
    enddo
    enddo
    enddo
 endif

 pbl_select: select case (trim(pbl_scheme))


    case("bl_tke_edmf")
! ------------------------------------------------------------------
! Copy updated prognostic TKE back to MPAS state scalars.
!
! Registry order in var_array scalars for this configuration is:
!   1=qv, 2=qc, 3=qr, 4=qi, 5=tke
! This path is only active for bl_tke_edmf, so other PBL schemes are unchanged.
! ------------------------------------------------------------------
       call mpas_pool_get_array(state,'scalars',scalars,1)

       do j = jts,jte
       do k = kts,kte
       do i = its,ite
          scalars(5,k,i) = max(1.0e-9_RKIND, tkepbl_p(i,k,j))
       enddo
       enddo
       enddo

    case("bl_ysu")
       call mpas_pool_get_array(diag_physics,'delta',delta  )
       call mpas_pool_get_array(diag_physics,'wstar' ,wstar )
       call mpas_pool_get_array(diag_physics,'exch_h',exch_h)

       do j = jts,jte
       do i = its,ite
          delta(i) = delta_p(i,j)
          wstar(i) = wstar_p(i,j)
       enddo
       do k = kts,kte
       do i = its,ite
          exch_h(k,i) = exch_p(i,k,j)
       enddo
       enddo
       enddo

    case("bl_mynn")
       call mpas_pool_get_array(diag_physics,'el_pbl'    ,el_pbl    )
       call mpas_pool_get_array(diag_physics,'cov'       ,cov       )
       call mpas_pool_get_array(diag_physics,'qke'       ,qke       )
       call mpas_pool_get_array(diag_physics,'qke_adv'   ,qke_adv   )
       call mpas_pool_get_array(diag_physics,'qsq'       ,qsq       )
       call mpas_pool_get_array(diag_physics,'tsq'       ,tsq       )
       call mpas_pool_get_array(diag_physics,'tke_pbl'   ,tke_pbl   )
       call mpas_pool_get_array(diag_physics,'sh3d'      ,sh3d      )
       call mpas_pool_get_array(diag_physics,'sm3d'      ,sm3d      )
       call mpas_pool_get_array(diag_physics,'dqke'      ,dqke      )
       call mpas_pool_get_array(diag_physics,'qbuoy'     ,qbuoy     )
       call mpas_pool_get_array(diag_physics,'qdiss'     ,qdiss     )
       call mpas_pool_get_array(diag_physics,'qshear'    ,qshear    )
       call mpas_pool_get_array(diag_physics,'qwt'       ,qwt       )
       call mpas_pool_get_array(diag_physics,'cldfrac_bl',cldfrac_bl)
       call mpas_pool_get_array(diag_physics,'qc_bl'     ,qc_bl     )
       call mpas_pool_get_array(diag_physics,'qi_bl'     ,qi_bl     )
       call mpas_pool_get_array(diag_physics,'edmf_a'    ,edmf_a    )
       call mpas_pool_get_array(diag_physics,'edmf_ent'  ,edmf_ent  )
       call mpas_pool_get_array(diag_physics,'edmf_qc'   ,edmf_qc   )
       call mpas_pool_get_array(diag_physics,'edmf_qt'   ,edmf_qt   )
       call mpas_pool_get_array(diag_physics,'edmf_thl'  ,edmf_thl  )
       call mpas_pool_get_array(diag_physics,'edmf_w'    ,edmf_w    )
       call mpas_pool_get_array(diag_physics,'sub_thl'   ,sub_thl   )
       call mpas_pool_get_array(diag_physics,'sub_qv'    ,sub_qv    )
       call mpas_pool_get_array(diag_physics,'det_thl'   ,det_thl   )
       call mpas_pool_get_array(diag_physics,'det_qv'    ,det_qv    )

       call mpas_pool_get_array(tend_physics,'rqsblten'  ,rqsblten  )

       do j = jts,jte
       do k = kts,kte
       do i = its,ite
          el_pbl(k,i)     = elpbl_p(i,k,j)
          cov(k,i)        = cov_p(i,k,j)
          qke(k,i)        = qke_p(i,k,j)
          qsq(k,i)        = qsq_p(i,k,j)
          tsq(k,i)        = tsq_p(i,k,j)
          sh3d(k,i)       = sh3d_p(i,k,j)
          sm3d(k,i)       = sm3d_p(i,k,j)
          tke_pbl(k,i)    = tkepbl_p(i,k,j)
          qke_adv(k,i)    = qkeadv_p(i,k,j)
          cldfrac_bl(k,i) = cldfrabl_p(i,k,j)
          qc_bl(k,i)      = qcbl_p(i,k,j)
          qi_bl(k,i)      = qibl_p(i,k,j)
          edmf_a(k,i)     = edmfa_p(i,k,j)
          edmf_ent(k,i)   = edmfent_p(i,k,j)
          edmf_qc(k,i)    = edmfqc_p(i,k,j)
          edmf_qt(k,i)    = edmfqt_p(i,k,j)
          edmf_thl(k,i)   = edmfthl_p(i,k,j)
          edmf_w(k,i)     = edmfw_p(i,k,j)
          sub_thl(k,i)    = subthl_p(i,k,j)
          sub_qv(k,i)     = subqv_p(i,k,j)
          det_thl(k,i)    = detthl_p(i,k,j)
          det_qv(k,i)     = detqv_p(i,k,j)
          dqke(k,i)       = dqke_p(i,k,j)
          qbuoy(k,i)      = qbuoy_p(i,k,j)
          qdiss(k,i)      = qdiss_p(i,k,j)
          qshear(k,i)     = qshear_p(i,k,j)
          qwt(k,i)        = qwt_p(i,k,j)

          rqsblten(k,i)   = rqsblten_p(i,k,j)
       enddo
       enddo
       enddo

       if(f_ni) then
          call mpas_pool_get_array(tend_physics,'rniblten',rniblten)
          do j = jts,jte
             do k = kts,kte
                do i = its,ite
                   rniblten(k,i) = rniblten_p(i,k,j)
                enddo
             enddo
          enddo
       endif
       if(f_nc .and. f_nifa .and. f_nwfa) then
          call mpas_pool_get_array(tend_physics,'rncblten'  ,rncblten  )
          call mpas_pool_get_array(tend_physics,'rnifablten',rnifablten)
          call mpas_pool_get_array(tend_physics,'rnwfablten',rnwfablten)
          do j = jts,jte
             do k = kts,kte
                do i = its,ite
                   rncblten(k,i)   = rncblten_p(i,k,j)
                   rnifablten(k,i) = rnifablten_p(i,k,j)
                   rnwfablten(k,i) = rnwfablten_p(i,k,j)
                enddo
             enddo
          enddo
       endif
       
    case("bl_mynnedmf")
       call mpas_pool_get_array(diag_physics,'maxmf'     ,maxmf     )
       call mpas_pool_get_array(diag_physics,'maxwidth'  ,maxwidth  )
       call mpas_pool_get_array(diag_physics,'ztop_plume',ztop_plume)
       call mpas_pool_get_array(diag_physics,'excess_h'  ,excess_h  )
       call mpas_pool_get_array(diag_physics,'excess_q'  ,excess_q  )
       call mpas_pool_get_array(diag_physics,'el_pbl'    ,el_pbl    )
       call mpas_pool_get_array(diag_physics,'cov'       ,cov       )
       call mpas_pool_get_array(diag_physics,'qke'       ,qke       )
       call mpas_pool_get_array(diag_physics,'qke_adv'   ,qke_adv   )
       call mpas_pool_get_array(diag_physics,'qsq'       ,qsq       )
       call mpas_pool_get_array(diag_physics,'tsq'       ,tsq       )
       call mpas_pool_get_array(diag_physics,'sh3d'      ,sh3d      )
       call mpas_pool_get_array(diag_physics,'sm3d'      ,sm3d      )
       call mpas_pool_get_array(diag_physics,'dqke'      ,dqke      )
       call mpas_pool_get_array(diag_physics,'qbuoy'     ,qbuoy     )
       call mpas_pool_get_array(diag_physics,'qdiss'     ,qdiss     )
       call mpas_pool_get_array(diag_physics,'qshear'    ,qshear    )
       call mpas_pool_get_array(diag_physics,'qwt'       ,qwt       )
       call mpas_pool_get_array(diag_physics,'cldfrac_bl',cldfrac_bl)
       call mpas_pool_get_array(diag_physics,'qc_bl'     ,qc_bl     )
       call mpas_pool_get_array(diag_physics,'qi_bl'     ,qi_bl     )
       call mpas_pool_get_array(diag_physics,'edmf_a'    ,edmf_a    )
       call mpas_pool_get_array(diag_physics,'edmf_ent'  ,edmf_ent  )
       call mpas_pool_get_array(diag_physics,'edmf_qc'   ,edmf_qc   )
       call mpas_pool_get_array(diag_physics,'edmf_qt'   ,edmf_qt   )
       call mpas_pool_get_array(diag_physics,'edmf_thl'  ,edmf_thl  )
       call mpas_pool_get_array(diag_physics,'edmf_w'    ,edmf_w    )
       call mpas_pool_get_array(diag_physics,'sub_thl'   ,sub_thl   )
       call mpas_pool_get_array(diag_physics,'sub_qv'    ,sub_qv    )
       call mpas_pool_get_array(diag_physics,'det_thl'   ,det_thl   )
       call mpas_pool_get_array(diag_physics,'det_qv'    ,det_qv    )

       call mpas_pool_get_array(tend_physics,'rqsblten'  ,rqsblten  )

       do j = jts,jte
       do i = its,ite
          maxmf(i)      = maxmfbl_p(i,j)
          maxwidth(i)   = maxwidthbl_p(i,j)
          ztop_plume(i) = zbl_plume_p(i,j)
          excess_h(i)   = excessh_p(i,j)
          excess_q(i)   = excessq_p(i,j)
       enddo
       enddo

       do j = jts,jte
       do k = kts,kte
       do i = its,ite
          el_pbl(k,i)     = elpbl_p(i,k,j)
          cov(k,i)        = cov_p(i,k,j)
          qke(k,i)        = qke_p(i,k,j)
          qsq(k,i)        = qsq_p(i,k,j)
          tsq(k,i)        = tsq_p(i,k,j)
          sh3d(k,i)       = sh3d_p(i,k,j)
          sm3d(k,i)       = sm3d_p(i,k,j)
          qke_adv(k,i)    = qkeadv_p(i,k,j)
          cldfrac_bl(k,i) = cldfrabl_p(i,k,j)
          qc_bl(k,i)      = qcbl_p(i,k,j)
          qi_bl(k,i)      = qibl_p(i,k,j)
          edmf_a(k,i)     = edmfa_p(i,k,j)
          edmf_ent(k,i)   = edmfent_p(i,k,j)
          edmf_qc(k,i)    = edmfqc_p(i,k,j)
          edmf_qt(k,i)    = edmfqt_p(i,k,j)
          edmf_thl(k,i)   = edmfthl_p(i,k,j)
          edmf_w(k,i)     = edmfw_p(i,k,j)
          sub_thl(k,i)    = subthl_p(i,k,j)
          sub_qv(k,i)     = subqv_p(i,k,j)
          det_thl(k,i)    = detthl_p(i,k,j)
          det_qv(k,i)     = detqv_p(i,k,j)
          dqke(k,i)       = dqke_p(i,k,j)
          qbuoy(k,i)      = qbuoy_p(i,k,j)
          qdiss(k,i)      = qdiss_p(i,k,j)
          qshear(k,i)     = qshear_p(i,k,j)
          qwt(k,i)        = qwt_p(i,k,j)

          rqsblten(k,i)   = rqsblten_p(i,k,j)
       enddo
       enddo
       enddo

       if(f_ni) then
          call mpas_pool_get_array(tend_physics,'rniblten',rniblten)
          do j = jts,jte
             do k = kts,kte
                do i = its,ite
                   rniblten(k,i) = rniblten_p(i,k,j)
                enddo
             enddo
          enddo
       endif
       if(f_nc .and. f_nifa .and. f_nwfa) then
          call mpas_pool_get_array(tend_physics,'rncblten'  ,rncblten  )
          call mpas_pool_get_array(tend_physics,'rnifablten',rnifablten)
          call mpas_pool_get_array(tend_physics,'rnwfablten',rnwfablten)
          do j = jts,jte
             do k = kts,kte
                do i = its,ite
                   rncblten(k,i)   = rncblten_p(i,k,j)
                   rnifablten(k,i) = rnifablten_p(i,k,j)
                   rnwfablten(k,i) = rnwfablten_p(i,k,j)
                enddo
             enddo
          enddo
       endif
       
    case("bl_myj")
       !inout
       call mpas_pool_get_array(diag_physics, 'thz0', thz0)
       call mpas_pool_get_array(diag_physics, 'uz0', uz0)
       call mpas_pool_get_array(diag_physics, 'vz0', vz0)
       call mpas_pool_get_array(diag_physics, 'qz0', qz0)
       call mpas_pool_get_array(diag_physics, 'akhs', akhs)
       call mpas_pool_get_array(diag_physics, 'akms', akms)
       call mpas_pool_get_array(diag_physics, 'tke_pbl', tke_pbl)
       call mpas_pool_get_array(diag_physics, 'qsfc', qsfc)
       call mpas_pool_get_array(diag_physics, 'ust', ust)
       call mpas_pool_get_array(diag_physics, 'znt', znt)
       !out
       call mpas_pool_get_array(diag_physics, 'ct', ct)
       call mpas_pool_get_array(diag_physics, 'mixht', mixht)
       call mpas_pool_get_array(diag_physics,'el_pbl'  ,el_pbl  )
       do j = jts, jte
       do k = kts,kte
       do i = its, ite
          if (k==kts) then
             !inout variables
             thz0(i) = thz0_p(i,j)
             qz0(i) = qz0_p(i,j)
             uz0(i) = uz0_p(i,j)
             vz0(i) = vz0_p(i,j)
             akms(i) = akms_p(i,j)
             akhs(i) = akhs_p(i,j)
             qsfc(i) = qsfc_p(i,j)
             ust(i) = ust_p(i,j)
             znt(i) = znt_p(i,j)
             !out variables
             ct(i) = ct_p(i,j)
             mixht(i) = mixht_p(i,j)
          endif
          tke_pbl(k,i)  = tkepbl_p(i,k,j)
          el_pbl(k,i)   = elpbl_p(i,k,j)
       enddo
       enddo
       enddo

    case default

 end select pbl_select

 end subroutine pbl_to_MPAS
 
!=================================================================================================================
 subroutine init_pbl(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

!local variables and pointers:
 character(len=StrKIND),pointer:: pbl_scheme
 character(len=StrKIND):: errmsg
 integer:: errflg

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_pbl_scheme',pbl_scheme)

 pbl_select: select case (trim(pbl_scheme))

    case("bl_mynn")
!      call mpas_log_write('--- enter subroutine bl_mynn_init:')
       call bl_mynn_init(cp,cpv,cice,cliq,ep_1,ep_2,gravity,karman,P0,R_d,R_v,svp1,svp2,svp3,svpt0, &
                         xlf,xls,xlv,errmsg,errflg)
!      call mpas_log_write('--- end subroutine bl_mynn_init:')

    case default

 end select pbl_select

 end subroutine init_pbl

!=================================================================================================================
 subroutine driver_pbl(itimestep,configs,state,mesh,sfc_input,diag_physics,tend_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(inout):: state

 integer,intent(in):: its,ite
 integer,intent(in):: itimestep

!inout arguments:
 type(mpas_pool_type),intent(inout):: sfc_input
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

!local pointers:
 logical,pointer:: config_do_DAcycling, &
                   config_do_restart,   &
                   bl_mynn_tkeadvect

 logical,pointer:: enh_mix

 character(len=StrKIND),pointer:: pbl_scheme

 integer,pointer:: bl_mynn_cloudpdf,    &
                   bl_mynn_mixlength,   &
                   bl_mynn_stfunc,      &
                   bl_mynn_topdown,     &
                   bl_mynn_scaleaware,  &
                   bl_mynn_dheat_opt,   &
                   bl_mynn_edmf,        &
                   bl_mynn_edmf_dd,     &
                   bl_mynn_edmf_mom,    &
                   bl_mynn_edmf_tke,    &
                   bl_mynn_edmf_output, &
                   bl_mynn_mixscalars,  &
                   bl_mynn_mixaerosols, &
                   bl_mynn_mixnumcon,   &
                   bl_mynn_cloudmix,    &
                   bl_mynn_mixqt,       &
                   bl_mynn_ess,         &
                   bl_mynn_tkebudget,   &
                   icloud_bl,           &
                   spp_pbl

 real(kind=RKIND),pointer:: bl_mynn_closure

!local variables:
 integer:: initflag
 integer:: i,k,j
 logical:: mix_chem

!CCPP-compliant flags:
 character(len=StrKIND):: errmsg
 integer:: errflg

!UFS TKE-EDMF adapter variables:
 type(mpas_satmedmfvdifq_config_type) :: tke_edmf_cfg
 character(len=256) :: tke_edmf_errmsg
 integer :: tke_edmf_errflg
 real(kind=RKIND), allocatable :: coszen_tkeedmf_p(:,:)
 real(kind=RKIND), allocatable :: lh_tkeedmf_p(:,:)
 real(kind=RKIND), allocatable :: stress_tkeedmf_p(:,:)
 real(kind=RKIND), allocatable :: zmid_tkeedmf_p(:,:,:)
 real(kind=RKIND), allocatable :: zint_tkeedmf_p(:,:,:)

! Rare-column TKE-EDMF safeguard. Use actual MPAS cell connectivity
! to replace only invalid PBL outputs with neighboring valid PBL outputs.
 integer, dimension(:), pointer :: nEdgesOnCell
 integer, dimension(:,:), pointer :: cellsOnCell
 logical, allocatable :: pbl_col_valid(:), pbl_col_bad(:)
 logical, allocatable :: pbl_press_valid(:), pbl_press_bad(:)
 real(kind=RKIND), allocatable :: p_mid_tkeedmf(:,:), p_int_tkeedmf(:,:)
 integer :: n, nbr, ngood, nbad_pbl, nbad_press, bad_reason
 integer, save :: pbl_repair_total_local = 0
 integer, save :: pbl_repair_window_local = 0
 integer, save :: pbl_pressure_total_local = 0
 integer, save :: pbl_pressure_window_local = 0
 integer, save :: nbad_nan_window  = 0
 integer, save :: nbad_hpbl_window = 0
 integer, save :: nbad_t_window    = 0
 integer, save :: nbad_qv_window   = 0
 integer, save :: nqv_limited_window = 0
 integer, save :: nqv_limited_total  = 0
 integer, save :: nqc_limited_window = 0
 integer, save :: nqc_limited_total  = 0
 integer, save :: nqi_limited_window = 0
 integer, save :: nqi_limited_total  = 0
 integer, save :: nbad_qc_window   = 0
 integer, save :: nbad_qi_window   = 0
 integer, save :: nbad_uv_window   = 0
 integer, save :: nbad_tke_window  = 0
 real(kind=RKIND) :: t_after, qv_after, qc_after, qi_after
 real(kind=RKIND) :: min_dp_press, p_lower, p_upper
 real(kind=RKIND) :: u_after, v_after
 real(kind=RKIND) :: sum_hpbl
 integer :: sum_kpbl

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('')
!call mpas_log_write('--- enter subroutine driver_pbl:')

!initialization of CCPP-compliant flags:
 errmsg = ' '
 errflg = 0

 mix_chem  = .false.

 call mpas_pool_get_config(configs,'config_do_DAcycling',config_do_DAcycling)
 call mpas_pool_get_config(configs,'config_do_restart'  ,config_do_restart  )
 call mpas_pool_get_config(configs,'config_pbl_scheme'  ,pbl_scheme         )

!copy MPAS arrays to local arrays:
 call pbl_from_MPAS(configs,state,mesh,sfc_input,diag_physics,tend_physics,its,ite)

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


    case("bl_tke_edmf")
       call mpas_timer_start('bl_tke_edmf')

       call mpas_pool_get_array(mesh,'nEdgesOnCell',nEdgesOnCell)
       call mpas_pool_get_array(mesh,'cellsOnCell' ,cellsOnCell )

       allocate(pbl_col_valid(ims:ime), pbl_col_bad(ims:ime))
       allocate(pbl_press_valid(ims:ime), pbl_press_bad(ims:ime))
       allocate(p_mid_tkeedmf(its:ite,kts:kte))
       allocate(p_int_tkeedmf(its:ite,kts:kte+1))
       pbl_col_valid   = .false.
       pbl_col_bad     = .false.
       pbl_press_valid = .false.
       pbl_press_bad   = .false.
       min_dp_press    = 1.0_RKIND

       ! Conservative initial settings.  Expose these through Registry later.
       tke_edmf_cfg%sa3dtke     = .false.
       tke_edmf_cfg%tte_edmf    = .false.
       tke_edmf_cfg%dspheat     = .true.
       tke_edmf_cfg%do_canopy   = .false.
       tke_edmf_cfg%cplaqm      = .false.
       tke_edmf_cfg%gen_tend    = .false.
       tke_edmf_cfg%ldiag3d     = .false.

       tke_edmf_cfg%sfc_rlm = 0
       tke_edmf_cfg%tc_pbl  = 0
       tke_edmf_cfg%use_lpt = 0

       tke_edmf_cfg%xkzm_m  = 1.0_RKIND
       tke_edmf_cfg%xkzm_h  = 1.0_RKIND
       tke_edmf_cfg%xkzm_s  = 1.0_RKIND
       tke_edmf_cfg%dspfac  = 1.0_RKIND
       tke_edmf_cfg%bl_upfr = 0.07_RKIND
       tke_edmf_cfg%bl_dnfr = 0.05_RKIND
       tke_edmf_cfg%rlmx    = 300.0_RKIND
       tke_edmf_cfg%elmx    = 300.0_RKIND

       tke_edmf_errmsg = ''
       tke_edmf_errflg = 0

       if (.not. allocated(coszen_tkeedmf_p)) then
          allocate(coszen_tkeedmf_p(ims:ime,jms:jme))
          allocate(lh_tkeedmf_p(ims:ime,jms:jme))
          allocate(stress_tkeedmf_p(ims:ime,jms:jme))
          allocate(zmid_tkeedmf_p(ims:ime,kms:kme,jms:jme))
          allocate(zint_tkeedmf_p(ims:ime,kms:kme+1,jms:jme))
       endif

       do j = jts,jte
       do i = its,ite
          coszen_tkeedmf_p(i,j) = 1.0_RKIND
          lh_tkeedmf_p(i,j) = qfx_p(i,j) * xlv
          stress_tkeedmf_p(i,j) = ust_p(i,j) * ust_p(i,j)

          ! Build geometric heights from MPAS layer thickness.
          ! zint is interface height and zmid is layer-center height.
          zint_tkeedmf_p(i,kts,j) = 0.0_RKIND
       enddo
       do k = kts,kte
       do i = its,ite
          zint_tkeedmf_p(i,k+1,j) = zint_tkeedmf_p(i,k,j)              &
                                  + dz_p(i,k,j)
          zmid_tkeedmf_p(i,k,j) = 0.5_RKIND *                         &
             (zint_tkeedmf_p(i,k,j) + zint_tkeedmf_p(i,k+1,j))
       enddo
       enddo
       enddo
       !---------------------------------------------------------------
       ! TEST: repair a pressure column before TKE-EDMF only when the
       ! interface pressure is non-finite, non-positive, or fails to
       ! decrease upward.  The repair uses actual MPAS cellsOnCell
       ! neighbors.  The MPAS pressure state itself is not overwritten;
       ! only these temporary arrays are passed to the PBL wrapper.
       !---------------------------------------------------------------
       p_mid_tkeedmf(:,:) = pres_hyd_p(its:ite,kts:kte,jts)
       p_int_tkeedmf(:,:) = pres2_hyd_p(its:ite,kts:kte+1,jts)

       ! Determine valid pressure profiles for owned and halo cells so
       ! that a bad owned cell can use valid physical neighbors.
       do i = ims, ime
          pbl_press_valid(i) = .true.
          do k = kts, kte
             if (.not. ieee_is_finite(pres_hyd_p(i,k,jts)) .or. &
                 .not. ieee_is_finite(pres2_hyd_p(i,k,jts)) .or. &
                 .not. ieee_is_finite(pres2_hyd_p(i,k+1,jts))) then
                pbl_press_valid(i) = .false.
                exit
             endif
             if (pres2_hyd_p(i,k,jts) <= 0.0_RKIND .or. &
                 pres2_hyd_p(i,k+1,jts) <= 0.0_RKIND .or. &
                 pres2_hyd_p(i,k+1,jts) >= pres2_hyd_p(i,k,jts) .or. &
                 pres_hyd_p(i,k,jts) >= pres2_hyd_p(i,k,jts) .or. &
                 pres_hyd_p(i,k,jts) <= pres2_hyd_p(i,k+1,jts)) then
                pbl_press_valid(i) = .false.
                exit
             endif
          enddo
       enddo

       do i = its, ite
          pbl_press_bad(i) = .not. pbl_press_valid(i)
       enddo
       nbad_press = count(pbl_press_bad(its:ite))
       pbl_pressure_total_local  = pbl_pressure_total_local  + nbad_press
       pbl_pressure_window_local = pbl_pressure_window_local + nbad_press

       do i = its, ite
          if (.not. pbl_press_bad(i)) cycle

          p_mid_tkeedmf(i,kts:kte)   = 0.0_RKIND
          p_int_tkeedmf(i,kts:kte+1) = 0.0_RKIND
          ngood = 0

          do n = 1, nEdgesOnCell(i)
             nbr = cellsOnCell(n,i)
             if (nbr < ims .or. nbr > ime) cycle
             if (.not. pbl_press_valid(nbr)) cycle

             ngood = ngood + 1
             p_mid_tkeedmf(i,kts:kte) = p_mid_tkeedmf(i,kts:kte) + &
                  pres_hyd_p(nbr,kts:kte,jts)
             p_int_tkeedmf(i,kts:kte+1) = p_int_tkeedmf(i,kts:kte+1) + &
                  pres2_hyd_p(nbr,kts:kte+1,jts)
          enddo

          if (ngood > 0) then
             p_mid_tkeedmf(i,kts:kte) = p_mid_tkeedmf(i,kts:kte) / &
                  real(ngood,RKIND)
             p_int_tkeedmf(i,kts:kte+1) = p_int_tkeedmf(i,kts:kte+1) / &
                  real(ngood,RKIND)
          else
             ! Emergency fallback when no valid neighbor is locally
             ! available.  Preserve the original surface interface and
             ! force a minimally monotonic profile only for this PBL call.
             p_int_tkeedmf(i,kts) = max(100.0_RKIND, &
                  pres2_hyd_p(i,kts,jts))
             do k = kts, kte
                p_int_tkeedmf(i,k+1) = min(p_int_tkeedmf(i,k)-min_dp_press, &
                     max(1.0_RKIND,pres2_hyd_p(i,k+1,jts)))
                p_mid_tkeedmf(i,k) = 0.5_RKIND * &
                     (p_int_tkeedmf(i,k)+p_int_tkeedmf(i,k+1))
             enddo
          endif

          ! Final strict monotonicity and midpoint check after averaging.
          do k = kts, kte
             p_lower = p_int_tkeedmf(i,k)
             p_upper = p_int_tkeedmf(i,k+1)
             if (p_upper >= p_lower-min_dp_press) then
                p_upper = max(1.0_RKIND,p_lower-min_dp_press)
                p_int_tkeedmf(i,k+1) = p_upper
             endif
             if (p_mid_tkeedmf(i,k) >= p_lower .or. &
                 p_mid_tkeedmf(i,k) <= p_upper .or. &
                 .not. ieee_is_finite(p_mid_tkeedmf(i,k))) then
                p_mid_tkeedmf(i,k) = 0.5_RKIND*(p_lower+p_upper)
             endif
          enddo
       enddo

       ! Quiet test: pressure repair remains active, but its routine
       ! per-step diagnostic is suppressed.  Only the post-PBL neighbor
       ! repair count is printed below.

       call mpas_call_satmedmfvdifq(                                  &
            (ite-its+1), (kte-kts+1), 4, dt_pbl,                      &
            zmid_tkeedmf_p(its:ite,kts:kte,jts),                                &
            zint_tkeedmf_p(its:ite,kts:kte+1,jts),                              &
            dx_p(its:ite,jts),                                        &
            u_p(its:ite,kts:kte,jts),                                 &
            v_p(its:ite,kts:kte,jts),                                 &
            t_p(its:ite,kts:kte,jts),                                 &
            qv_p(its:ite,kts:kte,jts),                                &
            qc_p(its:ite,kts:kte,jts),                                &
            qi_p(its:ite,kts:kte,jts),                                &
            tkepbl_p(its:ite,kts:kte,jts),                            &
            p_mid_tkeedmf(its:ite,kts:kte),                            &
            p_int_tkeedmf(its:ite,kts:kte+1),                         &
            pi_p(its:ite,kts:kte,jts),                                &
            rthraten_p(its:ite,kts:kte,jts),                          &
            exch_p(its:ite,kts:kte,jts),                              &
            coszen_tkeedmf_p(its:ite,jts),                            &
            tsk_p(its:ite,jts),                                       &
            hfx_p(its:ite,jts),                                       &
            lh_tkeedmf_p(its:ite,jts),                                &
            stress_tkeedmf_p(its:ite,jts),                            &
            znt_p(its:ite,jts),                                       &
            u10_p(its:ite,jts),                                       &
            v10_p(its:ite,jts),                                       &
            vegfra_p(its:ite,jts),                                    &
            br_p(its:ite,jts),                                        &
            psim_p(its:ite,jts),                                      &
            psih_p(its:ite,jts),                                      &
            hpbl_p(its:ite,jts),                                      &
            kpbl_p(its:ite,jts),                                      &
            rthblten_p(its:ite,kts:kte,jts),                          &
            rublten_p(its:ite,kts:kte,jts),                           &
            rvblten_p(its:ite,kts:kte,jts),                           &
            rqvblten_p(its:ite,kts:kte,jts),                          &
            rqcblten_p(its:ite,kts:kte,jts),                          &
            rqiblten_p(its:ite,kts:kte,jts),                          &
            tke_edmf_cfg, tke_edmf_errmsg, tke_edmf_errflg)

       !-----------------------------------------------------------------
       ! TEST: identify columns for which the returned TKE-EDMF tendency
       ! would produce a non-finite or clearly unphysical one-step state.
       ! Do not modify valid columns.  For a bad column, replace the full
       ! PBL tendency profile with the mean from valid cellsOnCell neighbors.
       ! This keeps the atmospheric state unchanged and repairs only the
       ! faulty PBL contribution.
       !-----------------------------------------------------------------
       if (tke_edmf_errflg == 0) then
          do i = its, ite
             pbl_col_valid(i) = .true.
             bad_reason = 0

             if (.not. ieee_is_finite(hpbl_p(i,jts))) then
                bad_reason = 1
             else if (hpbl_p(i,jts) < 0.0_RKIND .or. &
                      hpbl_p(i,jts) > 10000.0_RKIND) then
                bad_reason = 2
             endif

             if (bad_reason == 0) then
                do k = kts, kte
                   if (.not. ieee_is_finite(rthblten_p(i,k,jts)) .or. &
                       .not. ieee_is_finite(rqvblten_p(i,k,jts)) .or. &
                       .not. ieee_is_finite(rqcblten_p(i,k,jts)) .or. &
                       .not. ieee_is_finite(rqiblten_p(i,k,jts)) .or. &
                       .not. ieee_is_finite(rublten_p (i,k,jts)) .or. &
                       .not. ieee_is_finite(rvblten_p (i,k,jts)) .or. &
                       .not. ieee_is_finite(tkepbl_p  (i,k,jts))) then
                      bad_reason = 1
                      exit
                   endif

                   ! rthblten is potential-temperature tendency. Convert the
                   ! one-step impact to physical temperature with Exner pi_p.
                   t_after  = t_p(i,k,jts)  + dt_pbl*pi_p(i,k,jts)*rthblten_p(i,k,jts)
                   qv_after = qv_p(i,k,jts) + dt_pbl*rqvblten_p(i,k,jts)
                   qc_after = qc_p(i,k,jts) + dt_pbl*rqcblten_p(i,k,jts)
                   qi_after = qi_p(i,k,jts) + dt_pbl*rqiblten_p(i,k,jts)
                   u_after  = u_p(i,k,jts)  + dt_pbl*rublten_p(i,k,jts)
                   v_after  = v_p(i,k,jts)  + dt_pbl*rvblten_p(i,k,jts)

                   if (.not. ieee_is_finite(t_after)  .or. &
                       .not. ieee_is_finite(qv_after) .or. &
                       .not. ieee_is_finite(qc_after) .or. &
                       .not. ieee_is_finite(qi_after) .or. &
                       .not. ieee_is_finite(u_after)  .or. &
                       .not. ieee_is_finite(v_after)) then
                      bad_reason = 1
                      exit
                   endif

                   ! Remove harmless roundoff before deciding that the full
                   ! column requires neighbor repair.
                   if (tkepbl_p(i,k,jts) < 0.0_RKIND .and. &
                       tkepbl_p(i,k,jts) >= -1.0e-8_RKIND) then
                      tkepbl_p(i,k,jts) = 0.0_RKIND
                   endif

                   if (t_after < 140.0_RKIND .or. t_after > 400.0_RKIND) then
                      bad_reason = 3
                      exit
                   else if (qv_after < 0.0_RKIND) then
                      ! Positivity limiter: modify only the vapor tendency at
                      ! this level so that qv after the PBL update is exactly
                      ! zero.  Do not replace the entire PBL column.
                      rqvblten_p(i,k,jts) = -qv_p(i,k,jts) / &
                           max(dt_pbl, tiny(1.0_RKIND))
                      qv_after = qv_p(i,k,jts) + dt_pbl*rqvblten_p(i,k,jts)
                      nqv_limited_window = nqv_limited_window + 1
                      nqv_limited_total  = nqv_limited_total  + 1

                      ! Only a nonfinite or materially positive-overflow state
                      ! remains a reason for whole-column neighbor repair.
                      if (.not. ieee_is_finite(qv_after)) then
                         bad_reason = 1
                         exit
                      endif
                   else if (qv_after > 0.10_RKIND) then
                      bad_reason = 4
                      exit
                   else if (qc_after < 0.0_RKIND) then
                      ! Positivity limiter for cloud water: modify only the
                      ! cloud-water tendency at this level so qc_after = 0.
                      rqcblten_p(i,k,jts) = -qc_p(i,k,jts) / &
                           max(dt_pbl, tiny(1.0_RKIND))
                      qc_after = qc_p(i,k,jts) + dt_pbl*rqcblten_p(i,k,jts)
                      nqc_limited_window = nqc_limited_window + 1
                      nqc_limited_total  = nqc_limited_total  + 1
                      if (.not. ieee_is_finite(qc_after)) then
                         bad_reason = 1
                         exit
                      endif
                   else if (qi_after < 0.0_RKIND) then
                      ! Positivity limiter for cloud ice: modify only the
                      ! cloud-ice tendency at this level so qi_after = 0.
                      rqiblten_p(i,k,jts) = -qi_p(i,k,jts) / &
                           max(dt_pbl, tiny(1.0_RKIND))
                      qi_after = qi_p(i,k,jts) + dt_pbl*rqiblten_p(i,k,jts)
                      nqi_limited_window = nqi_limited_window + 1
                      nqi_limited_total  = nqi_limited_total  + 1
                      if (.not. ieee_is_finite(qi_after)) then
                         bad_reason = 1
                         exit
                      endif
                   else if (abs(u_after) > 250.0_RKIND .or. &
                            abs(v_after) > 250.0_RKIND) then
                      bad_reason = 7
                      exit
                   else if (tkepbl_p(i,k,jts) < -1.0e-8_RKIND .or. &
                            tkepbl_p(i,k,jts) > 1000.0_RKIND) then
                      bad_reason = 8
                      exit
                   endif
                enddo
             endif

             pbl_col_valid(i) = (bad_reason == 0)
             pbl_col_bad(i)   = .not. pbl_col_valid(i)

             ! Count only the first failure reason for each repaired column.
             select case (bad_reason)
             case (1)
                nbad_nan_window = nbad_nan_window + 1
             case (2)
                nbad_hpbl_window = nbad_hpbl_window + 1
             case (3)
                nbad_t_window = nbad_t_window + 1
             case (4)
                nbad_qv_window = nbad_qv_window + 1
             case (5)
                nbad_qc_window = nbad_qc_window + 1
             case (6)
                nbad_qi_window = nbad_qi_window + 1
             case (7)
                nbad_uv_window = nbad_uv_window + 1
             case (8)
                nbad_tke_window = nbad_tke_window + 1
             end select
          enddo

          nbad_pbl = count(pbl_col_bad(its:ite))
          pbl_repair_total_local  = pbl_repair_total_local  + nbad_pbl
          pbl_repair_window_local = pbl_repair_window_local + nbad_pbl

          do i = its, ite
             if (.not. pbl_col_bad(i)) cycle

             ! First zero the bad-column outputs.  Valid neighboring columns
             ! are accumulated below.
             rthblten_p(i,kts:kte,jts) = 0.0_RKIND
             rqvblten_p(i,kts:kte,jts) = 0.0_RKIND
             rqcblten_p(i,kts:kte,jts) = 0.0_RKIND
             rqiblten_p(i,kts:kte,jts) = 0.0_RKIND
             rublten_p (i,kts:kte,jts) = 0.0_RKIND
             rvblten_p (i,kts:kte,jts) = 0.0_RKIND
             tkepbl_p  (i,kts:kte,jts) = 0.0_RKIND
             sum_hpbl = 0.0_RKIND
             sum_kpbl = 0
             ngood = 0

             do n = 1, nEdgesOnCell(i)
                nbr = cellsOnCell(n,i)
                if (nbr < its .or. nbr > ite) cycle
                if (.not. pbl_col_valid(nbr)) cycle

                ngood = ngood + 1
                rthblten_p(i,kts:kte,jts) = rthblten_p(i,kts:kte,jts) + &
                                             rthblten_p(nbr,kts:kte,jts)
                rqvblten_p(i,kts:kte,jts) = rqvblten_p(i,kts:kte,jts) + &
                                             rqvblten_p(nbr,kts:kte,jts)
                rqcblten_p(i,kts:kte,jts) = rqcblten_p(i,kts:kte,jts) + &
                                             rqcblten_p(nbr,kts:kte,jts)
                rqiblten_p(i,kts:kte,jts) = rqiblten_p(i,kts:kte,jts) + &
                                             rqiblten_p(nbr,kts:kte,jts)
                rublten_p(i,kts:kte,jts)  = rublten_p(i,kts:kte,jts) + &
                                             rublten_p(nbr,kts:kte,jts)
                rvblten_p(i,kts:kte,jts)  = rvblten_p(i,kts:kte,jts) + &
                                             rvblten_p(nbr,kts:kte,jts)
                tkepbl_p(i,kts:kte,jts)   = tkepbl_p(i,kts:kte,jts) + &
                                             tkepbl_p(nbr,kts:kte,jts)
                sum_hpbl = sum_hpbl + hpbl_p(nbr,jts)
                sum_kpbl = sum_kpbl + kpbl_p(nbr,jts)
             enddo

             if (ngood > 0) then
                rthblten_p(i,kts:kte,jts) = rthblten_p(i,kts:kte,jts) / real(ngood,RKIND)
                rqvblten_p(i,kts:kte,jts) = rqvblten_p(i,kts:kte,jts) / real(ngood,RKIND)
                rqcblten_p(i,kts:kte,jts) = rqcblten_p(i,kts:kte,jts) / real(ngood,RKIND)
                rqiblten_p(i,kts:kte,jts) = rqiblten_p(i,kts:kte,jts) / real(ngood,RKIND)
                rublten_p(i,kts:kte,jts)  = rublten_p(i,kts:kte,jts)  / real(ngood,RKIND)
                rvblten_p(i,kts:kte,jts)  = rvblten_p(i,kts:kte,jts)  / real(ngood,RKIND)
                tkepbl_p(i,kts:kte,jts)   = max(1.0e-9_RKIND, &
                     tkepbl_p(i,kts:kte,jts) / real(ngood,RKIND))
                hpbl_p(i,jts) = sum_hpbl / real(ngood,RKIND)
                kpbl_p(i,jts) = nint(real(sum_kpbl,RKIND) / real(ngood,RKIND))
             else
                ! No valid owned neighbor was available (for example at a
                ! partition boundary).  Return zero PBL tendencies rather
                ! than sending an invalid update to the rest of physics.
                tkepbl_p(i,kts:kte,jts) = max(1.0e-9_RKIND, tkepbl_p(i,kts:kte,jts))
                hpbl_p(i,jts) = max(0.0_RKIND, hpbl_p(i,jts))
                kpbl_p(i,jts) = max(kts, min(kte, kpbl_p(i,jts)))
             endif
          enddo

          ! Print one compact LOCAL reason summary every 100 model steps.
          ! Each repaired column contributes to exactly one reason category.
          if (mod(itimestep,1000) == 0) then
             if (pbl_pressure_window_local + pbl_repair_window_local + &
                 nqv_limited_window + nqc_limited_window + &
                 nqi_limited_window > 0) then
                write(0,'(A,I7,16(A,I7))') &
                     'TKE-EDMF REPAIR REASONS step=', itimestep, &
                     ' press=', pbl_pressure_window_local, &
                     ' total=', pbl_repair_window_local, &
                     ' nan=', nbad_nan_window, &
                     ' hpbl=', nbad_hpbl_window, &
                     ' t=', nbad_t_window, &
                     ' qvbad=', nbad_qv_window, &
                     ' qvlim=', nqv_limited_window, &
                     ' qvlim_total=', nqv_limited_total, &
                     ' qclim=', nqc_limited_window, &
                     ' qclim_total=', nqc_limited_total, &
                     ' qilim=', nqi_limited_window, &
                     ' qilim_total=', nqi_limited_total, &
                     ' qcbad=', nbad_qc_window, &
                     ' qibad=', nbad_qi_window, &
                     ' uv=', nbad_uv_window, &
                     ' tke=', nbad_tke_window
             endif

             pbl_pressure_window_local = 0
             pbl_repair_window_local   = 0
             nbad_nan_window  = 0
             nbad_hpbl_window = 0
             nbad_t_window    = 0
             nbad_qv_window   = 0
             nqv_limited_window = 0
             nqc_limited_window = 0
             nqi_limited_window = 0
             nbad_qc_window   = 0
             nbad_qi_window   = 0
             nbad_uv_window   = 0
             nbad_tke_window  = 0
          endif
       endif

       if (allocated(pbl_col_valid))   deallocate(pbl_col_valid)
       if (allocated(pbl_col_bad))     deallocate(pbl_col_bad)
       if (allocated(pbl_press_valid)) deallocate(pbl_press_valid)
       if (allocated(pbl_press_bad))   deallocate(pbl_press_bad)
       if (allocated(p_mid_tkeedmf))   deallocate(p_mid_tkeedmf)
       if (allocated(p_int_tkeedmf))   deallocate(p_int_tkeedmf)

       if (tke_edmf_errflg /= 0) then
          write(*,*) 'FATAL in bl_tke_edmf: ', trim(tke_edmf_errmsg)
          call mpas_timer_stop('bl_tke_edmf')
          call mpas_timer_stop('PBL Scheme')
          stop 1
       endif

       call mpas_timer_stop('bl_tke_edmf')

    case("bl_mynn")
       call mpas_pool_get_config(configs,'config_mynn_cloudpdf'   ,bl_mynn_cloudpdf   )
       call mpas_pool_get_config(configs,'config_mynn_mixlength'  ,bl_mynn_mixlength  )
       call mpas_pool_get_config(configs,'config_mynn_stfunc'     ,bl_mynn_stfunc     )
       call mpas_pool_get_config(configs,'config_mynn_topdown'    ,bl_mynn_topdown    )
       call mpas_pool_get_config(configs,'config_mynn_scaleaware' ,bl_mynn_scaleaware )
       call mpas_pool_get_config(configs,'config_mynn_dheat_opt'  ,bl_mynn_dheat_opt  )
       call mpas_pool_get_config(configs,'config_mynn_edmf'       ,bl_mynn_edmf       )
       call mpas_pool_get_config(configs,'config_mynn_edmf_dd'    ,bl_mynn_edmf_dd    )
       call mpas_pool_get_config(configs,'config_mynn_edmf_mom'   ,bl_mynn_edmf_mom   )
       call mpas_pool_get_config(configs,'config_mynn_edmf_tke'   ,bl_mynn_edmf_tke   )
       call mpas_pool_get_config(configs,'config_mynn_edmf_output',bl_mynn_edmf_output)
       call mpas_pool_get_config(configs,'config_mynn_closure'    ,bl_mynn_closure    )
       call mpas_pool_get_config(configs,'config_mynn_mixscalars' ,bl_mynn_mixscalars )
       call mpas_pool_get_config(configs,'config_mynn_mixclouds'  ,bl_mynn_cloudmix   )
       call mpas_pool_get_config(configs,'config_mynn_mixqt'      ,bl_mynn_mixqt      )
       call mpas_pool_get_config(configs,'config_mynn_tkeadvect'  ,bl_mynn_tkeadvect  )
       call mpas_pool_get_config(configs,'config_mynn_tkebudget'  ,bl_mynn_tkebudget  )
       call mpas_pool_get_config(configs,'config_icloud_bl'       ,icloud_bl          )
       call mpas_pool_get_config(configs,'config_spp_pbl'         ,spp_pbl            )
!      call mpas_log_write(' ')
!      call mpas_log_write('--- enter subroutine mynn_bl_driver:')
!      call mpas_log_write('--- config_mynn_cloudpdf    = $i',intArgs=(/bl_mynn_cloudpdf/))
!      call mpas_log_write('--- config_mynn_mixlength   = $i',intArgs=(/bl_mynn_mixlength/))
!      call mpas_log_write('--- config_mynn_stfunc      = $i',intArgs=(/bl_mynn_stfunc/))
!      call mpas_log_write('--- config_mynn_topdown     = $i',intArgs=(/bl_mynn_topdown/))
!      call mpas_log_write('--- config_mynn_scaleaware  = $i',intArgs=(/bl_mynn_scaleaware/))
!      call mpas_log_write('--- config_mynn_dheat_opt   = $i',intArgs=(/bl_mynn_dheat_opt/))
!      call mpas_log_write('--- config_mynn_edmf        = $i',intArgs=(/bl_mynn_edmf/))
!      call mpas_log_write('--- config_mynn_edmf_dd     = $i',intArgs=(/bl_mynn_edmf_dd/))
!      call mpas_log_write('--- config_mynn_edmf_mom    = $i',intArgs=(/bl_mynn_edmf_mom/))
!      call mpas_log_write('--- config_mynn_edmf_tke    = $i',intArgs=(/bl_mynn_edmf_tke/))
!      call mpas_log_write('--- config_mynn_edmf_output = $i',intArgs=(/bl_mynn_edmf_output/))
!      call mpas_log_write('--- config_mynn_mixscalars  = $i',intArgs=(/bl_mynn_mixscalars/))
!      call mpas_log_write('--- config_mynn_mixclouds   = $i',intArgs=(/bl_mynn_cloudmix/))
!      call mpas_log_write('--- config_mynn_mixqt       = $i',intArgs=(/bl_mynn_mixqt/))
!      call mpas_log_write('--- config_mynn_tkeadvect   = $l',logicArgs=(/bl_mynn_tkeadvect/))
!      call mpas_log_write('--- config_mynn_tkebudget   = $i',intArgs=(/bl_mynn_tkebudget/))
!      call mpas_log_write('--- config_mynn_closure     = $r',realArgs=(/bl_mynn_closure/))
!      call mpas_log_write(' ')
!      call mpas_log_write('--- f_qc   = $l',logicArgs=(/f_qc/)  )
!      call mpas_log_write('--- f_qi   = $l',logicArgs=(/f_qi/)  )
!      call mpas_log_write('--- f_qs   = $l',logicArgs=(/f_qs/)  )
!      call mpas_log_write('--- f_qoz  = $l',logicArgs=(/f_qoz/) )
!      call mpas_log_write('--- f_nc   = $l',logicArgs=(/f_nc/)  )
!      call mpas_log_write('--- f_ni   = $l',logicArgs=(/f_ni/)  )
!      call mpas_log_write('--- f_nifa = $l',logicArgs=(/f_nifa/))
!      call mpas_log_write('--- f_nwfa = $l',logicArgs=(/f_nwfa/))
!      call mpas_log_write('--- f_nbca = $l',logicArgs=(/f_nbca/))

       call mpas_timer_start('bl_mynn')
       call mynn_bl_driver( &
                  f_qc       = f_qc         , f_qi       = f_qi         , f_qs        = f_qs            , &
                  f_qoz      = f_qoz        , f_nc       = f_nc         , f_ni        = f_ni            , &
                  f_nifa     = f_nifa       , f_nwfa     = f_nwfa       , f_nbca      = f_nbca          , &
                  icloud_bl  = icloud_bl    , delt       = dt_pbl       , dx          = dx_p            , &
                  xland      = xland_p      , ps         = psfc_p       , ts          = tsk_p           , &
                  qsfc       = qsfc_p       , ust        = ust_p        , ch          = ch_p            , &
                  hfx        = hfx_p        , qfx        = qfx_p        , rmol        = rmol_p          , &
                  wspd       = wspd_p       , znt        = znt_p        , uoce        = uoce_p          , &
                  voce       = voce_p       , dz         = dz_p         , u           = u_p             , &
                  v          = v_p          , w          = w_p          , th          = th_p            , &
                  tt         = t_p          , p          = pres_hyd_p   , exner       = pi_p            , &
                  rho        = rho_p        , qv         = qv_p         , qc          = qc_p            , &
                  qi         = qi_p         , qs         = qs_p         , nc          = nc_p            , &
                  ni         = ni_p         , nifa       = nifa_p       , nwfa        = nwfa_p          , &
                  rthraten   = rthraten_p   , pblh       = hpbl_p       , kpbl        = kpbl_p          , &
                  cldfra_bl  = cldfrabl_p   , qc_bl      = qcbl_p       , qi_bl       = qibl_p          , &
                  maxwidth   = maxwidthbl_p , maxmf      = maxmfbl_p    , ktop_plume  = kbl_plume_p     , &
                  ztop_plume = zbl_plume_p  , dqke       = dqke_p       , qke_adv     = qkeadv_p        , &
                  tsq        = tsq_p        , qsq        = qsq_p        , cov         = cov_p           , &
                  el_pbl     = elpbl_p      , rublten    = rublten_p    , rvblten     = rvblten_p       , &
                  rthblten   = rthblten_p   , rqvblten   = rqvblten_p   , rqcblten    = rqcblten_p      , &
                  rqiblten   = rqiblten_p   , rqsblten   = rqsblten_p   , rncblten    = rncblten_p      , &
                  rniblten   = rniblten_p   , rnifablten = rnifablten_p , rnwfablten  = rnwfablten_p    , &
                  edmf_a     = edmfa_p      , edmf_w     = edmfw_p      , edmf_qt     = edmfqt_p        , &
                  edmf_thl   = edmfthl_p    , edmf_ent   = edmfent_p    , edmf_qc     = edmfqc_p        , &
                  sub_thl    = subthl_p     , sub_sqv    = subqv_p      , det_thl     = detthl_p        , &
                  det_sqv    = detqv_p      , exch_h     = kzh_p        , exch_m      = kzm_p           , &
                  qke        = qke_p        , qwt        = qwt_p        , qshear      = qshear_p        , &
                  qbuoy      = qbuoy_p      , qdiss      = qdiss_p      , sh3d        = sh3d_p          , &
                  sm3d       = sm3d_p       , spp_pbl    = spp_pbl      , pattern_spp = pattern_spp_pbl , &
                  do_restart         = config_do_restart   ,                                              &
                  do_DAcycling       = config_do_DAcycling ,                                              &
                  initflag           = initflag            ,                                              &
                  bl_mynn_tkeadvect  = bl_mynn_tkeadvect   ,                                              &
                  bl_mynn_tkebudget  = bl_mynn_tkebudget   ,                                              &
                  bl_mynn_cloudpdf   = bl_mynn_cloudpdf    ,                                              &
                  bl_mynn_mixlength  = bl_mynn_mixlength   ,                                              &
                  bl_mynn_closure    = bl_mynn_closure     ,                                              &
                  bl_mynn_stfunc     = bl_mynn_stfunc      ,                                              &
                  bl_mynn_topdown    = bl_mynn_topdown     ,                                              &
                  bl_mynn_scaleaware = bl_mynn_scaleaware  ,                                              &
                  bl_mynn_dheat_opt  = bl_mynn_dheat_opt   ,                                              &
                  bl_mynn_edmf       = bl_mynn_edmf        ,                                              &
                  bl_mynn_edmf_dd    = bl_mynn_edmf_dd     ,                                              &
                  bl_mynn_edmf_mom   = bl_mynn_edmf_mom    ,                                              &
                  bl_mynn_edmf_tke   = bl_mynn_edmf_tke    ,                                              &
                  bl_mynn_output     = bl_mynn_edmf_output ,                                              &
                  bl_mynn_mixscalars = bl_mynn_mixscalars  ,                                              &
                  bl_mynn_cloudmix   = bl_mynn_cloudmix    ,                                              &
                  bl_mynn_mixqt      = bl_mynn_mixqt       ,                                              &
                  ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,                 &
                  ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme ,                 &
                  its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte ,                 &
                  errmsg = errmsg , errflg = errflg                                                       &
                          )
      call mpas_timer_stop('bl_mynn')
!     call mpas_log_write('--- exit subroutine mynn_bl_driver:')
!     call mpas_log_write(' ')
       
    case("bl_mynnedmf")
       call mpas_pool_get_config(configs,'config_mynn_cloudpdf'   ,bl_mynn_cloudpdf   )
       call mpas_pool_get_config(configs,'config_mynn_mixlength'  ,bl_mynn_mixlength  )
       call mpas_pool_get_config(configs,'config_mynn_stfunc'     ,bl_mynn_stfunc     )
       call mpas_pool_get_config(configs,'config_mynn_topdown'    ,bl_mynn_topdown    )
       call mpas_pool_get_config(configs,'config_mynn_scaleaware' ,bl_mynn_scaleaware )
       call mpas_pool_get_config(configs,'config_mynn_dheat_opt'  ,bl_mynn_dheat_opt  )
       call mpas_pool_get_config(configs,'config_mynn_edmf'       ,bl_mynn_edmf       )
       call mpas_pool_get_config(configs,'config_mynn_edmf_dd'    ,bl_mynn_edmf_dd    )
       call mpas_pool_get_config(configs,'config_mynn_edmf_mom'   ,bl_mynn_edmf_mom   )
       call mpas_pool_get_config(configs,'config_mynn_edmf_tke'   ,bl_mynn_edmf_tke   )
       call mpas_pool_get_config(configs,'config_mynn_edmf_output',bl_mynn_edmf_output)
       call mpas_pool_get_config(configs,'config_mynn_closure'    ,bl_mynn_closure    )
       call mpas_pool_get_config(configs,'config_mynn_mixscalars' ,bl_mynn_mixscalars )
       call mpas_pool_get_config(configs,'config_mynn_mixaerosols',bl_mynn_mixaerosols)
       call mpas_pool_get_config(configs,'config_mynn_mixnumcon'  ,bl_mynn_mixnumcon  )
       call mpas_pool_get_config(configs,'config_mynn_mixclouds'  ,bl_mynn_cloudmix   )
       call mpas_pool_get_config(configs,'config_mynn_mixqt'      ,bl_mynn_mixqt      )
       call mpas_pool_get_config(configs,'config_mynn_ess'        ,bl_mynn_ess        )
       call mpas_pool_get_config(configs,'config_mynn_tkeadvect'  ,bl_mynn_tkeadvect  )
       call mpas_pool_get_config(configs,'config_mynn_tkebudget'  ,bl_mynn_tkebudget  )
       call mpas_pool_get_config(configs,'config_icloud_bl'       ,icloud_bl          )
       call mpas_pool_get_config(configs,'config_spp_pbl'         ,spp_pbl            )
       call mpas_pool_get_config(configs,'config_mynn_enh_mix'    ,enh_mix            )

!      call mpas_log_write(' ')
!      call mpas_log_write('--- enter subroutine mynn_bl_driver:')
!      call mpas_log_write('--- config_mynn_cloudpdf    = $i',intArgs=(/bl_mynn_cloudpdf/))
!      call mpas_log_write('--- config_mynn_mixlength   = $i',intArgs=(/bl_mynn_mixlength/))
!      call mpas_log_write('--- config_mynn_stfunc      = $i',intArgs=(/bl_mynn_stfunc/))
!      call mpas_log_write('--- config_mynn_topdown     = $i',intArgs=(/bl_mynn_topdown/))
!      call mpas_log_write('--- config_mynn_scaleaware  = $i',intArgs=(/bl_mynn_scaleaware/))
!      call mpas_log_write('--- config_mynn_dheat_opt   = $i',intArgs=(/bl_mynn_dheat_opt/))
!      call mpas_log_write('--- config_mynn_edmf        = $i',intArgs=(/bl_mynn_edmf/))
!      call mpas_log_write('--- config_mynn_edmf_dd     = $i',intArgs=(/bl_mynn_edmf_dd/))
!      call mpas_log_write('--- config_mynn_edmf_mom    = $i',intArgs=(/bl_mynn_edmf_mom/))
!      call mpas_log_write('--- config_mynn_edmf_tke    = $i',intArgs=(/bl_mynn_edmf_tke/))
!      call mpas_log_write('--- config_mynn_edmf_output = $i',intArgs=(/bl_mynn_edmf_output/))
!      call mpas_log_write('--- config_mynn_mixscalars  = $i',intArgs=(/bl_mynn_mixscalars/))
!      call mpas_log_write('--- config_mynn_mixclouds   = $i',intArgs=(/bl_mynn_cloudmix/))
!      call mpas_log_write('--- config_mynn_mixqt       = $i',intArgs=(/bl_mynn_mixqt/))
!      call mpas_log_write('--- config_mynn_tkeadvect   = $l',logicArgs=(/bl_mynn_tkeadvect/))
!      call mpas_log_write('--- config_mynn_tkebudget   = $i',intArgs=(/bl_mynn_tkebudget/))
!      call mpas_log_write('--- config_mynn_closure     = $r',realArgs=(/bl_mynn_closure/))
!      call mpas_log_write(' ')
!      call mpas_log_write('--- f_qc   = $l',logicArgs=(/f_qc/)  )
!      call mpas_log_write('--- f_qi   = $l',logicArgs=(/f_qi/)  )
!      call mpas_log_write('--- f_qs   = $l',logicArgs=(/f_qs/)  )
!      call mpas_log_write('--- f_qoz  = $l',logicArgs=(/f_qoz/) )
!      call mpas_log_write('--- f_nc   = $l',logicArgs=(/f_nc/)  )
!      call mpas_log_write('--- f_ni   = $l',logicArgs=(/f_ni/)  )
!      call mpas_log_write('--- f_nifa = $l',logicArgs=(/f_nifa/))
!      call mpas_log_write('--- f_nwfa = $l',logicArgs=(/f_nwfa/))
!      call mpas_log_write('--- f_nbca = $l',logicArgs=(/f_nbca/))

       if (do_chemistry) mix_chem = .true.

       call mpas_timer_start('bl_mynnedmf')
       call mynnedmf_driver( &
                  f_qc       = f_qc         , f_qi       = f_qi        , f_qs        = f_qs             , &
                  f_qoz      = f_qoz        , f_nc       = f_nc        , f_ni        = f_ni             , &
                  f_nifa     = f_nifa       , f_nwfa     = f_nwfa      , f_nbca      = f_nbca           , &
                  icloud_bl  = icloud_bl    , delt       = dt_pbl      , dx          = dx_p             , &
                  xland      = xland_p      , ps         = psfc_p      , ts          = tsk_p            , &
                  qsfc       = qsfc_p       , ust        = ust_p       , ch          = ch_p             , &
                  hfx        = hfx_p        , qfx        = qfx_p       ,                                  &
                  wspd       = wspd_p       , znt        = znt_p       , uoce        = uoce_p           , &
                  voce       = voce_p       , dz         = dz_p        , u           = u_p              , &
                  v          = v_p          , w          = w_p         , th          = th_p             , &
                  tt         = t_p          , p          = pres_p      , exner       = pi_p             , &
                  rho        = rho_p        , qv         = qv_p        , qc          = qc_p             , &
                  qi         = qi_p         , qs         = qs_p        , ni          = ni_p             , &
                  nc         = nc_p         , nifa       = nifa_p      , nwfa        = nwfa_p           , &
                  rthraten   = rthraten_p   , pblh       = hpbl_p      , kpbl        = kpbl_p           , &
                  cldfra_bl  = cldfrabl_p   , qc_bl      = qcbl_p      , qi_bl       = qibl_p           , &
                  maxwidth   = maxwidthbl_p , maxmf      = maxmfbl_p   , ztop_plume  = zbl_plume_p      , &
                  excess_h   = excessh_p    , excess_q   = excessq_p   ,                                  &
                  qke        = qke_p        , qke_adv    = qkeadv_p    ,                                  &
                  tsq        = tsq_p        , qsq        = qsq_p       , cov         = cov_p            , &
                  el_pbl     = elpbl_p      , rublten    = rublten_p   , rvblten     = rvblten_p        , &
                  rthblten   = rthblten_p   , rqvblten   = rqvblten_p  , rqcblten    = rqcblten_p       , &
                  rqiblten   = rqiblten_p   , rqsblten   = rqsblten_p  , rniblten    = rniblten_p       , &
                  rncblten   = rncblten_p   , rnwfablten = rnwfablten_p, rnifablten  = rnifablten_p     , &
                  edmf_a     = edmfa_p      , edmf_w     = edmfw_p     , edmf_qt     = edmfqt_p         , &
                  edmf_thl   = edmfthl_p    , edmf_ent   = edmfent_p   , edmf_qc     = edmfqc_p         , &
                  sub_thl    = subthl_p     , sub_sqv    = subqv_p     , det_thl     = detthl_p         , &
                  det_sqv    = detqv_p      , exch_h     = kzh_p       , exch_m      = kzm_p            , &
                  dqke       = dqke_p       , qwt        = qwt_p       , qshear      = qshear_p         , &
                  qbuoy      = qbuoy_p      , qdiss      = qdiss_p     , sh3d        = sh3d_p           , &
                  sm3d       = sm3d_p       , spp_pbl    = spp_pbl     , pattern_spp = pattern_spp_pbl  , &
                  do_restart         = config_do_restart   ,                                              &
                  do_DAcycling       = config_do_DAcycling ,                                              &
                  initflag           = initflag            ,                                              &
                  bl_mynn_tkeadvect  = bl_mynn_tkeadvect   ,                                              &
                  bl_mynn_tkebudget  = bl_mynn_tkebudget   ,                                              &
                  bl_mynn_cloudpdf   = bl_mynn_cloudpdf    ,                                              &
                  bl_mynn_mixlength  = bl_mynn_mixlength   ,                                              &
                  bl_mynn_closure    = bl_mynn_closure     ,                                              &
                  bl_mynn_stfunc     = bl_mynn_stfunc      ,                                              &
                  bl_mynn_topdown    = bl_mynn_topdown     ,                                              &
                  bl_mynn_scaleaware = bl_mynn_scaleaware  ,                                              &
                  bl_mynn_dheat_opt  = bl_mynn_dheat_opt   ,                                              &
                  bl_mynn_edmf       = bl_mynn_edmf        ,                                              &
                  bl_mynn_edmf_dd    = bl_mynn_edmf_dd     ,                                              &
                  bl_mynn_edmf_mom   = bl_mynn_edmf_mom    ,                                              &
                  bl_mynn_edmf_tke   = bl_mynn_edmf_tke    ,                                              &
                  bl_mynn_output     = bl_mynn_edmf_output ,                                              &
                  bl_mynn_mixscalars = bl_mynn_mixscalars  ,                                              &
                  bl_mynn_mixaerosols= bl_mynn_mixaerosols ,                                              &
                  bl_mynn_mixnumcon  = bl_mynn_mixnumcon   ,                                              &
                  bl_mynn_cloudmix   = bl_mynn_cloudmix    ,                                              &
                  bl_mynn_mixqt      = bl_mynn_mixqt       ,                                              &
                  bl_mynn_ess        = bl_mynn_ess         ,                                              &
                  mix_chem           = mix_chem            ,                                              &
                  nchem              = num_chem            ,                                              &
                  ndvel              = num_chem            ,                                              &
                  chem3d             = chem_p              ,                                              &
                  vd3d               = ddvel_p             ,                                              &
                  settle3d           = tend_chem_settle_p  ,                                              &
                  enh_mix            = enh_mix             ,                                              &
                  frp_mean           = frp_out_p           ,                                              &
                  emis_ant_no        = aero_emis_for_enhmix_p ,                                           &
                  ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,                 &
                  ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme ,                 &
                  its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte ,                 &
                  errmsg = errmsg , errflg = errflg                                                       &
                          )
       call mpas_timer_stop('bl_mynnedmf')
!      call mpas_log_write('--- exit subroutine mynnedmf_driver:')
!      call mpas_log_write(' ')

    case("bl_myj")
       if (initflag==1) then
          do i = its,ite
          do k = kts,kte
          do j = jts,jte
             tkepbl_p(i,k,j) = epsq2
          enddo
          enddo
          enddo
       endif
       call mpas_timer_start('MYJ_pbl')
       call mpas_log_write('dt_pbl = $r',realArgs=(/dt_pbl/))
!       call mpas_log_write("BEFORE MYJPBL")
! call mpas_dmpar_min_real(domain % dminfo, minval(ht_p), scalar_min)
! call mpas_dmpar_min_real(domain % dminfo, maxval(ht_p), scalar_max)
!call mpas_log_write('MIN/MAX ht_p = $r, $r', realArgs=(/scalar_min,scalar_max/)
!
! call mpas_dmpar_min_real(domain % dminfo, minval(qsfc_p), scalar_min)
! call mpas_dmpar_min_real(domain % dminfo, maxval(qsfc_p), scalar_max)
!call mpas_log_write('MIN/MAX qsfc_p = $r, $r', realArgs=(/scalar_min,scalar_max/)
!
! call mpas_dmpar_min_real(domain % dminfo, minval(thz0_p), scalar_min)
! call mpas_dmpar_min_real(domain % dminfo, maxval(thz0_p), scalar_max)
!call mpas_log_write('MIN/MAX thz0_p = $r, $r',realArgs=(/scalar_min,scalar_max/)
!
! call mpas_dmpar_min_real(domain % dminfo, minval(qz0_p), scalar_min)
! call mpas_dmpar_min_real(domain % dminfo, maxval(qz0_p), scalar_max)
!call mpas_log_write('MIN/MAX qz0_p = $r, $r',realArgs=(/scalar_min,scalar_max/)
!
! call mpas_dmpar_min_real(domain % dminfo, minval(uz0_p), scalar_min)
! call mpas_dmpar_min_real(domain % dminfo, maxval(uz0_p), scalar_max)
!call mpas_log_write('MIN/MAX uz0_p = $r, $r', realArgs=(/scalar_min,scalar_max/)
!
! call mpas_dmpar_min_real(domain % dminfo, minval(vz0_p), scalar_min)
! call mpas_dmpar_min_real(domain % dminfo, maxval(vz0_p), scalar_max)
!call mpas_log_write('MIN/MAX vz0_p = $r, $r',realArgs=(/scalar_min,scalar_max/)
!
! call mpas_dmpar_min_real(domain % dminfo, minval(tkepbl_p), scalar_min)
! call mpas_dmpar_min_real(domain % dminfo, maxval(tkepbl_p), scalar_max)
!call mpas_log_write('MIN/MAX tkepbl_p = $r, $r',realArgs=(/scalar_min,scalar_max/)
!
! call mpas_dmpar_min_real(domain % dminfo, minval(kzh_p), scalar_min)
! call mpas_dmpar_min_real(domain % dminfo, maxval(kzh_p), scalar_max)
!call mpas_log_write('MIN/MAX kzh_p = $r, $r',realArgs=(/scalar_min,scalar_max/)
!
! call mpas_dmpar_min_real(domain % dminfo, minval(ust_p), scalar_min)
! call mpas_dmpar_min_real(domain % dminfo, maxval(ust_p), scalar_max)
!call mpas_log_write('MIN/MAX ust_p = $r, $r', realArgs=(/scalar_min,scalar_max/)
!
! call mpas_dmpar_min_real(domain % dminfo, minval(znt_p), scalar_min)
! call mpas_dmpar_min_real(domain % dminfo, maxval(znt_p), scalar_max)
!call mpas_log_write('MIN/MAX znt_p = $r, $r',realArgs=(/scalar_min,scalar_max/)
!
! call mpas_dmpar_min_real(domain % dminfo, minval(ct_p), scalar_min)
! call mpas_dmpar_min_real(domain % dminfo, maxval(ct_p), scalar_max)
!call mpas_log_write('MIN/MAX ct_p = $r, $r',realArgs=(/scalar_min,scalar_max/)
!
! call mpas_dmpar_min_real(domain % dminfo, minval(akhs_p), scalar_min)
! call mpas_dmpar_min_real(domain % dminfo, maxval(akhs_p), scalar_max)
!call mpas_log_write('MIN/MAX akhs_p = $r, $r', realArgs=(/scalar_min,scalar_max/)
!
! call mpas_dmpar_min_real(domain % dminfo, minval(akms_p), scalar_min)
! call mpas_dmpar_min_real(domain % dminfo, maxval(akms_p), scalar_max)
!call mpas_log_write('MIN/MAX akms_p = $r, $r', realArgs=(/scalar_min,scalar_max/)

       call myjpbl ( &
                dt       = dt_pbl    , stepbl   = 1          , ht       = ht_p       , &
                pmid     = pres_hyd_p, pint     = pres2_hyd_p, th       = th_p       , &
                t        = t_p       , exner    = pi_p       , qv       = qv_p       , &
                qcw      = qc_p      , qci      = qi_p       , qcs      = qs_p       , &
                qcr      = qr_p      , qcg      = qg_p       , u        = u_p        , &
                v        = v_p       , rho      = rho_p      , tsk      = tsk_p      , &
                qsfc     = qsfc_p    , chklowq  = chlowq_p   , thz0     = thz0_p     , &
                qz0      = qz0_p     , uz0      = uz0_p      , vz0      = vz0_p      , &
                lowlyr   = lowlyr_p  , xland    = xland_p    , sice     = xice_p     , &
                snow     = snow_p    , tke_myj  = tkepbl_p   , exch_h   = kzh_p      , &
                ustar    = ust_p     , znt      = znt_p      , el_myj   = elpbl_p    , &
                pblh     = hpbl_p    , kpbl     = kpbl_p     , ct       = ct_p       , &
                akhs     = akhs_p    , akms     = akms_p     , elflx    = lh_p       , &
                mixht    = mixht_p   , rublten  = rublten_p  , rvblten  = rvblten_p  , &
                rthblten = rthblten_p, rqvblten = rqvblten_p , rqcblten = rqcblten_p , &
                rqiblten = rqiblten_p, dz       = dz_p       ,                         &
                ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde, &
                ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme, &
                its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte  &
                            )
       call mpas_timer_stop('MYJ_pbl')
         

    case default

 end select pbl_select

!copy local arrays to MPAS grid:
 call pbl_to_MPAS(configs,state,diag_physics,tend_physics,its,ite)

!call mpas_log_write('--- end subroutine driver_pbl.')

 end subroutine driver_pbl

!=================================================================================================================
 end module mpas_atmphys_driver_pbl
!=================================================================================================================
