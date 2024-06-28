! Copyright (c) 2024 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_driver_seaice
 use mpas_kind_types
 use mpas_pool_routines,only: mpas_pool_get_array,mpas_pool_get_config,mpas_pool_type
 use mpas_log

 use mpas_atmphys_constants,only: rcp,gravity
 use mpas_atmphys_lsm_shared,only: correct_tsk_over_seaice
 use mpas_atmphys_vars
 use module_sf_noah_seaice_drv
 use module_sf_sfcdiags
 !RUC lsm ice module 
 use module_ruc_ice
 use module_sf_sfcdiags_ruclsm

 implicit none
 private
 public:: allocate_seaice,   &
          deallocate_seaice, &
          driver_seaice

 logical,parameter:: frpcpn   = .false.

 logical,parameter:: myjpbl   = .false.
 logical,parameter:: frpcpn_ruc  = .true.

!urban physics: MPAS does not plan to run the urban physics option.
 integer,parameter:: sf_urban_physics = 0 !activate urban canopy model (=0: no urban canopy)

 integer :: targetcell= 0 !99975 !81442 !107619

!MPAS driver for parameterization of surface processes over seaice points.
!Laura D. Fowler (laura@ucar.edu) / 2024-03-13.

!RUC ice model is added to the seaice driver, added variables for the call to RUC seaice model.
!Tanya Smirnova (Tanya.Smirnova@noaa.gov) / 2024-04-26

 contains


!=================================================================================================================
 subroutine allocate_seaice (config_lsm_scheme)
!=================================================================================================================

 character(len=StrKIND),pointer, intent(in) :: config_lsm_scheme

 if(.not.allocated(acsnom_p)    ) allocate(acsnom_p(ims:ime,jms:jme)    )
 if(.not.allocated(acsnow_p)    ) allocate(acsnow_p(ims:ime,jms:jme)    )
 if(.not.allocated(albsi_p)     ) allocate(albsi_p(ims:ime,jms:jme)     )
 if(.not.allocated(br_p)        ) allocate(br_p(ims:ime,jms:jme)        )
 if(.not.allocated(chs_p)       ) allocate(chs_p(ims:ime,jms:jme)       )
 if(.not.allocated(chs2_p)      ) allocate(chs2_p(ims:ime,jms:jme)      )
 if(.not.allocated(cpm_p)       ) allocate(cpm_p(ims:ime,jms:jme)       )
 if(.not.allocated(chklowq_p)   ) allocate(chklowq_p(ims:ime,jms:jme)   )
 if(.not.allocated(cqs2_p)      ) allocate(cqs2_p(ims:ime,jms:jme)      )
 if(.not.allocated(cqs_p)       ) allocate(cqs_p(ims:ime,jms:jme)       )
 if(.not.allocated(qgh_p)       ) allocate(qgh_p(ims:ime,jms:jme)       )
 if(.not.allocated(qsfc_p)      ) allocate(qsfc_p(ims:ime,jms:jme)      )
 if(.not.allocated(glw_p)       ) allocate(glw_p(ims:ime,jms:jme)       )
 if(.not.allocated(gsw_p)       ) allocate(gsw_p(ims:ime,jms:jme)       )
 if(.not.allocated(grdflx_p)    ) allocate(grdflx_p(ims:ime,jms:jme)    )
 if(.not.allocated(icedepth_p)  ) allocate(icedepth_p(ims:ime,jms:jme)  )
 if(.not.allocated(hfx_p)       ) allocate(hfx_p(ims:ime,jms:jme)       )
 if(.not.allocated(qfx_p)       ) allocate(qfx_p(ims:ime,jms:jme)       )
 if(.not.allocated(lh_p)        ) allocate(lh_p(ims:ime,jms:jme)        )
 if(.not.allocated(noahres_p)   ) allocate(noahres_p(ims:ime,jms:jme)   )
 if(.not.allocated(potevp_p)    ) allocate(potevp_p(ims:ime,jms:jme)    )
 if(.not.allocated(rainbl_p)    ) allocate(rainbl_p(ims:ime,jms:jme)    )
 if(.not.allocated(sfc_albbck_p)) allocate(sfc_albbck_p(ims:ime,jms:jme))
 if(.not.allocated(sfc_albedo_p)) allocate(sfc_albedo_p(ims:ime,jms:jme))
 if(.not.allocated(sfc_emiss_p) ) allocate(sfc_emiss_p(ims:ime,jms:jme) )
 if(.not.allocated(sfcrunoff_p) ) allocate(sfcrunoff_p(ims:ime,jms:jme) )
 if(.not.allocated(snoalb_p)    ) allocate(snoalb_p(ims:ime,jms:jme)    )
 if(.not.allocated(snow_p)      ) allocate(snow_p(ims:ime,jms:jme)      )
 if(.not.allocated(snowc_p)     ) allocate(snowc_p(ims:ime,jms:jme)     )
 if(.not.allocated(snowh_p)     ) allocate(snowh_p(ims:ime,jms:jme)     )
 if(.not.allocated(snopcx_p)    ) allocate(snopcx_p(ims:ime,jms:jme)    )
 if(.not.allocated(snowsi_p)    ) allocate(snowsi_p(ims:ime,jms:jme)    )
 if(.not.allocated(swdown_p)    ) allocate(swdown_p(ims:ime,jms:jme)    )
 if(.not.allocated(sr_p)        ) allocate(sr_p(ims:ime,jms:jme)        )
 if(.not.allocated(tsk_p)       ) allocate(tsk_p(ims:ime,jms:jme)       )
 if(.not.allocated(xice_p)      ) allocate(xice_p(ims:ime,jms:jme)      )
 if(.not.allocated(z0_p)        ) allocate(z0_p(ims:ime,jms:jme)        )
 if(.not.allocated(znt_p)       ) allocate(znt_p(ims:ime,jms:jme)       )
 if(.not.allocated(q2_p)        ) allocate(q2_p(ims:ime,jms:jme)        )
 if(.not.allocated(t2m_p)       ) allocate(t2m_p(ims:ime,jms:jme)       )
 if(.not.allocated(th2m_p)      ) allocate(th2m_p(ims:ime,jms:jme)      )

 if(.not.allocated(tsk_sea)     ) allocate(tsk_sea(ims:ime,jms:jme)     )
 if(.not.allocated(cpm_sea)     ) allocate(cpm_sea(ims:ime,jms:jme)     )
 if(.not.allocated(chs_sea)     ) allocate(chs_sea(ims:ime,jms:jme)     )
 if(.not.allocated(chs2_sea)    ) allocate(chs2_sea(ims:ime,jms:jme)    )
 if(.not.allocated(cqs2_sea)    ) allocate(cqs2_sea(ims:ime,jms:jme)    )
 if(.not.allocated(cqs_sea)     ) allocate(cqs_sea(ims:ime,jms:jme)     )
 if(.not.allocated(hfx_sea)     ) allocate(hfx_sea(ims:ime,jms:jme)     )
 if(.not.allocated(qfx_sea)     ) allocate(qfx_sea(ims:ime,jms:jme)     )
 if(.not.allocated(qgh_sea)     ) allocate(qgh_sea(ims:ime,jms:jme)     )
 if(.not.allocated(qsfc_sea)    ) allocate(qsfc_sea(ims:ime,jms:jme)    )
 if(.not.allocated(lh_sea)      ) allocate(lh_sea(ims:ime,jms:jme)      )
 if(.not.allocated(tsk_ice)     ) allocate(tsk_ice(ims:ime,jms:jme)     )
 if(.not.allocated(albsi_p)     ) allocate(albsi_p(ims:ime,jms:jme)     )
 if(.not.allocated(icedepth_p)  ) allocate(icedepth_p(ims:ime,jms:jme)  )
 if(.not.allocated(snowsi_p)    ) allocate(snowsi_p(ims:ime,jms:jme)    )

 if(.not.allocated(tslb_p)) allocate(tslb_p(ims:ime,1:num_soils,jms:jme))

  if(config_lsm_scheme=='sf_ruc') then
    if(.not.allocated(rhosnf_p)   )  allocate(rhosnf_p(ims:ime,jms:jme)    )
    if(.not.allocated(precipfr_p) )  allocate(precipfr_p(ims:ime,jms:jme)  )
    if(.not.allocated(qcg_p)      )  allocate(qcg_p(ims:ime,jms:jme)       )
    if(.not.allocated(qsg_p)      )  allocate(qsg_p(ims:ime,jms:jme)       )
    if(.not.allocated(qvg_p)      )  allocate(qvg_p(ims:ime,jms:jme)       )
    if(.not.allocated(dew_p)      )  allocate(dew_p(ims:ime,jms:jme)       )
    if(.not.allocated(soilt1_p)   )  allocate(soilt1_p(ims:ime,jms:jme)    )
    if(.not.allocated(tsnav_p)    )  allocate(tsnav_p(ims:ime,jms:jme)     )
    if(.not.allocated(acrunoff_p) )  allocate(acrunoff_p(ims:ime,jms:jme)  )
    if(.not.allocated(snowfallac_p)) allocate(snowfallac_p(ims:ime,jms:jme))
    if(.not.allocated(graupelncv_p)) allocate(graupelncv_p(ims:ime,jms:jme))
    if(.not.allocated(snowncv_p)  )  allocate(snowncv_p(ims:ime,jms:jme)   )
    if(.not.allocated(rainncv_p)  )  allocate(rainncv_p(ims:ime,jms:jme)   )
    if(.not.allocated(flqc_p)     )  allocate(flqc_p(ims:ime,jms:jme)      )
    if(.not.allocated(flhc_p)     )  allocate(flhc_p(ims:ime,jms:jme)      )
    if(.not.allocated(sfcevp_p)   )  allocate(sfcevp_p(ims:ime,jms:jme)    )
    if(.not.allocated(isltyp_p)   )  allocate(isltyp_p(ims:ime,jms:jme)    )
    if(.not.allocated(ivgtyp_p)   )  allocate(ivgtyp_p(ims:ime,jms:jme)    )
    if(.not.allocated(zs_p)       )  allocate(zs_p(1:num_soils)            )
    if(.not.allocated(tmn_p)      )  allocate(tmn_p(ims:ime,jms:jme)       )
    if(.not.allocated(xland_p)    )  allocate(xland_p(ims:ime,jms:jme)     )

 endif


 end subroutine allocate_seaice

!=================================================================================================================
 subroutine deallocate_seaice (config_lsm_scheme)
!=================================================================================================================

 character(len=StrKIND),pointer, intent(in) :: config_lsm_scheme

 if(allocated(acsnom_p)    ) deallocate(acsnom_p    )
 if(allocated(acsnow_p)    ) deallocate(acsnow_p    )
 if(allocated(albsi_p)     ) deallocate(albsi_p     )
 if(allocated(br_p)        ) deallocate(br_p        )
 if(allocated(cpm_p)       ) deallocate(cpm_p       )
 if(allocated(chs_p)       ) deallocate(chs_p       )
 if(allocated(chs2_p)      ) deallocate(chs2_p      )
 if(allocated(chklowq_p)   ) deallocate(chklowq_p   )
 if(allocated(cqs2_p)      ) deallocate(cqs2_p      )
 if(allocated(cqs_p)       ) deallocate(cqs_p       )
 if(allocated(qgh_p)       ) deallocate(qgh_p       )
 if(allocated(qsfc_p)      ) deallocate(qsfc_p      )
 if(allocated(glw_p)       ) deallocate(glw_p       )
 if(allocated(gsw_p)       ) deallocate(gsw_p       )
 if(allocated(grdflx_p)    ) deallocate(grdflx_p    )
 if(allocated(icedepth_p)  ) deallocate(icedepth_p  )
 if(allocated(isltyp_p)    ) deallocate(isltyp_p    )
 if(allocated(ivgtyp_p)    ) deallocate(ivgtyp_p    )
 if(allocated(hfx_p)       ) deallocate(hfx_p       )
 if(allocated(qfx_p)       ) deallocate(qfx_p       )
 if(allocated(lh_p)        ) deallocate(lh_p        )
 if(allocated(noahres_p)   ) deallocate(noahres_p   )
 if(allocated(potevp_p)    ) deallocate(potevp_p    )
 if(allocated(rainbl_p)    ) deallocate(rainbl_p    )
 if(allocated(sfc_albbck_p)) deallocate(sfc_albbck_p)
 if(allocated(sfc_albedo_p)) deallocate(sfc_albedo_p)
 if(allocated(sfc_emiss_p) ) deallocate(sfc_emiss_p )
 if(allocated(sfcrunoff_p) ) deallocate(sfcrunoff_p )
 if(allocated(snoalb_p)    ) deallocate(snoalb_p    )
 if(allocated(snow_p)      ) deallocate(snow_p      )
 if(allocated(snowc_p)     ) deallocate(snowc_p     )
 if(allocated(snowh_p)     ) deallocate(snowh_p     )
 if(allocated(snopcx_p)    ) deallocate(snopcx_p    )
 if(allocated(snowsi_p)    ) deallocate(snowsi_p    )
 if(allocated(swdown_p)    ) deallocate(swdown_p    )
 if(allocated(gsw_p)       ) deallocate(gsw_p       )
 if(allocated(sr_p)        ) deallocate(sr_p        )
 if(allocated(tsk_p)       ) deallocate(tsk_p       )
 if(allocated(xice_p)      ) deallocate(xice_p      )
 if(allocated(z0_p)        ) deallocate(z0_p        )
 if(allocated(znt_p)       ) deallocate(znt_p       )
 if(allocated(q2_p)        ) deallocate(q2_p        )
 if(allocated(t2m_p)       ) deallocate(t2m_p       )
 if(allocated(th2m_p)      ) deallocate(th2m_p      )

 if(allocated(cpm_sea)     ) deallocate(cpm_sea     )
 if(allocated(chs_sea)     ) deallocate(chs_sea     )
 if(allocated(chs2_sea)    ) deallocate(chs2_sea    )
 if(allocated(cqs2_sea)    ) deallocate(cqs2_sea    )
 if(allocated(cqs_sea)     ) deallocate(cqs_sea     )
 if(allocated(hfx_sea)     ) deallocate(hfx_sea     )
 if(allocated(qfx_sea)     ) deallocate(qfx_sea     )
 if(allocated(qgh_sea)     ) deallocate(qgh_sea     )
 if(allocated(qsfc_sea)    ) deallocate(qsfc_sea    )
 if(allocated(lh_sea)      ) deallocate(lh_sea      )
 if(allocated(tsk_sea)     ) deallocate(tsk_sea     )
 if(allocated(tsk_ice)     ) deallocate(tsk_ice     )
 if(allocated(albsi_p)     ) deallocate(albsi_p     )
 if(allocated(icedepth_p)  ) deallocate(icedepth_p  )
 if(allocated(snowsi_p)    ) deallocate(snowsi_p    )

 if(allocated(tslb_p)) deallocate(tslb_p)

 if (config_lsm_scheme == "sf_ruc") then
    if(allocated(rhosnf_p)    ) deallocate(rhosnf_p    )
    if(allocated(precipfr_p)  ) deallocate(precipfr_p  )
    if(allocated(qsg_p)       ) deallocate(qsg_p       )
    if(allocated(qvg_p)       ) deallocate(qvg_p       )
    if(allocated(qcg_p)       ) deallocate(qcg_p       )
    if(allocated(dew_p)       ) deallocate(dew_p       )
    if(allocated(soilt1_p)    ) deallocate(soilt1_p    )
    if(allocated(tsnav_p)     ) deallocate(tsnav_p     )
    if(allocated(acrunoff_p)  ) deallocate(acrunoff_p  )
    if(allocated(graupelncv_p)) deallocate(graupelncv_p)
    if(allocated(snowncv_p)   ) deallocate(snowncv_p   )
    if(allocated(rainncv_p)   ) deallocate(rainncv_p   )
    if(allocated(flqc_p)      ) deallocate(flqc_p      )
    if(allocated(flhc_p)      ) deallocate(flhc_p      )
    if(allocated(sfcevp_p)    ) deallocate(sfcevp_p    )
    if(allocated(zs_p)        ) deallocate(zs_p        )
    if(allocated(tmn_p)       ) deallocate(tmn_p       )
    if(allocated(xland_p)     ) deallocate(xland_p     )

 endif


 end subroutine deallocate_seaice

!=================================================================================================================
 subroutine seaice_from_MPAS(configs,mesh,diag_physics,sfc_input,its,ite)
!=================================================================================================================

!input and inout arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: sfc_input
 integer,intent(in):: its,ite

!local pointers:
 integer,pointer:: globalCells(:)
 character(len=StrKIND),pointer:: microp_scheme,    &
                                  convection_scheme, &
                                  lsm_scheme

 integer,dimension(:),pointer:: isltyp,ivgtyp

 real(kind=RKIND),dimension(:),pointer  :: acsnom,acsnow,cpm,chs,chs2,chklowq,cqs,cqs2,glw,     &
                                           grdflx,gsw,hfx,lh,noahres,potevp,qfx,qgh,qsfc,       &
                                           br,sfc_albedo,sfc_emiss,sfcrunoff,snopcx,sr,z0,znt
  real(kind=RKIND),dimension(:),pointer  :: rhosnf,precipfr,qsg,qvg,dew,soilt1,tsnav,acrunoff,  &
                                           snowfallac,graupelncv,snowncv,qcg,                   &
                                           flqc,flhc,sfcevp,mavail
 real(kind=RKIND),dimension(:),pointer  :: snoalb,sfc_albbck,snow,snowc,snowh,tmn,skintemp,xice,xland
 real(kind=RKIND),dimension(:),pointer  :: t2m,th2m,q2
 real(kind=RKIND),dimension(:),pointer  :: raincv,rainncv
 real(kind=RKIND),dimension(:,:),pointer:: tslb,dzs,zs

!local variables and arrays:
 integer:: i,j,n,nconfig_frac_seaice


!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('--- enter subroutine seaice_from_MPAS:')

 call mpas_pool_get_config(configs,'config_convection_scheme',convection_scheme)
 call mpas_pool_get_config(configs,'config_microp_scheme'    ,microp_scheme    )
 call mpas_pool_get_config(configs,'config_lsm_scheme'       ,lsm_scheme       )

 call mpas_pool_get_array(mesh,'indexToCellID',globalCells)

 call mpas_pool_get_array(diag_physics,'acsnom'           ,acsnom           )
 call mpas_pool_get_array(diag_physics,'acsnow'           ,acsnow           )
 call mpas_pool_get_array(diag_physics,'br'               ,br               )
 call mpas_pool_get_array(diag_physics,'chs'              ,chs              )
 call mpas_pool_get_array(diag_physics,'chs2'             ,chs2             )
 call mpas_pool_get_array(diag_physics,'cpm'              ,cpm              )
 call mpas_pool_get_array(diag_physics,'chklowq'          ,chklowq          )
 call mpas_pool_get_array(diag_physics,'cqs2'             ,cqs2             )
 call mpas_pool_get_array(diag_physics,'cqs'              ,cqs              )
 call mpas_pool_get_array(diag_physics,'qgh'              ,qgh              )
 call mpas_pool_get_array(diag_physics,'qsfc'             ,qsfc             )
 call mpas_pool_get_array(diag_physics,'glw'              ,glw              )
 call mpas_pool_get_array(diag_physics,'grdflx'           ,grdflx           )
 call mpas_pool_get_array(diag_physics,'gsw'              ,gsw              )
 call mpas_pool_get_array(diag_physics,'hfx'              ,hfx              )
 call mpas_pool_get_array(diag_physics,'qfx'              ,qfx              )
 call mpas_pool_get_array(diag_physics,'lh'               ,lh               )
 call mpas_pool_get_array(diag_physics,'noahres'          ,noahres          )
 call mpas_pool_get_array(diag_physics,'potevp'           ,potevp           )
 call mpas_pool_get_array(diag_physics,'sfc_albedo'       ,sfc_albedo       )
 call mpas_pool_get_array(diag_physics,'sfc_emiss'        ,sfc_emiss        )
 call mpas_pool_get_array(diag_physics,'sfcrunoff'        ,sfcrunoff        )
 call mpas_pool_get_array(diag_physics,'snopcx'           ,snopcx           )
 call mpas_pool_get_array(diag_physics,'z0'               ,z0               )
 call mpas_pool_get_array(diag_physics,'znt'              ,znt              )
 
 call mpas_pool_get_array(sfc_input,'snoalb'  ,snoalb  )
 call mpas_pool_get_array(sfc_input,'snow'    ,snow    )
 call mpas_pool_get_array(sfc_input,'snowc'   ,snowc   )
 call mpas_pool_get_array(sfc_input,'snowh'   ,snowh   )
 call mpas_pool_get_array(sfc_input,'skintemp',skintemp)
 call mpas_pool_get_array(sfc_input,'tslb'    ,tslb    )
 call mpas_pool_get_array(sfc_input,'xice'    ,xice    )

 if (lsm_scheme .eq. "sf_ruc") then
   call mpas_pool_get_array(diag_physics,'rhosnf'         ,rhosnf           )
   call mpas_pool_get_array(diag_physics,'precipfr'       ,precipfr         )
   call mpas_pool_get_array(diag_physics,'qsg'            ,qsg              )
   call mpas_pool_get_array(diag_physics,'qvg'            ,qvg              )
   call mpas_pool_get_array(diag_physics,'qcg'            ,qcg              )
   call mpas_pool_get_array(diag_physics,'dew'            ,dew              )
   call mpas_pool_get_array(diag_physics,'soilt1'         ,soilt1           )
   call mpas_pool_get_array(diag_physics,'tsnav'          ,tsnav            )
   call mpas_pool_get_array(diag_physics,'acrunoff'       ,acrunoff         )
   call mpas_pool_get_array(diag_physics,'snowfallac'     ,snowfallac       )
   call mpas_pool_get_array(diag_physics,'graupelncv'     ,graupelncv       )
   call mpas_pool_get_array(diag_physics,'snowncv'        ,snowncv          )
   call mpas_pool_get_array(diag_physics,'rainncv'        ,rainncv          )
   call mpas_pool_get_array(diag_physics,'flqc'           ,flqc             )
   call mpas_pool_get_array(diag_physics,'flhc'           ,flhc             )
   call mpas_pool_get_array(diag_physics,'sfcevp'         ,sfcevp           )
   call mpas_pool_get_array(diag_physics,'mavail'         ,mavail           )

   call mpas_pool_get_array(sfc_input,'zs'        ,zs        )
 endif

   call mpas_pool_get_array(sfc_input,'isltyp'    ,isltyp    )
   call mpas_pool_get_array(sfc_input,'ivgtyp'    ,ivgtyp    )
   call mpas_pool_get_array(sfc_input,'sfc_albbck',sfc_albbck)
   call mpas_pool_get_array(sfc_input,'tmn'       ,tmn       )
   call mpas_pool_get_array(sfc_input,'xland'     ,xland     )

 if (lsm_scheme == 'sf_ruc') then
    do n = 1,num_soils
       zs_p(n) = zs(n,its)
    enddo
 endif

 do j = jts,jte
    do i = its,ite
       !--- in variables:
       xice_p(i,j)      = xice(i)
       glw_p(i,j)       = glw(i)
       qgh_p(i,j)       = qgh(i)
       snoalb_p(i,j)    = snoalb(i)
       br_p(i,j)        = br(i)
       cpm_p(i,j)       = cpm(i)
       chs_p(i,j)       = chs(i)
       swdown_p(i,j)    = gsw(i)/(1._RKIND-sfc_albedo(i))
       gsw_p(i,j)       = gsw(i)
       if (lsm_scheme == 'sf_ruc') then
           chklowq_p(i,j)    = chklowq(i)
           isltyp_p(i,j)     = isltyp(i)
           ivgtyp_p(i,j)     = ivgtyp(i)
           sfc_emiss_p(i,j)  = sfc_emiss(i)
           sfc_albbck_p(i,j) = sfc_albbck(i)
           tmn_p(i,j)        = tmn(i)
           xland_p(i,j)      = xland(i)
           precipfr_p(i,j) = precipfr(i)
           graupelncv_p(i,j) = graupelncv(i)
           snowncv_p(i,j) = snowncv(i)
           rainncv_p(i,j) = rainncv(i)
           flqc_p(i,j) = flqc(i)
           flhc_p(i,j) = flhc(i)
       endif
       !--- inout variables:
       do n = 1,num_soils
          tslb_p(i,n,j) = tslb(n,i)
          !pattern_spp_lsm(i,n,j) = 0._RKIND
          !field_sf(i,n,j) = 0._RKIND
       enddo
       z0_p(i,j)        = z0(i)
       snow_p(i,j)      = snow(i)
       snowc_p(i,j)     = snowc(i)
       snowh_p(i,j)     = snowh(i)
       tsk_p(i,j)       = skintemp(i)
       cqs2_p(i,j)      = cqs2(i)
       cqs_p(i,j)       = cqs(i)
       acsnom_p(i,j)    = acsnom(i)
       acsnow_p(i,j)    = acsnow(i)
       sfcrunoff_p(i,j) = sfcrunoff(i)
       albsi_p(i,j)     = seaice_albedo_default
       snowsi_p(i,j)    = seaice_snowdepth_min
       icedepth_p(i,j)  = seaice_thickness_default
       !--- inout optional variables:
       potevp_p(i,j)    = potevp(i)
       snopcx_p(i,j)    = snopcx(i)
       if (lsm_scheme == 'sf_ruc') then
           rhosnf_p(i,j) = rhosnf(i)
           qcg_p(i,j) = qcg(i)
           qsg_p(i,j) = qsg(i)
           qvg_p(i,j) = qvg(i)
           soilt1_p(i,j) = soilt1(i)
           tsnav_p(i,j) = tsnav(i)
           snowfallac_p(i,j) = snowfallac(i)
           acrunoff_p(i,j) = acrunoff(i)
       endif
       !--- output variables:
       hfx_p(i,j)       = hfx(i)
       lh_p(i,j)        = lh(i)
       qfx_p(i,j)       = qfx(i)
       znt_p(i,j)       = znt(i)
       grdflx_p(i,j)    = grdflx(i)
       qsfc_p(i,j)      = qsfc(i)
       chs2_p(i,j)      = chs2(i)
       !--- output optional variables:
       noahres_p(i,j)   = noahres(i)
       if (lsm_scheme == 'sf_ruc') then
           dew_p(i,j) = dew(i)
           sfcevp_p(i,j) = sfcevp(i)
           if (globalcells(i)==targetcell) then
              print *,'from MPAS '
              print *,'globalcells(i),cqs_p,cqs2_p,chs_p,chs2_p',globalcells(i),cqs_p(i,j),cqs2_p(i,j),chs_p(i,j),chs2_p(i,j)
              print *,'flqc(i),mavail(i),rho_p(i,kts,j)',flqc(i),mavail(i),rho_p(i,kts,j)
           endif
       endif

       !modify the surface albedo and surface emissivity, and surface temperatures over sea-ice points:
       if(xice(i).ge.xice_threshold .and. xice(i).le.1._RKIND) then
          sfc_albedo_p(i,j) = (sfc_albedo(i) - 0.08_RKIND*(1._RKIND-xice(i))) / xice(i)
          sfc_emiss_p(i,j)  = (sfc_emiss(i) - 0.98_RKIND*(1._RKIND-xice(i))) / xice(i)
       else
          sfc_emiss_p(i,j)  = sfc_emiss(i)
          sfc_albedo_p(i,j) = sfc_albedo(i)
       endif
    enddo

    !calculate sea-surface and sea-ice temperatures over sea-ice grid cells:
    call correct_tsk_over_seaice(ims,ime,jms,jme,its,ite,jts,jte,xice_threshold,xice_p, &
                                 tsk_p,tsk_sea,tsk_ice)
    do i = its,ite
       tsk_p(i,j) = tsk_ice(i,j)
    enddo
 enddo

 do j = jts,jte
    do i = its,ite
       sr_p(i,j)     = 0._RKIND
       rainbl_p(i,j) = 0._RKIND
    enddo
    if(microp_scheme .ne. 'off') then
       call mpas_pool_get_array(diag_physics,'sr',sr)
       call mpas_pool_get_array(diag_physics,'rainncv',rainncv)
       do i = its,ite
          sr_p(i,j) = sr(i)
          rainbl_p(i,j) = rainbl_p(i,j) + rainncv(i)
       enddo
    endif
    if(convection_scheme .ne. 'off') then
       call mpas_pool_get_array(diag_physics,'raincv',raincv)
       do i = its,ite
          rainbl_p(i,j) = rainbl_p(i,j) + raincv(i)
       enddo
    endif
 enddo

!call mpas_log_write('--- end subroutine seaice_from_MPAS:')

 end subroutine seaice_from_MPAS

!=================================================================================================================
 subroutine seaice_to_MPAS(configs,mesh,diag_physics,sfc_input,its,ite)
!=================================================================================================================

!input and inout arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: sfc_input
 integer,intent(in):: its,ite

!local pointers:
 character(len=StrKIND),pointer:: microp_scheme, &
                                  lsm_scheme

 real(kind=RKIND),dimension(:),pointer:: acsnom,acsnow,cpm,chs,chs2,cqs,cqs2,qgh,qsfc, &
                                         grdflx,hfx, qfx,lh,noahres,snopcx,z0,znt,     &
                                         potevp,sfc_albedo,sfc_emiss,sfcrunoff
 real(kind=RKIND),dimension(:),pointer:: snow,snowc,snowh,skintemp,xice
 real(kind=RKIND),dimension(:),pointer:: t2m,th2m,q2
 real(kind=RKIND),dimension(:,:),pointer:: tslb

! RUC lsm variables
 real(kind=RKIND),dimension(:),pointer  :: dew,qcg,qsg,qvg,snowfallac,soilt1,tsnav

 integer,pointer:: globalCells(:)

!local variables and arrays:
 integer:: i,j,n

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('--- enter subroutine seaice_to_MPAS:')

 call mpas_pool_get_config(configs,'config_lsm_scheme',lsm_scheme)
 call mpas_pool_get_array(mesh,'indexToCellID',globalCells)

 if (lsm_scheme .eq. "sf_ruc") then
   call mpas_pool_get_array(diag_physics,'dew'       ,dew       )
   call mpas_pool_get_array(diag_physics,'qcg'       ,qcg       )
   call mpas_pool_get_array(diag_physics,'qsg'       ,qsg       )
   call mpas_pool_get_array(diag_physics,'qvg'       ,qvg       )
   call mpas_pool_get_array(diag_physics,'snowfallac',snowfallac)
   call mpas_pool_get_array(diag_physics,'soilt1'    ,soilt1    )
   call mpas_pool_get_array(diag_physics,'tsnav'     ,tsnav     )
 endif

 call mpas_pool_get_array(diag_physics,'acsnom'    ,acsnom    )
 call mpas_pool_get_array(diag_physics,'acsnow'    ,acsnow    )
 call mpas_pool_get_array(diag_physics,'chs'       ,chs       )
 call mpas_pool_get_array(diag_physics,'cpm'       ,cpm       )
 call mpas_pool_get_array(diag_physics,'chs2'      ,chs2      )
 call mpas_pool_get_array(diag_physics,'cqs2'      ,cqs2      )
 call mpas_pool_get_array(diag_physics,'cqs'       ,cqs       )
 call mpas_pool_get_array(diag_physics,'qgh'       ,qgh       )
 call mpas_pool_get_array(diag_physics,'qsfc'      ,qsfc      )
 call mpas_pool_get_array(diag_physics,'grdflx'    ,grdflx    )
 call mpas_pool_get_array(diag_physics,'hfx'       ,hfx       )
 call mpas_pool_get_array(diag_physics,'qfx'       ,qfx       )
 call mpas_pool_get_array(diag_physics,'lh'        ,lh        )
 call mpas_pool_get_array(diag_physics,'noahres'   ,noahres   )
 call mpas_pool_get_array(diag_physics,'potevp'    ,potevp    )
 call mpas_pool_get_array(diag_physics,'sfc_albedo',sfc_albedo)
 call mpas_pool_get_array(diag_physics,'sfc_emiss' ,sfc_emiss )
 call mpas_pool_get_array(diag_physics,'sfcrunoff' ,sfcrunoff )
 call mpas_pool_get_array(diag_physics,'snopcx'    ,snopcx    )
 call mpas_pool_get_array(diag_physics,'z0'        ,z0        )
 call mpas_pool_get_array(diag_physics,'znt'       ,znt       )
 call mpas_pool_get_array(diag_physics,'t2m'       ,t2m       )
 call mpas_pool_get_array(diag_physics,'th2m'      ,th2m      )
 call mpas_pool_get_array(diag_physics,'q2'        ,q2        )

 call mpas_pool_get_array(sfc_input,'snow'    ,snow    )
 call mpas_pool_get_array(sfc_input,'snowc'   ,snowc   )
 call mpas_pool_get_array(sfc_input,'snowh'   ,snowh   )
 call mpas_pool_get_array(sfc_input,'skintemp',skintemp)
 call mpas_pool_get_array(sfc_input,'tslb'    ,tslb    )
 call mpas_pool_get_array(sfc_input,'xice'    ,xice    )

!--- weigh local variables needed in the calculation of t2m, th2m, and q2 over seaice points:
 do j = jts,jte
    do i = its,ite
       if(xice_p(i,j).ge.xice_threshold .and. xice_p(i,j).le.1._RKIND) then
          cpm_p(i,j)  = xice_p(i,j)*cpm_p(i,j)  + (1._RKIND-xice_p(i,j))*cpm_sea(i,j)
          chs_p(i,j)  = xice_p(i,j)*chs_p(i,j)  + (1._RKIND-xice_p(i,j))*chs_sea(i,j)
          chs2_p(i,j) = xice_p(i,j)*chs2_p(i,j) + (1._RKIND-xice_p(i,j))*chs2_sea(i,j)
          cqs2_p(i,j) = xice_p(i,j)*cqs2_p(i,j) + (1._RKIND-xice_p(i,j))*cqs2_sea(i,j)
          cqs_p(i,j)  = xice_p(i,j)*cqs_p(i,j)  + (1._RKIND-xice_p(i,j))*cqs_sea(i,j)
          hfx_p(i,j)  = xice_p(i,j)*hfx_p(i,j)  + (1._RKIND-xice_p(i,j))*hfx_sea(i,j)
          lh_p(i,j)   = xice_p(i,j)*lh_p(i,j)   + (1._RKIND-xice_p(i,j))*lh_sea(i,j)
          qfx_p(i,j)  = xice_p(i,j)*qfx_p(i,j)  + (1._RKIND-xice_p(i,j))*qfx_sea(i,j)
          qgh_p(i,j)  = xice_p(i,j)*qgh_p(i,j)  + (1._RKIND-xice_p(i,j))*qgh_sea(i,j)
          qsfc_p(i,j) = xice_p(i,j)*qsfc_p(i,j) + (1._RKIND-xice_p(i,j))*qsfc_sea(i,j)
          tsk_p(i,j)  = xice_p(i,j)*tsk_p(i,j)  + (1._RKIND-xice_p(i,j))*tsk_sea(i,j)
          sfc_albedo_p(i,j) = xice_p(i,j)*sfc_albedo_p(i,j) + (1._RKIND-xice_p(i,j))*0.08_RKIND
          sfc_emiss_p(i,j)  = xice_p(i,j)*sfc_emiss_p(i,j)  + (1._RKIND-xice_p(i,j))*0.98_RKIND
       endif
    enddo
 enddo

 !call to 2-m diagnostics
 lsm_select: select case (trim(lsm_scheme))

    case("sf_noah")
       call sfcdiags( &
                     hfx   = hfx_p  , qfx     = qfx_p   , tsk  = tsk_p , qsfc = qsfc_p , chs = chs_p , &
                     chs2  = chs2_p , cqs2    = cqs2_p  , t2   = t2m_p , th2  = th2m_p , q2  = q2_p  , &
                     psfc  = psfc_p , t3d     = t_p     , qv3d = qv_p  , cp   = cp     , R_d = R_d   , &
                     rovcp = rcp    , ua_phys = ua_phys                                              , &
                     ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,           &
                     ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme ,           &
                     its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte             &
                    )
    case("sf_ruc")
       call mpas_log_write('--- call subroutine sfcdiags_ruclsm:')
       do j = jts,jte
          do i = its,ite
             if (globalcells(i)==targetcell) then
                print *,' before call to sfcdiags_ruclsm'
                print *,'globalcells(i),cqs_p,cqs2_p,chs_p,chs2_p',globalcells(i),cqs_p(i,j),cqs2_p(i,j),chs_p(i,j),chs2_p(i,j)
             endif
          enddo
       enddo
       call sfcdiags_ruclsm( &
                     hfx     = hfx_p  , qfx     = qfx_p   , tsk  = tsk_p       , qsfc = qsfc_p , chs = chs_p , &
                     chs2    = chs2_p , cqs2    = cqs2_p  , t2   = t2m_p       , th2  = th2m_p , q2  = q2_p  , &
                     psfc2d  = psfc_p , t3d  = t_p , qv3d  = qv_p , dz  = dz_p , cp   = cp     , R_d = R_d   , &
                     rovcp   = rcp    , rho3d   = rho_p   , p3d  = pres2_hyd_p , cqs = cqs_p,                  &
                     globalCells = globalCells,                                                                &
                     ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,                   &
                     ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme ,                   &
                     its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte                     &
                    )

    case default

 end select lsm_select

!--- update all variables for sea-ice points:
 do j = jts,jte
    do i = its,ite
       !--- inout variables:
       do n = 1,num_soils
          tslb(n,i) = tslb_p(i,n,j)
       enddo
       z0(i)        = z0_p(i,j)
       snow(i)      = snow_p(i,j)
       snowc(i)     = snowc_p(i,j)
       snowh(i)     = snowh_p(i,j)
       skintemp(i)  = tsk_p(i,j)
       acsnom(i)    = acsnom_p(i,j)
       acsnow(i)    = acsnow_p(i,j)
       sfcrunoff(i) = sfcrunoff_p(i,j)
       !--- inout optional variables:
       potevp(i)    = potevp_p(i,j)
       snopcx(i)    = snopcx_p(i,j)

       !--- output variables:
       znt(i)       = znt_p(i,j)
       grdflx(i)    = grdflx_p(i,j)
       !--- output optional variables:
       noahres(i)   = noahres_p(i,j)

       cpm(i)  = cpm_p(i,j)
       chs(i)  = chs_p(i,j)
       chs2(i) = chs2_p(i,j)
       cqs2(i) = cqs2_p(i,j)
       cqs(i)  = cqs_p(i,j)
       qsfc(i) = qsfc_p(i,j)
       qgh(i)  = qgh_p(i,j)
       hfx(i)  = hfx_p(i,j)
       qfx(i)  = qfx_p(i,j)
       lh(i)   = lh_p(i,j)
       sfc_albedo(i) = sfc_albedo_p(i,j)
       sfc_emiss(i)  = sfc_emiss_p(i,j)

       if(lsm_scheme .eq. "sf_ruc") then
          qcg(i) = qcg_p(i,j)
          qsg(i) = qsg_p(i,j)
          qvg(i) = qvg_p(i,j)
          dew(i) = dew_p(i,j)
          soilt1(i) = soilt1_p(i,j)
          tsnav(i) = tsnav_p(i,j)
          snowfallac(i) = snowfallac_p(i,j)
       endif

       !--- 2-meter diagnostics:
       q2(i)   = q2_p(i,j)
       t2m(i)  = t2m_p(i,j)
       th2m(i) = th2m_p(i,j)
    enddo
 enddo

!call mpas_log_write('--- end subroutine seaice_to_MPAS:')

 end subroutine seaice_to_MPAS

!=================================================================================================================
 subroutine driver_seaice(itimestep,configs,mesh,diag_physics,sfc_input,its,ite)
!=================================================================================================================

 use mpas_pool_routines,only: mpas_pool_get_array,mpas_pool_get_config,mpas_pool_type
!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs
 integer,intent(in):: its,ite
 integer,intent(in):: itimestep

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: sfc_input

!local pointers:
 integer:: i,j
 character(len=StrKIND),pointer:: lsm_scheme

 character(len=StrKIND),pointer:: mminlu
 integer,pointer:: isice,iswater
 integer,pointer:: spp_lsm_loc 
 integer,pointer:: globalCells(:)



!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
!call mpas_log_write('--- enter subroutine driver_seaice: xice_threshold = $r',realArgs=(/xice_threshold/))

 call mpas_pool_get_config(configs,'config_lsm_scheme',lsm_scheme)
 call mpas_pool_get_array(sfc_input,'mminlu',mminlu)
 call mpas_pool_get_array(sfc_input,'isice' ,isice )
 call mpas_pool_get_array(sfc_input,'iswater',iswater)
 call mpas_pool_get_array(diag_physics,'spp_lsm_loc', spp_lsm_loc)
 call mpas_pool_get_array(mesh,'indexToCellID',globalCells)


!copy MPAS arrays to local arrays:
 call seaice_from_MPAS(configs,mesh,diag_physics,sfc_input,its,ite)

!call to seaice model
 lsm_select: select case (trim(lsm_scheme))

    case("sf_noah")

       call seaice_noah( &
             dz8w      = dz_p        , p8w3d   = pres2_hyd_p  , t3d      = t_p      , &
             qv3d      = qv_p        , xice    = xice_p       , snoalb2d = snoalb_p , &
             glw       = glw_p       , swdown  = swdown_p     , rainbl   = rainbl_p , &
             sr        = sr_p        , qgh     = qgh_p        , tsk      = tsk_p    , &
             hfx       = hfx_p       , qfx     = qfx_p        , lh       = lh_p     , &
             grdflx    = grdflx_p    , potevp  = potevp_p     , qsfc     = qsfc_p   , &
             emiss     = sfc_emiss_p , albedo  = sfc_albedo_p , rib      = br_p     , &
             cqs2      = cqs2_p      , chs     = chs_p        , chs2     = chs2_p   , &
             z02d      = z0_p        , znt     = znt_p        , tslb     = tslb_p   , &
             snow      = snow_p      , snowc   = snowc_p      , snowh2d  = snowh_p  , &
             snopcx    = snopcx_p    , acsnow  = acsnow_p     , acsnom   = acsnom_p , &
             sfcrunoff = sfcrunoff_p , albsi   = albsi_p      , snowsi   = snowsi_p , &
             icedepth  = icedepth_p  , noahres = noahres_p    , dt       = dt_pbl   , &
             frpcpn    = frpcpn      ,                                                &
             seaice_albedo_opt        = seaice_albedo_opt        ,                    &
             seaice_albedo_default    = seaice_albedo_default    ,                    &
             seaice_thickness_opt     = seaice_thickness_opt     ,                    &
             seaice_thickness_default = seaice_thickness_default ,                    &
             seaice_snowdepth_opt     = seaice_snowdepth_opt     ,                    &
             seaice_snowdepth_max     = seaice_snowdepth_max     ,                    &
             seaice_snowdepth_min     = seaice_snowdepth_min     ,                    &
             xice_threshold           = xice_threshold           ,                    &
             num_soil_layers          = num_soils                ,                    &
             sf_urban_physics         = sf_urban_physics         ,                    &
             ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,  &
             ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme ,  &
             its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte    &
                 )
    case("sf_ruc")

       call mpas_log_write('--- call subroutine ruc_ice:')
       call ruc_ice( spp_lsm       = spp_lsm_loc   , iswater      = iswater     , chklowq   = chklowq_p,   &
                    rhosnf         = rhosnf_p      , precipfr     = precipfr_p  , chs       = chs_p,       &
                    qsg            = qsg_p         , qvg          = qvg_p       , dew       = dew_p,       &
                    soilt1         = soilt1_p      , tsnav        = tsnav_p     , acrunoff  = acrunoff_p,  &
                    snowfallac     = snowfallac_p  , stbolt       = stbolt      , qsfc      = qsfc_p,      &
                    graupelncv     = graupelncv_p  , snowncv      = snowncv_p   , rainncv   = rainncv_p,   &
                    qcg            = qcg_p         , flqc         = flqc_p      , flhc      = flhc_p,      &
                    dt             = dt_pbl        , ktau         = itimestep   , nsl       = num_soils,   &
                    rainbl         = rainbl_p      , snow         = snow_p      , snowh     = snowh_p,     &
                    snowc          = snowc_p       , frzfrac      = sr_p        , frpcpn    = frpcpn_ruc,  &
                    z3d            = dz_p          , p8w          = pres2_hyd_p , t3d       = t_p,         &
                    qv3d           = qv_p          , qc3d         = qc_p        , rho3d     = rho_p,       &
                    glw            = glw_p         , gsw          = gsw_p       , emiss     = sfc_emiss_p, &
                    alb            = sfc_albedo_p  , znt          = znt_p       , zs        = zs_p,        &
                    z0             = z0_p          , snoalb       = snoalb_p    , albbck    = sfc_albbck_p,&
                    tbot           = tmn_p         , ivgtyp       = ivgtyp_p    , isltyp    = isltyp_p,    &
                    xland          = xland_p       , isice        = isice       , xice      = xice_p,      &
                    xice_threshold = xice_threshold, cp           = cp          , rovcp     = rcp,         &
                    g0             = gravity       , lv           = xlv         , snom      = acsnom_p,    &
                    tso            = tslb_p        , soilt        = tsk_p       , hfx       = hfx_p,       &
                    qfx            = qfx_p         , lh           = lh_p        , sfcrunoff = sfcrunoff_p, &
                    grdflx         = grdflx_p      , acsnow       = acsnow_p    , sfcevp    = sfcevp_p,    &
                    myjpbl         = myjpbl        , mminlu       = mminlu      , globalCells=globalCells, &
                    ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,  &
                    ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme ,  &
                    its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte    &
                          )

    case default

 end select lsm_select

!copy local arrays to MPAS grid:
 call seaice_to_MPAS(configs,mesh,diag_physics,sfc_input,its,ite)

!call mpas_log_write('--- end subroutine driver_seaice:')

 end subroutine driver_seaice

!=================================================================================================================
 end module mpas_atmphys_driver_seaice
!=================================================================================================================
