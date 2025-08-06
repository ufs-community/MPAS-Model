! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_driver_smoke
 use mpas_kind_types
 use mpas_pool_routines
 use mpas_timer, only : mpas_timer_start, mpas_timer_stop
 use mpas_atmphys_manager, only: gmt,curr_julday,julday,year

 use mpas_atmphys_constants
 use mpas_atmphys_utilities
 use mpas_atmphys_vars

 use mpas_smoke_wrapper

 implicit none
 private
 public:: allocate_smoke,     &
          deallocate_smoke,   &
          init_smoke,         &
          driver_smoke

!MPAS driver for smoke.
!Haiqin Li (haiqin.li@noaa.gov).
!2024-09-25.

 contains

!=================================================================================================================
 subroutine allocate_smoke(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 
 character(len=StrKIND),pointer :: config_smoke_scheme
 character(len=StrKIND),pointer :: config_dust_scheme
 integer, pointer :: ebb_dcycle
 integer, pointer :: wetdep_ls_opt
 integer, pointer :: drydep_opt
 integer, pointer :: plumerise_opt
 logical, pointer :: add_fire_heat_flux, add_fire_moist_flux
 integer:: i,k,j

 call mpas_pool_get_config(configs,'config_smoke_scheme', config_smoke_scheme)
 call mpas_pool_get_config(configs,'config_dust_scheme', config_dust_scheme)
 call mpas_pool_get_config(configs,'ebb_dcycle', ebb_dcycle)
 call mpas_pool_get_config(configs,'wetdep_ls_opt',wetdep_ls_opt)
 call mpas_pool_get_config(configs,'drydep_opt',drydep_opt)
 call mpas_pool_get_config(configs,'plumerise_opt',plumerise_opt)
 call mpas_pool_get_config(configs,'add_fire_heat_flux',add_fire_heat_flux)
 call mpas_pool_get_config(configs,'add_fire_moist_flux',add_fire_moist_flux)


 if(.not.allocated(dx_p)        ) allocate(dx_p(ims:ime,jms:jme)                 )
 if(.not.allocated(area_p)      ) allocate(area_p(ims:ime,jms:jme)               )
 if(.not.allocated(xland_p)     ) allocate(xland_p(ims:ime,jms:jme)              )
 if(.not.allocated(u10_p)       ) allocate(u10_p(ims:ime,jms:jme)                )
 if(.not.allocated(v10_p)       ) allocate(v10_p(ims:ime,jms:jme)                )
 if(.not.allocated(ust_p)       ) allocate(ust_p(ims:ime,jms:jme)                )
 if(.not.allocated(xlat_p)      ) allocate(xlat_p(ims:ime,jms:jme)               )
 if(.not.allocated(xlon_p)      ) allocate(xlon_p(ims:ime,jms:jme)               )
 if(.not.allocated(vegfra_p)    ) allocate(vegfra_p(ims:ime,jms:jme)             )
 if(.not.allocated(ivgtyp_p)    ) allocate(ivgtyp_p(ims:ime,jms:jme)             )
 if(.not.allocated(landusef_p)  ) allocate(landusef_p(ims:ime,1:num_landcat,jms:jme))

 if(.not.allocated(tsk_p)       ) allocate(tsk_p(ims:ime,jms:jme)                )
 if(.not.allocated(hpbl_p)      ) allocate(hpbl_p(ims:ime,jms:jme)               )
 if(.not.allocated(t2m_p)       ) allocate(t2m_p(ims:ime,jms:jme)                )
! if(.not.allocated(w_p)         ) allocate(w_p(ims:ime,kms:kme,jms:jme)          )
 if(.not.allocated(isltyp_p)    ) allocate(isltyp_p(ims:ime,jms:jme)             )
 if(.not.allocated(swdown_p)    ) allocate(swdown_p(ims:ime,jms:jme)             )
 if(.not.allocated(z0_p)        ) allocate(z0_p(ims:ime,jms:jme)                 )
 if(.not.allocated(snowh_p)     ) allocate(snowh_p(ims:ime,jms:jme)              )
 if(.not.allocated(rmol_p)      ) allocate(rmol_p(ims:ime,jms:jme)               )
 if(.not.allocated(raincv_p)    ) allocate(raincv_p(ims:ime,jms:jme)             )
 if(.not.allocated(rainncv_p)   ) allocate(rainncv_p(ims:ime,jms:jme)            )
 if(.not.allocated(kpbl_p)      ) allocate(kpbl_p(ims:ime,jms:jme)               )
 if(.not.allocated(znt_p)       ) allocate(znt_p(ims:ime,jms:jme)                )
 if(.not.allocated(dpt2m_p)     ) allocate(dpt2m_p(ims:ime,jms:jme)              )
 if(.not.allocated(mavail_p)    )  allocate(mavail_p(ims:ime,jms:jme)            )

 if(.not.allocated(cldfrac_p)   ) allocate(cldfrac_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(ktop_deep_p) ) allocate(ktop_deep_p(ims:ime,jms:jme)        )

 if(.not.allocated(smois_p)     ) allocate(smois_p(ims:ime,1:num_soils,jms:jme)  )
 if(.not.allocated(tslb_p)      ) allocate(tslb_p(ims:ime,1:num_soils,jms:jme)   )

 if(.not.allocated(zgrid_p)     ) allocate(zgrid_p(ims:ime,kms:kme,jms:jme)      )

 if ( config_smoke_scheme .ne. 'off' ) then
    if(.not.allocated(frp_out_p)   )         allocate(frp_out_p(ims:ime,jms:jme))
    if(.not.allocated(fre_out_p)   )         allocate(fre_out_p(ims:ime,jms:jme))
    if(.not.allocated(hwp_p)       )         allocate(hwp_p(ims:ime,jms:jme))
    if(.not.allocated(coef_bb_dc_p))         allocate(coef_bb_dc_p(ims:ime,jms:jme))
    if (ebb_dcycle .eq. 1 ) then
       if(.not.allocated(frp_in_p)    )         allocate(frp_in_p(ims:ime,jms:jme))
       if(.not.allocated(fre_in_p)            )  allocate(fre_in_p(ims:ime,jms:jme))
    else
       if(.not.allocated(totprcp_prev24_p)    )  allocate(totprcp_prev24_p(ims:ime,jms:jme,24))
       if(.not.allocated(hwp_prev24_p)        )  allocate(hwp_prev24_p(ims:ime,jms:jme,24))
       if(.not.allocated(frp_prev24_p)        )  allocate(frp_prev24_p(ims:ime,jms:jme,24))
       if(.not.allocated(fre_prev24_p)        )  allocate(fre_prev24_p(ims:ime,jms:jme,24))
       if(.not.allocated(fire_end_hr_p)       )  allocate(fire_end_hr_p(ims:ime,jms:jme))
    endif
    if (add_fire_heat_flux) then 
       if(.not.allocated(hfx_bb_p)           )  allocate(hfx_bb_p(ims:ime,jms:jme))
       if(.not.allocated(frac_grid_burned_p) )  allocate(frac_grid_burned_p(ims:ime,jms:jme))
    endif
    if (add_fire_moist_flux) then
       if(.not.allocated(qfx_bb_p)           )  allocate(qfx_bb_p(ims:ime,jms:jme))
    endif
    if ( plumerise_opt .gt. 0 ) then
       if(.not.allocated(min_bb_plume_p))        allocate(min_bb_plume_p(ims:ime,jms:jme))
       if(.not.allocated(max_bb_plume_p))        allocate(max_bb_plume_p(ims:ime,jms:jme))
    endif
 endif
 if ( config_dust_scheme .ne. 'off' ) then
    if(.not.allocated(clayfrac_in_p))        allocate(clayfrac_in_p(ims:ime,jms:jme))
    if(.not.allocated(sandfrac_in_p))        allocate(sandfrac_in_p(ims:ime,jms:jme))
    if(.not.allocated(uthres_in_p))          allocate(uthres_in_p(ims:ime,jms:jme))
    if(.not.allocated(uthres_sg_in_p))       allocate(uthres_sg_in_p(ims:ime,jms:jme))
    if(.not.allocated(albedo_drag_p))        allocate(albedo_drag_p(ims:ime,jms:jme))
    if(.not.allocated(sep_in_p))             allocate(sep_in_p(ims:ime,jms:jme))
    if(.not.allocated(feff_p))               allocate(feff_p(ims:ime,jms:jme))
 endif

 if(.not.allocated(ddvel_p))              allocate(ddvel_p(ims:ime,jms:jme,1:num_chem))

 if ( wetdep_ls_opt .gt. 0 ) then
    if(.not.allocated(wetdep_resolved_p))    allocate(wetdep_resolved_p(ims:ime,jms:jme,1:num_chem))
 endif
 if ( drydep_opt .eq. 1 ) then
    if(.not.allocated(tend_chem_settle_p))   allocate(tend_chem_settle_p(ims:ime,kms:kme,jms:jme,1:num_chem))
 endif

 if(.not.allocated(chem_p))               allocate(chem_p(ims:ime,kms:kme,jms:jme,1:num_chem))
 if(.not.allocated(PM2_5_p))              allocate(PM2_5_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(PM10_p))               allocate(PM10_p(ims:ime,kms:kme,jms:jme))
 
 if ( num_e_bb_in .gt. 0 .and. config_smoke_scheme .ne. 'off' ) then
    if(.not.allocated(e_bb_in_p))           allocate(e_bb_in_p (ims:ime,1:kfire,jms:jme,1:num_e_bb_in) )
 endif
 if ( num_e_bb_out .gt. 0 .and. config_smoke_scheme .ne. 'off' ) then
    if(.not.allocated(e_bb_out_p))          allocate(e_bb_out_p  (ims:ime,kms:kme,jms:jme,1:num_e_bb_out)  )
 endif
 if ( num_e_dust_out .gt. 0 .and. config_dust_scheme .ne. 'off' ) then
    if(.not.allocated(e_dust_out_p))        allocate(e_dust_out_p(ims:ime,kms:kme,jms:jme,1:num_e_dust_out))
 endif

 end subroutine allocate_smoke

!=================================================================================================================
 subroutine deallocate_smoke(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

 if(allocated(dx_p)           ) deallocate(dx_p           )
 if(allocated(area_p)         ) deallocate(area_p         )
 if(allocated(xland_p)        ) deallocate(xland_p        )
 if(allocated(u10_p)          ) deallocate(u10_p          )
 if(allocated(v10_p)          ) deallocate(v10_p          )
 if(allocated(ust_p)          ) deallocate(ust_p          )
 if(allocated(xlat_p)         ) deallocate(xlat_p         )
 if(allocated(xlon_p)         ) deallocate(xlon_p         )
 if(allocated(vegfra_p)       ) deallocate(vegfra_p       )
 if(allocated(ivgtyp_p)       ) deallocate(ivgtyp_p       )
 if(allocated(landusef_p)     ) deallocate(landusef_p     )

 if(allocated(tsk_p)          ) deallocate(tsk_p          )
 if(allocated(hpbl_p)         ) deallocate(hpbl_p         )
 if(allocated(t2m_p)          ) deallocate(t2m_p          )
 if(allocated(isltyp_p)       ) deallocate(isltyp_p       )
 if(allocated(swdown_p)       ) deallocate(swdown_p       )
 if(allocated(z0_p)           ) deallocate(z0_p           )
 if(allocated(snowh_p)        ) deallocate(snowh_p        )
 if(allocated(rmol_p)         ) deallocate(rmol_p         )
 if(allocated(raincv_p)       ) deallocate(raincv_p       )
 if(allocated(rainncv_p)      ) deallocate(rainncv_p      )
 if(allocated(kpbl_p)         ) deallocate(kpbl_p         )
 if(allocated(znt_p)          ) deallocate(znt_p          )
 if(allocated(dpt2m_p)        ) deallocate(dpt2m_p        )
 if(allocated(mavail_p)       ) deallocate(mavail_p       )
 
 if(allocated(cldfrac_p)      ) deallocate(cldfrac_p      )
 if(allocated(ktop_deep_p)    ) deallocate(ktop_deep_p    )

 if(allocated(smois_p)        ) deallocate(smois_p        )
 if(allocated(tslb_p)         ) deallocate(tslb_p         )

 if(allocated(frp_in_p)       )       deallocate(frp_in_p     )
 if(allocated(frp_out_p)      )       deallocate(frp_out_p    )
 if(allocated(fre_in_p)       )       deallocate(fre_in_p     )
 if(allocated(fre_out_p)      )       deallocate(fre_out_p    )

 if(allocated(totprcp_prev24_p))   deallocate(totprcp_prev24_p)
 if(allocated(hwp_prev24_p))       deallocate(hwp_prev24_p)
 if(allocated(frp_prev24_p))       deallocate(frp_prev24_p)
 if(allocated(fre_prev24_p))       deallocate(fre_prev24_p)
 if(allocated(fire_end_hr_p))      deallocate(fire_end_hr_p)

 if(allocated(clayfrac_in_p)    )     deallocate(clayfrac_in_p    )
 if(allocated(sandfrac_in_p)    )     deallocate(sandfrac_in_p    )
 if(allocated(uthres_in_p)      )     deallocate(uthres_in_p      )
 if(allocated(uthres_sg_in_p)   )     deallocate(uthres_sg_in_p   )
 if(allocated(albedo_drag_p)    )     deallocate(albedo_drag_p    )
 if(allocated(feff_p)           )     deallocate(feff_p           )
 if(allocated(sep_in_p)         )     deallocate(sep_in_p         )

 if(allocated(ddvel_p)       )        deallocate(ddvel_p)
 if(allocated(wetdep_resolved_p)   )  deallocate(wetdep_resolved_p)
 if(allocated(tend_chem_settle_p))    deallocate(tend_chem_settle_p)

 if(allocated(chem_p)     )           deallocate(chem_p)
 if(allocated(PM2_5_p)    )           deallocate(PM2_5_p)
 if(allocated(PM10_p)     )           deallocate(PM10_p)

 if(allocated(e_bb_in_p))             deallocate(e_bb_in_p)

 if(allocated(e_bb_out_p))            deallocate(e_bb_out_p)
 if(allocated(e_dust_out_p))          deallocate(e_dust_out_p)
 
!-----------------------------------------------------------------------------------------------------------------

 end subroutine deallocate_smoke

 !=================================================================================================================
 subroutine smoke_from_MPAS(dt_dyn, time_lev, emis_input, state, configs,  &
                            mesh,sfc_input, diag_physics, tend_physics, its, ite       )
!=================================================================================================================

     !input arguments:
      type(mpas_pool_type),intent(in):: configs
      type(mpas_pool_type),intent(in):: mesh
      type(mpas_pool_type),intent(in):: diag_physics
      type(mpas_pool_type),intent(in):: sfc_input
      type(mpas_pool_type),intent(in):: tend_physics
      type(mpas_pool_type),intent(in):: state
      type(mpas_pool_type),intent(in):: emis_input
     
      integer,intent(in)             :: its,ite,time_lev
      real(kind=RKIND),intent(in)    :: dt_dyn
     
      !local pointers:
      integer,dimension(:),pointer               :: ktop_deep
      real(kind=RKIND),dimension(:),pointer      :: meshDensity
      real(kind=RKIND),dimension(:),pointer      :: xland
      real(kind=RKIND),pointer                   :: len_disp
      real(kind=RKIND),dimension(:),pointer      :: ust
      real(kind=RKIND),dimension(:),pointer      :: latCell,lonCell
      real(kind=RKIND),dimension(:),pointer      :: areaCell
      real(kind=RKIND),dimension(:),pointer      :: raincv,rainncv,mavail
      real(kind=RKIND),dimension(:),pointer      :: rmol,snowh,z0,t2m,dpt2m
      real(kind=RKIND),dimension(:),pointer      :: dewpoint_surface
      real(kind=RKIND),dimension(:),pointer      :: u10,v10,znt,hpbl,skintemp,vegfra
      real(kind=RKIND),dimension(:,:),pointer    :: cldfrac
      real(kind=RKIND),dimension(:,:),pointer    :: landusef
      real(kind=RKIND),dimension(:,:),pointer    :: smois,tslb
      real(kind=RKIND),dimension(:,:),pointer    :: zgrid
      real(kind=RKIND),dimension(:,:),pointer    :: ddvel
      real(kind=RKIND),dimension(:,:),pointer    :: wetdep_resolved
      real(kind=RKIND),dimension(:,:,:),pointer  :: tend_chem_settle
      real(kind=RKIND),dimension(:),pointer      :: frp_in, frp_out, fre_in, fre_out, hwp, coef_bb_dc, fire_end_hr
      real(kind=RKIND),dimension(:,:),pointer    :: totprcp_prev24, hwp_prev24, frp_prev24, fre_prev24
      real(kind=RKIND),dimension(:),pointer      :: hfx_bb, qfx_bb, frac_grid_burned
      integer,dimension(:),pointer               :: min_bb_plume, max_bb_plume
      real(kind=RKIND),dimension(:),pointer      :: sandfrac_in, clayfrac_in, uthres_in, uthres_sg_in, &
                                                    sep_in, albedo_drag, feff
      real(kind=RKIND),dimension(:,:,:),pointer  :: e_bb_in
      real(kind=RKIND),dimension(:,:,:),pointer  :: e_bb_out, e_dust_out
      integer,dimension(:),pointer               :: isltyp,ivgtyp,kpbl
     
      integer,pointer                            :: chemistry_start,chemistry_end

      integer,pointer                            :: ebb_dcycle
      character(len=StrKIND),pointer :: config_smoke_scheme
      character(len=StrKIND),pointer :: config_dust_scheme
      character(len=StrKIND),pointer :: config_convection_scheme
      integer, pointer :: wetdep_ls_opt
      integer, pointer :: drydep_opt
      integer, pointer :: plumerise_opt
      logical, pointer :: add_fire_heat_flux, add_fire_moist_flux
      
     
      real(kind=RKIND),dimension(:,:,:),pointer :: scalars
      real(kind=RKIND),dimension(:,:,:),pointer :: chem
     
      integer:: i,j,k,n,h
      
      call mpas_pool_get_array(state,'scalars',scalars,time_lev)
     
      call mpas_pool_get_dimension(state, 'chemistry_start', chemistry_start) 
      call mpas_pool_get_dimension(state, 'chemistry_end', chemistry_end)
      call mpas_pool_get_config(configs,'config_smoke_scheme', config_smoke_scheme)
      call mpas_pool_get_config(configs,'config_dust_scheme', config_dust_scheme)
      call mpas_pool_get_config(configs,'ebb_dcycle', ebb_dcycle)
      call mpas_pool_get_config(configs,'wetdep_ls_opt',wetdep_ls_opt)
      call mpas_pool_get_config(configs,'drydep_opt',drydep_opt)
      call mpas_pool_get_config(configs,'plumerise_opt',plumerise_opt)
      call mpas_pool_get_config(configs,'add_fire_heat_flux',add_fire_heat_flux)
      call mpas_pool_get_config(configs,'add_fire_moist_flux',add_fire_moist_flux)
      call mpas_pool_get_config(configs,'config_convection_scheme',config_convection_scheme)
     
      call mpas_pool_get_config(configs,'config_len_disp',len_disp)
      call mpas_pool_get_array(mesh,'meshDensity',meshDensity)
      call mpas_pool_get_array(mesh,'latCell',latCell)
      call mpas_pool_get_array(mesh,'lonCell',lonCell)
      call mpas_pool_get_array(mesh,'areaCell',areaCell)
      call mpas_pool_get_array(sfc_input,'xland',xland)
      call mpas_pool_get_array(diag_physics,'ust',ust)
     
      call mpas_pool_get_array(diag_physics,'u10',u10)
      call mpas_pool_get_array(diag_physics,'v10',v10)
      call mpas_pool_get_array(sfc_input,'skintemp',skintemp)
      call mpas_pool_get_array(diag_physics,'hpbl',hpbl)
      call mpas_pool_get_array(diag_physics,'kpbl',kpbl)
      call mpas_pool_get_array(diag_physics,'znt',znt)
      if ( config_convection_scheme .eq. 'cu_grell_freitas_li' ) then
         call mpas_pool_get_array(diag_physics,'ktop_deep',ktop_deep)
      endif
     
      call mpas_pool_get_array(mesh,'zgrid',zgrid)
      call mpas_pool_get_array(sfc_input,'smois',smois)
      call mpas_pool_get_array(sfc_input,'tslb',tslb)
      call mpas_pool_get_array(sfc_input,'isltyp',isltyp)
      call mpas_pool_get_array(sfc_input,'ivgtyp',ivgtyp)
      call mpas_pool_get_array(sfc_input,'vegfra',vegfra)
      call mpas_pool_get_array(sfc_input,'landusef',landusef)
      call mpas_pool_get_array(diag_physics, 'z0',z0)
      call mpas_pool_get_array(sfc_input,'snowh',snowh)
      call mpas_pool_get_array(diag_physics,'rmol',rmol)
      call mpas_pool_get_array(diag_physics,'raincv',raincv)
      call mpas_pool_get_array(diag_physics,'rainncv',rainncv)
     ! call mpas_pool_get_array(diag_physics,'dewpoint_surface',dewpoint_surface)
      call mpas_pool_get_array(diag_physics,'mavail',mavail)
      call mpas_pool_get_array(diag_physics,'cldfrac',cldfrac)
     
      call mpas_pool_get_array(diag_physics,'ddvel',ddvel)
      if ( wetdep_ls_opt .eq. 1 ) then 
         call mpas_pool_get_array(diag_physics,'wetdep_resolved',wetdep_resolved)
      endif
      call mpas_pool_get_array(tend_physics,'tend_chem_settle',tend_chem_settle)

      if ( config_dust_scheme .ne. 'off' ) then
         call mpas_pool_get_array(diag_physics, 'clayfrac_in',clayfrac_in)
         call mpas_pool_get_array(diag_physics, 'sandfrac_in',sandfrac_in)
         call mpas_pool_get_array(diag_physics, 'uthres_in',uthres_in)
         call mpas_pool_get_array(diag_physics, 'uthres_sg_in',uthres_sg_in)
         call mpas_pool_get_array(diag_physics, 'albedo_drag',albedo_drag) ! these have been updated to select
         call mpas_pool_get_array(diag_physics, 'feff',feff)               ! the correct month
         call mpas_pool_get_array(diag_physics, 'sep_in',sep_in)
      endif

      if (config_smoke_scheme .ne. 'off' .and. num_e_bb_in .gt. 0 ) then
         call mpas_pool_get_array(emis_input, 'e_bb_in', e_bb_in)
      endif

      if (config_smoke_scheme .ne. 'off' ) then 
         call mpas_pool_get_array(diag_physics,'frp_out',frp_out)
         call mpas_pool_get_array(diag_physics,'fre_out',fre_out)
         call mpas_pool_get_array(diag_physics,'hwp',hwp)
         call mpas_pool_get_array(diag_physics,'coef_bb_dc',coef_bb_dc)
         if ( ebb_dcycle .eq. 1 ) then
            call mpas_pool_get_array(diag_physics,'frp_in',frp_in)
            call mpas_pool_get_array(diag_physics,'fre_in',fre_in)
         else
            call mpas_pool_get_array(diag_physics,'totprcp_prev24',totprcp_prev24)
            call mpas_pool_get_array(diag_physics,'hwp_prev24',hwp_prev24)
            call mpas_pool_get_array(diag_physics,'frp_prev24',frp_prev24)
            call mpas_pool_get_array(diag_physics,'fre_prev24',fre_prev24)
            call mpas_pool_get_array(diag_physics,'fire_end_hr',fire_end_hr)
         endif
         if ( plumerise_opt .gt. 0 ) then
            call mpas_pool_get_array(diag_physics,'min_bb_plume',min_bb_plume)
            call mpas_pool_get_array(diag_physics,'max_bb_plume',max_bb_plume)
         endif
         if ( add_fire_heat_flux ) then
            call mpas_pool_get_array(diag_physics,'hfx_bb',hfx_bb)
            call mpas_pool_get_array(diag_physics,'frac_grid_burned',frac_grid_burned)
         endif
         if ( add_fire_moist_flux ) then
            call mpas_pool_get_array(diag_physics,'qfx_bb',qfx_bb)
         endif
      endif
      if ( num_e_bb_out .gt. 0 ) then
         call mpas_pool_get_array(diag_physics, 'e_bb_out',e_bb_out)
      endif
      if ( num_e_dust_out .gt. 0 ) then
         call mpas_pool_get_array(diag_physics, 'e_dust_out',e_dust_out)
      endif
     
      chem => scalars(chemistry_start:chemistry_end,:,:)
     
      do j = jts, jte
      do k = kts, kte
      do i = its, ite
         zgrid_p(i,k,j)      = zgrid(k,i)
         cldfrac_p(i,k,j)    = cldfrac(k,i)
         do n = 1, num_chem
            chem_p(i,k,j,n) = max(1.e-12_RKIND,chem(n,k,i))
            tend_chem_settle_p(i,k,j,n) = 0._RKIND
         enddo
      enddo
      enddo
      enddo
     
!
     if ( num_e_bb_in .gt. 0 ) then 
      do j = jts, jte
      do k = 1, kfire
      do i = its, ite
      do n = 1, num_e_bb_in
         e_bb_in_p(i,k,j,n) = e_bb_in(n,k,i)
      enddo
      enddo
      enddo
      enddo
     endif
!
      do j = jts,jte
      do n = 1,num_soils
      do i = its,ite
         smois_p(i,n,j)  = smois(n,i)
         tslb_p(i,n,j)   = tslb(n,i)
      enddo
      enddo
      enddo
     
      do j = jts,jte
      do i = its,ite
     ! mesh/grid/surf
         dx_p(i,j)      = len_disp / meshDensity(i)**0.25
         area_p(i,j)    = areaCell(i)
         xland_p(i,j)   = xland(i)
         ust_p(i,j)     = ust(i)
         xlat_p(i,j)    = latCell(i) / degrad
         xlon_p(i,j)    = lonCell(i) / degrad
         isltyp_p(i,j)  = isltyp(i)
         ivgtyp_p(i,j)  = ivgtyp(i)
         vegfra_p(i,j)  = vegfra(i)
     ! Met
         kpbl_p(i,j)    = kpbl(i)
         if ( config_convection_scheme .eq. 'cu_grell_freitas_li' ) then
            ktop_deep_p(i,j)=ktop_deep(i)
         else
            ktop_deep_p(i,j)=1
         endif
         u10_p(i,j)     = u10(i)
         v10_p(i,j)     = v10(i)
         znt_p(i,j)     = znt(i)
         tsk_p(i,j)     = skintemp(i)
         hpbl_p(i,j)    = hpbl(i)
         z0_p(i,j)      = z0(i)
         snowh_p(i,j)   = snowh(i)
         rmol_p(i,j)    = rmol(i)
         if ( config_convection_scheme .eq. 'cu_grell_freitas_li' ) then
            raincv_p(i,j)  = raincv(i)
         else
            raincv_p(i,j) = 0._RKIND
         endif
         rainncv_p(i,j) = rainncv(i)
         dpt2m_p(i,j)   = 280. !dewpoint_surface(i)
         mavail_p(i,j)  = mavail(i)
     !
     ! BB
     !
     if ( config_smoke_scheme .ne. 'off' ) then
         frp_out_p(i,j) = frp_out(i)
         fre_out_p(i,j) = fre_out(i)
         hwp_p(i,j)     = hwp(i)
         coef_bb_dc_p(i,j) = coef_bb_dc(i)
         if ( ebb_dcycle .eq. 1 ) then
            frp_in_p(i,j)  = frp_in(i)
            fre_in_p(i,j)  = fre_in(i)
         else
            do h = 1,24 
               totprcp_prev24_p(i,j,h) = totprcp_prev24(i,h)
               hwp_prev24_p(i,j,h)     = hwp_prev24(i,h)
               frp_prev24_p(i,j,h)     = frp_prev24(i,h)
               fre_prev24_p(i,j,h)     = fre_prev24(i,h)
               fire_end_hr_p(i,j)    = fire_end_hr(i)
            enddo
         endif
         if ( plumerise_opt .gt. 0 ) then   
            min_bb_plume_p(i,j) = 0
            max_bb_plume_p(i,j) = 0
         endif
         if ( add_fire_heat_flux ) then
            hfx_bb_p(i,j)       = 0._RKIND
            frac_grid_burned_p(i,j) = frac_grid_burned(i)
         endif
         if ( add_fire_moist_flux ) then
            qfx_bb_p(i,j)       = 0._RKIND
         endif
     endif
     !
     ! Dust
     !
     if ( config_dust_scheme .ne. 'off' ) then
         sandfrac_in_p(i,j)    = sandfrac_in(i)
         clayfrac_in_p(i,j)    = clayfrac_in(i)
         uthres_in_p(i,j)      = uthres_in(i)
         uthres_sg_in_p(i,j)   = uthres_sg_in(i)
         albedo_drag_p(i,j)    = albedo_drag(i)
         feff_p(i,j)           = feff(i)
         sep_in_p(i,j)         = sep_in(i)
     endif
     ! Chem diags
         do n= 1,num_chem
            ddvel_p(i,j,n) = 0._RKIND
            wetdep_resolved_p(i,j,n) = 0._RKIND
         enddo
     ! Land cat
         do n= 1,num_landcat
            landusef_p(i,n,j)= landusef(n,i)
         enddo
      enddo
      enddo
     
     ! Accumulated emissions or?
      do j = jts, jte
      do k = kts, kte
      do i = its, ite
         if (num_e_bb_out .gt. 0 ) then
         do n = 1, num_e_bb_out
            e_bb_out_p(i,k,j,n)  = e_bb_out(n,k,i)
         enddo
         endif
         if (num_e_dust_out .gt. 0 ) then
         do n = 1, num_e_dust_out
            e_dust_out_p(i,k,j,n) = e_dust_out(n,k,i)
         enddo
         endif
      enddo
      enddo
      enddo

 end subroutine smoke_from_MPAS

!=================================================================================================================
 subroutine smoke_to_MPAS(configs,time_lev,state,diag_physics,tend_physics,its,ite)
!=================================================================================================================
! This subroutine updates the MPAS arrays with the physics-oriented/working arrays 
! used in the smoke_driver. So anything that is modified by the driver needs to be 
! included here, primarily the scalars and any output diagnostics

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 integer,intent(in)             :: its,ite,time_lev

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics
 type(mpas_pool_type),intent(inout):: state

!local variables:
 !real(kind=RKIND),dimension(:,:)  ,pointer:: cldfrac
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars,chem
 real(kind=RKIND),dimension(:,:),  pointer:: ddvel, wetdep_resolved
 real(kind=RKIND),dimension(:,:,:),pointer:: tend_chem_settle
 real(kind=RKIND),dimension(:),    pointer:: frp_out, fre_out, hwp, coef_bb_dc
 real(kind=RKIND),dimension(:),    pointer:: hfx_bb, qfx_bb, frac_grid_burned
 real(kind=RKIND),dimension(:,:)  ,pointer:: PM2_5, PM10
 real(kind=RKIND),dimension(:,:,:),pointer:: e_bb_out, e_dust_out
 integer,dimension(:),pointer:: min_bb_plume, max_bb_plume
 
 integer,pointer :: chemistry_start,chemistry_end
 integer,pointer                :: ebb_dcycle
 character(len=StrKIND),pointer :: config_smoke_scheme
 character(len=StrKIND),pointer :: config_dust_scheme
 integer, pointer :: wetdep_ls_opt
 integer, pointer :: drydep_opt
 integer, pointer :: plumerise_opt
 logical, pointer :: add_fire_heat_flux, add_fire_moist_flux

 integer:: i,k,j,n

!local pointers:
 call mpas_pool_get_array(state,'scalars',scalars,time_lev)

 call mpas_pool_get_array(diag_physics,'PM2_5',PM2_5)
 call mpas_pool_get_array(diag_physics,'PM10',PM10)

 call mpas_pool_get_array(diag_physics,'ddvel',ddvel)
 call mpas_pool_get_array(diag_physics,'wetdep_resolved',wetdep_resolved)

 call mpas_pool_get_array(tend_physics,'tend_chem_settle',tend_chem_settle)

 call mpas_pool_get_dimension(state, 'chemistry_start', chemistry_start) 
 call mpas_pool_get_dimension(state, 'chemistry_end', chemistry_end) 
 call mpas_pool_get_config(configs,'config_smoke_scheme', config_smoke_scheme)
 call mpas_pool_get_config(configs,'config_dust_scheme', config_dust_scheme)
 call mpas_pool_get_config(configs,'ebb_dcycle', ebb_dcycle)
 call mpas_pool_get_config(configs,'wetdep_ls_opt',wetdep_ls_opt)
 call mpas_pool_get_config(configs,'drydep_opt',drydep_opt)
 call mpas_pool_get_config(configs,'plumerise_opt',plumerise_opt)
 call mpas_pool_get_config(configs,'add_fire_heat_flux',add_fire_heat_flux)
 call mpas_pool_get_config(configs,'add_fire_moist_flux',add_fire_moist_flux)

 if (num_e_bb_out .gt. 0 ) then
    call mpas_pool_get_array(diag_physics, 'e_bb_out',e_bb_out)
 endif
 if (num_e_dust_out .gt. 0 ) then
    call mpas_pool_get_array(diag_physics, 'e_dust_out',e_dust_out)
 endif
 if (config_smoke_scheme .ne. 'off' ) then
    call mpas_pool_get_array(diag_physics,'frp_out',frp_out)
    call mpas_pool_get_array(diag_physics,'fre_out',fre_out)
    call mpas_pool_get_array(diag_physics,'hwp',hwp)
    call mpas_pool_get_array(diag_physics,'coef_bb_dc',coef_bb_dc)
    if (plumerise_opt .gt. 0 ) then
       call mpas_pool_get_array(diag_physics,'min_bb_plume',min_bb_plume)
       call mpas_pool_get_array(diag_physics,'max_bb_plume',max_bb_plume)
    endif
    if (add_fire_heat_flux) then
       call mpas_pool_get_array(diag_physics,'hfx_bb',hfx_bb)
       call mpas_pool_get_array(diag_physics,'frac_grid_burned',frac_grid_burned)
    endif
    if (add_fire_moist_flux) then
       call mpas_pool_get_array(diag_physics,'qfx_bb',qfx_bb)
    endif
 endif

 chem => scalars(chemistry_start:chemistry_end,:,:)

 do j = jts,jte
 do k = kts,kte
 do i = its,ite
!
   PM2_5(k,i)     = PM2_5_p(i,k,j)
   PM10(k,i)      = PM10_p(i,k,j)
!
   do n = 1, num_chem
      chem(n,k,i) = chem_p(i,k,j,n)
      !TODO floating point exception
      !chem(n,k,i) = min(max(1.e-12_RKIND,chem_p(i,k,j,n)),5000._RKIND)
      tend_chem_settle(n,k,i) = tend_chem_settle_p(i,k,j,n)
   enddo
!
   if (num_e_bb_out .gt. 0 ) then
   do n = 1, num_e_bb_out
      e_bb_out(n,k,i) = e_bb_out_p(i,k,j,n)
   enddo
   endif
!
   if (num_e_dust_out .gt. 0 ) then
   do n = 1, num_e_dust_out
      e_dust_out(n,k,i) = e_dust_out_p(i,k,j,n)
   enddo
   endif
!
 enddo
 enddo
 enddo
 
 do j = jts,jte
 do i = its,ite
    if (config_smoke_scheme .ne. 'off') then
       frp_out(i)          = frp_out_p(i,j)
       fre_out(i)          = fre_out_p(i,j)
       hwp(i)              = hwp_p(i,j)
       coef_bb_dc(i)       = coef_bb_dc_p(i,j)
       if ( plumerise_opt .gt. 0 ) then
          min_bb_plume(i)     = min_bb_plume_p(i,j)
          max_bb_plume(i)     = max_bb_plume_p(i,j)
       endif
       if ( add_fire_heat_flux ) then
          hfx_bb(i)           = hfx_bb_p(i,j)
          frac_grid_burned(i) = frac_grid_burned_p(i,j)
       endif
       if ( add_fire_moist_flux ) then
          qfx_bb(i)           = qfx_bb_p(i,j)
       endif
    endif
    do n = 1,num_chem
      ddvel(n,i) = ddvel_p(i,j,n)
      if (wetdep_ls_opt .gt. 0) then
         wetdep_resolved(n,i) = wetdep_resolved_p(i,j,n)
      endif
    enddo
 enddo
 enddo

 end subroutine smoke_to_MPAS

!=================================================================================================================
 subroutine init_smoke(mesh,configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs

 end subroutine init_smoke

!=================================================================================================================
 subroutine driver_smoke(itimestep,time_lev,emis_input,state,configs,   &
                         mesh,sfc_input,diag_physics,tend_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in)   :: configs
 type(mpas_pool_type),intent(in)   :: mesh
 type(mpas_pool_type),intent(inout):: state
 type(mpas_pool_type),intent(in)   :: emis_input

 integer,intent(in):: its,ite
 integer,intent(in):: itimestep
 integer,intent(in):: time_lev

!inout arguments:
 type(mpas_pool_type),intent(inout):: sfc_input
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

 integer:: ktau,initflag
 real(kind=RKIND):: curr_secs, dt

 !local pointers:
 integer,pointer :: index_smoke_fine
 integer,pointer :: index_dust_fine, index_dust_coarse
 integer,pointer :: index_e_bb_in_smoke_fine
 integer,pointer :: index_e_bb_out_smoke_fine, index_e_dust_out_dust_fine, index_e_dust_out_dust_coarse
                  
 integer,pointer                :: chemistry_start
 !namelists
 logical,pointer                :: config_do_restart
 character(len=StrKIND),pointer :: config_smoke_scheme
 character(len=StrKIND),pointer :: config_dust_scheme
 logical                        :: do_mpas_smoke, do_mpas_dust
 integer,pointer                :: hwp_method
 real(kind=RKIND),pointer       :: hwp_alpha
 integer,pointer                :: wetdep_ls_opt
 real(kind=RKIND),pointer       :: wetdep_ls_alpha
 integer,pointer                :: plumerise_opt
 integer,pointer                :: plume_wind_eff
 real(kind=RKIND),pointer       :: plume_alpha
 real(kind=RKIND),pointer       :: bb_emis_scale_factor
 integer,pointer                :: ebb_dcycle
 integer,pointer                :: drydep_opt
 integer,pointer                :: pm_settling
 logical,pointer                :: add_fire_heat_flux
 logical,pointer                :: add_fire_moist_flux
 integer,pointer                :: plumerisefire_frq
 real(kind=RKIND),pointer       :: dust_alpha, dust_gamma
 real(kind=RKIND),pointer       :: dust_drylimit_factor, dust_moist_correction

 integer                        :: i,j,k

 do_mpas_smoke  = .false.
 do_mpas_dust   = .false.

!copy MPAS arrays to local arrays:
 call smoke_from_MPAS(dt_dyn, time_lev, emis_input, state, configs, mesh, &
                      sfc_input, diag_physics, tend_physics,its,ite)
 
 call mpas_pool_get_config(configs,'config_do_restart'       ,config_do_restart)
! Namelist: schemes
 call mpas_pool_get_config(configs,'config_smoke_scheme', config_smoke_scheme)
 call mpas_pool_get_config(configs,'config_dust_scheme', config_dust_scheme)
! Namelist: smoke
 call mpas_pool_get_config(configs,'hwp_method',hwp_method)
 call mpas_pool_get_config(configs,'hwp_alpha',hwp_alpha)
 call mpas_pool_get_config(configs,'plumerise_opt',plumerise_opt)
 call mpas_pool_get_config(configs,'plume_wind_eff',plume_wind_eff)
 call mpas_pool_get_config(configs,'plume_alpha',plume_alpha)
 call mpas_pool_get_config(configs,'bb_emis_scale_factor',bb_emis_scale_factor)
 call mpas_pool_get_config(configs,'ebb_dcycle',ebb_dcycle)
 call mpas_pool_get_config(configs,'add_fire_heat_flux',add_fire_heat_flux)
 call mpas_pool_get_config(configs,'add_fire_moist_flux',add_fire_moist_flux)
 call mpas_pool_get_config(configs,'plumerisefire_frq',plumerisefire_frq)
! Namelist: Wet/dry deposition
 call mpas_pool_get_config(configs,'wetdep_ls_opt',wetdep_ls_opt)
 call mpas_pool_get_config(configs,'wetdep_ls_alpha',wetdep_ls_alpha) 
 call mpas_pool_get_config(configs,'drydep_opt',drydep_opt)
 call mpas_pool_get_config(configs,'pm_settling',pm_settling)
! Namelist: dust
 call mpas_pool_get_config(configs,'dust_alpha', dust_alpha)
 call mpas_pool_get_config(configs,'dust_gamma', dust_gamma)
 call mpas_pool_get_config(configs,'dust_drylimit_factor' ,dust_drylimit_factor)
 call mpas_pool_get_config(configs,'dust_moist_correction', dust_moist_correction)
! Indexes: emissions in
 call mpas_pool_get_dimension(emis_input,'index_e_bb_in_smoke_fine',index_e_bb_in_smoke_fine)

! Indexes: emissions out
 call mpas_pool_get_dimension(diag_physics,'index_e_bb_out_smoke_fine',index_e_bb_out_smoke_fine)

 call mpas_pool_get_dimension(diag_physics,'index_e_dust_out_dust_fine',index_e_dust_out_dust_fine)
 call mpas_pool_get_dimension(diag_physics,'index_e_dust_out_dust_coarse',index_e_dust_out_dust_coarse)

! Tracer indexes
 call mpas_pool_get_dimension(state,'chemistry_start', chemistry_start) 
 call mpas_pool_get_dimension(state,'index_smoke_fine'  ,index_smoke_fine )
 call mpas_pool_get_dimension(state,'index_dust_fine'  ,index_dust_fine )
 call mpas_pool_get_dimension(state,'index_dust_coarse'  ,index_dust_coarse )

 initflag = 1
 if(config_do_restart .or. itimestep > 1) initflag = 0

 if ( config_smoke_scheme  .ne. 'off' .and. &
    (index_smoke_fine .gt. 0 )) then
     do_mpas_smoke = .true.
 endif
 if (config_dust_scheme   .ne. 'off'  .and. &
    (index_dust_fine .gt. 0 .or. index_dust_coarse .gt. 0 )) then
     do_mpas_dust  = .true.
 endif

 call mpas_timer_start('mpas_smoke')
 call mpas_smoke_driver(                                                                 &
           chem  = chem_p, num_chem = num_chem, chemistry_start = chemistry_start,       &
           kemit = kemit,  kfire = kfire, kbio = kbio, kvol = kvol,                      &
           index_smoke_fine = index_smoke_fine,                                          &  
           index_dust_fine  = index_dust_fine,  index_dust_coarse  = index_dust_coarse,  &
           index_e_bb_in_smoke_fine     = index_e_bb_in_smoke_fine,                      &
           index_e_bb_out_smoke_fine    = index_e_bb_out_smoke_fine,                     &
           index_e_dust_out_dust_fine   = index_e_dust_out_dust_fine,                    &
           index_e_dust_out_dust_coarse = index_e_dust_out_dust_coarse,                  &
           hwp     = hwp_p,                                                              &
           frp_in         = frp_in_p,            frp_out = frp_out_p,                    &
           fre_in         = fre_in_p,            fre_out = fre_out_p,                    &
           totprcp_prev24 = totprcp_prev24_p, hwp_prev24 = hwp_prev24_p,                 &
           frp_prev24     = frp_prev24_p,     fre_prev24 = fre_prev24_p,                 &
           hfx_bb         = hfx_bb_p,             qfx_bb = qfx_bb_p,                     &
           frac_grid_burned                              = frac_grid_burned_p,           &
           min_bb_plume   = min_bb_plume_p, max_bb_plume = max_bb_plume_p,               &
           coef_bb_dc     = coef_bb_dc_p,   e_bb_in       = e_bb_in_p,                   &
           sandfrac_in    = sandfrac_in_p, clayfrac_in    = clayfrac_in_p,               &
           uthres_in      = uthres_in_p,   uthres_sg_in   = uthres_sg_in_p,              &
           albedo_drag_in = albedo_drag_p, feff_in        = feff_p, sep_in = sep_in_p,   &
           e_dust_out     = e_dust_out_p,   &
           e_bb_out       = e_bb_out_p,  num_e_bb_in    = num_e_bb_in,                   &
           num_e_bb_out   = num_e_bb_out, num_e_dust_out = num_e_dust_out,               &
           ddvel = ddvel_p, wetdep_resolved = wetdep_resolved_p,                         &
           tend_chem_settle = tend_chem_settle_p,                                        &
           do_mpas_smoke = do_mpas_smoke, do_mpas_dust   = do_mpas_dust,                 &
           hwp_method = hwp_method, hwp_alpha = hwp_alpha, wetdep_ls_opt = wetdep_ls_opt,&
           wetdep_ls_alpha        = wetdep_ls_alpha, plumerise_opt = plumerise_opt,      &
           plume_wind_eff         = plume_wind_eff,  plume_alpha   = plume_alpha,        &
           bb_emis_scale_factor   = bb_emis_scale_factor,                                &
           ebb_dcycle             = ebb_dcycle,                                          &
           drydep_opt          = drydep_opt, pm_settling = pm_settling,                  &
           add_fire_heat_flux  = add_fire_heat_flux,                                     &
           add_fire_moist_flux = add_fire_moist_flux,                                    &
           plumerisefire_frq   = plumerisefire_frq,                                      &
           dust_alpha = dust_alpha, dust_gamma = dust_gamma,                             &
           dust_drylimit_factor  = dust_drylimit_factor,                                 &
           dust_moist_correction = dust_moist_correction,                                &
           ktau    = itimestep   , dt       = dt_dyn     , dxcell   = dx_p       ,       &
           area    = area_p      ,                                                       &
           xland   = xland_p     , u10      = u10_p      , v10      = v10_p      ,       &
           ust     = ust_p       , xlat     = xlat_p     , xlong    = xlon_p     ,       &
           tskin   = tsk_p       , pblh     = hpbl_p     , t2m      = t2m_p      ,       &
           p8w     = pres2_hyd_p , dz8w     = dz_p       , z_at_w   = zgrid_p    ,       &
           p_phy   = pres_hyd_p  , t_phy    = t_p        , u_phy    = u_p        ,       &
           v_phy   = v_p         , qv       = qv_p       , vvel     = w_p        ,       &
           pi_phy  = pi_p        , rho_phy  = rho_p      , kpbl     = kpbl_p     ,       &
           nsoil   = num_soils   , smois    = smois_p    , tslb     = tslb_p     ,       &
           ivgtyp  = ivgtyp_p    , isltyp   = isltyp_p   , nlcat    = num_landcat,       &
           swdown  = swdown_p    , z0       = z0_p       , snowh    = snowh_p    ,       &
           julian  = curr_julday , rmol     = rmol_p     , raincv   = raincv_p   ,       &
           rainncv = rainncv_p   , dpt2m    = dpt2m_p    , znt      = znt_p      ,       &                              
           mavail  = mavail_p    , g        = gravity    , vegfra   = vegfra_p   ,       &
           landusef = landusef_p , cldfrac  = cldfrac_p  , ktop_deep= ktop_deep_p,       &
           cp      = cp          , rd       = R_d        , gmt      = gmt        ,       &
           ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,       &
           ims = ims , ime = ime , jms = jms , jme = jme , kms = kds , kme = kme ,       &
           its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte         )
 call mpas_timer_stop('mpas_smoke')
 
 ! Update PM2.5/PM10
 
 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    PM2_5_p(i,k,j) = 0.0_RKIND
       if (index_smoke_fine .gt. 0)  PM2_5_p(i,k,j) = PM2_5_p(i,k,j) + chem_p(i,k,j,index_smoke_fine - chemistry_start + 1)
       if (index_dust_fine  .gt. 0)  PM2_5_p(i,k,j) = PM2_5_p(i,k,j) + chem_p(i,k,j,index_dust_fine  - chemistry_start + 1)
    PM2_5_p(i,k,j) = PM2_5_p(i,k,j) * rho_p(i,k,j)

    PM10_p(i,k,j) =  0.0_RKIND
     !if (index_smoke_coarse .gt. 0) PM10_p(i,k,j) = PM10_p(i,k,j) + chem_p(i,k,j,index_smoke_coarse - chemistry_start + 1)
      if (index_dust_coarse .gt.  0) PM10_p(i,k,j) = PM10_p(i,k,j) + chem_p(i,k,j,index_dust_coarse  - chemistry_start + 1)
    PM10_p(i,k,j) = PM10_p(i,k,j) *  rho_p(i,k,j)
 enddo
 enddo
 enddo



!copy local arrays to MPAS grid:
 call smoke_to_MPAS(configs,time_lev,state,diag_physics,tend_physics,its,ite)


 end subroutine driver_smoke

!=================================================================================================================
 end module mpas_atmphys_driver_smoke
!=================================================================================================================
