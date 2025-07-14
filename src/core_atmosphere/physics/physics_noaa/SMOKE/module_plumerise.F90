!>\file  module_plumerise.F90
!! This file is the fire plume rise driver.

 module module_plumerise

  use mpas_kind_types
!  real(RKIND),parameter :: p1000mb = 100000.  ! p at 1000mb (pascals)
!- Implementing the fire radiative power (FRP) methodology for biomass burning
!- emissions and convective energy estimation.
!- Saulo Freitas, Gabriel Pereira (INPE/UFJS, Brazil)
!- Ravan Ahmadov, Georg Grell (NOAA, USA)
!- The flag "plumerise_flag" defines the method:
!-    =1 => original method
!-    =2 => FRP based

CONTAINS
subroutine ebu_driver (      flam_frac,kfire,ebu_in,ebu,             &
                             !ebu_in_coarse,ebu_coarse,               &
                             theta_phy,q_vap,                        &   ! RAR: moist is replaced with q_vap, SRB: t_phy is repalced by theta_phy
                             rho_phy,vvel,u_phy,v_phy,pi_phy,        &   ! SRB: p_phy is replaced by pi_phy
                             wind_phy,                               &   ! SRB: added wind_phy
                             z_at_w,z,g,con_cp,con_rd,               &   ! scale_fire_emiss is part of config_flags
                             frp_inst, k_min, k_max,                 &   ! RAR:
                             wind_eff_opt,                           &
                             do_plumerise,                           &
                             kpbl_thetav, kpbl,                      &   ! SRB: added kpbl_thetav and kpbl
                             curr_secs,xlat, xlong , uspdavg2d,      &
                             hpbl2d,  alpha,                         &
                             frp_min, frp_wthreshold,zpbl_lim,uspd_lim,   &
                             ids,ide, jds,jde, kds,kde,              &
                             ims,ime, jms,jme, kms,kme,              &
                             its,ite, jts,jte, kts,kte,              & 
                             errmsg, errflg                          ) 

  use mpas_smoke_config
  use module_zero_plumegen_coms
  use module_smoke_plumerise
  IMPLICIT NONE

   REAL(RKIND), intent(in) :: frp_min, frp_wthreshold, zpbl_lim, uspd_lim
   integer, intent(in) :: kfire
    
   real(kind=RKIND), DIMENSION( ims:ime, jms:jme), INTENT(IN ) :: frp_inst         ! RAR: FRP array

   real(RKIND), DIMENSION(ims:ime,jms:jme), INTENT(IN) ::  xlat,xlong ! SRB

   integer,         dimension(ims:ime, jms:jme),     intent(in)    :: kpbl, kpbl_thetav

   logical, intent(in) :: do_plumerise

   character(*), intent(inout) :: errmsg
   integer, intent(inout) :: errflg
   INTEGER,      INTENT(IN   ) :: ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
   real(RKIND) :: curr_secs
   INTEGER,      INTENT(IN   ) :: wind_eff_opt
   REAL(RKIND), INTENT(IN)    :: alpha !  SRB: Enrainment constant for plumerise scheme
   real(kind=RKIND), DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT ) ::  ebu!, ebu_coarse
   real(kind=RKIND), INTENT(IN )  :: g, con_cp, con_rd
   real(kind=RKIND), DIMENSION( ims:ime, 1:kfire, jms:jme ), INTENT(IN )  :: ebu_in!, ebu_in_coarse
   real(kind=RKIND), DIMENSION( ims:ime, jms:jme ), INTENT(OUT ) :: flam_frac
   real(kind=RKIND), DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::   z,z_at_w,vvel,u_phy,v_phy,rho_phy,pi_phy,q_vap,theta_phy,wind_phy                     ! RAR, SRB

! Local variables...
      INTEGER :: nv, i, j, k,  kp1, kp2
      INTEGER :: icall
      INTEGER, DIMENSION(ims:ime, jms:jme), INTENT (INOUT) :: k_min, k_max      ! Min and max ver. levels for BB injection spread
      REAL, DIMENSION(ims:ime, jms:jme), INTENT (IN) :: uspdavg2d, hpbl2d ! SRB
      real(RKIND), dimension (kte) :: u_in ,v_in ,w_in ,theta_in ,pi_in, rho_phyin ,qv_in ,zmid, z_lev, uspd ! SRB
      real(kind=RKIND) :: dz_plume, cpor, con_rocp ! SRB

! MPI variables

        cpor    =con_cp/con_rd
        con_rocp=con_rd/con_cp

        if ( dbg_opt .and. (mod(int(curr_secs),1800) .eq. 0) ) then
            icall = 0
        endif

! RAR: setting to zero the ebu emissions at the levels k>1, this is necessary when the plumerise is called, so the emissions at k>1 are updated
       !do nv=1,num_ebu
          do j=jts,jte
            do k=kts,kte
               do i=its,ite
                 ebu(i,k,j)=0._RKIND
                 !ebu_coarse(i,k,j)=0._RKIND
               enddo
            enddo
          enddo
       !enddo

! RAR: new FRP based approach
! Haiqin: do_plumerise is added to the namelist options
check_pl:  IF (do_plumerise) THEN    ! if the namelist option is set for plumerise
       do j=jts,jte
          do i=its,ite

               do k=kts,kte
                  u_in(k)=  u_phy(i,k,j)
                  v_in(k)=  v_phy(i,k,j)
                  w_in(k)=  vvel(i,k,j)
                  qv_in(k)= q_vap(i,k,j)
                  pi_in(k)= pi_phy(i,k,j)
                  zmid(k)=  z(i,k,j)-z_at_w(i,kts,j)
                  z_lev(k)= z_at_w(i,k,j)-z_at_w(i,kts,j)
                  rho_phyin(k)= rho_phy(i,k,j)
                  theta_in(k)= theta_phy(i,k,j)
                  !uspd(k)= wind_phy(i,k,j) ! SRB
               enddo

! RAR: the plume rise calculation step:
               CALL plumerise(kte,1,1,1,1,1,1,                      &
                              u_in, v_in, w_in, theta_in ,pi_in,    &
                              rho_phyin, qv_in, zmid, z_lev,        &
                              wind_eff_opt,                         &
                              frp_inst(i,j), k_min(i,j),            & 
                              k_max(i,j), dbg_opt, g, con_cp,       &
                              con_rd, cpor, errmsg, errflg,         &
                              icall, xlat(i,j), xlong(i,j),         & 
                              curr_secs, alpha, frp_min )
               if(errflg/=0) return
              
               kp1= k_min(i,j)
               kp2= k_max(i,j)   

! SRB: Adding condition for overwriting plumerise levels               
               !uspdavg=SUM(uspd(kts:kpbl(i)))/kpbl(i) !Average wind speed within the boundary layer

! SRB: Adding output
               !uspdavg2(i,j) = uspdavg
               !hpbl_thetav2(i,j) = z_lev(kpbl(i))
               
               IF (frp_inst(i,j) .le. frp_min) THEN
                  !kp1=1
                  !kp2=2
                  flam_frac(i,j)= 0._RKIND 
               ELSE IF ( (frp_inst(i,j) .le. frp_wthreshold) .AND. ( uspdavg2d(i,1) .ge. uspd_lim ) .AND. & 
                       ( hpbl2d(i,1) .gt. zpbl_lim) .AND. (wind_eff_opt .eq. 1)) THEN
                  kp1=2
                  kp2=MAX(3,NINT(real(kpbl(i,j))/3._RKIND))
                  flam_frac(i,j)=0.85_RKIND
               ELSE
                  flam_frac(i,j)=0.9_RKIND  ! kp1,2 come from the plumerise scheme
               END IF
! SRB: End modification 

               ! RAR: emission distribution
               dz_plume= z_at_w(i,kp2,j) - z_at_w(i,kp1,j)
               do k=kp1,kp2-1
                     ebu(i,k,j)=flam_frac(i,j)*ebu_in(i,kts,j)*(z_at_w(i,k+1,j)-z_at_w(i,k,j))/dz_plume
                     !ebu_coarse(i,k,j)=flam_frac(i,j)*ebu_in_coarse(i,kts,j)*(z_at_w(i,k+1,j)-z_at_w(i,k,j))/dz_plume 
               enddo
               ebu(i,kts,j)= (1._RKIND - flam_frac(i,j))* ebu_in(i,kts,j)
               !ebu_coarse(i,kts,j)= (1._RKIND - flam_frac(i,j))* ebu_in_coarse(i,kts,j)
            !   ebu(i,kts,j) = ebu_in(i,j)

               ! For output diagnostic
               k_min(i,j) = kp1
               k_max(i,j) = kp2

               IF ( dbg_opt .and. (icall .le. n_dbg_lines)  .and. (frp_inst(i,j) .ge. frp_min) ) then
                 IF ( frp_inst(i,j) .ge. 2.e+10 ) then
                 END IF
                 icall = icall + 1
               END IF
!              endif check_frp
!              icall = icall + 1
            enddo
          enddo

        ENDIF check_pl

end subroutine ebu_driver

END module module_plumerise
