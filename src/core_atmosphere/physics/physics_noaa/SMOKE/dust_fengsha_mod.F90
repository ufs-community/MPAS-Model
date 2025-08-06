!>\file  dust_fengsha_mod.F90
!! This file contains the FENGSHA dust scheme.

module dust_fengsha_mod
!
!  This module developed by Barry Baker (NOAA ARL)
!  For serious questions contact barry.baker@noaa.gov
!
!  07/16/2019 - Adapted for NUOPC/GOCART, R. Montuoro
!  02/01/2020 - Adapted for FV3/CCPP, Haiqin Li

  use mpas_kind_types
  use dust_data_mod
  use mpas_smoke_init, only: p_dust_fine, p_dust_coarse

  implicit none

  private

  public :: gocart_dust_fengsha_driver

contains

  subroutine gocart_dust_fengsha_driver(dt,              &
       chem,rho_phy,smois,stemp,p8w,                     &
       isltyp,snowh,xland,area,g,                        &
       ust,znt,clay,sand,uthr,uthr_sg,                   &
       albedo_drag,feff,sep,                             & 
       e_dust_out, num_e_dust_out, &
       index_e_dust_out_dust_fine, index_e_dust_out_dust_coarse,   &
       num_emis_dust,num_chem,num_soil_layers,           &
       dust_alpha, dust_gamma, dust_drylimit_factor,     &
       dust_moist_correction,                            &
       ids,ide, jds,jde, kds,kde,                        &
       ims,ime, jms,jme, kms,kme,                        &
       its,ite, jts,jte, kts,kte                         )

    IMPLICIT NONE

    INTEGER,      INTENT(IN   ) ::                       &
         ids,ide, jds,jde, kds,kde,                      &
         ims,ime, jms,jme, kms,kme,                      &
         its,ite, jts,jte, kts,kte,                      &
         num_emis_dust,num_chem,num_soil_layers,         &
         num_e_dust_out, index_e_dust_out_dust_fine, index_e_dust_out_dust_coarse

    ! 2d input variables
    REAL(RKIND), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: sep     ! Sediment supply map
    REAL(RKIND), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: snowh   ! snow height (m)
    REAL(RKIND), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: xland   ! dominant land use type
    REAL(RKIND), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: area    ! area of grid cell
    REAL(RKIND), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: ust     ! friction velocity
    REAL(RKIND), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: znt     ! Surface Roughness length (m)
    REAL(RKIND), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: clay    ! Clay Fraction (-)
    REAL(RKIND), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: sand    ! Sand Fraction (-)
    REAL(RKIND), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: albedo_drag   ! Drag Partition (-)
    REAL(RKIND), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: feff          ! (New) Drag Partition
    REAL(RKIND), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: uthr    ! Dry Threshold Velocity (m/s)
    REAL(RKIND), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: uthr_sg ! (sg) Dry Threshold Velocity (m/s)
    REAL(RKIND), INTENT(IN) :: dust_alpha, dust_gamma, dust_drylimit_factor, dust_moist_correction

    INTEGER,         DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: isltyp  ! soil type

    ! 3d input variables
    REAL(RKIND), DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN) :: p8w
    REAL(RKIND), DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN) :: rho_phy
    REAL(RKIND), DIMENSION( ims:ime, kms:kme, jms:jme, 1:num_chem ), INTENT(INOUT) :: chem
    REAL(RKIND), DIMENSION( ims:ime, kms:kme, jms:jme,1:num_e_dust_out), INTENT(INOUT) :: e_dust_out
    REAL(RKIND), DIMENSION( ims:ime, 1:num_soil_layers, jms:jme ), INTENT(IN) :: smois, stemp

    !0d input variables 
    REAL(RKIND), INTENT(IN) :: dt ! time step
    REAL(RKIND), INTENT(IN) :: g  ! gravity (m/s**2)

    ! Local variables
    integer :: nmx,i,j,k,imx,jmx,lmx
    integer :: ilwi
    real(RKIND) :: airden ! air density
    REAL(RKIND) :: airmas ! dry air mass
    real(RKIND) :: dxy
    real(RKIND) :: conver,converi ! conversion values 
    real(RKIND) :: R ! local drag partition
    real(RKIND) :: ustar
    real(RKIND), DIMENSION (num_emis_dust) :: tc
    real(RKIND), DIMENSION (num_emis_dust) :: bems
    real(RKIND), DIMENSION (num_emis_dust) :: distribution
    real(RKIND), dimension (3) :: massfrac
    real(RKIND) :: erodtot
    real(RKIND) :: moist_volumetric, source_out_diag

    ! conversion values
    conver=1.e-9
    converi=1.e9

    ! Number of dust bins

    imx=1
    jmx=1
    lmx=1
    nmx=ndust

    k=kts
    do j=jts,jte
       do i=its,ite

          ! Don't do dust over water!!!
          ilwi=0
          if( (xland(i,j) - 1.5) .lt. 0.)then
             ilwi=1

             ! Total concentration at lowest model level. This is still hardcoded for 5 bins.

             !    if(config_flags%chem_opt == 2 .or. config_flags%chem_opt == 11 ) then
             !       tc(:)=1.e-16*conver
             !    else
             tc(1)=0._RKIND !chem(i,kts,j,p_dust_fine)*conver
             tc(2)=0._RKIND
             tc(3)=0._RKIND
             tc(4)=0._RKIND
             !tc(2)=chem(i,kts,j,p_dust_2)*conver
             !tc(3)=chem(i,kts,j,p_dust_3)*conver
             !tc(4)=chem(i,kts,j,p_dust_4)*conver
             tc(5)=0._RKIND !chem(i,kts,j,p_dust_coarse)*conver
             !    endif

             ! Air mass and density at lowest model level.

             airmas=-(p8w(i,kts+1,j)-p8w(i,kts,j))*area(i,j)/g
             airden=rho_phy(i,kts,j)
             ustar=ust(i,j)
             dxy=area(i,j)

             ! Mass fractions of clay, silt, and sand.
             massfrac(1)=clay(i,j)
             massfrac(2)=1._RKIND - (clay(i,j)+sand(i,j))
             massfrac(3)=sand(i,j)


             ! Total erodibility.
             
             erodtot = sep(i,j) ! SUM(erod(i,j,:))
             
             ! Don't allow roughness lengths greater than 20 cm to be lofted.
             ! This kludge accounts for land use types like urban areas and
             ! forests which would otherwise show up as high dust emitters.
             ! This is a placeholder for a more widely accepted kludge
             ! factor in the literature, which reduces lofting for rough areas.
             ! Forthcoming...

             IF (znt(i,j) .gt. 0.2_RKIND) then
                ilwi=0
             endif

             ! limit where there is lots of vegetation

             ! limit where there is snow on the ground
             if (snowh(i,j) .gt. 1.e-12_RKIND) then
                ilwi = 0
             endif

             ! Don't emit over frozen soil
             if (stemp(i,1,j) < 268._RKIND) then ! -5C
                ilwi = 0
             endif

             ! Do not allow areas with bedrock, lava, or land-ice to loft

             IF (isltyp(i,j) .eq. 15 .or. isltyp(i,j) .eq. 16. .or. &
                  isltyp(i,j) .eq. 18) then
                ilwi=0
             ENDIF
             IF (isltyp(i,j) .eq. 0)then
                ilwi=0
             endif
             if ( (ilwi == 0 ) .and. .not. ((isltyp(i,j) .eq. 15 .or. isltyp(i,j) .eq. 16. .or. &
                  isltyp(i,j) .eq. 18))) then
             endif
             if(ilwi == 0 ) cycle

             ! get drag partition
             ! FENGSHA uses the drag partition correction of MacKinnon et al 2004
             !     doi:10.1016/j.geomorph.2004.03.009
             if (dust_calcdrag .ne. 1) then
                call fengsha_drag(znt(i,j),R)
             else
                ! use the precalculated version derived from ASCAT; Prigent et al. (2012,2015)
                ! doi:10.1109/TGRS.2014.2338913 & doi:10.5194/amt-5-2703-2012
                ! pick only valid values
                if (feff(i,j) > 1.E-10_RKIND) then
                  R = real(feff(i,j), kind=RKIND)
                else
                  cycle
                endif
             endif

             ! soil moisture correction factor 
             moist_volumetric = dust_moist_correction * smois(i,2,j) 

             ! Call dust emission routine.
             bems(:) = 0._RKIND
             call source_dust(imx,jmx, lmx, nmx, dt, tc, ustar, massfrac, & 
                  erodtot, dxy, moist_volumetric, airden, airmas, bems, g, dust_alpha, dust_gamma, &
                  R, uthr(i,j),source_out_diag,dust_drylimit_factor)

             ! convert back to concentration
             
             chem(i,kts,j,p_dust_fine)  = chem(i,kts,j,p_dust_fine) + tc(1)*converi
             chem(i,kts,j,p_dust_coarse)= chem(i,kts,j,p_dust_coarse) + tc(5)*converi

             e_dust_out(i,kts,j,index_e_dust_out_dust_fine  ) = source_out_diag !bems(1)
             e_dust_out(i,kts,j,index_e_dust_out_dust_coarse) = bems(5)
             !chem(i,kts,j,p_dust_2)=tc(2)*converi
             !chem(i,kts,j,p_dust_3)=tc(3)*converi
             !chem(i,kts,j,p_dust_4)=tc(4)*converi

             ! For output diagnostics
!             emis_dust(i,1,j,p_edust1)=bems(1)
!             emis_dust(i,1,j,p_edust2)=bems(2)
!             emis_dust(i,1,j,p_edust3)=bems(3)
!             emis_dust(i,1,j,p_edust4)=bems(4)
!             emis_dust(i,1,j,p_edust5)=bems(5)
          endif
       enddo
    enddo
    !

  end subroutine gocart_dust_fengsha_driver


  subroutine source_dust(imx, jmx, lmx, nmx, dt1, tc, ustar, massfrac, &
                  erod, dxy, smois, airden, airmas, bems, g0, alpha, gamma, &
                  R, uthres,source_out_diag,dust_drylimit_factor)

    ! ****************************************************************************
    ! *  Evaluate the source of each dust particles size bin by soil emission
    ! *
    ! *  Input:
    ! *         EROD      Fraction of erodible grid cell                (-)
    ! *         smois     Volumetric  soil moisture                     (m3/m3)
    ! *         ALPHA     Constant to fudge the total emission of dust  (1/m)
    ! *         GAMMA     Tuning constant for erodibility               (-)
    ! *         DXY       Surface of each grid cell                     (m2)
    ! *         AIRMAS    Mass of air for each grid box                 (kg)
    ! *         AIRDEN    Density of air for each grid box              (kg/m3)
    ! *         USTAR     Friction velocity                             (m/s)
    ! *         DT1       Time step                                     (s)
    ! *         NMX       Number of dust bins                           (-)
    ! *         IMX       Number of I points                            (-)
    ! *         JMX       Number of J points                            (-)
    ! *         LMX       Number of L points                            (-)
    ! *         R         Drag Partition                                (-)
    ! *         UTHRES    FENGSHA Dry Threshold Velocities              (m/s)
    ! *
    ! *  Data:
    ! *         MASSFRAC  Fraction of mass in each of 3 soil classes    (-) (clay silt sand) 
    ! *         DEN_DUST  Dust density                                  (kg/m3)
    ! *         DEN_SALT  Saltation particle density                    (kg/m3)
    ! *         REFF_SALT Reference saltation particle diameter         (m)
    ! *         REFF_DUST Reference dust particle diameter              (m)
    ! *         LO_DUST   Lower diameter limits for dust bins           (m)
    ! *         UP_DUST   Upper diameter limits for dust bins           (m)
    ! *         FRAC_SALT Soil class mass fraction for saltation bins   (-)
    ! *
    ! *  Parameters:
    ! *         CMB       Constant of proportionality                   (-)
    ! *         MMD_DUST  Mass median diameter of dust                  (m)
    ! *         GSD_DUST  Geometric standard deviation of dust          (-)
    ! *         LAMBDA    Side crack propagation length                 (m)
    ! *         CV        Normalization constant                        (-)
    ! *         G0        Gravitational acceleration                    (m/s2)
    ! *
    ! *  Working:
    ! *         RHOA      Density of air in cgs                         (g/cm3)
    ! *         DS_REL    Saltation surface area distribution           (-)
    ! *         DLNDP     Dust bin width                                (-)
    ! *         EMIT      Total vertical mass flux                      (kg/m2/s)
    ! *         EMIT_VOL  Total vertical volume flux                    (m/s)
    ! *         DSRC      Mass of emitted dust               (kg/timestep/cell)
    ! *
    ! *  Output:
    ! *         TC        Total concentration of dust        (kg/kg/timestep/cell)
    ! *         BEMS      Source of each dust type           (kg/timestep/cell)
    ! *
    ! ****************************************************************************
    implicit none

    ! Input
    INTEGER,     INTENT(IN)    :: imx,jmx,lmx,nmx
    REAL(RKIND), INTENT(IN)    :: dt1
    REAL(RKIND), INTENT(IN)    :: ustar
    REAL(RKIND), INTENT(IN)    :: massfrac(3)
    REAL(RKIND), INTENT(IN)    :: erod
    REAL(RKIND), INTENT(IN)    :: dxy
    REAL(RKIND), INTENT(IN)    :: smois
    REAL(RKIND), INTENT(IN)    :: airden
    REAL(RKIND), INTENT(IN)    :: airmas
    REAL(RKIND), INTENT(IN)    :: g0
    REAL(RKIND), INTENT(IN)    :: alpha
    REAL(RKIND), INTENT(IN)    :: gamma
    REAL(RKIND), INTENT(IN)    :: R
    REAL(RKIND), INTENT(IN)    :: uthres
    REAL(RKIND), INTENT(IN)    :: dust_drylimit_factor

    ! Output
    REAL(RKIND), INTENT(INOUT) :: tc(nmx)

    REAL(RKIND), INTENT(OUT) :: source_out_diag

    ! Local Variables
    REAL(RKIND), INTENT(INOUT)   :: bems(nmx)
    
    REAL(RKIND) :: dvol(nmx)
    REAL(RKIND) :: distr_dust(nmx)
    REAL(RKIND) :: dlndp(nmx)
    REAL(RKIND) :: dsrc
    REAL(RKIND) :: dvol_tot
    REAL(RKIND) :: emit
    REAL(RKIND) :: emit_vol
    REAL(RKIND) :: rhoa
    REAL(RKIND) :: reason
    INTEGER   :: i, j, n

    ! Constant of proportionality from Marticorena et al, 1997 (unitless)
    ! Arguably more ~consistent~ fudge than alpha, which has many walnuts
    ! sprinkled throughout the literature. - GC

    REAL(RKIND), PARAMETER :: cmb=1.0_RKIND
    REAL(RKIND), PARAMETER :: kvhmax=2.0e-4_RKIND

    ! Parameters used in Kok distribution function. Advise not to play with
    ! these without the expressed written consent of someone who knows what
    ! they're doing. - GC

    REAL(RKIND), PARAMETER :: mmd_dust=3.4E-6_RKIND  ! median mass diameter (m)
    REAL(RKIND), PARAMETER :: gsd_dust=3.0_RKIND     ! geom. std deviation
    REAL(RKIND), PARAMETER :: lambda=12.0E-6_RKIND   ! crack propagation length (m)
    REAL(RKIND), PARAMETER :: cv=12.62E-6_RKIND      ! normalization constant
    REAL(RKIND), PARAMETER :: RHOSOIL=2650._RKIND


    ! calculate the total vertical dust flux 

    emit = 0._RKIND

    call DustEmissionFENGSHA(smois,massfrac(1),massfrac(3), massfrac(2), &
                                erod, R, airden, ustar, uthres, alpha, gamma, kvhmax, &
                                g0, RHOSOIL, emit,reason,dust_drylimit_factor)

    if ( isnan(emit) ) then
       write(*,*),'emit was NaN, inputs: smois,massfrac(1),massfrac(3), massfrac(2), &
                                erod, R, airden, ustar, uthres, alpha, gamma, kvhmax, &
                                g0, RHOSOIL',smois,massfrac(1),massfrac(3), massfrac(2), &
                                erod, R, airden, ustar, uthres, alpha, gamma, kvhmax, &
                                g0, RHOSOIL
    endif

    ! Now that we have the total dust emission, distribute into dust bins using
    ! lognormal distribution (Dr. Jasper Kok, in press), and
    ! calculate total mass emitted over the grid box over the timestep.
    !
    ! In calculating the Kok distribution, we assume upper and lower limits to each bin.
    ! For reff_dust=(/0.73D-6,1.4D-6,2.4D-6,4.5D-6,8.0D-6/) (default),
    ! lower limits were ASSUMED at lo_dust=(/0.1D-6,1.0D-6,1.8D-6,3.0D-6,6.0D-6/)
    ! upper limits were ASSUMED at up_dust=(/1.0D-6,1.8D-6,3.0D-6,6.0D-6,10.0D-6/)
    ! These may be changed within module_data_gocart_dust.F, but make sure it is
    ! consistent with reff_dust values.  These values were taken from the original
    ! GOCART bin configuration. We use them here to calculate dust bin width, dlndp.
    ! dVol is the volume distribution. You know...if you were wondering. GC

    dvol_tot=0._RKIND
    DO n=1,nmx
       dlndp(n)=LOG(up_dust(n)/lo_dust(n))
       dvol(n)=(2.0*reff_dust(n)/cv)*(1.+ERF(LOG(2.0*reff_dust(n)/mmd_dust)/(SQRT(2.)*LOG(gsd_dust))))*&
            EXP(-(2.0*reff_dust(n)/lambda)**3.0)*dlndp(n)
       dvol_tot=dvol_tot+dvol(n)
       ! Convert mass flux to volume flux
       !emit_vol=emit/den_dust(n) ! (m s^-1)
    END DO
    DO n=1,nmx
       distr_dust(n)=dvol(n)/dvol_tot
       !print *,"distr_dust(",n,")=",distr_dust(n)
    END DO

    ! Now distribute total vertical emission into dust bins and update concentration.

    DO n=1,nmx
       ! Calculate total mass emitted
       dsrc = emit*distr_dust(n)*dxy*dt1  ! (kg)
       IF (dsrc < 0._RKIND) dsrc = 0._RKIND
       
       ! Update dust mixing ratio at first model level.
       tc(n) = tc(n) + dsrc / airmas ! (kg/kg)
       !   bems(i,j,n) = dsrc  ! diagnostic
       !bems(i,j,n) = 1000.*dsrc/(dxy(j)*dt1) ! diagnostic (g/m2/s)
       bems(n) = 1.e+9_RKIND*dsrc/(dxy*dt1) ! diagnostic (ug/m2/s) !lzhang
       source_out_diag = reason !source_out_diag + dsrc
       
    END DO

    tc(1)=tc(1)+0.286_RKIND*tc(2)       ! This is just for RRFS-SD. DO NOT use in other models!!!
    tc(5)=0.714_RKIND*tc(2)+tc(3)+tc(4) ! This is just for RRFS-SD. DO NOT use in other models!!!
    tc(2) = 0._RKIND
    tc(3) = 0._RKIND
    tc(4) = 0._RKIND
    if ( isnan(tc(1)) ) then
       tc(1) = 0.0_RKIND
    endif
    if ( isnan(tc(5)) ) then
       tc(5) = 0.0_RKIND
    endif

    bems(1) = bems(1)+0.286_RKIND*bems(2)
    bems(5) = 0.714_RKIND*bems(2) + bems(3) + bems(4)
    bems(2) = 0._RKIND
    bems(3) = 0._RKIND
    bems(4) = 0._RKIND

  END SUBROUTINE source_dust


  subroutine fengsha_drag(z0,R)
    implicit none

    real(RKIND), intent(in) :: z0
    real(RKIND), intent(out) :: R
    real(RKIND), parameter :: z0s = 1.0e-04 !Surface roughness for ideal bare surface [m]
    ! ------------------------------------------------------------------------
    ! Function: Calculates the MacKinnon et al. 2004 Drag Partition Correction
    !
    !   R = 1.0 - log(z0 / z0s) / log( 0.7 * (12255./z0s) ** 0.8)
    !
    !--------------------------------------------------------------------------
    ! Drag partition correction. See MacKinnon et al. (2004),
    !     doi:10.1016/j.geomorph.2004.03.009
    R = 1.0 - log(z0 / z0s) / log( 0.7 * (12255./z0s) ** 0.8)

    ! Drag partition correction. See Marticorena et al. (1997),
    !     doi:10.1029/96JD02964
    !R = 1.0 - log(z0 / z0s) / log( 0.7 * (10./z0s) ** 0.8)

    return
  end subroutine fengsha_drag

  subroutine DustEmissionFENGSHA(slc, clay, sand, silt,  &
                                  sep, feff, airdens, ustar, uthrs, alpha, gamma, &
                                  kvhmax, grav, rhop, emissions,reason,dust_drylimit_factor)
    
    ! !USES:
    implicit NONE
    
! !INPUT PARAMETERS:
    REAL(RKIND), intent(in) :: slc      ! liquid water content of soil layer, volumetric fraction [1]
    REAL(RKIND), intent(in) :: clay     ! fractional clay content [1]
    REAL(RKIND), intent(in) :: sand     ! fractional sand content [1]
    REAL(RKIND), intent(in) :: silt     ! fractional silt content [1]
    REAL(RKIND), intent(in) :: sep      ! erosion map [1]
    REAL(RKIND), intent(in) :: feff    ! drag partition [1/m]
    REAL(RKIND), intent(in) :: airdens  ! air density at lowest level [kg/m^3]
    REAL(RKIND), intent(in) :: ustar    ! friction velocity [m/sec]
    REAL(RKIND), intent(in) :: uthrs    ! threshold velocity [m/2]
    REAL(RKIND), intent(in) :: alpha    ! scaling factor [1]
    REAL(RKIND), intent(in) :: gamma    ! scaling factor [1]
    REAL(RKIND), intent(in) :: kvhmax   ! max. vertical to horizontal mass flux ratio [1]
    REAL(RKIND), intent(in) :: grav     ! gravity [m/sec^2]
    REAL(RKIND), intent(in) :: rhop     ! soil class density [kg/m^3]
    real(RKIND), intent(in) :: dust_drylimit_factor
    real(RKIND), intent(out)     :: reason
    
    ! !OUTPUT PARAMETERS:
    REAL(RKIND), intent(inout) :: emissions ! binned surface emissions [kg/(m^2 sec)]
    
    ! !DESCRIPTION: Compute dust emissions using NOAA/ARL FENGSHA model
    !
    ! !REVISION HISTORY:
    !
    ! 22Feb2020 B.Baker/NOAA    - Original implementation
    ! 29Mar2021 R.Montuoro/NOAA - Refactored for process library
    ! 09Aug2022 B.Baker/NOAA    - Adapted for CCPP-Physics
    
    ! !Local Variables
    real(RKIND)                  :: alpha_grav
    real(RKIND)                  :: h
    real(RKIND)                  :: kvh
    real(RKIND)                  :: q
    real(RKIND)                  :: rustar
    real(RKIND)                  :: total_emissions
    real(RKIND)                  :: u_sum, u_thresh
    
!EOP
!-------------------------------------------------------------------------
!  Begin

!  Initialize emissions
!  --------------------
   emissions = 0._RKIND

!  Prepare scaling factor
!  ----------------------
   alpha_grav = alpha / grav

   ! Compute vertical-to-horizontal mass flux ratio
   ! ----------------------------------------------
   kvh = DustFluxV2HRatioMB95(clay, kvhmax)
   if (isnan(kvh) ) then
      write(*,*),'KVH was NaN, resetting to kvhmax'
      kvh = kvhmax
   endif

   ! Compute total emissions
   ! -----------------------
   emissions = alpha_grav * (sep ** gamma) * airdens * kvh

   !  Compute threshold wind friction velocity using drag partition
   !  -------------------------------------------------------------
   rustar = feff * ustar

   !  Now compute size-dependent total emission flux
   !  ----------------------------------------------

   if (dust_moist_opt .eq. 1) then

      ! Fecan moisture correction
      ! -------------------------
      h = moistureCorrectionFecan(slc, sand, clay, dust_drylimit_factor)
   else
      ! shao soil moisture correction 
      h = moistureCorrectionShao(slc)
   end if
   ! Adjust threshold
   ! ----------------
   u_thresh = uthrs * h
   
   u_sum = rustar + u_thresh
   
   ! Compute Horizontal Saltation Flux according to Eq (9) in Webb et al. (2020)
   ! ---------------------------------------------------------------------------
   q = max(0._RKIND, rustar - u_thresh) * u_sum * u_sum
   
   ! Distribute emissions to bins and convert to mass flux (kg s-1)
   reason = rustar - u_thresh
   if (q .gt. 0._RKIND .and. emissions .gt. 0._RKIND ) then
      write(*,*),'JLS  - emitting dust!'
   endif
   ! --------------------------------------------------------------
   emissions = emissions * q


 end subroutine DustEmissionFENGSHA
!-----------------------------------------------------------------
  real function soilMoistureConvertVol2Grav(vsoil, sandfrac)

! !USES:
    implicit NONE

! !INPUT PARAMETERS:
    REAL(RKIND), intent(in) :: vsoil       ! volumetric soil moisture fraction [1]
    REAL(RKIND), intent(in) :: sandfrac    ! fractional sand content [1]

! !DESCRIPTION: Convert soil moisture fraction from volumetric to gravimetric.
!
! !REVISION HISTORY:
!
!  02Apr2020, B.Baker/NOAA    - Original implementation
!  01Apr2020, R.Montuoro/NOAA - Adapted for GOCART process library

!  !Local Variables
    real :: vsat

!  !CONSTANTS:
    REAL(RKIND), parameter :: rhow = 1000._RKIND    ! density of water [kg m-3]
    REAL(RKIND), parameter :: rhop = 1700._RKIND    ! density of dry soil 
!EOP
!-------------------------------------------------------------------------
!  Begin...

!  Saturated volumetric water content (sand-dependent) ! [m3 m-3]
    vsat = 0.489_RKIND - 0.126_RKIND * sandfrac 
    

!  Gravimetric soil content
    soilMoistureConvertVol2Grav = 100._RKIND * (vsoil * rhow / rhop / ( 1._RKIND - vsat))

  end function soilMoistureConvertVol2Grav
!----------------------------------------------------------------
  real function moistureCorrectionFecan(slc, sand, clay, dust_drylimit_factor)

! !USES:
    implicit NONE

! !INPUT PARAMETERS:
    REAL(RKIND), intent(in) :: slc     ! liquid water content of top soil layer, volumetric fraction [1]
    REAL(RKIND), intent(in) :: sand    ! fractional sand content [1]
    REAL(RKIND), intent(in) :: clay    ! fractional clay content [1]
    REAL(RKIND), intent(in) :: dust_drylimit_factor

! !DESCRIPTION: Compute correction factor to account for Fecal soil moisture
!
! !REVISION HISTORY:
!
!  02Apr2020, B.Baker/NOAA    - Original implementation
!  01Apr2020, R.Montuoro/NOAA - Adapted for GOCART process library

!  !Local Variables
    real :: grvsoilm
    real :: drylimit

!EOP
!---------------------------------------------------------------
!  Begin...

!  Convert soil moisture from volumetric to gravimetric
    grvsoilm = soilMoistureConvertVol2Grav(slc, sand)

!  Compute fecan dry limit
    drylimit = dust_drylimit_factor * clay * (14.0 * clay + 17.0)

!  Compute soil moisture correction
    moistureCorrectionFecan = sqrt(1.0 + 1.21 * max(0., grvsoilm - drylimit)**0.68)

  end function moistureCorrectionFecan
!----------------------------------------------------------------
  real function moistureCorrectionShao(slc)

! !USES:
    implicit NONE

! !INPUT PARAMETERS:
    REAL(RKIND), intent(in) :: slc     ! liquid water content of top soil layer, volumetric fraction [1]

! !DESCRIPTION: Compute correction factor to account for Fecal soil moisture
!
! !REVISION HISTORY:
!
!  02Apr2020, B.Baker/NOAA    - Original implementation
!  01Apr2020, R.Montuoro/NOAA - Adapted for GOCART process library

!  !Local Variables
    real :: grvsoilm
    real :: drylimit

!EOP
!---------------------------------------------------------------
!  Begin...

    if (slc < 0.03) then
       moistureCorrectionShao = exp(22.7 * slc) 
    else
       moistureCorrectionShao = exp(95.3 * slc - 2.029)
    end if

  end function moistureCorrectionShao
!---------------------------------------------------------------
  real function DustFluxV2HRatioMB95(clay, kvhmax)

! !USES:
    implicit NONE

! !INPUT PARAMETERS:
    REAL(RKIND), intent(in) :: clay      ! fractional clay content [1]
    REAL(RKIND), intent(in) :: kvhmax    ! maximum flux ratio [1]

!  !CONSTANTS:
    REAL(RKIND), parameter :: clay_thresh = 0.2    ! clay fraction above which the maximum flux ratio is returned

! !DESCRIPTION: Computes the vertical-to-horizontal dust flux ratio according to
!               B.Marticorena, G.Bergametti, J.Geophys.Res., 100(D8), 164!               doi:10.1029/95JD00690
!
! !REVISION HISTORY:
!
! 22Feb2020 B.Baker/NOAA    - Original implementation
! 01Apr2021 R.Montuoro/NOAA - Adapted for GOCART process library
!
!EOP
!-------------------------------------------------------------------------
!  Begin...

    if (clay > clay_thresh) then
       DustFluxV2HRatioMB95 = kvhmax
    else
       DustFluxV2HRatioMB95 = 10.0_RKIND**(13.4_RKIND*clay-6.0_RKIND)
    end if

  end function DustFluxV2HRatioMB95
  
end module dust_fengsha_mod
