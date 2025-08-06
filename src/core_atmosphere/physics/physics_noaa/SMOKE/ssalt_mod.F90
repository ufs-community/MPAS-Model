!>\file  seas_mod.F90
!! This file contains the sea salt emission module.

module ssalt_mod

  use mpas_kind_types
  use seas_data_mod
  use seas_ngac_mod
  use mpas_smoke_init, only: p_ssalt_fine, p_ssalt_coarse

  implicit none

  integer, parameter :: CHEM_OPT_GOCART  = 300
  integer, parameter :: chem_opt  = 300

  ! -- NGAC parameters
  integer, parameter :: emission_scheme = 3    ! GEOSS 2012

  private

  public :: gocart_seasalt_driver

CONTAINS

  subroutine gocart_seasalt_driver(dt,alt,t_phy,u_phy,             &
         v_phy,num_chem,chem,rho_phy,dz8w,u10,v10,ustar,p8w,tsk,   &
         xland,xlat,xlong,area,g, &
             e_ss_out,num_e_ss_out,                                 &
             index_e_ss_out_ssalt_fine,                             &
             index_e_ss_out_ssalt_coarse,pi,                        &
         num_emis_seas,seas_opt,                                   &
         ids,ide, jds,jde, kds,kde,                                &
         ims,ime, jms,jme, kms,kme,                                &
         its,ite, jts,jte, kts,kte                                 )

     INTEGER,      INTENT(IN   ) :: num_emis_seas,num_chem,        &
                                    ids,ide, jds,jde, kds,kde,               &
                                    ims,ime, jms,jme, kms,kme,               &
                                    its,ite, jts,jte, kts,kte,seas_opt
     INTEGER, INTENT(IN) ::  num_e_ss_out,                                 &
             index_e_ss_out_ssalt_fine,                             &
             index_e_ss_out_ssalt_coarse
     REAL(RKIND), DIMENSION( ims:ime, kms:kme, jms:jme, 1:num_chem ),                 &
           INTENT(INOUT ) ::                                   chem
     REAL(RKIND), DIMENSION( ims:ime, kms:kme, jms:jme,1:num_e_ss_out),                    &
           INTENT(INOUT   ) ::                                                 &
           e_ss_out
     REAL(RKIND),  DIMENSION( ims:ime , jms:jme )                   ,               &
            INTENT(IN   ) ::                                                 &
                                                       u10,                  &
                                                       v10,                  &
                                                       ustar,tsk,            &
                                                       xland,                &
                                                       xlat,                 &
                                                       xlong,area
! JLS TODO FOR HAB
!     REAL(RKIND), DIMENSION( ims:ime, jms:jme ), 
!            INTENT(IN   ) :: bacteria_concentration
     REAL(RKIND),  DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
            INTENT(IN   ) ::                                                 &
                                                          alt,               &
                                                        t_phy,               &
                                                       dz8w,p8w,             &
                                                u_phy,v_phy,rho_phy

    REAL(RKIND), INTENT(IN   ) :: dt,g,pi
!
    integer, parameter :: p_seas_1=15
    integer, parameter :: p_seas_2=16
    integer, parameter :: p_seas_3=17
    integer, parameter :: p_seas_4=18
    integer, parameter :: p_seas_5=19

    integer, parameter :: p_eseas1=1
    integer, parameter :: p_eseas2=2
    integer, parameter :: p_eseas3=3
    integer, parameter :: p_eseas4=4
    integer, parameter :: p_eseas5=5
!
! local variables
!
    integer :: ipr,i,j,imx,jmx,lmx,n,rc
    integer,dimension (1,1) :: ilwi
    real(RKIND)               :: fsstemis, memissions, nemissions, tskin_c, ws10m
    real(RKIND) :: delp
    real(RKIND), DIMENSION (number_ss_bins) :: tc,bems
    real(RKIND), dimension (1,1) ::w10m,airmas,tskin
    real(RKIND), dimension (1) :: dxy

    real(RKIND), dimension(1,1,1) :: airmas1
    real(RKIND), dimension(1,1,1,number_ss_bins) :: tc1
    real(RKIND), dimension(1,1,number_ss_bins) :: bems1
    
!
! local parameters
!
    real(RKIND), parameter :: conver  = 1.e-9_RKIND
    real(RKIND), parameter :: converi = 1.e+9_RKIND
!
! number of dust bins
!
    imx=1
    jmx=1
    lmx=1

    ! -- original GOCART sea salt scheme
    do j = jts, jte
      do i = its, ite

        ! -- only use sea salt scheme over water
        if ( ( xland(i,j)-1.5) .gt. 0. ) then

          ! -- compute auxiliary variables
          delp = p8w(i,kts,j)-p8w(i,kts+1,j)
          if (dz8w(i,kts,j) < 12.) then
            w10m = sqrt(u_phy(i,kts,j)*u_phy(i,kts,j)+v_phy(i,kts,j)*v_phy(i,kts,j))
          else
            w10m = sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
          end if

          ilwi(1,1)=0
          tc = 0.
          tskin(1,1)=tsk(i,j)
          airmas(1,1)=area(i,j) * delp / g
          dxy(1)=area(i,j)
          ipr=0

          airmas1(1,1,1) = airmas(1,1)
          tc1(1,1,1,:) = tc
          bems1(1,1,:) = bems
          call source_ss( imx,jmx,lmx,number_ss_bins, dt, tc1, pi, ilwi, dxy, w10m, airmas1, bems1,ipr)
          tc   = tc1(1,1,1,:)
          bems = bems1(1,1,:)

          chem(i,kts,j,p_ssalt_fine)   = chem(i,kts,j,p_ssalt_fine) + &
                                          (tc(1) + 0.286*tc(2)) * converi
          chem(i,kts,j,p_ssalt_coarse) = chem(i,kts,j,p_ssalt_coarse) + &
                                          (0.714*tc(2) + tc(3) + tc(4)) * converi
          ! for output diagnostics
!           [ ug/m2/s ]
           e_ss_out(i,1,j,index_e_ss_out_ssalt_fine) = converi *( bems(1) + 0.286*bems(2) )
           e_ss_out(i,1,j,index_e_ss_out_ssalt_coarse) = converi *( 0.714*bems(2) + bems(3) + bems(4) )

! JLS TODO FOR HAB
!           if ( do_mpas_bact .and. p_bact_fine .gt. 0 ) then
!               call source_hab()
!               chem(i,kts,j,p_bact_fine) = chem(i,kts,j,p_bact_fine) + emis * converi
!           endif
 


        end if

      end do
    end do


  end subroutine gocart_seasalt_driver

  SUBROUTINE source_ss(imx,jmx,lmx,nmx, dt1, tc, pi,  &
                       ilwi, dxy, w10m, airmas, &
                       bems,ipr)

! ****************************************************************************
! *  Evaluate the source of each seasalt particles size classes  (kg/m3) 
! *  by soil emission.
! *  Input:
! *         SSALTDEN  Sea salt density                               (kg/m3)
! *         DXY       Surface of each grid cell                     (m2)
! *         NDT1      Time step                                     (s)
! *         W10m      Velocity at the anemometer level (10meters)   (m/s)
! *      
! *  Output:
! *         DSRC      Source of each sea salt bins       (kg/timestep/cell) 
! *
! *
! * Number flux density: Original formula by Monahan et al. (1986) adapted
! * by Sunling Gong (JGR 1997 (old) and GBC 2003 (new)).  The new version is
! * to better represent emission of sub-micron sea salt particles.
!
! * dFn/dr = c1*u10**c2/(r**A) * (1+c3*r**c4)*10**(c5*exp(-B**2))
! * where B = (b1 -log(r))/b2
! * see c_old, c_new, b_old, b_new below for the constants.
! * number fluxes are at 80% RH.
! *
! * To calculate the flux:
! * 1) Calculate dFn based on Monahan et al. (1986) and Gong (2003)
! * 2) Assume that wet radius r at 80% RH = dry radius r_d *frh
! * 3) Convert particles flux to mass flux :
! *    dFM/dr_d = 4/3*pi*rho_d*r_d^3 *(dr/dr_d) * dFn/dr
! *             = 4/3*pi*rho_d*r_d^3 * frh * dFn/dr
! *               where rho_p is particle density [kg/m3]
! *    The factor 1.e-18 is to convert in micro-meter r_d^3
! ****************************************************************************
   

    IMPLICIT NONE

    INTEGER, INTENT(IN)    :: nmx,imx,jmx,lmx,ipr
    INTEGER, INTENT(IN)    :: ilwi(imx,jmx)
    REAL(RKIND),    INTENT(IN)    :: dxy(jmx), w10m(imx,jmx), pi
    REAL(RKIND),    INTENT(IN)    :: airmas(imx,jmx,lmx)
    REAL(RKIND),    INTENT(INOUT) :: tc(imx,jmx,lmx,nmx)
    REAL(RKIND),    INTENT(OUT)   :: bems(imx,jmx,nmx)

    REAL(RKIND) :: c0(5), b0(2)
!  REAL(RKIND), PARAMETER :: c_old(5)=(/1.373, 3.41, 0.057, 1.05, 1.190/) 
!  REAL(RKIND), PARAMETER :: c_new(5)=(/1.373, 3.41, 0.057, 3.45, 1.607/)
    ! Change suggested by MC
    REAL(RKIND), PARAMETER :: c_old(5)=(/1.373, 3.2, 0.057, 1.05, 1.190/) 
    REAL(RKIND), PARAMETER :: c_new(5)=(/1.373, 3.2, 0.057, 3.45, 1.607/)
    REAL(RKIND), PARAMETER :: b_old(2)=(/0.380, 0.650/)
    REAL(RKIND), PARAMETER :: b_new(2)=(/0.433, 0.433/)
    REAL(RKIND), PARAMETER :: dr=5.0D-2 ! um   
    REAL(RKIND), PARAMETER :: theta=30.0
    ! Swelling coefficient frh (d rwet / d rd)
!!!  REAL(RKIND),    PARAMETER :: frh = 1.65
    REAL(RKIND),    PARAMETER :: frh = 2.d0
    LOGICAL, PARAMETER :: old=.TRUE., new=.FALSE.
    REAL(RKIND) :: rho_d, r0, r1, r, r_w, a, b, dfn, r_d, dfm, src
    INTEGER :: i, j, n, nr, ir
    REAL(RKIND) :: dt1,fudge_fac


    REAL(RKIND)    :: tcmw(nmx), ar(nmx), tcvv(nmx)
    REAL(RKIND)    :: ar_wetdep(nmx), kc(nmx)
    CHARACTER(LEN=20)     :: tcname(nmx), tcunits(nmx)
    LOGICAL               :: aerosol(nmx)


    REAL(RKIND) :: tc1(imx,jmx,lmx,nmx)
    REAL(RKIND), TARGET :: tcms(imx,jmx,lmx,nmx) ! tracer mass (kg; kgS for sulfur case)
    REAL(RKIND), TARGET :: tcgm(imx,jmx,lmx,nmx) ! g/m3

    !-----------------------------------------------------------------------  
    ! sea salt specific
    !-----------------------------------------------------------------------  
! REAL(RKIND), DIMENSION(nmx) :: ra, rb
! REAL(RKIND) :: ch_ss(nmx,12)

    !-----------------------------------------------------------------------  
    ! emissions (input)
    !-----------------------------------------------------------------------  
    REAL(RKIND) :: e_an(imx,jmx,2,nmx), e_bb(imx,jmx,nmx), &
            e_ac(imx,jmx,lmx,nmx)

    !-----------------------------------------------------------------------  
    ! diagnostics (budget)
    !-----------------------------------------------------------------------
!  ! tendencies per time step and process
!  REAL(RKIND), TARGET :: bems(imx,jmx,nmx), bdry(imx,jmx,nmx), bstl(imx,jmx,nmx)
!  REAL(RKIND), TARGET :: bwet(imx,jmx,nmx), bcnv(imx,jmx,nmx)!

!  ! integrated tendencies per process
!  REAL(RKIND), TARGET :: tems(imx,jmx,nmx), tstl(imx,jmx,nmx)
!  REAL(RKIND), TARGET :: tdry(imx,jmx,nmx), twet(imx,jmx,nmx), tcnv(imx,jmx,nmx)

    ! global mass balance per time step 
    REAL(RKIND) :: tmas0(nmx), tmas1(nmx)
    REAL(RKIND) :: dtems(nmx), dttrp(nmx), dtdif(nmx), dtcnv(nmx)
    REAL(RKIND) :: dtwet(nmx), dtdry(nmx), dtstl(nmx)
    REAL(RKIND) :: dtems2(nmx), dttrp2(nmx), dtdif2(nmx), dtcnv2(nmx)
    REAL(RKIND) :: dtwet2(nmx), dtdry2(nmx), dtstl2(nmx)

    ! detailed integrated budgets for individual emissions
    REAL(RKIND), TARGET :: ems_an(imx,jmx,nmx),    ems_bb(imx,jmx,nmx), ems_tp(imx,jmx)
    REAL(RKIND), TARGET :: ems_ac(imx,jmx,lmx,nmx)
    REAL(RKIND), TARGET :: ems_co(imx,jmx,nmx)

    ! executable statements
! decrease seasalt emissions (Colarco et al. 2010)
!
    !fudge_fac= 1. !.5
    !fudge_fac= .5 !lzhang
    fudge_fac= .25 !lzhang
!
    DO n = 1,nmx
       bems(:,:,n) = 0.0
       rho_d = den_seas(n)
       r0 = ra(n)*frh
       r1 = rb(n)*frh
       r = r0
       nr = INT((r1-r0)/dr+.001)
       DO ir = 1,nr
          r_w = r + dr*0.5
          r = r + dr
          IF (new) THEN
             a = 4.7*(1.0 + theta*r_w)**(-0.017*r_w**(-1.44))
             c0 = c_new
             b0 = b_new
          ELSE
             a = 3.0
             c0 = c_old
             b0 = b_old
          END IF
          !
          b = (b0(1) - LOG10(r_w))/b0(2)
          dfn = (c0(1)/r_w**a)*(1.0 + c0(3)*r_w**c0(4))* &
               10**(c0(5)*EXP(-(b**2)))
          
          r_d = r_w/frh*1.0D-6  ! um -> m
          dfm = 4.0/3.0*pi*r_d**3*rho_d*frh*dfn*dr*dt1 ! 3600 !dt1
          DO i = 1,imx
             DO j = 1,jmx
!              IF (water(i,j) > 0.0) THEN
                IF (ilwi(i,j) == 0) THEN
!                 src = dfm*dxy(j)*water(i,j)*w10m(i,j)**c0(2)
                   src = dfm*dxy(j)*w10m(i,j)**c0(2)
!                 src = ch_ss(n,dt(1)%mn)*dfm*dxy(j)*w10m(i,j)**c0(2)
                   tc(i,j,1,n) = tc(i,j,1,n) + fudge_fac*src/airmas(i,j,1)
                ELSE
                   src = 0.0
                END IF
                bems(i,j,n) = bems(i,j,n) + src*fudge_fac/(dxy(j)*dt1) !kg/m2/s
             END DO  ! i
          END DO ! j
       END DO ! ir
    END DO ! n

  END SUBROUTINE source_ss

end module ssalt_mod
