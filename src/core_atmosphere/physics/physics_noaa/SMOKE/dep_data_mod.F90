!>\file  dep_data_mod.F90
!! This file contains data for the dry deposition modules.
module dep_data_mod

  use mpas_kind_types

  integer, parameter :: nvegtype = 25
  real(RKIND), dimension(nvegtype), parameter :: &
        kpart = (/500., 500., 500., 500., 500., &
                  500., 500., 500., 500., 500., &
                  500., 500., 500., 500., 500., &
                  500., 500., 500., 500., 500., &
                  500., 500., 500., 500., 500.    /)
  real(RKIND), parameter :: max_dep_vel = 0.005                   ! m/s (may need to set per species)
  real(RKIND), parameter :: dep_ref_hgt = 2.0                     ! Meters 
  real(RKIND), parameter :: pi = 3.1415926536
!  3*PI
  REAL(RKIND), PARAMETER :: threepi=3.0*pi
  real(RKIND), parameter :: gravity =  9.81
! mean gravitational acceleration [ m/sec**2 ]
  REAL(RKIND), PARAMETER :: grav=9.80622
  real(RKIND), parameter :: boltzmann = 1.3807e-16
! universal gas constant [ J/mol-K ]
  REAL(RKIND), PARAMETER :: rgasuniv=8.314510
! Avogadro's Constant [ 1/mol ]
  REAL, PARAMETER :: avo=6.0221367E23
  ! Boltzmann's Constant [ J / K ]i\
  REAL(RKIND), PARAMETER :: boltz=rgasuniv/avo
  real(RKIND), parameter :: Cb = 2., Cim = 0.4, alpha = 0.8, Cin = 2.5, vv = 0.8
  real(RKIND), parameter :: A_for = 0.1 ! forest
  real(RKIND), parameter :: A_grs = 0.2 ! grass
  real(RKIND), parameter :: A_wat = 100. ! water
  real(RKIND), parameter :: eps0_for = 0.8*0.01 ! forest
  real(RKIND), parameter :: eps0_grs = 0.4*0.01 ! grass
  real(RKIND), parameter :: eps0_wat = 0.6*0.01 ! water

  REAL(RKIND), PARAMETER :: one3=1.0/3.0
  REAL(RKIND), PARAMETER :: two3=2.0/3.0
!  SQRT( 2 )
  REAL(RKIND), PARAMETER :: sqrt2=1.4142135623731
!  SQRT( PI )
  REAL(RKIND), PARAMETER :: sqrtpi=1.7724539
  REAL(RKIND) :: karman = 0.4                             ! von Karman constant
  REAL(RKIND), PARAMETER :: conmin= 1.E-16
  REAL(RKIND), PARAMETER :: pirs=3.14159265358979324
  REAL(RKIND), PARAMETER :: f6dpi=6.0/pirs
  REAL(RKIND), PARAMETER :: f6dpim9=1.0E-9*f6dpi
  REAL(RKIND), PARAMETER :: rhosmoke = 1.4E3
  REAL(RKIND), PARAMETER :: rhodust  = 2.6E3
  REAL(RKIND), PARAMETER :: smokefac=f6dpim9/rhosmoke
  REAL(RKIND), PARAMETER :: dustfac=f6dpim9/rhodust
!  starting standard surface temperature [ K ]
  REAL(RKIND), PARAMETER :: tss0=288.15
  REAL(RKIND), PARAMETER :: sigma1 = 1.8
  REAL(RKIND), PARAMETER :: mean_diameter1 = 4.e-8
  REAL(RKIND), PARAMETER :: fact_wfa = 1.e-9*6.0/pirs*exp(4.5*log(sigma1)**2)/mean_diameter1**3
  REAL(RKIND), PARAMETER :: sginia=2.00
!  initial sigma-G for nucleimode                 
  REAL(RKIND), PARAMETER :: sginin=1.70
! initial sigma-G for coarse mode               
  REAL(RKIND), PARAMETER :: sginic=2.5
!  starting standard surface pressure [ Pa ]  
  REAL(RKIND), PARAMETER :: pss0=101325.0
! lowest particle diameter ( m )   
  REAL(RKIND), PARAMETER :: dgmin=1.0E-09
! lowest particle density ( Kg/m**3 )
  REAL(RKIND), PARAMETER :: densmin=1.0E03
! index for Aitken mode number                  
  INTEGER, PARAMETER :: vdnnuc=1
! index for accumulation mode number            
  INTEGER, PARAMETER :: vdnacc=2
! index for coarse mode number                  
  INTEGER, PARAMETER :: vdncor=3
! index for Aitken mode mass                    
  INTEGER, PARAMETER :: vdmnuc=4
! index for accumulation mode                   
  INTEGER, PARAMETER :: vdmacc=5
! index for fine mode mass (Aitken + accumulation)
  INTEGER, PARAMETER :: vdmfine=6
! index for coarse mode mass                    
  INTEGER, PARAMETER :: vdmcor=7
! index for Aitken mode number                  
  INTEGER, PARAMETER :: vsnnuc=1
! index for Accumulation mode number            
  INTEGER, PARAMETER :: vsnacc=2
! index for coarse mode number                  
  INTEGER, PARAMETER :: vsncor=3
! index for Aitken mode mass                     
  INTEGER, PARAMETER :: vsmnuc=4
! index for accumulation mode mass              
  INTEGER, PARAMETER :: vsmacc=5
! index for coarse mass                         
  INTEGER, PARAMETER :: vsmcor=6
! coarse mode exp( log^2( sigmag )/8 )  
! nuclei        **4                    
      REAL(RKIND) :: esn04
! accumulation                         
      REAL(RKIND) :: esa04
      REAL(RKIND) :: esc04
! coarse                               
! nuclei        **5                    
      REAL(RKIND) :: esn05
      REAL(RKIND) :: esa05
! accumulation                         
! nuclei        **8                    
      REAL(RKIND) :: esn08
! accumulation                         
      REAL(RKIND) :: esa08
      REAL(RKIND) :: esc08
! coarse                               
! nuclei        **9                    
      REAL(RKIND) :: esn09
      REAL(RKIND) :: esa09
! accumulation                         
! nuclei        **12                   
      REAL(RKIND) :: esn12
! accumulation                         
      REAL(RKIND) :: esa12
      REAL(RKIND) :: esc12
! coarse mode                          
! nuclei        **16                   
      REAL(RKIND) :: esn16
! accumulation                         
      REAL(RKIND) :: esa16
      REAL(RKIND) :: esc16
! coarse                               
! nuclei        **20                   
      REAL(RKIND) :: esn20
! accumulation                         
      REAL(RKIND) :: esa20
      REAL(RKIND) :: esc20
! coarse                               
! nuclei        **25                   
      REAL(RKIND) :: esn25
      REAL(RKIND) :: esa25
! accumulation                         
! nuclei        **24                   
      REAL(RKIND) :: esn24
! accumulation                         
      REAL(RKIND) :: esa24
      REAL(RKIND) :: esc24
! coarse                               
! nuclei        **28                   
      REAL(RKIND) :: esn28
! accumulation                         
      REAL(RKIND) :: esa28
      REAL(RKIND) :: esc28
! coarse                               
! nuclei        **32                   
      REAL(RKIND) :: esn32
! accumulation                         
      REAL(RKIND) :: esa32
      REAL(RKIND) :: esc32
! coarese                              
! nuclei        **36                   
      REAL(RKIND) :: esn36
! accumulation                         
      REAL(RKIND) :: esa36
      REAL(RKIND) :: esc36
! coarse                               
! nuclei        **49                   
      REAL(RKIND) :: esn49
      REAL(RKIND) :: esa49
! accumulation                         
! nuclei        **52                   
      REAL(RKIND) :: esn52
      REAL(RKIND) :: esa52
! accumulation                         
! nuclei        **64                   
      REAL(RKIND) :: esn64
! accumulation                         
      REAL(RKIND) :: esa64
      REAL(RKIND) :: esc64
! coarse                               
      REAL(RKIND) :: esn100
! nuclei        **100                  
! nuclei        **(-20)                
      REAL(RKIND) :: esnm20
! accumulation                         
      REAL(RKIND) :: esam20
      REAL(RKIND) :: escm20
! coarse                               
! nuclei        **(-32)                
      REAL(RKIND) :: esnm32
! accumulation                         
      REAL(RKIND) :: esam32
      REAL(RKIND) :: escm32
!SAM 10/08 Gaussian quadrature constants for SOA_VBS deposition numerical
!integration
  INTEGER, PARAMETER ::  NGAUSdv= 7   ! Number of Gaussian Quadrature Points
  REAL(RKIND) :: xxlsgn, xxlsga, xxlsgc
  REAL(RKIND) :: Y_GQ(NGAUSdv), WGAUS(NGAUSdv)
end module dep_data_mod
