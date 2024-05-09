module bl_ugwp
use ccpp_kinds,only: kind_phys 
!===============================================================================
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: bl_ugwp_run
   PUBLIC :: bl_ugwp_init
   PUBLIC :: bl_ugwp_final
   PUBLIC :: bl_ugwp_timestep_init
   PUBLIC :: bl_ugwp_timestep_final

contains
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine bl_ugwp_run(sina, cosa,                                          &
                     rublten,rvblten,                                          &
                     dtaux3d,dtauy3d,                                          &
                     dtaux3d_ls,dtauy3d_ls,dtaux3d_bl,dtauy3d_bl,              &
                     dtaux3d_ss,dtauy3d_ss,dtaux3d_fd,dtauy3d_fd,              &
                     dusfcg,dvsfcg,                                            &
                     dusfc_ls,dvsfc_ls,dusfc_bl,dvsfc_bl,                      &
                     dusfc_ss,dvsfc_ss,dusfc_fd,dvsfc_fd,                      &
                     uproj, vproj,                                             &
                     t1, q1,                                                   &
                     prsi, prsl, prslk, zl, hpbl, kpbl, br1, xland1,           &
                     var, oc1,                                                 &
                     oa2d1, oa2d2,                                             &
                     oa2d3, oa2d4,                                             &
                     ol2d1, ol2d2,                                             &
                     ol2d3, ol2d4,                                             &
                     varss, oc1ss,                                             &
                     oa2d1ss, oa2d2ss,                                         &
                     oa2d3ss, oa2d4ss,                                         &
                     ol2d1ss, ol2d2ss,                                         &
                     ol2d3ss, ol2d4ss,                                         & 
                     g_, cp_, rd_, rv_, fv_, pi_,                              &
                     dxmeter, deltim,                                          &
                     its, ite, kte, kme,                                       &
                     errmsg, errflg                                            )
!-------------------------------------------------------------------------------
!  
!  abstract : 
!  this code handles the time tendencies of u v due to the effect of  mountain 
!  induced gravity wave drag from sub-grid scale orography. this routine 
!  not only treats the traditional upper-level wave breaking due to mountain 
!  variance (alpert 1988), but also the enhanced lower-tropospheric wave 
!  breaking due to mountain convexity and asymmetry (kim and arakawa 1995). 
!  thus, in addition to the terrain height data in a model grid gox, 
!  additional 10-2d topographic statistics files are needed, including 
!  orographic standard  deviation (var), convexity (oc1), asymmetry (oa4) 
!  and ol (ol4). these data sets are prepared based on the 30 sec usgs orography
!  hong (1999). the current scheme was implmented as in hong et al.(2008)
!
!  coded by song-you hong and young-joon kim and implemented by song-you hong
!
!  program history log:
!    2014-10-01  Hyun-Joo Choi (from KIAPS)  flow-blocking drag of kim and doyle
!                              with blocked height by dividing streamline theory
!    2017-04-06  Joseph Olson (from Gert-Jan Steeneveld) added small-scale
!                    orographic grabity wave drag:
!    2017-09-15  Joseph Olson, with some bug fixes from Michael Toy: added the
!                    topographic form drag of Beljaars et al. (2004, QJRMS)
!           Activation of each component is done by specifying the integer-parameters
!           (defined below) to 0: inactive or 1: active
!                    gsl_gwd_ls = 0 or 1: large-scale
!                    gsl_gwd_bl = 0 or 1: blocking drag 
!                    gsl_gwd_ss = 0 or 1: small-scale gravity wave drag
!                    gsl_gwd_fd = 0 or 1: topographic form drag
!
!  references:
!        hong et al. (2008), wea. and forecasting
!        kim and doyle (2005), Q. J. R. Meteor. Soc.
!        kim and arakawa (1995), j. atmos. sci.
!        alpet et al. (1988), NWP conference.
!        hong (1999), NCEP office note 424.
!        steeneveld et al (2008), JAMC
!        Tsiringakis et al. (2017), Q. J. R. Meteor. Soc.
!
!  notice : comparible or lower resolution orography files than model resolution
!           are desirable in preprocess (wps) to prevent weakening of the drag
!
!  input :                                                                
!    dudt, dvdt       - non-lin tendency for u and v wind component
!    uproj, vproj     - projection-relative U and V m/sec
!    u1, v1           - zonal and meridional wind m/sec  at t0-dt
!    t1               - temperature deg k at t0-dt
!    q1               - mixing ratio at t0-dt
!    deltim           - time step (s)                                      
!    del              - positive increment of pressure across layer (pa)
!    prslk, zl, prsl, prsi    - pressure and height variables
!    hpbl             - PBL height (m)
!    kpbl             - index level of PBL top
!    xland1           - land mask (1 for land, 2 for water)
!    oa4*, ol4*, var*, oc1* - orographic statistics 
!                                                                       
!  output :
!    dudt, dvdt - wind tendency due to gwdo
!    dtaux3d_*, dtauy3d_* - diagnosed orographic gwd
!    dusfcg_*, dvsfcg_*     - gw surface stress
!
!-------------------------------------------------------------------------------
   use ccpp_kinds, only: kind_phys
   implicit none
!
   integer, parameter                                                :: kts = 1
   integer                                           , intent(in   ) :: its, ite, kte, kme
   real(kind=kind_phys)                              , intent(in   ) :: g_, pi_, rd_, rv_, fv_,&
                                                                        cp_, deltim
   real(kind=kind_phys), dimension(its:ite)          , intent(in   ) :: dxmeter
   real(kind=kind_phys), dimension(its:ite,kts:kte)  , intent(inout) :: rublten, rvblten
  
   real(kind=kind_phys), dimension(its:ite,kts:kte)  , intent(  out) :: dtaux3d, dtauy3d 
   real(kind=kind_phys), dimension(its:ite,kts:kte)  , intent(  out) :: dtaux3d_ls, dtauy3d_ls, &
                                                                        dtaux3d_bl, dtauy3d_bl, &
                                                                        dtaux3d_ss, dtauy3d_ss, &
                                                                        dtaux3d_fd, dtauy3d_fd
  
   real(kind=kind_phys), dimension(its:ite)          , intent(  out) :: dusfcg, dvsfcg
   real(kind=kind_phys), dimension(its:ite)          , intent(  out) :: dusfc_ls, dvsfc_ls, &
                                                                        dusfc_bl, dvsfc_bl, &
                                                                        dusfc_ss, dvsfc_ss, &
                                                                        dusfc_fd, dvsfc_fd
   
   real(kind=kind_phys), dimension(its:ite)          , intent(in   ) :: sina, cosa
   real(kind=kind_phys), dimension(its:ite,kts:kte)  , intent(in   ) :: uproj, vproj
   real(kind=kind_phys), dimension(its:ite,kts:kte)  , intent(in   ) :: t1, q1, prslk, zl
   real(kind=kind_phys), dimension(its:ite)          , intent(in   ) :: hpbl, xland1, br1
   integer, dimension(its:ite)                       , intent(in   ) :: kpbl
!
   real(kind=kind_phys), dimension(its:ite,kts:kte)  , intent(in   ) :: prsl
   real(kind=kind_phys), dimension(its:ite,kts:kme)  , intent(in   ) :: prsi
!
   real(kind=kind_phys), dimension(its:ite)          , intent(in   ) :: var, oc1, &
                                                                        oa2d1, oa2d2, oa2d3, oa2d4, &
                                                                        ol2d1, ol2d2, ol2d3, ol2d4
   real(kind=kind_phys), dimension(its:ite)          , intent(in   ) :: varss, oc1ss, &
                                                                        oa2d1ss, oa2d2ss, &
                                                                        oa2d3ss, oa2d4ss, &
                                                                        ol2d1ss, ol2d2ss, &
                                                                        ol2d3ss, ol2d4ss
   character(len=*)                                  , intent(  out) :: errmsg
   integer                                           , intent(  out) :: errflg

! Variables for scale-awareness:
! Small-scale GWD + turbulent form drag
   real(kind=kind_phys), parameter :: dxmin_ss = 1000., dxmax_ss = 12000.  ! min,max range of tapering (m)
! Large-scale GWD
   real(kind=kind_phys), parameter :: dxmin_ls = 3000., dxmax_ls = 13000.  ! min,max range of tapering (m)
   real(kind=kind_phys), dimension(its:ite) :: ss_taper, ls_taper  ! large-scale tapering factors (-)
!
! Variables for limiting topographic standard deviation (var)
   real(kind=kind_phys), parameter      :: varmax_ss = 50.,      &  ! varmax_ss not used
                           varmax_fd = 500.,                     &
                           beta_fd = 0.2
!
   real(kind=kind_phys)                 :: var_temp, var_temp2
!
! GSL surface drag options to regulate specific components
! Each component is tapered off automatically as a function of dx, so best to 
! keep them activated (=1).
   integer, parameter ::         &
      gsl_gwd_ls      = 1,       & ! large-scale gravity wave drag
      gsl_gwd_bl      = 1,       & ! blocking drag 
      gsl_gwd_ss      = 1,       & ! small-scale gravity wave drag (Steeneveld et al. 2008)
      gsl_gwd_fd      = 1          ! form drag (Beljaars et al. 2004, QJRMS)

!
   real(kind=kind_phys), parameter      ::  ric     = 0.25    ! critical richardson number 
   real(kind=kind_phys), parameter      ::  dw2min  = 1.
   real(kind=kind_phys), parameter      ::  rimin   = -100.
   real(kind=kind_phys), parameter      ::  bnv2min = 1.0e-5
   real(kind=kind_phys), parameter      ::  efmin   = 0.0
   real(kind=kind_phys), parameter      ::  efmax   = 10.0
   real(kind=kind_phys), parameter      ::  xl      = 4.0e4  
   real(kind=kind_phys), parameter      ::  critac  = 1.0e-5
   real(kind=kind_phys), parameter      ::  gmax    = 1.    
   real(kind=kind_phys), parameter      ::  veleps  = 1.0                                                 
   real(kind=kind_phys), parameter      ::  factop  = 0.5
   real(kind=kind_phys), parameter      ::  frc     = 1.0      
   real(kind=kind_phys), parameter      ::  ce      = 0.8     
   real(kind=kind_phys), parameter      ::  cg      = 0.5    
   integer,parameter                    ::  kpblmin = 2
!
! local variables
!
   integer                              ::  kpblmax
   integer                              ::  latd,lond
   integer                              ::  i,k,lcap,lcapp1,nwd,idir,                          &
                                            klcap,kp1,ikount,kk
!
   real(kind=kind_phys)                 ::  fdir,cs,rcsks,                                     &
                                            wdir,ti,rdz,temp,tem2,dw2,shr2,bvf2,rdelks,        &
                                            wtkbj,tem,gfobnv,hd,fro,rim,temc,tem1,efact,       &
                                            temv
!
   real(kind=kind_phys), dimension(its:ite,kts:kte)   :: dudt, dvdt
   real(kind=kind_phys), dimension(its:ite,kts:kte)   :: dtaux2d_ls, dtauy2d_ls,  &
                                                         dtaux2d_bl, dtauy2d_bl
   real(kind=kind_phys), dimension(its:ite)           :: dusfc, dvsfc
   real(kind=kind_phys)                               :: dtaux, dtauy
   logical, dimension(its:ite)                        :: ldrag, icrilv, flag,kloop1
   logical                                            :: prop_test
   logical                                            :: ugwp_diags   ! flag for outputting optional diagnostics
   real(kind=kind_phys), dimension(its:ite)           :: coefm
!                                                                       
   real(kind=kind_phys), dimension(its:ite)           :: taub, xn, yn, ubar, vbar, fr,         &
                                                         ulow, rulow, bnv, rhobar,             &
                                                         oa, ol, oass, olss,                   &
                                                         dtfac, brvf, xlinv, delks,delks1,     &
                                                         zlowtop, govrth
   real(kind=kind_phys)                               :: cleff
   real(kind=kind_phys), dimension(its:ite,kts:kte+1) :: taup
   real(kind=kind_phys), dimension(its:ite,kts:kte-1) :: velco
   real(kind=kind_phys), dimension(its:ite,kts:kte)   :: bnv2, usqj, taud_ls, taud_bl,   &
                                                         rho, vtk, vtj
   real(kind=kind_phys), dimension(its:ite,kts:kte)   :: del 
   real(kind=kind_phys), dimension(its:ite,kts:kte)   :: u1, v1
   real(kind=kind_phys), dimension(its:ite,4)         :: oa4, ol4
   real(kind=kind_phys), dimension(its:ite,4)         :: oa4ss, ol4ss
!
   real(kind=kind_phys)                               :: inv_g  ! inverse of gravity
!
   integer, dimension(its:ite)                        :: kbl, klowtop
   integer, parameter                                 :: mdir=8
   integer, dimension(mdir)                           :: nwdir
   data nwdir/6,7,5,8,2,3,1,4/
!
   real(kind=kind_phys)                               :: xnbv, tauwavex0, tauwavey0
   real(kind=kind_phys), dimension(its:ite,kts:kte)   :: utendwave, vtendwave,        &
                                                         utendform, vtendform,        &
                                                         thx, thvx
   real(kind=kind_phys)                               :: a1, a2, wsp, H_efold, hpbl2, &
                                                         tvcon
   integer                                            :: kpbl2,kvar
!
! variables for flow-blocking drag
!
   real(kind=kind_phys), parameter                    :: frmax  = 10.
   real(kind=kind_phys), parameter                    :: olmin  = 1.0e-5
   real(kind=kind_phys), parameter                    :: odmin  = 0.1 
   real(kind=kind_phys), parameter                    :: odmax  = 10. 
!
   real(kind=kind_phys)                               :: fbdcd
   real(kind=kind_phys)                               :: zblk, tautem
   real(kind=kind_phys)                               :: pe, ke
   real(kind=kind_phys)                               :: fbdpe, fbdke 
   real(kind=kind_phys), dimension(its:ite)           :: delx, dely
   real(kind=kind_phys), dimension(its:ite,4)         :: dxy4, dxy4p
   real(kind=kind_phys), dimension(4)                 :: ol4p
   real(kind=kind_phys), dimension(its:ite)           :: dxy, dxyp, olp, od
   real(kind=kind_phys), dimension(its:ite,kts:kte+1) :: taufb
!
   integer, dimension(its:ite)                        :: komax
   integer                                            :: kblk
   real(kind=kind_phys)                               :: cd
!
! misc. coefficients
!
   real(kind=kind_phys), parameter ::      &
      clf_coeff = 3.4E+07, &  ! "effective length" coeff or large-scale GWD
                              ! Note:  Can be used as "tuning knob"
      a1_coeff = 0.00026615161, &  ! Coefficient for TOFD from Beljaars et al. (2004)
      a2_coeff = 0.005363, &       !  ""
      TOFD_coeff = 0.0759
!-------------------------------------------------------------------------------
!
! constants                                                         
!                                                                       
   lcap   = kte                                                         
   lcapp1 = lcap + 1                                                 
   fdir   = mdir / (2.0*pi_)
   inv_g  = 1.0/g_ 
!
! initialize CCPP error flag and message
!
    errmsg = ''
    errflg = 0
!
!--- calculate scale-aware tapering factors
!
do i = its,ite
   if ( dxmeter(i) .ge. dxmax_ls ) then
      ls_taper(i) = 1.
   else
      if ( dxmeter(i) .le. dxmin_ls) then
         ls_taper(i) = 0.
      else
         ls_taper(i) = 0.5 * ( sin(pi_*(dxmeter(i)-0.5*(dxmax_ls+dxmin_ls))/  &
                                   (dxmax_ls-dxmin_ls)) + 1. )
      end if
   end if
enddo

! Remove ss_tapering
ss_taper(:) = 1.

!
! calculate length of grid for flow-blocking drag
!
   delx(its:ite) = dxmeter(its:ite)
   dely(its:ite) = dxmeter(its:ite)
   dxy4(its:ite,1)  = delx(its:ite)
   dxy4(its:ite,2)  = dely(its:ite)
   dxy4(its:ite,3)  = sqrt(delx(its:ite)**2. + dely(its:ite)**2.)
   dxy4(its:ite,4)  = dxy4(its:ite,3)
   dxy4p(its:ite,1) = dxy4(its:ite,2)
   dxy4p(its:ite,2) = dxy4(its:ite,1)
   dxy4p(its:ite,3) = dxy4(its:ite,4)
   dxy4p(its:ite,4) = dxy4(its:ite,3)
!
! initialize arrays, array syntax is OK for OpenMP since these are local
!                                                                       
   ldrag   = .false. ; icrilv = .false. ; flag    = .true.
   ugwp_diags = .true.   ! Make this a namelist option eventually
!
   klowtop = 0       ; kbl    = 0
!
   dtaux   = 0.      ; dtauy  = 0.      ; xn      = 0.     ; yn      = 0.  
   ubar    = 0.      ; vbar   = 0.      ; rhobar  = 0.     ; ulow    = 0.
   oa      = 0.      ; ol     = 0.      ; oass    = 0.     ; olss    = 0.
   taub    = 0.
!
   usqj    = 0.      ; bnv2   = 0.      ; vtj     = 0.     ; vtk     = 0. 
   taup    = 0.      ; taud_ls = 0.     ; taud_bl = 0.
   dtaux3d = 0.      ; dtauy3d = 0.
   dtaux2d_ls = 0.   ; dtauy2d_ls = 0.  ; dtaux2d_bl = 0.  ; dtauy2d_bl = 0.
!
   dtfac   = 1.0     ; xlinv  = 1.0/xl
!
   komax = 0
   taufb = 0.0
!
   dudt = 0.0        ; dvdt = 0.0
!
   dusfc = 0.0       ; dvsfc = 0.0
   dusfcg = 0.0      ; dvsfcg = 0.0
if ( ugwp_diags ) then
   dusfc_ls = 0.0    ; dvsfc_ls = 0.0  ; dusfc_bl = 0.0  ; dvsfc_bl = 0.0
   dusfc_ss = 0.0    ; dvsfc_ss = 0.0  ; dusfc_fd = 0.0  ; dvsfc_fd = 0.0
endif
!
   do k = kts,kte
     do i = its,ite
       vtj(i,k) = t1(i,k)  * (1.+fv_*q1(i,k))
       vtk(i,k) = vtj(i,k) / prslk(i,k)

       !  Density (kg/m^3)

       rho(i,k) = 1./rd_ * prsl(i,k) / vtj(i,k)

       !  Delta p (positive) between interfaces levels (Pa)

       del(i,k) = prsi(i,k) - prsi(i,k+1)
       
       !  Earth-relative zonal and meridional winds (m/s)

       u1(i,k) = uproj(i,k)*cosa(i) - vproj(i,k)*sina(i)
       v1(i,k) = uproj(i,k)*sina(i) + vproj(i,k)*cosa(i)

     enddo
   enddo

!
   do i = its,ite
     zlowtop(i) = 2. * var(i)
   enddo
!
   do i = its,ite
     kloop1(i) = .true.
   enddo
!
   do k = kts+1,kte
     do i = its,ite
       if(zlowtop(i) .gt. 0.) then
       if (kloop1(i).and.zl(i,k)-zl(i,1).ge.zlowtop(i)) then
         klowtop(i) = k+1
         kloop1(i) = .false.
       endif
       endif
     enddo
   enddo
!
   kpblmax = kte
   do i = its,ite
     kbl(i)   = max(kpbl(i), klowtop(i))
     kbl(i)   = max(min(kbl(i),kpblmax),kpblmin)
   enddo
!
! determine the level of maximum orographic height
!
   ! komax(:) = kbl(:)
   komax(:) = klowtop(:) - 1    ! modification by NOAA/GSL March 2018
!
   do i = its,ite
     delks(i)  = 1.0 / (prsi(i,1) - prsi(i,kbl(i)))
     delks1(i) = 1.0 / (prsl(i,1) - prsl(i,kbl(i)))
   enddo
!
! compute low level averages within pbl
!
   do k = kts,kpblmax
     do i = its,ite
       if (k.lt.kbl(i)) then
         rcsks     = del(i,k) * delks(i)
         rdelks    = del(i,k)  * delks(i)
         ubar(i)   = ubar(i) + rcsks  * u1(i,k)      ! pbl  u  mean
         vbar(i)   = vbar(i) + rcsks  * v1(i,k)      ! pbl  v  mean
         rhobar(i)   = rhobar(i) + rdelks * rho(i,k)     ! pbl rho mean
       endif
     enddo
   enddo
!
! figure out low-level horizontal wind direction 
!
! nwd  1   2   3   4   5   6   7   8
! wd   w   s  sw  nw   e   n  ne  se
!
   do i = its,ite                      
     oa4(i,1) = oa2d1(i)
     oa4(i,2) = oa2d2(i)
     oa4(i,3) = oa2d3(i)
     oa4(i,4) = oa2d4(i)
     ol4(i,1) = ol2d1(i)
     ol4(i,2) = ol2d2(i)
     ol4(i,3) = ol2d3(i)
     ol4(i,4) = ol2d4(i)
     oa4ss(i,1) = oa2d1ss(i)
     oa4ss(i,2) = oa2d2ss(i)
     oa4ss(i,3) = oa2d3ss(i)
     oa4ss(i,4) = oa2d4ss(i)
     ol4ss(i,1) = ol2d1ss(i)
     ol4ss(i,2) = ol2d2ss(i)
     ol4ss(i,3) = ol2d3ss(i)
     ol4ss(i,4) = ol2d4ss(i)
     wdir  = atan2(ubar(i),vbar(i)) + pi_
     idir  = mod(nint(fdir*wdir),mdir) + 1
     nwd   = nwdir(idir)
     oa(i) = (1-2*int( (nwd-1)/4 )) * oa4(i,mod(nwd-1,4)+1)
     ol(i) = ol4(i,mod(nwd-1,4)+1) 
     ! Repeat for small-scale gwd
     oass(i)  = (1-2*int( (nwd-1)/4 )) * oa4ss(i,mod(nwd-1,4)+1)
     olss(i) = ol4ss(i,mod(nwd-1,4)+1) 
!
! compute orographic width along (ol) and perpendicular (olp) the wind direction
!
     ol4p(1) = ol4(i,2)
     ol4p(2) = ol4(i,1)
     ol4p(3) = ol4(i,4)
     ol4p(4) = ol4(i,3)
     olp(i)  = ol4p(mod(nwd-1,4)+1) 
!
! compute orographic direction (horizontal orographic aspect ratio)
!
     od(i) = olp(i)/max(ol(i),olmin)
     od(i) = min(od(i),odmax)
     od(i) = max(od(i),odmin)
!
! compute length of grid in the along(dxy) and cross(dxyp) wind directions
!
     dxy(i)  = dxy4(i,MOD(nwd-1,4)+1)
     dxyp(i) = dxy4p(i,MOD(nwd-1,4)+1)
   enddo
!
! END INITIALIZATION; BEGIN GWD CALCULATIONS:
!
if ( (gsl_gwd_ls .eq. 1).or.(gsl_gwd_bl .eq. 1) ) then

   do i = its,ite
      if ( ls_taper(i) .gt. 1.e-02 ) then  !====

!                                                                       
! saving richardson number in usqj for migwdi                       
!
         do k = kts,kte-1                                                
            ti        = 2.0 / (t1(i,k)+t1(i,k+1))                                
            rdz       = 1./(zl(i,k+1) - zl(i,k))
            tem1      = u1(i,k) - u1(i,k+1)
            tem2      = v1(i,k) - v1(i,k+1)   
            dw2       = tem1*tem1 + tem2*tem2
            shr2      = max(dw2,dw2min) * rdz * rdz
            bvf2      = g_*(g_/cp_+rdz*(vtj(i,k+1)-vtj(i,k))) * ti                
            usqj(i,k) = max(bvf2/shr2,rimin)                            
            bnv2(i,k) = 2.0*g_*rdz*(vtk(i,k+1)-vtk(i,k))/(vtk(i,k+1)+vtk(i,k))
            bnv2(i,k) = max( bnv2(i,k), bnv2min )                 
         enddo                                                             
!
! compute the "low level" or 1/3 wind magnitude (m/s)                
!                                                                                                                             
         ulow(i) = max(sqrt(ubar(i)*ubar(i) + vbar(i)*vbar(i)), 1.0)
         rulow(i) = 1./ulow(i)                                                          
!
         do k = kts,kte-1                                  
            velco(i,k) = 0.5 * ((u1(i,k)+u1(i,k+1)) * ubar(i)                  &
                              + (v1(i,k)+v1(i,k+1)) * vbar(i))                 
            velco(i,k) = velco(i,k) * rulow(i)                               
            if ((velco(i,k).lt.veleps) .and. (velco(i,k).gt.0.)) then
               velco(i,k) = veleps                                      
            endif                           
         enddo                                                             
!                                                                       
! no drag when critical level in the base layer                        
!                                                                                                                            
         ldrag(i) = velco(i,1).le.0.                                                                                              
!
!  no drag when velco.lt.0                                               
!                                                                       
         do k = kpblmin,kpblmax                             
            if (k .lt. kbl(i)) ldrag(i) = ldrag(i) .or. velco(i,k).le.0.                          
         enddo                                                             
!                                                                       
! the low level weighted average ri is stored in usqj(1,1; im)      
! the low level weighted average n**2 is stored in bnv2(1,1; im)    
! this is called bnvl2 in phy_gwd_alpert_sub not bnv2                           
! rdelks (del(k)/delks) vert ave factor so we can * instead of /    
!                                      
         wtkbj     = (prsl(i,1)-prsl(i,2)) * delks1(i)
         bnv2(i,1) = wtkbj * bnv2(i,1)                                
         usqj(i,1) = wtkbj * usqj(i,1)          
!
         do k = kpblmin,kpblmax            
            if (k .lt. kbl(i)) then
               rdelks    = (prsl(i,k)-prsl(i,k+1)) * delks1(i)
               bnv2(i,1) = bnv2(i,1) + bnv2(i,k) * rdelks
               usqj(i,1) = usqj(i,1) + usqj(i,k) * rdelks
            endif                                 
         enddo                                                             
!                    
         ldrag(i) = ldrag(i) .or. bnv2(i,1).le.0.0                         
         ldrag(i) = ldrag(i) .or. ulow(i).eq.1.0                           
         ldrag(i) = ldrag(i) .or. var(i) .le. 0.0                
!                                                                       
! set all ri low level values to the low level value          
!                                                                       
         do k = kpblmin,kpblmax                         
            if (k .lt. kbl(i)) usqj(i,k) = usqj(i,1)                 
         enddo                                                             
!
         if (.not.ldrag(i))   then   
            bnv(i) = sqrt( bnv2(i,1) )                                  
            fr(i) = bnv(i)  * rulow(i) * var(i) * od(i)
            fr(i) = min(fr(i),frmax)
            xn(i)  = ubar(i) * rulow(i)
            yn(i)  = vbar(i) * rulow(i)
         endif
!
! compute the base level stress and store it in taub
! calculate enhancement factor, number of mountains & aspect        
! ratio const. use simplified relationship between standard            
! deviation & critical hgt                                          
!                                      
         if (.not. ldrag(i))   then   
            efact    = (oa(i) + 2.) ** (ce*fr(i)/frc)                         
            efact    = min( max(efact,efmin), efmax )
            !!!!!!! cleff ("effective grid length") is highly tunable parameter
            ! Based on tuning by NOAA on the GFS, the value of cleff should
            ! be inversely proportional to the square root of the grid length
            ! (dxmeter).  The larger (smaller) the grid size the smaller (larger)
            ! the value of cleff.
            cleff    = clf_coeff / sqrt(dxmeter(i))                            
            coefm(i) = (1. + ol(i)) ** (oa(i)+1.)                   
            xlinv(i) = coefm(i) / cleff
            tem      = fr(i) * fr(i) * oc1(i)
            gfobnv   = gmax * tem / ((tem + cg)*bnv(i))   
               if ( gsl_gwd_ls .NE. 0 ) then
                  taub(i)  = xlinv(i) * rhobar(i) * ulow(i) * ulow(i)          &
                             * ulow(i) * gfobnv * efact          
               else     ! We've gotten what we need for the blocking scheme
                  taub(i) = 0.0
               end if
            else                                                          
            taub(i) = 0.0                                             
            xn(i)   = 0.0                                             
            yn(i)   = 0.0                                             
         endif                                                                               

      endif   ! if ( ls_taper(i) .gt. 1.e-02 )
   enddo   ! do i = its,ite

endif   ! (gsl_gwd_ls .eq. 1).or.(gsl_gwd_bl .eq. 1)

!=========================================================
! add small-scale wavedrag for stable boundary layer
!=========================================================
   xnbv=0.
   tauwavex0=0.
   tauwavey0=0.
   utendwave=0.
   vtendwave=0.
!
if ( gsl_gwd_ss .eq. 1 ) then
   do i = its,ite
      if ( ss_taper(i) .gt. 1.e-02 ) then

   !
   ! calculating potential temperature
   !
         do k = kts,kte
            thx(i,k) = t1(i,k)/prslk(i,k)
         enddo
   !
         do k = kts,kte
            tvcon = (1.+fv_*q1(i,k))
            thvx(i,k) = thx(i,k)*tvcon
         enddo
   !
         hpbl2 = hpbl(i)+10.
         kpbl2 = kpbl(i)
         !kvar = MIN(kpbl, k-level of var)
         kvar = 1
         do k=kts+1,MAX(kpbl(i),kts+1)
!            if (zl(i,k)>2.*var(i) .or. zl(i,k)>2*varmax) then
            if (zl(i,k)>300.) then
               kpbl2 = k
               if (k == kpbl(i)) then
                  hpbl2 = hpbl(i)+10.
               else
                  hpbl2 = zl(i,k)+10.
               endif
               exit
            endif
         enddo
         if ((xland1(i)-1.5).le.0. .and. 2.*varss(i).le.hpbl(i)) then
            if (br1(i).gt.0. .and. thvx(i,kpbl2)-thvx(i,kts) > 0.) then
              ! Modify xlinv to represent wave number of "typical" small-scale topography 
!              cleff_ss    = 3. * max(dx(i),cleff_ss)
!              cleff_ss    = 10. * max(dxmax_ss,cleff_ss)
!               cleff_ss    = 0.1 * 12000.
              xlinv(i) = 0.001*pi_   ! 2km horizontal wavelength
              govrth(i)=g_/(0.5*(thvx(i,kpbl2)+thvx(i,kts)))
              xnbv=sqrt(govrth(i)*(thvx(i,kpbl2)-thvx(i,kts))/hpbl2)
!
              ! check for possibility of vertical wave propagation
              ! (avoids division by zero if u1(i,kpbl2).eq.0.)
              if (u1(i,kpbl2).eq.0.) then
                 prop_test = .true.
              elseif (abs(xnbv/u1(i,kpbl2)).gt.xlinv(i)) then
                 prop_test = .true.
              else
                 prop_test = .false.
              endif
              if (prop_test) then
                ! Remove limit on varss
                var_temp = varss(i)
                ! Note:  This is a semi-implicit treatment of the time differencing
                var_temp2 = 0.5*xnbv*xlinv(i)*(2.*var_temp)**2*rho(i,kvar)  ! this is greater than zero
                tauwavex0=var_temp2*u1(i,kvar)/(1.+var_temp2*deltim)
                tauwavex0=tauwavex0*ss_taper(i)
              else
                tauwavex0=0.
              endif
!
              ! check for possibility of vertical wave propagation
              ! (avoids division by zero if v1(i,kpbl2).eq.0.)
              if (v1(i,kpbl2).eq.0.) then
                 prop_test = .true.
              elseif (abs(xnbv/v1(i,kpbl2)).gt.xlinv(i)) then
                 prop_test = .true.
              else
                 prop_test = .false.
              endif
              if (prop_test) then
                ! Remove limit on varss
                var_temp = varss(i)
                ! Note:  This is a semi-implicit treatment of the time differencing
                var_temp2 = 0.5*xnbv*xlinv(i)*(2.*var_temp)**2*rho(i,kvar)  ! this is greater than zero
                tauwavey0=var_temp2*v1(i,kvar)/(1.+var_temp2*deltim)
                tauwavey0=tauwavey0*ss_taper(i)
              else
                tauwavey0=0.
              endif

              do k=kts,kpbl(i) !MIN(kpbl2+1,km-1)
!original
                !utendwave(i,k)=-1.*tauwavex0*2.*max((1.-zl(i,k)/hpbl(i)),0.)/hpbl(i)
                !vtendwave(i,k)=-1.*tauwavey0*2.*max((1.-zl(i,k)/hpbl(i)),0.)/hpbl(i)
!new
                utendwave(i,k)=-1.*tauwavex0*2.*max((1.-zl(i,k)/hpbl2),0.)/hpbl2
                vtendwave(i,k)=-1.*tauwavey0*2.*max((1.-zl(i,k)/hpbl2),0.)/hpbl2
!mod-to be used in HRRRv3/RAPv4
                !utendwave(i,k)=-1.*tauwavex0 * max((1.-zl(i,k)/hpbl2),0.)**2
                !vtendwave(i,k)=-1.*tauwavey0 * max((1.-zl(i,k)/hpbl2),0.)**2
              enddo
            endif
         endif

         do k = kts,kte
            dudt(i,k)  = dudt(i,k) + utendwave(i,k)
            dvdt(i,k)  = dvdt(i,k) + vtendwave(i,k)
            dtaux3d(i,k) = dtaux3d(i,k) + utendwave(i,k)
            dtauy3d(i,k) = dtauy3d(i,k) + vtendwave(i,k)
            dusfc(i)   = dusfc(i) - inv_g * utendwave(i,k) * del(i,k)
            dvsfc(i)   = dvsfc(i) - inv_g * vtendwave(i,k) * del(i,k)
         enddo
         if ( ugwp_diags ) then
            do k = kts,kte
               dtaux3d_ss(i,k) = utendwave(i,k)
               dtauy3d_ss(i,k) = vtendwave(i,k)
               dusfc_ss(i) = dusfc_ss(i) + utendwave(i,k) * del(i,k)
               dvsfc_ss(i) = dvsfc_ss(i) + vtendwave(i,k) * del(i,k)
            enddo
         endif

      endif   ! if ( ss_taper(i) .gt. 1.e-02 ) then
   enddo   ! i = its,ite
endif  ! ( gsl_gwd_ss .eq. 1 )

if ( gsl_gwd_fd .eq. 1 ) THEN

   do i=its,its

      if ( ss_taper(i).gt.1.E-02 ) then

         utendform=0.
         vtendform=0.

         if ((xland1(i)-1.5) .le. 0.) then
            !(IH*kflt**n1)**-1 = (0.00102*0.00035**-1.9)**-1 = 0.00026615161
            var_temp = min(varss(i),varmax_fd) +                     &
                       max(0.,beta_fd*(varss(i)-varmax_fd))
            a1=a1_coeff*var_temp**2
!           a1=0.00026615161*MIN(varss(i),varmax)**2
!           a1=0.00026615161*(0.5*varss(i))**2
           ! k1**(n1-n2) = 0.003**(-1.9 - -2.8) = 0.003**0.9 = 0.005363
            a2=a1*a2_coeff
            ! Beljaars H_efold
            H_efold = 1500.
            do k=kts,kte
               wsp=sqrt(u1(i,k)**2 + v1(i,k)**2)
               ! alpha*beta*Cmd*Ccorr*2.109 = 12.*1.*0.005*0.6*2.109 = 0.0759
               var_temp = TOFD_coeff*EXP(-(zl(i,k)/H_efold)**1.5)*a2*          &
                                 zl(i,k)**(-1.2)*ss_taper(i) ! this is greater than zero
               !  Note:  This is a semi-implicit treatment of the time differencing
               !  per Beljaars et al. (2004, QJRMS)
               utendform(i,k) = - var_temp*wsp*u1(i,k)/(1. + var_temp*deltim*wsp)
               vtendform(i,k) = - var_temp*wsp*v1(i,k)/(1. + var_temp*deltim*wsp)
               !if(zl(i,k) > 4000.) exit
            enddo
         endif

         do k = kts,kte
            dudt(i,k)  = dudt(i,k) + utendform(i,k)
            dvdt(i,k)  = dvdt(i,k) + vtendform(i,k)
            dtaux3d(i,k) = dtaux3d(i,k) + utendform(i,k)
            dtauy3d(i,k) = dtauy3d(i,k) + vtendform(i,k)
            dusfc(i)   = dusfc(i) - inv_g * utendform(i,k) * del(i,k)
            dvsfc(i)   = dvsfc(i) - inv_g * vtendform(i,k) * del(i,k)
         enddo
         if ( ugwp_diags ) then
            do k = kts,kte
               dtaux3d_fd(i,k) = utendform(i,k)
               dtauy3d_fd(i,k) = vtendform(i,k)
               dusfc_fd(i) = dusfc_fd(i) + utendform(i,k) * del(i,k)
               dvsfc_fd(i) = dvsfc_fd(i) + vtendform(i,k) * del(i,k)
            enddo
         endif

      endif   ! if (ss_taper(i).gt.1.E-02)

   enddo  ! i=its,its

endif  ! if ( gsl_gwd_fd .eq. 1 )


if ( ugwp_diags ) then
   !  Finalize dusfc and dvsfc diagnostics for gsl small-scale drag components
   dusfc_ss(:) = -inv_g * dusfc_ss(:)
   dvsfc_ss(:) = -inv_g * dvsfc_ss(:)
   dusfc_fd(:) = -inv_g * dusfc_fd(:)
   dvsfc_fd(:) = -inv_g * dvsfc_fd(:)
endif



!=======================================================
! More for the large-scale gwd component
if ( gsl_gwd_ls.eq.1 ) then

   do i=its,ite
   
      if ( ls_taper(i).gt.1.e-02 ) then

!                                                                       
! now compute vertical structure of the stress.
!
         do k = kts,kpblmax
            if (k .le. kbl(i)) taup(i,k) = taub(i)
         enddo
!
         do k = kpblmin, kte-1                   ! vertical level k loop!
            kp1 = k + 1
!
! unstablelayer if ri < ric
! unstable layer if upper air vel comp along surf vel <=0 (crit lay)
! at (u-c)=0. crit layer exists and bit vector should be set (.le.)
!
            if (k .ge. kbl(i)) then
               icrilv(i) = icrilv(i) .or. ( usqj(i,k) .lt. ric)                &
                                     .or. (velco(i,k) .le. 0.0)
               brvf(i) = max(bnv2(i,k),bnv2min) ! brunt-vaisala frequency squared
               brvf(i) = sqrt(brvf(i))          ! brunt-vaisala frequency
            endif
!
            if (k .ge. kbl(i) .and. (.not. ldrag(i)))   then   
               if (.not.icrilv(i) .and. taup(i,k) .gt. 0.0 ) then
                  temv = 1.0 / velco(i,k)
                  tem1 = coefm(i)/dxy(i)*(rho(i,kp1)+rho(i,k))*brvf(i)*     &
                         velco(i,k)*0.5
                  hd   = sqrt(taup(i,k) / tem1)
                  fro  = brvf(i) * hd * temv
!
! rim is the  minimum-richardson number by shutts (1985)
!
                  tem2 = sqrt(usqj(i,k))
                  tem  = 1. + tem2 * fro
                  rim  = usqj(i,k) * (1.-fro) / (tem * tem)
!
! check stability to employ the 'saturation hypothesis'
! of lindzen (1981) except at tropospheric downstream regions
!
                  if (rim .le. ric) then  ! saturation hypothesis!
                     if ((oa(i) .le. 0.).or.(kp1 .ge. kpblmin )) then
                        temc = 2.0 + 1.0 / tem2
                        hd   = velco(i,k) * (2.*sqrt(temc)-temc) / brvf(i)
                        taup(i,kp1) = tem1 * hd * hd
                     endif
                  else                    ! no wavebreaking!
                     taup(i,kp1) = taup(i,k)
                  endif
               endif
            endif    
         enddo   ! k = kpblmin, kte-1
!
         if (lcap.lt.kte) then                                               
            do klcap = lcapp1,kte                                                                                      
               taup(i,klcap) = prsi(i,klcap) / prsi(i,lcap) * taup(i,lcap)                                                       
            enddo                                                          
         endif                                                             

      endif   ! if ( ls_taper(i).gt.1.e-02 )

   enddo   ! do i=its,ite

endif  ! if ( gsl_gwd_ls.eq.1 )


!===============================================================
!COMPUTE BLOCKING COMPONENT 
!===============================================================
if ( gsl_gwd_bl .eq. 1 ) then

   do i = its,ite

      if ( ls_taper(i).gt.1.e-02 ) then

         if (.not.ldrag(i)) then
!
!------- determine the height of flow-blocking layer
!
           kblk = 0
           pe = 0.0
           do k = kte, kpblmin, -1
             if(kblk.eq.0 .and. k.le.komax(i)) then
               pe = pe + bnv2(i,k)*(zl(i,komax(i))-zl(i,k))*del(i,k)/g_/rho(i,k)
               ke = 0.5*((u1(i,k))**2.+(v1(i,k))**2.)
!
!---------- apply flow-blocking drag when pe >= ke 
!
               if(pe.ge.ke) then
                 kblk = k
                 kblk = min(kblk,kbl(i))
                 zblk = zl(i,kblk)-zl(i,kts)
               endif
             endif
           enddo
           if(kblk.ne.0) then
!
!--------- compute flow-blocking stress
!
             cd = max(2.0-1.0/od(i),0.0)
             taufb(i,kts) = 0.5 * rhobar(i) * coefm(i) / max(dxmax_ls,dxy(i))**2 * cd * dxyp(i) &
                            * olp(i) * zblk * ulow(i)**2
             tautem = taufb(i,kts)/float(kblk-kts)
             do k = kts+1, kblk
               taufb(i,k) = taufb(i,k-1) - tautem
             enddo
!
!----------sum orographic GW stress and flow-blocking stress
!
          ! taup(i,:) = taup(i,:) + taufb(i,:)   ! Keep taup and taufb separate for now
           endif
!
        endif    ! if (.not.ldrag(i))

      endif   ! if ( ls_taper(i).gt.1.e-02 )

   enddo   ! i = its,ite

endif   ! if ( gsl_gwd_bl .eq. 1 ) -- end blocking drag

!===========================================================
do i = its,ite

   if ( (gsl_gwd_ls .eq. 1 .or. gsl_gwd_bl .eq. 1) .and. (ls_taper(i) .gt. 1.E-02) ) then
      !                                                                       
      !  calculate - (g)*d(tau)/d(pressure) and deceleration terms dtaux, dtauy
      !
      do k = kts,kte                                                       
       taud_ls(i,k) = 1. * (taup(i,k+1) - taup(i,k)) * g_ / del(i,k)
       taud_bl(i,k) = 1. * (taufb(i,k+1) - taufb(i,k)) * g_ / del(i,k)
      enddo                                                             
      !                                                                       
      !  limit de-acceleration (momentum deposition ) at top to 1/2 value 
      !  the idea is some stuff must go out the 'top'                     
      !                                                                       
      do klcap = lcap,kte                                               
         taud_ls(i,klcap) = taud_ls(i,klcap) * factop
         taud_bl(i,klcap) = taud_bl(i,klcap) * factop
      enddo                                                          
      !                                                                       
      !  if the gravity wave drag would force a critical line             
      !  in the lower ksmm1 layers during the next deltim timestep,     
      !  then only apply drag until that critical line is reached.        
      !                                                                       
      do k = kts,kpblmax-1
         if (k .le. kbl(i)) then
           if((taud_ls(i,k)+taud_bl(i,k)).ne.0.)                         &
              dtfac(i) = min(dtfac(i),abs(velco(i,k)                     &
                   /(deltim*(taud_ls(i,k)+taud_bl(i,k)))))
         endif
      enddo

      do k = kts,kte                                                       
         taud_ls(i,k)  = taud_ls(i,k) * dtfac(i) * ls_taper(i)
         taud_bl(i,k)  = taud_bl(i,k) * dtfac(i) * ls_taper(i)
         ! dtaux = taud(i,k) * xn(i)
         ! dtauy = taud(i,k) * yn(i)
         dtaux2d_ls(i,k) = taud_ls(i,k) * xn(i)
         dtauy2d_ls(i,k) = taud_ls(i,k) * yn(i)
         dtaux2d_bl(i,k) = taud_bl(i,k) * xn(i)
         dtauy2d_bl(i,k) = taud_bl(i,k) * yn(i)
         dudt(i,k)  = dtaux2d_ls(i,k) + dtaux2d_bl(i,k) + dudt(i,k)
         dvdt(i,k)  = dtauy2d_ls(i,k) + dtauy2d_bl(i,k) + dvdt(i,k)
         dtaux3d(i,k) = dtaux3d(i,k) + dtaux2d_ls(i,k) + dtaux2d_bl(i,k)
         dtauy3d(i,k) = dtauy3d(i,k) + dtauy2d_ls(i,k) + dtauy2d_bl(i,k)
      enddo                                                          

      if ( ugwp_diags ) then
         do k = kts,kte
            dtaux3d_ls(i,k) = dtaux2d_ls(i,k)
            dtauy3d_ls(i,k) = dtauy2d_ls(i,k)
            dtaux3d_bl(i,k) = dtaux2d_bl(i,k)
            dtauy3d_bl(i,k) = dtauy2d_bl(i,k)
            dusfc_ls(i)  = dusfc_ls(i) + dtaux2d_ls(i,k) * del(i,k)
            dvsfc_ls(i)  = dvsfc_ls(i) + dtauy2d_ls(i,k) * del(i,k)
            dusfc_bl(i)  = dusfc_bl(i) + dtaux2d_bl(i,k) * del(i,k)
            dvsfc_bl(i)  = dvsfc_bl(i) + dtauy2d_bl(i,k) * del(i,k) 
         enddo
      endif

   endif ! ( (gsl_gwd_ls .eq. 1 .or. gsl_gwd_bl .eq. 1) .and. (ls_taper(i) .gt. 1.E-02) )

enddo    ! i = its,ite

                                                  
if ( ugwp_diags ) then
!  Finalize dusfc and dvsfc diagnoses
   do i = its,ite
      dusfc_ls(i) = -inv_g * dusfc_ls(i)
      dvsfc_ls(i) = -inv_g * dvsfc_ls(i)
      dusfc_bl(i) = -inv_g * dusfc_bl(i)
      dvsfc_bl(i) = -inv_g * dvsfc_bl(i)
   enddo
endif


!
! rotate tendencies from zonal/meridional back to model grid
!
   do k = kts,kte
      do i = its,ite
         rublten(i,k) = rublten(i,k)+dudt(i,k)*cosa(i) + dvdt(i,k)*sina(i)
         rvblten(i,k) = rvblten(i,k)-dudt(i,k)*sina(i) + dvdt(i,k)*cosa(i)

         ! NOTE:  Is it desirable to rotate diagnostics "back" to (unstructured) model grid?
         !        It's probably better to keep them in zonal/meridional directions
         ! dtaux3d(i,k) = dtaux2d(i,k)*cosa(i) + dtauy2d(i,k)*sina(i)
         ! dtauy3d(i,k) =-dtaux2d(i,k)*sina(i) + dtauy2d(i,k)*cosa(i)
     enddo
   enddo
   ! NOTE:  Is it desirable to rotate diagnostics "back" to (unstructured) model grid?
   !        It's probably better to keep them in zonal/meridional directions
   ! do i = its,ite
   !    dusfcg(i) = dusfc(i)*cosa(i) + dvsfc(i)*sina(i)
   !    dvsfcg(i) =-dusfc(i)*sina(i) + dvsfc(i)*cosa(i)
   ! enddo



   return                                                            
   end subroutine bl_ugwp_run

!-------------------------------------------------------------------------------
  subroutine bl_ugwp_init (errmsg, errflg)

    character(len=*),        intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    ! This routine currently does nothing

    errmsg = ''
    errflg = 0

  end subroutine bl_ugwp_init

!-------------------------------------------------------------------------------
  subroutine bl_ugwp_final (errmsg, errflg)

    character(len=*),        intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    ! This routine currently does nothing

    errmsg = ''
    errflg = 0

  end subroutine bl_ugwp_final

!-------------------------------------------------------------------------------
  subroutine bl_ugwp_timestep_init (errmsg, errflg)

    character(len=*),        intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    ! This routine currently does nothing

    errmsg = ''
    errflg = 0

  end subroutine bl_ugwp_timestep_init

!-------------------------------------------------------------------------------
  subroutine bl_ugwp_timestep_final (errmsg, errflg)

    character(len=*),        intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    ! This routine currently does nothing

    errmsg = ''
    errflg = 0

  end subroutine bl_ugwp_timestep_final

!-------------------------------------------------------------------------------
end module bl_ugwp
