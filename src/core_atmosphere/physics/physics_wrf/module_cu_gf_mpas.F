!WRF:MODEL_LAYER:PHYSICS
!

MODULE module_cu_gf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     This convective parameterization is build to attempt     !
!     a smooth transition to cloud resolving scales as proposed!
!     by Arakawa et al (2011, ACP). It currently does not use  !
!     subsidencespreading as in G3. Difference and details     !
!     will be described in a forthcoming paper by              !
!     Grell and Freitas (2013). The parameterization also      !
!     offers options to couple with aerosols. Both, the smooth !
!     transition part as well as the aerosol coupling are      !
!     experimental. While the smooth transition part is turned !
!     on, nd has been tested dow to a resolution of about 3km  !
!     the aerosol coupling is turned off.                      !
!     More clean-up as well as a direct coupling to chemistry  !
!     will follow for V3.5.1                                   !
!                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS

!-----------------------------------------------------------------------------------------------------------------

!modifications to sourcecode for MPAS:
!   * replaced variable dx with variables dxCell and areaCell, to include the impact of the mean
!     distance between cells, and the area of individual cells when using variable-resolution meshes.
!   * replaced variable ichoice with variables ichoice_deep and ichoice_shallow to include separate
!     closures for deep and shallow convection.
!     Laura D. Fowler (laura@ucar.edu) / 2014-05-29.
!   * replaced .1 in the calculation of the initial radius with .2, prior to calculating frh.
!     Laura D. Fowler (laura@ucar.edu) / 2014-06-04.
!   * added choice for closures in the shallow convection scheme.
!     Laura D. Fowler (laura@ucar.edu) / 2104-06-19.
!   * renamed the original subroutine gfdrv to cu_grell_freitas, and clean-up the initialization of
!     local and inout variables before and after the calls to cup_gf and cup_gf_sh.
!     Laura D. Fowler (laura@ucar.edu) / 2016-03-30.
!     -> removed variables that are never computed from the argument list (for instance, all the apr_* arrays).
!     -> removed the variable ktop_deep which is a duplicate of the variable htop (ktop_deep was not properly
!        initialized in the original sourcecode).
!     -> corrected the declarations of the local arrays gswi,edti_out,massi_flx, and corrected the initialization
!        of those variables accordingly.
!     -> added the initialization of the in-cloud cloud water (gdc) and in-cloud cloud ice (gdc2) mixing ratios.
!        gdc and gdc2 were not initialized in the original sourcecode which led to non-zero values for gdc when
!        temperatures are colder than tcrit.
!     -> corrected the calculation of the moisture convergence. the conversion of w to omeg now uses the air
!        density rho interpolated to w levels.
!     Laura D. Fowler (laura@ucar.edu) / 2016-04-20.


 subroutine cu_grell_freitas( &
               itimestep,dt,dxcell,areacell                     &
              ,u,v,t,w,q,p,pi,rho,dz8w,p8w                      &
              ,xland,ht,hfx,qfx,gsw,rqvften,rthften             &
              ,rthblten,rqvblten,rthraten,kpbl,xlv,cp,g,r_v     &
              ,ichoice_deep,ichoice_shallow,ishallow_g3         &
              ,htop,hbot,k22_shallow,kbcon_shallow,ktop_shallow &
              ,xmb_total,xmb_shallow,raincv,pratec,gdc,gdc2     &
              ,rthcuten,rqvcuten,rqccuten,rqicuten              &
              ,rucuten,rvcuten                                  &
              ,pbl_scheme, maxmf, qc3d, qi3d                    &
              ,ims, ime, jms, jme, kms,kme                      &
              ,ids, ide, jds, jde, kds,kde                      &
              ,its, ite, jts, jte, kts,kte)

      use module_cu_gf_deep, only: cu_gf_deep_run,neg_check
      use module_cu_gf_sh, only: cu_gf_sh_run
 implicit none

!-----------------------------------------------------------------------------------------------------------------

!autoconv, 1=old c0, 2=berry c0
 integer, parameter:: autoconv      = 1
!aeroevap, 1=old,2=?, 3=average
 integer, parameter:: aeroevap      = 1
 integer, parameter:: use_excess    = 0
 integer, parameter:: use_excess_sh = 0

 real, parameter:: ccnclean = 250.
 real, parameter:: aodccn   = 0.1
 real, parameter:: beta     = 0.02
 real, parameter:: tcrit    = 258.

!-----------------------------------------------------------------------------------------------------------------

!intent arguments:
 integer,intent(in):: ids,ide,jds,jde,kds,kde, & 
                      ims,ime,jms,jme,kms,kme, & 
                      its,ite,jts,jte,kts,kte

 integer,intent(in):: ichoice_deep,ichoice_shallow,itimestep
 integer,intent(in):: ishallow_g3

 integer,dimension(ims:ime,jms:jme ),intent(in):: kpbl

 real,intent(in):: dt
 real,intent(in):: xlv,r_v,cp,g
 real,dimension(ims:ime,jms:jme),intent(in):: areaCell,dxCell
 real,dimension(ims:ime,jms:jme),intent(in):: hfx,qfx,gsw,ht,xland,maxmf

 real,dimension(ims:ime,kms:kme,jms:jme),intent(in):: u,v,w,p,pi,q,rho,t
 real,dimension(ims:ime,kms:kme,jms:jme),intent(in):: dz8w,p8w
 real,dimension(ims:ime,kms:kme,jms:jme),intent(in):: rqvblten,rthblten,rthraten
 real,dimension(ims:ime,kms:kme,jms:jme),intent(in),optional:: rthften,rqvften
 real,dimension(ims:ime,kms:kme,jms:jme),intent(inout):: qc3d,qi3d
 character(*),intent(in):: pbl_scheme


!inout arguments:
 integer,dimension(ims:ime,jms:jme),intent(inout):: k22_shallow,kbcon_shallow,ktop_shallow

 real,dimension(ims:ime,jms:jme),intent(inout):: hbot,htop,raincv,pratec,xmb_total,xmb_shallow
 real,dimension(ims:ime,kms:kme,jms:jme),intent(inout):: rthcuten,rqvcuten,rqccuten,rqicuten
 real,dimension(ims:ime,kms:kme,jms:jme),intent(inout):: rucuten,rvcuten
 real,dimension(ims:ime,kms:kme,jms:jme),intent(inout),optional:: gdc,gdc2

!local variables:
 character(len=50),dimension(its:ite):: ierrc,ierrcm
 character(len=50),dimension(its:ite):: ierrcs

 integer:: i,j,k,n
 integer:: ipr,jpr
 integer:: itf,jtf,ktf

 integer,dimension(its:ite):: ierr,ierrs,ierrm
 integer,dimension(its:ite):: kpbli
 integer,dimension(its:ite):: kbcon,ktop,k22s,k22,kbcons,ktops,jmin,kbconm,ktopm,k22m,jminm

 real:: dp,dq,pahfs,pgeoh,pqhfl,zkhvfl,zrho,zws,psum,clwtot
 real,dimension(its:ite):: area_loc,dx_loc
 real,dimension(its:ite):: xlandi,hfxi,qfxi
 real,dimension(its:ite):: xmb,xmbm,xmbs,xmb_dumm
 real,dimension(its:ite):: ccn
 real,dimension(its:ite):: cuten,psur,pret,pretm,prets,ter11,zqexec,ztexec,pmean,umean,vmean
 real,dimension(its:ite,kts:kte):: zo,t2d,q2d,po,p2d,us,vs,qc,qi,rhoi,tn,qo,tshall,qshall
 real,dimension(its:ite,kts:kte):: outt,outq,outqc,phh,cupclw,outu,outv
 real,dimension(its:ite,kts:kte):: outtm,outqm,outqcm,cupclwm,outum,outvm
 real,dimension(its:ite,kts:kte):: outts,outqs,outqcs,cupclws,outus,outvs,dhdt

 real,dimension(its:ite,jts:jte):: gswi,edti_out,massi_flx

 real,dimension (its:ite)         :: mconv
 real,dimension(its:ite,kts:kte)  :: omeg,qcheck
 !
 ! local variables required for new GF/C3 scheme, also for chem variables, should be input at some point...
 !
 integer, parameter:: nchem = 1, nranflag = 0, do_capsuppress = 0, dicycle = 0
 integer :: imid,ishallow,ktopmax,kbconmax
 real,dimension(its:ite,kts:kte,nchem)::chem3d
 real,dimension(nchem)::fscav
 real,dimension(its:ite,nchem)::wetdpc_deep,wetdpc_mid
 real,dimension (its:ite)         :: cap_suppress_j,rand_mom,rand_vmas
 integer,  dimension (its:ite) :: csum,csum_m
 real,dimension (its:ite,4)       :: rand_clos
 logical                          :: do_smoke_transport
 !
 ! local variables on output for diagnostics for new GF/C3 scheme
 !
 real,dimension (its:ite,10)       :: forcing,forcingm
 real,dimension (its:ite)          :: edto,edtd,edtm,frh_out,frhm
 real,dimension (its:ite,kts:kte)  :: zuo,zdo,zum,zus,zdm,zdd,cnvwt,cnvwtm,cnvwts
  
 itf = min(ite,ide-1)
 ktf = min(kte,kde-1)
 jtf = min(jte,jde-1)
!write(0,*) '--- its = ', its
!write(0,*) '--- kts = ', kts
!write(0,*) '--- jts = ', jts
!write(0,*) '--- itf = ', itf
!write(0,*) '--- ktf = ', ktf
!write(0,*) '--- jtf = ', jtf

 ipr = ite
 jpr = jte

 imid     = 1
 ishallow = ishallow_g3
 if (dx_loc(its)<6500.) imid = 0
 if (trim(pbl_scheme)=="bl_mynn") ishallow = 0

 do j = jts, jte
    do i = its, ite
       hbot(i,j)      = real(kte)
       htop(i,j)      = real(kts)
       xmb_total(i,j) = 0.
       raincv(i,j)    = 0.
       pratec(i,j)    = 0.

       !shallow convection:
       k22_shallow(i,j)   = 0
       kbcon_shallow(i,j) = 0
       ktop_shallow(i,j)  = 0
       xmb_shallow(i,j)   = 0.
    enddo
 enddo

!in-cloud cloud water and cloud ice mixing ratios:
 if(present(gdc)) then
    do j = jts, jte
       do k = jts, kte
          do i = its, ite
             gdc(i,k,j) = 0.
          enddo
       enddo
   enddo
 endif
 if(present(gdc2)) then
    do j = jts, jte
       do k = jts, kte
          do i = its, ite
             gdc2(i,k,j) = 0.
          enddo
       enddo
   enddo
 endif

 j_loop: do j = jts, jtf

!initialization of local variables:
    do i = its, itf
       edti_out(i,j)   = 0.
       gswi(i,j)       = gsw(i,j)
    enddo
 

       do k = kts, ktf
          do i = its, ite
             omeg(i,k) = 0.
          enddo
       enddo
       do i = its, itf
          mconv(i) = 0.
       enddo

    do i = its, itf
       ierrc(i)    = " "
       ierr(i)     = 0
       kbcon(i)    = 0
       kbconm(i)    = 0
       jmin(i)     = 0
       jminm(i)     = 0
       ktop(i)     = 0
       ktopm(i)     = 0
       xmb(i)      = 0.
       xmbm(i)      = 0.
       k22(i)     = 0
       k22m(i)     = 0

       !shallow convection:
       ierrcs(i)   = " "
       ierrs(i)    = 0
       kbcons(i)   = 0
       ktops(i)    = 0
       k22s(i)     = 0
       xmbs(i)     = 0.
       xmb_dumm(i) = 0.
    enddo

    do i = its, itf
       dx_loc(i)   = dxCell(i,j)
       area_loc(i) = areaCell(i,j)
       ter11(i)    = max(0.,ht(i,j))
       zo(i,kts)   = ter11(i) + 0.5*dz8w(i,1,j)
       do k = kts+1, ktf
         zo(i,k) = zo(i,k-1) + 0.5*(dz8w(i,k-1,j)+dz8w(i,k,j))
       enddo
       psur(i)     = p8w(i,1,j)*.01
       kpbli(i)    = kpbl(i,j)
       xlandi(i)   = xland(i,j)
       hfxi(i)     = hfx(i,j)
       qfxi(i)     = qfx(i,j)

       ccn(i)      = 1500.

       cuten(i)    = 0.
       umean(i)    = 0.
       vmean(i)    = 0.
       pmean(i)    = 0.
       pret(i)     = 0.
       pretm(i)     = 0.
       prets(i)     = 0.
       zqexec(i)   = 0.
       ztexec(i)   = 0.
    enddo

    do k = kts, ktf
       do i = its, itf
          us(i,k)      = u(i,k,j)
          vs(i,k)      = v(i,k,j)
          rhoi(i,k)    = rho(i,k,j)
          t2d(i,k)     = t(i,k,j)
          q2d(i,k)     = q(i,k,j)
          qc(i,k)     = max(0.,qc3d(i,k,j))
          qi(i,k)     = max(0.,qi3d(i,k,j))
          if(q2d(i,k) .lt. 1.e-08) q2d(i,k) = 1.e-08

          tn(i,k)      = t2d(i,k) + (rthften(i,k,j)+rthraten(i,k,j)+rthblten(i,k,j))*pi(i,k,j)*dt
          qo(i,k)      = q2d(i,k) + (rqvften(i,k,j)+rqvblten(i,k,j))*dt
          if(tn(i,k) .lt. 200.)   tn(i,k) = t2d(i,k)
          if(qo(i,k) .lt. 1.e-08) qo(i,k) = 1.e-08

          phh(i,k)     = p(i,k,j)
          po(i,k)      = phh(i,k)*.01
          p2d(i,k)     = po(i,k)

          cupclw(i,k)  = 0.
          outq(i,k)    = 0.
          outqm(i,k)    = 0.
          outqc(i,k)   = 0.
          outu(i,k)    = 0.
          outum(i,k)    = 0.
          outus(i,k)    = 0.
          outvm(i,k)    = 0.
          outvs(i,k)    = 0.
          outv(i,k)    = 0.

          !shallow convection:
          tshall(i,k)  = t2d(i,k) + rthblten(i,k,j)*pi(i,k,j)*dt
          qshall(i,k)  = q2d(i,k) + rqvblten(i,k,j)*dt
          dhdt(i,k)    = cp*rthblten(i,k,j)*pi(i,k,j) + xlv*rqvblten(i,k,j)

          cupclws(i,k) = 0.
          outqcs(i,k)  = 0.
          outqcm(i,k)  = 0.
          outqs(i,k)   = 0.
          outt(i,k)    = 0.
          outts(i,k)   = 0.
          outtm(i,k)   = 0.
       enddo
    enddo

    !calculation of the moisture convergence:
       do k = kts+1, ktf
          do i = its, itf
             omeg(i,k) = -g*0.5*(rho(i,k,j)+rho(i,k-1,j))*w(i,k,j)
          enddo
       enddo

       do k = kts+1, ktf
          do i = its, itf
             dq = (q2d(i,k)-q2d(i,k-1))
             mconv(i) = mconv(i) + omeg(i,k)*dq/g
          enddo
       enddo
       do i = its, itf
          if(mconv(i) .lt. 0.) mconv(i) = 0.
          if((dx_loc(i)<6500.).and.(trim(pbl_scheme)=="bl_mynn").and.(maxmf(i,jts).gt.0.))ierr(i)=555
       enddo

    if(use_excess.gt.0 .or. use_excess_sh.gt.0)then
       do i = its, itf
          zrho  = 100.*psur(i)/(287.04*(t2d(i,1)*(1.+0.608*q2d(i,1))))

          !- le and h fluxes 
          pahfs = -hfx(i,j) 
          pqhfl = -qfx(i,j)/xlv 
          !- buoyancy flux (h+le)
          zkhvfl = (pahfs/1004.64+0.608*t2d(i,1)*pqhfl)/zrho
          !- height of the 1st level
          pgeoh = zo(i,1)-ht(i,j) 
          !-convective-scale velocity w*
          zws = max(0.,0.001-1.5*0.41*zkhvfl*pgeoh/t2d(i,1))

          if(zws > tiny(pgeoh)) then
            !-convective-scale velocity w*
            zws = 1.2*zws**.3333
            !- temperature excess 
            ztexec(i)     = max(-1.5*pahfs/(zrho*zws*1004.64),0.0)
            !- moisture  excess
            zqexec(i)     = max(-1.5*pqhfl/(zrho*zws),0.)
          endif
        enddo
     endif  ! use_excess

     do k = kts+1, kte-1
        do i = its, itf
           if((p2d(i,1)-p2d(i,k)).gt.150. .and. p2d(i,k).gt.300.) then
               dp = -.5*(p2d(i,k+1)-p2d(i,k-1))
               umean(i) = umean(i) + us(i,k)*dp
               vmean(i) = vmean(i) + vs(i,k)*dp
               pmean(i) = pmean(i) + dp
            endif
        enddo
     enddo

!---- CALL CUMULUS PARAMETERIZATION:
!>Driver for the deep or congestus GF routine.
!! \section general_gf_deep Grell-Freitas Deep Convection General Algorithm
!
! following variables only for diagnostic output
!
         forcing(:,:)=0.
         forcingm(:,:)=0.
         edto(:)=0.
         edtd(:)=0.
         edtm(:)=0.
         cnvwt(:,:)=0.
         cnvwts(:,:)=0.
         cnvwtm(:,:)=0.
         zuo(:,:)=0.
         zdo(:,:)=0.
         zdm(:,:)=0.
         zdd(:,:)=0.
         zum(:,:)=0.
         zus(:,:)=0.
         frh_out(:)=0.
         frhm(:)=0.
!
! if tracers need transporting and/or scavenged (will need to define fscav)
!
         fscav(:)=0.
         wetdpc_deep(:,:)=0.
         wetdpc_mid(:,:)=0.
         chem3d(:,:,:)=0.
!
! if stochastics is included
!
         rand_mom(:)=0.
         rand_vmas(:)=0.
         rand_clos(:,:)=0.
!
! parameter to track mempory
!
         csum(:)=0.
         csum_m(:)=0.


!> - Call cu_gf_deep_run() for middle GF convection
      if(imid == 1)then
       call cu_gf_deep_run(        &
               itf,ktf,its,ite, kts,kte  &
              ,dicycle       &
              ,13            &
              ,ipr           &
              ,ccn           &
              ,ccnclean      &
              ,dt            &
              ,imid          &
              ,kpbli         &
              ,dhdt          &
              ,xlandi        &
              ,zo            &
              ,forcingm      &
              ,t2d           &
              ,q2d           &
              ,ter11         &
              ,tshall        &
              ,qshall        &
              ,p2d           &
              ,psur          &
              ,us            &
              ,vs            &
              ,qc            &
              ,qi            &
              ,rhoi          &
              ,hfxi          &
              ,qfxi          &
              ,dx_loc            &
              ,mconv         &
              ,omeg          &
              ,csum_m        &
              ,cnvwtm        &
              ,zum           &
              ,zdm           & ! hli
              ,zdd           &
              ,edtm          &
              ,edtd          & ! hli
              ,xmbm          &
              ,xmb_dumm      &
              ,xmbs          &
              ,pretm         &
              ,outum         &
              ,outvm         &
              ,outtm         &
              ,outqm         &
              ,outqcm        &
              ,kbconm        &
              ,ktopm         &
              ,cupclwm       &
              ,frhm          &
              ,ierrm         &
              ,ierrcm        &
              ,nchem         &
              ,fscav         &
              ,chem3d        &
              ,wetdpc_mid    &
              ,do_smoke_transport   &
!    the following should be set to zero if not available
              ,rand_mom      & ! for stochastics mom, if temporal and spatial patterns exist
              ,rand_vmas     & ! for stochastics vertmass, if temporal and spatial patterns exist
              ,rand_clos     & ! for stochastics closures, if temporal and spatial patterns exist
              ,nranflag      & ! flag to what you want perturbed
                               ! 1 = momentum transport
                               ! 2 = normalized vertical mass flux profile
                               ! 3 = closures
                               ! more is possible, talk to developer or
                               ! implement yourself. pattern is expected to be
                               ! betwee -1 and +1
              ,do_capsuppress,cap_suppress_j &
              ,k22m          &
              ,jminm)
      call neg_check('mid',j,dt,q2d,outqm,outtm,outum,outvm,outqcm,pretm   &
                     ,its,ite,kts,kte,itf,ktf,ktopm)
             do i=its,ite
               do k=kts,kte
                 qcheck(i,k)=q2d(i,k)+outqm(i,k)*dt
               enddo
             enddo
           endif
   call cu_gf_deep_run(        &
               itf,ktf,its,ite, kts,kte  &
              ,dicycle       &  ! diurnal cycle flag
              ,ichoice_deep  &  ! choice of closure, use "0" for ensemble average
              ,ipr           &  ! this flag can be used for debugging prints
              ,ccn           &  ! not well tested yet
              ,ccnclean      &
              ,dt            &  ! dt over which forcing is applied
              ,0             &  ! flag to turn on mid level convection
              ,kpbli          &  ! level of boundary layer height
              ,dhdt          &  ! boundary layer forcing (one closure for shallow)
              ,xlandi         &  ! land mask
              ,zo            &  ! heights above surface
              ,forcing       &  ! only diagnostic
              ,t2d             &  ! t before forcing
              ,q2d             &  ! q before forcing
              ,ter11            &  ! terrain
              ,tn            &  ! t including forcing
              ,qo            &  ! q including forcing
              ,po            &  ! pressure (mb)
              ,psur          &  ! surface pressure (mb)
              ,us            &  ! u on mass points
              ,vs            &  ! v on mass points
              ,qc            &
              ,qi            &
              ,rhoi          &  ! density
              ,hfxi          &  ! w/m2, positive upward
              ,qfxi          &  ! w/m2, positive upward
              ,dx_loc        &  ! dx is grid point dependent here
              ,mconv         &  ! integrated vertical advection of moisture
              ,omeg          &  ! omega (pa/s)
              ,csum          &  ! used to implement memory, set to zero if not avail
              ,cnvwt         &  ! gfs needs this
              ,zuo           &  ! nomalized updraft mass flux
              ,zdo           &  ! nomalized downdraft mass flux
              ,zum           &  ! nomalized downdraft mass flux from mid scheme
              ,edto          &  !
              ,edtm          &  !
              ,xmb           &  ! 
              ,xmbm          &  !
              ,xmbs          &  !
              ,pret          &  !
              ,outu          &  ! momentum tendencies at mass points
              ,outv          &  !
              ,outt          &  ! temperature tendencies
              ,outq          &  ! q tendencies
              ,outqc         &  ! ql/qice tendencies
              ,kbcon         &  ! lfc of parcel from k22
              ,ktop          &  ! cloud top
              ,cupclw        &  ! used for direct coupling to radiation, but with tuning factors
              ,frh_out       &  ! fractional coverage
              ,ierr          &  ! ierr flags are error flags, used for debugging
              ,ierrc         &  ! the following should be set to zero if not available
              ,nchem         &
              ,fscav         &
              ,chem3d        &
              ,wetdpc_deep   &
              ,do_smoke_transport   &
              ,rand_mom      &  ! for stochastics mom, if temporal and spatial patterns exist
              ,rand_vmas     &  ! for stochastics vertmass, if temporal and spatial patterns exist
              ,rand_clos     &  ! for stochastics closures, if temporal and spatial patterns exist
              ,nranflag      &  ! flag to what you want perturbed
                                !! 1 = momentum transport 
                                !! 2 = normalized vertical mass flux profile
                                !! 3 = closures
                                !! more is possible, talk to developer or
                                !! implement yourself. pattern is expected to be
                                !! betwee -1 and +1
              ,do_capsuppress,cap_suppress_j    &    !         
              ,k22                              &    !
              ,jmin)                         !

      call neg_check('deep',j,dt,qcheck,outq,outt,outu,outv,outqc,pret    &
                     ,its,ite,kts,kte,itf,ktf,ktop)
!    !... shallow convection:
    if(ishallow == 1 )then
!       call cup_gf_sh( &
!               xmbs,zo,outqcs,j,aaeq,t2d,q2d,ter11,tshall,qshall,p2d,      &
!               pret,p2d,outts,outqs,dt,itimestep,psur,us,vs,tcrit,ztexec,  &
!               zqexec,ccn,ccnclean,rhoi,dx_loc,area_loc,dhdt,kpbli,kbcons, &
!               ktops,cupclws,k22s,xlandi,gswi,tscl_kf,xlv,r_v,cp,g,        &
!               ichoice_shallow,0,0,ierrs,ierrcs,autoconv,itf,jtf,ktf,      & 
!               use_excess_sh,its,ite,jts,jte,kts,kte)
             !> - Call shallow: cu_gf_sh_run()
!
        do i = its, ite
          if(ierr(i).ne.0)ierrs(i)=999
       enddo
           call cu_gf_sh_run (us,vs,                                              &
! input variables, must be supplied
                          zo,t2d,q2d,ter11,tshall,qshall,p2d,psur,dhdt,kpbli,     &
                          rhoi,hfxi,qfxi,xlandi,3,tcrit,dt,         &
! input variables. ierr should be initialized to zero or larger than zero for
! turning off shallow convection for grid points
                          zus,xmbs,kbcons,ktops,k22s,ierrs,ierrcs,                &
! output tendencies
                          outts,outqs,outqcs,outus,outvs,cnvwts,prets,cupclws,     &
! dimesnional variables
                          itf,ktf,its,ite, kts,kte,ipr)


        do i = its, ite
          xmb_shallow(i,j)   = xmbm(i)
          k22_shallow(i,j)   = k22m(i)
          kbcon_shallow(i,j) = kbconm(i)
          ktop_shallow(i,j)  = ktopm(i)
       enddo
   endif

    do i = its, ite
       if(pret(i) .gt. 0. .or. pretm(i).gt.0. .or. prets(i).gt.0.) then
          xmb_total(i,j) = xmb(i)+xmbm(i)+xmbs(i)
          pratec(i,j)    = pret(i)+pretm(i)+prets(i)
          raincv(i,j)    = (pret(i)+pretm(i)+prets(i))*dt
          ktopmax=max(ktopm(i),ktop(i),ktops(i))
          kbconmax=max(kbconm(i),kbcon(i),kbcons(i))
          
          if(ktopmax > htop(i,j) ) htop(i,j) = ktopmax + .001
          if(kbconmax < hbot(i,j)) hbot(i,j) = kbconmax + .001
       endif
    enddo

    !... always save the tendencies of potential temperature, water vapor, cloud water, and cloud ice:
    do k = kts, ktf
       do i = its, itf
          rthcuten(i,k,j) = (outts(i,k) + outt(i,k)+outtm(i,k))/pi(i,k,j)
          rqvcuten(i,k,j) = outqs(i,k) + outq(i,k)+outqm(i,k)
          rucuten(i,k,j) = outu(i,k)+outum(i,k)+outus(i,k)
          rvcuten(i,k,j) = outv(i,k)+outvm(i,k)+outvs(i,k)

          if(t2d(i,k) .lt. tcrit) then
             rqccuten(i,k,j) = 0.
             rqicuten(i,k,j) = outqcs(i,k) + outqc(i,k) + outqcm(i,k)
             if(present(gdc2)) gdc2(i,k,j) = cupclws(i,k) + cupclw(i,k) + cupclwm(i,k)
          else
             rqicuten(i,k,j) = 0.
             rqccuten(i,k,j) = outqcs(i,k) + outqc(i,k) + outqcm(i,k)
             if(present(gdc)) gdc(i,k,j) = cupclws(i,k) + cupclw(i,k) + cupclwm(i,k)
          endif
       enddo
    enddo


 enddo j_loop

 end subroutine cu_grell_freitas

!-----------------------------------------------------------------------------------------------------------------

END MODULE module_cu_gf
