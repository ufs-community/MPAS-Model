!wrf:model_layer:physics
!
module module_sf_sfcdiags_ruclsm

       integer :: targetcell= 0 !99975 !81442 !107619

contains

   subroutine sfcdiags_ruclsm(hfx,qfx,tsk,qsfc,cqs,cqs2,chs,chs2,t2,th2,q2,  &
                     t3d,qv3d,dz,rho3d,p3d,psfc2d,                           &
                     cp,r_d,rovcp,globalcells,                               &
                     ids,ide, jds,jde, kds,kde,                              &
                     ims,ime, jms,jme, kms,kme,                              &        
                     its,ite, jts,jte, kts,kte                     )
!-------------------------------------------------------------------
      implicit none
!-------------------------------------------------------------------
      integer,  intent(in )   ::        ids,ide, jds,jde, kds,kde, &
                                        ims,ime, jms,jme, kms,kme, &
                                        its,ite, jts,jte, kts,kte
      integer,   intent(in   )     ::     globalcells(ims:ime)
      real,     dimension( ims:ime, jms:jme )                    , &
                intent(in)                  ::                hfx, &
                                                              qfx, &
                                                              tsk, &
                                                             qsfc
      real,     dimension( ims:ime, jms:jme )                    , &
                intent(inout)               ::                 q2, &
                                                              th2, &
                                                               t2
      real,     dimension( ims:ime, jms:jme )                    , &
                intent(in)                  ::                     &
                                                           psfc2d, &
                                                              chs, &
                                                              cqs, &
                                                             chs2, &
                                                             cqs2
      real,    dimension( ims:ime, kms:kme, jms:jme )            , &
               intent(in   )    ::                           qv3d, &
                                                               dz, &
                                                              t3d, &
                                                              p3d, &
                                                            rho3d

      real,     intent(in   )               ::       cp,r_D,rovcp
! local vars
      integer ::  i,j
      real    ::  rho, x2m, qlev1, tempc, qsat, p2m, qsfcprox, qsfcmr, &
                  psfc, dT, dQ, fh, fac, dz1
      real,     dimension( ims:ime, jms:jme ) :: t2_alt, q2_alt

      logical :: flux

      !flux = .true. ! used in RAP and HRRR
      flux = .false. ! used in UFS (RRFS)

      do j=jts,jte
        do i=its,ite
          rho = rho3d(i,1,j)
          psfc = psfc2d(i,j)

! 2-m Temperature - T2 
          if ( flux ) then
            if(chs2(I,J).lt.1.E-5) then
            ! when stable - let 2-m temperature be equal the first atm. level temp.
              th2(i,j) = t3d(i,1,j)*(1.E5/psfc)**rovcp 
            else
              th2(i,j) = tsk(i,j)*(1.E5/psfc)**rovcp - hfx(i,j)/(rho*cp*chs2(i,j))
            endif
            t2(I,J) = th2(I,J)*(1.E-5*psfc)**rovcp
          else
            t2(I,J) = tsk(i,j) - chs(I,J)/chs2(I,J)*(tsk(i,j) - t3d(i,1,j))
          endif ! flux method
          th2(i,j) = t2(i,j)*(1.E5/psfc)**rovcp

! 2-m Water vapor mixing ratio - Q2
          qlev1 = qv3d(i,1,j)
          ! saturation check
          tempc=t3d(i,1,j)-273.15
          if (tempc .le. 0.0) then
          ! over ice
             qsat = rsif(p3d(i,1,j), t3d(i,1,j))
          else
             qsat = rslf(p3d(i,1,j), t3d(i,1,j))
          endif
          !remove oversaturation at level 1
          qlev1 = min(qsat, qlev1)

          ! Compute QSFC proxy from QFX, qlev1 and CQS
          ! Use of QSFCprox is more accurate diagnostics for densely vegetated areas, like cropland in summer
          qsfcprox = qlev1+qfx(I,J)/(rho*cqs(I,J))
          qsfcmr = qsfc(i,j)/(1.-qsfc(i,j))

          if (globalcells(i)==targetcell) then          
            print *,' globalcells(i) =',globalcells(i)
            print *,'cqs,cqs2',cqs(i,j),cqs2(i,j)
            print *,'qsfc,qsfcmr,qsfcprox,qlev1',qsfc(i,j),qsfcmr,qsfcprox,qlev1
            print *,'(qsfcprox-qsfcmr)/qsfcmr =', (qsfcprox-qsfcmr)/qsfcmr
          endif
 
          if ( flux ) then
            if(cqs2(I,J).lt.1.E-5) then
            ! - under very stable conditions use first level for 2-m mixing ratio
              q2(I,J)=qlev1
            else
              x2m = qsfcprox - qfx(i,j)/(rho*cqs2(i,j))
              q2(i,j) = x2m
            endif
          else
          ! QFX is not used
            q2(I,J) = qsfcmr - cqs(I,J)/cqs2(i,j)*(qsfcmr - qlev1)
          endif  ! flux

          if (globalcells(i)==targetcell) then
             print *,' globalcells(i) =',globalcells(i)
             print *,'q2(i,j),qlev1', q2(i,j),qlev1
             print *,'(q2-qlev1)/qlev1 =', (q2(i,j)-qlev1)/qlev1
          endif

! Alternative logarithmic diagnostics:
          dt = t3d(i,1,j) - tsk(i,j)
          dq = qlev1 - qsfcmr
          dz1= 0.5 * dz(i,1,j)
          if (dt > 0.) then
          ! stable stratification        
             fh  = min(max(1.-dt/10.,0.01), 1.0)
             !for now, set zt = 0.05
             fac = log((2.  + .05)/(0.05 + fh))/ &
                   log((dz1 + .05)/(0.05 + fh))
             t2_alt(i,j) = tsk(i,j) + fac*(t3d(i,1,j) - tsk(i,j))
          else
             !no alternatives (yet) for unstable conditions
             t2_alt(i,j) = t2(i,j)
          endif
          if (dq > 0.) then
             fh  = min(max(1.-dq/0.003,0.01), 1.0)
             !for now, set zt = 0.05
             fac = log((2.  + .05)/(0.05 + fh))/ &
                   log((dz1 + .05)/(0.05 + fh))
             q2_alt(i,j) = qsfcmr + fac*(qlev1 - qsfcmr)
          else
             !no alternatives (yet) for unstable conditions
             q2_alt(i,j) = q2(i,j)
          endif
          if (globalcells(i)==targetcell) then
             print *,' globalcells(i) =',globalcells(i)
             print *,'dq,q2(i,j),q2_alt',dq, q2(i,j),q2_alt(i,j)
             print *,'dt,t2(i,j),t2_alt',dt, t2(i,j),t2_alt(i,j)
          endif

! Use alternative diagnostics for stable stratification
          t2(i,j) = t2_alt(i,j)
          q2(i,j) = q2_alt(i,j)

          ! check that T2 values lie in the range between TSK and T at the 1st level
          x2m     = MAX(MIN(tsk(i,j),t3d(i,1,j)) , t2(i,j))
          t2(i,j) = MIN(MAX(tsk(i,j),t3d(i,1,j)) , x2m)
          th2(i,j) = t2(i,j)*(1.e5/psfc)**rovcp

          ! check that q2 values lie between qsfcmr and qlev1
          x2m     = MAX(MIN(qsfcmr,qlev1) , q2(i,j))
          q2(i,j) = MIN(MAX(qsfcmr,qlev1) , x2m)
          ! saturation check
          tempc = t2(i,j)-273.15
          if (tempc .le. 0.0) then
          ! ice and supercooled water
            qsat = rsif(psfc, t2(i,j))
          else
          ! water
            qsat = rslf(psfc, t2(i,j))
          endif
          q2(i,j) = min(qsat, q2(i,j))
          if (globalcells(i)==targetcell) then
             print *,' globalcells(i) =',globalcells(i)
             print *,'FINAL - qsfc,qsfcmr,q2(i,j),qlev1', &
                              qsfc(i,j),qsfcmr,qsfcprox,q2(i,j),qlev1
             print *,'FINAL - tsk,t2(i,j),tlev1', &
                              tsk(i,j),t2(i,j),t3d(i,1,j)
          endif


        enddo
      enddo

  end subroutine sfcdiags_ruclsm

!- saturation functions are from Thompson microphysics scheme
      real function rslf(p,t)

      implicit none
      real, intent(in):: p, t
      real:: esl,x
      real, parameter:: c0= .611583699E03
      real, parameter:: c1= .444606896E02
      real, parameter:: c2= .143177157E01
      real, parameter:: c3= .264224321E-1
      real, parameter:: c4= .299291081E-3
      real, parameter:: c5= .203154182E-5
      real, parameter:: c6= .702620698E-8
      real, parameter:: c7= .379534310E-11
      real, parameter:: c8=-.321582393E-13

      x=max(-80.,t-273.16)

      esl=c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
      rslf=.622*esl/(p-esl)

      end function rslf
!
!    Alternative
!  ; Source: Murphy and Koop, Review of the vapour pressure of ice and
!             supercooled water for atmospheric applications, Q. J. R.
!             Meteorol. Soc (2005), 131, pp. 1539-1565.
!    Psat = EXP(54.842763 - 6763.22 / T - 4.210 * ALOG(T) + 0.000367 * T
!         + TANH(0.0415 * (T - 218.8)) * (53.878 - 1331.22
!         / T - 9.44523 * ALOG(T) + 0.014025 * T))
!
!+---+-----------------------------------------------------------------+
! this function calculates the ice saturation vapor mixing ratio as a
! function of temperature and pressure
!
      real function rsif(p,t)

      implicit none
      real, intent(in):: p, t
      real:: esi,x
      real, parameter:: c0= .609868993E03
      real, parameter:: c1= .499320233E02
      real, parameter:: c2= .184672631E01
      real, parameter:: c3= .402737184E-1
      real, parameter:: c4= .565392987E-3
      real, parameter:: c5= .521693933E-5
      real, parameter:: c6= .307839583E-7
      real, parameter:: c7= .105785160E-9
      real, parameter:: c8= .161444444E-12

      x=max(-80.,t-273.16)
      esi=c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
      rsif=.622*esi/(p-esi)

      end function rsif

end module module_sf_sfcdiags_ruclsm
