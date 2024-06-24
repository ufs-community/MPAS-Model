! Author: Tanya Smirnova (Tanya.Smirnova@noaa.gov)
! 15 march 2024 - RUC LSM is updated to the version used in RRFS version 1.
! Most changes are in the snow model.
!
#define iceruc_dbg_lvl 3000
!wrf:model_layer:physics
!
module module_ruc_ice

#if defined(mpas)
use mpas_atmphys_constants, g => gravity,rhowater=>rho_w, piconst => pii
use mpas_atmphys_utilities, only: physics_error_fatal,physics_message
use mpas_log, only: mpas_log_write
#define em_core 1
#define fatal_error(m) call physics_error_fatal( m )
#define wrf_at_debug_level(iceruc_dbg_lvl) .false.
#else
  use module_model_constants
  use module_wrf_error
#define fatal_error(m) call wrf_error_fatal( m )
#endif

   integer :: targetcell= 0 !45415


contains

!-----------------------------------------------------------------
   subroutine ruc_ice(spp_lsm,                              &
!-- in
#if (em_core==1)
              pattern_spp_lsm,field_sf,                     &
#endif
              dt,ktau,nsl,                                  &
#if (em_core==1)
              graupelncv,snowncv,rainncv,                   &
#endif
              zs,rainbl,frzfrac,frpcpn,rhosnf,precipfr,     &
              z3d,p8w,t3d,qv3d,qc3d,rho3d,                  & 
              glw,gsw,emiss,chklowq, chs,                   &
              flqc,flhc,alb,znt,z0,snoalb,albbck, mminlu,   &
              tbot,ivgtyp,isltyp,xland,myjpbl,              &
              iswater,isice,xice,xice_threshold,            &
!-- constants                   
              cp,rovcp,g0,lv,stbolt,                        &
!-- in/out                   
              tso,soilt,soilt1,tsnav,snow,snowh,snowc,      &
              qsfc,qsg,qvg,qcg,                             &
!-- out                   
              hfx,qfx,lh,dew,sfcrunoff,acrunoff,            &
              sfcevp,grdflx,snowfallac,acsnow,snom,         &
              globalcells,                                  &
              ids,ide, jds,jde, kds,kde,                    &
              ims,ime, jms,jme, kms,kme,                    &
              its,ite, jts,jte, kts,kte                     )
!-----------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------
!
! the RUC ice model is a part of RUC LSM described in:
!  Smirnova, T.G., J.M. Brown, and S.G. Benjamin, 1997:
!     performance of different soil model configurations in simulating
!     ground surface temperature and surface fluxes.
!     mon. wea. rev. 125, 1870-1884, 
!     https://doi.org/10.1175/1520-0493(1997)125%3C1870:PODSMC%3E2.0.CO;2 
!  Smirnova, T.G., J.M. Brown, and D. Kim, 2000: parameterization of
!     cold-season processes in the maps land-surface scheme.
!     j. geophys. res. 105, 4077-4086, https://doi.org/10.1029/1999JD901047
!  Smirnova, T. G., J. M. Brown, S. G. Benjamin, and J. S. Kenyon, 2016: 
!     Modifications to the Rapid Update Cycle land surface model (RUC LSM) 
!     available in the weather Research and forecasting model. 
!     Mon. Wea. Rev., 144, 1851–1865, https://doi.org/10.1175/MWR-D-15-0198.1  
!-----------------------------------------------------------------
!-- dt         time step (second)
!   ktau       number of time step
!   nsl        number of soil layers
!   nzs        number of levels in soil
!   zs         depth of soil levels (m)
!-- rainbl     accumulated rain in [mm] between the pbl calls
!-- rainncv    one time step grid scale precipitation (mm/step)
!-- snow       snow water equivalent [mm]
!-- frazfrac   fraction of frozen precipitation
!-- precipfr (mm) time step frozen precipitation
!-- snowc      flag indicating snow coverage (1 for snow cover)
!-- z3d        heights (m)
!-- p8w        3d pressure (pa)
!-- t3d        temperature (k)
!-- qv3d       3d water vapor mixing ratio (kg/kg)
!-- qc3d       3d cloud water mixing ratio (kg/kg)
!-- rho3d      3d air density (kg/m^3)
!-- glw        downward long wave flux at ground surface (w/m^2)
!-- gsw        absorbed short wave flux at ground surface (w/m^2)
!-- emiss      surface emissivity (between 0 and 1)
!   flqc       surface exchange coefficient for moisture (kg/m^2/s)
!   flhc       surface exchange coefficient for heat [w/m^2/s/degreek]
!   alb        surface albedo (between 0 and 1)
!   snoalb     maximum snow albedo (between 0 and 1)
!   albbck     snow-free albedo (between 0 and 1)
!   znt        roughness length [m]
!-- tbot       soil temperature at lower boundary (k)
!   ivgtyp     vegetation type
!   isltyp     soil type (16 classes)
!-- xland      land mask (1 for land, 2 for water)
!-- cp         heat capacity at constant pressure for dry air (j/kg/k)
!-- g0         acceleration due to gravity (m/s^2)
!-- lv         latent heat of melting (j/kg)
!-- stbolt     stefan-boltzmann constant (w/m^2/k^4)
!   tso        3d ice temp (k)
!-- soilt      ice surface temperature (k)
!-- hfx        upward heat flux at the surface (w/m^2)
!-- qfx        upward moisture flux at the surface (kg/m^2/s)
!-- lh         upward latent heat flux (w/m^2)
!   sfcrunoff  surface runoff [mm]
!   acrunoff   run-total surface runoff [mm]
!   sfcevp     total evaporation in [kg/m^2]
!   grdflx     soil heat flux (w/m^2: negative, if downward from surface)
!   snowfallac run-total snowfall accumulation [m]
!   acsnow     run-toral swe of snowfall [mm]
!-- chklowq    is either 0 or 1 (so far set equal to 1).
!--            used only in myjpbl.
!-- tice       sea ice temperture (c)
!-- rhosice    sea ice density (kg m^-3)
!-- capice     sea ice volumetric heat capacity (j/m^3/k)
!-- thdifice   sea ice thermal diffusivity (m^2/s)
!--
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-------------------------------------------------------------------------

   integer,     parameter       ::     nvegclas=24+3

   real,       intent(in   )    ::     dt
   logical,    intent(in   )    ::     myjpbl,frpcpn
   integer,    intent(in   )    ::     spp_lsm
   integer,    intent(in   )    ::     ktau, nsl, isice, iswater, &
                                       ims,ime, jms,jme, kms,kme, &
                                       ids,ide, jds,jde, kds,kde, &
                                       its,ite, jts,jte, kts,kte
   integer,   intent(in   )     ::     globalcells(ims:ime)

#if (em_core==1)
   real,    dimension( ims:ime, kms:kme, jms:jme ),optional::    pattern_spp_lsm
   real,    dimension( ims:ime, kms:kme, jms:jme ),optional::    field_sf
#endif
   real,    dimension( ims:ime, 1  :nsl, jms:jme )         ::    field_sf_loc

   real,    dimension( ims:ime, kms:kme, jms:jme )             , &
            intent(in   )    ::                            qv3d, &
                                                           qc3d, &
                                                            p8w, &
                                                          rho3d, &
                                                            t3d, &
                                                            z3d

   real,       dimension( ims:ime , jms:jme ),                   &
               intent(in   )    ::                       rainbl, &
                                                            glw, &
                                                            gsw, &
                                                         albbck, &
                                                           chs , &
                                                           xice, &
                                                          xland, &
                                                           tbot

   real,       dimension( ims:ime , jms:jme ),                   &
               intent(in   )    ::                    flqc,flhc


#if (em_core==1)
   real,       optional, dimension( ims:ime , jms:jme ),         &
               intent(in   )    ::                   graupelncv, &
                                                        snowncv, &
                                                        rainncv
#endif


   real,       dimension( 1:nsl), intent(in   )      ::      zs

   real,       dimension( ims:ime , jms:jme ),                   &
               intent(inout)    ::                      chklowq, &
                                                         snoalb, &
                                                           snow, &
                                                          snowh, &
                                                          snowc, &
                                                            alb, &
                                                            znt, &
                                                             z0, &
                                                          emiss

   real,       dimension( ims:ime , jms:jme ),                   &
               intent(in   )    ::                      frzfrac

   integer,    dimension( ims:ime , jms:jme ),                   &
               intent(in   )    ::                       ivgtyp, &
                                                         isltyp
   character(len=*), intent(in   )    ::                 mminlu

   real, intent(in   )          ::         cp,rovcp,g0,lv,stbolt,xice_threshold

   real,       dimension( ims:ime , 1:nsl, jms:jme )           , &
               intent(inout)    ::                       tso

   real,       dimension( ims:ime, jms:jme )                   , &
               intent(inout)    ::                        soilt, &
                                                            qvg, &
                                                            qcg, &
                                                            dew, &
                                                           qsfc, &
                                                            qsg, &
                                                         soilt1, &
                                                          tsnav
   real,       dimension( ims:ime, jms:jme )                   , &
               intent(inout)    ::                          hfx, &
                                                            qfx, &
                                                             lh, &
                                                         sfcevp, &
                                                      sfcrunoff, &
                                                       acrunoff, &
                                                         grdflx, &
                                                         acsnow, &
                                                           snom

   real,       dimension( ims:ime, jms:jme ), intent(out)     :: &
                                                         rhosnf, & ! rho of snowfall
                                                       precipfr, & ! time-step frozen precip
                                                     snowfallac

!-- local variables
   real,       dimension( its:ite, jts:jte )    ::               &
                                                        runoff1, &
                                                        runoff2, &
                                                         emissl, &
                                                         budget, &
                                                           zntl, &
                                                          smelt, &
                                                           snoh, &
                                                          snflx, &
                                                           edir, &
                                                             ec, &
                                                            ett, &
                                                         sublim, &
                                                           sflx, &
                                                            smf, &
                                                          evapl, &
                                                          prcpl, &
                                                         seaice, &
                                                        infiltr

!-- soil/snow properties
   real                                                          &
                             ::                           rhocs, &
                                                       rhonewsn, &
                                                          rhosn, &
                                                      rhosnfall, &
                                                           bclh, &
                                                            dqm, &
                                                           ksat, &
                                                           psis, &
                                                           qmin, &
                                                          qwrtz, &
                                                            ref, &
                                                           wilt, &
                                                       snowfrac, &
                                                          snhei, &
                                                           snwe

   real                                      ::              cn, &
                                                         sat,cw, &
                                                           c1sn, &
                                                           c2sn


   real,     dimension(1:nsl)                ::          zsmain, &
                                                         zshalf, &
                                                         dtdzs2

   real,     dimension(1:nsl)                ::          soilice,&
                                                         soiliqw           

   real,     dimension(1:2*(nsl-2))          ::           dtdzs

   real,     dimension(1:5001)               ::             tbq


   real,     dimension( 1:nsl )              ::           tso1d

   real                           ::                        rsm, &
                                                      snweprint, &
                                                     snheiprint

   real                           ::                     prcpms, &
                                                        newsnms, &
                                                      prcpncliq, &
                                                       prcpncfr, &
                                                      prcpculiq, &
                                                       prcpcufr, &
                                                           patm, &
                                                          patmb, &
                                                           tabs, &
                                                          qvatm, &
                                                          qcatm, &
                                                          q2sat, &
                                                         conflx, &
                                                            rho, &
                                                           qkms, &
                                                           tkms, &
                                                        snowrat, &
                                                       grauprat, &
                                                       graupamt, &
                                                         icerat, &
                                                          curat, &
                                                       infiltrp
   real      ::  cq,r61,r273,arp,brp,x,evs,eis

   integer   ::  iland,isoil,iforest

   integer   ::  i,j,k,nzs,nzs1,nddzs
   integer   ::  k1,l,k2,kp,km
   integer   ::  globalcellid
   character (len=132) :: message
   logical   :: print_flag = .false.
   real,dimension(ims:ime,1:nsl,jms:jme) :: rstoch
   real,dimension(1:nsl)::rstoch_temp,field_sf_temp
   real,dimension(ims:ime,jms:jme)::emisso,vegfrao,albo,snoalbo
   real,dimension(its:ite,jts:jte)::emisslo

!-----------------------------------------------------------------
   nzs=nsl
   nddzs=2*(nzs-2)

   rstoch=0.0
   field_sf_loc=0.0

#if (em_core==1)
   if (spp_lsm==1) then
      do j=jts,jte
         do i=its,ite
            do k=1,nsl
               rstoch(i,k,j) = pattern_spp_lsm(i,k,j)
               field_sf_loc(i,k,j)=field_sf(i,k,j)
            enddo
         enddo
      enddo
   endif
#endif
   do j=jts,jte
      do i=its,ite
      !-- check if ice fraction is higher than threshold
         if(xice(i,j).ge.xice_threshold) then
            seaice(i,j)=1.
         else
            seaice(i,j)=0.
         endif
      enddo
   enddo   


!---- table tbq is for resolution of balance equation in vilka
   cq=173.15-.05
   r273=1./273.15
   r61=6.1153*0.62198
   arp=77455.*41.9/461.525
   brp=64.*41.9/461.525

   do k=1,5001
      cq=cq+.05
      evs=exp(17.67*(cq-273.15)/(cq-29.65))
      eis=exp(22.514-6.15e3/cq)
      if(cq.ge.273.15) then
      ! tbq is in mb
         tbq(k) = r61*evs
      else
         tbq(k) = r61*eis
      endif
   end do

#if ( nmm_core == 1 )
   if(ktau+1.eq.1) then
#else
   if(ktau.eq.1) then
#endif
!--- initialize ice/snow variables at the first time step
      do j=jts,jte
         do i=its,ite

          if(seaice(i,j).gt.0.5)then
!--- Sea ice point
!--- initializing snow fraction, thereshold = 32 mm of snow water or ~100 mm of snow height
! call physics_message('!--- initializing inside snow temp if it is not defined')
            if((soilt1(i,j) .lt. 170.) .or. (soilt1(i,j) .gt.400.)) then
               if(snowc(i,j).gt.0.) then
                  soilt1(i,j)=0.5*(soilt(i,j)+tso(i,1,j))
                  if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
                     write ( message , fmt='(a,f8.3,2i6)' ) &
                     'temperature inside snow is initialized in ruc_ice ', soilt1(i,j),i,j
                     !call wrf_debug ( 0 , message )
                  endif
               else
                  soilt1(i,j) = tso(i,1,j)
               endif ! snowc
            endif ! soilt1
            !-- temperature inside snow is initialized
            tsnav(i,j) =0.5*(soilt(i,j)+tso(i,1,j))-273.15
            patmb=p8w(i,kms,j)*1.e-2
            qsg  (i,j) = qsn(soilt(i,j),tbq)/patmb
            if((qcg(i,j) < 0.) .or. (qcg(i,j) > 0.1)) then
               qcg  (i,j) = 0.
               if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
                   write ( message , fmt='(a,3f8.3,2i6)' ) &
                   'qcg is initialized in ruc_ice ', qvg(i,j),qsg(i,j),i,j
               endif
            endif ! qcg
            if((qvg(i,j) .le. 0.) .or. (qvg(i,j) .gt.0.1)) then
               qvg  (i,j) = qsg(i,j)
               if ( wrf_at_debug_level(lsmruc_dbg_lvl) ) then
                  write ( message , fmt='(a,3f8.3,2i6)' ) &
                  'qvg is initialized in rucice ', qvg(i,j)
               endif
            endif
            qsfc(i,j) = qvg(i,j)/(1.+qvg(i,j))
            smelt(i,j) = 0.
            snom (i,j) = 0.
            snowfallac(i,j) = 0.
            acsnow(i,j) = 0.
            precipfr(i,j) = 0.
            rhosnf(i,j) = -1.e3 ! non-zero flag
            snflx(i,j) = 0.
            dew  (i,j) = 0.
            zntl (i,j) = 0.
            runoff1(i,j) = 0.
            runoff2(i,j) = 0.
            sfcrunoff(i,j) = 0.
            acrunoff(i,j) = 0.
            emissl (i,j) = 0.

!---  for ruc lsm chklowq needed for myjpbl should
!     be 1 because is actual specific humidity at the surface, and
!     not the saturation value
            chklowq(i,j) = 1.
            infiltr(i,j) = 0.
            snoh  (i,j) = 0.
            edir  (i,j) = 0.
            ec    (i,j) = 0.
            ett   (i,j) = 0.
            sublim(i,j) = 0.
            sflx  (i,j) = 0.
            smf   (i,j) = 0.
            evapl (i,j) = 0.
            prcpl (i,j) = 0.
          endif ! sea ice point
         enddo
      enddo

      do k=1,nsl
         soilice(k)=1.
         soiliqw(k)=0.
      enddo
   endif
   !---  end of initialization
!-----------------------------------------------------------------

   prcpms = 0.
   newsnms = 0.
   prcpncliq = 0.
   prcpculiq = 0.
   prcpncfr = 0.
   prcpcufr = 0.


   !call mpas_log_write('--- in ruc_ice before main loop:')
   do j=jts,jte
      do i=its,ite
        if(seaice(i,j).gt.0.5)then
!--- sea ice point
         !if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
         if (globalcells(i)==targetcell) then
            print *,' sea ice point: ktau,globalcells(i),xice(i,j)', ktau,globalcells(i),xice(i,j)
            print *,' soilt,p8w',soilt(i,j),p8w(i,1,j)
            print *,' xland, qfx,hfx,chs,flhc from sfclay',xland(i,j),qfx(i,j),hfx(i,j),chs(i,j),flhc(i,j)
            print *,' gsw, glw =',gsw(i,j),glw(i,j)
            print *,' soilt, tso start of time step =',soilt(i,j),(tso(i,k,j),k=1,nsl)
            print *,' ivgtyp,isltyp,alb = ', ivgtyp(i,j),isltyp(i,j),alb(i,j)
            print *,' dt,rainbl =',dt,rainbl(i,j)
            print *,' snow, snowh, snowc =',snow(i,j), snowh(i,j), snowc(i,j)
         endif

         iland     = isice
         isoil     = isltyp(i,j)
         tabs      = t3d(i,kms,j)
         qvatm     = qv3d(i,kms,j)
         qcatm     = qc3d(i,kms,j)
         patm      = p8w(i,kms,j)*1.e-5
!-- z3d(1) is thickness between first full sigma level and the surface,
!-- but first mass level is at the half of the first sigma level
!-- (u and v are also at the half of first sigma level)
         conflx    = z3d(i,kms,j)*0.5
         rho       = rho3d(i,kms,j)
! -- initialize snow, graupel and ice fractions in frozen precip
         snowrat = 0.
         grauprat = 0.
         icerat = 0.
         curat = 0.

         if(frpcpn) then
#if (em_core==1)
            prcpncliq = rainncv(i,j)*(1.-frzfrac(i,j))
            prcpncfr = rainncv(i,j)*frzfrac(i,j)
!- apply the same frozen precipitation fraction to convective precip
!- 31 mar17 - add safety temperature check in case Thompson MP produces
!                 frozen precip at t > 273.
            if(frzfrac(i,j) > 0..and. tabs < 273.) then
               prcpculiq = max(0.,(rainbl(i,j)-rainncv(i,j))*(1.-frzfrac(i,j)))
               prcpcufr = max(0.,(rainbl(i,j)-rainncv(i,j))*frzfrac(i,j))
            else
               if(tabs < 273.) then
                  prcpcufr = max(0.,(rainbl(i,j)-rainncv(i,j)))
                  prcpculiq = 0.
               else
                  prcpcufr = 0.
                  prcpculiq = max(0.,(rainbl(i,j)-rainncv(i,j)))
               endif  ! tabs < 273.
            endif  ! frzfrac > 0.
!--- 1*e-3 is to convert from mm/s to m/s
            prcpms   = (prcpncliq + prcpculiq)/dt*1.e-3
            newsnms  = (prcpncfr + prcpcufr)/dt*1.e-3
 
            if ( present( graupelncv ) ) then
               graupamt = graupelncv(i,j)
            else
               graupamt = 0.
            endif

            if((prcpncfr + prcpcufr) > 0.) then
! -- calculate snow, graupel and ice fractions in falling frozen precip
               snowrat=min(1.,max(0.,snowncv(i,j)/(prcpncfr + prcpcufr)))
               grauprat=min(1.,max(0.,graupamt/(prcpncfr + prcpcufr)))
               icerat=min(1.,max(0.,(prcpncfr-snowncv(i,j)-graupamt) &
                                   /(prcpncfr + prcpcufr)))
               curat=min(1.,max(0.,(prcpcufr/(prcpncfr + prcpcufr))))
            endif
#else
            prcpms    = (rainbl(i,j)/dt*1.e-3)*(1-frzfrac(i,j))
            newsnms  = (rainbl(i,j)/dt*1.e-3)*frzfrac(i,j)
            if(newsnms == 0.) then
               snowrat = 0.
            else
               snowrat = min(1.,newsnms/(newsnms+prcpms))
            endif
#endif

         else  ! .not. frpcpn
            if (tabs.le.273.15) then
               prcpms    = 0.
               newsnms   = rainbl(i,j)/dt*1.e-3
!-- here no info about constituents of frozen precipitation,
!-- suppose it is all snow
               snowrat = 1.
            else
               prcpms    = rainbl(i,j)/dt*1.e-3
               newsnms   = 0.
            endif
         endif

! -- save time-step water equivalent of frozen precipitation in precipfr array to be used in
!    module_diagnostics
         precipfr(i,j) = newsnms * dt *1.e3

!--- convert exchange coeff qkms to [m/s]
         qkms=flqc(i,j)/rho
         tkms=flhc(i,j)/rho/(cp*(1.+0.84*qvatm))  ! mynnsfc uses cpm
!--- convert incoming snow and canwat from mm to m
         snwe=snow(i,j)*1.e-3
         snhei=snowh(i,j)

         snowfrac=snowc(i,j)
         rhosnfall=rhosnf(i,j)
!-----
         zsmain(1)=0.
         zshalf(1)=0.
         do k=2,nzs
            zsmain(k)= zs(k)
            zshalf(k)=0.5*(zsmain(k-1) + zsmain(k))
         enddo

!------------------------------------------------------------
!-----  ddzs and dsdz1 are for implicit solution of soil eqns.
!-------------------------------------------------------------
         nzs1=nzs-1
!-----
         if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
         !if (globalcells(i) == targetcell) then
            print *,' dt,nzs1, zsmain, zshalf --->', dt,nzs1,zsmain,zshalf
         endif

         do k=2,nzs1
            k1=2*k-3
            k2=k1+1
            x=dt/2./(zshalf(k+1)-zshalf(k))
            dtdzs(k1)=x/(zsmain(k)-zsmain(k-1))
            dtdzs2(k-1)=x
            dtdzs(k2)=x/(zsmain(k+1)-zsmain(k))
!           if (globalcells(i) == targetcell) then
!              print *,' k,k1,k2, dtdzs(k1), dtdzs(k2) ',k,k1,k2, dtdzs(k1), dtdzs(k2)
!           endif
         end do

         cw =4.183e6

!***********************************************************************
!--- constants for snow density calculations c1sn and c2sn

         c1sn=0.026
         c2sn=21.

!***********************************************************************

         rhonewsn = 200.
         if(snow(i,j).gt.0. .and. snowh(i,j).gt.0.) then
            rhosn = snow(i,j)/snowh(i,j)
         else
            rhosn = 300.
         endif

            iland = isice
            isoil = 16
            znt(i,j) = 0.011
            snoalb(i,j) = 0.75
            dqm = 1.
            ref = 1.
            qmin = 0.
            wilt = 0.
            emissl(i,j) = 0.98

            patmb=p8w(i,1,j)*1.e-2
            qvg  (i,j) = qsn(soilt(i,j),tbq)/patmb
            qsg  (i,j) = qvg(i,j)
            qsfc(i,j) = qvg(i,j)/(1.+qvg(i,j))

            do k=1,nzs
               tso(i,k,j) = min(271.4,tso(i,k,j))
               tso1d  (k) = tso(i,k,j)
            enddo

            !if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
            if (globalcells(i)==targetcell) then
               print *,'globalcells(i)= ',globalcells(i),i
               print *,'ice, tso1d,patm,tabs,qvatm,qcatm,rho,qvg',  &
                             tso1d,patm,tabs,qvatm,qcatm,rho,qvg(i,j)
               print *,'conflx =',conflx
            endif

            rstoch_temp = rstoch(i,:,j)
            field_sf_temp = field_sf_loc(i,:,j)
!-----------------------------------------------------------------
!-- solve energy budget

            call icetmp (spp_lsm,rstoch_temp,field_sf_temp,       &
!--- input variables
                 dt,ktau,conflx,i,j,nzs,nddzs,                    & 
                 iland,isoil,xland(i,j),ivgtyp(i,j),isltyp(i,j),  &
                 prcpms,newsnms,snowrat,grauprat,icerat,curat,    &
                 patm,tabs,qvatm,qcatm,rho,                       &
                 glw(i,j),gsw(i,j),emissl(i,j),                   &
                 qkms,tkms,alb(i,j),znt(i,j),                     &
                 snoalb(i,j),albbck(i,j),myjpbl,                  &
                 seaice(i,j),isice,                               &
                 zsmain,zshalf,dtdzs,dtdzs2,tbq,                  &
!--- constants
                 cp,rovcp,g0,lv,stbolt,cw,c1sn,c2sn,              &
!--- in/out
                 tso1d,soilt(i,j),soilt1(i,j),tsnav(i,j),         &
                 rhosn,rhonewsn,rhosnfall,                        &
                 snwe,snhei,snowfrac,                             &
!--- output variables
                 rsm,snweprint,snheiprint,                        &
                 dew(i,j),qvg(i,j),qsg(i,j),qcg(i,j),smelt(i,j),  &
                 snoh(i,j),snflx(i,j),snom(i,j),snowfallac(i,j),  &
                 acsnow(i,j),edir(i,j),ec(i,j),ett(i,j),qfx(i,j), &
                 lh(i,j),hfx(i,j),sflx(i,j),sublim(i,j),          &
                 evapl(i,j),prcpl(i,j),budget(i,j),runoff1(i,j),  &
                 runoff2(i,j),soilice,soiliqw,infiltrp,smf(i,j),  &
                 globalcells(i) )

            field_sf_loc(i,:,j) = field_sf_temp
!-----------------------------------------------------------------
! fill in field_sf to pass perturbed field of hydraulic cond. up to model driver and output
#if (em_core==1)
            if (spp_lsm==1) then
               do k=1,nsl
                  field_sf(i,k,j)=field_sf_loc(i,k,j)
               enddo
            endif
#endif


            do k=1,nzs
               tso(i,k,j) = tso1d(k)
            enddo

            z0       (i,j) = znt (i,j)
            patmb=p8w(i,1,j)*1.e-2
            q2sat=qsn(tabs,tbq)/patmb
            qsfc(i,j) = qvg(i,j)/(1.+qvg(i,j))

            if(snow(i,j)==0.) emissl(i,j) = 0.98
            emiss (i,j) = emissl(i,j)
! snow is in [mm], snwe is in [m]; canwat is in mm, canwatr is in m
            snow   (i,j) = snwe*1000.
            snowh  (i,j) = snhei

            if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
               print *,' land, i=,j=, qfx, hfx after icetmp', i,j,lh(i,j),hfx(i,j)
            endif
            grdflx (i,j) = -1. * sflx(i,j)

!--- snowc snow cover flag
            if(snowfrac > 0. .and. xice(i,j).ge.xice_threshold ) then
               snowfrac = snowfrac*xice(i,j)
            endif
 
            snowc(i,j)=snowfrac

!--- rhosnf - density of snowfall
            rhosnf(i,j)=rhosnfall

! accumulated moisture flux [kg/m^2]
            sfcevp (i,j) = sfcevp (i,j) + qfx (i,j) * dt

            if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
               print *,'land, i,j,tso1d,soilt - end of time step',         &
                              i,j,tso1d,soilt(i,j)
               print *,'land, qfx, hfx after icetmp', i,j,lh(i,j),hfx(i,j)
            endif

        endif ! sea ice point is done

      enddo
   enddo

!-----------------------------------------------------------------
   end subroutine ruc_ice
!-----------------------------------------------------------------



   subroutine icetmp (spp_lsm,rstochcol,fieldcol_sf,             &
!--- input variables
                delt,ktau,conflx,i,j,nzs,nddzs,                  &
                iland,isoil,xland,ivgtyp,isltyp,                 &
                prcpms,newsnms,snowrat,grauprat,icerat,curat,    &
                patm,tabs,qvatm,qcatm,rho,                       &
                glw,gsw,emiss,qkms,tkms,                         &
                alb,znt,alb_snow,alb_snow_free,                  &
                myj,seaice,isice,                                &
                zsmain,zshalf,dtdzs,dtdzs2,tbq,                  &
!--- constants
                cp,rovcp,g0,lv,stbolt,cw,c1sn,c2sn,              &
!--- in/out
                ts1d,soilt,soilt1,tsnav,                         &
                rhosn,rhonewsn,rhosnfall,                        &
                snwe,snhei,snowfrac,                             &
!--- output variables
                rsm,snweprint,snheiprint,                        &
                dew,qvg,qsg,qcg,                                 &
                smelt,snoh,snflx,snom,snowfallac,acsnow,         &
                edir1,ec1,ett1,eeta,qfx,hfx,s,sublim,            &
                evapl,prcpl,fltot,runoff1,runoff2,soilice,       &
                soiliqw,infiltr,smf,globalcellid)

!-----------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------

!--- input variables

   integer,  intent(in   )   ::  isice,i,j,ktau,nzs ,            &
                                 nddzs,globalcellid              !nddzs=2*(nzs-2)

   real,     intent(in   )   ::  delt,conflx
   real,     intent(in   )   ::  c1sn,c2sn
   logical,  intent(in   )   ::  myj
!--- 3-d atmospheric variables
   real                                                        , &
            intent(in   )    ::                            patm, &
                                                           tabs, &
                                                          qvatm, &
                                                          qcatm
   real                                                        , &
            intent(in   )    ::                             glw, &
                                                            gsw, &
                                                  alb_snow_free, &
                                                       alb_snow, &
                                                         seaice, &
                                                          xland, &
                                                            rho, &
                                                           qkms, &
                                                           tkms

   integer,   intent(in   )  ::                          ivgtyp, isltyp
!--- 2-d variables
   real                                                        , &
            intent(inout)    ::                           emiss, &
                                                       snowfrac, &
                                                            alb

   real,     intent(in   )   ::                                  &
                                                             cw, &
                                                             cp, &
                                                          rovcp, &
                                                             g0, &
                                                             lv, &
                                                         stbolt

   real,     dimension(1:nzs), intent(in)  ::            zsmain, &
                                                         zshalf, &
                                                         dtdzs2

   real,     dimension(1:nzs), intent(in)  ::          rstochcol
   real,     dimension(1:nzs), intent(inout) ::     fieldcol_sf


   real,     dimension(1:nddzs), intent(in)  ::           dtdzs

   real,     dimension(1:5001), intent(in)  ::              tbq


!--- input/output variables
!-------- 3-d soil moisture and temperature
   real,     dimension( 1:nzs )                                , &
             intent(inout)   ::                            ts1d

   real,  dimension(1:nzs), intent(inout)  ::           soilice, &
                                                        soiliqw


   integer, intent(inout)    ::                     iland,isoil
   integer                   ::                     ilands

!-------- 2-d variables
   real,      intent(inout)   ::                                 &
                                                            qvg, &
                                                            qsg, &
                                                            qcg, &
                                                          rhosn, &
                                                     snowfallac, &
                                                           snwe, &
                                                          snhei, &
                                                          soilt, &
                                                         soilt1, &
                                                          tsnav, &
                                                           snom, &
                                                            znt


   real,      intent(inout)   ::                            dew, &
                                                          edir1, &
                                                            ec1, &
                                                           ett1, &
                                                           eeta, &
                                                          evapl, &
                                                        infiltr, &
                                                       rhonewsn, &
                                                      rhosnfall, &
                                                        snowrat, &
                                                       grauprat, &
                                                         icerat, &
                                                          curat, &
                                                         sublim, &
                                                          prcpl, &
                                                            qfx, &
                                                            hfx, &
                                                          fltot, &
                                                            smf, &
                                                              s, &
                                                        runoff1, &
                                                        runoff2, &
                                                         acsnow, &
                                                          smelt, &
                                                           snoh, &
                                                          snflx

   real,     dimension(1:nzs)              ::                    &
                                                           tice, &
                                                        rhosice, &
                                                         capice, &
                                                       thdifice, &
                                                          ts1ds
!-------- 1-d variables
   real :: &
                                                            dews, &
                                                          edir1s, &
                                                            ec1s, &
                                                            csts, &
                                                           ett1s, &
                                                           eetas, &
                                                          evapls, &
                                                        infiltrs, &
                                                          prcpls, &
                                                            qvgs, &
                                                            qsgs, &
                                                            qcgs, &
                                                            qfxs, &
                                                            hfxs, &
                                                          fltots, &
                                                        runoff1s, &
                                                        runoff2s, &
                                                              ss, &
                                                          soilts




   real,  intent(inout)                     ::             rsm,  &
                                                      snweprint, &
                                                     snheiprint

   integer,   intent(in)                    ::     spp_lsm

!--- local variables

   integer ::  k,ilnb

   real    ::  bsn, xsn                                        , &
               rainf, snth, newsn, prcpms, newsnms             , &
               t3, upflux, xinet
   real    ::  snhei_crit, snhei_crit_newsn, keep_snow_albedo, snowfracnewsn
   real    ::  newsnowratio, dd1, snowfrac2, m

   real    ::  rhonewgr,rhonewice

   real    ::  rnet,gswnew,gswin,emissn,zntsn,emiss_snowfree
   real    ::  snow_mosaic, snfr
   real    ::  cice, albice, albsn, drip

!-----------------------------------------------------------------
   integer,   parameter      ::      ilsnow=99

   !if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
   if (globalcellid==targetcell) then
      print *,' in icetmp nzs,nddzs,snwe,rhosn,snom,ts1d ',nzs,nddzs,snwe,rhosn,snom,ts1d
   endif

   !-- snow fraction options
   !-- option 1: original formulation using critical snow depth to compute
   !-- snow fraction
   !-- option 2: the tanh formulation from Niu,G.-Y.,and Yang,Z.-L.
   !2007,JGR,doi:10.1029/2007jd008674.
   !-- option 3: the tanh formulation from Niu,G.-Y.,and Yang,Z.-L.
   !2007,JGR,doi:10.1029/2007jd008674.
   !   with vegetation dependent parameters from noah mp (personal
   !   communication with mike barlage)
   !-- snhei_crit is a threshold for fractional snow in isncovr_opt=1
   snhei_crit=0.01601*rhowater/rhosn
   snhei_crit_newsn=0.0005*rhowater/rhosn
   !--
   zntsn = 0.011

   snow_mosaic=0.
   snfr = 1.
   newsn=0.
   newsnowratio = 0.
   snowfracnewsn=0.
   rhonewsn = 100.
   if(snhei == 0.) snowfrac=0.
   smelt = 0.
   rainf = 0.
   rsm=0.
   infiltr=0.
   drip = 0.
   smf = 0.

!---initialize local arrays for sea ice
   do k=1,nzs
      tice(k) = 0.
      rhosice(k) = 0.
      cice = 0.
      capice(k) = 0.
      thdifice(k) = 0.
   enddo

   gswnew=gsw
   gswin=gsw/(1.-alb)
   albice=alb_snow_free
   albsn=alb_snow
   emissn = 0.98
   emiss_snowfree = 0.98

!--- sea ice properties
!--- N.N Zubov "Arctic Ice"
!--- no salinity dependence because we consider the ice pack
!--- to be old and to have low salinity (0.0002)
   do k=1,nzs
      tice(k) = ts1d(k) - 273.15
      rhosice(k) = 917.6/(1-0.000165*tice(k))
      cice = 2115.85 +7.7948*tice(k)
      capice(k) = cice*rhosice(k)
      thdifice(k) = 2.260872/capice(k)
   enddo
!-- sea ice alb dependence on ice temperature. when ice temperature is
!-- below critical value of -10c - no change to albedo.
!-- if temperature is higher that -10c then albedo is decreasing.
!-- the minimum albedo at t=0c for ice is 0.1 less.
   albice = min(alb_snow_free,max(alb_snow_free - 0.05,   &
                alb_snow_free - 0.1*(tice(1)+10.)/10. ))

   !if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
   if (globalcellid==targetcell) then
      print *,'alb_snow_free',alb_snow_free
      print *,'gsw,gswnew,glw,soilt,emiss,alb,albice,snwe',&
               gsw,gswnew,glw,soilt,emiss,alb,albice,snwe
   endif

   if(snhei.gt.0.0081*1.e3/rhosn) then
!*** update snow density for current temperature (Koren et al. 1999)
      bsn=delt/3600.*c1sn*exp(0.08*min(0.,tsnav)-c2sn*rhosn*1.e-3)
      if(bsn*snwe*100..gt.1.e-4) then
         xsn=rhosn*(exp(bsn*snwe*100.)-1.)/(bsn*snwe*100.)
         rhosn=min(max(58.8,xsn),500.) 
      endif
   endif

   !-- snow_mosaic from the previous time step 
   if(snowfrac < 0.75) snow_mosaic = 1.

   newsn=newsnms*delt

   if(newsn.gt.0.) then

      if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
         print *, 'there is new snow, newsn', newsn
      endif

      newsnowratio = min(1.,newsn/(snwe+newsn))

!--- density of hydrometeors as a function of temperature
      rhonewsn=min(125.,1000.0/max(8.,(17.*tanh((276.65-tabs)*0.15))))
      rhonewgr=min(500.,rhowater/max(2.,(3.5*tanh((274.15-tabs)*0.3333))))
      rhonewice=rhonewsn

!--- compute density of "snowfall" from weighted contribution
!                 of snow, graupel and ice fractions

      rhosnfall = min(500.,max(58.8,(rhonewsn*snowrat +  &  ! 13mar18-switch from 76.9 to 58.8
                  rhonewgr*grauprat + rhonewice*icerat + rhonewgr*curat)))

! from now on rhonewsn is the density of falling frozen precipitation
      rhonewsn=rhosnfall

!*** define average snow density of the snow pack considering
!*** the amount of fresh snow (eq. 9 in koren et al.(1999)
!*** without snow melt )
      xsn=(rhosn*snwe+rhonewsn*newsn)/                         &
          (snwe+newsn)
      rhosn=min(max(58.8,xsn),500.) ! 13mar18 - switch from 76.9 to 58.8

   endif ! end newsn > 0.

   if(prcpms.gt.0.) rainf = 1.

   drip = 0.

   if(newsn.gt.0.) then
      snwe=max(0.,snwe+newsn)
   endif
   snhei=snwe*rhowater/rhosn
   newsn=newsn*rhowater/rhonewsn

   if(snhei.gt.0.0) then
!--- snow on the ice
      iland=isice
      snowfrac=min(1.,snhei/(2.*snhei_crit))

      if(newsn > 0. ) then
         snowfracnewsn=min(1.,snowfallac*1.e-3/snhei_crit_newsn)
      endif

      if(snowfrac < 0.75) snow_mosaic = 1.

      keep_snow_albedo = 0.
      if (snowfracnewsn > 0.99 .and. rhosnfall < 450.) then
!--- new snow
         keep_snow_albedo = 1.
         snow_mosaic = 0.
      endif

      !if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
      if (globalcellid==targetcell) then
         print *,'snhei_crit,snowfrac,snhei_crit_newsn,snowfracnewsn', &
                  snhei_crit,snowfrac,snhei_crit_newsn,snowfracnewsn
      endif

!----- snow on ice
      if( snow_mosaic == 1.) then
         albsn=alb_snow
         emiss= emissn
      else
         albsn   = max(keep_snow_albedo*alb_snow,               &
                   min((albice + (alb_snow - albice) * snowfrac), alb_snow))
         emiss   = max(keep_snow_albedo*emissn,                 &
                   min((emiss_snowfree +                        &
           (emissn - emiss_snowfree) * snowfrac), emissn))
      endif

      !if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
      if (globalcellid==targetcell) then
         print *,'snow on ice snow_mosaic,albsn,emiss',albsn,emiss,snow_mosaic
      endif
!-- alb dependence on snow temperature. when snow temperature is
!-- below critical value of -10c - no change to albedo.
!-- if temperature is higher that -10c then albedo is decreasing.
      if(albsn.lt.alb_snow .or. keep_snow_albedo .eq.1.)then
         alb=albsn
      else
!-- change albedo when no fresh snow
         alb = min(albsn,max(albsn - 0.15*albsn*(soilt - 263.15)/  &
                  (273.15-263.15), albsn - 0.1))
      endif

      if (snow_mosaic==1.) then
!-- sea ice
!-- portion not covered with snow
!-- compute absorbed gsw for snow-free portion

         gswnew=gswin*(1.-albice)
!--------------
         t3      = stbolt*soilt*soilt*soilt
         upflux  = t3 *soilt
         xinet   = emiss_snowfree*(glw-upflux)
         rnet    = gswnew + xinet
         !if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
         if (globalcellid==targetcell) then
            print *,'fractional snow - snowfrac=',snowfrac
            print *,'snowfrac<1 gswin,gswnew -',gswin,gswnew,'soilt, rnet',soilt,rnet
         endif
         do k=1,nzs
            ts1ds(k) = ts1d(k)
         enddo
         soilts = soilt
         qvgs = qvg
         qsgs = qsg
         qcgs = qcg
         smelt=0.
         runoff1s=0.
         runoff2s=0.

         call sice(                                               &
!--- input variables
              i,j,iland,isoil,delt,ktau,conflx,nzs,nddzs,         &
              prcpms,rainf,patm,qvatm,qcatm,glw,gswnew,           &
              0.98,rnet,qkms,tkms,rho,myj,globalcellid,           &
!--- sea ice parameters
              tice,rhosice,capice,thdifice,                       &
              zsmain,zshalf,dtdzs,dtdzs2,tbq,                     &
!--- constants
              lv,cp,rovcp,cw,stbolt,tabs,                         &
!--- output variable
              ts1ds,dews,soilts,qvgs,qsgs,qcgs,                   &
              eetas,qfxs,hfxs,ss,evapls,prcpls,fltots             )

         edir1 = eeta*1.e-3
         ec1 = 0.
         ett1 = 0.
         runoff1 = prcpms
         runoff2 = 0.
         infiltr=0.

         !if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
         if (globalcellid==targetcell) then
            print *,'incoming gswnew snowfrac<1 -',gswnew
         endif

      endif ! snow_mosaic=1.

!--- recompute absorbed solar radiation and net radiation
!--- for updated value of snow albedo - alb
      gswnew=gswin*(1.-alb)
!--------------
      t3      = stbolt*soilt*soilt*soilt
      upflux  = t3 *soilt
      xinet   = emiss*(glw-upflux)
      rnet    = gswnew + xinet
      !if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
      if (globalcellid==targetcell) then
         print *,'rnet=',rnet
         print *,'snow - i,j,newsn,snwe,snhei,gsw,gswnew,glw,upflux,alb',&
                         i,j,newsn,snwe,snhei,gsw,gswnew,glw,upflux,alb
      endif

      if(snow_mosaic==1.)then
         snfr=1.
      else
         snfr=snowfrac
      endif

      call snowseaice (                                        &
           i,j,isoil,delt,ktau,conflx,nzs,nddzs,               &
           rhonewsn,snhei_crit,                                &  
           iland,prcpms,rainf,newsn,snhei,snwe,snfr,           &
           rhosn,patm,qvatm,qcatm,                             &
           glw,gswnew,emiss,rnet,                              &
           qkms,tkms,rho,myj,globalcellid,                     &
!--- sea ice parameters
           alb,znt,                                            &
           tice,rhosice,capice,thdifice,                       &
           zsmain,zshalf,dtdzs,dtdzs2,tbq,                     &
!--- constants
           lv,cp,rovcp,cw,stbolt,tabs,                         &
!--- output variables
           ilnb,snweprint,snheiprint,rsm,ts1d,                 &
           dew,soilt,soilt1,tsnav,qvg,qsg,qcg,                 &
           smelt,snoh,snflx,snom,eeta,                         &
           qfx,hfx,s,sublim,prcpl,fltot                        )

      edir1 = eeta*1.e-3
      ec1 = 0.
      ett1 = 0.
      runoff1 = smelt
      runoff2 = 0.
      infiltr=0.

      if (snow_mosaic==1.) then
! may 2014 - now combine snow covered and snow-free land fluxes, ice temp
! now combine fluxes for snow-free sea ice and snow-covered area
         do k=1,nzs
            ts1d(k) = ts1ds(k)*(1.-snowfrac) + ts1d(k)*snowfrac
         enddo
         dew = dews*(1.-snowfrac) + dew*snowfrac
         soilt = soilts*(1.-snowfrac) + soilt*snowfrac
         qvg = qvgs*(1.-snowfrac) + qvg*snowfrac
         qsg = qsgs*(1.-snowfrac) + qsg*snowfrac
         qcg = qcgs*(1.-snowfrac) + qcg*snowfrac
         eeta = eetas*(1.-snowfrac) + eeta*snowfrac
         qfx = qfxs*(1.-snowfrac) + qfx*snowfrac
         hfx = hfxs*(1.-snowfrac) + hfx*snowfrac
         s = ss*(1.-snowfrac) + s*snowfrac
         sublim = eeta
         prcpl = prcpls*(1.-snowfrac) + prcpl*snowfrac
         fltot = fltots*(1.-snowfrac) + fltot*snowfrac
!alb
         alb   = max(keep_snow_albedo*alb,              &
                 min((albice + (alb - alb_snow_free) * snowfrac), alb))

         emiss = max(keep_snow_albedo*emissn,           &
                 min((emiss_snowfree +                  &
               (emissn - emiss_snowfree) * snowfrac), emissn))

         runoff1 = runoff1s*(1.-snowfrac) + runoff1*snowfrac
         runoff2 = runoff2s*(1.-snowfrac) + runoff2*snowfrac
         if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
            print *,'soilt combined on ice', soilt
         endif

      endif ! snow_mosaic = 1.
 
      if(snhei.eq.0.) then
     !-- all snow is melted
        alb=alb_snow_free
        iland=isice
      else
      !-- snow on ice
         snowfrac=min(1.,snhei/(2.*snhei_crit))
      endif

!  run-total accumulated snow based on snowfall and snowmelt in [m]

      snowfallac = snowfallac + newsn * 1.e3    ! accumulated snow depth [mm], using variable snow den
      !snowfallac = snowfallac + max(0.,(newsn - rhowater/rhonewsn*smelt*delt*newsnowratio))
      acsnow = snowfallac
   else
!--- no snow
      snheiprint=0.
      snweprint=0.
      smelt=0.
      t3      = stbolt*soilt*soilt*soilt
      upflux  = t3 *soilt
      xinet   = emiss*(glw-upflux)

!--- if current ice albedo is not the same as from the previous time step, then
!--- update gsw, alb and rnet for surface energy budget
      if(alb.ne.albice) gswnew=gsw/(1.-alb)*(1.-albice)
      alb=albice
      rnet    = gswnew + xinet

      !if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
      if (globalcellid==targetcell) then
         print *,'globalcellid =',globalcellid
         print *,'no snow on the ground gswnew -',gswnew,'rnet=',rnet
      endif

      call sice(                                               &
!--- input variables
           i,j,iland,isoil,delt,ktau,conflx,nzs,nddzs,         &
           prcpms,rainf,patm,qvatm,qcatm,glw,gswnew,           &
           emiss,rnet,qkms,tkms,rho,myj,globalcellid,          &
!--- sea ice parameters
           tice,rhosice,capice,thdifice,                       &
           zsmain,zshalf,dtdzs,dtdzs2,tbq,                     &
!--- constants
           lv,cp,rovcp,cw,stbolt,tabs,                         &
!--- output variables
           ts1d,dew,soilt,qvg,qsg,qcg,                         &
           eeta,qfx,hfx,s,evapl,prcpl,fltot                    )

      edir1 = eeta*1.e-3
      ec1 = 0.
      ett1 = 0.
      runoff1 = prcpms
      runoff2 = 0.
      infiltr=0.

   endif

!---------------------------------------------------------------
   end subroutine icetmp
!---------------------------------------------------------------

   subroutine sice (                                            &
!--- input variables
              i,j,iland,isoil,delt,ktau,conflx,nzs,nddzs,       &
              prcpms,rainf,patm,qvatm,qcatm,glw,gsw,            &
              emiss,rnet,qkms,tkms,rho,myj,globalcellid,        &
!--- sea ice parameters
              tice,rhosice,capice,thdifice,                     &
              zsmain,zshalf,dtdzs,dtdzs2,tbq,                   &
!--- constants
              xlv,cp,rovcp,cw,stbolt,tabs,                      &
!--- output variables
              tso,dew,soilt,qvg,qsg,qcg,                        &
              eeta,qfx,hfx,s,evapl,prcpl,fltot                  &
                                                                )

!*****************************************************************
!   energy budget and  heat diffusion eqns. for
!   sea ice
!*************************************************************

   implicit none
!-----------------------------------------------------------------

!--- input variables

   integer,  intent(in   )   ::  ktau,nzs, nddzs
   integer,  intent(in   )   ::  i,j,iland,isoil,globalcellid
   real,     intent(in   )   ::  delt,conflx
   logical,  intent(in   )   ::  myj
!--- 3-d atmospheric variables
   real,                                                         &
            intent(in   )    ::                            patm, &
                                                          qvatm, &
                                                          qcatm
!--- 2-d variables
   real,                                                         &
            intent(in   )    ::                             glw, &
                                                            gsw, &
                                                          emiss, &
                                                            rho, &
                                                           qkms, &
                                                           tkms
!--- sea ice properties
   real,    dimension(1:nzs)                                   , &
            intent(in   )    ::                                  &
                                                           tice, &
                                                        rhosice, &
                                                         capice, &
                                                       thdifice


   real,     intent(in   )   ::                                  &
                                                             cw, &
                                                            xlv


   real,     dimension(1:nzs), intent(in)  ::            zsmain, &
                                                         zshalf, &
                                                         dtdzs2

   real,     dimension(1:nddzs), intent(in)  ::           dtdzs

   real,     dimension(1:5001), intent(in)  ::              tbq


!--- input/output variables
!----soil temperature
   real,     dimension( 1:nzs ),  intent(inout)   ::        tso
!-------- 2-d variables
   real,                                                         &
             intent(inout)   ::                             dew, &
                                                           eeta, &
                                                          evapl, &
                                                          prcpl, &
                                                            qvg, &
                                                            qsg, &
                                                            qcg, &
                                                           rnet, &
                                                            qfx, &
                                                            hfx, &
                                                              s, &
                                                          soilt

!--- local variables
   real    ::  x,x1,x2,x4,tn,denom
   real    ::  rainf,  prcpms                                  , &
               tabs, t3, upflux, xinet

   real    ::  cp,rovcp,g0,lv,stbolt,xlmelt,dzstop             , &
               epot,fltot,ft,fq,hft,ras,cvw

   real    ::  fkt,d1,d2,d9,d10,did,r211,r21,r22,r6,r7,d11     , &
               pi,h,fkq,r210,aa,bb,pp,q1,qs1,ts1,tq2,tx2       , &
               tdenom,qgold,snoh

   real    ::  aa1,rhcs, icemelt


   real,     dimension(1:nzs)  ::   cotso,rhtso

   integer ::  nzs1,nzs2,k,k1,kn,kk

!-----------------------------------------------------------------

!-- define constants
   xlmelt=3.35e+5
   cvw=cw

   prcpl=prcpms

   nzs1=nzs-1
   nzs2=nzs-2
   dzstop=1./(zsmain(2)-zsmain(1))
   ras=rho*1.e-3

   do k=1,nzs
      cotso(k)=0.
      rhtso(k)=0.
   enddo

   cotso(1)=0.
   rhtso(1)=tso(nzs)

   do 33 k=1,nzs2
      kn=nzs-k
      k1=2*kn-3
      x1=dtdzs(k1)*thdifice(kn-1)
      x2=dtdzs(k1+1)*thdifice(kn)
      ft=tso(kn)+x1*(tso(kn-1)-tso(kn))                             &
                -x2*(tso(kn)-tso(kn+1))
      denom=1.+x1+x2-x2*cotso(k)
      cotso(k+1)=x1/denom
      rhtso(k+1)=(ft+x2*rhtso(k))/denom
   33 continue

!************************************************************************
!--- the heat balance equation (Smirnova et al., 1996, eq. 21,26)
   rhcs=capice(1)
   h=1.
   fkt=tkms
   d1=cotso(nzs1)
   d2=rhtso(nzs1)
   tn=soilt
   d9=thdifice(1)*rhcs*dzstop
   d10=tkms*cp*rho
   r211=.5*conflx/delt
   r21=r211*cp*rho
   r22=.5/(thdifice(1)*delt*dzstop**2)
   r6=emiss *stbolt*.5*tn**4
   r7=r6/tn
   d11=rnet+r6
   tdenom=d9*(1.-d1+r22)+d10+r21+r7                              &
          +rainf*cvw*prcpms
   fkq=qkms*rho
   r210=r211*rho
   aa=xls*(fkq+r210)/tdenom
   bb=(d10*tabs+r21*tn+xls*(qvatm*fkq                            &
   +r210*qvg)+d11+d9*(d2+r22*tn)                                 &
   +rainf*cvw*prcpms*max(273.15,tabs)                            &
   )/tdenom
   aa1=aa
   pp=patm*1.e3
   aa1=aa1/pp
   !if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
   if (globalcellid==targetcell) then
       print *,' vilka-seaice1'
       print *,'d10,tabs,r21,tn,qvatm,fkq',                          &
                d10,tabs,r21,tn,qvatm,fkq
       print *,'rnet, emiss, stbolt, soilt',rnet, emiss, stbolt, soilt
       print *,'r210,qvg,d11,d9,d2,r22,rainf,cvw,prcpms,tdenom',     &
                r210,qvg,d11,d9,d2,r22,rainf,cvw,prcpms,tdenom
       print *,'tn,aa1,bb,pp,fkq,r210',                              &
                tn,aa1,bb,pp,fkq,r210
   endif
   qgold=qsg
   call vilka(tn,aa1,bb,pp,qs1,ts1,tbq,ktau,i,j,iland,isoil)
!--- it is saturation over sea ice
   qvg=qs1
   qsg=qs1
   tso(1)=min(271.4,ts1)
   qcg=0.
!--- sea ice melting is not included in this simple approach
!--- soilt - skin temperature
   soilt=tso(1)
!---- final solution for soil temperature - tso
   do k=2,nzs
      kk=nzs-k+1
      tso(k)=min(271.4,rhtso(kk)+cotso(kk)*tso(k-1))
   end do
!--- calculation of dew using new value of qsg or transp if no dew
   dew=0.

!--- the diagnostics of surface fluxes
   t3      = stbolt*tn*tn*tn
   upflux  = t3 *0.5*(tn+soilt)
   xinet   = emiss*(glw-upflux)
   hft=-tkms*cp*rho*(tabs-soilt)
   hfx=-tkms*cp*rho*(tabs-soilt)                                &
       *(p1000mb*0.00001/patm)**rovcp
   q1=-qkms*ras*(qvatm - qsg)
   if (q1.le.0.) then
! ---  condensation
      if(myj) then
!-- moisture flux for coupling with myj pbl
         eeta=-qkms*ras*(qvatm/(1.+qvatm) - qsg/(1.+qsg))*1.e3
         if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
            print *,'myj eeta',eeta
         endif
      else ! myj
!-- actual moisture flux from ruc lsm
         dew=qkms*(qvatm-qsg)
         eeta= - rho*dew
         if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
            print *,'ruc lsm eeta',eeta
         endif
      endif ! myj
      qfx= xls*eeta
      eeta= - rho*dew
   else
! ---  evaporation
      if(myj) then
!-- moisture flux for coupling with myj pbl
         eeta=-qkms*ras*(qvatm/(1.+qvatm) - qvg/(1.+qvg))*1.e3
         if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
            print *,'myj eeta',eeta
         endif
      else ! myj
! to convert from m s-1 to kg m-2 s-1: *rho water=1.e3************
!-- actual moisture flux from ruc lsm
         eeta = q1*1.e3
         if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
            print *,'ruc lsm eeta',eeta
         endif
      endif ! myj
      qfx= xls * eeta
      eeta = q1*1.e3
   endif ! q1<0

   evapl=eeta

   s=thdifice(1)*capice(1)*dzstop*(tso(1)-tso(2))
! heat storage in surface layer
   snoh=0.
! there is ice melt
   x= (cp*rho*r211+rhcs*zsmain(2)*0.5/delt)*(soilt-tn) +   &
       xls*rho*r211*(qsg-qgold)
   x=x &
! "heat" from rain
     -rainf*cvw*prcpms*(max(273.15,tabs)-soilt)

!-- excess energy spent on sea ice melt
   icemelt=rnet-xls*eeta -hft -s -x
   if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
      print *,'icemelt=',icemelt
   endif

   fltot=rnet-xls*eeta-hft-s-x-icemelt
   !if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
   if (globalcellid==targetcell) then
      print *,'sice - fltot,rnet,hft,qfx,s,snoh,x=', &
                      fltot,rnet,hft,xls*eeta,s,icemelt,x
   endif

!-------------------------------------------------------------------
   end subroutine sice
!-------------------------------------------------------------------

   subroutine snowseaice(                                       &
              i,j,isoil,delt,ktau,conflx,nzs,nddzs,             &
              rhonewsn,snhei_crit,                              & 
              iland,prcpms,rainf,newsnow,snhei,snwe,snowfrac,   &
              rhosn,patm,qvatm,qcatm,                           &
              glw,gsw,emiss,rnet,                               &
              qkms,tkms,rho,myj,globalcellid,                   &
!--- sea ice parameters
              alb,znt,                                          &
              tice,rhosice,capice,thdifice,                     &
              zsmain,zshalf,dtdzs,dtdzs2,tbq,                   &
!--- constants
              xlv,cp,rovcp,cw,stbolt,tabs,                      &
!--- output variables
              ilnb,snweprint,snheiprint,rsm,tso,                &
              dew,soilt,soilt1,tsnav,qvg,qsg,qcg,               &
              smelt,snoh,snflx,snom,eeta,                       &
              qfx,hfx,s,sublim,prcpl,fltot                      &
                                                                )
!***************************************************************
!   solving energy budget for snow on sea ice and heat diffusion
!   eqns. in snow and sea ice
!***************************************************************


   implicit none
!-------------------------------------------------------------------
!--- input variables

   integer,  intent(in   )   ::  ktau,nzs     ,                  &
                                 nddzs                         !nddzs=2*(nzs-2)
   integer,  intent(in   )   ::  i,j,isoil,globalcellid

   real,     intent(in   )   ::  delt,conflx,prcpms            , &
                                 rainf,newsnow,rhonewsn,         &
                                 snhei_crit
   real                      ::  rhonewcsn

   logical,  intent(in   )   ::  myj
!--- 3-d atmospheric variables
   real,                                                         &
            intent(in   )    ::                            patm, &
                                                          qvatm, &
                                                          qcatm
!--- 2-d variables
   real                                                        , &
            intent(in   )    ::                             glw, &
                                                            gsw, &
                                                            rho, &
                                                           qkms, &
                                                           tkms

!--- sea ice properties
   real,     dimension(1:nzs)                                  , &
            intent(in   )    ::                                  &
                                                           tice, &
                                                        rhosice, &
                                                         capice, &
                                                       thdifice

   real,     intent(in   )   ::                                  &
                                                             cw, &
                                                            xlv

   real,     dimension(1:nzs), intent(in)  ::            zsmain, &
                                                         zshalf, &
                                                         dtdzs2

   real,     dimension(1:nddzs), intent(in)  ::           dtdzs

   real,     dimension(1:5001), intent(in)  ::              tbq

!--- input/output variables
!-------- 3-d soil moisture and temperature
   real,     dimension(  1:nzs )                               , &
             intent(inout)   ::                             tso

   integer,  intent(inout)    ::                           iland


!-------- 2-d variables
   real                                                        , &
             intent(inout)   ::                             dew, &
                                                           eeta, &
                                                          rhosn, &
                                                         sublim, &
                                                          prcpl, &
                                                            alb, &
                                                          emiss, &
                                                            znt, &
                                                            qvg, &
                                                            qsg, &
                                                            qcg, &
                                                            qfx, &
                                                            hfx, &
                                                              s, &
                                                           snwe, &
                                                          snhei, &
                                                          smelt, &
                                                           snom, &
                                                           snoh, &
                                                          snflx, &
                                                          soilt, &
                                                         soilt1, &
                                                       snowfrac, &
                                                          tsnav

   integer, intent(inout)    ::                            ilnb

   real,     intent(out)                    ::              rsm, &
                                                      snweprint, &
                                                     snheiprint
!--- local variables


   integer ::  nzs1,nzs2,k,k1,kn,kk
   real    ::  x,x1,x2,dzstop,ft,tn,denom

   real    ::  snth, newsn                                     , &
               tabs, t3, upflux, xinet                         , &
               beta, snwepr,epdt,pp
   real    ::  cp,rovcp,g0,lv,xlvm,stbolt,xlmelt               , &
               epot,fltot,fq,hft,q1,ras,rhoice,ci,cvw          , &
               riw,deltsn,h

   real    ::  rhocsn,thdifsn,                                   &
               xsn,ddzsn,x1sn,d1sn,d2sn,d9sn,r22sn

   real    ::  cotsn,rhtsn,xsn1,ddzsn1,x1sn1,ftsnow,denomsn
   real    ::  fso,fsn,                                          &
               fkt,d1,d2,d9,d10,did,r211,r21,r22,r6,r7,d11,      &
               fkq,r210,aa,bb,qs1,ts1,tq2,tx2,                   &
               tdenom,aa1,rhcs,h1,tsob, snprim,                  &
               snodif,soh,tnold,qgold,snohgnew
   real,     dimension(1:nzs)  ::  cotso,rhtso

   real                   :: rnet,rsmfrac,soiltfrac,hsn,icemelt,rr
   integer                ::      nmelt


!-----------------------------------------------------------------
   xlmelt=3.35e+5
!-- heat of sublimation of water vapor
   xlvm=xlv+xlmelt

!--- deltsn - is the threshold for splitting the snow layer into 2 layers.
!--- with snow density 400 kg/m^3, this threshold is equal to 7.5 cm,
!--- equivalent to 0.03 m snwe. for other snow densities the threshold is
!--- computed using snwe=0.03 m and current snow density.
!--- snth - the threshold below which the snow layer is combined with
!--- the top sea ice layer. snth is computed using snwe=0.016 m, and
!--- equals 4 cm for snow density 400 kg/m^3.

   deltsn=0.05*1.e3/rhosn
   snth=0.01*1.e3/rhosn

! for 2-layer snow model when the snow depth is marginlly higher than deltsn,
! reset deltsn to half of snow depth.
   if(snhei.ge.deltsn+snth) then
      if(snhei-deltsn-snth.lt.snth) deltsn=0.5*(snhei-snth)
      if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
         print *,'deltsn ice is changed,deltsn,snhei,snth', &
                                   i,j, deltsn,snhei,snth
      endif
   endif

   rhoice=900.
   ci=rhoice*2100.
   ras=rho*1.e-3
   riw=rhoice*1.e-3

   xlmelt=3.35e+5
   rhocsn=2090.* rhosn
   rhonewcsn=2090.* rhonewsn
   thdifsn = 0.265/rhocsn
   ras=rho*1.e-3

   soiltfrac=soilt

   smelt=0.
   soh=0.
   snodif=0.
   snoh=0.
   snohgnew=0.
   rsm = 0.
   rsmfrac = 0.
   fsn=1.
   fso=0.
   cvw=cw

   nzs1=nzs-1
   nzs2=nzs-2

   qgold=qsg
   tnold=soilt
   dzstop=1./(zsmain(2)-zsmain(1))

   snweprint=0.
   snheiprint=0.
   prcpl=prcpms

!*** deltsn is the depth of the top layer of snow where
!*** there is a temperature gradient, the rest of the snow layer
!*** is considered to have constant temperature

   h=1.
   smelt=0.

   fq=qkms
   snhei=snwe*1.e3/rhosn
   snwepr=snwe

!  check if all snow can evaporate during dt
   beta=1.
   epot = -fq*(qvatm-qsg)
   epdt = epot * ras *delt
   if(epdt.gt.0. .and. snwepr.le.epdt) then
     beta=snwepr/max(1.e-8,epdt)
     snwe=0.
   endif

!******************************************************************************
!       coefficients for thomas algorithm for tso
!******************************************************************************

   cotso(1)=0.
   rhtso(1)=tso(nzs)
   do 33 k=1,nzs2
      kn=nzs-k
      k1=2*kn-3
      x1=dtdzs(k1)*thdifice(kn-1)
      x2=dtdzs(k1+1)*thdifice(kn)
      ft=tso(kn)+x1*(tso(kn-1)-tso(kn))                           &
            -x2*(tso(kn)-tso(kn+1))
      denom=1.+x1+x2-x2*cotso(k)
      cotso(k+1)=x1/denom
      rhtso(k+1)=(ft+x2*rhtso(k))/denom
 33 continue
!--- the nzs element in cotso and rhtso will be for snow
!--- there will be 2 layers in snow if it is deeper than deltsn+snth
   if(snhei.ge.snth) then
      if(snhei.le.deltsn+snth) then
!-- 1-layer snow model
         ilnb=1
         snprim=max(snth,snhei)
         soilt1=tso(1)
         tsob=tso(1)
         xsn = delt/2./(zshalf(2)+0.5*snprim)
         ddzsn = xsn / snprim
         x1sn = ddzsn * thdifsn
         x2 = dtdzs(1)*thdifice(1)
         ft = tso(1)+x1sn*(soilt-tso(1))                              &
              -x2*(tso(1)-tso(2))
         denom = 1. + x1sn + x2 -x2*cotso(nzs1)
         cotso(nzs)=x1sn/denom
         rhtso(nzs)=(ft+x2*rhtso(nzs1))/denom
         cotsn=cotso(nzs)
         rhtsn=rhtso(nzs)
!*** average temperature of snow pack (c)
         tsnav=0.5*(soilt+tso(1))                                     &
                     -273.15

      else
!-- 2 layers in snow, soilt1 is temperasture at deltsn depth
         ilnb=2
         snprim=deltsn
         tsob=soilt1
         xsn = delt/2./(0.5*snhei)
         xsn1= delt/2./(zshalf(2)+0.5*(snhei-deltsn))
         ddzsn = xsn / deltsn
         ddzsn1 = xsn1 / (snhei-deltsn)
         x1sn = ddzsn * thdifsn
         x1sn1 = ddzsn1 * thdifsn
         x2 = dtdzs(1)*thdifice(1)
         ft = tso(1)+x1sn1*(soilt1-tso(1))                            &
              -x2*(tso(1)-tso(2))
         denom = 1. + x1sn1 + x2 - x2*cotso(nzs1)
         cotso(nzs)=x1sn1/denom
         rhtso(nzs)=(ft+x2*rhtso(nzs1))/denom
         ftsnow = soilt1+x1sn*(soilt-soilt1)                          &
                 -x1sn1*(soilt1-tso(1))
         denomsn = 1. + x1sn + x1sn1 - x1sn1*cotso(nzs)
         cotsn=x1sn/denomsn
         rhtsn=(ftsnow+x1sn1*rhtso(nzs))/denomsn
!*** average temperature of snow pack (c)
         tsnav=0.5/snhei*((soilt+soilt1)*deltsn                       &
                        +(soilt1+tso(1))*(snhei-deltsn))                 &
                        -273.15
      endif
   endif

   if(snhei.lt.snth.and.snhei.gt.0.) then
!--- snow is too thin to be treated separately, therefore it
!--- is combined with the first sea ice layer.
      snprim=snhei+zsmain(2)
      fsn=snhei/snprim
      fso=1.-fsn
      soilt1=tso(1)
      tsob=tso(2)
      xsn = delt/2./((zshalf(3)-zsmain(2))+0.5*snprim)
      ddzsn = xsn /snprim
      x1sn = ddzsn * (fsn*thdifsn+fso*thdifice(1))
      x2=dtdzs(2)*thdifice(2)
      ft=tso(2)+x1sn*(soilt-tso(2))-                              &
                  x2*(tso(2)-tso(3))
      denom = 1. + x1sn + x2 - x2*cotso(nzs-2)
      cotso(nzs1) = x1sn/denom
      rhtso(nzs1)=(ft+x2*rhtso(nzs-2))/denom
      tsnav=0.5*(soilt+tso(1))                                    &
                  -273.15
      cotso(nzs)=cotso(nzs1)
      rhtso(nzs)=rhtso(nzs1)
      cotsn=cotso(nzs)
      rhtsn=rhtso(nzs)
   endif

!************************************************************************
!--- the heat balance equation
!18apr08 nmelt is the flag for melting, and snoh is heat of snow phase changes
   nmelt=0
   snoh=0.

   epot=-qkms*(qvatm-qsg)
   rhcs=capice(1)
   h=1.
   fkt=tkms
   d1=cotso(nzs1)
   d2=rhtso(nzs1)
   tn=soilt
   d9=thdifice(1)*rhcs*dzstop
   d10=tkms*cp*rho
   r211=.5*conflx/delt
   r21=r211*cp*rho
   r22=.5/(thdifice(1)*delt*dzstop**2)
   r6=emiss *stbolt*.5*tn**4
   r7=r6/tn
   d11=rnet+r6

   if(snhei.ge.snth) then
      if(snhei.le.deltsn+snth) then
!--- 1-layer snow
         d1sn = cotso(nzs)
         d2sn = rhtso(nzs)
      else
!--- 2-layer snow
         d1sn = cotsn
         d2sn = rhtsn
      endif
      d9sn= thdifsn*rhocsn / snprim
      r22sn = snprim*snprim*0.5/(thdifsn*delt)
   endif

   if(snhei.lt.snth.and.snhei.gt.0.) then
!--- thin snow is combined with sea ice
      d1sn = d1
      d2sn = d2
      d9sn = (fsn*thdifsn*rhocsn+fso*thdifice(1)*rhcs)/           &
              snprim
      r22sn = snprim*snprim*0.5                                   &
              /((fsn*thdifsn+fso*thdifice(1))*delt)
   endif

   if(snhei.eq.0.)then
!--- all snow is sublimated
      d9sn = d9
      r22sn = r22
      d1sn = d1
      d2sn = d2
   endif


!---- tdenom for snow
   tdenom = d9sn*(1.-d1sn +r22sn)+d10+r21+r7                    &
         +rainf*cvw*prcpms                                      &
         +rhonewcsn*newsnow/delt

   fkq=qkms*rho
   r210=r211*rho
   aa=xlvm*(beta*fkq+r210)/tdenom
   bb=(d10*tabs+r21*tn+xlvm*(qvatm*                             &
      (beta*fkq)                                                &
      +r210*qvg)+d11+d9sn*(d2sn+r22sn*tn)                        &
      +rainf*cvw*prcpms*max(273.15,tabs)                         &
      + rhonewcsn*newsnow/delt*min(273.15,tabs)                  &
      )/tdenom
   aa1=aa
   pp=patm*1.e3
   aa1=aa1/pp
!18apr08  - the iteration start point
 212 continue
   bb=bb-snoh/tdenom
   !if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
   if (globalcellid==targetcell) then
      print *,'vilka-snow on seaice'
      print *,'tn,aa1,bb,pp,fkq,r210',                          &
               tn,aa1,bb,pp,fkq,r210
      print *,'tabs,qvatm,tn,qvg=',tabs,qvatm,tn,qvg
   endif

   call vilka(tn,aa1,bb,pp,qs1,ts1,tbq,ktau,i,j,iland,isoil)
!--- it is saturation over snow
   qvg=qs1
   qsg=qs1
   qcg=0.

!--- soilt - skin temperature
   soilt=ts1

   !if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
   if (globalcellid==targetcell) then
      print *,' after vilka-snow on seaice'
      print *,' ts1,qs1: ', ts1,qs1
   endif
! solution for temperature at 7.5 cm depth and snow-seaice interface
   if(snhei.ge.snth) then
      if(snhei.gt.deltsn+snth) then
!-- 2-layer snow model
         soilt1=min(273.15,rhtsn+cotsn*soilt)
         tso(1)=min(271.4,(rhtso(nzs)+cotso(nzs)*soilt1))
         tsob=soilt1
      else
!-- 1 layer in snow
         tso(1)=min(271.4,(rhtso(nzs)+cotso(nzs)*soilt))
         soilt1=tso(1)
         tsob=tso(1)
      endif
   elseif  (snhei > 0. .and. snhei < snth) then
! blended
      tso(2)=min(271.4,(rhtso(nzs1)+cotso(nzs1)*soilt))
      tso(1)=min(271.4,(tso(2)+(soilt-tso(2))*fso))
      soilt1=tso(1)
      tsob=tso(2)
   else
! snow is melted
      tso(1)=min(271.4,soilt)
      soilt1=min(271.4,soilt)
      tsob=tso(1)
   endif
!---- final solution for tso in sea ice
   if (snhei > 0. .and. snhei < snth) then
! blended or snow is melted
      do k=3,nzs
         kk=nzs-k+1
         tso(k)=min(271.4,rhtso(kk)+cotso(kk)*tso(k-1))
      end do
   else
      do k=2,nzs
         kk=nzs-k+1
         tso(k)=min(271.4,rhtso(kk)+cotso(kk)*tso(k-1))
      end do
   endif

   if(nmelt.eq.1) go to 220

! if all snow can evaporate, then there is nothing to melt
   if(soilt.gt.273.15.and.snwepr-beta*epot*ras*delt.gt.0..and.snhei.gt.0.) then
!
      nmelt = 1
      soiltfrac=snowfrac*273.15+(1.-snowfrac)*min(271.4,soilt)

      qsg= qsn(soiltfrac,tbq)/pp
      t3      = stbolt*tnold*tnold*tnold
      upflux  = t3 * 0.5*(tnold+soiltfrac)
      xinet   = emiss*(glw-upflux)
      epot = -qkms*(qvatm-qsg)
      q1=epot*ras

      if (q1.le.0.) then
! ---  condensation
         dew=-epot
         qfx= xlvm*rho*dew
         eeta=qfx/xlvm
      else
! ---  evaporation
         eeta = q1 * beta *1.e3
! to convert from kg m-2 s-1 to m s-1: 1/rho water=1.e-3************
         qfx= - xlvm * eeta
      endif

      hfx=d10*(tabs-soiltfrac)

      if(snhei.ge.snth)then
         soh=thdifsn*rhocsn*(soiltfrac-tsob)/snprim
         snflx=soh
      else
         soh=(fsn*thdifsn*rhocsn+fso*thdifice(1)*rhcs)*                &
             (soiltfrac-tsob)/snprim
         snflx=soh
      endif
      x= (r21+d9sn*r22sn)*(soiltfrac-tnold) +                         &
          xlvm*r210*(qsg-qgold)
!-- snoh is energy flux of snow phase change
      snoh=rnet+qfx +hfx                                              &
               +rhonewcsn*newsnow/delt*(min(273.15,tabs)-soiltfrac)   &
               -soh-x+rainf*cvw*prcpms*                               &
               (max(273.15,tabs)-soiltfrac)

      !if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
      if (globalcellid==targetcell) then
         print *,'snowseaice melt i,j,snoh,rnet,qfx,hfx,soh,x',       &
                                  i,j,snoh,rnet,qfx,hfx,soh,x
         print *,'rainf*cvw*prcpms*(max(273.15,tabs)-soiltfrac)',     &
                  rainf*cvw*prcpms*(max(273.15,tabs)-soiltfrac)
      endif
      snoh=amax1(0.,snoh)
!-- smelt is speed of melting in m/s
      smelt= snoh /xlmelt*1.e-3
      smelt=amin1(smelt,snwepr/delt-beta*epot*ras)
      smelt=amax1(0.,smelt)

      if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
         print *,'1-smelt i,j',smelt,i,j
      endif
!18apr08 - egglston limit
      smelt= amin1 (smelt, 5.6e-8*max(1.,(soilt-273.15)))
      if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
         print *,'2-smelt i,j',smelt,i,j
      endif

! rr - potential melting
      rr=snwepr/delt-beta*epot*ras
      smelt=min(smelt,rr)
      if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
         print *,'3- smelt i,j,smelt,rr',i,j,smelt,rr
      endif
      snohgnew=smelt*xlmelt*1.e3
      snodif=amax1(0.,(snoh-snohgnew))
 
      snoh=snohgnew

      !if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
      if (globalcellid==targetcell) then
         print*,'soiltfrac,soilt,snohgnew,snodif=', &
             i,j,soiltfrac,soilt,snohgnew,snodif
         print *,'snoh,snodif',snoh,snodif
      endif

!*** from koren et al. (1999) 13% of snow melt stays in the snow pack
      rsmfrac=min(0.18,(max(0.08,snwepr/0.10*0.13)))
      if(snhei > 0.01) then
         rsm=rsmfrac*smelt*delt
      else
! do n ot keep melted water if snow depth is less that 1 cm
         rsm=0.
      endif
!18apr08 rsm is part of melted water that stays in snow as liquid
      smelt=amax1(0.,smelt-rsm/delt)

      !if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
      if (globalcellid==targetcell) then
         print *,'4-smelt i,j,smelt,rsm,snwepr,rsmfrac', &
                          i,j,smelt,rsm,snwepr,rsmfrac
      endif

!-- update liquid equivalent of snow depth
!-- for evaporation and snow melt
      snwe = amax1(0.,(snwepr-                                      &
                     (smelt+beta*epot*ras)*delt                        &
                                            ) )
      soilt=soiltfrac
!--- if there is no snow melting then just evaporation
!--- or condensation changes snwe
   else
      if(snhei.ne.0.) then
         epot=-qkms*(qvatm-qsg)
         snwe = amax1(0.,(snwepr-                               &
                beta*epot*ras*delt))
      endif
   endif

! no iteration for snow on sea ice, because it will produce
! skin temperature higher than it is possible with snow on sea ice
!      if(nmelt.eq.1) goto 212  ! second iteration
 220 continue

   if(smelt > 0..and.  rsm > 0.) then
      if(snwe.le.rsm) then
         !if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
         if (globalcellid==targetcell) then
            print *,'seaice snwe<rsm snwe,rsm,smelt*delt,epot*ras*delt,beta', &
                            snwe,rsm,smelt*delt,epot*ras*delt,beta
         endif
      else
!*** update snow density on effect of snow melt, melted
!*** from the top of the snow. 13% of melted water
!*** remains in the pack and changes its density.
!*** eq. 9 (with my correction) in koren et al. (1999)

         xsn=(rhosn*(snwe-rsm)+1.e3*rsm)/ snwe
         rhosn=min(max(58.8,xsn),500.)

         rhocsn=2090.* rhosn
         thdifsn = 0.265/rhocsn
      endif
   endif

   snweprint=snwe
!                                              &
!--- if vegfrac.ne.0. then some snow stays on the canopy
!--- and should be added to snwe for water conservation
! 4 nov 07                    +vegfrac*cst
   snheiprint=snweprint*1.e3 / rhosn

   !if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
   if (globalcellid==targetcell) then
      print *, 'snweprint : ',snweprint
      print *, 'd9sn,soilt,tsob : ', d9sn,soilt,tsob
   endif

   if(snhei.gt.0.) then
      if(ilnb.gt.1) then
         tsnav=0.5/snhei*((soilt+soilt1)*deltsn                     &
                   +(soilt1+tso(1))*(snhei-deltsn))                 &
                   -273.15
      else
         tsnav=0.5*(soilt+tso(1)) - 273.15
      endif
   endif
!--- recalculation of dew using new value of qsg
   dew=0.
   pp=patm*1.e3
   qsg= qsn(soilt,tbq)/pp
   epot = -fq*(qvatm-qsg)
   if(epot.lt.0.) then
! sublimation
      dew=-epot
   endif

   snom=snom+smelt*delt*1.e3

!--- the diagnostics of surface fluxes

   t3      = stbolt*tnold*tnold*tnold
   upflux  = t3 *0.5*(soilt+tnold)
   xinet   = emiss*(glw-upflux)
   hft=-tkms*cp*rho*(tabs-soilt)
   hfx=-tkms*cp*rho*(tabs-soilt)                        &
            *(p1000mb*0.00001/patm)**rovcp
   q1 = - fq*ras* (qvatm - qsg)
   if (q1.lt.0.) then
! ---  condensation
      if(myj) then
!-- moisture flux for coupling with myj pbl
         eeta=-qkms*ras*(qvatm/(1.+qvatm) - qsg/(1.+qsg))*1.e3
      else ! myj
!-- actual moisture flux from ruc lsm
         dew=qkms*(qvatm-qsg)
         eeta= - rho*dew
      endif ! myj
      qfx= xlvm*eeta
      eeta= - rho*dew
      sublim = eeta
   else
! ---  evaporation
      if(myj) then
!-- moisture flux for coupling with myj pbl
         eeta=-qkms*ras*beta*(qvatm/(1.+qvatm) - qvg/(1.+qvg))*1.e3
      else ! myj
! to convert from m s-1 to kg m-2 s-1: *rho water=1.e3************
!-- actual moisture flux from ruc lsm
         eeta = q1*beta*1.e3
      endif ! myj
      qfx= xlvm * eeta
      eeta = q1*beta*1.e3
      sublim = eeta
   endif

   icemelt=0.
   if(snhei.ge.snth)then
      s=thdifsn*rhocsn*(soilt-tsob)/snprim
      snflx=s
   elseif(snhei.lt.snth.and.snhei.gt.0.) then
      s=(fsn*thdifsn*rhocsn+fso*thdifice(1)*rhcs)*                &
              (soilt-tsob)/snprim
      snflx=s
      if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
          print *,'snow is thin, snflx',i,j,snflx
      endif
   else
      snflx=d9sn*(soilt-tsob)
      if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
         print *,'snow is gone, snflx',i,j,snflx
      endif
   endif

   snhei=snwe *1.e3 / rhosn

   if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
      print *,'snhei,snoh',i,j,snhei,snoh
   endif
!
   x= (r21+d9sn*r22sn)*(soilt-tnold) +                            &
       xlvm*r210*(qsg-qgold)
   !if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
   if (globalcellid==targetcell) then
      print *,'snowseaice storage ',i,j,x
      print *,'r21,d9sn,r22sn,soiltfrac,tnold,qsg,qgold,snprim',  &
               r21,d9sn,r22sn,soiltfrac,tnold,qsg,qgold,snprim
   endif
   x=x &
   -rhonewcsn*newsnow/delt*(min(273.15,tabs)-soilt)               &
   -rainf*cvw*prcpms*(max(273.15,tabs)-soilt)

! -- excess energy is spent on ice melt
   icemelt = rnet-hft-xlvm*eeta-s-snoh-x
   if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
      print *,'snowseaice icemelt=',icemelt
   endif

   fltot=rnet-hft-xlvm*eeta-s-snoh-x-icemelt
   !if ( wrf_at_debug_level(iceruc_dbg_lvl) ) then
   if (globalcellid==targetcell) then
      print *,'i,j,snhei,qsg,soilt,soilt1,tso,tabs,qvatm', &
               i,j,snhei,qsg,soilt,soilt1,tso,tabs,qvatm
      print *,'snowseaice - fltot,rnet,hft,qfx,s,snoh,icemelt,snodif,x,soilt=' &
                   ,fltot,rnet,hft,xlvm*eeta,s,snoh,icemelt,snodif,x,soilt
   endif
!-- restore sea-ice parameters if snow is less than threshold
   if(snhei.eq.0.)  then
      tsnav=soilt-273.15
      emiss=0.98
      znt=0.011
      alb=0.55
   endif

!------------------------------------------------------------------------
   end subroutine snowseaice
!------------------------------------------------------------------------

   subroutine vilka(tn,d1,d2,pp,qs,ts,tt,nstep,ii,j,iland,isoil)
!--------------------------------------------------------------
!--- vilka finds the solution of energy budget at the surface
!--- using table t,qs computed from clausius-klapeiron
!--------------------------------------------------------------
   real,     dimension(1:5001),  intent(in   )   ::  tt
   real,     intent(in  )   ::  tn,d1,d2,pp
   integer,  intent(in  )   ::  nstep,ii,j,iland,isoil

   real,     intent(out  )  ::  qs, ts

   real    ::  f1,t1,t2,rn
   integer ::  i,i1

       i=(tn-1.7315e2)/.05+1
       t1=173.1+float(i)*.05
       f1=t1+d1*tt(i)-d2
       i1=i-f1/(.05+d1*(tt(i+1)-tt(i)))
       i=i1
       if(i.gt.5000.or.i.lt.1) goto 1
  10   i1=i
       t1=173.1+float(i)*.05
       f1=t1+d1*tt(i)-d2
       rn=f1/(.05+d1*(tt(i+1)-tt(i)))
       i=i-int(rn)
       if(i.gt.5000.or.i.lt.1) goto 1
       if(i1.ne.i) goto 10
       ts=t1-.05*rn
       qs=(tt(i)+(tt(i)-tt(i+1))*rn)/pp
       goto 20
!   1   print *,'crash in surface energy budget - stop'
!   1   print *,'     avost in vilka     table index= ',i
    1   print *,d1,d2,pp,tt(i),tt(i+1)
    !   write(0,*),'i,j=',ii,j,'lu_index = ',iland, 'psfc[hpa] = ',pp, 'tsfc = ',tn,&
    !             'tti=',tt(i),'ttip1=',tt(i+1),'d1=',d1
       fatal_error ('  crash in surface energy budget  ' )
   20  continue
!----------------------------------------------------------------
   end subroutine vilka
!----------------------------------------------------------------   

   function qsn(tn,t)
!----------------------------------------------------------------   
   real,     dimension(1:5001),  intent(in   )   ::  t
   real,     intent(in  )   ::  tn

   real    qsn, r,r1,r2
   integer i

       r=(tn-173.15)/.05+1.
       i=int(r)
       if(i.ge.1) goto 10
       i=1
       r=1.
  10   if(i.le.5000) goto 20
       i=5000
       r=5001.
  20   r1=t(i)
       r2=r-i
       qsn=(t(i+1)-r1)*r2 + r1
!-----------------------------------------------------------------------
   end function qsn
!------------------------------------------------------------------------

end module module_ruc_ice
