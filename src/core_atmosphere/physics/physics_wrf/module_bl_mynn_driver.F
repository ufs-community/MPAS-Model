!=================================================================================================================
 module module_bl_mynn_driver
 use mpas_kind_types,only: kind_phys => RKIND
 use mpas_log

 use module_bl_mynn,only: mynnedmf
 use bl_mynn_post,only: bl_mynn_post_run
 use bl_mynn_pre,only: bl_mynn_pre_run

 implicit none
 private
 public:: mynn_bl_driver


 contains


!=================================================================================================================
 subroutine mynn_bl_driver &
                 (ids               , ide               , jds                , jde                , &
                  kds               , kde               , ims                , ime                , &
                  jms               , jme               , kms                , kme                , &
                  its               , ite               , jts                , jte                , &
                  kts               , kte               , f_qc               , f_qi               , &
                  f_qs              , f_qoz             , f_nc               , f_ni               , &
                  f_nifa            , f_nwfa            , f_nbca             , initflag           , &
                  do_restart        , do_DAcycling      , icloud_bl          , delt               , &
                  dx                , xland             , ps                 , ts                 , &
                  qsfc              , ust               , ch                 , hfx                , &
                  qfx               , rmol              , wspd               , znt                , &
                  uoce              , voce              , dz                 , u                  , &
                  v                 , w                 , th                 , tt                 , &
                  p                 , exner             , rho                , qv                 , &
                  qc                , qi                , qs                 , nc                 , &
                  ni                , nifa              , nwfa               , nbca               , &
                  qoz               , rthraten          , pblh               , kpbl               , &
                  cldfra_bl         , qc_bl             , qi_bl              , maxwidth           , &
                  maxmf                                 , ztop_plume         , qke                , &
                  qke_adv           , tsq               , qsq                , cov                , &
                  el_pbl            , rublten           , rvblten            , rthblten           , &
                  rqvblten          , rqcblten          , rqiblten           , rqsblten           , &
                  rncblten          , rniblten          , rnifablten         , rnwfablten         , &
                  rnbcablten        , rqozblten         , edmf_a             , edmf_w             , &
                  edmf_qt           , edmf_thl          , edmf_ent           , edmf_qc            , &
                  sub_thl           , sub_sqv           , det_thl            , det_sqv            , &
                  exch_h            , exch_m            , dqke               , qwt                , &
                  qshear            , qbuoy             , qdiss              , sh3d               , &
                  sm3d              , spp_pbl           , pattern_spp        ,                      &
                  bl_mynn_tkeadvect , bl_mynn_tkebudget , bl_mynn_cloudpdf   , bl_mynn_mixlength  , &
                  bl_mynn_closure   , bl_mynn_stfunc    , bl_mynn_topdown    , bl_mynn_scaleaware , &
                  bl_mynn_dheat_opt , bl_mynn_edmf      , bl_mynn_edmf_dd    , bl_mynn_edmf_mom   , &
                  bl_mynn_edmf_tke  , bl_mynn_output    , bl_mynn_mixscalars , bl_mynn_cloudmix   , &
                  bl_mynn_mixqt     , errmsg            , errflg                                    &
#if(WRF_CHEM == 1)
                 ,mix_chem   , nchem        , kdvel       , ndvel        , chem3d        , vd3d   , &
                  frp_mean   , emis_ant_no                                                          &
#endif
               )

!=================================================================================================================

!--- input arguments:
 logical,intent(in):: &
    f_qc,               &! if true,the physics package includes the cloud liquid water mixing ratio.
    f_qi,               &! if true,the physics package includes the cloud ice mixing ratio.
    f_qs,               &! if true,the physics package includes the snow mixing ratio.
    f_qoz,              &! if true,the physics package includes the ozone mixing ratio.
    f_nc,               &! if true,the physics package includes the cloud liquid water number concentration.
    f_ni,               &! if true,the physics package includes the cloud ice number concentration.
    f_nifa,             &! if true,the physics package includes the "ice-friendly" aerosol number concentration.
    f_nwfa,             &! if true,the physics package includes the "water-friendly" aerosol number concentration.
    f_nbca               ! if true,the physics package includes the number concentration of black carbon.

 logical,intent(in):: &
    bl_mynn_tkeadvect    !

 logical,intent(in):: &
    do_restart,         &!
    do_DAcycling         !

 integer,intent(in):: &
    ids,ide,jds,jde,kds,kde, &
    ims,ime,jms,jme,kms,kme, &
    its,ite,jts,jte,kts,kte

 integer,intent(in):: &
    bl_mynn_cloudpdf,   &!
    bl_mynn_mixlength,  &!
    bl_mynn_stfunc,     &!
    bl_mynn_topdown,    &!
    bl_mynn_scaleaware, &!
    bl_mynn_dheat_opt,  &!
    bl_mynn_edmf,       &!
    bl_mynn_edmf_dd,    &!
    bl_mynn_edmf_mom,   &!
    bl_mynn_edmf_tke,   &!
    bl_mynn_output,     &!
    bl_mynn_mixscalars, &!
    bl_mynn_cloudmix,   &!
    bl_mynn_mixqt,      &!
    bl_mynn_tkebudget    !
 
 integer,intent(in):: &
    initflag,           &!
    icloud_bl,          &!
    spp_pbl              !

 real(kind=kind_phys),intent(in):: &
    bl_mynn_closure

 real(kind=kind_phys),intent(in):: &
    delt                 !

 real(kind=kind_phys),intent(in),dimension(ims:ime,jms:jme):: &
    dx,                 &!
    xland,              &!
    ps,                 &!
    ts,                 &!
    qsfc,               &!
    ust,                &!
    ch,                 &!
    hfx,                &!
    qfx,                &!
    rmol,               &!
    wspd,               &!
    uoce,               &!
    voce,               &!
    znt                  !

 real(kind=kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme):: &
    dz,      &!
    u,       &!
    w,       &!
    v,       &!
    th,      &!
    tt,      &!
    p,       &!
    exner,   &!
    rho,     &!
    qv,      &!
    rthraten  !

 real(kind=kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    qc,      &!
    qi,      &!
    qs,      &!
    qoz,     &!
    nc,      &!
    ni,      &!
    nifa,    &!
    nwfa,    &!
    nbca

 real(kind=kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    pattern_spp   !


!--- inout arguments:
 integer,intent(inout),dimension(ims:ime,jms:jme):: &
    kpbl

 real(kind=kind_phys),intent(inout),dimension(ims:ime,jms:jme):: &
    pblh          !

 real(kind=kind_phys),intent(inout),dimension(ims:ime,kms:kme,jms:jme):: &
    cldfra_bl,   &!
    qc_bl,       &!
    qi_bl         !

 real(kind=kind_phys),intent(inout),dimension(ims:ime,kms:kme,jms:jme):: &
    el_pbl,      &!
    qke,         &!
    qke_adv,     &!
    cov,         &!
    qsq,         &!
    tsq,         &!
    sh3d,        &!
    sm3d

 real(kind=kind_phys),intent(inout),dimension(ims:ime,kms:kme,jms:jme):: &
    rublten,     &!
    rvblten,     &!
    rthblten,    &!
    rqvblten      !

 real(kind=kind_phys),intent(inout),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    rqcblten,    &!
    rqiblten,    &!
    rqsblten,    &!
    rqozblten,   &!
    rncblten,    &!
    rniblten,    &!
    rnifablten,  &!
    rnwfablten,  &!
    rnbcablten    !

 real(kind=kind_phys),intent(inout),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    edmf_a,      &!
    edmf_w,      &!
    edmf_qt,     &!
    edmf_thl,    &!
    edmf_ent,    &!
    edmf_qc,     &!
    sub_thl,     &!
    sub_sqv,     &!
    det_thl,     &!
    det_sqv       !


!--- output arguments:
 character(len=*),intent(out):: &
    errmsg        ! output error message (-).

 integer,intent(out):: &
    errflg        ! output error flag (-).

 real(kind=kind_phys),intent(out),dimension(ims:ime,jms:jme):: &
    maxwidth,    &!
    maxmf,       &!
    ztop_plume

 real(kind=kind_phys),intent(out),dimension(ims:ime,kms:kme,jms:jme):: &
    exch_h,      &!
    exch_m        !

 real(kind=kind_phys),intent(out),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    dqke,        &!
    qwt,         &!
    qshear,      &!
    qbuoy,       &!
    qdiss         !

#if(WRF_CHEM == 1)
!--- input arguments for PBL and free-tropospheric mixing of chemical species:
 logical,intent(in):: mix_chem
 integer,intent(in):: kdvel,nchem,ndvel

 real(kind=kind_phys),intent(in),dimension(ims:ime,jms:jme):: frp_mean,ems_ant_no
 real(kind=kind_phys),intent(in),dimension(ims:ime,kdvel,jms:jme,ndvel):: vd3d
 real(kind=kind_phys),intent(inout),dimension(ims:ime,kms:kme,jms:jme,nchem):: chem3
 logical, parameter :: &
    rrfs_sd    =.false.,  &
    smoke_dbg  =.false.,  &
    enh_mix    =.false.
#else
 logical, parameter :: &
    mix_chem   =.false.,  &
    enh_mix    =.false.,  &
    rrfs_sd    =.false.,  &
    smoke_dbg  =.false.
 integer, parameter :: nchem=2, ndvel=2, kdvel=1, num_vert_mix = 1
#endif
 real(kind=kind_phys):: frp_v,emisant_no_v
 real(kind=kind_phys),dimension(kdvel,ndvel):: vd_v
 real(kind=kind_phys),dimension(kts:kte,nchem):: chem_v

!local variables and arrays:
 logical:: mynn_edmf_l,mynn_edmf_dd_l,mynn_edmf_mom_l,mynn_edmf_tke_l
 logical:: mynn_mixscalars_l,mynn_mixclouds_l,mynn_mixqt_l
 logical:: mynn_tkebudget_l
 logical:: mynn_output_l,mynn_dheatopt_l,mynn_scaleaware_l,mynn_topdown_l

 integer:: i,k,j

 integer:: dheat_opt
 integer:: kpbl_v

 real(kind=kind_phys):: denom

 real(kind=kind_phys):: &
    dx_v,xland_v,ps_v,ts_v,qsfc_v,ust_v,ch_v,hfx_v,qfx_v, &
    rmol_v,wspd_v,uoce_v,voce_v,znt_v

 real(kind=kind_phys),dimension(kts:kte):: &
    dz_v,u_v,v_v,th_v,tt_v,p_v,exner_v,rho_v,qv_v,rthraten_v

 real(kind=kind_phys),dimension(kts:kme):: &
    w_v

 real(kind=kind_phys),dimension(kts:kte):: &
    qc_v,qi_v,qs_v,nc_v,ni_v,nifa_v,nwfa_v,nbca_v,qoz_v

 real(kind=kind_phys),dimension(kts:kte):: &
    pattern_spp_v

 real(kind=kind_phys):: &
    pblh_v

 real(kind=kind_phys),dimension(kts:kte):: &
    cldfrabl_v,qcbl_v,qibl_v,elpbl_v,qke_v,qkeadv_v,cov_v,qsq_v,tsq_v,sh3d_v,sm3d_v

 real(kind=kind_phys),dimension(kts:kte):: &
    rublten_v,rvblten_v,rthblten_v,rqvblten_v,rqcblten_v,rqiblten_v,rqsblten_v, &
    rncblten_v,rniblten_v,rnifablten_v,rnwfablten_v,rnbcablten_v,rqozblten_v

 real(kind=kind_phys),dimension(kts:kte):: &
    edmfa_v,edmfw_v,edmfqt_v,edmfthl_v,edmfent_v,edmfqc_v, &
    subthl_v,subsqv_v,detthl_v,detsqv_v

 real(kind=kind_phys):: &
    maxwidth_v,maxmf_v,ztopplume_v

 real(kind=kind_phys),dimension(kts:kte):: &
    exchh_v,exchm_v,dqke_v,qwt_v,qshear_v,qbuoy_v,qdiss_v

 real(kind=kind_phys),dimension(kts:kte):: &
    sqv_v,sqc_v,sqi_v,sqs_v

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
!call mpas_log_write('--- enter subroutine mynn_bl_driver:')

 errmsg = " "
 errflg = 0

 mynn_edmf_l     = .false.
 mynn_edmf_dd_l  = .false.
 mynn_edmf_mom_l = .false.
 mynn_edmf_tke_l = .false.
 if(bl_mynn_edmf     == 1) mynn_edmf_l     = .true.
 if(bl_mynn_edmf_dd  == 1) mynn_edmf_dd_l  = .true.
 if(bl_mynn_edmf_mom == 1) mynn_edmf_mom_l = .true.
 if(bl_mynn_edmf_tke == 1) mynn_edmf_tke_l = .true.

 mynn_mixscalars_l = .false.
 mynn_mixclouds_l  = .false.
 mynn_mixqt_l      = .false.
 if(bl_mynn_mixscalars == 1) mynn_mixscalars_l = .true.
 if(bl_mynn_cloudmix   == 1) mynn_mixclouds_l  = .true.
 if(bl_mynn_mixqt      == 1) mynn_mixqt_l       = .true.

 mynn_tkebudget_l = .false.
 if(bl_mynn_tkebudget == 1) mynn_tkebudget_l = .true.

 mynn_output_l     = .false.
 mynn_dheatopt_l   = .false.
 mynn_scaleaware_l = .false.
 mynn_topdown_l    = .false.
 if(bl_mynn_output     == 1) mynn_output_l     = .true.
 if(bl_mynn_dheat_opt  == 1) mynn_dheatopt_l   = .true.
 if(bl_mynn_scaleaware == 1) mynn_scaleaware_l = .true.
 if(bl_mynn_topdown    == 1) mynn_topdown_l    = .true.

 dheat_opt = bl_mynn_dheat_opt

 do j = jts,jte
 do i = its,ite
     
    !--- input arguments
    dx_v    = dx(i,j)
    xland_v = xland(i,j)
    ps_v    = ps(i,j)
    ts_v    = ts(i,j)
    qsfc_v  = qsfc(i,j)
    ust_v   = ust(i,j)
    ch_v    = ch(i,j)
    hfx_v   = hfx(i,j)
    qfx_v   = qfx(i,j)
    rmol_v  = rmol(i,j)
    wspd_v  = wspd(i,j)
    uoce_v  = uoce(i,j)
    voce_v  = voce(i,j)
    znt_v   = znt(i,j)

    do k = kts,kte
       dz_v(k)       = dz(i,k,j)
       u_v(k)        = u(i,k,j)
       v_v(k)        = v(i,k,j)
       w_v(k)        = w(i,k,j)
       th_v(k)       = th(i,k,j)
       tt_v(k)       = tt(i,k,j)
       p_v(k)        = p(i,k,j)
       exner_v(k)    = exner(i,k,j)
       rho_v(k)      = rho(i,k,j)
       qv_v(k)       = qv(i,k,j)
       rthraten_v(k) = rthraten(i,k,j)
    enddo
    w_v(kte+1) = w(i,kte+1,j)

    !--- input arguments for cloud mixing ratios and number concentrations; input argument
    !    for the ozone mixing ratio; input arguments for aerosols from the aerosol-aware
    !    Thompson cloud microphysics:
    do k = kts,kte
       qc_v(k)   = 0._kind_phys
       qi_v(k)   = 0._kind_phys
       qs_v(k)   = 0._kind_phys
       qoz_v(k)  = 0._kind_phys
       nc_v(k)   = 0._kind_phys
       ni_v(k)   = 0._kind_phys
       nifa_v(k) = 0._kind_phys
       nwfa_v(k) = 0._kind_phys
       nbca_v(k) = 0._kind_phys
    enddo
    if(f_qc .and. present(qc)) then
       do k = kts,kte
          qc_v(k) = qc(i,k,j)
       enddo
    endif
    if(f_qi .and. present(qi)) then
       do k = kts,kte
          qi_v(k) = qi(i,k,j)
       enddo
    endif
    if(f_qs .and. present(qs)) then
       do k = kts,kte
          qs_v(k) = qs(i,k,j)
       enddo
    endif
    if(f_nc .and. present(nc)) then
       do k = kts,kte
          nc_v(k) = nc(i,k,j)
       enddo
    endif
    if(f_ni .and. present(ni)) then
       do k = kts,kte
          ni_v(k) = ni(i,k,j)
       enddo
    endif
    if(f_nifa .and. present(nifa)) then
       do k = kts,kte
          nifa_v(k) = nifa(i,k,j)
       enddo
    endif
    if(f_nwfa .and. present(nwfa)) then
       do k = kts,kte
          nwfa_v(k) = nwfa(i,k,j)
       enddo
    endif
    if(f_nbca .and. present(nbca)) then
       do k = kts,kte
          nbca_v(k) = nbca(i,k,j)
       enddo
    endif
    if(f_qoz .and. present(qoz)) then
       do k = kts,kte
          qoz_v(k) = qoz(i,k,j)
       enddo
    endif

    !--- conversion from mixing ratios to specific contents:
    call bl_mynn_pre_run(its,ite,kte,f_qc,f_qi,f_qs,qv_v,qc_v,qi_v,qs_v,sqv_v,sqc_v, &
                         sqi_v,sqs_v,errmsg,errflg)

    !--- initialization of the stochastic forcing in the PBL:
    if(spp_pbl > 0 .and. present(pattern_spp)) then
       do k = kts,kte
          pattern_spp_v(k) = pattern_spp(i,k,j)
       enddo
    else
       do k = kts,kte
          pattern_spp_v(k) = 0._kind_phys
       enddo
    endif

    !--- inout arguments:
    pblh_v = pblh(i,j)
    kpbl_v = kpbl(i,j)

    do k = kts,kte
       cldfrabl_v(k) = cldfra_bl(i,k,j)
       qcbl_v(k)     = qc_bl(i,k,j)
       qibl_v(k)     = qi_bl(i,k,j)
    enddo

    do k = kts,kte
       elpbl_v(k)  = el_pbl(i,k,j)
       qke_v(k)    = qke(i,k,j)
       qkeadv_v(k) = qke_adv(i,k,j)
       cov_v(k)    = cov(i,k,j)
       tsq_v(k)    = tsq(i,k,j)   
       qsq_v(k)    = qsq(i,k,j)
       sh3d_v(k)   = sh3d(i,k,j)
       sm3d_v(k)   = sm3d(i,k,j)
    enddo

#if(WRF_CHEM == 1)
    do ic = 1,nchem
       do k = kts,kte
          chem_v(k,ic) = chem3d(i,k,j,ic)
       enddo
    enddo
    do ic = 1,ndvel
       do k = 1,kdvel
          vd_v(k,ic) = vd3d(i,k,j,ic)
       enddo
    enddo
    frp_v        = frp_mean(i,j)
    emisant_no_v = emis_ant_no(i,j)
#else
    chem_v       = 0.0
    vd_v         = 0.0
    frp_v        = 0.0
    emisant_no_v = 0.0
#endif

    do k = kts,kte
       rqcblten_v(k)   = 0._kind_phys
       rqiblten_v(k)   = 0._kind_phys
       rqsblten_v(k)   = 0._kind_phys
       rqozblten_v(k)  = 0._kind_phys
       rncblten_v(k)   = 0._kind_phys
       rniblten_v(k)   = 0._kind_phys
       rnifablten_v(k) = 0._kind_phys
       rnwfablten_v(k) = 0._kind_phys
       rnbcablten_v(k) = 0._kind_phys
    enddo

    call mynnedmf( &
            i               = i             , j           = j             ,                              &
            initflag        = initflag      , restart     = do_restart    , cycling     = do_DAcycling , &
            delt            = delt          , dz1         = dz_v          , dx          = dx_v         , &
            znt             = znt_v         , u1          = u_v           , v1          = v_v          , &
            w1              = w_v           , th1         = th_v          , sqv1        = sqv_v        , &
            sqc1            = sqc_v         , sqi1        = sqi_v         , sqs1        = sqs_v        , &
            qnc1            = nc_v          , qni1        = ni_v          , qnwfa1      = nwfa_v       , &
            qnifa1          = nifa_v        , qnbca1      = nbca_v        , ozone1      = qoz_v        , &
            p1              = p_v           , ex1         = exner_v       , rho1        = rho_v        , &
            tk1             = tt_v          , xland       = xland_v       , ts          = ts_v         , &
            qsfc            = qsfc_v        , ps          = ps_v          , ust         = ust_v        , &
            ch              = ch_v          , hfx         = hfx_v         , qfx         = qfx_v        , &
            rmol            = rmol_v        , wspd        = wspd_v        , uoce        = uoce_v       , &
            voce            = voce_v        , qke1        = qke_v         , qke_adv1    = qkeadv_v     , &
            tsq1            = tsq_v         , qsq1        = qsq_v         , cov1        = cov_v        , &
            rthraten1       = rthraten_v    , du1         = rublten_v     , dv1         = rvblten_v    , &
            dth1            = rthblten_v    , dqv1        = rqvblten_v    , dqc1        = rqcblten_v   , &
            dqi1            = rqiblten_v    , dqs1        = rqsblten_v    , dqnc1       = rncblten_v   , &
            dqni1           = rniblten_v    , dqnwfa1     = rnwfablten_v  , dqnifa1     = rnifablten_v , &
            dqnbca1         = rnbcablten_v  , dozone1     = rqozblten_v   , kh1         = exchh_v      , &
            km1             = exchm_v       , pblh        = pblh_v        , kpbl        = kpbl_v       , &
            el1             = elpbl_v       , dqke1       = dqke_v        , qwt1        = qwt_v        , &
            qshear1         = qshear_v      , qbuoy1      = qbuoy_v       , qdiss1      = qdiss_v      , &
            sh1             = sh3d_v        , sm1         = sm3d_v        , qc_bl1      = qcbl_v       , &
            qi_bl1          = qibl_v        , cldfra_bl1  = cldfrabl_v    , &
            edmf_a1         = edmfa_v       , edmf_w1     = edmfw_v       , edmf_qt1    = edmfqt_v     , &
            edmf_thl1       = edmfthl_v     , edmf_ent1   = edmfent_v     , edmf_qc1    = edmfqc_v     , &
            sub_thl1        = subthl_v      , sub_sqv1    = subsqv_v      , det_thl1    = detthl_v     , &
            det_sqv1        = detsqv_v      ,                                                            &
            maxwidth        = maxwidth_v    , maxmf       = maxmf_v       , ztop_plume  = ztopplume_v  , &
            flag_qc         = f_qc          , flag_qi     = f_qi          , flag_qs     = f_qs         , &
            flag_ozone      = f_qoz         , flag_qnc    = f_nc          , flag_qni    = f_ni         , &
            flag_qnwfa      = f_nwfa        , flag_qnifa  = f_nifa        , flag_qnbca  = f_nbca       , &
            pattern_spp_pbl1= pattern_spp_v                                                              &
!#if(WRF_CHEM == 1)
           ,mix_chem  = mix_chem  , enh_mix = enh_mix , rrfs_sd     = rrfs_sd        , &
            smoke_dbg = smoke_dbg , nchem   = nchem   , kdvel       = kdvel          , &
            ndvel     = ndvel     , chem    = chem_v  , emis_ant_no = emisant_no_v   , &
            frp       = frp_v     , vdep    = vd_v                                     &
!#endif
           ,bl_mynn_tkeadvect  = bl_mynn_tkeadvect    , &
            tke_budget         = bl_mynn_tkebudget    , &
            bl_mynn_cloudpdf   = bl_mynn_cloudpdf     , &
            bl_mynn_mixlength  = bl_mynn_mixlength    , &
!            bl_mynn_stfunc     = bl_mynn_stfunc       , &
!            bl_mynn_dheatopt   = mynn_dheatopt_l      , &
!            bl_mynn_scaleaware = mynn_scaleaware_l    , &
!            bl_mynn_topdown    = mynn_topdown_l       , &
            closure            = bl_mynn_closure      , &
            bl_mynn_edmf       = bl_mynn_edmf         , &
!            bl_mynn_edmf_dd    = mynn_edmf_dd         , &
            bl_mynn_edmf_mom   = bl_mynn_edmf_mom     , &
            bl_mynn_edmf_tke   = bl_mynn_edmf_tke     , &
            bl_mynn_mixscalars = bl_mynn_mixscalars   , &
            bl_mynn_output     = bl_mynn_output       , &
            bl_mynn_cloudmix   = bl_mynn_cloudmix     , &
            bl_mynn_mixqt      = bl_mynn_mixqt        , &
            icloud_bl          = icloud_bl            , &
            spp_pbl            = spp_pbl              , &
            kts = kts , kte = kte , errmsg = errmsg , errflg = errflg )


    !--- conversion of tendencies in terms of specific contents to in terms of mixing ratios:
    call  bl_mynn_post_run(its,ite,kte,f_qc,f_qi,f_qs,delt,qv_v,qc_v,qi_v,qs_v,rqvblten_v,rqcblten_v, &
                           rqiblten_v,rqsblten_v,errmsg,errflg)

    !--- inout arguments:
    pblh(i,j)  = pblh_v
    kpbl(i,j)  = kpbl_v
    do k = kts,kte
       cldfra_bl(i,k,j) = cldfrabl_v(k)
       qc_bl(i,k,j)     = qcbl_v(k)
       qi_bl(i,k,j)     = qibl_v(k)
    enddo

    do k = kts,kte
       el_pbl(i,k,j)  = elpbl_v(k)
       qke(i,k,j)     = qke_v(k)
       qke_adv(i,k,j) = qkeadv_v(k)
       cov(i,k,j)     = cov_v(k)
       tsq(i,k,j)     = tsq_v(k)
       qsq(i,k,j)     = qsq_v(k)
       sh3d(i,k,j)    = sh3d_v(k)
       sm3d(i,k,j)    = sm3d_v(k)
    enddo

    !--- inout tendencies:
    do k = kts,kte
       rublten(i,k,j)    = rublten_v(k) 
       rvblten(i,k,j)    = rvblten_v(k) 
       rthblten(i,k,j)   = rthblten_v(k) 
       rqvblten(i,k,j)   = rqvblten_v(k) 
    enddo
    if(f_qc .and. present(rqcblten)) then
       do k = kts,kte
          rqcblten(i,k,j) = rqcblten_v(k) 
       enddo
    endif
    if(f_qi .and. present(rqiblten)) then
       do k = kts,kte
          rqiblten(i,k,j) = rqiblten_v(k) 
       enddo
    endif
    if(f_qs .and. present(rqsblten)) then
       do k = kts,kte
          rqsblten(i,k,j) = rqsblten_v(k)
       enddo
    endif
    if(f_qoz .and. present(rqozblten)) then
       do k = kts,kte
          rqozblten(i,k,j) = rqozblten_v(k) 
       enddo
    endif
    if(f_nc .and. present(rncblten)) then
       do k = kts,kte
          rncblten(i,k,j) = rncblten_v(k) 
       enddo
    endif
    if(f_ni .and. present(rniblten)) then
       do k = kts,kte
          rniblten(i,k,j) = rniblten_v(k) 
       enddo
    endif
    if(f_nifa .and. present(rnifablten)) then
       do k = kts,kte
          rnifablten(i,k,j) = rnifablten_v(k) 
       enddo
    endif
    if(f_nwfa .and. present(rnwfablten)) then
       do k = kts,kte
          rnwfablten(i,k,j) = rnwfablten_v(k) 
       enddo
    endif
    if(f_nbca .and. present(rnbcablten)) then
       do k = kts,kte
          rnbcablten(i,k,j) = rnbcablten_v(k) 
       enddo
    endif

    do k = kts,kte
       edmf_a(i,k,j)   = edmfa_v(k)
       edmf_w(i,k,j)   = edmfw_v(k)
       edmf_qt(i,k,j)  = edmfqt_v(k)
       edmf_thl(i,k,j) = edmfthl_v(k)
       edmf_ent(i,k,j) = edmfent_v(k)
       edmf_qc(i,k,j)  = edmfqc_v(k)
       sub_thl(i,k,j)  = subthl_v(k)
       sub_sqv(i,k,j)  = subsqv_v(k)
       det_thl(i,k,j)  = detthl_v(k)
       det_sqv(i,k,j)  = detsqv_v(k)
    enddo

    !--- output arguments:
    maxwidth(i,j)   = maxwidth_v
    maxmf(i,j)      = maxmf_v
    ztop_plume(i,j) = ztopplume_v

    do k = kts,kte
       exch_h(i,k,j) = exchh_v(k)
       exch_m(i,k,j) = exchm_v(k)
    enddo

    if(present(qwt)   .and. present(qbuoy) .and. present(qshear) .and. &
       present(qdiss) .and. present(dqke)) then
       do k = kts,kte
          dqke(i,k,j)   = dqke_v(k)
          qwt(i,k,j)    = qwt_v(k)
          qshear(i,k,j) = qshear_v(k)
          qbuoy(i,k,j)  = qbuoy_v(k)
          qdiss(i,k,j)  = qdiss_v(k)
       enddo
    endif

#if (WRF_CHEM == 1)
    if (mix_chem) then
       do ic = 1,nchem
          do k = kts,kte
             chem3d(i,k,j,ic) = max(1.e-12, chem(k,ic))
          enddo
       enddo
    endif
#endif

 enddo !i
 enddo !j

!call mpas_log_write('--- end subroutine mynn_bl_driver:')

 end subroutine mynn_bl_driver

!=================================================================================================================
 end module module_bl_mynn_driver
!=================================================================================================================

