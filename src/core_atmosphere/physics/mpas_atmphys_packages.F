! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_packages
 use mpas_kind_types
 use mpas_derived_types,only : mpas_pool_type,mpas_io_context_type,MPAS_LOG_ERR, MPAS_streamInfo_type
 use mpas_pool_routines,only : mpas_pool_get_config,mpas_pool_get_package
 use mpas_log,only : mpas_log_write

 implicit none
 private
 public:: atmphys_setup_packages

!mpas_atmphys_packages contains the definitions of all physics packages.
!Laura D. Fowler (laura@ucar.edu) / 2016-03-10.


 contains


!=================================================================================================================
 function atmphys_setup_packages(configs,streamInfo,packages,iocontext) result(ierr)
!=================================================================================================================

!inout arguments:
 type (mpas_pool_type), intent(inout) :: configs
 type (MPAS_streamInfo_type), intent(in) :: streamInfo    ! SH (Mar-17-2024): Not used in this function. 
 type (mpas_pool_type), intent(inout) :: packages
 type (mpas_io_context_type), intent(inout) :: iocontext

!local variables:
 character(len=StrKIND),pointer:: config_microp_scheme
 character(len=StrKIND),pointer:: config_nssl_moments
 character(len=StrKIND),pointer:: config_convection_scheme
 character(len=StrKIND),pointer:: config_pbl_scheme
 character(len=StrKIND),pointer:: config_lsm_scheme
 logical,pointer:: config_tempo_hailaware, config_tempo_aerosolaware
 logical,pointer:: mp_kessler_in,mp_thompson_in,mp_thompson_aers_in,mp_wsm6_in,mp_nssl2m_in,nssl3m_in
 logical,pointer:: mp_tempo_in, tempo_hailaware_in, tempo_aerosolaware_in
 logical,pointer:: cu_grell_freitas_in,cu_kain_fritsch_in,cu_ntiedtke_in
 logical,pointer:: bl_mynn_in,bl_mynnedmf_in,bl_ysu_in,bl_myj_in
 logical,pointer:: lsm_noah_in,lsm_ruc_in
 logical,pointer:: sf_noahmp_in

 integer :: ierr

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('')
!call mpas_log_write('--- enter subroutine atmphys_setup_packages:')

 ierr = 0

 call mpas_log_write('----- Setting up package variables -----')
 call mpas_log_write('')

!--- initialization of all packages for parameterizations of cloud microphysics:

 call mpas_pool_get_config(configs,'config_microp_scheme',config_microp_scheme)

 call mpas_pool_get_config(configs,'config_nssl_moments',config_nssl_moments)
 call mpas_pool_get_config(configs,'config_tempo_hailaware',config_tempo_hailaware)

 nullify(mp_kessler_in)
 call mpas_pool_get_package(packages,'mp_kessler_inActive',mp_kessler_in)

 nullify(mp_thompson_in)
 call mpas_pool_get_package(packages,'mp_thompson_inActive',mp_thompson_in)

 nullify(mp_thompson_aers_in)
 call mpas_pool_get_package(packages,'mp_thompson_aers_inActive',mp_thompson_aers_in)

 nullify(mp_tempo_in)
 call mpas_pool_get_package(packages,'mp_tempo_inActive',mp_tempo_in)

 nullify(tempo_aerosolaware_in)
 call mpas_pool_get_package(packages,'tempo_aerosolaware_inActive',tempo_aerosolaware_in)
 
 nullify(tempo_hailaware_in)
 call mpas_pool_get_package(packages,'tempo_hailaware_inActive',tempo_hailaware_in)
 
 nullify(mp_wsm6_in)
 call mpas_pool_get_package(packages,'mp_wsm6_inActive',mp_wsm6_in)

 nullify(mp_nssl2m_in)
 call mpas_pool_get_package(packages,'mp_nssl2m_inActive',mp_nssl2m_in)

 nullify(nssl3m_in)
 call mpas_pool_get_package(packages,'nssl3m_inActive',nssl3m_in)

 if(.not.associated(mp_kessler_in      )  .or. &
    .not.associated(mp_thompson_in     ) .or. &
    .not.associated(mp_thompson_aers_in) .or. &
    .not.associated(mp_tempo_in     ) .or. &
    .not.associated(mp_nssl2m_in) .or. &
    .not.associated(mp_wsm6_in)) then
    call mpas_log_write('====================================================================================',messageType=MPAS_LOG_ERR)
    call mpas_log_write('* Error while setting up packages for cloud microphysics options in atmosphere core.',messageType=MPAS_LOG_ERR)
    call mpas_log_write('====================================================================================',messageType=MPAS_LOG_ERR)
    ierr = 1
    return
 endif

 mp_kessler_in       = .false.
 mp_thompson_in      = .false.
 mp_thompson_aers_in = .false.
 mp_tempo_in         = .false.
 tempo_aerosolaware_in  = .false.
 tempo_hailaware_in  = .false.
 mp_wsm6_in          = .false.
 mp_nssl2m_in            = .false.
 nssl3m_in               = .false.

 if(config_microp_scheme == 'mp_kessler') then
    mp_kessler_in = .true.
 elseif(config_microp_scheme == 'mp_thompson') then
    mp_thompson_in = .true.
 elseif(config_microp_scheme == 'mp_thompson_aerosols') then
    mp_thompson_aers_in = .true.
 elseif(config_microp_scheme == 'mp_wsm6') then
    mp_wsm6_in = .true.
 elseif(config_microp_scheme == 'mp_tempo') then
    mp_tempo_in = .true.
    if (config_tempo_aerosolaware) then
       tempo_aerosolaware_in = .true.
    endif
    if (config_tempo_hailaware) then
       tempo_hailaware_in = .true.
    endif
 elseif(config_microp_scheme == 'mp_nssl2m') then
    mp_nssl2m_in = .true.
    IF ( config_nssl_moments == 'nssl3m' ) THEN
      nssl3m_in = .true.
    ENDIF
 endif

 call mpas_log_write('    mp_kessler_in           = $l', logicArgs=(/mp_kessler_in/))
 call mpas_log_write('    mp_thompson_in          = $l', logicArgs=(/mp_thompson_in/))
 call mpas_log_write('    mp_thompson_aers_in     = $l', logicArgs=(/mp_thompson_aers_in/))
 call mpas_log_write('    mp_tempo_in             = $l', logicArgs=(/mp_tempo_in/))
 call mpas_log_write('    tempo_aerosolaware_in   = $l', logicArgs=(/tempo_aerosolaware_in/))
 call mpas_log_write('    tempo_hailaware_in      = $l', logicArgs=(/tempo_hailaware_in/))
 call mpas_log_write('    mp_wsm6_in              = $l', logicArgs=(/mp_wsm6_in/))
 call mpas_log_write('    mp_nssl2m_in            = $l', logicArgs=(/mp_nssl2m_in/))
 call mpas_log_write('    nssl3m_in               = $l', logicArgs=(/nssl3m_in/))

!--- initialization of all packages for parameterizations of convection:

 call mpas_pool_get_config(configs,'config_convection_scheme',config_convection_scheme)

 nullify(cu_grell_freitas_in)
 call mpas_pool_get_package(packages,'cu_grell_freitas_inActive',cu_grell_freitas_in)

 nullify(cu_kain_fritsch_in)
 call mpas_pool_get_package(packages,'cu_kain_fritsch_inActive',cu_kain_fritsch_in)

 nullify(cu_ntiedtke_in)
 call mpas_pool_get_package(packages,'cu_ntiedtke_inActive',cu_ntiedtke_in)

 if(.not.associated(cu_grell_freitas_in) .or. &
    .not.associated(cu_kain_fritsch_in)  .or. &
    .not.associated(cu_ntiedtke_in)    ) then
    call mpas_log_write('====================================================================================',messageType=MPAS_LOG_ERR)
    call mpas_log_write('* Error while setting up packages for convection options in atmosphere core.',        messageType=MPAS_LOG_ERR)
    call mpas_log_write('====================================================================================',messageType=MPAS_LOG_ERR)
    ierr = 1
    return
 endif

 cu_grell_freitas_in = .false.
 cu_kain_fritsch_in  = .false.
 cu_ntiedtke_in      = .false.

 if(config_convection_scheme=='cu_grell_freitas') then
    cu_grell_freitas_in = .true.
 elseif(config_convection_scheme == 'cu_kain_fritsch') then
    cu_kain_fritsch_in = .true.
 elseif(config_convection_scheme == 'cu_tiedtke' .or. &
        config_convection_scheme == 'cu_ntiedtke') then
    cu_ntiedtke_in = .true.
 endif

 call mpas_log_write('    cu_grell_freitas_in     = $l', logicArgs=(/cu_grell_freitas_in/))
 call mpas_log_write('    cu_kain_fritsch_in      = $l', logicArgs=(/cu_kain_fritsch_in/))
 call mpas_log_write('    cu_ntiedtke_in          = $l', logicArgs=(/cu_ntiedtke_in/))

!--- initialization of all packages for parameterizations of surface layer and planetary boundary layer:

 call mpas_pool_get_config(configs,'config_pbl_scheme',config_pbl_scheme)

 nullify(bl_mynn_in)
 call mpas_pool_get_package(packages,'bl_mynn_inActive',bl_mynn_in)

 nullify(bl_mynnedmf_in)
 call mpas_pool_get_package(packages,'bl_mynnedmf_inActive',bl_mynnedmf_in)
      
 nullify(bl_ysu_in)
 call mpas_pool_get_package(packages,'bl_ysu_inActive',bl_ysu_in)

 nullify(bl_myj_in)
 call mpas_pool_get_package(packages,'bl_myj_inActive',bl_myj_in)

 if(.not.associated(bl_mynn_in) .or. &
    .not.associated(bl_mynnedmf_in) .or. &
    .not.associated(bl_ysu_in) .or. &
    .not.associated(bl_myj_in)) then
    call mpas_log_write('====================================================================================',messageType=MPAS_LOG_ERR)
    call mpas_log_write('* Error while setting up packages for planetary layer  options in atmosphere core.',  messageType=MPAS_LOG_ERR)
    call mpas_log_write('====================================================================================',messageType=MPAS_LOG_ERR)
    ierr = 1
    return
 endif

 bl_mynn_in = .false.
 bl_mynnedmf_in = .false. 
 bl_ysu_in  = .false.
 bl_myj_in = .false.
 if(config_pbl_scheme=='bl_mynn') then
    bl_mynn_in = .true.
 elseif(config_pbl_scheme == 'bl_mynnedmf') then
    bl_mynnedmf_in = .true.
 elseif(config_pbl_scheme == 'bl_ysu') then
    bl_ysu_in = .true.
 elseif(config_pbl_scheme == 'bl_myj') then
    bl_myj_in = .true.
 endif

 call mpas_log_write('    bl_mynn_in              = $l', logicArgs=(/bl_mynn_in/))
 call mpas_log_write('    bl_mynnedmf_in          = $l', logicArgs=(/bl_mynnedmf_in/))
 call mpas_log_write('    bl_ysu_in               = $l', logicArgs=(/bl_ysu_in/))
 call mpas_log_write('    bl_myj_in               = $l', logicArgs=(/bl_myj_in/))
 call mpas_log_write('')

!--- initialization of all packages for parameterizations of the land surface:

  call mpas_pool_get_config(configs,'config_lsm_scheme', config_lsm_scheme)
  
  nullify(lsm_noah_in)
  call mpas_pool_get_package(packages,'lsm_noah_inActive', lsm_noah_in)

  nullify(lsm_ruc_in)
  call mpas_pool_get_package(packages,'lsm_ruc_inActive', lsm_ruc_in)

 if(.not.associated(lsm_noah_in) .or. &
    .not.associated(lsm_ruc_in)) then
    call mpas_log_write('================================================================================',messageType=MPAS_LOG_ERR)
    call mpas_log_write('* Error while setting up packages for land surface options in atmosphere core.',  messageType=MPAS_LOG_ERR)
    call mpas_log_write('================================================================================',messageType=MPAS_LOG_ERR)
    ierr = 1
    return
 endif

 lsm_noah_in = .false.
 lsm_ruc_in = .false.
 
 if(config_lsm_scheme=='sf_noah') then
   lsm_noah_in = .true.
 elseif(config_lsm_scheme=='sf_ruc') then
   lsm_ruc_in = .true.
 endif

 call mpas_log_write('    lsm_noah_in              = $l', logicArgs=(/lsm_noah_in/))
 call mpas_log_write('    lsm_ruc_in               = $l', logicArgs=(/lsm_ruc_in/))
 call mpas_log_write('')

!--- initialization of all packages for parameterizations of land surface processes:

 call mpas_pool_get_config(configs,'config_lsm_scheme',config_lsm_scheme)

 nullify(sf_noahmp_in)
 call mpas_pool_get_package(packages,'sf_noahmp_inActive',sf_noahmp_in)

 if(.not.associated(sf_noahmp_in)) then
    call mpas_log_write('====================================================================================',messageType=MPAS_LOG_ERR)
    call mpas_log_write('* Error while setting up packages for land surface options in atmosphere core.'    ,  messageType=MPAS_LOG_ERR)
    call mpas_log_write('====================================================================================',messageType=MPAS_LOG_ERR)
    ierr = 1
    return
 endif

 if(config_lsm_scheme=='sf_noahmp') then
    sf_noahmp_in = .true.
 endif

 call mpas_log_write('    sf_noahmp_in            = $l', logicArgs=(/sf_noahmp_in/))
 call mpas_log_write('')


 end function atmphys_setup_packages

!=================================================================================================================
 end module mpas_atmphys_packages
!=================================================================================================================



