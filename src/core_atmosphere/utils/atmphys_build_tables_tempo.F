! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module atmphys_build_tables_tempo
 use module_mp_tempo, only : tempo_init
 use module_mp_tempo_utils, only : qi_aut_qs, freezeH2O, qr_acr_qs, qr_acr_qg
 use module_mp_tempo_params, only : dimNRHG, build_hail_aware_table, &
      tcg_racg, tmr_racg, tcr_gacr, tnr_racg, tnr_gacr, &
      tcs_racs1, tmr_racs1, tcs_racs2, tmr_racs2, tcr_sacr1, tms_sacr1, tcr_sacr2, tms_sacr2, &
      tnr_racs1, tnr_racs2, tnr_sacr1, tnr_sacr2, &
      tpi_qrfz, tni_qrfz, tpg_qrfz, tnr_qrfz, tpi_qcfz, tni_qcfz, &
      tpi_ide, tps_iaus, tni_iaus
 
 implicit none
 private
 public:: build_tables_tempo_mp

!builds the files containing the look-up tables for the TEMPO cloud microphysics scheme.


 contains


!=================================================================================================================
 subroutine build_tables_tempo_mp
!=================================================================================================================

 use mpas_io_units, only : mpas_new_unit, mpas_release_unit

!local variables:
 logical, parameter:: l_mp_tables = .false.
 integer:: istatus
 integer:: mp_unit

!-----------------------------------------------------------------------------------------------------------------
!--- partial initialization before building the look-up tables:
!--- hail_aware_flag defined in module_mp_tempo_params

 call tempo_init(l_mp_tables=l_mp_tables, hail_aware_flag=build_hail_aware_table)

 call mpas_new_unit(mp_unit, unformatted = .true.)

!--- building look-up table for rain collecting graupel:
 write(0,*)
 if (build_hail_aware_table) then
    write(0,*) '--- building MP_TEMPO_HAILAWARE_QRacrQG_DATA.DBL'
    open(unit=mp_unit,file='MP_TEMPO_HAILAWARE_QRacrQG_DATA.DBL',form='unformatted',status='new',iostat=istatus)
    if (istatus /= 0) then
       call print_parallel_mesg('MP_TEMPO_HAILAWARE_QRacrQG_DATA.DBL')
       return
    end if
 else
    write(0,*) '--- building MP_TEMPO_QRacrQG_DATA.DBL'
    open(unit=mp_unit,file='MP_TEMPO_QRacrQG_DATA.DBL',form='unformatted',status='new',iostat=istatus)
    if (istatus /= 0) then
       call print_parallel_mesg('MP_TEMPO_QRacrQG_DATA.DBL')
       return
    end if
 endif

 call qr_acr_qg(dimNRHG)
 write(mp_unit) tcg_racg
 write(mp_unit) tmr_racg
 write(mp_unit) tcr_gacr
! write(mp_unit) tmg_gacr
 write(mp_unit) tnr_racg
 write(mp_unit) tnr_gacr
 close(unit=mp_unit)

!--- building look-up table for rain collecting snow:
 write(0,*)
 write(0,*) '--- building MP_TEMPO_QRacrQS_DATA.DBL'
 open(unit=mp_unit,file='MP_TEMPO_QRacrQS_DATA.DBL',form='unformatted',status='new',iostat=istatus)
 if (istatus /= 0) then
    call print_parallel_mesg('MP_TEMPO_QRacrQS_DATA.DBL')
    return
 end if
 call qr_acr_qs
 write(mp_unit)tcs_racs1
 write(mp_unit)tmr_racs1
 write(mp_unit)tcs_racs2
 write(mp_unit)tmr_racs2
 write(mp_unit)tcr_sacr1
 write(mp_unit)tms_sacr1
 write(mp_unit)tcr_sacr2
 write(mp_unit)tms_sacr2
 write(mp_unit)tnr_racs1
 write(mp_unit)tnr_racs2
 write(mp_unit)tnr_sacr1
 write(mp_unit)tnr_sacr2
 close(unit=mp_unit)

!--- building look-up table for freezing of cloud droplets:
 write(0,*)
 write(0,*) '--- building MP_TEMPO_freezeH2O_DATA.DBL'
 open(unit=mp_unit,file='MP_TEMPO_freezeH2O_DATA.DBL',form='unformatted',status='new',iostat=istatus)
 if (istatus /= 0) then
    call print_parallel_mesg('MP_TEMPO_freezeH2O_DATA.DBL')
    return
 end if
 call freezeH2O
 write(mp_unit) tpi_qrfz
 write(mp_unit) tni_qrfz
 write(mp_unit) tpg_qrfz
 write(mp_unit) tnr_qrfz
 write(mp_unit) tpi_qcfz
 write(mp_unit) tni_qcfz
 close(unit=mp_unit)
 
!--- building look-up table for autoconversion of cloud ice to snow:
 write(0,*)
 write(0,*) '--- building MP_TEMPO_QIautQS_DATA.DBL'
 open(unit=mp_unit,file='MP_TEMPO_QIautQS_DATA.DBL',form='unformatted',status='new',iostat=istatus)
 if (istatus /= 0) then
    call print_parallel_mesg('MP_TEMPO_QIautQS_DATA.DBL')
    return
 end if
 call qi_aut_qs
 write(mp_unit) tpi_ide
 write(mp_unit) tps_iaus
 write(mp_unit) tni_iaus
 close(unit=mp_unit)
 call mpas_release_unit(mp_unit)

 write(0,*)
 write(0,*) 'Finished building all tables.'
 write(0,*)
 write(0,*) '*******************************************************************************'
 write(0,*)
 if (build_hail_aware_table) then
    write(0,*) '  MP_TEMPO_HAILAWARE_QRacrQG_DATA.DBL'
 else
    write(0,*) '  MP_TEMPO_QRacrQG_DATA.DBL'
 endif
 write(0,*) '  MP_TEMPO_QRacrQS_DATA.DBL'
 write(0,*) '  MP_TEMPO_freezeH2O_DATA.DBL'
 write(0,*) '  MP_TEMPO_QIautQS_DATA.DBL'
 write(0,*)
 write(0,*) '*******************************************************************************'

 end subroutine build_tables_tempo_mp


!=================================================================================================================
 subroutine print_parallel_mesg(filename)
!=================================================================================================================

 character(len=*), intent(in) :: filename

 write(0,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 write(0,*) '! Error encountered while trying to create new file '//trim(filename)
 write(0,*) '! '
 write(0,*) '! Please ensure that this file does not exist before running ''build_tables'','
 write(0,*) '! and ensure that ''build_tables'' is *NOT* run in parallel.'
 write(0,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

 end subroutine print_parallel_mesg

!=================================================================================================================
 end module atmphys_build_tables_tempo
!=================================================================================================================
