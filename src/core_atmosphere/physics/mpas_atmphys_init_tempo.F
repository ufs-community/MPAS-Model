! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
#define DM_BCAST_MACRO(A) call mpas_dmpar_bcast_real4s(dminfo,size(A),A)

!=================================================================================================================
 module mpas_atmphys_init_tempo
 use mpas_dmpar
 use mpas_kind_types
 use mpas_log
 use mpas_pool_routines

 use mpas_atmphys_utilities
 use module_mp_tempo_params, only: ntb_arc,ntb_arw,ntb_art,ntb_arr,ntb_ark,tnccn_act, &
      naCCN1,naCCN0,naIN0,naIN1,nwfa_default,aero_max
                              
 implicit none
 private
 public:: init_tempo_aerosols_forMPAS

!MPAS main initialization of the TEMPO parameterization of cloud microphysics with nucleation of cloud
!droplets based on distributions of CCNs and INs (aerosol-aware parameterization).

 
 contains


!=================================================================================================================
 subroutine init_tempo_aerosols_forMPAS(do_restart,dminfo,mesh,state,time_lev,diag_physics)
!=================================================================================================================

!input variables:
 type(dm_info),intent(in):: dminfo
 type(mpas_pool_type),intent(in):: mesh
 logical,intent(in):: do_restart
 integer,intent(in):: time_lev

!inout variables:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: state

!local variables and pointers:
 integer,pointer:: nCellsSolve,nVertLevels
 integer,pointer:: index_nifa,index_nwfa

 real(kind=RKIND),dimension(:),pointer    :: areaCell
 real(kind=RKIND),dimension(:),pointer    :: nifa2d,nwfa2d
 real(kind=RKIND),dimension(:,:),pointer  :: zgrid,zz
 real(kind=RKIND),dimension(:,:),pointer  :: rho_zz,nifa,nwfa
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

 character(len=StrKIND):: mess

 integer:: iCell, k

 real(kind=RKIND):: max_test
 real(kind=RKIND):: airmass, deltaz
 real(kind=RKIND):: h_01
 real(kind=RKIND):: niIN3,niCCN3
 real(kind=RKIND):: nifa_max,nifa_min,global_nifa_max,global_nifa_min
 real(kind=RKIND):: nwfa_max,nwfa_min,global_nwfa_max,global_nwfa_min
 real(kind=RKIND),dimension(:,:),allocatable:: hgt

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('--- enter subroutine init_tempo_aerosols_forMPAS:')

!... read a static file containing CCN activation of aerosols. The data were created from a parcel model by
!... Feingold & Heymsfield with further changes by Eidhammer and Kriedenweis.
 call table_ccnAct(dminfo)
 call mpas_log_write('--- end read table_ccnAct:')

!... if do_restart is true, then we do not need to check the initialization of nwfa, nifa, and nwfa2d. If false,
!    then, we proceed with the initialization:
 if(do_restart) return

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_array(mesh,'areaCell',areaCell)
 call mpas_pool_get_array(mesh,'zgrid'   ,zgrid   )
 call mpas_pool_get_array(mesh,'zz'      ,zz      )

 call mpas_pool_get_array(diag_physics,'nifa2d',nifa2d)
 call mpas_pool_get_array(diag_physics,'nwfa2d',nwfa2d)

 call mpas_pool_get_dimension(state,'index_nifa'  ,index_nifa  )
 call mpas_pool_get_dimension(state,'index_nwfa'  ,index_nwfa  )

 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 nifa   => scalars(index_nifa,:,:)
 nwfa   => scalars(index_nwfa,:,:)

 call mpas_pool_get_array(state,'rho_zz',rho_zz,time_lev)

 if(.not.allocated(hgt)) allocate(hgt(1:nVertLevels,1:nCellsSolve))
 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       hgt(k,iCell) = 0.5_RKIND * (zgrid(k,iCell)+zgrid(k+1,iCell))
    enddo
 enddo

!... initialize the distribution of hygroscopic ("water friendly") aerosols if not already initialized using
!    GOCART data:
 global_nwfa_min = 0._RKIND
 global_nwfa_max = 0._RKIND
 nwfa_min = minval(nwfa(:,1:nCellsSolve))
 nwfa_max = maxval(nwfa(:,1:nCellsSolve))
 call mpas_dmpar_min_real(dminfo,nwfa_min,global_nwfa_min)
 call mpas_dmpar_max_real(dminfo,nwfa_max,global_nwfa_max)
 call mpas_log_write('--- global_nwfa_min = $r',realArgs=(/global_nwfa_min/))
 call mpas_log_write('--- global_nwfa_max = $r',realArgs=(/global_nwfa_max/))

 if(global_nwfa_min == 0._RKIND .and. global_nwfa_max == 0._RKIND) then
    call mpas_log_write('--- initialize nwfa using an exponential distribution of CCN as a function of height.')
    do iCell = 1, nCellsSolve
       if(hgt(1,iCell).le.1000.0) then
          h_01 = 0.8
       elseif(hgt(1,iCell).ge.2500.0) then
          h_01 = 0.01
       else
          h_01 = 0.8*cos(hgt(1,iCell)*0.001 - 1.0)
       endif
       niCCN3 = -1.0*ALOG(naCCN1/naCCN0)/h_01
       nwfa(1,iCell) = naCCN1+naCCN0*exp(-((hgt(2,iCell)-hgt(1,iCell))/1000.)*niCCN3)
       do k = 2, nVertLevels
          nwfa(k,iCell) = naCCN1+naCCN0*exp(-((hgt(k,iCell)-hgt(1,iCell))/1000.)*niCCN3)
       enddo
    enddo
 else
    call mpas_log_write('--- initialize nwfa using the climatological GOCART data.')
 endif

!... initialize the distribution of nonhygroscopic ("ice friendly") aerosols if not already initialized using
!    GOCART data:
 global_nifa_min = 0._RKIND
 global_nifa_max = 0._RKIND
 nifa_min = minval(nifa(:,1:nCellsSolve))
 nifa_max = maxval(nifa(:,1:nCellsSolve))
 call mpas_dmpar_min_real(dminfo,nifa_min,global_nifa_min)
 call mpas_dmpar_max_real(dminfo,nifa_max,global_nifa_max)
 call mpas_log_write('--- global_nifa_min = $r',realArgs=(/global_nifa_min/))
 call mpas_log_write('--- global_nifa_max = $r',realArgs=(/global_nifa_max/))

 if(global_nifa_min == 0._RKIND .and. global_nifa_max == 0._RKIND) then
    call mpas_log_write('--- initialize nifa using an exponential distribution of IN as a function of height.')
    do iCell = 1, nCellsSolve
       if(hgt(1,iCell).le.1000.0) then
          h_01 = 0.8
       elseif(hgt(1,iCell).ge.2500.0) then
          h_01 = 0.01
       else
          h_01 = 0.8*cos(hgt(1,iCell)*0.001 - 1.0)
       endif
       niIN3 = -1.0*ALOG(naIN1/naIN0)/h_01
       nifa(1,iCell) = naIN1+naIN0*exp(-((hgt(2,iCell)-hgt(1,iCell))/1000.)*niIN3)
       do k = 2, nVertLevels
          nifa(k,iCell) = naIN1+naIN0*exp(-((hgt(k,iCell)-hgt(1,iCell))/1000.)*niIN3)
       enddo
    enddo
 else
    call mpas_log_write('--- initialize nifa using the climatological GOCART data.')
 endif

!... scale the lowest level aerosol data into an emissions rate.  This is very far from ideal, but
!... need higher emissions where larger amount of (climo) existing and lesser emissions where there
!... exists fewer to begin as a first-order simplistic approach.  Later, proper connection to emission
!... inventory would be better, but, for now, scale like this:
!... where: Nwfa=50 per cc, emit 0.875E4 aerosols per second per grid box unit
!... that was tested as ~(20kmx20kmx50m = 2.E10 m**3).

 k = 1
 do iCell = 1, nCellsSolve
    ! airmass = rho_zz(k,iCell)*zz(k,iCell)
    ! airmass = airmass*(zgrid(k+1,iCell)-zgrid(k,iCell))*areaCell(iCell) ! (in kg)
    ! nwfa2d(iCell) = nwfa(k,iCell)*0.000196*airmass*0.5e-10

    deltaz = zgrid(k+1,iCell)-zgrid(k,iCell)
    if (deltaz < 9.0) then
       deltaz = 9.0
    endif
    nwfa2d(iCell) = max(nwfa_default, min(aero_max, nwfa(k,iCell))) * &
         0.000196 * 50.0 / deltaz
    nifa2d(iCell) = 0._RKIND
!   call mpas_log_write('$i $r $r $r',intArgs=(/iCell/),realArgs=(/airmass,nwfa2d(iCell),nifa2d(iCell)/))
 enddo

!... deallocate local arrays:
 if(allocated(hgt)) deallocate(hgt)

!call mpas_log_write('--- end subroutine init_tempo_aerosols_forMPAS.')

 end subroutine init_tempo_aerosols_forMPAS

!=================================================================================================================
 subroutine table_ccnAct(dminfo)
!=================================================================================================================

!input variables:
 type(dm_info),intent(in):: dminfo

!local variables:
 logical:: opened
 integer:: ccn_unit,i,istat
 character(len=StrKIND):: errmess
!-----------------------------------------------------------------------------------------------------------------

 if(.not.allocated(tnccn_act)) allocate(tnccn_act(ntb_arc,ntb_arw,ntb_art,ntb_arr,ntb_ark))

!get a unit to open binary file:
 istat = -999
 if(dminfo % my_proc_id == IO_NODE) then
    do i = 10,99
       inquire(i,opened = opened,iostat=istat)
       if(.not. opened ) then
          ccn_unit = i
          exit
       endif
    enddo
    if(istat /= 0) &
       call physics_error_fatal('mpas_atmphys_init_microphysics table_ccnAct: Can not '// &
                                'find unused fortran unit to read in lookup table.' )
 endif

!distribute unit to other processors:
 call mpas_dmpar_bcast_int(dminfo,ccn_unit)

!open binary file:
 istat = -999
 if(dminfo % my_proc_id == IO_NODE) then
    open(ccn_unit,file='CCN_ACTIVATE_DATA',form='UNFORMATTED',status='OLD',iostat=istat)
    if(istat /= 0) then
       write(errmess,'(A,I4)') 'mpas_atmphys_init_microphysics table_ccnAct:: '// &
                               'error opening CCN_ACTIVATE_DATA on unit', ccn_unit
       call physics_error_fatal(errmess)
    endif
 endif

!read and broadcast data to all nodes:
 istat = -999
 if(dminfo % my_proc_id == IO_NODE) then
    read(ccn_unit,iostat=istat) tnccn_act
    if(istat /= 0) then
       write(errmess,'(A,I4)') 'mpas_atmphys_init_microphysics table_ccnAct:: '// &
                               'error reading tnccn_act on unit', ccn_unit
       call physics_error_fatal(errmess)
    endif
 endif

 DM_BCAST_MACRO(tnccn_act)

 end subroutine table_ccnAct

!=================================================================================================================
 end module mpas_atmphys_init_tempo
!=================================================================================================================
