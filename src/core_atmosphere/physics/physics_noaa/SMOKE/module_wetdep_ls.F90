!>\file module_wetdep_ls.F90
!! This file contains aerosol wet deposition module.

module module_wetdep_ls
  use mpas_kind_types

contains

subroutine wetdep_ls(dt,var,rain,moist,                                         &
                     rho,nchem,dz8w,vvel,p_phy,alpha,                           &
                     wetdep_resolved,                                           &
                     ids,ide, jds,jde, kds,kde,                                 &
                     ims,ime, jms,jme, kms,kme,                                 &
                     its,ite, jts,jte, kts,kte                                  )
   implicit none

   integer,      intent(in) :: nchem,                                           &
                               ids,ide, jds,jde, kds,kde,                       &
                               ims,ime, jms,jme, kms,kme,                       &
                               its,ite, jts,jte, kts,kte
   real(RKIND), intent(in) :: dt
   real(RKIND), intent(in) :: alpha
   real(RKIND), dimension( ims:ime, kms:kme, jms:jme),intent(in) :: moist, p_phy
   real(RKIND), dimension( ims:ime, kms:kme, jms:jme),intent(in) :: rho,dz8w,vvel        
   real(RKIND), dimension( ims:ime, jms:jme),         intent(in) :: rain
   real(RKIND), dimension( ims:ime, kms:kme, jms:jme,1:nchem),intent(inout) :: var        
   real(RKIND), dimension( ims:ime, jms:jme, 1:nchem ),       intent(out)   ::  wetdep_resolved

! Local Variables
   real(RKIND), dimension( its:ite, jts:jte) :: var_sum,var_rmv
   real(RKIND), dimension( its:ite, kts:kte, jts:jte) :: var_rmvl
   real(RKIND), dimension( its:ite, jts:jte) :: frc,var_sum_clw,rain_clw     
   real(RKIND) :: dvar,factor,clsum,factor_diag
   integer :: nv,i,j,k,km,kb,kbeg

    do nv=1,nchem
      do i=its,ite
       do j=jts,jte
        var_sum_clw(i,j)=0.
        var_sum(i,j)=0.
        var_rmvl(i,:,j)=0.
        frc(i,j)=0.
        rain_clw(i,j)=0.
        if(rain(i,j).gt.1.e-10)then
! convert rain back to rate
!
           rain_clw(i,j)=rain(i,j)/dt
! total cloud water
!
           do k=1,kte-1
              dvar=max(0.,moist(i,k,j)*rho(i,k,j)*vvel(i,k,j)*dz8w(i,k,j))
              var_sum_clw(i,j)=var_sum_clw(i,j)+dvar
              var_sum(i,j)=var_sum(i,j)+var(i,k,j,nv)*rho(i,k,j)
           enddo
           if(var_sum(i,j).gt.1.e-10 .and. var_sum_clw(i,j).gt.1.e-10 ) then
!          assuming that frc is onstant, it is my conversion factor 
!          (just like in convec. parameterization)
              frc(i,j)=rain_clw(i,j)/var_sum_clw(i,j)
              frc(i,j)=max(1.e-6,min(frc(i,j),.005))
           endif
        endif
       enddo
       enddo
!
! get rid of it
!
       do i=its,ite
       do j=jts,jte
       if(rain(i,j).gt.1.e-10 .and. var_sum(i,j).gt.1.e-10 .and. var_sum_clw(i,j).gt.1.e-10)then
         do k=kts,kte-2
          if(var(i,k,j,nv).gt.1.e-16 .and. moist(i,k,j).gt.0.)then
            factor = max(0.,frc(i,j)*rho(i,k,j)*dz8w(i,k,j)*vvel(i,k,j))
            dvar=alpha*factor/(1+factor)*var(i,k,j,nv)
            factor_diag = ( dvar * (p_phy(i,k+1,j) - p_phy(i,k,j))*100. / 9.8 / dt )
! Accumulate diags
            wetdep_resolved(i,j,nv) = wetdep_resolved(i,j,nv) + factor_diag
            var(i,k,j,nv)=max(1.e-16,var(i,k,j,nv)-dvar)
          endif
         enddo
       endif
       enddo
       enddo
      enddo ! nv
end subroutine wetdep_ls 
end module module_wetdep_ls
