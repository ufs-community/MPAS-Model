!=================================================================================================================
 module bl_mynn_post
 use ccpp_kind_types,only: kind_phys

 implicit none
 private
 public:: bl_mynn_post_init,     &
          bl_mynn_post_finalize, &
          bl_mynn_post_run


 contains


!=================================================================================================================
!>\section arg_table_bl_mynn_post_init
!!\html\include bl_mynn_post_init.html
!!
 subroutine bl_mynn_post_init(errmsg,errflg)
!=================================================================================================================

!--- output arguments:
 character(len=*),intent(out):: &
    errmsg      ! output error message (-).

 integer,intent(out):: &
    errflg      ! output error flag (-).

!-----------------------------------------------------------------------------------------------------------------

!--- output error flag and message:
 errflg = 0
 errmsg = " "

 end subroutine bl_mynn_post_init

!=================================================================================================================
!>\section arg_table_bl_mynn_post_finalize
!!\html\include bl_mynn_post_finalize.html
!!
 subroutine bl_mynn_post_finalize(errmsg,errflg)
!=================================================================================================================

!--- output arguments:
 character(len=*),intent(out):: &
    errmsg      ! output error message (-).

 integer,intent(out):: &
    errflg      ! output error flag (-).

!-----------------------------------------------------------------------------------------------------------------

!--- output error flag and message:
 errflg = 0
 errmsg = " "

 end subroutine bl_mynn_post_finalize

!=================================================================================================================
 subroutine bl_mynn_post_run(its,ite,kte,f_qc,f_qi,f_qs,delt,qv,qc,qi,qs,dqv,dqc,dqi,dqs,errmsg,errflg)
!=================================================================================================================

!--- input arguments:
 logical,intent(in):: &
    f_qc, &! if true,the physics package includes the cloud liquid water mixing ratio.
    f_qi, &! if true,the physics package includes the cloud ice mixing ratio.
    f_qs   ! if true,the physics package includes the snow mixing ratio.

 integer,intent(in):: its,ite
 integer,intent(in):: kte

 real(kind=kind_phys),intent(in):: &
    delt   !

 real(kind=kind_phys),intent(in),dimension(1:kte):: &
    qv,   &!
    qc,   &!
    qi,   &!
    qs     !


!--- inout arguments:
 real(kind=kind_phys),intent(inout),dimension(1:kte):: &
    dqv,  &!
    dqc,  &!
    dqi,  &!
    dqs    !


!--- output arguments:
 character(len=*),intent(out):: errmsg
 integer,intent(out):: errflg


!--- local variables:
 integer:: k,kts
 real(kind=kind_phys):: rq,sq,tem
 real(kind=kind_phys),dimension(1:kte):: sqv,sqc,sqi,sqs

!-----------------------------------------------------------------------------------------------------------------

!--- initialization:
 kts = 1

!---
 do k = kts,kte
    sq = qv(k)/(1.+qv(k))      !conversion of qv at time-step n from mixing ratio to specific humidity.
    sqv(k) = sq + dqv(k)*delt  !calculation of specific humidity at time-step n+1.
    rq = sqv(k)/(1.-sqv(k))    !conversion of qv at time-step n+1 from specific humidity to mixing ratio.
    dqv(k) = (rq - qv(k))/delt !calculation of the tendency.
 enddo

 if (f_qc) then
    do k = kts,kte
       sq = qc(k)/(1.+qv(k))
       sqc(k) = sq + dqc(k)*delt
       rq  = sqc(k)*(1.+sqv(k))
       dqc(k) = (rq - qc(k))/delt
    enddo
 endif

 if (f_qi) then
    do k = kts,kte
       sq = qi(k)/(1.+qv(k))
       sqi(k) = sq + dqi(k)*delt
       rq = sqi(k)*(1.+sqv(k))
       dqi(k) = (rq - qi(k))/delt
    enddo
 endif

 if (f_qs) then
    do k = kts,kte
       sq = qs(k)/(1.+qv(k))
       sqs(k) = sq + dqs(k)*delt
       rq = sqs(k)*(1.+sqv(k))
       dqs(k) = (rq - qs(k))/delt
    enddo
 endif

!--- output error flag and message:
 errmsg = " "
 errflg = 0

 end subroutine bl_mynn_post_run

!=================================================================================================================
 end module bl_mynn_post
!=================================================================================================================
