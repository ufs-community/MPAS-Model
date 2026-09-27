module mpas_gfs_pressure_repair
 use mpas_kind_types
 use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
 implicit none
 private
 public :: repair_gfs_pressure_columns

contains

 subroutine repair_gfs_pressure_columns(p_mid_native,p_int_native,       &
                                        p_mid_gfs,p_int_gfs,             &
                                        nEdgesOnCell,cellsOnCell,        &
                                        its,ite,kts,kte)
   integer,intent(in) :: its,ite,kts,kte
   integer,intent(in) :: nEdgesOnCell(:)
   integer,intent(in) :: cellsOnCell(:,:)
   real(kind=RKIND),intent(in)  :: p_mid_native(its:ite,kts:kte)
   real(kind=RKIND),intent(in)  :: p_int_native(its:ite,kts:kte+1)
   real(kind=RKIND),intent(out) :: p_mid_gfs(its:ite,kts:kte)
   real(kind=RKIND),intent(out) :: p_int_gfs(its:ite,kts:kte+1)

   logical,allocatable :: press_valid(:)
   integer :: i,k,n,nbr,ngood
   real(kind=RKIND) :: min_dp_press,p_lower,p_upper

   min_dp_press=1.0_RKIND
   allocate(press_valid(its:ite))

   p_mid_gfs(:,:) = p_mid_native(:,:)
   p_int_gfs(:,:) = p_int_native(:,:)
   press_valid(:)=.true.

   ! Same validity test used by the previous common-interface repair.
   do i=its,ite
      do k=kts,kte
         if (.not.ieee_is_finite(p_mid_native(i,k)) .or. &
             .not.ieee_is_finite(p_int_native(i,k)) .or. &
             .not.ieee_is_finite(p_int_native(i,k+1))) then
            press_valid(i)=.false.
            exit
         endif
         if (p_int_native(i,k) <= 0.0_RKIND .or. &
             p_int_native(i,k+1) <= 0.0_RKIND .or. &
             p_int_native(i,k+1) >= p_int_native(i,k) .or. &
             p_mid_native(i,k) >= p_int_native(i,k) .or. &
             p_mid_native(i,k) <= p_int_native(i,k+1)) then
            press_valid(i)=.false.
            exit
         endif
      enddo
   enddo

   do i=its,ite
      if (press_valid(i)) cycle

      p_mid_gfs(i,kts:kte)=0.0_RKIND
      p_int_gfs(i,kts:kte+1)=0.0_RKIND
      ngood=0

      ! Preserve previous behavior: use only valid owned neighbors.
      do n=1,nEdgesOnCell(i)
         nbr=cellsOnCell(n,i)
         if (nbr < its .or. nbr > ite) cycle
         if (.not.press_valid(nbr)) cycle
         ngood=ngood+1
         p_mid_gfs(i,kts:kte)=p_mid_gfs(i,kts:kte)+p_mid_native(nbr,kts:kte)
         p_int_gfs(i,kts:kte+1)=p_int_gfs(i,kts:kte+1)+p_int_native(nbr,kts:kte+1)
      enddo

      if (ngood > 0) then
         p_mid_gfs(i,kts:kte)=p_mid_gfs(i,kts:kte)/real(ngood,RKIND)
         p_int_gfs(i,kts:kte+1)=p_int_gfs(i,kts:kte+1)/real(ngood,RKIND)
      else
         if (ieee_is_finite(p_int_native(i,kts)) .and. &
             p_int_native(i,kts) > 0.0_RKIND) then
            p_int_gfs(i,kts)=p_int_native(i,kts)
         else
            p_int_gfs(i,kts)=100000.0_RKIND
         endif

         do k=kts,kte
            p_lower=p_int_gfs(i,k)
            if (ieee_is_finite(p_int_native(i,k+1)) .and. &
                p_int_native(i,k+1) > 0.0_RKIND .and. &
                p_int_native(i,k+1) < p_lower) then
               p_upper=p_int_native(i,k+1)
            else if (p_lower > 2.0_RKIND*min_dp_press) then
               p_upper=p_lower-min_dp_press
            else
               p_upper=0.5_RKIND*p_lower
            endif

            if (p_lower > 2.0_RKIND*min_dp_press) then
               if (p_lower-p_upper < min_dp_press) p_upper=p_lower-min_dp_press
            else if (p_upper >= p_lower .or. p_upper <= 0.0_RKIND) then
               p_upper=0.5_RKIND*p_lower
            endif

            p_int_gfs(i,k+1)=p_upper
            p_mid_gfs(i,k)=0.5_RKIND*(p_lower+p_upper)
         enddo
      endif

      ! Same final guarantee pass as the previous common repair.
      if (.not.ieee_is_finite(p_int_gfs(i,kts)) .or. &
          p_int_gfs(i,kts) <= 0.0_RKIND) p_int_gfs(i,kts)=100000.0_RKIND

      do k=kts,kte
         p_lower=p_int_gfs(i,k)
         p_upper=p_int_gfs(i,k+1)

         if (.not.ieee_is_finite(p_lower) .or. p_lower <= 0.0_RKIND) then
            if (k == kts) then
               p_lower=100000.0_RKIND
            else
               p_lower=max(2.0_RKIND*tiny(1.0_RKIND), &
                           0.5_RKIND*p_int_gfs(i,k-1))
            endif
            p_int_gfs(i,k)=p_lower
         endif

         if (.not.ieee_is_finite(p_upper) .or. p_upper <= 0.0_RKIND .or. &
             p_upper >= p_lower) then
            if (p_lower > 2.0_RKIND*min_dp_press) then
               p_upper=p_lower-min_dp_press
            else
               p_upper=0.5_RKIND*p_lower
            endif
         else if (p_lower > 2.0_RKIND*min_dp_press .and. &
                  p_lower-p_upper < min_dp_press) then
            p_upper=p_lower-min_dp_press
         endif
         p_int_gfs(i,k+1)=p_upper

         if (.not.ieee_is_finite(p_mid_gfs(i,k)) .or. &
             p_mid_gfs(i,k) >= p_lower .or. &
             p_mid_gfs(i,k) <= p_upper) then
            p_mid_gfs(i,k)=0.5_RKIND*(p_lower+p_upper)
         endif
      enddo
   enddo

   deallocate(press_valid)
 end subroutine repair_gfs_pressure_columns
end module mpas_gfs_pressure_repair
