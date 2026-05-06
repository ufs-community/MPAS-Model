      module canopy_utils_mod

      use gfs_tke_edmf_machine, only : kind_phys
      implicit none

      contains

      function IntegrateTrapezoid(x, y)
      real(kind=kind_phys), intent(in) :: x(:)
      real(kind=kind_phys), intent(in) :: y(size(x))
      real(kind=kind_phys)             :: IntegrateTrapezoid
      integer                          :: n

      n = size(x)

      IntegrateTrapezoid = sum((y(2:n) + y(1:n-1)) *
     &                   (x(2:n) - x(1:n-1))) / 2.0_kind_phys

      end function IntegrateTrapezoid

      end module canopy_utils_mod
