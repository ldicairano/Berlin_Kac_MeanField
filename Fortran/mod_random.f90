module mod_random
  use mod_kinds, only: dp
  implicit none
contains
  subroutine gauss01(x)
    real(dp), intent(out) :: x
    real(dp) :: u1, u2, twopi
    twopi = 2.0_dp*acos(-1.0_dp)
    call random_number(u1)
    call random_number(u2)
    if (u1 < 1.0e-16_dp) u1 = 1.0e-16_dp
    x = sqrt(-2.0_dp*log(u1)) * cos(twopi*u2)
  end subroutine gauss01
end module mod_random
