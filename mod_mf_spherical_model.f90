module mod_mf_spherical_model
  use mod_kinds, only: dp
  use mod_types, only: SimParams
  implicit none
contains

  pure subroutine mf_force_potential(par, Ssum, Fconst, pot, mag)
    type(SimParams), intent(in) :: par
    real(dp), intent(in)  :: Ssum
    real(dp), intent(out) :: Fconst, pot, mag
    real(dp) :: m

    m = Ssum / real(par%N,dp)
    mag = m
    Fconst = par%J * m
    pot = -0.5_dp * par%J * real(par%N,dp) * m*m
  end subroutine mf_force_potential

end module mod_mf_spherical_model
