module mod_types
  use mod_kinds, only: dp
  implicit none

  type :: SimParams
    integer  :: N = 0
    integer  :: n_samp = 0
    real(dp) :: en_in = 0.0_dp
    real(dp) :: Etot  = 0.0_dp
    real(dp) :: J     = 1.0_dp
    real(dp) :: dt    = 0.01_dp
  end type SimParams

  type :: SimState
    real(dp), allocatable :: q(:), p(:)
    real(dp) :: Ssum = 0.0_dp   ! sum_i q_i (MF cache)
  end type SimState

  type :: SimObs
    ! running sums
    real(dp) :: av_en=0.0_dp
    real(dp) :: av_kin=0.0_dp, av_pot=0.0_dp
    real(dp) :: av_kin_1=0.0_dp, av_kin_2=0.0_dp, av_kin_3=0.0_dp
    real(dp) :: av_mag=0.0_dp

    ! current averages
    real(dp) :: c_en=0.0_dp
    real(dp) :: c_kin=0.0_dp, c_pot=0.0_dp
    real(dp) :: c_kin_1=0.0_dp, c_kin_2=0.0_dp, c_kin_3=0.0_dp
    real(dp) :: c_mag=0.0_dp

    ! derived thermodynamic quantities
    real(dp) :: beta=0.0_dp, ders_2=0.0_dp, ders_3=0.0_dp
  end type SimObs

contains
  pure real(dp) function pi_dp()
    pi_dp = 3.1415926535897932384626433832795_dp
  end function pi_dp

end module mod_types
