module mod_io_observables
  use mod_kinds, only: dp
  use mod_types, only: SimParams, SimObs
  implicit none
contains

  subroutine reset_obs(obs)
    type(SimObs), intent(inout) :: obs
    obs = SimObs()  ! azzera tutto (default init)
  end subroutine reset_obs

  subroutine update_obs(par, obs, n_MC, ekin, pot, mag)
    type(SimParams), intent(in)    :: par
    type(SimObs),    intent(inout) :: obs
    integer, intent(in) :: n_MC
    real(dp), intent(in) :: ekin, pot, mag

    real(dp) :: Neff, kin_per_dof
    real(dp) :: a1, a2, a3

    Neff = real(par%N -1, dp)     ! dof cinetici effettivi
    kin_per_dof = ekin / Neff

    ! accumuli
    obs%av_mag  = obs%av_mag  + abs(mag)

    obs%av_pot  = obs%av_pot + pot / real(par%N,dp)
    obs%av_kin  = obs%av_kin + ekin / real(par%N,dp)

    obs%av_kin_1 = obs%av_kin_1 + 1.0_dp/kin_per_dof
    obs%av_kin_2 = obs%av_kin_2 + 1.0_dp/(kin_per_dof*kin_per_dof)
    obs%av_kin_3 = obs%av_kin_3 + 1.0_dp/(kin_per_dof*kin_per_dof*kin_per_dof)

    obs%av_en   = obs%av_en + (ekin + pot) / real(par%N,dp)

    ! medie correnti
    obs%c_mag  = obs%av_mag  / real(n_MC,dp)

    obs%c_pot = obs%av_pot / real(n_MC,dp)
    obs%c_kin = obs%av_kin / real(n_MC,dp)

    obs%c_kin_1 = obs%av_kin_1 / real(n_MC,dp)
    obs%c_kin_2 = obs%av_kin_2 / real(n_MC,dp)
    obs%c_kin_3 = obs%av_kin_3 / real(n_MC,dp)

    obs%c_en = obs%av_en / real(n_MC,dp)

    ! beta e derivate nello stesso stile del tuo vecchio codice:
    ! beta = (1/2 - 1/Neff) < 1/(K/Neff) >
    a1 = (0.5_dp - 1.0_dp/Neff)
    a2 = (0.5_dp - 2.0_dp/Neff)
    a3 = (0.5_dp - 3.0_dp/Neff)

    obs%beta = a1 * obs%c_kin_1

    obs%ders_2 = Neff * ( a1*a2*obs%c_kin_2 - (a1*a1)*(obs%c_kin_1*obs%c_kin_1) )

    obs%ders_3 = Neff*Neff * ( &
         a1*a2*a3*obs%c_kin_3 &
       - 3.0_dp*(a1*a1)*a2*obs%c_kin_2*obs%c_kin_1 &
       + 2.0_dp*(a1*a1*a1)*(obs%c_kin_1*obs%c_kin_1*obs%c_kin_1) )

  end subroutine update_obs

subroutine append_observables_legacy(par, obs, time_tot, n_MC, n_realiz, n_rest)
  use mod_kinds, only: dp
  use mod_types, only: SimParams, SimObs
  implicit none
  type(SimParams), intent(in) :: par
  type(SimObs),    intent(in) :: obs
  real(dp), intent(in) :: time_tot
  integer, intent(in) :: n_MC, n_realiz, n_rest
  character(len=256) :: fname
  integer :: u

  write(fname,'("obs_",I0,"_",I0,"_",I0,"_",I0,".dat")') par%n_samp, n_realiz, n_rest, par%N
  open(newunit=u, file=fname, status="unknown", action="write")

  ! time_tot, n_MC, en_in, c_kin, c_pot, c_kin_1, c_kin_2, c_kin_3,
  ! beta, ders_2, ders_3, c_mag, c_mag2, c_mag4
  write(u,'(ES14.6,1X,I10,1X,ES14.6,1X,ES14.6,1X,ES14.6,1X,ES14.6,1X,ES14.6,1X,ES14.6,1X,ES14.6,1X,&
&  ES14.6,1X,ES14.6,1X,ES14.6,1X,ES14.6,1X,ES14.6,1X,ES14.6,1X,ES14.6)') &
       time_tot, n_MC, par%en_in, obs%c_en, &
       obs%c_kin, obs%c_pot, obs%c_kin_1, obs%c_kin_2, obs%c_kin_3, &
       obs%beta, obs%ders_2, obs%ders_3, obs%c_mag

  close(u)
end subroutine append_observables_legacy


end module mod_io_observables
