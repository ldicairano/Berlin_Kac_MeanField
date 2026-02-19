module mod_io_restart
  use mod_kinds, only: dp
  use mod_types, only: SimParams, SimState, SimObs
  implicit none
contains

  subroutine save_restart(par, st, obs, n_MC, time_tot, n_realiz, n_samp, n_rest)
    type(SimParams), intent(in) :: par
    type(SimState),  intent(in) :: st
    type(SimObs),    intent(in) :: obs
    integer, intent(in) :: n_MC, n_realiz, n_samp, n_rest
    real(dp), intent(in) :: time_tot

    character(len=256) :: fname
    integer :: u

    write(fname,'("restart_",I0,"_",I0,"_",I0,"_",I0,".bin")') n_samp, n_realiz, n_rest, par%N
    fname = trim(adjustl(fname))

    open(newunit=u, file=fname, form="unformatted", status="replace", action="write")

    ! --- params (for sanity check on reload)
    write(u) par%N, par%J, par%dt, par%en_in, par%Etot

    ! --- progress
    write(u) n_MC, time_tot

    ! --- state
    write(u) st%Ssum
    write(u) st%q
    write(u) st%p

    ! --- observables (store running sums + derived/current fields)
    write(u) obs%av_en, obs%av_kin, obs%av_pot
    write(u) obs%av_kin_1, obs%av_kin_2, obs%av_kin_3
    write(u) obs%av_mag

    write(u) obs%c_en, obs%c_kin, obs%c_pot
    write(u) obs%c_kin_1, obs%c_kin_2, obs%c_kin_3
    write(u) obs%c_mag

    write(u) obs%beta, obs%ders_2, obs%ders_3

    close(u)
  end subroutine save_restart


  subroutine load_restart(par, st, obs, n_MC, time_tot, n_realiz, n_samp, n_rest)
    type(SimParams), intent(inout) :: par
    type(SimState),  intent(inout) :: st
    type(SimObs),    intent(inout) :: obs
    integer, intent(out) :: n_MC
    real(dp), intent(out) :: time_tot
    integer, intent(in) :: n_realiz, n_samp, n_rest

    character(len=256) :: fname
    integer :: u, Nfile
    real(dp) :: Jfile, dtfile, enfile, Etfile

    write(fname,'("restart_",I0,"_",I0,"_",I0,"_",I0,".bin")') n_samp, n_realiz, n_rest, par%N
    fname = trim(adjustl(fname))

    open(newunit=u, file=fname, form="unformatted", status="old", action="read")

    read(u) Nfile, Jfile, dtfile, enfile, Etfile

    ! check/restore params (safe)
    par%N     = Nfile
    par%J     = Jfile
    par%dt    = dtfile
    par%en_in = enfile
    par%Etot  = Etfile

    read(u) n_MC, time_tot

    read(u) st%Ssum
    read(u) st%q
    read(u) st%p

    read(u) obs%av_en, obs%av_kin, obs%av_pot
    read(u) obs%av_kin_1, obs%av_kin_2, obs%av_kin_3
    read(u) obs%av_mag

    read(u) obs%c_en, obs%c_kin, obs%c_pot
    read(u) obs%c_kin_1, obs%c_kin_2, obs%c_kin_3
    read(u) obs%c_mag

    read(u) obs%beta, obs%ders_2, obs%ders_3

    close(u)
  end subroutine load_restart

end module mod_io_restart
