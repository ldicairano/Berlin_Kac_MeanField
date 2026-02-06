program main_berlin_kac_rattle
  use mod_kinds, only: dp
  use mod_types, only: SimParams, SimState, SimObs
  use mod_input, only: InputData, read_input_file
  use mod_init_spherical, only: init_microcanonical_mf
  use mod_rattle_spherical, only: rattle_step_mf
  use mod_io_restart, only: save_restart, load_restart
  use mod_io_observables, only: reset_obs, update_obs, append_observables_legacy
  implicit none


  !                        ██████╗ ██╗  ██╗    ███╗   ███╗ ███████╗                         
  !                        ██╔══██╗██║ ██╔╝    ████╗ ████║ ██╔════╝                         
  !                        ██████╔╝█████╔╝     ██╔████╔██║ █████╗                           
  !                        ██╔══██╗██╔═██╗     ██║╚██╔╝██║ ██╔══╝                           
  !                        ██████╔╝██║  ██╗    ██║ ╚═╝ ██║ ██║                              
  !                        ╚═════╝ ╚═╝  ╚═╝    ╚═╝     ╚═╝ ╚═╝                              
  !                                                                                
  !  //============================================================================\\ 
  !  |                           BERLIN-KAC  (BK)                                  |
  !  |                     MEAN-FIELD  SPHERICAL  MODEL                            |
  !  |                     MICROCANONICAL  +  RATTLE                               |
  !  \\===========================================================================// 
                                                                                

  type(InputData) :: inp
  type(SimParams) :: par
  type(SimState)  :: st
  type(SimObs)    :: obs

  character(len=256) :: input_file, arg
  integer :: nargs
  integer :: n, n_samp, n_realiz, n_rest
  integer :: n_steps, n_jump
  integer :: n_MC, n_camp, n_initial
  real(dp) :: cluster_time
  real(dp) :: time_tot
  real(dp) :: ekin, pot, mag, Etot_now
  real(dp) :: t_in, t_out, t_in_step, t_out_step

  ! ---------------------------
  ! CLI parsing
  ! Usage:
  !   ./BKMF_RATTLE -i input.inp n n_samp n_realiz n_rest
  ! or
  !   ./BKMF_RATTLE n_samp n_realiz n_rest    (uses input.inp)
  ! ---------------------------
  input_file = "input.inp"
  nargs = command_argument_count()

  if (nargs < 3) then
    write(*,*) "Usage:"
    write(*,*) "  ./BKMF_RATTLE -i input.inp n n_samp n_realiz n_rest"
    stop
  end if

  ! optional -i <file> in first two args
  call get_command_argument(1, arg)
  if (trim(arg) == "-i") then
    if (nargs < 6) then
      write(*,*) "ERROR: with -i you must provide: -i input.inp n_samp n_realiz n_rest"
      stop
    end if
    call get_command_argument(2, input_file)
  end if

  ! last 3 args are always n_samp, n_realiz, n_rest
  call get_command_argument(nargs-3, arg); read(arg,*) n
  call get_command_argument(nargs-2, arg); read(arg,*) n_samp
  call get_command_argument(nargs-1, arg); read(arg,*) n_realiz
  call get_command_argument(nargs,   arg); read(arg,*) n_rest

  ! ---------------------------
  ! Read input.inp (J,n_steps,n_jump,cluster_time,dt)
  ! ---------------------------
  call read_input_file(trim(input_file), inp)

  cluster_time = inp%cluster_time
  n_steps      = inp%n_steps
  n_jump       = inp%n_jump

  ! map into SimParams
  par%N  = N
  par%J  = inp%J
  par%dt = inp%dt

  par%n_samp = n_samp
  par%en_in  = real(n_samp,dp)*0.01_dp
  par%Etot   = real(par%N,dp)*par%en_in

  allocate(st%q(par%N), st%p(par%N))
  st%Ssum = 0.0_dp

  ! ---------------------------
  ! Init / Restart
  ! Semantics: n_rest labels output segment.
  ! If n_rest>1, load from previous segment (n_rest-1) and continue.
  ! ---------------------------
  time_tot = 0.0_dp
  call cpu_time(t_in)

  if (n_rest > 1) then
    call load_restart(par, st, obs, n_initial, time_tot, n_realiz, par%n_samp, n_rest-1)
    n_initial = n_initial + 1
  else
    call reset_obs(obs)
    call init_microcanonical_mf(par, st, ekin, pot, mag)
    n_initial = 1
    time_tot  = 0.0_dp
  end if

  call cpu_time(t_out)
  time_tot = time_tot + abs(t_out - t_in)

  ! ---------------------------
  ! Main loop
  ! ---------------------------
  do n_MC = n_initial, n_steps

    call cpu_time(t_in_step)
    do n_camp = 1, n_jump
      call rattle_step_mf(par, st, ekin, pot, mag, Etot_now)
    end do
    call cpu_time(t_out_step)
    time_tot = time_tot + abs(t_out_step - t_in_step)

    if (mod(n_MC,1000) == 0) then
      call append_observables_legacy(par, obs, time_tot, n_MC, n_realiz, n_rest)
      call save_restart(par, st, obs, n_MC, time_tot, n_realiz, n_samp, n_rest)
    end if

    if (time_tot >= 0.95_dp*cluster_time) exit

    call update_obs(par, obs, n_MC, ekin, pot, mag)
  end do

  ! final dump
  call append_observables_legacy(par, obs, time_tot, n_MC-1, n_realiz, n_rest)
  call save_restart(par, st, obs, n_MC, time_tot, n_realiz, n_samp, n_rest)

  deallocate(st%q, st%p)

end program main_berlin_kac_rattle
