module mod_input
  use mod_kinds, only: dp
  implicit none

  type :: InputData
    real(dp) :: cluster_time = 172800.0_dp
    real(dp) :: J = 1.0_dp
    integer  :: n_steps = 0
    integer  :: n_jump  = 0
    real(dp) :: dt = 0.01_dp
  end type InputData

contains

  subroutine read_input_file(filename, inp)
    character(len=*), intent(in) :: filename
    type(InputData), intent(inout) :: inp
    integer :: u, ios

    ! local vars for NAMELIST
    real(dp) :: cluster_time, J, dt
    integer  :: n_steps, n_jump

    namelist /system/ cluster_time
    namelist /model/  J
    namelist /run/    n_steps, n_jump, dt

    ! defaults
    cluster_time = inp%cluster_time
    J = inp%J
    n_steps = inp%n_steps
    n_jump  = inp%n_jump
    dt = inp%dt

    open(newunit=u, file=filename, status="old", action="read", iostat=ios)
    if (ios /= 0) then
      write(*,*) "ERROR: cannot open input file: ", trim(filename)
      stop
    end if

    rewind(u); read(u, nml=system, iostat=ios)
    rewind(u); read(u, nml=model,  iostat=ios)
    rewind(u); read(u, nml=run,    iostat=ios)
    close(u)

    inp%cluster_time = cluster_time
    inp%J = J
    inp%n_steps = n_steps
    inp%n_jump  = n_jump
    inp%dt = dt

    ! sanity
    if (inp%n_steps <= 0) then
      write(*,*) "ERROR: n_steps must be set in &run"
      stop
    end if
    if (inp%n_jump <= 0) then
      write(*,*) "ERROR: n_jump must be set in &run"
      stop
    end if
  end subroutine read_input_file

end module mod_input
