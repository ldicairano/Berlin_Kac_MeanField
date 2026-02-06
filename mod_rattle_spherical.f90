module mod_rattle_spherical
  use mod_kinds, only: dp
  use mod_types, only: SimParams, SimState
  use mod_mf_spherical_model, only: mf_force_potential
  implicit none
contains


subroutine rattle_step_mf(par, st, ekin, pot, mag, Etot_now)
  use mod_kinds, only: dp
  use mod_types, only: SimParams, SimState
  use mod_mf_spherical_model, only: mf_force_potential
  implicit none
  type(SimParams), intent(in)    :: par
  type(SimState),  intent(inout) :: st
  real(dp), intent(out) :: ekin, pot, mag, Etot_now

  integer :: i, N
  real(dp) :: dt, Nreal, invN
  real(dp) :: F1, F2
  real(dp) :: S0, S1
  real(dp) :: Psum, P2
  real(dp) :: A, C, D
  real(dp) :: alpha, beta, gamma, disc, sdisc
  real(dp) :: lam_n, lam_np1, lam_lim, r1, r2
  real(dp) :: qp_half, qF2
  real(dp) :: p2_final

  N     = par%N
  Nreal = real(N,dp)
  invN  = 1.0_dp / Nreal
  dt    = par%dt

  associate(q => st%q, p => st%p)

    ! ---------- Force at time n (uniform scalar) ----------
    S0 = st%Ssum
    call mf_force_potential(par, S0, F1, pot, mag)

    ! ---------- LOOP 1: compute Psum and P2 at time n ----------
    Psum = 0.0_dp
    P2   = 0.0_dp
    do i = 1, N
      Psum = Psum + p(i)
      P2   = P2   + p(i)*p(i)
    end do

    ! Scalars for quadratic enforcing |q_{n+1}|^2 = N
    ! A = q·F = F1 * sum(q) = F1 * S0
    A = F1 * S0
    ! C = p·F = F1 * sum(p) = F1 * Psum
    C = F1 * Psum
    ! D = F·F = N * F1^2
    D = Nreal * F1*F1

    alpha = (dt*dt/4.0_dp) * Nreal
    beta  = -Nreal - (dt*dt/2.0_dp) * A
    gamma = A + P2 + dt*C + (dt*dt/4.0_dp)*D

    disc = beta*beta - 4.0_dp*alpha*gamma
    if (disc < 0.0_dp) disc = 0.0_dp
    sdisc = sqrt(disc)

    r1 = (-beta + sdisc) / (2.0_dp*alpha)
    r2 = (-beta - sdisc) / (2.0_dp*alpha)

    ! Root selection continuous with dt->0: lambda -> (A+P2)/N
    lam_lim = (A + P2) * invN
    if (abs(r1 - lam_lim) < abs(r2 - lam_lim)) then
      lam_n = r1
    else
      lam_n = r2
    end if

    ! ---------- LOOP 2: half-kick + drift, accumulate S1 and qp_half ----------
    S1      = 0.0_dp
    qp_half = 0.0_dp
    do i = 1, N
      p(i) = p(i) + 0.5_dp*dt*( F1 - lam_n*q(i) )
      q(i) = q(i) + dt*p(i)
      S1      = S1 + q(i)
      qp_half = qp_half + q(i)*p(i)
    end do
    st%Ssum = S1

    ! ---------- Force at time n+1 ----------
    call mf_force_potential(par, S1, F2, pot, mag)

    ! lambda_{n+1} from velocity constraint q_{n+1}·p_{n+1}=0
    qF2    = F2 * S1
    lam_np1 = ( qF2 + (2.0_dp/dt)*qp_half ) * invN

    ! ---------- LOOP 3: second half-kick and compute kinetic energy ----------
    p2_final = 0.0_dp
    do i = 1, N
      p(i) = p(i) + 0.5_dp*dt*( F2 - lam_np1*q(i) )
      p2_final = p2_final + p(i)*p(i)
    end do

    ekin = 0.5_dp * p2_final
    Etot_now = ekin + pot

  end associate
end subroutine rattle_step_mf




end module mod_rattle_spherical
