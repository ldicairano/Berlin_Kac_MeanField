module mod_init_spherical
  use mod_kinds,  only: dp
  use mod_types,  only: SimParams, SimState
  use mod_random, only: gauss01
  use mod_mf_spherical_model, only: mf_force_potential
  implicit none
contains

  subroutine init_microcanonical_mf(par, st, ekin, pot, mag)
    type(SimParams), intent(in)    :: par
    type(SimState),  intent(inout) :: st
    real(dp), intent(out) :: ekin, pot, mag

    integer :: i, N, it
    real(dp) :: Nreal, invN, sqrtN
    real(dp) :: eps, m0, m_try, Ktarget, K0, scale
    real(dp) :: x, qp, c
    real(dp) :: Fdummy
    real(dp) :: u(par%N), v(par%N)

    N     = par%N
    Nreal = real(N,dp)
    invN  = 1.0_dp / Nreal
    sqrtN = sqrt(Nreal)
    eps   = par%en_in

    ! ------------------------------------------------------------
    ! Choose a robust starting magnetization guess:
    ! For MF spherical (J>0): below eps_c=J/2 equilibrium has m^2=(J/2 - eps)/J.
    ! Above eps_c: m ~ 0.
    ! ------------------------------------------------------------
    if (eps < 0.5_dp*par%J) then
      m0 = sqrt( max(0.0_dp, (0.5_dp*par%J - eps) / par%J ) )
    else
      m0 = 0.0_dp
    end if
    m0 = min(m0, 0.999999_dp)

    ! ------------------------------------------------------------
    ! Build u (uniform direction) and a random v orthogonal to u
    ! ------------------------------------------------------------
    do i=1,N
      u(i) = 1.0_dp/sqrtN
    end do

    ! If m0=0, we can just sample a random point on the sphere as before
    ! but keep it in the same framework for consistency.
    m_try = m0

    do it = 1, 50
      call random_orthonormal_to_u(N, u, v)  ! v unit, v·u=0

      ! q = sqrtN( m_try u + sqrt(1-m_try^2) v )
      do i=1,N
        st%q(i) = sqrtN * ( m_try*u(i) + sqrt(max(0.0_dp,1.0_dp-m_try*m_try))*v(i) )
      end do

      st%Ssum = 0.0_dp
      do i=1,N
        st%Ssum = st%Ssum + st%q(i)
      end do

      call mf_force_potential(par, st%Ssum, Fdummy, pot, mag)

      Ktarget = par%Etot - pot

      if (Ktarget >= 0.0_dp) exit

      ! If Ktarget < 0, increase |m| (makes V more negative => increases Ktarget)
      m_try = min(0.999999_dp, max(m_try + 0.05_dp, 0.05_dp))
      if (it == 50) then
        write(*,*) "ERROR init: Etot - V(q) < 0 even after increasing m."
        write(*,*) "Etot=",par%Etot," pot=",pot," eps=",eps," J=",par%J
        stop
      end if
    end do

    ! ------------------------------------------------------------
    ! Initialize p ~ Gaussian, project to tangent space (q·p=0),
    ! then scale to match Ktarget.
    ! ------------------------------------------------------------
    do i=1,N
      call gauss01(x)
      st%p(i) = x
    end do

    qp = 0.0_dp
    do i=1,N
      qp = qp + st%q(i)*st%p(i)
    end do
    c = qp / ( Nreal )   ! since q·q = N
    do i=1,N
      st%p(i) = st%p(i) - c*st%q(i)
    end do

    K0 = 0.0_dp
    do i=1,N
      K0 = K0 + 0.5_dp*st%p(i)*st%p(i)
    end do
    if (K0 <= 0.0_dp) then
      write(*,*) "ERROR init: K0 <= 0 after projection."
      stop
    end if

    scale = sqrt( Ktarget / K0 )
    do i=1,N
      st%p(i) = scale*st%p(i)
    end do

    ekin = 0.0_dp
    do i=1,N
      ekin = ekin + 0.5_dp*st%p(i)*st%p(i)
    end do

  contains

  subroutine random_orthonormal_to_u(N, u, v)
    use mod_kinds,  only: dp
    use mod_random, only: gauss01
    implicit none
    integer, intent(in) :: N
    real(dp), intent(in)  :: u(N)
    real(dp), intent(out) :: v(N)
  
    integer :: i, it
    real(dp) :: r(N), dotru, normv, x
  
    do it = 1, 50
      ! random r
      do i=1,N
        call gauss01(x)
        r(i) = x
      end do
  
      ! v = r - (r·u) u
      dotru = 0.0_dp
      do i=1,N
        dotru = dotru + r(i)*u(i)
      end do
      do i=1,N
        v(i) = r(i) - dotru*u(i)
      end do
  
      ! normalize v
      normv = 0.0_dp
      do i=1,N
        normv = normv + v(i)*v(i)
      end do
      normv = sqrt(normv)
  
      if (normv > 1.0d-300) then
        do i=1,N
          v(i) = v(i)/normv
        end do
        return
      end if
    end do
  
    write(*,*) "ERROR: failed to generate vector orthonormal to u (numerical issue)."
    stop
  end subroutine random_orthonormal_to_u
  

  end subroutine init_microcanonical_mf

end module mod_init_spherical
