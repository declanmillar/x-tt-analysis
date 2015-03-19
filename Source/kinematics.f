module kinematics

  use configuration, only: n_fb_asymmetries

  implicit none

  real :: m3,m4,m5,m6,m7,m8
  real :: s
  real :: sigma
  real, parameter :: pi = 3.14159265358979323846d0
  real :: unit_conv ! GeV^{-2} to nb
  parameter (unit_conv=0.38937966d9)
  real :: Xsec_polar(20,-1:1,-1:1),error_polar(20,-1:1,-1:1)
  real :: Xsec_FB(n_fb_asymmetries,20,-1:1),error_FB(n_fb_asymmetries,20,-1:1)
  integer :: npoints

end module kinematics

function longitudinal_neutrino_momentum(p_l,pT_nu)

  use mathematics, only: solve_quadratic
  use quantum_field_theory, only: rm_W

  implicit none

  ! finds the longitudinal neutrino momentum for semi-hadronic decay
  ! assuming all particles are massless

  real :: longitudinal_neutrino_momentum
  real :: p_l(0:3), pT_nu(1:2), p_l0
  real :: pz_nu
  complex, dimension(2) :: pz_nu_solutions
  real :: a,b,c,k

!   ! recalculate lepton energy in zero mass approximation
!   p_l0 = sqrt(p_l(1)*p_l(1) + p_l(2)*p_l(2) + p_l(3)*p_l(3))

!   ! check this matches the 4-vector energy
!   if ( abs(p_l0 - p_l(0)) > 1.d-12) print *, "p_l0 doesn't match"

  k = rm_W*rm_W/2 + p_l(1)*pT_nu(1) + p_l(2)*pT_nu(2)

  a = p_l(1)*p_l(1) + p_l(2)*p_l(2)

  b = -2*k*p_l(3)

  c = (pT_nu(1)*pT_nu(1) + pT_nu(2)*pT_nu(2))*p_l(0)*p_l(0) - k*k

  pz_nu_solutions = solve_quadratic(a, b, c)

  pz_nu = 1.d0

  longitudinal_neutrino_momentum = pz_nu

  return

end function longitudinal_neutrino_momentum
