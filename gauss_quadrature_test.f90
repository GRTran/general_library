program gauss_quadrature_test
  use VarPrecision
  use BlaspackInterface
  use gauss_quadrature
  implicit none

  integer :: number_quad_points
  real(8), allocatable :: abcsicca(:)
  real(8), allocatable :: weights(:)
  real(8) :: integral_weight
  integer :: status
  integer :: i
  real :: numerical_integral, exact_integral

  number_quad_points = 3

  allocate( abcsicca(number_quad_points), weights(number_quad_points), stat=status )
  if( status /= 0 ) stop 'Cannot allocate quadrature result array, not enough ram'


  !! integral of weight function w(x) dx
  integral_weight = 2._wp



  ! determine the quadrature rules depending on the integral weight and 3 term recurrence relation
  call determine_quadrature_rules( number_quad_points, legendre, abcsicca, weights )

  ! calculate the numerical integral of the function
  numerical_integral = 0._wp
  do i = 1, number_quad_points
    numerical_integral = numerical_integral + cos(abcsicca(i))*weights(i)
  enddo

  exact_integral = sin(1._wp) - sin(-1._wp)

  write(*,*) numerical_integral, exact_integral

contains



end program
