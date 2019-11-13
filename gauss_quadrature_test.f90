program gauss_quadrature_test
  use VarPrecision
  use BlaspackInterface
  use gauss_quadrature
  use HermiteNISP
  implicit none

  integer :: number_quad_points
  real(8), allocatable :: abcsicca(:)
  real(8), allocatable :: weights(:)
  integer :: status
  integer :: i

  ! nisp for hermite polynomial to represent a normal-distribution
  type(hermitePCClass) :: hermite_pc
  real(wp), allocatable :: input_random_numbers(:)
  real(wp), allocatable :: pdf(:)
  real(wp), allocatable :: pc_sample_responses(:)
  real(wp), allocatable :: monte_carlo_sample_responses(:)
  real(wp) :: mean
  real(wp) :: variance
  real(wp) :: uniform_rand_num1, uniform_rand_num2
  integer :: poly_order
  integer :: num_samples


  number_quad_points = 10
  poly_order = 2

  mean = 0.5_wp
  variance = 4._wp

  num_samples = 5000

  allocate( pc_sample_responses(num_samples), monte_carlo_sample_responses(num_samples), &
            input_random_numbers(num_samples), pdf(num_samples), stat=status)
  if( status /= 0 ) stop 'Cannot allocate sample data result array, not enough ram'

  allocate( abcsicca(number_quad_points), weights(number_quad_points), stat=status )
  if( status /= 0 ) stop 'Cannot allocate quadrature result array, not enough ram'


  ! determine the quadrature rules depending on the integral weight and 3 term recurrence relation
  call determine_quadrature_rules( number_quad_points, hermite, abcsicca, weights )

  ! test the polynomial chaos expansion of the function that depends on a random variable
  call hermite_pc%initialise(poly_order, mean, variance, abcsicca, weights, myfunc )
  write(*,*) hermite_pc%response_poly_coeffs
  write(*,*) hermite_pc%nisp_evaluation(1._wp)

  do i = 1, num_samples
    ! generating a normal distribution using the acceptance/rejection algorithm
    do while (uniform_rand_num2 < (uniform_rand_num1 - 1)**2/2._wp)
      call RANDOM_NUMBER(uniform_rand_num1)
      call RANDOM_NUMBER(uniform_rand_num2)
      uniform_rand_num1 = -log(uniform_rand_num1)
      uniform_rand_num2 = -log(uniform_rand_num2)
    enddo
    call RANDOM_NUMBER(uniform_rand_num2)
    if( uniform_rand_num2 <= 0.5_wp ) then
      uniform_rand_num1 = abs(uniform_rand_num1)
    else
      uniform_rand_num1 = -abs(uniform_rand_num1)
    endif
    ! calculating x and pdf
    input_random_numbers(i) = sqrt(variance)*uniform_rand_num1 + mean
    pdf(i) = (1._wp/sqrt(2._wp*pi*variance)) * exp(- (input_random_numbers(i)-mean)**2/(2._wp*variance))

    pc_sample_responses(i) = hermite_pc%nisp_evaluation(input_random_numbers(i))
    monte_carlo_sample_responses(i) = cos(input_random_numbers(i))
  enddo

  ! write information to file for plotting
  open(145,file='random_inputs.dat')
  do i = 1, size(input_random_numbers)
    write(145,*) input_random_numbers(i), pdf(i), pc_sample_responses(i), monte_carlo_sample_responses(i)
  enddo

  ! deallocate variables
  deallocate( pc_sample_responses, monte_carlo_sample_responses, input_random_numbers, pdf )
  deallocate( abcsicca, weights )

contains

  real(wp) function myfunc( x )
    real(wp), intent(in) :: x
    myfunc = cos(x)
  end function



end program
