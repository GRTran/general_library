module gauss_quadrature
  use VarPrecision
  use BlaspackInterface
  implicit none


contains

  subroutine determine_quadrature_rules( number_quad_points, relation_func, abcsicca, weights )
    !! determines the gauss quadrature points and associated weighting using the Golub Welsch algorithm. The matrix
    !! Ax = bx is solved for the eigenvalues b which in turn are used to find eigenvectors x_i. The eigenvectors
    !! are then used to determine the weights, the abcsicca are the eigenvalues or roots of the n-th order polynomial.
    !! The order of the polynomial is the same as the number of quadrature points.
    integer, intent(in) :: number_quad_points
    real(8), intent(out) :: abcsicca(:)
    real(8), intent(out) :: weights(:)

    integer :: i
    integer :: order_polynomial
    integer :: status
    real(8) :: alpha, beta_right, beta_left
    real(8), allocatable :: J(:,:)
    real(8), allocatable :: eigenvalues(:)
    real(8), allocatable :: eigenvectors(:,:)
    real(8) :: recurrence_relation(3)

    interface
      function relation_func( i, coeff ) result( three_term_weight )
        integer, intent(in) :: i
        integer, intent(in) :: coeff
        real(8) :: step
        real(8) :: three_term_weight
      end function
    end interface

    order_polynomial = number_quad_points + 1

    allocate( eigenvalues(order_polynomial), J(order_polynomial, order_polynomial), &
              eigenvectors(order_polynomial, order_polynomial), stat=status)
    if( status /= 0 ) stop 'Error cannot allocate temporary gauss quadrature array'

    ! form the tridiagonal symmetric matrix
    J = 0._wp
    alpha = -relation_func(1, 2)/relation_func(1, 1)
    beta_right = sqrt(relation_func(2,3)/(relation_func(1,1)*relation_func(2,1)))
    J(1,1) = alpha
    J(1,2) = beta_right
    do i = 2, order_polynomial - 1
      alpha = -relation_func(i, 2)/relation_func(i, 1)
      beta_left = beta_right
      beta_right = sqrt(relation_func(i+1,3)/(relation_func(i,1)*relation_func(i+1,1)))
      J(i,i) = alpha
      J(i,i+1) = beta_right
      J(i,i-1) = beta_left
    enddo
    alpha = -relation_func(order_polynomial, 2)/relation_func(order_polynomial, 1)
    J(order_polynomial,order_polynomial) = alpha
    J(order_polynomial-1, order_polynomial) = beta_right

    ! solve the tridiagonal matrix for all of the eigenvalues and corresponding eigenvectors
    call lapack_eigenvalues_eigenvectors( J, eigenvalues, eigenvectors )
    ! calculate the eigenvalues
    abcsicca = eigenvalues(1:number_quad_points)
    ! calculate the corresponding weights using the first component of the eigenvector
    weights = eigenvectors(1,1:size(weights))**2*relation_func(0,4)

    deallocate(eigenvalues, J, eigenvectors, stat=status)
    if( status /= 0 ) stop 'Error cannot deallocate temporary gauss quadrature array'
  end subroutine

  function legendre( i, coeff ) result( three_term_weight )
    integer, intent(in) :: i
    integer, intent(in) :: coeff
    real(wp) :: step
    real(wp) :: three_term_weight

    step = i
    if( coeff == 1 ) then
      three_term_weight = (2._wp*(step-1._wp)+1._wp)/(step)
    elseif( coeff == 2 ) then
      three_term_weight = 0._wp
    elseif ( coeff == 3 ) then
      three_term_weight = (step-1._wp)/(step)
    elseif ( coeff == 4 ) then
      three_term_weight = 2._wp
    endif
  end function

  function hermite( i, coeff ) result( three_term_weight )
    integer, intent(in) :: i
    integer, intent(in) :: coeff
    real(wp) :: step
    real(wp) :: three_term_weight

    step = i
    if( coeff == 1 ) then
      three_term_weight = 2._wp
    elseif( coeff == 2 ) then
      three_term_weight = 0._wp
    elseif ( coeff == 3 ) then
      three_term_weight = 2._wp*(step-1._wp)
    elseif ( coeff == 4 ) then
      three_term_weight = sqrt(pi)
    endif
  end function

  subroutine gauss_quad_test(tolerance, number_quad_points, quadrature_method)
    real, intent(in) :: tolerance
    integer, intent(in) :: number_quad_points
    character(len=*), intent(in) :: quadrature_method
    real(8), allocatable :: abcsicca(:)
    real(8), allocatable :: weights(:)
    integer :: status
    integer :: i
    real :: numerical_integral, exact_integral

    allocate( abcsicca(number_quad_points), weights(number_quad_points), stat=status )
    if( status /= 0 ) stop 'Cannot allocate quadrature result array, not enough ram'
    ! determine the quadrature rules depending on the integral weight and 3 term recurrence relation
    if( quadrature_method == 'legendre' ) then
      call determine_quadrature_rules( number_quad_points, legendre, abcsicca, weights )
    endif
    ! calculate the numerical integral of the function
    numerical_integral = 0._wp
    do i = 1, number_quad_points
      numerical_integral = numerical_integral + cos(abcsicca(i))*weights(i)
    enddo
    exact_integral = sin(1._wp) - sin(-1._wp)
    if( abs(numerical_integral - exact_integral) < tolerance ) then
      write(*,*) 'Success: Quadrature Test 1 of 1 passed', abs(numerical_integral - exact_integral)
    else
      write(*,*) 'Fail: Quadrature Test 0 of 1 passed', abs(numerical_integral - exact_integral)
    endif
    deallocate( abcsicca, weights, stat=status )
    if( status /= 0 ) stop 'Cannot deallocate quadrature result array, not enough ram'
  end subroutine

end module
