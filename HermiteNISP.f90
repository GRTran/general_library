module HermiteNISP
  use gauss_quadrature
  implicit none

  type hermitePCClass
    real(wp), allocatable :: response_poly_coeffs(:)
    integer               :: poly_order
  contains
    procedure, public, pass :: initialise => initialise_hermite_pc
    procedure, public, pass :: nisp_evaluation => hermite_nisp_evaluation
  endtype

contains

  subroutine initialise_hermite_pc( this, poly_order, mean, variance, abscissa, weights, response_func )
    class(hermitePCClass), intent(out) :: this
    real(wp), intent(in) :: abscissa(:)
    real(wp), intent(in) :: weights(:)
    real(wp), intent(in) :: mean
    real(wp), intent(in) :: variance
    integer , intent(in) :: poly_order

    integer :: i, j
    integer :: status

    interface
      real(wp) function response_func( x )
        import wp
        real(wp), intent(in) :: x
      end function
    end interface

    this%poly_order = poly_order
    allocate( this%response_poly_coeffs(poly_order+1), stat=status )
    if( status /= 0 ) stop 'Error allocating response polynomial coefficients for Hermite polynomial, not enough RAM'

    this%response_poly_coeffs = 0._wp

    ! perform a loop to determine the polynomial response coefficient R_i for all coefficients in total order
    do j = 0, poly_order
      ! perform gauss-hermite integral with f(x) = R(x)*He(x) of f(x)w(x) where w(x) is the hermite weight function
      ! adjust result for a non-standard normal distribution
      ! write(*,*) hermite_expansion(sqrt(2._wp)*abscissa(1),j)
      do i = 1, size(abscissa)
        this%response_poly_coeffs(j+1) = this%response_poly_coeffs(j+1) + hermite_expansion(sqrt(2._wp)*abscissa(i),j)*response_func((mean+sqrt(variance)*sqrt(2._wp)*abscissa(i)))*weights(i)
      enddo
      write(*,*) this%response_poly_coeffs(j+1) , j, factorial(j)
      this%response_poly_coeffs(j+1) = this%response_poly_coeffs(j+1) * sqrt(2._wp) / (sqrt(2._wp*pi)*factorial(j))
      ! write(*,*) this%response_poly_coeffs(j+1) , j, factorial(j)
    enddo

  end subroutine

  real(wp) function hermite_nisp_evaluation( this, x )
    class(hermitePCClass), intent(in) :: this
    real(wp)             , intent(in) :: x

    integer :: i

    hermite_nisp_evaluation = 0._wp
    do i = 0, this%poly_order
      hermite_nisp_evaluation = hermite_nisp_evaluation + this%response_poly_coeffs(i+1) * hermite_expansion(x,i)
    enddo

  end function

  real(wp) function hermite_expansion(x, order)
    real(wp), intent(in) :: x
    integer , intent(inout) :: order

    integer :: i
    real(wp) :: temp


    hermite_expansion = 0._wp
    do i = 0, order/2
      ! write(*,*) 'it',i
      temp = ((-1)**i / (real(factorial(i)) * real(factorial(order-2*i)))) * (x**(order-2._wp*i)/2._wp**i)
      ! write(*,*) 'temp val', temp
      hermite_expansion = hermite_expansion + temp
    enddo
    hermite_expansion = hermite_expansion*real(factorial(order))
  end function

  integer function factorial(n)
    integer :: n

    integer :: i

    factorial = 1
    do i = 1, n
      factorial = factorial*i
    enddo
  end function

  subroutine hermite_nisp_test()
    type(hermitePCClass) :: test_pc

    real(wp) :: x
    integer :: order

    x = 2._wp
    order = 2
    write(*,*) hermite_expansion(x, order), x**2 - 1._wp
  end subroutine

end module
