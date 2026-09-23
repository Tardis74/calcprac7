! Модуль с тестовыми функциями.
module test_funcs_mod
    use precision_mod
    implicit none
    private
    public :: quadratic, rosenbrock, himmelblau

contains

    function quadratic(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        f = (x(1) - 1.0_dp)**2 + 2.0_dp * (x(2) + 3.0_dp)**2 + &
            0.5_dp * (x(3) - 2.0_dp)**2
    end function quadratic

    function rosenbrock(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        f = (1.0_dp - x(1))**2 + 100.0_dp * (x(2) - x(1)**2)**2
    end function rosenbrock

    function himmelblau(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        f = (x(1)**2 + x(2) - 11.0_dp)**2 + &
            (x(1) + x(2)**2 - 7.0_dp)**2
    end function himmelblau

end module test_funcs_mod
