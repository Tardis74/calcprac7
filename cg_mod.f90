module cg_mod
    use precision_mod
    implicit none
    private
    public :: cg_minimize

    abstract interface
        function func_iface(x) result(f)
            import :: dp
            real(dp), intent(in) :: x(:)
            real(dp) :: f
        end function func_iface
    end interface

    real(dp), parameter :: c1 = 1.0e-4_dp
    real(dp), parameter :: c2 = 0.9_dp
    integer,  parameter :: max_ls   = 30
    integer,  parameter :: max_zoom = 30

contains

	function cg_minimize(f, x0, eps, max_iter, restart_period, restart_gamma, &
                         beta_method, status, iterations, ftol) result(xmin)
        procedure(func_iface) :: f
        real(dp), intent(in) :: x0(:)
        real(dp), intent(in) :: eps
        integer, intent(in),  optional :: max_iter
        integer, intent(in),  optional :: restart_period
        real(dp), intent(in),  optional :: restart_gamma
        integer, intent(in),  optional :: beta_method
        integer, intent(out), optional :: status
        integer, intent(out), optional :: iterations
        real(dp), intent(in),  optional :: ftol

        real(dp), allocatable :: xmin(:)

        integer :: n, maxit, period, method, iter
        real(dp) :: tol, gamma, beta, alpha, gnorm2
        real(dp) :: tol_f, f_prev, f_cur, f_best
        real(dp), allocatable :: x(:), g(:), p(:), g_prev(:)
        real(dp), allocatable :: x_new(:), g_new(:)
        real(dp), allocatable :: x_best(:)              ! лучшая найденная точка
        logical :: ok, restart, converged

        n = size(x0)
        maxit = 1000
        if (present(max_iter)) maxit = max_iter
        period  = n
        if (present(restart_period)) period = restart_period
        gamma   = 0.2_dp
        if (present(restart_gamma)) gamma = restart_gamma
        method  = 2
        if (present(beta_method)) method = beta_method
        tol = eps
        tol_f = 1.0e-12_dp
        if (present(ftol)) tol_f = ftol

        allocate(x(n), g(n), p(n), g_prev(n), x_new(n), g_new(n), x_best(n))
        x = x0
        call numerical_gradient(f, x, g)
        p = -g
        iter = 0
        gnorm2 = dot_product(g, g)

        ! "лучшее" значение по функции.
        f_cur = f(x)
        f_best = f_cur
        x_best = x

        converged = .false.

        do while (.not. converged .and. iter < maxit)
            if (dot_product(g, p) >= 0.0_dp) p = -g

            g_prev = g
            f_prev = f_cur

            call line_search(f, x, p, g, alpha, x_new, g_new, f_cur, ok)
            if (.not. ok) then
                p = -g
                call line_search(f, x, p, g, alpha, x_new, g_new, f_cur, ok)
                if (.not. ok) exit
            end if

            x = x_new
            g = g_new
            iter = iter + 1
            gnorm2 = dot_product(g, g)

            ! Обновление лучшей точки, если функция уменьшилась.
            if (f_cur < f_best) then
                f_best = f_cur
                x_best = x
            end if

            ! Критерии остановки
            if (sqrt(gnorm2) <= tol) then
                converged = .true.
                exit
            end if
            if (abs(f_cur - f_prev) <= tol_f * (1.0_dp + abs(f_prev))) then
                ! Отсутствие дальнейшего убывания
                converged = .true.
                exit
            end if

            ! Вычисление beta
            if (method == 1) then
                beta = gnorm2 / dot_product(g_prev, g_prev)
            else
                beta = dot_product(g, g - g_prev) / dot_product(g_prev, g_prev)
                if (beta < 0.0_dp) beta = 0.0_dp
            end if

            ! Перезапуск
            restart = .false.
            if (period > 0) then
                if (mod(iter, period) == 0) restart = .true.
            end if
            if (dot_product(g, g_prev) > gamma * gnorm2) restart = .true.

            if (restart) then
                p = -g
            else
                p = -g + beta * p
            end if
        end do

        ! Лучшая найденная точка.
        allocate(xmin(n))
        xmin = x_best

        if (present(status)) then
            if (converged) then
                status = 0
            else
                status = 1
            end if
        end if
        if (present(iterations)) iterations = iter
    end function cg_minimize

    ! Линейный поиск
    subroutine line_search(f, x, p, g, alpha, x_new, g_new, f_new, success)
        procedure(func_iface) :: f
        real(dp), intent(in)  :: x(:), p(:), g(:)
        real(dp), intent(out) :: alpha, x_new(:), g_new(:), f_new
        logical, intent(out) :: success

        real(dp) :: phi0, derphi0, alpha_prev, phi_prev, alpha_i, phi_i, derphi_i
        real(dp), parameter :: alpha_max = 10.0_dp
        integer :: it

        phi0 = f(x)
        derphi0 = dot_product(g, p)
        if (derphi0 >= 0.0_dp) then
            success = .false.; alpha = 0.0_dp; x_new = x; g_new = g; f_new = phi0
            return
        end if

        alpha_prev = 0.0_dp
        phi_prev = phi0
        alpha_i = 1.0_dp
        success = .false.

        do it = 1, max_ls
            x_new = x + alpha_i * p
            phi_i = f(x_new)
            if (phi_i > phi0 + c1 * alpha_i * derphi0 .or. &
                (it > 1 .and. phi_i >= phi_prev)) then
                call zoom(f, x, p, phi0, derphi0, alpha_prev, alpha_i, &
                          alpha, x_new, g_new, f_new, success)
                return
            end if
            call numerical_gradient(f, x_new, g_new)
            derphi_i = dot_product(g_new, p)
            if (abs(derphi_i) <= -c2 * derphi0) then
                alpha = alpha_i; f_new = phi_i; success = .true.; return
            end if
            if (derphi_i >= 0.0_dp) then
                call zoom(f, x, p, phi0, derphi0, alpha_i, alpha_prev, &
                          alpha, x_new, g_new, f_new, success)
                return
            end if
            alpha_prev = alpha_i
            phi_prev   = phi_i
            alpha_i    = min(alpha_max, alpha_i * 2.0_dp)
        end do

        success = .false.
        alpha = alpha_prev
        x_new = x + alpha * p
        call numerical_gradient(f, x_new, g_new)
        f_new = f(x_new)
    end subroutine line_search

    subroutine zoom(f, x, p, phi0, derphi0, alpha_lo_in, alpha_hi_in, &
                    alpha_out, x_new, g_new, f_new, success)
        procedure(func_iface) :: f
        real(dp), intent(in)  :: x(:), p(:), phi0, derphi0, alpha_lo_in, alpha_hi_in
        real(dp), intent(out) :: alpha_out, x_new(:), g_new(:), f_new
        logical,  intent(out) :: success

        real(dp) :: alpha_lo, alpha_hi, alpha, phi_lo, phi, derphi
        integer  :: it
        real(dp) :: xt(size(x))

        alpha_lo = alpha_lo_in
        alpha_hi = alpha_hi_in
        xt       = x + alpha_lo * p
        phi_lo   = f(xt)
        success  = .false.

        do it = 1, max_zoom
            alpha = 0.5_dp * (alpha_lo + alpha_hi)
            xt    = x + alpha * p
            phi   = f(xt)

            if (phi > phi0 + c1 * alpha * derphi0 .or. phi >= phi_lo) then
                alpha_hi = alpha
            else
                call numerical_gradient(f, xt, g_new)
                derphi = dot_product(g_new, p)
                if (abs(derphi) <= -c2 * derphi0) then
                    alpha_out = alpha; x_new = xt; f_new = phi
                    success   = .true.; return
                end if
                if (derphi * (alpha_hi - alpha_lo) >= 0.0_dp) alpha_hi = alpha_lo
                alpha_lo = alpha
                phi_lo   = phi
            end if
            if (abs(alpha_hi - alpha_lo) < 1.0e-12_dp) exit
        end do

        alpha = alpha_lo
        x_new = x + alpha * p
        call numerical_gradient(f, x_new, g_new)
        f_new = f(x_new)
        alpha_out = alpha
        success = (f_new < phi0)
    end subroutine zoom

    ! Численный градиент
    subroutine numerical_gradient(f, x, g)
        procedure(func_iface) :: f
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: g(:)
        integer :: i, n
        real(dp) :: xi, h, fp, fm
        real(dp) :: xt(size(x))

        n = size(x)
        xt = x
        do i = 1, n
            xi = x(i)
            h  = epsilon(1.0_dp)**(1.0_dp / 3.0_dp) * max(1.0_dp, abs(xi))
            xt(i) = xi + h;  fp = f(xt)
            xt(i) = xi - h;  fm = f(xt)
            xt(i) = xi
            g(i)  = (fp - fm) / (2.0_dp * h)
        end do
    end subroutine numerical_gradient

end module cg_mod
