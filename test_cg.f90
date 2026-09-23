program test_cg
    use precision_mod
    use cg_mod
    use test_funcs_mod
    implicit none

    real(dp), allocatable :: x0(:), xmin(:)
    integer  :: status, beta_method, max_iter, restart_period, iter, u
    real(dp) :: eps, restart_gamma, ftol

    open(newunit=u, file='result.txt', status='replace', action='write')

    !Квадратичная
    eps = 1.0e-8_dp
    max_iter = 2000
    restart_period = 5
    restart_gamma = 0.2_dp
    beta_method = 2
    ftol = 1.0e-14_dp
    x0 = [0.0_dp, 0.0_dp, 0.0_dp]

    call report(u, "Тест 1: квадратичная функция", x0, eps)
    xmin = cg_minimize(quadratic, x0, eps, max_iter, restart_period, &
                       restart_gamma, beta_method, status, iter, ftol)
    call report_result(u, xmin, quadratic(xmin), status, iter)

    !Розенброк
    eps = 1.0e-7_dp
    max_iter = 50000
    restart_period = 0
    restart_gamma = 0.1_dp
    beta_method = 2
    ftol = 1.0e-14_dp
    x0 = [-1.2_dp, 1.0_dp]

    call report(u, "Тест 2: функция Розенброка", x0, eps)
    xmin = cg_minimize(rosenbrock, x0, eps, max_iter, restart_period, &
                       restart_gamma, beta_method, status, iter, ftol)
    call report_result(u, xmin, rosenbrock(xmin), status, iter)

    !Химмельблау
    eps = 1.0e-8_dp
    max_iter = 2000
    restart_period = 5
    restart_gamma = 0.2_dp
    beta_method = 2
    ftol = 1.0e-14_dp
    x0 = [1.0_dp, 1.0_dp]

    call report(u, "Тест 3: функция Химмельблау", x0, eps)
    xmin = cg_minimize(himmelblau, x0, eps, max_iter, restart_period, &
                       restart_gamma, beta_method, status, iter, ftol)
    call report_result(u, xmin, himmelblau(xmin), status, iter)

    close(u)

contains

    subroutine report(u, title, x0, eps)
        integer,  intent(in) :: u
        character(*), intent(in) :: title
        real(dp), intent(in) :: x0(:), eps
        character(len=64) :: buf

        write(buf,'(A)') "=== " // trim(title) // " ==="
        write(*,'(A)') trim(buf); write(u,'(A)') trim(buf)
        write(buf,'(A,ES12.4)') "  eps  = ", eps
        write(*,'(A)') trim(buf); write(u,'(A)') trim(buf)
        write(buf,'(A,*(F10.4,1X))') "  x0   = ", x0
        write(*,'(A)') trim(buf); write(u,'(A)') trim(buf)
    end subroutine report

    subroutine report_result(u, xmin, fmin, status, iter)
        integer,  intent(in) :: u, status, iter
        real(dp), intent(in) :: xmin(:), fmin
        character(len=80) :: buf

        write(buf,'(A,*(ES15.6,1X))') "  xmin    = ", xmin
        write(*,'(A)') trim(buf); write(u,'(A)') trim(buf)
        write(buf,'(A,ES15.6)') "  f(xmin) = ", fmin
        write(*,'(A)') trim(buf); write(u,'(A)') trim(buf)
        write(buf,'(A,I0)') "  status  = ", status
        write(*,'(A)') trim(buf); write(u,'(A)') trim(buf)
        write(buf,'(A,I0)') "  iter    = ", iter
        write(*,'(A)') trim(buf); write(u,'(A)') trim(buf)
        write(*,'(A)') ""; write(u,'(A)') ""
    end subroutine report_result

end program test_cg
