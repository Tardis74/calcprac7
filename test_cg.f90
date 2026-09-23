program test_cg
    use precision_mod
    use cg_mod,       only: cg_minimize, func_iface
    use test_funcs_mod
    implicit none

    integer :: u
    real(dp) :: eps, restart_gamma, ftol
    integer :: max_iter, restart_period

    open(newunit=u, file='result.txt', status='replace', action='write')

    restart_gamma  = 0.2_dp
    ftol = 1.0e-14_dp

    ! Тест 1: квадратичная функция
    call compare_methods(u, quadratic, &
        "Тест 1: квадратичная функция", &
        [0.0_dp, 0.0_dp, 0.0_dp], &
        eps=1.0e-8_dp, max_iter=2000, restart_period=5, &
        restart_gamma=restart_gamma, ftol=ftol)

    ! Тест 2: функция Розенброка
    call compare_methods(u, rosenbrock, &
        "Тест 2: функция Розенброка", &
        [-1.2_dp, 1.0_dp], &
        eps=1.0e-7_dp, max_iter=50000, restart_period=0, &
        restart_gamma=0.1_dp, ftol=ftol)

    ! Тест 3: функция Химмельблау
    call compare_methods(u, himmelblau, &
        "Тест 3: функция Химмельблау", &
        [1.0_dp, 1.0_dp], &
        eps=1.0e-8_dp, max_iter=2000, restart_period=5, &
        restart_gamma=restart_gamma, ftol=ftol)

    close(u)

contains

    subroutine compare_methods(u, f, title, x0, eps, max_iter, &
                               restart_period, restart_gamma, ftol)
        integer, intent(in) :: u
        procedure(func_iface)    :: f
        character(*), intent(in) :: title
        real(dp), intent(in) :: x0(:)
        real(dp), intent(in) :: eps, restart_gamma, ftol
        integer, intent(in) :: max_iter, restart_period

        real(dp), allocatable :: x_fr(:), x_pr(:)
        integer :: st_fr, st_pr, it_fr, it_pr
        real(dp) :: f_fr, f_pr, g_fr, g_pr
        character(len=128) :: buf

        write(buf,'(A)') "=== " // trim(title) // " ==="
        call emit(u, buf)
        write(buf,'(A,ES12.4)') "  eps    = ", eps
        call emit(u, buf)
        write(buf,'(A,*(F10.4,1X))') "  x0     = ", x0
        call emit(u, buf)
        call emit(u, "")

        ! FR (beta_method = 1)
        x_fr = cg_minimize(f, x0, eps, max_iter, restart_period, &
                           restart_gamma, 1, st_fr, it_fr, ftol, g_fr)
        f_fr = f(x_fr)
        call print_method(u, "FR ", x_fr, f_fr, g_fr, st_fr, it_fr)

        ! PR (beta_method = 2)
        x_pr = cg_minimize(f, x0, eps, max_iter, restart_period, &
                           restart_gamma, 2, st_pr, it_pr, ftol, g_pr)
        f_pr = f(x_pr)
        call print_method(u, "PR+", x_pr, f_pr, g_pr, st_pr, it_pr)

        ! Итог сравнения
        write(buf,'(A)') "  Итог: " // trim(verdict( &
            st_fr, it_fr, f_fr, st_pr, it_pr, f_pr))
        call emit(u, buf)
        call emit(u, "")
    end subroutine compare_methods

    subroutine print_method(u, tag, x, fval, gnorm, status, iter)
        integer, intent(in) :: u, status, iter
        character(*), intent(in) :: tag
        real(dp), intent(in) :: x(:), fval, gnorm
        character(len=128) :: buf, prefix

        write(prefix,'(A)') "  [" // trim(tag) // "] "
        write(buf,'(A,A,*(ES15.6,1X))') trim(prefix), "xmin    = ", x
        call emit(u, buf)
        write(buf,'(A,A,ES15.6)') trim(prefix), "f(xmin) = ", fval
        call emit(u, buf)
        write(buf,'(A,A,ES15.6)') trim(prefix), "||g||   = ", gnorm
        call emit(u, buf)
        write(buf,'(A,A,I0)')    trim(prefix), "status  = ", status
        call emit(u, buf)
        write(buf,'(A,A,I0)')    trim(prefix), "iter    = ", iter
        call emit(u, buf)
        call emit(u, "")
    end subroutine print_method

    subroutine emit(u, line)
        integer, intent(in) :: u
        character(*), intent(in) :: line
        write(*,'(A)') trim(line)
        write(u,'(A)') trim(line)
    end subroutine emit

    function verdict(st_fr, it_fr, f_fr, st_pr, it_pr, f_pr) result(s)
        integer, intent(in) :: st_fr, it_fr, st_pr, it_pr
        real(dp), intent(in) :: f_fr, f_pr
        character(len=96) :: s

        if (st_fr == 0 .and. st_pr == 0) then
            if (it_fr < it_pr) then
                write(s,'(A,I0,A,I0,A)') "FR лучше (", it_fr, &
                    " итераций против ", it_pr, ")"
            else if (it_pr < it_fr) then
                write(s,'(A,I0,A,I0,A)') "PR+ лучше (", it_pr, &
                    " итераций против ", it_fr, ")"
            else
                write(s,'(A,I0,A)') "оба метода одинаковы (", it_fr, " итераций)"
            end if
        else if (st_fr == 0) then
            s = "сошёлся только FR"
        else if (st_pr == 0) then
            s = "сошёлся только PR+"
        else
            write(s,'(A,ES10.2,A,ES10.2)') "ни один не сошёлся; f_FR = ", &
                f_fr, ", f_PR = ", f_pr
        end if
    end function verdict

end program test_cg
