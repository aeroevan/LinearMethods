program test
    use linearmethods
    implicit none
    real(kind=wp), allocatable :: A(:,:), b(:), x(:)
    integer, parameter :: n = 500

    allocate(A(n,n), b(n), x(n))

    call init_random_seed()
    call random_number(A)
    call random_number(b)

    x = solve_lu(A,b)

    !print *, x
end program test
