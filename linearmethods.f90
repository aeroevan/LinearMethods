module linearmethods
    use mpi
    implicit none
    integer, parameter :: sp = kind(1.0), dp = kind(1.0d0)
    ! Working precision
    integer, parameter :: wp = dp
contains
    ! {{{ Givens QR routines...
    ! {{{ givens(a, b, c, s)
    pure subroutine givens(a, b, c, s)
        real(kind=wp), intent(in) :: a, b
        real(kind=wp), intent(out) :: c, s
        real(kind=wp) :: tau

        if (b == 0.0_wp) then
            c = 1.0_wp
            s = 0.0_wp
        else if (abs(b) > abs(a)) then
            tau = -a/b
            s = 1.0_wp/(sqrt(1.0_wp+tau*tau))
            c = s*tau
        else
            tau = -b/a
            c = 1.0_wp/(sqrt(1.0_wp+tau*tau))
            s = c*tau
        end if
    end subroutine givens
    ! }}}
    ! {{{ givens_row(A, c, s)
    subroutine givens_row(A, c, s)
        real(kind=wp), intent(inout) :: A(:,:)
        real(kind=wp), intent(in) :: c, s
        real(kind=wp) :: tau, sigma
        integer :: k, q

        q = size(A,2)

        do k = 1, q
            tau = A(1,k)
            sigma = A(2,k)
            A(1,k) = c*tau - s*sigma
            A(2,k) = s*tau + c*sigma
        end do
    end subroutine givens_row
    ! }}}
    ! {{{ givens_qr(A, Q, R)
    subroutine old_givens_qr(A, Q, R)
        real(kind=wp), intent(in) :: A(:,:)
        real(kind=wp), intent(out) :: Q(size(A,1), size(A,2)), R(size(A,1), size(A,2))
        integer :: n, m, i, j, T
        real(kind=wp) :: c, s
        logical :: updated

        n = size(A,1)
        m = size(A,2)
        
        Q = 0.0_wp
        do i = 1, n
            Q(i,i) = 1.0_wp
        end do

        updated = .true.
        T = 1

        R = A

        do while(updated)
            updated = .false.
            !$OMP parallel do private(c,s)
            do i = n, 2, -1
                do j = 1, i-1
                    if (i-2*j==n-1-T) then
                        updated = .true.
                        call givens(R(i-1,j), R(i,j), c, s)
                        ! Only need to update j:m here (rest should be zero)
                        call givens_row(R(i-1:i,j:m), c, s)
                        call givens_row(Q(i-1:i,1:m), c, s)
                    end if
                end do
            end do
            !$OMP end parallel do
            T = T + 1
        end do
        Q = transpose(Q)

        !print *, "A:"
        !call print_matrix(A)
        !print *, "Q*R:"
        !call print_matrix(matmul(Q,R))
    end subroutine old_givens_qr
    ! }}}
    subroutine givens_qr(A, Q, R)
        real(kind=wp), intent(in) :: A(:,:)
        real(kind=wp), intent(out) :: Q(size(A,1), size(A,2)), R(size(A,1), size(A,2))
        integer :: n, m
        ! MPI stuff
        integer :: id, ierr

        n = size(A,1)
        m = size(A,2)
        call MPI_BCAST(m, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        call MPI_COMM_RANK(MPI_COMM_WORLD, id, ierr)
        if (id == 0) then
            call main_givens_qr(A, Q, R)
        else
            call worker_givens_qr(m)
        end if
    end subroutine givens_qr
    !{{{ givens_schedule(n, T, is, js)
    subroutine givens_schedule(n, T, is, js)
        integer, intent(in) :: n, T
        integer, allocatable, intent(out) :: is(:), js(:)
        integer, allocatable :: tmp(:)
        integer :: i, j, s
        if (allocated(is)) then
            deallocate(is)
        end if
        if (allocated(js)) then
            deallocate(js)
        end if
        do i = n, 2, -1
            do j = 1, i-1
                if (i-2*j==n-1-T) then
                    if (allocated(is)) then
                        s = size(is,1)
                        allocate(tmp(s+1))
                        tmp(1:s) = is
                        tmp(s+1) = i
                        ! move_alloc is f2003...
                        call move_alloc(from=tmp, to=is)
                    else
                        allocate(is(1))
                        is(1) = i
                    end if
                    if (allocated(js)) then
                        s = size(js,1)
                        allocate(tmp(s+1))
                        tmp(1:s) = js
                        tmp(s+1) = j
                        call move_alloc(from=tmp, to=js)
                    else
                        allocate(js(1))
                        js(1) = j
                    end if
                end if
            end do
        end do
    end subroutine givens_schedule
    ! }}}
    !{{{ main_givens_qr(A, Q, R)
    subroutine main_givens_qr(A, Q, R)
        real(kind=wp), intent(in) :: A(:,:)
        real(kind=wp), intent(out) :: Q(size(A,1), size(A,2)), R(size(A,1), size(A,2))
        integer :: n, m, i, j, T, number_sent, number_to_send, number_received, status(MPI_STATUS_SIZE)
        integer, allocatable :: is(:), js(:)
        ! MPI stuff
        integer :: ierr, nproc, nsend, source, dest

        n = size(A,1)
        m = size(A,2)
        Q = 0.0_wp
        forall(i=1:n) Q(i,i) = 1.0_wp

        R = A

        T = 1

        call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
        nproc = nproc - 1

        do while(.true.)
            number_sent = 0
            number_received = 0
            call givens_schedule(n, T, is, js)
            if (allocated(is)) then
                number_to_send = size(is,1)
                nsend = min(nproc, number_to_send)
                do i = 1, nsend
                    !print *, i, number_to_send, "Send"
                    call MPI_SEND(is(i), 1, MPI_INTEGER, i, i, MPI_COMM_WORLD, ierr)
                    call MPI_SEND(js(i), 1, MPI_INTEGER, i, i, MPI_COMM_WORLD, ierr)
                    if (wp == sp) then
                        call MPI_SEND(R(is(i)-1:is(i),js(i):m), 2*(m-js(i)+1), MPI_REAL, i, i, MPI_COMM_WORLD, ierr)
                    else
                        call MPI_SEND(R(is(i)-1:is(i),js(i):m), 2*(m-js(i)+1), MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, ierr)
                    end if
                    if (wp == sp) then
                        call MPI_SEND(Q(is(i)-1:is(i),1:m), 2*m, MPI_REAL, i, i, MPI_COMM_WORLD, ierr)
                    else
                        call MPI_SEND(Q(is(i)-1:is(i),1:m), 2*m, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, ierr)
                    end if
                end do
                number_sent = nsend
                i = 0
                do while (number_received < number_to_send)
                    i = i + 1
                    if (wp == sp) then
                        call MPI_RECV(R(is(i)-1:is(i),js(i):m), 2*(m-js(i)+1), &
                            MPI_REAL, MPI_ANY_SOURCE, i, MPI_COMM_WORLD, status, ierr)
                    else
                        call MPI_RECV(R(is(i)-1:is(i),js(i):m), 2*(m-js(i)+1), &
                            MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, i, MPI_COMM_WORLD, status, ierr)
                    end if
                    if (wp == sp) then
                        call MPI_RECV(Q(is(i)-1:is(i),1:m), 2*m, &
                            MPI_REAL, MPI_ANY_SOURCE, i, MPI_COMM_WORLD, status, ierr)
                    else
                        call MPI_RECV(Q(is(i)-1:is(i),1:m), 2*m, &
                            MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, i, MPI_COMM_WORLD, status, ierr)
                    end if
                    source = status(MPI_SOURCE)
                    dest = source
                    number_received = number_received + 1
                    if (number_sent < number_to_send) then
                        number_sent = number_sent + 1
                        j = number_sent
                        call MPI_SEND(is(j), 1, MPI_INTEGER, dest, j, MPI_COMM_WORLD, ierr)
                        call MPI_SEND(js(j), 1, MPI_INTEGER, dest, j, MPI_COMM_WORLD, ierr)
                        if (wp == sp) then
                            call MPI_SEND(R(is(j)-1:is(j),js(j):m), 2*(m-js(j)+1), &
                                MPI_REAL, dest, j, MPI_COMM_WORLD, ierr)
                        else
                            call MPI_SEND(R(is(j)-1:is(j),js(j):m), 2*(m-js(j)+1), &
                                MPI_DOUBLE_PRECISION, dest, j, MPI_COMM_WORLD, ierr)
                        end if
                        if (wp == sp) then
                            call MPI_SEND(Q(is(j)-1:is(j),1:m), 2*m, &
                                MPI_REAL, dest, j, MPI_COMM_WORLD, ierr)
                        else
                            call MPI_SEND(Q(is(j)-1:is(j),1:m), 2*m, &
                                MPI_DOUBLE_PRECISION, dest, j, MPI_COMM_WORLD, ierr)
                        end if
                    end if
                end do
            else
                exit
            end if
            deallocate(is, js)
            T = T + 1
        end do
        do i = 1, nproc
            call MPI_SEND(0, 1, MPI_INTEGER, i, i, MPI_COMM_WORLD, ierr)
        end do

    end subroutine main_givens_qr
    ! }}}
    subroutine worker_givens_qr(m)
        integer, intent(in) :: m
        integer :: i, j, ierr, source, tag, status(MPI_STATUS_SIZE)
        real(kind=wp) :: c, s
        real(kind=wp), allocatable :: R(:,:), Q(:,:)

        do
            call MPI_RECV(i, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
            if (i == 0) then
                call MPI_FINALIZE(ierr)
                stop
            end if
            call MPI_RECV(j, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
            allocate(R(2,j:m), Q(2,m))
            if (wp == sp) then
                call MPI_RECV(R(1:2,j:m), 2*(m-j+1), &
                    MPI_REAL, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
            else
                call MPI_RECV(R(1:2,j:m), 2*(m-j+1), &
                    MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
            end if
            if (wp == sp) then
                call MPI_RECV(Q(1:2,1:m), 2*m, &
                    MPI_REAL, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
            else
                call MPI_RECV(Q(1:2,1:m), 2*m, &
                    MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
            end if
            source = status(MPI_SOURCE)
            tag = status(MPI_TAG)

            call givens(R(1,j), R(2,j), c, s)
            call givens_row(R(1:2,j:m), c, s)
            call givens_row(Q(1:2,1:m), c, s)

            if (wp == sp) then
                call MPI_SEND(R(1:2,j:m), 2*(m-j+1), MPI_REAL, source, tag, MPI_COMM_WORLD, ierr)
            else
                call MPI_SEND(R(1:2,j:m), 2*(m-j+1), MPI_DOUBLE_PRECISION, source, tag, MPI_COMM_WORLD, ierr)
            end if
            if (wp == sp) then
                call MPI_SEND(Q(1:2,1:m), 2*m, MPI_REAL, source, tag, MPI_COMM_WORLD, ierr)
            else
                call MPI_SEND(Q(1:2,1:m), 2*m, MPI_DOUBLE_PRECISION, source, tag, MPI_COMM_WORLD, ierr)
            end if
            deallocate(R, Q)
        end do

    end subroutine worker_givens_qr
    ! {{{ solve_qr(A, b) result(x)
    function solve_qr(A, b) result(x)
        ! Solve a linear system using QR decomposition
        real(kind=wp) :: A(:,:), b(:), x(size(b))
        real(kind=wp) :: Q(size(A,1), size(A,2)), R(size(A,1), size(A,2))
        integer :: l, m, n
        n = size(A, 1)
        m = size(A, 2)
        if (n /= m) then
            print *, "n /= m, A is not square"
            stop
        end if
        l = size(b, 1)
        if (n /= l) then
            print *, "n /= l, A and b are of different sizes"
            stop
        end if

        call givens_qr(A, Q, R)

        ! Ax = QRx = b => Rx = Q'b
        Q = transpose(Q)
        b = matmul(Q, b)
        x = backsolve_r(R, b)
    end function solve_qr
    ! }}}
    ! {{{ backsolve_r(R, b) result(x)
    function backsolve_r(R, b) result(x)
        ! Backsolve right triangular systems as is common with QR
        ! decomposition
        real(kind=wp) :: R(:, :), b(:), x(size(b))
        integer :: i, j, n
        n = size(R, 1)

        ! Solve Q x = b
        do j = n, 1, -1
            x(j) = b(j) / R(j, j)
            !$OMP parallel do
            do i = j-1, 1, -1
                b(i) = b(i) - R(i, j)*x(j)
            end do
            !$OMP end parallel do
        end do
    end function backsolve_r
    ! }}}
    ! }}}
    ! {{{ LU Routines
    ! {{{ solve_lu(A, b) result(x)
    function solve_lu(A, b) result(x)
        ! Solve a linear system using LU decomposition
        real(kind=wp) :: A(:,:), b(:), x(size(b))
        integer :: l, m, n
        n = size(A, 1)
        m = size(A, 2)
        if (n /= m) then
            print *, "n /= m, A is not square"
            stop
        end if
        l = size(b, 1)
        if (n /= l) then
            print *, "n /= l, A and b are of different sizes"
            stop
        end if

        call decomp_lu(A)
        x = backsolve_lu(A, b)
    end function solve_lu
    ! }}}
    ! {{{ backsolve_lu(LU, b) result(x)
    function backsolve_lu(LU, b) result(x)
        ! Backsolve the two triangular systems as is common with LU
        ! decomposition
        real(kind=wp) :: LU(:, :), b(:), x(size(b))
        integer :: i, j, n
        n = size(LU, 1)

        ! Solve L y = b
        do j = 1, n
            x(j) = b(j)
            do i = j+1, n
                b(i) = b(i) - LU(i, j)*x(j)
            end do
        end do
        ! Now b is actually y...
        ! Solve U x = y
        do j = n, 1, -1
            x(j) = b(j) / LU(j, j)
            !$OMP parallel do
            do i = j-1, 1, -1
                b(i) = b(i) - LU(i, j)*x(j)
            end do
            !$OMP end parallel do
        end do
    end function backsolve_lu
    ! }}}
    ! {{{ decomp_lu(A)
    subroutine decomp_lu(A)
        real(kind=wp), intent(inout) :: A(:,:)
        integer :: i, j, k, m, n

        n = size(A, 1)
        m = size(A, 2)
        if (n /= m) then
            print *, "n /= m, A is not square"
            stop
        end if

        do k = 1, n-1
            ! Column normilization
            A(k+1:n, k) = A(k+1:n, k) / A(k, k)
            do i = k+1, n
                do j = k+1, n
                    A(i, j) = A(i, j) - A(i, k)*A(k, j)
                end do
            end do
        end do
    end subroutine decomp_lu
    ! }}}
    ! }}}
    subroutine print_matrix(A)
        real(kind=wp) :: A(:, :)
        integer :: i, n

        n = size(A,1)
        do i=1, n
            print *, A(i,:)
        end do
    end subroutine print_matrix
    subroutine init_random_seed()
        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed
        
        call random_seed(size=n)
        allocate(seed(n))

        call system_clock(count=clock)

        seed = clock + 37 * [(i -1, i=1, n)]

        call random_seed(put = seed)
        deallocate(seed)
    end subroutine init_random_seed
end module linearmethods
! vim: set foldmethod=marker :
