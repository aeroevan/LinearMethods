program prob10
    use linearmethods
    implicit none
    integer :: Ai(6,6), bi(6)
    real(kind=wp) :: Ar(6,6), br(6), xr(6)

    integer :: i, j

    ! Initialize A and b first way...
    do i=1, 6
        do j=1, 6
            Ai(i,j) = i**j
        end do
    end do
    Ar = Ai
    bi(1) = 1
    do i=2, 6
        bi(i) = bi(i-1) + i**5
    end do
    br = bi

    xr = solve_qr(Ar, br)
    print *, xr
    ! Solution:
    !    33.453449      -23.775845       3.2380152       3.4227219     -0.59532136      0.26479077
    ! Initialize A and b first way again
    do i=1, 6
        do j=1, 6
            Ai(i,j) = i**j
        end do
    end do
    bi(1) = 1
    do i=2, 6
        bi(i) = bi(i-1) + i**5
    end do
    ! Reverse rows...
    do i = 1, 6
        Ar(i,:) = Ai(7-i,:)
        br(i) = bi(7-i)
    end do
    xr = solve_qr(Ar, br)
    print *, xr

contains
end program prob10
