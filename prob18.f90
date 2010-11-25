program prob18
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
    Ar = real(Ai)
    bi(1) = 1
    do i=2, 6
        bi(i) = bi(i-1) + i**5
    end do
    br = real(bi)

    xr = solve_lu(Ar, br)
    print *, xr
    ! Solution:
    !   1.49011612E-07 -8.33334923E-02   0.0000000      0.41666666      0.50000000      0.16666667    

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

    xr = solve_lu(Ar, br)
    print *, xr
    ! Solution:
    !   1.89185143E-04 -8.33648667E-02   0.0000000      0.41666666      0.50000000      0.16666667    
end program prob18
