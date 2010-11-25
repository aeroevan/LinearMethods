module prob9
    implicit none
contains
    subroutine givens_rotation(A, i, j, Q)
        double precision, intent(inout) :: A(:,:)
        integer, intent(in) :: i, j
        double precision, intent(out) :: Q(size(A,1),size(A,2))
        integer :: n, k
        double precision :: invdenom, ctheta, stheta, rowi(size(A,1)), rowj(size(A,1))

        n = size(A,1)

        invdenom = 1.0/sqrt(A(i,i)*A(i,i) + A(j,i)*A(j,i))
        ctheta = A(i,i)*invdenom
        stheta = A(j,i)*invdenom

        Q = 0.
        do k=1, n
            Q(k,k) = 1.0
        end do
        Q(i,i) = ctheta
        Q(j,j) = ctheta
        Q(i,j) = stheta
        Q(j,i) = -stheta
        A = matmul(Q,A)
        !A(j,i) = 0.0d0
    end subroutine givens_rotation
    subroutine givens_qr(A, Q, R)
        double precision, intent(in) :: A(:,:)                    ! R is triangular...
        double precision, intent(out) :: Q(size(A,1), size(A,2)), R(size(A,1), size(A,2))
        double precision :: Qtmp(size(A,1), size(A,2))
        integer :: i, n, j, k, T
        logical :: updated

        n = size(A,1)

        ! Initialize Q
        Q = 0.0d0
        do k=1, n
            Q(k,k) = 1.0d0
        end do

        ! Initialize R
        R = A

        ! I'm sure there's a better way to handle the schedule... but this works
        T = 1
        updated = .true.
        do while(updated)
            updated = .false.
            do j=1, n-1
                do k=n, j, -1
                    if (k-2*j==n-1-T) then
                        call givens_rotation(R, j, k, Qtmp)
                        Q = matmul(Qtmp, Q)
                        updated = .true.
                    end if
                end do
            end do
            T = T + 1
        end do
        Q = transpose(Q)
    end subroutine givens_qr
    subroutine print_matrix(A)
        double precision :: A(:, :)
        integer :: i, n

        n = size(A,1)
        do i=1, n
            print *, A(i,:)
        end do
    end subroutine print_matrix
end module prob9

!program prob9main
!    use prob9
!    implicit none
!    double precision :: A(3,3), Q(3,3), R(3,3)
!
!    A(1,1) = 12.
!    A(1,2) = -51.
!    A(1,3) = 4.
!
!    A(2,1) = 6.
!    A(2,2) = 167.
!    A(2,3) = -68.
!
!    A(3,1) = -4.
!    A(3,2) = 24.
!    A(3,3) = -41.
!    print *, "A:"
!    call print_matrix(A)
!    call givens_qr(A, Q, R)
!    print *, "Q:"
!    call print_matrix(Q)
!    print *, "R:"
!    call print_matrix(R)
!    A = matmul(Q,R)
!    print *, "A:"
!    call print_matrix(A)
!
!end program prob9main

! Output of test program prob9main:
! A:
!   12.000000      -51.000000       4.0000000    
!   6.0000000       167.00000      -68.000000    
!  -4.0000000       24.000000      -41.000000    
! Q:
!  0.85714281     -0.39428577      0.33142859    
!  0.42857146      0.90285718     -3.42857055E-02
! -0.28571427      0.17142859      0.94285715    
! R:
!   14.000000       21.000008      -14.000004    
! -4.76493739E-07   175.00002      -70.000008    
!  1.80946937E-08   0.0000000      -35.000000    
! A:
!   11.999999      -51.000008       4.0000029    
!   6.0000000       167.00002      -68.000015    
!  -3.9999998       24.000004      -41.000000
