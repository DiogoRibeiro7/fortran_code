program LAPACKExample
    implicit none
    
    ! Declare variables
    integer, parameter :: n = 3, nrhs = 1, lda = n, ldb = n, info = 0
    integer :: ipiv(n)
    real*8 :: A(lda,n), b(ldb), work(n)
    
    ! Initialize the coefficient matrix A and the right-hand side vector b
    A = reshape([2.0, 1.0, -1.0, -3.0, -1.0, 2.0, -2.0, 1.0, 2.0], [lda, n])
    b = [8.0, -11.0, -3.0]
    
    ! Solve the system of linear equations using LAPACK routine dgesv
    call dgesv(n, nrhs, A, lda, ipiv, b, ldb, info)
    
    ! Check if the solution was successful
    if (info == 0) then
        print *, "Solution:"
        print *, b
    else
        print *, "Failed to solve the system."
    end if
end program LAPACKExample


