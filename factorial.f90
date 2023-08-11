function Factorial(n) result(fact)
    integer, intent(in) :: n
    real(kind=8) :: fact
    integer :: i
    
    fact = 1.0d0  ! Using double precision literal
    do i = 1, n
        fact = fact * i
    end do
end function Factorial

