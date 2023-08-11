function Factorial(n) result(fact)
    integer, intent(in) :: n
    integer :: fact, i
    
    fact = 1
    i = 1
    do while (i <= n)
        fact = fact * i
        i = i + 1
    end do
end function Factorial
