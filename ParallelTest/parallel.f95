program parallel
    implicit none
    !$OMP PARALLEL num_threads(10)
        write(*,*) "Hello"
    !$OMP END PARALLEL
end program parallel
