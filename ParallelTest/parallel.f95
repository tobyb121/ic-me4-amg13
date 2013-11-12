program parallel
    use sparse
    implicit none
    type(sparse_matrix) sA
    sA%columns=10
    sA%rows=10
    call initialise_sparse(sA,28)
    !$OMP PARALLEL num_threads(10)
        write(*,*) "Hello"
    !$OMP END PARALLEL
    call free_matrix(sA)
end program parallel
