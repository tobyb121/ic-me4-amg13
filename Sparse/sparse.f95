module sparse
    implicit none

    type sparse_matrix
        integer :: rows
        integer :: columns
        integer,private :: nnz
        integer,private :: nnz_max
        integer,  allocatable :: rowIndexes(:)
        integer,  allocatable :: columnIndexes(:)
        real, allocatable :: values(:)
    end type sparse_matrix
contains
    subroutine initialise_sparse(A,nnz_max)
        type(sparse_matrix),intent(inout) :: A
        integer,intent(in) :: nnz_max
        integer :: i
        A%nnz_max=nnz_max
        A%nnz=0
        allocate(A%columnIndexes(A%columns))
        allocate(A%rowIndexes(nnz_max))
        allocate(A%values(nnz_max))
        do i=1,A%columns
            A%columnIndexes(i)=i
        end do
        do i=1,nnz_max
            A%rowIndexes(i)=1
        end do
        do i=1,nnz_max
            A%values(i)=0
        end do
    end subroutine initialise_sparse

    subroutine free_matrix(A)
        type(sparse_matrix),intent(inout) :: A
        deallocate(A%columnIndexes)
        deallocate(A%rowIndexes)
        deallocate(A%values)
    end subroutine free_matrix

    logical function contains(A,i,j)
        type(sparse_matrix),intent(in) :: A
        integer,intent(in) :: i
        integer,intent(in) :: j
        integer :: n,max_n
        logical :: result = .false.
        if (j<A%columns) then
            max_n=A%columnIndexes(j+1)-1
        else
            max_n=A%nnz
        endif
        do n=A%columnIndexes(j),max_n
            if (A%rowIndexes(n)==i) then
                result=.true.
            endif
        end do
        contains=result
    end function contains

    logical function assign(A,i,j,value)
        type(sparse_matrix),intent(inout) :: A
        integer,intent(in) :: i
        integer,intent(in) :: j
        real,intent(in) :: value
        integer :: n,max_n,m
        assign=.false.
        if (j<A%columns) then
            max_n=A%columnIndexes(j+1)-1
        else
            max_n=A%nnz
        endif
        if(contains(A,i,j)) then
            do n=A%columnIndexes(j),max_n
                if (A%rowIndexes(n)==i) then
                    A%values(n)=value
                endif
            end do
            assign=.true.
        else if(A%nnz<A%nnz_max) then
            do n=A%columnIndexes(j),max_n
                if(A%rowIndexes(n)>i) then
                    exit
                endif
            end do
            do m = A%nnz,n,-1
                        A%values(m+1)=A%values(m)
                        A%rowIndexes(m+1)=A%rowIndexes(m)
                    end do
                    do m=j+1,A%columns
                        A%columnIndexes(m)=A%columnIndexes(m+1)
                    end do
                    A%rowIndexes(n)=i
                    A%values(n)=value
        endif

    end function assign

end module sparse

