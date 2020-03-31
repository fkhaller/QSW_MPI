module linked_list

    use, intrinsic :: iso_fortran_env, only: dp => real64

    type :: element
        integer :: row_index
        integer :: col_index
        complex(dp) :: value
        type(element), pointer :: previous
        type(element), pointer :: next
    end type element

    type :: COO
        integer :: rows
        integer :: columns
        integer :: nnz = 0
        type(element), pointer :: current
        type(element), pointer :: first
        type(element), pointer :: last
    end type COO

    contains

    subroutine deallocate_coo(A)

        type(COO), intent(inout) :: A
        type(element), pointer :: current, next

        current => A%first
        next => A%first%next

        do while (associated(next))
            deallocate(current)
            current => next
            next => current%next
        enddo

        deallocate(current)

    end subroutine deallocate_coo

    subroutine initialize_list(row_index, col_index, value, list)

        integer, intent(in) :: row_index
        integer, intent(in) :: col_index
        complex(dp), intent(in) :: value
        type(element), pointer, intent(inout) :: list

        allocate(list)

        list = element( row_index, &
                        col_index, &
                        value, &
                        null(), &
                        null())

    end subroutine initialize_list

    subroutine append_element(row_index, col_index, value, elements)

        integer, intent(in) :: row_index
        integer, intent(in) :: col_index
        complex(dp), intent(in) :: value
        type(element), pointer, intent(inout) :: elements

        if (.not. associated(elements)) then
            call initialize_list(row_index, col_index, value, elements)
            return
        endif

        allocate(elements%next)

        elements%next = element(row_index, &
                                col_index, &
                                value, &
                                elements, &
                                null())

        elements => elements%next

    end subroutine append_element

    subroutine append_non_zero(row_index, col_index, value, A)

        integer, intent(in):: row_index
        integer, intent(in):: col_index
        complex(dp), intent(in) :: value
        type(COO), intent(inout) :: A


        if (A%nnz == 0) then
            call initialize_list(row_index, col_index, value, A%current)
            A%nnz = 1
            A%first => A%current
            A%last => A%current
            return
        endif

        A%current => A%last

        allocate(A%current%next)

        A%current%next = element(   row_index, &
                                    col_index, &
                                    value, &
                                    A%current, &
                                    null())


        A%nnz = A%nnz + 1
        A%current => A%current%next
        A%last => A%current

    end subroutine append_non_zero

    subroutine delete_non_zero(A)

        type(COO), intent(inout) :: A
        type(element), pointer :: non_zero

        non_zero => A%current

        A%current%previous%next => A%current%next
        A%current%next%previous => A%current%previous

        A%current => A%current%next

        deallocate(non_zero)

        A%nnz = A%nnz - 1

    end subroutine delete_non_zero

    subroutine print_list(A)

        type(COO), intent(inout) :: A
        integer :: i, n

        A%current => A%first
        n = 1

        do i = 1, A%nnz
            write(*,*) n, A%current%row_index, A%current%col_index, A%current%value
            A%current => A%current%next
            n = n + 1
        enddo

    end subroutine print_list

    subroutine sum_coo(A, B, C)

        type(coo), intent(inout) :: A
        type(coo), intent(inout) :: B
        type(coo), intent(out) :: C

        A%current => A%first
        B%current => B%first

        do while(associated(A%current) .and. associated(B%current))
            if ((A%current%col_index == B%current%col_index) .and. &
                (A%current%row_index == B%current%row_index)) then
                call append_non_zero(A%current%row_index, A%current%col_index, A%current%value + B%current%value, C)
                !write(*,*) "add", A%current%row_index, B%current%row_index
                A%current => A%current%next
                B%current => B%current%next
            elseif ((A%current%col_index < B%current%col_index) .and. &
                (A%current%row_index == B%current%row_index)) then
                call append_non_zero(A%current%row_index, A%current%col_index, A%current%value, C)
                !write(*,*) "A", A%current%row_index, B%current%row_index
                A%current => A%current%next
            elseif ((A%current%col_index > B%current%col_index) .and. &
                (A%current%row_index == B%current%row_index)) then
                call append_non_zero(B%current%row_index, B%current%col_index, B%current%value, C)
                !write(*,*) "B", A%current%row_index, B%current%row_index
                B%current => B%current%next
            elseif (A%current%row_index < B%current%row_index) then
                do while (A%current%row_index < B%current%row_index)
                    call append_non_zero(A%current%row_index, A%current%col_index, A%current%value, C)
                    !write(*,*) "C", A%current%row_index, B%current%row_index
                    A%current => A%current%next
                    if (.not. associated(A%current)) exit
                enddo
            elseif (A%current%row_index > B%current%row_index) then
                do while (A%current%row_index > B%current%row_index)
                    call append_non_zero(B%current%row_index, B%current%col_index, B%current%value, C)
                    !write(*,*) "D", A%current%row_index, B%current%row_index
                    B%current => B%current%next
                    if (.not. associated(B%current)) exit
                enddo
            endif
            if (.not. associated(B%current)) then
                do while (associated(A%current))
                    call append_non_zero(A%current%row_index, A%current%col_index, A%current%value, C)
                    !write(*,*) "E", A%current%row_index, B%current%row_index
                    A%current => A%current%next
                enddo
            elseif (.not. associated(A%current)) then
                do while (associated(B%current))
                    call append_non_zero(B%current%row_index, B%current%col_index, B%current%value, C)
                    !write(*,*) "F", A%current%row_index, B%current%row_index
                    B%current => B%current%next
                enddo
            endif
        enddo

    end subroutine sum_coo

end module linked_list

module utilities

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use :: sparse

    implicit none
    contains

    subroutine prefix_sum(array)

        integer, dimension(:), intent(inout) :: array
        integer :: i

        do i = 2, size(array,1)
            array(i) = array(i) + array(i-1)
        enddo

    end subroutine prefix_sum

    subroutine transpose_csr(A, A_D)

        type(CSR), intent(in) :: A
        type(CSR), intent(out) :: A_D

        integer, dimension(:), allocatable :: A_D_row_indexes
        integer :: nnz

        integer :: A_D_row_ind, A_D_col_ind
        complex(dp) val

        integer :: i, j

        nnz = size(A%col_indexes,1)

        A_D%rows = A%rows
        A_D%columns = A%columns

        allocate(A_D%row_starts(A%rows + 1))
        allocate(A_D%col_indexes(nnz))
        allocate(A_D%values(nnz))

        A_D%row_starts = 0
        A_D%row_starts(1) = 1
        A_D%col_indexes = 0
        A_D%values = A%values

        allocate(A_D_row_indexes(nnz))

        A_D_row_indexes = A%col_indexes

        do i = 1, A%rows
            do j = A%row_starts(i), A%row_starts(i + 1) - 1
                A_D%col_indexes(j) = i
                A_D%row_starts(A%col_indexes(j) + 1) = A_D%row_starts(A%col_indexes(j) + 1) + 1
            enddo
        enddo

        call prefix_sum(A_D%row_starts)

        do i = 2, nnz

            A_D_row_ind = A_D_row_indexes(i)
            A_D_col_ind = A_D%col_indexes(i)
            val = A_D%values(i)

            j = i - 1

            do while (j >= 1)

                if (A_D_row_indexes(j) <= A_D_row_ind) exit
                    A_D_row_indexes(j + 1) = A_D_row_indexes(j)
                    A_D%col_indexes(j + 1) = A_D%col_indexes(j)
                    A_D%values(j + 1) = A_D%values(j)
                    j = j - 1
            enddo

            A_D_row_indexes(j + 1) = A_D_row_ind
            A_D%col_indexes(j + 1) = A_D_col_ind
            A_D%values(j + 1) = val

        enddo

        deallocate(A_D_row_indexes)

    end subroutine transpose_csr

end module utilities

program test

    use linked_list
    use sparse
    use utilities

    implicit none

    type(CSR) :: H, L
    integer :: i, j, k, ll

    type(COO) :: HSO1, HSO2, LSO1, LSO2, LSO3

    type(COO) :: SO

    !LSO1
    integer :: row_block, col_block, L_row, L_nz, row, row_nz

    integer :: cnt, L_check
    logical :: hit, last, lapped
    integer, dimension(:), allocatable :: temp_row_starts

    !LSO2
    type(CSR) :: L_T
    type(COO) :: LSO2_elements
    complex(dp) :: value
    integer :: row_index, col_index

    !LSO3
    type(COO) :: LSO3_elements
    type(element), pointer :: row_elements

    type(COO) :: SO_temp_1, SO_temp_2

    H%rows = 4
    H%columns = 4

    allocate(H%row_starts(H%rows + 1))
    allocate(H%col_indexes(12))
    allocate(H%values(12))

    H%row_starts = [1,4,7,10,13]
    H%col_indexes = [1,2,4,1,2,3,2,3,4,1,3,4]
    H%values = [2,-1,-1,-1,2,-1,-1,2,-1,-1,-1,2]

    L%rows = 4
    L%columns = 4

    allocate(L%row_starts(L%rows + 1))
    allocate(L%col_indexes(12))
    allocate(L%values(12))

    L%row_starts = [1,4,7,10,13]
    L%col_indexes = [1,2,4,1,2,3,2,3,4,1,3,4]
    L%values =  [1.8139880718188666, 0.4023998995620808, 0.7235171145993213, &
                0.17664548926530443, 1.1378510756236193, 0.38146933019585005,&
                0.3584470100181627, 0.46245073679237114, 0.7197312284825397, &
                0.9743159835366508, 0.39866263852355677, 1.2538626845607173]

    ! LSO1

    LSO1%rows = L%rows**2
    LSO1%columns = L%columns**2

    allocate(temp_row_starts(LSO1%rows + 1))
    temp_row_starts(1) = 1
    temp_row_starts = 0
    do i = 1, LSO1%rows
        L_row = i - (ceiling(i/float(L%rows)) - 1)*L%rows
        L_nz = L%row_starts(L_row + 1) - L%row_starts(L_row)
        do j = 1, L%rows
            do k = L%row_starts(L_row), L%row_starts(L_row + 1) - 1
                call append_non_zero(i, (j - 1)*L%rows + L%col_indexes(k), L%values(k), LSO1)
                temp_row_starts(i + 1) = temp_row_starts(i  + 1) + 1
            enddo
        enddo
    enddo

    call prefix_sum(temp_row_starts)

    LSO1%current => LSO1%first
    col_block = ceiling(LSO1%current%col_index/float(L%columns))

    do while (associated(LSO1%current))

        row_block = ceiling(LSO1%current%row_index/float(L%rows))
        row_nz = temp_row_starts(LSO1%current%row_index + 1) &
                - temp_row_starts(LSO1%current%row_index)
        i = 0
        j = L%row_starts(row_block)
        hit = .false.
        lapped = .false.

        do while(i < row_nz)

            do while ((L%col_indexes(j) == col_block).and.(i < row_nz))
                LSO1%current%value = LSO1%current%value * conjg(L%values(j))
                if (associated(LSO1%current%next)) then
                    LSO1%current => LSO1%current%next
                    col_block = ceiling(LSO1%current%col_index/float(L%columns))
                else
                    LSO1%current => LSO1%current%next
                    exit
                endif
                hit = .true.
                i = i + 1
            enddo

            do while ((L%col_indexes(j) > col_block) .and. (i < row_nz))
                call delete_non_zero(LSO1)
                if (associated(LSO1%current)) then
                    col_block = ceiling(LSO1%current%col_index/float(L%columns))
                else
                    LSO1%current => LSO1%current%next
                    exit
                endif
                hit = .true.
                i = i + 1
            enddo

            do while ((L%row_starts(row_block + 1) < j) .and. (i < row_nz))
                call delete_non_zero(LSO1)
                if (associated(LSO1%current)) then
                    col_block = ceiling(LSO1%current%col_index/float(L%columns))
                else
                    LSO1%current => LSO1%current%next
                    exit
                endif
                hit = .true.
                i = i + 1
            enddo

            if (.not. hit) then
                j = j + 1
            endif

            hit = .false.
            if (.not. associated(LSO1%current)) exit

        enddo
        write(*,*) j
    enddo

    deallocate(temp_row_starts)

    call print_list(LSO1)
    write(*,*)

    open (unit = 10, file = "LSO1")

    LSO1%current => LSO1%first

    do while (associated(LSO1%current))
        write(10,*) LSO1%current%col_index, LSO1%current%row_index, real(LSO1%current%value)
        LSO1%current => LSO1%current%next
    enddo

    close (10)

    ! LSO2

    call transpose_csr(L, L_T)

    do i = 1, L%rows
        do j = 1, L%rows
            value = 0
            do k = L_T%row_starts(i), L_T%row_starts(i + 1) - 1
                do ll = L_T%row_starts(j), L_T%row_starts(j + 1) - 1
                    if(L_T%col_indexes(ll) == L_T%col_indexes(k)) then
                        value = value + conjg(L_T%values(k))*L_T%values(ll)
                    endif
                    if(L_T%col_indexes(ll) > L_T%col_indexes(k)) exit
                enddo
            enddo
            call append_non_zero(i,j,value,LSO2_elements)
        enddo
    enddo

    do i = 1, L%rows
        LSO2_elements%current => LSO2_elements%first
        do while (associated(LSO2_elements%current))
            row_index = (i-1)*L%rows + LSO2_elements%current%row_index
            col_index = (i-1)*L%columns + LSO2_elements%current%col_index
            call append_non_zero(row_index, col_index, &
                -0.5d0 * LSO2_elements%current%value, LSO2)
           LSO2_elements%current => LSO2_elements%current%next
        enddo
    enddo

    call deallocate_coo(LSO2_elements)

    !call print_list(LSO2)
    write(*,*)

    !LSO3

    do i = 1, L%rows
        do j = 1, L%rows
            value = 0
            do k = L_T%row_starts(i), L_T%row_starts(i + 1) - 1
                do ll = L_T%row_starts(j), L_T%row_starts(j + 1) - 1
                    if(L_T%col_indexes(ll) == L_T%col_indexes(k)) then
                        value = value + L_T%values(k)*conjg(L_T%values(ll))
                    endif
                    if(L_T%col_indexes(ll) > L_T%col_indexes(k)) exit
                enddo
            enddo
            call append_non_zero(i,j,value,LSO3_elements)
        enddo
    enddo

    LSO3_elements%current => LSO3_elements%first

    do i = 1, L%rows
        cnt = 0
        row_elements => LSO3_elements%current
        do while(associated(LSO3_elements%current) .and. (LSO3_elements%current%row_index == i))
            row_index = (LSO3_elements%current%row_index - 1) * L%rows + 1
            col_index = (LSO3_elements%current%col_index - 1) * L%columns + 1
            value = -0.5d0 * LSO3_elements%current%value
            call append_non_zero(row_index, col_index, value, LSO3)
            write(*,*) row_index, col_index, -0.5d0 * LSO3_elements%current%value
            cnt = cnt + 1
            LSO3_elements%current => LSO3_elements%current%next
        enddo
        if (cnt /= 0) then
            do j = 2, L%rows
            LSO3_elements%current => row_elements
            row_index = (LSO3_elements%current%row_index - 1) * L%rows + j
                do k = 1, cnt
                    col_index = (LSO3_elements%current%col_index - 1) * L%columns + j
                    value = -0.5d0 * LSO3_elements%current%value
                    call append_non_zero(row_index, col_index, value, LSO3)
                    LSO3_elements%current => LSO3_elements%current%next
                enddo
            enddo
        endif
    enddo

    call deallocate_coo(LSO3_elements)

    call print_list(LSO3)
    write(*,*)

    call sum_coo(LSO1, LSO2, SO_temp_1)
    call print_list(SO_temp_1)
    call sum_coo(LSO3, SO_temp_1, SO_temp_2)
    write(*,*)

    call print_list(SO_temp_2)

    ! HSO1



end program test
