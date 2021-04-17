program main
	implicit none
	include 'mpif.h'
	! Переменные для MPI
	integer :: iErr         
    	integer :: rank
    	integer :: nProcess
    	integer, dimension(MPI_STATUS_SIZE) :: status
	
	!!! Переменные
	! A - исходная матрица коэф., B - правая часть, X - вектор решения, Ax(A*x), Rk(невязка)
	real, allocatable :: A(:,:), B(:), X(:), Ax(:), Rk(:), temp(:)
	integer :: length, i, j
	real :: timeStart, timeStop, tau, eps=0.000001, vectors_scalar_multiplication
	logical :: check_result

	!!! Body
	write (*, *) "Начало программы"
	
	! Выделяем память
	open(1,file="matrix.f90")
	read(1,*) length
	allocate(A(length,length))
	allocate(B(length))
	allocate(X(length))
	allocate(Ax(length))
	allocate(Rk(length))
	allocate(temp(length))
	
	! Считываем данные
	do i=1,length
		read(1,*)(A(i,j), j=1,length), B(i)
	end do

	! Вывод матрицы
	write(*,*) "Размерность матрицы = ", length
	write(*,*) ""
	
	call show_matrix(A, length, "Исходна матрица А:")	
	call show_vector(A, length, "Исходная матрица B:")
	
	call cpu_time(timeStart)
	
	call MPI_INIT(iErr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcess, iErr)
	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, iErr)
	
	! Задаем начальные значения
	do i=1,length
		! Преобладание диагональных элементов для достаточного условия сходимости
		X(i)=B(i)/A(i,i)
	end do
	
	! Вычисляем вектор невязок (Ax-b)
	call multiplication_matrix_vector(A, X, length, Ax)
	call vectors_difference(Ax, B, length, Rk)
	
	do while(check_result(Rk, length, eps))
		call multiplication_matrix_vector(A, Rk, length, Ax)
		tau = vectors_scalar_multiplication(Rk, Rk, length) / vectors_scalar_multiplication(Ax, Rk, length)
		call multiplication_num_vector(Rk, tau, length, temp)
		call multiplication_num_vector(Ax, tau, length, Ax)
		call vectors_difference(X, temp, length, X)
		call vectors_difference(Rk, Ax, length, Rk)
	end do
	call MPI_FINALIZE(iErr)
	
	! Результаты
	call show_vector(X, length, "Результирующий вектор X:")
		
	! Чистим память
	deallocate(A)
	deallocate(B)
	deallocate(X)
	deallocate(Ax)
	deallocate(Rk)
	deallocate(temp)

	call cpu_time(timeStop)

	write(*,*) "Время работы: ", timeStop - timeStart

	write(*,*) "Программа завершена"
end

! Разность векторов
subroutine vectors_difference(A, B, length, result_vector)
	implicit none
	integer :: i
	integer, intent(in) :: length	
	real, intent(in out) :: A(length), B(length), result_vector(length)
	do i=1,length
		result_vector(i) = A(i) - B(i)
	end do
end subroutine vectors_difference

! Умножение матрицы на вектор
subroutine multiplication_matrix_vector(matrix, vector, length, result_vector)
	implicit none
	integer :: i, j
	integer, intent(in) :: length	
	real, intent(in) :: matrix(length,length), vector(length)
	real, intent(out) :: result_vector(length)
	do i=1,length
		result_vector(i)=0
		do j=1,length
			result_vector(i) = result_vector(i) + matrix(i,j) * vector(j)
		end do
	end do
end subroutine multiplication_matrix_vector

! Умножение числа на вектор
subroutine multiplication_num_vector(vector, num, length, result_vector)
	implicit none
	integer :: i
	integer, intent(in) :: length	
	real, intent(in out) :: vector(length), num, result_vector(length)
	do i=1,length
		result_vector(i) = vector(i) * num
	end do
end subroutine multiplication_num_vector

! Вывод вектора
subroutine show_vector(vector, length, text)
	implicit none
	integer :: i
	integer, intent(in) :: length
	real, intent(in) :: vector(length)
	character(len = 48) :: text
	write(*,*) text
	do i=1,length
		write(*,20) i, vector(i)
	end do
	write(*,*) ""
	20 format(1x,i4,20f8.2)
end subroutine

! Вывод матрицы
subroutine show_matrix(matrix, length, text)
	implicit none
	integer :: i, j
	integer, intent(in) :: length
	real, intent(in) :: matrix(length, length)
	character(len = 40) :: text
	write(*,*) text
	do i=1,length
		write(*,20) i, (matrix(i,j),j=1,length)
	end do
	write(*,*) ""
	20 format(1x,i4,20f8.2)
end subroutine

! Скалярное произведение векторов
real function vectors_scalar_multiplication(A, B, length) result(res)
	implicit none
	integer, intent(in) :: length
	real, intent(in) :: A(length), B(length)
	integer :: i
	res = 0
	do i=1,length
		res = res + A(i) * B(i)
	end do
	return
end function vectors_scalar_multiplication

! Проверка сходимости
logical function check_result(A, length, eps) result(res)
	implicit none
	real, intent(in) :: A(length), eps
	integer, intent(in) :: length
	integer :: i
	do i=1,length
		res = ABS(A(i)) > eps
			if (.not. res) then
				exit
			end if
	end do
	return
end function check_result


