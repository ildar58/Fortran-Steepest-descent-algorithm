program main
	implicit none
	!!! Переменные
	! A - исходная матрица коэф., B - правая часть, X - вектор решения, Ax(A*x), Rk(невязка)
	real, allocatable :: A(:,:), B(:), X(:), Ax(:), Rk(:), temp(:)
	integer :: n, i, j, iter
	! SM(Rk,Rk), SMD(A*Rk, Rk) - скалярные произведения
	real :: timeStart, timeStop, tau, eps=0.000001, vectors_scalar_multiplication
	logical :: check_result

	!!! Body
	write (*, *) "Начало программы"
	
	! Выделяем память
	open(1,file="matrix.f90")
	read(1,*) n
	allocate(A(n,n))
	allocate(B(n))
	allocate(X(n))
	allocate(Ax(n))
	allocate(Rk(n))
	allocate(temp(n))
	
	! Считываем данные
	do i=1,n
		read(1,*)(A(i,j), j=1,n), B(i)
	end do

	! Вывод матрицы
	write(*,*) "Размерность матрицы = ", n
	write(*,*) ""
	
	call show_matrix(A, n, "Исходна матрица А:")	
	call show_vector(A, n, "Исходная матрица B:")
	
	call cpu_time(timeStart)
	
	! Задаем начальные значения
	do i=1,n
		! Преобладание диагональных элементов для достаточного условия сходимости
		X(i)=B(i)/A(i,i)
	end do
	iter = 0
	
	! Вычисляем вектор невязок (Ax-b)
	call multiplication_matrix_vector(A, X, n, Ax)
	call vectors_difference(Ax, B, n, Rk)
	
	do while(check_result(Rk, n, eps))
		call multiplication_matrix_vector(A, Rk, n, Ax)
		tau = vectors_scalar_multiplication(Rk, Rk, n) / vectors_scalar_multiplication(Ax, Rk, n)
		call multiplication_num_vector(Rk, tau, n, temp)
		call multiplication_num_vector(Ax, tau, n, Ax)
		call vectors_difference(X, temp, n, X)
		call vectors_difference(Rk, Ax, n, Rk)
	end do
	
	! Результаты
	call show_vector(X, n, "Результирующий вектор X:")
		
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
subroutine vectors_difference(A, B, n, result_vector)
	implicit none
	integer :: i
	integer, intent(in) :: n	
	real, intent(in out) :: A(n), B(n), result_vector(n)
	do i=1,n
		result_vector(i) = A(i) - B(i)
	end do
end subroutine vectors_difference

! Умножение матрицы на вектор
subroutine multiplication_matrix_vector(matrix, vector, n, result_vector)
	implicit none
	integer :: i, j
	integer, intent(in) :: n	
	real, intent(in) :: matrix(n,n), vector(n)
	real, intent(out) :: result_vector(n)
	do i=1,n
		result_vector(i)=0
		do j=1,n
			result_vector(i) = result_vector(i) + matrix(i,j) * vector(j)
		end do
	end do
end subroutine multiplication_matrix_vector

! Умножение числа на вектор
subroutine multiplication_num_vector(vector, num, n, result_vector)
	implicit none
	integer :: i
	integer, intent(in) :: n	
	real, intent(in out) :: vector(n), num, result_vector(n)
	do i=1,n
		result_vector(i) = vector(i) * num
	end do
end subroutine multiplication_num_vector

! Вывод вектора
subroutine show_vector(vector, n, text)
	implicit none
	integer :: i
	integer, intent(in) :: n
	real, intent(in) :: vector(n)
	character(len = 48) :: text
	write(*,*) text
	do i=1,n
		write(*,20) i, vector(i)
	end do
	write(*,*) ""
	20 format(1x,i4,20f8.2)
end subroutine

! Вывод матрицы
subroutine show_matrix(matrix, n, text)
	implicit none
	integer :: i, j
	integer, intent(in) :: n
	real, intent(in) :: matrix(n, n)
	character(len = 40) :: text
	write(*,*) text
	do i=1,n
		write(*,20) i, (matrix(i,j),j=1,n)
	end do
	write(*,*) ""
	20 format(1x,i4,20f8.2)
end subroutine

! Скалярное произведение векторов
real function vectors_scalar_multiplication(A, B, n) result(res)
	implicit none
	integer, intent(in) :: n
	real, intent(in) :: A(n), B(n)
	integer :: i
	res = 0
	do i=1,n
		res = res + A(i) * B(i)
	end do
	return
end function vectors_scalar_multiplication

! Проверка сходимости
logical function check_result(A, n, eps) result(res)
	implicit none
	real, intent(in) :: A(n), eps
	integer, intent(in) :: n
	integer :: i
	do i=1,n
		res = ABS(A(i)) > eps
			if (.not. res) then
				exit
			end if
	end do
	return
end function check_result


