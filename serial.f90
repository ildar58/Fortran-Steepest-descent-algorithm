program main
	use procedures

	implicit none
	!!! Переменные
	! A - исходная матрица коэф., B - правая часть, X - вектор решения, Ax(A*x), Rk(невязка)
	integer :: length = 4
	integer :: i, j
	real*8, allocatable :: A(:, :), B(:), X(:), Ax(:), Rk(:), temp(:)
	real*8 :: timeStart, timeStop, tau, eps=0.000001

	!!! Body
	write (*, *) "Начало программы"
	
	!! Выделяем память
	allocate(A(length,length))
	allocate(B(length))
	allocate(X(length))
	allocate(Ax(length))
	allocate(Rk(length))
	allocate(temp(length))
	
	! Генерируем данные
	call fill_matrix(A, length)
	call fill_vector(B, length)
	
	! Вывод матрицы
	write(*,*) "Размерность матрицы = ", length
	write(*,*) ""
	
	write(*,*) "Матрица A:"
	call show_matrix(A)
	write(*,*) "Матрица B:"	
	call show_vector(B)
	
	call cpu_time(timeStart)
	
	! Задаем начальные значения
	do i=1,length
		! Преобладание диагональных элементов для достаточного условия сходимости
		X(i)=B(i) / A(i,i)
	end do
	
	! Вычисляем вектор невязок (Ax-b)
	call multiplication_matrix_vector(A, X, Ax)
	call vectors_difference(Ax, B, Rk)
	
	do while(check_result(Rk, eps))
		call multiplication_matrix_vector(A, Rk, Ax)
		tau = vectors_scalar_multiplication(Rk, Rk) / vectors_scalar_multiplication(Ax, Rk)
		call multiplication_num_vector(Rk, tau, temp)
		call multiplication_num_vector(Ax, tau, Ax)
		call vectors_difference(X, temp, X)
		call vectors_difference(Rk, Ax, Rk)
	end do
	
	! Результаты
	write(*,*) "Матрица X:"
	call show_vector(X)
		
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


