!DOC (https://www.sgu.ru/sites/default/files/textdocsfiles/2014/12/29/lin.alg.pdf), стр. 9

program main
	implicit none
	!!! Переменные
	! A - исходная матрица коэф., B - правая часть, X - вектор решения, Ax(A*x), Rk(невязка)
	real, allocatable :: A(:,:), B(:), X(:), Ax(:), Rk(:), temp(:)
	integer :: n, i, j, iter
	! SM(Rk,Rk), SMD(A*Rk, Rk) - скалярные произведения
	real :: timeStart, timeStop, tau, SM, SMD, eps=0.001, vv_mult
	logical :: check

	!!! Body
	write (*, *) "Начало программы"
	
	call cpu_time(timeStart)
	
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
	
	call show_m(A, n, "Исходна матрица А:")	
	call show_v(A, n, "Исходная матрица B:")
	
	! Задаем начальные значения
	do i=1,n
		! Преобладание диагональных элементов для достаточного условия сходимости
		X(i)=B(i)/A(i,i)
	end do
	iter = 0
	
	! Вычисляем вектор невязок (Ax-b)
	call mv_mult(A, X, n, Ax)
	call vv_substr(Ax, B, n, Rk)
	
	do while(check(Rk, n, eps))
		! Вычисляем произведение частного скалярных величин на вектор
		call mv_mult(A, Rk, n, Ax)
		SM = vv_mult(Rk, Rk, n)
		SMD = vv_mult(Ax, Rk, n)
		tau = SM / SMD
		call vn_mult(Rk, tau, n, temp)
		call vn_mult(Ax, tau, n, Ax)
		call vv_substr(X, temp, n, X)
		call vv_substr(Rk, Ax, n, Rk)
		iter = iter + 1
	end do
	
	! Результаты
	call show_v(X, n, "Результирующий вектор X:")
		
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
subroutine vv_substr(A, B, n, C)
	implicit none
	integer :: i
	integer, intent(in) :: n	
	real, intent(in out) :: A(n), B(n), C(n)
	do i=1,n
		C(i) = A(i) - B(i)
	end do
end subroutine vv_substr

! Умножение матрицы на вектор
subroutine mv_mult(A, B, n, C)
	implicit none
	integer :: i, j
	integer, intent(in) :: n	
	real, intent(in) :: A(n,n), B(n)
	real, intent(out) :: C(n)
	do i=1,n
		C(i)=0
		do j=1,n
			C(i) = C(i) + A(i,j) * B(j)
		end do
	end do
end subroutine mv_mult

! Умножение числа на вектор
subroutine vn_mult(A, num, n, B)
	implicit none
	integer :: i
	integer, intent(in) :: n	
	real, intent(in out) :: A(n), num, B(n)
	do i=1,n
		B(i) = A(i) * num
	end do
end subroutine vn_mult

! Вывод вектора
subroutine show_v(A, n, text)
	implicit none
	integer :: i
	integer, intent(in) :: n
	real, intent(in) :: A(n)
	character(len = 48) :: text
	write(*,*) text
	do i=1,n
		write(*,20) i, A(i)
	end do
	write(*,*) ""
	20 format(1x,i4,20f8.2)
end subroutine

! Вывод матрицы
subroutine show_m(A, n, text1)
	implicit none
	integer :: i, j
	integer, intent(in) :: n
	real, intent(in) :: A(n, n)
	character(len = 40) :: text1
	write(*,*) text1
	do i=1,n
		write(*,20) i, (A(i,j),j=1,n)
	end do
	write(*,*) ""
	20 format(1x,i4,20f8.2)
end subroutine

! Скалярное произведение векторов
real function vv_mult(A, B, n) result(res)
	implicit none
	integer, intent(in) :: n
	real, intent(in) :: A(n), B(n)
	integer :: i
	res = 0
	do i=1,n
		res = res + A(i) * B(i)
	end do
	return
end function vv_mult

! Проверка сходимости
logical function check(A, n, eps) result(res)
	implicit none
	real, intent(in) :: A(n), eps
	integer, intent(in) :: n
	integer :: i
	do i=1,n
		res = A(i)<eps
			if (.not. res) then
				exit
			end if
	end do
	return
end function check


