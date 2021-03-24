program main
	implicit none
	! Variables
	real, allocatable :: A(:,:), B(:), X(:)
	integer :: n, i, j
	real :: timeStart, timeStop

	! Body
	write (*, *) "Начало программы"
	
	! Выделяем память и считываем данные
	open(1,file="matrix.f90")
	read(1,*) n
	allocate(A(n,n))
	allocate(B(n))
	allocate(X(n))
	do i=1,n
		read(1,*)(A(i,j), j=1,n), B(i)
	enddo

	! Вывод матрицы
	write(*,*) "Размерность матрицы = ", n
	write(*,*) ""
	write(*,*) "Исходная матрица A:"
	do i=1,n
		write(*,21) i, (A(i,j),j=1,n)
	enddo
	write(*,*) ""
	write(*,*) "Исходная матрица B:"
	do i=1,n
		write(*,21) i, B(i)
	enddo
	write(*,*) ""

	21 format(1x,i4,20f8.2)

	! Чистим память
	deallocate(A)
	deallocate(B)
	deallocate(X)

	write(*,*) "Программа завершена"
end

subroutine cpu_time(time)
	
end
