	program main
		use procedures
		implicit none
		include 'mpif.h'
		
		!!! Переменные
		! A - исходная матрица коэф., B - правая часть, X - вектор решения, Ax(A*x), Rk(невязка)
		integer :: length = 20, lenghtThread
		integer :: i, j, iErr, nProcess, rank, status(MPI_STATUS_SIZE)
		real*8, allocatable :: A(:, :), B(:), X(:), Ax(:), Rk(:), temp(:), Athread(:,:), Bthread(:), Rkthread(:), Axthread(:)
		real*8 :: timeStart, timeStop, tau, eps=0.000001

		call MPI_INIT(iErr)
		call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcess, iErr)
		call MPI_COMM_RANK(MPI_COMM_WORLD, rank, iErr)
		

		!! Выделяем память
		allocate(A(length,length))
		allocate(Athread(length/4,length))
		allocate(X(length))
		allocate(B(length))
		allocate(Bthread(length/4))
		allocate(Ax(length))
		allocate(Axthread(length/4))
		allocate(Rk(length))
		allocate(Rkthread(length/4))
		allocate(temp(length/4))
	! Генерируем данные
		call fill_matrix(A, length)
		call fill_vector(B, length)
		do i=1,length
			X(i)=(0)
		end do

		call cpu_time(timeStart)
		if (rank == 0) then
			Athread = A(1:length/4, 1:length)
			call MPI_SEND(A(length/4+1:length/2,1:length), (length/4)*length*2, MPI_REAL, 1, 0, MPI_COMM_WORLD, iErr)
			call MPI_SEND(A(length/2+1:(length/4)*3,1:length), (length/4)*length*2, MPI_REAL, 2, 0, MPI_COMM_WORLD, iErr)
			call MPI_SEND(A(length/4*3+1:length,1:length), (length/4)*length*2, MPI_REAL, 3, 0, MPI_COMM_WORLD, iErr)

			Bthread = B(1:length/4)
			call MPI_SEND(B(length/4+1:length/2), (length/2), MPI_REAL, 1, 0, MPI_COMM_WORLD, iErr)
			call MPI_SEND(B(length/2+1:(length/4)*3), (length/2), MPI_REAL, 2, 0, MPI_COMM_WORLD, iErr)
			call MPI_SEND(B(length/4*3+1:length), (length/2), MPI_REAL, 3, 0, MPI_COMM_WORLD, iErr)
		else
			call MPI_RECV(Athread, (length/4)*length*2, MPI_REAL, 0, 0, MPI_COMM_WORLD, status,iErr)
			call MPI_RECV(Bthread, (length/4)*length*2, MPI_REAL, 0, 0, MPI_COMM_WORLD, status,iErr)
		end if

		! Вычисляем вектор невязок (Ax-b)
		call multiplication_matrix_vector(Athread, X, Axthread)
		call vectors_difference(Axthread, Bthread, Rkthread)
		!call MPI_GATHER(Rkthread, (length/4)*length*2, MPI_REAL, Rk, length*length*2, MPI_REAL, 0, MPI_COMM_WORLD, iErr)
		if (rank == 0) then
			Rk(1:length/4) = Rkthread

			call MPI_RECV(Rkthread, (length/2), MPI_REAL, 1, 1, MPI_COMM_WORLD, status,iErr)
			Rk(length/4+1:length/2) = Rkthread

			call MPI_RECV(Rkthread, (length/2), MPI_REAL, 2, 2, MPI_COMM_WORLD, status,iErr)
			Rk(length/2+1:(length/4)*3) = Rkthread

			call MPI_RECV(Rkthread, (length/2), MPI_REAL, 3, 3, MPI_COMM_WORLD, status,iErr)
			Rk(length/4*3+1:length) = Rkthread

			call show_vector(Rk)

			call multiplication_matrix_vector(A, X, Ax)
			call vectors_difference(Ax, B, Rk)

			call show_vector(Rk)
		else
			call MPI_SEND(Rkthread, (length/2), MPI_REAL, 0, rank, MPI_COMM_WORLD, iErr)
		end if


		! Результаты
		!write(*,*) "Матрица X:"
		!call show_vector(X)
			
		! Чистим память
		deallocate(X)
		deallocate(Ax)
		deallocate(Rk)
		deallocate(temp)
	
		!call cpu_time(timeStop)
	
		!write(*,*) "Время работы: ", timeStop - timeStart
	
		!write(*,*) "Программа завершена"
		call MPI_FINALIZE(iErr)
	end program main
	
	
	