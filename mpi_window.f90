	program main
		use procedures
		implicit none
		include 'mpif.h'
		
		!!! Переменные
		! A - исходная матрица коэф., B - правая часть, X - вектор решения, Ax(A*x), Rk(невязка)
		integer :: length = 20
		integer :: i, j, iErr, nProcess, rank, status(MPI_STATUS_SIZE), iter=0
		real*8, allocatable :: A_BUF(:, :), B_BUF(:)
		integer (KIND=MPI_ADDRESS_KIND) lb, extent, target_disp
		real*8, allocatable :: A(:, :), B(:), X(:), Ax(:), Rk(:), temp(:), Athread(:,:), Bthread(:), Rkthread(:), Axthread(:)
		real*8, allocatable :: tempThread(:)
		real*8 :: timeStart, timeStop, tau, eps=0.000001
		integer :: winA, winB

		call MPI_INIT(iErr)
		call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcess, iErr)
		call MPI_COMM_RANK(MPI_COMM_WORLD, rank, iErr)
	
		
		call MPI_TYPE_GET_EXTENT(MPI_REAL, lb, extent, ierr)
		allocate(A_BUF(length,length))
		allocate(B_BUF(length))
		
		call MPI_WIN_CREATE(A_BUF, length*length*extent, extent, MPI_INFO_NULL, MPI_COMM_WORLD, winA, iErr)
		call MPI_WIN_CREATE(B_BUF, length*extent, extent, MPI_INFO_NULL, MPI_COMM_WORLD, winB, iErr)
		call MPI_WIN_FENCE(0, winA, iErr)
		call MPI_WIN_FENCE(0, winB, iErr)

		!! Выделяем память
		allocate(A(length,length))
		allocate(B(length))
		allocate(Athread(length/4,length))
		allocate(X(length))
		allocate(Bthread(length/4))
		allocate(Ax(length))
		allocate(Axthread(length/4))
		allocate(Rk(length))
		allocate(Rkthread(length/4))
		allocate(temp(length))
		allocate(tempThread(length/4))

		call cpu_time(timeStart)
	! Генерируем данные
		call fill_matrix(A, length)
		call fill_vector(B, length)
		do i=1,length
			X(i)= 0
		end do

		call cpu_time(timeStart)
		if (rank == 0) then
			Athread = A(1:length/4, 1:length)

			call MPI_Put(A(length/4+1:length/2,1:length), (length/4)*length*2, MPI_REAL, 1, &
			2*length*length, (length/4)*length*2, MPI_REAL, winA, iErr)
			call MPI_Put(A(length/2+1:(length/4)*3,1:length), (length/4)*length*2, MPI_REAL, 2,&
			 3*length*length, (length/4)*length*2, MPI_REAL, winA, iErr)
			call MPI_Put(A(length/4*3+1:length,1:length), (length/4)*length*2, MPI_REAL, 3, &
			4*length*length, (length/4)*length*2, MPI_REAL, winA, iErr)
			Bthread = B(1:length/4)

			call MPI_Put(B(length/4+1:length/2), (length/2), MPI_REAL, 1, 2*length, (length/2), MPI_REAL, winB, iErr)
			call MPI_Put(B(length/2+1:(length/4)*3), (length/2), MPI_REAL, 2, 3*length, (length/2), MPI_REAL, winB, iErr)
			call MPI_Put(B(length/4*3+1:length), (length/2), MPI_REAL, 3, 4*length, (length/2), MPI_REAL, winB, iErr)
		else
			!call MPI_RECV(Athread, (length/4)*length*2, MPI_REAL, 0, 0, MPI_COMM_WORLD, status,iErr)
			call MPI_Get(Athread, (length/4)*length*2, MPI_REAL, 0, rank*length*length, (length/4)*length*2, MPI_REAL, winA, iErr)
			call MPI_WIN_FENCE(0, winA, iErr)
			!call MPI_RECV(Bthread, (length/4)*length*2, MPI_REAL, 0, 0, MPI_COMM_WORLD, status,iErr)
			call MPI_Get(Bthread, length/2, MPI_REAL, 0, rank*length, length/2, MPI_REAL, winB, iErr)
			call MPI_WIN_FENCE(0, winB, iErr)
		end if
		if (rank /= 0) then
			call show_matrix(Athread, length, length/4)
		end if

		! Вычисляем вектор невязок (Ax-b)
		call multiplication_matrix_vector(Athread, X, Axthread)
		call vectors_difference(Axthread, Bthread, Rkthread)
		if (rank == 0) then
			Rk(1:length/4) = Rkthread

			call MPI_RECV(Rkthread, (length/2), MPI_REAL, 1, 1, MPI_COMM_WORLD, status,iErr)
			Rk(length/4+1:length/2) = Rkthread

			call MPI_RECV(Rkthread, (length/2), MPI_REAL, 2, 2, MPI_COMM_WORLD, status,iErr)
			Rk(length/2+1:(length/4)*3) = Rkthread

			call MPI_RECV(Rkthread, (length/2), MPI_REAL, 3, 3, MPI_COMM_WORLD, status,iErr)
			Rk(length/4*3+1:length) = Rkthread
		else
			call MPI_SEND(Rkthread, (length/2), MPI_REAL, 0, rank, MPI_COMM_WORLD, iErr)
		end if

		call MPI_BCAST(Rk, length*2, MPI_REAL, 0, MPI_COMM_WORLD, iErr)
		call multiplication_matrix_vector(Athread, Rk, Axthread)
		call MPI_BARRIER(MPI_COMM_WORLD,iErr)
		!-----------------------------------------------------------------------------------------
		do while(check_result(Rk, eps))
			call multiplication_matrix_vector(Athread, Rk, Axthread)
			call MPI_BARRIER(MPI_COMM_WORLD,iErr)

			if (rank == 0) then
				Ax(1:length/4) = Axthread
	
				call MPI_RECV(Axthread, (length/2), MPI_REAL, 1, 1, MPI_COMM_WORLD, status,iErr)
				Ax(length/4+1:length/2) = Axthread
	
				call MPI_RECV(Axthread, (length/2), MPI_REAL, 2, 2, MPI_COMM_WORLD, status,iErr)
				Ax(length/2+1:(length/4)*3) = Axthread
	
				call MPI_RECV(Axthread, (length/2), MPI_REAL, 3, 3, MPI_COMM_WORLD, status,iErr)
				Ax(length/4*3+1:length) = Axthread
			else
				call MPI_SEND(Axthread, (length/2), MPI_REAL, 0, rank, MPI_COMM_WORLD, iErr)
			end if
			
			if (rank == 0) then
				tau = vectors_scalar_multiplication(Rk, Rk) / vectors_scalar_multiplication(Ax, Rk)
			endif
			call MPI_BARRIER(MPI_COMM_WORLD,iErr)
			call MPI_BCAST(tau, 2, MPI_REAL, 0, MPI_COMM_WORLD, iErr)

			call multiplication_num_vector(Rkthread, tau, tempThread)
			call multiplication_num_vector(Axthread, tau, Axthread)

			call MPI_BARRIER(MPI_COMM_WORLD,iErr)
			if (rank == 0) then
				temp(1:length/4) = tempThread
	
				call MPI_RECV(tempThread, (length/2), MPI_REAL, 1, 1, MPI_COMM_WORLD, status,iErr)
				temp(length/4+1:length/2) = tempThread
	
				call MPI_RECV(tempThread, (length/2), MPI_REAL, 2, 2, MPI_COMM_WORLD, status,iErr)
				temp(length/2+1:(length/4)*3) = tempThread
	
				call MPI_RECV(tempThread, (length/2), MPI_REAL, 3, 3, MPI_COMM_WORLD, status,iErr)
				temp(length/4*3+1:length) = tempThread
			else
				call MPI_SEND(tempThread, (length/2), MPI_REAL, 0, rank, MPI_COMM_WORLD, iErr)
			end if

			if (rank == 0) then
				call vectors_difference(X, temp, X)
			endif
			call MPI_BARRIER(MPI_COMM_WORLD,iErr)
			call MPI_BCAST(X, length*2, MPI_REAL, 0, MPI_COMM_WORLD, iErr)

			call vectors_difference(Rkthread, Axthread, Rkthread)

			call MPI_BARRIER(MPI_COMM_WORLD,iErr)
			if (rank == 0) then
				Rk(1:length/4) = Rkthread
	
				call MPI_RECV(Rkthread, (length/2), MPI_REAL, 1, 1, MPI_COMM_WORLD, status,iErr)
				Rk(length/4+1:length/2) = Rkthread
	
				call MPI_RECV(Rkthread, (length/2), MPI_REAL, 2, 2, MPI_COMM_WORLD, status,iErr)
				Rk(length/2+1:(length/4)*3) = Rkthread
	
				call MPI_RECV(Rkthread, (length/2), MPI_REAL, 3, 3, MPI_COMM_WORLD, status,iErr)
				Rk(length/4*3+1:length) = Rkthread
			else
				call MPI_SEND(Rkthread, (length/2), MPI_REAL, 0, rank, MPI_COMM_WORLD, iErr)
			end if

			call MPI_BCAST(Rk, length*2, MPI_REAL, 0, MPI_COMM_WORLD, iErr)
			call MPI_BARRIER(MPI_COMM_WORLD, iErr)
			iter = iter + 1
		end do
		!-----------------------------------------------------------------------------------------
		

		if (rank == 0) then
			! Результаты
			write(*,*) "Матрица X:"
			call show_vector(X)
			call show_vector(Rk)
			write(*,*) iter
		end if
		
			
		! Чистим память
		deallocate(X)
		deallocate(Ax)
		deallocate(Rk)
		deallocate(temp)

		call cpu_time(timeStop)

		write(*,*) "Время работы: ", timeStop - timeStart
	
		!call cpu_time(timeStop)
	
		!write(*,*) "Время работы: ", timeStop - timeStart
	
		!write(*,*) "Программа завершена"
		call MPI_Win_free(winA, iErr)
		call MPI_Win_free(winB, iErr)
		call MPI_FINALIZE(iErr)
	end program main
	
	
	
