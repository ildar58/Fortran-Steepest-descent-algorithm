program ver_MPI
implicit none
include 'mpif.h'
	! Declorations for MPI
	integer :: iErr         
    integer :: rank
    integer :: nProcess
    integer, dimension(MPI_STATUS_SIZE) :: status
    ! Variables
    integer, allocatable :: values(:,:,:)
    integer, allocatable :: newValues(:,:,:)
    integer :: xSize, ySize, zSize, x, y, z, maxValue, anMaxValue, i, zPartBegin, zPartEnd, zPartSize
    integer :: cube(3,3,3)
    integer :: xMask(3,3,3), yMask(3,3,3), zMask(3,3,3)
    real :: gX, gY, gZ, k, timeStart, timeStop
    CHARACTER(100) :: numChar

    CALL GET_COMMAND_ARGUMENT(1,numChar)
    READ(numChar,*)xSize
    CALL GET_COMMAND_ARGUMENT(2,numChar)
    READ(numChar,*)ySize 
    CALL GET_COMMAND_ARGUMENT(3,numChar)
    READ(numChar,*)zSize

    print *, 'xSize', xSize
    print *, 'ySize', ySize
    print *, 'zSize', zSize
    
    call MPI_INIT(iErr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcess, iErr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, iErr)
    
!!!! Разделяем для разделения данных
    ! Читаем маски
    open (unit = 2,file = "mask_x.txt")
    open (unit = 3,file = "mask_y.txt")
    open (unit = 4,file = "mask_z.txt")
    do x = 1, 3
        do y = 1, 3
            do z = 1, 3
              read (2,*) xMask(x,y,z)
              read (3,*) yMask(x,y,z)  
              read (4,*) zMask(x,y,z)  
            end do
        end do
    end do  
    close(1)
    close(2)
    close(3)
    print *, '- Process',rank, 'The masks were read.'

    zPartSize  = (zSize/nProcess)

    ! выделяем память, генерируем данные, отправляем данные
    if (rank == 0) then
        ! Выделяем память и генерируем данные
        allocate (values(xSize+2,ySize+2,zSize+2))
        allocate (newValues(xSize,ySize,zSize))
		
		print *, '- Process',rank, 'Generating and sending data...'
		values = 0
		do z = 2, zSize+1
			do y = 2, ySize+1
				do x = 2, xSize+1
                    !call random_number(k)
                    values(x,y,z) = mod(x*y*z,256)
				end do
			end do
		end do
        print *, '- Process',rank, 'Data was generated.'

        CALL CPU_TIME(timeStart)

        ! Отправляем данные
        if (nProcess > 1) then
            do i = 1, nProcess - 1
                call MPI_SSend(&
                   values(:,:,(i)*zPartSize+1 : (i+1)*zPartSize+2),&
                   (xSize+2)*(ySize+2)*(zPartSize+2),MPI_INT, i, 100, MPI_COMM_WORLD, iErr)
            end do
            print *, '- Process',rank, 'Data was sent.'
        endif

	! выделяем память, принимаем
	else
        allocate (values(xSize+2,ySize+2,zPartSize+2))
        allocate (newValues(xSize,ySize,zPartSize))

        call MPI_Recv(&
            values,&
            (xSize+2)*(ySize+2)*(zPartSize+2),MPI_INT,0,100,MPI_COMM_WORLD,status,iErr)
      
        print *, '- Process',rank, 'The data obtained.'
        
    end if
    
    ! Обрабатываем
    print *, '- Process',rank, 'Data processing...'
    newValues = 0
    maxValue = -1

    do z = 2, zPartSize+1
        do y = 2, ySize+1
            do x = 2, xSize+1
                cube = values(x-1:x+1, y-1:y+1,z-1:z+1)
                gX = real(sum(cube*xMask))
                gY = real(sum(cube*yMask)) 
                gZ = real(sum(cube*zMask))
                newValues(x-1,y-1,z-1) = int(sqrt(gX**2+gY**2+gZ**2))
                           
                if (newValues(x-1,y-1,z-1) > maxValue) then
                    maxValue = newValues(x-1,y-1,z-1) 
                end if
            end do
        end do
    end do
    print *, '- Process',rank,'The data was processed.'

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!!!! Разделяем для соединения данных
    if (rank == 0) then
        if (nProcess > 1) then
            do i = 1, nProcess - 1

                call MPI_Recv(&
                newValues(:,:,(i)*zPartSize+1 :(i+1)*zPartSize ),&
                zPartSize*xSize*ySize, MPI_INT, i, 200, MPI_COMM_WORLD,status,iErr)

                call MPI_Recv(anMaxValue, 1, MPI_INT, i, 300, MPI_COMM_WORLD,status,iErr)

                if (anMaxValue > maxValue) then
                    maxValue = anMaxValue
                end if
            end do
            print *, '- Process',rank, 'The results obtained.'
        end if
    else
        call MPI_SSend(newValues,zPartSize*xSize*ySize,MPI_INT,0,200,MPI_COMM_WORLD,iErr)
        call MPI_SSend(maxValue ,1                    ,MPI_INT,0,300,MPI_COMM_WORLD,iErr)
        print *, '- Process',rank, 'Results were sent.'
    end if

call MPI_FINALIZE(iErr)   
!!!! Разделяем, т.к. всё остальное доделается в 0 процессе
    if (rank == 0) then
        CALL CPU_TIME(timeStop)
        print *, '- Process',rank, 'The data was processed. Time:'
        print *, '- Process',rank, timeStop-timeStart
        write (*,*) maxValue
 
        print *, '- Process',rank, 'Postprocess...'     
        k = 255 / real(maxValue)

        ! Приводим ответ к 0..255    
        do x = 1, xSize
            do y = 1, ySize
                do z = 1, zSize        
                    newValues(x,y,z) = newValues(x,y,z)*k
                end do
            end do
        end do    
        
        print *, '- Process',rank, 'Finish!'
        write (*,*) 'Bye.'
    end if

    deallocate(values)
    deallocate(newvalues)
end program ver_MPI

