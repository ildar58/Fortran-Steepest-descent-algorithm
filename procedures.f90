module procedures
	implicit none
	public
	
	contains
	
	! Разность векторов
	subroutine vectors_difference(A, B, result_vector)
		integer :: i, length
		real*8, dimension(:), intent(in) :: A, B
		real*8, dimension(:), intent(out) :: result_vector
		
		length = size(A)
		
		do i=1,length
			result_vector(i) = A(i) - B(i)
		end do
	end subroutine vectors_difference
	
	! Умножение матрицы на вектор
	subroutine multiplication_matrix_vector(matrix, vector, result_vector)
		integer :: i, j, length_m, length_v
		real*8, dimension(:,:), intent(in) :: matrix
		real*8, dimension(:), intent(in) :: vector
		real*8, dimension(:), intent(out) :: result_vector
		
		length_m = size(matrix, 1)
		length_v = size(vector, 1)
		
		do i=1,length_m
			result_vector(i) = 0
			do j=1,length_v
				result_vector(i) = result_vector(i) + matrix(i,j) * vector(j)
			end do
		end do
	end subroutine multiplication_matrix_vector
	
	! Умножение числа на вектор
	subroutine multiplication_num_vector(vector, num, result_vector)
		integer :: i, length
		real*8, intent(in) :: num
		real*8, dimension(:), intent(in) :: vector
		real*8, dimension(:), intent(out) :: result_vector
		
		length = size(vector)
		
		do i=1,length
			result_vector(i) = vector(i) * num
		end do
	end subroutine multiplication_num_vector
	
	! Вывод вектора
	subroutine show_vector(vector)
		integer :: i, length
		real*8, dimension(:), intent(in) :: vector
		
		length = size(vector)
		
		do i=1,length
			write(*,20) i, vector(i)
		end do
		
		write(*,*) ""
		20 format(1x,i4,20f8.2)
	end subroutine
	
	! Вывод матрицы
	subroutine show_matrix(matrix, lengthA, lengthB)
		integer :: i, j, length, lengthA, lengthB
		real*8, dimension(:,:), intent(in) :: matrix

		do i=1,lengthA
			write(*,20) i, (matrix(i,j),j=1,lengthB)
		end do
		
		write(*,*) ""
		20 format(1x,i4,20f8.2)
	end subroutine
	
	! Скалярное произведение векторов
	real*8 function vectors_scalar_multiplication(A, B) result(res)
		integer :: i, length
		real*8, dimension(:), intent(in) :: A, B
		
		length = size(A)
		res = 0
		
		do i=1,length
			res = res + A(i) * B(i)
		end do
		
		return
	end function vectors_scalar_multiplication
	
	! Проверка сходимости
	logical function check_result(vector, eps) result(res)
		real*8, dimension(:), intent(in) :: vector
		real*8, intent(in) :: eps
		integer :: i, length
		
		length = size(vector)
		
		do i=1,length
			res = ABS(vector(i)) > eps
				if (.not. res) then
					exit
				end if
		end do
		return
	end function check_result
	
	subroutine fill_matrix(matrix, length)
		integer :: i, j, k
		integer, intent(in) :: length
		real*8, dimension(length, length), intent(in out) :: matrix
		call random_seed()
		do i=1, length
			do j=1, length
				!call random_number(matrix(i,j))
				matrix(i,j) = 0.1
			end do
		end do
		do i=1, length
			do j=i, length
					if (i == j) then
						do k = 1,length
							matrix(i,j) = matrix(i,j) + matrix(i,k)
						end do
					else
				end if
			end do
		end do
	end subroutine
	
	subroutine fill_vector(vector, length)
		integer :: i, j
		integer, intent(in) :: length
		real*8, dimension(length), intent(in out) :: vector
		do i=1, length
			vector(i) = i*5.42
		end do		
	end subroutine
end module procedures
