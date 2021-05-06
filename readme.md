Компиляция gfortran -o serial procedures.f90 serial.f90 для Serial
/usr/bin/mpif90 -o mpi procedures.f90 mpi.f90  для MPI (запуск в 4х потоках ) /usr/bin/mpirun -n 4 ./mpi