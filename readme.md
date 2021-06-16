For serial compilation use gfortran -o serial procedures.f90 serial.f90
* To run use ./serial

For MPI compilation use /usr/bin/mpif90 -o mpi procedures.f90 mpi.f90  
* To run in 4 threads use  /usr/bin/mpirun -n 4 ./mpi
