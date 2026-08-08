C     Report MPI rank placement in single-node and multi-node tests.
      program mpi_identity
      include 'mpif.h'
      integer ierr
      integer rank
      integer world_size
      integer name_length
      character*(MPI_MAX_PROCESSOR_NAME) processor_name

      call MPI_Init(ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
      call MPI_Comm_size(MPI_COMM_WORLD,world_size,ierr)
      call MPI_Get_processor_name(processor_name,name_length,ierr)

      if (world_size .lt. 2) then
         call MPI_Abort(MPI_COMM_WORLD,3,ierr)
      endif

      write(*,100) rank,world_size,processor_name(1:name_length)
  100 format('rank ',I6,' of ',I6,' on ',A)

      call MPI_Finalize(ierr)
      end
