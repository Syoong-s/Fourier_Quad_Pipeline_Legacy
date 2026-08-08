C     Verify MPI Fortran startup and a reference LAPACK solve.
      program mpi_lapack_smoke
      include 'mpif.h'
      integer ierr
      integer info
      integer ipiv(1)
      double precision matrix(1,1)
      double precision rhs(1)

      call MPI_Init(ierr)
      matrix(1,1) = 2.0d0
      rhs(1) = 4.0d0
      call dgesv(1,1,matrix,1,ipiv,rhs,1,info)

      if (info .ne. 0) then
         call MPI_Abort(MPI_COMM_WORLD,1,ierr)
      endif
      if (dabs(rhs(1)-2.0d0) .gt. 1.0d-12) then
         call MPI_Abort(MPI_COMM_WORLD,2,ierr)
      endif

      call MPI_Finalize(ierr)
      end
