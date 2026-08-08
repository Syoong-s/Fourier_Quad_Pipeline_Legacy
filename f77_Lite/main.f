      program main
      use mpi
      include 'para.inc'

      integer ierr, my_id, num_procs
      common /MPIpar/my_id, num_procs

      character*(strl) EXPO_LIST,root_dir,filename
      character*(strl) EXPO_FILE(NMAX_EXPO)
      integer N_EXPO
      common /filename_pass/ EXPO_FILE,N_EXPO

      real expo_para(6,NMAX_EXPO),expo_para_t(6,NMAX_EXPO)
      common /expo_para_pass/ expo_para
      integer iexpo,i,j,rng_seed

      external pre_process,proc_astrometry,proc_source,proc_comb
      external proc_FourierT_st1,proc_FourierT_st2
      external proc_PSF,proc_shear,proc_info

      !----- MPI initialization ---------
      call MPI_Init( ierr )
      call MPI_comm_rank(MPI_cOMM_WORLD, my_id, ierr )
      call MPI_comm_size(MPI_cOMM_WORLD, num_procs, ierr )
c ==========================================
c Function: Initialize one RNG stream per MPI rank
c Method: Mix the local clock and rank before the first RNG call
c ==========================================
      call initialize_ran1_seed(my_id,num_procs,rng_seed)
      write(*,*) 'RNG_SEED rank seed: ',my_id,rng_seed
      ! my_id = 0,1,...,num_procs-1
      call MPI_BARRIER(MPI_cOMM_WORLD, ierr ) ! synchronize all nodes

      call getarg(1,EXPO_LIST)
      ! ----------------------------------------------------------------------------
      if (my_id.eq.0) call initialize(EXPO_LIST)
      call MPI_Bcast(N_EXPO,1,mpi_int,0,MPI_cOMM_WORLD,ierr)
      call MPI_Bcast(EXPO_FILE,NMAX_EXPO*strl,mpi_character,0
     .,MPI_cOMM_WORLD,ierr)

      call MPI_BARRIER(MPI_cOMM_WORLD,ierr)
      ! ----------------------------------------------------------------------------
      if (mod(PROcESS_stage,2).eq.0)
     .call mpi_distribute(N_EXPO,pre_process,'Pre-process...')
      call MPI_BARRIER(MPI_cOMM_WORLD,ierr)

      if (mod(PROcESS_stage,3).eq.0)
     .call mpi_distribute(N_EXPO,proc_astrometry,'Astrometry...')
      call MPI_BARRIER(MPI_cOMM_WORLD,ierr)

      if (mod(PROcESS_stage,5).eq.0)
     .call mpi_distribute(N_EXPO,proc_source,'Sources ...')
      call MPI_BARRIER(MPI_cOMM_WORLD,ierr)

      if (mod(PROcESS_stage,7).eq.0)
     .call mpi_distribute(N_EXPO,proc_FourierT_st1,'FFT st1...')
      call MPI_BARRIER(MPI_cOMM_WORLD,ierr)

      if (mod(PROcESS_stage,11).eq.0)
     .call mpi_distribute(N_EXPO,proc_PSF,'PSF ...')
      call MPI_BARRIER(MPI_cOMM_WORLD,ierr)

      if (mod(PROcESS_stage,13).eq.0)
     .call mpi_distribute(N_EXPO,proc_FourierT_st2,'FFT st2...')
      call MPI_BARRIER(MPI_cOMM_WORLD,ierr)


      if (mod(PROcESS_stage,17).eq.0)
     .call mpi_distribute(N_EXPO,proc_shear,'Shear ...')
      call MPI_BARRIER(MPI_cOMM_WORLD,ierr)

      do iexpo=1,NMAX_EXPO
        do i=1,6
          expo_para(i,iexpo)=0.
        enddo
      enddo

      if (mod(PROcESS_stage,19).eq.0)
     .call mpi_distribute(N_EXPO,proc_info,'Info ...')

      call MPI_BARRIER(MPI_cOMM_WORLD,ierr)
      call MPI_AllReduce(expo_para,expo_para_t,NMAX_EXPO*6
     .,mpi_REAL,MPI_SUM,MPI_cOMM_WORLD,ierr)

      do iexpo=1,NMAX_EXPO
        do i=1,6
          expo_para(i,iexpo)=expo_para_t(i,iexpo)
        enddo
      enddo

      if (my_id.eq.0) then
        call get_dir(EXPO_LIST,root_dir,1)
        filename=trim(root_dir)//'/expo_info.dat'
        open(unit=10,file=filename,status='replace')
        rewind 10
        write(10,*) 'N-valid-chip PSF-FWHM(arcsec) chi_d-stars'    
     .,' nstar-per-chip cRVAL1 cRVAL2 expo_name '
        do i=1,N_EXPO
          write(10,*) (expo_para(j,i),j=1,6),trim(EXPO_FILE(i))
        enddo
        close(10)
      endif

      if (mod(PROcESS_stage,23).eq.0)
     .call mpi_distribute(N_EXPO,proc_comb,'combine ...')
      call MPI_BARRIER(MPI_cOMM_WORLD,ierr)

c      !----- MPI finalization -----------------------------------------------
      call MPI_BARRIER(MPI_cOMM_WORLD,ierr)
      call MPI_Finalize (ierr)

      stop
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initialize(EXPO_LIST)
      implicit none
      include 'para.inc'

      character*(strl) EXPO_LIST
      character*(strl) EXPO_FILE(NMAX_EXPO)
      integer N_EXPO
      common /filename_pass/ EXPO_FILE,N_EXPO

      character*(strl) expo_name
      integer ierror,nchip

      N_EXPO=0

      open(unit=10,file=EXPO_LIST,status='old',iostat=ierror)
      rewind 10
      if (ierror.ne.0) then
        write(*,*) 'EXPO_LIST reading error!!'
        stop
      endif
      do while (ierror.ge.0)
        read(10,*,iostat=ierror) expo_name,nchip
        if (ierror.lt.0) cycle
        N_EXPO=N_EXPO+1
        EXPO_FILE(N_EXPO)=trim(expo_name)
      enddo
      close(10)
      write(*,*) 'Total number of EXPOSURE: ',N_EXPO


      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
