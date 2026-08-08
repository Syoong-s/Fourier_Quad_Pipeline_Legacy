      module psf_storage_mod

      implicit none
      save

      include 'para.inc'
      include 'cust_para.inc'

      double precision, allocatable :: global_components(:,:,:) 
      double precision, allocatable :: global_mean_psf(:,:)
      real, allocatable :: global_poly_coefs(:,:,:,:,:)
      
      logical :: is_data_loaded = .false.

      contains
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init_and_load_all_psf(dir_output, my_rank)
      use mpi
      character*(*), intent(in) :: dir_output
      integer, intent(in) :: my_rank  
      
      integer :: i_ccd, bx, by, u, j, k, ierror
      character(len=strl) :: filename
      character(len=2)    :: c_chip
      character(len=1)    :: c_bx, c_by
      
      integer :: mpi_err
      integer :: total_count

      if (is_data_loaded) return 


      if (my_rank .eq. 0) write(*,*) 'Allocating memory on all ranks...'
      
      allocate(global_components(Camera_ccd_num,nsns,n_pcs))
      allocate(global_mean_psf(Camera_ccd_num,nsns))
      allocate(global_poly_coefs(Camera_ccd_num, 2, 2, n_pcs, npp6th))


      global_components = 0.0d0
      global_mean_psf = 0.0d0
      global_poly_coefs = 0.0


      if (my_rank .eq. 0) then
         write(*,*) 'Rank 0 is reading files from disk...'

         do i_ccd = 1, Camera_ccd_num
            if (i_ccd .eq. 2 .or. i_ccd .eq. 61) cycle
             write(c_chip, '(I2.2)') i_ccd
             filename = trim(dir_output)//'/dat_pcs/'
     .                //'pcs_ccd'//trim(c_chip)//'.dat'
             
             open(unit=30, file=filename, status='old', action='read',
     .            iostat=ierror)
             if (ierror .eq. 0) then
                rewind 30
                do k = 1, nsns
                  read(30,*) (global_components(i_ccd,k,j), j=1,n_pcs),
     .                        global_mean_psf(i_ccd,k)
                enddo
                close(30)
                if (global_components(i_ccd,1,1) .lt. -1.0d20) then
                  write(*,*) 'CCD',i_ccd,' has bad PCS data.'
                endif
             else
                global_components(i_ccd,1,1) = -1.0d30
                write(*,*) 'CCD',i_ccd,' PCS data read error.'
                stop
             endif


             do bx = 1, nblocks
               do by = 1, nblocks
                write(c_bx, '(I1.1)') bx
                write(c_by, '(I1.1)') by
                filename = trim(dir_output)//'/dat_pcs/'//'coeff_ccd'// 
     .             trim(c_chip)//'_'//trim(c_bx)//trim(c_by)//'.dat'
                open(unit=20, file=filename, status='old', 
     .                action='read', iostat=ierror)
                if (ierror .eq. 0) then
                  rewind 20
                  do u = 1, n_pcs
         read(20,*) (global_poly_coefs(i_ccd,bx,by,u,j),j=1, npp6th)
                      if (isnan(global_poly_coefs(i_ccd,bx,by,u,1))) 
     .global_poly_coefs(i_ccd,bx,by,1,1) = -1.0e30
                  enddo
                  close(20)
      if (global_poly_coefs(i_ccd,bx,by,1,1) .lt. -1.0e20) then
        write(*,*) 'CCD',i_ccd,'field',bx,by,' has bad polynomial data.'
                  endif
                else
                    global_poly_coefs(i_ccd,bx,by,1,1) = -1.0e30
                    write(*,*) 'CCD',i_ccd,'field',bx,by,' polynomial'
     .//' data read error.'
                    stop
                 endif
               enddo
             enddo
         enddo
         write(*,*) 'Rank 0 finished reading. Starting Broadcast...'
      endif
c     -------------------------------------------------------------------
      total_count = Camera_ccd_num * nsns * n_pcs
      call MPI_BCAST(global_components, total_count
     ., MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_err)


      total_count = Camera_ccd_num * nsns
      call MPI_BCAST(global_mean_psf, total_count
     ., MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_err)


      total_count = Camera_ccd_num * 2 * 2 * n_pcs * npp6th
      call MPI_BCAST(global_poly_coefs, total_count, MPI_REAL, 0, 
     .               MPI_COMM_WORLD, mpi_err)

      is_data_loaded = .true.
      
      call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

      if (my_rank .eq. 0) write(*,*) 'Broadcast finished. All ready.'

      end subroutine init_and_load_all_psf
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine free_psf_memory()
        if (allocated(global_components)) deallocate(global_components)
        if (allocated(global_mean_psf))   deallocate(global_mean_psf)
        if (allocated(global_poly_coefs)) deallocate(global_poly_coefs)
        is_data_loaded = .false.
      end subroutine free_psf_memory

      end module psf_storage_mod