      subroutine proc_astrometry(iexpo)
      implicit none
      include 'para.inc'

      integer iexpo,nchip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT

      call get_image_list(iexpo,IMAGE_FILE,nchip,DIR_OUTPUT)

      call chip_process_astrometry(IMAGE_FILE,nchip,DIR_OUTPUT)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine chip_process_astrometry(IMAGE_FILE,nchip,DIR_OUTPUT)
      implicit none
      include 'para.inc'

      integer nchip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT
      integer proc_error,ichip
      character PREFIX*(strl),filename*(strl)
      double precision cRPIX(2),cD(2,2),PU(2,npd),cRVAL(2)
      double precision cRPIX_b(2),cD_b(2,2),PU_b(2,npd),cRVAL_b(2)
      integer n,i,j,k
      double precision ra,dec,x,y,ra2,dec2

      call get_astrometry(IMAGE_FILE,nchip,DIR_OUTPUT)

      call get_PREFIX_expo(IMAGE_FILE(1),PREFIX)
      filename=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX)//'_check.dat'

      open(unit=50,file=trim(filename),status='replace')
      rewind 50

      do ichip=1,nchip
        proc_error=0

        call get_PREFIX_expo(IMAGE_FILE(1),PREFIX)
        filename=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX)//'.head'

        call read_astrometry_para(filename,ichip    
     .,cRPIX,cD,cRVAL,PU,npd,proc_error)

        if (proc_error.eq.0) then

          call get_PREFIX(IMAGE_FILE(ichip),PREFIX)
          filename=trim(DIR_OUTPUT)//'/stamps/'
     .//trim(PREFIX)//'_norm.fits'
          call update_para(filename,cRPIX,cD)

          filename=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX)//'_astro.dat'
          open(unit=40,file=trim(filename),status='old')
          rewind 40
          read(40,*)
          read(40,*)
          read(40,*) n,j,k
          do i=1,n
            read(40,*) ra,dec,x,y
            call coordinate_transfer_PU(ra2,dec2,x,y,1
     .,cRPIX,cD,cRVAL,PU,npd)
            write(50,*) ra,dec,ra2,dec2
          enddo
          close(40)
        endif
      enddo

      close(50)


      return
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_astrometry(IMAGE_FILE,nchip,DIR_OUTPUT)
      implicit none
      include 'para.inc'

      integer nchip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT
      integer ichip

      character*(strl) PREFIX,filename
      integer nss(nchip),n_user,n_ref,i,k
      integer nss_max
      parameter (nss_max=10000)
      double precision ra2(nchip,nss_max),dec2(nchip,nss_max)
      double precision x2(nchip,nss_max),y2(nchip,nss_max)

      double precision cRPIX2(nchip,2),cD2(nchip,2,2)
      double precision PU(2,npd),cRVAL2(2)
      integer valid(nchip),tot_valid,tot_source

      tot_valid=0
      tot_source=0

      do ichip=1,nchip

        call get_PREFIX(IMAGE_FILE(ichip),PREFIX)
        filename=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX)//'_astro.dat'

        open(unit=10,file=filename,status='old')
        rewind 10
        read(10,*) cRPIX2(ichip,1),cRPIX2(ichip,2)
     .,cRVAL2(1),cRVAL2(2)
        read(10,*) cD2(ichip,1,1),cD2(ichip,1,2)
     .,cD2(ichip,2,1),cD2(ichip,2,2)
        read(10,*) nss(ichip),n_user,n_ref
        do i=1,nss(ichip)
          read(10,*) ra2(ichip,i),dec2(ichip,i)
     .,x2(ichip,i),y2(ichip,i)
        enddo
        close(10)
        if (nss(ichip).ge.10) then
          valid(ichip)=1
          tot_valid=tot_valid+1
          tot_source=tot_source+nss(ichip)
        else
          valid(ichip)=0
        endif
      enddo


      if (tot_source.ge.(npd+tot_valid*3)*3) then
        call measure_astrometry_global(nss_max,nss,nchip,ra2
     .,dec2,x2,y2,cRPIX2,cD2,cRVAL2,PU,npd,valid)
        call check_astrometry_global(nss_max,nss,nchip,ra2
     .,dec2,x2,y2,cRPIX2,cD2,cRVAL2,PU,npd,valid)
      else
        do i=1,nchip
          valid(i)=0
        enddo
      endif


      call get_PREFIX_expo(IMAGE_FILE(1),PREFIX)
      filename=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX)//'.head'
      open(unit=10,file=filename,status='replace')
      rewind 10
      write(10,*) cRVAL2(1),cRVAL2(2)
      do i=1,npd
        write(10,*) PU(1,i),PU(2,i)
      enddo
      do k=1,nchip
        write(10,*) k,valid(k),cRPIX2(k,1),cRPIX2(k,2)
     .,cD2(k,1,1),cD2(k,1,2),cD2(k,2,1),cD2(k,2,2)
      enddo
      close(10)


      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_astrometry_para(filename,ichip
     .,cRPIX,cD,cRVAL,PU,npd,proc_error)
      implicit none

      integer ichip,npd,proc_error
      character filename*(*)
      double precision cRPIX(2),cD(2,2),PU(2,npd),cRVAL(2)
      integer valid,j,k,i

      if (proc_error.eq.1) return

      open(unit=11,file=filename,status='old')
      rewind 11
      read(11,*) cRVAL(1),cRVAL(2)
      do i=1,npd
        read(11,*) PU(1,i),PU(2,i)
      enddo
      do k=1,ichip-1
        read(11,*)
      enddo
      read(11,*) j,valid,cRPIX(1),cRPIX(2)
     .,cD(1,1),cD(1,2),cD(2,1),cD(2,2)

      close(11)

      if (valid.eq.0) proc_error=1

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
