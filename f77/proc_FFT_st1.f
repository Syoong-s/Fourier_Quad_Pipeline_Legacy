      subroutine proc_FourierT_st1(iexpo)
      implicit none
      include 'para.inc'

      integer iexpo,nchip,ichip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT

      call get_image_list(iexpo,IMAGE_FILE,nchip,DIR_OUTPUT)

      do ichip=1,nchip
        call chip_process_Fourier_T_st1(IMAGE_FILE(ichip),DIR_OUTPUT)
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc\
      subroutine chip_process_Fourier_T_st1(IMAGE_FILE,DIR_OUTPUT)
      implicit none
      include 'para.inc'

      character*(strl) IMAGE_FILE,DIR_OUTPUT

      integer ierror
      integer i,j,u,v,nn1,nn2
      character*(strl) PREFIX,filename

      integer nsource
      real source_coll(ngal_max,ns,ns)
      real noise_coll(ngal_max,ns,ns)
      real power_coll(ngal_max,ns,ns)
      real power_ori_coll(ngal_max,ns,ns)

      real source_p(ns,ns),source(ns,ns)
      real noise_p(ns,ns),noise(ns,ns)

      real aa(npara)

c-------------------------------------------------------------------
      if (ext_PSF.eq.1) return
      call get_PREFIX(IMAGE_FILE,PREFIX)
      PREFIX=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)
c-------------------------------------------------------------------
      nsource=0
      filename=trim(PREFIX)//'_star_can_info.dat'
      open(unit=10,file=filename,status='old',iostat=ierror)
      rewind 10
      if (ierror.ne.0) then
        write(*,*) filename
        write(*,*) 'Error / FFT1 star_can_info catalog file error!!'
        read(*,*)
      endif
      read(10,*)
c      read(10,*) 'ig xp yp SNR'
      do while (ierror.ge.0)
        read(10,*,iostat=ierror) (aa(i),i=1,4)
        if (ierror.lt.0) cycle
        nsource=nsource+1
      enddo
      close(10)

      if (nsource.gt.0) then 
        nn1=ns*len_s
        nn2=ns*(int(nsource/len_s)+1)
        filename=trim(PREFIX)//'_star_can.fits'
        call read_stamps(ngal_max,1,nsource,ns,ns
     .,source_coll,nn1,nn2,filename)

        filename=trim(PREFIX)//'_star_can_noise.fits'
        call read_stamps(ngal_max,1,nsource,ns,ns
     .,noise_coll,nn1,nn2,filename)

        do i=1,nsource
          do u=1,ns
            do v=1,ns
              source(u,v)=source_coll(i,u,v)
              noise(u,v)=noise_coll(i,u,v)
            enddo
          enddo
          call get_power(ns,ns,source,source_p,star_smooth)
          call get_power(ns,ns,noise,noise_p,star_smooth)
          call process_powers(ns,source_p,noise_p)
          call regularize_power(ns,ns,source_p,star_smooth)
          do u=1,ns
            do v=1,ns
              power_coll(i,u,v)=source_p(u,v)
            enddo
          enddo
c        call get_power(ns,ns,source,source_p,0)
c        call get_power(ns,ns,noise,noise_p,0)
c        call process_powers(ns,source_p,noise_p)
c        call regularize_power(ns,ns,source_p,0)
c        do u=1,ns
c          do v=1,ns
c            power_ori_coll(i,u,v)=source_p(u,v)
c          enddo
c        enddo
        enddo

        filename=trim(PREFIX)//'_star_can_power.fits'
        nn1=ns*len_s
        nn2=ns*(int(nsource/len_s)+1)
        call write_stamps(ngal_max,1,nsource,ns,ns
     .,power_coll,nn1,nn2,filename)
      endif

c      filename=trim(PREFIX)//'_star_ori_power.fits'
c      nn1=ns*len_s
c      nn2=ns*(int(nsource/len_s)+1)
c      call write_stamps(ngal_max,1,nsource,ns,ns
c     .,power_ori_coll,nn1,nn2,filename)

c----------------------------------------------------
      ! write(*,*) trim(PREFIX),nsource



      return
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
