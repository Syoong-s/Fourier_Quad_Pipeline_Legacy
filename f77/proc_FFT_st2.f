      subroutine proc_FourierT_st2(iexpo)
      implicit none
      include 'para.inc'
      include 'cust_para.inc'

      integer iexpo,nchip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT

      call get_image_list(iexpo,IMAGE_FILE,nchip,DIR_OUTPUT)

      call expo_Fourier_T_st2(nchip,IMAGE_FILE,DIR_OUTPUT,chipnx,chipny)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ==========================================
c Function: Transform source stamps for one exposure
c Method: Read split source products and write source power by chip
c ==========================================
      subroutine expo_Fourier_T_st2(nchip,IMAGE_FILE,DIR_OUTPUT
     .                                                  ,chipnx,chipny)
      implicit none
      include 'para.inc'
      include 'path_layout.inc'

      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT

      integer ichip,nchip,chipnx,chipny
      integer ierror
      integer i,j,u,v,nn1,nn2
      character*(strl) PREFIX,filename

      integer nsource
      real source_para(ngal_max,npara)
      real source_coll(ngal_max,ns,ns)
      real noise_coll(ngal_max,ns,ns)
      real power_coll(ngal_max,ns,ns)

      real source_p(ns,ns),source(ns,ns)
      real noise_p(ns,ns),noise(ns,ns)

      real aa(npara),flux_alt,temp

      real pc
      common /pc_pass/ pc

      do ichip=1,nchip    
c------------------------------------------------------------------
        nsource=0
        call fq_chip_product_path(IMAGE_FILE(ichip),DIR_OUTPUT,
     .  DIR_SRC_INFO,
     .    '_source_info.dat',filename)
        open(unit=10,file=filename,status='old',iostat=ierror)
        rewind 10
        if (ierror.ne.0) then
          write(*,*) filename
          write(*,*) 'Error / FFT2 source_info catalog file error!!'
          read(*,*)
        endif
        read(10,*)
c      read(10,*) 'ig xc yc sigma peak imax jmax'
c     .,'half_light_flux half_light_area flag'
        do while (ierror.ge.0)
          read(10,*,iostat=ierror) (aa(i),i=1,iflag)
          if (ierror.lt.0) cycle
          nsource=nsource+1
          do i=1,iflag
            source_para(nsource,i)=aa(i)
          enddo
        enddo
        close(10)
        if (nsource.eq.0) goto 991

        nn1=ns*len_g
        nn2=ns*(int(nsource/len_g)+1)
        call fq_chip_product_path(IMAGE_FILE(ichip),DIR_OUTPUT,
     .  DIR_SRC,
     .    '_source.fits',filename)
        call read_stamps(ngal_max,1,nsource,ns,ns
     .,source_coll,nn1,nn2,filename)

        call fq_chip_product_path(IMAGE_FILE(ichip),DIR_OUTPUT,
     .  DIR_NOISE,
     .    '_noise.fits',filename)
        call read_stamps(ngal_max,1,nsource,ns,ns
     .,noise_coll,nn1,nn2,filename)


        do i=1,nsource
          do u=1,ns
            do v=1,ns
              source(u,v)=source_coll(i,u,v)
              noise(u,v)=noise_coll(i,u,v)
            enddo
          enddo

          call get_power(ns,ns,source,source_p,2)
          source_para(i,11)=sqrt(max(pc,source_p(ns_2+1,ns_2+1)))
          source_para(i,12)=source_para(i,11)/source_para(i,4)*ns

          call get_power(ns,ns,source,source_p,gal_smooth)
          call get_power(ns,ns,noise,noise_p,gal_smooth)
          call process_powers(ns,source_p,noise_p)
          do u=1,ns
            do v=1,ns
              power_coll(i,u,v)=source_p(u,v)
            enddo
          enddo
        enddo

991     call fq_chip_product_path(IMAGE_FILE(ichip),DIR_OUTPUT,
     .  DIR_SRC_INFO,
     .    '_source_info.dat',filename)
        open(unit=10,file=filename,status='replace')
        rewind 10
        write(10,*) 'ig xp yp sigma peak imax jmax '
     .,'half_light_flux half_light_area flag flux2 SNR_F'
        if (nsource.gt.0) then
          do i=1,nsource
            write(10,*) (source_para(i,j),j=1,iSNR_F)
          enddo
          close(10)

          nn1=ns*len_g
          nn2=ns*(int(nsource/len_g)+1)
          call fq_chip_product_path(IMAGE_FILE(ichip),DIR_OUTPUT,
     .  DIR_SRC_P,
     .      '_source_p.fits',filename)
          call write_stamps(ngal_max,1,nsource,ns,ns
     .,power_coll,nn1,nn2,filename)
        else
          close(10)
        endif
      enddo

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
