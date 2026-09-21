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
c ==========================================
c Function: Refine and validate exposure astrometry
c Method: Consume aligned chip data and write head and check products
c ==========================================
      subroutine chip_process_astrometry(IMAGE_FILE,nchip,DIR_OUTPUT)
      implicit none
      include 'para.inc'
      include 'path_layout.inc'

      integer nchip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT
      integer proc_error,ichip
      character PREFIX*(strl),filename*(strl)
      double precision cRPIX(2),cD(2,2),PU(2,npd),cRVAL(2)
      double precision cRPIX_b(2),cD_b(2,2),PU_b(2,npd),cRVAL_b(2)
      integer n,i,j,k
      double precision ra,dec,x,y,ra2,dec2

      if (ASTROMETRY_trivial.eq.1) then
        call get_astrometry_trivial(IMAGE_FILE,nchip,DIR_OUTPUT)
      else
        call get_astrometry(IMAGE_FILE,nchip,DIR_OUTPUT)
c        if (ext_cat.eq.1) call get_astrometry_gal(IMAGE_FILE,nchip
c     .,DIR_OUTPUT)
      endif

      if (ASTROMETRY_trivial.eq.1) return

      call fq_expo_product_path(IMAGE_FILE(1),DIR_OUTPUT,
     .  DIR_ASTRO_CHECK,
     .  '_check.dat',filename)

      open(unit=50,file=trim(filename),status='replace')
      rewind 50

      do ichip=1,nchip
        proc_error=0

        call fq_expo_product_path(IMAGE_FILE(1),DIR_OUTPUT,
     .  DIR_ASTRO_HEAD,
     .    '.head',filename)

        call read_astrometry_para(filename,ichip
     .,cRPIX,cD,cRVAL,PU,npd,proc_error)

        if (proc_error.eq.0) then

          call fq_chip_product_path(IMAGE_FILE(ichip),DIR_OUTPUT,
     .  DIR_NORM,
     .      '_norm.fits',filename)
          call update_para(filename,cRPIX,cD)

          call fq_chip_product_path(IMAGE_FILE(ichip),DIR_OUTPUT,
     .      DIR_ASTRO_DATA,'_astro.dat',filename)
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
c ==========================================
c Function: Assemble the trivial exposure astrometry head
c Method: Read aligned chip astrometry and write the exposure product
c ==========================================
      subroutine get_astrometry_trivial(IMAGE_FILE,nchip,DIR_OUTPUT)
      implicit none
      include 'para.inc'
      include 'path_layout.inc'

      integer nchip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT

      character PREFIX*(strl),filename*(strl)
      integer i,valid,ichip,k

      double precision cRPIX2(nchip,2),cD2(nchip,2,2)
      double precision PU(2,npd),cRVAL2(2)

      do ichip=1,nchip

        call fq_chip_product_path(IMAGE_FILE(ichip),DIR_OUTPUT,
     .    DIR_ASTRO_DATA,'_astro.dat',filename)

        open(unit=10,file=filename,status='old')
        rewind 10
        read(10,*) cRPIX2(ichip,1),cRPIX2(ichip,2)
     .,cRVAL2(1),cRVAL2(2)
        read(10,*) cD2(ichip,1,1),cD2(ichip,1,2)
     .,cD2(ichip,2,1),cD2(ichip,2,2)
        close(10)

      enddo

      do i=1,npd
        PU(1,i)=0d0
        PU(2,i)=0d0
      enddo

      valid=1

      call fq_expo_product_path(IMAGE_FILE(1),DIR_OUTPUT,
     .  DIR_ASTRO_HEAD,
     .  '.head',filename)
      open(unit=10,file=filename,status='replace')
      rewind 10
      write(10,*) cRVAL2(1),cRVAL2(2)
      do i=1,npd
        write(10,*) PU(1,i),PU(2,i)
      enddo
      do k=1,nchip
        write(10,*) k,valid,cRPIX2(k,1),cRPIX2(k,2)
     .,cD2(k,1,1),cD2(k,1,2),cD2(k,2,1),cD2(k,2,2)
      enddo
      close(10)


      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ==========================================
c Function: Fit the exposure astrometry solution
c Method: Read aligned chip matches and write the exposure head product
c ==========================================
      subroutine get_astrometry(IMAGE_FILE,nchip,DIR_OUTPUT)
      implicit none
      include 'para.inc'
      include 'path_layout.inc'

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

        call fq_chip_product_path(IMAGE_FILE(ichip),DIR_OUTPUT,
     .    DIR_ASTRO_DATA,'_astro.dat',filename)

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


      call fq_expo_product_path(IMAGE_FILE(1),DIR_OUTPUT,
     .  DIR_ASTRO_HEAD,
     .  '.head',filename)
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
      subroutine gen_astrometry_data_gal(cat_gal,nx,ny,npx,npy
     .,map,sigmap,cRPIX,cD,cRVAL,filename,proc_error,saturation_thresh)
      implicit none

      real astrometry_shift_ratio
      parameter (astrometry_shift_ratio=0.3)

      character*(*) cat_gal,filename
      integer nx,ny,npx,npy,ierror,n_user,n_ref,proc_error
      real map(npx,npy),sigmap(npx,npy),flux_min,saturation_thresh

      integer i,j,u,v,k
      double precision ra,dra,diffra,dec(2),a,d,flux,mag,tmp
      integer nss_max,nss
      parameter (nss_max=10000)
      double precision xs(nss_max),ys(nss_max)
      double precision xr(nss_max),yr(nss_max)
      double precision ra_r(nss_max),dec_r(nss_max)
      integer box(nss_max),astrometry_shift_range
      double precision ra2(nss_max),dec2(nss_max)
      double precision x2(nss_max),y2(nss_max)

      double precision cRPIX(2),cD(2,2),cRVAL(2)


      if (proc_error.eq.1) then
        nss=0
        n_user=0
        n_ref=0
        goto 40
      endif


c Get the ranges of ra & dec:
      call get_ra_dec_range(nx,ny,ra,dec,dra,cRPIX,cD,cRVAL
     .,astrometry_shift_ratio)

c------------------------------------------------------------------
c Get the positions of the stars in the reference catalog:

      n_ref=0
      open(unit=20,file=trim(cat_gal),status='old',iostat=ierror)
      rewind 20
      if (ierror.ne.0) then
        write(*,'(A)') cat_gal
        write(*,*) ierror
        write(*,*) 'catalog file error!!'
        read(*,*)
      endif
      read(20,*)
      do while (ierror.ge.0)
        read(20,*,iostat=ierror) a,d
        if (ierror.lt.0) cycle
        if (abs(diffra(a,ra)).gt.dra*0.5) cycle
        if (d.lt.dec(1).or.d.gt.dec(2)) cycle
        n_ref=n_ref+1
        ra_r(n_ref)=a
        dec_r(n_ref)=d
        call coordinate_transfer_simple(a,d,xr(n_ref),yr(n_ref),-1
     .,cRPIX,cD,cRVAL)
      enddo
      close(20)


      call get_astrometry_catalog(nx,ny,npx,npy,map
     .,sigmap,n_ref,nss_max,nss,xs,ys,saturation_thresh*1.5)


      n_user=nss

c Match the point sources in the "user" and "ref" catalogs:
      astrometry_shift_range=int(max(nx,ny)*astrometry_shift_ratio)

      call pattern_matching(nss_max,n_ref,xr,yr,nss_max,n_user,xs,ys
     .,astrometry_shift_range,box)

c Get the Astrometric calibration parameters:


      nss=0
      do i=1,n_ref
        if (box(i).eq.0) cycle
        nss=nss+1
        ra2(nss)=ra_r(i)
        dec2(nss)=dec_r(i)
        j=box(i)
        x2(nss)=xs(j)
        y2(nss)=ys(j)
      enddo

      write(*,*) nss,n_ref,n_user,trim(filename)

40    open(unit=10,file=trim(filename),status='replace')
      rewind 10
      write(10,*) cRPIX(1),cRPIX(2),cRVAL(1),cRVAL(2)
      write(10,*) cD(1,1),cD(1,2),cD(2,1),cD(2,2)
      write(10,*) nss,n_user,n_ref
      do i=1,nss
        write(10,*) ra2(i),dec2(i),x2(i),y2(i)
      enddo
      close(10)


      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ==========================================
c Function: Fit the optional galaxy-based astrometry solution
c Method: Read aligned chip matches and write the alternate head product
c ==========================================
      subroutine get_astrometry_gal(IMAGE_FILE,nchip,DIR_OUTPUT)
      implicit none
      include 'para.inc'
      include 'path_layout.inc'

      integer nchip,npchip
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

        call fq_chip_product_path(IMAGE_FILE(ichip),DIR_OUTPUT,
     .    DIR_ASTRO_DATA,'_astro_gal.dat',filename)

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
      else
        do i=1,nchip
          valid(i)=0
        enddo
      endif


      call fq_expo_product_path(IMAGE_FILE(1),DIR_OUTPUT,
     .  DIR_ASTRO_HEAD,
     .  '_gal.head',filename)
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
