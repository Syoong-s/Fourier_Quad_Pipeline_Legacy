      subroutine proc_info(iexpo)
      implicit none
      include 'para.inc'

      integer iexpo,nchip,i
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT
      real para(6)

      real expo_para(6,NMAX_EXPO)
      common /expo_para_pass/ expo_para

      call get_image_list(iexpo,IMAGE_FILE,nchip,DIR_OUTPUT)
      call get_expo_info(IMAGE_FILE,nchip,DIR_OUTPUT,para)

      do i=1,6
        expo_para(i,iexpo)=para(i)
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ==========================================
c Function: Aggregate exposure-level diagnostics
c Method: Read each chip WCS with a freshly initialized status flag
c ==========================================
      subroutine get_expo_info(IMAGE_FILE,nchip,DIR_OUTPUT,para)
      implicit none
      include 'para.inc'

      integer nchip
      real para(6)
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT
      character*(strl) PREFIX,fstar,fastro,fexpo
      integer ierror,i,j,ichip,nstar,nvalid
      double precision cRPIX(2),cD(2,2),cRVAL(2),PU(2,npd)
      real FWHM,e1,e2,chi_d
      real FWHM_AVE,chi_d_AVE,nstar_AVE,cRVAL1,cRVAL2

      call get_PREFIX_expo(IMAGE_FILE(1),PREFIX)
      fstar=trim(DIR_OUTPUT)
     .//'/stamps/'//trim(PREFIX)//'_star_info_expo.dat'
      fastro=trim(DIR_OUTPUT)
     .//'/astrometry/'//trim(PREFIX)//'.head'
      fexpo=trim(DIR_OUTPUT)
     .//'/result/'//trim(PREFIX)//'_expo_info.dat'

      open(unit=10,file=fexpo,status='replace')
      rewind 10
      write(10,*) '# ichip nstar FWHM e1 e2 chi_d '
     .,' cRPIX_1 cRPIX_2 cD_11 cD_12 cD_21 cD_22'
      open(unit=20,file=fstar,status='old')
      rewind 20
      read(20,*)

      nvalid=0
      FWHM_AVE=0.
      chi_d_AVE=0.
      nstar_AVE=0.
      cRVAL1=0.
      cRVAL2=0.

      do ichip=1,nchip
        read(20,*) i,nstar,FWHM,e1,e2,chi_d
        if (nstar.eq.0) then
          write(10,*) ichip,0,-99.,-99.,-99.,-99.,-99.
     .,-99.,-99.,-99.,-99.,-99.
          cycle
        endif
        nvalid=nvalid+1
        FWHM_AVE=FWHM_AVE+FWHM
        chi_d_AVE=chi_d_AVE+chi_d
        nstar_AVE=nstar_AVE+nstar
        ierror=0
        call read_astrometry_para(fastro,ichip         
     .,cRPIX,cD,cRVAL,PU,npd,ierror)
        cRVAL1=cRVAL(1)
        cRVAL2=cRVAL(2)
        write(10,*) ichip,nstar,FWHM,e1,e2,chi_d
     .,cRPIX(1),cRPIX(2),cD(1,1),cD(1,2),cD(2,1),cD(2,2)
      enddo
      close(20)
      close(10)

      if (nvalid.gt.0) then
        FWHM_AVE=FWHM_AVE/nvalid
        chi_d_AVE=chi_d_AVE/nvalid
        nstar_AVE=nstar_AVE/nvalid
      endif

      para(1)=nvalid
      para(2)=FWHM_AVE
      para(3)=chi_d_AVE
      para(4)=nstar_AVE
      para(5)=cRVAL1
      para(6)=cRVAL2

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
