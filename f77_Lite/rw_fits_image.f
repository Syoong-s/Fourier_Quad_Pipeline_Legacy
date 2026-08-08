      SUBROUTINE read_CCDNUM(filename,ccd_num)
! read the size of the 2D fits image from a fits file.
	  IMPLICIT NONE

      INTEGER status,unit,readwrite,blocksize,nfound
      INTEGER ccd_num
      CHARACTER filename*(*)
      character comment*100
      
      status=0
      
      CALL ftgiou(unit,status)

      readwrite=0
      CALL ftopen(unit,filename,readwrite,blocksize,status)
      CALL ftgkyj(unit,'CCDNUM',ccd_num,comment,status)

      CALL ftclos(unit, status)
      CALL ftfiou(unit, status)

      IF (status .gt. 0) CALL printerror(status)
      return
      END 
!*************************************************************************
      SUBROUTINE read_para(filename,nx,ny)
! read the size of the 2D fits image from a fits file.
	    IMPLICIT NONE

      INTEGER status,unit,readwrite,blocksize,nfound
      INTEGER naxes(2)
      INTEGER nx,ny
      CHARACTER filename*(*)
      
      status=0
      
      CALL ftgiou(unit,status)

      readwrite=0
      CALL ftopen(unit,filename,readwrite,blocksize,status)
      CALL ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)

      if (nfound .ne. 2)then
          print *,'READIMAGE failed to read the NAXISn keywords of:'
     .,filename	 
          return
      endif

      nx=naxes(1)
      ny=naxes(2)

      CALL ftclos(unit, status)
      CALL ftfiou(unit, status)

      IF (status .gt. 0) CALL printerror(status)
      return
      END 
!*************************************************************************
      subroutine read_stamps(np,nstart,n,nx,ny,stamp,n1,n2
     .,filename)
      implicit none

      integer np,nstart,n,nx,ny,n1,n2
      integer NMAX
      parameter (NMAX=7000)
      real stamp(np,nx,ny),large_stamp(NMAX,NMAX)
      character filename*(*)
      integer i,j,offx,offy,k

      if (n2.ge.NMAX) then
        write(*,*) '(Read Stamp)Large_stamp is too small!!'
        stop
      endif
      call readimage(filename,n1,n2,NMAX,NMAX,large_stamp)

      if (large_stamp(1,1) .lt. (-99990.0)) then
        write(*,*) 'Error opening large file:',filename
        stop
      endif

      offx=0	
      offy=0
      
      do k=nstart,n
        do i=1,nx
          do j=1,ny
            stamp(k,i,j)=large_stamp(i+offx,j+offy)
          enddo
        enddo
        offx=offx+nx
        if (offx+nx.gt.n1) then
          offx=0
          offy=offy+ny
        endif
      enddo


      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine write_stamps_2(np,n,nx,ny,stamp,opt,val,n1,n2
     .,filename)
      implicit none

      integer np,n,nx,ny,n1,n2
      integer NMAX,opt(np),val
      parameter (NMAX=7000)
      real stamp(np,nx,ny),large_stamp(NMAX,NMAX)
      character filename*(*)
      integer i,j,offx,offy,k

      do i=1,n1
        do j=1,n2
          large_stamp(i,j)=0.
        enddo
      enddo

      offx=0	
      offy=0
      
      do k=1,n
              if (offy+ny.gt.n2) then
                write(*,*) 'large_stamp is too small!!'
                read(*,*)
              endif
              if (opt(k).ne.val) cycle
        do i=1,nx
          do j=1,ny
            large_stamp(i+offx,j+offy)=stamp(k,i,j)
          enddo
        enddo
        offx=offx+nx
        if (offx+nx.gt.n1) then
          offx=0
          offy=offy+ny
        endif
      enddo

            call writeimage(filename,n1,n2,NMAX,NMAX,large_stamp)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine write_stamps(np,nstart,n,nx,ny,stamp,n1,n2
     .,filename)
      implicit none

      integer np,nstart,n,nx,ny,n1,n2
      integer NMAX
      parameter (NMAX=7000)
      real stamp(np,nx,ny),large_stamp(NMAX,NMAX)
      character filename*(*)
      integer i,j,offx,offy,k

      if (n2.gt.NMAX) then
        write(*,*) '(Write Stamp)Large_stamp is too small!!'
        stop
      endif

      do i=1,n1
        do j=1,n2
          large_stamp(i,j)=0.
        enddo
      enddo

      offx=0	
      offy=0
      
      do k=nstart,n
        if (offy+ny.gt.n2) then
          write(*,*) 'write_stamps ERROR!!'
          stop
        endif
        do i=1,nx
          do j=1,ny
            large_stamp(i+offx,j+offy)=stamp(k,i,j)
          enddo
        enddo
        offx=offx+nx
        if (offx+nx.gt.n1) then
          offx=0
          offy=offy+ny
        endif
      enddo

      call writeimage(filename,n1,n2,NMAX,NMAX,large_stamp)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This is a simple module contains subroutine that can read and write 2D fits images
! with no extensions to the primary array. It is stored in single precision.
!IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE readimage_para(filename,nx,ny,npx,npy,array
     .,CRPIX,CD,CRVAL)
! read a 2D fits images from a fits file.
	  IMPLICIT NONE

      INTEGER status,unit,readwrite,blocksize,nfound
      INTEGER group,firstpix,nbuffer,npixels,i,status2,j
      INTEGER naxes(2)
      INTEGER nx,ny,npx,npy
      double precision CRPIX(2),CD(2,2),CRVAL(2),temp(2)
      REAL array(npx,npy)
      REAL nullval,anyf
      LOGICAL anynull
      CHARACTER filename*(*),comment*100
      
!  The STATUS parameter must always be initialized.
      status=0
      nullval=0.
      
!  Get an unused Logical Unit Number to use to open the FITS file.
      CALL ftgiou(unit,status)

      readwrite=0
      CALL ftopen(unit,filename,readwrite,blocksize,status)
      IF (status .gt. 0) THEN
!         仅输出一行：包含状态码与文件名。使用 TRIM 消除文件名后的多余空格
          PRINT *, 'Error Status=', status, ' | File=', TRIM(filename)
          array(1,1)=-99999.0
          
!         清理内存并关闭可能残留的句柄
          CALL ftclos(unit, status)
          CALL ftfiou(unit, status)
          
!         将 status 强制置为 0，防止末尾的 CALL printerror 打印长篇堆栈信息
          status = 0
          GOTO 999 
      END IF


!  Determine the coordinate parameters of the image.
      CALL ftgknd(unit,'CRPIX',1,2,temp,nfound,status)
!  Check that it found both CRPIX1 and CRPIX2 keywords.
      IF (nfound .ne. 2 .or. status .gt. 0) THEN
          PRINT *, 'Error Status=', status, ' | Missing CRPIXn | File='
     ., TRIM(filename)
          array(1,1) = -99999.0
          
!         【关键】必须在退出前关闭文件并释放 Unit，防止句柄泄漏
          CALL ftclos(unit, status)
          CALL ftfiou(unit, status)
          
!         强制清零 status 并跳转到末尾，静默 CFITSIO 的多行报错堆栈
          status = 0
          GOTO 999 
      END IF
      CRPIX(1)=temp(1)
      CRPIX(2)=temp(2)

      CALL ftgknd(unit,'CRVAL',1,2,temp,nfound,status)
!  Check that it found both CRVAL1 and CRVAL2 keywords.
      IF (nfound .ne. 2 .or. status .gt. 0) THEN
          PRINT *, 'Error Status=', status, ' | Missing CRVAL | File='
     ., TRIM(filename)
          array(1,1) = -99999.0
          
!         【关键】必须在退出前关闭文件并释放 Unit，防止句柄泄漏
          CALL ftclos(unit, status)
          CALL ftfiou(unit, status)
          
!         强制清零 status 并跳转到末尾，静默 CFITSIO 的多行报错堆栈
          status = 0
          GOTO 999 
      END IF
      CRVAL(1)=temp(1)
      CRVAL(2)=temp(2)

!      CALL ftgknd(unit,'CDELT',1,2,temp,nfound,status)
!  Check that it found both CDELT1 and CDELT2 keywords.
!      if (nfound .ne. 2)then
!          print *,'READIMAGE failed to read the CDELTn keywords.'
!          return
!      end if
!      CDELT(1)=temp(1)
!      CDELT(2)=temp(2)

      CALL ftgknd(unit,'CD1_',1,2,temp,nfound,status)
!  Check that it found both CD1_1 and CD1_2 keywords.
      IF (nfound .ne. 2 .or. status .gt. 0) THEN
          PRINT *, 'Error Status=', status, ' | Missing CD1_ | File='
     ., TRIM(filename)
          array(1,1) = -99999.0
          
!         【关键】必须在退出前关闭文件并释放 Unit，防止句柄泄漏
          CALL ftclos(unit, status)
          CALL ftfiou(unit, status)
          
!         强制清零 status 并跳转到末尾，静默 CFITSIO 的多行报错堆栈
          status = 0
          GOTO 999 
      END IF
      CD(1,1)=temp(1)
      CD(1,2)=temp(2)

      CALL ftgknd(unit,'CD2_',1,2,temp,nfound,status)
!  Check that it found both CD2_1 and CD2_2 keywords.
      IF (nfound .ne. 2 .or. status .gt. 0) THEN
          PRINT *, 'Error Status=', status, ' | Missing CD2_ | File='
     ., TRIM(filename)
          array(1,1) = -99999.0
          
!         【关键】必须在退出前关闭文件并释放 Unit，防止句柄泄漏
          CALL ftclos(unit, status)
          CALL ftfiou(unit, status)
          
!         强制清零 status 并跳转到末尾，静默 CFITSIO 的多行报错堆栈
          status = 0
          GOTO 999 
      END IF
      CD(2,1)=temp(1)
      CD(2,2)=temp(2)

c      CALL ftgkyd(unit,'ZP',ZP,comment,status)

c      CALL ftgkyd(unit,'EXPTIME',EXPTIME,comment,status)
      

!  Determine the size of the image.
      CALL ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)

!  Check that it found both NAXIS1 and NAXIS2 keywords.
      IF (nfound .ne. 2 .or. status .gt. 0) THEN
          PRINT *, 'Error Status=', status, ' | Missing NAXIS | File='
     ., TRIM(filename)
          array(1,1) = -99999.0
          
!         【关键】必须在退出前关闭文件并释放 Unit，防止句柄泄漏
          CALL ftclos(unit, status)
          CALL ftfiou(unit, status)
          
!         强制清零 status 并跳转到末尾，静默 CFITSIO 的多行报错堆栈
          status = 0
          GOTO 999 
      END IF

      group=1
      nx=naxes(1)
      ny=naxes(2)

      status2=0
	  
      CALL FTG2DE(unit,group,nullval,npx,nx,ny,array,anyf,status)     

!  The FITS file must always be closed before exiting the program. 
!  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      CALL ftclos(unit, status)
      CALL ftfiou(unit, status)
999   CONTINUE
!  Check for any error, and if so print out error messages.
!  The PRINTERROR subroutine is listed near the end of this file.
      IF (status .gt. 0) CALL printerror(status)
      return
      END 
! *************************************************************************
      SUBROUTINE update_para(filename,CRPIX,CD)
      IMPLICIT NONE

      INTEGER status,unit,readwrite,blocksize
      double precision CRPIX(2),CD(2,2)
      CHARACTER filename*(*),comment*100
      
      status=0
      comment='replaced'
      CALL ftgiou(unit,status)
      readwrite=1
      CALL ftopen(unit,filename,readwrite,blocksize,status)
  
      call ftukyd(unit,'CRPIX1',CRPIX(1),9,comment,status)
      call ftukyd(unit,'CRPIX2',CRPIX(2),9,comment,status)

      call ftukyd(unit,'CD1_1',CD(1,1),9,comment,status)
      call ftukyd(unit,'CD1_2',CD(1,2),9,comment,status)
      call ftukyd(unit,'CD2_1',CD(2,1),9,comment,status)
      call ftukyd(unit,'CD2_2',CD(2,2),9,comment,status)
	

      CALL ftclos(unit, status)
      CALL ftfiou(unit, status)

      IF (status .gt. 0) CALL printerror(status)
      return
      END 
! *************************************************************************
      SUBROUTINE writeimage_copyhdu(file1,filename,nx,ny,npx,npy,array)

!  Create a FITS primary array containing a 2-D image
      IMPLICIT NONE
      
      INTEGER status,unit,blocksize,bitpix,naxis,readwrite
      INTEGER i,j,group,unit1
      INTEGER nx,ny,npx,npy
      INTEGER naxes(2)
      REAL array(npx,npy)
      CHARACTER filename*(*),file1*(*)      
      LOGICAL simple,extend

      status=0      

!  Delete the file if it already exists, so we can then recreate it.
!  The deletefile subroutine is listed at the end of this file.
      CALL deletefile(filename,status)

!  Get an unused Logical Unit Number to use to open the FITS file.
!  This routine is not required;  programmers can choose any unused
!  unit number to open the file.
      CALL ftgiou(unit,status)

!  Create the new empty FITS file.  The blocksize parameter is a
!  historical artifact and the value is ignored by FITSIO.
      blocksize=1
      CALL ftinit(unit,filename,blocksize,status)

      status=0
      CALL ftgiou(unit1,status)
      readwrite=0
      CALL ftopen(unit1,file1,readwrite,blocksize,status)
      CALL ftcphd(unit1,unit,status)

      CALL ftmkyj(unit,'bitpix',-32,'new',status)	
c      CALL ftmkyd(unit,'bscale',1d0,'new',status)
c      CALL ftmkyd(unit,'bzero',0d0,'new',status)	

!  Initialize parameters about the FITS image.
!  The size of the image is given by the NAXES values. 
!  The EXTEND = TRUE parameter indicates that the FITS file
!  may contain extensions following the primary array.

c	simple=.true.
c      bitpix=-32
c      naxis=2
c      naxes(1)=nx
c      naxes(2)=ny
c      extend=.false.

!  Write the required header keywords to the file
c      CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

!  Write the array to the FITS file.
!  The last letter of the subroutine name defines the datatype of the
!  array argument; in this case the 'J' indicates that the array has an
!  integer*4 datatype. ('I' = I*2, 'E' = Real*4, 'D' = Real*8).
!  The 2D array is treated as a single 1-D array with NAXIS1 * NAXIS2
!  total number of pixels.  GROUP is seldom used parameter that should
!  almost always be set = 1.
      group=1
      CALL FTP2DE(unit,group,npx,nx,ny,array,status)

!  Write another optional keyword to the header
!  The keyword record will look like this in the FITS file:
!
!  EXPOSURE=                 1500 / Total Exposure Time
!
!     CALL ftpkyj(unit,'EXPOSURE',1500,'Total Exposure Time',status)


!  The FITS file must always be closed before exiting the program. 
!  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      CALL ftclos(unit, status)
      CALL ftfiou(unit, status)
      CALL ftclos(unit1, status)
      CALL ftfiou(unit1, status)

!  Check for any errors, and if so print out error messages.
!  The PRINTERROR subroutine is listed near the end of this file.
      IF (status .gt. 0) CALL printerror(status)
      return
      END
! *************************************************************************
      SUBROUTINE writeimage(filename,nx,ny,npx,npy,array)

!  Create a FITS primary array containing a 2-D image
      IMPLICIT NONE
      
      INTEGER status,unit,blocksize,bitpix,naxis
      INTEGER i,j,group
      INTEGER nx,ny,npx,npy
      INTEGER naxes(2)
      REAL array(npx,npy)
      CHARACTER filename*(*)      
      LOGICAL simple,extend
      
!  The STATUS parameter must be initialized before using FITSIO.  A
!  positive value of STATUS is returned whenever a serious error occurs.
!  FITSIO uses an `inherited status' convention, which means that if a
!  subroutine is called with a positive input value of STATUS, then the
!  subroutine will exit immediately, preserving the status value. For 
!  simplicity, this program only checks the status value at the end of 
!  the program, but it is usually better practice to check the status 
!  value more frequently.

      status=0

!  Delete the file if it already exists, so we can then recreate it.
!  The deletefile subroutine is listed at the end of this file.
      CALL deletefile(filename,status)

!  Get an unused Logical Unit Number to use to open the FITS file.
!  This routine is not required;  programmers can choose any unused
!  unit number to open the file.
      CALL ftgiou(unit,status)

!  Create the new empty FITS file.  The blocksize parameter is a
!  historical artifact and the value is ignored by FITSIO.
      blocksize=1
      CALL ftinit(unit,'!'//trim(filename),blocksize,status)

!  Initialize parameters about the FITS image.
!  The size of the image is given by the NAXES values. 
!  The EXTEND = TRUE parameter indicates that the FITS file
!  may contain extensions following the primary array.
      simple=.true.
      bitpix=-32
      naxis=2
      naxes(1)=nx
      naxes(2)=ny
      extend=.false.

!  Write the required header keywords to the file
      CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

!  Write the array to the FITS file.
!  The last letter of the subroutine name defines the datatype of the
!  array argument; in this case the 'J' indicates that the array has an
!  integer*4 datatype. ('I' = I*2, 'E' = Real*4, 'D' = Real*8).
!  The 2D array is treated as a single 1-D array with NAXIS1 * NAXIS2
!  total number of pixels.  GROUP is seldom used parameter that should
!  almost always be set = 1.
      group=1
      CALL FTP2DE(unit,group,npx,nx,ny,array,status)

!  Write another optional keyword to the header
!  The keyword record will look like this in the FITS file:
!
!  EXPOSURE=                 1500 / Total Exposure Time
!
!     CALL ftpkyj(unit,'EXPOSURE',1500,'Total Exposure Time',status)

!  The FITS file must always be closed before exiting the program. 
!  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      CALL ftclos(unit, status)
      CALL ftfiou(unit, status)

!  Check for any errors, and if so print out error messages.
!  The PRINTERROR subroutine is listed near the end of this file.
      IF (status .gt. 0) CALL printerror(status)
      return
      END
! *************************************************************************
      SUBROUTINE readimage(filename,nx,ny,npx,npy,array)
! read a 2D fits images from a fits file.
	  IMPLICIT NONE

      INTEGER status,unit,readwrite,blocksize,nfound
      INTEGER group,firstpix,nbuffer,npixels,i,status2
      INTEGER naxes(2)
      INTEGER nx,ny,npx,npy
      REAL array(npx,npy)
      REAL nullval,anyf
      LOGICAL anynull
      CHARACTER filename*(*)
      
!  The STATUS parameter must always be initialized.
      status=0
      nullval=0.
      
!  Get an unused Logical Unit Number to use to open the FITS file.
      CALL ftgiou(unit,status)

      readwrite=0
      CALL ftopen(unit,filename,readwrite,blocksize,status)
      IF (status .gt. 0) THEN
!         仅输出一行：包含状态码与文件名。使用 TRIM 消除文件名后的多余空格
          PRINT *, 'Error Status=', status, ' | File=', TRIM(filename)
          array(1,1)=-99999.0
          
!         清理内存并关闭可能残留的句柄
          CALL ftclos(unit, status)
          CALL ftfiou(unit, status)
          
!         【关键】将 status 强制置为 0，防止末尾的 CALL printerror 打印长篇堆栈信息
          status = 0
          GOTO 999 
      END IF

!  Determine the size of the image.
      CALL ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)

!  Check that it found both NAXIS1 and NAXIS2 keywords.
      IF (nfound .ne. 2 .or. status .gt. 0) THEN
          PRINT *, 'Error Status=', status, ' | Missing NAXIS | File='
     ., TRIM(filename)
          array(1,1) = -99999.0
          
!         【关键】必须在退出前关闭文件并释放 Unit，防止句柄泄漏
          CALL ftclos(unit, status)
          CALL ftfiou(unit, status)
          
!         强制清零 status 并跳转到末尾，静默 CFITSIO 的多行报错堆栈
          status = 0
          GOTO 999 
      END IF

      group=1
      nx=naxes(1)
      ny=naxes(2)

      status2=0
	  
      CALL FTG2DE(unit,group,nullval,npx,nx,ny,array,anyf,status)

!  The FITS file must always be closed before exiting the program. 
!  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      CALL ftclos(unit, status)
      CALL ftfiou(unit, status)
999   CONTINUE
!  Check for any error, and if so print out error messages.
!  The PRINTERROR subroutine is listed near the end of this file.
      IF (status .gt. 0) CALL printerror(status)
      return
      END 
! *************************************************************************
      SUBROUTINE printerror(status)

!  This subroutine prints out the descriptive text corresponding to the
!  error status value and prints out the contents of the internal
!  error message stack generated by FITSIO whenever an error occurs.
      IMPLICIT NONE
      INTEGER status
      CHARACTER errtext*30
      CHARACTER errmessage*100

!  Check if status is OK (no error); if so, simply return
      IF (status .le. 0) RETURN

!  The FTGERR subroutine returns a descriptive 30-character text string that
!  corresponds to the integer error status number.  A complete list of all
!  the error numbers can be found in the back of the FITSIO User's Guide.
      CALL ftgerr(status,errtext)
      WRITE(*,*) 'FITSIO Error Status =',status,': ',errtext

!  FITSIO usually generates an internal stack of error messages whenever
!  an error occurs.  These messages provide much more information on the
!  cause of the problem than can be provided by the single integer error
!  status value.  The FTGMSG subroutine retrieves the oldest message from
!  the stack and shifts any remaining messages on the stack down one
!  position.  FTGMSG is called repeatedly until a blank message is
!  returned, which indicates that the stack is empty.  Each error message
!  may be up to 100 characters in length.  Another subroutine, called
!  FTCMSG, is available to simply clear the whole error message stack in
!  cases where one is not interested in the contents.
      CALL ftgmsg(errmessage)
      DO WHILE (errmessage .ne. ' ')
          WRITE(*,*) errmessage
          CALL ftgmsg(errmessage)
      END DO
      return 
      END 
! *************************************************************************
      SUBROUTINE deletefile(filename,status)

!  A simple little routine to delete a FITS file
      IMPLICIT NONE
      
      INTEGER status,unit,blocksize
      CHARACTER*(*) filename

!  Simply return if status is greater than zero
      IF (status .gt. 0) RETURN

!  Get an unused Logical Unit Number to use to open the FITS file
      CALL ftgiou(unit,status)

!  Try to open the file, to see if it exists
      CALL ftopen(unit,filename,1,blocksize,status)

      IF (status .eq. 0) THEN
!         file was opened;  so now delete it 
          CALL ftdelt(unit,status)
      ELSE IF (status .eq. 103) THEN
!         file doesn't exist, so just reset status to zero and clear errors
          status=0
          CALL ftcmsg
      ELSE
!         there was some other error opening the file; delete the file anyway
          status=0
          CALL ftcmsg
          CALL ftdelt(unit,status)
      END IF

!  Free the unit number for later reuse
      CALL ftfiou(unit, status)
      return
      END 
!*****************************************************************************
      subroutine write_copyhdu(infilename,outfilename
     .,nx,ny,npx,npy,array)

C     copy the 1st and 3rd HDUs from the input file to a new FITS file

      integer status,inunit,outunit,readwrite,blocksize,morekeys,hdutype
      character*(*) infilename,outfilename

      INTEGER group
      INTEGER nx,ny,npx,npy
      REAL array(npx,npy)

1     status=0

C     Delete the file if it already exists, so we can then recreate it
2     call deletefile(outfilename,status)

C     Get  unused Logical Unit Numbers to use to open the FITS files
3     call ftgiou(inunit,status)
      call ftgiou(outunit,status)

C     open the input FITS file, with readonly access
      readwrite=0
4     call ftopen(inunit,infilename,readwrite,blocksize,status)

C     create the new empty FITS file with the standard block size
      blocksize=1
5     call ftinit(outunit,outfilename,blocksize,status)

C     copy the primary array from the input file to the output file
      morekeys=0
6     call ftcopy(inunit,outunit,morekeys,status)

      group=1
      CALL FTP2DE(outunit,group,npx,nx,ny,array,status)

C     append/create a new empty extension on the end of the output file
c 7    call ftcrhd(outunit,status)
C     skip to the 3rd extension in the input file
c 8    call ftmahd(inunit,3,hdutype,status)
C     copy this extension from the input file to the output file
c 9    call ftcopy(inunit,outunit,morekeys,status)  

C     close the FITS file and free the unit numbers
10    call ftclos(inunit, status)
      call ftclos(outunit, status)
11    call ftfiou(-1, status)

C     check for any error, and if so print out error messages
12    if (status .gt. 0)call printerror(status)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine write_serial_stamps(filename,nx,ny,npx,npy,array
     .,sign_order)

C     copy a series of stamp images to a .fits file.

      integer status,unit,readwrite,blocksize,morekeys,hdutype
      save status,unit
      character filename*(*)

      INTEGER group,sign_order
      INTEGER nx,ny,npx,npy
      REAL array(npx,npy)

      INTEGER bitpix,naxis
      INTEGER i,j
      INTEGER naxes(2)
c      LOGICAL simple,extend
	
      character item_name*8,item_explanation*24
      integer item_value
      common /item_pass/ item_name,item_value,item_explanation

        bitpix=-32
        naxis=2
        group=1

      if (sign_order.eq.0) then
        status=0
        call deletefile(filename,status)
        call ftgiou(unit,status)
        blocksize=1
5      call ftinit(unit,filename,blocksize,status)

      elseif (sign_order.eq.1) then
        naxes(1)=nx
        naxes(2)=ny
        CALL ftphps(unit,bitpix,naxis,naxes,status)
        CALL FTP2DE(unit,group,npx,nx,ny,array,status)

      elseif (sign_order.eq.2) then
        call ftcrhd(unit,status)
        naxes(1)=nx
        naxes(2)=ny
        CALL ftphps(unit,bitpix,naxis,naxes,status)
        CALL FTP2DE(unit,group,npx,nx,ny,array,status)

      elseif (sign_order.eq.3) then
        call ftpkyj(unit,item_name,item_value,item_explanation,status)

      elseif (sign_order.eq.-1) then
        call ftclos(unit, status)
        call ftfiou(-1, status)

      endif

C     check for any error, and if so print out error messages
12    if (status .gt. 0)call printerror(status)
      return
      end
