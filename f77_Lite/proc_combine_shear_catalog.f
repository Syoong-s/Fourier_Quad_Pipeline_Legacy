c ==========================================
c Function: Combine one exposure's calibrated shear catalog
c Method: Pass its one-based exposure index into the Stage-9 combiner
c ==========================================
      subroutine proc_comb(iexpo)
      implicit none
      include 'para.inc'

      integer iexpo,nchip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT
      real chi2
      real expo_para(6,NMAX_EXPO)
      common /expo_para_pass/ expo_para

      call get_image_list(iexpo,IMAGE_FILE,nchip,DIR_OUTPUT)

      chi2=expo_para(3,iexpo)
      call combine_expo_catalog(nchip,IMAGE_FILE,DIR_OUTPUT,iexpo,chi2)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ==========================================
c Function: Write the final C++-compatible exposure catalog
c Method: Join aligned rows with exposure and physical CCD identifiers
c ==========================================
      subroutine combine_expo_catalog(nchip,IMAGE_FILE
     .,DIR_OUTPUT,iexpo,chi2)
      implicit none
      include 'para.inc'
      include 'path_layout.inc'

      integer nchip,iexpo
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT

      real chi2
      character*(strl) PREFIX,filename
      character*1000 cat_content,cat_list1,cat_list2
      integer ichip,ierror,u,i,m,n,chip_index
      real cat(iparity),g1c,g2c


      call fq_expo_product_path(IMAGE_FILE(1),DIR_OUTPUT,
     .  DIR_RESULT,
     .  '_all.cat',filename)
      open(unit=20,file=filename,status='replace',iostat=ierror)
      rewind 20

      do ichip=1,nchip
        call fq_chip_product_path(IMAGE_FILE(ichip),DIR_OUTPUT,
     .  DIR_CAT_ORIG,
     .    '_orig.cat',filename)
        open(unit=15,file=filename,status='old',iostat=ierror)
        rewind 15
        if (ierror.ne.0) then
          write(*,*) trim(filename),' is missing!'
          stop
        endif
        read(15,'(A)') cat_list2
        read(15,'(A)',iostat=ierror)  cat_content
        if (ierror.lt.0) then
          close(15)
          cycle
        else
          close(15)
          exit
        endif
      enddo

      n=0
      m=0
      do ichip=1,nchip
        call get_chip_id(IMAGE_FILE(ichip),chip_index)
        call get_PREFIX(IMAGE_FILE(ichip),PREFIX)
        call fq_chip_product_path(IMAGE_FILE(ichip),DIR_OUTPUT,
     .  DIR_SHEAR,
     .    '_shear.dat',filename)
        open(unit=10,file=filename,status='old',iostat=ierror)
        rewind 10
        if (ierror.ne.0) then
          write(*,*) trim(filename),' is missing!'
          stop
        endif
        read(10,'(A)') cat_list1

        call fq_chip_product_path(IMAGE_FILE(ichip),DIR_OUTPUT,
     .  DIR_CAT_ORIG,
     .    '_orig.cat',filename)
        open(unit=15,file=filename,status='old',iostat=ierror)
        rewind 15
        if (ierror.ne.0) then
          write(*,*) trim(filename),' is missing!'
          stop
        endif
        read(15,*)
        if (ichip.eq.1) then
          write(20,*) trim(cat_list2),' EXPO_NUM ccD_NUM ',
     .      trim(cat_list1),' Chi2'
          if (chi2.gt.chi2_thresh) then
            close(10)
            close(15)
            close(20)
            write(*,*) trim(PREFIX)//' contains no valid sources!'
            return
          endif
        endif
        do while (ierror.ge.0)
          read(10,*,iostat=ierror) (cat(u),u=1,iparity)
          if (ierror.lt.0) cycle
          read(15,'(A)') cat_content
          if (cat(i_imax).ge.ns.or.cat(i_jmax).ge.ns) then
            m=m+1
            cycle
          endif
          if (cat(1) .lt. -900.0) then
            m=m+1
            cycle
          endif
          n=n+1
          ! g1c=cat(igf1)+g1_c
          ! g2c=cat(igf2)+g2_c
          g1c = 0.
          g2c = 0.
          cat(ig1)=cat(ig1)-g1c*cat(ide)+g1c*cat(ih1)+g2c*cat(ih2)
          cat(ig2)=cat(ig2)-g2c*cat(ide)+g1c*cat(ih2)-g2c*cat(ih1)
          write(20,*) trim(cat_content),iexpo,chip_index
     .,(cat(u),u=1,iparity),chi2
        enddo
        close(15)
        close(10)
      enddo
      write(*,*) trim(PREFIX),n,m
      close(20)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
