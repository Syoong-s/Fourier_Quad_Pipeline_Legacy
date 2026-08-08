      subroutine proc_source(iexpo)
      implicit none
      include 'para.inc'

      integer iexpo,nchip,ichip,i
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT,FLAT_FILE
      character*2 chip_name
      call get_image_list(iexpo,IMAGE_FILE,nchip,DIR_OUTPUT)

      do ichip=1,nchip
c        call get_flatname(IMAGE_FILE(ichip),FLAT_FILE)
        ! call get_chip_id(IMAGE_FILE(ichip),i)
        ! write(chip_name,'(I2.2)') i
    !     FLAT_FILE=trim(FLAT_PATH)//'/flat_'
    !  .//trim(chip_name)//'_weight.fits'
        call chip_process_source(IMAGE_FILE,ichip,DIR_OUTPUT,FLAT_FILE)
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_flatname(imagefile,flatfile)
      implicit none
      include 'para.inc'

      character*(*) imagefile,flatfile

      integer i,n,id
      character*2 cc
      character b
      character*4 time,tttt
      integer ym

20    format(I2.2)
30    format(I4.4)

      call get_chip_id(imagefile,id)
      write(cc,20) id

      n=len(trim(imagefile))
      do i=n,1,-1
        if (imagefile(i:i).eq.'/') exit
      enddo
      b=imagefile(i+23:i+23)
      time=imagefile(i+5:i+8)

      read(time,*) ym
      if (b.eq.'z') then
        if (ym.lt.1308) then
          tttt='1303'
        elseif (ym.lt.1403) then
          tttt='1308'
        elseif (ym.lt.1408) then
          tttt='1403'
        elseif (ym.lt.1411) then
          tttt='1408'
        elseif (ym.lt.1412) then
          tttt='1411'
        elseif (ym.lt.1501) then
          tttt='1412'
        elseif (ym.lt.1502) then
          tttt='1501'
        elseif (ym.lt.1508) then
          tttt='1502'
        elseif (ym.lt.1511) then
          tttt='1508'
        elseif (ym.lt.1512) then
          tttt='1511'
        elseif (ym.lt.1601) then
          tttt='1512'
        elseif (ym.lt.1606) then
          tttt='1601'
        elseif (ym.lt.1608) then
          tttt='1606'
        elseif (ym.lt.1609) then
          tttt='1608'
        elseif (ym.lt.1610) then
          tttt='1609'
        elseif (ym.lt.1611) then
          tttt='1610'
        elseif (ym.lt.1612) then
          tttt='1611'
        elseif (ym.lt.1701) then
          tttt='1612'
        elseif (ym.lt.1703) then
          tttt='1701'
        elseif (ym.lt.1704) then
          tttt='1703'
        elseif (ym.lt.1707) then
          tttt='1704'
        elseif (ym.lt.1709) then
          tttt='1707'
        elseif (ym.lt.1711) then
          tttt='1709'
        elseif (ym.lt.1801) then
          tttt='1711'
        elseif (ym.lt.1805) then
          tttt='1801'
        elseif (ym.lt.1809) then
          tttt='1805'
        elseif (ym.lt.1810) then
          tttt='1809'
        elseif (ym.lt.1811) then
          tttt='1810'
        elseif (ym.lt.1901) then
          tttt='1811'
        else
          tttt='1901'
        endif

      elseif (b.eq.'r') then
        if (ym.lt.1403) then
          tttt='1304'
        elseif (ym.lt.1409) then
          tttt='1403'
        elseif (ym.lt.1410) then
          tttt='1409'
        elseif (ym.lt.1411) then
          tttt='1410'
        elseif (ym.lt.1412) then
          tttt='1411'
        elseif (ym.lt.1501) then
          tttt='1412'
        elseif (ym.lt.1503) then
          tttt='1501'
        elseif (ym.lt.1506) then
          tttt='1503'
        elseif (ym.lt.1511) then
          tttt='1506'
        elseif (ym.lt.1512) then
          tttt='1511'
        elseif (ym.lt.1601) then
          tttt='1512'
        elseif (ym.lt.1602) then
          tttt='1601'
        elseif (ym.lt.1606) then
          tttt='1602'
        elseif (ym.lt.1609) then
          tttt='1606'
        elseif (ym.lt.1610) then
          tttt='1609'
        elseif (ym.lt.1611) then
          tttt='1610'
        elseif (ym.lt.1612) then
          tttt='1611'
        elseif (ym.lt.1701) then
          tttt='1612'
        elseif (ym.lt.1703) then
          tttt='1701'
        elseif (ym.lt.1704) then
          tttt='1703'
        elseif (ym.lt.1707) then
          tttt='1704'
        elseif (ym.lt.1708) then
          tttt='1707'
        elseif (ym.lt.1709) then
          tttt='1708'
        elseif (ym.lt.1710) then
          tttt='1709'
        elseif (ym.lt.1801) then
          tttt='1710'
        elseif (ym.lt.1803) then
          tttt='1801'
        elseif (ym.lt.1805) then
          tttt='1803'
        elseif (ym.lt.1808) then
          tttt='1805'
        elseif (ym.lt.1809) then
          tttt='1808'
        elseif (ym.lt.1810) then
          tttt='1809'
        elseif (ym.lt.1811) then
          tttt='1810'
        elseif (ym.lt.1901) then
          tttt='1811'
        else
          tttt='1901'
        endif

      elseif (b.eq.'g') then
        if (ym.lt.1310) then
          tttt='1308'
        elseif (ym.lt.1401) then
          tttt='1310'
        elseif (ym.lt.1409) then
          tttt='1401'
        elseif (ym.lt.1410) then
          tttt='1409'
        elseif (ym.lt.1412) then
          tttt='1410'
        elseif (ym.lt.1501) then
          tttt='1412'
        elseif (ym.lt.1502) then
          tttt='1501'
        elseif (ym.lt.1506) then
          tttt='1502'
        elseif (ym.lt.1511) then
          tttt='1506'
        elseif (ym.lt.1601) then
          tttt='1511'
        elseif (ym.lt.1602) then
          tttt='1601'
        elseif (ym.lt.1606) then
          tttt='1602'
        elseif (ym.lt.1609) then
          tttt='1606'
        elseif (ym.lt.1610) then
          tttt='1609'
        elseif (ym.lt.1611) then
          tttt='1610'
        elseif (ym.lt.1612) then
          tttt='1611'
        elseif (ym.lt.1701) then
          tttt='1612'
        elseif (ym.lt.1703) then
          tttt='1701'
        elseif (ym.lt.1704) then
          tttt='1703'
        elseif (ym.lt.1707) then
          tttt='1704'
        elseif (ym.lt.1708) then
          tttt='1707'
        elseif (ym.lt.1709) then
          tttt='1708'
        elseif (ym.lt.1710) then
          tttt='1709'
        elseif (ym.lt.1801) then
          tttt='1710'
        elseif (ym.lt.1803) then
          tttt='1801'
        elseif (ym.lt.1805) then
          tttt='1803'
        elseif (ym.lt.1808) then
          tttt='1805'
        elseif (ym.lt.1809) then
          tttt='1808'
        elseif (ym.lt.1810) then
          tttt='1809'
        elseif (ym.lt.1811) then
          tttt='1810'
        elseif (ym.lt.1901) then
          tttt='1811'
        else
          tttt='1901'
        endif

      endif

      flatfile=trim(FLAT_PATH)//'/expo_'//b//'/'
     .//b//tttt//'/flat_'//cc//'_smooth.fits'


      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine chip_process_source(IMAGE_FILE,ichip,DIR_OUTPUT
     .,FLAT_FILE)
      implicit none
      include 'para.inc'

      integer ichip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT,FLAT_FILE
      character*(strl) catfile,PREFIX,PREFIX_head,filename
      character*(strl) sortfile(27)
      integer sortnum
      integer proc_error

      integer nx,ny
      real array(npx,npy),sigmap(npx,npy),normap(npx,npy)
      real flat(npx,npy)

      integer weight(npx,npy)

      integer i,j,u,v,ii
      double precision cRPIX(2),cD(2,2),PU(2,npd),cRVAL(2)

      integer ngal,nstar,nxc

      real sigabc(2,3)


      proc_error=0
      nstar=0
      ngal=0

      if (include_FLAT.eq.1)                    
     . call readimage(FLAT_FILE,nx,ny,npx,npy,flat)
      call readimage(IMAGE_FILE(ichip),nx,ny,npx,npy,array)
      call get_PREFIX(IMAGE_FILE(ichip),PREFIX)

      PREFIX=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)
      filename=trim(PREFIX)//'_norm.fits'
      call readimage(filename,nx,ny,npx,npy,normap)
      if (normap(1,1).ge.0. .or. normap(1,1).lt.(-99990.)
     .    .or. normap(1,1).ne.normap(1,1)) then
        proc_error=1
        write(*,*) 'Error / proc_source rejected failed norm chip',
     .    trim(IMAGE_FILE(ichip))
        return
      endif
      do i=1,ccD_split
        do j=1,3
          sigabc(i,j)=normap(1+i,j)
        enddo
      enddo
      do i=1,nx
        do j=1,ny
          weight(i,j)=1
          if (normap(i,j).lt.(-900.)) weight(i,j)=0
          if (include_FLAT.eq.1) then
            ! if (abs(flat(i,j)).le.flat_thresh) then
              array(i,j)=array(i,j)*flat(i,j)
            ! else
              ! weight(i,j)=0
            ! endif
          endif
        enddo
      enddo

      nxc=nx/2
      if (ccD_split.eq.2) then
        do i=1,nxc
          ii=i+nxc
          do j=1,ny
            sigmap(i,j)=sigabc(1,1)+sigabc(1,2)*i+sigabc(1,3)*j
            sigmap(i,j)=sqrt(0.5*sigmap(i,j))
            sigmap(ii,j)=sigabc(2,1)+sigabc(2,2)*ii+sigabc(2,3)*j
            sigmap(ii,j)=sqrt(0.5*sigmap(ii,j))
          enddo
        enddo
      else
        do i=1,nx
          do j=1,ny
            sigmap(i,j)=sigabc(1,1)+sigabc(1,2)*i+sigabc(1,3)*j
            sigmap(i,j)=sqrt(0.5*sigmap(i,j))
          enddo
        enddo
      endif

c------------------------------------------------------
      call get_expo_catalog(PREFIX,nx,ny,sigmap,weight,normap   
     .,proc_error)

      if (ext_cat.eq.0) then

        call gen_source_catalog(PREFIX,nx,ny,array,weight
     .,ngal,proc_error)

        if (ext_PSF.eq.0) call gen_star_candidate(PREFIX,nstar
     .,proc_error)

      else

        call get_PREFIX_expo(IMAGE_FILE(1),PREFIX_head)
        filename=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX_head)//'.head'
        call read_astrometry_para(filename,ichip    
     .,cRPIX,cD,cRVAL,PU,npd,proc_error)

        catfile=SOURcE_cAT

        if (proc_error.eq.0)
     . call generate_gal_cat_file_name(cRVAL,catfile,sortfile,sortnum)

        if (deblending.eq.1) call de_blending(sortfile,sortnum,nx,ny  
     .,weight,cRPIX,cD,cRVAL,PU,proc_error)


        call gen_source_ext_catalog(sortfile,sortnum               
     .,PREFIX,nx,ny,array,weight,sigmap,cRPIX,cD,cRVAL,PU
     .,ngal,proc_error)

        if (ext_PSF.eq.0) call gen_star_candidate_direct(PREFIX    
     .,nx,ny,array,weight,nstar,proc_error)

      endif
      
      if (proc_error.ne.0) write(*,*) 'Error / proc_source'
     ., trim(IMAGE_FILE(ichip)),proc_error,nstar,ngal

      return
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine de_blending(sortfile,sortnum,nx,ny,weight
     .,cRPIX,cD,cRVAL,PU,proc_error)
      IMPLICIT NONE
      include 'para.inc'

      character*(strl) sortfile(27)
      integer n,sortnum
      double precision cRPIX(2),cD(2,2),cRVAL(2)
      double precision PU(2,npd),z1
      double precision xx,yy,ra,dec,z,dra,ra_c,dec_bound(2),diffra
      real astrometry_shift_ratio
      parameter (astrometry_shift_ratio=0.2)

      integer nx,ny,proc_error,ierr
      integer weight(npx,npy)
      integer i,j,ix,iy,iz,old_w,new_w
      double precision flags_ft,flags_fg,flags_gold,ext_mash
      double precision mag_g,mag_r,mag_i,mag_z,mag_y
      double precision magerr_g,magerr_r,magerr_i,magerr_z,magerr_y
      double precision zp,zperr

      if (proc_error.eq.1) return

      call get_ra_dec_range_fine(nx,ny,ra_c,dec_bound,dra  
     .,cRPIX,cD,cRVAL,PU,npd,astrometry_shift_ratio)

      do n=1,sortnum
        open(unit=10,file=trim(sortfile(n)),status='old',action='read'
     .,iostat=ierr)
        if (ierr.ne.0) then
    !       write(*,*) trim(sortfile(n))
    !  .,'  catalog file error in deblending!!'
          cycle
        endif
        rewind 10
        read(10,*)

        do while (ierr.ge.0)
          read(10,*,iostat=ierr) flags_ft,flags_fg,flags_gold,ext_mash,
     .ra,dec,mag_g,magerr_g,mag_r,magerr_r,
     .mag_i,magerr_i,mag_z,magerr_z,mag_y,magerr_y,zp,zperr
          if (ierr.lt.0) cycle
          if (zp.lt.0 .or. zp.gt.5.) cycle

          z=zp
          if (abs(diffra(ra,ra_c)).gt.dra*0.5) cycle
          if (dec.lt.dec_bound(1).or.dec.gt.dec_bound(2)) cycle
          call coordinate_transfer_PU(ra,dec,xx,yy,-1
     .,cRPIX,cD,cRVAL,PU,npd)
          ix=int(xx+0.5)
          iy=int(yy+0.5)
          if (ix.lt.1 .or. ix.gt.nx .or. iy.lt.1 .or. iy.gt.ny) cycle
          iz=int(z*1000.+10)

          if (weight(ix,iy).lt.2) then
            cycle
          elseif (weight(ix,iy).gt.2) then
            z1=0.001*(weight(ix,iy)-10.)
            if (abs(z-z1).gt.dz_thresh) then
              new_w=0
              old_w=weight(ix,iy)
              call fill_patch(nx,ny,weight,ix,iy,old_w,new_w)
            endif
          else
            old_w=2
            new_w=iz
            call fill_patch(nx,ny,weight,ix,iy,old_w,new_w)
          endif
        enddo
        close(10)
      enddo

      do i=1,nx
        do j=1,ny
          if (weight(i,j).gt.2) weight(i,j)=2
        enddo
      enddo

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fill_patch(nx,ny,map,ix,iy,old_v,new_v)
      implicit none
      include 'para.inc'

      integer nx,ny,ix,iy,old_v,new_v
      integer map(npx,npy)
      integer changed,xmin,xmax,ymin,ymax,i,j,u,v

      if (old_v.eq.new_v) return

      map(ix,iy)=new_v
      changed=1
      xmin=ix
      xmax=xmin
      ymin=iy
      ymax=ymin

      do while (changed.eq.1)
        changed=0
        do i=xmin,xmax
          do j=ymin,ymax
            if (map(i,j).eq.new_v) then
              map(i,j)=-new_v
              do u=max(i-1,1),min(i+1,nx)
                do v=max(j-1,1),min(j+1,ny)
                  if (map(u,v).eq.old_v) then
                    map(u,v)=new_v
                    changed=1
                    xmin=min(xmin,u)
                    ymin=min(ymin,v)
                    xmax=max(xmax,u)
                    ymax=max(ymax,v)
                  endif
                enddo
              enddo
            endif
          enddo
        enddo
      enddo

      do i=xmin,xmax
        do j=ymin,ymax
          if (map(i,j).eq.-new_v) map(i,j)=new_v
        enddo
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_expo_catalog(PREFIX,nx,ny,sigmap
     .,weight,normap,ierror)
      implicit none
      include 'para.inc'

      character*(strl) catname,PREFIX
      integer nx,ny,ierror
      real sigmap(npx,npy),normap(npx,npy)
      integer weight(npx,npy),mark(npx,npy)

      integer nb,nbb,toobig,buffer(area_max,2)
      integer i,j,ix,iy,jx,jy,u,v,k1,k2,k
      real xp,yp,r2,rmax,sig,temp,peak,thresh,xc,yc,SNR
      real total_flux,half_light_flux
      integer half_light_area,total_area


      if (ierror.eq.1) return

      do i=1,nx
        do j=1,ny
          if (normap(i,j).ge.source_thresh.and.weight(i,j).ge.1) then
            mark(i,j)=1
          else
            mark(i,j)=0
          endif
        enddo
      enddo

      catname=trim(PREFIX)//'.cat'

      open(unit=10,file=catname,status='replace')
      rewind 10

      write(10,*) ' xp ',' yp '
     .,' total_area ',' half_light_area ',' sig '
     .,' sig_normed_total_flux ',' sig_normed_half_light_flux '
     .,' sig_normed_peak ',' rmax '

      do i=1,nx
        do j=1,ny
          if (mark(i,j).eq.1) then
            nbb=0
            nb=1
            buffer(nb,1)=i
            buffer(nb,2)=j
            mark(i,j)=-1
            toobig=0
            do while (nb.gt.nbb)
              k1=nbb+1
              k2=nb
              nbb=nb
              do k=k1,k2
                ix=buffer(k,1)
                iy=buffer(k,2)
                do u=max(ix-1,1),min(ix+1,nx)
                  do v=max(iy-1,1),min(iy+1,ny)
                    if (mark(u,v).eq.1) then
                      nb=nb+1
                      buffer(nb,1)=u
                      buffer(nb,2)=v
                      mark(u,v)=-1
                      if (nb.eq.area_max) then
                        toobig=1
                        goto 20
                      endif
                    elseif (mark(u,v).gt.1) then
                      toobig=1
                      goto 20
                    endif
                  enddo
                enddo
              enddo
            enddo
20          if (toobig.eq.1) then
              do k=1,nb
                mark(buffer(k,1),buffer(k,2))=area_max
                weight(buffer(k,1),buffer(k,2))=2
              enddo
            else
              if (nb.ge.area_thresh) then
                xc=0.
                yc=0.
                total_flux=0.
                sig=0.
                peak=-100000.
                do k=1,nb
                  ix=buffer(k,1)
                  iy=buffer(k,2)
                  mark(ix,iy)=nb
                  weight(ix,iy)=2

                  xc=xc+ix*normap(ix,iy)
                  yc=yc+iy*normap(ix,iy)
                  total_flux=total_flux+normap(ix,iy)
                  sig=sig+sigmap(ix,iy)

                  if (normap(ix,iy).gt.peak) then
                    peak=normap(ix,iy)
                    xp=ix
                    yp=iy
                  endif
                enddo
                xc=xc/total_flux
                yc=yc/total_flux
                total_area=nb
                sig=sig/total_area

                thresh=peak*0.5
                rmax=0
                half_light_area=0
                half_light_flux=0.
                do k=1,nb
                  ix=buffer(k,1)
                  iy=buffer(k,2)
                  r2=(ix-xc)**2+(iy-yc)**2
                  rmax=max(rmax,r2)
                  if (normap(ix,iy).ge.thresh) then
                    half_light_area=half_light_area+1
                    half_light_flux=half_light_flux+normap(ix,iy)
                  endif
                enddo
                rmax=sqrt(rmax)
                if (peak.ge.core_thresh) then
                  write(10,*) xp,yp,total_area,half_light_area
     .,sig,total_flux,half_light_flux,peak,rmax
                endif
              else
                do k=1,nb
                  ix=buffer(k,1)
                  iy=buffer(k,2)
                  mark(ix,iy)=nb
                enddo
              endif
            endif
          endif
        enddo
      enddo


      close(10)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gen_source_catalog(PREFIX,nx,ny,array
     .,weight,ngal,proc_error)
      implicit none
      include 'para.inc'

      character*(strl) catname,filename,PREFIX
      integer ierror,flag,proc_error,ngal

      integer nx,ny
      real array(npx,npy)
      integer weight(npx,npy)

      integer i,j,x1,x2,y1,y2,u,v,nn,nn1,nn2,ig,ix,iy,cc,tempi

      real source_para(ngal_max,npara)
      real source_collect(ngal_max,ns,ns)
      real noise_collect(ngal_max,ns,ns)
      real source(ns,ns),noise(ns,ns)

      real SNR,temp,flux_alt,flux2

      real xp,yp,sig,total_flux,half_light_flux,peak,rf
      integer total_area,half_light_area
      common /stamp_pass/ xp,yp,total_area,half_light_area,sig
     .,total_flux,half_light_flux,peak,rf

      integer imax,jmax
      common /defect_pass/ imax,jmax

      ngal=0

      if (proc_error.eq.1) goto 40

      catname=trim(PREFIX)//'.cat'

      open(unit=10,file=catname,status='old',iostat=ierror)
      rewind 10

      if (ierror.ne.0) then
        write(*,*) catname
        stop 'Error / gen_source_catalog catalog file error!!'
      endif
      read(10,*)

      ig=0
      do while (ierror.ge.0)

        read(10,*,iostat=ierror) xp,yp,total_area
     .,half_light_area,sig,total_flux,half_light_flux,peak,rf

        if (ierror.lt.0) cycle

        ig=ig+1
        call find_noise(flag,noise,nx,ny,array,weight)
        if (flag.lt.0) cycle

        call check_source(flag,source,nx,ny,array,weight)
        if (flag.lt.0) cycle

        ngal=ngal+1

        do u=1,ns
          do v=1,ns
            source_collect(ngal,u,v)=source(u,v)
            noise_collect(ngal,u,v)=noise(u,v)
          enddo
        enddo

        source_para(ngal,1)=ig
        source_para(ngal,2)=xp
        source_para(ngal,3)=yp
        source_para(ngal,4)=sig
        source_para(ngal,5)=peak
        source_para(ngal,6)=imax
        source_para(ngal,7)=jmax
        source_para(ngal,8)=half_light_flux
        source_para(ngal,9)=half_light_area
        source_para(ngal,10)=flag

        if (ngal.ge.ngal_max) exit

      enddo
      close(10)


40    nn1=ns*len_g
      nn2=ns*(int(ngal/len_g)+1)
      filename=trim(PREFIX)//'_source.fits'
      call write_stamps(ngal_max,1,ngal,ns,ns
     .,source_collect,nn1,nn2,filename)

      filename=trim(PREFIX)//'_noise.fits'
      call write_stamps(ngal_max,1,ngal,ns,ns
     .,noise_collect,nn1,nn2,filename)

      filename=trim(PREFIX)//'_source_info.dat'
      open(unit=10,file=filename,status='replace')
      rewind 10
      write(10,*) 'ig xp yp sigma peak imax jmax '
     .,'half_light_flux half_light_area flag'
      do i=1,ngal
        write(10,*) (source_para(i,j),j=1,iflag)
      enddo
      close(10)


      return
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gen_source_ext_catalog(sortfile,sortnum,PREFIX
     .,nx,ny,array,weight,sigmap,cRPIX,cD,cRVAL,PU,ngal,proc_error)
      implicit none
      include 'para.inc'

c      The purpose of this subroutine to extract sources from an existing external catalog.

      character*(strl) sortfile(27),filename,PREFIX
      integer n,sortnum
      integer ierror,flag,ngal,proc_error

      character*500 cat_content,cat_list

      integer nx,ny
      real array(npx,npy),sigmap(npx,npy)
      integer weight(npx,npy)

      integer i,j,u,v,nn1,nn2,ig,igal

      integer sid(ngal_max)

      real source_para(ngal_max,npara)
      real source_collect(ngal_max,ns,ns)
      real noise_collect(ngal_max,ns,ns)
      real source(ns,ns),noise(ns,ns)

      real xp,yp,sig,total_flux,half_light_flux,peak,rf
      integer total_area,half_light_area
      common /stamp_pass/ xp,yp,total_area,half_light_area,sig
     .,total_flux,half_light_flux,peak,rf

      real SNR,temp

      double precision cRPIX(2),cD(2,2),cRVAL(2)
      double precision PU(2,npd)
      double precision xx,yy,ra,dec,dra,ra_c,dec_bound(2),diffra
      real astrometry_shift_ratio
      parameter (astrometry_shift_ratio=0.2)
      real dummy1,dummy2,dummy3,dummy4

      integer imax,jmax
      common /defect_pass/ imax,jmax

      integer orighead

      ngal=0
      if (proc_error.eq.1) goto 40

      call get_ra_dec_range_fine(nx,ny,ra_c,dec_bound,dra
     .,cRPIX,cD,cRVAL,PU,npd,astrometry_shift_ratio)

      ig=0
      do n=1,sortnum
        open(unit=10,file=trim(sortfile(n)),status='old',action='read'
     .,iostat=ierror)
        if (ierror.ne.0) then
    !       write(*,*) trim(sortfile(n))
    !  .,'  catalog file error in source_ext!!'
          cycle
        endif
        rewind 10
        read(10,*)

        do while (ierror.ge.0)
          read(10,*,iostat=ierror) dummy1,dummy2,dummy3,dummy4,ra,dec
          if (ierror.lt.0) cycle
          ig=ig+1

          if (abs(diffra(ra,ra_c)).gt.dra*0.5) cycle
          if (dec.lt.dec_bound(1).or.dec.gt.dec_bound(2)) cycle

          call coordinate_transfer_PU(ra,dec,xx,yy,-1
     .,cRPIX,cD,cRVAL,PU,npd)

          xp=xx
          yp=yy

          if (xp-nl_2.lt.chip_margin.or.xp+nl_2.gt.nx-chip_margin
     . .or. yp-nl_2.lt.chip_margin.or.yp+nl_2.gt.ny-chip_margin) cycle

          sig=sigmap(int(xp+0.5),int(yp+0.5))

          call find_noise(flag,noise,nx,ny,array,weight)
          if (flag.lt.0) cycle
          call check_source(flag,source,nx,ny,array,weight)
          if (flag.lt.0) cycle

          ngal=ngal+1
          do u=1,ns
            do v=1,ns
              source_collect(ngal,u,v)=source(u,v)
              noise_collect(ngal,u,v)=noise(u,v)
            enddo
          enddo
          source_para(ngal,1)=ig
          source_para(ngal,2)=xp
          source_para(ngal,3)=yp
          source_para(ngal,4)=sig
          source_para(ngal,5)=peak
          source_para(ngal,6)=imax
          source_para(ngal,7)=jmax
          source_para(ngal,8)=half_light_flux
          source_para(ngal,9)=half_light_area
          source_para(ngal,10)=flag

          sid(ngal)=ig

          if (ngal.ge.ngal_max) then
            close(10)
            goto 40
          endif
        enddo
        close(10)
      enddo

40    if (ngal.gt.0) then
        nn1=ns*len_g
        nn2=ns*(int(ngal/len_g)+1)
        filename=trim(PREFIX)//'_source.fits'
        call write_stamps(ngal_max,1,ngal,ns,ns
     .,source_collect,nn1,nn2,filename)

        filename=trim(PREFIX)//'_noise.fits'
        call write_stamps(ngal_max,1,ngal,ns,ns
     .,noise_collect,nn1,nn2,filename)
      endif

      filename=trim(PREFIX)//'_source_info.dat'
      open(unit=10,file=filename,status='replace')
      rewind 10
      write(10,*) 'ig xp yp sigma peak imax '
     .,'jmax half_light_flux half_light_area flag'
      do i=1,ngal
        write(10,*) (source_para(i,j),j=1,iflag)
      enddo
      close(10)

      orighead = 0
      filename=trim(PREFIX)//'_orig.cat'
      open(unit=15,file=filename,status='replace')
      rewind 15
      if (proc_error.eq.1 .or. ngal.eq.0) then
        write(15,*) 'No sources!!'
      else
        i=0
        igal=1
        ig=sid(igal)
        do n=1,sortnum
        open(unit=10,file=trim(sortfile(n)),status='old',action='read'
     .,iostat=ierror)
          if (ierror.ne.0) then
            cycle
          endif
          rewind 10
          read(10,'(A)',iostat=ierror) cat_list
          if (orighead.eq.0) then
            write(15,'(A)') trim(cat_list)
            orighead = 1
          endif
          do while (ierror.ge.0 .and. igal.le.ngal)
            read(10,'(A)',iostat=ierror) cat_content
            if (ierror.lt.0) cycle
            i=i+1
            if (i.eq.ig) then
              write(15,'(A)') trim(cat_content)
              igal=igal+1
              if (igal.le.ngal) ig=sid(igal)
            endif
          enddo
          close(10)
        enddo
      endif

      close(15)

      return
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine find_noise(flag,stamps,nx,ny,array,weight)
      implicit none
      include 'para.inc'

      integer x1,x2,y1,y2,flag
      real stamps(ns,ns),stampl(nl,nl)
      integer weights(ns,ns),weightl(nl,nl)

      integer nx,ny,nsteph
      real array(npx,npy)
      integer weight(npx,npy)

      integer i,j,ix1,iy1,ix2,iy2,nd,u,v,jx,jy,xmid
      real npeak,pmax,temp,roughness

      real xp,yp,sig,total_flux,half_light_flux,peak,rf
      real rough
      integer total_area,half_light_area
      common /stamp_pass/ xp,yp,total_area,half_light_area
     .,sig,total_flux,half_light_flux,peak,rf

      flag=0

      x1=int(xp+0.5)-nl_2-1
      y1=int(yp+0.5)-nl_2-1

      temp=100000.
      pmax=temp
      xmid=nx/2

      nsteph=5
      do i=1,nsteph
        do j=1,nsteph
          if ((i.eq.1).or.(i.eq.nsteph).or.(j.eq.1).or.(j.eq.nsteph))
     . then
            ix1=x1+nl*(i-(nsteph+1)/2)
            iy1=y1+nl*(j-(nsteph+1)/2)
            ix2=ix1+nl-1
            iy2=iy1+nl-1
            if (ix1.lt.chip_edge_margin.or.ix2.gt.nx-chip_edge_margin
     . .or. iy1.lt.chip_edge_margin
     . .or. iy2.gt.ny-chip_edge_margin) cycle

            npeak=-temp
            nd=0
            do u=ix1,ix2
              do v=iy1,iy2
                if (array(u,v).gt.npeak) npeak=array(u,v)
                if (weight(u,v).eq.0) nd=nd+1
              enddo
            enddo
            if (nd.le.nl.and.npeak.lt.pmax) then
              pmax=npeak
              jx=ix1
              jy=iy1
            endif
          endif
        enddo
      enddo

      if (pmax.eq.temp) then
        flag=-1
        return
      endif

c//////////To decorate the defects \\\\\\\\\\\\\\

      do i=jx,jx+nl-1
        u=i-jx+1
        do j=jy,jy+nl-1
          v=j-jy+1
          stampl(u,v)=array(i,j)
          weightl(u,v)=weight(i,j)
        enddo
      enddo

      call flatten_stamp_2D(ns,nl,stampl,weightl,flag)

      if (flag.lt.0) return

      call mark_noise(nl,stampl,weightl,sig,source_thresh
     .,core_thresh)

      x1=nl_2-ns_2
      y1=nl_2-ns_2

      do i=1,ns
        do j=1,ns
          stamps(i,j)=stampl(x1+i,y1+j)
          weights(i,j)=weightl(x1+i,y1+j)
        enddo
      enddo


c/////////To decorate the image\\\\\\\\\\\\\\\\

      call decorate_stamp(ns,sig,weights,stamps)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine check_source(flag,stamps,nx,ny,array,weight)
      implicit none
      include 'para.inc'

      integer x,y,x1,x2,y1,y2,flag,xmid
      real stamps(ns,ns),stampl(nl,nl)
      integer weights(ns,ns),weightl(nl,nl)

      integer nx,ny
      real array(npx,npy),roughness,vec(ns)
      integer weight(npx,npy)

      integer i,j,k,changed,u,v,ic,ixc,iyc,boundx(2),boundy(2)
      integer offx,offy,xcenter,ycenter

      real xp,yp,sig,total_flux,half_light_flux,peak,rf
      real rough
      integer total_area,half_light_area
      common /stamp_pass/ xp,yp,total_area,half_light_area
     .,sig,total_flux,half_light_flux,peak,rf

      integer imax,jmax
      common /defect_pass/ imax,jmax

      real radius,temp

      flag=0

      x=int(xp+0.5)
      y=int(yp+0.5)


      if (x-nl_2.lt.chip_edge_margin.or.x+nl_2.gt.nx-chip_edge_margin
     . .or. y-nl_2.lt.chip_edge_margin
     . .or. y+nl_2.gt.ny-chip_edge_margin) then
        flag=-1
        return
      endif


      x1=x-nl_2-1
      y1=y-nl_2-1

      do i=1,nl
        do j=1,nl
          stampl(i,j)=array(i+x1,j+y1)
          weightl(i,j)=weight(i+x1,j+y1)
        enddo
      enddo

      call flatten_stamp_2D(ns,nl,stampl,weightl,flag)

      if (flag.lt.0) return

      call mark_source(nl,stampl,weightl,sig,source_thresh
     .,core_thresh,boundx,boundy,total_flux,total_area,peak
     .,half_light_flux,half_light_area,flag,radius,xcenter,ycenter)

      if (flag.lt.0) return

      if (peak.gt.saturation_thresh/sig) then
        flag=-1
        return
      endif

      temp=ns_2-flag_thresh
      if (radius.ge.temp) then
        flag=-1
        return
      endif

      x1=xcenter-ns_2
      y1=ycenter-ns_2

      if (x1.lt.1 .or. x1+ns-1.gt.nl
     . .or. y1.lt.1 .or. y1+ns-1.gt.nl) then
        flag=-1
        return
      endif

      imax=0
      do i=1,ns
        u=0
        do j=1,ns
          stamps(i,j)=stampl(x1+i-1,y1+j-1)
          weights(i,j)=weightl(x1+i-1,y1+j-1)
          if (weights(i,j).eq.0) u=u+1
        enddo
        if (u.gt.imax) imax=u
      enddo

      jmax=0
      do j=1,ns
        u=0
        do i=1,ns
          if (weights(i,j).eq.0) u=u+1
        enddo
        if (u.gt.jmax) jmax=u
      enddo


c/////////To decorate the image\\\\\\\\\\\\\\\\

      call decorate_stamp(ns,sig,weights,stamps)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gen_star_candidate(PREFIX,nstar,proc_error)
      implicit none
      include 'para.inc'

      character*(strl) filename,PREFIX
      integer ierror,nstar,proc_error,ngal

      real source_para(ngal_max,npara)
      real source_collect(ngal_max,ns,ns)
      real noise_collect(ngal_max,ns,ns)
      real star_para(nstar_max,npara),aa(npara)
      real SNR,temp
      integer i,j,u,v,nn1,nn2


      nstar=0

      if (proc_error.eq.1) goto 40

      ngal=0

      filename=trim(PREFIX)//'_source_info.dat'
      open(unit=10,file=filename,status='old',iostat=ierror)
      rewind 10
      if (ierror.ne.0) then
        write(*,*) filename
        stop 'Error / gen_star_candidate catalog file error!!'
      endif
      read(10,*)
c      read(10,*) 'ig xp yp sigma peak imax jmax half_light_flux'
c     .,' half_light_area flag'
      do while (ierror.ge.0)
        read(10,*,iostat=ierror) (aa(i),i=1,iflag)
        if (ierror.lt.0) cycle
        ngal=ngal+1
        do i=1,iflag
          source_para(ngal,i)=aa(i)
        enddo
        SNR=aa(ih_flux)/sqrt(aa(ih_area))
        source_para(ngal,11)=SNR

      enddo
      close(10)

      if (ngal.gt.0) then
        nn1=ns*len_g
        nn2=ns*(int(ngal/len_g)+1)

        filename=trim(PREFIX)//'_noise.fits'
        call read_stamps(ngal_max,1,ngal,ns,ns
     .,noise_collect,nn1,nn2,filename)

        filename=trim(PREFIX)//'_source.fits'
        call read_stamps(ngal_max,1,ngal,ns,ns
     .,source_collect,nn1,nn2,filename)
      endif

      do i=1,ngal

        if (source_para(i,11).lt.SNR_PSF) cycle

        nstar=nstar+1
        do u=1,ns
          do v=1,ns
            source_collect(nstar,u,v)=source_collect(i,u,v)
            noise_collect(nstar,u,v)=noise_collect(i,u,v)
          enddo
        enddo
        star_para(nstar,1)=source_para(i,1)
        star_para(nstar,2)=source_para(i,2)
        star_para(nstar,3)=source_para(i,3)
        star_para(nstar,4)=source_para(i,11)

        if (nstar.ge.nstar_max) exit

      enddo

40    filename=trim(PREFIX)//'_star_can_info.dat'
      open(unit=20,file=filename,status='replace')
      rewind 20
      write(20,*) 'ig xp yp SNR'
      do i=1,nstar
        write(20,*) (star_para(i,j),j=1,4)
      enddo
      close(20)

      if (nstar.gt.0) then
        filename=trim(PREFIX)//'_star_can.fits'
        nn1=ns*len_s
        nn2=ns*(int(nstar/len_s)+1)
        call write_stamps(ngal_max,1,nstar,ns,ns
     .,source_collect,nn1,nn2,filename)

        filename=trim(PREFIX)//'_star_can_noise.fits'
        nn1=ns*len_s
        nn2=ns*(int(nstar/len_s)+1)
        call write_stamps(ngal_max,1,nstar,ns,ns
     .,noise_collect,nn1,nn2,filename)
      endif


      return
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gen_star_candidate_direct(PREFIX,nx,ny,array,weight
     .,nstar,proc_error)
      implicit none
      include 'para.inc'

      character*(strl) catname,filename,PREFIX
      integer ierror,flag,proc_error,nstar

      integer nx,ny
      real array(npx,npy)
      integer weight(npx,npy)

      integer i,j,u,v,nn1,nn2,ig

      real star_para(nstar_max,npara),aa(npara)
      real source_coll(ngal_max,ns,ns)
      real noise_coll(ngal_max,ns,ns)

      real source(ns,ns),noise(ns,ns)

      real SNR,temp

      real xp,yp,sig,total_flux,half_light_flux,peak,rf
      integer total_area,half_light_area
      common /stamp_pass/ xp,yp,total_area,half_light_area,sig
     .,total_flux,half_light_flux,peak,rf

      nstar=0

      if (proc_error.eq.1) goto 40

      catname=trim(PREFIX)//'.cat'

      open(unit=10,file=catname,status='old',iostat=ierror)
      rewind 10

      if (ierror.ne.0) then
        write(*,*) catname
        stop 'Error / gen_star_candidate_direct catalog file error!!'
      endif
      read(10,*)

      do while (ierror.ge.0)

        read(10,*,iostat=ierror) xp,yp,total_area
     .,half_light_area,sig,total_flux,half_light_flux,peak,rf

        if (ierror.lt.0) cycle

        temp=half_light_area
        SNR=half_light_flux/sqrt(temp)
        if (SNR.lt.SNR_PSF*0.5) cycle

        call find_noise(flag,noise,nx,ny,array,weight)
        if (flag.lt.0) cycle

        call check_source(flag,source,nx,ny,array,weight)
        if (flag.lt.0) cycle

        temp=half_light_area
        SNR=half_light_flux/sqrt(temp)

        if (SNR.lt.SNR_PSF) cycle

        nstar=nstar+1
        do u=1,ns
          do v=1,ns
            source_coll(nstar,u,v)=source(u,v)
            noise_coll(nstar,u,v)=noise(u,v)
          enddo
        enddo
        star_para(nstar,1)=nstar
        star_para(nstar,2)=xp
        star_para(nstar,3)=yp
        star_para(nstar,4)=SNR

        if (nstar.ge.nstar_max) exit

      enddo
      close(10)

40    filename=trim(PREFIX)//'_star_can_info.dat'
      open(unit=20,file=filename,status='replace')
      rewind 20
      write(20,*) 'ig xp yp SNR'
      do i=1,nstar
        write(20,*) (star_para(i,j),j=1,4)
      enddo
      close(20)

      if (nstar.gt.0) then
        filename=trim(PREFIX)//'_star_can.fits'
        nn1=ns*len_s
        nn2=ns*(int(nstar/len_s)+1)
        call write_stamps(ngal_max,1,nstar,ns,ns
     .,source_coll,nn1,nn2,filename)

        filename=trim(PREFIX)//'_star_can_noise.fits'
        nn1=ns*len_s
        nn2=ns*(int(nstar/len_s)+1)
        call write_stamps(ngal_max,1,nstar,ns,ns
     .,noise_coll,nn1,nn2,filename)
      endif


      return
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
