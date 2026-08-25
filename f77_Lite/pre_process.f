      subroutine pre_process(iexpo)
      implicit none
      include 'para.inc'

      integer iexpo,nchip,ichip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT
      integer i

      call get_image_list(iexpo,IMAGE_FILE,nchip,DIR_OUTPUT)

      do ichip=1,nchip
        call get_chip_id(IMAGE_FILE(ichip),i)

        call chip_pre_process(IMAGE_FILE(ichip),DIR_OUTPUT,i)
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine chip_pre_process(IMAGE_FILE,DIR_OUTPUT,cid)
      implicit none
      include 'para.inc'
      include 'sig_para.inc'

      character*(*) IMAGE_FILE,DIR_OUTPUT
      integer proc_error
      character*(strl) PREFIX,filename,catfile,PREFIX_e

      integer nx,ny
      real array(npx,npy),normap(npx,npy)
      integer weight(npx,npy)

      integer i,j,u,v,nxc
      double precision cRPIX(2),cD(2,2),cRVAL(2)

      real sigabc(2,3)
      integer order
      common /sig_pass/ sigabc,order

      real aa,bb,cc

      character*(strl) MASK_FILE
      integer nxx , nyy
      real flat_weight(npx,npy)

      integer cid
      character*(2) c_chip


      proc_error=0
      do i=1,2
        do j=1,3
          sigabc(i,j)=0.
        enddo
      enddo

      call readimage_para(IMAGE_FILE
     .,nx,ny,npx,npy,array,cRPIX,cD,cRVAL)
      if (array(1,1).lt.(-99990.0)) then
        write(*,*) 'Error / Pre_proc cant find image file!'
        proc_error=1
      endif

      do i=1,nx
        do j=1,ny
          weight(i,j)=1
          if (array(i,j).gt.saturation_thresh) weight(i,j)=0
          normap(i,j)=array(i,j)
        enddo
      enddo

c ==========================================
c Function: Apply the DQ mask before astrometry and defect merging
c Method: Keep set_sig independent of DQ, then reject every nonzero
c         DQ pixel from all later preprocessing stages.
c ==========================================
      if (proc_error.eq.0) then
        write(c_chip,'(I2)') cid
        call get_PREFIX_expo(IMAGE_FILE,PREFIX_e)
        MASK_FILE=trim(DIR_OUTPUT)//'/dqmask/'//trim(PREFIX_e)
     .  //'_'//trim(adjustl(c_chip))//'.fits'
        call readimage(MASK_FILE,nxx,nyy,npx,npy,flat_weight)
        if (flat_weight(1,1).lt.(-99990.0)) then
          write(*,*) 'Error / cant find mask file!'
          proc_error=1
        elseif (nxx.ne.nx .or. nyy.ne.ny) then
          write(*,*) 'Error / wrong size of DQ file!'
          proc_error=1
        endif
        if (proc_error.eq.0) then
          do i=1,nx
            do j=1,ny
              if (abs(flat_weight(i,j)).gt.1e-7) then
                weight(i,j)=0
              endif
            enddo
          enddo
        endif
      endif
c------------------------------------------------------

      nxc=nx/2

      if (ccD_split.eq.2) then

        call set_background(1,nxc,1,ny,npx,npy,normap
     .,blocksize,nct,ncx,proc_error)
        call set_background(nxc+1,nx,1,ny,npx,npy,normap
     .,blocksize,nct,ncx,proc_error)
        call set_sig(1,nxc,1,ny,npx,npy,normap,weight
     .,aa,bb,cc,proc_error,sig_scale)
        if (proc_error.eq.0) then
          sigabc(1,1)=aa
          sigabc(1,2)=bb
          sigabc(1,3)=cc
        endif
        if (proc_error.eq.0) then
          call set_sig(nxc+1,nx,1,ny,npx,npy,normap,weight
     .    ,aa,bb,cc,proc_error,sig_scale)
        endif
        if (proc_error.eq.0) then
          sigabc(2,1)=aa
          sigabc(2,2)=bb
          sigabc(2,3)=cc
        endif
      else
        call set_background(1,nx,1,ny,npx,npy,normap        
     .,blocksize,nct,ncx,proc_error)
        call set_sig(1,nx,1,ny,npx,npy,normap,weight
     .,aa,bb,cc,proc_error,sig_scale)
        if (proc_error.eq.0) then
          sigabc(1,1)=aa
          sigabc(1,2)=bb
          sigabc(1,3)=cc
        endif
      endif


c--------------------------------------------------------------
      call get_PREFIX(IMAGE_FILE,PREFIX)
      filename=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX)//'_astro.dat'
      catfile=ASTROMETRY_cAT
      call generate_gaia_file_name(cRVAL,catfile,proc_error)
      call gen_astrometry_data(catfile,nx,ny,npx,npy
     .,normap,weight,cRPIX,cD,cRVAL,filename,proc_error)

      call locate_defects(nx,ny,npx,npy,array,normap
     .,weight,area_max,area_thresh,proc_error)
      call merge_defects(nx,ny,npx,npy,weight,normap
     .,area_max,source_thresh,area_thresh,proc_error)

c------------------------------------------------------------------

      if (proc_error.eq.0) then
c Serialize the final combined saturation, DQ, and detected-defect mask.
        do i=1,nx
          do j=1,ny
            if (weight(i,j).eq.0) normap(i,j)=-1000.
          enddo
        enddo
      else
        do i=1,nx
          do j=1,ny
            normap(i,j)=-1000.
          enddo
        enddo
      endif

      if (proc_error.eq.0) then
        normap(1,1)=-1.
      else
        normap(1,1)=1.
      endif

      do i=1,ccD_split
        do j=1,3
          normap(1+i,j)=sigabc(i,j)
        enddo
      enddo

      filename=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)
     .//'_norm.fits'
      call writeimage_copyhdu(IMAGE_FILE,filename
     .,nx,ny,npx,npy,normap)

c-----------------------------------------------------------
    !   if (proc_error.eq.0) then
    !     write(*,*) 'Status of processing ',trim(IMAGE_FILE)
    !  .,': OK.'
    !   else
    !     write(*,*) 'Status of processing ',trim(IMAGE_FILE)
    !  .,': ERROR!'
    !   endif
    !   if (proc_error.ne.0) 
    !  .write(*,*) 'Error / pre-processing ', trim(IMAGE_FILE)

c      pause

      return
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine set_background(nx1,nx2,ny1,ny2,npx,npy,image
     .,blocksize,nct,ncx,ierror)
      implicit none

c The purpose of this code is to achieve a fine adjustment of the background.

c input & output:
      integer nx1,nx2,ny1,ny2,npx,npy,ierror,blocksize,nct,ncx
      real image(npx,npy)

c parameters:
      integer npp,nfit
      parameter (npp=1000)
      parameter (nfit=5000)

c local variables:
      real pix(npp),arr(npp,3),c(nct)
      real mean_sam(nfit,3),sam(nfit),tempo(nfit),tmp_fit(nfit,3)
      integer indx(npp),nsam,nsam1,changed
      integer ix,iy,i,j,k,nbx,nby,xmin,xmax,ymin,ymax,nct_min
      real i2,j2,ij,mean,sig,aa,bb,cc,med,x,y
      real bound1,bound2
c sub-region affine normalization (numerical_fix F1):
      real xmid,ymid,xhalf_1,yhalf_1
c uses:
      real ran1,func_val

      if (ierror.eq.1) return

c a rough flattening of the field first:
      call flatten_chip(nx1,nx2,ny1,ny2,npx,npy,image,4,2,ierror)

c      return

      if (ierror.eq.1) return

c numerical_fix F1: the legacy code used x=i*ratio with
c ratio=1/max(dx,dy) and the ABSOLUTE pixel index, so the
c second amplifier landed on
c x in [0.275,0.470] where {1,x,x^2} are nearly collinear
c (cond(A^T A)=2.0e9, eps32*cond=235 -> no significant digit left).
c Mapping each sub-region onto [-1,1]x[-1,1] leaves the fitted surface
c unchanged in exact arithmetic (the tensor-product monomial span is
c invariant under an affine change of x and of y) and drops the
c condition number to 5.3e2.
      xmid=0.5*(nx1+nx2)
      ymid=0.5*(ny1+ny2)
      xhalf_1=2./real(max(nx2-nx1,1))
      yhalf_1=2./real(max(ny2-ny1,1))

c divide the chip into nbx*nby blocks, each with sidelength of "blocksize":
      nbx=max((nx2-nx1)/blocksize,1)
      nby=max((ny2-ny1)/blocksize,1)

c get the median of each block:
      nsam=0
      do i=1,nbx
        xmin=(i-1)*blocksize+nx1
        xmax=min(xmin+blocksize,nx2)
        do j=1,nby
          ymin=(j-1)*blocksize+ny1
          ymax=min(ymin+blocksize,ny2)

          do k=1,npp
            ix=int(ran1()*(xmax-xmin)+xmin)
            iy=int(ran1()*(ymax-ymin)+ymin)
            arr(k,1)=(ix-xmid)*xhalf_1
            arr(k,2)=(iy-ymid)*yhalf_1
            arr(k,3)=image(ix,iy)
            pix(k)=arr(k,3)
          enddo

          call indexx(npp,npp,pix,indx)

          k=indx(npp/2)
          nsam=nsam+1

          mean_sam(nsam,1)=arr(k,1)
          mean_sam(nsam,2)=arr(k,2)
          mean_sam(nsam,3)=arr(k,3)
        enddo
      enddo

c remove the blocks that are outliers:

      nct_min=nct*3/2

      nsam1=nsam
      changed=1
      do while (changed.eq.1 .and. nsam1.ge.nct_min)
        call find_slope_2D(nfit,nsam1,mean_sam,aa,bb,cc)

        do i=1,nsam1
          sam(i)=mean_sam(i,3)-aa-bb*mean_sam(i,1)-cc*mean_sam(i,2)
          tempo(i)=sam(i)
        enddo
        call sort(nsam1,nfit,tempo)
        mean=tempo(nsam1/2)
        sig=0.5*(tempo(nsam1*5/6)-tempo(nsam1/6))
        nsam=0
        changed=0
        do i=1,nsam1
          if (abs(sam(i)-mean).lt.3.*sig) then
            nsam=nsam+1
            mean_sam(nsam,1)=mean_sam(i,1)
            mean_sam(nsam,2)=mean_sam(i,2)
            mean_sam(nsam,3)=mean_sam(i,3)
          else
            changed=1
          endif
        enddo
        nsam1=nsam
      enddo

      if (nsam1.lt.nct_min) then
        write(*,*) 'Background not stable enough!',nsam1
        ierror=1
        return
      endif

      call fit_2D(nfit,nsam1,mean_sam,nct,ncx,c)

      do i=nx1,nx2
        do j=ny1,ny2

c numerical_fix F1: evaluate on the same normalized frame used to fit.
          x=(i-xmid)*xhalf_1
          y=(j-ymid)*yhalf_1

          image(i,j)=image(i,j)-func_val(x,y,nct,ncx,c)

        enddo
      enddo


      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine flatten_chip(nx1,nx2,ny1,ny2,npx,npy,array
     .,nct,ncx,ierror)
      implicit none

c The purpose of this code is to achieve a rough flattening of the background.

c input & output:
      integer nx1,nx2,ny1,ny2,npx,npy,ierror,nct,ncx
      real array(npx,npy)

c parameter:
      integer npp
      parameter (npp=1000)

c local variables:
      real pix(npp),arr(npp,3),arr2(npp,3),c(nct)
      real aa,bb,cc,arr_min,arr_max,x,y
      integer ix,iy,i,j,nr
c sub-region affine normalization (numerical_fix F1):
      real xmid,ymid,xhalf_1,yhalf_1
c Uses:
      real ran1,func_val

      if (ierror.eq.1) return

c numerical_fix F1: same defect and same cure as in set_background --
c the legacy x=i*ratio kept the absolute pixel index, so the second
c amplifier sat far from the origin and the {1,x,y,xy} basis degraded.
c The affine map onto [-1,1]x[-1,1] preserves the fitted surface.
      xmid=0.5*(nx1+nx2)
      ymid=0.5*(ny1+ny2)
      xhalf_1=2./real(max(nx2-nx1,1))
      yhalf_1=2./real(max(ny2-ny1,1))

c pick "npp" random positions:
      do i=1,npp
        ix=int(ran1()*(nx2-nx1)+nx1)
        iy=int(ran1()*(ny2-ny1)+ny1)
        pix(i)=array(ix,iy)
        arr2(i,1)=(ix-xmid)*xhalf_1
        arr2(i,2)=(iy-ymid)*yhalf_1
        arr2(i,3)=pix(i)
      enddo
      call sort(npp,npp,pix)

      if (pix(1).eq.pix(npp)) then
        ierror=1
        return
      endif

c remove the pixels of extreme values:
      arr_min=pix(npp/3)
      arr_max=pix(2*npp/3)
      nr=0
      do i=1,npp
        if (arr2(i,3).ge.arr_min.and.arr2(i,3).le.arr_max) then
          nr=nr+1
          arr(nr,1)=arr2(i,1)
          arr(nr,2)=arr2(i,2)
          arr(nr,3)=arr2(i,3)
        endif
      enddo

      if (nr.lt.npp/10) then
        ierror=1
        return
      endif

      call fit_2D(npp,nr,arr,nct,ncx,c)

      do i=nx1,nx2
        do j=ny1,ny2
c numerical_fix F1: evaluate on the same normalized frame used to fit.
          x=(i-xmid)*xhalf_1
          y=(j-ymid)*yhalf_1
          array(i,j)=array(i,j)-func_val(x,y,nct,ncx,c)
        enddo
      enddo


      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine set_sig(nx1,nx2,ny1,ny2,npx,npy,image
     .,weight,aa,bb,cc,ierror,sig_scale_in)
      implicit none

c ==========================================
c Function: Estimate, validate, and apply one amplifier noise plane
c Method: Find the image-density mode, clip a private symmetric bar,
c         refine it once with a provisional plane, and fit all surviving
c         triples. Validate and normalize this amplifier immediately;
c         no DQ image is read or required.
c ==========================================
      integer nx1,nx2,ny1,ny2,npx,npy,ierror
      integer weight(npx,npy)
      real image(npx,npy),aa,bb,cc,sig_scale_in

      include 'sig_para.inc'

      integer i,j,ib,jb,nbx,nby,nblock,ibin,imax,iq
      integer xmin,xmax,ymin,ymax,nval,ntri,width,height
      integer nhist,nbelow,iter,use_local,npossible,nuse,nmask
      integer score,scoremax,scorem,scorep
      integer hist(sig_hist_nbin)
      logical*1 src(npx,npy)
      double precision blkval(sig_block_max)
      double precision bmeds(sig_max_blocks),bsigs(sig_max_blocks)
      double precision bmed,pmed,slocal,dval
      double precision dx,dy,pixd
      double precision center_seed,sigma_seed,center,sigma0
      double precision hlo,hhi,hstep,denom,delta,target,cum,frac
      double precision qval,rawa,rawb,rawc,aad,bbd,ccd,tmin,tmax
      logical sig_finite_d
      external sig_finite_d

      aa=0.
      bb=0.
      cc=0.
      if (ierror.ne.0) return

      if (nx2.le.nx1 .or. ny2.le.ny1) then
        write(*,*) 'Error / set_sig invalid amplifier geometry'
        ierror=1
        return
      endif
      if (.not.sig_finite_d(dble(sig_scale_in))
     .    .or. sig_scale_in.le.0.) then
        write(*,*) 'Error / set_sig invalid sig_scale'
        ierror=1
        return
      endif
      if (sig_hist_nbin.lt.4 .or. mod(sig_hist_nbin,2).ne.0
     .    .or. sig_clip_niter.lt.1) then
        write(*,*) 'Error / set_sig invalid configured iteration'
        ierror=1
        return
      endif

c First pass: block medians seed the amp-wide mode and noise width.
      width=nx2-nx1+1
      height=ny2-ny1+1
      nbx=(width+sig_blocksize-1)/sig_blocksize
      nby=(height+sig_blocksize-1)/sig_blocksize
      if (nbx*nby.gt.sig_max_blocks) then
        write(*,*) 'Error / set_sig block table overflow',nbx*nby
        ierror=1
        return
      endif
      nblock=0

      do ib=1,nbx
        xmin=nx1+(ib-1)*width/nbx
        xmax=nx1+ib*width/nbx-1
        do jb=1,nby
          ymin=ny1+(jb-1)*height/nby
          ymax=ny1+jb*height/nby-1

          nval=0
          do i=xmin,xmax
            do j=ymin,ymax
              if (weight(i,j).le.0) cycle
              dval=dble(image(i,j))
              if (.not.sig_finite_d(dval)) then
                write(*,*) 'Error / set_sig nonfinite image pixel'
                ierror=1
                return
              endif
              nval=nval+1
              if (nval.gt.sig_block_max) then
                write(*,*) 'Error / set_sig block buffer overflow'
                ierror=1
                return
              endif
              blkval(nval)=dval
            enddo
          enddo
          if (nval.lt.sig_min_block_pixels) then
            write(*,*) 'Error / set_sig too few block pixels',nval
            ierror=1
            return
          endif
          call sort_doub(nval,sig_block_max,blkval)
          if (mod(nval,2).eq.0) then
            bmed=0.5d0*(blkval(nval/2)+blkval(nval/2+1))
          else
            bmed=blkval(nval/2+1)
          endif

          ntri=0
          do i=xmin,xmax-1
            do j=ymin,ymax-1
              if (weight(i,j).le.0) cycle
              if (weight(i+1,j).le.0) cycle
              if (weight(i,j+1).le.0) cycle
              dval=dble(image(i,j))
              dx=dval-dble(image(i+1,j))
              dy=dval-dble(image(i,j+1))
              pixd=0.5d0*(dx*dx+dy*dy)
              if (.not.sig_finite_d(pixd)) then
                write(*,*) 'Error / set_sig nonfinite block pix'
                ierror=1
                return
              endif
              ntri=ntri+1
              if (ntri.gt.sig_block_max) then
                write(*,*) 'Error / set_sig triple buffer overflow'
                ierror=1
                return
              endif
              blkval(ntri)=pixd
            enddo
          enddo
          if (ntri.lt.sig_min_block_triples) then
            write(*,*) 'Error / set_sig too few block triples',ntri
            ierror=1
            return
          endif
          call sort_doub(ntri,sig_block_max,blkval)
          if (mod(ntri,2).eq.0) then
            pmed=0.5d0*(blkval(ntri/2)+blkval(ntri/2+1))
          else
            pmed=blkval(ntri/2+1)
          endif
          if (.not.sig_finite_d(pmed)
     .        .or. pmed.le.dble(sig_plane_min)) then
            write(*,*) 'Error / set_sig invalid block median'
            ierror=1
            return
          endif

          slocal=sqrt(pmed/dble(sig_median_ratio))
          nblock=nblock+1
          bmeds(nblock)=bmed
          bsigs(nblock)=slocal
        enddo
      enddo

      if (nblock.lt.sig_min_blocks) then
        write(*,*) 'Error / set_sig too few valid blocks',nblock
        ierror=1
        return
      endif

c Median block statistics give deterministic mode seeds.
      call sort_doub(nblock,sig_max_blocks,bmeds)
      call sort_doub(nblock,sig_max_blocks,bsigs)
      if (mod(nblock,2).eq.0) then
        center_seed=0.5d0*(bmeds(nblock/2)+bmeds(nblock/2+1))
        sigma_seed=0.5d0*(bsigs(nblock/2)+bsigs(nblock/2+1))
      else
        center_seed=bmeds(nblock/2+1)
        sigma_seed=bsigs(nblock/2+1)
      endif
      if (.not.sig_finite_d(center_seed)
     .    .or. .not.sig_finite_d(sigma_seed)
     .    .or. sigma_seed.le.0d0) then
        write(*,*) 'Error / set_sig invalid mode seeds'
        ierror=1
        return
      endif

c Locate the highest-density brightness bin and interpolate its peak.
      do ibin=1,sig_hist_nbin
        hist(ibin)=0
      enddo
      hlo=center_seed-dble(sig_hist_range)*sigma_seed
      hhi=center_seed+dble(sig_hist_range)*sigma_seed
      hstep=(hhi-hlo)/dble(sig_hist_nbin)
      if (.not.sig_finite_d(hstep) .or. hstep.le.0d0) then
        write(*,*) 'Error / set_sig invalid histogram range'
        ierror=1
        return
      endif
      nhist=0
      do i=nx1,nx2
        do j=ny1,ny2
          if (weight(i,j).le.0) cycle
          dval=dble(image(i,j))
          if (.not.sig_finite_d(dval)) then
            write(*,*) 'Error / set_sig nonfinite mode pixel'
            ierror=1
            return
          endif
          if (dval.lt.hlo .or. dval.ge.hhi) cycle
          ibin=int((dval-hlo)/hstep)+1
          ibin=max(1,min(sig_hist_nbin,ibin))
          hist(ibin)=hist(ibin)+1
          nhist=nhist+1
        enddo
      enddo
c Smooth five adjacent bins before choosing the density maximum. This
c keeps the mode definition but prevents a single noisy narrow bin from
c moving the centre on otherwise uniform Gaussian data.
      imax=3
      scoremax=-1
      do ibin=3,sig_hist_nbin-2
        score=hist(ibin-2)+2*hist(ibin-1)+3*hist(ibin)
     .       +2*hist(ibin+1)+hist(ibin+2)
        if (score.gt.scoremax) then
          scoremax=score
          imax=ibin
        endif
      enddo
      if (nhist.lt.sig_min_lower_count
     .    .or. hist(imax).lt.sig_min_mode_count
     .    .or. imax.le.3 .or. imax.ge.sig_hist_nbin-2) then
        write(*,*) 'Error / set_sig unresolved brightness mode'
        ierror=1
        return
      endif
      scorem=hist(imax-3)+2*hist(imax-2)+3*hist(imax-1)
     .      +2*hist(imax)+hist(imax+1)
      score=hist(imax-2)+2*hist(imax-1)+3*hist(imax)
     .     +2*hist(imax+1)+hist(imax+2)
      scorep=hist(imax-1)+2*hist(imax)+3*hist(imax+1)
     .      +2*hist(imax+2)+hist(imax+3)
      denom=dble(scorem-2*score+scorep)
      delta=0d0
      if (denom.lt.0d0) then
        delta=0.5d0*dble(scorem-scorep)/denom
        delta=max(-0.5d0,min(0.5d0,delta))
      endif
      center=hlo+(dble(imax)-0.5d0+delta)*hstep

c Recenter a second histogram so the lower half ends at the mode.
      do ibin=1,sig_hist_nbin
        hist(ibin)=0
      enddo
      hlo=center-dble(sig_hist_range)*sigma_seed
      hhi=center+dble(sig_hist_range)*sigma_seed
      hstep=(hhi-hlo)/dble(sig_hist_nbin)
      do i=nx1,nx2
        do j=ny1,ny2
          if (weight(i,j).le.0) cycle
          dval=dble(image(i,j))
          if (dval.lt.hlo .or. dval.ge.hhi) cycle
          ibin=int((dval-hlo)/hstep)+1
          ibin=max(1,min(sig_hist_nbin,ibin))
          hist(ibin)=hist(ibin)+1
        enddo
      enddo
      nbelow=0
      do ibin=1,sig_hist_nbin/2
        nbelow=nbelow+hist(ibin)
      enddo
      if (nbelow.lt.sig_min_lower_count) then
        write(*,*) 'Error / set_sig too few lower-half pixels',nbelow
        ierror=1
        return
      endif
      target=dble(sig_lower_quantile)*dble(nbelow)
      cum=0d0
      iq=0
      qval=0d0
      do ibin=1,sig_hist_nbin/2
        if (iq.eq.0 .and. hist(ibin).gt.0
     .      .and. cum+dble(hist(ibin)).ge.target) then
          frac=(target-cum)/dble(hist(ibin))
          frac=max(0d0,min(1d0,frac))
          qval=hlo+(dble(ibin-1)+frac)*hstep
          iq=ibin
        endif
        cum=cum+dble(hist(ibin))
      enddo
      sigma0=center-qval
      if (iq.eq.0 .or. .not.sig_finite_d(center)
     .    .or. .not.sig_finite_d(sigma0)
     .    .or. sigma0.le.sqrt(dble(sig_plane_min))) then
        write(*,*) 'Error / set_sig invalid mode width'
        ierror=1
        return
      endif

c Alternate private image clipping and all-survivor triple fitting.
      rawa=0d0
      rawb=0d0
      rawc=0d0
      do iter=1,sig_clip_niter
        use_local=0
        if (iter.gt.1) use_local=1
        nmask=0
        call build_sig_private_mask(nx1,nx2,ny1,ny2,npx,npy
     .  ,image,weight,center,sigma0,use_local,rawa,rawb,rawc
     .  ,src,nmask,ierror)
        if (ierror.ne.0) then
          write(*,*) 'F6_DIAG mask',nx1,nx2,iter,center,sigma0
     .               ,nmask
          return
        endif
        npossible=0
        nuse=0
        call fit_sig_masked_plane(nx1,nx2,ny1,ny2,npx,npy
     .  ,image,weight,src,rawa,rawb,rawc,npossible,nuse,ierror)
        if (ierror.ne.0) then
          write(*,*) 'F6_DIAG fit',nx1,nx2,iter,center,sigma0
     .               ,nmask,nuse,npossible
          return
        endif
      enddo

c Apply the required publication scale only to the final raw fit.
      aad=rawa*dble(sig_scale_in)
      bbd=rawb*dble(sig_scale_in)
      ccd=rawc*dble(sig_scale_in)
      if (.not.sig_finite_d(aad)
     .    .or. .not.sig_finite_d(bbd)
     .    .or. .not.sig_finite_d(ccd)) then
        write(*,*) 'Error / set_sig nonfinite coefficients'
        ierror=1
        return
      endif

      call validate_sig_plane_d(nx1,nx2,ny1,ny2,aad,bbd,ccd
     .                         ,tmin,tmax,ierror)
      if (ierror.ne.0) return

      aa=real(aad)
      bb=real(bbd)
      cc=real(ccd)
      if (.not.sig_finite_d(dble(aa))
     .    .or. .not.sig_finite_d(dble(bb))
     .    .or. .not.sig_finite_d(dble(cc))) then
        aa=0.
        bb=0.
        cc=0.
        ierror=1
      endif
      if (ierror.eq.0) then
        call validate_sig_plane(nx1,nx2,ny1,ny2,npx,npy
     .  ,aa,bb,cc,ierror)
      endif
      if (ierror.eq.0) then
        call apply_sig_plane(nx1,nx2,ny1,ny2,npx,npy,image
     .  ,aa,bb,cc,ierror)
      endif
      if (ierror.ne.0) then
        aa=0.
        bb=0.
        cc=0.
      endif

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine build_sig_private_mask(nx1,nx2,ny1,ny2,npx,npy
     .,image,weight,center,sigma0,use_local,aa,bb,cc,src,nmask
     .,ierror)
      implicit none

c ==========================================
c Function: Build the image-only private F6 brightness mask
c Method: Symmetrically clip around the amp mode, optionally scale the
c         width by a prior plane, and dilate without changing weight.
c ==========================================
      include 'sig_para.inc'
      integer nx1,nx2,ny1,ny2,npx,npy,use_local,nmask,ierror
      integer weight(npx,npy)
      real image(npx,npy)
      logical*1 src(npx,npy)
      integer i,j,ii,jj,i1,i2,j1,j2,nbase
      double precision center,sigma0,aa,bb,cc
      double precision pcenter,plane,sigloc,dval,threshold
      logical sig_finite_d
      external sig_finite_d

      if (ierror.ne.0) return
      do i=nx1,nx2
        do j=ny1,ny2
          src(i,j)=.false.
        enddo
      enddo

      pcenter=1d0
      if (use_local.ne.0) then
        pcenter=aa+bb*0.5d0*dble(nx1+nx2)
     .            +cc*0.5d0*dble(ny1+ny2)
        if (.not.sig_finite_d(pcenter)
     .      .or. pcenter.le.dble(sig_plane_min)) then
          write(*,*) 'Error / set_sig invalid local mask plane'
          ierror=1
          return
        endif
      endif

      nbase=0
      do i=nx1,nx2
        do j=ny1,ny2
          if (weight(i,j).le.0) cycle
          nbase=nbase+1
          dval=dble(image(i,j))
          if (.not.sig_finite_d(dval)) then
            write(*,*) 'Error / set_sig nonfinite mask pixel'
            ierror=1
            return
          endif
          sigloc=sigma0
          if (use_local.ne.0) then
            plane=aa+bb*dble(i)+cc*dble(j)
            if (.not.sig_finite_d(plane)
     .          .or. plane.le.dble(sig_plane_min)) then
              write(*,*) 'Error / set_sig invalid local mask value'
              ierror=1
              return
            endif
            sigloc=sigma0*sqrt(plane/pcenter)
          endif
          threshold=dble(sig_clip_k)*sigloc
          if (abs(dval-center).le.threshold) cycle
          i1=max(nx1,i-sig_rdil)
          i2=min(nx2,i+sig_rdil)
          j1=max(ny1,j-sig_rdil)
          j2=min(ny2,j+sig_rdil)
          do ii=i1,i2
            do jj=j1,j2
              src(ii,jj)=.true.
            enddo
          enddo
        enddo
      enddo
      if (nbase.lt.sig_min_fit_triples) then
        write(*,*) 'Error / set_sig too few base pixels',nbase
        ierror=1
        return
      endif

      nmask=0
      do i=nx1,nx2
        do j=ny1,ny2
          if (weight(i,j).gt.0 .and. src(i,j)) nmask=nmask+1
        enddo
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fit_sig_masked_plane(nx1,nx2,ny1,ny2,npx,npy
     .,image,weight,src,aa,bb,cc,npossible,nuse,ierror)
      implicit none

c ==========================================
c Function: Fit the raw F6 2*sigma^2 plane from private-mask survivors
c Method: Accumulate every eligible triple in double precision, Jacobi
c         scale the normal equations, and solve by checked Cholesky.
c ==========================================
      include 'sig_para.inc'
      integer nx1,nx2,ny1,ny2,npx,npy,npossible,nuse,ierror
      integer weight(npx,npy)
      real image(npx,npy)
      logical*1 src(npx,npy)
      integer i,j,ii,jj
      double precision aa,bb,cc,dval,dx,dy,pixd
      double precision xmid,ymid,xhalf_1,yhalf_1,xn,yn
      double precision bas(3),cm(3,3),cv(3),dsc(3),sol(3)
      double precision tmin,tmax
      logical sig_finite_d
      external sig_finite_d

      if (ierror.ne.0) return
      xmid=0.5d0*dble(nx1+nx2)
      ymid=0.5d0*dble(ny1+ny2)
      xhalf_1=2d0/dble(nx2-nx1)
      yhalf_1=2d0/dble(ny2-ny1)
      do i=1,3
        cv(i)=0d0
        do j=1,3
          cm(i,j)=0d0
        enddo
      enddo
      npossible=0
      nuse=0

      do i=nx1,nx2-1
        do j=ny1,ny2-1
          if (weight(i,j).le.0) cycle
          if (weight(i+1,j).le.0) cycle
          if (weight(i,j+1).le.0) cycle
          npossible=npossible+1
          if (src(i,j)) cycle
          if (src(i+1,j)) cycle
          if (src(i,j+1)) cycle
          dval=dble(image(i,j))
          dx=dval-dble(image(i+1,j))
          dy=dval-dble(image(i,j+1))
          pixd=0.5d0*(dx*dx+dy*dy)
          if (.not.sig_finite_d(pixd)) then
            write(*,*) 'Error / set_sig nonfinite final pix'
            ierror=1
            return
          endif
          xn=(dble(i)-xmid)*xhalf_1
          yn=(dble(j)-ymid)*yhalf_1
          bas(1)=1d0
          bas(2)=xn
          bas(3)=yn
          do ii=1,3
            cv(ii)=cv(ii)+pixd*bas(ii)
            do jj=1,3
              cm(ii,jj)=cm(ii,jj)+bas(ii)*bas(jj)
            enddo
          enddo
          nuse=nuse+1
        enddo
      enddo

      if (npossible.lt.sig_min_fit_triples
     .    .or. nuse.lt.sig_min_fit_triples) then
        write(*,*) 'Error / set_sig too few final triples',nuse
        ierror=1
        return
      endif
      if (dble(nuse)/dble(npossible)
     .    .lt.dble(sig_min_fit_frac)) then
        write(*,*) 'Error / set_sig private mask too large',nuse
        ierror=1
        return
      endif

      do i=1,3
        if (.not.sig_finite_d(cm(i,i))
     .      .or. cm(i,i).le.0d0) then
          write(*,*) 'Error / set_sig invalid normal matrix'
          ierror=1
          return
        endif
        dsc(i)=1d0/sqrt(cm(i,i))
      enddo
      do i=1,3
        cv(i)=cv(i)*dsc(i)
        do j=1,3
          cm(i,j)=cm(i,j)*dsc(i)*dsc(j)
          if (.not.sig_finite_d(cm(i,j))) then
            write(*,*) 'Error / set_sig nonfinite normal matrix'
            ierror=1
            return
          endif
        enddo
      enddo
      call solve_sig_plane3(cm,cv,sol,ierror)
      if (ierror.ne.0) return

      bb=sol(2)*dsc(2)*xhalf_1
      cc=sol(3)*dsc(3)*yhalf_1
      aa=sol(1)*dsc(1)-bb*xmid-cc*ymid
      if (.not.sig_finite_d(aa) .or. .not.sig_finite_d(bb)
     .    .or. .not.sig_finite_d(cc)) then
        write(*,*) 'Error / set_sig nonfinite raw coefficients'
        ierror=1
        return
      endif
      call validate_sig_plane_d(nx1,nx2,ny1,ny2,aa,bb,cc
     .                         ,tmin,tmax,ierror)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine validate_sig_plane_d(nx1,nx2,ny1,ny2,aa,bb,cc
     .,tmin,tmax,ierror)
      implicit none

c ==========================================
c Function: Validate a double-precision linear F6 plane on an amplifier
c Method: Use the four rectangle corners for finite, positive extrema
c         and require the configured maximum-to-minimum ratio.
c ==========================================
      include 'sig_para.inc'
      integer nx1,nx2,ny1,ny2,ierror,i
      double precision aa,bb,cc,tmin,tmax,corner(4)
      logical sig_finite_d
      external sig_finite_d

      if (ierror.ne.0) return
      corner(1)=aa+bb*dble(nx1)+cc*dble(ny1)
      corner(2)=aa+bb*dble(nx2)+cc*dble(ny1)
      corner(3)=aa+bb*dble(nx1)+cc*dble(ny2)
      corner(4)=aa+bb*dble(nx2)+cc*dble(ny2)
      tmin=corner(1)
      tmax=corner(1)
      do i=1,4
        if (.not.sig_finite_d(corner(i))) then
          write(*,*) 'Error / set_sig nonfinite plane corner'
          ierror=1
          return
        endif
        tmin=min(tmin,corner(i))
        tmax=max(tmax,corner(i))
      enddo
      if (tmin.le.dble(sig_plane_min)) then
        write(*,*) 'Error / set_sig nonpositive noise plane',tmin,tmax
        ierror=1
        return
      endif
      if (tmax/tmin.gt.dble(sig_max_plane_ratio)) then
        write(*,*) 'Error / set_sig excessive plane variation',tmin,tmax
        ierror=1
      endif
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine validate_sig_plane(nx1,nx2,ny1,ny2,npx,npy
     .,aa,bb,cc,ierror)
      implicit none

c ==========================================
c Function: Preflight one real-coefficient F6 plane before image
c           mutation
c Method: Promote the coefficients and reuse the estimator's exact
c         corner validation before either amplifier is committed.
c ==========================================
      integer nx1,nx2,ny1,ny2,npx,npy,ierror
      real aa,bb,cc
      double precision tmin,tmax

      if (npx.lt.1 .or. npy.lt.1) then
        write(*,*) 'Error / set_sig invalid image bounds'
        ierror=1
        return
      endif
      call validate_sig_plane_d(nx1,nx2,ny1,ny2,dble(aa),dble(bb)
     .                         ,dble(cc),tmin,tmax,ierror)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function sig_finite_d(x)
      implicit none

c ==========================================
c Function: Check a double-precision F6 value without an IEEE module
c Method: Reject NaN through x.ne.x and reject infinities or overflowed
c         values with a conservative GCC-4.8-compatible finite bound.
c ==========================================
      double precision x,sig_finite_limit
      parameter (sig_finite_limit=1.797693134862315d308)

      sig_finite_d=.false.
      if (x.ne.x) return
      if (abs(x).gt.sig_finite_limit) return
      sig_finite_d=.true.
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine solve_sig_plane3(a,b,x,ierror)
      implicit none

c Function: Solve the scaled 3 by 3 F6 normal equations with status.
c Method: Cholesky rejects singular/nonfinite pivots; unlike the legacy
c         LU helper, it never pauses for input or fabricates one.
      include 'sig_para.inc'
      integer ierror,i,j,k
      double precision a(3,3),b(3),x(3),l(3,3),y(3),s
      logical sig_finite_d
      external sig_finite_d

      do i=1,3
        x(i)=0d0
        y(i)=0d0
        do j=1,3
          l(i,j)=0d0
        enddo
      enddo

      do i=1,3
        do j=1,i
          s=a(i,j)
          do k=1,j-1
            s=s-l(i,k)*l(j,k)
          enddo
          if (.not.sig_finite_d(s)) then
            write(*,*) 'Error / set_sig nonfinite Cholesky term'
            ierror=1
            return
          endif
          if (i.eq.j) then
            if (s.le.sig_pivot_min) then
              write(*,*) 'Error / set_sig singular normal matrix'
              ierror=1
              return
            endif
            l(i,j)=sqrt(s)
          else
            if (l(j,j).le.sig_pivot_min) then
              write(*,*) 'Error / set_sig invalid Cholesky pivot'
              ierror=1
              return
            endif
            l(i,j)=s/l(j,j)
          endif
        enddo
      enddo

      do i=1,3
        s=b(i)
        do k=1,i-1
          s=s-l(i,k)*y(k)
        enddo
        y(i)=s/l(i,i)
        if (.not.sig_finite_d(y(i))) then
          write(*,*) 'Error / set_sig nonfinite forward solve'
          ierror=1
          return
        endif
      enddo

      do i=3,1,-1
        s=y(i)
        do k=i+1,3
          s=s-l(k,i)*x(k)
        enddo
        x(i)=s/l(i,i)
        if (.not.sig_finite_d(x(i))) then
          write(*,*) 'Error / set_sig nonfinite backward solve'
          ierror=1
          return
        endif
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine apply_sig_plane(nx1,nx2,ny1,ny2,npx,npy,image
     .,aa,bb,cc,ierror)
      implicit none

c ==========================================
c Function: Apply a validated F6 noise plane to one amplifier.
c Method: After set_sig validation, divide the image by sqrt(plane/2)
c         in double precision.
c ==========================================
      integer nx1,nx2,ny1,ny2,npx,npy,ierror,i,j
      real image(npx,npy),aa,bb,cc
      double precision plane

      if (ierror.ne.0) return
      do i=nx1,nx2
        do j=ny1,ny2
          plane=dble(aa)+dble(bb)*dble(i)+dble(cc)*dble(j)
          image(i,j)=real(dble(image(i,j))/sqrt(0.5d0*plane))
        enddo
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine merge_defects(nx,ny,npx,npy
     .,weight,normap,area_max,source_thresh,area_thresh
     .,ierror)
      implicit none

      integer nx,ny,npx,npy,area_max,area_thresh,ierror
      real normap(npx,npy),source_thresh
      integer weight(npx,npy),mark(nx,ny)

      integer nb,nbb,toobig,buffer(area_max,2)
      integer i,j,ix,iy,jx,jy,u,v,k1,k2,k


      if (ierror.eq.1) return

      do i=1,nx
        do j=1,ny
          if (normap(i,j).ge.source_thresh.and.weight(i,j).eq.1) then
            mark(i,j)=1
          else
            mark(i,j)=0
          endif
        enddo
      enddo

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
                    elseif (mark(u,v).gt.1
     . .or. (mark(u,v).eq.0 .and. weight(u,v).eq.0)) then
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
                weight(buffer(k,1),buffer(k,2))=0
              enddo
            else
              do k=1,nb
                mark(buffer(k,1),buffer(k,2))=nb
              enddo
            endif
          endif
        enddo
      enddo


      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine locate_defects(nx,ny,npx,npy,array,normap,weight
     .,area_max,area_thresh,ierror)
      implicit none

      integer nx,ny,npx,npy,ierror,area_max,area_thresh
      real normap(npx,npy),map(nx,ny),array(npx,npy)
      integer weight(npx,npy)
      real diffx(nx,ny),diffy(nx,ny)
      integer i,j,ix,iy
      character filename*100

      integer margin
      parameter (margin=10)
      real defect_halo_thresh
      parameter (defect_halo_thresh=1.)
      integer y_smooth,x_smooth
      parameter (y_smooth=200)
      parameter (x_smooth=100)

      real sig,med,sigx,medx,sigy,medy

      real loga,iden
      external loga,iden

      do j=1,ny
        do i=nx/2-margin,nx/2+margin
          weight(i,j)=0
        enddo
        do i=1,margin
          weight(i,j)=0
          weight(nx+1-i,j)=0
        enddo
      enddo
      do i=1,nx
        do j=1,margin
          weight(i,j)=0
          weight(i,ny+1-j)=0
        enddo
      enddo

      if (ierror.eq.1) return

      do i=1,nx
        do j=1,ny
c          map(i,j)=loga(normap(i,j),1)
c          map(i,j)=log(normap(i,j)
          map(i,j)=loga(array(i,j),1)
        enddo
      enddo


      call remove_continuous(nx,ny,nx,ny,map,iden,4)
c      filename='map1.fits'
c      call writeimage(filename,nx,ny,nx,ny,map)
c      call get_sig_med(nx,ny,map,sig,med)
c      do i=1,nx
c        do j=1,ny
c          if (abs(map(i,j)).gt.5.*sig) weight(i,j)=0
c        enddo
c      enddo
      do i=1,nx
        ix=mod(i,nx)+1
        do j=1,ny
          iy=mod(j,ny)+1
          diffx(i,j)=map(i,j)-map(ix,j)
          diffy(i,j)=map(i,j)-map(i,iy)
        enddo
      enddo
      call get_sig_med(nx,ny,diffx,sigx,medx)
      do i=1,nx
        do j=1,ny
          if (abs(diffx(i,j)).gt.8.*sigx) weight(i,j)=0
        enddo
      enddo
      call get_sig_med(nx,ny,diffy,sigy,medy)
      do i=1,nx
        do j=1,ny
          if (abs(diffy(i,j)).gt.8.*sigy) weight(i,j)=0
        enddo
      enddo

      call mask_source_regions(nx,ny,npx,npy
     .,weight,normap,area_max,defect_halo_thresh*2.,area_thresh)

      call detect_stripes(nx,ny,npx,npy,normap,weight
     .,x_smooth,y_smooth)

      call detect_artificial_stripes(nx,ny,npx,npy,weight
     .,diffx,diffy,sigx,sigy,medx,medy)

      do i=1,nx
        do j=1,ny
          if (weight(i,j).gt.1) weight(i,j)=1
        enddo
      enddo

      call detect_stellar_halo(nx,ny,npx,npy,normap,weight
     .,area_max,defect_halo_thresh)

      call detect_dent(nx,ny,npx,npy,normap,weight
     .,area_max,defect_halo_thresh)



      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine detect_artificial_stripes(nx,ny,npx,npy,weight
     .,diffx,diffy,sigx,sigy,medx,medy)
      implicit none

      integer nx,ny,npx,npy
      integer weight(npx,npy)
      real diffx(nx,ny),diffy(nx,ny)
      real entropy(nx,ny),sig,med,sigx,sigy,medx,medy

      integer i,j

      call get_entropy(nx,ny,diffx,sigx,medx,2,entropy)
      call get_sig_med(nx,ny,entropy,sig,med)
c      call writeimage('entropy1.fits',nx,ny,nx,ny,entropy)
      do i=1,nx
        do j=1,ny
          if (weight(i,j).eq.1 .and. abs(entropy(i,j)-med).gt.10.*sig)
     . weight(i,j)=0
        enddo
      enddo

      call get_entropy(nx,ny,diffy,sigy,medy,2,entropy)
      call get_sig_med(nx,ny,entropy,sig,med)
c      call writeimage('entropy2.fits',nx,ny,nx,ny,entropy)
      do i=1,nx
        do j=1,ny
          if (weight(i,j).eq.1 .and. abs(entropy(i,j)-med).gt.10.*sig)
     . weight(i,j)=0
        enddo
      enddo

c      pause

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mask_source_regions(nx,ny,npx,npy
     .,weight,normap,area_max,source_thresh,area_thresh)
      implicit none

      integer nx,ny,npx,npy,area_max,area_thresh
      real normap(npx,npy),source_thresh
      integer weight(npx,npy),mark(nx,ny)

      integer nb,nbb,toobig,buffer(area_max,2)
      integer i,j,ix,iy,jx,jy,u,v,k1,k2,k

      do i=1,nx
        do j=1,ny
          if (normap(i,j).ge.source_thresh.and.weight(i,j).eq.1) then
            mark(i,j)=1
          else
            mark(i,j)=0
          endif
        enddo
      enddo

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
                do k=1,nb
                  ix=buffer(k,1)
                  iy=buffer(k,2)
                  mark(ix,iy)=nb
                  weight(ix,iy)=2
                enddo
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


      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ==========================================
c Function: Detect stripe outliers including partial edge blocks
c Method: Use ceiling block counts and sum/sqrt(nvalid), so full,
c         partial, and DQ-masked blocks share one noise scale.
c ==========================================
      subroutine detect_stripes(nx,ny,npx,npy,normap,weight
     .,x_smooth,y_smooth)
      implicit none

      integer nx,ny,npx,npy
      integer y_smooth,x_smooth
      real normap(npx,npy)
      real ymap(nx,(ny+y_smooth-1)/y_smooth)
      real xmap((nx+x_smooth-1)/x_smooth,ny)
      integer weight(npx,npy)
      integer yvalid(nx,(ny+y_smooth-1)/y_smooth)
      integer xvalid((nx+x_smooth-1)/x_smooth,ny)

      integer i,j,startj,endj,jj,numy,numx
      integer starti,endi,ii,nvalid
      real sig,med

      numy=(ny+y_smooth-1)/y_smooth
      do i=1,nx
        do j=1,numy
          ymap(i,j)=0.
          yvalid(i,j)=0
          nvalid=0
          startj=(j-1)*y_smooth+1
          endj=min(j*y_smooth,ny)
          do jj=startj,endj
            if (weight(i,jj).eq.1) then
              ymap(i,j)=ymap(i,j)+normap(i,jj)
              nvalid=nvalid+1
            endif
          enddo
          if (nvalid.gt.0) then
            ymap(i,j)=ymap(i,j)/sqrt(real(nvalid))
            yvalid(i,j)=1
          endif
        enddo
      enddo
      call get_stripe_sig_med(nx*numy,ymap,yvalid,sig,med)
      do i=1,nx
        do j=1,numy
          if (yvalid(i,j).eq.1 .and. sig.gt.0.
     .        .and. abs(ymap(i,j)-med).gt.sig*4.) then
            startj=(j-1)*y_smooth+1
            endj=min(j*y_smooth,ny)
            do jj=startj,endj
              weight(i,jj)=0
            enddo
          endif
        enddo
      enddo

      numx=(nx+x_smooth-1)/x_smooth
      do j=1,ny
        do i=1,numx
          xmap(i,j)=0.
          xvalid(i,j)=0
          nvalid=0
          starti=(i-1)*x_smooth+1
          endi=min(i*x_smooth,nx)
          do ii=starti,endi
            if (weight(ii,j).eq.1) then
              xmap(i,j)=xmap(i,j)+normap(ii,j)
              nvalid=nvalid+1
            endif
          enddo
          if (nvalid.gt.0) then
            xmap(i,j)=xmap(i,j)/sqrt(real(nvalid))
            xvalid(i,j)=1
          endif
        enddo
      enddo
      call get_stripe_sig_med(numx*ny,xmap,xvalid,sig,med)
      do j=1,ny
        do i=1,numx
          if (xvalid(i,j).eq.1 .and. sig.gt.0.
     .        .and. abs(xmap(i,j)-med).gt.sig*4.) then
            starti=(i-1)*x_smooth+1
            endi=min(i*x_smooth,nx)
            do ii=starti,endi
              weight(ii,j)=0
            enddo
          endif
        enddo
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ==========================================
c Function: Estimate stripe median and width from valid blocks
c Method: Preserve the 1000-sample robust estimator while drawing only
c         from blocks that contain at least one usable pixel.
c ==========================================
      subroutine get_stripe_sig_med(np,map,valid,sig,med)
      implicit none

      integer np,valid(np)
      real map(np),sig,med

      integer npp
      parameter (npp=1000)
      integer i,k,nuse
      real values(np),pix(npp),ran1

      nuse=0
      do i=1,np
        if (valid(i).eq.1) then
          nuse=nuse+1
          values(nuse)=map(i)
        endif
      enddo

      if (nuse.eq.0) then
        sig=-1.
        med=0.
        return
      endif

      do k=1,npp
        i=int(ran1()*real(nuse))+1
        pix(k)=values(i)
      enddo

      call sort(npp,npp,pix)
      sig=0.5*(pix(5*npp/6)-pix(npp/6))
      med=pix(npp/2)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine detect_stellar_halo(nx,ny,npx,npy,normap,weight
     .,npmax,defect_halo_thresh)
      implicit none

c The purpose of this code is to mark the parts affected by stellar halos.

c input & output:
      integer nx,ny,npx,npy,npmax
      real normap(npx,npy),defect_halo_thresh
      integer weight(npx,npy)
      real smoothed(nx,ny)

c local variables:
      integer dmark(nx,ny),buffer(npmax,2)
      integer nb,nbb,toobig,i,j,ix,iy,jx,jy,u,v,k1,k2,k


      do i=1,nx
        do j=1,ny
          smoothed(i,j)=normap(i,j)
        enddo
      enddo

      call smooth_image55(nx,ny,smoothed,1)


      do i=1,nx
        do j=1,ny
          if (smoothed(i,j).ge.defect_halo_thresh) then
            dmark(i,j)=1
          else
            dmark(i,j)=0
          endif
        enddo
      enddo

      do i=1,nx
        do j=1,ny
          if (dmark(i,j).eq.1) then
            nbb=0
            nb=1
            buffer(nb,1)=i
            buffer(nb,2)=j
            dmark(i,j)=-1
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
                    if (dmark(u,v).eq.1) then
                      nb=nb+1
                      buffer(nb,1)=u
                      buffer(nb,2)=v
                      dmark(u,v)=-1
                      if (nb.eq.npmax) then
                        toobig=1
                        goto 20
                      endif
                    elseif (dmark(u,v).gt.1) then
                      toobig=1
                      goto 20
                    endif
                  enddo
                enddo
              enddo
            enddo
20          if (toobig.eq.1) then
              do k=1,nb
                dmark(buffer(k,1),buffer(k,2))=npmax
                weight(buffer(k,1),buffer(k,2))=0
              enddo
            else
              do k=1,nb
                dmark(buffer(k,1),buffer(k,2))=nb
              enddo
            endif
          endif
        enddo
      enddo



      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine detect_dent(nx,ny,npx,npy,normap,weight
     .,npmax,defect_halo_thresh)
      implicit none

c input & output:
      integer nx,ny,npx,npy,npmax
      real normap(npx,npy),defect_halo_thresh
      integer weight(npx,npy)

c local variables:
      integer dmark(nx,ny),buffer(npmax,2)
      integer nb,nbb,toobig,i,j,ix,iy,jx,jy,u,v,k1,k2,k


      do i=1,nx
        do j=1,ny
          if (normap(i,j).le.-defect_halo_thresh) then
            dmark(i,j)=1
          else
            dmark(i,j)=0
          endif
        enddo
      enddo

      do i=1,nx
        do j=1,ny
          if (dmark(i,j).eq.1) then
            nbb=0
            nb=1
            buffer(nb,1)=i
            buffer(nb,2)=j
            dmark(i,j)=-1
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
                    if (dmark(u,v).eq.1) then
                      nb=nb+1
                      buffer(nb,1)=u
                      buffer(nb,2)=v
                      dmark(u,v)=-1
                      if (nb.eq.npmax) then
                        toobig=1
                        goto 20
                      endif
                    elseif (dmark(u,v).gt.1) then
                      toobig=1
                      goto 20
                    endif
                  enddo
                enddo
              enddo
            enddo
20          if (toobig.eq.1) then
              do k=1,nb
                dmark(buffer(k,1),buffer(k,2))=npmax
                weight(buffer(k,1),buffer(k,2))=0
              enddo
            else
              do k=1,nb
                dmark(buffer(k,1),buffer(k,2))=nb
              enddo
            endif
          endif
        enddo
      enddo



      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
