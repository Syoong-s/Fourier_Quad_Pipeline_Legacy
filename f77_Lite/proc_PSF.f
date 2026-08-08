      subroutine proc_PSF(iexpo)
      implicit none
      include 'para.inc'

      integer iexpo,nchip,ichip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT
      integer nc
      real p_chip(NMAX_cHIP,4)

      call get_image_list(iexpo,IMAGE_FILE,nchip,DIR_OUTPUT)

      call read_in_candidates(nchip,IMAGE_FILE,DIR_OUTPUT,nc,p_chip)

      call star_selection(nchip)

      call plot_star_expo(nchip,IMAGE_FILE,DIR_OUTPUT)

      call plot_stars(nchip,IMAGE_FILE,DIR_OUTPUT
     .,nc,p_chip)

      call make_PSF_local_fit(nchip,IMAGE_FILE,DIR_OUTPUT)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_in_candidates(nchip,IMAGE_FILE,DIR_OUTPUT
     .,nc,p_chip)
      implicit none
      include 'para.inc'

      integer nchip,nc
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT
      character*(strl) PREFIX,filename,headname

      integer nstar(NMAX_cHIP)
      double precision star_para(NMAX_cHIP,nstar_max,npara)
      common /star_info_pass/ star_para,nstar

      integer ierror,ntot,i,nn1,nn2,u,v,k,j
      double precision cRPIX(2),cD(2,2),PU(2,npd),cRVAL(2)
      double precision x,y,xx,yy,step
      real p_chip(NMAX_cHIP,4),aa(npara),source_p(ns,ns)
      real star(nstar_max,ns,ns),ee(2),size,temp,FWHM
      real map1(ns,ns),map2(ns,ns)

      real chi_d(NMAX_cHIP,nstar_max,nstar_max)
      common /chi_d_pass/ chi_d

      call get_PREFIX_expo(IMAGE_FILE(1),PREFIX)
      headname=trim(DIR_OUTPUT)//'/astrometry/'//trim(PREFIX)//'.head'

      nc=0
      do k=1,nchip

        nstar(k)=0

        ierror=0
        call read_astrometry_para(headname,k,cRPIX,cD,cRVAL,PU
     .,npd,ierror)

        if (ierror.eq.1) cycle
        nc=nc+1
        x=1d0
        y=1d0
        call xy_to_xxyy(x,y,xx,yy,cRPIX,cD)
        p_chip(nc,1)=xx
        p_chip(nc,2)=yy
        x=2046d0
        y=4094d0
        call xy_to_xxyy(x,y,xx,yy,cRPIX,cD)
        p_chip(nc,3)=xx
        p_chip(nc,4)=yy

        call get_PREFIX(IMAGE_FILE(k),PREFIX)
        PREFIX=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)
        filename=trim(PREFIX)//'_star_can_info.dat'
        open(unit=10,file=filename,status='old',iostat=ierror)
        rewind 10
        if (ierror.ne.0) then
          write(*,*) filename
          stop 'Error / PSF star_can_info catalog file error!!'
        endif
        read(10,*)
c        read(10,*) 'ig xp yp SNR'
        do while (ierror.ge.0)
          read(10,*,iostat=ierror) (aa(i),i=1,4)
          if (ierror.lt.0) cycle
          nstar(k)=nstar(k)+1
          do i=1,4
            star_para(k,nstar(k),i)=aa(i)
          enddo
        enddo
        close(10)

        if (nstar(k).gt.0) then
          nn1=ns*len_s
          nn2=ns*(int(nstar(k)/len_s)+1)
          filename=trim(PREFIX)//'_star_can_power.fits'
          call read_stamps(nstar_max,1,nstar(k),ns,ns,star
     .,nn1,nn2,filename)

          do i=1,nstar(k)
            star_para(k,i,5)=1.
            do u=1,ns
              do v=1,ns
                source_p(u,v)=star(i,u,v)
              enddo
            enddo

            x=star_para(k,i,2)
            y=star_para(k,i,3)
            call xy_to_xxyy(x,y,xx,yy,cRPIX,cD)
            star_para(k,i,6)=xx
            star_para(k,i,7)=yy
            call get_power_all(ns,ns,source_p,ee,size,0.02)
            star_para(k,i,8)=size
            star_para(k,i,9)=ee(1)
            star_para(k,i,10)=ee(2)
            call get_PSF_FWHM(source_p,FWHM)
            star_para(k,i,11)=FWHM

            temp=0.
            do u=1,ns
              do v=1,ns
                temp=temp+source_p(u,v)
              enddo
            enddo
            star_para(k,i,12)=1./temp
          enddo

          do i=1,nstar_max
            do j=1,nstar_max
              chi_d(k,i,j)=0.
            enddo
          enddo

          do i=1,nstar(k)-1
            do j=i+1,nstar(k)
              do u=1,ns
                do v=1,ns
                  map1(u,v)=star(i,u,v)*star_para(k,i,12)
                  map2(u,v)=star(j,u,v)*star_para(k,j,12)
                enddo
              enddo
              call ana_chi2(ns,map1,map2,temp)
              temp=sqrt(temp)
              chi_d(k,i,j)=temp
              chi_d(k,j,i)=temp
            enddo
          enddo
        endif

c        write(*,*) 'reading star cans:',ntot
      enddo

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine star_selection(nchip)
      implicit none
      include 'para.inc'

      integer nchip
      integer nstar(NMAX_cHIP)
      double precision star_para(NMAX_cHIP,nstar_max,npara)
      common /star_info_pass/ star_para,nstar

      integer npsam
      parameter (npsam=NMAX_cHIP*nstar_max*nstar_max)

      integer i,j,k,u,v,ntot,w,n,u2,v2
      real tmp(npsam),temp
      real thresh,peak,sig,chimin(NMAX_cHIP,nstar_max)
      real chi_d(NMAX_cHIP,nstar_max,nstar_max)
      common /chi_d_pass/ chi_d

      integer id(nstar_max),group_size(nstar_max)
      integer group_id(nstar_max),max_group_id
      integer max_gsize,max_gid,max2_gsize,max2_gid

      ntot=0
      do k=1,nchip
        ntot=ntot+nstar(k)
        do i=1,nstar_max
          chimin(k,i)=1000.
        enddo
      enddo
      if (ntot.lt.nstar_min*2) then
        do k=1,nchip
          do i=1,nstar(k)
            star_para(k,i,5)=-1
          enddo
        enddo
        return
      endif


      ntot=0
      do k=1,nchip
        do i=1,nstar(k)
          ntot=ntot+1
          tmp(ntot)=star_para(k,i,8)
        enddo
      enddo
      call sort(ntot,npsam,tmp)
      thresh=tmp((ntot*2)/3)

!----Determine the threshold for chi^2 ------------------------------

      ntot=0
      do k=1,nchip
        do i=1,nstar(k)-1
          do j=i+1,nstar(k)
            temp=chi_d(k,i,j)
            chimin(k,i)=min(temp,chimin(k,i))
            chimin(k,j)=min(temp,chimin(k,j))
            if (star_para(k,i,8).lt.thresh) cycle
            if (star_para(k,j,8).lt.thresh) cycle
            ntot=ntot+1
            tmp(ntot)=temp
          enddo
        enddo
      enddo

      call get_peak_width_low_side(npsam,ntot,tmp,peak,sig)

      thresh=peak+4.*sig

      do k=1,nchip
        ntot=0
        do i=1,nstar(k)
          if (chimin(k,i).gt.thresh) then
            star_para(k,i,5)=-1
            cycle
          endif
          ntot=ntot+1
        enddo
        if (ntot.lt.nstar_min_local) then
          do i=1,nstar(k)
            star_para(k,i,5)=-1
          enddo
          cycle
        endif

        ntot=0
        do i=1,nstar(k)
          if (star_para(k,i,5).lt.0) cycle
          ntot=ntot+1
          id(ntot)=i
        enddo
        do i=1,ntot
          group_id(i)=0
          group_size(i)=0
        enddo
        max_group_id=0
        do i=1,ntot
          if (group_id(i).eq.0) then
            max_group_id=max_group_id+1
            group_id(i)=max_group_id
          endif
          do j=i+1,ntot
            if (group_id(j).eq.group_id(i)) cycle
            if (chi_d(k,id(i),id(j)).gt.thresh) cycle
            if (group_id(j).eq.0) then
              group_id(j)=group_id(i)
            else
              u=group_id(j)
              v=group_id(i)
              do w=1,ntot
                if (group_id(w).eq.u) group_id(w)=v
              enddo
            endif
          enddo
        enddo
!-----------------Find the largest group ---------------------------
        do i=1,ntot
          u=group_id(i)
          group_size(u)=group_size(u)+1
        enddo

        if (nstar(k).gt.0) then
          max_gsize=group_size(1)
          max_gid=1
          max2_gsize=0
          max2_gid=0

          do i=2,max_group_id
            if (group_size(i).gt.max_gsize) then
              max2_gsize=max_gsize
              max2_gid=max_gid
              max_gsize=group_size(i)
              max_gid=i
            elseif (group_size(i).gt.max2_gsize) then
              max2_gsize=group_size(i)
              max2_gid=i
            endif
          enddo
        endif
c        if (max2_gsize.lt.nstar_min_local/2) then
          do i=1,ntot
            j=id(i)
            if (group_id(i).ne.max_gid) star_para(k,j,5)=-1
          enddo
c        else
c          do i=1,ntot
c            j=id(i)
c            if (group_id(i).ne.max_gid.and.group_id(i).ne.max2_gid)
c    .star_para(k,j,5)=-1
c          enddo
c        endif

      enddo

      do k=1,nchip
        n=0
        do i=1,nstar(k)
          if (star_para(k,i,5).lt.0) cycle
          n=n+1
        enddo
        if (n.lt.nstar_min_local) then
          do i=1,nstar(k)
            star_para(k,i,5)=-1
          enddo
        endif
      enddo


      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine plot_stars(nchip,IMAGE_FILE
     .,DIR_OUTPUT,nc,p_chip)
      implicit none
      include 'para.inc'

      integer nchip,nc
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT
      character*(strl) PREFIX,filename
      real p_chip(NMAX_cHIP,4)

      integer nstar(NMAX_cHIP)
      double precision star_para(NMAX_cHIP,nstar_max,npara)
      common /star_info_pass/ star_para,nstar

      integer nm
      parameter (nm=1000)
      real PSFmap(nm,nm),sk(NMAX_cHIP*nstar_max,5),source_p(ns,ns)
      integer nmax_stamp,opt(NMAX_cHIP*nstar_max)
      parameter (nmax_stamp=5000)

      integer i,j,k,u,v,ntot,w,nums,nn1,nn2
      real FWHM,FWHM_ave,e1_ave,e2_ave,chi_d_ave
      real chi_d(NMAX_cHIP,nstar_max,nstar_max)
      common /chi_d_pass/ chi_d


      call get_PREFIX_expo(IMAGE_FILE(1),PREFIX)
      PREFIX=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)

      filename=trim(PREFIX)//'_star_info_expo.dat'
      open(unit=10,file=filename,status='replace')
      rewind 10
      write(10,*) '# ichip nstar FWHM e1 e2 chi_d'

      ntot=0

      do k=1,nchip

        FWHM_ave=0.
        e1_ave=0.
        e2_ave=0.
        chi_d_ave=0.
        nums=0

        do i=1,nstar(k)
          if (star_para(k,i,5).le.0) cycle
          nums=nums+1
          ntot=ntot+1
          sk(ntot,1)=star_para(k,i,6)
          sk(ntot,2)=star_para(k,i,7)
          sk(ntot,3)=star_para(k,i,8)
          sk(ntot,4)=star_para(k,i,9)
          sk(ntot,5)=star_para(k,i,10)
          FWHM=star_para(k,i,11)

          FWHM_ave=FWHM_ave+FWHM
          e1_ave=e1_ave+star_para(k,i,9)
          e2_ave=e2_ave+star_para(k,i,10)
          if (nums.ge.2) chi_d_ave=chi_d_ave+chi_d(k,i,j)
          j=i
        enddo

        if (nums.ge.nstar_min_local) then
          FWHM_ave=FWHM_ave/nums
          e1_ave=e1_ave/nums
          e2_ave=e2_ave/nums
          chi_d_ave=chi_d_ave/(nums-1.)
          write(10,*) k,nums,FWHM_ave,e1_ave,e2_ave,chi_d_ave
        else
          write(10,*) k,0,-99.,-99.,-99.,-99.
        endif

      enddo

      ! write(*,*) trim(PREFIX),' total no. of stars:',ntot
      close(10)


      call draw_shear_expo(nm,PSFmap,nchip,NMAX_cHIP,nc,p_chip
     .,NMAX_cHIP*nstar_max,ntot,sk,200.,1.)
      filename=trim(PREFIX)//'_PSF_source.fits'
      call writeimage(filename,nm,nm,nm,nm,PSFmap)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine make_PSF_local_fit(nchip,IMAGE_FILE,DIR_OUTPUT)
      implicit none
      include 'para.inc'
      include 'cust_para.inc'

      integer nchip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT
      character*(strl) PREFIX,filename,PREFIX_e

      integer nstar(NMAX_cHIP)
      double precision star_para(NMAX_cHIP,nstar_max,npara)
      common /star_info_pass/ star_para,nstar

      integer ntot,nc,k,i,j,w,u,v,nums
      double precision posi(nstar_max,2)
      double precision xx,yy,PSF_coe_l(ns,ns,npl+1)
      real px,py,sshape(nstar_max,3),msshape(nstar_max,3)
      real model0(ns,ns),model(ns,ns),poly_ave,poly_std
      real ee(2),size,star(nstar_max,ns,ns)

      integer nn1,nn2
      real,allocatable :: poly_cochi2(:)

      call get_PREFIX_expo(IMAGE_FILE(1),PREFIX_e)
      filename=trim(DIR_OUTPUT)//'/result/'//trim(PREFIX_e)//
     .'_star_comp_expo.dat'
      open(unit=90,file=filename,status='replace')
      rewind 90

      ntot=0
      do k=1,nchip
        nums=0
        call get_PREFIX(IMAGE_FILE(k),PREFIX)

        if (nstar(k).gt.0) then
          nn1=ns*len_s
          nn2=ns*(int(nstar(k)/len_s)+1)
          filename=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)
     .//'_star_can_power.fits'
          call read_stamps(nstar_max,1,nstar(k),ns,ns,star
     .,nn1,nn2,filename)
        endif
        
        filename=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)//
     .'_PSF_coe_local.dat'
        open(unit=10,file=filename,status='replace')
        rewind 10

        do i=1,nstar(k)
          if (star_para(k,i,5).lt.0) cycle
          ntot=ntot+1
          nums=nums+1
          posi(nums,1)=star_para(k,i,2)
          posi(nums,2)=star_para(k,i,3)
          sshape(nums,1)=star_para(k,i,8)
          sshape(nums,2)=star_para(k,i,9)
          sshape(nums,3)=star_para(k,i,10)
          do u=1,ns
            do v=1,ns
              star(nums,u,v)=star(i,u,v)
            enddo
          enddo
        enddo

        if (nums.ge.nstar_min_local) then
          write(90,*) k,nums,1
          call ITP_norm_PSF(nums,nstar_max,star,posi,ns
     .,npl,nplx,chipnx,chipny,PSF_coe_l)
          allocate(poly_cochi2(nums))
          do i=1,nums
            xx=2.0 * (posi(i,1) / dble(chipnx)) - 1.0
            yy=2.0 * (posi(i,2) / dble(chipny)) - 1.0
            call get_PSF_model(ns,npl,nplx,PSF_coe_l,xx,yy,model,model0)
            call ana_chi2_simple(ns,model,model0,poly_cochi2(i))
            ! poly_cochi2(i) = MINVAL(model)
            call get_power_all(ns,ns,model,ee,size,0.02)
            msshape(i,1)=size
            msshape(i,2)=ee(1)
            msshape(i,3)=ee(2)
            px=posi(i,1)
            py=posi(i,2)
            write(90,*) px,py,(sshape(i,u),u=1,3),(msshape(i,v),v=1,3)
          enddo
          call get_array_ave_std(poly_cochi2,nums,poly_ave,poly_std)
          if (allocated(poly_cochi2)) deallocate(poly_cochi2)
          write(10,*) nums,1,poly_ave,poly_std
          do i=1,ns
            do j=1,ns
              write(10,*) (PSF_coe_l(i,j,u),u=1,(npl+1))
            enddo
          enddo

        else
          write(10,*) nums,-1,-1,-1
          write(90,*) k,nums,-1
        endif
        close(10)

      enddo

      close(90)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_PSF_model(ns,npp,nppx,PSF_coe,xx,yy,modelp,model0)
      implicit none

      integer ns,npp,nppx
      double precision xx,yy,PSF_coe(ns,ns,npp+1)
      real x,y,func_val,func_val_2
      real coep(npp),coe0(1),modelp(ns,ns),model0(ns,ns)
      integer i,j,k

      x=xx
      y=yy

      do i=1,ns
        do j=1,ns
          do k=1,npp
            coep(k)=PSF_coe(i,j,k)
          enddo
          coe0(1)=PSF_coe(i,j,npp+1)
          model0(i,j)=func_val_2(x,y,1,coe0)
          modelp(i,j)=func_val_2(x,y,npp,coep)
          if (isnan(modelp(i,j))) then
            modelp(1,1) = modelp(i,j)
            return
          endif
        enddo
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_PSF_FWHM(power,FWHM)
      implicit none
      include 'para.inc'

      real power(ns,ns),thresh
      integer i,j
      real area,FWHM,beta

      thresh=power(ns/2+1,ns/2+1)*exp(-1.)

      area=-1e-5
      do i=1,ns
        do j=1,ns
          if (power(i,j).ge.thresh) area=area+1.
        enddo
      enddo
      if (area.le.0.) then
        FWHM=0.0
        return
      endif
      beta=ns/(2.*pi)/sqrt(area/pi)
      FWHM=beta*2.*sqrt(2.*log(2.))*pixel_size

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_power_all(nx,ny,power,e,size,thresh_ratio)
      implicit none

      integer nx,ny
      real e(2),thresh_ratio,power(nx,ny),size
      integer area

c      integer i,j,cx,cy,mark(nx,ny),stack(nx*ny,2),area0
c      integer u,v,x,y,tempi
c      real kx,ky,thresh,norm

      call get_power_area(nx,ny,power,area,thresh_ratio)
      size=area
      call get_power_e(nx,ny,power,e,thresh_ratio)

c      cx=nx/2+1
c      cy=ny/2+1
c      thresh=power(cx,cy)*thresh_ratio
c      e(1)=0.
c      e(2)=0.
c      norm=0.
c      size=0.

c      do i=1,nx
c        do j=1,ny
c          mark(i,j)=0
c        enddo
c      enddo

c      area0=0
c      mark(cx,cy)=1
c      area=1
c      stack(area,1)=cx
c      stack(area,2)=cy

c      do while (area.gt.area0)
c        tempi=area
c        do i=area0+1,tempi
c          x=stack(i,1)
c          y=stack(i,2)
c          do u=max(x-1,1),min(x+1,nx)
c            do v=max(y-1,1),min(y+1,ny)
c              if (mark(u,v).eq.0.and.power(u,v).ge.thresh) then
c                mark(u,v)=1
c                area=area+1
c                stack(area,1)=u
c                stack(area,2)=v
c
c                kx=u-cx
c                ky=v-cy
c                e(1)=e(1)+power(u,v)*(kx*kx-ky*ky)
c                e(2)=e(2)+power(u,v)*2.*kx*ky
c                norm=norm+power(u,v)*(kx*kx+ky*ky)
c                size=size+power(u,v)

c              endif
c            enddo
c          enddo
c        enddo
c        area0=tempi
c      enddo
c
c      e(1)=e(1)/norm
c      e(2)=e(2)/norm
c      size=norm/size

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_power_area(nx,ny,power,area,thresh_ratio)
      implicit none

      integer nx,ny,area
      real thresh_ratio
      real power(nx,ny),thresh,de,rs_2,flx,temp,dx,dy,dx2,dy2,r
      integer i,j,cx,cy,mark(nx,ny),stack(nx*ny,2),area0
      integer u,v,x,y,tempi

      cx=nx/2+1
      cy=ny/2+1

      thresh=power(cx,cy)*thresh_ratio

      do i=1,nx
        do j=1,ny
          mark(i,j)=0
        enddo
      enddo

      area0=0
      mark(cx,cy)=1
      area=1
      stack(area,1)=cx
      stack(area,2)=cy

      do while (area.gt.area0)
        tempi=area
        do i=area0+1,tempi
          x=stack(i,1)
          y=stack(i,2)
          do u=max(x-1,1),min(x+1,nx)
            do v=max(y-1,1),min(y+1,ny)
              if (mark(u,v).eq.0 .and. power(u,v).ge.thresh) then
                mark(u,v)=1
                area=area+1
                stack(area,1)=u
                stack(area,2)=v
              endif
            enddo
          enddo
        enddo
        area0=tempi
      enddo

      area=(area-1)/2

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_power_e(nx,ny,power,e,thresh_ratio)
      implicit none

      integer nx,ny,area
      real e(2),norm,thresh_ratio,power(nx,ny),thresh
      integer i,j,cx,cy,mark(nx,ny),stack(nx*ny,2),area0
      integer u,v,x,y,tempi
      real kx,ky

      cx=nx/2+1
      cy=ny/2+1

      thresh=power(cx,cy)*thresh_ratio

      e(1)=0.
      e(2)=0.
      norm=0.

      do i=1,nx
        do j=1,ny
          mark(i,j)=0
        enddo
      enddo

      area0=0
      mark(cx,cy)=1
      area=1
      stack(area,1)=cx
      stack(area,2)=cy

      do while (area.gt.area0)
        tempi=area
        do i=area0+1,tempi
          x=stack(i,1)
          y=stack(i,2)
          do u=max(x-1,1),min(x+1,nx)
            do v=max(y-1,1),min(y+1,ny)
              if (mark(u,v).eq.0 .and. power(u,v).ge.thresh) then
                mark(u,v)=1
                area=area+1
                stack(area,1)=u
                stack(area,2)=v

                kx=u-cx
                ky=v-cy
                e(1)=e(1)+power(u,v)*(kx*kx-ky*ky)
                e(2)=e(2)+power(u,v)*2.*kx*ky
                norm=norm+power(u,v)*(kx*kx+ky*ky)

              endif
            enddo
          enddo
        enddo
        area0=tempi
      enddo

      e(1)=e(1)/norm
      e(2)=e(2)/norm

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine plot_star_expo(nchip,IMAGE_FILE,DIR_OUTPUT)
      implicit none
      include 'para.inc'

      integer nchip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT
      character*(strl) PREFIX,filename

      integer nstar(NMAX_cHIP)
      double precision star_para(NMAX_cHIP,nstar_max,npara)
      real star_test(NMAX_cHIP*nstar_max,ns,ns)
      common /star_info_pass/ star_para,nstar

      integer nmax_stamp,opt(NMAX_cHIP*nstar_max)
      parameter (nmax_stamp=5000)

      integer ntot ,ichip, w, nn1, nn2, start
      integer k, i, u, v
      real star(nstar_max,ns,ns)

      real chi_d(NMAX_cHIP,nstar_max,nstar_max)
      common /chi_d_pass/ chi_d

      do i=1,nstar_max*NMAX_cHIP
        opt(i)=0
      enddo

      ntot = 0
      w = 0
      start = 0
      do ichip=1,nchip
        if (nstar(ichip).eq.0) cycle
        nn1=ns*len_s
        nn2=ns*(int(nstar(ichip)/len_s)+1)

        call get_PREFIX(IMAGE_FILE(ichip),PREFIX)
        PREFIX=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)
        filename=trim(PREFIX)//'_star_can_power.fits'
        call read_stamps(nstar_max,1,nstar(ichip),ns,ns
     .,star,nn1,nn2,filename)
        do i=1,nstar(ichip)
          do u=1,ns
            do v=1,ns
              star_test(start+i,u,v)=star(i,u,v)
            enddo
          enddo
          w = start + i
          if (star_para(ichip,i,5).le.0) cycle
          ntot = ntot + 1
          if (ntot.lt.nmax_stamp) opt(w)=1
        enddo
        start = start + nstar(ichip)
      enddo

      call get_PREFIX_expo(IMAGE_FILE(1),PREFIX)
      PREFIX=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)

      if (ntot.gt.0) then 
        nn1=ns*len_sam
        nn2=ns*(int(min(ntot,nmax_stamp)/len_sam)+1)
        filename=trim(PREFIX)//'_star_power_expo.fits'
        call write_stamps_2(nstar_max*NMAX_cHIP,w,ns,ns,star_test
     .,opt,1,nn1,nn2,filename)
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

