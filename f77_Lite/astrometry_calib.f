c ==========================================
c Function: Synced astrometry routine gen_astrometry_data
c Method: Preserve the current full-f77 numerical implementation.
c ==========================================
      subroutine gen_astrometry_data(cat_standard,nx,ny,npx,npy
     .,map,weight,cRPIX,cD,cRVAL,filename,proc_error)
      implicit none

      real astrometry_shift_ratio
      parameter (astrometry_shift_ratio=0.2)

      character cat_standard*(*),filename*(*)
      integer nx,ny,npx,npy,ierror,n_user,n_ref,proc_error
      real map(npx,npy),flux_min
      integer weight(npx,npy)
      integer i,j,u,v,k
      double precision ra,dra,diffra,dec(2),a,d,flux,mag
      integer nss_max,nss,n_user_max
      parameter (nss_max=10000)
      parameter (n_user_max=200)

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
      open(unit=20,file=trim(cat_standard)
     .,status='old',iostat=ierror)
      rewind 20
      if (ierror.ne.0) then
        write(*,'(A)') cat_standard
        write(*,*) ierror
        write(*,*) 'Error / gen_astrometry_data catalog file error!!'
        proc_error = 1
        nss=0
        n_user=0
        n_ref=0
        goto 40
      endif
      read(20,*)
      do while (ierror.ge.0)
        read(20,*,iostat=ierror) a,d
        if (ierror.lt.0) cycle
        if (abs(diffra(a,ra)).gt.dra*0.5) cycle
        if (d.lt.dec(1).or.d.gt.dec(2)) cycle
        n_ref=n_ref+1
        if (n_ref.gt.nss_max) then
          write(*,*) 'n_ref is too large!!',filename
          close(20)
          nss=0
          n_user=0
          n_ref=0
          goto 40
        endif
        ra_r(n_ref)=a
        dec_r(n_ref)=d
        call coordinate_transfer_simple(a,d,xr(n_ref),yr(n_ref),-1
     .,cRPIX,cD,cRVAL)
      enddo
      close(20)


      call get_astrometry_catalog(nx,ny,npx,npy,map
     .,weight,n_user_max,nss_max,nss,xs,ys)


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

      ! write(*,*) nss,n_user,n_ref,trim(filename)

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
c Function: Synced astrometry routine get_astrometry_catalog
c Method: Preserve the current full-f77 numerical implementation.
c ==========================================
      subroutine get_astrometry_catalog(nx,ny,npx,npy,image
     .,weight,n_user_max,ns_max,ns,xs,ys)
      implicit none

      integer nx,ny,npx,npy,n_user_max
      real image(npx,npy),flux
      integer mark(nx,ny),weight(npx,npy)

      integer sizemax
      parameter (sizemax=10000)
      integer nb,nbb,buffer(sizemax,2)
      integer i,j,ix,iy,jx,jy,u,v,k1,k2,k,xp,yp
      real temp,peak,xc,yc

      integer ns_max,ns,toobig,nsb
      double precision xs(ns_max),ys(ns_max)
      double precision xsb(ns_max),ysb(ns_max)

      real flux_array(ns_max)
      integer order(ns_max)

      integer area,area_limit
      parameter (area_limit=400)


      do i=1,nx
        do j=1,ny
          if (image(i,j).ge.5. .and. weight(i,j).gt.0) then
            mark(i,j)=1
          else
            mark(i,j)=0
          endif
        enddo
      enddo

      nsb=0

      do i=1,nx
        do j=1,ny
          if (mark(i,j).eq.1) then
            nbb=0
            nb=1
            buffer(nb,1)=i
            buffer(nb,2)=j
            mark(i,j)=-1
            xp=i
            yp=j
            peak=image(i,j)
            flux=0.
            toobig=0
            do while (nb.gt.nbb)
              k1=nbb+1
              k2=nb
              nbb=nb
              do k=k1,k2
                ix=buffer(k,1)
                iy=buffer(k,2)
                do u=max(ix-3,1),min(ix+3,nx)
                  do v=max(iy-3,1),min(iy+3,ny)
                    if (mark(u,v).eq.1) then
                      nb=nb+1
                      buffer(nb,1)=u
                      buffer(nb,2)=v
                      mark(u,v)=-1
                      flux=flux+image(u,v)
                      if (image(u,v).gt.peak) then
                        peak=image(u,v)
                        xp=u
                        yp=v
                      endif
                      if (nb.eq.sizemax) then
                        toobig=1
                        goto 20
                      endif
                    elseif (mark(u,v).eq.2) then
                      toobig=1
                      goto 20
                    endif
                  enddo
                enddo
              enddo
            enddo
20          if (toobig.eq.1) then
              do k=1,nb
                ix=buffer(k,1)
                iy=buffer(k,2)
                mark(ix,iy)=2
              enddo
            else
              area=0
              flux=0.
              xc=0.
              yc=0.
              temp=peak*0.5
              do k=1,nb
                ix=buffer(k,1)
                iy=buffer(k,2)
                if (image(ix,iy).ge.temp) then
                  area=area+1
                  flux=flux+image(ix,iy)
                  xc=xc+image(ix,iy)*ix
                  yc=yc+image(ix,iy)*iy
                endif
              enddo
              xc=xc/flux
              yc=yc/flux

              if (area.le.area_limit) then
                nsb=nsb+1
                xsb(nsb)=xc
                ysb(nsb)=yc
                flux_array(nsb)=flux
                order(nsb)=nsb
                if (nsb.eq.ns_max) then
                  ns=0
                  return
                endif
              endif
            endif
          endif
        enddo
      enddo

      call indexx(nsb,ns_max,flux_array,order)

      ns=0
      i=nsb
      do while (ns.lt.n_user_max.and.i.ge.1)
        ns=ns+1
        xs(ns)=xsb(order(i))
        ys(ns)=ysb(order(i))
c        write(*,*) ns,flux_array(order(i)),xs(ns),ys(ns)
        i=i-1
      enddo
c      pause

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ==========================================
c Function: Synced astrometry routine pattern_matching
c Method: Preserve the current full-f77 numerical implementation.
c ==========================================
      subroutine pattern_matching(np0,n0,x0,y0,np1,n1,x1,y1
     .,shift_range,box_final)
      implicit none

      integer np0,np1,n0,n1,shift_range,n_match
      double precision x0(np0),y0(np0),x1(np1),y1(np1)

      integer radius,r_fof
      parameter (radius=40)
      parameter (r_fof=1)

      integer box(n0,0:n1),i,j,k
      integer box_final(np0),u,v,dx,dy,xp,yp

      real mark1(-shift_range:shift_range,-shift_range:shift_range)
      real mark2(-shift_range:shift_range,-shift_range:shift_range)
      integer mark3(-shift_range:shift_range,-shift_range:shift_range)

      integer nn,dpl,changed,xmin,xmax,ymin,ymax
      real thresh,peak,temp,norm

      real map(2*shift_range+1,2*shift_range+1)
      real arr((2*shift_range+1)*(2*shift_range+1))

c      character filename*100
      double precision const
      double precision xx(np0,2),xxt(np0,2),coe(2,3)

      integer crowdy
      parameter (crowdy=1)

      nn=shift_range
      do i=-nn,nn
        do j=-nn,nn
          mark1(i,j)=0.
          mark2(i,j)=0.
          mark3(i,j)=0
        enddo
      enddo

      do i=1,n0
        box(i,0)=0
        do j=1,n1
          dx=int(x1(j)-x0(i)+0.5)
          dy=int(y1(j)-y0(i)+0.5)
          if (dx.ge.-nn.and.dx.le.nn.and.dy.ge.-nn.and.dy.le.nn) then
            mark1(dx,dy)=mark1(dx,dy)+1.
            box(i,0)=box(i,0)+1
            box(i,box(i,0))=j
          endif
        enddo
      enddo

      do i=-nn,nn
        do j=-nn,nn
          mark2(i,j)=mark1(i,j)
        enddo
      enddo

c      do i=1,2*nn+1
c        do j=1,2*nn+1
c          map(i,j)=mark1(i-nn-1,j-nn-1)
c        enddo
c      enddo
c      filename='./pattern.fits'
c      call writeimage(filename,2*nn+1,2*nn+1,2*nn+1,2*nn+1,map)
c      pause


      if (radius.gt.0) then
        const=(radius*0.25d0)**(-2)
        do i=-nn,nn
          do j=-nn,nn
            if (mark1(i,j).gt.0) then
              do u=max(i-radius,-shift_range),min(i+radius,shift_range)
                do v=max(j-radius,-shift_range),
     .min(j+radius,shift_range)
                  mark2(u,v)=mark2(u,v)
     .+mark1(i,j)/(((i-u)**2+(j-v)**2)*const+1.)
                enddo
              enddo
            endif
          enddo
        enddo

c        do i=1,2*nn+1
c          do j=1,2*nn+1
c            map(i,j)=mark2(i-nn-1,j-nn-1)
c          enddo
c        enddo
c        filename='./pattern_smoothed.fits'
c        call writeimage(filename,2*nn+1,2*nn+1,2*nn+1,2*nn+1,map)
      endif

      peak=0.
      xp=0
      yp=0
      u=0
      do i=-nn+radius,nn-radius
        do j=-nn+radius,nn-radius
          u=u+1
          arr(u)=mark2(i,j)
          if (mark2(i,j).gt.peak) then
            peak=mark2(i,j)
            xp=i
            yp=j
          endif
        enddo
      enddo

      call sort(u,(2*nn+1)**2,arr)

c      thresh=int(arr(u/2))*3+1

      thresh=(arr(u/2)+peak)*0.5

      changed=1
      xmin=xp
      xmax=xp
      ymin=yp
      ymax=yp
      mark3(xp,yp)=1
      do while (changed.eq.1)
        changed=0
        do i=xmin,xmax
          do j=ymin,ymax
            if (mark3(i,j).eq.1) then
              do u=i-r_fof,i+r_fof
                do v=j-r_fof,j+r_fof
                  if (mark3(u,v).eq.0.and.mark2(u,v).gt.thresh) then
                    mark3(u,v)=1
                    xmin=min(xmin,u)
                    xmax=max(xmax,u)
                    ymin=min(ymin,v)
                    ymax=max(ymax,v)
                    changed=1
                  endif
                enddo
              enddo
            endif
          enddo
        enddo
      enddo


      do i=1,n0
        box_final(i)=0
        do j=1,box(i,0)
          k=box(i,j)

          dx=int(x1(k)-x0(i)+0.5)
          dy=int(y1(k)-y0(i)+0.5)
          if (dx.ge.-nn.and.dx.le.nn.and.dy.ge.-nn.and.dy.le.nn) then
            if (mark3(dx,dy).eq.1) then
              if (box_final(i).eq.0) then
                box_final(i)=k
              else
                box_final(i)=0
                goto 30
              endif
            endif
          endif
        enddo
30    enddo

      n_match=0
      do i=1,n0
        if (box_final(i).eq.0) cycle
        dpl=0
        do j=i+1,n0
          if (box_final(j).eq.0) cycle
          if (box_final(i).eq.box_final(j)) then
            box_final(j)=0
            dpl=1
          endif
        enddo
        if (dpl.eq.1) then
          box_final(i)=0
          cycle
        endif
        n_match=n_match+1
        xx(n_match,1)=x0(i)
        xx(n_match,2)=y0(i)
        xxt(n_match,1)=x1(box_final(i))
        xxt(n_match,2)=y1(box_final(i))

      enddo

      if (n_match.lt.6) then
        do i=1,n0
          box_final(i)=0
        enddo
        n_match=0
        return
      endif

      call fit_linear_2D(n_match,np0,xx,xxt,coe)

      do i=-nn,nn
        do j=-nn,nn
          mark1(i,j)=0.
          mark3(i,j)=0
        enddo
      enddo

      do i=1,n0
        do j=1,n1
          dx=int(x1(j)-(x0(i)*coe(1,1)+y0(i)*coe(1,2)+coe(1,3))+0.5)
          dy=int(y1(j)-(x0(i)*coe(2,1)+y0(i)*coe(2,2)+coe(2,3))+0.5)
          if (dx.ge.-nn.and.dx.le.nn.and.dy.ge.-nn.and.dy.le.nn) then
            mark1(dx,dy)=mark1(dx,dy)+1.
          endif
        enddo
      enddo

c      do i=1,2*nn+1
c        do j=1,2*nn+1
c          map(i,j)=mark1(i-nn-1,j-nn-1)
c        enddo
c      enddo
c      filename='./pattern_corrected.fits'
c      call writeimage(filename,2*nn+1,2*nn+1,2*nn+1,2*nn+1,map)
c      pause

      if (crowdy.eq.1) then
c--------------------------------------------------------------
c        FOR HIGH SOURcE DENSITY:
        u=0
        do i=-nn,nn
          do j=-nn,nn
            u=u+1
            arr(u)=mark1(i,j)
          enddo
        enddo
        call sort(u,u,arr)

        peak=mark1(0,0)
        thresh=(arr(u/2)+peak)*0.5
c---------------------------------------------------------------
      else
c--------------------------------------------------------------
c        FOR LOW SOURcE DENSITY:
        thresh=0.5
c--------------------------------------------------------------
      endif

      changed=1
      xmin=0
      xmax=0
      ymin=0
      ymax=0
      mark3(0,0)=1
      do while (changed.eq.1)
        changed=0
        do i=xmin,xmax
          do j=ymin,ymax
            if (mark3(i,j).eq.1) then
              do u=i-r_fof,i+r_fof
                do v=j-r_fof,j+r_fof
                  if (mark3(u,v).eq.0 .and. mark1(u,v).gt.thresh) then
                    mark3(u,v)=1
                    xmin=min(xmin,u)
                    xmax=max(xmax,u)
                    ymin=min(ymin,v)
                    ymax=max(ymax,v)
                    changed=1
                  endif
                enddo
              enddo
            endif
          enddo
        enddo
      enddo
cc--------------------------------------------------------------------------
      do i=1,n0
        box_final(i)=0
      enddo
      n_match=0
      do i=1,n0
        do j=1,n1
          dx=int(x1(j)-(x0(i)*coe(1,1)+y0(i)*coe(1,2)+coe(1,3))+0.5)
          dy=int(y1(j)-(x0(i)*coe(2,1)+y0(i)*coe(2,2)+coe(2,3))+0.5)
          if (dx.ge.-nn.and.dx.le.nn.and.dy.ge.-nn.and.dy.le.nn) then
            if (mark3(dx,dy).eq.1) then
              n_match=n_match+1
              box_final(i)=j
              goto 50
            endif
          endif
        enddo
50    enddo


      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ==========================================
c Function: Synced astrometry routine check_astrometry_global
c Method: Preserve the current full-f77 numerical implementation.
c ==========================================
      subroutine check_astrometry_global(np,n,nc,ra,dec,x,y
     .,cRPIX,cD,cRVAL,PU,npd,valid)
      implicit none

      integer np,nc,npd,ic
      integer valid(nc),n(nc),i
      double precision x(nc,np),y(nc,np),ra(nc,np),dec(nc,np)
      double precision xx(np),yy(np),a(np),d(np),ra_c,dra,dec_c,ddec
      double precision cRPIX(nc,2),cD(nc,2,2),cRVAL(2),PU(2,npd)
      double precision tolerate_shift,aa,dd,crp(2),cdd(2,2),diffra
      parameter (tolerate_shift=0.1d0)

      do ic=1,nc

        if (valid(ic).eq.0) cycle

        do i=1,n(ic)
          a(i)=ra(ic,i)
          d(i)=dec(ic,i)
        enddo

        call get_ra_dec_bound(np,n(ic),a,d,ra_c,dra,dec_c,ddec)

        crp(1)=cRPIX(ic,1)
        crp(2)=cRPIX(ic,2)

        cdd(1,1)=cD(ic,1,1)
        cdd(1,2)=cD(ic,1,2)
        cdd(2,1)=cD(ic,2,1)
        cdd(2,2)=cD(ic,2,2)

        do i=1,n(ic)

          call coordinate_transfer_PU(aa,dd,x(ic,i),y(ic,i)
     .,1,crp,cdd,cRVAL,PU,npd)
          if (abs(diffra(aa,ra_c)).gt.dra*(0.5+tolerate_shift)
     . .or. abs(dd-dec_c).gt.ddec*(0.5+tolerate_shift)
     . .or. isnan(aa).or.isnan(dd)) then
            valid(ic)=0
            exit
          endif
        enddo

      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ==========================================
c Function: Synced astrometry routine measure_astrometry_global
c Method: Preserve the current full-f77 numerical implementation.
c ==========================================
      subroutine measure_astrometry_global(np,n,nc,ra,dec,x,y
     .,cRPIX,cD,cRVAL,PU,npd,valid)
      implicit none

      integer np,nc,npd,ic
      integer valid(nc),n(nc)
      double precision x(nc,np),y(nc,np),ra(nc,np),dec(nc,np)
      double precision xi(nc,np),eta(nc,np)
      double precision cRPIX(nc,2),cD(nc,2,2),cRVAL(2),PU(2,npd)

      integer i,j,k,l,u,v,nn,ntot,order,px,py


      double precision vec(npd+nc*3),vec1(npd+nc*3),vec2(npd+nc*3)
      double precision bvec1(npd+nc*3),bvec2(npd+nc*3)
      double precision matx(npd+nc*3,npd+nc*3)
      double precision r
c numerical_fix F4 workspace (matx_1 is no longer needed, see below):
      double precision dscal(npd+nc*3),dsign
      integer indx(npd+nc*3)


      ntot=npd+nc*3

      do u=1,ntot
        do v=1,ntot
          matx(u,v)=0d0
        enddo
        bvec1(u)=0d0
        bvec2(u)=0d0
      enddo

      ic=0
      do k=1,nc
        if (valid(k).eq.0) cycle
        ic=ic+1
        do i=1,n(k)
          call ra_dec_to_xi_eta(ra(k,i),dec(k,i),xi(k,i),eta(k,i)
     .,cRVAL(1),cRVAL(2))

          do j=1,ntot
            vec(j)=0d0
          enddo


c          r=sqrt(xi(k,i)**2+eta(k,i)**2)
          px=0
          py=1
          order=1
          nn=0

          do while (nn.lt.npd)
            if (py.eq.order) then
c              if (mod(order,2).eq.1) then
c                nn=nn+1
c                vec(nn)=r**order
c                if (nn.eq.npd) cycle
c              endif
              order=order+1
              px=order
              py=0
            else
              px=px-1
              py=py+1
            endif
            nn=nn+1
            vec(nn)=xi(k,i)**px*eta(k,i)**py
          enddo



          j=npd+(ic-1)*3
          vec(j+1)=x(k,i)
          vec(j+2)=y(k,i)
          vec(j+3)=1d0

          j=npd+ic*3

          do u=1,j
            do v=1,j
              matx(u,v)=matx(u,v)+vec(u)*vec(v)
            enddo
            bvec1(u)=bvec1(u)+xi(k,i)*vec(u)
            bvec2(u)=bvec2(u)+eta(k,i)*vec(u)
          enddo

        enddo
      enddo

      k=npd+ic*3

c numerical_fix F4: the legacy path called matrix_inverse_doub, which
c (a) "normalized" the matrix by dividing every entry by one global
c scalar max|matx(i,j)| -- that leaves the condition number completely
c unchanged and does nothing about the column imbalance -- and then
c (b) built an EXPLICIT inverse and multiplied it into the right hand
c side. The design columns span 17 orders of magnitude (the PU terms
c xi**px*eta**py fall to ~1e-14 at 7th order while the per-chip affine
c columns carry pixel coordinates ~4e3), giving cond(A)=3.1e9 and
c cond(A^T A)=8.9e20, far past the double precision limit of 4.5e15.
c The published chip1 solution therefore missed its own least squares
c optimum by 0.089 arcsec: with a per-chip constant column present the
c optimal residual mean is identically zero, yet the released solution
c had deta=+0.0851 arcsec (19 sigma).
c
c Two changes, both exact in exact arithmetic:
c   1. symmetric Jacobi scaling matx -> D*matx*D, D=diag(1/sqrt(N_ii)),
c      which equilibrates the columns that the scalar norm could not;
c   2. one LU factorization reused for both right hand sides through
c      backsubstitution, instead of an explicit inverse.
c An error decomposition on the real data showed step 2 alone is what
c matters: the same NR LU driven by backsubstitution lands within
c 1e-5 arcsec, while LAPACK applied as inv(N)*b still produces the full
c 0.094 arcsec error. Step 1 additionally guarantees aamax>=1 inside
c ludcmp_doub, so its blocking read(*,*) branch cannot be reached.

      do i=1,k
        if (matx(i,i).gt.0d0) then
          dscal(i)=1d0/sqrt(matx(i,i))
        else
          dscal(i)=1d0
        endif
      enddo

      do i=1,k
        do j=1,k
          matx(i,j)=matx(i,j)*dscal(i)*dscal(j)
        enddo
        vec1(i)=bvec1(i)*dscal(i)
        vec2(i)=bvec2(i)*dscal(i)
      enddo

      call ludcmp_doub(matx,k,ntot,indx,dsign)
      call lubksb_doub(matx,k,ntot,indx,vec1)
      call lubksb_doub(matx,k,ntot,indx,vec2)

      do i=1,k
        vec1(i)=vec1(i)*dscal(i)
        vec2(i)=vec2(i)*dscal(i)
      enddo


      px=0
      py=1
      order=1
      nn=0
      do while (nn.lt.npd)
        if (py.eq.order) then
c          if (mod(order,2).eq.1) then
c            nn=nn+1
c            PU(1,nn)=vec1(nn)
c            PU(2,nn)=vec2(nn)
c            if (nn.eq.npd) cycle
c          endif
          order=order+1
          px=order
          py=0
        else
          px=px-1
          py=py+1
        endif
        nn=nn+1
        PU(1,nn)=vec1(nn)
        PU(2,nn)=vec2(nn+order-py*2)
      enddo


      ic=0
      do k=1,nc
        if (valid(k).eq.0) cycle
        ic=ic+1
        j=npd+(ic-1)*3
        cD(k,1,1)=vec1(j+1)
        cD(k,1,2)=vec1(j+2)

        cD(k,2,1)=vec2(j+1)
        cD(k,2,2)=vec2(j+2)

c        cD(k,1,1)*cRPIX(k,1)+cD(k,1,2)*cRPIX(k,2)=-vec1(j+3)
c        cD(k,2,1)*cRPIX(k,1)+cD(k,2,2)*cRPIX(k,2)=-vec2(j+3)

        cRPIX(k,1)=-(vec1(j+3)*cD(k,2,2)-vec2(j+3)*cD(k,1,2))
     ./(cD(k,1,1)*cD(k,2,2)-cD(k,2,1)*cD(k,1,2))
        cRPIX(k,2)=-(vec2(j+3)*cD(k,1,1)-vec1(j+3)*cD(k,2,1))
     ./(cD(k,1,1)*cD(k,2,2)-cD(k,2,1)*cD(k,1,2))

c        write(*,*) k,cD(k,1,1),cD(k,1,2),cD(k,2,1),cD(k,2,2)
c     .,cRPIX(k,1),cRPIX(k,2)
      enddo
c      pause

c r=sqrt(xi(i)^2+eta(i)^2)

c chi^2=\sum_i{cD11(k)*x(i)+cD12(k)*y(i)+S1(k)-xi(i)+PU(1,1)*r
c +PU(1,2)*xi(i)^2+PU(1,3)*xi(i)*eta(i)+PU(1,4)*eta(i)^2
c +PU(1,5)*xi(i)^3+PU(1,6)*xi(i)^2*eta(i)+PU(1,7)*xi(i)*eta(i)^2
c +PU(1,8)*eta(i)^3+PU(1,9)*r^3+... ...}^2

c      +\sum_i{cD21(k)*x(i)+cD22(k)*y(i)+S2(k)-eta(i)+PU(2,1)*r
c +PU(2,2)*eta(i)^2+PU(2,3)*xi(i)*eta(i)+PU(2,4)*xi(i)^2
c +PU(2,5)*eta(i)^3+PU(2,6)*eta(i)^2*xi(i)+PU(2,7)*eta(i)*xi(i)^2
c +PU(2,8)*xi(i)^3+PU(2,9)*r^3+... ...}^2



      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ==========================================
c Function: Synced astrometry routine mapping_PU
c Method: Preserve the current full-f77 numerical implementation.
c ==========================================
      subroutine mapping_PU(xx,yy,xi,eta,npd,PU,direc)
      implicit none

      integer npd,direc
      double precision xi,eta,xx,yy,PU(2,npd),dxi,deta,r

      integer n,px,py,order,i


       ! 计算场畸变后ccd平面 - 天球坐标的变换
      if (direc.eq.1) then
        xi=xx
        eta=yy
c        r=sqrt(xi**2+eta**2)

        do i=1,3
          dxi=0d0
          deta=0d0

          px=0
          py=1
          order=1
          n=0
          do while (n.lt.npd)
            if (py.eq.order) then
c              if (mod(order,2).eq.1) then
c                n=n+1
c                dxi=dxi+PU(1,n)*r**order
c                deta=deta+PU(2,n)*r**order
c                if (n.eq.npd) cycle
c              endif
              order=order+1
              px=order
              py=0
            else
              px=px-1
              py=py+1
            endif
            n=n+1
            dxi=dxi+PU(1,n)*xi**px*eta**py
            deta=deta+PU(2,n)*eta**px*xi**py
          enddo

          xi=xx+dxi
          eta=yy+deta
c          r=sqrt(xi**2+eta**2)

        enddo

      else

        xx=xi
        yy=eta
c        r=sqrt(xi**2+eta**2)

        px=0
        py=1
        order=1
        n=0
        do while (n.lt.npd)
          if (py.eq.order) then
c            if (mod(order,2).eq.1) then
c              n=n+1
c              xx=xx-PU(1,n)*r**order
c              yy=yy-PU(2,n)*r**order
c              if (n.eq.npd) cycle
c            endif
            order=order+1
            px=order
            py=0
          else
            px=px-1
            py=py+1
          endif
          n=n+1
          xx=xx-PU(1,n)*xi**px*eta**py
          yy=yy-PU(2,n)*eta**px*xi**py
        enddo

      endif

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ==========================================
c Function: Synced astrometry routine ra_dec_to_xi_eta
c Method: Preserve the current full-f77 numerical implementation.
c ==========================================
      subroutine ra_dec_to_xi_eta(ra,dec,xi,eta,cRVAL1,cRVAL2)
      implicit none

      double precision ra,dec,xi,eta,cRVAL1,cRVAL2,diffra

      double precision da,dd,const1,cosda,tandc,tandd,tanda,x,y
      double precision pi
      parameter (pi=3.1415926d0)

      const1=pi/180d0
      tandc=tan(cRVAL2*const1)

      da=diffra(ra,cRVAL1)*const1
      dd=(dec-cRVAL2)*const1
      tandd=tan(dd)
      cosda=cos(da)
      tanda=tan(da)

      y=tandc*(cosda-1.)-(1.+cosda*tandc**2)*tandd
      y=y/(tandd*tandc*(cosda-1.)-(cosda+tandc**2))
      x=tanda*(cos(cRVAL2*const1)*(1.-y*tandc))

      xi=x/const1
      eta=y/const1

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ==========================================
c Function: Synced astrometry routine diffra
c Method: Preserve the current full-f77 numerical implementation.
c ==========================================
      function diffra(ra1,ra2)
      implicit none

      double precision ra1,ra2,diffra

      diffra=ra1-ra2
      if (diffra.lt.-180d0) then
        diffra=diffra+360d0
      elseif (diffra.gt.180d0) then
        diffra=diffra-360d0
      endif

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ==========================================
c Function: Synced astrometry routine sumra
c Method: Preserve the current full-f77 numerical implementation.
c ==========================================
      function sumra(dra,ra)
      implicit none

      double precision dra,ra,sumra

      sumra=ra+dra
      if (sumra.ge.360d0) then
        sumra=sumra-360d0
      elseif (sumra.lt.0d0) then
        sumra=sumra+360d0
      endif

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ==========================================
c Function: Synced astrometry routine xy_to_xxyy
c Method: Preserve the current full-f77 numerical implementation.
c ==========================================
      subroutine xy_to_xxyy(x,y,xx,yy,cRPIX,cD)
      implicit none

      double precision x,y,xx,yy
      double precision cRPIX(2),cD(2,2)

      xx=cD(1,1)*(x-cRPIX(1))+cD(1,2)*(y-cRPIX(2))
      yy=cD(2,1)*(x-cRPIX(1))+cD(2,2)*(y-cRPIX(2))

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ==========================================
c Function: Synced astrometry routine coordinate_transfer_PU
c Method: Preserve the current full-f77 numerical implementation.
c ==========================================
      subroutine coordinate_transfer_PU(a,d,x,y,direc
     .,cRPIX,cD,cRVAL,PU,npd)
      implicit none

      integer direc,npd
      double precision a,d,x,y,xx,yy,xxx,yyy,rr,dxx,dyy
      double precision cRPIX(2),cD(2,2),cD_1(2,2),cRVAL(2)
      double precision PU(2,npd)

      double precision da,dd,const1,const2,ds,xi,eta,diffra
      double precision cosda,tandc,tandd,tanda,temp,sumra
      double precision pi
      parameter (pi=3.1415926d0)

c      integer tmp_sig
c      common /temp_pass/ tmp_sig

      const1=pi/180d0
      tandc=tan(cRVAL(2)*const1)

      if (direc.eq.1) then                            ! direc 是什么？
        xx=cD(1,1)*(x-cRPIX(1))+cD(1,2)*(y-cRPIX(2))  ! 转换到天球坐标
        yy=cD(2,1)*(x-cRPIX(1))+cD(2,2)*(y-cRPIX(2))

        call mapping_PU(xx,yy,xi,eta,npd,PU,1)

        xxx=xi*const1
        yyy=eta*const1

        da=xxx/(cos(cRVAL(2)*const1)*(1.-yyy*tandc))
        da=da-da*da*da*0.33333333333
        a=sumra(da/const1,cRVAL(1))
        cosda=1.-da*da*0.5+da**4/24.

        dd=(yyy*(cosda+tandc**2)+tandc*(cosda-1.))
     ./(yyy*tandc*(cosda-1.)+1.+cosda*tandc**2)
        dd=dd-dd*dd*dd*0.3333333333
        d=dd/const1+cRVAL(2)

      else

        da=diffra(a,cRVAL(1))*const1
        dd=(d-cRVAL(2))*const1
        tandd=tan(dd)
        cosda=cos(da)
        tanda=tan(da)

        yy=tandc*(cosda-1.)-(1.+cosda*tandc**2)*tandd
        yy=yy/(tandd*tandc*(cosda-1.)-(cosda+tandc**2))
        xx=tanda*(cos(cRVAL(2)*const1)*(1.-yy*tandc))

        xi=xx/const1
        eta=yy/const1

c        if (tmp_sig.eq.1.and.direc.eq.-1) then
c          write(*,*) da/const1,dd/const1,xi,eta
c        endif

        call mapping_PU(xx,yy,xi,eta,npd,PU,2)

        temp=cD(1,1)*cD(2,2)-cD(1,2)*cD(2,1)
        temp=1./temp

        cD_1(1,1)=cD(2,2)*temp
        cD_1(2,2)=cD(1,1)*temp
        cD_1(1,2)=-cD(1,2)*temp
        cD_1(2,1)=-cD(2,1)*temp

        x=xx*cD_1(1,1)+yy*cD_1(1,2)+cRPIX(1)
        y=xx*cD_1(2,1)+yy*cD_1(2,2)+cRPIX(2)

      endif

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ==========================================
c Function: Synced astrometry routine field_distortion_PU
c Method: Preserve the current full-f77 numerical implementation.
c ==========================================
      subroutine field_distortion_PU(x,y,npd
     .,PU,cD,cRPIX,g1,g2,cos2,sin2,parity)
      implicit none

      integer npd
      double precision xi,eta,x,y,xx,yy,g1,g2,cos2,sin2,aa,bb
      double precision cRPIX(2),cD(2,2),PU(2,npd),cos1,sin1
      integer parity

      double precision mat(2,2),dm(2,2),temp,det,sqrt_det
      double precision dx_dxx,dx_dyy,dy_dxx,dy_dyy
      double precision dxx_dxi,dxx_deta,dyy_dxi,dyy_deta
      double precision r,dr_dxi,dr_deta
      integer px,py,n,order

      ! error will occur if both x==cRPIX(1) and y==cRPIX(2) !
      xx=cD(1,1)*(x-cRPIX(1))+cD(1,2)*(y-cRPIX(2))
      yy=cD(2,1)*(x-cRPIX(1))+cD(2,2)*(y-cRPIX(2))

      call mapping_PU(xx,yy,xi,eta,npd,PU,1)
 
      temp=1d0/(cD(1,1)*cD(2,2)-cD(1,2)*cD(2,1))

      dx_dxx=cD(2,2)*temp
      dx_dyy=-cD(1,2)*temp
      dy_dxx=-cD(2,1)*temp
      dy_dyy=cD(1,1)*temp

      dxx_dxi=1d0
      dyy_deta=1d0
      dxx_deta=0d0
      dyy_dxi=0d0

c      r=sqrt(xi**2+eta**2)
c      if (r.gt.0d0) then
c        dr_dxi=xi/r
c        dr_deta=eta/r
c      else
c        dr_dxi=0d0
c        dr_deta=0d0
c      endif

      px=0
      py=1
      order=1
      n=0
      do while (n.lt.npd)
        if (py.eq.order) then
c          if (mod(order,2).eq.1) then
c            n=n+1
c            dxx_dxi=dxx_dxi-PU(1,n)*order*r**(order-1)*dr_dxi
c            dxx_deta=dxx_deta-PU(1,n)*order*r**(order-1)*dr_deta
c            dyy_deta=dyy_deta-PU(2,n)*order*r**(order-1)*dr_deta
c            dyy_dxi=dyy_dxi-PU(2,n)*order*r**(order-1)*dr_dxi
c            if (n.eq.npd) cycle
c          endif
          order=order+1
          px=order
          py=0
        else
          px=px-1
          py=py+1
        endif
        n=n+1
        dxx_dxi=dxx_dxi-PU(1,n)*px*xi**(px-1)*eta**py
        dxx_deta=dxx_deta-PU(1,n)*xi**px*py*eta**(py-1)
        dyy_deta=dyy_deta-PU(2,n)*px*eta**(px-1)*xi**py
        dyy_dxi=dyy_dxi-PU(2,n)*eta**px*py*xi**(py-1)

      enddo


      mat(1,1)=dx_dxx*dxx_dxi+dx_dyy*dyy_dxi
      mat(1,2)=dx_dxx*dxx_deta+dx_dyy*dyy_deta
      mat(2,1)=dy_dxx*dxx_dxi+dy_dyy*dyy_dxi
      mat(2,2)=dy_dxx*dxx_deta+dy_dyy*dyy_deta

      det=mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)
      temp=1d0/det

      dm(1,1)=mat(2,2)*temp
      dm(1,2)=-mat(1,2)*temp
      dm(2,1)=-mat(2,1)*temp
      dm(2,2)=mat(1,1)*temp

      parity=1
      if (det.lt.0) then
        dm(1,1)=-dm(1,1)
        dm(1,2)=-dm(1,2)
        parity=-1
      endif

      sqrt_det=sqrt(dm(1,1)*dm(2,2)-dm(1,2)*dm(2,1))

      dm(1,1)=dm(1,1)/sqrt_det
      dm(1,2)=dm(1,2)/sqrt_det
      dm(2,1)=dm(2,1)/sqrt_det
      dm(2,2)=dm(2,2)/sqrt_det

      cos1=0.5d0*(dm(1,1)+dm(2,2))
      sin1=0.5d0*(dm(1,2)-dm(2,1))

      aa=-0.5d0*(dm(1,2)+dm(2,1))
      bb=0.5d0*(dm(2,2)-dm(1,1))

      g1=aa*sin1+bb*cos1
      g2=aa*cos1-bb*sin1

      cos2=cos1*cos1-sin1*sin1
      sin2=2d0*sin1*cos1

      if (parity.eq.-1) g2=-g2

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ==========================================
c Function: Synced astrometry routine get_ra_dec_bound
c Method: Preserve the current full-f77 numerical implementation.
c ==========================================
      subroutine get_ra_dec_bound(np,n,a,d,ra,dra,dec,ddec)
      implicit none

      integer np,n
      double precision ra,dec,dra,ddec
      double precision a(np),d(np),ra1,ra2,dec1,dec2
      double precision tmp(n)
      integer i

      dec1=d(1)
      dec2=d(1)

      do i=1,n
        dec1=min(dec1,d(i))
        dec2=max(dec2,d(i))
      enddo
      dec=(dec1+dec2)*0.5d0
      ddec=dec2-dec1

      ra1=a(1)
      ra2=a(1)
      do i=1,n
        ra1=min(ra1,a(i))
        ra2=max(ra2,a(i))
      enddo

      ra=0.5*(ra1+ra2)
      dra=ra2-ra1

      if (dra.gt.180d0) then

        do i=1,n
          if (a(i).lt.180d0) then
            tmp(i)=a(i)+360d0
          else
            tmp(i)=a(i)
          endif
        enddo

        ra1=tmp(1)
        ra2=tmp(1)
        do i=1,n
          ra1=min(ra1,tmp(i))
          ra2=max(ra2,tmp(i))
        enddo

        dra=ra2-ra1
        ra=0.5*(ra1+ra2)
        if (ra.ge.360d0) ra=ra-360d0

      endif


      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ==========================================
c Function: Synced astrometry routine get_ra_dec_range_fine
c Method: Preserve the current full-f77 numerical implementation.
c ==========================================
      subroutine get_ra_dec_range_fine(nx,ny,ra,dec,dra
     .,cRPIX,cD,cRVAL,PU,npd,astrometry_shift_ratio)
      implicit none

      integer nx,ny,npd
      double precision ra,dec(2),dra,dec_c,ddec
      double precision x,y,a(4),d(4)

      double precision cRPIX(2),cD(2,2),cRVAL(2)
      double precision PU(2,npd)

      real astrometry_shift_ratio

      x=1d0
      y=1d0
      call coordinate_transfer_PU(a(1),d(1),x,y,1,cRPIX,cD,cRVAL,PU,npd)
      x=nx
      y=1d0
      call coordinate_transfer_PU(a(2),d(2),x,y,1,cRPIX,cD,cRVAL,PU,npd)
      x=nx
      y=ny
      call coordinate_transfer_PU(a(3),d(3),x,y,1,cRPIX,cD,cRVAL,PU,npd)
      x=1d0
      y=ny
      call coordinate_transfer_PU(a(4),d(4),x,y,1,cRPIX,cD,cRVAL,PU,npd)

      call get_ra_dec_bound(4,4,a,d,ra,dra,dec_c,ddec)

      dec(1)=dec_c-0.5*ddec
      dec(2)=dec_c+0.5*ddec

      ddec=ddec*astrometry_shift_ratio

      dec(1)=dec(1)-ddec
      dec(2)=dec(2)+ddec

      dra=dra*(1d0+2d0*astrometry_shift_ratio)


      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ==========================================
c Function: Synced astrometry routine get_ra_dec_range
c Method: Preserve the current full-f77 numerical implementation.
c ==========================================
      subroutine get_ra_dec_range(nx,ny,ra,dec,dra
     .,cRPIX,cD,cRVAL,astrometry_shift_ratio)
      implicit none

      integer nx,ny
      double precision ra,dec(2),dra,dec_c,ddec
      double precision x,y,a(4),d(4)

      double precision cRPIX(2),cD(2,2),cRVAL(2)
      real astrometry_shift_ratio

      x=1d0
      y=1d0
      call coordinate_transfer_simple(a(1),d(1),x,y,1,cRPIX,cD,cRVAL)
      x=nx
      y=1d0
      call coordinate_transfer_simple(a(2),d(2),x,y,1,cRPIX,cD,cRVAL)
      x=nx
      y=ny
      call coordinate_transfer_simple(a(3),d(3),x,y,1,cRPIX,cD,cRVAL)
      x=1d0
      y=ny
      call coordinate_transfer_simple(a(4),d(4),x,y,1,cRPIX,cD,cRVAL)

      call get_ra_dec_bound(4,4,a,d,ra,dra,dec_c,ddec)

      dec(1)=dec_c-0.5*ddec
      dec(2)=dec_c+0.5*ddec

      ddec=ddec*astrometry_shift_ratio

      dec(1)=dec(1)-ddec
      dec(2)=dec(2)+ddec

      dra=dra*(1d0+2d0*astrometry_shift_ratio)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ==========================================
c Function: Synced astrometry routine coordinate_transfer_simple
c Method: Preserve the current full-f77 numerical implementation.
c ==========================================
      subroutine coordinate_transfer_simple(a,d,x,y,direc
     .,cRPIX,cD,cRVAL)
      implicit none

      integer direc
      double precision a,d,x,y,xx,yy,xxx,yyy
      double precision cRPIX(2),cD(2,2),cD_1(2,2),cRVAL(2)

      double precision da,dd,const1,const2,ds,diffra,sumra
      double precision cosda,tandc,tandd,tanda,temp
      double precision pi
      parameter (pi=3.1415926d0)

      const1=pi/180d0
      tandc=tan(cRVAL(2)*const1)

      if (direc.eq.1) then
        xx=cD(1,1)*(x-cRPIX(1))+cD(1,2)*(y-cRPIX(2))
        yy=cD(2,1)*(x-cRPIX(1))+cD(2,2)*(y-cRPIX(2))

        xxx=xx*const1
        yyy=yy*const1

        da=xxx/(cos(cRVAL(2)*const1)*(1.-yyy*tandc))
        da=da-da*da*da*0.33333333333
        a=sumra(da/const1,cRVAL(1))
        cosda=1.-da*da*0.5+da**4/24.

        dd=(yyy*(cosda+tandc**2)+tandc*(cosda-1.))
     ./(yyy*tandc*(cosda-1.)+1.+cosda*tandc**2)
        dd=dd-dd*dd*dd*0.3333333333
        d=dd/const1+cRVAL(2)

      else

        da=diffra(a,cRVAL(1))*const1
        dd=(d-cRVAL(2))*const1
        tandd=tan(dd)
        cosda=cos(da)
        tanda=tan(da)

        yy=tandc*(cosda-1.)-(1.+cosda*tandc**2)*tandd
        yy=yy/(tandd*tandc*(cosda-1.)-(cosda+tandc**2))
        xx=tanda*(cos(cRVAL(2)*const1)*(1.-yy*tandc))

        xx=xx/const1
        yy=yy/const1

        temp=cD(1,1)*cD(2,2)-cD(1,2)*cD(2,1)
        temp=1./temp

        cD_1(1,1)=cD(2,2)*temp
        cD_1(2,2)=cD(1,1)*temp
        cD_1(1,2)=-cD(1,2)*temp
        cD_1(2,1)=-cD(2,1)*temp

        x=xx*cD_1(1,1)+yy*cD_1(1,2)+cRPIX(1)
        y=xx*cD_1(2,1)+yy*cD_1(2,2)+cRPIX(2)

      endif

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

