      subroutine chip_psf_recons(nexpo)
      use mpi
      use psf_storage_mod, only: init_and_load_all_psf
      implicit none
      include 'para.inc'
      include 'cust_para.inc'

      integer my_id,num_procs
      common /MPIpar/ my_id,num_procs

      integer nexpo
      integer ierr
      character*(strl) DIR_OUTPUT, IMA1(NMAX_cHIP)
      integer nc1

      external chip_res_pca_fit,Plot_residuals_v2

      call get_image_list(1,IMA1,nc1,DIR_OUTPUT)

      call mpi_forcov(procs_pn,work_pn,Camera_ccd_num,chip_res_pca_fit
     .,'fitting residual...',nexpo)
      if (my_id.eq.0) then
        write(*,*) 'PSF PCA fitting completed for all chips.'
      endif

      call MPI_BARRIER(MPI_cOMM_WORLD, ierr ) 
      call init_and_load_all_psf(DIR_OUTPUT,my_id)
      call MPI_BARRIER(MPI_cOMM_WORLD, ierr ) 

      call mpi_distribute(nexpo,Plot_residuals_v2
     .,'Mapping Modified Residuals...')
      call MPI_BARRIER(MPI_cOMM_WORLD, ierr ) 
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine chip_res_pca_fit(ichip,nexpo)
      implicit none
      include 'para.inc'
      include 'cust_para.inc'

      integer u,v,ierror,ntot,i,nc
      integer nn1,nn2

      integer ichip,nexpo
      integer nstar,valid
      integer j,k,info,ibstar
      integer fit_num

      character*(strl) DIR_OUTPUT, IMAGE_FILE(NMAX_cHIP)
      character*(strl) PREFIX,PREFIX_e,filename
      character*(2) c_chip
      character*(1) c_bx, c_by

      real px,py
      integer valid_num
      real xc(nmax_star_pchip),yc(nmax_star_pchip)
      real,allocatable :: psf_residual(:,:,:)
      real res_slice(nsns)
      real coe_fit(npp6th)
      real lxd,uxd,lyd,uyd,nxc,nyc

      integer, parameter :: block_size = 2000
      integer buf_cnt
      double precision term_corr
      double precision, allocatable :: cov_arr(:,:)
      double precision, allocatable :: mean_arr(:)
      double precision, allocatable :: block_dble(:,:)

      double precision,allocatable :: eigval(:)
      double precision,allocatable :: eigvec(:,:)
      double precision, allocatable :: components(:,:)

      integer curr_star_idx
      real,allocatable :: xfit(:),yfit(:),zfit(:)
      integer, allocatable :: idx_list(:)
      double precision,allocatable :: temp(:),coeff(:,:)

100   format(I2.2)
200   format(I1.1)

      if (ichip.eq.2 .or. ichip.eq.61) return

      allocate(block_dble(block_size, nsns))
      allocate(cov_arr(nsns,nsns))
      allocate(mean_arr(nsns))
      !--------------------- Read star residuals ------------------------!
      ntot=0
      mean_arr = 0.0d0
      cov_arr = 0.0d0
      buf_cnt = 0
      do i=1,nexpo

        call get_image_list(i,IMAGE_FILE,nc,DIR_OUTPUT)
        call get_PREFIX_expo(IMAGE_FILE(1),PREFIX_e)
        write(c_chip,'(I2)') ichip
        filename=trim(DIR_OUTPUT)//'/starxy/'//trim(PREFIX_e)//
     .'_'//trim(adjustl(c_chip))//'_star_xy.dat'
        open(unit=10,file=filename,status='old',action='read'
     ,,iostat=ierror)
        if (ierror.ne.0) then
            cycle
        endif
        rewind 10
        read(10,*) nstar, valid
        if (nstar .le. 0 .or. valid .lt. 0) then
          close(10)
          cycle
        endif

        allocate(psf_residual(nstar,ns,ns))

        nn1=ns*len_s
        nn2=ns*(int(nstar/len_s)+1)
        filename=trim(DIR_OUTPUT)//'/fits_psfresi/'//trim(PREFIX_e)
     .//'_'//trim(adjustl(c_chip))//'_psf_p_resi.fits' 
        call read_stamps(nstar,1,nstar,ns,ns,psf_residual
     .,nn1,nn2,filename)

        valid_num = 0
        do while (ierror.ge.0)
          read(10,*,iostat=ierror) px,py
          if (ierror.lt.0) cycle
          valid_num = valid_num + 1
          if (px .lt. 0.0 .or. py .lt. 0.0) cycle
          ntot = ntot + 1
          buf_cnt = buf_cnt + 1

          do u=1,ns
            do v=1,ns
        block_dble(buf_cnt,(u-1)*ns+v) = psf_residual(valid_num,u,v)
            enddo
          enddo

          if (buf_cnt .eq. block_size) then
            call accumulate_block(block_size, nsns, block_dble,
     .                               mean_arr, cov_arr)
             if (mean_arr(32*ns+33).lt.1d-30) 
     .write(*,*) 'mean sums error'
               buf_cnt = 0
            endif
        enddo
        close(10)
        deallocate(psf_residual)
      enddo

      if (buf_cnt .gt. 0) then
        call accumulate_block(buf_cnt, nsns, block_dble(1:buf_cnt,:),
     .                         mean_arr, cov_arr)
      endif

      write(*,*) 'chip ',ichip,' total stars: ',ntot
      if (allocated(block_dble)) deallocate(block_dble)

      !--------------------- PCA decomponent ----------------------------!
      if (ntot.eq.0) then
        write(*,*) 'No valid stars found for this ccd ',ichip
        if (allocated(cov_arr)) deallocate(cov_arr)
        if (allocated(mean_arr)) deallocate(mean_arr)
        return
      else
        mean_arr = mean_arr / dble(ntot)
        do j=1,nsns
          do k=1,nsns
            term_corr = dble(ntot) * mean_arr(j) * mean_arr(k)
            cov_arr(j,k) = (cov_arr(j,k) - term_corr) / dble(ntot - 1)
          enddo
        enddo
      endif

      allocate(eigval(nsns))
      allocate(eigvec(nsns, nsns))
      allocate(components(nsns, n_pcs))

      call get_eigval_vec(nsns,cov_arr,eigval,eigvec,info)
      if (info.ne.0) then
        components(1,1) = -1.0d30
        write(*,*) 'Error: PCA failed for chip :',ichip
      else                             
        do j = 1, n_pcs
          do k = 1, nsns
          components(k,j) = eigvec(k, nsns+1-j)
          enddo
        enddo
      endif
      
      call get_image_list(1,IMAGE_FILE,nc,DIR_OUTPUT)
      write(c_chip,100) ichip
      filename = trim(DIR_OUTPUT)//'/dat_pcs/'//'pcs_ccd'//trim(c_chip)
     .//'.dat'
      open(unit=30,file=filename,status='replace',iostat=ierror)
      rewind 30
      do k = 1,nsns
        write(30,*) (components(k,j) ,j=1,n_pcs), dble(mean_arr(k))
      enddo
      close(30)
      write(*,*) 'PCA finished chip ...',ichip

      deallocate(eigval,eigvec)
      deallocate(cov_arr)

      !---------- Fit PCA coefficients as functions of position -----!
      allocate(temp(n_pcs))
      allocate(coeff(ntot,n_pcs))

      ntot = 0
      do i=1,nexpo
        call get_image_list(i,IMAGE_FILE,nc,DIR_OUTPUT)
        call get_PREFIX_expo(IMAGE_FILE(1),PREFIX_e)

        write(c_chip,'(I2)') ichip
        filename=trim(DIR_OUTPUT)//'/starxy/'//trim(PREFIX_e)//
     .            '_'//trim(adjustl(c_chip))//'_star_xy.dat'
        open(unit=10,file=filename,status='old',action='read'
     .,iostat=ierror)
        if (ierror.ne.0) then
          cycle
        endif
        rewind 10
        read(10,*) nstar, valid
        if (nstar .le. 0 .or. valid .lt. 0) then
          close(10)
          cycle
        endif
         
        allocate(psf_residual(nstar,ns,ns))
        nn1=ns*len_s
        nn2=ns*(int(nstar/len_s)+1)
        filename=trim(DIR_OUTPUT)//'/fits_psfresi/'//trim(PREFIX_e)
     .            //'_'//trim(adjustl(c_chip))//'_psf_p_resi.fits'
         
        call read_stamps(nstar,1,nstar,ns,ns,psf_residual,
     .                    nn1,nn2,filename)

        valid_num = 0
        do while (ierror.ge.0)
          read(10,*,iostat=ierror) px,py
          if (ierror.lt.0) cycle
          valid_num = valid_num + 1
          if (px .lt. 0.0 .or. py .lt. 0.0) cycle
          ntot = ntot + 1

          xc(ntot) = px
          yc(ntot) = py
          
          do u=1,ns
            do v=1,ns
              res_slice((u-1)*ns+v) = psf_residual(valid_num,u,v) -
     .                                    real(mean_arr((u-1)*ns+v))
            enddo
          enddo
            
          temp(1:n_pcs) = matmul(dble(res_slice), components(:,1:n_pcs))
          coeff(ntot, :) = temp(:)
        enddo
        close(10)
        deallocate(psf_residual)
      enddo

      if (allocated(mean_arr)) deallocate(mean_arr)
      allocate(xfit(ntot))
      allocate(yfit(ntot))
      allocate(zfit(ntot))
      allocate(idx_list(ntot))

      call get_image_list(1,IMAGE_FILE,nc,DIR_OUTPUT)
      write(c_chip,100) ichip
      nxc = real(chipnx) / nblocks
      nyc = real(chipny) / nblocks
      do j=1,nblocks
        lxd = (j-1.0)*nxc
        uxd = j*nxc
        do k=1,nblocks
          lyd = (k-1.0)*nyc
          uyd = k*nyc
          fit_num = 0
          if (components(1,1) .lt. -1.0d20) then
            write(*,*) 'PCA decomposition failed for chip ',ichip
            goto 99
          endif
          do ibstar=1,ntot
            if (xc(ibstar).lt.lxd.or.xc(ibstar).gt.uxd) cycle
            if (yc(ibstar).lt.lyd.or.yc(ibstar).gt.uyd) cycle
            fit_num = fit_num + 1
            idx_list(fit_num) = ibstar
            xfit(fit_num) = 2.0*(xc(ibstar) - lxd) / (uxd - lxd) -1.0
            yfit(fit_num) = 2.0*(yc(ibstar) - lyd) / (uyd - lyd) -1.0
          enddo

99        write(c_bx,200) j
          write(c_by,200) k
      filename = trim(DIR_OUTPUT)//'/dat_pcs/'//'coeff_ccd'
     .//trim(c_chip)//'_'//trim(c_bx)//trim(c_by)//'.dat'
          open(unit=20,file=filename,status='replace',iostat=ierror)
          rewind 20
          if (fit_num .le. (npp6th+10)) then
            write(*,*) 'Not enough stars for fitting in chip ',ichip,j,k
            write(20,*) (-1.0e30,i=1,npp6th)
            close(20)
            cycle
          endif
          do u=1,n_pcs
            do i=1,fit_num
                curr_star_idx = idx_list(i)  
                zfit(i) = real(coeff(curr_star_idx, u))
             enddo
      call interpolate_6th(fit_num,xfit(1:fit_num),yfit(1:fit_num)
     .,zfit(1:fit_num),npp6th,coe_fit)
            if (isnan(coe_fit(1))) 
     .write(*,*) 'NaN in 6th fitting coef.',ichip,j,k
            write(20,*) (coe_fit(i),i=1,npp6th)
          enddo
          close(20)

        enddo
      enddo

      deallocate(xfit,yfit,zfit,temp,coeff,idx_list)
      deallocate(components)


      return
      end
c ------------------------------------------------------------------
      subroutine accumulate_block(n, m, blk, sum_x, sum_xxt)
      implicit none
      integer, intent(in) :: n, m
      double precision, intent(in) :: blk(n, m)
      double precision, intent(inout) :: sum_x(m)
      double precision, intent(inout) :: sum_xxt(m, m)
      integer j


      do j = 1, m
         sum_x(j) = sum_x(j) + sum(blk(1:n, j))
      enddo

      sum_xxt = sum_xxt + matmul(transpose(blk(1:n, :)), blk(1:n, :))

      return
      end subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_eigval_vec(d, cov, eigval, eigvec, info)
        implicit none
      integer d, lwork, liwork, info, query_info
      double precision cov(d,d), eigval(d), eigvec(d,d)
      double precision, allocatable :: work(:)
      integer, allocatable :: iwork(:)
      double precision work_query
      integer iwork_query

      eigvec(:,:) = cov(:,:)

      ! Query optimal work and iwork sizes
      lwork = -1
      liwork = -1
      call dsyevd('V', 'U', d, eigvec, d, eigval
     ., work_query, lwork, iwork_query, liwork, query_info)

      if (query_info .ne. 0) then
          write(*,*) 'LAPACK workspace query failed...'
          lwork = 26 * d
          liwork = 10 * d
      else
          lwork = int(work_query) + 1
          liwork = iwork_query + 1
      endif

      allocate(work(lwork))
      allocate(iwork(liwork))

      ! Actual diagonalization
      call dsyevd('V', 'U', d, eigvec, d, eigval
     ., work, lwork, iwork, liwork, info)

      if (info .ne. 0) then
          write(*,*) 'Eigen decomposition failed, info=', info
      endif

      deallocate(work)
      deallocate(iwork)
      return
      end subroutine get_eigval_vec
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine interpolate_6th(nsam,x,y,z,nc,coef)
      implicit none

c      fit a 6D function with f(x,y)=\Sigma_{i=1}^{nc}func_i(x,y)c_i

      integer nsam,nc
      real x(nsam),y(nsam),z(nsam),coef(nc)
      real coe(nc,nc),coe_1(nc,nc),vec(nc),temp(nc)
      real xi,yi,zi
      real fit_func_2
      external fit_func_2
      integer i,j,u,v

      coe=0.
      vec=0.

      do i=1,nsam
        xi=x(i)
        yi=y(i)
        zi=z(i)
        do j=1,nc
          temp(j)=fit_func_2(xi,yi,j)
        enddo
        do u=1,nc
          do v=1,nc
            coe(u,v)=coe(u,v)+temp(u)*temp(v)
          enddo
          vec(u)=vec(u)+zi*temp(u)
        enddo
      enddo

      call matrix_inverse(coe,nc,coe_1)

      do j=1,nc
        coef(j)=0.
        do i=1,nc
          coef(j)=coef(j)+coe_1(j,i)*vec(i)
        enddo
        if (isnan(coef(j))) then
          coef(1) = coef(j)
          return
        endif 
      enddo

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mpi_forcov(ppn,work_pn,num_job,job,message,nexpo)
      use mpi
      implicit none

      integer my_id,num_procs
      common /MPIpar/ my_id,num_procs

      external job

      character*(*) message

      integer i,j,k,complete,ppn,work_pn,num_job,nexpo
      integer source,tag,ierr,status(mpi_status_size)

      integer id_innode


      call MPI_BARRIER(MPI_cOMM_WORLD, ierr ) ! synchronize all nodes
      
      i=0
      j=0
      if (my_id.eq.0) then
        i=1
        j=num_job
      endif

      complete = 1
      id_innode = mod(my_id,ppn)
      if (id_innode .gt. 0 .and. id_innode .le. work_pn) complete = 0
      if (my_id .eq. 0) complete = 0

      do while (complete.eq.0)
        if (my_id.ne.0) then
          call MPI_SEND(i,1,mpi_int,0,0,MPI_cOMM_WORLD,ierr)
          call MPI_REcV(i,1,mpi_int,0,0,MPI_cOMM_WORLD,status,ierr)
          if (i.eq.0) then
            complete=1
          else
            call job(i,nexpo)
          endif
        else
          call MPI_REcV(k,1,mpi_int,MPI_ANY_SOURcE
     .,MPI_ANY_TAG,MPI_cOMM_WORLD,status,ierr)
          tag=status(MPI_TAG)
          source=status(MPI_SOURcE)
          call MPI_SEND(i,1,mpi_int,source,tag,MPI_cOMM_WORLD,ierr)
          if (k.gt.0) then
            j=j-1
          endif
          write(*,*) source,i,j,trim(message)
          if (i.ne.0) i=i+1
          if (i.gt.num_job) i=0
          if (j.eq.0) complete=1
        endif
      enddo

      call MPI_BARRIER(MPI_cOMM_WORLD, ierr )


      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_PSF_model_hierarchical(i_ccd,x,y,refactor
     .,local_coe,psf_model)
      use psf_storage_mod, only: global_components, 
     .                           global_mean_psf,
     .                           global_poly_coefs,
     .                           is_data_loaded
      implicit none
      include 'para.inc'
      include 'cust_para.inc'


      double precision x, y ,xx_norm, yy_norm     
      double precision local_coe(ns,ns,npl+1)
      integer i_ccd
      real psf_model(ns,ns)
      real refactor

      integer k, j, u, v

      integer bx, by
      real nxc, nyc        
      real lxd, uxd, lyd, uyd
      real x_norm, y_norm     

      real vec_b(npp6th)
      real coeff_val(n_pcs)
      real psf_layer1(ns,ns), psf_layer2(ns,ns)

      real fit_func_2         
      external fit_func_2


      if (.not. is_data_loaded) then
          write(*,*) 'Error: PSF data not loaded! Call init first.'
          stop
      endif

      xx_norm = 2.0*(x / dble(chipnx)) -1.0
      yy_norm = 2.0*(y / dble(chipny)) -1.0
      call get_PSF_model(ns,npl,nplx,local_coe,xx_norm,yy_norm
     .,psf_layer1,psf_layer2)
      if (isnan(psf_layer1(1,1))) then
        write(*,*) 'Error in getting PSF model layer1 for chip ',i_ccd
        psf_model = psf_layer1
        return
      endif

      bx = 0
      by = 0
      nxc = real(chipnx) / nblocks
      nyc = real(chipny) / nblocks
      do j=1,nblocks
        lxd = (j-1.0)*nxc
        uxd = j*nxc
        do k=1,nblocks
          lyd = (k-1.0)*nyc
          uyd = k*nyc
            if (real(x).lt.lxd.or.real(x).gt.uxd) cycle
            if (real(y).lt.lyd.or.real(y).gt.uyd) cycle
            bx = j
            by = k
            x_norm = 2.0*real(x - lxd) / (uxd - lxd) -1.0
            y_norm = 2.0*real(y - lyd) / (uyd - lyd) -1.0
            goto 100
        enddo
      enddo

      if (bx.eq.0.or.by.eq.0) then 
        write(*,*) 'Error: cannot find the block for source ',i_ccd
        psf_model = psf_layer1
        return 
      endif

100   if (global_components(i_ccd,1,1) .lt. -1.0d20) then
        !  write(*,*) 'PCA decomposition failed for chip ', i_ccd
         psf_model = psf_layer1
         return
      endif

      do j = 1, npp6th
        vec_b(j) = fit_func_2(x_norm, y_norm, j)
      enddo

      if (global_poly_coefs(i_ccd,bx,by,1,1) .lt. -1.0e20) then
          ! write(*,*) 'Errors in fitting in chip',i_ccd,bx,by
          psf_model = psf_layer1
          return
      endif

      do u = 1, n_pcs
        coeff_val(u) = 0.0
        do j = 1, npp6th
          coeff_val(u) = coeff_val(u) + 
     .       global_poly_coefs(i_ccd,bx,by,u,j) * vec_b(j)
        enddo
      enddo
      
      do u = 1,ns
        do v = 1,ns
          k = (u-1)*ns + v
          psf_layer2(u, v) = real(global_mean_psf(i_ccd,k))
          do j = 1, n_pcs
             psf_layer2(u, v) = psf_layer2(u, v) + 
     .            real(global_components(i_ccd, k, j)) * coeff_val(j)
          enddo
        enddo
      enddo

      call PSF_unscale(psf_layer2, refactor)

      psf_model = psf_layer1 + psf_layer2

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine PSF_rescale(resi,refactor)
      implicit none
      include 'para.inc'

      real resi(ns,ns),temp(ns,ns)
      integer i,j
      real refactor

      integer lx, ly
      real xc, yc, rx, ry
      real d1, d2, d3, d4
      real w1, w2, w3, w4, wtot
      real val1, val2, val3, val4
      
      
      xc = real(ns)/2.0 + 1.0
      yc = real(ns)/2.0 + 1.0
      do j=1,ns
        do i=1,ns
           temp(i,j) = 0.0
        enddo
      enddo

      do j=1,ns
        do i=1,ns
          rx = (real(i) - xc) * refactor + xc
          ry = (real(j) - yc) * refactor + yc
          lx = floor(rx)
          ly = floor(ry)    
          
          val1 = 0.0
          val2 = 0.0
          val3 = 0.0
          val4 = 0.0
          IF (lx .le. ns .and. ly .le. ns 
     .       .and. lx .ge. 1 .and. ly .ge. 1) THEN
              val1 = resi(lx, ly)
          ENDIF
          IF ((lx+1) .le. ns .and. ly .le. ns 
     .       .and. (lx+1) .ge. 1 .and. ly .ge. 1) THEN
              val2 = resi(lx+1, ly)
          ENDIF
          IF (lx .le. ns .and. (ly+1) .le. ns 
     .       .and. lx .ge. 1 .and. (ly+1) .ge. 1) THEN
              val3 = resi(lx, ly+1)
          ENDIF
          IF ((lx+1) .le. ns .and. (ly+1) .le. ns 
     .       .and. (lx+1) .ge. 1 .and. (ly+1) .ge. 1) THEN
              val4 = resi(lx+1, ly+1)
          ENDIF

    !       d1 = sqrt((rx-real(lx))**2 + (ry-real(ly))**2)
    !       d2 = sqrt((rx-real(lx+1))**2 + (ry-real(ly))**2)
    !       d3 = sqrt((rx-real(lx))**2 + (ry-real(ly+1))**2)
    !       d4 = sqrt((rx-real(lx+1))**2 + (ry-real(ly+1))**2)

    !       w1 = 1.0 / (d1)
    !       w2 = 1.0 / (d2)
    !       w3 = 1.0 / (d3)
    !       w4 = 1.0 / (d4)
    !       wtot = w1 + w2 + w3 + w4
    !       rescalemodel(i,j) = (val1*w1 + val2*w2 + val3*w3 + val4*w4)
    !  ./ wtot

          d1 = real(lx+1)-rx
          d2 = rx - real(lx)
          d3 = real(ly+1)-ry
          d4 = ry - real(ly)
          
          w1 = val1 * d1 * d3
          w2 = val2 * d2 * d3
          w3 = val3 * d1 * d4
          w4 = val4 * d2 * d4

          temp(i,j) = w1 + w2 + w3 + w4

        enddo
      enddo

      do i = 1,ns
        do j = 1,ns
          resi(i,j) = temp(i,j)
        enddo
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine PSF_unscale(model_in, scale_factor)
      implicit none
      include 'para.inc'

      real model_in(ns,ns)    
      real psf_out(ns,ns)       
      real scale_factor         
      
      integer i, j, lx, ly
      real xc, yc, rx, ry
      real val1, val2, val3, val4
      real d1, d2, d3, d4, wtot
      real w1, w2, w3, w4      
      real inv_scale


      inv_scale = 1.0 / scale_factor
      
      xc = real(ns)/2.0 + 1.0
      yc = real(ns)/2.0 + 1.0
      do j = 1, ns
        do i = 1, ns
            psf_out(i,j) = 0.0
        enddo
      enddo

      do j = 1, ns
        do i = 1, ns
          rx = (real(i) - xc) * inv_scale + xc
          ry = (real(j) - yc) * inv_scale + yc
          lx = floor(rx)
          ly = floor(ry)   

          val1 = 0.0
          val2 = 0.0
          val3 = 0.0
          val4 = 0.0
          IF (lx .le. ns .and. ly .le. ns 
     .       .and. lx .ge. 1 .and. ly .ge. 1) THEN
              val1 = model_in(lx, ly)
          ENDIF
          IF ((lx+1) .le. ns .and. ly .le. ns 
     .       .and. (lx+1) .ge. 1 .and. ly .ge. 1) THEN
              val2 = model_in(lx+1, ly)
          ENDIF
          IF (lx .le. ns .and. (ly+1) .le. ns 
     .       .and. lx .ge. 1 .and. (ly+1) .ge. 1) THEN
              val3 = model_in(lx, ly+1)
          ENDIF
          IF ((lx+1) .le. ns .and. (ly+1) .le. ns 
     .       .and. (lx+1) .ge. 1 .and. (ly+1) .ge. 1) THEN
              val4 = model_in(lx+1, ly+1)
          ENDIF

          ! d1 = sqrt((rx-real(lx))**2   + (ry-real(ly))**2)
          ! d2 = sqrt((rx-real(lx+1))**2 + (ry-real(ly))**2)
          ! d3 = sqrt((rx-real(lx))**2   + (ry-real(ly+1))**2)
          ! d4 = sqrt((rx-real(lx+1))**2 + (ry-real(ly+1))**2)

          ! w1 = 1.0 / (d1 + tiny_val)
          ! w2 = 1.0 / (d2 + tiny_val)
          ! w3 = 1.0 / (d3 + tiny_val)
          ! w4 = 1.0 / (d4 + tiny_val)

          ! wtot = w1 + w2 + w3 + w4
          
          ! psf_out(i,j) = (val1*w1 + val2*w2 + val3*w3 + val4*w4) / wtot

          d1 = real(lx+1)-rx
          d2 = rx - real(lx)
          d3 = real(ly+1)-ry
          d4 = ry - real(ly)
          
          w1 = val1 * d1 * d3
          w2 = val2 * d2 * d3
          w3 = val3 * d1 * d4
          w4 = val4 * d2 * d4

          psf_out(i,j) = w1 + w2 + w3 + w4

        enddo
      enddo

      do i = 1,ns
        do j = 1,ns
          model_in(i,j) = psf_out(i,j)
        enddo
      enddo
      return
      end  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Plot_residuals_v2(iexpo)
      implicit none
      include 'para.inc'
      include 'cust_para.inc'

      integer iexpo,nchip,ierror
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT
      character*(strl) PREFIX,filename,PREFIXc

      integer ichip, nstar, valid,chip_index,chip_circle,u,v
      integer undefined_flag
      real psf_para(8)
      double precision px,py

      integer nstar_coe,status,proc_error,i,j,k
      double precision local_coe(ns,ns,npl+1)
      real size , ee(2), psf_model(ns,ns)
      real res_factor

      call get_image_list(iexpo,IMAGE_FILE,nchip,DIR_OUTPUT)
      call get_PREFIX_expo(IMAGE_FILE(1),PREFIX)

      filename=trim(DIR_OUTPUT)//'/rescale/'//trim(PREFIX)//
     .'_factor.dat'
      open(unit=91,file=filename,status='old',iostat=ierror)
      if (ierror.ne.0) then
        write(*,*) 'cannot find rescale factor file' 
        stop
      endif
      rewind 91
      read(91,*) res_factor
      close(91)
      
      filename=trim(DIR_OUTPUT)//'/dat_starcomp/'//trim(PREFIX)//
     .'_star_comp_expo_v2.dat'
      open(unit=11,file=trim(filename),status='replace')
      rewind 11

      filename=trim(DIR_OUTPUT)//'/result/'//trim(PREFIX)//
     .'_star_comp_expo.dat'
      open(unit=10,file=filename,status='old',action='read'
     .,iostat=ierror)
      if (ierror.ne.0) then
          write(*,*) filename, ' does not exist!!'
          return
      endif
      rewind 10

      do while (ierror.ge.0)
        read(10,*,iostat=ierror) ichip, nstar, valid
        if (ierror.lt.0) cycle
        proc_error = 0
        
50      if (valid.lt.0) then
          write(11,*) ichip, nstar, valid
          do chip_circle=1,nstar
            read(10,*,iostat=ierror) (psf_para(u),u=1,8)
            write(11,*) (-999.0 ,u=1,8)
          enddo
          cycle
        endif

        call get_PREFIX(IMAGE_FILE(ichip),PREFIXc)
        filename=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIXc)
     . //'_PSF_coe_local.dat'
        open(unit=13,file=filename,status='old',iostat=ierror)
        rewind 13
        read(13,*) nstar_coe,status
        if (status.eq.-1) then
          proc_error=1
        else
          do i=1,ns
            do j=1,ns
              read(13,*) (local_coe(i,j,k),k=1,npl+1)
            enddo
          enddo
        endif
        close(13)
        if (proc_error.eq.1) then
          valid = -1
          goto 50
        endif

        call get_chip_id(IMAGE_FILE(ichip),chip_index)
        write(11,*) ichip, nstar, valid
        do chip_circle=1,nstar
          read(10,*,iostat=ierror) (psf_para(u),u=1,8)
          if (ierror.lt.0) then
              write(*,*) 'Error in reading file:',filename
              return
          endif

          undefined_flag = 0
          do u=1,8
            if (isnan(psf_para(u))) then
                ! NaN check（NaN ont equal to itself）
                write(*,*) 'NaN detected!',filename
                undefined_flag = 1
                exit
            endif
            if (abs(psf_para(u)).gt.1.0d300) then
                ! check Inf
                write(*,*) 'Inf detected!',filename
                undefined_flag = 1
                exit
            endif
          enddo
          if (undefined_flag.eq.1) then
            write(11,*) (psf_para(u),u=1,8)

          else
            px = psf_para(1)
            py = psf_para(2)
            call get_PSF_model_hierarchical(chip_index,px,py
     .,res_factor,local_coe,psf_model)
            if (isnan(psf_model(1,1))) then
              write(11,*) -999.0,-999.0,(psf_model(1,1), u=1,6)
            else
              call get_power_all(ns,ns,psf_model,ee,size,0.02)
              write(11,*) (psf_para(u),u=1,5),size,(ee(v),v=1,2)
            endif

          endif
        enddo
      enddo
      close(10)
      close(11)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ITP_norm_PSF(nsam,npsam,image,posi,ns
     .,npp,nppx,nx,ny,PSF_coe)
      implicit none

      integer ns,npp,nppx,nsam,npsam,nx,ny
      double precision posi(npsam,2),PSF_coe(ns,ns,npp+1)
      real image(npsam,ns,ns)
      integer i,j,k
      real arr(nsam,3),coep(npp),coe0(1)

      do i=1,ns
        do j=1,ns
          do k=1,nsam
            arr(k,1)=2.0*(posi(k,1) / dble(nx)) - 1.0
            arr(k,2)=2.0*(posi(k,2) / dble(ny)) - 1.0
            arr(k,3)=image(k,i,j)
          enddo
          call fit_2D_2(nsam,nsam,arr,1,coe0)
          call fit_2D_2(nsam,nsam,arr,npp,coep)
          do k=1,npp
            PSF_coe(i,j,k)=coep(k)
          enddo
          PSF_coe(i,j,npp+1)=coe0(1)
        enddo
      enddo

      return
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      subroutine ITP_norm_PSF_cov(nsam,npsam,image,posi,ns
     .,npp,nppx,nx,ny,PSF_coe,sigmarr,comat)
      implicit none

      integer ns,npp,nppx,nsam,npsam,nx,ny
      double precision posi(npsam,2),PSF_coe(ns,ns,npp+1)
      real image(npsam,ns,ns)
      integer i,j,k
      real arr(nsam,3),coep(npp),coe0(1)
      real comat(npp,npp)
      real sigma, sigmarr(ns,ns)

      do i=1,ns
        do j=1,ns
          do k=1,nsam
            arr(k,1)=2.0*(posi(k,1) / dble(nx)) - 1.0
            arr(k,2)=2.0*(posi(k,2) / dble(ny)) - 1.0
            arr(k,3)=image(k,i,j)
          enddo
          call fit_2D_2(nsam,nsam,arr,1,coe0)
          call fit_2D_2_cov(nsam,nsam,arr,npp,coep,sigma,comat)
          do k=1,npp
            PSF_coe(i,j,k)=coep(k)
          enddo
          PSF_coe(i,j,npp+1)=coe0(1)
        enddo
      enddo

      return
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! covariance matrix    
      subroutine fit_2D_2_cov(n,np,arr,nc,c,sigma2,coe_1)
      implicit none

c      fit a 2D function with f(x,y)=\Sigma_{i=1}^{nc}func_i(x,y)c_i

      integer n,np,nc
      real arr(np,3),c(nc)
      real coe(nc,nc),coe_1(nc,nc),vec(nc),temp(nc)
      real x,y,f,smax
      real fit_func_2
      external fit_func_2
      integer i,j,u,v
      real sigma2, ssr, y_pred

      coe=0.
      vec=0.

      do i=1,n
        x=arr(i,1)
        y=arr(i,2)
        f=arr(i,3)
        do j=1,nc
          temp(j)=fit_func_2(x,y,j)
        enddo
        do u=1,nc
          do v=1,nc
            coe(u,v)=coe(u,v)+temp(u)*temp(v)
          enddo
          vec(u)=vec(u)+f*temp(u)
        enddo
      enddo

      call matrix_inverse(coe,nc,coe_1)

      do j=1,nc
        c(j)=0.
        do i=1,nc
          c(j)=c(j)+coe_1(j,i)*vec(i)
        enddo
      enddo

      ssr = 0.
      do i=1,n
        x = arr(i,1)
        y = arr(i,2)
        f = arr(i,3)
        
        ! 1. 计算拟合模型在 (x,y) 处的预测值
        y_pred = 0.
        do j=1,nc
          y_pred = y_pred + c(j) * fit_func_2(x,y,j)
        enddo
        
        ! 2. 累加残差平方和
        ssr = ssr + (f - y_pred)**2
      enddo

      ! 3. 估算本底噪声方差 (除以自由度 n - nc)
      sigma2 = ssr / real(n - nc)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
