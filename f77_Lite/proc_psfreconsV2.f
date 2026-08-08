c ==========================================
c f77_lite: what remains of the original proc_psfreconsV2.f.
c The PSF_Ms=1 multi-scale / PCA PSF-reconstruction stage (chip_psf_recons,
c chip_res_pca_fit, accumulate_block, get_eigval_vec, interpolate_6th,
c mpi_forcov, get_PSF_model_hierarchical, PSF_rescale, PSF_unscale,
c Plot_residuals_v2) has been deleted together with the never-called
c ITP_norm_PSF_cov / fit_2D_2_cov pair.  Only ITP_norm_PSF survives: it is
c the local 2-D polynomial PSF interpolator that make_PSF_local_fit
c (PSF_type=1) calls.
c ==========================================
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
