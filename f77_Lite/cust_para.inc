      integer chipnx ,chipny , Camera_ccd_num
      parameter (chipnx=2046,chipny=4094)
      parameter (Camera_ccd_num=62)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c f77_lite: the PCA-decomposition parameters (rescale_size, procs_pn,
c work_pn, nblocks, n_pcs, npp6th, nmax_star_pchip) belonged exclusively
c to the PSF_Ms=1 multi-scale PSF reconstruction, which is gone.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
