      integer chipnx ,chipny , Camera_ccd_num
      parameter (chipnx=2046,chipny=4094)
      parameter (Camera_ccd_num=62)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     parameters for pca decomposition
      
      real rescale_size
      parameter (rescale_size = 1.2)

      integer work_pn , procs_pn , n_pcs
      integer nblocks , npp6th
      integer nmax_star_pchip

      parameter (procs_pn = 40)
      parameter (work_pn = 10)
      parameter (nblocks = 2)

      parameter (n_pcs = 100)
      parameter (npp6th = 28)
      parameter (nmax_star_pchip = 1000000)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
