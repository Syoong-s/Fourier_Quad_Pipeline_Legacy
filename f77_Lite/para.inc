      real pi,arc_convert
      parameter (pi=3.1415926)
      parameter (arc_convert=pi/180.)

      integer npx,npy,strl
      parameter (npx=3000)
      parameter (npy=5000)
      parameter (strl=150)

c ==========================================
c f77_lite: the following build-time branch selectors of the full
c pipeline are FROZEN and their unused branches have been deleted from
c the sources.  They no longer exist as parameters; the code implements
c exactly one behaviour each:
c   ASTROMETRY_trivial = 0  -> Gaia-based astrometry only
c   include_FLAT       = 0  -> no super-flat multiplication
c   include_Mask       = 2  -> per-chip DQ mask from DIR_OUTPUT/dqmask
c   ext_cat            = 1  -> external (SOURCE_CAT) source catalogue
c   ext_PSF            = 0  -> PSF measured from the stars of the frame
c   deblending         = 1  -> de-blending always applied
c   PSF_type           = 1  -> local polynomial PSF fit
c   PSF_Ms             = 0  -> no multi-scale / PCA PSF reconstruction
c Still selectable: PROCESS_stage, CCD_split, gal_smooth, star_smooth.
c ==========================================
      integer PROCESS_stage
      parameter (PROCESS_stage=2*3*5*7*11*13*17*19*23)


      character*200 ASTROMETRY_CAT,SOURCE_CAT
      parameter (ASTROMETRY_CAT
     .='/lustre/home/acct-phyzj/phyzj/jzhang/gaia/gaia_cat_sorted')

c  parameter (SOURCE_CAT='/lustre/home/acct-phyzj/share/dr8_rongpu/gal_cat')

c  parameter (SOURCE_CAT.='/lustre/home/acct-phyzj/phyzj/DECALS/gal_cat')

      parameter (SOURCE_CAT='/lustre/home/acct-phyzj/share/DES'
     .//'/testy/des_y6_cat')

      integer CCD_split,nct,ncx
      parameter (CCD_split=2)
      parameter (nct=12)
      parameter (ncx=3)

      integer nstar_min,psf_order,npo,npox,npl,nplx,nstar_min_local
      parameter (psf_order=8)
      parameter (npo=64)
      parameter (npox=8)
      parameter (nstar_min=npo*3/2)
      parameter (npl=10)
      parameter (nplx=2)
      parameter (nstar_min_local=16)
c	    parameter (nstar_min_local=20)

      integer ns,chip_margin,ns_2,nl,nl_2,flag_thresh,chip_edge_margin
      integer nsns
      parameter (ns=64)
      parameter (nsns=ns*ns)
      parameter (chip_margin=8)
      parameter (ns_2=ns/2)
      parameter (nl_2=ns_2+chip_margin)
      parameter (nl=nl_2*2)
      parameter (flag_thresh=3)
c	    parameter (chip_edge_margin=100)
      parameter (chip_edge_margin=chip_margin)

      real dz_thresh
      parameter (dz_thresh=0.1)

      integer len_g,len_s,nstar_max,ngal_max,npara,len_sam
      parameter (len_g=40)
      parameter (len_s=15)
      parameter (ngal_max=4000)
      parameter (nstar_max=2000)
      parameter (npara=25)
      parameter (len_sam=50)

      integer npd
      parameter (npd=33)	

c 1+3+4+1+5+6+1+7+8+1+9+10+1+11+12+1+13+14+1

      integer blocksize
      parameter (blocksize=200)

      real source_thresh,core_thresh
      parameter (source_thresh=2.)
      parameter (core_thresh=4.)
      integer area_max,area_thresh
      parameter (area_max=ns*ns)
      parameter (area_thresh=6)

      integer gal_smooth,star_smooth
      parameter (gal_smooth=2)
      parameter (star_smooth=2)

      real saturation_thresh,SNR_PSF
      parameter (SNR_PSF=100)
      parameter (saturation_thresh=25000)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real pixel_size
      parameter (pixel_size=0.2628)
c in unit of arcsec
	

      integer iSNR_F,igf1,igf2,ira,idec,ih_flux,iPSF
      integer i_imax,ig1,ig2,ide,ih1,ih2,ipeak,iflag
      integer i_jmax,ih_area,istar,isig
      integer icos2,isin2,iparity
	
      parameter (isig=4)
      parameter (istar=5)
      parameter (ipeak=5)
      parameter (i_imax=6)
      parameter (i_jmax=7)
      parameter (ih_flux=8)
      parameter (ih_area=9)
      parameter (iflag=10)
      parameter (iPSF=11)
      parameter (iSNR_F=12)
      parameter (ira=13)
      parameter (idec=14)
      parameter (igf1=15)
      parameter (igf2=16)
      parameter (ig1=17)
      parameter (ig2=18)
      parameter (ide=19)
      parameter (ih1=20)
      parameter (ih2=21)
      parameter (icos2=22)
      parameter (isin2=23)
      parameter (iparity=24)
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer NMAX_EXPO,NMAX_CHIP
      parameter (NMAX_EXPO=25000)
      parameter (NMAX_CHIP=62)
	
      real g1_c,g2_c
c g_band:
      parameter (g1_c=-0.001)
      parameter (g2_c=-0.0003)
c r_band:
c        parameter (g1_c=0.0004)
c	parameter (g2_c=-0.0)
c z_band: 
c        parameter (g1_c=0.)
c	parameter (g2_c=-0.)

      real chi2_thresh
      parameter (chi2_thresh=0.01)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



