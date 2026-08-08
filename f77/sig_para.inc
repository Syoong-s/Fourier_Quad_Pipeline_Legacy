c ==========================================
c Configuration: numerical_fix F6 noise-plane estimator
c Method: Find the per-amplifier background mode, construct an
c         image-only private symmetric mask, then fit all surviving
c         pixel triples. No DQ mask is read or required by F6.
c ==========================================
c
c The fitted response is
c   pix=0.5*((I-Ix)^2+(I-Iy)^2),
c whose expectation is 2*sigma^2 for independent, equal-variance noise.
c Block medians provide only robust seeds for the brightness mode and
c initial noise width. The final plane always uses all eligible triples;
c pix itself is never value-trimmed.
c
c The private mask is local to set_sig. It may read the pre-F6 base
c weight (for example a flat-field validity flag), but never changes it.
c In particular set_sig never reads, constructs, or requires a DQ mask.
c
c Blocks are partitioned nearly evenly with target side sig_blocksize.
c The two compile-time bounds protect the fixed-form work arrays.
      integer sig_blocksize,sig_block_max,sig_max_blocks
      integer sig_min_block_pixels,sig_min_block_triples
      integer sig_min_blocks,sig_hist_nbin,sig_min_mode_count
      integer sig_min_lower_count,sig_min_fit_triples,sig_rdil
      integer sig_clip_niter
      real sig_hist_range,sig_clip_k,sig_lower_quantile
      real sig_median_ratio,sig_min_fit_frac
      real sig_plane_min,sig_max_plane_ratio
      double precision sig_pivot_min
      parameter (sig_blocksize=200)
      parameter (sig_block_max=sig_blocksize*sig_blocksize)
      parameter (sig_max_blocks=2048)
      parameter (sig_min_block_pixels=1000)
      parameter (sig_min_block_triples=1000)
      parameter (sig_min_blocks=4)
      parameter (sig_hist_nbin=256)
      parameter (sig_hist_range=6.0)
      parameter (sig_min_mode_count=500)
      parameter (sig_min_lower_count=1000)
c Q_0.3173105 of the conditional lower half is center minus one sigma.
      parameter (sig_lower_quantile=0.3173105)
      parameter (sig_clip_k=3.0)
      parameter (sig_rdil=2)
c Two mask/fit passes: global sigma0, then one local-plane refinement.
      parameter (sig_clip_niter=2)
      parameter (sig_min_fit_triples=1000)
      parameter (sig_min_fit_frac=0.20)
c median(pix)/sigma^2 for 1.5*chi2_1+0.5*chi2_1.
      parameter (sig_median_ratio=1.2678405)
      parameter (sig_plane_min=1.e-8)
c Plane-value ratio 4 permits a factor-two sigma change across an amp.
      parameter (sig_max_plane_ratio=4.0)
      parameter (sig_pivot_min=1.d-8)
c
c sig_scale converts the final masked-mean fit into the published
c 2*sigma^2 convention and remains a required set_sig argument.
c The B correction below was calibrated by verify_setsig from 24
c deterministic Gaussian frames using the exact symmetric image
c clip, dilation, triple selection, and local-refinement implementation.
c
c C selects the output convention:
c   C = 0.809486 (stage one) records the historical effective scale.
c   C = 1 (stage two, ACTIVE) makes the reported sigma an honest sigma.
c The 2026-07-29 bounded production-threshold scan retained C=1, so
c sig_scale_s2 stays 1.027786. sig_scale_s1 is a historical reference
c and is not co-recalibrated when the production threshold is frozen.
      real sig_scale,sig_scale_s1,sig_scale_s2
      parameter (sig_scale_s1=0.673475)
      parameter (sig_scale_s2=1.027786)
      parameter (sig_scale=sig_scale_s2)
