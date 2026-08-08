# Source provenance and integrity

The Docker build downloads each dependency from its upstream project over HTTPS
and verifies it against `checksums.sha256` before extraction.

| Component | Version | Upstream archive |
| --- | --- | --- |
| GCC | 4.8.5 | <https://ftp.gnu.org/gnu/gcc/gcc-4.8.5/gcc-4.8.5.tar.gz> |
| GMP | 4.3.2 | <https://ftp.gnu.org/gnu/gmp/gmp-4.3.2.tar.bz2> |
| MPFR | 2.4.2 | <https://ftp.gnu.org/gnu/mpfr/mpfr-2.4.2.tar.bz2> |
| MPC | 0.8.1 | <https://gcc.gnu.org/pub/gcc/infrastructure/mpc-0.8.1.tar.gz> |
| MPICH | 4.1.2 | <https://www.mpich.org/static/downloads/4.1.2/mpich-4.1.2.tar.gz> |
| CFITSIO | 4.3.1 | <https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-4.3.1.tar.gz> |
| LAPACK | 3.8.0 | <https://www.netlib.org/lapack/lapack-3.8.0.tar.gz> |

The base image is pinned by digest. The build additionally asserts Rocky Linux
8.10's GCC 8.5.0 bootstrap compiler and glibc 2.28 before compiling the required
GNU 4.8.5 toolchain.

`patches/gcc-4.8.5-el8-compat.patch` contains the compatibility changes required
to bootstrap GCC 4.8.5 against the Rocky Linux 8 userspace. It is applied only
after the original GCC archive passes its checksum.
