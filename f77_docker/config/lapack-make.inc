####################################################################
# LAPACK 3.8.0 configuration matching the pilogin build.
# BLAS and LAPACK are built with GNU Fortran 4.8.5.
####################################################################

SHELL = /bin/sh

CC     = gcc
CFLAGS = -O3

FORTRAN = gfortran
OPTS    = -O2 -frecursive
DRVOPTS = $(OPTS)
NOOPT   = -O0 -frecursive

LOADER   = gfortran
LOADOPTS =

ARCH      = ar
ARCHFLAGS = cr
RANLIB    = ranlib

TIMER = INT_ETIME

BLASLIB      = ../../libblas.a
CBLASLIB     = ../../libcblas.a
LAPACKLIB    = liblapack.a
TMGLIB       = libtmglib.a
LAPACKELIB   = liblapacke.a
