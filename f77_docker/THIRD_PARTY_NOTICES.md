# Third-party notices

This repository does not redistribute the dependency source archives. During a
build, the Dockerfile downloads GCC, GMP, MPFR, MPC, MPICH, CFITSIO, and LAPACK
from the upstream locations recorded in `SOURCES.md`.

Each downloaded project remains subject to its own license. Review the license
files included in each upstream source distribution before redistributing a
built image.

The GCC compatibility patch modifies files from the GCC source distribution and
is therefore provided under the same applicable GNU project license terms as
those files. The repository's MIT license does not replace those terms.
