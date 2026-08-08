ARG BASE_IMAGE=quay.io/rockylinux/rockylinux:8.10@sha256:e8a49c5403b687db05d4d67333fa45808fbe74f36e683cec7abb1f7d0f2338c6

FROM ${BASE_IMAGE} AS base-builder

ARG BUILD_JOBS=4
ARG TARGETARCH

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN test -z "${TARGETARCH}" || test "${TARGETARCH}" = "amd64"

RUN dnf -y --allowerasing --enablerepo=powertools \
        --setopt=install_weak_deps=False install \
        autoconf \
        automake \
        binutils \
        bison \
        bzip2 \
        ca-certificates \
        curl \
        curl-devel \
        diffutils \
        file \
        findutils \
        flex \
        gcc \
        gcc-c++ \
        gcc-gfortran \
        gzip \
        hostname \
        libatomic \
        libibverbs \
        libibverbs-devel \
        libstdc++-static \
        libtool \
        librdmacm \
        librdmacm-devel \
        m4 \
        make \
        numactl-devel \
        patch \
        perl \
        python3 \
        tar \
        texinfo \
        which \
        xz \
        zlib-devel \
    && test "$(gcc -dumpfullversion)" = "8.5.0" \
    && test "$(getconf GNU_LIBC_VERSION)" = "glibc 2.28" \
    && dnf clean all

WORKDIR /build

FROM base-builder AS source-fetcher

ARG GCC_URL=https://ftp.gnu.org/gnu/gcc/gcc-4.8.5/gcc-4.8.5.tar.gz
ARG GMP_URL=https://ftp.gnu.org/gnu/gmp/gmp-4.3.2.tar.bz2
ARG MPFR_URL=https://ftp.gnu.org/gnu/mpfr/mpfr-2.4.2.tar.bz2
ARG MPC_URL=https://gcc.gnu.org/pub/gcc/infrastructure/mpc-0.8.1.tar.gz
ARG MPICH_URL=https://www.mpich.org/static/downloads/4.1.2/mpich-4.1.2.tar.gz
ARG CFITSIO_URL=https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-4.3.1.tar.gz
ARG LAPACK_URL=https://www.netlib.org/lapack/lapack-3.8.0.tar.gz

COPY checksums.sha256 /build/dist/checksums.sha256

RUN curl --fail --location --retry 5 --connect-timeout 30 \
        "${GCC_URL}" --output /build/dist/gcc-4.8.5.tar.gz \
    && curl --fail --location --retry 5 --connect-timeout 30 \
        "${GMP_URL}" --output /build/dist/gmp-4.3.2.tar.bz2 \
    && curl --fail --location --retry 5 --connect-timeout 30 \
        "${MPFR_URL}" --output /build/dist/mpfr-2.4.2.tar.bz2 \
    && curl --fail --location --retry 5 --connect-timeout 30 \
        "${MPC_URL}" --output /build/dist/mpc-0.8.1.tar.gz \
    && curl --fail --location --retry 5 --connect-timeout 30 \
        "${MPICH_URL}" --output /build/dist/mpich-4.1.2.tar.gz \
    && curl --fail --location --retry 5 --connect-timeout 30 \
        "${CFITSIO_URL}" --output /build/dist/cfitsio-4.3.1.tar.gz \
    && curl --fail --location --retry 5 --connect-timeout 30 \
        "${LAPACK_URL}" --output /build/dist/lapack-3.8.0.tar.gz \
    && cd /build/dist \
    && sha256sum --check checksums.sha256

FROM source-fetcher AS gcc-builder

ARG BUILD_JOBS=4

COPY patches/gcc-4.8.5-el8-compat.patch /build/patches/gcc-4.8.5-el8-compat.patch

RUN tar -xzf /build/dist/gcc-4.8.5.tar.gz -C /build \
    && tar -xjf /build/dist/gmp-4.3.2.tar.bz2 -C /build/gcc-4.8.5 \
    && tar -xjf /build/dist/mpfr-2.4.2.tar.bz2 -C /build/gcc-4.8.5 \
    && tar -xzf /build/dist/mpc-0.8.1.tar.gz -C /build/gcc-4.8.5 \
    && ln -s gmp-4.3.2 /build/gcc-4.8.5/gmp \
    && ln -s mpfr-2.4.2 /build/gcc-4.8.5/mpfr \
    && ln -s mpc-0.8.1 /build/gcc-4.8.5/mpc \
    && patch --batch --forward --strip=1 \
        --directory=/build/gcc-4.8.5 \
        < /build/patches/gcc-4.8.5-el8-compat.patch

RUN mkdir -p /build/gcc-4.8.5-build \
    && cd /build/gcc-4.8.5-build \
    && /build/gcc-4.8.5/configure \
        --prefix=/opt/gcc-4.8.5 \
        --disable-multilib \
        --enable-languages=c,c++,fortran \
        --disable-nls \
        --with-gnu-ld \
        --with-ld=/usr/bin/ld \
        --with-gnu-as \
        --with-as=/usr/bin/as \
        --enable-bootstrap \
        --without-isl \
    && { make -j"${BUILD_JOBS}" \
        || { \
            find /build/gcc-4.8.5-build -path '*/libsanitizer/config.log' \
                -exec tail -n 240 {} \; 2>/dev/null; \
            exit 1; \
        }; \
    } \
    && make install \
    && test "$(/opt/gcc-4.8.5/bin/gfortran -dumpversion)" = "4.8.5"

FROM gcc-builder AS stack-builder

ARG BUILD_JOBS=4

ENV PATH=/opt/gcc-4.8.5/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV LD_LIBRARY_PATH=/opt/gcc-4.8.5/lib64:/opt/gcc-4.8.5/lib

RUN tar -xzf /build/dist/mpich-4.1.2.tar.gz -C /build \
    && cd /build/mpich-4.1.2 \
    && ./configure --prefix=/opt/f77stack \
    && make -j"${BUILD_JOBS}" \
    && make install \
    && /opt/f77stack/bin/mpichversion \
    && /opt/f77stack/bin/mpiexec -info \
        | grep -F "Launchers available:" | grep -F "slurm" \
    && /opt/f77stack/bin/mpiexec -info \
        | grep -F "Resource management kernels available:" | grep -F "slurm"

RUN tar -xzf /build/dist/cfitsio-4.3.1.tar.gz -C /build \
    && cd /build/cfitsio-4.3.1 \
    && CC=gcc ./configure --prefix=/opt/f77stack \
    && make -j"${BUILD_JOBS}" \
    && make install \
    && test -e /opt/f77stack/lib/libcfitsio.so.10

RUN tar -xzf /build/dist/lapack-3.8.0.tar.gz -C /build

COPY config/lapack-make.inc /build/lapack-3.8.0/make.inc

RUN cd /build/lapack-3.8.0 \
    && make -j"${BUILD_JOBS}" blaslib \
    && make -j"${BUILD_JOBS}" lapacklib \
    && install -m 0644 libblas.a /opt/f77stack/lib/libblas.a \
    && install -m 0644 liblapack.a /opt/f77stack/lib/liblapack.a \
    && test -s /opt/f77stack/lib/libblas.a \
    && test -s /opt/f77stack/lib/liblapack.a

RUN test "$(/opt/gcc-4.8.5/bin/gcc -dumpversion)" = "4.8.5" \
    && test "$(/opt/gcc-4.8.5/bin/gfortran -dumpversion)" = "4.8.5" \
    && /opt/f77stack/bin/mpichversion | grep -F "MPICH Version:" | grep -F "4.1.2" \
    && /opt/f77stack/bin/mpiexec -info \
        | grep -F "Launchers available:" | grep -F "slurm" \
    && /opt/f77stack/bin/mpiexec -info \
        | grep -F "Resource management kernels available:" | grep -F "slurm" \
    && grep -F "Version: 4.3.1" /opt/f77stack/lib/pkgconfig/cfitsio.pc \
    && test -e /opt/f77stack/lib/libmpi.so.12 \
    && test -e /opt/f77stack/lib/libmpifort.so.12 \
    && test -e /opt/f77stack/lib/libcfitsio.so.10 \
    && test -e /opt/f77stack/lib/liblapack.a \
    && test -e /opt/f77stack/lib/libblas.a

FROM ${BASE_IMAGE} AS runtime-dev

ARG F77_UID=1000
ARG F77_GID=1000

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN dnf -y --allowerasing --setopt=install_weak_deps=False install \
        binutils \
        file \
        findutils \
        glibc-devel \
        hostname \
        libatomic \
        libcurl \
        libibverbs \
        librdmacm \
        make \
        openssh-clients \
        shadow-utils \
        which \
        zlib \
    && test "$(getconf GNU_LIBC_VERSION)" = "glibc 2.28" \
    && dnf clean all \
    && groupadd --gid "${F77_GID}" f77 \
    && useradd --uid "${F77_UID}" --gid "${F77_GID}" \
        --create-home --shell /bin/bash f77 \
    && mkdir -p \
        /workspace/f77 \
        /data/catalogs/AstroDir \
        /data/catalogs/ExtSrcDir \
        /data/calib/FlatDir \
        /data/DataProcess \
    && chown -R f77:f77 /workspace/f77 /data

COPY --from=stack-builder /opt/gcc-4.8.5 /opt/gcc-4.8.5
COPY --from=stack-builder /opt/f77stack /opt/f77stack
COPY config/f77stack.conf /etc/ld.so.conf.d/f77stack.conf

ENV PATH=/opt/f77stack/bin:/opt/gcc-4.8.5/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV LD_LIBRARY_PATH=/opt/f77stack/lib:/opt/gcc-4.8.5/lib64:/opt/gcc-4.8.5/lib
ENV LIBRARY_PATH=/opt/f77stack/lib
ENV CPATH=/opt/f77stack/include
ENV PKG_CONFIG_PATH=/opt/f77stack/lib/pkgconfig
ENV LAPACK_LIB_DIR=/opt/f77stack/lib
ENV CFITSIO_LIB_DIR=/opt/f77stack/lib

RUN ldconfig \
    && test "$(gcc -dumpversion)" = "4.8.5" \
    && test "$(gfortran -dumpversion)" = "4.8.5" \
    && mpichversion | grep -F "MPICH Version:" | grep -F "4.1.2" \
    && mpiexec -info \
        | grep -F "Launchers available:" | grep -F "slurm" \
    && mpiexec -info \
        | grep -F "Resource management kernels available:" | grep -F "slurm" \
    && grep -F "Version: 4.3.1" /opt/f77stack/lib/pkgconfig/cfitsio.pc \
    && mpiexec -n 2 /bin/true

LABEL org.opencontainers.image.title="f77pipeline development runtime"
LABEL org.opencontainers.image.description="GNU 4.8.5, MPICH 4.1.2, CFITSIO 4.3.1 and LAPACK 3.8.0"
LABEL f77pipeline.gcc.version="4.8.5"
LABEL f77pipeline.mpich.version="4.1.2"
LABEL f77pipeline.cfitsio.version="4.3.1"
LABEL f77pipeline.lapack.version="3.8.0"

USER f77
WORKDIR /workspace/f77

CMD ["/bin/bash"]
