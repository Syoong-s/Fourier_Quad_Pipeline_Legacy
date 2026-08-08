# f77pipeline 在无 root HPC 上运行的完整中文教程

本教程面向第一次接触 Docker、Apptainer、Slurm 和 MPI 的用户。按照顺序操作后，你将能够：

1. 在自己的 Linux 工作站或 WSL2 中构建 Docker 镜像；
2. 把镜像转换成 HPC 使用的 Apptainer SIF；
3. 在不把 f77pipeline 源码放进镜像的前提下编译源码；
4. 挂载星表、平场、EXPO_LIST 和 FITS 数据；
5. 通过 Slurm 在多个计算节点上运行 MPI；
6. 避免覆盖原始数据；
7. 判断一次运行是否真正成功。

本项目已经在 pilogin 集群完成以下实测：

- Apptainer 1.5.2 可以从 Docker archive 生成 SIF；
- GNU Fortran 4.8.5、MPICH 4.1.2、CFITSIO 4.3.1、LAPACK 3.8.0 均正确；
- 4 个 MPI rank 可以分布到 2 个计算节点；
- 完整 f77pipeline 从 Pre-process 到 combine 全部完成；
- 测试处理 1 个 exposure、5 个 FITS chip，耗时约 4 分钟；
- stderr 为空，最终合并星表成功生成。
- 加载 GCC 12.3.0/OpenMPI 4.1.6 后，宿主 OpenMPI 基线跨节点成功；
- 同一新版模块环境下，`srun --mpi=pmi2` 可以直接启动容器内 MPICH；
- PMI2 模式完整回归耗时 3 分 52 秒，所有处理阶段和结果均正常。

不同集群的 Slurm 和网络配置可能不同，所以即使使用同一镜像，也必须先做本教程中的 MPI 烟雾测试。

## 1. 先理解几个基本概念

### 1.1 Dockerfile

`Dockerfile` 是镜像的构建说明书。它描述基础系统、编译器、MPI 和科学计算库如何安装。

本项目的 Dockerfile 会构建：

- GNU GCC/GFortran 4.8.5；
- MPICH 4.1.2；
- CFITSIO 4.3.1；
- LAPACK 3.8.0；
- reference BLAS。

### 1.2 Docker 镜像

Docker 镜像可以理解为一个只读的软件环境快照。本项目的镜像只包含编译器和库，不包含：

- f77pipeline 源码；
- `para.inc`；
- EXPO_LIST；
- FITS 图像；
- Gaia 或其他星表；
- 平场文件；
- 处理结果。

因此修改 f77pipeline 源码后不需要重新构建镜像，只需要重新挂载源码并编译。

### 1.3 容器

容器是镜像的一次运行实例。容器不是虚拟机，不会模拟另一台物理计算机。它仍然使用宿主机的 CPU、内核、网络和文件系统挂载。

本地 Docker Compose 中的服务名 `FourierQuad-F77` 可以修改；Docker
容器名也可以自行指定。但 HPC 使用 `apptainer exec` 时通常不会创建一个
需要命名和长期保持运行的容器，每个 MPI rank 只是临时执行一次 SIF。

### 1.4 Apptainer 与 SIF

大多数 HPC 不允许普通用户运行 Docker daemon，也不会授予 root 权限。Apptainer 专门用于这种环境。

SIF 是 Apptainer 常用的只读镜像文件，例如：

```text
f77pipeline-gnu4.8.5.sif
```

一个 SIF 文件可以被多个计算节点同时读取，不需要 root，也不需要启动后台服务。

Docker 镜像 tag 和 SIF 文件名都可以自定义，只要后续命令以及
`F77_SIF` 使用同一个名称。

Apptainer 默认使用提交作业的真实集群 UID/GID 访问挂载目录，因此 HPC
配置不需要 Docker Compose 使用的 `HOST_UID` 和 `HOST_GID`。能否读写文件
仍由你的集群账户权限决定。

### 1.5 bind 挂载

bind 挂载把宿主机目录映射到容器内目录。例如：

```text
宿主机：/shared/catalogs/gaia
容器内：/data/catalogs/AstroDir
权限：只读
```

程序在容器内读取 `/data/catalogs/AstroDir`，实际读取的是宿主机 `/shared/catalogs/gaia`。

### 1.6 Slurm、节点和 MPI rank

- Slurm：集群作业调度器；
- 节点：一台计算服务器；
- MPI rank：一个 MPI 进程；
- `--nodes=2`：申请 2 个节点；
- `--ntasks=4`：总共启动 4 个 rank；
- `--ntasks-per-node=2`：每节点放 2 个 rank；
- `--cpus-per-task=1`：每个 rank 使用 1 个 CPU 核。

## 2. 整体运行流程

完整流程如下：

```text
本地 Linux/WSL2
  │
  ├─ Dockerfile 构建 Docker 镜像
  │
  ├─ 验证镜像版本和动态库
  │
  └─ 发布 OCI 镜像，或者 docker save 生成归档
          │
          ▼
HPC 共享文件系统
  │
  ├─ Apptainer 把 OCI/Docker archive 转换成 SIF
  ├─ 挂载源码、星表、平场和数据
  ├─ 使用容器内 mpif77 编译源码副本
  └─ 已验证的 MPI/Slurm 启动器启动多个 apptainer exec
          │
          ▼
多个计算节点上的 f77pipeline MPI rank
```

本项目支持两种已经在 pilogin 验证的启动模式。

模式 A 是 host-MPICH hybrid：

1. Slurm 分配计算节点；
2. 集群 MPICH 4.1.2 的 `mpiexec` 启动进程；
3. 每个进程执行一次 `apptainer exec`；
4. 应用使用 SIF 内相同版本的 MPICH。

模式 B 是 pilogin 当前推荐的 Slurm PMI2：

1. Slurm 分配计算节点；
2. 作业读取 `f77pipeline.env`，并加载其中配置的 GCC 12.3.0/OpenMPI 4.1.6 模块；
3. 不使用宿主 OpenMPI 的 `mpirun` 启动 MPICH 应用；
4. `srun --mpi=pmi2` 为每个 rank 启动一个容器；
5. 容器内 MPICH 4.1.2 通过 Slurm PMI2 完成初始化。

宿主 GCC 12.3.0 不会替换容器内编译器。f77pipeline 仍由 SIF 内的
GFortran 4.8.5 编译，运行时仍链接容器内 MPICH、CFITSIO 和 LAPACK。

后文命令分为两类：

- “本地执行”表示在安装了 Docker 的 Linux/WSL2 计算机上执行；
- “HPC 执行”表示先用 SSH 登录集群，再在登录节点提交 Slurm 作业。

不要在 HPC 登录节点上运行 Docker 构建，也不要在本地计算机上执行
`sbatch`。

## 3. 运行前需要准备什么

### 3.1 本地计算机

需要：

- Linux，或者 Windows 11 的 WSL2；
- Docker；
- Git；
- 能够通过 SSH 登录 HPC；
- 足够的磁盘空间。首次构建会下载和编译多个源码包。

检查 Docker：

```bash
docker info
```

检查 Git：

```bash
git --version
```

### 3.2 HPC

需要：

- Apptainer 或 Singularity；
- Slurm；
- 与镜像兼容的宿主 MPICH，或者经过实测兼容的 Slurm PMI 插件；
- 所有计算节点都能访问的共享文件系统；
- 源码目录可写，因为编译结果要保存到源码副本；
- 数据输出目录可写；
- 星表、原始 FITS 和原始 dqmask 可以只读。

登录集群后检查：

```bash
command -v apptainer
apptainer version
module load gcc/12.3.0 openmpi/4.1.6-gcc-12.3.0
gfortran --version
mpirun --version
srun --mpi=list true
```

pilogin 实测环境为：

| 项目 | 版本 |
| --- | --- |
| Apptainer | 1.5.2 |
| Slurm | 25.11.2 |
| 宿主模块 GFortran | 12.3.0 |
| 宿主模块 MPI | OpenMPI 4.1.6 |
| Slurm MPI plugin | `pmi2` |
| 容器 GFortran | 4.8.5 |
| 容器 MPI | MPICH 4.1.2，ABI 15:1:3，`ch4:ofi` |

如果集群的 `mpiexec` 属于 Open MPI，不要用它直接启动当前 MPICH
镜像。应先验证 `srun` PMI 模式，或者构建与宿主 OpenMPI 匹配的镜像变体。

本项目的 HPC 包装脚本在本地使用 Bash 5.3.9 完成语法验证；生产集群需要
Bash 4 或更高版本，以支持 `HPC_MODULES` 的索引数组语法。编译器和科学
计算库仍全部来自 SIF。

## 4. 规划目录

下面使用示例路径。请把 `/shared/...` 替换成你在集群上的真实共享路径。

```text
/shared/project/f77pipeline/
├── code/                             # 挂载并编译的 f77 源码副本
├── runner/                           # 单独上传的运行目录
│   ├── f77pipeline.slurm
│   ├── mpi-smoke-test.slurm
│   ├── run-apptainer.sh
│   ├── f77pipeline.env
│   └── tests/                        # MPI 烟雾测试源码
├── images/
│   ├── f77pipeline.docker.tar
│   └── f77pipeline-gnu4.8.5.sif
├── apptainer-cache/
├── apptainer-tmp/
└── logs/

/shared/data/cpp_test/
├── f2019/                            # 原始数据，只读使用
└── apptainer-f77-test/               # 隔离测试输出
```

所有这些路径必须在每个计算节点上可见。不要把 SIF、源码或数据放在只对登录节点可见的临时目录。

先登录 HPC，并创建项目使用的目录：

```bash
ssh pilogin

mkdir -p \
    /shared/project/f77pipeline/code \
    /shared/project/f77pipeline/images \
    /shared/project/f77pipeline/apptainer-cache \
    /shared/project/f77pipeline/apptainer-tmp \
    /shared/project/f77pipeline/logs
```

集群只需要完整的 `runner/`，不需要上传整个 GitHub 项目。后文会从本地
项目上传该目录。所有集群运行命令都应先进入：

```bash
cd /shared/project/f77pipeline/runner
```

## 5. 在本地构建 Docker 镜像

在本地计算机进入项目：

```bash
cd Fourier_Quad_Pipeline/f77_docker
```

构建：

```bash
docker build \
    --platform linux/amd64 \
    --build-arg BUILD_JOBS=4 \
    -t f77pipeline-dev:gnu4.8.5 \
    .
```

第一次构建会编译 GCC、MPICH、CFITSIO 和 LAPACK，耗时可能较长。`BUILD_JOBS` 可根据本机 CPU 和内存调整。

验证：

```bash
bash scripts/verify-image.sh f77pipeline-dev:gnu4.8.5
```

成功时最后应显示：

```text
Image verification passed: f77pipeline-dev:gnu4.8.5
```

检查镜像：

```bash
docker image inspect f77pipeline-dev:gnu4.8.5
```

## 6. 把镜像交给 HPC

有两种方式，只需要选择一种。

无论使用哪种镜像传输方式，先在本地项目根目录上传完整的 `runner/`：

```bash
scp -r runner pilogin:/shared/project/f77pipeline/
```

然后登录集群并进入运行目录：

```bash
cd /shared/project/f77pipeline/runner
```

### 6.1 方式 A：通过 GHCR 或其他 OCI registry

这是正式发布时推荐的方法。

将镜像推送到 GHCR 后，在集群创建并编辑配置：

```bash
cp f77pipeline.env.example f77pipeline.env
vi f77pipeline.env
```

把文件中的对应两行改为：

```bash
OCI_IMAGE_URI=ghcr.io/OWNER/REPOSITORY@sha256:IMAGE_DIGEST
F77_SIF=/shared/project/f77pipeline/images/f77pipeline-gnu4.8.5.sif
```

第 13 节会填写该文件中的其余路径。

生产运行应使用 digest，不要只依赖可变的 tag。

然后在集群执行：

```bash
bash pull-sif.sh
```

该操作不需要 root。

### 6.2 方式 B：使用 Docker archive

这是 pilogin 实际验证过的方法，适合镜像尚未发布到 registry 的情况。

在本地生成归档：

```bash
docker save \
    --output f77pipeline-dev-gnu4.8.5.docker.tar \
    f77pipeline-dev:gnu4.8.5
```

计算校验值：

```bash
sha256sum f77pipeline-dev-gnu4.8.5.docker.tar
```

传到集群：

```bash
scp f77pipeline-dev-gnu4.8.5.docker.tar \
    pilogin:/shared/project/f77pipeline/images/
```

在集群再次检查：

```bash
sha256sum \
    /shared/project/f77pipeline/images/f77pipeline-dev-gnu4.8.5.docker.tar
```

两端 SHA-256 必须完全相同。

镜像转换需要解包和压缩，不要在登录节点直接运行。创建：

```text
/shared/project/f77pipeline/build-sif.slurm
```

内容如下：

```bash
#!/usr/bin/env bash

#SBATCH --job-name=f77-sif-build
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/shared/project/f77pipeline/logs/%x-%j.out
#SBATCH --error=/shared/project/f77pipeline/logs/%x-%j.err

set -euo pipefail

export APPTAINER_TMPDIR=/shared/project/f77pipeline/apptainer-tmp
export APPTAINER_CACHEDIR=/shared/project/f77pipeline/apptainer-cache

apptainer build --disable-cache \
    /shared/project/f77pipeline/images/f77pipeline-gnu4.8.5.sif \
    docker-archive:/shared/project/f77pipeline/images/f77pipeline-dev-gnu4.8.5.docker.tar
```

注意 `docker-archive:` 后面只有一个冒号，不是 `docker-archive://`。

提交：

```bash
sbatch /shared/project/f77pipeline/build-sif.slurm
```

查看：

```bash
squeue -u "$USER"
```

完成后：

```bash
sacct -j JOB_ID \
    --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS
```

必须看到：

```text
State=COMPLETED
ExitCode=0:0
```

## 7. 验证 SIF

`pull-sif.sh` 和下面这些登录节点上的独立验证命令不会读取并加载
`HPC_MODULES`，因为 `run-apptainer.sh` 在 MPI 作业中还会被每个 rank
调用。如果 Apptainer/Singularity 本身由 module 提供，请在执行独立命令前
手动加载容器运行时模块，或把 `APPTAINER_BIN` 设为其可执行文件路径。
正式 pipeline 和 MPI 烟雾测试会自行加载 env 中的 `HPC_MODULES`。

查看镜像标签：

```bash
apptainer inspect \
    /shared/project/f77pipeline/images/f77pipeline-gnu4.8.5.sif
```

检查编译器：

```bash
apptainer exec --no-home \
    /shared/project/f77pipeline/images/f77pipeline-gnu4.8.5.sif \
    gfortran --version
```

检查 MPI：

```bash
apptainer exec --no-home \
    /shared/project/f77pipeline/images/f77pipeline-gnu4.8.5.sif \
    mpichversion
```

检查 MPI 动态库：

```bash
apptainer exec --no-home \
    /shared/project/f77pipeline/images/f77pipeline-gnu4.8.5.sif \
    ldd /opt/f77stack/lib/libmpi.so
```

预期版本：

```text
GNU Fortran 4.8.5
MPICH 4.1.2
MPICH ABI 15:1:3
MPICH Device ch4:ofi
```

## 8. 准备 f77 源码副本

镜像内没有 f77pipeline 源码。不要直接修改公共或生产源码目录，先创建副本：

```bash
cp -a \
    /shared/original/f77-source/. \
    /shared/project/f77pipeline/code/
```

以后只修改：

```text
/shared/project/f77pipeline/code
```

如果不同 pipeline 阶段需要不同 `para.inc`，应为每个阶段创建独立源码副本，避免两个作业同时执行 `make clean`。

## 9. 理解 para.inc 路径

`para.inc` 中至少要关注：

- `ASTROMETRY_CAT`
- `SOURCE_CAT`
- `FLAT_PATH`
- `PROCESS_stage`

推荐使用稳定的容器路径：

```text
ASTROMETRY_CAT=/data/catalogs/AstroDir
SOURCE_CAT=/data/catalogs/ExtSrcDir
FLAT_PATH=/data/calib/FlatDir
```

运行时再映射：

```text
/shared/catalogs/gaia        → /data/catalogs/AstroDir:ro
/shared/catalogs/source      → /data/catalogs/ExtSrcDir:ro
/shared/calibration/flat     → /data/calib/FlatDir:ro
```

也可以保留 `para.inc` 中原有的宿主绝对路径，然后把宿主目录挂载到容器内相同的绝对路径。pilogin 实测采用过这种方式。

无论选择哪种方式，`para.inc` 中的字符串必须与 bind 的容器目标路径完全相同。修改 `para.inc` 后必须重新编译。

## 10. 编译源码副本

本节展示的是默认容器路径下的手工诊断命令。如果在
`f77pipeline.env` 中修改了 `F77_SOURCE_CONTAINER` 或其他 bind 目标，请在
手工命令中使用相同的值；日常 Slurm 作业会自动读取 env，不需要修改运行
脚本中的路径。

不要让 Makefile 链接集群上的旧库路径。即使 Makefile 仍保留原默认值，也可以在命令行覆盖：

```bash
apptainer exec \
    --no-home \
    --bind /shared/project/f77pipeline/code:/workspace/f77:rw \
    --pwd /workspace/f77 \
    /shared/project/f77pipeline/images/f77pipeline-gnu4.8.5.sif \
    make -C /workspace/f77 clean
```

编译：

```bash
apptainer exec \
    --no-home \
    --bind /shared/project/f77pipeline/code:/workspace/f77:rw \
    --pwd /workspace/f77 \
    /shared/project/f77pipeline/images/f77pipeline-gnu4.8.5.sif \
    make -C /workspace/f77 \
    LAPACK_LIB_DIR=/opt/f77stack/lib \
    CFITSIO_LIB_DIR=/opt/f77stack/lib \
    CFITSIO_LIB=/opt/f77stack/lib/libcfitsio.so
```

关键点：

- 编译器是镜像内的 `mpif77`；
- LAPACK、BLAS 和 CFITSIO 都来自 `/opt/f77stack/lib`；
- 生成的 `Fourier_Quad_Pipe` 保存在宿主源码副本中；
- SIF 本身保持只读。

正式编译应放在 Slurm 作业中。登录节点上的上述命令只适合作为命令结构说明或非常短的检查。

## 11. 为什么不能直接对原始数据运行

当前 f77pipeline 会从 FITS 路径向上两级推导 `DIR_OUTPUT`。

例如输入：

```text
/shared/data/cpp_test/f2019/science/example_1.fits
```

程序会把输出根目录推导为：

```text
/shared/data/cpp_test/f2019
```

这意味着直接使用原始 FITS 路径可能在原数据树中创建或覆盖：

- `astrometry/`
- `stamps/`
- `result/`
- `rescale/`
- `starxy/`
- `fits_psfresi/`
- `dat_pcs/`
- `dat_starcomp/`
- `expo_info.dat`

因此推荐建立独立的运行数据树，用符号链接引用原始 FITS 和 dqmask。

## 12. 建立安全的数据镜像

先定义两个概念：

```text
原始数据：/shared/data/cpp_test
测试输出：/shared/data/cpp_test/apptainer-f77-test
容器路径：/data/DataProcess
```

创建目录：

```bash
mkdir -p \
    /shared/data/cpp_test/apptainer-f77-test/f2019/science \
    /shared/data/cpp_test/apptainer-f77-test/f2019/dqmask \
    /shared/data/cpp_test/apptainer-f77-test/f2019/stamps \
    /shared/data/cpp_test/apptainer-f77-test/f2019/astrometry \
    /shared/data/cpp_test/apptainer-f77-test/f2019/result \
    /shared/data/cpp_test/apptainer-f77-test/f2019/rescale \
    /shared/data/cpp_test/apptainer-f77-test/f2019/starxy \
    /shared/data/cpp_test/apptainer-f77-test/f2019/fits_psfresi \
    /shared/data/cpp_test/apptainer-f77-test/f2019/dat_pcs \
    /shared/data/cpp_test/apptainer-f77-test/f2019/dat_starcomp
```

先查看原始 chip list，确认其中只包含预期 FITS：

```bash
sed -n '1,100p' \
    /shared/data/cpp_test/f2019/stamps/EXPOSURE_NAME.list
```

为每个 FITS 创建链接。以下仅以一个 chip 为例：

```bash
ln -s \
    /shared/data/cpp_test/f2019/science/EXPOSURE_1.fits \
    /shared/data/cpp_test/apptainer-f77-test/f2019/science/EXPOSURE_1.fits
```

对应 dqmask：

```bash
ln -s \
    /shared/data/cpp_test/f2019/dqmask/EXPOSURE_1.fits \
    /shared/data/cpp_test/apptainer-f77-test/f2019/dqmask/EXPOSURE_1.fits
```

每个 chip 都要建立 science 和 dqmask 两个链接。

新的 chip list 必须使用容器路径：

```text
/data/DataProcess/f2019/science/EXPOSURE_1.fits
/data/DataProcess/f2019/science/EXPOSURE_2.fits
```

保存为：

```text
/shared/data/cpp_test/apptainer-f77-test/f2019/stamps/EXPOSURE_NAME.list
```

新的顶层 EXPO_LIST 也必须使用容器路径：

```text
"/data/DataProcess/f2019/stamps/EXPOSURE_NAME.list"     NCHIP
```

将 `EXPOSURE_NAME` 换成实际 exposure 名称，将 `NCHIP` 换成该列表中的
实际 chip 数，例如 5；不要把占位文字原样保留。

保存为：

```text
/shared/data/cpp_test/apptainer-f77-test/expo_fsingle.list
```

检查所有链接：

```bash
find /shared/data/cpp_test/apptainer-f77-test \
    -type l -printf '%p -> %l\n'
```

检查列表：

```bash
sed -n '1,20p' \
    /shared/data/cpp_test/apptainer-f77-test/expo_fsingle.list

sed -n '1,100p' \
    /shared/data/cpp_test/apptainer-f77-test/f2019/stamps/EXPOSURE_NAME.list
```

## 13. 配置 runner/f77pipeline.env

在 pilogin 使用新版模块时，复制已经验证的 PMI2 配置；如果文件已经
存在，不要再次复制覆盖：

```bash
test -e f77pipeline.env || \
    cp f77pipeline.pilogin-openmpi.env.example \
       f77pipeline.env
```

如果沿用旧版 `f77pipeline.env`，必须补充 `F77_SOURCE_CONTAINER`、
`F77_EXPO_LIST_CONTAINER` 和 Bash 数组 `HPC_MODULES`。脚本会拒绝缺少这些
变量的旧配置，避免悄悄退回脚本内的重复默认值。

用文本编辑器打开当前 `runner/` 中的 `f77pipeline.env`，完整内容应类似：

```bash
OCI_IMAGE_URI=ghcr.io/OWNER/REPOSITORY@sha256:IMAGE_DIGEST
F77_SIF=/shared/project/f77pipeline/images/f77pipeline-gnu4.8.5.sif

F77_SOURCE_HOST=/shared/project/f77pipeline/code
F77_SOURCE_CONTAINER=/workspace/f77

ASTROMETRY_CAT_HOST=/shared/catalogs/gaia
ASTROMETRY_CAT_CONTAINER=/data/catalogs/AstroDir

SOURCE_CAT_HOST=/shared/catalogs/source
SOURCE_CAT_CONTAINER=/data/catalogs/ExtSrcDir

FLAT_PATH_HOST=/shared/calibration/flat
FLAT_PATH_CONTAINER=/data/calib/FlatDir

PROCESS_DATA_HOST=/shared/data/cpp_test/apptainer-f77-test
PROCESS_DATA_CONTAINER=/data/DataProcess
F77_EXPO_LIST_CONTAINER="${PROCESS_DATA_CONTAINER%/}/expo_list.list"

HPC_SHARED_SCRATCH_HOST=/shared/data/cpp_test/apptainer-f77-test

APPTAINER_BIN=apptainer

HPC_EXTRA_BINDS=/shared/data/cpp_test/f2019/science:/shared/data/cpp_test/f2019/science:ro,/shared/data/cpp_test/f2019/dqmask:/shared/data/cpp_test/f2019/dqmask:ro

HPC_MODULES=(gcc/12.3.0 openmpi/4.1.6-gcc-12.3.0)
SITE_ENV_SCRIPT=

MPI_LAUNCH_MODE=srun
MPI_LAUNCHER=
SLURM_MPI_TYPE=pmi2
HPC_SCRUB_OPENMPI_ENV=1

F77_BUILD_JOBS=4
F77_MAKE_CLEAN=1
F77_EXECUTABLE="${F77_SOURCE_CONTAINER%/}/Fourier_Quad_Pipe"
```

这里设置 `HPC_EXTRA_BINDS`，是因为第 12 节的符号链接目标仍是原始 science
和 dqmask 的宿主绝对路径；把这两个原始目录挂载到容器内同一路径后，链接
才能解析。若隔离树中存放的是实际文件而不是符号链接，可以留空。原始目录
始终使用 `:ro`。

`PROCESS_DATA_CONTAINER` 可以自定义，例如 `/data/data_process`。一旦
修改，派生的 `F77_EXPO_LIST_CONTAINER` 会自动使用新的容器目录；EXPO_LIST
和 chip list 内部记录的路径仍必须是容器可见路径。Slurm 和 Apptainer
脚本不需要同步修改 `--pwd`。该目录不要求预先存在于 SIF 中，Apptainer
bind 会提供它。

`F77_SOURCE_CONTAINER` 同时控制源码 bind 目标、编译工作目录和默认可执行
文件位置。`HPC_MODULES` 必须使用 Bash 数组语法；不需要加载宿主模块时写
成 `HPC_MODULES=()`。pipeline 和 MPI 烟雾测试会读取同一个模块列表。

`HPC_SCRUB_OPENMPI_ENV=1` 只移除宿主模块强制设置的 OpenMPI transport
参数，保留 `PMI_*` 和 `SLURM_*`。不要在 PMI2 运行命令上使用
`apptainer --cleanenv`，否则可能删除 MPI 初始化需要的环境变量。

如果改回 host-MPICH 模式，再使用通用的
`f77pipeline.env.example`，并把 `MPI_LAUNCHER` 设置为兼容的宿主 MPICH
绝对路径。

## 14. pilogin 必须处理的 Slurm 差异

### 14.1 低核数多节点任务需要 exclusive

pilogin 会拒绝普通的：

```text
2 个节点 × 每节点 2 个任务
```

错误通常是：

```text
Please add --exclusive or set --ntasks-per-node=40 for multi-node jobs.
```

如果测试时只想每节点运行 2 个 rank，需要：

```bash
#SBATCH --exclusive
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
```

这会预留两个完整节点，但实际只启动 4 个 MPI rank。不要把 `--ntasks-per-node` 假写成 40，同时又只申请 4 个总任务；站点会认为参数不合理。

### 14.2 统一 Slurm CPU 环境变量

实测曾出现：

```text
srun: fatal: cpus-per-task set by two different environment variables
SLURM_CPUS_PER_TASK=2 != SLURM_TRES_PER_TASK=cpu=1
```

在调用 `srun` 或宿主 `mpiexec` 前显式统一：

```bash
export SLURM_CPUS_PER_TASK=1
export SLURM_TRES_PER_TASK=cpu=1
```

这两个值必须与：

```bash
#SBATCH --cpus-per-task=1
```

一致。

新版 pilogin Slurm 包装器已经设置这两个变量。

### 14.3 从 f77pipeline.env 加载新版模块

pilogin 配置模板包含：

```bash
HPC_MODULES=(gcc/12.3.0 openmpi/4.1.6-gcc-12.3.0)
```

通用 pipeline 和烟雾测试脚本读取 env 后，会在所有编译、Apptainer 和
MPI 命令之前执行 `module load "${HPC_MODULES[@]}"`。模块加载失败时作业
立即退出。`SITE_ENV_SCRIPT` 只用于某些集群在 `module` 命令可用前所需的
额外初始化；不要在其中切换成另一套 MPI。

## 15. 先做 MPI 烟雾测试

不要直接运行真实数据。新版模块环境使用专用脚本：

```bash
bash -n mpi-smoke-test-pilogin-openmpi.slurm
sbatch mpi-smoke-test-pilogin-openmpi.slurm
```

脚本已经配置：

- 2 个节点；
- 4 个 rank；
- 每节点 2 个 rank；
- `--exclusive`；
- `--cpus-per-task=1`；
- `srun --mpi=pmi2`；
- 从 `f77pipeline.env` 加载 GCC 12.3.0/OpenMPI 4.1.6；
- 容器内使用 GFortran 4.8.5/MPICH 4.1.2 编译测试。

如站点要求 `--account`、`--qos` 或指定日志目录，先复制站点副本再编辑：

```bash
cp mpi-smoke-test-pilogin-openmpi.slurm \
   mpi-smoke-test.site.slurm
```

成功标准：

- 状态是 `COMPLETED`；
- 退出码是 `0:0`；
- 输出中能看到两个不同节点名；
- rank 不都是 0；
- MPI/LAPACK smoke test 成功；
- stderr 没有 OFI、PMI、Hydra、RDMA 或 Apptainer 错误。

## 16. 创建 pilogin 运行脚本

新版模块环境使用：

```bash
bash -n f77pipeline-pilogin-openmpi.slurm
```

该包装器负责修正 Slurm CPU 变量，再调用通用 `f77pipeline.slurm`。
通用脚本从 `f77pipeline.env` 加载模块，并完成镜像检查、源码编译和 MPI
运行。包装器默认申请
2 个节点、4 个 rank、每节点 2 个 rank、每节点 16 GiB、最长 2 小时。

需要增加 `--account`、`--qos` 或修改日志目录时，复制站点副本：

```bash
cp f77pipeline-pilogin-openmpi.slurm \
   f77pipeline.site.slurm
```

通用模板编译时会执行 Makefile。为确保链接镜像内的库，Makefile 应使用：

```text
LAPACK_LIB_DIR=/opt/f77stack/lib
CFITSIO_LIB_DIR=/opt/f77stack/lib
CFITSIO_LIB=/opt/f77stack/lib/libcfitsio.so
```

镜像已经导出 `LAPACK_LIB_DIR`、`CFITSIO_LIB_DIR` 和 `CFITSIO_LIB`，
因此 Makefile 中使用 `?=` 的旧默认值会被容器环境覆盖。无需加载宿主
CFITSIO、LAPACK 或 MPI 编译器。

## 17. 提交真实作业

提交：

```bash
sbatch f77pipeline-pilogin-openmpi.slurm
```

未提供外部参数时，程序使用 `F77_EXPO_LIST_CONTAINER`。临时运行其他列表
时，可显式覆盖：

```bash
sbatch f77pipeline-pilogin-openmpi.slurm \
    /data/DataProcess/expo_fsingle.list
```

如果创建了站点副本，则提交 `f77pipeline.site.slurm`。EXPO_LIST 必须
使用容器内路径；若配置采用同路径 bind，也可以使用在容器中可见的宿主
绝对路径。

查看队列：

```bash
squeue -u "$USER"
```

查看状态：

```bash
sacct -j JOB_ID \
    --format=JobID,JobName,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS,NodeList
```

实时查看 stdout：

```bash
tail -f /shared/project/f77pipeline/logs/f77pipeline-JOB_ID.out
```

查看 stderr：

```bash
cat /shared/project/f77pipeline/logs/f77pipeline-JOB_ID.err
```

正常日志应依次出现：

```text
RNG_SEED rank seed
Total number of EXPOSURE
Pre-process...
Astrometry...
Sources ...
FFT st1...
PSF ...
FFT st2...
Shear ...
Info ...
combine ...
```

最终应看到：

```text
State=COMPLETED
ExitCode=0:0
```

## 18. 如何判断处理结果正确生成

检查总体容量：

```bash
du -sh /shared/data/cpp_test/apptainer-f77-test
```

统计普通文件：

```bash
find /shared/data/cpp_test/apptainer-f77-test \
    -type f -printf '.' | wc -c
```

查看 exposure 汇总：

```bash
sed -n '1,20p' \
    /shared/data/cpp_test/apptainer-f77-test/expo_info.dat
```

查看 result：

```bash
find /shared/data/cpp_test/apptainer-f77-test/f2019/result \
    -maxdepth 1 -type f -printf '%s %f\n'
```

至少应检查：

- 每个有效 chip 的 shear 文件；
- exposure 级 `_expo_info.dat`；
- `_star_comp_expo.dat`；
- 合并后的 `_all.cat`；
- stdout 包含全部阶段；
- stderr 为空；
- Slurm 退出码是 `0:0`。

pilogin 实测结果为：

- 2 个计算节点；
- 4 个 rank；
- 1 个 exposure；
- 5 个有效 chip；
- 78 个普通结果文件；
- 约 575 MiB；
- host-MPICH 模式完整运行约 3 分 57 秒；
- GCC 12.3/OpenMPI 模块加 `srun --mpi=pmi2` 模式约 3 分 52 秒；
- MPI job step 峰值内存约 5.4 GiB；
- 合并星表约 4.7 MiB。

## 19. 日常修改源码后的工作流

只修改 f77 源码、`para.inc` 或处理阶段时：

1. 修改共享文件系统中的源码副本；
2. 不重建 Docker 镜像；
3. 重新提交编译作业；
4. 确认新可执行文件生成；
5. 使用新的隔离输出目录；
6. 提交 pipeline。

只有修改以下内容时才需要重新构建镜像：

- Dockerfile；
- GCC/GFortran 版本；
- MPICH 版本或 configure 参数；
- CFITSIO 版本；
- LAPACK/BLAS 版本；
- Rocky Linux 系统依赖；
- `/opt/f77stack` 构建方式。

## 20. 常见错误

### 20.1 `Please add --exclusive`

原因：站点不允许低任务数的非独占多节点作业。

处理：

```bash
#SBATCH --exclusive
```

同时保证：

```text
nodes × ntasks-per-node = ntasks
```

### 20.2 `SLURM_CPUS_PER_TASK != SLURM_TRES_PER_TASK`

原因：Slurm 新旧 CPU 环境变量冲突，MPI launcher 拒绝启动。

处理：

```bash
export SLURM_CPUS_PER_TASK=1
export SLURM_TRES_PER_TASK=cpu=1
```

如果 stdout 为 0 字节且没有 RNG seed，说明应用 rank 没有启动。host-MPICH
模式通常出现 `hydra_pmi_proxy` step；PMI2 模式通常显示 `srun` 启动的
Apptainer step。

### 20.3 看到多个独立的 rank 0

原因：MPI launcher 与容器内 MPI/PMI 不兼容，或者每个容器都作为独立单进程启动。

处理：

- host-MPICH 模式使用兼容的 MPICH 4.1.2 `mpiexec`；
- 宿主只有 OpenMPI 时，不要用它的 `mpirun` 启动 MPICH 程序；
- pilogin 新版模块使用已验证的 `srun --mpi=pmi2`；
- 先提交 `mpi-smoke-test-pilogin-openmpi.slurm`。

### 20.4 宿主是 OpenMPI，容器是 MPICH

OpenMPI 和 MPICH 不是可互换的运行时。宿主默认 GCC 版本不同并不要求
重建镜像，但宿主 `mpirun` 不能直接启动当前 MPICH 应用。

pilogin 的验证方案是：

```text
module load GCC 12.3.0 + OpenMPI 4.1.6
Slurm srun --mpi=pmi2
Apptainer
容器 MPICH 4.1.2 应用
```

该模式必须保留 PMI/Slurm 环境变量。`run-apptainer.sh` 只清除宿主
OpenMPI 的 `OMPI_MCA_mtl` 和 `OMPI_MCA_osc`，不会使用
`--cleanenv` 启动 MPI rank。

### 20.5 找不到 EXPO_LIST

检查：

- 未传入外部参数时，`F77_EXPO_LIST_CONTAINER` 指向正确的容器路径；
- 传给程序的是容器内路径；
- `PROCESS_DATA_HOST` 已挂载到 `PROCESS_DATA_CONTAINER`；
- EXPO_LIST 内的 chip list 也是容器路径；
- chip list 内的 FITS 路径在容器中存在。

### 20.6 程序写到了原始数据目录

原因：FITS 路径决定 `DIR_OUTPUT`。

处理：

- 为每次测试建立独立数据树；
- chip list 使用隔离树的容器路径；
- 原始 FITS 和 dqmask 只读挂载。

### 20.7 链接到了集群库，而不是镜像库

检查：

```bash
apptainer exec \
    --bind /shared/project/f77pipeline/code:/workspace/f77:ro \
    /shared/project/f77pipeline/images/f77pipeline-gnu4.8.5.sif \
    ldd /workspace/f77/Fourier_Quad_Pipe
```

应看到：

```text
libcfitsio.so.10 => /opt/f77stack/lib/...
libmpifort.so.12 => /opt/f77stack/lib/...
libmpi.so.12 => /opt/f77stack/lib/...
libgfortran.so.3 => /opt/gcc-4.8.5/lib64/...
```

### 20.8 OFI、libfabric、InfiniBand 或 RDMA 错误

先确认单节点和双节点 smoke test。不要直接把宿主整个 `/lib64` 覆盖到容器。

只有管理员或计算节点证据明确要求时，才通过 `HPC_EXTRA_BINDS` 添加具体 provider 文件或目录的只读挂载。

pilogin 的宿主 rdma-core 与镜像版本不同，但 SONAME 兼容，实际双节点完整 pipeline 已成功运行。

### 20.9 SIF 是只读的

这是正常行为。不要把源码或输出写入 SIF：

- 源码目录以 `rw` 挂载；
- 处理输出目录以 `rw` 挂载；
- 星表、原始 FITS、dqmask 和平场以 `ro` 挂载。

## 21. 查看、取消和清理

查看自己的任务：

```bash
squeue -u "$USER"
```

查看历史：

```bash
sacct -j JOB_ID
```

取消任务：

```bash
scancel JOB_ID
```

取消前先确认 Job ID 属于自己。取消不可恢复。

清理时先检查：

```bash
du -sh /shared/data/cpp_test/apptainer-f77-test
find /shared/data/cpp_test/apptainer-f77-test \
    -maxdepth 2 -printf '%y %p\n'
```

确认结果不再需要后，再由你自行删除明确的隔离测试目录。不要对共享目录、用户 home 或包含原始数据的父目录执行递归删除。

Docker archive 在 SIF 验证成功后可以清理，但生产环境建议保留：

- Dockerfile；
- 镜像 digest；
- SIF SHA-256；
- Slurm 脚本；
- `f77pipeline.env`；
- `site-env.sh`；
- stdout/stderr；
- f77 源码提交版本或校验值。

## 22. 最短日常检查清单

每次提交前确认：

- [ ] SIF 存在且可读；
- [ ] 容器内 `gfortran` 是 4.8.5；
- [ ] host-MPICH 版本匹配，或 Slurm PMI2 模式已经过烟雾测试；
- [ ] `HPC_MODULES` 包含当前集群需要的编译器、MPI 和 Apptainer 模块；
- [ ] 源码是独立可写副本；
- [ ] Makefile 链接 `/opt/f77stack/lib`；
- [ ] `para.inc` 路径与 bind 目标完全一致；
- [ ] EXPO_LIST 和 chip list 使用容器路径；
- [ ] 原始 FITS、dqmask 和星表为只读挂载；
- [ ] 结果写入独立目录；
- [ ] Slurm 资源参数互相一致；
- [ ] pilogin 低核多节点任务使用 `--exclusive`；
- [ ] 两个 Slurm CPU 环境变量都为 1；
- [ ] 先通过双节点 smoke test；
- [ ] stdout/stderr 写入共享日志目录；
- [ ] 作业结束后检查 `COMPLETED` 和 `0:0`。

## 23. 相关文件

- [`Dockerfile`](../Dockerfile)：软件环境构建说明；
- [`f77pipeline.env.example`](f77pipeline.env.example)：HPC 路径、模块、默认 EXPO_LIST 和 MPI 配置模板；
- [`f77pipeline.pilogin-openmpi.env.example`](f77pipeline.pilogin-openmpi.env.example)：pilogin 新版模块与 PMI2 配置模板；
- [`run-apptainer.sh`](run-apptainer.sh)：Apptainer bind 包装器；
- [`pull-sif.sh`](pull-sif.sh)：从 OCI registry 获取 SIF；
- [`inspect-cluster-mpi.sh`](inspect-cluster-mpi.sh)：集群 MPI 环境审计；
- [`mpi-smoke-test.slurm`](mpi-smoke-test.slurm)：双节点 MPI 测试模板；
- [`mpi-smoke-test-pilogin-openmpi.slurm`](mpi-smoke-test-pilogin-openmpi.slurm)：pilogin 资源设置与 PMI2 烟雾测试入口；
- [`f77pipeline.slurm`](f77pipeline.slurm)：通用 pipeline 作业模板；
- [`f77pipeline-pilogin-openmpi.slurm`](f77pipeline-pilogin-openmpi.slurm)：pilogin 资源设置与 PMI2 pipeline 入口；
- [`PILOGIN-AUDIT.md`](PILOGIN-AUDIT.md)：pilogin 与镜像兼容性记录；
- [`README.md`](README.md)：英文 HPC 说明。
