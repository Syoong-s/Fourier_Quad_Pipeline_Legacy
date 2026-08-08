# 使用说明

镜像和 HPC 适配的设计演进、关键决策与实测结论见
[DEVELOPMENT-SUMMARY-CN.md](DEVELOPMENT-SUMMARY-CN.md)。

进入克隆后的发行版目录：

```bash
cd Fourier_Quad_Pipeline/f77_docker
```

## 1. 构建镜像

确认 Docker 正常：

```bash
docker info
```

开始构建：

```bash
docker build \
  --platform linux/amd64 \
  --build-arg BUILD_JOBS=4 \
  -t f77pipeline-dev:gnu4.8.5 \
  .
```

`BUILD_JOBS` 可以按内存和 CPU 调整，例如 `8`。第一次会下载并编译 GCC、MPICH、CFITSIO 和 LAPACK，耗时较长。

构建后验证：

```bash
bash scripts/verify-image.sh f77pipeline-dev:gnu4.8.5
```

验证成功时最后会显示：

```text
Image verification passed: f77pipeline-dev:gnu4.8.5
```

## 2. 配置日常挂载目录

推荐使用 Docker Compose。先创建配置文件：

```bash
cp .env.example .env
```

编辑 `.env`：

```bash
nano .env
```

主要配置示例：

```dotenv
IMAGE_NAME=f77pipeline-dev:gnu4.8.5
BUILD_JOBS=4
HOST_UID=1000
HOST_GID=1000

F77_SOURCE_HOST=/实际的f77源码目录

ASTROMETRY_CAT_HOST=/实际的/测天星表目录
ASTROMETRY_CAT_CONTAINER=/data/catalogs/AstroDir

SOURCE_CAT_HOST=/实际的/源星表目录
SOURCE_CAT_CONTAINER=/data/catalogs/ExtSrcDir

FLAT_PATH_HOST=/实际的/平场目录
FLAT_PATH_CONTAINER=/data/calib/FlatDir

PROCESS_DATA_HOST=/实际的数据处理目录
PROCESS_DATA_CONTAINER=/data/DataProcess
```

获取自己的 UID 和 GID：

```bash
id -u
id -g
```

把结果分别填入 `HOST_UID` 和 `HOST_GID`。

特别注意，下面三个容器内路径必须和 `para.inc` 中编译进去的值完全一致：

- `ASTROMETRY_CAT_CONTAINER` 对应 `ASTROMETRY_CAT`
- `SOURCE_CAT_CONTAINER` 对应 `SOURCE_CAT`
- `FLAT_PATH_CONTAINER` 对应 `FLAT_PATH`

例如 `para.inc` 中是：

```text
ASTROMETRY_CAT = '/data/catalogs/AstroDir'
```

那么 `.env` 里也必须是：

```dotenv
ASTROMETRY_CAT_CONTAINER=/data/catalogs/AstroDir
```

`PROCESS_DATA_HOST` 应指向包含 `EXPO_LIST` 和 FITS 图像的数据处理目录。

## 3. 启动开发容器

进入交互式容器：

```bash
docker compose run --rm FourierQuad-F77
```

进入后：

- f77pipeline 源码位于 `/workspace/f77`
- 数据处理目录位于 `$PROCESS_DATA_CONTAINER`
- 库统一位于 `/opt/f77stack/lib`
- 头文件位于 `/opt/f77stack/include`

检查环境：

```bash
gfortran --version
mpif77 -show
mpichversion
```

## 4. 编译 f77pipeline

在容器中：

```bash
cd /workspace/f77
make clean
make -j4
```

*注意：如果当前 Makefile 仍然包含集群上的绝对库路径，需要改为镜像中的统一目录：*

```text
-L/opt/f77stack/lib
-I/opt/f77stack/include
```

常见链接参数可以使用：

```text
-lcfitsio -llapack -lblas
```

MPI Fortran 编译器应使用：

```text
mpif77
```

或：

```text
mpifort
```

源码目录是宿主机挂载进来的，因此：

- 在容器中修改源码，会直接反映到宿主机。
- 在容器中重新编译，生成的目标文件和程序也会保存在宿主机的 f77 目录。
- 删除容器不会删除源码和处理数据。

## 5. 运行 pipeline

进入数据处理目录：

```bash
cd "$PROCESS_DATA_CONTAINER"
```

确认曝光列表：

```bash
head -n 10 EXPO_LIST
```

然后按照 pipeline 输入参数规范运行程序：

```bash
mpiexec -n (ntask) /workspace/f77/Fourier_Quad_Pipe EXPO_LIST
```

具体可执行文件名称和参数以 f77 目录当前 Makefile及各阶段脚本为准。

## 6. 日常推荐流程

每次调整代码后：

```bash
cd Fourier_Quad_Pipeline/f77_docker
docker compose run --rm FourierQuad-F77
```

在容器里：

```bash
cd /workspace/f77
make clean
make -j4

cd "$PROCESS_DATA_CONTAINER"
mpiexec -n (ntask) /workspace/f77/Fourier_Quad_Pipe EXPO_LIST
```

退出容器：

```bash
exit
```

因为使用了 `--rm`，退出后临时容器会自动删除，但源码、编译结果、星表和处理数据都保留在宿主机。

## 7. 什么时候需要重新构建镜像

修改以下内容时才需要重建：

- Dockerfile
- GCC、MPICH、CFITSIO 或 LAPACK 版本
- GCC 兼容补丁
- 镜像系统依赖
- `/opt/f77stack` 的构建配置

重建：

```bash
docker compose build
```

仅修改 f77pipeline 源码、`para.inc` 或处理阶段代码时，不需要重建镜像，只需进入容器重新运行 `make`。

## 8. 无 root、多节点 HPC

Docker Compose 只用于本地开发。在 Slurm 集群上使用 Apptainer 或
Singularity，把公开的 Docker/OCI 镜像转换为只读 SIF，并由宿主机
MPICH 或经过验证的 `srun` PMI 模式为每个 MPI rank 启动一个容器进程。

HPC 环境不需要 `HOST_UID`、`HOST_GID`、Compose 服务名或容器名。
Apptainer 会以提交作业的真实集群用户身份运行。

在本地项目根目录把完整的 `runner/` 上传到集群：

```bash
scp -r runner pilogin:/shared/project/f77pipeline/
```

登录集群，进入 `runner/` 后再执行审计和配置：

```bash
cd /shared/project/f77pipeline/runner
bash inspect-cluster-mpi.sh
cp f77pipeline.pilogin-openmpi.env.example f77pipeline.env
```

该文件统一配置 HPC 使用的宿主/容器路径、宿主模块、MPI 启动方式和默认
EXPO_LIST。Slurm 的节点、进程、内存和时限仍在 `.slurm` 文件中设置。
其他集群应改为复制通用的 `f77pipeline.env.example`，再按该集群实测结果
配置模块和 MPI 启动方式。
默认容器路径沿用当前 Compose 配置：

- `/workspace/f77`
- `/data/catalogs/AstroDir`
- `/data/catalogs/ExtSrcDir`
- `/data/calib/FlatDir`
- `/data/DataProcess`

获得 SIF 后，继续在当前 `runner/` 目录验证并提交：

```bash
bash run-apptainer.sh --check
sbatch mpi-smoke-test-pilogin-openmpi.slurm
sbatch f77pipeline-pilogin-openmpi.slurm
```

通用作业脚本会加载 `f77pipeline.env` 的
`HPC_MODULES=(gcc/12.3.0 openmpi/4.1.6-gcc-12.3.0)`，并通过
`srun --mpi=pmi2` 启动容器内 MPICH 应用。不要使用宿主 OpenMPI 的
`mpirun` 直接启动 MPICH 程序。

MPI 版本相同不代表一定能跨节点运行，还必须核对 MPICH configure 参数、
Slurm PMI、libfabric provider 和 InfiniBand/RDMA 用户态库。面向初学者的
完整中文流程、pilogin 实测参数和故障排查见
[runner/README-CN.md](runner/README-CN.md)，英文参考见
[runner/README.md](runner/README.md)。

本项目已经验证两条 pilogin 路径：兼容宿主 MPICH 的 `mpiexec` hybrid
模式，以及加载 GCC 12.3/OpenMPI 4.1.6 后使用 Slurm PMI2 启动容器
MPICH 4.1.2。两种模式均通过双节点 4-rank 烟雾测试和完整 pipeline；
当前模块化环境推荐第二种。迁移到其他集群时仍需重新测试。

审计同时发现 pilogin 使用 KOS5 的 `libibverbs/librdmacm 37.2`，Rocky
镜像当前是相同 SONAME 的 rdma-core 48.0。先在计算节点使用默认配置测试；
如果管理员或测试结果要求使用宿主 provider，再通过
`HPC_EXTRA_BINDS` 精确添加只读映射，不要覆盖整个容器 `/lib64`。

完整英文说明见 [README.md](README.md)。
