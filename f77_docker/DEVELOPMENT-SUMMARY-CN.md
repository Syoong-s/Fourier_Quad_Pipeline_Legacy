# f77pipeline Docker 与 Apptainer 适配开发总结

更新日期：2026-07-30

## 1. 项目目标

本项目的目标是为 Fourier_Quad F77 pipeline 建立一套可复现的编译与运行
环境，并同时支持：

- 本地 Linux/WSL2 下使用 Docker Compose 反复修改、编译和调试 F77 源码；
- 无 root 权限的 Slurm HPC 集群使用 Apptainer/Singularity 运行；
- 多节点、多进程 MPI 任务；
- pipeline 不同阶段修改 `para.inc` 或源码后重新编译；
- 公开发布 Dockerfile，而不公开私有数据、集群路径或第三方源码压缩包。

开发过程中始终遵循两个边界：

1. 镜像只提供编译器和依赖环境，不包含 f77pipeline 源码与观测数据；
2. 集群上的生产源码和只读输入不直接修改，测试使用独立源码副本和明确授权的
   处理目录。

## 2. 目标软件栈

镜像以 pilogin 原始 F77 运行环境为基准，固定以下版本：

| 组件 | 版本或配置 |
| --- | --- |
| 目标架构 | Linux `x86-64` / `linux/amd64` |
| 基础系统 | Rocky Linux 8.10，镜像 digest 固定 |
| glibc | 2.28 |
| GCC、G++、GFortran | GNU 4.8.5 |
| MPICH | 4.1.2 |
| MPICH ABI | 15:1:3 |
| MPICH device | `ch4:ofi` |
| CFITSIO | 4.3.1 |
| LAPACK | 3.8.0 |
| BLAS | LAPACK 3.8.0 自带 reference BLAS |

最终统一安装位置为：

```text
/opt/gcc-4.8.5
/opt/f77stack/bin
/opt/f77stack/include
/opt/f77stack/lib
```

MPICH、CFITSIO、LAPACK 和 BLAS 不再分散在不同目录。F77 Makefile 可以统一
使用 `/opt/f77stack/lib`，镜像也导出了 `LIBRARY_PATH`、`CPATH`、
`LAPACK_LIB_DIR` 和 `CFITSIO_LIB_DIR`。当前 Makefile 会根据
`CFITSIO_LIB_DIR` 派生 `CFITSIO_LIB`。

## 3. Docker 镜像开发过程

### 3.1 集群环境审计

首先对 pilogin 的编译器、MPICH、glibc、Slurm PMI、Apptainer、libfabric
和 RDMA 用户态库进行了只读审计。审计确认原始环境的核心组合是：

- GFortran 4.8.5；
- MPICH 4.1.2，ABI 15:1:3；
- `ch4:ofi`；
- Hydra 支持 Slurm launcher 和 resource manager；
- glibc 2.28。

MPICH 不能只比较版本号。configure 选项、ABI、device、process manager、
Hydra、PMI 和网络 provider 都会影响多节点兼容性。审计结果记录在
[runner/PILOGIN-AUDIT.md](runner/PILOGIN-AUDIT.md)。

### 3.2 从“本地源码包构建”转向公开可复现构建

早期验证参考了集群上已有的第三方源码包和解压目录，以确认版本和构建方式。
公开版随后改为：

1. 从 GCC、MPICH、HEASARC 和 Netlib 官方地址下载源码；
2. 下载后立即用仓库中的 [checksums.sha256](checksums.sha256) 校验；
3. 校验成功后才解压和编译；
4. Dockerfile 不从本地 `sources/` 复制任何压缩包；
5. 最终镜像不保留第三方源码树或源码压缩包。

源码来源和版本记录在 [SOURCES.md](SOURCES.md)，许可说明记录在
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)。

### 3.3 在 Rocky 8 上构建 GNU 4.8.5

Rocky Linux 8.10 使用 GCC 8.5.0 作为 bootstrap 编译器。GNU 4.8.5 年代较
早，直接在 EL8 用户态构建需要兼容补丁，因此项目保留了
[patches/gcc-4.8.5-el8-compat.patch](patches/gcc-4.8.5-el8-compat.patch)。

Dockerfile 使用 multi-stage build：

1. `base-builder` 安装构建工具并验证 bootstrap GCC 8.5.0 和 glibc 2.28；
2. `source-fetcher` 下载并校验全部源码；
3. `gcc-builder` 构建 GCC/G++/GFortran 4.8.5；
4. `stack-builder` 使用 GNU 4.8.5 构建 MPICH、CFITSIO、LAPACK 和 BLAS；
5. `runtime-dev` 只复制编译器和安装结果，并保留运行期重新编译 pipeline
   所需的 `make`、头文件和基础开发环境。

它是“开发运行镜像”，不是只包含可执行文件的最小运行镜像，因为 pipeline
各阶段需要修改源码后在容器内重新编译。

### 3.4 MPICH 构建调整

MPICH 最终使用接近集群的简化配置：

```text
./configure --prefix=/opt/f77stack
```

这让 MPICH 自动选择与 pilogin 原始构建一致的编译器标志、`ch4:ofi`、
embedded libfabric、embedded hwloc 和 Hydra Slurm 支持。安装前后均检查
`mpichversion` 与 `mpiexec -info`，避免只得到版本相同但功能不同的 MPI。

### 3.5 镜像内不放 pipeline 源码

f77pipeline 源码没有 `COPY` 进镜像。原因是：

- 不同处理阶段需要修改 `para.inc`；
- F77 源码会持续更新；
- 每次修改源码都重建整套 GCC/MPICH 镜像成本过高；
- 编译产物应保留在宿主共享目录，便于后续任务复用。

运行时将宿主 `code/` 挂载到容器 `/workspace/f77`，并以读写方式编译。
镜像本身保持通用和不可变。

## 4. 本地 Docker Compose 适配

[compose.yaml](compose.yaml) 负责把本地目录映射为 pipeline 使用的稳定容器
路径：

| 内容 | 默认容器路径 | 权限 |
| --- | --- | --- |
| F77 源码 | `/workspace/f77` | 读写 |
| 测天星表 | `/data/catalogs/AstroDir` | 只读 |
| 源星表 | `/data/catalogs/ExtSrcDir` | 只读 |
| 平场文件 | `/data/calib/FlatDir` | 只读 |
| EXPO_LIST、FITS 和处理输出 | `/data/DataProcess` | 读写 |

`para.inc` 中的 `ASTROMETRY_CAT`、`SOURCE_CAT` 和 `FLAT_PATH` 必须与对应
容器路径完全一致。`PROCESS_DATA_CONTAINER` 可以自定义，默认 EXPO_LIST
路径会从它派生。

本地 Docker 使用 `HOST_UID` 和 `HOST_GID`，使挂载目录中新生成的目标文件
属于宿主用户。该机制只服务 Docker；Apptainer 天然以提交作业的集群用户
身份运行，不需要 UID/GID 映射、Compose 服务名或固定容器名。

本地开发入口为：

```bash
cd f77_docker
cp .env.example .env
docker compose build
docker compose run --rm FourierQuad-F77
```

## 5. 从 Docker 到 Apptainer

### 5.1 为什么 HPC 不直接使用 Docker

多数 HPC 用户没有 root 权限，也不能访问 Docker daemon。适配方案是：

1. 在外部 Linux 构建并发布 OCI/Docker 镜像；
2. 集群普通用户用 Apptainer/Singularity 转换为只读 SIF；
3. Slurm 为每个 MPI rank 启动一个 `apptainer exec`；
4. 源码、星表、平场和处理数据继续通过 bind mount 提供。

[runner/pull-sif.sh](runner/pull-sif.sh) 从公开 OCI 地址生成 SIF，并拒绝覆盖
已有文件。生产环境推荐在 `OCI_IMAGE_URI` 中使用镜像 digest，而不是可变
tag。

### 5.2 集群目录重构

最终集群运行布局为：

```text
f77pipeline/
├── code/                 # 可修改、可编译的 F77 源码副本
├── runner/               # Slurm、Apptainer、env 和烟雾测试
├── images/               # SIF 或临时 Docker archive
├── apptainer-cache/
├── apptainer-tmp/
└── logs/
```

集群不需要上传完整 GitHub 项目，只需完整上传 `runner/`，再准备 SIF、
`code/` 和数据目录。MPI 测试源码已放入 `runner/tests/`，因此烟雾测试不
依赖仓库根目录。

用户进入 `runner/` 后直接提交 Slurm 文件。作业脚本通过
`SLURM_SUBMIT_DIR` 定位自身，不会额外拼接目录层级。
`F77_RUNNER_DIR` 和 `F77_RUNNER_ENV_FILE` 只是内部高级覆盖入口，正常使用
不需要写入 `f77pipeline.env`。

### 5.3 单一 env 配置

HPC 路径、容器路径、模块和 MPI 选择集中在一个
`runner/f77pipeline.env`。其中包括：

- SIF 和 OCI 地址；
- 源码、星表、平场、处理数据的宿主与容器路径；
- 默认 `F77_EXPO_LIST_CONTAINER`；
- pipeline 可执行文件；
- `HPC_MODULES`；
- `MPI_LAUNCH_MODE`、`MPI_LAUNCHER` 和 `SLURM_MPI_TYPE`；
- 可选的 RDMA/provider bind 和环境变量；
- `make clean` 与并行编译设置。

节点数、rank 数、内存、时限、partition、account 和 QoS 仍属于 Slurm 资源
请求，保留在 `.slurm` 文件或 `sbatch` 命令行中。

通用模板默认 `HPC_MODULES=()`，表示日常环境不需要主动加载模块。只有站点
确实要求时才在 env 中填写模块。pilogin 的新版环境模板使用：

```text
HPC_MODULES=(gcc/12.3.0 openmpi/4.1.6-gcc-12.3.0)
```

这只是对“宿主默认环境不是原 MPICH/GCC 组合”的适配，不会替换 SIF 内的
GNU 4.8.5、MPICH 4.1.2、CFITSIO 4.3.1 或 LAPACK 3.8.0。

## 6. MPI 启动模型

项目验证了两条不同的多节点路线。

### 6.1 宿主 MPICH hybrid 模式

当宿主和容器具有兼容的 MPICH 4.1.2 构建时：

1. Slurm 分配节点；
2. 宿主 MPICH `mpiexec` 通过 Hydra Slurm 集成启动每个 rank；
3. 每个 rank 执行一次 `run-apptainer.sh`；
4. 应用在 SIF 内加载匹配的 MPICH 4.1.2。

这是 `MPI_LAUNCH_MODE=mpiexec` 模式。它要求宿主 MPICH 与容器 MPICH 的
ABI、device、process manager 和网络能力尽量一致。

### 6.2 宿主 OpenMPI 模块与 Slurm PMI2 模式

pilogin 统一加载 GCC 12.3.0 和 OpenMPI 4.1.6 后，不能使用宿主 OpenMPI
`mpirun` 启动容器内链接 MPICH 的程序。正确做法是：

1. 批处理脚本根据 env 加载站点模块；
2. Slurm 使用 `srun --mpi=pmi2` 直接启动每个 Apptainer rank；
3. 容器内 MPICH 4.1.2 通过 Slurm PMI2 初始化；
4. runner 只清除可能污染 MPICH 的 `OMPI_MCA_mtl` 和 `OMPI_MCA_osc`；
5. `SLURM_*` 和 `PMI_*` 变量继续保留。

这条路线的成功不表示 OpenMPI 与 MPICH ABI 兼容，也不允许混用宿主
OpenMPI `mpirun` 和容器 MPICH 应用。

pilogin 的低 rank 多节点测试还确认需要：

- 申请 `--exclusive`；
- 统一设置 `SLURM_CPUS_PER_TASK=1`；
- 统一设置 `SLURM_TRES_PER_TASK=cpu=1`。

这些站点资源修正保存在专用 pilogin Slurm wrapper 中，通用脚本仍由 env
决定 MPI 模式。

## 7. 烟雾测试设计

真实数据运行前必须先通过
[runner/mpi-smoke-test.slurm](runner/mpi-smoke-test.slurm)。它会在共享临时
目录内完成：

1. 使用容器 `mpif77` 编译 `runner/tests/mpi_identity.f`；
2. 使用容器 `mpif77` 和 `/opt/f77stack/lib` 编译
   `runner/tests/mpi_lapack_smoke.f`；
3. 跨 Slurm allocation 启动 MPI identity 测试；
4. 输出每个 rank 的编号、world size 和节点名；
5. 运行一个 `DGESV` 线性方程求解，验证 MPI Fortran 初始化以及 LAPACK/
   BLAS 链接和数值结果；
6. 作业退出时只清理由当前 Job 创建的唯一临时目录。

该测试能够发现以下常见问题：

- 所有进程都错误地成为 rank 0；
- 实际没有跨节点；
- host launcher 与 container MPI 不兼容；
- Slurm PMI 类型错误；
- LAPACK/BLAS 链接失败；
- 共享目录在不同节点不可见；
- OFI、libfabric 或 RDMA provider 初始化失败。

## 8. 验证结果

### 8.1 Docker 镜像验证

[scripts/verify-image.sh](scripts/verify-image.sh) 检查：

- GCC/GFortran 4.8.5；
- glibc 2.28；
- MPICH 4.1.2、ABI 15:1:3 和 `ch4:ofi`；
- Hydra Slurm launcher/resource manager；
- CFITSIO 4.3.1；
- MPICH、CFITSIO、LAPACK、BLAS 的实际库文件与动态链接；
- 两 rank MPI 启动；
- MPI/LAPACK Fortran 和 CFITSIO C 程序的真实编译与运行；
- 最终镜像中不存在 pipeline 源码和第三方源码包。

### 8.2 pilogin 计算节点验证

实测完成了：

- 两节点、四 rank、每节点两个 rank 的 MPI identity 测试；
- MPI/LAPACK 烟雾测试；
- 一次代表性真实数据任务，从预处理到目录合并完成一个 exposure 的五个有效
  chip；
- GCC 12.3.0/OpenMPI 4.1.6 模块环境下的 Slurm PMI2 路线；
- 完整 pipeline 零退出码、空 stderr，并重新生成合并目录。

模块环境下的一次完整运行记录为约 3 分 52 秒 MPI runtime、约 5.4 GiB
峰值内存。它是特定数据和资源配置下的验证记录，不应视为通用性能基准。

其他集群必须重新执行计算节点烟雾测试。相同版本号或相同 SONAME 只能作为
兼容前提，不能替代对 PMI、CPU affinity、网络设备、provider 和跨节点通信
的实际验证。

## 9. 自动检查与公开发布

[scripts/check-public-repo.sh](scripts/check-public-repo.sh) 在发布前检查：

- 是否误放第三方源码归档或超大文件；
- 是否包含私有本地/集群路径；
- Dockerfile 是否从本地源码目录复制依赖；
- Docker Compose 和两个 HPC env 示例的挂载默认值是否一致；
- HPC env 是否满足单一配置契约；
- `runner/` 是否自包含；
- shell 文件是否保持 LF；
- 固定格式 Fortran 测试是否超过 72 列；
- 运行脚本是否保留可执行权限。

Docker 项目已合并到主仓库的 `f77_docker/` 子目录。仓库根部的
[`../.github/workflows/f77pipeline-container.yml`](../.github/workflows/f77pipeline-container.yml)
在手动触发或推送 `v*` tag 时：

1. 运行公开仓库检查；
2. 使用 `f77_docker/` 作为 Docker build context；
3. 构建 `linux/amd64` 镜像；
4. 使用 GitHub `GITHUB_TOKEN` 推送到 GHCR；
5. 生成 `gnu4.8.5`、版本 tag 和 commit SHA tag；
6. 使用 GitHub Actions cache 加速后续源码编译。

本次文档与仓库集成自检使用 Bash 5.3.9、Docker 29.6.2 和 Docker Compose
v5.3.1。已通过 10 个 runner/scripts 文件的 Bash 语法检查、公开仓库检查、
Compose 配置解析、GitHub Actions YAML 解析和 Docker build 静态检查。
生产集群只要求 Bash 4 或更高版本，并以集群提供的 Slurm 和
Apptainer/Singularity 为准。

## 10. 当前项目边界

当前实现已经解决版本复现、源码外置、统一库目录、本地调试、无 root SIF、
Slurm 多 rank、两种 MPI 启动路线、集中 env 和公开发布问题，但仍有以下
边界：

- 镜像只正式验证 `linux/amd64`；
- pilogin 的成功结论不能直接推广到另一套 Slurm/PMI/RDMA 环境；
- GHCR 构建工作流目前依赖 Dockerfile 内置断言和公开仓库检查，完整
  `verify-image.sh` 仍应在本地或发布后执行；
- 同一个源码副本不应被多个作业同时 `make clean` 和重编译；
- 修改 `para.inc` 的阶段应使用独立 `code/` 副本；
- 高速网络出现问题时应基于管理员建议精确 bind provider 文件，不能覆盖
  整个容器 `/lib64`；
- 生产任务应固定 OCI digest，避免可变 tag 导致环境漂移。

## 11. 开发与生产入口

本地开发：

```bash
cd f77_docker
docker compose build
bash scripts/verify-image.sh f77pipeline-dev:gnu4.8.5
docker compose run --rm FourierQuad-F77
```

集群生产：

```bash
cd /shared/project/f77pipeline/runner
bash pull-sif.sh
bash run-apptainer.sh --check
sbatch mpi-smoke-test.slurm
sbatch f77pipeline.slurm
```

pilogin 新版模块环境分别使用
`mpi-smoke-test-pilogin-openmpi.slurm` 和
`f77pipeline-pilogin-openmpi.slurm`。

更详细的操作步骤见 [runner/README-CN.md](runner/README-CN.md)，英文配置
参考见 [runner/README.md](runner/README.md)。
