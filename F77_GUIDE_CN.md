# Fourier_Quad F77 流水线指南

Fortran 77（`f77` / `f77_Lite`）流水线的源码结构、参数、编译、Docker 环境与 HPC runner。项目总览与 C++ 流水线见
[`README_CN.md`](README_CN.md) 与 [`CPP_GUIDE_CN.md`](CPP_GUIDE_CN.md)。

> **English**: see [F77_GUIDE.md](F77_GUIDE.md)

---

## 源码结构

### 源码目录

#### `f77/` - 完整 Fortran 77 流水线

| 文件 | 说明 |
|---|---|
| `main.f` | 主程序入口。初始化 MPI，读取曝光列表，分发流水线各阶段。 |
| `para.inc` | 主参数文件。控制图像尺寸、阶段选择（`PROCESS_stage`）、星表路径、PSF 阶数、stamp 尺寸、阈值及星表列索引。 |
| `path_layout.inc` | 完整 `f77/` 流水线使用的 WFST 兼容产品目录约定。 |
| `cust_para.inc` | 自定义参数：CCD 几何尺寸与 PCA PSF 分解参数。 |
| `sig_para.inc` | F6 mode-bar 噪声平面估计器参数。 |
| `pre_process.f` | **阶段 1**：平场、掩膜处理、背景/噪声估计。 |
| `proc_astrometry.f` | **阶段 2**：基于 Gaia 参考星表的天体测量校准。 |
| `proc_source.f` | **阶段 3**：源检测与 stamp 提取。 |
| `proc_FFT_st1.f` | **阶段 4**：星系 stamp 的第一阶段傅里叶变换。 |
| `proc_PSF.f` | **阶段 5**：基于恒星 stamp 的 PSF 建模（局域多项式拟合）。 |
| `proc_psfreconsV2.f` | PSF PCA 重建（多尺度，仅 `PSF_Ms=1` 时启用）。 |
| `00_psf_module.f` | Fortran 90 模块，用于全局 PSF PCA 存储。 |
| `proc_FFT_st2.f` | **阶段 6**：第二阶段傅里叶变换。 |
| `proc_shear.f` | **阶段 7**：Fourier\_Quad 剪切估计。 |
| `proc_info.f` | **阶段 8**：单次曝光统计信息收集。 |
| `proc_combine_shear_catalog.f` | **阶段 9**：剪切星表合并与标定。 |
| `astrometry_calib.f` | 天体测量校准工具例程。 |
| `ex_star.f` | 恒星提取与分类。 |
| `rw_fits_image.f` | FITS 图像读写例程（基于 CFITSIO）。 |
| `universal.f` | 通用工具函数。 |
| `univ_imag_proc.f` | 图像处理工具。 |
| `press.f` | Numerical Recipes 例程（排序、统计、插值）。 |
| `mpi_routines.f` | MPI 分发辅助函数（`mpi_distribute`）。 |
| `FFTPACK.f` | FFTPACK 5.1 快速傅里叶变换库。 |
| `Makefile` | 构建文件。使用 `mpif77`，链接 LAPACK、BLAS 和 CFITSIO。 |

#### `f77_Lite/` - 精简 Fortran 77 流水线

文件集与 `f77/` 接近，但不含 `00_psf_module.f`（PCA PSF 重建已移除）。8 个编译期
开关均冻结为生产取值，死代码分支已物理删除。

## 流水线阶段

流水线共 9 个阶段。阶段执行由 `PROCESS_stage` 参数（定义于 `para.inc` /
`LensingConfig.hpp`）控制，该参数为素因数之积。当 `PROCESS_stage` 能被某阶段的
素数整除时，该阶段执行。默认值 `2 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23` 启用全部
阶段。

| 阶段 | 素数 | 函数 | 说明 |
|---|---|---|---|
| 1 | 2 | `pre_process` / `PreProcess` | 读取 FITS 图像，平场与掩膜校正，背景噪声估计（F6 mode-bar 估计器）。 |
| 2 | 3 | `proc_astrometry` / `Astrometry` | 基于 Gaia 参考星表的天体测量校准，WCS 拟合。 |
| 3 | 5 | `proc_source` / `SourceExtractor` | 源检测、去混叠、stamp 提取。 |
| 4 | 7 | `proc_FFT_st1` / `FourierTransformSt1` | 星系 stamp 第一阶段傅里叶变换。 |
| 5 | 11 | `proc_PSF` / `PSFModel` | 基于恒星 stamp 的 PSF 建模（局域多项式拟合）。可选 PCA 重建（`PSF_Ms=1`）。 |
| 6 | 13 | `proc_FFT_st2` / `FourierTransformSt2` | 第二阶段傅里叶变换。 |
| 7 | 17 | `proc_shear` / `ShearMeasurement` | 基于四阶傅里叶矩的 Fourier\_Quad 剪切估计。 |
| 8 | 19 | `proc_info` / `ExposureInfo` | 收集单次曝光统计（PSF FWHM、恒星数等）。 |
| 9 | 23 | `proc_combine_shear_catalog` / `CatalogCombiner` | 跨曝光合并剪切星表并应用标定校正。 |

如需跳过某阶段，将 `PROCESS_stage` 除以该阶段的素因数即可。例如设置
`PROCESS_stage = 2 * 3 * 5 * 7 * 11 * 13 * 17 * 19`（省去 23）可跳过星表合并。

---


## 参数配置

### Fortran 77（`para.inc`、`cust_para.inc`、`sig_para.inc`）

所有编译期参数分布于三个 include 文件：

- **`para.inc`** - 图像尺寸（`npx`、`npy`）、阶段控制（`PROCESS_stage`）、星表
  路径、PSF 参数、stamp 尺寸、检测阈值及星表列索引。
- **`cust_para.inc`** - CCD 几何（`chipnx`、`chipny`、`Camera_ccd_num`）及 PCA
  PSF 分解参数。
- **`sig_para.inc`** - F6 mode-bar 噪声平面估计器参数。

> `ASTROMETRY_CAT`、`SOURCE_CAT`、`FLAT_PATH` 等路径必须与 Docker `.env` 文件中
> 定义的容器挂载路径一致。

`SOURCE_CAT_TILE_PREFIX` 控制外部源星表 tile 的 basename。默认值为 `extern_`，
与 WFST initializer 和 `process_extcat` 保持一致。


## 数据集结构

```text
<dataset>/
├── science/<exposure>/<chip>.fits
├── dqmask/<exposure>/<exposure>_<CCDNUM>.fits
├── expolists/<exposure>.list
├── astrometry/{Head,dat_Chk,dat_Astro}/...
├── stamps/
│   ├── {Norm,cat_Orig,dat_StarCanInfo,fits_StarCan}/<exposure>/...
│   ├── {fits_StarCanN,fits_StarCanP,dat_SrcInfo,fits_Src}/<exposure>/...
│   ├── {fits_Noise,fits_SrcP,dat_PsfFit,fits_PsfLocal}/<exposure>/...
│   ├── {dat_Shear,dat_StarXY,fits_PsfResi}/<exposure>/...
│   └── {dat_StarInfo,fits_StarP,fits_PsfSrc,dat_ExpoInfo}/...
└── result/<exposure>_all.cat
```

## 源码编译

### Fortran 77

前置条件：`gfortran`、`mpif77`（MPICH）、CFITSIO、LAPACK/BLAS。

```bash
cd f77          # 或 f77_Lite
# 编辑 para.inc 设置星表路径和参数
make
# 可执行文件：./Fourier_Quad_Pipe
```

Makefile 支持变量覆盖：

```bash
make LAPACK_LIB_DIR=/path/to/lapack/lib \
     CFITSIO_LIB_DIR=/path/to/cfitsio/lib
```

本次 I/O 迁移的本地验证环境为 GNU Fortran 15.2.0、Open MPI 5.0.10、
CFITSIO 4.6.4 与 LAPACK/BLAS 3.11.0。现代 GNU Fortran 编译仓库既有 FFTPACK
隐式参数接口时需要以下兼容覆盖：

```bash
export SCIENCE_PREFIX=/path/to/scientific-stack
make FC=mpifort \
     FFLAGS='-mcmodel=medium -w -fallow-argument-mismatch' \
     LAPACK_LIB_DIR="$SCIENCE_PREFIX/lib" \
     CFITSIO_LIB_DIR="$SCIENCE_PREFIX/lib"
```

生产集群应加载站点提供的 MPI/GNU Fortran、CFITSIO、LAPACK/BLAS modules，或使用
仓库已记录的容器工具链；通用编译与运行说明保持相对路径，不依赖本机绝对路径。


### 运行流水线

```bash
mpirun -np <N> ./Fourier_Quad_Pipe <EXPO_LIST>                 # Fortran
```

---


## Docker 环境


Docker 环境提供可复现的构建工具链，无需手动安装编译器和科学库。


### Docker 目录

#### `f77_docker/`

构建与 pilogin 集群工具链匹配的开发镜像：

| 组件 | 版本 |
|---|---|
| 基础镜像 | Rocky Linux 8.10 (glibc 2.28) |
| GCC / G++ / GFortran | 4.8.5（源码编译） |
| MPICH | 4.1.2 (ch4:ofi 设备，Slurm 启动器) |
| CFITSIO | 4.3.1 |
| LAPACK + 参考 BLAS | 3.8.0 |

关键文件：`Dockerfile`、`compose.yaml`、`.env.example`、
`scripts/verify-image.sh`、`scripts/check-public-repo.sh`、`runner/`（HPC 部署
脚本）、`config/`、`patches/`、`checksums.sha256`、`SOURCES.md`、
`THIRD_PARTY_NOTICES.md`。


### 快速开始

```bash
cd f77_docker
cp .env.example .env
# 编辑 .env：设置 F77_SOURCE_HOST、星表路径和 PROCESS_DATA_HOST

docker compose build
docker compose run --rm FourierQuad-F77
```

进入容器后：

```bash
cd /workspace/f77
make
mpirun -np 4 ./Fourier_Quad_Pipe /data/DataProcess/expo_list.list
```


### 验证镜像

每个 Docker 目录附带验证脚本：

```bash
bash f77_docker/scripts/verify-image.sh f77pipeline-dev:gnu4.8.5
```

Docker 环境详细文档请参见：
- [`f77_docker/README.md`](f77_docker/README.md) / [`f77_docker/README-CN.md`](f77_docker/README-CN.md)

---


## HPC 部署

Docker 环境包含 `runner/` 目录，提供 Slurm/Apptainer 部署脚本。典型流程：

1. 拉取 GHCR 镜像并转换为 SIF：
   ```bash
   bash f77_docker/runner/pull-sif.sh   
   ```

2. 配置环境：
   ```bash
   cp f77_docker/runner/f77pipeline.env.example f77_docker/runner/f77pipeline.env
   # 编辑集群路径
   ```

3. 审计 MPI 兼容性（只读）：
   ```bash
   bash f77_docker/runner/inspect-cluster-mpi.sh
   ```

4. 提交流水线：
   ```bash
   sbatch f77_docker/runner/f77pipeline.slurm
   ```

HPC runner 详细文档请参见：
- [`f77_docker/runner/README.md`](f77_docker/runner/README.md) / [`f77_docker/runner/README-CN.md`](f77_docker/runner/README-CN.md)

---
