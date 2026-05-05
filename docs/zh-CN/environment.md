# cuGMEC 环境配置与编译

本文档说明 cuGMEC 所需依赖的安装方式、环境变量设置和编译命令。目标环境为 Linux GPU 工作站或 GPU 集群。

## 1. 依赖总览

| 组件 | 要求 | 下载地址 |
|---|---|---|
| CUDA Toolkit | CUDA 12.x | <https://developer.nvidia.com/cuda-toolkit> |
| NCCL | 与 CUDA 12.x 匹配 | <https://developer.nvidia.com/nccl> |
| cuDSS | 与 CUDA 12.x 匹配 | <https://docs.nvidia.com/cuda/cudss/> |
| Intel MPI | Linux Offline Installer | <https://www.intel.cn/content/www/cn/zh/developer/tools/oneapi/mpi-library-download.html> |

本文示例版本为 CUDA 12.6、NCCL 2.23.4、cuDSS 0.4.0.2、Intel MPI 2021.14.1.7。实际使用时可更新版本，但需要保证 CUDA、NCCL、cuDSS、GPU 驱动和 Intel MPI 彼此兼容。

## 2. 安装 CUDA 并设置环境变量

根据 Linux 发行版从 CUDA Toolkit 官网下载 installer，并按 NVIDIA 官方文档安装。以下以 CUDA 12.6 为例：

<pre style="background:#f6f8fa;border-radius:4px;padding:12px;overflow:auto;color:#000;"><code><span style="color:#008000;font-style:italic;"># 编辑环境变量文件</span>
<span style="color:#0000ff;">vim</span> ~/.bashrc

<span style="color:#008000;font-style:italic;"># 在 ~/.bashrc 中加入，路径改成实际的 CUDA 目录</span>
<span style="color:#0000ff;">export</span> CUDA_HOME=/usr/local/cuda-12.6
<span style="color:#0000ff;">export</span> PATH=${CUDA_HOME}/bin${PATH:+:${PATH}}
<span style="color:#0000ff;">export</span> LD_LIBRARY_PATH=${CUDA_HOME}/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}

<span style="color:#008000;font-style:italic;"># 使环境变量生效并检查安装</span>
<span style="color:#0000ff;">source</span> ~/.bashrc
<span style="color:#0000ff;">nvcc</span> --version</code></pre>

## 3. 安装 NCCL 并设置环境变量

从 NCCL 官网下载与 CUDA 版本匹配的 local installer 或 archive package，并复制到目标集群。以下以 NCCL 2.23.4 + CUDA 12.6 为例：

<pre style="background:#f6f8fa;border-radius:4px;padding:12px;overflow:auto;color:#000;"><code><span style="color:#008000;font-style:italic;"># 解压并重命名</span>
<span style="color:#0000ff;">tar xvf</span> nccl_2.23.4-1+cuda12.6_x86_64.txz
<span style="color:#0000ff;">mv</span> nccl_2.23.4-1+cuda12.6_x86_64 nccl

<span style="color:#008000;font-style:italic;"># 编辑环境变量文件</span>
<span style="color:#0000ff;">vim</span> ~/.bashrc

<span style="color:#008000;font-style:italic;"># 在 ~/.bashrc 中加入，路径改成实际的 NCCL 目录</span>
<span style="color:#0000ff;">export</span> NCCL_HOME=/path/to/nccl
<span style="color:#0000ff;">export</span> LD_LIBRARY_PATH=${NCCL_HOME}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}

<span style="color:#008000;font-style:italic;"># 使环境变量生效</span>
<span style="color:#0000ff;">source</span> ~/.bashrc</code></pre>

## 4. 安装 cuDSS 并设置环境变量

从 cuDSS 官网下载与 CUDA 12.x 匹配的 Linux archive package，并复制到目标集群。以下以 cuDSS 0.4.0.2 为例：

<pre style="background:#f6f8fa;border-radius:4px;padding:12px;overflow:auto;color:#000;"><code><span style="color:#008000;font-style:italic;"># 解压并重命名</span>
<span style="color:#0000ff;">tar xvf</span> libcudss-linux-x86_64-0.4.0.2_cuda12-archive.tar.xz
<span style="color:#0000ff;">mv</span> libcudss-linux-x86_64-0.4.0.2_cuda12-archive cudss

<span style="color:#008000;font-style:italic;"># 编辑环境变量文件</span>
<span style="color:#0000ff;">vim</span> ~/.bashrc

<span style="color:#008000;font-style:italic;"># 在 ~/.bashrc 中加入，路径改成实际的 cuDSS 目录</span>
<span style="color:#0000ff;">export</span> CUDSS_HOME=/path/to/cudss
<span style="color:#0000ff;">export</span> LD_LIBRARY_PATH=${CUDSS_HOME}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}

<span style="color:#008000;font-style:italic;"># 使环境变量生效</span>
<span style="color:#0000ff;">source</span> ~/.bashrc</code></pre>

## 5. 安装 Intel MPI 并设置环境变量

在 Intel MPI 下载页面选择 Linux offline installer，并将 `.sh` 安装文件复制到目标集群。以下以 Intel MPI 2021.14.1.7 为例：

<pre style="background:#f6f8fa;border-radius:4px;padding:12px;overflow:auto;color:#000;"><code><span style="color:#008000;font-style:italic;"># 运行安装程序</span>
<span style="color:#0000ff;">sh</span> intel-mpi-2021.14.1.7_offline.sh

<span style="color:#008000;font-style:italic;"># 编辑环境变量文件</span>
<span style="color:#0000ff;">vim</span> ~/.bashrc

<span style="color:#008000;font-style:italic;"># 在 ~/.bashrc 中加入，路径改成实际的 Intel oneAPI 目录</span>
<span style="color:#0000ff;">source</span> /path/to/intel/oneapi/setvars.sh

<span style="color:#008000;font-style:italic;"># 使环境变量生效并检查安装</span>
<span style="color:#0000ff;">source</span> ~/.bashrc
<span style="color:#0000ff;">mpirun</span> --version</code></pre>

## 6. 编译

当前 `Makefile` 使用以下变量：

| 变量 | 默认值 | 说明 |
|---|---|---|
| `NVCC` | `nvcc` | CUDA 编译器 |
| `ARCH` | `sm_80` | GPU 目标架构 |
| `NCCL_HOME` | 空 | NCCL 安装目录 |
| `CUDSS_HOME` | 空 | cuDSS 安装目录 |

编译前需要确认 CUDA、NCCL、cuDSS 和 Intel MPI 已经配置环境变量并生效。如果已经按前面的步骤在 `~/.bashrc` 中设置了 `NCCL_HOME` 和 `CUDSS_HOME`，并执行过 `source ~/.bashrc`，这里不需要再手动填写路径，直接 `make` 即可。

<pre style="background:#f6f8fa;border-radius:4px;padding:12px;overflow:auto;color:#000;"><code><span style="color:#0000ff;">make</span></code></pre>

编译成功后，项目根目录下会生成 `cuGMEC` 可执行文件。

常见 GPU 架构示例：

| GPU | `ARCH` |
|---|---|
| NVIDIA A100 | `sm_80` |
| NVIDIA A800 | `sm_80` |
| NVIDIA RTX 4090 | `sm_89` |
| NVIDIA H100 | `sm_90` |
