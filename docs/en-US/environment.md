# cuGMEC Environment Setup and Build

This document describes how to install the dependencies required by cuGMEC, configure environment variables, and build the executable. The target environment is a Linux GPU workstation or GPU cluster.

## 1. Dependency Overview

| Component | Requirement | Download |
|---|---|---|
| CUDA Toolkit | CUDA 12.x | <https://developer.nvidia.com/cuda-toolkit> |
| NCCL | Compatible with CUDA 12.x | <https://developer.nvidia.com/nccl> |
| cuDSS | Compatible with CUDA 12.x | <https://docs.nvidia.com/cuda/cudss/> |
| Intel MPI | Linux Offline Installer | <https://www.intel.cn/content/www/cn/zh/developer/tools/oneapi/mpi-library-download.html> |

The example versions used in this document are CUDA 12.6, NCCL 2.23.4, cuDSS 0.4.0.2, and Intel MPI 2021.14.1.7. Newer versions may be used, but CUDA, NCCL, cuDSS, the GPU driver, and Intel MPI must be mutually compatible.

## 2. Install CUDA and Set Environment Variables

Download the installer for your Linux distribution from the CUDA Toolkit website, and install it according to the official NVIDIA documentation. The following example uses CUDA 12.6:

<pre style="background:#f6f8fa;border-radius:4px;padding:12px;overflow:auto;color:#000;"><code><span style="color:#008000;font-style:italic;"># Edit the environment variable file</span>
<span style="color:#0000ff;">vim</span> ~/.bashrc

<span style="color:#008000;font-style:italic;"># Add the following lines to ~/.bashrc. Replace the path with the actual CUDA directory.</span>
<span style="color:#0000ff;">export</span> CUDA_HOME=/usr/local/cuda-12.6
<span style="color:#0000ff;">export</span> PATH=${CUDA_HOME}/bin${PATH:+:${PATH}}
<span style="color:#0000ff;">export</span> LD_LIBRARY_PATH=${CUDA_HOME}/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}

<span style="color:#008000;font-style:italic;"># Apply the environment variables and check the installation</span>
<span style="color:#0000ff;">source</span> ~/.bashrc
<span style="color:#0000ff;">nvcc</span> --version</code></pre>

## 3. Install NCCL and Set Environment Variables

Download a local installer or archive package from the NCCL website that matches your CUDA version, and copy it to the target cluster. The following example uses NCCL 2.23.4 for CUDA 12.6:

<pre style="background:#f6f8fa;border-radius:4px;padding:12px;overflow:auto;color:#000;"><code><span style="color:#008000;font-style:italic;"># Extract and rename</span>
<span style="color:#0000ff;">tar xvf</span> nccl_2.23.4-1+cuda12.6_x86_64.txz
<span style="color:#0000ff;">mv</span> nccl_2.23.4-1+cuda12.6_x86_64 nccl

<span style="color:#008000;font-style:italic;"># Edit the environment variable file</span>
<span style="color:#0000ff;">vim</span> ~/.bashrc

<span style="color:#008000;font-style:italic;"># Add the following lines to ~/.bashrc. Replace the path with the actual NCCL directory.</span>
<span style="color:#0000ff;">export</span> NCCL_HOME=/path/to/nccl
<span style="color:#0000ff;">export</span> LD_LIBRARY_PATH=${NCCL_HOME}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}

<span style="color:#008000;font-style:italic;"># Apply the environment variables</span>
<span style="color:#0000ff;">source</span> ~/.bashrc</code></pre>

## 4. Install cuDSS and Set Environment Variables

Download a Linux archive package from the cuDSS website that matches CUDA 12.x, and copy it to the target cluster. The following example uses cuDSS 0.4.0.2:

<pre style="background:#f6f8fa;border-radius:4px;padding:12px;overflow:auto;color:#000;"><code><span style="color:#008000;font-style:italic;"># Extract and rename</span>
<span style="color:#0000ff;">tar xvf</span> libcudss-linux-x86_64-0.4.0.2_cuda12-archive.tar.xz
<span style="color:#0000ff;">mv</span> libcudss-linux-x86_64-0.4.0.2_cuda12-archive cudss

<span style="color:#008000;font-style:italic;"># Edit the environment variable file</span>
<span style="color:#0000ff;">vim</span> ~/.bashrc

<span style="color:#008000;font-style:italic;"># Add the following lines to ~/.bashrc. Replace the path with the actual cuDSS directory.</span>
<span style="color:#0000ff;">export</span> CUDSS_HOME=/path/to/cudss
<span style="color:#0000ff;">export</span> LD_LIBRARY_PATH=${CUDSS_HOME}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}

<span style="color:#008000;font-style:italic;"># Apply the environment variables</span>
<span style="color:#0000ff;">source</span> ~/.bashrc</code></pre>

## 5. Install Intel MPI and Set Environment Variables

On the Intel MPI download page, select the Linux offline installer and copy the `.sh` installer to the target cluster. The following example uses Intel MPI 2021.14.1.7:

<pre style="background:#f6f8fa;border-radius:4px;padding:12px;overflow:auto;color:#000;"><code><span style="color:#008000;font-style:italic;"># Run the installer</span>
<span style="color:#0000ff;">sh</span> intel-mpi-2021.14.1.7_offline.sh

<span style="color:#008000;font-style:italic;"># Edit the environment variable file</span>
<span style="color:#0000ff;">vim</span> ~/.bashrc

<span style="color:#008000;font-style:italic;"># Add the following line to ~/.bashrc. Replace the path with the actual Intel oneAPI directory.</span>
<span style="color:#0000ff;">source</span> /path/to/intel/oneapi/setvars.sh

<span style="color:#008000;font-style:italic;"># Apply the environment variables and check the installation</span>
<span style="color:#0000ff;">source</span> ~/.bashrc
<span style="color:#0000ff;">mpirun</span> --version</code></pre>

## 6. Build

The current `Makefile` uses the following variables:

| Variable | Default | Description |
|---|---|---|
| `NVCC` | `nvcc` | CUDA compiler |
| `ARCH` | `sm_80` | Target GPU architecture |
| `NCCL_HOME` | Empty | NCCL installation directory |
| `CUDSS_HOME` | Empty | cuDSS installation directory |

Before building, make sure the environment variables for CUDA, NCCL, cuDSS, and Intel MPI have been configured and applied. If `NCCL_HOME` and `CUDSS_HOME` have already been set in `~/.bashrc` as shown above, and `source ~/.bashrc` has been run, no paths need to be specified manually here. Run `make` directly.

<pre style="background:#f6f8fa;border-radius:4px;padding:12px;overflow:auto;color:#000;"><code><span style="color:#0000ff;">make</span></code></pre>

After a successful build, the `cuGMEC` executable will be generated in the project root directory.

Common GPU architecture examples:

| GPU | `ARCH` |
|---|---|
| NVIDIA A100 | `sm_80` |
| NVIDIA A800 | `sm_80` |
| NVIDIA RTX 4090 | `sm_89` |
| NVIDIA H100 | `sm_90` |
