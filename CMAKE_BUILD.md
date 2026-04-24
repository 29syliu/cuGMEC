# cuGMEC CMake Build Guide

## Build Requirements

- CMake >= 3.18
- CUDA Toolkit
- MPI
- NCCL
- cuDSS
- OpenMP

## Build Instructions

### 1. Basic Build (Using Default Paths)

```bash
mkdir build
cd build
cmake ..
cmake --build . --parallel $(nproc)
```

### 2. Specify Custom Library Paths

If your NCCL or cuDSS libraries are not in default paths, specify custom paths:

```bash
mkdir build
cd build
cmake .. \
    -DNCCL_HOME=/path/to/nccl \
    -DCUDSS_HOME=/path/to/cudss
cmake --build . --parallel $(nproc)
```

### 3. Using Build Script (Recommended)

```bash
chmod +x build.sh
NCCL_HOME=/path/to/nccl CUDSS_HOME=/path/to/cudss ./build.sh
```

## Build Output

The compiled executable is located at `build/bin/cuGMEC`

## Additional Build Options

### Debug Mode
```bash
cmake .. -DCMAKE_BUILD_TYPE=Debug
```

### Specify CUDA Architecture
Edit `CMAKE_CUDA_ARCHITECTURES` in CMakeLists.txt, for example:
```cmake
set(CMAKE_CUDA_ARCHITECTURES 80)  # For A100
set(CMAKE_CUDA_ARCHITECTURES 70)  # For V100
```

## Clean Build

```bash
rm -rf build
```

## Troubleshooting

### NCCL/cuDSS Not Found
Make sure you correctly set the `NCCL_HOME` and `CUDSS_HOME` environment variables or CMake options.

### MPI Compilation Errors
Ensure MPI (e.g., OpenMPI or MPICH) is properly installed and environment variables are correctly configured.

### CUDA Architecture Mismatch
Check if the `CMAKE_CUDA_ARCHITECTURES` setting in CMakeLists.txt matches your GPU.

## Original Compilation Command Reference

Original compilation command:
```bash
nvcc -I/cudss/include -L/cudss/lib -I/nccl/include -L/nccl/lib \
     -l mpi -l nccl -l cufft -l cudss -l cusparse \
     -Xcompiler -fopenmp cuGMEC.cu
```

Now fully reconfigured with CMake, supporting modular management and cross-platform compilation.
