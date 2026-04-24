#!/bin/bash
# Build script for cuGMEC with CMake

# Set default paths if not already set
NCCL_HOME=${NCCL_HOME:-/nccl}
CUDSS_HOME=${CUDSS_HOME:-/cudss}
BUILD_DIR=${BUILD_DIR:-build}

# Create build directory
mkdir -p $BUILD_DIR
cd $BUILD_DIR

# Configure with CMake
cmake .. \
    -DNCCL_HOME=$NCCL_HOME \
    -DCUDSS_HOME=$CUDSS_HOME \
    -DCMAKE_BUILD_TYPE=Release

# Build
cmake --build . --parallel $(nproc)

echo "Build complete! Executable is at: $BUILD_DIR/bin/cuGMEC"
