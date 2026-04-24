cmake_minimum_required(VERSION 3.18)
project(cuGMEC CUDA CXX)

# Set CUDA standard
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set output directory for executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Find required packages
find_package(CUDA REQUIRED)
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
find_package(nlohmann_json QUIET)

if(NOT nlohmann_json_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        nlohmann_json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG v3.11.3
    )
    FetchContent_MakeAvailable(nlohmann_json)
endif()

# Optional: Find NCCL and cuDSS if they are installed in non-standard locations
# Users can specify paths via cmake options:
# -DNCCL_HOME=/path/to/nccl
# -DCUDSS_HOME=/path/to/cudss

# NCCL Configuration
if(NOT NCCL_HOME)
    set(NCCL_HOME "/nccl" CACHE STRING "Path to NCCL installation")
endif()

find_library(NCCL_LIB NAMES nccl
    HINTS ${NCCL_HOME}/lib
    NO_DEFAULT_PATH)

if(NOT NCCL_LIB)
    find_library(NCCL_LIB NAMES nccl)
endif()

find_path(NCCL_INCLUDE NAMES nccl.h
    HINTS ${NCCL_HOME}/include
    NO_DEFAULT_PATH)

if(NOT NCCL_INCLUDE)
    find_path(NCCL_INCLUDE NAMES nccl.h)
endif()

# cuDSS Configuration
if(NOT CUDSS_HOME)
    set(CUDSS_HOME "/cudss" CACHE STRING "Path to cuDSS installation")
endif()

find_library(CUDSS_LIB NAMES cudss
    HINTS ${CUDSS_HOME}/lib
    NO_DEFAULT_PATH)

if(NOT CUDSS_LIB)
    find_library(CUDSS_LIB NAMES cudss)
endif()

find_path(CUDSS_INCLUDE NAMES cudss.h
    HINTS ${CUDSS_HOME}/include
    NO_DEFAULT_PATH)

if(NOT CUDSS_INCLUDE)
    find_path(CUDSS_INCLUDE NAMES cudss.h)
endif()

# Verify required libraries are found
if(NOT NCCL_LIB)
    message(FATAL_ERROR "NCCL library not found. Please set -DNCCL_HOME=/path/to/nccl")
endif()

if(NOT NCCL_INCLUDE)
    message(FATAL_ERROR "NCCL include directory not found. Please set -DNCCL_HOME=/path/to/nccl")
endif()

if(NOT CUDSS_LIB)
    message(FATAL_ERROR "cuDSS library not found. Please set -CUDSS_HOME=/path/to/cudss")
endif()

if(NOT CUDSS_INCLUDE)
    message(FATAL_ERROR "cuDSS include directory not found. Please set -CUDSS_HOME=/path/to/cudss")
endif()

# Source files
set(SOURCES
    src/cuGMEC.cu
)

# Create executable
add_executable(cuGMEC ${SOURCES})

# Include directories
target_include_directories(cuGMEC PRIVATE
    src
    lib/BSI/BSplineInterpolation
    ${CUDA_INCLUDE_DIRS}
    ${NCCL_INCLUDE}
    ${CUDSS_INCLUDE}
    ${MPI_CXX_INCLUDE_DIRS}
)

# Link libraries
target_link_libraries(cuGMEC PRIVATE
    ${CUDA_LIBRARIES}
    ${NCCL_LIB}
    ${CUDSS_LIB}
    MPI::MPI_CXX
    OpenMP::OpenMP_CXX
    nlohmann_json::nlohmann_json
    cuda
    cufft
    cusparse
    cublas
)

# CUDA-specific compiler flags
# Automatically detect GPU architecture or allow manual override
set(CMAKE_CUDA_ARCHITECTURES 75 CACHE STRING "CUDA architectures to target")  # Default: sm_75

# C++ compiler flags for CUDA host compilation
set_target_properties(cuGMEC PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
)

# Print configuration info
message(STATUS "======== cuGMEC Configuration ========")
message(STATUS "CUDA_TOOLKIT_ROOT_DIR: ${CUDA_TOOLKIT_ROOT_DIR}")
message(STATUS "NCCL_HOME: ${NCCL_HOME}")
message(STATUS "NCCL_LIB: ${NCCL_LIB}")
message(STATUS "NCCL_INCLUDE: ${NCCL_INCLUDE}")
message(STATUS "CUDSS_HOME: ${CUDSS_HOME}")
message(STATUS "CUDSS_LIB: ${CUDSS_LIB}")
message(STATUS "CUDSS_INCLUDE: ${CUDSS_INCLUDE}")
message(STATUS "MPI_CXX_COMPILER: ${MPI_CXX_COMPILER}")
message(STATUS "========================================")
