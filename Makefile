# cuGMEC build

NVCC       ?= nvcc
ARCH       ?= sm_80
CUDSS_HOME ?=
NCCL_HOME  ?=

SRC := src/cuGMEC_main.cu
BIN := cuGMEC

INC_FLAGS := $(if $(CUDSS_HOME),-I$(CUDSS_HOME)/include) \
             $(if $(NCCL_HOME),-I$(NCCL_HOME)/include)
LIB_FLAGS := $(if $(CUDSS_HOME),-L$(CUDSS_HOME)/lib) \
             $(if $(NCCL_HOME),-L$(NCCL_HOME)/lib)

LIBS      := -lmpi -lnccl -lcufft -lcudss -lcusparse
NVCCFLAGS := -arch=$(ARCH) -std=c++20 --expt-relaxed-constexpr -Xcompiler -fopenmp

DEPS := $(wildcard src/*.h) $(wildcard src/kernels/*.cuh)

.PHONY: all clean

all: $(BIN)

$(BIN): $(SRC) $(DEPS)
	$(NVCC) $(INC_FLAGS) $(LIB_FLAGS) $(NVCCFLAGS) $(SRC) $(LIBS) -o $(BIN)

clean:
	rm -f $(BIN)
