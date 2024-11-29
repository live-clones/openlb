# You may need to add to the PATH variable! (but better configure in your terminal depending on machine)
#export PATH := /usr/local/cuda-12.6/bin:$(PATH)
#export LIBRARY_PATH := $(LIBRARY_PATH):/usr/local/cuda-12.6/lib64
export MPI_ROOT := ${MPI_ROOT}#/home/hbrs/openmpi#/usr/lib/x86_64-linux-gnu/openmpi
# export PATH := ${MPI_ROOT}/include:$(PATH)

# Include config.mk environment (optional)
-include $(OLB_ROOT)/config.mk
CXXFLAGS		+= -diag-suppress 20012
# LIBS = -L/lib/x86_64-linux-gnu/ -lhwloc  # libhwloc.so.15:$(LD_LIBRARY_PATH)

# Select mixed compilation mode if separate CUDA compiler is given
ifdef CUDA_CXX
include $(OLB_ROOT)/default.mixed.mk
# otherwise use single CXX for all of OpenLB
else
include $(OLB_ROOT)/default.single.mk
endif
