# You may need to add to the PATH variable! (but better configure in your terminal depending on machine)
export PATH := /usr/local/cuda-12.1/bin:$(PATH)
export PATH := ${HOME}/openmpi/bin:$(PATH)

# Include config.mk environment (optional)
-include $(OLB_ROOT)/config.mk#/cpu_gcc_openmpi.mk#/config/gpu_only.mk#gpu_hybrid_mixed.mk#
# LIBS = -L/lib/x86_64-linux-gnu/ -lhwloc  # libhwloc.so.15:$(LD_LIBRARY_PATH)

# Select mixed compilation mode if separate CUDA compiler is given
ifdef CUDA_CXX
include $(OLB_ROOT)/default.mixed.mk
# otherwise use single CXX for all of OpenLB
else
include $(OLB_ROOT)/default.single.mk
endif
