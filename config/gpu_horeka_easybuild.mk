# Example build config for OpenLB using mixed compilation of CUDA with OpenMPI
#
# Adapted for use of easybuild modules on the HoreKa supercomputer at KIT
#
# Usage:
#  - Copy this file to OpenLB root as `config.mk`
#  - Load modules
#      module purge
#      module use /software/easybuild/modules/all/
#      module load make/4.4.1-GCCcore-13.3.0 
#      module load OpenMPI/5.0.7-GCC-14.2.0-CUDA-12.8.0
#  - Run `make clean; make`
#  - Switch to example directory, e.g. `examples/turbulence/nozzle3d`
#  - Run `make`
#  - Use `mpirun --map-by ppr:2:socket:pe=19 --bind-to core bash -c 'export CUDA_VISIBLE_DEVICES=${OMPI_COMM_WORLD_LOCAL_RANK}; ./nozzle3d' for launch

CXX             := mpic++
CC              := gcc

CXXFLAGS        := -O3 -Wall -march=native -mtune=native
CXXFLAGS        += -std=c++20

PARALLEL_MODE   := MPI

PLATFORMS       := CPU_SISD GPU_CUDA

CUDA_CXX        := nvcc
CUDA_CXXFLAGS   := -O3 -std=c++20

CUDA_ARCH       := 80

FLOATING_POINT_TYPE := float

USE_EMBEDDED_DEPENDENCIES := ON
