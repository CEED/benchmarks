# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
# reserved. See files LICENSE and NOTICE for details.
#
# This file is part of CEED, a collection of benchmarks, miniapps, software
# libraries and APIs for efficient high-order finite element and spectral
# element discretizations for exascale applications. For more information and
# source code availability see http://github.com/ceed.
#
# The CEED research is supported by the Exascale Computing Project (17-SC-20-SC)
# a collaborative effort of two U.S. Department of Energy organizations (Office
# of Science and the National Nuclear Security Administration) responsible for
# the planning and preparation of a capable exascale ecosystem, including
# software, applications, hardware, advanced system engineering and early
# testbed platforms, in support of the nation's exascale computing imperative.

function setup_bigmem()
{
   ARCH=$(uname -m)
   if [[ "$ARCH" == "x86_64" ]]; then
     CFLAGS+=""
   elif [[ "$ARCH" == "ppc64" ]]; then
     CFLAGS+=" -m64"
   fi
}

function setup_intelmpi()
{
   MPICC=mpiicc
   MPICXX=mpiicpc
   MPIFC=mpiifort
   MPIF77=mpiifort
   # OpenMPI
   export OMPI_CC="$CC"
   export OMPI_CXX="$CXX"
   export OMPI_FC="$FC"
   # or MPICH
   export MPICH_CC="$CC"
   export MPICH_CXX="$CXX"
   export MPICH_FC="$FC"
}

function setup_mpi()
{
   MPICC=mpicc
   MPICXX=mpicxx
   MPIFC=mpifort
   MPIF77=mpif77
   # OpenMPI
   export OMPI_CC="$CC"
   export OMPI_CXX="$CXX"
   export OMPI_FC="$FC"
   # or MPICH
   export MPICH_CC="$CC"
   export MPICH_CXX="$CXX"
   export MPICH_FC="$FC"
}

function setup_intel()
{
   CC=icc
   CXX=icpc
   FC=ifort

   setup_intelmpi

   CFLAGS="-O3"
   setup_bigmem
   FFLAGS="$CFLAGS"

   # The following options assume GCC:
   # TEST_EXTRA_CFLAGS="-march=native --param max-completely-peel-times=3"
   # "-std=c++11 -pedantic -Wall -fdump-tree-optimized-blocks"

   NATIVE_CFLAG="-xHost"
}

function setup_gcc()
{
   CC=gcc
   CXX=g++
   FC=gfortran

   setup_mpi

   CFLAGS="-O3"
   CXXFLAGS="-O3"
   setup_bigmem
   FFLAGS="$CFLAGS"

   # The following options assume GCC:
   TEST_EXTRA_CFLAGS="-march=native --param max-completely-peel-times=3"
   # "-std=c++11 -pedantic -Wall -fdump-tree-optimized-blocks"

   NATIVE_CFLAG="-march=native"
}

function setup_gcc_no_peel()
{
   CC=gcc
   CXX=g++
   FC=gfortran

   setup_mpi

   CFLAGS="-O3"
   setup_bigmem
   FFLAGS="$CFLAGS"

   # The following options assume GCC:
   TEST_EXTRA_CFLAGS="-march=native"
   # "-std=c++11 -pedantic -Wall -fdump-tree-optimized-blocks"

   NATIVE_CFLAG="-march=native"
}


function setup_clang()
{
   CC=clang
   CXX=clang++
   # Use gfortran
   FC=gfortran

   setup_mpi

   CFLAGS="-O3"
   setup_bigmem
   FFLAGS="$CFLAGS"

   TEST_EXTRA_CFLAGS="-march=native -fcolor-diagnostics -fvectorize"
   TEST_EXTRA_CFLAGS+=" -fslp-vectorize -fslp-vectorize-aggressive"
   TEST_EXTRA_CFLAGS+=" -ffp-contract=fast"

   NATIVE_CFLAG="-march=native"
}


function set_mpi_options()
{
   # Run all tasks on the same node.
   num_proc_node=${num_proc_run}
   compose_mpi_run_command
}


search_file_list LAPACK_LIB \
   "/usr/lib64/atlas/libsatlas.so" "/usr/lib64/libopenblas.so"

valid_compilers="gcc clang intel"
# Number of processors to use for building packages and tests:
num_proc_build=${num_proc_build:-8}
# Default number of processors and processors per node for running tests:
num_proc_run=${num_proc_run:-4}
num_proc_node=${num_proc_run}
# Total memory per node:
memory_per_node=8

# Optional (default): MPIEXEC (mpirun), MPIEXEC_OPTS (), MPIEXEC_NP (-np)

cuda_home=${CUDA_HOME:-/usr/local/cuda}
cuda_path=${cuda_home}/bin
cuda_arch=${cuda_arch:-sm_80}
