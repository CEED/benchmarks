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

function setup_openmpi()
{
   MPICC=mpicc
   MPICXX=mpicxx
   MPIFC=mpifort
   MPIF77=mpif77
   #OMP
   export OMPI_CC="$CC"
   export OMPI_CXX="$CXX"
   export OMPI_FC="$FC"
   export OMPI_F77="$FC"
}


function setup_clang()
{
   # Homebrew open-mpi; default is to use apple clang
   CC=clang
   CXX=clang++
   # Use gfortran
   FC=gfortran

   setup_openmpi

   CFLAGS="-O3"
   FFLAGS="$CFLAGS"
   TEST_EXTRA_CFLAGS="-march=native -fcolor-diagnostics -fvectorize"
   TEST_EXTRA_CFLAGS+=" -fslp-vectorize -fslp-vectorize-aggressive"
   TEST_EXTRA_CFLAGS+=" -ffp-contract=fast"
   # "-std=c++11 -pedantic -Wall"
   NATIVE_CFLAG="-march=native"

   NEK5K_BIGMEM="no" # used by: package-builders/nek5000.sh
}


function setup_gcc_6()
{
   # Homebrew open-mpi with gcc v6
   CC=gcc-6
   CXX=g++-6
   FC=gfortran-6

   setup_openmpi

   CFLAGS="-O3"
   FFLAGS="$CFLAGS"
   OCCA_CXXFLAGS="-O3 -march=native"
   TEST_EXTRA_CFLAGS="-march=native --param max-completely-peel-times=3"
   # "-std=c++11 -pedantic -Wall -fdump-tree-optimized-blocks"
   NATIVE_CFLAG="-march=native"

   NEK5K_BIGMEM="no" # used by: package-builders/nek5000.sh
}


function setup_gcc_7()
{
   # Homebrew open-mpi with gcc v7
   CC=gcc-7
   CXX=g++-7
   FC=gfortran-7
   OCCA_CXX="$CXX"

   setup_openmpi

   CFLAGS="-O3"
   FFLAGS="$CFLAGS"
   OCCA_CXXFLAGS="-O3 -march=native"
   TEST_EXTRA_CFLAGS="-march=native --param max-completely-peel-times=3"
   # "-std=c++11 -pedantic -Wall -fdump-tree-optimized-blocks"
   NATIVE_CFLAG="-march=native"

   NEK5K_BIGMEM="no" # used by: package-builders/nek5000.sh
}


function setup_gcc_7a()
{
   # Homebrew open-mpi with gcc v7 without "--param max-completely-peel-times=3"
   setup_gcc_7
   TEST_EXTRA_CFLAGS="-march=native"
}


function setup_gcc_8()
{
   # Homebrew open-mpi with gcc v8
   CC=gcc-8
   CXX=g++-8
   FC=gfortran-8
   OCCA_CXX="$CXX"

   setup_openmpi

   CFLAGS="-O3"
   FFLAGS="$CFLAGS"
   OCCA_CXXFLAGS="-O3 -march=native"
   TEST_EXTRA_CFLAGS="-march=native --param max-completely-peel-times=3"
   # "-std=c++11 -pedantic -Wall -fdump-tree-optimized-blocks"
   NATIVE_CFLAG="-march=native"

   NEK5K_BIGMEM="no" # used by: package-builders/nek5000.sh
}


function set_mpi_options()
{
   num_proc_node=${num_proc_run}
   compose_mpi_run_command
}


# MFEM_EXTRA_CONFIG="MFEM_TIMER_TYPE=4"
GLVIS_EXTRA_CONFIG="GLVIS_MULTISAMPLE=4 GLVIS_MS_LINEWIDTH=0.01"
search_file_list LAPACK_LIB \
   "$(brew --prefix)/opt/openblas/lib/libopenblas.dylib" || \
LAPACK_LIB="-framework Accelerate"

valid_compilers="clang gcc_6 gcc_7 gcc_7a gcc_8"
num_proc_detect="$(getconf _NPROCESSORS_ONLN)"
num_proc_build=${num_proc_build:-$num_proc_detect}
num_proc_run=${num_proc_run:-2}
num_proc_node=${num_proc_run}
memory_per_node=8

# Optional (default): MPIEXEC (mpirun), MPIEXEC_OPTS (), MPIEXEC_NP (-np)
