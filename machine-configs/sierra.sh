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

function setup_xlc()
{
   MPICC=mpixlc
   MPICXX=mpixlC
   MPIF77=mpixlf
   mpi_info_flag="-qversion=verbose"

   CFLAGS="-O3 -qarch=auto -qcache=auto -qhot -qtune=auto"
   # CFLAGS+=" -qipa=threads"

   FFLAGS="$CFLAGS"

   TEST_EXTRA_CFLAGS="-O5"
   # TEST_EXTRA_CFLAGS+=" -std=c++11 -qreport"

   add_to_path after PATH "$cuda_path"
   CUFLAGS="-O3"
}


function setup_gcc()
{
   MPICC=mpigcc
   MPICXX=mpig++
   MPIF77=mpigfortran
   mpi_info_flag="--version"

   CFLAGS="-O3 -mcpu=native -mtune=native"
   FFLAGS="$CFLAGS"
   TEST_EXTRA_CFLAGS=""
   # TEST_EXTRA_CFLAGS+=" -std=c++11 -fdump-tree-optimized-blocks"

   add_to_path after PATH "$cuda_path"
   CUFLAGS="-O3"
}


function setup_clang()
{
   MPICC=mpiclang
   MPICXX=mpiclang++
   # MPIF77=?
   mpi_info_flag="--version"

   CFLAGS="-O3 -mcpu=native -mtune=native"
   # FFLAGS=?
   TEST_EXTRA_CFLAGS="-fcolor-diagnostics -fvectorize"
   TEST_EXTRA_CFLAGS+=" -fslp-vectorize -fslp-vectorize-aggressive"
   TEST_EXTRA_CFLAGS+=" -ffp-contract=fast"

   add_to_path after PATH "$cuda_path"
   CUFLAGS="-O3"
}


function set_mpi_options()
{
   local account="${account:-guests}"
   local partition="${partition:-pbatch}"
   # Time limit in minutes
   local TIME_LIMIT=30
   local BSUB_OPTS="-q ${partition} -G ${account} -nnodes ${num_nodes}"
   BSUB_OPTS="${BSUB_OPTS} -W ${TIME_LIMIT}"
   MPIEXEC_OPTS="-T $num_proc_node"
   MPIEXEC_OPTS="${BSUB_OPTS} lrun ${MPIEXEC_OPTS}"
   compose_mpi_run_command
}


valid_compilers="xlc gcc clang"
num_proc_build=${num_proc_build:-16}
num_proc_run=${num_proc_run:-4}
num_proc_node=${num_proc_node:-4}
memory_per_node=256

# Optional (default): MPIEXEC (mpirun), MPIEXEC_OPTS (), MPIEXEC_NP (-np)
bind_sh=
MPIEXEC="bsub"
MPIEXEC_NP="-n"

cuda_path=${CUDA_HOME:-/usr/local/cuda}/bin
