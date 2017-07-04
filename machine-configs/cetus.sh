# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory. LLNL-CODE-XXXXXX. All Rights
# reserved. See file LICENSE for details.
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

function setup_gcc()
{
   # GCC 4.7.2
   MPICC=mpigcc-4.7.2-fastmpi
   MPICXX=mpig++-4.7.2-fastmpi
   MPIF77=mpigfortran-4.7.2-fastmpi

   CFLAGS="-O3 -mcpu=a2 -mtune=a2"
   FFLAGS="$CFLAGS"
   TEST_EXTRA_CFLAGS=""

   NEK5K_EXTRA_PPLIST=""
}

MFEM_EXTRA_CONFIG="MFEM_TIMER_TYPE=0"

valid_compilers="gcc"
num_proc_build=${num_proc_build:-16}
num_proc_run=${num_proc_run:-16}
num_proc_node=${num_proc_node:-16}
memory_per_node=16
node_virt_mem_lim=16

# Optional (default): MPIEXEC (mpirun), MPIEXEC_OPTS (), MPIEXEC_NP (-np)
MPIEXEC=qsub
MPIEXEC_NP=-n
