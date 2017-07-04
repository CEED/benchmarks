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

function setup_mpi()
{
   MPICC=mpicc
   MPICXX=mpicxx
   MPIF77=mpif77

   CFLAGS="-O3 -mcpu=a2 -mtune=a2"
   FFLAGS="$CFLAGS"
   TEST_EXTRA_CFLAGS=""
   NEK5K_EXTRA_PPLIST=""
}

function set_mpi_options()
{
   local account="${account:-CEED_ECPAD}"
   local queue="${queue:-default}"
   local mode="${mode:-c32}"
   local runtime="${runtime:-0:60:00}"

   MPIEXEC_OPTS="-I --mode ${mode} -A ${account} -q ${queue}"
   MPIEXEC_OPTS+=" -t ${runtime}"

   compose_mpi_run_command
}

MFEM_EXTRA_CONFIG="MFEM_TIMER_TYPE=0"

valid_compilers="gcc"
num_proc_build=${num_proc_build:-32}
num_proc_run=${num_proc_run:-32}
num_proc_node=${num_proc_node:-32}
memory_per_node=16
node_virt_mem_lim=16

# Optional (default): MPIEXEC (mpirun), MPIEXEC_OPTS (), MPIEXEC_NP (-np)
MPIEXEC=qsub
MPIEXEC_NP=-n
