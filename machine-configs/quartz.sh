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

# Configuration for LLNL's Quartz system

function setup_intel()
{
   module load intel/16.0.3
   module load mvapich2/2.2
   MPICC=mpicc
   MPICXX=mpicxx
   MPIF77=mpif77
   mpi_info_flag="-v"

   CFLAGS="-O3"
   FFLAGS="-O3"
   TEST_EXTRA_CFLAGS="-xHost"
   NATIVE_CFLAG="-xHost"

   NEK5K_EXTRA_PPLIST=""
}


function setup_gcc()
{
   module load gcc/7.1.0
   module load mvapich2/2.2
   MPICC=mpicc
   MPICXX=mpicxx
   MPIF77=mpif77

   CFLAGS="-O3"
   FFLAGS="$CFLAGS"
   TEST_EXTRA_CFLAGS="-march=native --param max-completely-peel-times=3"
   # TEST_EXTRA_CFLAGS+=" -std=c++11 -fdump-tree-optimized-blocks"
   NATIVE_CFLAG="-march=native"

   NEK5K_EXTRA_PPLIST=""
}


function setup_gcc_b()
{
   setup_gcc
   TEST_EXTRA_CFLAGS=""
}


function set_mpi_options()
{
   # Pools     Max nodes/job     Max runtime
   # pdebug             8(*)      30 minutes
   # pbatch          1200         24 hours
   #
   # (*) on a PER USER basis

   local account="${account:-ceed}"
   local partition="${partition:-pdebug}"
   MPIEXEC_OPTS="-A ${account} -p ${partition}"
   MPIEXEC_OPTS+=" --ntasks-per-node $num_proc_node"
   if [[ "$num_proc_node" -gt "36" ]]; then
      MPIEXEC_OPTS+=" --overcommit"
   fi
   compose_mpi_run_command
}


MFEM_EXTRA_CONFIG=""

valid_compilers="intel gcc gcc_b"
num_proc_build=${num_proc_build:-36}
num_proc_run=${num_proc_run:-36}
num_proc_node=${num_proc_node:-36}
memory_per_node=128
# node_virt_mem_lim=

# Optional (default): MPIEXEC (mpirun), MPIEXEC_OPTS (), MPIEXEC_NP (-np)
MPIEXEC=srun
MPIEXEC_NP=-n
