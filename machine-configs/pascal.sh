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

# Configuration for LLNL's Pascal system

function setup_intel()
{
   module load intel/16.0.3
   module load mvapich2/2.2

   CC=icc
   CXX=icpc

   MPICC=mpicc
   MPICXX=mpicxx
   MPIF77=mpif77
   mpi_info_flag="-v"

   CFLAGS="-O3"
   FFLAGS="-O3"
   TEST_EXTRA_CFLAGS="-xHost"
   NATIVE_CFLAG="-xHost"

   NEK5K_EXTRA_PPLIST=""

   module load cuda/9.2.88
   cuda_home=${CUDA_HOME:-/usr/local/cuda}
   cuda_path=${cuda_home}/bin
   CUFLAGS="-O3"
}


function setup_gcc()
{
   module load gcc/7.3.0
   module load mvapich2/2.2

   CC=gcc
   CXX=g++

   MPICC=mpicc
   MPICXX=mpicxx
   MPIF77=mpif77

   CFLAGS="-O3"
   FFLAGS="$CFLAGS"
   TEST_EXTRA_CFLAGS="-march=native --param max-completely-peel-times=3"
   # TEST_EXTRA_CFLAGS+=" -std=c++11 -fdump-tree-optimized-blocks"
   NATIVE_CFLAG="-march=native"

   NEK5K_EXTRA_PPLIST=""

   module load cuda/9.2.88
   cuda_home=${CUDA_HOME:-/usr/local/cuda}
   cuda_path=${cuda_home}/bin
   CUFLAGS="-O3"
}


function set_mpi_options()
{
   # Pools     Max nodes/job     Max runtime
   # pbatch          32          24 hours

   if [[ -z "${SLURM_JOB_CPUS_PER_NODE}" ]]; then
      local account="${account:-ceed}"
      local partition="${partition:-pbatch}"
      MPIEXEC_OPTS="-A ${account} -p ${partition}"
      MPIEXEC_OPTS+=" --ntasks-per-node $num_proc_node"
      if [[ "$num_proc_node" -gt "36" ]]; then
         MPIEXEC_OPTS+=" --overcommit"
      fi
   else
      local job_num_nodes=$SLURM_NNODES
      if (( job_num_nodes < num_nodes )); then
         echo "Insufficient number of nodes in the job allocation:"
         echo "   ($job_num_nodes < $num_nodes)"
         exit 1
      fi
      MPIEXEC_OPTS="--ntasks-per-node $num_proc_node"
      if [[ "$num_proc_node" -gt "36" ]]; then
         MPIEXEC_OPTS+=" --overcommit"
      fi
   fi
   compose_mpi_run_command
}


MFEM_EXTRA_CONFIG=""

valid_compilers="intel gcc"
num_proc_build=${num_proc_build:-36}
num_proc_run=${num_proc_run:-36}
num_proc_node=${num_proc_node:-36}
memory_per_node=256
# node_virt_mem_lim=

# Optional (default): MPIEXEC (mpirun), MPIEXEC_OPTS (), MPIEXEC_NP (-np)
MPIEXEC=srun
MPIEXEC_NP=-n
