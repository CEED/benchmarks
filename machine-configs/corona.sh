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

function setup_hip()
{
    module load opt
    module load rocm/2.7
   
    
    CXX=hipcc

   MPICC=mpicc
   MPICXX=mpicxx
   MPIF77=mpif77

   CFLAGS="-O3"
   FFLAGS="$CFLAGS"
   
   hip_home=${HIP_HOME:-/opt/rocm/hip}
   # /opt/rocm/lib64 /opt/rocm/hsa/lib
   hip_path=${hip_home}/bin
   hsa_lib=${hip_home}/../hsa/lib
   hip_lib=${hip_home}/lib64
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

valid_compilers="hip"
num_proc_build=${num_proc_build:-48}
num_proc_run=${num_proc_run:-48}
num_proc_node=${num_proc_node:-48}
memory_per_node=256
# node_virt_mem_lim=

# Optional (default): MPIEXEC (mpirun), MPIEXEC_OPTS (), MPIEXEC_NP (-np)
MPIEXEC=srun
MPIEXEC_NP=-n

hip_arch=gfx900
