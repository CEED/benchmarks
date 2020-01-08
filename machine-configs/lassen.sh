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

   CFLAGS="-O2 -qmaxmem=-1 -qarch=auto -qcache=auto -qtune=auto"
   # CFLAGS="-O3 -qarch=auto -qcache=auto -qhot -qtune=auto"
   # CFLAGS+=" -qipa=threads"

   FFLAGS="$CFLAGS"

   TEST_EXTRA_CFLAGS="-O5"
   # TEST_EXTRA_CFLAGS+=" -std=c++11 -qreport"

   add_to_path after PATH "$cuda_path"
   CUFLAGS="-O3"

   # OpenMP
   omp_flag="-qsmp=omp"
}


function setup_gcc()
{
   module load gcc/7.3.1

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

   # OpenMP
   omp_flag="-fopenmp"
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
   if [[ -z "$gpu_direct" ]]; then
      LRUN_GPU_OPT=
   else
      LRUN_GPU_OPT="-M -gpu"
   fi

   # Autodetect if running inside a job
   if [[ -z "${LSB_JOBID}" ]]; then
      local account="${account:-ceed}"
      local partition="${partition:-pbatch}"
      # Time limit in minutes
      local TIME_LIMIT=${time_limit:-30}
      local BSUB_OPTS="-I -q ${partition} -G ${account} -nnodes ${num_nodes}"
      BSUB_OPTS="${BSUB_OPTS} -W ${TIME_LIMIT}"
      MPIEXEC_OPTS="-N $num_nodes"
      MPIEXEC_OPTS="${BSUB_OPTS} lrun ${MPIEXEC_OPTS}"
      MPIEXEC_POST_OPTS="${LRUN_GPU_OPT}"
      MPIEXEC="bsub"
      MPIEXEC_NP="-n"
   else
      # LSB_DJOB_NUMPROC=num
      #   - The number of processors (slots) allocated to the job.
      # LSB_MCPU_HOSTS="hostA num_processors1 hostB num_processors2..."
      # On lassen, the LSB_MCPU_HOSTS seem to include the "launch" node with
      # 1 processor. Also, LSB_DJOB_NUMPROC includes the 1 processor of the
      # "launch" node.
      local job_nodes_list=($LSB_MCPU_HOSTS)
      local job_num_nodes=$(( ${#job_nodes_list[@]} / 2 - 1 ))
      if (( job_num_nodes < num_nodes )); then
         echo "Insufficient number of nodes in the job allocation:"
         echo "   ($job_num_nodes < $num_nodes)"
         exit 1
      fi
      # Note: The following check is disabled to allow oversubscription.
      # if (( LSB_DJOB_NUMPROC-1 < num_proc_run )); then
      #    echo "Insufficient number of processors in the job allocation:"
      #    echo "   ($LSB_DJOB_NUMPROC - 1 < $num_proc_run)"
      #    exit 1
      # fi
      MPIEXEC_OPTS="-N $num_nodes"
      MPIEXEC_POST_OPTS="${LRUN_GPU_OPT}"
      MPIEXEC="lrun"
      MPIEXEC_NP="-n"
   fi
   compose_mpi_run_command
}


valid_compilers="xlc gcc clang"
num_proc_build=${num_proc_build:-16}
num_proc_run=${num_proc_run:-4}
num_proc_node=${num_proc_node:-4}
memory_per_node=256

# Optional (default): MPIEXEC (mpirun), MPIEXEC_OPTS (), MPIEXEC_NP (-np)
bind_sh=
MPIEXEC=
MPIEXEC_NP=

cuda_home=${CUDA_HOME:-/usr/local/cuda}
cuda_path=${cuda_home}/bin
cuda_arch=sm_70
