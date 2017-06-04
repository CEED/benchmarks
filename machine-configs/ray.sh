# This file is part of CEED. For more details, see exascaleproject.org.


function setup_xlc()
{
   MPICC=mpixlc
   MPICXX=mpixlC
   mpi_info_flag="-qversion=verbose"

   CFLAGS="-O3 -qarch=auto -qcache=auto -qhot -qtune=auto"
   TEST_EXTRA_CFLAGS="-O5"
   # TEST_EXTRA_CFLAGS+=" -std=c++11 -qreport"

   add_to_path PATH "$cuda_path"
}


function setup_gcc()
{
   # GCC 4.9.3 is the default GCC as of 2017-05-28.

   MPICC=mpigcc
   MPICXX=mpig++
   mpi_info_flag="--version"

   CFLAGS="-O3 -mcpu=native -mtune=native"
   TEST_EXTRA_CFLAGS=""
   # TEST_EXTRA_CFLAGS+=" -std=c++11 -fdump-tree-optimized-blocks"

   add_to_path PATH "$cuda_path"
}


function setup_clang()
{
   # Clang coral compilers: clang-coral-2017.05.19 as of 2017-05-29.

   MPICC=mpiclang
   MPICXX=mpiclang++
   mpi_info_flag="--version"

   CFLAGS="-O3 -mcpu=native -mtune=native"
   TEST_EXTRA_CFLAGS="-fcolor-diagnostics -fvectorize"
   TEST_EXTRA_CFLAGS+=" -fslp-vectorize -fslp-vectorize-aggressive"
   TEST_EXTRA_CFLAGS+=" -ffp-contract=fast"

   add_to_path PATH "$cuda_path"
}


function set_mpi_options()
{
   MPIEXEC_OPTS="-npernode $num_proc_node"
   # Q: Is bind.sh equivalent to the mpirun option '-bind-to core'?
   # MPIEXEC_OPTS="-bind-to core $MPIEXEC_OPTS"
   # MPIEXEC_OPTS+=" -report-bindings"
   # echo "$LSB_HOSTS" > hostfile
   # MPIEXEC_OPTS+=" -hostfile hostfile"
   if [ "$num_proc_node" -gt "20" ]; then
      MPIEXEC_OPTS="-oversubscribe $MPIEXEC_OPTS"
   fi
   local BSUB_OPTS="-q pbatch -G guests -n $((num_nodes*20)) -I"
   MPIEXEC_OPTS="$BSUB_OPTS mpirun $MPIEXEC_OPTS"
   compose_mpi_run_command
}


valid_compilers="xlc gcc clang"
num_proc_build=10
num_proc_run=${num_proc_run:-20}
num_proc_node=${num_proc_node:-20}
memory_per_node=256

# Optional (default): MPIEXEC (mpirun), MPIEXEC_OPTS (), MPIEXEC_NP (-np)
bind_sh=mpibind
# bind_sh=bind.sh
# MPIEXEC=jsrun (under dev)
MPIEXEC="bsub"
MPIEXEC_NP="-n"

cuda_path=/usr/local/cuda-8.0/bin