# This file is part of CEED. For more details, see exascaleproject.org.


function setup_xlc()
{
   MPICC=mpixlc_r-fastmpi
   MPICXX=mpixlcxx_r-fastmpi
   mpi_info_flag="-qversion=verbose"

   # CFLAGS=""
   # CFLAGS="-O2 -qmaxmem=-1"
   # CFLAGS="-O3 -qstrict -qdebug=forcesqrt -qdebug=nsel -qmaxmem=-1"
   CFLAGS="-O3 -qnounwind -qsuppress=1540-1088:1540-1090:1540-1101"
   # TEST_EXTRA_CFLAGS=""
   TEST_EXTRA_CFLAGS="-O5 -qnounwind -qstrict"
   TEST_EXTRA_CFLAGS+=" -qsuppress=1540-1088:1540-1090:1540-1101"
   # TEST_EXTRA_CFLAGS+=" -qnoeh"
   # TEST_EXTRA_CFLAGS+=" -qreport -qlistopt -qlist -qskipsrc=hide -qsource"
}


function setup_gcc()
{
   # GCC 4.7.2
   MPICC=mpigcc-4.7.2-fastmpi
   MPICXX=mpig++-4.7.2-fastmpi

   CFLAGS="-O3 -mcpu=a2 -mtune=a2"
   TEST_EXTRA_CFLAGS=""
   # TEST_EXTRA_CFLAGS+=" -std=c++11 -fdump-tree-optimized-blocks"
}


function setup_gcc_b()
{
   # GCC 4.7.2
   MPICC=mpigcc-4.7.2-fastmpi
   MPICXX=mpig++-4.7.2-fastmpi

   CFLAGS="-O3 -mcpu=a2 -mtune=a2"
   TEST_EXTRA_CFLAGS="--param max-completely-peel-times=3"
   # TEST_EXTRA_CFLAGS+=" -std=c++11 -fdump-tree-optimized-blocks"
}


function setup_clang()
{
   # clang 3.7.0
   MPICC=mpiclang-fastmpi
   MPICXX=mpiclang++-fastmpi

   CFLAGS="-O3 -mcpu=a2 -mtune=a2"
   TEST_EXTRA_CFLAGS="-fcolor-diagnostics -fvectorize -fslp-vectorize"
   TEST_EXTRA_CFLAGS+=" -fslp-vectorize-aggressive -ffp-contract=fast"
   # "-std=c++11 -pedantic -Wall"
}


function set_mpi_options()
{
   # pdebug (<= 1K nodes & <= 1h)
   # psmall (<= 1K nodes & <= 12h)
   # pbatch (> 1K nodes & <= 8K nodes & <= 12h)
   MPIEXEC_OPTS="-A ceed -p pdebug"
   # MPIEXEC_OPTS="-A ceed -p psmall"
   MPIEXEC_OPTS+=" --ntasks-per-node $num_proc_node"
   if [[ "$num_proc_node" -gt "16" ]]; then
      MPIEXEC_OPTS+=" --overcommit"
   fi
   compose_mpi_run_command
}


MFEM_EXTRA_CONFIG="MFEM_TIMER_TYPE=0"

valid_compilers="xlc gcc gcc_b clang"
num_proc_build=16
num_proc_run=${num_proc_run:-16}
num_proc_node=${num_proc_node:-16}
memory_per_node=16

# Optional (default): MPIEXEC (mpirun), MPIEXEC_OPTS (), MPIEXEC_NP (-np)
MPIEXEC=srun
MPIEXEC_NP=-n
