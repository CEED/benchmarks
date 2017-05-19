# This file is part of CEED. For more details, see exascaleproject.org.


function setup_clang()
{
   # Homebrew open-mpi; default is to use apple clang
   mpi_cc=mpicc
   mpi_cxx=mpicxx

   CFLAGS="-O3"
   TEST_EXTRA_CFLAGS="-march=native -fcolor-diagnostics -fvectorize"
   TEST_EXTRA_CFLAGS+=" -fslp-vectorize -fslp-vectorize-aggressive"
   TEST_EXTRA_CFLAGS+=" -ffp-contract=fast"
   # "-std=c++11 -pedantic -Wall"
}


function setup_gcc_6()
{
   # Homebrew open-mpi with gcc v6
   mpi_cc=mpicc
   mpi_cxx=mpicxx
   export OMPI_CC=gcc-6
   export OMPI_CXX=g++-6

   CFLAGS="-O3"
   TEST_EXTRA_CFLAGS="-march=native --param max-completely-peel-times=3"
   # "-std=c++11 -pedantic -Wall -fdump-tree-optimized-blocks"
}


function set_mpi_options()
{
   num_proc_node=${num_proc_run}
}


GLVIS_EXTRA_CONFIG="GLVIS_MULTISAMPLE=4 GLVIS_MS_LINEWIDTH=0.01"

valid_compilers="clang gcc_6"
num_proc_build=4
num_proc_run=${num_proc_run:-2}
num_proc_node=${num_proc_run}
memory_per_node=8

# Optional (default): MPIEXEC (mpirun), MPIEXEC_OPTS (), MPIEXEC_NP (-np)
