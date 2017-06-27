# This file is part of CEED. For more details, see exascaleproject.org.

function setup_bigmem()
{
   ARCH=$(uname -m)
   if [[ "$ARCH" == "x86_64" ]]; then
     CFLAGS+=" -mcmodel=medium"
   elif [[ "$ARCH" == "ppc64" ]]; then
     CFLAGS+=" -m64"
   fi
}

function setup_intelmpi()
{
   MPICC=mpiicc
   MPICXX=mpiicpc
   MPIFC=mpiifort
   MPIF77=mpiifort
   # OpenMPI
   export OMPI_CC="$CC"
   export OMPI_CXX="$CXX"
   export OMPI_FC="$FC"
   # or MPICH
   export MPICH_CC="$CC"
   export MPICH_CXX="$CXX"
   export MPICH_FC="$FC"
}

function setup_mpi()
{
   MPICC=mpicc
   MPICXX=mpicxx
   MPIFC=mpifort
   MPIF77=mpif77
   # OpenMPI
   export OMPI_CC="$CC"
   export OMPI_CXX="$CXX"
   export OMPI_FC="$FC"
   # or MPICH
   export MPICH_CC="$CC"
   export MPICH_CXX="$CXX"
   export MPICH_FC="$FC"
}

function setup_intel()
{
   CC=icc
   CXX=icpc
   FC=ifort

   setup_intelmpi

   CFLAGS="-O3"
   setup_bigmem
   FFLAGS="$CFLAGS"

   # The following options assume GCC:
   TEST_EXTRA_CFLAGS="-march=native --param max-completely-peel-times=3"
   # "-std=c++11 -pedantic -Wall -fdump-tree-optimized-blocks"
}

function setup_gcc()
{
   CC=gcc
   CXX=g++
   FC=gfortran

   setup_mpi

   CFLAGS="-O3"
   setup_bigmem
   FFLAGS="$CFLAGS"

   # The following options assume GCC:
   TEST_EXTRA_CFLAGS="-march=native --param max-completely-peel-times=3"
   # "-std=c++11 -pedantic -Wall -fdump-tree-optimized-blocks"
}


function setup_clang()
{
   CC=clang
   CXX=clang++
   # Use gfortran
   FC=gfortran

   setup_mpi

   CFLAGS="-O3"
   setup_bigmem
   FFLAGS="$CFLAGS"

   TEST_EXTRA_CFLAGS="-march=native -fcolor-diagnostics -fvectorize"
   TEST_EXTRA_CFLAGS+=" -fslp-vectorize -fslp-vectorize-aggressive"
   TEST_EXTRA_CFLAGS+=" -ffp-contract=fast"
}


function set_mpi_options()
{
   # Run all tasks on the same node.
   num_proc_node=${num_proc_run}
   compose_mpi_run_command
}


search_file_list LAPACK_LIB \
   "/usr/lib64/atlas/libsatlas.so" "/usr/lib64/libopenblas.so"

valid_compilers="gcc clang intel"
# Number of processors to use for building packages and tests:
num_proc_build=8
# Default number of processors and processors per node for running tests:
num_proc_run=${num_proc_run:-4}
num_proc_node=${num_proc_run}
# Total memory per node:
memory_per_node=8

# Optional (default): MPIEXEC (mpirun), MPIEXEC_OPTS (), MPIEXEC_NP (-np)
