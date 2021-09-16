# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
0;95;0c0;95;0c# the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
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

function setup_mix()
{
    echo "${cyan}FUGAKU MIX setup${none}"
    CXX_HOME=$(dirname $(which g++))/..
    CXX=g++
    MPICC=gcc
    MPICXX=g++
    MPI_HOME=$(dirname $(which mpiFCC)/..)
    INCFLAGS=-I${MPI_HOME}/include/mpi/fujitsu
    LDFLAGS="-L${MPI_HOME}/lib64 -lmpi $CXX_HOME/lib64/libstdc++.a"
    CFLAGS="-O3 -march=native -mtune=native $INCFLAGS"
    CXX11FLAG="--std=c++11"
    export MFEM_CPPFLAGS=$INCFLAGS
    TEST_EXTRA_CFLAGS="-ffast-math -msve-vector-bits=512 --param max-completely-peel-times=3"
}

function setup_gcc11()
{
    echo "${cyan}FUGAKU GCC11 setup${none}"
    CXX=/home/ra010009/a04177/usr/local/gcc/11.2.0/bin/g++
    MPICC=/home/ra010009/a04177/usr/local/gcc/11.2.0/bin/gcc
    MPICXX=/home/ra010009/a04177/usr/local/gcc/11.2.0/bin/g++
    MPI_HOME=/home/ra010009/a04177/usr/local/openmpi-4.1.1
    INCFLAGS=-I${MPI_HOME}/include
    LDFLAGS="-L${MPI_HOME}/lib -lmpi /home/ra010009/a04177/usr/local/gcc/11.2.0/lib64/libstdc++.a"
    CFLAGS="-O3 -ffast-math -march=native -msve-vector-bits=512 $INCFLAGS"
    CXX11FLAG="--std=c++11"
    export MFEM_CPPFLAGS=$INCFLAGS
    TEST_EXTRA_CFLAGS="--param max-completely-peel-times=8"
}

function setup_gcc10()
{
    echo "${cyan}FUGAKU GCC10 setup${none}"
    CXX=/home/ra010009/a04177/usr/local/gcc/10.2.0/bin/g++
    MPICC=/home/ra010009/a04177/usr/local/gcc/10.2.0/bin/gcc
    MPICXX=/home/ra010009/a04177/usr/local/gcc/10.2.0/bin/g++
    MPI_HOME=/home/ra010009/a04177/usr/local/openmpi/4.1.0
    INCFLAGS=-I${MPI_HOME}/include
    LDFLAGS="-L${MPI_HOME}/lib -lmpi /home/ra010009/a04177/usr/local/gcc/10.2.0/lib64/libstdc++.a"
    CFLAGS="-O3 -march=armv8.2-a+sve -msve-vector-bits=512 $INCFLAGS"
    CXX11FLAG="--std=c++11"
    export MFEM_CPPFLAGS=$INCFLAGS
    TEST_EXTRA_CFLAGS="--param max-completely-peel-times=8"
}


function setup_fujitsu()
{
    echo "${cyan}FUGAKU setup${none}"
    CXX=FCC
    MPICC=mpifcc
    MPICXX=mpiFCC
    CFLAGS="-Kfast,simd_reg_size=512"
    CXX11FLAG="--std=c++11"
}

function set_mpi_options()
{
    echo "${cyan}MPI setup${none}"
    num_proc_node=${num_proc_run}
    #mpi_run="pjsub -j -L node=1 -L elapse=2:00 --mpi proc=${num_proc_run} "
    #mpi_run+="--mpi max-proc-per-node=48 "
    compose_mpi_run_command
}

valid_compilers="mix gcc10 gcc11 fujitsu"

num_proc_detect="$(getconf _NPROCESSORS_ONLN)"
num_proc_build=${num_proc_build:-$num_proc_detect}
num_proc_run=${num_proc_run:-1}
memory_per_node=256

# Optional (default): MPIEXEC (mpirun), MPIEXEC_OPTS (), MPIEXEC_NP (-np)
MPIEXEC=mpirun
MPIEXEC_NP=-n
