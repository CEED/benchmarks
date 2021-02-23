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

function setup_gcc10()
{
    echo "${cyan}FUGAKU setup${none}"
    source /vol0001/apps/oss/spack/share/spack/setup-env.sh
    spack compiler find
    module load gcc-10.1.0-gcc-8.3.1-2wrm3no
    CXX=g++
    MPICC=gcc
    MPICXX=g++
    MPI_HOME=$(dirname $(which mpiFCC))/..
    INCFLAGS="-I${MPI_HOME}/include/mpi/fujitsu"
    LDFLAGS="-L${MPI_HOME}/lib64 -lmpi -lmpi_cxx"
    CFLAGS="-O3 -march=armv8.2-a+sve $INCFLAGS"
    CXX11FLAG="--std=c++11"
    export MFEM_CPPFLAGS=$INCFLAGS
}

function setup_fujitsu()
{
    echo "${cyan}FUGAKU setup${none}"
    CXX=FCC
    MPICC=mpifcc
    MPICXX=mpiFCC
    CFLAGS="-O3 -march=armv8.2-a+sve"
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

valid_compilers="gcc10 fujitsu"

num_proc_detect="$(getconf _NPROCESSORS_ONLN)"
num_proc_build=${num_proc_build:-$num_proc_detect}
num_proc_run=${num_proc_run:-1}
memory_per_node=256

# Optional (default): MPIEXEC (mpirun), MPIEXEC_OPTS (), MPIEXEC_NP (-np)
MPIEXEC=mpirun
MPIEXEC_NP=-n
