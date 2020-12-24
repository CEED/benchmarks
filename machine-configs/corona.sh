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

# Configuration for LLNL's Corona system, HSA Agents:
#     - 8 Agents: AMD EPYC 7401 24-Core Processor
#     - 4 Agents: AMD/ATI Vega 10 [Radeon Instinct MI25] (rev 01) - gfx900
#     - 4 Agents: AMD/ATI Vega 20 [Radeon Instinct MI60] (rev 01) - gfx906

hip_arch=${hip_arch:-gfx906}

function setup_hip()
{
    echo "${cyan}HIP setup${none}"
    module load opt
    module load rocm/3.1
    module load gcc/8.1.0
    module load mvapich2/2.3
    CXX=hipcc
    MPICC=mpicc
    MPICXX=mpicxx
    MPI_HOME=$(dirname $(which $MPICC))/..
    export MFEM_CPPFLAGS="-I${MPI_HOME}/include"
    export LDFLAGS="-L${MPI_HOME}/lib -lmpi"
    CFLAGS="-O3"
    CXX11FLAG="--std=c++11 -fno-gpu-rdc"
    hip_home=${HIP_HOME:-/opt/rocm/hip}
    hip_path=${hip_home}/bin
    hip_lib=${hip_home}/lib64
}

function setup_gcc8()
{
    # Load additional modules from 'opt'; needed for the 'rocm' module
    module load opt
    # ROCm 3.8.0 is using /opt/rh/devtoolset-8 for linking with libstdc++, so
    # we use the module gcc/8.3.1 which seems to be the same as the version in
    # /opt/rh/devtoolset-8
    module load gcc/8.3.1
    module load mvapich2/2.3

    # ROCm and HIP
    # module load rocm/3.8.0
    module load rocm/3.10.0

    MPICC=mpicc
    MPICXX=mpicxx
    CFLAGS="-O3"
    # CXX11FLAG="--std=c++11 -fno-gpu-rdc"
    # Used e.g. by libCEED:
    NATIVE_CFLAG="-march=native"

    hip_path=$(dirname $(which hipcc))
    hip_home=$(dirname ${hip_path})
    hip_lib=${hip_home}/lib64
}

function set_mpi_options()
{
    echo "${cyan}MPI setup${none}"
    num_proc_node=${num_proc_run}
    compose_mpi_run_command
}

valid_compilers="hip gcc8"

num_proc_detect="$(getconf _NPROCESSORS_ONLN)"
num_proc_build=${num_proc_build:-$num_proc_detect}
num_proc_run=${num_proc_run:-1}
memory_per_node=256

# Optional (default): MPIEXEC (mpirun), MPIEXEC_OPTS (), MPIEXEC_NP (-np)
MPIEXEC=srun
MPIEXEC_NP=-n
