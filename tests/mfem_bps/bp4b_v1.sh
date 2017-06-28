# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory. LLNL-CODE-XXXXXX. All Rights
# reserved. See file LICENSE for details.
#
# This file is part of CEED, a collection of benchmarks, miniapps, software
# libraries and APIs for efficient high-order finite element and spectral
# element discretizations for exascale applications. For more information and
# source code availability see http://github.com/ceed.
#
# The CEED research was supported by the Exascale Computing Project
# (17-SC-20-SC), a collaborative effort of two U.S. Department of Energy
# organizations (Office of Science and the National Nuclear Security
# Administration) responsible for the planning and preparation of a capable
# exascale ecosystem, including software, applications, hardware, advanced
# system engineering and early testbed platforms, in support of the nation's
# exascale computing imperative.


if [[ -z "$root_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   return 1
fi

# problem: 0 - diffusion, 1 - mass
problem=0
vdim=3
# Ordering::byVDIM or Ordering::byNODES
vec_layout="Ordering::byNODES"
source ${root_dir}/tests/mfem_bps/bp1_v1.sh
