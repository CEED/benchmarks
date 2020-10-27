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

# Pseudo-package that allows CUDA support to be enabled/disabled for other
# packages that can use CUDA.

if [[ -z "$pkg_sources_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   return 1
fi
if [[ -z "$OUT_DIR" ]]; then
   echo "The variable 'OUT_DIR' is not set. Stop."
   return 1
fi

pkg_src_dir="cuda"
pkg_bld_dir="$OUT_DIR/cuda"
if [[ -n "$cuda_home" ]]; then
   CUDA_ENABLED="YES"
else
   CUDA_ENABLED=""
fi
pkg="CUDA"


function build_package()
{
   if [[ -n "$CUDA_ENABLED" ]]; then
      echo "CUDA is enabled. [ cuda_home = $cuda_home ]"
   else
      echo "Error: cannot enable CUDA: 'cuda_home' is not configured. Stop."
      return 1
   fi
}
