# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory. LLNL-CODE-XXXXXX. All Rights
# reserved. See file LICENSE for details.
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
# testbed platforms, in support of the nation’s exascale computing imperative.

# Clone and build Acrotensor.

if [[ -z "$pkg_sources_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   return 1
fi
if [[ -z "$OUT_DIR" ]]; then
   echo "The variable 'OUT_DIR' is not set. Stop."
   return 1
fi
pkg_src_dir="acrotensor"
ACROTENSOR_SOURCE_DIR="$pkg_sources_dir/$pkg_src_dir"
pkg_bld_dir="$OUT_DIR/acrotensor"
ACROTENSOR_DIR="$pkg_bld_dir"
CUDA_DIR="/usr/local/cuda"
# TODO How should compiler/option dependencies be handled?
ACROTENSOR_CXXFLAGS="-x cu --std=c++11 -DACRO_HAVE_CUDA --compiler-options=\"-fPIC\""
pkg="acrotensor"
llnluser=$(whoami)

function acrotensor_clone()
{
   pkg_repo_list=(
      "https://${llnluser}@lc.llnl.gov/bitbucket/scm/~fisher47/acrotensor.git"
      "ssh://git@cz-bitbucket.llnl.gov:7999/~fisher47/acrotensor.git")
   pkg_git_branch="master"
   cd "$pkg_sources_dir" || return 1
   if [[ -d "$pkg_src_dir" ]]; then
      update_git_package
      return
   fi
   for pkg_repo in "${pkg_repo_list[@]}"; do
      echo "Cloning $pkg from $pkg_repo ..."
      git clone "$pkg_repo" "$pkg_src_dir" && \
      cd "$pkg_src_dir" && \
      git checkout "$pkg_git_branch" && \
      return 0
   done
   echo "Could not successfully clone $pkg. Stop."
   return 1
}


function acrotensor_build()
{
   if package_build_is_good; then
      echo "Using successfully built $pkg from OUT_DIR."
      return 0
   elif [[ ! -d "$pkg_bld_dir" ]]; then
      cd "$OUT_DIR" && git clone "$ACROTENSOR_SOURCE_DIR" || {
         echo "Cloning $ACROTENSOR_SOURCE_DIR to OUT_DIR failed. Stop."
         return 1
      }
   fi
   echo "Building $pkg, sending output to ${pkg_bld_dir}_build.log ..." && {
      cd "$pkg_bld_dir" && \
      make config && \
      make \
         CXX=nvcc \
         DEBUG=NO \
         UTILCXX=nvcc \
         CXX_OPT="$ACROTENSOR_CXXFLAGS" \
         CXX_DEBUG="$ACROTENSOR_CXXFLAGS -g" \
         $ACROTENSOR_EXTRA_CONFIG \
         -j $num_proc_build
   } &> "${pkg_bld_dir}_build.log" || {
      echo " ... building $pkg FAILED, see log for details."
      return 1
   }
   echo "Build succesful."
   : > "${pkg_bld_dir}_build_successful"
}


function build_package()
{
   acrotensor_clone && get_package_git_version && acrotensor_build
}
