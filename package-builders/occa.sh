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

# Clone and build OCCA.

if [[ -z "$pkg_sources_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   return 1
fi
if [[ -z "$OUT_DIR" ]]; then
   echo "The variable 'OUT_DIR' is not set. Stop."
   return 1
fi
pkg_src_dir="occa"
OCCA_SOURCE_DIR="$pkg_sources_dir/$pkg_src_dir"
pkg_bld_dir="$OUT_DIR/occa"
OCCA_DIR="$pkg_bld_dir"
# The OCCA branch can be set on the command line of go.sh too.
occa_branch="${occa_branch:-master}"
OCCA_BRANCH="${occa_branch}"
pkg="OCCA (branch ${occa_branch})"
pkg_var_prefix="occa_"


function occa_clone()
{
   pkg_repo_list=("git@github.com:libocca/occa.git"
                  "https://github.com/libocca/occa.git")
   pkg_git_branch="${occa_branch}"
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


function occa_build()
{
   if package_build_is_good; then
      echo "Using successfully built $pkg from OUT_DIR."
      return 0
   elif [[ ! -d "$pkg_bld_dir" ]]; then
      cd "$OUT_DIR" && git clone "$OCCA_SOURCE_DIR" || {
         echo "Cloning $OCCA_SOURCE_DIR to OUT_DIR failed. Stop."
         return 1
      }
   fi
   echo "Building $pkg, sending output to ${pkg_bld_dir}_build.log ..." && {
      cd "$pkg_bld_dir" && \
      make \
         CXX="$MPICXX" \
         CXXFLAGS="$CFLAGS" \
         $OCCA_EXTRA_CONFIG \
         -j $num_proc_build
   } &> "${pkg_bld_dir}_build.log" || {
      echo " ... building $pkg FAILED, see log for details."
      return 1
   }
   echo "Build successful."
   print_variables "$pkg_var_prefix" \
      OCCA_BRANCH \
      > "${pkg_bld_dir}_build_successful"
}


function build_package()
{
   occa_clone && get_package_git_version && occa_build || return 1

   add_to_path before PATH "${OCCA_DIR}/bin"
   add_to_path before LD_LIBRARY_PATH "${OCCA_DIR}/lib"
   add_to_path before DYLD_LIBRARY_PATH "${OCCA_DIR}/lib"
   OCCA_CACHE_DIR="${OCCA_DIR}/cache"
   mkdir -p "${OCCA_CACHE_DIR}"
   OCCA_CXX="${OCCA_CXX:-$MPICXX}"
   OCCA_CXXFLAGS="${OCCA_CXXFLAGS:-$CFLAGS}"
   OCCA_CUDA_COMPILER_FLAGS="${CUFLAGS}"
   export PATH LD_LIBRARY_PATH DYLD_LIBRARY_PATH
   export OCCA_DIR OCCA_CACHE_DIR OCCA_CXX OCCA_CXXFLAGS
   export OCCA_CUDA_COMPILER_FLAGS
   # If OCCA_VERBOSE is set (e.g. on the go.sh command line) export it:
   [[ -n "$OCCA_VERBOSE" ]] && export OCCA_VERBOSE
   return 0
}
