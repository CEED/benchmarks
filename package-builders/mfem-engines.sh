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

# Clone MFEM and build the parallel version.

if [[ -z "$pkg_sources_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   return 1
fi
if [[ -z "$OUT_DIR" ]]; then
   echo "The variable 'OUT_DIR' is not set. Stop."
   return 1
fi
pkg_src_dir="mfem"
MFEM_SOURCE_DIR="$pkg_sources_dir/$pkg_src_dir"
pkg_bld_dir="$OUT_DIR/mfem-engines"
MFEM_DIR="$pkg_bld_dir"
# The mfem_branch can be set at the command line of the go.sh call
mfem_branch="${mfem_branch:-engines-dev}"
MFEM_BRANCH="${mfem_branch}"
pkg_var_prefix="mfem_engines_"
pkg="MFEM engines (branch ${mfem_branch})"


function mfem_branch_clone()
{
   pkg_repo_list=("git@github.com:mfem/mfem.git"
                  "https://github.com/mfem/mfem.git")
   pkg_git_branch="${mfem_branch}"
   cd "$pkg_sources_dir" || return 1
   if [[ -d "$pkg_src_dir" ]]; then
      update_git_package
      return
   fi
   for pkg_repo in "${pkg_repo_list[@]}"; do
      echo "Cloning $pkg from $pkg_repo ..."
      git clone "$pkg_repo" "$pkg_src_dir" && return 0
   done
   echo "Could not successfully clone $pkg. Stop."
   return 1
}


function mfem_branch_build()
{
   if package_build_is_good; then
      echo "Using successfully built $pkg from OUT_DIR."
      return 0
   elif [[ ! -d "$pkg_bld_dir" ]]; then
      mkdir -p "$pkg_bld_dir"
   fi
   if [[ -z "$HYPRE_DIR" ]]; then
      echo "The required variable 'HYPRE_DIR' is not set. Stop."
      return 1
   fi
   if [[ -z "$METIS_DIR" ]]; then
      echo "The required variable 'METIS_DIR' is not set. Stop."
      return 1
   fi
   if [[ -z "$OCCA_DIR" ]]; then
      echo "The required variable 'OCCA_DIR' is not set. Stop."
      return 1
   fi
   local METIS_5="NO"
   [[ "$METIS_VERSION" = "5" ]] && METIS_5="YES"
   local SUNDIALS_MAKE_OPTS=()
   if [[ -n "$SUNDIALS_DIR" ]]; then
      SUNDIALS_MAKE_OPTS=(
         "MFEM_USE_SUNDIALS=YES"
         "SUNDIALS_DIR=$SUNDIALS_DIR")
   else
      echo "${magenta}Warning: Building $pkg without SUNDIALS ...${none}"
   fi
   local GECKO_MAKE_OPTS=()
   if [[ -n "$GECKO_DIR" ]]; then
      GECKO_MAKE_OPTS=(
         "MFEM_USE_GECKO=YES"
         "GECKO_DIR=$GECKO_DIR")
   else
      echo "${magenta}Warning: Building $pkg without Gecko ...${none}"
   fi
   echo "Building $pkg, sending output to ${pkg_bld_dir}_build.log ..." && {
      local num_nodes=1  # for 'make check' or 'make test'
      set_mpi_options    # for 'make check' or 'make test'
      cd "$pkg_bld_dir" && \
      make config \
         -f "$MFEM_SOURCE_DIR/makefile" \
         MFEM_USE_MPI=YES \
         MFEM_USE_BACKENDS=YES \
         MFEM_USE_OCCA=YES \
         $MFEM_EXTRA_CONFIG \
         MPICXX="$MPICXX" \
         CXXFLAGS="$CFLAGS" \
         HYPRE_DIR="$HYPRE_DIR/src/hypre" \
         METIS_DIR="$METIS_DIR" \
         OCCA_DIR="$OCCA_DIR" \
         MFEM_USE_METIS_5="$METIS_5" \
         "${SUNDIALS_MAKE_OPTS[@]}" \
         "${GECKO_MAKE_OPTS[@]}" \
         LDFLAGS="${LDFLAGS[*]}" \
         MFEM_MPIEXEC="${MPIEXEC:-mpirun}" \
         MFEM_MPIEXEC_NP="${MPIEXEC_OPTS} ${MPIEXEC_NP:--np}" && \
      make -j $num_proc_build
   } &> "${pkg_bld_dir}_build.log" || {
      echo " ... building $pkg FAILED, see log for details."
      return 1
   }
   echo "Build successful."
   print_variables "$pkg_var_prefix" \
      HYPRE_DIR METIS_DIR METIS_VERSION SUNDIALS_DIR GECKO_DIR OCCA_DIR \
      MFEM_BRANCH \
      > "${pkg_bld_dir}_build_successful"
}


function build_package()
{
   mfem_branch_clone && get_package_git_version && mfem_branch_build
}
