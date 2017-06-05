# This file is part of CEED. For more details, see exascaleproject.org.

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
pkg="OCCA"


function occa_clone()
{
   pkg_repo_list=("git@github.com:libocca/occa.git"
                  "https://github.com/libocca/occa.git")
   pkg_git_branch="1.0"
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
   : > "${pkg_bld_dir}_build_successful"
}


function build_package()
{
   occa_clone && occa_build || return 1

   add_to_path PATH "${OCCA_DIR}/bin"
   add_to_path LD_LIBRARY_PATH "${OCCA_DIR}/lib"
   add_to_path DYLD_LIBRARY_PATH "${OCCA_DIR}/lib"
   OCCA_CXX="${OCCA_CXX:-$MPICXX}"
   OCCA_CXXFLAGS="${OCCA_CXXFLAGS:-$CFLAGS}"
   OCCA_CUDA_COMPILER_FLAGS="${CUFLAGS}"
   export PATH LD_LIBRARY_PATH DYLD_LIBRARY_PATH OCCA_CXX OCCA_CXXFLAGS
   export OCCA_CUDA_COMPILER_FLAGS
}
