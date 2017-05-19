# This file is part of CEED. For more details, see exascaleproject.org.

# Clone and build hypre.

if [[ -z "$pkg_sources_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   exit 1
fi
if [[ -z "$OUT_DIR" ]]; then
   echo "The variable 'OUT_DIR' is not set. Stop."
   exit 1
fi
pkg_src_dir="hypre"
HYPRE_SOURCE_DIR="$pkg_sources_dir/$pkg_src_dir"
pkg_bld_dir="$OUT_DIR/hypre"
HYPRE_DIR="$pkg_bld_dir"


function hypre_clone()
{
   local pkg="hypre"
   pkg_repo_list=("git@github.com:LLNL/hypre.git"
                  "https://github.com/LLNL/hypre.git")
   cd "$pkg_sources_dir" || return 1
   [[ -d "$pkg_src_dir" ]] && return 0
   for pkg_repo in "${pkg_repo_list[@]}"; do
      echo "Cloning $pkg from $pkg_repo ..."
      git clone "$pkg_repo" "$pkg_src_dir" && return 0
   done
   echo "Could not successfully clone $pkg. Stop."
   exit 1
}


function hypre_build()
{
   local pkg="hypre"
   if [[ ! -d "$pkg_bld_dir" ]]; then
      cd "$OUT_DIR" && git clone "$HYPRE_SOURCE_DIR" || {
         echo "Cloning $HYPRE_SOURCE_DIR to OUT_DIR failed. Stop."
         exit 1
      }
   elif [[ -e "${pkg_bld_dir}_build_successful" ]]; then
      echo "Using successfully built $pkg from OUT_DIR."
      return 0
   fi
   echo "Building $pkg, sending output to ${pkg_bld_dir}_build.log ..." && {
      cd "$pkg_bld_dir/src" && \
      if [[ -e config/Makefile.config ]]; then
         make distclean
      fi && \
      ./configure \
         CC="$mpi_cc" \
         CXX="$mpi_cxx" \
         CFLAGS="$CFLAGS" \
         CXXFLAGS="$CFLAGS" \
         $HYPRE_EXTRA_CONFIG \
         --disable-fortran \
         --without-fei && \
      make -j $num_proc_build
   } &> "${pkg_bld_dir}_build.log" || {
      echo " ... building $pkg FAILED, see log for details."
      exit 1
   }
   echo "Build succesful."
   : > "${pkg_bld_dir}_build_successful"
}


hypre_clone && hypre_build
