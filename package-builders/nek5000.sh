# This file is part of CEED. For more details, see exascaleproject.org.

# Clone Nek5000.

if [[ -z "$pkg_sources_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   return 1
fi
if [[ -z "$OUT_DIR" ]]; then
   echo "The variable 'OUT_DIR' is not set. Stop."
   return 1
fi

pkg_src_dir="Nek5000"
NEK5K_SOURCE_DIR="$pkg_sources_dir/$pkg_src_dir"
pkg_bld_dir="$OUT_DIR/$pkg_src_dir"
NEK5K_DIR="$pkg_bld_dir"
pkg="Nek5000"

function nek5k_clone()
{
   pkg_repo_list=("git@github.com:Nek5000/Nek5000.git"
                  "https://github.com/Nek5000/Nek5000.git")
   pkg_git_branch="master"
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

function nek5k_build()
{
   if [[ ! -d "$pkg_bld_dir" ]]; then
      cd "$OUT_DIR" && git clone "$NEK5K_SOURCE_DIR" || {
         echo "Cloning $NEK5K_SOURCE_DIR to OUT_DIR failed. Stop."
         return 1
      }
   elif [[ -e "${pkg_bld_dir}_build_successful" ]]; then
      echo "Using successfully built $pkg from OUT_DIR."
      return 0
   fi
   [[ -n "$FC" && -n "$CC" ]] || {
      echo "$pkg requires both Fortran and C compilers."
      echo "Please update the used config/compiler settings. Stop."
      return 1
   }
   echo "Building $pkg, sending output to ${pkg_bld_dir}_build.log ..." && {
      ## Just build the requited tools: genbox and genmap
      cd "$pkg_bld_dir/tools" && {
         # replace set F77 and CC in 'maketools'
         [[ -e "maketools.orig" ]] || cp -p maketools maketools.orig
         sed -e "s/^F77=.*$/F77=\"$FC $FFLAGS\"/" \
             -e "s/^CC=.*$/CC=\"$CC $CFLAGS\"/" \
             maketools.orig > maketools
      } && \
      ./maketools genmap && \
      ./maketools genbox
   } &> "${pkg_bld_dir}_build.log" || {
      echo " ... building $pkg FAILED, see log for details."
      return 1
   }
   echo "Build succesful."
   : > "${pkg_bld_dir}_build_successful"
}

function build_package()
{
   nek5k_clone && nek5k_build
}
