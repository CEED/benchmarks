# This file is part of CEED. For more details, see exascaleproject.org.

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
   pkg_repo_list=("https://${llnluser}@lc.llnl.gov/bitbucket/scm/~fisher47/acrotensor.git")
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
   if [[ ! -d "$pkg_bld_dir" ]]; then
      cd "$OUT_DIR" && git clone "$ACROTENSOR_SOURCE_DIR" || {
         echo "Cloning $ACROTENSOR_SOURCE_DIR to OUT_DIR failed. Stop."
         return 1
      }
   elif [[ -e "${pkg_bld_dir}_build_successful" ]]; then
      echo "Using successfully built $pkg from OUT_DIR."
      return 0
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
   acrotensor_clone && acrotensor_build
}
