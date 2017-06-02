# This file is part of CEED. For more details, see exascaleproject.org.

# Download and build SUNDIALS.

if [[ -z "$pkg_sources_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   return 1
fi
if [[ -z "$OUT_DIR" ]]; then
   echo "The variable 'OUT_DIR' is not set. Stop."
   return 1
fi
SUNDIALS_VERSION="${SUNDIALS_VERSION:-2.7.0}"
pkg_src_url_base="https://computation.llnl.gov/projects/sundials/download"
pkg_src_dir="sundials-$SUNDIALS_VERSION"
pkg_src_file="${pkg_src_dir}.tar.gz"
pkg_src_url="$pkg_src_url_base/$pkg_src_file"
SUNDIALS_SOURCE_DIR="$pkg_sources_dir/$pkg_src_dir"
SUNDIALS_SOURCE_FILE="$pkg_sources_dir/$pkg_src_file"
pkg_bld_dir="$OUT_DIR/$pkg_src_dir"
SUNDIALS_DIR="$pkg_bld_dir"
pkg="SUNDIALS v$SUNDIALS_VERSION"


function sundials_download()
{
   # Download and extract SUNDIALS in $pkg_sources_dir for out-of-source build.
   cd "$pkg_sources_dir" || return 1
   [[ -d "$pkg_src_dir" ]] && return 0
   if [[ ! -e "$pkg_src_file" ]]; then
      echo "Downloading $pkg from $pkg_src_url ..."
      {
         command -v wget &> /dev/null && wget "$pkg_src_url"
      } || {
         command -v curl &> /dev/null && curl -O "$pkg_src_url"
      } || {
         echo "Could not download $pkg. Stop."
         return 1
      }
   fi
   echo "Extracting $pkg_src_file ..."
   tar zxf "$pkg_src_file" || {
      echo "Error extracting \"$pkg_src_file\". Stop."
      return 1
   }
}


function sundials_build()
{
   # Use an out-of-source build.
   if package_build_is_good; then
      echo "Using successfully built $pkg from OUT_DIR."
      return 0
   elif [[ ! -d "$pkg_bld_dir" ]]; then
      mkdir -p "$pkg_bld_dir" || return 1
   fi
   local HYPRE_CMAKE_OPTS=("-DHYPRE_ENABLE=OFF")
   if [[ -n "$HYPRE_DIR" ]]; then
      HYPRE_CMAKE_OPTS=(
         "-DHYPRE_ENABLE=ON"
         "-DHYPRE_INCLUDE_DIR=$HYPRE_DIR/src/hypre/include"
         "-DHYPRE_LIBRARY_DIR=$HYPRE_DIR/src/hypre/lib")
   else
      echo "${magenta}Warning: Building $pkg without HYPRE ...${none}"
   fi
   local LAPACK_CMAKE_OPTS=("-DLAPACK_ENABLE=OFF")
   if [[ -n "$LAPACK_LIB" ]]; then
      LAPACK_CMAKE_OPTS=(
         "-DLAPACK_ENABLE=ON"
         "-DLAPACK_LIBRARIES=\"$LAPACK_LIB\"")
      array_union LDFLAGS "$LAPACK_LIB"
   else
      echo "${magenta}Warning: Building $pkg without LAPACK ...${none}"
   fi
   local SUITESPARSE_CMAKE_OPTS=("-DKLU_ENABLE=OFF")
   if [[ -n "$SUITESPARSE_DIR" ]]; then
      SUITESPARSE_CMAKE_OPTS=(
         "-DKLU_ENABLE=ON"
         "-DKLU_INCLUDE_DIR=$SUITESPARSE_DIR/include"
         "-DKLU_LIBRARY_DIR=$SUITESPARSE_DIR/lib")
   else
      echo "${magenta}Warning: Building $pkg without SuiteSparse/KLU ...${none}"
   fi
   echo "Building $pkg, sending output to ${pkg_bld_dir}_build.log ..." && \
   {
      cd "$pkg_bld_dir" && \
      rm -rf ./* && \
      cmake "$SUNDIALS_SOURCE_DIR" \
         -DCMAKE_C_COMPILER="$MPICC" \
         -DCMAKE_C_FLAGS_RELEASE="$CFLAGS" \
         -DCMAKE_BUILD_TYPE="Release" \
         -DBUILD_SHARED_LIBS="ON" \
         -DCMAKE_INSTALL_NAME_DIR="$pkg_bld_dir/lib" \
         -DCMAKE_INSTALL_RPATH="$pkg_bld_dir/lib" \
         -DCMAKE_INSTALL_RPATH_USE_LINK_PATH="ON" \
         -DCMAKE_INSTALL_PREFIX="$pkg_bld_dir" \
         -DEXAMPLES_ENABLE="OFF" \
         -DMPI_ENABLE="ON" \
         "${HYPRE_CMAKE_OPTS[@]}" \
         "${LAPACK_CMAKE_OPTS[@]}" \
         "${SUITESPARSE_CMAKE_OPTS[@]}" \
         -DCMAKE_VERBOSE_MAKEFILE="OFF" && \
      make -j $num_proc_build && \
      make install
   } &> "${pkg_bld_dir}_build.log" || {
      echo " ... building $pkg FAILED, see log for details."
      return 1
   }
   echo "Build succesful."
   : > "${pkg_bld_dir}_build_successful"
}


function build_package()
{
   sundials_download && sundials_build
}
