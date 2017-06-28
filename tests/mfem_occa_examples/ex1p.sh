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
# testbed platforms, in support of the nation's exascale computing imperative.


function build_and_run_tests()
{
   local dev_info= occa_verbose_opt=
   $dry_run cd "$MFEM_DIR/examples/occa" && \
   $dry_run make ex1p && {

      local occa_verbose="0"

      dev_info="mode: 'Serial'"

      # dev_info="mode: 'CUDA', deviceID: 0"

      # dev_info="mode: 'OpenMP'"
      # $dry_run export OMP_NUM_THREADS="1"

      if [[ "$occa_verbose" = "1" ]]; then
         export OCCA_VERBOSE=1
         occa_verbose_opt="--occa-verbose"
      else
         occa_verbose_opt="--no-occa-verbose"
      fi

      set_mpi_options

      # $dry_run $mpi_run occa info

      $dry_run $mpi_run ./ex1p \
         --mesh ../../data/fichera.mesh \
         --order 3 \
         --preconditioner none \
         --device-info "$dev_info" \
         "$occa_verbose_opt" \
         --no-visualization

   }
}


test_required_packages="metis hypre occa mfem-occa"
