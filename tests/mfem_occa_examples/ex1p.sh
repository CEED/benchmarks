# This file is part of CEED. For more details, see exascaleproject.org.


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
