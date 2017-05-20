# This file is part of CEED. For more details, see exascaleproject.org.


function build_and_run_tests()
{
   $dry_run cd "$MFEM_DIR/examples/occa" && \
   $dry_run make ex1 && {

      $dry_run occa info

      $dry_run ./ex1 \
         --mesh ../../data/fichera.mesh \
         --order 3 \
         --preconditioner none \
         --device-info "mode: 'Serial'" \
         --no-occa-verbose \
         --no-visualization

   } || return 1
}


test_required_packages="metis hypre occa mfem-occa"
