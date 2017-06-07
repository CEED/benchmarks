function configure_tests()
{
    case_orders=( 1 2 3 4 5 6 7 8 )
    case_refs=( 6 5 4 4 4 3 3 3 )
}

function setup_occa_test()
{
    # This test name
    local test_name=dtp_occa

    # Export OCCA_CACHE_DIR
    export OCCA_CACHE_DIR=$test_exe_dir/.occa

    # configuration from dtp.sh.inl
    configure_tests

    # Build the executable
    $dry_run cd "$test_dir"

    $dry_run make "$test_exe_dir/$test_name" "MFEM_DIR=$MFEM_DIR" "BLD=$test_exe_dir/"
    [ -z $dry_run ] && {
        occa info && occa clear -ay && make cache-kernels "MFEM_DIR=$MFEM_DIR" "MFEM_SOURCE_DIR=$MFEM_SOURCE_DIR"
    }
    $dry_run cp $MFEM_DIR/data/fichera.mesh $test_exe_dir/
}
