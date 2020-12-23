#Using Block Jacobi Smoother
 ./go.sh -c lassen -m xlc -r tests/mfem_amgx/ex1p_amgx_solver.sh num_proc_build=36  >> amgx_solver.txt &&
 ./go.sh -c lassen -m xlc -r tests/mfem_amgx/ex1p_amgx_precon.sh num_proc_build=36  >> amgx_precon.txt &&
 ./go.sh -c lassen -m xlc -r tests/mfem_amgx/ex1p_amgx_agg_precon.sh num_proc_build=36  >> amgx_agg_precon.txt &&
 ./go.sh -c lassen -m xlc -r tests/mfem_amgx/ex1p_amgx_agg_solver.sh num_proc_build=36  >> amgx_agg_solver.txt
