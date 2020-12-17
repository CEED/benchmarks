#Using Block JACOBI SMOOTHER
 ./go.sh -c lassen -m xlc -r tests/mfem_amgx/ex1p_amgx_precon.sh num_proc_build=36  >> amgx_precon &&
 ./go.sh -c lassen -m xlc -r tests/mfem_amgx/ex1p_amgx_solver.sh num_proc_build=36  >> amgx_solver &&
 ./go.sh -c lassen -m xlc -r tests/mfem_amgx/ex1p_amgx_agg_precon.sh num_proc_build=36  >> amgx_agg_precon &&
 ./go.sh -c lassen -m xlc -r tests/mfem_amgx/ex1p_amgx_agg_solver.sh num_proc_build=36  >> amgx_agg_solver &&

#Using L1 JACOBI SMOOTHER
 ./go.sh -c lassen -m xlc -r tests/mfem_amgx/ex1p_amgx_precon_l1_jacobi.sh num_proc_build=36  >> amgx_precon_l1_jacobi &&
 ./go.sh -c lassen -m xlc -r tests/mfem_amgx/ex1p_amgx_solver_l1_jacobi.sh num_proc_build=36  >> amgx_solver_l1_jacobi &&
 ./go.sh -c lassen -m xlc -r tests/mfem_amgx/ex1p_amgx_agg_precon_l1_jacobi.sh num_proc_build=36  >> amgx_agg_precon_l1_jacobi &&
 ./go.sh -c lassen -m xlc -r tests/mfem_amgx/ex1p_amgx_agg_solver_l1_jacobi.sh num_proc_build=36  >> amgx_agg_solver_l1_jacobi
