#Using L1 Jacobi Smoother
./go.sh -c lassen -m xlc -r tests/mfem_amgx/ex1p_amgx_solver_l1_jacobi.sh num_proc_build=36  >> amgx_solver_l1_jacobi.txt &&
./go.sh -c lassen -m xlc -r tests/mfem_amgx/ex1p_amgx_precon_l1_jacobi.sh num_proc_build=36  >> amgx_precon_l1_jacobi.txt
#./go.sh -c lassen -m xlc -r tests/mfem_amgx/ex1p_amgx_agg_precon_l1_jacobi.sh num_proc_build=36  >> amgx_agg_precon_l1_jacobi.txt &
#./go.sh -c lassen -m xlc -r tests/mfem_amgx/ex1p_amgx_agg_solver_l1_jacobi.sh num_proc_build=36  >> amgx_agg_solver_l1_jacobi.txt &
