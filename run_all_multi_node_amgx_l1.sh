#Using L1 Jacobi Smoother

# 4 GPUS   
./go.sh -c lassen -m xlc -r tests/mfem_amgx/ex1p_amgx_solver_l1_jacobi.sh -n 4 >> amgx_solver_l1_jacobi_4.txt &&
./go.sh -c lassen -m xlc -r tests/mfem_amgx/ex1p_amgx_precon_l1_jacobi.sh -n 4  >> amgx_precon_l1_jacobi_4.txt &&

#16 GPUS
./go.sh -c lassen -m xlc -r tests/mfem_amgx/ex1p_amgx_solver_l1_jacobi.sh -n 16 >> amgx_solver_l1_jacobi_16.txt &&
./go.sh -c lassen -m xlc -r tests/mfem_amgx/ex1p_amgx_precon_l1_jacobi.sh -n 16 >> amgx_precon_l1_jacobi_16.txt

#32 GPUS
#./go.sh -c lassen -m xlc -r tests/mfem_amgx/ex1p_amgx_solver_l1_jacobi.sh -n 32 >> amgx_solver_l1_jacobi_32.txt &&
#./go.sh -c lassen -m xlc -r tests/mfem_amgx/ex1p_amgx_precon_l1_jacobi.sh -n 32 >> amgx_precon_l1_jacobi_32.txt &&

#64 GPUS
#./go.sh -c lassen -m xlc -r tests/mfem_amgx/ex1p_amgx_solver_l1_jacobi.sh -n 64 >> amgx_solver_l1_jacobi_64.txt &&
#./go.sh -c lassen -m xlc -r tests/mfem_amgx/ex1p_amgx_precon_l1_jacobi.sh -n 64 >> amgx_precon_l1_jacobi_64.txt 
