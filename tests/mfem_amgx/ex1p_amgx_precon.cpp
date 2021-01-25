//                       MFEM Example 1 - Parallel Version
//
// Compile with: make ex1p
//
// Description:  This example code demonstrates the use of MFEM to define a
//               simple finite element discretization of the Laplace problem
//               -Delta u = 1 with homogeneous Dirichlet boundary conditions.
//               Specifically, we discretize using a FE space of the specified
//               order, or if order < 1 using an isoparametric/isogeometric
//               space (i.e. quadratic for quadratic curvilinear mesh, NURBS for
//               NURBS mesh, etc.)

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

Mesh *make_mesh(int myid, int num_procs, int dim, int level,
                int &par_ref_levels, Array<int> &nxyz);

int main(int argc, char *argv[])
{
   // 1. Initialize MPI.
   int num_procs, myid;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);

   // 2. Parse command-line options.
   int dim = 3;
   int level = 0;
   int order = 1;
   const char *device_config = "cuda";

   OptionsParser args(argc, argv);
   args.AddOption(&dim, "-dim", "--mesh-dimension",
                  "Solve 2D or 3D problem.");
   args.AddOption(&level, "-l", "--refinement-level",
                  "Set the problem size: 2^level mesh elements per processor.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
   args.Parse();
   if (!args.Good())
   {
      if (myid == 0)
      {
         args.PrintUsage(cout);
      }
      MPI_Finalize();
      return 1;
   }
   if (myid == 0)
   {
      args.PrintOptions(cout);
   }

   // 3. Enable hardware devices such as GPUs, and programming models such as
   //    CUDA, OCCA, RAJA and OpenMP based on command line options.
   Device device(device_config);
   if (myid == 0) { device.Print(); }

   // 4. Read the (serial) mesh from the given mesh file on all processors.  We
   //    can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
   //    and volume meshes with the same code.
   int par_ref_levels;
   Array<int> nxyz;
   Mesh *mesh = make_mesh(myid, num_procs, dim, level, par_ref_levels, nxyz);

   // 5. Refine the serial mesh on all processors to increase the resolution. In
   //    this example we do 'ref_levels' of uniform refinement. We choose
   //    'ref_levels' to be the largest number that gives a final mesh with no
   //    more than 10,000 elements.
   {
     // Disabled
   }

   // 6. Define a parallel mesh by a partitioning of the serial mesh. Refine
   //    this mesh further in parallel to increase the resolution. Once the
   //    parallel mesh is defined, the serial mesh can be deleted.
   int *partitioning = nxyz.Size() ? mesh->CartesianPartitioning(nxyz) : NULL;
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh, partitioning);
   delete [] partitioning;
   delete mesh;
   for (int l = 0; l < par_ref_levels; l++)
   {
      pmesh->UniformRefinement();
   }
   // pmesh->PrintInfo();
   long global_ne = pmesh->ReduceInt(pmesh->GetNE());
   if (myid == 0)
   {
      cout << "Total number of elements: " << global_ne << endl;
   }

   // 7. Define a parallel finite element space on the parallel mesh. Here we
   //    use continuous Lagrange finite elements of the specified order. If
   //    order < 1, we instead use an isoparametric/isogeometric space.
   FiniteElementCollection *fec;
   MFEM_VERIFY(order > 0, "invalid 'order': " << order);
   fec = new H1_FECollection(order, dim);
   ParFiniteElementSpace *fespace = new ParFiniteElementSpace(pmesh, fec);
   HYPRE_Int size = fespace->GlobalTrueVSize();
   if (myid == 0)
   {
      cout << "Number of finite element unknowns: " << size << endl;
   }

   // 8. Determine the list of true (i.e. parallel conforming) essential
   //    boundary dofs. In this example, the boundary conditions are defined
   //    by marking all the boundary attributes from the mesh as essential
   //    (Dirichlet) and converting them to a list of true dofs.
   Array<int> ess_tdof_list;
   if (pmesh->bdr_attributes.Size())
   {
      Array<int> ess_bdr(pmesh->bdr_attributes.Max());
      ess_bdr = 1;
      fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }

   // 9. Set up the parallel linear form b(.) which corresponds to the
   //    right-hand side of the FEM linear system, which in this case is
   //    (1,phi_i) where phi_i are the basis functions in fespace.
   ParLinearForm *b = new ParLinearForm(fespace);
   ConstantCoefficient one(1.0);
   b->AddDomainIntegrator(new DomainLFIntegrator(one));
   b->Assemble();

   // 10. Define the solution vector x as a parallel finite element grid function
   //     corresponding to fespace. Initialize x with initial guess of zero,
   //     which satisfies the boundary conditions.
   ParGridFunction x(fespace);
   x = 0.0;

   // 11. Set up the parallel bilinear form a(.,.) on the finite element space
   //     corresponding to the Laplacian operator -Delta, by adding the Diffusion
   //     domain integrator.
   ParBilinearForm *a = new ParBilinearForm(fespace);
   a->AddDomainIntegrator(new DiffusionIntegrator(one));

   // 12. Assemble the parallel bilinear form and the corresponding linear
   //     system, applying any necessary transformations such as: parallel
   //     assembly, eliminating boundary conditions, applying conforming
   //     constraints for non-conforming AMR, static condensation, etc.
   a->Assemble();

   OperatorPtr A;
   Vector B, X;
   a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);

   // 13. Solve the linear system A X = B.
   //     * With full assembly, use the preconditioner from AmgX.

   bool amgx_verbose = true;

   string amgx_config = "{\n"
                    " \"config_version\": 2, \n"
                    " \"solver\": { \n"
                    "   \"solver\": \"AMG\", \n"
                    "     \"smoother\": { \n"
                    "     \"scope\": \"jacobi\", \n"
                    "     \"solver\": \"BLOCK_JACOBI\", \n"
                    "     \"relaxation_factor\": 0.7 \n"
                    "       }, \n"
                    "   \"presweeps\": 1, \n"
                    "   \"interpolator\": \"D2\", \n"
                    "   \"max_row_sum\" : 0.9, \n"
                    "   \"strength_threshold\" : 0.25, \n"
                    "   \"postsweeps\": 1, \n"
                    "   \"max_iters\": 1, \n"
                    "   \"cycle\": \"V\"";
   if (amgx_verbose)
   {
     amgx_config = amgx_config + ",\n"
                   "   \"obtain_timings\": 1 \n";
   }
   else
   {
     amgx_config = amgx_config + "\n";
   }
   amgx_config = amgx_config + " }\n" + "}\n";

   AmgXSolver *prec = new AmgXSolver;
   prec->SetConvergenceCheck(false);
   prec->ReadParameters(amgx_config, AmgXSolver::INTERNAL);
   prec->InitExclusiveGPU(MPI_COMM_WORLD);

   CGSolver cg(MPI_COMM_WORLD);
   const int max_cg_iter = 200;
   const int cg_print_level = 3;
   cg.SetRelTol(1e-12);
   cg.SetMaxIter(max_cg_iter);
   cg.SetPrintLevel(cg_print_level);
   if (prec) { cg.SetPreconditioner(*prec); }
   cg.SetOperator(*A);

   // Warm-up CG solve (in case of JIT to avoid timing it)
   {
      Vector Xtmp(X);
      cg.SetMaxIter(2);
      cg.SetPrintLevel(-1);
      cg.Mult(B, Xtmp);
      cg.SetMaxIter(max_cg_iter);
      cg.SetPrintLevel(cg_print_level);
   }

   // Sync all ranks
   MPI_Barrier(pmesh->GetComm());
   
   // Start CG timing.
   tic_toc.Clear();
   
   // Start & Stop CG timing.
   tic_toc.Start();
   cg.Mult(B, X);
   tic_toc.Stop();
   double rt_min, rt_max, my_rt;
   my_rt = tic_toc.RealTime();
   MPI_Reduce(&my_rt, &rt_min, 1, MPI_DOUBLE, MPI_MIN, 0, pmesh->GetComm());
   MPI_Reduce(&my_rt, &rt_max, 1, MPI_DOUBLE, MPI_MAX, 0, pmesh->GetComm());
   
   // Print timing results.
   if (myid == 0)
   {
      int cg_iter = cg.GetNumIterations();
      // Note: In the pcg algorithm, the number of operator Mult() calls is
      //       N_iter and the number of preconditioner Mult() calls is N_iter+1.
      cout << '\n'
           << "Total iterations "<< cg_iter;
      cout << '\n'
           << "Total CG time:    " << rt_max << " (" << rt_min << ") sec."
           << endl;
      cout << "Time per CG step: "
           << rt_max / cg_iter << " ("
           << rt_min / cg_iter << ") sec." << endl;
      cout << "\n\"DOFs/sec\" in CG: "
           << 1e-6*size*cg_iter/rt_max << " ("
           << 1e-6*size*cg_iter/rt_min << ") million.\n"
           << endl;
   }
   delete prec;

   // 14. Recover the parallel grid function corresponding to X. This is the
   //     local finite element solution on each processor.
   a->RecoverFEMSolution(X, *b, x);

   // 15. Save the refined mesh and the solution in parallel. This output can
   //     be viewed later using GLVis: "glvis -np <np> -m mesh -g sol".
   {
      //disabled
   }

   // 16. Send the solution by socket to a GLVis server.
   //if (visualization)
   {
     //disabled
   }

   // 17. Free the used memory.
   delete a;
   delete b;
   delete fespace;
   if (order > 0) { delete fec; }
   delete pmesh;

   MPI_Finalize();

   return 0;
}


Mesh *make_mesh(int myid, int num_procs, int dim, int level,
                int &par_ref_levels, Array<int> &nxyz)
{
   int log_p = (int)floor(log((double)num_procs)/log(2.0) + 0.5);
   MFEM_VERIFY((1 << log_p) == num_procs,
               "number of processor is not a power of 2: " << num_procs);
   MFEM_VERIFY(dim == 3, "dim = " << dim << " is NOT implemented!");

   // Determine processor decomposition.
   int s[3];
   s[0] = log_p/3 + (log_p%3 > 0 ? 1 : 0);
   s[1] = log_p/3 + (log_p%3 > 1 ? 1 : 0);
   s[2] = log_p/3;
   nxyz.SetSize(dim);
   nxyz[0] = 1 << s[0];
   nxyz[1] = 1 << s[1];
   nxyz[2] = 1 << s[2];

   // Determine mesh size.
   int ser_level = level%3;
   par_ref_levels = level/3;
   int log_n = log_p + ser_level;
   int t[3];
   t[0] = log_n/3 + (log_n%3 > 0 ? 1 : 0);
   t[1] = log_n/3 + (log_n%3 > 1 ? 1 : 0);
   t[2] = log_n/3;

   // Create the Mesh.
   const bool gen_edges = true;
   const bool sfc_ordering = true;
   Mesh *mesh = new Mesh(1 << t[0], 1 << t[1], 1 << t[2],
                         Element::HEXAHEDRON, gen_edges,
                         1.0, 1.0, 1.0, sfc_ordering);
   if (myid == 0)
   {
      cout << "Processor partitioning: ";
      nxyz.Print(cout, dim);

      // Mesh dimensions AFTER parallel refinement:
      cout << "Mesh dimensions: "
           << (1 << (t[0]+par_ref_levels)) << ' '
           << (1 << (t[1]+par_ref_levels)) << ' '
           << (1 << (t[2]+par_ref_levels)) << endl;
   }

   return mesh;
}
