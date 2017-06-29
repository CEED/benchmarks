#include "mfem-performance.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

#ifndef PROBLEM
#define PROBLEM 1
#endif
#ifndef SOL_P
#define SOL_P 3
#endif
#ifndef IR_ORDER
#define IR_ORDER 2*(SOL_P+3)-1
#endif

// Define template parameters for optimized build.
const int            dim      = 3;
const Geometry::Type geom     = Geometry::CUBE; // mesh elements  (default: hex)
const int            mesh_p   = 1;              // mesh curvature (default: 1)
const int            sol_p    = SOL_P;      // solution order (default: 3)
const int            rdim     = Geometry::Constants<geom>::Dimension;
const int            ir_order = IR_ORDER;

// Static mesh type
typedef H1_FiniteElement<geom,mesh_p>         mesh_fe_t;
typedef H1_FiniteElementSpace<mesh_fe_t>      mesh_fes_t;
typedef VectorLayout<Ordering::byNODES,dim>   mesh_layout_t;
typedef VectorLayout<Ordering::byNODES,1>     scal_layout_t;
typedef VectorLayout<Ordering::byVDIM,dim>    vec_layout_t;
typedef TMesh<mesh_fes_t,mesh_layout_t>       mesh_t;

// Static solution finite element space type
typedef H1_FiniteElement<geom,sol_p>          sol_fe_t;
typedef H1_FiniteElementSpace<sol_fe_t>       sol_fes_t;

// Static quadrature, coefficient and integrator types
typedef TIntegrationRule<geom,ir_order>       int_rule_t;
typedef TConstantCoefficient<>                coeff_t;
typedef TIntegrator<coeff_t,TMassKernel>      mass_integ_t;
typedef TIntegrator<coeff_t,TDiffusionKernel> diffusion_integ_t;

// Static bilinear form type, combining the above types
#if PROBLEM == 1
typedef TBilinearForm<mesh_t,sol_fes_t,int_rule_t,mass_integ_t,scal_layout_t> HPCBilinearForm;
#elif PROBLEM == 2
typedef TBilinearForm<mesh_t,sol_fes_t,int_rule_t,diffusion_integ_t,scal_layout_t> HPCBilinearForm;
#elif PROBLEM == 3
typedef TBilinearForm<mesh_t,sol_fes_t,int_rule_t,mass_integ_t,vec_layout_t> HPCBilinearForm;
#elif PROBLEM == 4
typedef TBilinearForm<mesh_t,sol_fes_t,int_rule_t,diffusion_integ_t,vec_layout_t> HPCBilinearForm;
#else
#error "Invalid bake-off problem."
#endif

int main(int argc, char *argv[])
{
   // Initialize MPI.
   int num_procs, myid;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);

   // Parse command-line options.
   int order = sol_p;
   const char *pc = "none";
   bool visualization = 1;
   int num_procs_x = num_procs;
   int num_procs_y = 1;
   int num_procs_z = 1;
   int el_per_proc_x = 4;
   int el_per_proc_y = 4;
   int el_per_proc_z = 4;
   const bool vec = PROBLEM == 3 || PROBLEM == 4;

   OptionsParser args(argc, argv);
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree) or -1 for"
                  " isoparametric space.");
   args.AddOption(&pc, "-pc", "--preconditioner",
                  "Preconditioner: "
                  "lor - low-order-refined (matrix-free) AMG, "
                  "ho - high-order (assembled) AMG, "
                  "jacobi, "
                  "lumpedmass, "
                  "none.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&num_procs_x, "-nx", "--num-procs-x",
                  "Number of MPI ranks in x-dimension.");
   args.AddOption(&num_procs_y, "-ny", "--num-procs-y",
                  "Number of MPI ranks in y-dimension.");
   args.AddOption(&num_procs_z, "-nz", "--num-procs-z",
                  "Number of MPI ranks in z-dimension.");
   args.AddOption(&el_per_proc_x, "-ex", "--num-el-per-proc-x",
                  "Number of elements per MPI rank in x-dimension.");
   args.AddOption(&el_per_proc_y, "-ey", "--num-el-per-proc-y",
                  "Number of elements per MPI rank in y-dimension.");
   args.AddOption(&el_per_proc_z, "-ez", "--num-el-per-proc-z",
                  "Number of elements per MPI rank in z-dimension.");
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

   if (num_procs_x * num_procs_y * num_procs_z != num_procs)
   {
      mfem_error("Invalid dimensions for MPI ranks");
      return -1;
   }

   enum PCType { NONE, LOR, HO, JACOBI, LUMPEDMASS };
   PCType pc_choice;
   if (!strcmp(pc, "ho"))          { pc_choice = HO; }
   else if (!strcmp(pc, "lor"))    { pc_choice = LOR; }
   else if (!strcmp(pc, "jacobi")) { pc_choice = JACOBI; }
   else if (!strcmp(pc, "lumpedmass")) { pc_choice = LUMPEDMASS; }
   else if (!strcmp(pc, "none"))   { pc_choice = NONE; }
   else
   {
      mfem_error("Invalid preconditioner specified");
      return 3;
   }


   // Generate serial mesh
   Mesh *mesh;
   const int nx = num_procs_x * el_per_proc_x;
   const int ny = num_procs_y * el_per_proc_y;
   const int nz = num_procs_z * el_per_proc_z;
   mesh = new Mesh(nx, ny, nz, Element::HEXAHEDRON, 1,
                   1.0, 1.0, 1.0);

   // Check if the generated mesh matches the optimized version
   if (myid == 0)
   {
       cout << "High-performance version using integration rule with "
            << int_rule_t::qpts << " points ..." << endl;
   }
   if (!mesh_t::MatchesGeometry(*mesh))
   {
       if (myid == 0)
       {
           cout << "The given mesh does not match the optimized 'geom' parameter.\n"
                << "Recompile with suitable 'geom' value." << endl;
       }
       delete mesh;
       MPI_Finalize();
       return 4;
   }
   else if (!mesh_t::MatchesNodes(*mesh))
   {
       if (myid == 0)
       {
           cout << "Switching the mesh curvature to match the "
                << "optimized value (order " << mesh_p << ") ..." << endl;
       }
       mesh->SetCurvature(mesh_p, false, dim, Ordering::byNODES);
   }

   // Define a parallel mesh by a partitioning of the serial
   // mesh. Refine this mesh further in parallel to increase the
   // resolution. Once the parallel mesh is defined, the serial mesh
   // can be deleted.
   if (myid == 0)
   {
      cout << "Initializing parallel mesh ..." << endl;
   }
   int *partitioning = new int[nx*ny*nz];
   for (int k=0; k<nz; ++k)
   {
      const int pz = k / el_per_proc_z;
      for (int j=0; j<ny; ++j)
      {
         const int py = j / el_per_proc_y;
         for (int i=0; i<nx; ++i)
         {
            const int px = i / el_per_proc_x;
            partitioning[i+j*nx+k*nx*ny] = px + py * num_procs_x + pz * num_procs_x * num_procs_y;
         }
      }
   }
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh, partitioning);
   delete mesh;
   if (pmesh->MeshGenerator() & 1) // simplex mesh
   {
      MFEM_VERIFY(pc_choice != LOR,
                  "triangle and tet meshes do not support the LOR "
                  "preconditioner yet");
   }

   // Define a parallel finite element space on the parallel mesh
   const int basis = BasisType::GaussLobatto;
   FiniteElementCollection *fec = new H1_FECollection(order, dim, basis);
   ParFiniteElementSpace *fespace
     = new ParFiniteElementSpace(pmesh, fec,
                                 vec ? dim : 1,
                                 vec ? Ordering::byVDIM : Ordering::byNODES);
   HYPRE_Int size = fespace->GlobalTrueVSize();
   if (myid == 0)
   {
      cout << "Number of finite element unknowns: " << size << endl;
   }
   ParMesh *pmesh_lor = NULL;
   FiniteElementCollection *fec_lor = NULL;
   ParFiniteElementSpace *fespace_lor = NULL;
   if (pc_choice == LOR)
   {
      pmesh_lor = new ParMesh(pmesh, order, basis);
      fec_lor = new H1_FECollection(1, dim);
      fespace_lor
        = new ParFiniteElementSpace(pmesh_lor,
                                    fec_lor,
                                    vec ? dim : 1,
                                    vec ? Ordering::byVDIM : Ordering::byNODES);
   }

   // Check if the optimized version matches the given space
   if (!sol_fes_t::Matches(*fespace))
   {
      if (myid == 0)
      {
         cout << "The given order does not match the optimized parameter.\n"
              << "Recompile with suitable 'sol_p' value." << endl;
      }
      delete fespace;
      delete fec;
      delete mesh;
      MPI_Finalize();
      return 5;
   }

   // Determine the list of true (i.e. parallel conforming) essential
   // boundary dofs
   Array<int> ess_tdof_list;
   if (pmesh->bdr_attributes.Size())
   {
      Array<int> ess_bdr(pmesh->bdr_attributes.Max());
      ess_bdr = 0;
      fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }

   // Define the solution vector x and RHS vector b
   // Note: subtract mean value if solving stiffness matrix problem
   ParGridFunction x0(fespace), x(fespace), b(fespace), ones(fespace);
   ones = 1.0;
   x0.Randomize();
#if PROBLEM == 2
   x0 -= (x0 * ones) / x0.Size();
#elif PROBLEM == 4
   for (int d=0; d<dim; ++d)
   {
      const int ndofs = fespace->GetNDofs();
      double mean = 0;
      for (int i=d*ndofs; i<(d+1)*ndofs; ++i)
      {
         mean += x0(i);
      }
      MPI_Allreduce(MPI_IN_PLACE, &mean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      mean /= (x0.Size() / dim);
      for (int i=d*ndofs; i<(d+1)*ndofs; ++i)
      {
         x0(i) -= mean;
      }
   }
#endif
   x = x0;
   b = -1.0;

   // Set up bilinear form for preconditioner
   ParBilinearForm *a_pc = NULL;
   if (pc_choice == LOR)
   {
      a_pc = new ParBilinearForm(fespace_lor);
   }
   if (pc_choice == HO || pc_choice == JACOBI)
   {
      a_pc = new ParBilinearForm(fespace);
   }

   // 13. Assemble the parallel bilinear form and the corresponding linear
   //     system, applying any necessary transformations such as: parallel
   //     assembly, eliminating boundary conditions, applying conforming
   //     constraints for non-conforming AMR, static condensation, etc.
   if (myid == 0)
   {
      cout << "Assembling the matrix ..." << flush;
   }
   tic_toc.Clear();
   tic_toc.Start();

   // High-performance assembly/evaluation using the templated operator type
   Operator *a = NULL;
#if PROBLEM == 1 || PROBLEM == 3
   a = new HPCBilinearForm(mass_integ_t(coeff_t(1.0)), *fespace);
#else
   a = new HPCBilinearForm(diffusion_integ_t(coeff_t(1.0)), *fespace);
#endif
   ((HPCBilinearForm*) a)->Assemble();

   tic_toc.Stop();
   if (myid == 0)
   {
      cout << " done, " << tic_toc.RealTime() << "s." << endl;
   }

   // Apply operator matrix
   if (myid == 0)
   {
      cout << "Applying the matrix ..." << flush;
   }
   tic_toc.Clear();
   tic_toc.Start();
   a->Mult(x, b);
   tic_toc.Stop();
   if (myid == 0)
   {
      cout << " done, " << tic_toc.RealTime() << "s." << endl;
   }
   x = 0.0;

   Operator *a_oper = NULL;
   Vector B, X;
   a->FormLinearSystem(ess_tdof_list, x, b, a_oper, X, B);
   HYPRE_Int glob_size = fespace->GlobalTrueVSize();
   if (myid == 0)
   {
       cout << "Size of linear system: " << glob_size << endl;
   }

   // Setup the matrix used for preconditioning
   if (myid == 0)
   {
      cout << "Assembling the preconditioning matrix ..." << flush;
   }
   tic_toc.Clear();
   tic_toc.Start();

   HypreParMatrix *A_pc = NULL;
   if (pc_choice == LOR)
   {
      // TODO: assemble the matrix using the performance code
      ConstantCoefficient one(1.0);
#if PROBLEM == 1
      a_pc->AddDomainIntegrator(new MassIntegrator(one));
#elif PROBLEM == 2
      a_pc->AddDomainIntegrator(new DiffusionIntegrator(one));
#elif PROBLEM == 3
      a_pc->AddDomainIntegrator(new VectorMassIntegrator(one));
#elif PROBLEM == 4
      a_pc->AddDomainIntegrator(new VectorDiffusionIntegrator(one));
#endif
      A_pc = new HypreParMatrix();
      a_pc->UsePrecomputedSparsity();
      a_pc->Assemble();
      a_pc->FormSystemMatrix(ess_tdof_list, *A_pc);        
   }
   else if (pc_choice == HO || pc_choice == JACOBI)
   {
      a_pc->UsePrecomputedSparsity();
      ((HPCBilinearForm*) a)->AssembleBilinearForm(*a_pc);
      A_pc = new HypreParMatrix();
      a_pc->FormSystemMatrix(ess_tdof_list, *A_pc);
   }
   else if (pc_choice == LUMPEDMASS)
   {
      ParGridFunction lumped_mass_diag(fespace);
      a->Mult(ones, lumped_mass_diag);
      HypreParVector* lumped_mass_vec
        = lumped_mass_diag.ParallelAssemble();
      const int local_size = lumped_mass_vec->Size();
      int* I = new int[local_size+1];
      int* J = new int[local_size];
      const int my_col_start = lumped_mass_vec->Partitioning()[0];
      for (int i=0; i<local_size; ++i)
      {
         I[i] = i;
         J[i] = my_col_start + i;
      }
      I[local_size] = local_size;
      A_pc = new HypreParMatrix(lumped_mass_vec->GetComm(),
                                lumped_mass_vec->Size(),
                                lumped_mass_vec->GlobalSize(),
                                lumped_mass_vec->GlobalSize(),
                                I,
                                J,
                                lumped_mass_vec->GetData(),
                                lumped_mass_vec->Partitioning(),
                                lumped_mass_vec->Partitioning());
      delete lumped_mass_vec;
      delete[] I;
      delete[] J;
   }
   tic_toc.Stop();
   if (myid == 0)
   {
      cout << " done, " << tic_toc.RealTime() << "s." << endl;
   }

   // Solve with CG or PCG, depending if the matrix A_pc is available
   CGSolver *pcg;
   pcg = new CGSolver(MPI_COMM_WORLD);
   pcg->SetRelTol(1e-6);
   pcg->SetMaxIter(1000);
   pcg->SetPrintLevel(1);

   HypreSolver* pc_oper = NULL;
   pcg->SetOperator(*a_oper);
   if (pc_choice == HO || pc_choice == LOR)
   {
      pc_oper = new HypreBoomerAMG(*A_pc);
      pcg->SetPreconditioner(*pc_oper);
   }
   if (pc_choice == JACOBI || pc_choice == LUMPEDMASS)
   {
      pc_oper = new HypreDiagScale(*A_pc);
      pcg->SetPreconditioner(*pc_oper);
   }

   tic_toc.Clear();
   tic_toc.Start();

   pcg->Mult(B, X);

   tic_toc.Stop();
   delete pc_oper;

   if (myid == 0)
   {
      // Note: In the pcg algorithm, the number of operator Mult() calls is
      //       N_iter and the number of preconditioner Mult() calls is N_iter+1.
      cout << "Time in CG: " << tic_toc.RealTime() << "s." << endl;
      cout << "Time per CG step: "
           << tic_toc.RealTime() / pcg->GetNumIterations() << "s." << endl;
   }

   // Check relative error in solution
   a->RecoverFEMSolution(X, b, x);
#if 0
#if PROBLEM == 2
   x -= (x * ones) / x.Size();
#elif PROBLEM == 4
   for (int d=0; d<dim; ++d)
   {
      const int ndofs = fespace->GetNDofs();
      double mean = 0;
      for (int i=d*ndofs; i<(d+1)*ndofs; ++i)
      {
         mean += x(i);
      }
      MPI_Allreduce(MPI_IN_PLACE, &mean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      mean /= (x.Size() / dim);
      for (int i=d*ndofs; i<(d+1)*ndofs; ++i)
      {
         x(i) -= mean;
      }
   }
#endif
   ConstantCoefficient zero(0.0);
   const double norm_x = x0.ComputeL2Error(zero);
   x -= x0;
   const double norm_err = x.ComputeL2Error(zero);
   if (myid == 0)
   {
      cout << "|| x - x0 ||_2 / || x0 ||_2 = " << norm_err << "/" << norm_x << "=" << norm_err/norm_x << endl;
   }
#endif

   // 17. Send the solution by socket to a GLVis server.
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock << "parallel " << num_procs << " " << myid << "\n";
      sol_sock.precision(8);
      sol_sock << "solution\n" << *pmesh << x << "\n";
      sol_sock << "keys maaAcvvv" << endl;
   }

   // 18. Free the used memory.
   delete a;
   if (A_pc) { delete A_pc; }
   if (a_pc) { delete a_pc; }
   delete fespace;
   delete fespace_lor;
   delete fec_lor;
   delete pmesh_lor;
   delete fec;
   delete pmesh;
   delete[] partitioning;
   delete pcg;

   MPI_Finalize();

   return 0;
}
