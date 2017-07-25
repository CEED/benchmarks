// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
// reserved. See files LICENSE and NOTICE for details.
//
// This file is part of CEED, a collection of benchmarks, miniapps, software
// libraries and APIs for efficient high-order finite element and spectral
// element discretizations for exascale applications. For more information and
// source code availability see http://github.com/ceed.
//
// The CEED research is supported by the Exascale Computing Project
// (17-SC-20-SC), a collaborative effort of two U.S. Department of Energy
// organizations (Office of Science and the National Nuclear Security
// Administration) responsible for the planning and preparation of a capable
// exascale ecosystem, including software, applications, hardware, advanced
// system engineering and early testbed platforms, in support of the nation's
// exascale computing imperative.

//==============================================================================
//                  MFEM Bake-off Problems 1, 2, 3, and 4
//                                Version 2
//
// Compile with: see ../../README.md
//
// Sample runs:  see ../../README.md
//
// Description:  These benchmarks (CEED Bake-off Problems BP1-4) test the
//               performance of high-order mass and stiffness matrix operator
//               evaluation with "partial assembly" algorithms.
//
//               Code is based on MFEM's HPC ex1, http://mfem.org/performance.
//
//               More details about CEED's bake-off problems can be found at
//               http://ceed.exascaleproject.org/bps.
//==============================================================================

#include "mfem-performance.hpp"
#include <fstream>
#include <iostream>
#include <vector>
#include <list>

using namespace std;
using namespace mfem;

#ifndef PROBLEM
#define PROBLEM 1
#endif

#ifndef GEOM
#define GEOM Geometry::CUBE
#endif

#ifndef MESH_P
#define MESH_P 1
#endif

#ifndef SOL_P
#define SOL_P 3
#endif

#ifndef IR_ORDER
#define IR_ORDER 2*SOL_P+3
#endif

// Define template parameters for optimized build.
const Geometry::Type geom     = GEOM;
const int            mesh_p   = MESH_P;
const int            sol_p    = SOL_P;
const int            ir_order = IR_ORDER;
const int            dim      = Geometry::Constants<geom>::Dimension;
const bool           vec      = PROBLEM == 2 || PROBLEM == 4;

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
typedef TBilinearForm<mesh_t,sol_fes_t,int_rule_t,mass_integ_t,vec_layout_t> HPCBilinearForm;
#elif PROBLEM == 3
typedef TBilinearForm<mesh_t,sol_fes_t,int_rule_t,diffusion_integ_t,scal_layout_t> HPCBilinearForm;
#elif PROBLEM == 4
typedef TBilinearForm<mesh_t,sol_fes_t,int_rule_t,diffusion_integ_t,vec_layout_t> HPCBilinearForm;
#else
#error "Invalid bake-off problem."
#endif

// Naive factorization of number into roughly balanced factors
vector<int> balanced_factorization(int number, int num_factors) {

   // Return trivial factorization if possible
   if (num_factors <= 1)
   {
      return vector<int>(1, number);
   }
   if (number == 0)
   {
      return vector<int>(num_factors, 0);
   }
   if (number == 1)
   {
      return vector<int>(num_factors, 1);
   }
   if (number == -1)
   {
      vector<int> f(num_factors, 1);
      f[0] = -1;
      return f;
   }

   // Initial list of factors
   list<int> factors(num_factors-1, 1);
   factors.push_back(number > 0 ? number : -number);

   // Iterate algorithm until list converges
   for(int iter = 0; iter < number; ++iter)
   {

      // Smallest and largest factors
      int smallest = factors.front();
      int largest = factors.back();

      // Find smallest prime factor of largest factor
      int p = 2;
      while (largest % p != 0) { ++p; }

      // Transfer prime factor from largest to smallest list entry
      if (largest / p <= smallest) { break; }
      else
      {
         factors.pop_front();
         factors.pop_back();
         factors.insert(lower_bound(factors.begin(),
                                    factors.end(),
                                    smallest * p),
                        smallest * p);
         factors.insert(lower_bound(factors.begin(),
                                    factors.end(),
                                    largest / p),
                        largest / p);
      }
      
   }

   // Negate first factor is initial number is negative
   if (number < 0) { factors.front() *= -1; }

   // Return factorization
   return vector<int>(factors.begin(), factors.end());

}

int main(int argc, char *argv[])
{
   // Initialize MPI.
   int num_procs, myid;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);

   // Initialize timers
   double my_rt_start, my_rt, rt_min, rt_max;

   // Parse command-line options.
   int el_per_proc = 1;
   double tol = 1e-6;
   int max_iters = 500;
   const char *pc = "none";
   bool essential_bcs = true;
   bool visualization = 1;

   OptionsParser args(argc, argv);
   args.AddOption(&el_per_proc, "-e", "--num-el-per-proc",
                  "Number of elements per MPI rank.");
   args.AddOption(&tol, "-tol", "--pcg-rel-tol",
                  "Relative tolerance for PCG convergence.");
   args.AddOption(&max_iters, "-i", "--max-iters",
                  "Maximum number of PCG iterations.");
   args.AddOption(&pc, "-pc", "--preconditioner",
                  "Preconditioner: "
                  "lor - low-order-refined (matrix-free) AMG, "
                  "ho - high-order (assembled) AMG, "
                  "jacobi, "
                  "lumpedmass, "
                  "none.");
   args.AddOption(&essential_bcs, "-ess-bc", "--essential-bcs",
                  "-nat-bc", "--natural-bcs",
                  "Essential or natural boundary conditions.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
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

   // Factorize number of processes and elements
   vector<int> num_procs_dims = balanced_factorization(num_procs, dim);
   int unrefined_el_per_proc = el_per_proc;
   int refinement_levels = 0;
   int refinement_factor = 1;
   for (int d = 0; d < dim; ++d) { refinement_factor *= 2; }
   while (unrefined_el_per_proc % refinement_factor == 0)
   {
      unrefined_el_per_proc /= refinement_factor;
      refinement_levels++;
   }
   vector<int> unrefined_el_per_proc_dims
     = balanced_factorization(unrefined_el_per_proc, dim);
   reverse(unrefined_el_per_proc_dims.begin(),
           unrefined_el_per_proc_dims.end());

   // Generate serial mesh
   Mesh *mesh;
   switch (dim)
   {
   case 1:
      mesh = new Mesh(num_procs * unrefined_el_per_proc, 1.0);
      break;
   case 2:
      mesh = new Mesh(num_procs_dims[0] * unrefined_el_per_proc_dims[0],
                      num_procs_dims[1] * unrefined_el_per_proc_dims[1],
                      Element::QUADRILATERAL, 1, 1.0, 1.0);
      break;
   case 3:
      mesh = new Mesh(num_procs_dims[0] * unrefined_el_per_proc_dims[0],
                      num_procs_dims[1] * unrefined_el_per_proc_dims[1],
                      num_procs_dims[2] * unrefined_el_per_proc_dims[2],
                      Element::HEXAHEDRON, 1, 1.0, 1.0, 1.0);
      break;
   default:
      mfem_error("Invalid number of dimensions");
      return -1;
   }

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

   // Define a parallel mesh by a partitioning of the serial mesh.
   // Once the parallel mesh is defined, the serial mesh can be
   // deleted.
   if (myid == 0)
   {
      cout << "Initializing parallel mesh ..." << endl;
   }
   int *partitioning = mesh->CartesianPartitioning(num_procs_dims.data());
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh, partitioning);
   delete[] partitioning;
   delete mesh;
   for (int l = 0; l < refinement_levels; ++l)
   {
      if (myid == 0)
      {
         cout << "Parallel refinement: level " << l << " -> level " << l+1
              << " ..." << flush;
      }
      pmesh->UniformRefinement();
      MPI_Barrier(MPI_COMM_WORLD);
      if (myid == 0)
      {
         cout << " done." << endl;
      }
   }
   if (pmesh->MeshGenerator() & 1) // simplex mesh
   {
      MFEM_VERIFY(pc_choice != LOR,
                  "triangle and tet meshes do not support the LOR "
                  "preconditioner yet");
   }

   // Define a parallel finite element space on the parallel mesh
   const int basis = BasisType::GaussLobatto;
   FiniteElementCollection *fec = new H1_FECollection(sol_p, dim, basis);
   ParFiniteElementSpace *fespace
     = new ParFiniteElementSpace(pmesh, fec,
                                 vec ? dim : 1,
                                 vec ? Ordering::byVDIM : Ordering::byNODES);
   HYPRE_Int size = fespace->GlobalTrueVSize();
   int local_size = fespace->GetTrueVSize();
   if (myid == 0)
   {
      cout << "Number of finite element unknowns: " << size << endl;
   }
   ParMesh *pmesh_lor = NULL;
   FiniteElementCollection *fec_lor = NULL;
   ParFiniteElementSpace *fespace_lor = NULL;
   if (pc_choice == LOR)
   {
      pmesh_lor = new ParMesh(pmesh, sol_p, basis);
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
      delete pmesh;
      MPI_Finalize();
      return 5;
   }

   // Determine the list of true (i.e. parallel conforming) essential
   // boundary dofs
   Array<int> ess_tdof_list;
   if (pmesh->bdr_attributes.Size())
   {
      Array<int> ess_bdr(pmesh->bdr_attributes.Max());
      ess_bdr = essential_bcs ? 1 : 0;
      fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }

   // Define the solution vector x and RHS vector b
   // Note: subtract mean value if solving stiffness matrix problem
   ParGridFunction x0(fespace), x(fespace), b(fespace);
   HypreParVector* x0_tdofs = fespace->NewTrueDofVector();
   x0_tdofs->Randomize(0);
   x0_tdofs->SetSubVector(ess_tdof_list, 0.0);
   x0.SetFromTrueDofs(*x0_tdofs);
   delete x0_tdofs;
#if PROBLEM == 3
   if (!essential_bcs)
   {
      double mean = x0.Sum() / size;
      MPI_Allreduce(MPI_IN_PLACE, &mean, 1, MPI_DOUBLE, MPI_SUM, pmesh->GetComm());
      x0 -= mean;
   }
#elif PROBLEM == 4
   if (!essential_bcs)
   {
      const int npts = local_size / dim;
      for (int d=0; d<dim; ++d)
      {
         double mean = 0;
         for (int i=d*npts; i<(d+1)*npts; ++i)
         {
            mean += x0(i);
         }
         MPI_Allreduce(MPI_IN_PLACE, &mean, 1, MPI_DOUBLE, MPI_SUM, pmesh->GetComm());
         mean /= size / dim;
         for (int i=d*npts; i<(d+1)*npts; ++i)
         {
            x0(i) -= mean;
         }
      }
   }
#endif
   x = x0;

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

   // High-performance assembly/evaluation using the templated operator type
   HPCBilinearForm *a = NULL;
   if (myid == 0)
   {
      cout << "Assembling the local matrix ..." << flush;
   }
#ifdef USE_MPI_WTIME
   my_rt_start = MPI_Wtime();
#else
   tic_toc.Clear();
   tic_toc.Start();
#endif
#if PROBLEM == 1 || PROBLEM == 2
   a = new HPCBilinearForm(mass_integ_t(coeff_t(1.0)), *fespace);
#elif PROBLEM == 3 || PROBLEM == 4
   a = new HPCBilinearForm(diffusion_integ_t(coeff_t(1.0)), *fespace);
#endif
   a->Assemble();
#ifdef USE_MPI_WTIME
   my_rt = MPI_Wtime() - my_rt_start;
#else
   tic_toc.Stop();
   my_rt = tic_toc.RealTime();
#endif
   MPI_Reduce(&my_rt, &rt_min, 1, MPI_DOUBLE, MPI_MIN, 0, pmesh->GetComm());
   MPI_Reduce(&my_rt, &rt_max, 1, MPI_DOUBLE, MPI_MAX, 0, pmesh->GetComm());
   if (myid == 0)
   {
      cout << " done, " << rt_max << " (" << rt_min << ") s." << endl;
      cout << "\n\"DOFs/sec\" in local assembly: "
           << 1e-6*size/rt_max << " ("
           << 1e-6*size/rt_min << ") million.\n" << endl;
   }

   // Apply operator matrix
   if (myid == 0)
   {
      cout << "Applying the matrix ..." << flush;
   }
#ifdef USE_MPI_WTIME
   my_rt_start = MPI_Wtime();
#else
   tic_toc.Clear();
   tic_toc.Start();
#endif
   a->Mult(x, b);
#ifdef USE_MPI_WTIME
   my_rt = MPI_Wtime() - my_rt_start;
#else
   tic_toc.Stop();
   my_rt = tic_toc.RealTime();
#endif
   MPI_Reduce(&my_rt, &rt_min, 1, MPI_DOUBLE, MPI_MIN, 0, pmesh->GetComm());
   MPI_Reduce(&my_rt, &rt_max, 1, MPI_DOUBLE, MPI_MAX, 0, pmesh->GetComm());
   if (myid == 0)
   {
      cout << " done, " << rt_max << " (" << rt_min << ") s." << endl;
      cout << "\n\"DOFs/sec\" in matrix multiplication: "
           << 1e-6*size/rt_max << " ("
           << 1e-6*size/rt_min << ") million.\n" << endl;
   }
   x = 0.0;

   if (myid == 0)
   {
      cout << "FormLinearSystem() ..." << flush;
   }
   Operator *a_oper = NULL;
   Vector B, X;
#ifdef USE_MPI_WTIME
   my_rt_start = MPI_Wtime();
#else
   tic_toc.Clear();
   tic_toc.Start();
#endif
   a->FormLinearSystem(ess_tdof_list, x, b, a_oper, X, B);
#ifdef USE_MPI_WTIME
   my_rt = MPI_Wtime() - my_rt_start;
#else
   tic_toc.Stop();
   double rt_min, rt_max, my_rt;
   my_rt = tic_toc.RealTime();
#endif
   MPI_Reduce(&my_rt, &rt_min, 1, MPI_DOUBLE, MPI_MIN, 0, pmesh->GetComm());
   MPI_Reduce(&my_rt, &rt_max, 1, MPI_DOUBLE, MPI_MAX, 0, pmesh->GetComm());
   if (myid == 0)
   {
      cout << " done, " << rt_max << " (" << rt_min << ") s." << endl;
      cout << "\n\"DOFs/sec\" in FormLinearSystem(): "
           << 1e-6*size/rt_max << " ("
           << 1e-6*size/rt_min << ") million.\n" << endl;
   }

   // Setup the matrix used for preconditioning
   if (myid == 0)
   {
      cout << "Assembling the preconditioning matrix ..." << flush;
   }
#ifdef USE_MPI_WTIME
   my_rt_start = MPI_Wtime();
#else
   tic_toc.Clear();
   tic_toc.Start();
#endif

   HypreParMatrix *A_pc = NULL;
   if (pc_choice == LOR)
   {
      // TODO: assemble the matrix using the performance code
      ConstantCoefficient one(1.0);
#if PROBLEM == 1
      a_pc->AddDomainIntegrator(new MassIntegrator(one));
#elif PROBLEM == 2
      a_pc->AddDomainIntegrator(new VectorMassIntegrator(one));
#elif PROBLEM == 3
      a_pc->AddDomainIntegrator(new DiffusionIntegrator(one));
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
      a->AssembleBilinearForm(*a_pc);
      A_pc = new HypreParMatrix();
      a_pc->FormSystemMatrix(ess_tdof_list, *A_pc);
   }
   else if (pc_choice == LUMPEDMASS)
   {
      ParGridFunction lumped_mass_diag(fespace), ones(fespace);
      ones = 1.0;
      a->Mult(ones, lumped_mass_diag);
      HypreParVector* lumped_mass_vec
        = lumped_mass_diag.ParallelAssemble();
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
#ifdef USE_MPI_WTIME
   my_rt = MPI_Wtime() - my_rt_start;
#else
   tic_toc.Stop();
   my_rt = tic_toc.RealTime();
#endif
   MPI_Reduce(&my_rt, &rt_min, 1, MPI_DOUBLE, MPI_MIN, 0, pmesh->GetComm());
   MPI_Reduce(&my_rt, &rt_max, 1, MPI_DOUBLE, MPI_MAX, 0, pmesh->GetComm());
   if (myid == 0)
   {
      cout << " done, " << rt_max << "s." << endl;
   }

   // Solve with CG or PCG, depending if the matrix A_pc is available
   CGSolver *pcg;
   pcg = new CGSolver(pmesh->GetComm());
   pcg->SetRelTol(tol);
   pcg->SetMaxIter(max_iters);
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

#ifdef USE_MPI_WTIME
   my_rt_start = MPI_Wtime();
#else
   tic_toc.Clear();
   tic_toc.Start();
#endif

   pcg->Mult(B, X);

#ifdef USE_MPI_WTIME
   my_rt = MPI_Wtime() - my_rt_start;
#else
   tic_toc.Stop();
   my_rt = tic_toc.RealTime();
#endif
   delete pc_oper;

   MPI_Reduce(&my_rt, &rt_min, 1, MPI_DOUBLE, MPI_MIN, 0, pmesh->GetComm());
   MPI_Reduce(&my_rt, &rt_max, 1, MPI_DOUBLE, MPI_MAX, 0, pmesh->GetComm());
   if (myid == 0)
   {
      // Note: In the pcg algorithm, the number of operator Mult() calls is
      //       N_iter and the number of preconditioner Mult() calls is N_iter+1.
      cout << "Total CG time:    " << rt_max << " (" << rt_min << ") sec."
           << endl;
      cout << "Time per CG step: "
           << rt_max / pcg->GetNumIterations() << " ("
           << rt_min / pcg->GetNumIterations() << ") sec." << endl;
      cout << "\n\"DOFs/sec\" in CG: "
           << 1e-6*size*pcg->GetNumIterations()/rt_max << " ("
           << 1e-6*size*pcg->GetNumIterations()/rt_min << ") million.\n"
           << endl;
   }

   // Check relative error in solution
   a->RecoverFEMSolution(X, b, x);
#if PROBLEM == 3
   if (!essential_bcs)
   {
      double mean = x.Sum() / size;
      MPI_Allreduce(MPI_IN_PLACE, &mean, 1, MPI_DOUBLE, MPI_SUM, pmesh->GetComm());
      x -= mean;
   }
#elif PROBLEM == 4
   if (!essential_bcs)
   {
      const int npts = local_size / dim;
      for (int d=0; d<dim; ++d)
      {
         double mean = 0;
         for (int i=d*npts; i<(d+1)*npts; ++i)
         {
            mean += x(i);
         }
         MPI_Allreduce(MPI_IN_PLACE, &mean, 1, MPI_DOUBLE, MPI_SUM, pmesh->GetComm());
         mean /= size / dim;
         for (int i=d*npts; i<(d+1)*npts; ++i)
         {
            x(i) -= mean;
         }
      }
   }
#endif
   ConstantCoefficient zero_coefficient(0.0);
#if PROBLEM == 2 || PROBLEM == 4
   VectorGridFunctionCoefficient x0_coefficient(&x0);
#else
   GridFunctionCoefficient x0_coefficient(&x0);
#endif
   const double norm_x0 = x0.ComputeL2Error(zero_coefficient);
   const double norm_err = x.ComputeL2Error(x0_coefficient);
   if (myid == 0)
   {
      cout << "Relative error: " << norm_err / norm_x0 << endl;
   }

   // Send the solution by socket to a GLVis server.
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

   // Free the used memory.
   delete a;
   if (A_pc) { delete A_pc; }
   if (a_pc) { delete a_pc; }
   delete fespace;
   delete fespace_lor;
   delete fec_lor;
   delete pmesh_lor;
   delete fec;
   delete pmesh;
   delete pcg;

   MPI_Finalize();

   return 0;
}
