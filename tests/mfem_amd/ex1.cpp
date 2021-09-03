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

//
//   CEED Bake-off PROBLEM, based on MFEM Example 1p
//

#include "mfem.hpp"
#include "mfem/linalg/solvers_device.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

Mesh *make_mesh(int myid, int num_procs, int dim, int level,
                int &par_ref_levels, Array<int> &nxyz);

template <typename INTEGRATOR>
void bk2_vector_pa_integrator(int order, Mesh *mesh)
{
   constexpr int dim = 3;

   H1_FECollection fec(order, dim);
   FiniteElementSpace fes(mesh, &fec, dim);

   GridFunction x(&fes), y_fa(&fes), y_pa(&fes);
   x.Randomize(1);

   BilinearForm blf_fa(&fes);
   blf_fa.AddDomainIntegrator(new INTEGRATOR);
   blf_fa.Assemble();
   blf_fa.Finalize();
   blf_fa.Mult(x, y_fa);

   BilinearForm blf_pa(&fes);
   blf_pa.SetAssemblyLevel(AssemblyLevel::PARTIAL);
   blf_pa.AddDomainIntegrator(new INTEGRATOR);
   blf_pa.Assemble();
   blf_pa.Mult(x, y_pa);

   y_pa -= y_fa;
   const double dot = y_pa*y_pa;
   const double epsilon = numeric_limits<double>::epsilon();
   MFEM_VERIFY(fabs(dot) < 10.*epsilon, "Error dot: "<<dot);

   const int dofs = fes.GetVSize();
   constexpr int iter = 64;

   tic_toc.Clear();
   for (int i=0; i<iter; i++)
   {
      MFEM_DEVICE_SYNC;
      tic_toc.Start();
      blf_pa.Mult(x, y_pa);
      MFEM_DEVICE_SYNC;
      tic_toc.Stop();
   }
   const double real_time = tic_toc.RealTime();
   const double mdofs = ((1e-6 * dofs) * iter) / real_time;
   //mfem::out << "\033[38;5;87m";
   //mfem::out << "[VectorMassIntegrator] ";
   mfem::out << "\"DOFs/sec\" in CG: " << mdofs << " million.";
   //mfem::out << "\033[m" << std::endl;
}

int main(int argc, char *argv[])
{
   // 1. Initialize MPI.
   const int num_procs = 1, myid = 0;

   // 2. Parse command-line options.
   int ni = 32;
   int dim = 3;
   int level = 0;
   int order = 1;
   int problem = 0;
   const char *device_config = "cpu";

   OptionsParser args(argc, argv);
   args.AddOption(&ni, "-i", "--ni", "Number of outer iterations.");
   args.AddOption(&dim, "-dim", "--mesh-dimension",
                  "Solve 2D or 3D problem.");
   args.AddOption(&level, "-l", "--refinement-level",
                  "Set the problem size: 2^level mesh elements per processor.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&problem, "-p", "--problem", "Problem 0:Mass, 1:Diffusion.");
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   // 2. Enable hardware devices such as GPUs, and programming models such as
   //    CUDA, OCCA, RAJA and OpenMP based on command line options.
   Device device(device_config);
   device.Print();

   // 4. Read the (serial) mesh from the given mesh file on all processors.  We
   //    can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
   //    and volume meshes with the same code.
   int par_ref_levels;
   Array<int> nxyz;
   Mesh *mesh = make_mesh(myid, num_procs, dim, level, par_ref_levels, nxyz);

   // 6. Define a parallel mesh by a partitioning of the serial mesh. Refine
   //    this mesh further in parallel to increase the resolution. Once the
   //    parallel mesh is defined, the serial mesh can be deleted.
   for (int l = 0; l < par_ref_levels; l++)
   {
      mesh->UniformRefinement();
   }

   long global_ne = mesh->ReduceInt(mesh->GetNE());
   cout << "Total number of elements: " << global_ne << endl;

   if (problem==1)
   {
      bk2_vector_pa_integrator<VectorMassIntegrator>(order,mesh);
      return 0;
   }

   // 7. Define a parallel finite element space on the parallel mesh. Here we
   //    use continuous Lagrange finite elements of the specified order. If
   //    order < 1, we instead use an isoparametric/isogeometric space.
   FiniteElementCollection *fec;
   MFEM_VERIFY(order > 0, "invalid 'order': " << order);
   fec = new H1_FECollection(order, dim);
   FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);
   int size = fespace->GetTrueVSize();
   cout << "Number of finite element unknowns: " << size << endl;

   // 8. Determine the list of true (i.e. parallel conforming) essential
   //    boundary dofs. In this example, the boundary conditions are defined
   //    by marking all the boundary attributes from the mesh as essential
   //    (Dirichlet) and converting them to a list of true dofs.
   Array<int> ess_tdof_list;
   if (mesh->bdr_attributes.Size())
   {
      Array<int> ess_bdr(mesh->bdr_attributes.Max());
      ess_bdr = 1;
      fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }

   // 9. Set up the parallel linear form b(.) which corresponds to the
   //    right-hand side of the FEM linear system, which in this case is
   //    (1,phi_i) where phi_i are the basis functions in fespace.
   LinearForm *b = new LinearForm(fespace);
   b->SetAssemblyLevel(LinearAssemblyLevel::PARTIAL);
   ConstantCoefficient one(1.0);
   b->AddDomainIntegrator(new DomainLFIntegrator(one));
   b->Assemble();

   // 10. Define the solution vector x as a parallel finite element grid function
   //     corresponding to fespace. Initialize x with initial guess of zero,
   //     which satisfies the boundary conditions.
   GridFunction x(fespace);
   x = 0.0;

   // 11. Set up the parallel bilinear form a(.,.) on the finite element space
   //     corresponding to the Laplacian operator -Delta, by adding the Diffusion
   //     domain integrator.
   BilinearForm *a = new BilinearForm(fespace);
   a->SetAssemblyLevel(AssemblyLevel::PARTIAL);
   assert(problem == 0);
   MassIntegrator *mass = new MassIntegrator(one);
   a->AddDomainIntegrator(mass);

   // 12. Assemble the parallel bilinear form and the corresponding linear
   //     system, applying any necessary transformations such as: parallel
   //     assembly, eliminating boundary conditions, applying conforming
   //     constraints for non-conforming AMR, static condensation, etc.
   a->Assemble();

   OperatorPtr A;
   Vector B, X;
   a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);

   // 13. Solve the linear system A X = B.
   //     * With full assembly, use the BoomerAMG preconditioner from hypre.
   //     * With partial assembly, use no preconditioner, for now.
   const int max_cg_iter = 200;
   const int cg_print_level = -1;
   const double rtol = 1e-12;
   const double atol = 0.0;

   CGSolver cg;
   cg.SetRelTol(rtol);
   cg.SetAbsTol(sqrt(atol));
   cg.SetOperator(*A);
   cg.iterative_mode = false;

   // Warm-up CG solve (in case of JIT to avoid timing it)
   {
      Vector Y(X);
      cg.SetMaxIter(2);
      cg.SetPrintLevel(-1);
      cg.Mult(B, Y);
      MFEM_DEVICE_SYNC;
   }

   cg.SetMaxIter(max_cg_iter);
   cg.SetPrintLevel(cg_print_level);

   const int dofs = A->Height();
   double max_mdofs = 0.0;

   const int kernel_version = Device::KernelsVersion();

   for (int i=0; i<ni; i++)
   {
      int cg_iter = 1;
      double real_time = 0.0;

      if (kernel_version != 5)
      {
         tic_toc.Clear();
         // Start & Stop CG timing.
         MFEM_DEVICE_SYNC;
         tic_toc.Start();
         cg.Mult(B, X);
         MFEM_DEVICE_SYNC;
         tic_toc.Stop();

         real_time = tic_toc.RealTime();
         cg_iter = cg.GetNumIterations();
      }
      else // kernel_version == 5
      {
         // benchmark this problem through DeviceCG
         DeviceCG(*fespace, *A, B, X, ess_tdof_list, mass->GetPAData(),
                  cg_iter, real_time, cg_print_level, max_cg_iter, rtol, atol);
      }
      const double mdofs = ((1e-6 * dofs) * cg_iter) / real_time;
      /*mfem::out << "\033[33m";
      mfem::out << "\"DOFs/sec\" in CG: " << mdofs << " million.";
      mfem::out << "\033[m" << std::endl;*/
      max_mdofs = fmax(max_mdofs, mdofs);
   }
   mfem::out << "\033[32m";
   mfem::out << "\"DOFs/sec\" in CG: " << max_mdofs << " million.";
   mfem::out << "\033[m" << std::endl;

   // 17. Free the used memory.
   delete a;
   delete b;
   delete fespace;
   delete fec;

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
