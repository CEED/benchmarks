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

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <memory>

using namespace std;
using namespace mfem;

Mesh *make_mesh(int myid, int num_procs, int dim, int level,
                int &par_ref_levels, Array<int> &nxyz);

// Prescribed periodic pressure
double p_periodic(const Vector &x)
{
   return cos(5.0*M_PI*x(0)) * sin(2.0*M_PI*x(1)); // + cos(3.0*x(2));
}

int main(int argc, char *argv[])
{
   // **************************************************************************
   // Initialize MPI.
   int num_procs, myid;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);

   // **************************************************************************
   // Parse command-line options.
   int dim = 3;
   int level = 0;
   int order = 1;
   const char *device_config = "cpu";

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
   const int vel_order = order;
   const int pres_order = order - 1;

   // **************************************************************************
   // Enable hardware devices such as GPUs, and programming models such as
   // CUDA, OCCA, RAJA and OpenMP based on command line options.
   Device device(device_config);
   if (myid == 0) { device.Print(); }

   // **************************************************************************
   // Read the (serial) mesh from the given mesh file on all processors.  We
   // can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
   // and volume meshes with the same code.
   Array<int> nxyz;
   int par_ref_levels;
   Mesh *mesh = make_mesh(myid, num_procs, dim, level, par_ref_levels, nxyz);
   //mesh->SetCurvature(order, false, -1, Ordering::byNODES);

   // **************************************************************************
   // Define a parallel mesh by a partitioning of the serial mesh. Refine
   // this mesh further in parallel to increase the resolution. Once the
   //  parallel mesh is defined, the serial mesh can be deleted.
   int *partitioning = nxyz.Size() ? mesh->CartesianPartitioning(nxyz) : NULL;
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh, partitioning);
   delete [] partitioning;
   delete mesh;
   for (int l = 0; l < par_ref_levels; l++)
   {
      pmesh->UniformRefinement();
   }
   const long global_ne = pmesh->ReduceInt(pmesh->GetNE());
   if (myid == 0)
   {
      cout << "Total number of elements: " << global_ne << endl;
   }

   // **************************************************************************
   // Define a finite element space on the mesh. Here we use continuous
   // Lagrange finite elements of the specified order.
   // Define vector FE space for velocity and scalar FE space for pressure
   H1_FECollection vel_fec{vel_order, dim};
   H1_FECollection pres_fec{pres_order};
   ParFiniteElementSpace vel_fes{pmesh, &vel_fec, dim};
   ParFiniteElementSpace pres_fes{pmesh, &pres_fec};
   HYPRE_Int vel_size = vel_fes.GlobalTrueVSize();
   HYPRE_Int pres_size = pres_fes.GlobalTrueVSize();
   const HYPRE_Int size = vel_size + pres_size;
   if (myid == 0)
   {
      cout << "Number of finite element unknowns: " << size << endl;
   }
   std::cout << "Number of velocity unknowns: "
             << vel_fes.GetTrueVSize() << std::endl;
   std::cout << "Number of pressure unknowns: "
             << pres_fes.GetTrueVSize() << std::endl;

   // **************************************************************************
   Array<int> ess_trial_tdof_list, ess_test_tdof_list;
   const int nattr = pmesh->bdr_attributes.Max();
   Array<int> ess_bdr(nattr);
   if (nattr >= 1)
   {
      ess_bdr = 0;
      ess_bdr[0] = 1;
      pres_fes.GetEssentialTrueDofs(ess_bdr, ess_trial_tdof_list);
   }
   if (nattr >= 2)
   {
      ess_bdr = 0;
      ess_bdr[1] = 1;
      vel_fes.GetEssentialTrueDofs(ess_bdr, ess_test_tdof_list);
   }

   FunctionCoefficient pcoeff{p_periodic};
   ParGridFunction vel{&vel_fes};
   ParGridFunction pressure{&pres_fes};
   pressure.Randomize(3);

   Vector v_tdofs(vel_fes.GetTrueVSize());
   Vector p_tdofs(pres_fes.GetTrueVSize());
   vel.ParallelProject(v_tdofs);
   pressure.ParallelProject(p_tdofs);

   // **************************************************************************
   // Set up the bilinear form d(.,.) on the finite element space
   // corresponding to the vector divergence form.
   ParMixedBilinearForm g(&pres_fes, &vel_fes);
   g.SetAssemblyLevel(AssemblyLevel::PARTIAL);
   g.AddDomainIntegrator(new GradientIntegrator);
   g.Assemble();
   OperatorHandle G;
   g.FormRectangularSystemMatrix(ess_trial_tdof_list, ess_test_tdof_list, G);
   G->Mult(p_tdofs, v_tdofs);
   //vel.Distribute(v_tdofs);

   // **************************************************************************
#if 0
   ParGridFunction vel2 {&vel_fes};
   {
      // Also do non-pa version
      ParMixedBilinearForm gform2(&pres_fes, &vel_fes);
      gform2.AddDomainIntegrator(new GradientIntegrator);
      gform2.Assemble();
      OperatorHandle G2;
      gform2.FormRectangularSystemMatrix(ess_trial_tdof_list, ess_test_tdof_list, G2);
      G2->Mult(p_tdofs, v_tdofs);
      vel2.Distribute(v_tdofs);
   }
   // Setup all integration rules for any element type
   int order_quad = max(2, 2*order+1);
   const IntegrationRule *irs[Geometry::NumGeom];
   for (int i=0; i < Geometry::NumGeom; ++i)
   {
      irs[i] = &(IntRules.Get(i, order_quad));
   }
   ConstantCoefficient zero{0.0};

   double norm_u = vel.ComputeL2Error(zero, irs);
   std::cout << "||G*p_{test}|| = " << norm_u << std::endl;
   MPI_Barrier(MPI_COMM_WORLD);
   vel -= vel2;
   vel.HostReadWrite();
   double error_v = GlobalLpNorm(infinity(), vel.Normlinf(), MPI_COMM_WORLD);
   std::cout << "||G*p_{pa} - G*p_{non-pa}|| = " << error_v << std::endl;
   std::cout << "Setup time = " << assembly_time << std::endl;
   std::cout << "Mult time = " << mult_time << std::endl;
#endif

   // **************************************************************************
   tic_toc.Clear();
   tic_toc.Start();
   constexpr int N = 2000;
   for (int k=0; k<N; k++)
   {
      G->Mult(p_tdofs, v_tdofs);
      MFEM_DEVICE_SYNC;
   }
   tic_toc.Stop();

   // **************************************************************************
   double rt_min, rt_max, my_rt;
   my_rt = tic_toc.RealTime();
   MPI_Reduce(&my_rt, &rt_min, 1, MPI_DOUBLE, MPI_MIN, 0, pmesh->GetComm());
   MPI_Reduce(&my_rt, &rt_max, 1, MPI_DOUBLE, MPI_MAX, 0, pmesh->GetComm());
   if (myid == 0)
   {
      int cg_iter = N;
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

   // **************************************************************************
   delete pmesh;
   MPI_Finalize();
   return 0;
}

// *****************************************************************************
Mesh *make_mesh(int myid, int num_procs, int dim, int level,
                int &par_ref_levels, Array<int> &nxyz)
{
   int log_p = (int)floor(log((double)num_procs)/log(2.0) + 0.5);
   MFEM_VERIFY((1 << log_p) == num_procs,
               "number of processor is not a power of 2: " << num_procs)
   MFEM_VERIFY(dim == 3, "dim = " << dim << " is NOT implemented!")

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
