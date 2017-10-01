//                             MFEM Bake-Off Problems
//
// Compile with: See bps.sh
//
// Sample runs:  See bps.sh
//
// Description: This example code demonstrates how to solve the
//              bake-off problems using the proposed operator-based
//              interface proposed in the pa-oper-dev and associated
//              branches. These tests are probably not very efficient
//              at the moment, but the purpose of creating this is to
//              build on top of this and improve performance over time
//              in the various implementations.
//

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{
   MPI_Session mpi(argc, argv);

   // Parse command-line options.
   const char *mesh_file = "inline-hex.mesh";
   int bakeoff = 1;
   int ser_ref_levels = -1;
   int par_ref_levels = 1;
   int order = 1;
   const char *basis_type = "G"; // Gauss-Lobatto
   bool matrix_free = true;
   int ir_order = -1;
   int max_iter = 500;
   int print_level = 0;
   bool write_solution = false;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&bakeoff, "-bp", "--bakeoff",
                  "Bake-off problem to solve.");
   args.AddOption(&ser_ref_levels, "-rs", "--refine-serial",
                  "Number of serial refinements.");
   args.AddOption(&par_ref_levels, "-rp", "--refine-parallel",
                  "Number of parallel refinements.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree) or -1 for"
                  " isoparametric space.");
   args.AddOption(&basis_type, "-b", "--basis-type",
                  "Basis: G - Gauss-Lobatto, P - Positive, U - Uniform");
   args.AddOption(&matrix_free, "-mf", "--matrix-free", "-asm", "--assembly",
                  "Use matrix-free evaluation or efficient matrix assembly.");
   args.AddOption(&ir_order, "-qo", "--quadrature-order",
                  "Quadrature order (2q-1 = 2p+3 by default if negative)");
   args.AddOption(&max_iter, "-mi", "--max-iter",
                  "Maximum number of iterations.");
   args.AddOption(&print_level, "-pl", "--print-level",
                  "Print level for CG.");
   args.AddOption(&write_solution, "-ws", "--write-solution",
                  "-no-ws", "--no-write-solution",
                  "Enable or disable solution and mesh output.");

   args.Parse();
   if (!args.Good())
   {
      if (mpi.Root()) { args.PrintUsage(cout); }
      return 1;
   }
   if (mpi.Root()) { args.PrintOptions(cout); }

   // See class BasisType in fem/fe_coll.hpp for available basis types
   int basis = BasisType::GetType(basis_type[0]);
   if (mpi.Root()) { cout << "Using " << BasisType::Name(basis) << " basis ..." << endl; }

   // Read the serial mesh from the given mesh file on all processors.
   // Refine the mesh in serial to increase the resolution.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   for (int lev = 0; lev < ser_ref_levels; lev++) { mesh->UniformRefinement(); }

   // Parallel partitioning of the mesh.
   ParMesh *pmesh = NULL;
   const int dim = mesh->Dimension();
   const int num_tasks = mpi.WorldSize();
   const int partitions = floor(pow(num_tasks, 1.0 / dim) + 1e-2);
   int *nxyz = new int[dim];
   int product = 1;
   for (int d = 0; d < dim; d++)
   {
      nxyz[d] = partitions;
      product *= partitions;
   }
   if (product == num_tasks)
   {
      int *partitioning = mesh->CartesianPartitioning(nxyz);
      pmesh = new ParMesh(MPI_COMM_WORLD, *mesh, partitioning);
      delete partitioning;
   }
   else
   {
      if (mpi.Root())
      {
         cout << "Non-Cartesian partitioning through METIS will be used.\n";
#ifndef MFEM_USE_METIS
         cout << "MFEM was built without METIS. "
              << "Adjust the number of tasks to use a Cartesian split." << endl;
#endif
      }
#ifndef MFEM_USE_METIS
      return 1;
#endif
      pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   }
   delete [] nxyz;
   delete mesh;
   for (int l = 0; l < par_ref_levels; l++) { pmesh->UniformRefinement(); }

   pmesh->PrintInfo(cout);

   // Define a finite element space on the mesh.
   H1_FECollection fec(order, dim);
   ParFiniteElementSpace fespace(pmesh, &fec, (!bakeoff%2) ? dim : 1);
   const HYPRE_Int size = fespace.GlobalTrueVSize();

   if (mpi.Root()) { cout << "Number of finite element unknowns: " << size << endl; }

   // Determine the list of true (i.e. conforming) essential boundary dofs.
   Array<int> ess_tdof_list;
   if (pmesh->bdr_attributes.Size())
   {
      Array<int> ess_bdr(pmesh->bdr_attributes.Max());
      ess_bdr = 1;
      fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }

   // Set up the linear form b(.) which corresponds to the right-hand
   // side of the FEM linear system, which in this case is (1,phi_i)
   // where phi_i are the basis functions in the finite element
   // fespace.
   LinearForm b(&fespace);
   ConstantCoefficient one(1.0);
   b.AddDomainIntegrator(new DomainLFIntegrator(one));
   b.Assemble();

   // Define the solution vector x as a finite element grid function
   // corresponding to fespace. Initialize x with initial guess of
   // zero, which satisfies the boundary conditions.
   ParGridFunction x(&fespace);
   x = 0.0;

   // If the quadrature order was not passed, set the default defined
   // for all bake-off problems.
   if (ir_order < 0) ir_order = 2 * order + 3;

   if (mpi.Root())
   {
      cout << "Using integration rule with "
           << std::pow((ir_order + 1)/2, dim) << " points ..." << endl;
   }

   // Set up the bilinear form a(.,.) on the finite element space
   // corresponding to the bake-off problem definition.
   if (mpi.Root()) { cout << "Assembling the local matrix ..." << flush; }

   BilinearForm *a = NULL;

#ifdef USE_MPI_WTIME
   double my_rt_start = MPI_Wtime();
#else
   tic_toc.Clear();
   tic_toc.Start();
#endif
   if (matrix_free)
   {
      BilinearFormOperator *a_oper = new BilinearFormOperator(&fespace);
      BilinearFormIntegrator *integ = NULL;
      if (bakeoff/2)
      {
         integ = new PADiffusionIntegrator(&fespace, ir_order);
      }
      else
      {
         integ = new PAMassIntegrator(&fespace, ir_order);
      }
      a_oper->AddDomainIntegrator(integ);
      a = a_oper;
   }
   else
   {
      BilinearForm *a_mat = new BilinearForm(&fespace);
      BilinearFormIntegrator *integ = NULL;
      switch (bakeoff)
      {
      case 1:
         integ = new MassIntegrator(); break;
      case 2:
         integ = new VectorMassIntegrator(); break;
      case 3:
         integ = new DiffusionIntegrator(); break;
      case 4:
         integ = new VectorDiffusionIntegrator(); break;
      default:
         mfem_error("Undefined bake-off problem");
      }
      integ->SetIntRule(
         &IntRules.Get(fespace.GetFE(0)->GetGeomType(), ir_order)
         );
      a_mat->AddDomainIntegrator(integ);
      a_mat->Assemble();
      a = a_mat;
   }
#ifdef USE_MPI_WTIME
   double rt_min, rt_max, my_rt;
   my_rt = MPI_Wtime() - my_rt_start;
#else
   tic_toc.Stop();
   double rt_min, rt_max, my_rt;
   my_rt = tic_toc.RealTime();
#endif
   MPI_Reduce(&my_rt, &rt_min, 1, MPI_DOUBLE, MPI_MIN, 0, pmesh->GetComm());
   MPI_Reduce(&my_rt, &rt_max, 1, MPI_DOUBLE, MPI_MAX, 0, pmesh->GetComm());
   if (mpi.Root())
   {
      cout << " done, " << rt_max << " (" << rt_min << ") s." << endl;
      cout << "\n\"DOFs/sec\" in assembly: "
           << 1e-6*size/rt_max << " ("
           << 1e-6*size/rt_min << ") million.\n" << endl;
   }

   SparseMatrix A_sp;
   Operator *A;
   Vector B, X;
   if (mpi.Root())
   {
      cout << "FormLinearSystem() ..." << endl;
   }
#ifdef USE_MPI_WTIME
   my_rt_start = MPI_Wtime();
#else
   tic_toc.Clear();
   tic_toc.Start();
#endif

   if (matrix_free)
   {
      BilinearFormOperator *a_oper = dynamic_cast<BilinearFormOperator *>(a);
      a_oper->FormLinearSystem(ess_tdof_list, x, b, A, X, B);
   }
   else
   {
      a->FormLinearSystem(ess_tdof_list, x, b, A_sp, X, B);
      A = &A_sp;
   }

#ifdef USE_MPI_WTIME
   my_rt = MPI_Wtime() - my_rt_start;
#else
   tic_toc.Stop();
   my_rt = tic_toc.RealTime();
#endif
   MPI_Reduce(&my_rt, &rt_min, 1, MPI_DOUBLE, MPI_MIN, 0, pmesh->GetComm());
   MPI_Reduce(&my_rt, &rt_max, 1, MPI_DOUBLE, MPI_MAX, 0, pmesh->GetComm());
   if (mpi.Root())
   {
      cout << "FormLinearSystem() ... done, " << rt_max << " (" << rt_min
           << ") s." << endl;
      cout << "\n\"DOFs/sec\" in FormLinearSystem(): "
           << 1e-6*size/rt_max << " ("
           << 1e-6*size/rt_min << ") million.\n" << endl;
   }

   if (mpi.Root()) { cout << "Size of linear system: " << A->Height() << endl; }

   // Solve for solution X.
   CGSolver cg(pmesh->GetComm());
   cg.SetRelTol(1e-6);
   cg.SetMaxIter(max_iter);
   cg.SetPrintLevel(print_level);
   cg.SetOperator(*A);

#ifdef USE_MPI_WTIME
   my_rt_start = MPI_Wtime();
#else
   tic_toc.Clear();
   tic_toc.Start();
#endif

   cg.Mult(B, X);

#ifdef USE_MPI_WTIME
   my_rt = MPI_Wtime() - my_rt_start;
#else
   tic_toc.Stop();
   my_rt = tic_toc.RealTime();
#endif
   MPI_Reduce(&my_rt, &rt_min, 1, MPI_DOUBLE, MPI_MIN, 0, pmesh->GetComm());
   MPI_Reduce(&my_rt, &rt_max, 1, MPI_DOUBLE, MPI_MAX, 0, pmesh->GetComm());
   if (mpi.Root())
   {
      // Note: In the pcg algorithm, the number of operator Mult() calls is
      //       N_iter and the number of preconditioner Mult() calls is N_iter+1.
      cout << "Total CG time:    " << rt_max << " (" << rt_min << ") sec."
           << endl;
      cout << "Time per CG step: "
           << rt_max / cg.GetNumIterations() << " ("
           << rt_min / cg.GetNumIterations() << ") sec." << endl;
      cout << "\n\"DOFs/sec\" in CG: "
           << 1e-6*size*cg.GetNumIterations()/rt_max << " ("
           << 1e-6*size*cg.GetNumIterations()/rt_min << ") million.\n"
           << endl;
   }

   // Recover the solution X as a finite element grid function x.
   a->RecoverFEMSolution(X, b, x);

   // Save the refined mesh and the solution (plot with GLVis).
   if (write_solution)
   {
      ostringstream mesh_name, sol_name;
      mesh_name << "mesh." << setfill('0') << setw(6) << mpi.WorldRank();
      sol_name << "sol." << setfill('0') << setw(6) << mpi.WorldRank();

      ofstream mesh_ofs(mesh_name.str().c_str());
      mesh_ofs.precision(8);
      mesh_ofs << pmesh;

      ofstream sol_ofs(sol_name.str().c_str());
      sol_ofs.precision(8);
      sol_ofs << x;
   }

   // Free the used memory.
   delete pmesh;

   return 0;
}
