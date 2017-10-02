//                                MFEM Example 0
//
// Compile with: make ex0
//
// Sample runs:  ex0 -m ../../data/square-disc.mesh
//               ex0 -m ../../data/star.mesh
//               ex0 -m ../../data/escher.mesh
//               ex0 -m ../../data/fichera.mesh
//               ex0 -m ../../data/square-disc-p2.vtk -o 2
//               ex0 -m ../../data/square-disc-p3.mesh -o 3
//               ex0 -m ../../data/square-disc-nurbs.mesh -o -1
//               ex0 -m ../../data/disc-nurbs.mesh -o -1
//               ex0 -m ../../data/pipe-nurbs.mesh -o -1
//               ex0 -m ../../data/star-surf.mesh
//               ex0 -m ../../data/square-disc-surf.mesh
//               ex0 -m ../../data/inline-segment.mesh
//               ex0 -m ../../data/amr-quad.mesh
//               ex0 -m ../../data/amr-hex.mesh
//               ex0 -m ../../data/fichera-amr.mesh
//               ex0 -m ../../data/mobius-strip.mesh
//               ex0 -m ../../data/mobius-strip.mesh -o -1 -sc
//
// Description:  This example code demonstrates the use of MFEM to define a
//               simple finite element discretization of the Laplace problem
//               -Delta u = 1 with homogeneous Dirichlet boundary conditions.
//               Specifically, we discretize using a FE space of the specified
//               order, or if order < 1 using an isoparametric/isogeometric
//               space (i.e. quadratic for quadratic curvilinear mesh, NURBS for
//               NURBS mesh, etc.)
//
//               The example highlights the use of mesh refinement, finite
//               element grid functions, as well as linear and bilinear forms
//               corresponding to the left-hand side and right-hand side of the
//               discrete linear system. We also cover the explicit elimination
//               of essential boundary conditions, static condensation, and the
//               optional connection to the GLVis tool for visualization.

#include <fstream>
#include <iostream>
#include <cmath>

#include "mfem.hpp"
#include "occa.hpp"

using namespace std;
using namespace mfem;

#ifndef MFEM_USE_ACROTENSOR
typedef OccaMassIntegrator AcroMassIntegrator;
#endif

double InitialCondition(const Vector& x)
{
  double r2 = 0;
  for (int i = 0; i < x.Size(); i++)
    r2 += (x(i)-0.5)*(x(i)-0.5);
  return std::exp(-r2);
}

void VectorInitialCondition(const Vector& x, Vector& v)
{
  double r2 = 0;
  for (int i = 0; i < x.Size(); i++)
    r2 += (x(i)-0.5)*(x(i)-0.5);
  v(0) = std::exp(-r2);
  if (x.Size() > 1)
    v(1) = std::exp(-r2);
  if (x.Size() > 2)
    v(2) = std::exp(r2);
}


int main(int argc, char *argv[])
{
  // 1. Parse command-line options.
  const char *mesh_file = "../../data/inline-hex.mesh";
  const char *basis_type = "G";
  int order = 3;
  int ref_levels = 2;
  int cg_iter = 100;
  int quad_add = 0;
  const char *device_info = "mode: 'Serial'";
  bool occa_verbose = false;
  bool use_acrotensor = false;
  bool visualization = 0;

  OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh",
                 "Mesh file to use.");
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree) or -1 for"
                 " isoparametric space.");
  args.AddOption(&ref_levels, "-r", "--ref-levels",
                 "Number of uniform mesh refinements"
                 " base mesh.");
  args.AddOption(&cg_iter, "-i", "--cg-iter",
                 "Number of uniform mesh refinements"
                 " base mesh.");
  args.AddOption(&quad_add, "-q", "--quad-add",
                 "Additional quadrature order: 2p+3 + quad_add"
                 " base quadrature: 2p+3.");
  args.AddOption(&device_info, "-d", "--device-info",
                 "Device information to run example on (default: \"mode: 'Serial'\").");
  args.AddOption(&occa_verbose,
                 "-ov", "--occa-verbose",
                 "--no-ov", "--no-occa-verbose",
                 "Print verbose information about OCCA kernel compilation.");
  args.AddOption(&use_acrotensor,
                 "-ac", "--use-acro",
                 "--no-ac", "--no-acro",
                 "Use Acrotensor.");
  args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                 "--no-visualization",
                 "Enable or disable GLVis visualization.");
  args.Parse();
  if (!args.Good())
    {
      args.PrintUsage(cout);
      return 1;
    }
  args.PrintOptions(cout);

#ifndef MFEM_USE_ACROTENSOR
  if (use_acrotensor) {
    cout << "MFEM not compiled with Acrotensor, reverting to OCCA\n";
    use_acrotensor = false;
  }
#endif

  // Set the OCCA device to run example in
  occa::setDevice(device_info);

  // Load cached kernels
  occa::loadKernels();
  occa::loadKernels("mfem");

  // Set as the background device
  occa::settings()["verboseCompilation"] = occa_verbose;

  // See class BasisType in fem/fe_coll.hpp for available basis types
  int basis = BasisType::GetType(basis_type[0]);
  cout << "Using " << BasisType::Name(basis) << " basis ..." << endl;

  // 2. Read the mesh from the given mesh file. We can handle triangular,
  //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
  //    the same code.
  Mesh *mesh = new Mesh(mesh_file, 1, 1);
  int dim = mesh->Dimension();

  // 3. Refine the mesh to increase the resolution. In this example we do
  //    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
  //    largest number that gives a final mesh with no more than 50,000
  //    elements.
  {
    for (int l = 0; l < ref_levels; l++)
      {
        mesh->UniformRefinement();
      }
  }

  // 4. Define a finite element space on the mesh. Here we use continuous
  //    Lagrange finite elements of the specified order. If order < 1, we
  //    instead use an isoparametric/isogeometric space.
#if 1
  FiniteElementCollection *fec = new H1_FECollection(order, dim, basis);
#else
  FiniteElementCollection *fec = new L2_FECollection(order, dim, basis);
#endif

  OccaFiniteElementSpace *ofespace = new OccaFiniteElementSpace(mesh, fec, dim);
  FiniteElementSpace *fespace = ofespace->GetFESpace();
  const int size = fespace->GetTrueVSize();

  cout << "Number of finite element unknowns: "
       << fespace->GetTrueVSize() << endl;

  // 6. Set up the linear form b(.) which corresponds to the right-hand side of
  //    the FEM linear system, which in this case is (1,phi_i) where phi_i are
  //    the basis functions in the finite element fespace.
  LinearForm lf_b(fespace);
#if 0
  ConstantCoefficient u_0(1.0);
#else
  VectorFunctionCoefficient u_0(dim, VectorInitialCondition);
#endif
  lf_b.AddDomainIntegrator(new VectorDomainLFIntegrator(u_0));
  lf_b.Assemble();
  OccaVector b(lf_b);

  // 7. Define the solution vector x as a finite element grid function
  //    corresponding to fespace. Initialize x with initial guess of zero,
  //    which satisfies the boundary conditions.
  OccaGridFunction x(ofespace);
  x = 0.0;

  // 8. Set up the bilinear form a(.,.) on the finite element space
  //    corresponding to the Laplacian operator -Delta, by adding the Diffusion
  //    domain integrator.
  cout << "Assembling the local matrix ..." << endl;
  tic_toc.Clear();
  tic_toc.Start();
  OccaBilinearForm *a = new OccaBilinearForm(ofespace);

  if (use_acrotensor) {
    mfem_error("Not yet implemented");
  } else {
    const IntegrationRule *ir = &(IntRules.Get(ofespace->GetFESpace()->GetFE(0)->GetGeomType(), 2 * order + 3 + quad_add));
    OccaIntegrator *di = new OccaVectorMassIntegrator(1.0, ir);
    a->AddDomainIntegrator(di);

    const int points = di->GetIntegrationRule().GetNPoints();
    std::cout << " using integration rule with " << points << " points" << std::endl;
  }

  a->Assemble();
  tic_toc.Stop();
  const double rt_assemb = tic_toc.RealTime();
  cout << " done, " << rt_assemb << "s." << endl;
  cout << "\n\"DOFs/sec\" in assembly: " << 1e-6*size/rt_assemb << " million.\n" << endl;

  Operator *A;
  OccaVector B, X;
  Array<int> ess_tdof_list;

  cout << "FormLinearSystem() ..." << flush;
  tic_toc.Clear();
  tic_toc.Start();

  a->FormLinearSystem(ess_tdof_list, x, b, A, X, B, 1);

  tic_toc.Stop();
  const double rt_form = tic_toc.RealTime();
  cout << " done, " << rt_form << "s." << endl;
  cout << "\n\"DOFs/sec\" in FormLinearSystem(): " << 1e-6*size/rt_form << " million.\n" << endl;

  cout << "Size of linear system: " << A->Height() << endl;

  cout << "Running CG ...\n" << flush;
  tic_toc.Clear();
  tic_toc.Start();

  CG(*A, B, X, 0, cg_iter, 0.0, 0.0);

  occa::finish();
  tic_toc.Stop();
  const double rt_solve = tic_toc.RealTime();
  cout << " done, " << rt_solve << "s." << endl;
  cout << "\n\"DOFs/sec\" in CG: " << 1e-6*size*cg_iter/rt_solve << " million.\n" << endl;

  // 11. Recover the solution as a finite element grid function.
  a->RecoverFEMSolution(X, b, x);

  // 12. Save the refined mesh and the solution. This output can be viewed later
  //     using GLVis: "glvis -m refined.mesh -g sol.gf".
  ofstream mesh_ofs("refined.mesh");
  mesh_ofs.precision(8);
  mesh->Print(mesh_ofs);
  ofstream sol_ofs("sol.gf");
  // Reuse GridFunction's Save
  GridFunction gf_x(fespace);
  gf_x = x;
  gf_x.Save(sol_ofs);

  // 13. Send the solution by socket to a GLVis server.
  if (visualization)
    {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock.precision(8);
      sol_sock << "solution\n" << *mesh << gf_x << flush;
    }

  // 14. Free the used memory.
  delete a;
  delete fespace;
  if (order > 0) { delete fec; }
  delete mesh;

  return 0;
}
