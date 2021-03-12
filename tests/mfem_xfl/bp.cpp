// Copyright (c) 2010-2021, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.
#include "mfem.hpp"

#include "general/forall.hpp"
#include "linalg/kernels.hpp"

// Kernels addons //////////////////////////////////////////////////////////////
#ifndef MFEM_USE_MPI
#define HYPRE_Int int
typedef int MPI_Session;
#define ParMesh Mesh
#define ParFESpace FESpace
#define GetParMesh GetMesh
#define GlobalTrueVSize GetVSize
#define ParLinearForm LinearForm
#define ParBilinearForm BilinearForm
#define ParGridFunction GridFunction
#define ParFiniteElementSpace FiniteElementSpace
#define PFesGetParMeshGetComm(...)
#define PFesGetParMeshGetComm0(...) 0
#define MPI_Init(...)
#define MPI_Barrier(...)
#define MPI_Finalize()
#define MPI_Comm_size(...)
#define MPI_Comm_rank(...)
#define MPI_Allreduce(src,dst,...) *dst = *src
#define MPI_Reduce(src, dst, n, T,...) *dst = *src
#define NewParMesh(pmesh, mesh, partitioning) pmesh = mesh
#else
#define NewParMesh(pmesh, mesh, partitioning) \
    pmesh = new ParMesh(MPI_COMM_WORLD, *mesh, partitioning);\
    delete mesh;
#endif

using namespace mfem;

namespace mfem
{

mfem::Mesh *CreateMeshEx7(int order);

using FE = mfem::FiniteElement;
using QI = mfem::QuadratureInterpolator;

// Kernels addons //////////////////////////////////////////////////////////////
namespace kernels
{

template<typename T> MFEM_HOST_DEVICE inline
void HouseholderReflect(T *A, const T *v,
                        const T b, const int m, const int n,
                        const int row, const int col)
{
   for (int j = 0; j < n; j++)
   {
      T w = A[0*row + j*col];
      for (int i = 1; i < m; i++) { w += v[i] * A[i*row + j*col]; }
      A[0*row + j*col] -= b * w;
      for (int i = 1; i < m; i++) { A[i*row + j*col] -= b * w * v[i]; }
   }
}

template<int Q1D, typename T> MFEM_HOST_DEVICE inline
void HouseholderApplyQ(T *A, const T *Q, const T *tau,
                       const int k, const int row, const int col)
{
   T v[Q1D];
   for (int ii=0; ii<k; ii++)
   {
      const int i = k-1-ii;
      for (int j = i+1; j < Q1D; j++) { v[j] = Q[j*k+i]; }
      // Apply Householder reflector (I - tau v v^T) coG^T
      HouseholderReflect(&A[i*row], &v[i], tau[i], Q1D-i, Q1D, row, col);
   }
}

template<int D1D, int Q1D, typename T> MFEM_HOST_DEVICE inline
void QRFactorization(T *mat, T *tau)
{
   T v[Q1D];
   DeviceMatrix B(mat, D1D, Q1D);
   for (int i = 0; i < D1D; i++)
   {
      // Calculate Householder vector, magnitude
      T sigma = 0.0;
      v[i] = B(i,i);
      for (int j = i + 1; j < Q1D; j++)
      {
         v[j] = B(i,j);
         sigma += v[j] * v[j];
      }
      T norm = std::sqrt(v[i]*v[i] + sigma); // norm of v[i:m]
      T Rii = -copysign(norm, v[i]);
      v[i] -= Rii;
      // norm of v[i:m] after modification above and scaling below
      //   norm = sqrt(v[i]*v[i] + sigma) / v[i];
      //   tau = 2 / (norm*norm)
      tau[i] = 2 * v[i]*v[i] / (v[i]*v[i] + sigma);
      for (int j=i+1; j<Q1D; j++) { v[j] /= v[i]; }
      // Apply Householder reflector to lower right panel
      HouseholderReflect(&mat[i*D1D+i+1], &v[i], tau[i],
                         Q1D-i, D1D-i-1, D1D, 1);
      // Save v
      B(i,i) = Rii;
      for (int j=i+1; j<Q1D; j++) { B(i,j) = v[j]; }
   }
}

template<int D1D, int Q1D, typename T = double> MFEM_HOST_DEVICE inline
void GetCollocatedGrad(DeviceTensor<2,const T> b,
                       DeviceTensor<2,const T> g,
                       DeviceTensor<2,T> CoG)
{
   T tau[Q1D];
   T B1d[Q1D*D1D];
   T G1d[Q1D*D1D];
   DeviceMatrix B(B1d, D1D, Q1D);
   DeviceMatrix G(G1d, D1D, Q1D);

   for (int d = 0; d < D1D; d++)
   {
      for (int q = 0; q < Q1D; q++)
      {
         B(d,q) = b(q,d);
         G(d,q) = g(q,d);
      }
   }
   QRFactorization<D1D,Q1D>(B1d, tau);
   // Apply Rinv, colograd1d = grad1d Rinv
   for (int i = 0; i < Q1D; i++)
   {
      CoG(0,i) = G(0,i)/B(0,0);
      for (int j = 1; j < D1D; j++)
      {
         CoG(j,i) = G(j,i);
         for (int k = 0; k < j; k++) { CoG(j,i) -= B(j,k)*CoG(k,i); }
         CoG(j,i) /= B(j,j);
      }
      for (int j = D1D; j < Q1D; j++) { CoG(j,i) = 0.0; }
   }
   // Apply Qtranspose, colograd = colograd Qtranspose
   HouseholderApplyQ<Q1D>((T*)CoG, B1d, tau, D1D, 1, Q1D);
}

}  // namespace kernels

// XFL addons //////////////////////////////////////////////////////////////////
namespace xfl
{

struct XElementRestriction : public ElementRestriction
{
   XElementRestriction(const ParFiniteElementSpace *fes,
                       ElementDofOrdering ordering)
      : ElementRestriction(*fes, ordering) { }
   const Array<int> &GatherMap() const { return gatherMap; }
};

/** ****************************************************************************
 * @brief The Operator class
 **************************************************************************** */
template <int DIM> class Operator;

/** ****************************************************************************
 * @brief The 2D Operator class
 **************************************************************************** */
template <>
class Operator<2> : public mfem::Operator
{
protected:
   static constexpr int DIM = 2;
   static constexpr int NBZ = 1;

   const DofToQuad::Mode mode = DofToQuad::TENSOR;
   const int flags = GeometricFactors::JACOBIANS |
                     GeometricFactors::COORDINATES |
                     GeometricFactors::DETERMINANTS;
   const ElementDofOrdering e_ordering = ElementDofOrdering::LEXICOGRAPHIC;

   mfem::ParMesh *mesh;
   const ParFiniteElementSpace *pfes;
   const GridFunction *nodes;
   const mfem::FiniteElementSpace *nfes;
   const int p, q;
   const xfl::XElementRestriction ER;
   const mfem::Operator *NR;
   const Geometry::Type type;
   const IntegrationRule &ir;
   const GeometricFactors *geom;
   const DofToQuad *maps;
   const QuadratureInterpolator *nqi;
   const int SDIM, VDIM, NDOFS, NE, NQ, D1D, Q1D;
   mutable Vector val_xq, grad_xq;
   Vector J0, dx;
   const mfem::Operator *P, *R;

public:
   Operator(const ParFiniteElementSpace *pfes)
      : mfem::Operator(pfes->GetVSize()),
        mesh(pfes->GetParMesh()),
        pfes(pfes),
        nodes((mesh->EnsureNodes(), mesh->GetNodes())),
        nfes(nodes->FESpace()),
        p(pfes->GetFE(0)->GetOrder()),
        q(2 * p + mesh->GetElementTransformation(0)->OrderW()),
        ER(pfes, e_ordering),
        NR(nfes->GetElementRestriction(e_ordering)),
        type(mesh->GetElementBaseGeometry(0)),
        ir(IntRules.Get(type, q)),
        geom(mesh->GetGeometricFactors(ir, flags, mode)),
        maps(&pfes->GetFE(0)->GetDofToQuad(ir, mode)),
        nqi(nfes->GetQuadratureInterpolator(ir, mode)),
        SDIM(mesh->SpaceDimension()),
        VDIM(pfes->GetVDim()),
        NDOFS(pfes->GetNDofs()),
        NE(mesh->GetNE()),
        NQ(ir.GetNPoints()),
        D1D(pfes->GetFE(0)->GetOrder() + 1),
        Q1D(IntRules.Get(Geometry::SEGMENT, ir.GetOrder()).GetNPoints()),
        val_xq(NQ * VDIM * NE),
        grad_xq(NQ * VDIM * DIM * NE),
        J0(SDIM * DIM * NQ * NE),
        dx(NQ * NE, MemoryType::HOST_32),
        P(pfes->GetProlongationMatrix()),
        R(pfes->GetRestrictionMatrix())
   {
      MFEM_VERIFY(DIM == 2, "");
      MFEM_VERIFY(VDIM == 1, "");
      MFEM_VERIFY(SDIM == DIM, "");
      MFEM_VERIFY(NQ == Q1D * Q1D, "");
      MFEM_VERIFY(DIM == mesh->Dimension(), "");
      nqi->SetOutputLayout(QVectorLayout::byVDIM);
      const FiniteElement *fe = nfes->GetFE(0);
      const int vdim = nfes->GetVDim();
      const int nd = fe->GetDof();
      Vector Enodes(vdim * nd * NE);
      NR->Mult(*nodes, Enodes);
      nqi->Derivatives(Enodes, J0);
   }

   virtual void Mult(const mfem::Vector &, mfem::Vector &) const {}

   virtual void QMult(const mfem::Vector &, mfem::Vector &) const {}

   virtual const mfem::Operator *GetProlongation() const { return P; }

   virtual const mfem::Operator *GetRestriction() const { return R; }
};

/** ****************************************************************************
 * @brief The 3D Operator class
 **************************************************************************** */
template <>
class Operator<3> : public mfem::Operator
{
protected:
   static constexpr int DIM = 3;

   const DofToQuad::Mode mode = DofToQuad::TENSOR;
   const int flags = GeometricFactors::JACOBIANS |
                     GeometricFactors::COORDINATES |
                     GeometricFactors::DETERMINANTS;
   const ElementDofOrdering e_ordering = ElementDofOrdering::LEXICOGRAPHIC;

   mfem::ParMesh *mesh;
   const ParFiniteElementSpace *pfes;
   const GridFunction *nodes;
   const mfem::FiniteElementSpace *nfes;
   const int p, q;
   const xfl::XElementRestriction ER;
   const mfem::Operator *NR;
   const Geometry::Type type;
   const IntegrationRule &ir;
   const GeometricFactors *geom;
   const DofToQuad *maps;
   const QuadratureInterpolator *qi, *nqi;
   const int SDIM, VDIM, NDOFS, NE, NQ, D1D, Q1D;
   mutable Vector val_xq, grad_xq;
   Vector J0, dx;
   const mfem::Operator *P, *R;
   mutable Array<double> CoG;

public:
   Operator(const ParFiniteElementSpace *pfes)
      : mfem::Operator(pfes->GetVSize()),
        mesh(pfes->GetParMesh()),
        pfes(pfes),
        nodes((mesh->EnsureNodes(), mesh->GetNodes())),
        nfes(nodes->FESpace()),
        p(pfes->GetFE(0)->GetOrder()),
        q(2 * p + mesh->GetElementTransformation(0)->OrderW()),
        ER(pfes, e_ordering),
        NR(nfes->GetElementRestriction(e_ordering)),
        type(mesh->GetElementBaseGeometry(0)),
        ir(IntRules.Get(type, q)),
        geom(mesh->GetGeometricFactors(ir, flags, mode)),
        maps(&pfes->GetFE(0)->GetDofToQuad(ir, mode)),
        qi(pfes->GetQuadratureInterpolator(ir, mode)),
        nqi(nfes->GetQuadratureInterpolator(ir, mode)),
        SDIM(mesh->SpaceDimension()),
        VDIM(pfes->GetVDim()),
        NDOFS(pfes->GetNDofs()),
        NE(mesh->GetNE()),
        NQ(ir.GetNPoints()),
        D1D(pfes->GetFE(0)->GetOrder() + 1),
        Q1D(IntRules.Get(Geometry::SEGMENT, ir.GetOrder()).GetNPoints()),
        val_xq(NQ * VDIM * NE),
        grad_xq(NQ * VDIM * DIM * NE),
        J0(SDIM * DIM * NQ * NE),
        dx(NQ * NE),
        P(pfes->GetProlongationMatrix()),
        R(pfes->GetRestrictionMatrix())
   {
      MFEM_VERIFY(DIM == 3, "");
      MFEM_VERIFY(VDIM == 1, "");
      MFEM_VERIFY(SDIM == DIM, "");
      MFEM_VERIFY(NQ == Q1D * Q1D * Q1D, "");
      MFEM_VERIFY(DIM == mesh->Dimension(), "");
      qi->SetOutputLayout(QVectorLayout::byVDIM);
      nqi->SetOutputLayout(QVectorLayout::byVDIM);
      const FiniteElement *fe = nfes->GetFE(0);
      const int vdim = nfes->GetVDim();
      const int nd = fe->GetDof();
      Vector Enodes(vdim * nd * NE);
      NR->Mult(*nodes, Enodes);
      nqi->Derivatives(Enodes, J0);
   }

   virtual void Mult(const mfem::Vector &, mfem::Vector &) const {}

   virtual void QMult(const mfem::Vector &, mfem::Vector &) const {}

   virtual const mfem::Operator *GetProlongation() const { return P; }

   virtual const mfem::Operator *GetRestriction() const { return R; }
};

/** ****************************************************************************
 * @brief The Problem struct
 ******************************************************************************/
struct Problem
{
   mfem::Operator *QM {nullptr};
   mfem::ParLinearForm &b;
   Problem(mfem::ParLinearForm &b, mfem::Operator *QM) : QM(QM), b(b) {}
   ~Problem() { delete QM; }
};

/** ****************************************************************************
 * @brief The QForm class
 ******************************************************************************/
class QForm
{
public:
   const char *qs;
   mfem::Operator *QM;
   mfem::ParLinearForm *b = nullptr;
   mfem::ConstantCoefficient *constant_coeff = nullptr;
   mfem::FunctionCoefficient *function_coeff = nullptr;
   mfem::ParFiniteElementSpace *pfes;

public:

   QForm(ParFiniteElementSpace *pfes, const char *qs, mfem::Operator *QM)
      : qs(qs), QM(QM), pfes(pfes) { }

   ~QForm() { }

   // Create problem
   Problem *operator==(QForm &rhs)
   {
      assert(!b);
      mfem::ParLinearForm *b = new mfem::ParLinearForm(rhs.ParFESpace());
      assert(b);
      if (!rhs.ConstantCoeff() && !rhs.FunctionCoeff())
      {
         ConstantCoefficient *cst = new ConstantCoefficient(1.0);
         b->AddDomainIntegrator(new DomainLFIntegrator(*cst));
      }
      else if (rhs.ConstantCoeff())
      {
         ConstantCoefficient *cst = rhs.ConstantCoeff();
         b->AddDomainIntegrator(new DomainLFIntegrator(*cst));
      }
      else if (rhs.FunctionCoeff())
      {
         FunctionCoefficient *func = rhs.FunctionCoeff();
         b->AddDomainIntegrator(new DomainLFIntegrator(*func));
      }
      else
      {
         assert(false);
      }

      return new Problem(*b, QM);
   }

   // + operator on QForms
   QForm &operator+(QForm &rhs)
   {
      assert(false);  // not supported
      return *this;
   }

   mfem::ParFiniteElementSpace *ParFESpace() { return pfes; }
   mfem::ConstantCoefficient *ConstantCoeff() const { return constant_coeff; }
   mfem::FunctionCoefficient *FunctionCoeff() const { return function_coeff; }
};

/** ****************************************************************************
 * @brief Function class
 ******************************************************************************/
class Function : public ParGridFunction
{
public:
   Function(ParFiniteElementSpace *pfes) : ParGridFunction(pfes)
   {
      assert(pfes);
      assert(pfes->GlobalTrueVSize() > 0);
   }
   void operator=(double value) { ParGridFunction::operator=(value); }
   int geometric_dimension() { return fes->GetMesh()->SpaceDimension(); }
   ParFiniteElementSpace *ParFESpace() { return ParGridFunction::ParFESpace(); }
   const ParFiniteElementSpace *ParFESpace() const { return ParGridFunction::ParFESpace(); }
   ConstantCoefficient *ConstantCoeff() const { return nullptr; }
   FunctionCoefficient *FunctionCoeff() const { return nullptr; }
};

/** ****************************************************************************
 * @brief TrialFunction class
 ******************************************************************************/
class TrialFunction : public Function
{
public:
   TrialFunction(ParFiniteElementSpace *pfes) : Function(pfes) { }
   ~TrialFunction() { }
};

/** ****************************************************************************
 * @brief TestFunction class
 ******************************************************************************/
class TestFunction : public Function
{
public:
   TestFunction(ParFiniteElementSpace *pfes) : Function(pfes) { }
   TestFunction(const TestFunction &) = default;
   ~TestFunction() { }
};

/** ****************************************************************************
 * @brief Constant class
 ******************************************************************************/
class Constant
{
   const double value = 0.0;
   ConstantCoefficient *cst = nullptr;

public:
   Constant(double val) : value(val), cst(new ConstantCoefficient(val)) { }
   ~Constant() { /*delete cst;*/ }
   ParFiniteElementSpace *ParFESpace() const { return nullptr; }
   double Value() const { return value; }
   double Value() { return value; }
   double operator*(TestFunction &v) { return 0.0; }
   ConstantCoefficient *ConstantCoeff() const { return cst; }
   FunctionCoefficient *FunctionCoeff() const { return nullptr; }
   operator const double *() const { return nullptr; }  // qf eval
};

/** ****************************************************************************
 * @brief Expressions
 ******************************************************************************/
class Expression
{
   FunctionCoefficient *fct = nullptr;

public:
   Expression(std::function<double(const Vector &)> F)
      : fct(new FunctionCoefficient(F)) { }
   ~Expression() { /*delete fct;*/ }
   ParFiniteElementSpace *ParFESpace() const { return nullptr; }  // qf args
   ConstantCoefficient *ConstantCoeff() const { return nullptr; }
   FunctionCoefficient *FunctionCoeff() const { return fct; }
   operator const double *() const { return nullptr; }  // qf eval
   // double operator *(TestFunction &v) { return 0.0;} // qf eval
};

/** ****************************************************************************
 * @brief Mesh
 ******************************************************************************/
static mfem::ParMesh *MeshToPMesh(mfem::Mesh *mesh)
{
   int num_procs = 1;
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   int n111[3] {1, 1, 1};
   int n211[3] {2, 1, 1};
   int n221[3] {2, 2, 1};
   int n222[3] {2, 2, 2};
   int n422[3] {4, 2, 2};
   int n442[3] {4, 4, 2};
   int n444[3] {4, 4, 4};
   int *nxyz = (num_procs == 1 ? n111 :
                num_procs == 2 ? n211 :
                num_procs == 4 ? n221 :
                num_procs == 8 ? n222 :
                num_procs == 16 ? n422 :
                num_procs == 32 ? n442 :
                num_procs == 64 ? n444 : nullptr);
   assert(nxyz);
   const int mesh_p = 1;
   mesh->SetCurvature(mesh_p, false, -1, Ordering::byNODES);
   int *partitioning = mesh->CartesianPartitioning(nxyz);
   ParMesh *pmesh = nullptr;
   NewParMesh(pmesh, mesh, partitioning);
   return pmesh;
}

mfem::ParMesh &Mesh(const char *mesh_file)  // and & for mesh
{
   return *MeshToPMesh(new mfem::Mesh(mesh_file, 1, 1));
}

#ifdef MFEM_USE_MPI
mfem::ParMesh &Mesh(mfem::Mesh *mesh) { return *MeshToPMesh(mesh); }

mfem::ParMesh &Mesh(mfem::ParMesh *pmesh) { return *pmesh; }
#else
mfem::Mesh &Mesh(mfem::Mesh *mesh) { return *mesh; }
#endif

mfem::ParMesh &UnitSquareMesh(int nx, int ny)
{
   const double sx = 1.0, sy = 1.0;
   Element::Type quad = Element::Type::QUADRILATERAL;
   const bool edges = false, sfc = true;
   mfem::Mesh *mesh =
      new mfem::Mesh(nx, ny, quad, edges, sx, sy, sfc);
   return *MeshToPMesh(mesh);
}

mfem::ParMesh &UnitHexMesh(int nx, int ny, int nz)
{
   Element::Type hex = Element::Type::HEXAHEDRON;
   const bool edges = false, sfc = true;
   const double sx = 1.0, sy = 1.0, sz = 1.0;
   mfem::Mesh *mesh =
      new mfem::Mesh(nx, ny, nz, hex, edges, sx, sy, sz, sfc);
   return *MeshToPMesh(mesh);
}

/** ****************************************************************************
 * @brief Device
 ******************************************************************************/
mfem::Device Device(const char *device_config) { return {device_config}; }

/** ****************************************************************************
 * @brief FiniteElement
 ******************************************************************************/
FiniteElementCollection *FiniteElement(std::string family, int type, int p)
{
   MFEM_VERIFY(family == "Lagrange", "Unsupported family!");
   MFEM_VERIFY(type == Element::Type::QUADRILATERAL ||
               type == Element::Type::HEXAHEDRON, "Unsupported type!");
   const int dim = (type == Element::Type::QUADRILATERAL) ? 2 :
                   (type == Element::Type::HEXAHEDRON)  ? 3 : 0;
   const int btype = BasisType::GaussLobatto;
   return new H1_FECollection(p, dim, btype);
}

/** ****************************************************************************
 * @brief Function Spaces
 ******************************************************************************/
class FunctionSpace : public ParFiniteElementSpace {};

ParFiniteElementSpace *FunctionSpace(mfem::ParMesh *pmesh,
                                     std::string family,
                                     int p)
{
   assert(false);
   const int dim = pmesh->Dimension();
   MFEM_VERIFY(family == "P", "Unsupported FE!");
   FiniteElementCollection *fec = new H1_FECollection(p, dim);
   return new ParFiniteElementSpace(pmesh, fec);
}

ParFiniteElementSpace *FunctionSpace(mfem::ParMesh &pmesh, std::string f, int p)
{
   assert(false);
   return FunctionSpace(&pmesh, f, p);
}

ParFiniteElementSpace *FunctionSpace(mfem::ParMesh &pmesh,
                                     FiniteElementCollection *fec)
{
   const int vdim = 1;
   const Ordering::Type ordering = Ordering::byNODES;
   ParFiniteElementSpace *pfes =
      new ParFiniteElementSpace(&pmesh, fec, vdim, ordering);
   return pfes;
}

ParFiniteElementSpace *FunctionSpace(mfem::ParMesh &pmesh,
                                     FiniteElementCollection *fec,
                                     const int vdim)
{
   assert(false);
   return new ParFiniteElementSpace(&pmesh, fec, vdim);
}

/** ****************************************************************************
 * @brief Vector Function Space
 ******************************************************************************/
ParFiniteElementSpace *VectorFunctionSpace(mfem::ParMesh *pmesh,
                                           std::string family,
                                           const int p)
{
   const int dim = pmesh->Dimension();
   MFEM_VERIFY(family == "P", "Unsupported FE!");
   FiniteElementCollection *fec = new H1_FECollection(p, dim);
   return new ParFiniteElementSpace(pmesh, fec, dim);
}

ParFiniteElementSpace *VectorFunctionSpace(mfem::ParMesh &pmesh,
                                           std::string family,
                                           const int p)
{
   return VectorFunctionSpace(&pmesh, family, p);
}

ParFiniteElementSpace *VectorFunctionSpace(mfem::ParMesh &pmesh,
                                           FiniteElementCollection *fec)
{
   return new ParFiniteElementSpace(&pmesh, fec, pmesh.Dimension());
}

/** ****************************************************************************
 * @brief Boundary Conditions
 ******************************************************************************/
Array<int> DirichletBC(mfem::ParFiniteElementSpace *pfes)
{
   assert(pfes);
   Array<int> ess_tdof_list;
   mfem::ParMesh *pmesh = pfes->GetParMesh();
   if (pmesh->bdr_attributes.Size())
   {
      Array<int> ess_bdr(pmesh->bdr_attributes.Max());
      ess_bdr = 1;
      pfes->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }
   return ess_tdof_list;
}

/** ****************************************************************************
 * @brief Math namespace
 ******************************************************************************/
namespace math
{

Constant Pow(Function &gf, double exp)
{
   return Constant(gf.Vector::Normlp(exp));
}

}  // namespace math

/** ****************************************************************************
 * @brief solve with boundary conditions
 ******************************************************************************/
int solve(xfl::Problem *pb, xfl::Function &x, Array<int> ess_tdof_list)
{
   assert(x.FESpace());
   ParFiniteElementSpace *fes = x.ParFESpace();
   MFEM_VERIFY(UsesTensorBasis(*fes), "FE Space must Use Tensor Basis!");

   Vector B, X;
   pb->b.Assemble();
   mfem::Operator *A = nullptr;

   mfem::Operator &op = *(pb->QM);
   op.FormLinearSystem(ess_tdof_list, x, pb->b, A, X, B);
   CG(*A, B, X, 1, 400, 1e-12, 0.0);
   op.RecoverFEMSolution(X, pb->b, x);
   x.HostReadWrite();

   delete pb;
   return 0;
}

/// solve with empty boundary conditions
int solve(xfl::Problem *pb, xfl::Function &x)
{
   Array<int> empty_tdof_list;
   return solve(pb, x, empty_tdof_list);
}

/** ****************************************************************************
 * @brief benchmark this prblem with boundary conditions
 ******************************************************************************/
int benchmark(xfl::Problem *pb, xfl::Function &x, Array<int> ess_tdof_list,
              const double rtol, const int max_it, const int print_lvl)
{
   assert(x.ParFESpace());
   ParFiniteElementSpace *pfes = x.ParFESpace();
   assert(pfes->GlobalTrueVSize() > 0);
   MFEM_VERIFY(UsesTensorBasis(*pfes), "FE Space must Use Tensor Basis!");

   mfem::ParLinearForm &b = pb->b;
   b.Assemble();

   Vector B, X;
   mfem::Operator *A;
   mfem::Operator &a = *(pb->QM);
   a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

#ifndef MFEM_USE_MPI
   CGSolver cg;
#else
   CGSolver cg(MPI_COMM_WORLD);
#endif

   cg.SetRelTol(rtol);
   cg.SetOperator(*A);

   // Warm-up CG solve (in case of JIT to avoid timing it)
   {
      Vector Y(X);
      cg.SetMaxIter(2);
      cg.SetPrintLevel(-1);
      cg.Mult(B, Y);
   }

   // benchmark this problem
   {
      tic_toc.Clear();
      cg.SetMaxIter(max_it);
      cg.SetPrintLevel(print_lvl);
      {
         tic_toc.Start();
         cg.Mult(B, X);
         tic_toc.Stop();
      }
   }

   // MFEM_VERIFY(cg.GetConverged(), "CG did not converged!");
   MFEM_VERIFY(cg.GetNumIterations() <= max_it, "");
   a.RecoverFEMSolution(X, b, x);
   x.HostReadWrite();

   int myid = 0;
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);

   const double rt = tic_toc.RealTime();
   double rt_min, rt_max;
   MPI_Reduce(&rt, &rt_min, 1, MPI_DOUBLE, MPI_MIN, 0,
              pfes->GetParMesh()->GetComm());
   MPI_Reduce(&rt, &rt_max, 1, MPI_DOUBLE, MPI_MAX, 0,
              pfes->GetParMesh()->GetComm());

   const HYPRE_Int dofs = pfes->GlobalTrueVSize();
   const int cg_iter = cg.GetNumIterations();
   const double mdofs_max = ((1e-6 * dofs) * cg_iter) / rt_max;
   const double mdofs_min = ((1e-6 * dofs) * cg_iter) / rt_min;

   if (myid == 0)
   {
      std::cout << "Number of finite element unknowns: " << dofs <<  std::endl;
      std::cout << "Total CG time:    " << rt_max << " (" << rt_min << ") sec."
                << std::endl;
      std::cout << "Time per CG step: "
                << rt_max / cg_iter << " ("
                << rt_min / cg_iter << ") sec." << std::endl;
      std::cout << "\033[32m";
      std::cout << "\"DOFs/sec\" in CG: " << mdofs_max << " ("
                << mdofs_min << ") million.";
      std::cout << "\033[m" << std::endl;
   }
   delete pb;
   return 0;
}

int benchmark(xfl::Problem *pb, xfl::Function &x, Array<int> ess_tdof_list)
{
   return benchmark(pb, x, ess_tdof_list, 1e-12, 200, -1);
}

/** ****************************************************************************
 * @brief plot the x gridfunction
 ******************************************************************************/
int plot(xfl::Function &x)
{
   int num_procs = 1, myid = 0;
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);

   ParFiniteElementSpace *fes = x.ParFESpace();
   assert(fes);
   mfem::ParMesh *pmesh = fes->GetParMesh();
   assert(pmesh);
   char vishost[] = "localhost";
   int visport = 19916;
   socketstream sol_sock(vishost, visport);
   sol_sock << "parallel " << num_procs << " " << myid << "\n";
   sol_sock.precision(8);
   sol_sock << "solution\n" << *pmesh << x << std::flush;
   return 0;
}

/** ****************************************************************************
 * @brief plot the mesh
 ******************************************************************************/
int plot(mfem::ParMesh *pmesh)
{
   int num_procs = 1, myid = 0;
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);

   char vishost[] = "localhost";
   int visport = 19916;
   socketstream sol_sock(vishost, visport);
   sol_sock << "parallel " << num_procs << " " << myid << "\n";
   sol_sock.precision(8);
   sol_sock << "mesh\n" << *pmesh << std::flush;
   return 0;
}

/** ****************************************************************************
 * @brief save the x gridfunction
 ******************************************************************************/
int save(xfl::Function &x, const char *filename)
{
   int myid = 0;
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   std::ostringstream sol_name;
   sol_name << filename << "." << std::setfill('0') << std::setw(6) << myid;

   std::ofstream sol_ofs(sol_name.str().c_str());
   sol_ofs.precision(8);
   x.Save(sol_ofs);
   return 0;
}

/** ****************************************************************************
 * @brief save the x gridfunction
 ******************************************************************************/
int save(mfem::ParMesh &mesh, const char *filename)
{
   int myid = 0;
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   std::ostringstream mesh_name;
   mesh_name << filename << "." << std::setfill('0') << std::setw(6) << myid;

   std::ofstream mesh_ofs(mesh_name.str().c_str());
   mesh_ofs.precision(8);
   mesh.Print(mesh_ofs);
   return 0;
}

//constexpr int point = Element::Type::POINT;
//constexpr int segment = Element::Type::SEGMENT;
//constexpr int triangle = Element::Type::TRIANGLE;
constexpr int quadrilateral = Element::Type::QUADRILATERAL;
//constexpr int tetrahedron = Element::Type::TETRAHEDRON;
constexpr int hexahedron = Element::Type::HEXAHEDRON;
//constexpr int wedge = Element::Type::WEDGE;

}  // namespace xfl

template <typename... Args>
void print(const char *fmt, Args... args)
{
   std::cout << std::flush;
   std::printf(fmt, args...);
   std::cout << std::endl;
}

inline bool UsesTensorBasis(const FiniteElementSpace *fes)
{
   return mfem::UsesTensorBasis(*fes);
}

int sym(int u) { return u; }
int dot(int u, int v) { return u * v; }

}  // namespace mfem


// SIMD - Headers and Sizes
#include "linalg/simd.hpp"
using Real = AutoSIMDTraits<double,double>::vreal_t;
#define SIMD_SIZE (MFEM_SIMD_BYTES/sizeof(double))
#define SMS SIMD_SIZE

// Pragma - UNROLL & ALIGN
#define MFEM_ALIGN(T) alignas(alignof(T)) T
#define PRAGMA(X) _Pragma(#X)
#ifdef __clang__
#define UNROLL(N) PRAGMA(unroll(N))
#elif defined(__FUJITSU)
#define UNROLL(N) PRAGMA(loop unroll N)
#elif defined(__GNUC__) || defined(__GNUG__)
#define UNROLL(N) PRAGMA(GCC unroll(N))
#else
#define UNROLL(N)
#endif

// Blocking
#define BLK_SZ 1
#define FOREACH_BLOCK(i,N) for(int i##i=0; i##i<N+BLK_SZ-1; i##i+=BLK_SZ)
#define FOREACH_INNER(i) for(int i=i##i; i<i##i+BLK_SZ; i++)

// Foreach Thread
#define FOREACH_THREAD(i,k,N) for(int i=0; i<N; i+=1)
#define SYNC_THREADS


template<int DIM, int DX0, int DX1, int Q1D> inline static
void KSetup1(const int ndofs, const int vdim, const int NE,
             const double * __restrict__ J0,
             const double * __restrict__ w,
             double * __restrict__ dx)
{
   assert(vdim == 1);
   const auto J = Reshape(J0, DIM,DIM, Q1D,Q1D,Q1D, NE);
   const auto W = Reshape(w, Q1D,Q1D,Q1D);
   auto DX = Reshape(dx, Q1D,Q1D,Q1D, 6);
   const int e = 0;
   {
      MFEM_FOREACH_THREAD(qz,z,Q1D)
      {
         MFEM_FOREACH_THREAD(qy,y,Q1D)
         {
            MFEM_FOREACH_THREAD(qx,x,Q1D)
            {
               const double irw = W(qx,qy,qz);
               const double *Jtr = &J(0,0,qx,qy,qz, e);
               const double detJ = kernels::Det<DIM>(Jtr);
               const double wd = irw * detJ;
               double Jrt[DIM*DIM];
               kernels::CalcInverse<DIM>(Jtr, Jrt);
               double A[DX0*DX1];
               const double D[DX0*DX1] = {wd,0,0,0,wd,0,0,0,wd};
               kernels::MultABt(DIM,DIM,DIM,D,Jrt,A);
               double B[DX0*DX1];
               kernels::Mult(DIM,DIM,DIM,A,Jrt,B);
               DX(qx,qy,qz,0) = B[0+DX0*0];
               DX(qx,qy,qz,1) = B[1+DX0*0];
               DX(qx,qy,qz,2) = B[2+DX0*0];
               DX(qx,qy,qz,3) = B[1+DX0*1];
               DX(qx,qy,qz,4) = B[2+DX0*1];
               DX(qx,qy,qz,5) = B[2+DX0*2];
            }
         }
      }
      MFEM_SYNC_THREAD;
   }
}

template<int DIM, int DX0, int DX1, int D1D, int Q1D> inline static
void KMult1(const int ndofs, const int vdim, const int NE,
            const double * __restrict__ B,
            const double * __restrict__ G,
            const int * __restrict__ map,
            const double * __restrict__ dx,
            const double * __restrict__ xd,
            double * __restrict__ yd)
{
   // kernel operations: u,G,*,D,*,v,G,T,*,Ye
   const auto b = Reshape(B, Q1D,D1D);
   const auto g = Reshape(G, Q1D,Q1D);
   const auto DX = Reshape(dx, Q1D,Q1D,Q1D, 6);
   const auto MAP = Reshape(map, D1D,D1D,D1D, NE);
   const auto XD = Reshape(xd, ndofs);
   auto YD = Reshape(yd, ndofs);
   MFEM_ALIGN(Real) s_Iq[Q1D][Q1D][Q1D];
   MFEM_ALIGN(double) s_D[Q1D][Q1D];
   MFEM_ALIGN(double) s_I[Q1D][D1D];
   MFEM_VERIFY((NE % SIMD_SIZE) == 0, "NE vs SIMD_SIZE error!");
   for (int e = 0; e < NE; e+=SMS)
   {
      MFEM_ALIGN(Real) r_qt[Q1D][Q1D];
      MFEM_ALIGN(Real) r_q[Q1D][Q1D][Q1D];
      MFEM_ALIGN(Real) r_Aq[Q1D][Q1D][Q1D];

      // Scatter X
      FOREACH_BLOCK(j,Q1D)
      {
         FOREACH_INNER(j)
         {
            if (j>=Q1D) { continue; }
            FOREACH_THREAD(i,x,Q1D)
            {
               s_D[j][i] = g(i,j);
               if (i<D1D) { s_I[j][i] = b(j,i); }
               if (i<D1D && j<D1D)
               {
                  UNROLL(6)
                  for (int k = 0; k < D1D; k++)
                  {
                     MFEM_ALIGN(Real) vXD;
                     UNROLL(SMS)
                     for (size_t v = 0; v < SMS; v++)
                     {
                        const int gid = MAP(i, j, k, e + v);
                        const int idx = gid >= 0 ? gid : -1 - gid;
                        vXD[v] = XD(idx);
                     }
                     r_q[j][i][k] = vXD;
                  }
               }
            }
         }
      } SYNC_THREADS;

      // Grad1X
      FOREACH_THREAD(b,y,D1D)
      {
         FOREACH_THREAD(a,x,D1D)
         {
            UNROLL(7)
            for (int k=0; k<Q1D; ++k)
            {
               MFEM_ALIGN(Real) res; res = 0.0;
               UNROLL(6)
               for (int c=0; c<D1D; ++c)
               {
                  res.fma(s_I[k][c], r_q[b][a][c]);
               }
               s_Iq[k][b][a] = res;
            }
         }
      } SYNC_THREADS;

      // Grad1Y
      FOREACH_THREAD(k,y,Q1D)
      {
         FOREACH_THREAD(a,x,D1D)
         {
            for (int b=0; b<D1D; ++b)
            {
               r_Aq[k][a][b] = s_Iq[k][b][a];
            }
            UNROLL(7)
            for (int j=0; j<Q1D; ++j)
            {
               MFEM_ALIGN(Real) res; res = 0;
               UNROLL(6)
               for (int b=0; b<D1D; ++b)
               {
                  res.fma(s_I[j][b], r_Aq[k][a][b]);
               }
               s_Iq[k][j][a] = res;
            }
         }
      } SYNC_THREADS;

      // Grad1Z
      FOREACH_THREAD(k,y,Q1D)
      {
         FOREACH_THREAD(j,x,Q1D)
         {
            for (int a=0; a<D1D; ++a)
            {
               r_Aq[k][j][a] = s_Iq[k][j][a];
            }
            UNROLL(7)
            for (int i=0; i<Q1D; ++i)
            {
               MFEM_ALIGN(Real) res; res = 0;
               UNROLL(6)
               for (int a=0; a<D1D; ++a)
               {
                  res.fma(s_I[i][a], r_Aq[k][j][a]);
               }
               s_Iq[k][j][i] = res;
            }
         }
      } SYNC_THREADS;

      // Flush
      FOREACH_THREAD(j,y,Q1D)
      {
         FOREACH_THREAD(i,x,Q1D)
         {
            UNROLL(7)
            for (int k = 0; k < Q1D; k++) { r_Aq[j][i][k] = 0.0; }
         }
      } SYNC_THREADS;

      // Q-Function
      UNROLL(7)
      for (int k = 0; k < Q1D; k++)
      {
         SYNC_THREADS;
         MFEM_ALIGN(Real) r_Gqr[Q1D][Q1D];
         MFEM_ALIGN(Real) r_Gqs[Q1D][Q1D];
         FOREACH_THREAD(j,y,Q1D)
         {
            FOREACH_BLOCK(i,Q1D)
            {
               FOREACH_INNER(i)
               {
                  if (i>=Q1D) { continue; }
                  MFEM_ALIGN(Real) qr, qs, qt; qr = 0.0; qs = 0.0; qt = 0.0;
                  UNROLL(7)
                  for (int m = 0; m < Q1D; m++)
                  {
                     const double Dim = s_D[i][m];
                     const double Djm = s_D[j][m];
                     const double Dkm = s_D[k][m];
                     qr.fma(Dim, s_Iq[k][j][m]);
                     qs.fma(Djm, s_Iq[k][m][i]);
                     qt.fma(Dkm, s_Iq[m][j][i]);
                  }
                  const double D00 = DX(i,j,k,0);
                  const double D01 = DX(i,j,k,1);
                  const double D02 = DX(i,j,k,2);
                  const double D10 = D01;
                  const double D11 = DX(i,j,k,3);
                  const double D12 = DX(i,j,k,4);
                  const double D20 = D02;
                  const double D21 = D12;
                  const double D22 = DX(i,j,k,5);
                  r_Gqr[j][i] = 0.0;
                  r_Gqr[j][i].fma(D00,qr);
                  r_Gqr[j][i].fma(D10,qs);
                  r_Gqr[j][i].fma(D20,qt);
                  r_Gqs[j][i] = 0.0;
                  r_Gqs[j][i].fma(D01,qr);
                  r_Gqs[j][i].fma(D11,qs);
                  r_Gqs[j][i].fma(D21,qt);
                  r_qt[j][i] = 0.0;
                  r_qt[j][i].fma(D02,qr);
                  r_qt[j][i].fma(D12,qs);
                  r_qt[j][i].fma(D22,qt);
               }
            }
         } // BLOCK
         SYNC_THREADS;
         FOREACH_THREAD(j,y,Q1D)
         {
            FOREACH_THREAD(i,x,Q1D)
            {
               MFEM_ALIGN(Real) Aqtmp; Aqtmp = 0.0;
               UNROLL(7)
               for (int m = 0; m < Q1D; m++)
               {
                  const double Dmi = s_D[m][i];
                  const double Dmj = s_D[m][j];
                  const double Dkm = s_D[k][m];
                  Aqtmp.fma(Dmi, r_Gqr[j][m]);
                  Aqtmp.fma(Dmj, r_Gqs[m][i]);
                  r_Aq[j][i][m].fma(Dkm, r_qt[j][i]);
               }
               r_Aq[j][i][k] += Aqtmp;
            }
         } SYNC_THREADS;
      }
      // GradZT
      FOREACH_THREAD(j,y,Q1D)
      {
         FOREACH_THREAD(i,x,Q1D)
         {
            UNROLL(6)
            for (int c=0; c<D1D; ++c)
            {
               MFEM_ALIGN(Real) res; res = 0;
               UNROLL(7)
               for (int k=0; k<Q1D; ++k)
               {
                  res.fma(s_I[k][c], r_Aq[j][i][k]);
               }
               s_Iq[c][j][i] = res;
            }
         }
      } SYNC_THREADS;
      // GradYT
      FOREACH_THREAD(c,y,D1D)
      {
         FOREACH_THREAD(i,x,Q1D)
         {
            UNROLL(7)
            for (int j=0; j<Q1D; ++j)
            {
               r_Aq[c][i][j] = s_Iq[c][j][i];
            }
            UNROLL(6)
            for (int b=0; b<D1D; ++b)
            {
               MFEM_ALIGN(Real) res; res = 0;
               UNROLL(7)
               for (int j=0; j<Q1D; ++j)
               {
                  res.fma(s_I[j][b], r_Aq[c][i][j]);
               }
               s_Iq[c][b][i] = res;
            }
         }
      } SYNC_THREADS;
      // GradXT
      FOREACH_THREAD(c,y,D1D)
      {
         FOREACH_THREAD(b,x,D1D)
         {
            UNROLL(7)
            for (int i=0; i<Q1D; ++i)
            {
               r_Aq[c][b][i] = s_Iq[c][b][i];
            }
            UNROLL(6)
            for (int a=0; a<D1D; ++a)
            {
               MFEM_ALIGN(Real) res; res = 0;
               UNROLL(7)
               for (int i=0; i<Q1D; ++i)
               {
                  res.fma(s_I[i][a], r_Aq[c][b][i]);
               }
               s_Iq[c][b][a] = res;
               s_Iq[c][b][a] = res;
            }
         }
      } SYNC_THREADS;
      // Gather
      FOREACH_THREAD(j,y,D1D)
      {
         FOREACH_THREAD(i,x,D1D)
         {
            UNROLL(6)
            for (int k = 0; k < D1D; k++)
            {
               UNROLL(SMS)
               for (size_t v = 0; v < SMS; v++)
               {
                  const int gid = MAP(i, j, k, e + v);
                  const int idx = gid >= 0 ? gid : -1 - gid;
                  YD(idx) += (s_Iq[k][j][i])[v];
               }
            }
         }
      }
   } // MFEM_FORALL_2D
} // KMult1

#ifndef GEOM
#define GEOM Geometry::CUBE
#endif

#ifndef MESH_P
#define MESH_P 1
#endif

#ifndef SOL_P
#define SOL_P 5
#endif

#ifndef IR_ORDER
#define IR_ORDER 0
#endif

#ifndef IR_TYPE
#define IR_TYPE 0
#endif

#ifndef PROBLEM
#define PROBLEM 0
#endif

#ifndef VDIM
#define VDIM 1
#endif

int main(int argc, char* argv[])
{
   int status = 0;
   int num_procs = 1, myid = 0;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);

   assert(VDIM==1); assert(MESH_P==1); assert(IR_TYPE==0); assert(IR_ORDER==0);
   assert(PROBLEM==0); assert(GEOM==Geometry::CUBE);
   const char *mesh_file = "../../data/hex-01x01x01.mesh"; int ser_ref_levels = 4;
   int par_ref_levels = 0; Array<int> nxyz; int order = SOL_P;
   const char *basis_type = "G"; bool static_cond = false; const char *pc = "none";
   bool perf = true; bool matrix_free = true; int max_iter = 50;
   bool visualization = false; OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
   args.AddOption(&ser_ref_levels, "-rs", "--refine-serial",
                  "Number of times to refine the mesh uniformly in serial.");
   args.AddOption(&par_ref_levels, "-rp", "--refine-parallel",
                  "Number of times to refine the mesh uniformly in parallel.");
   args.AddOption(&nxyz, "-c", "--cartesian-partitioning",
                  "Use Cartesian partitioning.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree) or -1 for" " isoparametric space.");
   args.AddOption(&basis_type, "-b", "--basis-type",
                  "Basis: G - Gauss-Lobatto, P - Positive, U - Uniform");
   args.AddOption(&perf, "-perf", "--hpc-version", "-std", "--standard-version",
                  "Enable high-performance, tensor-based, assembly/evaluation.");
   args.AddOption(&matrix_free, "-mf", "--matrix-free", "-asm", "--assembly",
                  "Use matrix-free evaluation or efficient matrix assembly in "
                  "the high-performance version.");
   args.AddOption(&pc, "-pc", "--preconditioner",
                  "Preconditioner: lor - low-order-refined (matrix-free) AMG, "
                  "ho - high-order (assembled) AMG, none.");
   args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
                  "--no-static-condensation", "Enable static condensation.");
   args.AddOption(&max_iter, "-mi", "--max-iter", "Maximum number of iterations.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization", "Enable or disable GLVis visualization."); args.Parse();
   if (!args.Good()) { if (myid == 0) { args.PrintUsage(std::cout); } return 1; } if (
      myid == 0) { args.PrintOptions(std::cout); } assert(SOL_P == order);
   ParMesh *pmesh = nullptr; { Mesh *mesh = new Mesh(mesh_file, 1, 1); int dim = mesh->Dimension(); { int ref_levels = (int)floor(log(10000./mesh->GetNE())/log(2.)/dim); ref_levels = (ser_ref_levels != -1) ? ser_ref_levels : ref_levels; for (int l = 0; l < ref_levels; l++) { if (myid == 0) { std::cout << "Serial refinement: level " << l << " -> level " << l+1 << " ..." << std::flush; } mesh->UniformRefinement(); MPI_Barrier(MPI_COMM_WORLD); if (myid == 0) { std::cout << " done." << std::endl; } } } MFEM_VERIFY(nxyz.Size() == 0 || nxyz.Size() == mesh->SpaceDimension(), "Expected " << mesh->SpaceDimension() << " integers with the " "option --cartesian-partitioning."); int *partitioning = nxyz.Size() ? mesh->CartesianPartitioning(nxyz) : NULL; NewParMesh(pmesh, mesh, partitioning); delete [] partitioning; { for (int l = 0; l < par_ref_levels; l++) { if (myid == 0) { std::cout << "Parallel refinement: level " << l << " -> level " << l+1 << " ..." << std::flush; } pmesh->UniformRefinement(); MPI_Barrier(MPI_COMM_WORLD); if (myid == 0) { std::cout << " done." << std::endl; } } } pmesh->PrintInfo(std::cout); }

   const int p = 5;
   const int dim = 3;
   auto &mesh = xfl::Mesh(pmesh);
   const int el = (dim==2)?xfl::quadrilateral:xfl::hexahedron;
   FiniteElementCollection *fe = xfl::FiniteElement("Lagrange", el, p);
   mfem::ParFiniteElementSpace *fes = xfl::FunctionSpace(mesh, fe);
   const Array<int> bc = xfl::DirichletBC(fes);
   xfl::Function x = xfl::Function(fes);
   x = 0.0;
   xfl::TrialFunction u = xfl::TrialFunction(fes);
   xfl::TestFunction v = xfl::TestFunction(fes);
   auto b = [&]()
   {
      constexpr const char *qs0 = "v";
      // var:[v], ops:[Xe,v,]
      // Test FES: 'fes':v (Eval)
      ParFiniteElementSpace *fes0 = fes;
      mfem::Operator *QM0 = nullptr;
      xfl::QForm QForm0(fes0, qs0, QM0);
      return QForm0;
   }();
   auto a = [&]()
   {
      constexpr const char *qs1 = "dot(grad(u), grad(v))";
      // var:[u,v], ops:[Xe,u,G,*,D,*,v,G,T,*,Ye]
      // Trial FES: 'fes':u (Grad)
      // Test FES: 'fes':v (Grad)
      ParFiniteElementSpace *fes1 = fes;
      struct QMult1: public xfl::Operator<3>
      {
         QMult1(const ParFiniteElementSpace *fes): xfl::Operator<3>(fes) { Setup(); }
         ~QMult1() { }
         void Setup()
         {
            int myid = 0; MPI_Comm_rank(MPI_COMM_WORLD, &myid);
            if (myid == 0)
            {
               std::cout << "XFL(SIMD_" << SIMD_SIZE
                         <<") version using integration rule with 343 points ...\n";
               std::cout << "D1D:6, Q1D:7\n";
            }
            dx.SetSize(NQ*NE*3*3, Device::GetDeviceMemoryType()); // DX shape: 3x3
            KSetup1<3,3,3,7>(NDOFS, VDIM, NE, J0.Read(), ir.GetWeights().Read(),
                             dx.Write());
            CoG.SetSize(7*7);
            // Compute the collocated gradient d2q->CoG
            kernels::GetCollocatedGrad<6,7>(
               ConstDeviceMatrix(maps->B.HostRead(),7,6),
               ConstDeviceMatrix(maps->G.HostRead(),7,6),
               DeviceMatrix(CoG.HostReadWrite(),7,7));
         }
         void Mult(const mfem::Vector &x, mfem::Vector &y) const
         {
            y = 0.0;
            KMult1<3,3,3,6,7>(NDOFS /*0*/,VDIM /*0*/, NE /*0*/,maps->B.Read(), CoG.Read(),
                              ER.GatherMap().Read(), dx.Read(), x.Read(), y.ReadWrite());
         }
      }; // QMult struct
      QMult1 *QM1 = new QMult1(fes);
      xfl::QForm QForm1(fes1, qs1, QM1);
      return QForm1;
   }();
   status |= xfl::benchmark(a==b, x, bc, 0, 50, 3);
   const bool glvis = true;
   if (glvis) { xfl::plot(x); }
   assert(SOL_P == p);
   MPI_Finalize();
   return status;
}

