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

/// Use optimized CG
#define MFEM_OPTIMIZED_CG

void OptimizedCG(const FiniteElementSpace &fes,
                 const mfem::Operator &A,
                 const Vector &B,
                 Vector &X,
                 const Array<int> &ess_tdof_list,
                 const Vector &pa_data,
                 int &cg_iter, double &real_time,
                 int print_iter = 0, int max_num_iter = 1000,
                 double RTOLERANCE = 1e-12, double ATOLERANCE = 1e-24);

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

   XElementRestriction(const FiniteElementSpace &fes,
                       ElementDofOrdering ordering)
      : ElementRestriction(fes, ordering) { }

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
        geom(mesh->GetGeometricFactors(ir, flags)),
        maps(&pfes->GetFE(0)->GetDofToQuad(ir, mode)),
        nqi(nfes->GetQuadratureInterpolator(ir)),
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
        geom(mesh->GetGeometricFactors(ir, flags)),
        maps(&pfes->GetFE(0)->GetDofToQuad(ir, mode)),
        qi(pfes->GetQuadratureInterpolator(ir)),
        nqi(nfes->GetQuadratureInterpolator(ir)),
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

   Vector& GetData() { return dx; }
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
   const Vector &dx;

public:

   QForm(ParFiniteElementSpace *pfes,
         const char *qs,
         mfem::Operator *QM,
         const Vector &dx = Vector())
      : qs(qs), QM(QM), pfes(pfes), dx(dx) { }

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
      MFEM_CONTRACT_VAR(rhs);
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
   double operator*(TestFunction &v) { MFEM_CONTRACT_VAR(v); return 0.0; }
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
   Element::Type quad = Element::Type::QUADRILATERAL;
   mfem::Mesh mesh = Mesh::MakeCartesian2D(nx, ny, quad);
   return *MeshToPMesh(&mesh);
}

mfem::ParMesh &UnitHexMesh(int nx, int ny, int nz)
{
   Element::Type hex = Element::Type::HEXAHEDRON;
   mfem::Mesh mesh = Mesh::MakeCartesian3D(nx, ny, nz, hex);
   return *MeshToPMesh(&mesh);
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
int benchmark(xfl::Problem *pb,
              xfl::Function &x,
              Array<int> ess_tdof_list,
              const Vector &pa_data,
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

#ifndef MFEM_OPTIMIZED_CG
#ifndef MFEM_USE_MPI
   CGSolver cg;
#else
   CGSolver cg(MPI_COMM_WORLD);
#endif


   cg.SetRelTol(rtol);
   cg.SetAbsTol(0.0);
   cg.SetOperator(*A);
   cg.iterative_mode = false;

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
   const HYPRE_Int dofs = pfes->GlobalTrueVSize();
   const int cg_iter = cg.GetNumIterations();
#else
   int final_iter;
   double real_time;
   OptimizedCG(*pfes,
               *A, B, X,
               ess_tdof_list, pa_data,
               final_iter, real_time,
               print_lvl, max_it, rtol, 0.0);
   const HYPRE_Int dofs = pfes->GlobalTrueVSize();
   const double mdofs = ((1e-6 * dofs) * max_it) / real_time;
   mfem::out << "\"DOFs/sec\" in CG: " << mdofs << " million.\n";
   const int cg_iter = max_it;
   assert(false);
#endif

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

int benchmark(xfl::Problem *pb,
              xfl::Function &x,
              const Vector &pa_data,
              Array<int> ess_tdof_list)
{
   return benchmark(pb, x, ess_tdof_list, pa_data, 1e-12, 200, -1);
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

// Foreach Thread
#define FOREACH_THREAD(i,k,N) \
    PRAGMA(unroll(N))\
    for(int i=0; i<N; ++i)
#define SYNC_THREADS

#define D1D (SOL_P + 1)
#define Q1D (SOL_P + 2)

template<int DIM, int DX0, int DX1> inline static
void KSetup1(const int NE,
             const double * __restrict__ J0,
             const double * __restrict__ w,
             double * __restrict__ dx)
{
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

template<int DIM, int DX0, int DX1> inline static
void KMult1(const int ndofs,
            const int NE,
            const double * __restrict__ B_,
            const double * __restrict__ G_,
            const int * __restrict__ map,
            const double * __restrict__ dx,
            const double * __restrict__ xd,
            double * __restrict__ yd)
{
   // kernel operations: u,G,*,D,*,v,G,T,*,Ye
   const auto B = Reshape(B_, Q1D,D1D);
   const auto G = Reshape(G_, Q1D,Q1D);
   const auto DX = Reshape(dx, Q1D,Q1D,Q1D, 6);
   const auto MAP = Reshape(map, D1D,D1D,D1D, NE);
   const auto XD = Reshape(xd, ndofs);
   auto YD = Reshape(yd, ndofs);

   MFEM_ALIGN(Real) s_q[Q1D][Q1D][Q1D];
   MFEM_ALIGN(double) s_G[Q1D][Q1D];
   MFEM_ALIGN(double) s_B[Q1D][D1D];

   MFEM_ALIGN(Real) s_Gqr[Q1D][D1D];
   MFEM_ALIGN(Real) s_Gqs[Q1D][D1D];

   MFEM_ALIGN(Real) l_q[Q1D];
   MFEM_ALIGN(Real) r_qt[Q1D]/**/[Q1D]/**/;
   MFEM_ALIGN(Real) r_q[Q1D]/**/[Q1D][Q1D]/**/;

   for (int e = 0; e < NE; e+=SMS)
   {
      // Scatter X
      FOREACH_THREAD(b,y,Q1D)
      {
         FOREACH_THREAD(a,x,Q1D)
         {
            s_G[b][a] = G(b,a);
            if (a<D1D) { s_B[b][a] = B(b,a); }
            if (a<D1D && b<D1D)
            {
               UNROLL(D1D)
               for (int c=0; c<D1D; ++c)
               {
                  MFEM_ALIGN(Real) vXD;
                  UNROLL(SMS)
                  for (size_t v = 0; v < SMS; v++)
                  {
                     const int gid = MAP(a,b,a,e+v);
                     vXD[v] = XD(gid);
                  }
                  s_q[c][b][a] = vXD;
               }
            }
         }
      } SYNC_THREADS;

      // Grad1Z
      FOREACH_THREAD(b,y,D1D)
      {
         FOREACH_THREAD(a,x,D1D)
         {
            UNROLL(Q1D)
            for (int k=0; k<Q1D; ++k)
            {
               MFEM_ALIGN(Real) u; u = 0.0;
               UNROLL(D1D)
               for (int c=0; c<D1D; ++c) { u.fma(s_B[k][c], s_q[c][b][a]); }
               s_q[k][b][a] = u;
            }
         }
      } SYNC_THREADS;

      // Grad1Y
      FOREACH_THREAD(k,y,Q1D)
      {
         FOREACH_THREAD(a,x,D1D)
         {
            UNROLL(D1D)
            for (int b=0; b<D1D; ++b) { l_q[b] = s_q[k][b][a]; }
            UNROLL(Q1D)
            for (int j=0; j<Q1D; ++j)
            {
               MFEM_ALIGN(Real) u; u = 0.0;
               UNROLL(D1D)
               for (int b=0; b<D1D; ++b) { u.fma(s_B[j][b], l_q[b]); }
               s_q[k][j][a] = u;
            }
         }
      } SYNC_THREADS;

      // Grad1X
      FOREACH_THREAD(k,y,Q1D)
      {
         FOREACH_THREAD(j,x,Q1D)
         {
            for (int a=0; a<D1D; ++a) { l_q[a] = s_q[k][j][a]; }
            UNROLL(Q1D)
            for (int i=0; i<Q1D; ++i)
            {
               MFEM_ALIGN(Real) u; u = 0.0;
               UNROLL(D1D)
               for (int a=0; a<D1D; ++a) { u.fma(s_B[i][a], l_q[a]); }
               s_q[k][j][i] = u;
            }
         }
      } SYNC_THREADS;

      // Flush
      FOREACH_THREAD(j,y,Q1D)
      {
         FOREACH_THREAD(i,x,Q1D)
         {
            UNROLL(Q1D)
            for (int k = 0; k < Q1D; ++k) { r_q[i][j][k] = 0.0; }
         }
      } SYNC_THREADS;

      // Q-Function
      UNROLL(Q1D)
      for (int k = 0; k < Q1D; ++k)
      {
         SYNC_THREADS;
         FOREACH_THREAD(j,y,Q1D)
         {
            FOREACH_THREAD(i,x,Q1D)
            {
               MFEM_ALIGN(Real) qr, qs, qt;
               qr = 0.0; qs = 0.0; qt = 0.0;
               UNROLL(Q1D)
               for (int m = 0; m < Q1D; ++m)
               {
                  const double Dim = s_G[i][m];
                  const double Djm = s_G[j][m];
                  const double Dkm = s_G[k][m];
                  qr.fma(Dim, s_q[k][j][m]);
                  qs.fma(Djm, s_q[k][m][i]);
                  qt.fma(Dkm, s_q[m][j][i]);
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

               s_Gqr[j][i] = 0.0;
               s_Gqr[j][i].fma(D00,qr);
               s_Gqr[j][i].fma(D10,qs);
               s_Gqr[j][i].fma(D20,qt);

               s_Gqs[j][i] = 0.0;
               s_Gqs[j][i].fma(D01,qr);
               s_Gqs[j][i].fma(D11,qs);
               s_Gqs[j][i].fma(D21,qt);

               r_qt[j][i] = 0.0;
               r_qt[j][i].fma(D02,qr);
               r_qt[j][i].fma(D12,qs);
               r_qt[j][i].fma(D22,qt);
            }
         }
         SYNC_THREADS;

         FOREACH_THREAD(j,y,Q1D)
         {
            FOREACH_THREAD(i,x,Q1D)
            {
               MFEM_ALIGN(Real) u; u = 0.0;
               UNROLL(Q1D)
               for (int m = 0; m < Q1D; ++m)
               {
                  const double Dmi = s_G[m][i];
                  const double Dmj = s_G[m][j];
                  const double Dkm = s_G[k][m];
                  u.fma(Dmi, s_Gqr[j][m]);
                  u.fma(Dmj, s_Gqs[m][i]);
                  r_q[j][i][m].fma(Dkm, r_qt[j][i]);
               }
               r_q[i][j][k] += u;
            }
         } SYNC_THREADS;
      }

      // GradZT
      FOREACH_THREAD(j,y,Q1D)
      {
         FOREACH_THREAD(i,x,Q1D)
         {
            UNROLL(D1D)
            for (int c=0; c<D1D; ++c)
            {
               MFEM_ALIGN(Real) u; u = 0.0;
               UNROLL(Q1D)
               for (int k=0; k<Q1D; ++k) { u.fma(s_B[k][c], r_q[i][j][k]); }
               s_q[c][j][i] = u;
            }
         }
      } SYNC_THREADS;

      // GradYT
      FOREACH_THREAD(c,y,D1D)
      {
         FOREACH_THREAD(i,x,Q1D)
         {
            UNROLL(Q1D)
            for (int j=0; j<Q1D; ++j) { l_q[j] = s_q[c][j][i]; }
            UNROLL(D1D)
            for (int b=0; b<D1D; ++b)
            {
               MFEM_ALIGN(Real) u; u = 0.0;
               UNROLL(Q1D)
               for (int j=0; j<Q1D; ++j) { u.fma(s_B[j][b], l_q[j]); }
               s_q[c][b][i] = u;
            }
         }
      } SYNC_THREADS;

      // GradXT
      FOREACH_THREAD(c,y,D1D)
      {
         FOREACH_THREAD(b,x,D1D)
         {
            UNROLL(Q1D)
            for (int i=0; i<Q1D; ++i) { l_q[i] = s_q[c][b][i]; }
            UNROLL(D1D)
            for (int a=0; a<D1D; ++a)
            {
               MFEM_ALIGN(Real) u; u = 0.0;
               UNROLL(Q1D)
               for (int i=0; i<Q1D; ++i) { u.fma(s_B[i][a], l_q[i]); }
               s_q[c][b][a] = u;
            }
         }
      } SYNC_THREADS;

      // Gather
      FOREACH_THREAD(b,y,D1D)
      {
         FOREACH_THREAD(a,x,D1D)
         {
            UNROLL(D1D)
            for (int c = 0; c < D1D; c++)
            {
               UNROLL(SMS)
               for (size_t v = 0; v < SMS; v++)
               {
                  const int gid = MAP(a,b,c,e+v);
                  YD(gid) += (s_q[c][b][a])[v];
               }
            }
         }
      }
   } // MFEM_FORALL
} // KMult1
#undef D1D
#undef Q1D


/// Optimized CG solver

#define MFEM_FORALL_GRID_1D(i,N) for(int i=0; i<N; ++i)

/// OptDot
template<int Q1D> static
double OptDot(const double * __restrict__ A,
              const double * __restrict__ B,
              double *__restrict__ dot_result,
              const int N)
{
   *dot_result = 0.0;
   MFEM_FORALL_GRID_1D(i,N) { *dot_result += A[i] * B[i]; }
   return *dot_result;
}

/// OperMult
template<int D1D, int Q1D, int NBZ> static
void OperMult(const int N,
              const int csz,
              const int *__restrict__ ess_tdof_list,
              const int NE,
              const DeviceTensor<4,const int> MAP,
              const DeviceTensor<2,const double> B,
              const DeviceTensor<4,const double> Q,
              const double *__restrict__ d,
              double *__restrict__ w,
              double *__restrict__ z)
{
   MFEM_FORALL_GRID_1D(i,N) { w[i] = d[i]; z[i] = 0.0; }
   MFEM_FORALL_GRID_1D(i,csz) { w[ess_tdof_list[i]] = 0.0; }
   // kPAMassApply3D<D1D,Q1D,NBZ>(NE,MAP,B,Q,w,z);
#warning here
   /*
   KMult1<3,3,3>(N, NE,
                 maps->B.Read(),
                 CoG.Read(),
                 MAP,
                 pa_data, w, z);
                 */
   MFEM_FORALL_GRID_1D(i,csz) { z[ess_tdof_list[i]] = d[ess_tdof_list[i]]; }
}

/// Optimized CG Solver
template<int D1D, int Q1D, int NBZ>  static
void OptCGSolver(const int N,
                 /* CG Solver */
                 const double rel_tol,
                 const double abs_tol,
                 double *__restrict__ dot_result,
                 int max_iter,
                 int print_level,
                 double *__restrict__ r,
                 double *__restrict__ d,
                 double *__restrict__ z,
                 double *__restrict__ w,
                 const double *__restrict__ b,
                 double *__restrict__ x,
                 /* ConstrainedOperator */
                 const int csz,
                 const int *__restrict__ etdl,
                 /* SmRgPAMassApply3D */
                 const int NE,
                 const DeviceTensor<4,const int> MAP,
                 const DeviceTensor<2,const double> B,
                 const DeviceTensor<4,const double> Q)
{
   MFEM_CONTRACT_VAR(rel_tol);
   MFEM_CONTRACT_VAR(abs_tol);
   int iter;
   double betanom;

   MFEM_FORALL_GRID_1D(i,N)
   {
      x[i] = 0.0;  // x = 0.0;
      d[i] = r[i] = b[i]; // d = r = b;
   }

   double nom = OptDot<Q1D>(d,r,dot_result, N);

   if (print_level > 0)
   {
      printf("   Iteration : %3d  (B r, r) = %.5e\n", 0, nom);
   }

   OperMult<D1D,Q1D,NBZ>(N,csz,etdl,NE,MAP,B,Q,d,w,z); // 1ms order 6

   double den = OptDot<Q1D>(z,d,dot_result,N); // .2

   for (iter = 1; iter <= max_iter; ++iter)
   {
      const double alpha = nom / den;
      MFEM_FORALL_GRID_1D(i,N)
      {
         x[i] += alpha * d[i]; //  x = x + alpha d
         r[i] -= alpha * z[i]; //  r = r - alpha A d
      }
      betanom = OptDot<Q1D>(r,r,dot_result,N);
      const double beta = betanom / nom;
      MFEM_FORALL_GRID_1D(i,N) { d[i] = r[i] + beta * d[i]; }
      OperMult<D1D,Q1D,NBZ>(N,csz,etdl,NE,MAP,B,Q,d,w,z);
      den = OptDot<Q1D>(d,z,dot_result,N);
      nom = betanom;
   }

   if (print_level > 0)
   {
      printf("   Iteration : %3d  (B r, r) = %.5e\n", iter-1, betanom);
   }
}

/// Optimized Conjugate gradient method
class OptimizedCGSolver : public IterativeSolver
{
protected:
   const FiniteElementSpace &fes;
   const Array<int> &ess_tdof_list;
   const Vector &pa_data;
   mfem::Mesh *mesh;
   const mfem::FiniteElement *el;
   ElementTransformation *Tr;
   const int qorder;
   const IntegrationRule &ir;
   const DofToQuad &maps;
   const unsigned int D1D, Q1D, ND, NE;
   mutable Vector r, d, z, w;

   void UpdateVectors()
   {
      MemoryType mt = GetMemoryType(oper->GetMemoryClass());
      r.SetSize(width, mt); r.UseDevice(true);
      d.SetSize(width, mt); d.UseDevice(true);
      z.SetSize(width, mt); z.UseDevice(true);
      w.SetSize(width, mt); w.UseDevice(true);
   }

public:
   OptimizedCGSolver(const FiniteElementSpace &fes,
                     const Array<int> &ess_tdof_list,
                     const Vector &pa_data):
      fes(fes),
      ess_tdof_list(ess_tdof_list),
      pa_data(pa_data),
      mesh(fes.GetMesh()),
      el(fes.GetFE(0)),
      Tr(mesh->GetElementTransformation(0)),
      qorder(el->GetOrder() + el->GetOrder() + Tr->OrderW()),
      ir(IntRules.Get(el->GetGeomType(), qorder)),
      maps(el->GetDofToQuad(ir, DofToQuad::TENSOR)),
      D1D(maps.ndof),
      Q1D(maps.nqpt),
      ND(fes.GetNDofs()),
      NE(mesh->GetNE()) { /**/ }

   virtual void SetOperator(const Operator &op)
   {
      IterativeSolver::SetOperator(op);
      UpdateVectors();
   }

   virtual void Mult(const Vector &b, Vector &x) const
   {
      assert(!monitor);
      assert(!iterative_mode);
      assert(!prec);

      final_iter = max_iter;
      OptLaunch(b, x);
   }

private:
   void OptLaunch(const Vector &b, Vector &x) const
   {
      const int N = r.Size();

      double *R = r.ReadWrite();
      double *D = d.ReadWrite();
      double *Z = z.ReadWrite();
      double *W = w.ReadWrite();

      const double *B = b.Read();
      double *X = x.ReadWrite();

      const int csz = ess_tdof_list.Size();
      const int *ETDL = ess_tdof_list.Read();

      constexpr ElementDofOrdering ordering = ElementDofOrdering::LEXICOGRAPHIC;
      //const Operator *ERop = fes.GetElementRestriction(ordering);
      XElementRestriction XER(fes, ordering);
      const int *map = XER.GatherMap().Read();
      const DeviceTensor<4,const int> &MAP = Reshape(map, D1D,D1D,D1D, NE);
      const ConstDeviceMatrix &Bqd = Reshape(maps.B.Read(), Q1D, D1D);
      const DeviceTensor<4,const double> &Dq =
         Reshape(pa_data.Read(), Q1D,Q1D,Q1D,NE);

      void (*Kernel)(const int N,
                     /* CG Solver */
                     const double rel_tol,
                     const double abs_tol,
                     double *dot_result,
                     int max_iter,
                     int print_level,
                     double *r,
                     double *d,
                     double *z,
                     double *w,
                     const double *b,
                     double *x,
                     /* ConstrainedOperator */
                     const int csz,
                     const int *ess_tdof_list,
                     /* SmRgPAMassApply3D */
                     const int NE,
                     const DeviceTensor<4,const int> MAP,
                     const DeviceTensor<2,const double> Bqd,
                     const DeviceTensor<4,const double> Dq) = nullptr;

      const int id = (D1D << 4) | Q1D;
      switch (id) // orders 1~6
      {
         case 0x23: Kernel=OptCGSolver<2,3,32>; break; // 1
         case 0x34: Kernel=OptCGSolver<3,4,16>; break; // 2
         case 0x45: Kernel=OptCGSolver<4,5, 4>; break; // 3
         case 0x56: Kernel=OptCGSolver<5,6, 2>; break; // 4
         case 0x67: Kernel=OptCGSolver<6,7, 2>; break; // 5
         case 0x78: Kernel=OptCGSolver<7,8, 2>; break; // 6
         default: MFEM_ABORT("Unknown kernel 0x" << std::hex << id);
      }


      static Vector dot_result(1, Device::GetDeviceMemoryType());
      dot_result.UseDevice(true);
      double *device_dot = dot_result.Write();

      assert(Kernel);
      Kernel(N,
             rel_tol, abs_tol,
             device_dot,
             max_iter, print_level,
             R,D,Z,W,B,X,
             csz,ETDL,
             NE, MAP, Bqd, Dq);
   }
};

/// Conjugate gradient method on device. (tolerances are squared)
void OptimizedCG(const FiniteElementSpace &fes,
                 const mfem::Operator &A,
                 const Vector &B,
                 Vector &X,
                 const Array<int> &ess_tdof_list,
                 const Vector &pa_data,
                 int &cg_iter, double &real_time,
                 int print_iter = 0, int max_num_iter = 1000,
                 double RTOLERANCE = 1e-12, double ATOLERANCE = 1e-24)
{
   OptimizedCGSolver cg(fes, ess_tdof_list, pa_data);
   cg.SetRelTol(RTOLERANCE);
   cg.SetAbsTol(sqrt(ATOLERANCE));
   cg.SetOperator(A);
   cg.iterative_mode = false;

   // Warm-up CG solve (in case of JIT to avoid timing it)
   {
      Vector Y(X);
      cg.SetMaxIter(2);
      cg.SetPrintLevel(-1);
      cg.Mult(B, Y);
      MFEM_DEVICE_SYNC;
   }

   {
      tic_toc.Clear();
      cg.SetPrintLevel(print_iter);
      cg.SetMaxIter(max_num_iter);
      {
         MFEM_DEVICE_SYNC;
         tic_toc.Start();
         cg.Mult(B, X);
         MFEM_DEVICE_SYNC;
         tic_toc.Stop();
      }
   }
   real_time = tic_toc.RealTime();
   cg_iter = cg.GetNumIterations();
   assert(cg_iter==max_num_iter);
}


/// MAIN
int main(int argc, char* argv[])
{
   assert(MESH_P==1);
   assert(IR_TYPE==0);
   assert(IR_ORDER==0);
   assert(PROBLEM==0);
   assert(GEOM==Geometry::CUBE);

   int status = 0;
   int num_procs = 1, myid = 0;

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);

#ifdef MFEM_OPTIMIZED_CG
   assert(num_procs==1); // serial for now
#endif

   const char *mesh_file = "../../data/hex-01x01x01.mesh";
   int ser_ref_levels = 4;
   int par_ref_levels = 0;
   Array<int> nxyz; int order = SOL_P;
   int max_iter = 50;
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
   args.AddOption(&max_iter, "-mi", "--max-iter", "Maximum number of iterations.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization", "Enable or disable GLVis visualization."); args.Parse();

   if (!args.Good()) { if (myid == 0) { args.PrintUsage(std::cout); } return 1; }

   if (myid == 0) { args.PrintOptions(std::cout); }


   ParMesh *pmesh = nullptr;
   {
      Mesh *mesh = new Mesh(mesh_file, 1, 1);
      const int dim = mesh->Dimension();
      {
         int ref_levels = (int)floor(log(10000./mesh->GetNE())/log(2.)/dim);
         ref_levels = (ser_ref_levels != -1) ? ser_ref_levels : ref_levels;
         for (int l = 0; l < ref_levels; l++)
         {
            if (myid == 0) { std::cout << "Serial refinement: level " << l << " -> level " << l+1 << " ..." << std::flush; } mesh->UniformRefinement();
            MPI_Barrier(MPI_COMM_WORLD);
            if (myid == 0) { std::cout << " done." << std::endl; }
         }
      }

      MFEM_VERIFY(nxyz.Size() == 0 ||
                  nxyz.Size() == mesh->SpaceDimension(),
                  "Expected " << mesh->SpaceDimension() << " integers with the "
                  "option --cartesian-partitioning.");
      int *partitioning = nxyz.Size() ? mesh->CartesianPartitioning(nxyz) : NULL;
      NewParMesh(pmesh, mesh, partitioning);
      delete [] partitioning;

      {
         for (int l = 0; l < par_ref_levels; l++)
         {
            if (myid == 0)
            {
               std::cout << "Parallel refinement: level " << l
                         << " -> level " << l+1 << " ..."
                         << std::flush;
            }
            pmesh->UniformRefinement();
            MPI_Barrier(MPI_COMM_WORLD);
            if (myid == 0) { std::cout << " done." << std::endl; }
         }
      }
      pmesh->PrintInfo(std::cout);
   }

   const int dim = 3;
   const int p = SOL_P;
   assert(SOL_P == order);
   auto &mesh = xfl::Mesh(pmesh);

   const int NE = mesh.GetNE();
   //int NE_min;
   //MPI_Reduce(&NE, &NE_min, 1, MPI_INT, MPI_MIN, 0, mesh.GetComm());
   assert((NE % SIMD_SIZE) == 0);

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
      //  Test FES: 'fes':v (Grad)
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
                         <<") version using integration rule with "
                         << (Q1D*Q1D*Q1D) <<" points ...\n";
               std::cout << "D1D:"<<D1D<<", Q1D:"<<Q1D<<"\n";
            }
            dx.SetSize(NQ*3*3, Device::GetDeviceMemoryType());
            KSetup1<3,3,3>(NE, J0.Read(), ir.GetWeights().Read(), dx.Write());
            CoG.SetSize(Q1D*Q1D);
            // Compute the collocated gradient d2q->CoG
            kernels::GetCollocatedGrad<D1D,Q1D>(
               ConstDeviceMatrix(maps->B.HostRead(),Q1D,D1D),
               ConstDeviceMatrix(maps->G.HostRead(),Q1D,D1D),
               DeviceMatrix(CoG.HostReadWrite(),Q1D,Q1D));
         }
         void Mult(const mfem::Vector &x, mfem::Vector &y) const
         {
            y = 0.0;
            KMult1<3,3,3>(NDOFS, NE,
                          maps->B.Read(), CoG.Read(),
                          ER.GatherMap().Read(),
                          dx.Read(), x.Read(), y.ReadWrite());
         }
      }; // QMult struct
      QMult1 *QM1 = new QMult1(fes);
      xfl::QForm QForm1(fes1, qs1, QM1, QM1->GetData());
      return QForm1;
   }();
   status |= xfl::benchmark(a==b, x, bc, a.dx, 0, max_iter, 3);
   if (visualization) { xfl::plot(x); }
   MPI_Finalize();
   return status;
}

