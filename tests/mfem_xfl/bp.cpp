// Copyright (c) 2010-2020, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>
#include <typeindex>
#include <utility>
#include <vector>
#include <thread>

#include "mfem.hpp"
#define dbg(...)

#include "general/forall.hpp"
#include "linalg/kernels.hpp"


#include "linalg/simd.hpp"
#define SIMD_SIZE (MFEM_SIMD_BYTES / sizeof(double))

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
#define IR_ORDER 0
#endif

#ifndef IR_TYPE
#define IR_TYPE 0
#endif

#ifndef PROBLEM
#define PROBLEM 0
#endif

#define MFEM_ALIGN alignas(32)

using namespace mfem;

namespace mfem
{

using FE = FiniteElement;
using QI = QuadratureInterpolator;
using Real = AutoSIMDTraits<double, double>::vreal_t;
using ConstDeviceMatrix = DeviceTensor<2,const double>;

// XFL addons //////////////////////////////////////////////////////////////////
namespace xfl
{

class XElementRestriction : public ElementRestriction
{
   const ParFiniteElementSpace &fes;

public:
   XElementRestriction(const ParFiniteElementSpace *fes, ElementDofOrdering edo)
      : ElementRestriction(*fes, edo), fes(*fes) {}

   const Array<int> &ScatterMap() const { return offsets; }
   const Array<int> &ScatterIdx() const { return indices; }
   const Array<int> &GatherMap() const { return gatherMap; }

   void Mult(const Vector &x, Vector &y) const
   {
      const int ndof = fes.GetFE(0)->GetDof();
      const int ndofs = fes.GetNDofs();

      const auto d_x = Reshape(x.Read(), ndofs);
      const auto d_j = gatherMap.Read();

      auto d_y = Reshape(y.Write(), ndof, ne);

      MFEM_FORALL(i, ndof * ne, { d_y(i % ndof, i / ndof) = d_x(d_j[i]); });
   }

   void MultTranspose(const Vector &x, Vector &y) const
   {
      const int nd = fes.GetFE(0)->GetDof();

      const auto d_offsets = offsets.Read();
      const auto d_indices = indices.Read();

      const auto d_x = Reshape(x.Read(), nd, ne);
      auto d_y = Reshape(y.Write(), ndofs);
      MFEM_FORALL(i, ndofs,
      {
         double dofValue = 0.0;
         const int offset = d_offsets[i];
         const int nextOffset = d_offsets[i + 1];
         for (int j = offset; j < nextOffset; ++j)
         {
            const bool plus = d_indices[j] >= 0;
            const int idx_j = plus ? d_indices[j] : -1 - d_indices[j];
            const double value = d_x(idx_j % nd, idx_j / nd);
            dofValue += (plus) ? value : -value;
         }
         d_y(i) = dofValue;
      });
   }
};

/** ****************************************************************************
 * @brief The Operator class
 **************************************************************************** */
template <int DIM>
class Operator;

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
        geom(mesh->GetGeometricFactors(ir, flags)),//, mode)),
        maps(&pfes->GetFE(0)->GetDofToQuad(ir, mode)),
        qi(pfes->GetQuadratureInterpolator(ir)),//, mode)),
        nqi(nfes->GetQuadratureInterpolator(ir)),//, mode)),
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
      dbg("VDIM:%d, SDIM:%d", VDIM, SDIM);
      dbg("p:%d q:%d OrderW:%d, D1D:%d, Q1D:%d", p, q,
          mesh->GetElementTransformation(0)->OrderW(), D1D, Q1D);
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
   mfem::Operator *QM{nullptr};
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
   // Constructor
   QForm(mfem::ParFiniteElementSpace *pfes, const char *qs, mfem::Operator *QM)
      : qs(qs), QM(QM), pfes(pfes) {}

   ~QForm() {}

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
      else {  assert(false); }
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
   const ParFiniteElementSpace *ParFESpace() const
   {
      return ParGridFunction::ParFESpace();
   }
   ConstantCoefficient *ConstantCoeff() const { return nullptr; }
   FunctionCoefficient *FunctionCoeff() const { return nullptr; }
};

/** ****************************************************************************
 * @brief TrialFunction class
 ******************************************************************************/
class TrialFunction : public Function
{
public:
   TrialFunction(ParFiniteElementSpace *pfes) : Function(pfes) {}
   ~TrialFunction() {}
};

/** ****************************************************************************
 * @brief TestFunction class
 ******************************************************************************/
class TestFunction : public Function
{
public:
   TestFunction(ParFiniteElementSpace *pfes) : Function(pfes) {}
   TestFunction(const TestFunction &) = default;
   ~TestFunction() {}
};

/** ****************************************************************************
 * @brief Constant class
 ******************************************************************************/
class Constant
{
   const double value = 0.0;
   ConstantCoefficient *cst = nullptr;

public:
   Constant(double val) : value(val), cst(new ConstantCoefficient(val)) {}
   ~Constant() { dbg(); /*delete cst;*/ }
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
      : fct(new FunctionCoefficient(F)) {}
   ~Expression()   /*delete fct;*/
   {
   }
   ParFiniteElementSpace *ParFESpace() const { return nullptr; }  // qf args
   ConstantCoefficient *ConstantCoeff() const { return nullptr; }
   FunctionCoefficient *FunctionCoeff() const { return fct; }
   operator const double *() const { return nullptr; }  // qf eval
   // double operator *(TestFunction &v) { dbg(); return 0.0;} // qf eval
};

/** ****************************************************************************
 * @brief Mesh
 ******************************************************************************/
mfem::ParMesh &Mesh(mfem::ParMesh *pmesh) { return *pmesh; }

/** ****************************************************************************
 * @brief FiniteElement
 ******************************************************************************/
FiniteElementCollection *FiniteElement(std::string family, int type, int p)
{
   MFEM_VERIFY(family == "Lagrange", "Unsupported family!");
   MFEM_VERIFY(
      type == Element::Type::QUADRILATERAL || type == Element::Type::HEXAHEDRON,
      "Unsupported type!");
   const int dim = (type == Element::Type::QUADRILATERAL) ? 2
                   : (type == Element::Type::HEXAHEDRON)  ? 3
                   : 0;
   const int btype = BasisType::GaussLobatto;
   return new H1_FECollection(p, dim, btype);
}

/** ****************************************************************************
 * @brief Function Spaces
 ******************************************************************************/
class FunctionSpace : public ParFiniteElementSpace {};

ParFiniteElementSpace *FunctionSpace(mfem::ParMesh *pmesh, std::string family,
                                     int p)
{
   assert(false);
   const int dim = pmesh->Dimension();
   MFEM_VERIFY(family == "P", "Unsupported FE!");
   FiniteElementCollection *fec = new H1_FECollection(p, dim);
   return new ParFiniteElementSpace(pmesh, fec);
}

ParFiniteElementSpace *FunctionSpace(mfem::ParMesh &pmesh, std::string f,
                                     int p)
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
                                           std::string family, const int p)
{
   const int dim = pmesh->Dimension();
   MFEM_VERIFY(family == "P", "Unsupported FE!");
   FiniteElementCollection *fec = new H1_FECollection(p, dim);
   return new ParFiniteElementSpace(pmesh, fec, dim);
}

ParFiniteElementSpace *VectorFunctionSpace(mfem::ParMesh &pmesh,
                                           std::string family, const int p)
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

   dbg("second solve with the inline QMult");
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

   CGSolver cg(MPI_COMM_WORLD);
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

   int myid;
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   MPI_Comm comm = pfes->GetParMesh()->GetComm();

   const double rt = tic_toc.RealTime();
   double rt_min, rt_max;
   MPI_Reduce(&rt, &rt_min, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
   MPI_Reduce(&rt, &rt_max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

   const HYPRE_Int dofs = pfes->GlobalTrueVSize();
   const int cg_iter = cg.GetNumIterations();
   const double mdofs_max = ((1e-6 * dofs) * cg_iter) / rt_max;
   const double mdofs_min = ((1e-6 * dofs) * cg_iter) / rt_min;

   if (myid == 0)
   {
      std::cout << "Number of finite element unknowns: " << dofs << std::endl;
      std::cout << "Total CG time:    " << rt_max << " (" << rt_min << ") sec."
                << std::endl;
      std::cout << "Time per CG step: " << rt_max / cg_iter << " ("
                << rt_min / cg_iter << ") sec." << std::endl;
      std::cout << "\"DOFs/sec\" in CG: " << mdofs_max << " (" << mdofs_min
                << ") million.";
      std::cout << std::endl;
   }
   delete pb;
   return 0;
}

int benchmark(xfl::Problem *pb, xfl::Function &x, Array<int> ess_tdof_list)
{
   return benchmark(pb, x, ess_tdof_list, 1e-12, 200, -1);
}

constexpr int point = Element::Type::POINT;
constexpr int segment = Element::Type::SEGMENT;
constexpr int triangle = Element::Type::TRIANGLE;
constexpr int quadrilateral = Element::Type::QUADRILATERAL;
constexpr int tetrahedron = Element::Type::TETRAHEDRON;
constexpr int hexahedron = Element::Type::HEXAHEDRON;
constexpr int wedge = Element::Type::WEDGE;

}  // namespace xfl

inline bool UsesTensorBasis(const FiniteElementSpace *fes)
{
   return mfem::UsesTensorBasis(*fes);
}

int sym(int u) { return u; }
int dot(int u, int v) { return u * v; }

}  // namespace mfem


template <int DIM, int DX0, int DX1>
inline static void KSetup1(const int ndofs, const int vdim, const int NE,
                           const double *__restrict__ J0,
                           const double *__restrict__ w,
                           double *__restrict__ dx)
{
   assert(vdim == 1);
   static constexpr int Q1D = (SOL_P + 2);

   // kernel operations: u,G,*,D,*,v,G,T,*,Ye

   const auto J = Reshape(J0, DIM, DIM, Q1D, Q1D, Q1D, NE);
   const auto W = Reshape(w, Q1D, Q1D, Q1D);
   auto DX = Reshape(dx, DX0, DX1, Q1D, Q1D, Q1D, NE);

   MFEM_FORALL_3D(e, NE, Q1D, Q1D, Q1D,
   {
      MFEM_FOREACH_THREAD(qz, z, Q1D)
      {
         MFEM_FOREACH_THREAD(qy, y, Q1D)
         {
            MFEM_FOREACH_THREAD(qx, x, Q1D)
            {
               const double irw = W(qx, qy, qz);
               const double *Jtr = &J(0, 0, qx, qy, qz, e);
               const double detJ = kernels::Det<DIM>(Jtr);
               const double wd = irw * detJ;
               double Jrt[DIM * DIM];
               kernels::CalcInverse<DIM>(Jtr, Jrt);
               double A[DX0 * DX1];
               double D[DX0 * DX1] = {wd, 0, 0, 0, wd, 0, 0, 0, wd};
               kernels::MultABt(DIM, DIM, DIM, D, Jrt, A);
               kernels::Mult(DIM, DIM, DIM, A, Jrt, &DX(0, 0, qx, qy, qz, e));
            }
         }
      }
      MFEM_SYNC_THREAD;
   });
}

template<int DIM, int DX0, int DX1> inline static
void KMult1(const int ndofs, const int vdim, const int NE,
            const double * __restrict__ B,
            const double * __restrict__ G,
            const int * __restrict__ map,
            const double * __restrict__ dx,
            const double * __restrict__ xd,
            double * __restrict__ yd)
{

   static constexpr int D1D = 4;
   static constexpr int MD1 = 4;
   static constexpr int Q1D = 5;
   static constexpr int MQ1 = 5;
   static constexpr int SMS = SIMD_SIZE;

   // kernel operations: u,G,*,D,*,v,G,T,*,Ye

   assert(vdim == 1);

   const auto b = Reshape(B, Q1D, D1D);
   const auto g = Reshape(G, Q1D, Q1D);
   const auto DX = Reshape(dx, DX0,DX1, Q1D,Q1D,Q1D, NE);
   const auto MAP = Reshape(map, D1D,D1D,D1D, NE);
   const auto XD = Reshape(xd, ndofs);
   auto YD = Reshape(yd, ndofs);
   // kernel operations: u,G,*,D,*,v,G,T,*,Ye

   MFEM_VERIFY((NE % SIMD_SIZE) == 0, "NE vs SIMD_SIZE error!");

   for (int e = 0; e < NE; e+=SIMD_SIZE)
   {
      MFEM_SHARED Real s_Iq[MQ1][MQ1][MQ1];
      MFEM_SHARED double s_D[MQ1][MQ1];
      MFEM_SHARED double s_I[MQ1][MD1];
      MFEM_SHARED Real s_Gqr[MQ1][MQ1];
      MFEM_SHARED Real s_Gqs[MQ1][MQ1];

      Real r_qt[MQ1][MQ1];
      Real r_q[MQ1][MQ1][MQ1];
      Real r_Aq[MQ1][MQ1][MQ1];

      // Load X
      MFEM_FOREACH_THREAD(j,y,Q1D)
      {
         MFEM_FOREACH_THREAD(i,x,Q1D)
         {
            s_D[j][i] = g(i,j);
            if (i<D1D) { s_I[j][i] = b(j,i); }
            if (i<D1D && j<D1D)
            {
#pragma unroll D1D
               for (int k = 0; k < D1D; k++)
               {
                  Real vXD;
                  for (int v = 0; v < SMS; v++)
                  {
                     const int gid = MAP(i, j, k, e + v);
                     const int idx = gid >= 0 ? gid : -1 - gid;
                     vXD[v] = XD(idx);
                  }
                  r_q[j][i][k] = vXD;
               }
            }
         }
      } MFEM_SYNC_THREAD;

      // Grad1X
      MFEM_FOREACH_THREAD(b,y,Q1D)
      {
         MFEM_FOREACH_THREAD(a,x,Q1D)
         {
            if (a<D1D && b<D1D)
            {
#pragma unroll Q1D
               for (int k=0; k<Q1D; ++k)
               {
                  Real res; res = 0.0;
#pragma unroll D1D
                  for (int c=0; c<D1D; ++c)
                  {
                     res += s_I[k][c]*r_q[b][a][c];
                  }
                  s_Iq[k][b][a] = res;
               }
            }
         }
      } MFEM_SYNC_THREAD;

      // Grad1Y
      MFEM_FOREACH_THREAD(k,y,Q1D)
      {
         MFEM_FOREACH_THREAD(a,x,Q1D)
         {
            if (a<D1D)
            {
               for (int b=0; b<D1D; ++b)
               {
                  r_Aq[k][a][b] = s_Iq[k][b][a];
               }
#pragma unroll Q1D
               for (int j=0; j<Q1D; ++j)
               {
                  Real res; res = 0;
#pragma unroll D1D
                  for (int b=0; b<D1D; ++b)
                  {
                     res += s_I[j][b]*r_Aq[k][a][b];
                  }
                  s_Iq[k][j][a] = res;
               }
            }
         }
      } MFEM_SYNC_THREAD;

      // Grad1Z
      MFEM_FOREACH_THREAD(k,y,Q1D)
      {
         MFEM_FOREACH_THREAD(j,x,Q1D)
         {
            for (int a=0; a<D1D; ++a)
            {
               r_Aq[k][j][a] = s_Iq[k][j][a];
            }
#pragma unroll Q1D
            for (int i=0; i<Q1D; ++i)
            {
               Real res; res = 0;
#pragma unroll D1D
               for (int a=0; a<D1D; ++a)
               {
                  res += s_I[i][a]*r_Aq[k][j][a];
               }
               s_Iq[k][j][i] = res;
            }
         }
      } MFEM_SYNC_THREAD;

      // Flush
      MFEM_FOREACH_THREAD(j,y,Q1D)
      {
         MFEM_FOREACH_THREAD(i,x,Q1D)
         {
#pragma unroll Q1D
            for (int k = 0; k < Q1D; k++) { r_Aq[j][i][k] = 0.0; }
         }
      } MFEM_SYNC_THREAD;

      // Q-Function
#pragma unroll Q1D
      for (int k = 0; k < Q1D; k++)
      {
         MFEM_SYNC_THREAD;
         MFEM_FOREACH_THREAD(j,y,Q1D)
         {
            MFEM_FOREACH_THREAD(i,x,Q1D)
            {
               Real qr, qs; qr = 0.0; qs = 0.0;
               r_qt[j][i] = 0.0;
#pragma unroll Q1D
               for (int m = 0; m < Q1D; m++)
               {
                  const double Dim = s_D[i][m];
                  const double Djm = s_D[j][m];
                  const double Dkm = s_D[k][m];
                  qr += Dim*s_Iq[k][j][m];
                  qs += Djm*s_Iq[k][m][i];
                  r_qt[j][i] += Dkm*s_Iq[m][j][i];
               }
               const Real qt = r_qt[j][i];
               const Real u[DX0] = {qr, qs, qt}; Real v[DX0];
               kernels::Mult(DX0,DX1,&DX(0,0,i,j,k,e),u,v);
               s_Gqr[j][i] = v[0];
               s_Gqs[j][i] = v[1];
               r_qt[j][i]  = v[2];
            }
         }
         MFEM_SYNC_THREAD;
         MFEM_FOREACH_THREAD(j,y,Q1D)
         {
            MFEM_FOREACH_THREAD(i,x,Q1D)
            {
               Real Aqtmp; Aqtmp = 0.0;
#pragma unroll Q1D
               for (int m = 0; m < Q1D; m++)
               {
                  const double Dmi = s_D[m][i];
                  const double Dmj = s_D[m][j];
                  const double Dkm = s_D[k][m];
                  Aqtmp += Dmi*s_Gqr[j][m];
                  Aqtmp += Dmj*s_Gqs[m][i];
                  r_Aq[j][i][m] += Dkm*r_qt[j][i];
               }
               r_Aq[j][i][k] += Aqtmp;
            }
         } MFEM_SYNC_THREAD;
      }
      // GradZT
      MFEM_FOREACH_THREAD(j,y,Q1D)
      {
         MFEM_FOREACH_THREAD(i,x,Q1D)
         {
#pragma unroll D1D
            for (int c=0; c<D1D; ++c)
            {
               Real res; res = 0;
#pragma unroll Q1D
               for (int k=0; k<Q1D; ++k)
               {
                  res += s_I[k][c]*r_Aq[j][i][k];
               }
               s_Iq[c][j][i] = res;
            }
         }
      } MFEM_SYNC_THREAD;
      // GradYT
      MFEM_FOREACH_THREAD(c,y,Q1D)
      {
         MFEM_FOREACH_THREAD(i,x,Q1D)
         {
            if (c<D1D)
            {
#pragma unroll Q1D
               for (int j=0; j<Q1D; ++j)
               {
                  r_Aq[c][i][j] = s_Iq[c][j][i];
               }
#pragma unroll D1D
               for (int b=0; b<D1D; ++b)
               {
                  Real res; res = 0;
#pragma unroll Q1D
                  for (int j=0; j<Q1D; ++j)
                  {
                     res += s_I[j][b]*r_Aq[c][i][j];
                  }
                  s_Iq[c][b][i] = res;
               }
            }
         }
      } MFEM_SYNC_THREAD;
      // GradXT
      MFEM_FOREACH_THREAD(c,y,Q1D)
      {
         MFEM_FOREACH_THREAD(b,x,Q1D)
         {
            if (b<D1D && c<D1D)
            {
#pragma unroll Q1D
               for (int i=0; i<Q1D; ++i)
               {
                  r_Aq[c][b][i] = s_Iq[c][b][i];
               }
#pragma unroll D1D
               for (int a=0; a<D1D; ++a)
               {
                  Real res; res = 0;
#pragma unroll Q1D
                  for (int i=0; i<Q1D; ++i)
                  {
                     res += s_I[i][a]*r_Aq[c][b][i];
                  }
                  s_Iq[c][b][a] = res;
               }
            }
         }
      } MFEM_SYNC_THREAD;
      // Scatter
      MFEM_FOREACH_THREAD(j,y,Q1D)
      {
         MFEM_FOREACH_THREAD(i,x,Q1D)
         {
            if (i<D1D && j<D1D)
            {
#pragma unroll D1D
               for (int k = 0; k < D1D; k++)
               {
                  for (int v = 0; v < SMS; v++)
                  {
                     const int gid = MAP(i, j, k, e + v);
                     const int idx = gid >= 0 ? gid : -1 - gid;
                     AtomicAdd(YD(idx), (s_Iq[k][j][i])[v]);
                  }
               }
            }
         }
      }
   } // MFEM_FORALL_2D
} // KMult1

int main(int argc, char *argv[])
{
   int status = 0;
   int num_procs, myid;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);

   assert(MESH_P == 1);
   assert(IR_TYPE == 0);
   assert(IR_ORDER == 0);
   assert(PROBLEM == 0);
   assert(GEOM == Geometry::CUBE);

   /*
      alignas(32) double Z[4];
      printf("sizeof(Z):%d, "
             "alignof(Z):%d, "
             "MFEM_ALIGN_BYTES:%d, "
             "MFEM_ALIGN_SIZE(4,double):%d"
             "\n",
             sizeof(Z),
             alignof(Z),
             MFEM_ALIGN_BYTES,
             MFEM_ALIGN_SIZE(4,double));*/

   const char *mesh_file = "../../data/hex-01x01x01.mesh";
   int ser_ref_levels = 3;
   int par_ref_levels = 1;
   Array<int> nxyz;
   int order = SOL_P;
   const char *basis_type = "G";
   bool static_cond = false;
   const char *pc = "none";
   bool perf = true;
   bool matrix_free = true;
   int max_iter = 50;
   bool visualization = 0;
   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
   args.AddOption(&ser_ref_levels, "-rs", "--refine-serial",
                  "Number of times to refine the mesh uniformly in serial.");
   args.AddOption(&par_ref_levels, "-rp", "--refine-parallel",
                  "Number of times to refine the mesh uniformly in parallel.");
   args.AddOption(&nxyz, "-c", "--cartesian-partitioning",
                  "Use Cartesian partitioning.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree) or -1 for"
                  " isoparametric space.");
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
   args.AddOption(&max_iter, "-mi", "--max-iter",
                  "Maximum number of iterations.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.Parse();
   if (!args.Good())
   {
      if (myid == 0)
      {
         args.PrintUsage(std::cout);
      }
      return 1;
   }
   if (myid == 0)
   {
      args.PrintOptions(std::cout);
   }
   ParMesh *pmesh = nullptr;
   {
      Mesh *mesh = new Mesh(mesh_file, 1, 1);
      int dim = mesh->Dimension();
      {
         int ref_levels = (int)floor(log(10000. / mesh->GetNE()) / log(2.) / dim);
         ref_levels = (ser_ref_levels != -1) ? ser_ref_levels : ref_levels;
         for (int l = 0; l < ref_levels; l++)
         {
            if (myid == 0)
            {
               std::cout << "Serial refinement: level " << l << " -> level " << l + 1
                         << " ..." << std::flush;
            }
            mesh->UniformRefinement();
            MPI_Barrier(MPI_COMM_WORLD);
            if (myid == 0)
            {
               std::cout << " done." << std::endl;
            }
         }
      }
      MFEM_VERIFY(nxyz.Size() == 0 || nxyz.Size() == mesh->SpaceDimension(),
                  "Expected " << mesh->SpaceDimension()
                  << " integers with the "
                  "option --cartesian-partitioning.");
      int *partitioning = nxyz.Size() ? mesh->CartesianPartitioning(nxyz) : NULL;
      pmesh = new ParMesh(MPI_COMM_WORLD, *mesh, partitioning);
      delete[] partitioning;
      delete mesh;
      {
         for (int l = 0; l < par_ref_levels; l++)
         {
            if (myid == 0)
            {
               std::cout << "Parallel refinement: level " << l << " -> level "
                         << l + 1 << " ..." << std::flush;
            }
            pmesh->UniformRefinement();
            MPI_Barrier(MPI_COMM_WORLD);
            if (myid == 0)
            {
               std::cout << " done." << std::endl;
            }
         }
      }
      pmesh->PrintInfo(std::cout);
   }

   static constexpr int Q1D = (SOL_P + 2);
   if (myid == 0)
   {
      std::cout << "XFL(SIMD_" << SIMD_SIZE << ") "
                << "version using integration rule with " << (Q1D * Q1D * Q1D)
                << " points ...\n";
   }

   const int p = SOL_P;
   const int dim = 3;
   auto &mesh = xfl::Mesh(pmesh);
   const int el = (dim == 2) ? xfl::quadrilateral : xfl::hexahedron;
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
      struct QMult1 : public xfl::Operator<3>
      {
         QMult1(const ParFiniteElementSpace *fes) : xfl::Operator<3>(fes)
         {
            dbg();
            Setup();
         }
         ~QMult1() { dbg(); }
         void Setup()
         {
            dx.SetSize(NQ * NE * 3 * 3,
                       Device::GetDeviceMemoryType());  // DX shape: 3x3
            KSetup1<3, 3, 3>(NDOFS, VDIM, NE, J0.Read(), ir.GetWeights().Read(),
                             dx.Write());
         }
         void Mult(const mfem::Vector &x, mfem::Vector &y) const
         {
            y = 0.0;
            KMult1<3, 3, 3>(NDOFS, VDIM, NE, maps->B.Read(), maps->G.Read(),
                            ER.GatherMap().Read(), dx.Read(), x.Read(),
                            y.ReadWrite());
         }
      };  // QMult struct
      QMult1 *QM1 = new QMult1(fes);
      xfl::QForm QForm1(fes1, qs1, QM1);
      return QForm1;
   }();
   status |= xfl::benchmark(a == b, x, bc, 0, 50, -1);
   MPI_Finalize();
   return status;
}
