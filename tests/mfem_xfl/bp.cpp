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

#include "mfem.hpp"
#define dbg(...)

#include "general/forall.hpp"
#include "linalg/kernels.hpp"

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

using namespace mfem;

namespace mfem {

using FE = mfem::FiniteElement;
using QI = mfem::QuadratureInterpolator;

// Kernels addons //////////////////////////////////////////////////////////////
namespace kernels {

/// Load B1d & G1d matrices into shared memory
template <int MD1, int MQ1>
MFEM_HOST_DEVICE inline void LoadBG(const int D1D, const int Q1D,
                                    const ConstDeviceMatrix b,
                                    const ConstDeviceMatrix g,
                                    double sBG[2][MQ1 * MD1]) {
  const int tidz = MFEM_THREAD_ID(z);
  DeviceMatrix B(sBG[0], MD1, MQ1);
  DeviceMatrix G(sBG[1], MD1, MQ1);

  if (tidz == 0) {
    MFEM_FOREACH_THREAD(d, y, D1D) {
      MFEM_FOREACH_THREAD(q, x, Q1D) {
        B(d, q) = b(q, d);
        G(d, q) = g(q, d);
      }
    }
  }
  MFEM_SYNC_THREAD;
}

/// Load Bt1d & Gt1d matrices into shared memory
template <int MD1, int MQ1, typename T>
MFEM_HOST_DEVICE inline void LoadBGt(const int D1D, const int Q1D,
                                     const DeviceTensor<2, const double> b,
                                     const DeviceTensor<2, const double> g,
                                     T sBG[2][MQ1 * MD1]) {
  const int tidz = MFEM_THREAD_ID(z);
  DeviceTensor<2, T> Bt(sBG[0], MQ1, MD1);
  DeviceTensor<2, T> Gt(sBG[1], MQ1, MD1);

  if (tidz == 0) {
    MFEM_FOREACH_THREAD(d, y, D1D) {
      MFEM_FOREACH_THREAD(q, x, Q1D) {
        Bt(q, d) = b(q, d);
        Gt(q, d) = g(q, d);
      }
    }
  }
  MFEM_SYNC_THREAD;
}

template <int MD1, typename T = double, int SMS = 1>
MFEM_HOST_DEVICE inline void LoadXDGather(
    const int e, const int D1D, const DeviceTensor<4, const int> MAP,
    const DeviceTensor<1, const double> xd, T sm[MD1 * MD1 * MD1]) {
  DeviceTensor<3, T> X(sm, MD1, MD1, MD1);

  MFEM_FOREACH_THREAD(dz, z, D1D) {
    MFEM_FOREACH_THREAD(dy, y, D1D) {
      MFEM_FOREACH_THREAD(dx, x, D1D) {
        // Gather
        T XD;
        for (int i = 0; i < SMS; i++) {
          const int gid = MAP(dx, dy, dz, e + i);
          const int j = gid >= 0 ? gid : -1 - gid;
          XD[i] = xd(j);
        }
        X(dx, dy, dz) = XD;
      }
    }
  }
  MFEM_SYNC_THREAD;
}

///////////////////////////////////////////////////////////////////////////////
/// Push 3D Scalar Evaluation

/// 3D Scalar Gradient, 1/3
template <int MD1, int MQ1, typename T>
MFEM_HOST_DEVICE inline void Grad1X(const int D1D, const int Q1D,
                                    const double (*__restrict__ BG)[MQ1 * MD1],
                                    const T(*__restrict__ DDD),
                                    T (*__restrict__ DDQ)[MD1 * MD1 * MQ1]) {
  DeviceTensor<2, const double> B(BG[0], MD1, MQ1);
  DeviceTensor<2, const double> G(BG[1], MD1, MQ1);
  DeviceTensor<3, const T> X(DDD, MD1, MD1, MD1);
  DeviceTensor<3, T> XB(DDQ[0], MQ1, MD1, MD1);
  DeviceTensor<3, T> XG(DDQ[1], MQ1, MD1, MD1);

  MFEM_FOREACH_THREAD(dz, z, D1D) {
    MFEM_FOREACH_THREAD(dy, y, D1D) {
      MFEM_FOREACH_THREAD(qx, x, Q1D) {
        T u;
        u = 0.0;
        T v;
        v = 0.0;
        for (int dx = 0; dx < D1D; ++dx) {
          const T xx = X(dx, dy, dz);
          const double Bx = B(dx, qx);
          const double Gx = G(dx, qx);
          u += Bx * xx;
          v += Gx * xx;
        }
        XB(qx, dy, dz) = u;
        XG(qx, dy, dz) = v;
      }
    }
  }
  MFEM_SYNC_THREAD;
}

/// 3D Scalar Gradient, 2/3
template <int MD1, int MQ1, typename T>
MFEM_HOST_DEVICE inline void Grad1Y(
    const int D1D, const int Q1D, const double (*__restrict__ BG)[MQ1 * MD1],
    const T (*__restrict__ DDQ)[MD1 * MD1 * MQ1],
    T (*__restrict__ DQQ)[MD1 * MQ1 * MQ1]) {
  DeviceTensor<2, const double> B(BG[0], MD1, MQ1);
  DeviceTensor<2, const double> G(BG[1], MD1, MQ1);
  DeviceTensor<3, const T> XB(DDQ[0], MQ1, MD1, MD1);
  DeviceTensor<3, const T> XG(DDQ[1], MQ1, MD1, MD1);
  DeviceTensor<3, T> XBB(DQQ[0], MQ1, MQ1, MD1);
  DeviceTensor<3, T> XBG(DQQ[1], MQ1, MQ1, MD1);
  DeviceTensor<3, T> XGB(DQQ[2], MQ1, MQ1, MD1);

  MFEM_FOREACH_THREAD(dz, z, D1D) {
    MFEM_FOREACH_THREAD(qy, y, Q1D) {
      MFEM_FOREACH_THREAD(qx, x, Q1D) {
        T u;
        u = 0.0;
        T v;
        v = 0.0;
        T w;
        w = 0.0;
        for (int dy = 0; dy < D1D; ++dy) {
          const double By = B(dy, qy);
          const double Gy = G(dy, qy);
          u += XB(qx, dy, dz) * By;
          v += XG(qx, dy, dz) * By;
          w += XB(qx, dy, dz) * Gy;
        }
        XBB(qx, qy, dz) = u;
        XBG(qx, qy, dz) = v;
        XGB(qx, qy, dz) = w;
      }
    }
  }
  MFEM_SYNC_THREAD;
}

/// 3D Scalar Gradient, 3/3
template <int MD1, int MQ1, typename T = double>
MFEM_HOST_DEVICE inline void Grad1Z(
    const int D1D, const int Q1D, const double (*__restrict__ BG)[MQ1 * MD1],
    const T (*__restrict__ DQQ)[MD1 * MQ1 * MQ1],
    T (*__restrict__ QQQ)[MQ1 * MQ1 * MQ1]) {
  DeviceTensor<2, const double> B(BG[0], MD1, MQ1);
  DeviceTensor<2, const double> G(BG[1], MD1, MQ1);
  DeviceTensor<3, const T> XBB(DQQ[0], MQ1, MQ1, MD1);
  DeviceTensor<3, const T> XBG(DQQ[1], MQ1, MQ1, MD1);
  DeviceTensor<3, const T> XGB(DQQ[2], MQ1, MQ1, MD1);
  DeviceTensor<3, T> XBBG(QQQ[0], MQ1, MQ1, MQ1);
  DeviceTensor<3, T> XBGB(QQQ[1], MQ1, MQ1, MQ1);
  DeviceTensor<3, T> XGBB(QQQ[2], MQ1, MQ1, MQ1);

  MFEM_FOREACH_THREAD(qz, z, Q1D) {
    MFEM_FOREACH_THREAD(qy, y, Q1D) {
      MFEM_FOREACH_THREAD(qx, x, Q1D) {
        T u;
        u = 0.0;
        T v;
        v = 0.0;
        T w;
        w = 0.0;
        for (int dz = 0; dz < D1D; ++dz) {
          const double Bz = B(dz, qz);
          const double Gz = G(dz, qz);
          u += XBG(qx, qy, dz) * Bz;
          v += XGB(qx, qy, dz) * Bz;
          w += XBB(qx, qy, dz) * Gz;
        }
        XBBG(qx, qy, qz) = u;
        XBGB(qx, qy, qz) = v;
        XGBB(qx, qy, qz) = w;
      }
    }
  }
  MFEM_SYNC_THREAD;
}

/// Pull 3D Scalar Gradient
template <int MQ1, typename T>
MFEM_HOST_DEVICE inline void PullGrad1(
    const int x, const int y, const int z,
    const T(*__restrict__ QQQ) /*[3]*/[MQ1 * MQ1 * MQ1], T *__restrict__ A) {
  DeviceTensor<3, const T> BBG(QQQ[0], MQ1, MQ1, MQ1);
  DeviceTensor<3, const T> BGB(QQQ[1], MQ1, MQ1, MQ1);
  DeviceTensor<3, const T> GBB(QQQ[2], MQ1, MQ1, MQ1);

  A[0] = BBG(x, y, z);
  A[1] = BGB(x, y, z);
  A[2] = GBB(x, y, z);
}

/// Push 3D Scalar Gradient
template <int MQ1, typename T>
MFEM_HOST_DEVICE inline void PushGrad1(
    const int x, const int y, const int z, const T *__restrict__ A,
    T(*__restrict__ QQQ) /*[3]*/[MQ1 * MQ1 * MQ1]) {
  DeviceTensor<3, T> BBG(QQQ[0], MQ1, MQ1, MQ1);
  DeviceTensor<3, T> BGB(QQQ[1], MQ1, MQ1, MQ1);
  DeviceTensor<3, T> GBB(QQQ[2], MQ1, MQ1, MQ1);

  BBG(x, y, z) = A[0];
  BGB(x, y, z) = A[1];
  GBB(x, y, z) = A[2];
}

/// 3D Transposed Scalar Gradient, 1/3
template <int MD1, int MQ1, typename T>
MFEM_HOST_DEVICE inline void Grad1Zt(
    const int D1D, const int Q1D,
    const double(*__restrict__ BG) /*[2]*/[MQ1 * MD1],
    const T(*__restrict__ QQQ) /*[3]*/[MQ1 * MQ1 * MQ1],
    T(*__restrict__ DQQ) /*[3]*/[MD1 * MQ1 * MQ1]) {
  DeviceTensor<2, const double> Bt(BG[0], MQ1, MD1);
  DeviceTensor<2, const double> Gt(BG[1], MQ1, MD1);
  DeviceTensor<3, const T> XBBG(QQQ[0], MQ1, MQ1, MQ1);
  DeviceTensor<3, const T> XBGB(QQQ[1], MQ1, MQ1, MQ1);
  DeviceTensor<3, const T> XGBB(QQQ[2], MQ1, MQ1, MQ1);
  DeviceTensor<3, T> XBB(DQQ[0], MQ1, MQ1, MD1);
  DeviceTensor<3, T> XBG(DQQ[1], MQ1, MQ1, MD1);
  DeviceTensor<3, T> XGB(DQQ[2], MQ1, MQ1, MD1);

  MFEM_FOREACH_THREAD(qz, z, Q1D) {
    MFEM_FOREACH_THREAD(qy, y, Q1D) {
      MFEM_FOREACH_THREAD(dx, x, D1D) {
        T u;
        u = 0.0;
        T v;
        v = 0.0;
        T w;
        w = 0.0;
        for (int qx = 0; qx < Q1D; ++qx) {
          const double Btx = Bt(qx, dx);
          const double Gtx = Gt(qx, dx);
          u += XBBG(qx, qy, qz) * Gtx;
          v += XBGB(qx, qy, qz) * Btx;
          w += XGBB(qx, qy, qz) * Btx;
        }
        XBB(qz, qy, dx) = u;
        XBG(qz, qy, dx) = v;
        XGB(qz, qy, dx) = w;
      }
    }
  }
  MFEM_SYNC_THREAD;
}

/// 3D Transposed Scalar Gradient, 2/3
template <int MD1, int MQ1, typename T>
MFEM_HOST_DEVICE inline void Grad1Yt(
    const int D1D, const int Q1D,
    const double(*__restrict__ BG) /*[2]*/[MQ1 * MD1],
    const T(*__restrict__ DQQ) /*[3]*/[MD1 * MQ1 * MQ1],
    T DDQ[3][MD1 * MD1 * MQ1]) {
  DeviceTensor<2, const double> Bt(BG[0], MQ1, MD1);
  DeviceTensor<2, const double> Gt(BG[1], MQ1, MD1);
  DeviceTensor<3, const T> XBB(DQQ[0], MQ1, MQ1, MD1);
  DeviceTensor<3, const T> XBG(DQQ[1], MQ1, MQ1, MD1);
  DeviceTensor<3, const T> XGB(DQQ[2], MQ1, MQ1, MD1);
  DeviceTensor<3, T> XB(DDQ[0], MQ1, MD1, MD1);
  DeviceTensor<3, T> XG(DDQ[1], MQ1, MD1, MD1);
  DeviceTensor<3, T> XC(DDQ[2], MQ1, MD1, MD1);

  MFEM_FOREACH_THREAD(qz, z, Q1D) {
    MFEM_FOREACH_THREAD(dy, y, D1D) {
      MFEM_FOREACH_THREAD(dx, x, D1D) {
        T u;
        u = 0.0;
        T v;
        v = 0.0;
        T w;
        w = 0.0;
        for (int qy = 0; qy < Q1D; ++qy) {
          const double Bty = Bt(qy, dy);
          const double Gty = Gt(qy, dy);
          u += XBB(qz, qy, dx) * Bty;
          v += XBG(qz, qy, dx) * Gty;
          w += XGB(qz, qy, dx) * Bty;
        }
        XB(qz, dy, dx) = u;
        XG(qz, dy, dx) = v;
        XC(qz, dy, dx) = w;
      }
    }
  }
  MFEM_SYNC_THREAD;
}

/// 3D Transposed Gradient, 3/3
template <int MD1, int MQ1, typename T, int SMS>
MFEM_HOST_DEVICE inline void Grad1XtDScatter(
    const int D1D, const int Q1D,
    const double(*__restrict__ BG) /*[2]*/[MQ1 * MD1],
    const T(*__restrict__ DDQ) /*[3]*/[MD1 * MD1 * MQ1],
    const DeviceTensor<4, const int> MAP, DeviceTensor<1, double> YD,
    const int e) {
  DeviceTensor<2, const double> Bt(BG[0], MQ1, MD1);
  DeviceTensor<2, const double> Gt(BG[1], MQ1, MD1);
  DeviceTensor<3, const T> XB(DDQ[0], MQ1, MD1, MD1);
  DeviceTensor<3, const T> XG(DDQ[1], MQ1, MD1, MD1);
  DeviceTensor<3, const T> XC(DDQ[2], MQ1, MD1, MD1);

  MFEM_FOREACH_THREAD(dz, z, D1D) {
    MFEM_FOREACH_THREAD(dy, y, D1D) {
      MFEM_FOREACH_THREAD(dx, x, D1D) {
        T u;
        u = 0.0;
        T v;
        v = 0.0;
        T w;
        w = 0.0;
        for (int qz = 0; qz < Q1D; ++qz) {
          const double Btz = Bt(qz, dz);
          const double Gtz = Gt(qz, dz);
          u += XB(qz, dy, dx) * Btz;
          v += XG(qz, dy, dx) * Btz;
          w += XC(qz, dy, dx) * Gtz;
        }
        const T value = u + v + w;
        for (int i = 0; i < SMS; i++) {
          const int gid = MAP(dx, dy, dz, e + i);
          const int j = gid >= 0 ? gid : -1 - gid;
          AtomicAdd(YD(j), value[i]);
        }
      }
    }
  }
}

}  // namespace kernels

// XFL addons //////////////////////////////////////////////////////////////////
namespace xfl {

class XElementRestriction : public ElementRestriction {
  const ParFiniteElementSpace &fes;

 public:
  XElementRestriction(const ParFiniteElementSpace *fes, ElementDofOrdering edo)
      : ElementRestriction(*fes, edo), fes(*fes) {}

  const Array<int> &ScatterMap() const { return offsets; }
  const Array<int> &ScatterIdx() const { return indices; }
  const Array<int> &GatherMap() const { return gatherMap; }

  void Mult(const Vector &x, Vector &y) const {
    const int ndof = fes.GetFE(0)->GetDof();
    const int ndofs = fes.GetNDofs();

    const auto d_x = Reshape(x.Read(), ndofs);
    const auto d_j = gatherMap.Read();

    auto d_y = Reshape(y.Write(), ndof, ne);

    MFEM_FORALL(i, ndof * ne, { d_y(i % ndof, i / ndof) = d_x(d_j[i]); });
  }

  void MultTranspose(const Vector &x, Vector &y) const {
    const int nd = fes.GetFE(0)->GetDof();

    const auto d_offsets = offsets.Read();
    const auto d_indices = indices.Read();

    const auto d_x = Reshape(x.Read(), nd, ne);
    auto d_y = Reshape(y.Write(), ndofs);
    MFEM_FORALL(i, ndofs, {
      double dofValue = 0.0;
      const int offset = d_offsets[i];
      const int nextOffset = d_offsets[i + 1];
      for (int j = offset; j < nextOffset; ++j) {
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
class Operator<3> : public mfem::Operator {
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
        R(pfes->GetRestrictionMatrix()) {
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
struct Problem {
  mfem::Operator *QM{nullptr};
  mfem::ParLinearForm &b;
  Problem(mfem::ParLinearForm &b, mfem::Operator *QM) : QM(QM), b(b) {}
  ~Problem() {
    dbg();
    delete QM;
  }
};

/** ****************************************************************************
 * @brief The QForm class
 ******************************************************************************/
class QForm {
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
  Problem *operator==(QForm &rhs) {
    assert(!b);
    mfem::ParLinearForm *b = new mfem::ParLinearForm(rhs.ParFESpace());
    assert(b);
    if (!rhs.ConstantCoeff() && !rhs.FunctionCoeff()) {
      ConstantCoefficient *cst = new ConstantCoefficient(1.0);
      b->AddDomainIntegrator(new DomainLFIntegrator(*cst));
    } else if (rhs.ConstantCoeff()) {
      ConstantCoefficient *cst = rhs.ConstantCoeff();
      b->AddDomainIntegrator(new DomainLFIntegrator(*cst));
    } else if (rhs.FunctionCoeff()) {
      FunctionCoefficient *func = rhs.FunctionCoeff();
      b->AddDomainIntegrator(new DomainLFIntegrator(*func));
    } else {
      assert(false);
    }

    return new Problem(*b, QM);
  }

  // + operator on QForms
  QForm &operator+(QForm &rhs) {
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
class Function : public ParGridFunction {
 public:
  Function(ParFiniteElementSpace *pfes) : ParGridFunction(pfes) {
    assert(pfes);
    assert(pfes->GlobalTrueVSize() > 0);
  }
  void operator=(double value) { ParGridFunction::operator=(value); }
  int geometric_dimension() { return fes->GetMesh()->SpaceDimension(); }
  ParFiniteElementSpace *ParFESpace() { return ParGridFunction::ParFESpace(); }
  const ParFiniteElementSpace *ParFESpace() const {
    return ParGridFunction::ParFESpace();
  }
  ConstantCoefficient *ConstantCoeff() const { return nullptr; }
  FunctionCoefficient *FunctionCoeff() const { return nullptr; }
};

/** ****************************************************************************
 * @brief TrialFunction class
 ******************************************************************************/
class TrialFunction : public Function {
 public:
  TrialFunction(ParFiniteElementSpace *pfes) : Function(pfes) {}
  ~TrialFunction() {}
};

/** ****************************************************************************
 * @brief TestFunction class
 ******************************************************************************/
class TestFunction : public Function {
 public:
  TestFunction(ParFiniteElementSpace *pfes) : Function(pfes) {}
  TestFunction(const TestFunction &) = default;
  ~TestFunction() {}
};

/** ****************************************************************************
 * @brief Constant class
 ******************************************************************************/
class Constant {
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
class Expression {
  FunctionCoefficient *fct = nullptr;

 public:
  Expression(std::function<double(const Vector &)> F)
      : fct(new FunctionCoefficient(F)) {}
  ~Expression() { /*delete fct;*/
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
FiniteElementCollection *FiniteElement(std::string family, int type, int p) {
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
                                     int p) {
  assert(false);
  const int dim = pmesh->Dimension();
  MFEM_VERIFY(family == "P", "Unsupported FE!");
  FiniteElementCollection *fec = new H1_FECollection(p, dim);
  return new ParFiniteElementSpace(pmesh, fec);
}

ParFiniteElementSpace *FunctionSpace(mfem::ParMesh &pmesh, std::string f,
                                     int p) {
  assert(false);
  return FunctionSpace(&pmesh, f, p);
}

ParFiniteElementSpace *FunctionSpace(mfem::ParMesh &pmesh,
                                     FiniteElementCollection *fec) {
  const int vdim = 1;
  const Ordering::Type ordering = Ordering::byNODES;
  ParFiniteElementSpace *pfes =
      new ParFiniteElementSpace(&pmesh, fec, vdim, ordering);
  return pfes;
}

ParFiniteElementSpace *FunctionSpace(mfem::ParMesh &pmesh,
                                     FiniteElementCollection *fec,
                                     const int vdim) {
  assert(false);
  return new ParFiniteElementSpace(&pmesh, fec, vdim);
}

/** ****************************************************************************
 * @brief Vector Function Space
 ******************************************************************************/
ParFiniteElementSpace *VectorFunctionSpace(mfem::ParMesh *pmesh,
                                           std::string family, const int p) {
  const int dim = pmesh->Dimension();
  MFEM_VERIFY(family == "P", "Unsupported FE!");
  FiniteElementCollection *fec = new H1_FECollection(p, dim);
  return new ParFiniteElementSpace(pmesh, fec, dim);
}

ParFiniteElementSpace *VectorFunctionSpace(mfem::ParMesh &pmesh,
                                           std::string family, const int p) {
  return VectorFunctionSpace(&pmesh, family, p);
}

ParFiniteElementSpace *VectorFunctionSpace(mfem::ParMesh &pmesh,
                                           FiniteElementCollection *fec) {
  return new ParFiniteElementSpace(&pmesh, fec, pmesh.Dimension());
}

/** ****************************************************************************
 * @brief Boundary Conditions
 ******************************************************************************/
Array<int> DirichletBC(mfem::ParFiniteElementSpace *pfes) {
  assert(pfes);
  Array<int> ess_tdof_list;
  mfem::ParMesh *pmesh = pfes->GetParMesh();
  if (pmesh->bdr_attributes.Size()) {
    Array<int> ess_bdr(pmesh->bdr_attributes.Max());
    ess_bdr = 1;
    pfes->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
  }
  return ess_tdof_list;
}

/** ****************************************************************************
 * @brief solve with boundary conditions
 ******************************************************************************/
int solve(xfl::Problem *pb, xfl::Function &x, Array<int> ess_tdof_list) {
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
int solve(xfl::Problem *pb, xfl::Function &x) {
  Array<int> empty_tdof_list;
  return solve(pb, x, empty_tdof_list);
}

/** ****************************************************************************
 * @brief benchmark this prblem with boundary conditions
 ******************************************************************************/
int benchmark(xfl::Problem *pb, xfl::Function &x, Array<int> ess_tdof_list,
              const double rtol, const int max_it, const int print_lvl) {
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

  if (myid == 0) {
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

int benchmark(xfl::Problem *pb, xfl::Function &x, Array<int> ess_tdof_list) {
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

inline bool UsesTensorBasis(const FiniteElementSpace *fes) {
  return mfem::UsesTensorBasis(*fes);
}

int sym(int u) { return u; }
int dot(int u, int v) { return u * v; }

}  // namespace mfem

#include <thread>

#include "linalg/simd.hpp"
using Real = AutoSIMDTraits<double, double>::vreal_t;
#define SIMD_SIZE (MFEM_SIMD_BYTES / sizeof(double))

template <int DIM, int DX0, int DX1>
inline static void KSetup1(const int ndofs, const int vdim, const int NE,
                           const double *__restrict__ J0,
                           const double *__restrict__ w,
                           double *__restrict__ dx) {
  assert(vdim == 1);
  static constexpr int Q1D = (SOL_P + 2);

  // kernel operations: u,G,*,D,*,v,G,T,*,Ye

  const auto J = Reshape(J0, DIM, DIM, Q1D, Q1D, Q1D, NE);
  const auto W = Reshape(w, Q1D, Q1D, Q1D);
  auto DX = Reshape(dx, DX0, DX1, Q1D, Q1D, Q1D, NE);

  MFEM_FORALL_3D(e, NE, Q1D, Q1D, Q1D, {
    MFEM_FOREACH_THREAD(qz, z, Q1D) {
      MFEM_FOREACH_THREAD(qy, y, Q1D) {
        MFEM_FOREACH_THREAD(qx, x, Q1D) {
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

template <int DIM, int DX0, int DX1>
inline static void KMult1(const int ndofs, const int vdim, const int NE,
                          const double *__restrict__ B,
                          const double *__restrict__ G,
                          const int *__restrict__ map,
                          const double *__restrict__ dx,
                          const double *__restrict__ xd,
                          double *__restrict__ yd) {
  static constexpr int D1D = (SOL_P + 1);
  static constexpr int MD1 = D1D;
  static constexpr int Q1D = (SOL_P + 2);
  static constexpr int MQ1 = Q1D;
  static constexpr int SMS = SIMD_SIZE;

  // kernel operations: u,G,*,D,*,v,G,T,*,Ye

  assert(vdim == 1);

  const auto b = Reshape(B, Q1D, D1D);
  const auto g = Reshape(G, Q1D, D1D);
  const auto DX = Reshape(dx, DX0, DX1, Q1D, Q1D, Q1D, NE);
  const auto MAP = Reshape(map, D1D, D1D, D1D, NE);
  const auto XD = Reshape(xd, ndofs);
  auto YD = Reshape(yd, ndofs);

  double BG[2][MQ1 * MD1];
  kernels::LoadBG<MD1, MQ1>(D1D, Q1D, b, g, BG);

  double BGt[2][MQ1 * MD1];
  kernels::LoadBGt<MD1, MQ1>(D1D, Q1D, b, g, BGt);

  //if ((NE % SIMD_SIZE) != 0) { return; }
  MFEM_VERIFY((NE % SIMD_SIZE) == 0, "NE vs SIMD_SIZE error!")

#ifndef MFEM_USE_THREADS
  //#pragma omp parallel for
  int BATCH_SIZE = 32;
  while ((NE % BATCH_SIZE) != 0) {
    BATCH_SIZE >>= 1;
  }
  while (((NE / BATCH_SIZE) % SIMD_SIZE) != 0) {
    BATCH_SIZE >>= 1;
  }
  MFEM_VERIFY((NE % BATCH_SIZE) == 0, "NE vs BATCH_SIZE error!")
  MFEM_VERIFY(((NE / BATCH_SIZE) % SIMD_SIZE) == 0,
              "NE/BATCH_SIZE vs SIMD_SIZE error!")
  for (size_t eb = 0; eb < (NE / (BATCH_SIZE * SIMD_SIZE)); eb += 1) {
    for (size_t e = eb * BATCH_SIZE * SIMD_SIZE;
         e < (eb + 1) * BATCH_SIZE * SIMD_SIZE; e += SIMD_SIZE) {
#else
  std::vector<std::thread> threads;
  static const unsigned int num_threads = std::thread::hardware_concurrency();
  dbg("NE:%d, num_threads:%d", NE, num_threads);
  MFEM_VERIFY((NE % num_threads) == 0, "NE vs #Threads error")
  int e0 = 0;
  const int NEB = NE / num_threads;
  int NE0 = NEB;
  for (unsigned int tid = 0; tid < num_threads; ++tid) {
      threads.push_back(std::thread(
                           [&](const int tid, const int e0, const int NE0)
      {
      // printf("[#%d] e0:%d, NE0:%d", tid, e0, NE0);
      for (size_t e = e0; e < NE0; e += SIMD_SIZE) {
#endif
      Real DDD[MD1 * MD1 * MD1];
      kernels::LoadXDGather<MD1, Real, SMS>(e, D1D, MAP, XD, DDD);
      // kernel operations: u,G,*,D,*,v,G,T,*,Ye
      // [push] trial u:2
      Real sm0[3][MQ1 * MQ1 * MQ1];
      Real sm1[3][MQ1 * MQ1 * MQ1];
      Real(*DDQ)[MD1 * MD1 * MQ1] = (Real(*)[MD1 * MD1 * MQ1])(sm0);
      Real(*DQQ)[MD1 * MQ1 * MQ1] = (Real(*)[MD1 * MQ1 * MQ1])(sm1);
      Real(*QQQ)[MQ1 * MQ1 * MQ1] = (Real(*)[MQ1 * MQ1 * MQ1])(sm0);
      // Grad(u)
      kernels::Grad1X<MD1, MQ1>(D1D, Q1D, BG, DDD, DDQ);
      kernels::Grad1Y<MD1, MQ1>(D1D, Q1D, BG, DDQ, DQQ);
      kernels::Grad1Z<MD1, MQ1>(D1D, Q1D, BG, DQQ, QQQ);
      // [ pop] u
      for (int qz = 0; qz < Q1D; qz++) {
        for (int qy = 0; qy < Q1D; qy++) {
          for (int qx = 0; qx < Q1D; qx++) {
            Real u[DX0], v[DX0];
            kernels::PullGrad1<MQ1>(qx, qy, qz, QQQ, u);
            kernels::Mult(DX0, DX1, &DX(0, 0, qx, qy, qz, e), u, v);
            kernels::PushGrad1<MQ1>(qx, qy, qz, v, QQQ);
            // [push] test v:2
          }
        }
      }
      // Grad(v)
      kernels::Grad1Zt<MD1, MQ1>(D1D, Q1D, BGt, QQQ, DQQ);
      kernels::Grad1Yt<MD1, MQ1>(D1D, Q1D, BGt, DQQ, DDQ);
      kernels::Grad1XtDScatter<MD1, MQ1, Real, SMS>(D1D, Q1D, BGt, DDQ, MAP, YD,
                                                    e);
      // [ pop] v
    }
  }  // Element for loop
#ifdef MFEM_USE_THREADS
}, tid, e0, NE0)); // lambda & thread vector push_back
e0 += NEB;
NE0 += NEB;
}  // Thread for loop
for (auto &thr : threads) {
  thr.join();
}
#endif
}  // KMult1

int main(int argc, char *argv[]) {
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
  if (!args.Good()) {
    if (myid == 0) {
      args.PrintUsage(std::cout);
    }
    return 1;
  }
  if (myid == 0) {
    args.PrintOptions(std::cout);
  }
  ParMesh *pmesh = nullptr;
  {
    Mesh *mesh = new Mesh(mesh_file, 1, 1);
    int dim = mesh->Dimension();
    {
      int ref_levels = (int)floor(log(10000. / mesh->GetNE()) / log(2.) / dim);
      ref_levels = (ser_ref_levels != -1) ? ser_ref_levels : ref_levels;
      for (int l = 0; l < ref_levels; l++) {
        if (myid == 0) {
          std::cout << "Serial refinement: level " << l << " -> level " << l + 1
                    << " ..." << std::flush;
        }
        mesh->UniformRefinement();
        MPI_Barrier(MPI_COMM_WORLD);
        if (myid == 0) {
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
      for (int l = 0; l < par_ref_levels; l++) {
        if (myid == 0) {
          std::cout << "Parallel refinement: level " << l << " -> level "
                    << l + 1 << " ..." << std::flush;
        }
        pmesh->UniformRefinement();
        MPI_Barrier(MPI_COMM_WORLD);
        if (myid == 0) {
          std::cout << " done." << std::endl;
        }
      }
    }
    pmesh->PrintInfo(std::cout);
  }

  static constexpr int Q1D = (SOL_P + 2);
  if (myid == 0) {
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
  auto b = [&]() {
    constexpr const char *qs0 = "v";
    // var:[v], ops:[Xe,v,]
    // Test FES: 'fes':v (Eval)
    ParFiniteElementSpace *fes0 = fes;
    mfem::Operator *QM0 = nullptr;
    xfl::QForm QForm0(fes0, qs0, QM0);
    return QForm0;
  }();
  auto a = [&]() {
    constexpr const char *qs1 = "dot(grad(u), grad(v))";
    // var:[u,v], ops:[Xe,u,G,*,D,*,v,G,T,*,Ye]
    // Trial FES: 'fes':u (Grad)
    // Test FES: 'fes':v (Grad)
    ParFiniteElementSpace *fes1 = fes;
    struct QMult1 : public xfl::Operator<3> {
      QMult1(const ParFiniteElementSpace *fes) : xfl::Operator<3>(fes) {
        dbg();
        Setup();
      }
      ~QMult1() { dbg(); }
      void Setup() {
        dx.SetSize(NQ * NE * 3 * 3,
                   Device::GetDeviceMemoryType());  // DX shape: 3x3
        KSetup1<3, 3, 3>(NDOFS, VDIM, NE, J0.Read(), ir.GetWeights().Read(),
                         dx.Write());
      }
      void Mult(const mfem::Vector &x, mfem::Vector &y) const {
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
