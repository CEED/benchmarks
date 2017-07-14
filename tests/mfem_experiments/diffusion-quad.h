/*
  Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
  the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
  reserved. See files LICENSE and NOTICE for details.

  This file is part of CEED, a collection of benchmarks, miniapps, software
  libraries and APIs for efficient high-order finite element and spectral
  element discretizations for exascale applications. For more information and
  source code availability see http://github.com/ceed.

  The CEED research is supported by the Exascale Computing Project (17-SC-20-SC)
  a collaborative effort of two U.S. Department of Energy organizations (Office
  of Science and the National Nuclear Security Administration) responsible for
  the planning and preparation of a capable exascale ecosystem, including
  software, applications, hardware, advanced system engineering and early
  testbed platforms, in support of the nation's exascale computing imperative.
*/

void add_mult_diffusion_quad(
   int ndof_1d,      /* number of 1D dofs (points) */
   int nqpt_1d,      /* number of 1D quadrature points */
   int nelem,        /* number of elements */
   double *D,        /* nqpt_1d x nqpt_1d x 3 x nelem; (3) -> (xx,xy,yy) */
   double *B1d,      /* nqpt_1d x ndof_1d dense matrix, column-major layout */
   double *B1d_t,    /* trasnspose of B1d */
   double *G1d,      /* nqpt_1d x ndof_1d dense matrix, column-major layout */
   double *G1d_t,    /* trasnspose of B1d */
   int *dof_offsets, /* array of size ndofs_1d x ndofs_1d x nelem representing a
                        boolean P */
   double *x,        /* input vector */
   double *y         /* result, input-output vector */
);
