/*
  Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
  the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
  reserved. See file LICENSE for details.

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

void add_mult_mass_quad(
   int ndof_1d,      /* number of 1D dofs (points) */
   int nqpt_1d,      /* number of 1D quadrature points */
   int nelem,        /* number of elements */
   double *D,        /* nqpt_1d x nqpt_1d x nelem */
   double *B1d,      /* nqpt_1d x ndof_1d dense matrix, column-major layout */
   double *B1d_t,    /* trasnspose of B1d */
   int *dof_offsets, /* array of size ndofs_1d x ndofs_1d x nelem representing a
                        boolean P */
   double *x,        /* input vector */
   double *y         /* result, input-output vector */
)
{
   int i, j, k1, k2, k3;
   int n = ndof_1d, m = nqpt_1d;
   int ndofs = n*n, nqpts = m*m;
   /* variable-lenght arrays for simplicity */
   double x_loc[ndofs], t[m*n], x_qpt[nqpts];

   for (i = 0; i < nelem; i++)
   {
      /* Action of P */
      for (j = 0; j < ndofs; j++)
      {
         x_loc[j] = x[dof_offsets[j+ndofs*i]];
      }

      /* Action of B */

      /* B1d contraction: (m x n) x (n x n) -> (m x n) */
      /*                    B1d   x  x_loc  ->    t    */
      /* Loop variables:  k2   k3   k3   k1    k2   k1 */
      for (k1 = 0; k1 < n; k1++)
      {
         for (k2 = 0; k2 < m; k2++)
         {
            t[k2+m*k1] = 0.0;
            for (k3 = 0; k3 < n; k3++)
            {
               t[k2+m*k1] += B1d[k2+m*k3] * x_loc[k3+n*k1];
            }
         }
      }

      /* B1d contraction: (m x n) x (n x m) -> (m x m) */
      /*                     t   x   B1d^T  ->  x_qpt  */
      /* Loop variables:  k2   k3   k3   k1    k2   k1 */
      for (k1 = 0; k1 < m; k1++)
      {
         for (k2 = 0; k2 < m; k2++)
         {
            x_qpt[k2+m*k1] = 0.0;
            for (k3 = 0; k3 < n; k3++)
            {
               x_qpt[k2+m*k1] += t[k2+m*k3] * B1d[k1+m*k3];
            }
         }
      }

      /* Action of D */
      for (j = 0; j < nqpts; j++)
      {
         x_qpt[j] *= D[j + nqpts*i];
      }

      /* Action of B^T */

      /* B1d contraction: (m x m) x (m x n) -> (m x n) */
      /*                   x_qpt  x   B1d   ->    t    */
      /* Loop variables:  k2   k3   k3   k1    k2   k1 */
      for (k1 = 0; k1 < n; k1++)
      {
         for (k2 = 0; k2 < m; k2++)
         {
            t[k2+m*k1] = 0.0;
            for (k3 = 0; k3 < m; k3++)
            {
               t[k2+m*k1] += x_qpt[k2+m*k3] * B1d[k3+m*k1];
            }
         }
      }

      /* B1d contraction: (n x m) x (m x n) -> (n x n) */
      /*                   B1d^T  x    t    ->  x_loc  */
      /* Loop variables:  k2   k3   k3   k1    k2   k1 */
      for (k1 = 0; k1 < n; k1++)
      {
         for (k2 = 0; k2 < n; k2++)
         {
            x_loc[k2+n*k1] = 0.0;
            for (k3 = 0; k3 < m; k3++)
            {
               x_loc[k2+n*k1] += B1d_t[k2+n*k3] * t[k3+m*k1];
            }
         }
      }

      /* Action of P^T */
      for (j = 0; j < ndofs; j++)
      {
         y[dof_offsets[j+ndofs*i]] += x_loc[j];
      }
   }
}
