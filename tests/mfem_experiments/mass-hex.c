// This file is part of CEED. For more details, see exascaleproject.org.

void add_mult_mass_hex(
   int ndof_1d,      /* number of 1D dofs (points) */
   int nqpt_1d,      /* number of 1D quadrature points */
   int nelem,        /* number of elements */
   double *D,        /* (nqpt_1d)^3 x nelem */
   double *B1d,      /* nqpt_1d x ndof_1d dense matrix, column-major layout */
   double *B1d_t,    /* trasnspose of B1d */
   int *dof_offsets, /* array of size (ndofs_1d)^3 x nelem representing a
                        boolean P */
   double *x,        /* input vector */
   double *y         /* result, input-output vector */
)
{
   int i, j, k1, k2, k3, kz;
   int n = ndof_1d, m = nqpt_1d;
   int nn = n*n, ndofs = n*nn, mm = m*m, nqpts = m*mm;
   /* variable-lenght arrays, for simplicity */
#if 0
   double x_loc[ndofs], t1[m*nn], t2[mm*n], x_qpt[nqpts];
#else
   double x_qpt[nqpts], t2[mm*n], *t1 = x_qpt, *x_loc = t2;
   if (m < n)
   {
      fprintf(stderr, "\n"
              "add_mult_mass_hex: the case m < n is not implemented."
              " abort.\n");
      abort();
   }
#endif

   for (i = 0; i < nelem; i++)
   {
      /* Action of P */
      for (j = 0; j < ndofs; j++)
      {
         x_loc[j] = x[dof_offsets[j + ndofs*i]];
      }

      /* Action of B */

      /* B1d contraction: (m x n) x (n x (n x n)) -> (m x (n x n)) */
      /*                    B1d   x    x_loc      ->      t1       */
      /* Loop variables:  k2   k3   k3     k1        k2     k1     */
      for (k1 = 0; k1 < nn; k1++)
      {
         for (k2 = 0; k2 < m; k2++)
         {
            t1[k2+m*k1] = 0.0;
            for (k3 = 0; k3 < n; k3++)
            {
               t1[k2+m*k1] += B1d[k2+m*k3] * x_loc[k3+n*k1];
            }
         }
      }

      /* B1d contraction: (m x n)   x (n x m) ->  (m x m)   */
      /*                 t1[:,:,kz] x  B1d^T  -> t2[:,:,kz] */
      /* Loop variables:  k2   k3     k3   k1     k2   k1      */
      for (kz = 0; kz < n; kz++)
      {
         for (k1 = 0; k1 < m; k1++)
         {
            for (k2 = 0; k2 < m; k2++)
            {
               t2[k2+m*(k1+m*kz)] = 0.0;
               for (k3 = 0; k3 < n; k3++)
               {
                  t2[k2+m*(k1+m*kz)] += t1[k2+m*(k3+n*kz)] * B1d[k1+m*k3];
               }
            }
         }
      }

      /* B1d contraction: ((m x m) x n) x (n x m) -> ((m x m) x m) */
      /*                        t2         B1d^T  ->     x_qpt     */
      /* Loop variables:     k2      k3   k3   k1       k2      k1 */
      for (k1 = 0; k1 < m; k1++)
      {
         for (k2 = 0; k2 < mm; k2++)
         {
            x_qpt[k2+mm*k1] = 0.0;
            for (k3 = 0; k3 < n; k3++)
            {
               x_qpt[k2+mm*k1] += t2[k2+mm*k3] * B1d[k1+m*k3];
            }
         }
      }

      /* Action of D */
      for (j = 0; j < nqpts; j++)
      {
         x_qpt[j] *= D[j + nqpts*i];
      }

      /* Action of B^T */

      /* B1d contraction: ((m x m) x m) x (m x n) -> ((m x m) x n) */
      /*                      x_qpt     x   B1d   ->       t2      */
      /* Loop variables:     k2      k3   k3   k1       k2      k1 */
      for (k1 = 0; k1 < n; k1++)
      {
         for (k2 = 0; k2 < mm; k2++)
         {
            t2[k2+mm*k1] = 0.0;
            for (k3 = 0; k3 < m; k3++)
            {
               t2[k2+mm*k1] += x_qpt[k2+mm*k3] * B1d[k3+m*k1];
            }
         }
      }

      /* B1d contraction: (m x m) x (m x n) ->  (m x n)   */
      /*                 t2[:,:,kz]   B1d   -> t1[:,:,kz] */
      /* Loop variables:  k2   k3   k3   k1     k2   k1   */
      for (kz = 0; kz < n; kz++)
      {
         for (k1 = 0; k1 < n; k1++)
         {
            for (k2 = 0; k2 < m; k2++)
            {
               t1[k2+m*(k1+n*kz)] = 0.0;
               for (k3 = 0; k3 < m; k3++)
               {
                  t1[k2+m*(k1+n*kz)] += t2[k2+m*(k3+m*kz)] * B1d[k3+m*k1];
               }
            }
         }
      }

      /* B1d contraction: (n x m) x (m x (n x n)) -> (n x (n x n)) */
      /*                   B1d^T  x      t1       ->    x_loc      */
      /* Loop variables:  k2   k3   k3     k1        k2     k1     */
      for (k1 = 0; k1 < nn; k1++)
      {
         for (k2 = 0; k2 < n; k2++)
         {
            x_loc[k2+n*k1] = 0.0;
            for (k3 = 0; k3 < m; k3++)
            {
               /* x_loc[k2+n*k1] += B1d[k3+m*k2] * t1[k3+m*k1]; */
               x_loc[k2+n*k1] += B1d_t[k2+n*k3] * t1[k3+m*k1];
            }
         }
      }

      /* Action of P^T */
      for (j = 0; j < ndofs; j++)
      {
         y[dof_offsets[j + ndofs*i]] += x_loc[j];
      }
   }
}
