/* This file is part of CEED. For more details, see exascaleproject.org. */

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
)
{
    int i, j, k1, k2, k3;
    int n = ndof_1d, m = nqpt_1d;
    int ndofs = n*n, nqpts = m*m;
    /* variable-lenght arrays for simplicity */
    double x_loc[ndofs], t[m*n], x_qpt[nqpts*2], vx, vy;

    for (i = 0; i < nelem; i++)
    {
       /* Action of P */
       for (j = 0; j < ndofs; j++)
       {
          x_loc[j] = x[dof_offsets[j+ndofs*i]];
       }

       /* Action of G */

       /* B1d contraction: (n x n) x (n x m) -> (n x m) */
       /*                   x_loc  x  B1d^T  ->    t    */
       /* Loop variables:  k2   k3   k3   k1    k2   k1 */
       for (k1 = 0; k1 < m; k1++)
       {
          for (k2 = 0; k2 < n; k2++)
          {
             t[k2+n*k1] = 0.0;
             for (k3 = 0; k3 < n; k3++)
             {
                t[k2+n*k1] += x_loc[k2+n*k3] * B1d[k1+m*k3];
             }
          }
       }

       /* G1d contraction: (m x n) x (n x m) ->     (m x m)  */
       /*                    G1d   x    t    -> x_qpt[:,:,0] */
       /* Loop variables:  k2   k3   k3   k1        k2   k1  */
       for (k1 = 0; k1 < m; k1++)
       {
          for (k2 = 0; k2 < m; k2++)
          {
             x_qpt[k2+m*k1] = 0.0;
             for (k3 = 0; k3 < n; k3++)
             {
                x_qpt[k2+m*k1] += G1d[k2+m*k3] * t[k3+n*k1];
             }
          }
       }

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

       /* G1d contraction: (m x n) x (n x m) ->     (m x m)  */
       /*                     t   x   G1d^T  -> x_qpt[:,:,1] */
       /* Loop variables:  k2   k3   k3   k1        k2   k1  */
       for (k1 = 0; k1 < m; k1++)
       {
          for (k2 = 0; k2 < m; k2++)
          {
             x_qpt[k2+m*k1+nqpts] = 0.0;
             for (k3 = 0; k3 < n; k3++)
             {
                x_qpt[k2+m*k1+nqpts] += t[k2+m*k3] * G1d[k1+m*k3];
             }
          }
       }

       /* Action of D */
       for (j = 0; j < nqpts; j++)
       {
          vx = x_qpt[j];
          vy = x_qpt[j+nqpts];
          x_qpt[j      ] = D[j+nqpts*(  3*i)] * vx + D[j+nqpts*(1+3*i)] * vy;
          x_qpt[j+nqpts] = D[j+nqpts*(1+3*i)] * vx + D[j+nqpts*(2+3*i)] * vy;
       }

       /* Action of G^T */

       /* G1d contraction: (m x m)  x (m x n) -> (m x n) */
       /*              x_qpt[:,:,1] x   G1d   ->    t    */
       /* Loop variables:  k2   k3    k3   k1    k2   k1 */
       for (k1 = 0; k1 < n; k1++)
       {
          for (k2 = 0; k2 < m; k2++)
          {
             t[k2+m*k1] = 0.0;
             for (k3 = 0; k3 < m; k3++)
             {
                t[k2+m*k1] += x_qpt[k2+m*k3+nqpts] * G1d[k3+m*k1];
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

       /* G1d contraction: (n x m) x     (m x m)  -> (n x m) */
       /*                   G1d^T  x x_qpt[:,:,0] ->    t    */
       /* Loop variables:  k2   k3       k3   k1     k2   k1 */
       for (k1 = 0; k1 < m; k1++)
       {
          for (k2 = 0; k2 < n; k2++)
          {
             t[k2+n*k1] = 0.0;
             for (k3 = 0; k3 < m; k3++)
             {
                t[k2+n*k1] += G1d_t[k2+n*k3] * x_qpt[k3+m*k1];
             }
          }
       }

       /* B1d contraction: (n x m) x (m x n) -> (n x n) */
       /*          x_loc +    t    x   B1d   ->  x_loc  */
       /* Loop variables:  k2   k3   k3   k1    k2   k1 */
       for (k1 = 0; k1 < n; k1++)
       {
          for (k2 = 0; k2 < n; k2++)
          {
             for (k3 = 0; k3 < m; k3++)
             {
                x_loc[k2+n*k1] += t[k2+n*k3] * B1d[k3+m*k1];
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