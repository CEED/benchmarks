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
);
