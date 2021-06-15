#ifndef MYNUMERICS_H
#define MYNUMERICS_H


// Function that returns array filled with commuters from a cell
void commutersFrom(double **commuters, int dimension, int l, double *a);
// Function that returns array filled with commuters to a cell
void commutersTo(double **commuters, int dimension, int l, double *a);
// Funtion to calculate the effective population in every cell
void effective_population(double **commuters, double *population, int dimension, double *array);
// Function to calculate the effective infected in every cell
void effective_infected(double **commuters, double *population, int dimension, double *infected, double *Ieff);

// typedef for function to be integrated
typedef
int ode_func(double, double[], double[], double**, double*, double);
// Runge Kutta 4 for this system
void rk4_step(double t, double delta_t, double y[], ode_func func, int dim, double **commuters, double *population, double params);


//----------------------- LINEAR ALGEBRA --------------------------/*

/* Gauss elimination algorithm that transforms a nxm matrix into a triangular matrix */
int gauss(int n, int m, double A[n*m], int pivoting);
/* Returns the determinant of a nxn matrix A */
double det(int n, double A[n*n], int pivoting);
/* Calculates the solution x of linear equation Ax=b */
void lgs_solve(int n, double A[n*n],double b[n], double x[n], int pivoting);
/* Calculates the inverse of nxn matrix A */
void inverse(int n, double A[n*n], double A_inverse[n*n]);
/* Calculates the matrix product of nxm matrix A and mxp matrix B */
void matrix_product(int n, int m, int p, double A[n*m], double B[m*p], double AB[n*p]);
/* Calculates the transpose of nxm matrix A */
void transpose(int n, int m, double A[n*m], double A_transpose[m*n]);


#endif