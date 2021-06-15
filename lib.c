#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>


// ############################# COMMUTERS FROM ############################

/**
 * @brief Function to calculate the commuters going out of cell l
 * 
 * @param commuters The matrix with all the commuters
 * @param dimension Integer telling how big the system is
 * @param l         Integer saying from what cell the commuters come from
 * @param a         Matrix to be written in
 */
void commutersFrom(double **commuters, int dimension, int l, double *a){
    

    // fill
    for (int i = 0; i < dimension; i++){
        a[i] = commuters[i][l];
    }

}


// ############################### COMMUTERS TO ################################


/**
 * @brief Function to calculate the commuters going to cell l
 * 
 * @param commuters The matrix with all the commuters
 * @param dimension Integer telling how big the system is
 * @param l         Integer saying to what cell the commuters are going
 * @return double*  Array with the commuters going to cell l
 * @param a         Matrix to be written in
 */
void commutersTo(double **commuters, int dimension, int l, double *a){

    // fill
    for (int i = 0; i < dimension; i++){
        a[i] = commuters[l][i];
    }

}

// ########################## EFFECTIVE POPULATION #######################################

/**
 * @brief function to calculate the effective population in a cell. It is defined as initial population subtracted with the commuters commuting out of the cell.
 * 
 * @param commuters Matrix (38 * 38) filled with commuters
 * @param population Array with population
 * @return double*  Array with effective population
 * @param array         Matrix to be written in
 */

void effective_population(double **commuters, double *population, int dimension, double *array){

    // iterating over every cell and calculating effective population for every cell
    for (int i = 0; i < dimension; i++){

        // setting the value of return array to the population of the respective cell
        array[i] = population[i];

        // iterating over commuters and subtracting them from the population
        for (int k = 0; k < dimension; k++){
            array[i] -= commuters[k][i];
        }

    }
}

// ############################### EFFECTIVE INFECTED #############################


/**
 * @brief Function to calculate the effective infected for every cell in the system
 * 
 * @param commuters   Matrix with commuters
 * @param population  Array with population
 * @param dimension   Dimension of the system
 * @param infected    Array with the infected
 * @return double*    Array with effective infected
 * @param Ieff         Matrix to be written in
 */

void effective_infected(double **commuters, double *population, int dimension, double *infected, double* Ieff){

    // iterating over every cell to fill return array and calculating the effective infected
    for (int i = 0; i < dimension; i++){

        // setting array values
        Ieff[i] = infected[i];

        // setting two help variables 
        double cfrom = 0;
        double cto = 0;

        // loop for adjusting parameters
        for (int k = 0; k < dimension; k++){
            cfrom -= commuters[k][i];
            cto += commuters[i][k]*infected[k];
        }

        // calculating the effective infected
        Ieff[i] += 1.0/population[i]*(cfrom*infected[i] + cto);
    }
    
}


// ################################# TYPEDEF FUNC #####################################

typedef
int ode_func(double, const double[], double[], double**, double*, double);



// ######################################## RK4 #####################################


/**
 * @brief Function to calculate a step in the differential equation using Runge-Kutta 4
 * 
 * @param t Current time
 * @param delta_t Time step
 * @param y State array
 * @param func Function to calculate the derivatives
 * @param dim System size of differential equation
 * @param commuters 2D array with commuters
 * @param population 1D array with population of the cells
 * @param params parameters for the function to be integrated
 */

// TODO expand the rk4 to fit variaton in parameters such as alpha, beta etc.
void rk4_step(double t, double delta_t, double y[], ode_func func, int dim, double **commuters, double *population, double params){

	// initializing arrays
	double *support = (double*)malloc(sizeof(double)*dim);
	double *k1 = (double*)malloc(sizeof(double)*dim);
	double *k2 = (double*)malloc(sizeof(double)*dim);
	double *k3 = (double*)malloc(sizeof(double)*dim);
	double *k4 = (double*)malloc(sizeof(double)*dim);
		
	// filling k1
	func(t, y, k1, commuters, population, params);
	
	// multiplying by delta_t
	for (int k = 0; k < dim; k++) k1[k] *= delta_t;
	
	// filling support
	for (int i = 0; i < dim; i++) support[i] = y[i] + k1[i]/2.0;
	
	// filling k2
	func(t + delta_t/2.0, support, k2, commuters, population, params);
	
	// multiplying by delta_t
	for (int k = 0; k < dim; k++) k2[k] *= delta_t;
	
	// filling support
	for (int i = 0; i < dim; i++) support[i] = y[i] + k2[i]/2.0;
	
	// filling k3
	func(t + delta_t/2.0, support, k3, commuters, population, params);
	
	// multiplying by delta_t
	for (int k = 0; k < dim; k++) k3[k] *= delta_t;
	
	// filling support
	for (int i = 0; i < dim; i++) support[i] = y[i] + k3[i];
	
	// filling k4
	func(t + delta_t, support, k4, commuters, population, params);
	
	// multiplying by delta_t
	for (int k = 0; k < dim; k++) k4[k] *= delta_t;
	
	// updating y
	for (int i = 0; i < dim; i++){
		y[i] += 1.0/6.0*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
	}
	
	// free the mallocs
	free(support);
	free(k1);
	free(k2);
	free(k3);
	free(k4);
	return;

}

/* ##################################### LINEAR ALGEBRA ####################################### */

/**
 * @brief Gauss elimination algorithm that transforms a nxm matrix A into a triangular matrix
 * 
 * @param n number of rows in matrix A
 * @param m number of columns in matrix A
 * @param A 1D array containing the matrix whose (i,j)-element has been mapped to index (i*m +j)
 * @param pivoting if 1: algorithm chooses highest pivot element in each step
 * 
 * @return int	sign of determinant of resulting matrix (always 1 when pivoting is not activated)
 */
int gauss(int n, int m, double A[n*m], int pivoting)
{
	int sign = 1;
	
	for (int i = 0; i < n-1; i++)
	{	
		if(pivoting == 1)
		{
			for(int k = i+1; k < n ; k++)
			{
				if(fabs(A[i*m + i]) < fabs(A[k*m + i]))
				{
					for(int j = 0; j < m; j++)
					{
						double temp = A[i*m + j];
						A[i*m + j] = A[k*m + j];
						A[k*m + j] = temp;
					}
					
					sign = (-1)*sign;
				}
			}
		}
		
		for (int k = i+1; k < n; k++)
		{
			double c = A[k*m + i] / A[i*m + i];
			for (int j = 0; j < m; j++)
			{
				A[k*m + j] = A[k*m + j] - c * A[i*m + j];
			}
		}
	}
	
	return sign;
}

/**
 * @brief Returns the determinant of a nxn matrix A
 * 
 * @param n number of rows/columns in square matrix A
 * @param A 1D array containing the matrix whose (i,j)-element has been mapped to index (i*n +j)
 * @param pivoting if 1: algorithm chooses biggest pivot element pivot element in each step of Gaussian elimination
 * 
 * @return double	determinant of A
 */
double det(int n, double A[n*n], int pivoting)
{
	double B[n*n];
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			B[i*n + j] = A[i*n + j];
		}
	}
	
	int sign = gauss(n, n, B, pivoting);
	
	double product = 1.0;
	for (int i = 0; i < n; i++)
	{
		product = product * B[i*n + i];
	}
	
	return sign * product;	
}

/**
 * @brief Calculates the solution x of linear equation Ax=b
 * 
 * @param n number of rows/columns in square matrix A
 * @param A 1D array containing the matrix whose (i,j)-element has been mapped to index (i*n +j)
 * @param b 1D array containing the components of vector b
 * @param x 1D array that specifies the destination of solution vector x
 * @param pivoting if 1: algorithm chooses biggest pivot element in each step of Gaussian elimination
 * 
 */
void lgs_solve(int n, double A[n*n], double b[n], double x[n], int pivoting)
{
	double C[n*(n + 1)];
	for(int i = 0; i < n; i++)
	{
		for(int j = 0;  j < n; j++)
		{
			C[i*(n+1) + j] = A[i*n + j];
		}
		C[i*(n+1) + n] = b[i];
	}
	
	gauss(n, n+1, C, pivoting);
	
	for(int i = n-1; i >= 0; i--)
	{
		double sum = 0;
		for(int j = i+1; j < n; j++)
		{
			sum += C[i*(n+1) + j] * x[j];
		}
		x[i] = (C[i*(n+1) + n] - sum) / C[i*(n+1) + i];
	}
}

/**
 * @brief Calculates the inverse of nxn matrix A
 * 
 * @param n number of rows/columns in square matrix A
 * @param A 1D array containing the matrix whose (i,j)-element has been mapped to index (i*n +j)
 * @param A_inverse 1D array that specifies the destination of inverse whose (i,j)-element is mapped to index (i*n +j)
 * 
 */
void inverse(int n, double A[n*n], double A_inverse[n*n])
{
	if (det(n, A, 1) == 0)
	{
		printf("How dare you trying to invert a matrix with vanishing determinant!");
	}
	else
	{
		for(int j = 0; j < n; j++)
		{
			double b[n];
			double x[n];
			
			for(int i = 0; i < n; i++)
			{
				b[i] = 0;
			}
			b[j] = 1;
			
			lgs_solve(n, A, b, x, 1);
			
			for(int i = 0; i < n; i++)
			{
				A_inverse[i*n + j] = x[i];
			}
			
		}
	}
}

/**
 * @brief Calculates the matrix product of nxm matrix A and mxp matrix B
 * 
 * @param n number of rows in matrix A
 * @param m number of columns in matrix A / rows in matrix B
 * @param p number of columns in matrix B
 * @param A 1D array containing the matrix A whose (i,j)-element has been mapped to index (i*m +j)
 * @param B 1D array containing the matrix B whose (i,j)-element has been mapped to index (i*p +j)
 * @param AB 1D array that specifies the destination of matrix product whose (i,j)-element is mapped to index (i*p +j)
 * 
 */
void matrix_product(int n, int m, int p, double A[n*m], double B[m*p], double AB[n*p])
{
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < p; j++)
		{
			AB[i*p + j] = 0;
			
			for(int k = 0; k < m; k++)
			{
				AB[i*p + j] += A[i*m + k] * B[k*p + j];
			}
		}
	}
	
}

/**
 * @brief Calculates the transpose of nxm matrix A
 * 
 * @param n number of rows in matrix A
 * @param m number of columns in matrix A
 * @param A 1D array containing the matrix A whose (i,j)-element has been mapped to index (i*m +j)
 * @param A_transpose 1D array that specifies the destination of transpose whose (i,j)-element is mapped to index (i*n +j)
 * 
 */
void transpose(int n, int m, double A[n*m], double A_transpose[m*n])
{
	for(int i = 0; i < m; i++)
	{
		for(int j = 0; j < n; j++)
		{
			A_transpose[i*n +j] = A[j*m +i];
		}
	}
}
