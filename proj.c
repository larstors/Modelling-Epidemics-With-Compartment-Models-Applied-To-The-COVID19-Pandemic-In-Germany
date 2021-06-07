#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<tgmath.h>


#include"lib.h"

// FUNTION INITIALIZATION


int function_of_system(double t, double y[], double f[], double **commuters, double *population);

void read_commuters(double **commu, int N);

void fill_pop_array();

void fill_y_array(double* y);



// initializing parameters. For more explanation of these see the Latex/PDF 

const double alpha = 0.2; // rate of infected
const double beta = 1.0/14.0; // recovery rate
const int dimension = 12; // dimension of the system
const double p = 0.0264; // death rate
// TODO find a source or a better estimate for the value (in case it isn't good enough already)
const double t0 = 10.0/24.0; // ratio of "commuters not home" (in hours) to hours a day.


// ####################################################################
// #                            MAIN                                  #
// ####################################################################

int main(){
    // size of matrix
    int N = 38;
    // creating matrix for commuters. Double pointer 
    double **commu = (double**)malloc(sizeof(double) * N);
    // set the double pointer on mallocs
    for (int i = 0; i < N; i++){
        commu[i] = (double*) malloc(sizeof(double) * N);
    }


    // filling the population vector
    double *population = malloc(sizeof(double) * dimension); // array with population of system
    fill_pop_array(population);

    

    // filling the commuter matrix
    read_commuters(commu, N);


    // solving the differential equation using rk4
    // the time period is from July 24. to November 1. 

    // time of starting. Note that the unit is days
    double t_start = 0.0;
    // The simmulation spans 100 days
    double t_end = 100.0;

    // choosing a delta_t
    // TODO maybe adjust this for better resolution/results
    double delta_t = 0.001; // this time step corresponds to about 1.4 minutes (someone please check the math on this)

    // dimension of differential equation (as there are 4 differential equations per cell)
    int dim_deq = 4 * dimension;

    // array with initial conditions
    double y[dim_deq];
    fill_y_array(y);


    // file for saving the solution of calculations
    FILE *sol = fopen("DiffEq.txt", "w");

    // loop to solve the differential equation
    while (t_start < t_end){
        // applying the rk4 step function to calculate the next step in y
        rk4_step(t_start, delta_t, y, function_of_system, dim_deq, commu, population);

        // adding time to the file
        fprintf(sol, "&lf,", t_start);

        // adding all the y values to file
        for (int i = 0; i < dim_deq; i++){
            fprintf(sol, "%lf,", y[i]);
        }

        // new line in file
        fprintf(sol, "\n");


        // increasing time by delta_t
        t_start += delta_t;

    }


    // free the mallocs
    for (int i = 0; i < N; i++){
        free(commu[i]);
    }
    free(commu);
    
    free(population);

    // closing file
    fclose(sol);

    return 0;

}


// ####################################################################
// #                            FUNCTIONS                             #
// ####################################################################


// ################################### COMMUTERS ###################################

/**
 * @brief function to read the txt file with the commuters
 * 
 * @param com Array to fill with the commuters
 * @param N The size of the matrix that will be read. 
 */
void read_commuters(double **commu, int N){

    // Commuters
    // loading txt with commuters
    FILE *com = fopen("Pendler.txt", "r");

    // filling array
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            // Break condition if it can't read in value
            if(!fscanf(com, "%lf\t", &commu[i][j])){
                printf("The matrix didn’t read properly");
                break;
            } 
        }
    }

    fclose(com);
    return;
}

// ##################################### FUNCTION OF SYSTEM ######################################################


/**
 * @brief The function of the system. The system is in this case is modeled by the approach of constant coefficients, i.e. the weighted mean of the standard
 *        and the dynamic model
 * 
 * @param t time at which it is to be calculated
 * @param y array contaning the system quantities [S,...S, I,...I, R,...R, D,...D]
 * @param f array which is to be filled with derivatives
 * @param commuters matrix with commuters
 * @param population array with population
 * @return int 
 */
int function_of_system(double t, double y[], double f[], double **commuters, double *population){
    // arrays for different quantities
    // effective infected
    double *Ieff = (double*)malloc(sizeof(double) * dimension);
    // S
    double *So = (double*)malloc(sizeof(double) * dimension);
    // I
    double *Io = (double*)malloc(sizeof(double) * dimension);
    // R
    double *Ro = (double*)malloc(sizeof(double) * dimension);
    // D
    double *Do = (double*)malloc(sizeof(double) * dimension);
    // effective population
    double *Neff = (double*)malloc(sizeof(double) * dimension);
    
    // filling the arrays
    for (int i = 0; i < dimension; i++){
        So[i] = y[i];
        Io[i] = y[i + dimension];
        Ro[i] = y[i + 2*dimension];
        Do[i] = y[i + 3*dimension];
        Neff[i] = effective_population(commuters, population, dimension)[i];
    }

    // fill the effective infected
    for (int i = 0; i < dimension; i++){
        Ieff[i] = effective_infected(commuters, population, dimension, Io)[i];
    }

    //TODO finish the function
    for (int i = 0; i < dimension; i++){
        // sum for last term in derivative. For detail see PDF/LaTeX
        double sum = 0;
        for (int k = 0; k < dimension; k++){
            sum += commutersFrom(commuters, dimension, i)[k] * Ieff[k];
        }

        // dSdt
        f[i] = - (1 - t0) * alpha * So[i] * Io[i] - t0/population[i] * So[i] * (Neff[i] * Ieff[i] - sum);

        // dIdt
        f[i + dimension] = - f[i] - beta*Io[i];

        // dRdt
        f[i + 2 * dimension] = beta * (1 - p) * Io[i];

        // dDdt
        f[i + 3 * dimension] = beta * p * Io[i];
    }


    // freeing the malloc
    free(Ieff);
    free(So);
    free(Io);
    free(Ro);
    free(Do);
    free(Neff);
    return 0;
}

// ######################################### POPULATION ###################################################
/**
 * @brief Function to read in the population from txt file and save in the array "population"
 * 
 * @param population Array to save the population in
 */
void fill_pop_array(double* population)
{
    FILE* popdata = fopen("Internal Data/popdata38.txt", "r");
    for (int j = 0; j < dimension; ++j)
        fscanf(popdata, "%lf", &population[j]);
    fclose(popdata);
}

// ######################################### INITIAL CONDITIONS ###################################################
/**
 * @brief Function to read in the initial condition (SIRD distribution) from txt file and save in the array "y"
 * 
 * @param y Array to save the initial conditions in
 */
void fill_y_array(double* y)
{
    double *initial_unsorted = malloc(4 * dimension * sizeof(double));
    FILE* initial = fopen("Internal Data/initial_data.txt", "r");
    for (int j = 0; j < 4 * dimension; ++j)
        fscanf(initial, "%lf", &initial_unsorted[j]);
    fclose(initial);
    for (int k = 0; k < dimension; ++k)
    {
        for (int l = 0; l < 4; ++l)
        {
            y[dimension * l + k] = initial_unsorted[4 * k + l];
        }
    }
    free(initial_unsorted);
        
}

