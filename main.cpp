#include "boundary_val.hpp"
#include "helper.hpp"
#include "init.hpp"
#include "sor.hpp"
#include "utilities.hpp"
#include "uvp.hpp"
#include "visual.hpp"
#include <cstdio>
#include <fcntl.h>
#include <iostream>
#include <unistd.h>

/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed. Use the predefined matrix<typename>
 * type and give initial values in the constructor.
 * - perform the main loop
 * - at the end: destroy any memory allocated and print some useful statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two-dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop, the following steps are required (for some of the
 * operations, a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */

int main(int argn, char **args) {

// Reading command line flags and clear output dir
#pragma region
  std::string szFileName;
  std::string problemName;
  std::string vtkFileName;
  std::string logFileName;

  // parse flags entered by user to set the problem type
  ParseCommandLineOptions(argn, args, szFileName, vtkFileName, problemName);
  // clearing output directory for the vtk files, to be sure of the output data
  clear_output_dir();
#pragma endregion

// Statically allocated no need for dynamic allocations for simulation
// parameters
#pragma region
  // Solvers Parameters,
  int itermax;
  double eps;
  double alpha;
  double omg;
  double tau;

  // Grid Data
  double xlength;
  double ylength;
  double dx;
  double dy;
  int imax;
  int jmax;

  // Fluid Data, reynolds number, initial conditions for velocity(x and y
  // direction) and pressure
  double Re;
  double UI;
  double VI;
  double PI;

  // Forces Data, gravitational forces
  double GX;
  double GY;

  // Time Stepping parameters
  double t_end;
  double dt;
  double dt_value;

  // Loop parameters
  double t = 0.0;       // time variable for time grid points
  double t_print = 0.0; // time variable for printing output files
  double res = 1.0;     // residual value for the SOR iteration
  int n = 0;            // counter for the variables value at each step
  int it = 0;           // iteration index for the SOR solver

  // Testing data print parameters
  int step = 1;
  int file_descriptor;
  // Solver required matrices
  matrix<double> F;
  matrix<double> G;
  matrix<double> RS;
  // Grid matrices
  matrix<double> U;
  matrix<double> V;
  matrix<double> P;
#pragma endregion

// Reading parameters
#pragma region
  printf("\x1B[96m Reading Simulation parameters from the data file\033[0m\n");
  // Reading parameters to the static local variables from the parameters file
  // provided by the user
  read_parameters(szFileName, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength,
                  &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau,
                  &itermax, &eps, &dt_value);
#pragma endregion

// Grid initliazation and matrices memory allocation
#pragma region
  // Initialize Grid
  Grid domain = Grid(imax, jmax, boundary_size, PI, UI, VI);

  // Allocating memory for matrices
  F.resize(domain.imaxb(), std::vector<double>(domain.jmaxb(), 0.0));
  G.resize(domain.imaxb(), std::vector<double>(domain.jmaxb(), 0.0));
  RS.resize(domain.imaxb(), std::vector<double>(domain.jmaxb(), 0.0));
  // Allocating memory for matrices
  U.resize(domain.imaxb(), std::vector<double>(domain.jmaxb(), 0.0));
  V.resize(domain.imaxb(), std::vector<double>(domain.jmaxb(), 0.0));
  P.resize(domain.imaxb(), std::vector<double>(domain.jmaxb(), 0.0));
#pragma endregion

  // Setting the boudnary conditions of the grid
  boundaryvalues(imax, jmax, domain);

// Printing Problem Data for user
#pragma region
  printf("\x1B[35m Solving %s with size %.2fx%.2f, Grid: %dx%d Reynolds: "
         "%.2f\033[0m\n",
         problemName.c_str(), xlength, ylength, imax, jmax, Re);
  printf(
      "\x1B[35m Solver Parameters: omg: %.2f alpha: %.2f dt: %.2f  \033[0m\n",
      omg, alpha, dt);
#pragma endregion

  std::ofstream fout(Temp_log_dt_s);

// Solver
#pragma region
  // Main Loop
  while (t < t_end) {

    printf("\x1B[37m t = %.3f n = %d dt = %.5f\033[0m\n", t, n, dt);

    // Calculating suitable step size
    calculate_dt(Re, tau, &dt, dx, dy, imax, jmax, domain);

    // Set boundary values for u and v
    boundaryvalues(imax, jmax, domain);

    // Compute F n and G n
    calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, domain, F, G);

    // Compute the right-hand side rhs of the pressure equation
    calculate_rs(dt, dx, dy, imax, jmax, F, G, RS);

    it = 0;
    res = 1.0;
    // Solve system using SOR
    while (it < itermax && res > eps) {
      sor(omg, dx, dy, imax, jmax, domain, RS, &res);
      it++;
    }

    if (it == itermax)
      printf("\x1B[33m Warning SOR exited due to maximum iteration reached at "
             "step: %d\033[0m\n",
             n);

    if (t >= t_print) {
      t_print += dt_value;
      VTKHelper::printVTKFile(domain, dx, dy, vtkFileName, Temp_Results_s, n);
    }

    calculate_uv(dt, dx, dy, imax, jmax, domain, F, G);

    // printing to log files
    if (t > 45 && (step <= 2)) {
      if (step == 1) {
        fout << std::setprecision(4) << dt;
        fout << std::endl;

        logFileName = Temp_log_F_s + std::to_string(step);
        write_matrix(logFileName, F, domain.imaxb(), domain.jmaxb());

        logFileName = Temp_log_G_s + std::to_string(step);
        write_matrix(logFileName, G, domain.imaxb(), domain.jmaxb());

        logFileName = Temp_log_RS_s + std::to_string(step);
        write_matrix_less_precision(logFileName, RS, domain.imaxb(),
                                    domain.jmaxb());

        logFileName = Temp_log_P_s + std::to_string(step);
        domain.pressure(P);
        write_matrix(logFileName, P, domain.imaxb(), domain.jmaxb());

        logFileName = Temp_log_U_s + std::to_string(step);
        domain.velocity(U, velocity_type::U);
        write_matrix(logFileName, U, domain.imaxb(), domain.jmaxb());

        logFileName = Temp_log_V_s + std::to_string(step);
        domain.velocity(V, velocity_type::V);
        write_matrix(logFileName, V, domain.imaxb(), domain.jmaxb());
        step++;
      } else {

        // fout << std::setprecision(10) <<  dt;
        // fout << std::endl;

        logFileName = Temp_log_F_s + std::to_string(step);
        write_matrix_less_precision(logFileName, F, domain.imaxb(),
                                    domain.jmaxb());
        logFileName =
            Temp_log_F_s + std::to_string(step) + std::to_string(step);
        write_matrix(logFileName, F, domain.imaxb(), domain.jmaxb());

        logFileName = Temp_log_G_s + std::to_string(step);
        write_matrix_less_precision(logFileName, G, domain.imaxb(),
                                    domain.jmaxb());
        logFileName =
            Temp_log_G_s + std::to_string(step) + std::to_string(step);
        write_matrix(logFileName, G, domain.imaxb(), domain.jmaxb());

        logFileName = Temp_log_RS_s + std::to_string(step);
        write_matrix(logFileName, RS, domain.imaxb(), domain.jmaxb());

        logFileName = Temp_log_P_s + std::to_string(step);
        domain.pressure(P);
        write_matrix_less_precision(logFileName, P, domain.imaxb(),
                                    domain.jmaxb());

        logFileName =
            Temp_log_P_s + std::to_string(step) + std::to_string(step);
        domain.pressure(P);
        write_matrix(logFileName, P, domain.imaxb(), domain.jmaxb());

        logFileName = Temp_log_U_s + std::to_string(step);
        domain.velocity(U, velocity_type::U);
        write_matrix_less_precision(logFileName, U, domain.imaxb(),
                                    domain.jmaxb());

        logFileName = Temp_log_V_s + std::to_string(step);
        domain.velocity(V, velocity_type::V);
        write_matrix_less_precision(logFileName, V, domain.imaxb(),
                                    domain.jmaxb());
        step++;
      }
    }

    t = t + dt;
    n = n + 1;
  }

  printf("\x1B[32m Simulation Finished Succesfully in %d time steps\033[0m\n",
         n);

  // Output values for visualization
  VTKHelper::printVTKFile(domain, dx, dy, vtkFileName, Temp_Results_s, n);

  printf("\x1B[32m vtk file written successfully\033[0m\n");
#pragma endregion

  return 0;
}
