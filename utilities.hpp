#include <openssl/md5.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <string.h>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <sstream>
#include "datastructures.hpp"


/** This File is for the utilities, which includes 
 * methods for reading input flags as well the data
 * paths
*/
/*------------Directory paths for the Folder structure in the project-------------*/
// Paths are relative to the build path of the solver executable not to tests
#pragma region 
// Directory of Testing Data, log files for matrices
const std::string Ref_Test_Data_s = "../Ref_Test_Data/";
const std::string Temp_Test_Data_s = "../Temp_Test_Data/";
const std::string Cmp_Test_Data_s = "../Cmp_Test_Data/";
// Directory of Input Data
const std::string Sim_Data_s = "../Simulation_Data/";
const std::string Deb_Data_s = "../Debugging_Data/";
// Directory of Results in vtk format
const std::string Ref_Results_s = "../Ref_Results/";
const std::string Temp_Results_s = "../Temp_Results/";
#pragma endregion
                                            
/*--------------------- Log file paths-------------------------------------------*/
#pragma region 
const std::string Temp_log_U_s = std::string(Temp_Test_Data_s) + "log_U";
const std::string Temp_log_V_s = std::string(Temp_Test_Data_s) + "log_V";
const std::string Temp_log_F_s = std::string(Temp_Test_Data_s) + "log_F";
const std::string Temp_log_G_s = std::string(Temp_Test_Data_s) + "log_G";
const std::string Temp_log_P_s = std::string(Temp_Test_Data_s) + "log_P";
const std::string Temp_log_RS_s = std::string(Temp_Test_Data_s) + "log_RS";
const std::string Temp_log_dt_s = std::string(Temp_Test_Data_s) + "log_dt";
#pragma endregion


/*Directories for the input data and for saving output vtk file*/
const std::string sim_dir = "../";

/** Files names for data and geometry for different example problems
 */
const std::string Lid_driven_cavity_s = "Cavity";


/*help string for the user for the command line flags*/
const std::string help_message = std::string("usage ./sim <problem> \n") +
                                             " problem        -c For Lid driven cavity\n";
                                            
/*-----Boundary sizes-----*/
const int boundary_size = 1;
const int boundary_size_mpi = 2;
/*------------------------*/

/*----------Command line & directories handling methods----------------*/
#pragma region 
/* Method to parse command line input flags*/
void ParseCommandLineOptions(int argc, char *argv[],std::string &sceFileName,std::string &vtkFileName,std::string &problemName);
/* Method to clear output dir*/
void clear_output_dir();
#pragma endregion

/*-----------Printing methods-----------------------------------------*/
#pragma region 
void print_matrix(matrix<double> m,int xdim,int ydim);
#pragma endregion


/*----------I/P & O/P streams----------------*/
#pragma region 
void read_matrix(std::string fileName,matrix<double>& m,int xdim,int ydim);

void write_matrix(std::string fileName,matrix<double>& m,int xdim,int ydim);

void write_matrix_less_precision(std::string fileName,matrix<double>& m,int xdim,int ydim);
#pragma endregion

/*----------computing the hash for log files-------------*/
#pragma region 
std::string get_md5hash( const std::string fname);
#pragma endregion
