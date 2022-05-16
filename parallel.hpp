#include "mpi.h"
#include "grid.hpp"

/** This operation initializes necessary values for the domain decomposition
 * into subdomains over the processes in the communicator 
 *
 * @param iproc      maximum nuber of process in the i index
 * @param jproc      maximum nuber of process in the i index
 * @param imax       number of cells in x-direction
 * @param jmax       number of cells in Y-direction
 * @param myrank     rank of current process
 * @param il         index i to start from left boundary for current subdomain
 * @param ir         index i to end at righ boundary for current subdomain
 * @param jb         index j to start from left boundary for current subdo
 * @param jt         index j to end at righ boundary for current subdomain
 * @param rank_l     rank of left neighbor process
 * @param rank_r     rank of right neighbor process
 * @param rank_b     rank of bottom neighbor process
 * @param rank_t     rank of top neighbor process
 * @param omg_i      i index for the current process
 * @param omg_j      j index for the current process
 * @param num_proc   number of processes of the commuicator(where domain is splitted to)
 */

void Program_Message (char *txt);
/* produces a stderr text output  */


void Programm_Sync (char *txt);
/* produces a stderr textoutput and synchronize all processes */


void Programm_Stop (char *txt);

void init_parallel(int iproc,int jproc,int imax,int jmax,int *myrank,
                   int *il,int *ir,int *jb,int *jt,
                   int *rank_l,int *rank_r,int *rank_b,int *rank_t,
                   int *omg_i,int *omg_j,int num_proc);

/* this method exchanges pressure values between processes that treat
 adjacent sub-domains*/
/* Called within SOR before actually calculating the residual */
void pressure_comm (double **P,
                    int il, int ir,
                    int jb,int jt,
                    int rank_l,int rank_r,
                    int rank_b,int rank_t,
                    double *bufSend,double *bufRecv, 
                    MPI_Status *status, int chunk);


void pressure_comm(Grid& grid,
                  int imax,
                  int jmax,
                  int rank_l,
                  int rank_r,
                  int rank_b,
                  int rank_t,
                  double *bufSend,
                  double *bufRecv,
                  matrix<double> &P);


void uv_comm(Grid& grid,
             int imax,
             int jmax,
             int rank_l,
             int rank_r,
             int rank_b,
             int rank_t,
             int myrank,
             double *bufSend,
             double *bufRecv,
             matrix<double> &U,
             matrix<double> &V);
/* this method exchanges the velocity values U und V between the processes
treating adjacent subdomains*/
/* Called at the end of calculate_uv */
void uv_comm   (double **U,double **V,
                int il,int ir, 
                int jb,int jt,
                int rank_l,int rank_r,
                int rank_b,int rank_t,
                double *bufSend, double *bufRecv, 
                MPI_Status *status, int chunk);
