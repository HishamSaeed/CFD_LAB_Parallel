#include "parallel.hpp"
#include <mpi.h>
#include <iomanip>



void init_parallel(int iproc,
                   int jproc,
                   int imax,
                   int jmax,
                   int *myrank,
                   int *il,
                   int *ir,
                   int *jb,
                   int *jt,
                   int *rank_l,
                   int *rank_r,
                   int *rank_b,
                   int *rank_t,
                   int *omg_i,
                   int *omg_j,
                   int num_proc)
{
    int process;
    int root = 0;
    if(*myrank == 0){
        for(process = num_proc-1;process >= 0;process--)
        {
            *omg_i = (process % iproc) + 1;
            *omg_j = ((process+1-*omg_i)/iproc)+1;

            *il = (*omg_i-1)*(imax/iproc) + 1;
            *ir = (*omg_i!=iproc)?((*omg_i)*(imax/iproc)):imax;

            *jb = (*omg_j-1)*(jmax/jproc) + 1;
            *jt = (*omg_j!=jproc)?((*omg_j)*(jmax/jproc)):jmax;

            if(*il == 1)      *rank_l = MPI_PROC_NULL;
            else              *rank_l = process - 1;

            if(*ir == imax)   *rank_r = MPI_PROC_NULL;
            else              *rank_r = process + 1;


            if(*jb == 1)      *rank_b = MPI_PROC_NULL;
            else              *rank_b = process - iproc;

            if(*jt == jmax)   *rank_t = MPI_PROC_NULL;
            else              *rank_t = process + iproc;

            if(process == 0 ) continue;

            // Sending computed values to other processes
            MPI_Send(omg_i, 1, MPI_INT,process,0,MPI_COMM_WORLD);
            MPI_Send(omg_j, 1, MPI_INT,process,0,MPI_COMM_WORLD);
            MPI_Send(il, 1, MPI_INT, process,0,MPI_COMM_WORLD);
            MPI_Send(ir, 1, MPI_INT, process,0,MPI_COMM_WORLD);
            MPI_Send(jb, 1, MPI_INT, process,0,MPI_COMM_WORLD);
            MPI_Send(jt, 1, MPI_INT, process,0,MPI_COMM_WORLD);
            MPI_Send(rank_l, 1, MPI_INT, process, 0, MPI_COMM_WORLD);
            MPI_Send(rank_r, 1, MPI_INT, process, 0, MPI_COMM_WORLD);
            MPI_Send(rank_b, 1, MPI_INT, process, 0, MPI_COMM_WORLD);
            MPI_Send(rank_t, 1, MPI_INT, process, 0, MPI_COMM_WORLD);
        }
    }
    else{
        MPI_Recv(omg_i, 1, MPI_INT, root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(omg_j, 1, MPI_INT, root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(il, 1, MPI_INT, root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(ir, 1, MPI_INT, root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(jb, 1, MPI_INT, root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(jt, 1, MPI_INT, root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(rank_l, 1, MPI_INT, root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(rank_r, 1, MPI_INT, root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(rank_b, 1, MPI_INT, root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(rank_t, 1, MPI_INT, root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    printf("Thread_id: %d omg_ij: %d%d \nil: %d, ir: %d, jb: %d, jt: %d \n \n",*myrank,*omg_i,*omg_j, *il,*ir,*jb,*jt);
    printf("Thread_id: %d rank_l: %d, rank_r: %d, rank_b: %d, rank_t: %d \n \n",*myrank, *rank_l,*rank_r,*rank_b,*rank_t);

    
}



void pressure_comm(Grid& grid,
                  int imax,
                  int jmax,
                  int rank_l,
                  int rank_r,
                  int rank_b,
                  int rank_t,
                  double *bufSend,
                  double *bufRecv,
                  matrix<double> &P)
{

  


  
  /*---------------------------Sending from right to left--------------------------------------------------------*/
  // Sending to left
  if (rank_l != MPI_PROC_NULL)
	{
		for(int j=2; j<=jmax+1; j++)
		{
			bufSend[j-2] = P[2][j];
		}
		MPI_Send( bufSend, jmax, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD );
	}

  // Receive from right
  if (rank_r != MPI_PROC_NULL) 
	{
		MPI_Recv( bufRecv, jmax, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
		for(int j=2; j<=jmax+1; j++)
		{
			P[imax+2][j] = bufRecv[j-2];
		}
 
  }

  /*---------------------------Sending from left to right--------------------------------------------------------*/
  //Sending to right
  if (rank_r != MPI_PROC_NULL)
  {
		for(int j=2; j<=jmax+2; j++)
		{
			bufSend[j-2] = P[imax+1][j];
		}
		MPI_Send( bufSend, jmax, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD );
	}

  // Receive from left
	if (rank_l != MPI_PROC_NULL)
	{
		MPI_Recv(bufRecv, jmax, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
		for(int j=2; j<=jmax+1; j++)
		{
			P[1][j] = bufRecv[j-2];
		}
  
	}

  /*---------------------------Sending from bottom to top--------------------------------------------------------*/
  // Send to top
  if (rank_t != MPI_PROC_NULL)
	{
		for(int i=2; i<=imax+1; i++)
		{
			bufSend[i-2] = P[i][jmax+1];
		}
		MPI_Send( bufSend, imax, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD );
	}
  // Receive from bottom
	if (rank_b != MPI_PROC_NULL)
	{
		MPI_Recv( bufRecv, imax, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
		for(int i=2; i<=imax+1; i++)
		{
			// P[i][jmax+2] = bufRecv[i-2];
      P[i][1] = bufRecv[i-2];
		}
   
  }

  /*---------------------------Sending from bottom to top--------------------------------------------------------*/
  // Sending to bottom
   if (rank_b != MPI_PROC_NULL)
	{
		for(int i=2; i<=imax+1; i++)
		{
			bufSend[i-2] = P[i][2];
		}
		MPI_Send( bufSend, imax, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD );
	}
  // Receive from top
	if (rank_t != MPI_PROC_NULL) // Receive from the top
	{
		MPI_Recv( bufRecv, imax, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
		for(int i=2; i<=imax+1; i++)
		{
			P[i][jmax+2] = bufRecv[i-2];
		}
 
	}

  
}



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
             matrix<double> &V)
{

  if(rank_l != MPI_PROC_NULL && rank_r != MPI_PROC_NULL)
  {
    /*---------------------From Left to right---------------------*/

    /*-----------------------Velocity U---------------------*/
    
    /*-------------------------------------------------------*/
    for (int j = 2; j <= jmax+1; ++j)
    {
      bufSend[j-2] = U[imax][j];
    }
    MPI_Sendrecv(bufSend, jmax, MPI_DOUBLE, rank_r, 1,
                 bufRecv, jmax, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    for (int j = 2; j <= jmax+1; ++j)
    {
      U[0][j] = bufRecv[j-2];
    }

    /*-------------------------------------------------------*/
    for (int j = 2; j <= jmax+1; ++j)
    {
      bufSend[j-2] = U[imax+1][j];
    }
    MPI_Sendrecv(bufSend, jmax, MPI_DOUBLE, rank_r, 1,
                 bufRecv, jmax, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    for (int j = 2; j <= jmax+1; ++j)
    {
      U[1][j] = bufRecv[j-2];
    }
    /*-------------------------------------------------------*/

    /*-----------------------Vertival Velocities V---------------------*/
    
    /*-------------------------------------------------------*/
    for (int j = 2; j <= jmax+1; ++j)
    {
      bufSend[j-2] = V[imax][j];
    }
    MPI_Sendrecv(bufSend, jmax, MPI_DOUBLE, rank_r, 1,
                 bufRecv, jmax, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    for (int j = 2; j <= jmax+1; ++j)
    {
      V[0][j] = bufRecv[j-2];
    }

    /*-------------------------------------------------------*/
    for (int j = 2; j <= jmax+1; ++j)
    {
      bufSend[j-2] = V[imax+1][j];
    }
    MPI_Sendrecv(bufSend, jmax, MPI_DOUBLE, rank_r, 1,
                 bufRecv, jmax, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    for (int j = 2; j <= jmax+1; ++j)
    {
      V[1][j] = bufRecv[j-2];
    }
    /*-------------------------------------------------------*/



    /*---------------------From right to left---------------------*/

    /*-----------------------Velocity U---------------------*/
    
    /*-------------------------------------------------------*/
    for (int j = 2; j <= jmax+1; ++j)
    {
      bufSend[j-2] = U[2][j];
    }
    MPI_Sendrecv(bufSend, jmax, MPI_DOUBLE, rank_l, 1,
                 bufRecv, jmax, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    for (int j = 2; j <= jmax+1; ++j)
    {
      U[imax+2][j] = bufRecv[j-2];
    }

    /*-------------------------------------------------------*/
    for (int j = 2; j <= jmax+1; ++j)
    {
      bufSend[j-2] = U[3][j];
    }
    MPI_Sendrecv(bufSend, jmax, MPI_DOUBLE, rank_l, 1,
                 bufRecv, jmax, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    for (int j = 2; j <= jmax+1; ++j)
    {
      U[imax+3][j] = bufRecv[j-2];
    }
    /*-------------------------------------------------------*/

    /*-----------------------Velocity V---------------------*/
    
    /*-------------------------------------------------------*/
    for (int j = 2; j <= jmax+1; ++j)
    {
      bufSend[j-2] = V[2][j];
    }
    MPI_Sendrecv(bufSend, jmax, MPI_DOUBLE, rank_l, 1,
                 bufRecv, jmax, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    for (int j = 2; j <= jmax+1; ++j)
    {
      V[imax+2][j] = bufRecv[j-2];
    }

    /*-------------------------------------------------------*/
    for (int j = 2; j <= jmax+1; ++j)
    {
      bufSend[j-2] = V[3][j];
    }
    MPI_Sendrecv(bufSend, jmax, MPI_DOUBLE, rank_l, 1,
                 bufRecv, jmax, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    for (int j = 2; j <= jmax+1; ++j)
    {
      V[imax+3][j] = bufRecv[j-2];
    }
    /*-------------------------------------------------------*/
  }
  else if(rank_l == MPI_PROC_NULL && rank_r != MPI_PROC_NULL)
  {
    /*---------------------From Left to right---------------------*/

    /*---------------------- Velocity U -------------------- */
    
    /*-------------------------------------------------------*/
    for (int j = 2; j <= jmax+1; ++j)
    {
      bufSend[j-2] = U[imax][j];
    }
    MPI_Send(bufSend, jmax, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD);
    
    /*-------------------------------------------------------*/
    for (int j = 2; j <= jmax+1; ++j)
    {
      bufSend[j-2] = U[imax+1][j];
    }
    MPI_Send(bufSend, jmax, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD);

    /*---------------------- Velocity V -------------------- */
    
    /*-------------------------------------------------------*/
    for (int j = 2; j <= jmax+1; ++j)
    {
      bufSend[j-2] = V[imax][j];
    }
    MPI_Send(bufSend, jmax, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD);
    
    /*-------------------------------------------------------*/
    for (int j = 2; j <= jmax+1; ++j)
    {
      bufSend[j-2] = V[imax+1][j];
    }
    MPI_Send(bufSend, jmax, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD);

    /*---------------------From right to left---------------------*/

    /*----------------------- Velocity U --------------------------*/
    
    /*------------------------------------------------------------*/
    MPI_Recv(bufRecv, jmax, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int j = 2; j <= jmax+1; ++j)
    {
      U[imax+2][j] = bufRecv[j-2];
    }
    /*-----------------------------------------------------------*/
    MPI_Recv(bufRecv, jmax, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int j = 2; j <= jmax+1; ++j)
    {
      U[imax+3][j] = bufRecv[j-2];
    }
    /*----------------------- Velocity V --------------------------*/
    
    /*------------------------------------------------------------*/
    MPI_Recv(bufRecv, jmax, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int j = 2; j <= jmax+1; ++j)
    {
      V[imax+2][j] = bufRecv[j-2];
    }
    /*-----------------------------------------------------------*/
    MPI_Recv(bufRecv, jmax, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int j = 2; j <= jmax+1; ++j)
    {
      V[imax+3][j] = bufRecv[j-2];
    }


  }
  else if(rank_l != MPI_PROC_NULL && rank_r == MPI_PROC_NULL)
  {
    /*---------------------From Left to right---------------------*/

    /*----------------------- Velocity U --------------------------*/
    
    MPI_Recv(bufRecv, jmax, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int j = 2; j <= jmax+1; ++j)
    {
      U[0][j] = bufRecv[j-2];
    }
    /*---------------------------------------------------*/
    MPI_Recv(bufRecv, jmax, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int j = 2; j <= jmax+1; ++j)
    {
      U[1][j] = bufRecv[j-2];
    }

    /*----------------------- Velocity V --------------------------*/
    
    MPI_Recv(bufRecv, jmax, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int j = 2; j <= jmax+1; ++j)
    {
      V[0][j] = bufRecv[j-2];
    }
    /*---------------------------------------------------*/
    MPI_Recv(bufRecv, jmax, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int j = 2; j <= jmax+1; ++j)
    {
      V[1][j] = bufRecv[j-2];
    }

    /*---------------------From right to left---------------------*/

    /*----------------------- Velocity U --------------------------*/
      // Send First Ghost layer
    for (int j = 2; j <= jmax+1; ++j)
    {
      bufSend[j-2] = U[2][j];
    }
    MPI_Send(bufSend, jmax, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD);
    // Send Second Ghost layer
    for (int j = 2; j <= jmax+1; ++j)
    {
      bufSend[j-2] = U[3][j];
    }
    MPI_Send(bufSend, jmax, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD);

    /*----------------------- Velocity V --------------------------*/
      // Send First Ghost layer
    for (int j = 2; j <= jmax+1; ++j)
    {
      bufSend[j-2] = V[2][j];
    }
    MPI_Send(bufSend, jmax, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD);
    // Send Second Ghost layer
    for (int j = 2; j <= jmax+1; ++j)
    {
      bufSend[j-2] = V[3][j];
    }
    MPI_Send(bufSend, jmax, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD);

  }



  if(rank_b != MPI_PROC_NULL && rank_t != MPI_PROC_NULL)
   {
     /*------------From Top to bottom--------------------*/

     /*------------Horizontal Velocities U ---------------*/

     /*---------------------------------------------------*/
     for(int i = 2; i <= imax+1; i++)
     {
       bufSend[i-2] = U[i][2];
     }
     MPI_Sendrecv(bufSend, imax, MPI_DOUBLE, rank_b, 1,
     bufRecv, imax, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     for (int i = 2; i <=imax+1; i++)
     {
       U[i][jmax+2] = bufRecv[i-2];
     }
     /*---------------------------------------------------*/
     for(int i = 2; i <= imax+1; i++)
     {
       bufSend[i-2] = U[i][3];
     }
     MPI_Sendrecv(bufSend, imax, MPI_DOUBLE, rank_b, 1,
     bufRecv, imax, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     for (int i = 2; i <=imax+1; i++)
     {
       U[i][jmax+3] = bufRecv[i-2];
     }
     /*---------------------------------------------------*/

     /*------------Vertical Velocities V ---------------*/

     /*---------------------------------------------------*/
     for(int i = 2; i <= imax+1; i++)
     {
       bufSend[i-2] = V[i][2];
     }
     MPI_Sendrecv(bufSend, imax, MPI_DOUBLE, rank_b, 1,
     bufRecv, imax, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     for (int i = 2; i <=imax+1; i++)
     {
       V[i][jmax+2] = bufRecv[i-2];
     }
     /*---------------------------------------------------*/
     for(int i = 2; i <= imax+1; i++)
     {
       bufSend[i-2] = V[i][3];
     }
     MPI_Sendrecv(bufSend, imax, MPI_DOUBLE, rank_b, 1,
     bufRecv, imax, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     for (int i = 2; i <=imax+1; i++)
     {
       V[i][jmax+3] = bufRecv[i-2];
     }
     /*---------------------------------------------------*/


     /*------------From bottom to Top--------------------*/

     /*------------Horizontal Velocities U ---------------*/

     /*---------------------------------------------------*/
     for(int i = 2; i <= imax+1; i++)
     {
       bufSend[i-2] = U[i][jmax];
     }
     MPI_Sendrecv(bufSend, imax, MPI_DOUBLE, rank_t, 1,
     bufRecv, imax, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     for (int i = 2; i <=imax+1; i++)
     {
       U[i][0] = bufRecv[i-2];
     }
     /*---------------------------------------------------*/
     for(int i = 2; i <= imax+1; i++)
     {
       bufSend[i-2] = U[i][jmax+1];
     }
     MPI_Sendrecv(bufSend, imax, MPI_DOUBLE, rank_t, 1,
     bufRecv, imax, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     for (int i = 2; i <=imax+1; i++)
     {
       U[i][1] = bufRecv[i-2];
     }

     /*------------Vertical Velocities V ---------------*/

     /*---------------------------------------------------*/
     for(int i = 2; i <= imax+1; i++)
     {
       bufSend[i-2] = V[i][jmax];
     }
     MPI_Sendrecv(bufSend, imax, MPI_DOUBLE, rank_t, 1,
     bufRecv, imax, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     for (int i = 2; i <=imax+1; i++)
     {
       V[i][0] = bufRecv[i-2];
     }
     /*---------------------------------------------------*/
     for(int i = 2; i <= imax+1; i++)
     {
       bufSend[i-2] = V[i][jmax+1];
     }
     MPI_Sendrecv(bufSend, imax, MPI_DOUBLE, rank_t, 1,
     bufRecv, imax, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     for (int i = 2; i <=imax+1; i++)
     {
       V[i][1] = bufRecv[i-2];
     } 


    // grid.set_velocity(U,velocity_type::U);
    // grid.set_velocity(V,velocity_type::V);

   }
  else if(rank_t == MPI_PROC_NULL && rank_b != MPI_PROC_NULL)
   {

    /*------------From Top to bottom--------------------*/

    /*------------Horizontal Velocities U ---------------*/

    /*---------------------------------------------------*/
     
    for(int i = 2; i <= imax+1; i++)
    {
      bufSend[i-2] = U[i][2];
    }
    MPI_Send(bufSend, imax, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD);

    /*---------------------------------------------------*/
    for(int i = 2; i <= imax+1; i++)
    {
      bufSend[i-2] = U[i][3];
    }
    MPI_Send(bufSend, imax, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD);
    /*---------------------------------------------------*/

    /*------------Vertical Velocities V -----------------*/

    /*---------------------------------------------------*/
    for(int i = 2; i <= imax+1; i++)
    {
      bufSend[i-2] = V[i][2];
    }
    MPI_Send(bufSend, imax, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD);

    /*---------------------------------------------------*/
    for(int i = 2; i <= imax+1; i++)
    {
      bufSend[i-2] = V[i][3];
    }
    MPI_Send(bufSend, imax, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD);


    /*------------From bottom to Top---------------------*/

    /*------------Horizontal Velocities U ---------------*/

    /*---------------------------------------------------*/
    MPI_Recv(bufRecv, imax, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i = 2; i <=imax+1; i++)
    {
      U[i][0] = bufRecv[i-2];
    }

    /*---------------------------------------------------*/
    MPI_Recv(bufRecv, imax, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i = 2; i <=imax+1; i++)
    {
      U[i][1] = bufRecv[i-2];
    }
    /*---------------------------------------------------*/

    /*------------Vertical Velocities V -----------------*/

    /*---------------------------------------------------*/
    MPI_Recv(bufRecv, imax, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i = 2; i <=imax+1; i++)
    {
      V[i][0] = bufRecv[i-2];
    }

    /*---------------------------------------------------*/
    MPI_Recv(bufRecv, imax, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i = 2; i <=imax+1; i++)
    {
      V[i][1] = bufRecv[i-2];
    }



    // grid.set_velocity(U,velocity_type::U);
    // grid.set_velocity(V,velocity_type::V);

   }
  else if(rank_t != MPI_PROC_NULL && rank_b == MPI_PROC_NULL)
   {

     /*------------From Top to bottom--------------------*/

     /*------------Horizontal Velocities U ---------------*/

     /*---------------------------------------------------*/

    MPI_Recv(bufRecv, imax, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i = 2; i <=imax+1; i++)
    {
      U[i][jmax+2] = bufRecv[i-2];
    }

    /*---------------------------------------------------*/
    MPI_Recv(bufRecv, imax, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i = 2; i <=imax+1; i++)
    {
      U[i][jmax+3] = bufRecv[i-2];
    }
    
    /*------------Vertical Velocities V ----------------*/

    /*--------------------------------------------------*/
    MPI_Recv(bufRecv, imax, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i = 2; i <=imax+1; i++)
    {
      V[i][jmax+2] = bufRecv[i-2];
    }

    /*---------------------------------------------------*/
    MPI_Recv(bufRecv, imax, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i = 2; i <=imax+1; i++)
    {
      V[i][jmax+3] = bufRecv[i-2];
    }


    /*------------From bottom to Top--------------------*/

    /*------------Horizontal Velocities U --------------*/

    /*--------------------------------------------------*/
    for(int i = 2; i <= imax+1; i++)
    {
      bufSend[i-2] = U[i][jmax];
    }
    MPI_Send(bufSend, imax, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD);

    /*--------------------------------------------------*/
    for(int i = 2; i <= imax+1; i++)
    {
      bufSend[i-2] = U[i][jmax+1];
    }
    MPI_Send(bufSend, imax, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD);

    /*------------Vertical Velocities V ----------------*/

    /*--------------------------------------------------*/
    for(int i = 2; i <= imax+1; i++)
    {
      bufSend[i-2] = V[i][jmax];
    }
    MPI_Send(bufSend, imax, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD);

    /*--------------------------------------------------*/
    for(int i = 2; i <= imax+1; i++)
    {
      bufSend[i-2] = V[i][jmax+1];
    }
    MPI_Send(bufSend, imax, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD);



    // grid.set_velocity(U,velocity_type::U);
    // grid.set_velocity(V,velocity_type::V);

   }


// free(U);
// free(V);

}            


