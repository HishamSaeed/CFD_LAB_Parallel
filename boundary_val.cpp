#include "boundary_val.hpp"
#include "datastructures.hpp"
#include "grid.hpp"


void boundaryvalues(int imax, int jmax, Grid& grid) {

    // Setting the boundary values for the vertical boundaries

    // Setting the inner cells boundary conditions First
    // right Wall boundary conditions
    for(int j = 1; j <= jmax; j++)
    {
      // Horizontal Velocities
      grid.cell(imax,j).velocity(velocity_type::U) = 0.0; 
    }

    for(int i = 1; i <= imax; i++)
    {
      // Vertical velocities
      grid.cell(i,jmax).velocity(velocity_type::V) = 0.0;
    }

    // Left and Right wall
    for(int j = 1; j <= jmax; j++)
    {
          //left Wall
        {
          // Vertical Velocities
          grid.cell(0,j).velocity(velocity_type::V) = (-1.0) * (grid.cell(1,j).velocity(velocity_type::V));
          
          // Horizontal Velocities
          grid.cell(0,j).velocity(velocity_type::U) = 0.0;
        }    
        //right wall
        {
          // Vertical Velocities     
          grid.cell(imax+1,j).velocity(velocity_type::V) = (-1.0) * grid.cell(imax,j).velocity(velocity_type::V); 
        }
    }
    // Top and Bottom wall
    for(int i = 1; i <= imax; i++)
    {
           // bottom wall
       {
         // Vertical velocities
         grid.cell(i,0).velocity(velocity_type::V) = 0.0;

         // Horizontal Velocities
         grid.cell(i,0).velocity(velocity_type::U) = (-1.0) * grid.cell(i,1).velocity(velocity_type::U);
       }  
       // moving wall
       {            
           // Horizontal Velocitiy     
           grid.cell(i,jmax+1).velocity(velocity_type::U) = 2.0 - grid.cell(i,jmax).velocity(velocity_type::U);
       } 
    }

    // Setting Boundary values for right wall
    // for(int j = 1; j <= jmax; j++)
    // {
    //   // Vertical Velocities     
    //   grid.cell(imax+1,j).velocity(velocity_type::V) = (-1.0) * grid.cell(imax,j).velocity(velocity_type::V);
    //   // Horizontal Velocities
    //   grid.cell(imax,j).velocity(velocity_type::U) = 0.0;  
    // }

    // // Setting Boundary values for Top wall
    // for(int i = 1; i <= imax; i++)
    // {
    //   // Vertical velocities
    //   grid.cell(i,jmax).velocity(velocity_type::V) = 0.0;
    //   // Horizontal Velocitiy     
    //   grid.cell(i,jmax+1).velocity(velocity_type::U) = 2.0 - grid.cell(i,jmax).velocity(velocity_type::U);
    // }

    // // Setting Boundary values for Left wall
    // for(int j = 1;j <= jmax; j++)
    // {
    //   // Vertical Velocities
    //   grid.cell(0,j).velocity(velocity_type::V) = (-1.0) * (grid.cell(1,j).velocity(velocity_type::V));      
    //   // Horizontal Velocities
    //   grid.cell(0,j).velocity(velocity_type::U) = 0.0;
    // }

    // // Setting Boundary values for Bottom wall
    // for(int i = 1; i <= imax; i++)
    // {
    //   // Vertical velocities
    //   grid.cell(i,0).velocity(velocity_type::V) = 0.0;
    //   // Horizontal Velocities
    //   grid.cell(i,0).velocity(velocity_type::U) = (-1.0) * grid.cell(i,1).velocity(velocity_type::U);
    // }



    // for(int j = 1; j <= jmax; j++){
    //     //left Wall
    //     {
    //       // Vertical Velocities
    //       grid.cell(0,j).velocity(velocity_type::V) = (-1.0) * (grid.cell(1,j).velocity(velocity_type::V));
          
    //       // Horizontal Velocities
    //       grid.cell(0,j).velocity(velocity_type::U) = 0.0;
    //     }    
    //     //right wall
    //     {
    //       // Vertical Velocities     
    //       grid.cell(imax+1,j).velocity(velocity_type::V) = (-1.0) * grid.cell(imax,j).velocity(velocity_type::V);

    //       // Horizontal Velocities
    //       grid.cell(imax,j).velocity(velocity_type::U) = 0.0;   
    //     }
    // }


    // // Setting the boundary values for the horizontal boundaries
    // for(int i =1; i <= imax; i++){
    //    // bottom wall
    //    {
    //      // Vertical velocities
    //      grid.cell(i,0).velocity(velocity_type::V) = 0.0;

    //      // Horizontal Velocities
    //      grid.cell(i,0).velocity(velocity_type::U) = (-1.0) * grid.cell(i,1).velocity(velocity_type::U);
    //    }  
    //    // moving wall
    //    {
    //        // Vertical velocities
    //        grid.cell(i,jmax).velocity(velocity_type::V) = 0.0;
            
    //        // Horizontal Velocitiy     
    //        grid.cell(i,jmax+1).velocity(velocity_type::U) = 2.0 - grid.cell(i,jmax).velocity(velocity_type::U);
    //    }   
    // }
}


