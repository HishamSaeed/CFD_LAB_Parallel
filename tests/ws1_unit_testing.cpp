// Definition of sample test
#define CATCH_CONFIG_RUNNER

#include "catch.hpp"
#include "../init.hpp"
#include "../grid.hpp"
#include "../utilities.hpp"
#include "../boundary_val.hpp"
#include "../uvp.hpp"
#include "../sor.hpp"
#include "test_utilities.hpp"

int main( int argc, char* argv[] ) {
  clear_output_dir_test();  
  int result = Catch::Session().run( argc, argv );

  return result;
}


TEST_CASE("Test grid parameters", "[read_grid_parameters]"){
    
    // Dummy Grid Paramters
    #pragma region 
    
    int imax = 10;
    int jmax = 10;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;

    #pragma endregion
    // Initialize Dummy Grid
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI);
    // Test Dummy Grid
    #pragma region
    REQUIRE ( domain.boundary_size() == boundary_size );
    REQUIRE ( domain.imax() == imax );
    REQUIRE ( domain.jmax() == jmax );
    REQUIRE ( domain.imaxb() == (jmax + (2*boundary_size)) );
    REQUIRE ( domain.jmaxb() == (jmax + (2*boundary_size)) );

    for (size_t i = 0; i < domain.imaxb(); i++)
    {
        for (size_t j = 0; j < domain.jmaxb(); j++)
        {
            REQUIRE ( domain.cell(i, j).velocity(velocity_type::U) == UI );
            REQUIRE ( domain.cell(i, j).velocity(velocity_type::V) == VI );
            REQUIRE ( domain.cell(i, j).pressure() == PI );
        }
        
    }

    for (size_t i = 0; i < domain.imaxb(); i++)
    {
        REQUIRE ( domain.cell(i, 0).border(border_position::BOTTOM) == true );
        REQUIRE ( domain.cell(i, domain.jmaxb()-1).border(border_position::TOP) == true );
    }
    
    for (size_t j = 0; j < domain.jmaxb(); j++)
    {
        REQUIRE ( domain.cell(0, j).border(border_position::LEFT) == true );
        REQUIRE ( domain.cell(domain.imaxb()-1, j).border(border_position::RIGHT) == true );
    }
    #pragma endregion

}

TEST_CASE( "Test read parameters", "[read_parameters] [!mayfail]" ) {
 
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

  // Fluid Data, reynolds number, initial conditions for velocity(x and y direction) and pressure
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
  #pragma endregion

  read_parameters(Lid_Driven_Cavity_t,&Re,&UI,&VI,&PI,&GX,&GY,&t_end,&xlength,&ylength,&dt,&dx,&dy,&imax,&jmax,&alpha,&omg,&tau,&itermax,&eps,&dt_value);


  #pragma region   
  CHECK( xlength == 1 );
  CHECK( ylength == 1 );
  CHECK( imax == 50 );
  CHECK( jmax == 50 );
  CHECK( dt == 0.05 );
  CHECK( t_end == 50.0 );
  CHECK( tau == 0.5 );
  CHECK( dt_value == 0.5 );  
  CHECK( itermax == 100 );
  CHECK( eps == 0.001 );
  CHECK( omg == 1.7 );
  CHECK( alpha == 0.5 );
  CHECK( Re == 100 );
  CHECK( GX == 0.0 );
  CHECK( GY == 0.0 );
  CHECK( PI == 0.0 );
  CHECK( UI == 0.0 );
  CHECK( VI == 0.0 );
  #pragma endregion  
}

TEST_CASE("Test boundary values", "[boundaryvalues] [!mayfail]")
{
    // temp variables
    double ghost_cell;
    double inner_cell;

    // tolerance
    double eps = 0.00001;

    // Dummy Grid
    double PI = 0.0;
    double UI = 3.0;
    double VI = 4.0;

    int imax = 10;
    int jmax = 10;

    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI) ;

    boundaryvalues(imax,jmax,domain);

    // Testing upper and lower boundary
    for(int i = 1; i <= imax; i++)
    {
        ghost_cell = domain.cell(i,0).velocity(velocity_type::U);
        inner_cell = domain.cell(i,1).velocity(velocity_type::U);
        REQUIRE( fabs(ghost_cell + inner_cell) < eps);

        ghost_cell = domain.cell(i,0).velocity(velocity_type::V);
        REQUIRE( ghost_cell == 0.0);

        ghost_cell = domain.cell(i,imax+1).velocity(velocity_type::U);
        inner_cell = domain.cell(i,imax).velocity(velocity_type::U);
        REQUIRE( fabs(ghost_cell + inner_cell) == 2);

        ghost_cell = domain.cell(i,imax).velocity(velocity_type::V);
        REQUIRE( ghost_cell == 0.0);
    }

    // Testing left and right boundary
    for(int j = 1; j <= jmax; j++)
    {
        ghost_cell = domain.cell(0,j).velocity(velocity_type::V);
        inner_cell = domain.cell(1,j).velocity(velocity_type::V);
        REQUIRE( fabs(ghost_cell + inner_cell) < eps);

        ghost_cell = domain.cell(0,j).velocity(velocity_type::U);
        REQUIRE( ghost_cell == 0.0);

        ghost_cell = domain.cell(jmax+1,j).velocity(velocity_type::V);
        inner_cell = domain.cell(jmax,j).velocity(velocity_type::V);
        REQUIRE( fabs(ghost_cell + inner_cell) < eps);

        ghost_cell = domain.cell(jmax,j).velocity(velocity_type::U);
        REQUIRE( ghost_cell == 0.0);
    }
}

TEST_CASE("Test calculate dt ", "[calculate_dt]")
{
    // Dummy Grid Parameters
    #pragma region 
    double PI;
    double UI;
    double VI;

    double Re = 100;

    int imax = 10;
    int jmax = 10;

    double dx = 0.1;
    double dy = 0.1;

    double tau = 0.5;
    double dt = 0.5;
    #pragma endregion
    // Initialzie Dummy Grid
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI) ;
    // Data Structure Declaration
    #pragma region 
    matrix<double> U;
    matrix<double> V;
    matrix<double> P;
    #pragma endregion

    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();

     
    /** Dummy Scenario 1
     * dx/umax = inf
     * dy/vmax = inf
     * (Re/2)*(1/dx) + (1/dy) = 0.25
     * result tau * 0.25 = 0.125
     * */
    #pragma region 
    PI = 0.0;
    UI = 0.0;
    VI = 0.0;
    U.resize(imaxb,std::vector<double>(jmaxb,UI));
    V.resize(imaxb,std::vector<double>(jmaxb,VI));
    P.resize(imaxb,std::vector<double>(jmaxb,PI));

    domain.set_velocity(U,velocity_type::U);
    domain.set_velocity(V,velocity_type::V);
    domain.set_pressure(P);
    calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,domain);

    REQUIRE( dt == 0.125);
    
    P.clear();
    U.clear();
    V.clear();
    #pragma endregion

    /** Dummy Scenario 2
     * dx/umax = inf
     * dy/vmax = 0.025
     * (Re/2)*(1/dx) + (1/dy) = 0.25
     * result tau * 0.025 = 0.0125
     * */
    #pragma region 
    PI = 0.0;
    UI = 0.0;
    VI = 4.0;
    P.resize(imaxb,std::vector<double>(jmaxb,PI));
    U.resize(imaxb,std::vector<double>(jmaxb,UI));
    V.resize(imaxb,std::vector<double>(jmaxb,VI));
    
    domain.set_pressure(P);
    domain.set_velocity(U,velocity_type::U);
    domain.set_velocity(V,velocity_type::V);  

    calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,domain);

    REQUIRE( dt == 0.0125);

    P.clear();
    U.clear();
    V.clear();
    #pragma endregion

    /** Dummy Scenario 3
     * dx/umax = 0.05
     * dy/vmax = inf
     * (Re/2)*(1/dx) + (1/dy) = 0.25
     * result tau * 0.05 = 0.025
     * */
    #pragma region 
    PI = 0.0;
    UI = 2.0;
    VI = 0.0;
    U.resize(imaxb,std::vector<double>(jmaxb,UI));
    V.resize(imaxb,std::vector<double>(jmaxb,VI));
    P.resize(imaxb,std::vector<double>(jmaxb,PI));

    domain.set_velocity(U,velocity_type::U);
    domain.set_velocity(V,velocity_type::V);
    domain.set_pressure(P);
    calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,domain);

    REQUIRE( dt == 0.025);

    P.clear();
    U.clear();
    V.clear();
    #pragma endregion

    /* Real Scenario*/
    #pragma region 
    double dt_ref;

    // File and directory paths
    std::string logFileName;

    std::ifstream fin(ref_log_dt_t);
    fin >> dt_ref;

    U.resize(imaxb,std::vector<double>(jmaxb,UI));
    V.resize(imaxb,std::vector<double>(jmaxb,VI));
    
    // Reading log files in data structures
    #pragma region 
    logFileName = ref_log_U_t + "1";
    read_matrix(logFileName.c_str(),U,imaxb,jmaxb);
    logFileName = ref_log_V_t + "1";
    read_matrix(logFileName.c_str(),V,imaxb,jmaxb);
    #pragma endregion

    domain.set_velocity(U, velocity_type::U);
    domain.set_velocity(V, velocity_type::V);

    calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,domain);

    std::ofstream fout(cmp_log_dt_t);
    fout << std::setprecision(4) << dt;
    fout << std::endl;

    const std::string hashdt = get_md5hash(ref_log_dt_t);
    const std::string hash_cmpdt = get_md5hash(cmp_log_dt_t);

    REQUIRE (hashdt == hash_cmpdt);

    U.clear();
    V.clear();
    #pragma endregion
}

TEST_CASE( "Test F and G matrix", "[calculate_FG]" )
{
    std::cout << "  " << std::endl;
    std::cout<< "================================================================="<<std::endl;
    std::cout<< "---------------------Testing calculate_fg()----------------------"<<std::endl;
    std::cout<< "================================================================="<<std::endl;
    std::cout << "  " << std::endl;
    std::string szFileName = "../../cavity100.dat";

    // Dummy Grid Paramters
    #pragma region 
    
    int imax = 10;
    int jmax = 10;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;

    #pragma endregion
    // Initialize Dummy Grid
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI);
    // Dummy Scenario Parameters
    #pragma region 
    double Re = 100;

    double dx = 0.1;
    double dy = 0.1;

    double GX = 0.0;
    double GY = 0.0;

    double alpha = 0.5;
    double dt = 0.0803324;
    double eps = 0.001;

    #pragma endregion


    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();

    // File and directory paths
    std::string logFileName;

    // Data Structure Declaration
    #pragma region 
    // Data structures for first time step
    matrix<double> F1;
    matrix<double> G1;
    matrix<double> RS1;
    matrix<double> P1;
    matrix<double> U1;
    matrix<double> V1;

    // Data structures for second time step
    matrix<double> F2;
    matrix<double> G2;
    matrix<double> RS2;
    matrix<double> P2;
    matrix<double> U2;
    matrix<double> V2;
    #pragma endregion

    // Data Structure Allocation
    #pragma region 
    // Allocating memory for matrices
    F1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    P1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    U1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V1.resize(imaxb,std::vector<double>(jmaxb,0.0));

    // Allocating memory for matrices
    F2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    P2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    U2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion

    // Reading log files in data structures
    #pragma region 
    logFileName = ref_log_F_t + "1";
    read_matrix(logFileName.c_str(),F1,imaxb,jmaxb);
    logFileName = ref_log_G_t + "1";
    read_matrix(logFileName.c_str(),G1,imaxb,jmaxb);
    logFileName = ref_log_RS_t + "1";
    read_matrix(logFileName.c_str(),RS1,imaxb,jmaxb);
    logFileName = ref_log_P_t + "1";
    read_matrix(logFileName.c_str(),P1,imaxb,jmaxb);
    logFileName = ref_log_U_t + "1";
    read_matrix(logFileName.c_str(),U1,imaxb,jmaxb);
    logFileName = ref_log_V_t + "1";
    read_matrix(logFileName.c_str(),V1,imaxb,jmaxb);
    #pragma endregion


    domain.set_velocity(U1, velocity_type::U);
    domain.set_velocity(V1, velocity_type::V);
    domain.set_pressure(P1);

    // std::cout << "-------------" << std::endl;
    logFileName = ref_log_F_t + "2";
    read_matrix(logFileName.c_str(),F2,imaxb,jmaxb);
    const std::string hashF2 = get_md5hash(logFileName);

    logFileName = ref_log_G_t + "2";
    read_matrix(logFileName.c_str(),G2,imaxb,jmaxb);
    const std::string hashG2 = get_md5hash(logFileName);

    const double normf2 = frobeniusNorm(F2, imaxb, jmaxb);
    const double normg2 = frobeniusNorm(G2, imaxb, jmaxb);
    
    // Dummy Scenario
    #pragma region
    

    boundaryvalues(imax, jmax, domain);
    calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, domain, F2, G2);

    #pragma endregion
    


    write_matrix_less_precision(cmp_log_F_t, F2, imaxb, jmaxb);
    write_matrix_less_precision(cmp_log_G_t, G2, imaxb, jmaxb);

    const double norm_cmpf2 = frobeniusNorm(F2, imaxb, jmaxb);
    const double norm_cmpg2 = frobeniusNorm(G2, imaxb, jmaxb);

    const std::string hash_cmpF2 = get_md5hash(cmp_log_F_t);
    const std::string hash_cmpG2 = get_md5hash(cmp_log_G_t);

    #pragma region

    CHECKED_IF (hashF2 == hash_cmpF2)
    {
        std::cout<<hashF2<<" == "<<hash_cmpF2<<std::endl;   
    }
    CHECKED_ELSE (hashF2 == hash_cmpF2) 
    {
        norm_diff_check(normf2, norm_cmpf2, eps, "F");
    }

    CHECKED_IF (hashG2 == hash_cmpG2)
    {
        std::cout<<hashG2<<" == "<<hash_cmpG2<<std::endl; 
    }
    CHECKED_ELSE (hashG2 == hash_cmpG2) 
    {
        norm_diff_check(normg2, norm_cmpg2, eps, "G");
    }

    #pragma endregion
    
}

TEST_CASE( "Test rhs of PPE matrix", "[calculate_rs]" )
{
    std::cout << "  " << std::endl;
    std::cout<< "================================================================="<<std::endl;
    std::cout<< "---------------------Testing calculate_rhs()----------------------"<<std::endl;
    std::cout<< "================================================================="<<std::endl;
    std::cout << "  " << std::endl;
    std::string szFileName = "../../cavity100.dat";


    // Dummy Grid Paramters
    #pragma region 
    
    int imax = 10;
    int jmax = 10;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;

    #pragma endregion
    // Initialize Dummy Grid
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI);
    // Dummy Scenario Parameters
    #pragma region 
    double Re = 100;

    double dx = 0.1;
    double dy = 0.1;

    double GX = 0.0;
    double GY = 0.0;

    double alpha = 0.5;
    // double dt = 0.0803324;
    double dt = 0.08033235832; 
    double eps = 0.001;
    #pragma endregion
    
    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();



    // File and directory paths
    std::string logFileName;

    // Data structure Declaration
    #pragma region 
    // Data structures for first time step
    matrix<double> F1;
    matrix<double> G1;
    matrix<double> RS1;
    matrix<double> P1;
    matrix<double> U1;
    matrix<double> V1;

    // Data structures for second time step
    matrix<double> F2;
    matrix<double> G2;
    matrix<double> RS2;
    matrix<double> P2;
    matrix<double> U2;
    matrix<double> V2;
    #pragma endregion

    // Data Structure initialization
    #pragma region 
    // Allocating memory for matrices
    F1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion

    // Reading log files in data structure
    #pragma region 
    logFileName = ref_log_F_t + "1";
    read_matrix(logFileName.c_str(),F1,imaxb,jmaxb);
    logFileName = ref_log_G_t + "1";
    read_matrix(logFileName.c_str(),G1,imaxb,jmaxb);
    logFileName = ref_log_RS_t + "1";
    #pragma endregion    


    logFileName = ref_log_RS_t + "1";
    read_matrix(logFileName.c_str(),RS1,imaxb,jmaxb);
    const std::string hashRS1 = get_md5hash(logFileName);

    const double normrs1 = frobeniusNorm(RS1, imaxb, jmaxb);
    
    #pragma region
    calculate_rs(dt, dx, dy, imax, jmax, F1, G1, RS1);

    write_matrix_less_precision(cmp_log_RS_t, RS1, imaxb, jmaxb);

    const std::string hash_cmpRS1 = get_md5hash(cmp_log_RS_t);

    const double norm_cmprs1 = frobeniusNorm(RS1, imaxb, jmaxb);

    CHECKED_IF (hashRS1 == hash_cmpRS1)
    {
        std::cout<<hashRS1<<" == "<<hash_cmpRS1<<std::endl; 
    }
    CHECKED_ELSE (hashRS1 == hash_cmpRS1) 
    {
        norm_diff_check(normrs1, norm_cmprs1, eps, "RS");
    }
    #pragma endregion
    
}

TEST_CASE( "Test SOR", "[sor_iteration]" )
{
    std::cout << "  " << std::endl;
    std::cout<< "================================================================="<<std::endl;
    std::cout<< "---------------------------Testing sor()-------------------------"<<std::endl;
    std::cout<< "================================================================="<<std::endl;
    std::cout << "  " << std::endl;
    std::string szFileName = "../../cavity100.dat";


    // Dummy Grid Paramters
    #pragma region 
    
    int imax = 10;
    int jmax = 10;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;

    #pragma endregion
    // Initialize Dummy Grid
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI);
    // Dummy Scenario Parameters
    #pragma region 
    double Re = 100;

    double dx = 0.1;
    double dy = 0.1;

    double GX = 0.0;
    double GY = 0.0;

    double alpha = 0.5;
    double omg = 1.7;
    double eps = 0.001;
    // double dt = 0.5;
    double dt = 0.08033235832; 

    int itermax = 100;

    #pragma endregion
    
    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();


    // File and directory paths
    std::string logFileName;

    // Data structure Declaration
    #pragma region 
    // Data structures for first time step
    matrix<double> F1;
    matrix<double> G1;
    matrix<double> RS1;
    matrix<double> P1;
    matrix<double> U1;
    matrix<double> V1;

    // Data structures for second time step
    matrix<double> F2;
    matrix<double> G2;
    matrix<double> RS2;
    matrix<double> P2;
    matrix<double> U2;
    matrix<double> V2;
    #pragma endregion

    // Data structure Initialization
    #pragma region 
    // Allocating memory for matrices
    F1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    P1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    U1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V1.resize(imaxb,std::vector<double>(jmaxb,0.0));

    // Allocating memory for matrices
    F2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    P2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    U2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion
    

    logFileName = ref_log_RS_t + "1";
    read_matrix(logFileName.c_str(),RS1,imaxb,jmaxb);
    logFileName = ref_log_P_t + "1";
    read_matrix(logFileName.c_str(),P1,imaxb,jmaxb);
    logFileName = ref_log_U_t + "1";
    read_matrix(logFileName.c_str(),U1,imaxb,jmaxb);
    logFileName = ref_log_V_t + "1";
    read_matrix(logFileName.c_str(),V1,imaxb,jmaxb);

    logFileName = ref_log_RS_t + "2";
    read_matrix(logFileName.c_str(),RS2,imaxb,jmaxb);
    logFileName = ref_log_P_t + "2";
    read_matrix(logFileName.c_str(),P2,imaxb,jmaxb);
    const std::string hashP2 = get_md5hash(logFileName);

    
    domain.set_velocity(U1, velocity_type::U);
    domain.set_velocity(V1, velocity_type::V);
    domain.set_pressure(P1);

    const double normp2 = frobeniusNorm(P2, imaxb, jmaxb);

    // boundaryvalues(imax, jmax, domain);
    double res = 1.0;
    int it = 0;
    while(it < itermax && res > eps){
      sor(omg,dx,dy,imax,jmax,domain,RS2,&res);
      it++;
    }

    domain.pressure(P2);
    write_matrix_less_precision(cmp_log_P_t, P2, imaxb, jmaxb);

    const double norm_cmpp2 = frobeniusNorm(P2, imaxb, jmaxb);

    const std::string hash_cmpP2 = get_md5hash(cmp_log_P_t);

    CHECKED_IF (hashP2 == hash_cmpP2)
    {
        std::cout<<hashP2<<" == "<<hash_cmpP2<<std::endl; 
    }
    CHECKED_ELSE (hashP2 == hash_cmpP2) 
    {
        norm_diff_check(normp2, norm_cmpp2, eps, "P");
    }  

}

TEST_CASE("Test velocity (U and v) values", "[calculate_uv]")
{
    std::cout << "  " << std::endl;
    std::cout<< "================================================================="<<std::endl;
    std::cout<< "---------------------Testing calculate_uv()----------------------"<<std::endl;
    std::cout<< "================================================================="<<std::endl;
    std::cout << "  " << std::endl;
    std::string szFileName = "../../cavity100.dat";

    // Dummy Grid Paramters
    #pragma region 
    
    int imax = 10;
    int jmax = 10;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;

    #pragma endregion
    // Initialize Dummy Grid
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI);
    // Dummy Scenario Parameters
    #pragma region 
    double Re = 100;

    double dx = 0.1;
    double dy = 0.1;

    double GX = 0.0;
    double GY = 0.0;

    double alpha = 0.5;
    double omg = 1.7;
    double eps = 0.001;
    double dt = 0.08033235832; 
    
    int itermax = 100;
    #pragma endregion 

    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();

    // File and directory paths
    std::string logFileName;

    // Data Structure Declaration
    #pragma region 
    // Data structures for first time step
    matrix<double> U1;
    matrix<double> V1;

    // Data structures for first time step
    matrix<double> F2;
    matrix<double> G2;
    matrix<double> P2;
    matrix<double> U2;
    matrix<double> V2;
    #pragma endregion

    // Data Structure Initialization
    #pragma region 
    // Allocating memory for matrices
    U1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V1.resize(imaxb,std::vector<double>(jmaxb,0.0));

    // Allocating memory for matrices
    F2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    P2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    U2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion

    // Reading log files in data structures
    #pragma region 
    logFileName = ref_log_U_t + "1";
    read_matrix(logFileName.c_str(),U1,imaxb,jmaxb);
    logFileName = ref_log_V_t + "1";
    read_matrix(logFileName.c_str(),V1,imaxb,jmaxb);
    logFileName = ref_log_F_t + "22";
    read_matrix(logFileName.c_str(),F2,imaxb,jmaxb);
    logFileName = ref_log_G_t + "22";
    read_matrix(logFileName.c_str(),G2,imaxb,jmaxb);
    logFileName = ref_log_P_t + "22";
    read_matrix(logFileName.c_str(),P2,imaxb,jmaxb);
    #pragma endregion

    // Calculating ref files hasing
    #pragma region 
    logFileName = ref_log_U_t + "2";
    read_matrix(logFileName.c_str(),U2,imaxb,jmaxb);
    const std::string hashU2 = get_md5hash(logFileName);
    logFileName = ref_log_V_t + "2";
    read_matrix(logFileName.c_str(),V2,imaxb,jmaxb);
    const std::string hashV2 = get_md5hash(logFileName);
    #pragma endregion

    const double normu2 = frobeniusNorm(U2, imaxb, jmaxb);
    const double normv2 = frobeniusNorm(V2, imaxb, jmaxb);

    // Setting the Grid values to reference values
    #pragma region 
    domain.set_velocity(U1, velocity_type::U);
    domain.set_velocity(V1, velocity_type::V);
    domain.set_pressure(P2);
    #pragma endregion

    boundaryvalues(imax, jmax, domain);
    calculate_uv(dt, dx, dy, imax, jmax, domain, F2, G2);

    domain.velocity(U2, velocity_type::U);
    domain.velocity(V2, velocity_type::V);


    write_matrix_less_precision(cmp_log_U_t, U2, imaxb, jmaxb);
    write_matrix_less_precision(cmp_log_V_t, V2, imaxb, jmaxb);

    const double norm_cmpu2 = frobeniusNorm(U2, imaxb, jmaxb);
    const double norm_cmpv2 = frobeniusNorm(V2, imaxb, jmaxb);


    const std::string hash_cmpU2 = get_md5hash(cmp_log_U_t);
    const std::string hash_cmpV2 = get_md5hash(cmp_log_V_t);

    CHECKED_IF (hashU2 == hash_cmpU2)
    {
        std::cout<<hashU2<<" == "<<hash_cmpU2<<std::endl; 
    }
    CHECKED_ELSE (hashU2 == hash_cmpU2) 
    {
        norm_diff_check(normu2, norm_cmpu2, eps, "U");
    }

    CHECKED_IF (hashV2 == hash_cmpV2)
    {
        std::cout<<hashV2<<" == "<<hash_cmpV2<<std::endl; 
    }
    CHECKED_ELSE (hashV2 == hash_cmpV2) 
    {
        norm_diff_check(normv2, norm_cmpv2, eps, "V");
    }


}


TEST_CASE( "Test F & G integration with rhs", "[integrate F&G and RHS]" )
{
    std::cout << "  " << std::endl;
    std::cout<< "================================================================="<<std::endl;
    std::cout<< "---------------Testing integration of F&G and RHS----------------"<<std::endl;
    std::cout<< "================================================================="<<std::endl;
    std::cout << "  " << std::endl;
    std::string szFileName = "../../cavity100.dat";

    // Dummy Grid Paramters
    #pragma region 
    
    int imax = 10;
    int jmax = 10;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;

    #pragma endregion
    // Initialize Dummy Grid
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI);
    // Dummy Scenario Parameters
    #pragma region 
    double Re = 100;

    double dx = 0.1;
    double dy = 0.1;

    double GX = 0.0;
    double GY = 0.0;

    double alpha = 0.5;
    double dt = 0.0803324;

    double eps = 0.001;

    #pragma endregion


    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();

    // File and directory paths
    std::string logFileName;

    // Data Structure Declaration
    #pragma region 
    // Data structures for first time step
    matrix<double> F1;
    matrix<double> G1;
    matrix<double> RS1;
    matrix<double> P1;
    matrix<double> U1;
    matrix<double> V1;

    // Data structures for second time step
    matrix<double> F2;
    matrix<double> G2;
    matrix<double> RS2;
    matrix<double> P2;
    matrix<double> U2;
    matrix<double> V2;
    #pragma endregion

    // Data Structure Allocation
    #pragma region 
    // Allocating memory for matrices
    F1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    P1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    U1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V1.resize(imaxb,std::vector<double>(jmaxb,0.0));

    // Allocating memory for matrices
    F2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    P2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    U2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion

    // Reading log files in data structures
    #pragma region 
    logFileName = ref_log_F_t + "1";
    read_matrix(logFileName.c_str(),F1,imaxb,jmaxb);
    logFileName = ref_log_G_t + "1";
    read_matrix(logFileName.c_str(),G1,imaxb,jmaxb);
    logFileName = ref_log_RS_t + "1";
    read_matrix(logFileName.c_str(),RS1,imaxb,jmaxb);
    logFileName = ref_log_P_t + "1";
    read_matrix(logFileName.c_str(),P1,imaxb,jmaxb);
    logFileName = ref_log_U_t + "1";
    read_matrix(logFileName.c_str(),U1,imaxb,jmaxb);
    logFileName = ref_log_V_t + "1";
    read_matrix(logFileName.c_str(),V1,imaxb,jmaxb);
    #pragma endregion


    domain.set_velocity(U1, velocity_type::U);
    domain.set_velocity(V1, velocity_type::V);
    domain.set_pressure(P1);


    logFileName = ref_log_RS_t + "2";
    read_matrix(logFileName.c_str(),RS2,imaxb,jmaxb);
    
    logFileName = ref_log_RS_t + "2_integ";
    write_matrix_less_precision(logFileName, RS2, imaxb, jmaxb);

    const double normrs2 = frobeniusNorm(RS2, imaxb, jmaxb);

    const std::string hashRS2 = get_md5hash(logFileName);
    
    // Dummy Scenario
    #pragma region
    

    boundaryvalues(imax, jmax, domain);
    calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, domain, F2, G2);
    calculate_rs(dt, dx, dy, imax, jmax, F2, G2, RS2);

    #pragma endregion
    write_matrix_less_precision(cmp_log_RS_t, RS2, imaxb, jmaxb);

    const double norm_cmprs2 = frobeniusNorm(RS2, imaxb, jmaxb);

    const std::string hash_cmpRS2 = get_md5hash(cmp_log_RS_t);
    
    #pragma region
    
    CHECKED_IF (hashRS2 == hash_cmpRS2)
    {
        std::cout<<hashRS2<<" == "<<hash_cmpRS2<<std::endl; 
    }
    CHECKED_ELSE (hashRS2 == hash_cmpRS2) 
    {
        norm_diff_check(normrs2, norm_cmprs2, eps, "RS");
    }

    #pragma endregion
    
}


TEST_CASE( "Test F & G integration with SOR", "[integrate F&G and SOR]" )
{
    std::cout << "  " << std::endl;
    std::cout<< "================================================================="<<std::endl;
    std::cout<< "---------------Testing integration of F&G and SOR----------------"<<std::endl;
    std::cout<< "================================================================="<<std::endl;
    std::cout << "  " << std::endl;
    std::string szFileName = "../../cavity100.dat";

    // Dummy Grid Paramters
    #pragma region 
    
    int imax = 10;
    int jmax = 10;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;

    #pragma endregion
    // Initialize Dummy Grid
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI);
    // Dummy Scenario Parameters
    #pragma region 
    double Re = 100;

    double dx = 0.1;
    double dy = 0.1;

    double GX = 0.0;
    double GY = 0.0;

    double alpha = 0.5;
    double dt = 0.0803324;

    int itermax = 100;
    double eps = 0.001;
    double omg = 1.7;

    #pragma endregion


    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();

    // File and directory paths
    std::string logFileName;

    // Data Structure Declaration
    #pragma region 
    // Data structures for first time step
    matrix<double> F1;
    matrix<double> G1;
    matrix<double> RS1;
    matrix<double> P1;
    matrix<double> U1;
    matrix<double> V1;

    // Data structures for second time step
    matrix<double> F2;
    matrix<double> G2;
    matrix<double> RS2;
    matrix<double> P2;
    matrix<double> U2;
    matrix<double> V2;
    #pragma endregion

    // Data Structure Allocation
    #pragma region 
    // Allocating memory for matrices
    F1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    P1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    U1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V1.resize(imaxb,std::vector<double>(jmaxb,0.0));

    // Allocating memory for matrices
    F2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    P2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    U2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion

    // Reading log files in data structures
    #pragma region 
    logFileName = ref_log_F_t + "1";
    read_matrix(logFileName.c_str(),F1,imaxb,jmaxb);
    logFileName = ref_log_G_t + "1";
    read_matrix(logFileName.c_str(),G1,imaxb,jmaxb);
    logFileName = ref_log_RS_t + "1";
    read_matrix(logFileName.c_str(),RS1,imaxb,jmaxb);
    logFileName = ref_log_P_t + "1";
    read_matrix(logFileName.c_str(),P1,imaxb,jmaxb);
    logFileName = ref_log_U_t + "1";
    read_matrix(logFileName.c_str(),U1,imaxb,jmaxb);
    logFileName = ref_log_V_t + "1";
    read_matrix(logFileName.c_str(),V1,imaxb,jmaxb);
    #pragma endregion


    domain.set_velocity(U1, velocity_type::U);
    domain.set_velocity(V1, velocity_type::V);
    domain.set_pressure(P1);


    logFileName = ref_log_P_t + "2";
    const std::string hashP2 = get_md5hash(logFileName);
    
    // Dummy Scenario
    #pragma region
    
    const double normp2 = frobeniusNorm(P2, imaxb, jmaxb);

    boundaryvalues(imax, jmax, domain);
    calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, domain, F2, G2);
    calculate_rs(dt, dx, dy, imax, jmax, F2, G2, RS2);
    
    double res = 1.0;
    int it = 0;
    while(it < itermax && res > eps){
      sor(omg,dx,dy,imax,jmax,domain,RS2,&res);
      it++;
    }

    #pragma endregion
    
    domain.pressure(P2);
    write_matrix_less_precision(cmp_log_P_t, P2, imaxb, jmaxb);

    const double norm_cmpp2 = frobeniusNorm(P2, imaxb, jmaxb);

    const std::string hash_cmpP2 = get_md5hash(cmp_log_P_t);

    #pragma region

    CHECKED_IF (hashP2 == hash_cmpP2)
    {
        std::cout<<hashP2<<" == "<<hash_cmpP2<<std::endl; 
    }
    CHECKED_ELSE (hashP2 == hash_cmpP2) 
    {
        norm_diff_check(normp2, norm_cmpp2, eps, "P");
    }

    #pragma endregion
    
}

TEST_CASE( "Test F & G integration with U&V", "[integrate F&G and calculate_U&V]" )
{
    std::cout << "  " << std::endl;
    std::cout<< "================================================================="<<std::endl;
    std::cout<< "-----------Testing integration of F&G and calculate_uv-----------"<<std::endl;
    std::cout<< "================================================================="<<std::endl;
    std::cout << "  " << std::endl;
    std::string szFileName = "../../cavity100.dat";

    // Dummy Grid Paramters
    #pragma region 
    
    int imax = 10;
    int jmax = 10;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;

    #pragma endregion
    // Initialize Dummy Grid
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI);
    // Dummy Scenario Parameters
    #pragma region 
    double Re = 100;

    double dx = 0.1;
    double dy = 0.1;

    double GX = 0.0;
    double GY = 0.0;

    double alpha = 0.5;
    double dt = 0.0803324;

    int itermax = 100;
    double eps = 0.001;
    double omg = 1.7;

    #pragma endregion


    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();

    // File and directory paths
    std::string logFileName;

    // Data Structure Declaration
    #pragma region 
    // Data structures for first time step
    matrix<double> F1;
    matrix<double> G1;
    matrix<double> RS1;
    matrix<double> P1;
    matrix<double> U1;
    matrix<double> V1;

    // Data structures for second time step
    matrix<double> F2;
    matrix<double> G2;
    matrix<double> RS2;
    matrix<double> P2;
    matrix<double> U2;
    matrix<double> V2;
    #pragma endregion

    // Data Structure Allocation
    #pragma region 
    // Allocating memory for matrices
    F1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    P1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    U1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V1.resize(imaxb,std::vector<double>(jmaxb,0.0));

    // Allocating memory for matrices
    F2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    P2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    U2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion

    // Reading log files in data structures
    #pragma region 
    logFileName = ref_log_F_t + "1";
    read_matrix(logFileName.c_str(),F1,imaxb,jmaxb);
    logFileName = ref_log_G_t + "1";
    read_matrix(logFileName.c_str(),G1,imaxb,jmaxb);
    logFileName = ref_log_RS_t + "1";
    read_matrix(logFileName.c_str(),RS1,imaxb,jmaxb);
    logFileName = ref_log_P_t + "1";
    read_matrix(logFileName.c_str(),P1,imaxb,jmaxb);
    logFileName = ref_log_U_t + "1";
    read_matrix(logFileName.c_str(),U1,imaxb,jmaxb);
    logFileName = ref_log_V_t + "1";
    read_matrix(logFileName.c_str(),V1,imaxb,jmaxb);
    #pragma endregion


    domain.set_velocity(U1, velocity_type::U);
    domain.set_velocity(V1, velocity_type::V);
    domain.set_pressure(P1);


    logFileName = ref_log_U_t + "2";
    const std::string hashU2 = get_md5hash(logFileName);
    logFileName = ref_log_V_t + "2";
    const std::string hashV2 = get_md5hash(logFileName);
    
    // Dummy Scenario
    #pragma region
    
    const double normu2 = frobeniusNorm(U2, imaxb, jmaxb);
    const double normv2 = frobeniusNorm(V2, imaxb, jmaxb);

    boundaryvalues(imax, jmax, domain);
    calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, domain, F2, G2);
    calculate_rs(dt, dx, dy, imax, jmax, F2, G2, RS2);
    
    double res = 1.0;
    int it = 0;
    while(it < itermax && res > eps){
      sor(omg,dx,dy,imax,jmax,domain,RS2,&res);
      it++;
    }
    calculate_uv(dt, dx, dy, imax, jmax, domain, F2, G2);

    domain.pressure(P2);
    domain.velocity(U2, velocity_type::U);
    domain.velocity(V2, velocity_type::V);

    #pragma endregion

    write_matrix_less_precision(cmp_log_U_t, U2, imaxb, jmaxb);
    write_matrix_less_precision(cmp_log_V_t, V2, imaxb, jmaxb);

    const double norm_cmpu2 = frobeniusNorm(U2, imaxb, jmaxb);
    const double norm_cmpv2 = frobeniusNorm(V2, imaxb, jmaxb);

    const std::string hash_cmpU2 = get_md5hash(cmp_log_U_t);
    const std::string hash_cmpV2 = get_md5hash(cmp_log_V_t);

    #pragma region

    CHECKED_IF (hashU2 == hash_cmpU2)
    {
        std::cout<<hashU2<<" == "<<hash_cmpU2<<std::endl; 
    }
    CHECKED_ELSE (hashU2 == hash_cmpU2) 
    {
        norm_diff_check(normu2, norm_cmpu2, eps, "U");
    }

    CHECKED_IF (hashV2 == hash_cmpV2)
    {
        std::cout<<hashV2<<" == "<<hash_cmpV2<<std::endl; 
    }
    CHECKED_ELSE (hashV2 == hash_cmpV2) 
    {
        norm_diff_check(normv2, norm_cmpv2, eps, "V");
    }

    #pragma endregion
    
}

TEST_CASE( "Test whole solver", "[single_iteration_test]" )
{
    std::cout << "  " << std::endl;
    std::cout<< "================================================================="<<std::endl;
    std::cout<< "---------------Testing single iteration of a solver--------------"<<std::endl;
    std::cout<< "================================================================="<<std::endl;
    std::cout << "  " << std::endl;
    std::string szFileName = "../../cavity100.dat";

    // Dummy Grid Paramters
    #pragma region 
    
    int imax = 10;
    int jmax = 10;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;

    #pragma endregion
    // Initialize Dummy Grid
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI);
    // Dummy Scenario Parameters
    #pragma region 
    double Re = 100;

    double dx = 0.1;
    double dy = 0.1;

    double GX = 0.0;
    double GY = 0.0;

    double alpha = 0.5;
    double dt = 0.0803324;

    int itermax = 100;
    double eps = 0.001;
    double omg = 1.7;

    #pragma endregion


    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();

    // File and directory paths
    std::string logFileName;

    // Data Structure Declaration
    #pragma region 
    // Data structures for first time step
    matrix<double> F1;
    matrix<double> G1;
    matrix<double> RS1;
    matrix<double> P1;
    matrix<double> U1;
    matrix<double> V1;

    // Data structures for second time step
    matrix<double> F2;
    matrix<double> G2;
    matrix<double> RS2;
    matrix<double> P2;
    matrix<double> U2;
    matrix<double> V2;
    #pragma endregion

    // Data Structure Allocation
    #pragma region 
    // Allocating memory for matrices
    F1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    P1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    U1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V1.resize(imaxb,std::vector<double>(jmaxb,0.0));

    // Allocating memory for matrices
    F2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    P2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    U2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion

    // Reading log files in data structures
    #pragma region 
    logFileName = ref_log_F_t + "1";
    read_matrix(logFileName.c_str(),F1,imaxb,jmaxb);
    logFileName = ref_log_G_t + "1";
    read_matrix(logFileName.c_str(),G1,imaxb,jmaxb);
    logFileName = ref_log_RS_t + "1";
    read_matrix(logFileName.c_str(),RS1,imaxb,jmaxb);
    logFileName = ref_log_P_t + "1";
    read_matrix(logFileName.c_str(),P1,imaxb,jmaxb);
    logFileName = ref_log_U_t + "1";
    read_matrix(logFileName.c_str(),U1,imaxb,jmaxb);
    logFileName = ref_log_V_t + "1";
    read_matrix(logFileName.c_str(),V1,imaxb,jmaxb);
    #pragma endregion


    domain.set_velocity(U1, velocity_type::U);
    domain.set_velocity(V1, velocity_type::V);
    domain.set_pressure(P1);

    logFileName = ref_log_F_t + "2";
    read_matrix(logFileName.c_str(),F2,imaxb,jmaxb);
    const std::string hashF2 = get_md5hash(logFileName);
    logFileName = ref_log_G_t + "2";
    read_matrix(logFileName.c_str(),G2,imaxb,jmaxb);
    const std::string hashG2 = get_md5hash(logFileName);
    logFileName = ref_log_RS_t + "2";
    read_matrix(logFileName.c_str(),RS2,imaxb,jmaxb);
    const std::string hashRS2 = get_md5hash(logFileName);
    logFileName = ref_log_P_t + "2";
    read_matrix(logFileName.c_str(),P2,imaxb,jmaxb);
    const std::string hashP2 = get_md5hash(logFileName);
    logFileName = ref_log_U_t + "2";
    read_matrix(logFileName.c_str(),U2,imaxb,jmaxb);
    const std::string hashU2 = get_md5hash(logFileName);
    logFileName = ref_log_V_t + "2";
    read_matrix(logFileName.c_str(),V2,imaxb,jmaxb);
    const std::string hashV2 = get_md5hash(logFileName);
    
    const double normf2 = frobeniusNorm(F2, imaxb, jmaxb);
    const double normg2 = frobeniusNorm(G2, imaxb, jmaxb);
    const double normrs2 = frobeniusNorm(RS2, imaxb, jmaxb);
    const double normp2 = frobeniusNorm(P2, imaxb, jmaxb);
    const double normu2 = frobeniusNorm(U2, imaxb, jmaxb);
    const double normv2 = frobeniusNorm(V2, imaxb, jmaxb);

    // Dummy Scenario
    #pragma region
    

    boundaryvalues(imax, jmax, domain);
    calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, domain, F2, G2);
    calculate_rs(dt, dx, dy, imax, jmax, F2, G2, RS2);
    
    double res = 1.0;
    int it = 0;
    while(it < itermax && res > eps){
      sor(omg,dx,dy,imax,jmax,domain,RS2,&res);
      it++;
    }
    calculate_uv(dt, dx, dy, imax, jmax, domain, F2, G2);

    domain.pressure(P2);
    domain.velocity(U2, velocity_type::U);
    domain.velocity(V2, velocity_type::V);

    #pragma endregion

    write_matrix(cmp_log_F_t, F2, imaxb, jmaxb);
    write_matrix(cmp_log_G_t, G2, imaxb, jmaxb);
    write_matrix(cmp_log_RS_t, RS2, imaxb, jmaxb);
    write_matrix(cmp_log_U_t, U2, imaxb, jmaxb);
    write_matrix(cmp_log_V_t, V2, imaxb, jmaxb);
    write_matrix(cmp_log_P_t, P2, imaxb, jmaxb);

    const double norm_cmpf2 = frobeniusNorm(F2, imaxb, jmaxb);
    const double norm_cmpg2 = frobeniusNorm(G2, imaxb, jmaxb);
    const double norm_cmprs2 = frobeniusNorm(RS2, imaxb, jmaxb);
    const double norm_cmpp2 = frobeniusNorm(P2, imaxb, jmaxb);
    const double norm_cmpu2 = frobeniusNorm(U2, imaxb, jmaxb);
    const double norm_cmpv2 = frobeniusNorm(V2, imaxb, jmaxb);

    const std::string hash_cmpF2 = get_md5hash(cmp_log_F_t);
    const std::string hash_cmpG2 = get_md5hash(cmp_log_G_t);
    const std::string hash_cmpRS2 = get_md5hash(cmp_log_RS_t);
    const std::string hash_cmpP2 = get_md5hash(cmp_log_P_t);
    const std::string hash_cmpU2 = get_md5hash(cmp_log_U_t);
    const std::string hash_cmpV2 = get_md5hash(cmp_log_V_t);

    #pragma region

    CHECKED_IF (hashF2 == hash_cmpF2)
    {
        std::cout<<hashF2<<" == "<<hash_cmpF2<<std::endl;   
    }
    CHECKED_ELSE (hashF2 == hash_cmpF2) 
    {
        norm_diff_check(normf2, norm_cmpf2, eps, "F");
    }

    CHECKED_IF (hashG2 == hash_cmpG2)
    {
        std::cout<<hashG2<<" == "<<hash_cmpG2<<std::endl; 
    }
    CHECKED_ELSE (hashG2 == hash_cmpG2) 
    {
        norm_diff_check(normg2, norm_cmpg2, eps, "G");
    }

    CHECKED_IF (hashRS2 == hash_cmpRS2)
    {
        std::cout<<hashRS2<<" == "<<hash_cmpRS2<<std::endl; 
    }
    CHECKED_ELSE (hashRS2 == hash_cmpRS2) 
    {
        norm_diff_check(normrs2, norm_cmprs2, eps, "RS");
    }

    CHECKED_IF (hashP2 == hash_cmpP2)
    {
        std::cout<<hashP2<<" == "<<hash_cmpP2<<std::endl; 
    }
    CHECKED_ELSE (hashP2 == hash_cmpP2) 
    {
        norm_diff_check(normp2, norm_cmpp2, eps, "P");
    }

    CHECKED_IF (hashU2 == hash_cmpU2)
    {
        std::cout<<hashU2<<" == "<<hash_cmpU2<<std::endl; 
    }
    CHECKED_ELSE (hashU2 == hash_cmpU2) 
    {
        norm_diff_check(normu2, norm_cmpu2, eps, "U");
    }

    CHECKED_IF (hashV2 == hash_cmpV2)
    {
        std::cout<<hashV2<<" == "<<hash_cmpV2<<std::endl; 
    }
    CHECKED_ELSE (hashV2 == hash_cmpV2) 
    {
        norm_diff_check(normv2, norm_cmpv2, eps, "V");
    }

    #pragma endregion
    
}



