#define CATCH_CONFIG_RUNNER

#include "catch.hpp"
#include "../init.hpp"
#include "../grid.hpp"
#include "../utilities.hpp"
#include "../boundary_val.hpp"
#include "../uvp.hpp"
#include "../sor.hpp"
#include "test_utilities.hpp"
#include "../parallel.hpp"

int main( int argc, char* argv[] ) {
  clear_output_dir_test();

  MPI_Init(&argc, &argv);
  int result = Catch::Session().run( argc, argv );
  MPI_Finalize();

  return result;
}

SCENARIO("Parallel Environment Testing", "[init_parallel]") {
    
    int num_proc = 0;
    int myrank = 0;


    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    #pragma region

    // Parallel data
    int iproc = 0;
    int jproc = 0;
    int il = 0;
    int ir = 0;
    int jb = 0;
    int jt = 0;
    int rank_l = 0;
    int rank_r = 0;
    int rank_b = 0;
    int rank_t = 0;
    int omg_i = 0;
    int omg_j = 0;

    #pragma endregion

    int imax = 300;
    int jmax = 300;
    iproc = 2;
    jproc = 1;
  

    init_parallel(iproc,jproc,imax,jmax,&myrank,&il,&ir,&jb,&jt,&rank_l,
                  &rank_r,&rank_b,&rank_t,&omg_i,&omg_j,num_proc);

    REQUIRE ( num_proc == (iproc*jproc) );
    if (myrank == 0)
    {
        CHECK ( il == 1 );
        CHECK ( ir == 150 );
        CHECK ( jb == 1 );
        CHECK ( jt == 300 );
        CHECK ( omg_i == 1 );
        CHECK ( omg_j == 1 );
        CHECK ( rank_l == -2 );
        CHECK ( rank_r == 1 );
        CHECK ( rank_t == -2 );
        CHECK ( rank_b == -2 );

    }
    if (myrank == 1)
    {
        CHECK ( il == 151 );
        CHECK ( ir == 300 );
        CHECK ( jb == 1 );
        CHECK ( jt == 300 );
        CHECK ( omg_i == 2 );
        CHECK ( omg_j == 1 );
        CHECK ( rank_l == 0 );
        CHECK ( rank_r == -2 );
        CHECK ( rank_t == -2 );
        CHECK ( rank_b == -2 );
    }
    // if (myrank == 2)
    // {
    //     CHECK ( il == 1 );
    //     CHECK ( ir == 150 );
    //     CHECK ( jb == 101 );
    //     CHECK ( jt == 200 );
    //     CHECK ( omg_i == 1 );
    //     CHECK ( omg_j == 2 );
    //     CHECK ( rank_l == -2 );
    //     CHECK ( rank_r == 3 );
    //     CHECK ( rank_t == 4 );
    //     CHECK ( rank_b == 0 );
    // }
    // if (myrank == 3)
    // {
    //     CHECK ( il == 151 );
    //     CHECK ( ir == 300 );
    //     CHECK ( jb == 101 );
    //     CHECK ( jt == 200 );
    //     CHECK ( omg_i == 2 );
    //     CHECK ( omg_j == 2 );
    //     CHECK ( rank_l == 2 );
    //     CHECK ( rank_r == -2 );
    //     CHECK ( rank_t == 5 );
    //     CHECK ( rank_b == 1 );
    // }
    // if (myrank == 4)
    // {
    //     CHECK ( il == 1 );
    //     CHECK ( ir == 150 );
    //     CHECK ( jb == 201 );
    //     CHECK ( jt == 300 );
    //     CHECK ( omg_i == 1 );
    //     CHECK ( omg_j == 3 );
    //     CHECK ( rank_l == -2 );
    //     CHECK ( rank_r == 5 );
    //     CHECK ( rank_t == -2 );
    //     CHECK ( rank_b == 2 );
    // }
    // if (myrank == 5)
    // {
    //     CHECK ( il == 151 );
    //     CHECK ( ir == 300 );
    //     CHECK ( jb == 201 );
    //     CHECK ( jt == 300 );
    //     CHECK ( omg_i == 2 );
    //     CHECK ( omg_j == 3 );
    //     CHECK ( rank_l == 4 );
    //     CHECK ( rank_r == -2 );
    //     CHECK ( rank_t == -2 );
    //     CHECK ( rank_b == 3 );
    // }

}

SCENARIO("U&V Communication Testing", "[uv_comm]") {

    int num_proc = 0;
    int myrank = 0;


    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    #pragma region

    // Parallel data
    int iproc = 0;
    int jproc = 0;
    int il = 0;
    int ir = 0;
    int jb = 0;
    int jt = 0;
    int rank_l = 0;
    int rank_r = 0;
    int rank_b = 0;
    int rank_t = 0;
    int omg_i = 0;
    int omg_j = 0;

    matrix<double> U;
    matrix<double> V;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;

    #pragma endregion

    int imax = 300;
    int jmax = 300;
    iproc = 2;
    jproc = 1;

    init_parallel(iproc,jproc,imax,jmax,&myrank,&il,&ir,&jb,&jt,&rank_l,
                  &rank_r,&rank_b,&rank_t,&omg_i,&omg_j,num_proc);

    imax = ir-il+1;
    jmax = jt-jb+1;

    Grid domain = Grid(imax, jmax, boundary_size_mpi, PI, UI, VI) ;

    U.resize(domain.imaxb(), std::vector<double>(domain.jmaxb(), 0.0));
    V.resize(domain.imaxb(), std::vector<double>(domain.jmaxb(), 0.0));

    int buffsize = (domain.imaxb())>(domain.jmaxb())?(domain.imaxb()):(domain.jmaxb());
    double *bufSend = (double*)malloc(sizeof(double)*2*buffsize);
    double *bufRecv = (double*)malloc(sizeof(double)*2*buffsize);

    std::string logFileName;
    logFileName = ref_log_U_t + std::to_string(myrank) + "_before";
    read_matrix(logFileName, U, domain.imaxb(), domain.jmaxb());
    logFileName = ref_log_V_t + std::to_string(myrank) + "_before";
    read_matrix(logFileName, V, domain.imaxb(), domain.jmaxb());

    uv_comm(domain,imax,jmax,rank_l,rank_r,rank_b,rank_t,myrank,bufSend,bufRecv,U,V);

    logFileName = cmp_log_U_t + std::to_string(myrank);
    write_matrix(logFileName, U, domain.imaxb(), domain.jmaxb());
    logFileName = cmp_log_V_t + std::to_string(myrank);
    write_matrix(logFileName, V, domain.imaxb(), domain.jmaxb());

    if (myrank == 0 )
    {
        logFileName = ref_log_U_t + std::to_string(myrank) + "_after";
        const std::string u0_hash = get_md5hash(logFileName);
        logFileName = ref_log_V_t + std::to_string(myrank) + "_after";
        const std::string v0_hash = get_md5hash(logFileName);

        logFileName = cmp_log_U_t + std::to_string(myrank);
        const std::string u0_hash_cmp = get_md5hash(logFileName);
        logFileName = cmp_log_V_t + std::to_string(myrank);
        const std::string v0_hash_cmp = get_md5hash(logFileName);

        CHECK ( u0_hash == u0_hash_cmp );
        CHECK ( v0_hash == v0_hash_cmp );
    }

    if (myrank == 1 )
    {
        logFileName = ref_log_U_t + std::to_string(myrank) + "_after";
        const std::string u1_hash = get_md5hash(logFileName);
        logFileName = ref_log_V_t + std::to_string(myrank) + "_after";
        const std::string v1_hash = get_md5hash(logFileName);

        logFileName = cmp_log_U_t + std::to_string(myrank);
        const std::string u1_hash_cmp = get_md5hash(logFileName);
        logFileName = cmp_log_V_t + std::to_string(myrank);
        const std::string v1_hash_cmp = get_md5hash(logFileName);

        CHECK ( u1_hash == u1_hash_cmp );
        CHECK ( v1_hash == v1_hash_cmp );
    }

    // if (myrank == 2 )
    // {
    //     logFileName = ref_log_U_t + std::to_string(myrank) + "_after";
    //     const std::string u2_hash = get_md5hash(logFileName);
    //     logFileName = ref_log_V_t + std::to_string(myrank) + "_after";
    //     const std::string v2_hash = get_md5hash(logFileName);

    //     logFileName = cmp_log_U_t + std::to_string(myrank);
    //     const std::string u2_hash_cmp = get_md5hash(logFileName);
    //     logFileName = cmp_log_V_t + std::to_string(myrank);
    //     const std::string v2_hash_cmp = get_md5hash(logFileName);

    //     CHECK ( u2_hash == u2_hash_cmp );
    //     CHECK ( v2_hash == v2_hash_cmp );
    // }

    // if (myrank == 3 )
    // {
    //     logFileName = ref_log_U_t + std::to_string(myrank) + "_after";
    //     const std::string u3_hash = get_md5hash(logFileName);
    //     logFileName = ref_log_V_t + std::to_string(myrank) + "_after";
    //     const std::string v3_hash = get_md5hash(logFileName);

    //     logFileName = cmp_log_U_t + std::to_string(myrank);
    //     const std::string u3_hash_cmp = get_md5hash(logFileName);
    //     logFileName = cmp_log_V_t + std::to_string(myrank);
    //     const std::string v3_hash_cmp = get_md5hash(logFileName);

    //     CHECK ( u3_hash == u3_hash_cmp );
    //     CHECK ( v3_hash == v3_hash_cmp );
    // }

    // if (myrank == 4 )
    // {
    //     logFileName = ref_log_U_t + std::to_string(myrank) + "_after";
    //     const std::string u4_hash = get_md5hash(logFileName);
    //     logFileName = ref_log_V_t + std::to_string(myrank) + "_after";
    //     const std::string v4_hash = get_md5hash(logFileName);

    //     logFileName = cmp_log_U_t + std::to_string(myrank);
    //     const std::string u4_hash_cmp = get_md5hash(logFileName);
    //     logFileName = cmp_log_V_t + std::to_string(myrank);
    //     const std::string v4_hash_cmp = get_md5hash(logFileName);

    //     CHECK ( u4_hash == u4_hash_cmp );
    //     CHECK ( v4_hash == v4_hash_cmp );
    // }

    // if (myrank == 5 )
    // {
    //     logFileName = ref_log_U_t + std::to_string(myrank) + "_after";
    //     const std::string u5_hash = get_md5hash(logFileName);
    //     logFileName = ref_log_V_t + std::to_string(myrank) + "_after";
    //     const std::string v5_hash = get_md5hash(logFileName);

    //     logFileName = cmp_log_U_t + std::to_string(myrank);
    //     const std::string u5_hash_cmp = get_md5hash(logFileName);
    //     logFileName = cmp_log_V_t + std::to_string(myrank);
    //     const std::string v5_hash_cmp = get_md5hash(logFileName);

    //     CHECK ( u5_hash == u5_hash_cmp );
    //     CHECK ( v5_hash == v5_hash_cmp );
    // }
}


SCENARIO("P Communication Testing", "[pressure_comm]")
{
    int num_proc = 0;
    int myrank = 0;

    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    #pragma region

    // Parallel data
    int iproc = 0;
    int jproc = 0;
    int il = 0;
    int ir = 0;
    int jb = 0;
    int jt = 0;
    int rank_l = 0;
    int rank_r = 0;
    int rank_b = 0;
    int rank_t = 0;
    int omg_i = 0;
    int omg_j = 0;

    matrix<double> P;
    
    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;

    #pragma endregion

    int imax = 300;
    int jmax = 300;
    iproc = 2;
    jproc = 1;

    init_parallel(iproc,jproc,imax,jmax,&myrank,&il,&ir,&jb,&jt,&rank_l,
                  &rank_r,&rank_b,&rank_t,&omg_i,&omg_j,num_proc);
    
    imax = ir-il+1;
    jmax = jt-jb+1;

    Grid domain = Grid(imax, jmax, boundary_size_mpi, PI, UI, VI);

    P.resize(domain.imaxb(),std::vector<double>(domain.jmaxb(),PI));

    int buffsize = (domain.imaxb())>(domain.jmaxb())?(domain.imaxb()):(domain.jmaxb());
    double *bufSend = (double*)malloc(sizeof(double)*2*buffsize);
    double *bufRecv = (double*)malloc(sizeof(double)*2*buffsize);

    std::string logFileName;
    logFileName = ref_log_P_t + std::to_string(myrank) + "_before";
    read_matrix(logFileName, P, domain.imaxb(), domain.jmaxb());

    pressure_comm(domain,imax,jmax,rank_l,rank_r,rank_b,rank_t,bufSend,bufRecv,P);

    logFileName = cmp_log_P_t + std::to_string(myrank);
    write_matrix(logFileName, P, domain.imaxb(), domain.jmaxb());

    if (myrank == 0)
    {
        logFileName = ref_log_P_t + std::to_string(myrank) + "_after";
        const std::string p0_hash = get_md5hash(logFileName);

        logFileName = cmp_log_P_t + std::to_string(myrank);
        const std::string p0_hash_cmp = get_md5hash(logFileName);

        CHECK (p0_hash == p0_hash_cmp);
    }

    if (myrank == 1)
    {
        logFileName = ref_log_P_t + std::to_string(myrank) + "_after";
        const std::string p1_hash = get_md5hash(logFileName);

        logFileName = cmp_log_P_t + std::to_string(myrank);
        const std::string p1_hash_cmp = get_md5hash(logFileName);

        CHECK (p1_hash == p1_hash_cmp);
    }

    // if (myrank == 2)
    // {
    //     logFileName = ref_log_P_t + std::to_string(myrank) + "_after";
    //     const std::string p2_hash = get_md5hash(logFileName);

    //     logFileName = cmp_log_P_t + std::to_string(myrank);
    //     const std::string p2_hash_cmp = get_md5hash(logFileName);

    //     CHECK (p2_hash == p2_hash_cmp);
    // }

    // if (myrank == 3)
    // {
    //     logFileName = ref_log_P_t + std::to_string(myrank) + "_after";
    //     const std::string p3_hash = get_md5hash(logFileName);

    //     logFileName = cmp_log_P_t + std::to_string(myrank);
    //     const std::string p3_hash_cmp = get_md5hash(logFileName);

    //     CHECK (p3_hash == p3_hash_cmp);
    // }

    // if (myrank == 4)
    // {
    //     logFileName = ref_log_P_t + std::to_string(myrank) + "_after";
    //     const std::string p4_hash = get_md5hash(logFileName);

    //     logFileName = cmp_log_P_t + std::to_string(myrank);
    //     const std::string p4_hash_cmp = get_md5hash(logFileName);

    //     CHECK (p4_hash == p4_hash_cmp);
    // }

    // if (myrank == 5)
    // {
    //     logFileName = ref_log_P_t + std::to_string(myrank) + "_after";
    //     const std::string p5_hash = get_md5hash(logFileName);

    //     logFileName = cmp_log_P_t + std::to_string(myrank);
    //     const std::string p5_hash_cmp = get_md5hash(logFileName);

    //     CHECK (p5_hash == p5_hash_cmp);
    // }


}