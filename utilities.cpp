#include "utilities.hpp"


static void ParseError(const char *message)
{
      printf("\x1B[31m %s\033[0m\n", message);
      printf("\x1B[37m %s\033[0m\n", help_message.c_str());
      exit(-1);
}


/*----------Command line & directories handling methods----------------*/
#pragma region 
void ParseCommandLineOptions(int argc, char *argv[],std::string &sceFileName,std::string &vtkFileName,std::string &problemName)
{   
    // Error message
    std::string message;
    // Input file directory
    std::string dir = Sim_Data_s;

    // Check if no flags entered
    if(argc == 1){
        message = std::string("ERROR:No example problem chosen\n");
        ParseError(message.c_str());
    }

    if(strcmp(argv[1], "-c") == 0)
    {
        sceFileName = Lid_driven_cavity_s + ".dat";
        vtkFileName = Lid_driven_cavity_s;
        problemName = Lid_driven_cavity_s;
        
    }
    else 
    {
        message = "ERROR:Example problem does not exist";
        ParseError(message.c_str());
    }

    if(argc > 2)
    {
        int i = 2;
        while(i < argc){
            if(strcmp(argv[i], "-d") == 0)
            {
                dir = Deb_Data_s;
            }
            
            i++;
        }
    } 


    // Adding directory to file name
    sceFileName = dir + sceFileName;

    std::ifstream sceTest(sceFileName);
    // Checking if the file exists
    if(sceTest.fail())
    {
        message = "ERROR:file path" + sceFileName + "does not exist";
        ParseError(message.c_str());
    }
}

void clear_output_dir(){
    
    struct stat st = {0};
    
	if (stat(Temp_Results_s.c_str(), &st) == -1)
    {
    	mkdir(Temp_Results_s.c_str(), 0700);
	}
    else
    {
        system(("exec rm -r " + Temp_Results_s + "*").c_str());
    }

    if (stat(Temp_Test_Data_s.c_str(), &st) == -1)
    {
    	mkdir(Temp_Test_Data_s.c_str(), 0700);
	}
    else
    {
        system(("exec rm -r " + Temp_Test_Data_s + "*").c_str());
    }

}
#pragma endregion

/*----------I/P & O/P streams----------------*/
#pragma region 
void read_matrix(std::string fileName,matrix<double>& m,int xdim,int ydim)
{
    std::ifstream fin(fileName);

    for(int y = ydim-1 ; y >= 0 ;y--) {
        for(int x = 0;x<xdim ; x++){                
                fin >> std::fixed >> std::showpos >> std::setprecision(7) >> m.at(x).at(y);
                // fin >> m.at(x).at(y);                        
            }
           
        }
}

void write_matrix(std::string fileName,matrix<double>& m,int xdim,int ydim)
{
    std::ofstream fout(fileName);

    for(int y = ydim-1 ; y >= 0 ;y--)
    {
        for(int x = 0;x<xdim ; x++)
        {           
                fout <<  std::fixed << std::showpos << std::setprecision(15) << m.at(x).at(y) << " ";
                // fout  << m.at(x).at(y) << " ";
                         
        }     
        fout << std::endl;      
    }
}

void write_matrix_less_precision(std::string fileName,matrix<double>& m,int xdim,int ydim)
{
    std::ofstream fout(fileName);

    for(int y = ydim-1 ; y >= 0 ;y--)
    {
        for(int x = 0;x<xdim ; x++)
        {           
                fout <<  std::fixed << std::showpos << std::setprecision(5) << m.at(x).at(y) << " ";
                // fout  << m.at(x).at(y) << " ";
                         
        }     
        fout << std::endl;      
    }
}
#pragma endregion

/*----------Printing data structures----------------*/
#pragma region 
void print_matrix(matrix<double> m,int xdim,int ydim)
{
    for(int y = ydim-1 ; y >= 0 ;y--) {
        for(int x = 0;x<xdim ; x++){
                std::cout <<  std::fixed << std::showpos << std::setprecision(7) << m.at(x).at(y) << " ";            
            }
            // Print new line
            std::cout << std::endl;
        }
}
#pragma endregion

/*computing the hash for log files*/
#pragma region 
std::string get_md5hash( const std::string fname)
{
    std::ifstream mySource;
    mySource.open(fname, std::ios_base::binary);
    mySource.seekg(0, std::ios_base::end);
    int BUFFSIZE = mySource.tellg();
    mySource.close();

    unsigned char digest[MD5_DIGEST_LENGTH];

    std::stringstream ss;
    std::string md5string;

    char buffer[BUFFSIZE];

    std::ifstream ifs(fname, std::ios::in | std::ifstream::binary);
    // FILE * pFile;
    // char * buffer;
    // long lSize;
    // size_t result;

    // pFile = fopen ( fname.c_str() , "rb" );
    // if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
    // fseek (pFile , 0 , SEEK_END);
    // lSize = ftell (pFile);
    // rewind (pFile);

    // // allocate memory to contain the whole file:
    // buffer = (char*) malloc (sizeof(char)*lSize);
    // if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

    MD5_CTX md5Context;

    MD5_Init(&md5Context);

    // copy the file into the buffer:
    // result = fread (buffer,1,lSize,pFile);
    // if (result != lSize) {fputs ("Reading error",stderr); exit (3);}
    // MD5_Update(&md5Context, &buffer, lSize); 
    //std::cout<<buffer<<std::endl;


    // std::stringstream oss;
    //oss << ifs.rdbuf();

    //const std::string& file_str = oss.str();
    //const char* cstr = file_str.c_str();
    while (ifs.good()) 
    {
        ifs.read(buffer, BUFFSIZE);
        MD5_Update(&md5Context, &buffer, ifs.gcount()); 
    }
    
    if (ifs)
    { 
        std::cout<< fname <<" couldn't be opened\n"; 
        exit(1); 
    }
    // std::cout<<"for file"<<fname<<"\n"<<cstr<<"\n";
    
    // oss.seekg(0, std::ios::end);
    // int size = oss.tellg();

    ifs.close();
    // fclose (pFile);
    // free (buffer);

    int res = MD5_Final(digest, &md5Context);

    if( res == 0 ) // hash failed
        return {};   // or raise an exception

    // set up stringstream format
    ss << std::hex << std::uppercase << std::setfill('0');


    for(unsigned char uc: digest)
    ss << std::setw(2) << (int)uc;


    md5string = ss.str();

    return md5string;
}
#pragma endregion
