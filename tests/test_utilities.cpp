#include "test_utilities.hpp"

/*----------Output Directories Handling----------------*/
#pragma region 
void clear_output_dir_test()
{
    struct stat st = {0};

    if (stat(Cmp_Test_Data_t.c_str(), &st) == -1)
    {
    	mkdir(Cmp_Test_Data_t.c_str(), 0700);
	}
    else
    {
        system(("exec rm -r " + Cmp_Test_Data_t + "*").c_str());
    }
}
#pragma endregion

double frobeniusNorm(matrix<double> mat, int imax, int jmax) 
{ 
  
    // To store the sum of squares of the 
    // elements of the given matrix 
    int sumSq = 0; 
    for (int i = 0; i < imax; i++) { 
        for (int j = 0; j < jmax; j++) { 
            sumSq += pow(mat[i][j], 2); 
        } 
    } 
  
    // Return the square root of 
    // the sum of squares 
    float res = sqrt(sumSq); 
    return res; 
} 

void norm_diff_check(const double norm_given, const double norm_calc, 
                     const double thresh, const std::string mat)
{
    
    //normcheck
    if (abs(norm_given - norm_calc) < thresh)
    {
        std::cout<<"Norm difference is negligible between computed and given "<< mat <<std::endl;
    }
    else
    {
        std::cout<<"Norm difference is greater then threshold between given and calculated "<< mat <<" check respective unit"<<std::endl;
    }
    

}