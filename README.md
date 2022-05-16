# CFD Lab code skeleton

Code Skeleton for the CFD Lab taught at TUM Informatics

This repository contains:

* an example input file (`cavity100.dat`), with `100` stemming from the `itermax`
* the headers
* the files with the respective method stubs

Please [fork this repository and add us as collaborators](https://gitlab.lrz.de/tum-i05/public/cfdlabcodeskeleton/-/wikis/home).

Find more information in the [documentation](https://tum-i05.pages.gitlab.lrz.de/public/cfdlabcodeskeleton/).

## Software Requirements

* VTK 7 or higher
* GCC 9 (optional) 

### GCC version

You can get you current version of GCC by running:

```shell
g++ -v
```

### Defining your GCC version

If you have GCC 9 or newer, you can set in the `CMakeLists.txt` file:

```cmake
set(gpp9 True)
```

If you have a version lower than 9, then you don't have to modify the `CMakeLists.txt` file.

This will affect how we are using the C++ filesystem library, which is available already in GCC 7 as an experimental feature.

### Setup of VTK and GCC 9 (Ubuntu **20.04**)

```
apt-get update &&
apt-get upgrade -y &&
apt-get install -y build-essential cmake libvtk7-dev libfmt-dev
```

### Setup of VTK and GCC 9 (Ubuntu **18.04**)

If you want, you can upgrade your compiler version to have access to more recent C++ features.
This is, however, optional.

```
apt-get update &&
apt-get install -y software-properties-common &&
add-apt-repository -y ppa:ubuntu-toolchain-r/test &&
apt-get upgrade -y &&
apt-get install -y build-essential cmake libvtk7-dev libfmt-dev gcc-9 g++-9
apt-get install -y gcc-9 g++-9
```

## Using CMake

CMake is a C++ build system generator, which simplifies the building process compared e.g. to a system-specific Makefile. The CMake configuration is defined in the `CMakeList.txt` file.

In order to build your code with CMake, you can follow this (quite common) procedure:

1. Create a build directory: `mkdir build`
2. Get inside it: `cd build`
3. Configure and generate the build system: `cmake ..` (Note the two dots, this means that the `CmakeLists.txt` File is in the folder above)
4. Build your code: `make` (build the executable)

### Troubleshooting: VTK not found

You might run into a problem where the VTK library is not found. To fix this, you can try the following steps:

1. Find the installation path of your VTK library 
2. Define this path as an environment variable, as e.g. `export VTK_DIR=".../lib/cmake/vtk-8.2"`
3. Start in a clean build folder
4. Run `cmake ..` again

### Set a different GCC version

If you have multiple compiler versions installed you can set the GCC version which should be used by `cmake` like this:

```shell
export CXX=`which g++-7`
```

Make sure to use a backtick (\`) to get the `which` command executed. Afterwards, you can run `cmake ..`.

## Testing 

The folder `tests` contains example files for simple unit testing with the [Catch2](https://github.com/catchorg/Catch2) unit testing framework. The tests can be compiled using `cmake` as well.

After building, run:

```
ctest --verbose
```

With the `verbose` option you can get more details for failing tests.

## Documentation 

The generate the documentation, run:

```
pip3 install --user mkdocs mkdocs-material
```

and then serve it to be able to access it through your browser:

```
mkdocs serve
```
## Drive

Necessary files containing the results of the simulation with different parameters as well as documents disscussing these results
exist on a google drive please refer to the link: https://drive.google.com/drive/folders/1SrLzJQuSy84o-_T1Ql7MRzbVWiN6-YvY
the drive will have folder for each worksheet, the worksheet folder contains results as vtk file and some log files as txt files

### Execution

To run the simulation follow the shown steps:
1. Clone repository to your local machine
2. Create build directory `mkdir build` then change directory to build `cd build`
3. Run `cmake ..` to create make file inside build directory, the `CMakeLists.txt` in the main project directory
4. Run `make` in build directory and you will find executable created `sim`
5. To Run the simulation a problem has to be defined; running the program `./sim <problem> <optional-flag>`, 
        `<problem>` is a flag to choose example problem
        Current supported problem
            -Lid driven cavity        -c
        `<optional-flags>` flags for developping features(e.g. run on small grid to debug)
        Current supported optional feature
            debug grid                -d
6. If you enter wrong command, you should see help printed for the flags to set

## Folder Structure

To Familiraize with the folder in the repository, there exist 3 types of folders:
1. Input Data
2. VTK Results
3. Log Data for Testing

For Input Data, 2 directories
1. Simulation_Data: These are the dat and pgm files, for the solver with the actual grid in the example problem with full discretization
                    Not modified or created by solver nor tests, only developer modify it
2. Debugging_Data: These are the dat and pgm files, for the solver with the debug grid in the example problem with coarse discretization
                       Not modified or created by solver nor tests, only developer modify it

For VTK Results, 2 directories
1. Ref_Results: These are the vtk files of the simulation as reference after making sure the solver making correct result, 
                Not modified or created by solver nor tests, only developer modify it
2. Temp_Results: These are the vtk files of the simulation that solver writes during simulation,it does not exist on repository, the solver
                 creates this directory if it does not exist, or erase all files in it if it exists to avoid confusion with old results.
                 This directory where you should read vtk files by paraview to check results after each simulation 
                 Created and modified by solver but not tests

For Log Data for Testing, 3 directories
1. Ref_Test_Data: These are the log files of the simulation as reference after making sure the solver making correct result,these log files are
                  used for comparison in testing phase. 
                  Not modified or created by solver nor tests, only developer modify it, even that tests read data from it
2. Temp_Test_Data: These are the log files of the simulation that solver writes during simulation,it does not exist on repository, the solver
                   creates this directory if it does not exist, or erase all files in it if it exists to avoid confusion with old results.
                   Created and modified by solver but not tests
3. Cmp_Test_Data: These are the log files of the simulation that tests writes during testing,it does not exist on repository, the tests
                  creates this directory if it does not exist, or erase all files in it if it exists to avoid confusion with old results.
                  This directory where you should read vtk files by paraview to check results after each simulation 
                  Created and modified by tests but not solver

## Skeleton Modification

The code skeleton did not change, only 2 methods were added to the grid class (grid.hpp grid.cpp), these 2 methods created to return 2D array
instead of matrix containg pressure and velocity. The methods are written to copy data to 2D array required for the function that writes vtk file

void velocity(double** vec, velocity_type type); line 22 grid.hpp
void pressure(double** vec); line 28 grid.hpp

## Worksheet 1 Results

![GitHub Logo](/Final_Results/WS_1/pressure.png)
Format: ![Alt Text](url)