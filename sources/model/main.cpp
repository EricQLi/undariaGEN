// -------------------------------------------------------------------------- //
// UndariaGEN - A spatially-explicit agent-based model of invasive algal 
//            species in the marine environment
// Copyright (C) 2005-2015  James Thomas Murphy, Ray Walshe
//
// This file is a part of UndariaGEN
//
// UndariaGEN is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// UndariaGEN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// -------------------------------------------------------------------------- //


// -------------------------------------------------------------------------- //
// TO RUN PROGRAM:
//    ./bin/coastgen ./input/inFile1.in inFile2.in inFile3.in inFile4.in
//    inFile1.in = input parameters
//    inFile2.in = map of substrate types (imageJ xy coords)
//    (optional) inFile3.in = map of coastal elevations (imageJ pixel values)
//    (optional) inFile4.in = map of Undaria locations (imageJ xy coords)
// e.g. to run (x4 processors):
//   - mpiexec -n 4 ./bin/coastgen ./input/test_coastgen.in ./input/port_edit514_482.in
// -------------------------------------------------------------------------- //

// main.cpp

// -------------------------------------------------------------------------- //
// Libraries                                                                  //
// -------------------------------------------------------------------------- //

#include "mpi.h"
#include "../engine/myException.hpp"
#include "bait.hpp"
#include "parameters.hpp"
#include <exception>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <string>
#include <sstream>
using namespace std;


// JM - Functions for inputing parameters
static void inputDefaultParams(Parameters *);
static void inputFileParams(Parameters *, double *, int);
static char* nameOutputFile(time_t, int, string, string);
static char* nameOutputFileLB(time_t, int, string, string);
static void endTimekeeping(Bait *, time_t, int);
static void inputMapData(int *, int *, ifstream *, ifstream *, Parameters *);
static void inputMapDataMPI(int *, ifstream *, Parameters *, int, int);
static void inputAgentLocations(int *, ifstream *, Parameters *);
static void inputAgentLocationsMPI(int *, ifstream *, Parameters *, int, int);
static void inputSubstrateMap(int *, ifstream *, Parameters *);
static void inputSubstrateMapMPI(int *, ifstream *, Parameters *, int, int);



// -------------------------------------------------------------------------- //
// Main function                                                              //
// -------------------------------------------------------------------------- //

int main( int argc, char **argv )
{
   try
   {
      Bait               simulation      ;
      Parameters         params          ;
      const int numParameters = 48;
      int myRank;
      int mpiSize;
      int numOutFilesLB = 10, numOutFilesLB2=6, numOutFilesMG=10;
      ifstream inputFile, inputFile2, inputFile3, inputFile4, inputFile5, inputFile6;
      ofstream **outputFilesLB=NULL;
      ofstream **outputFilesLB2=NULL;
      ofstream **outputFilesMG=NULL;
      ofstream outputFile;
      ofstream outputFile2;
      int *xy_elevations;
      int *agentLocations;
      int *substrateLocations;

      MPI_Init(&argc, &argv);
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);   
      MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

      // JM - timekeeping code
      time_t start = time(NULL);

      // Processor 0 deals with time-keeping.
      if (myRank == 0) 
         cout << "Start time = " << start << endl;
         

      //************************************************************************//
      
      //******************** 1. Create Output Files ****************************//
      
      // Processor 0 writes results output files
      if (myRank == 0)
      {
         outputFilesLB = new ofstream*[numOutFilesLB];
         outputFilesLB2 = new ofstream*[numOutFilesLB2];
         outputFilesMG = new ofstream*[numOutFilesMG];
      
         // Output file "mgRes<date>_<myRank>.out" created.
         char *filenamePtr = nameOutputFile(start, myRank, argv[1], "");
         char *filenamePtr2 = nameOutputFile(start, myRank, argv[1], "_diag1");

         outputFile.open(filenamePtr, ios::out);
         outputFile2.open(filenamePtr2, ios::out);
        
         // create array of output files for Lattice Boltzmann results 
         // (see Plate::configLBE()) 
         for (int i=0; i<numOutFilesLB; i++)
         {
            string s;
            std::stringstream out;
            out << (60+i);
            s = out.str();
            char *filenamePtrLB = nameOutputFileLB( start, myRank, argv[1], s);
            outputFilesLB[i] = new ofstream(filenamePtrLB, ios::out);
         }
         // create array of output files for Lattice Boltzmann results (see Plate::profileLBE()) 
         for (int i=1; i<=numOutFilesLB2; i++)
         {
            string s;
            std::stringstream out;
            out << i;
            s = out.str();
            char *filenamePtrLB = nameOutputFileLB( start, myRank, argv[1], s);
            outputFilesLB2[i-1] = new ofstream(filenamePtrLB, ios::out);
         }
         // create array of output files for GNUPlot map of nutrient levels
         for (int i=0; i<numOutFilesMG; i++)
         {
            string s;
            std::stringstream out;
            out << (70+i);
            s = out.str();
            char *filenamePtrMG = nameOutputFileLB( start, myRank, argv[1], s);
            outputFilesMG[i] = new ofstream(filenamePtrMG, ios::out);
         }  


         if (!outputFile || !outputFile2)
         {
            cout << "Could not open output files " << endl;
            //system ("PAUSE");
            return -1;
         }
         
      }  // end if (myRank == 0)
      
      
      //************************************************************************//


      // ********************* 2. Process Input Files **************************//

      // JM - code for receiving input file as command-line argument
      // Input file contains the simulation parameters. 
      
      // If no input file specified then exit.
      if (argc < 2)
      {
         cout << "No input files specified - exiting.\n" << endl;
         //inputDefaultParams(&params);
         return -1;
      }
     
      cout << "argc = " << argc  << ", " << argv[1] << endl;
      inputFile.open(argv[1], ios::in);
      if (!inputFile)
      {
         cout << "Could not open input file 1: " << "\"" << argv[1] 
            << "\"" << endl;
         //system ("PAUSE");
         return -1;
      }
      
      // Input file with wind rose data
      if (argc > 2)
      {
         inputFile2.open(argv[2], ios::in);
         
         if (!inputFile2)
         {
            cout << "Could not open input file 2 (substrate map): " << "\"" << argv[2] 
            << "\"" << endl;
            //system ("PAUSE");
            return -1;
         }
      }
      /*   
      // Input files with map layer data (elevation data and agent locations)
      if (argc > 4)
      {
         inputFile5.open(argv[3], ios::in);
         inputFile6.open(argv[5], ios::in);
         
         if (!inputFile5 || !inputFile6)
         {
            cout << "Could not open input file 5 (elevations) and/or 6 (agent locations): " << "\"" << argv[5] 
               << "\", \"" << argv[6] << "\"" << endl;
            //system ("PAUSE");
            return -1;
         }         
      }*/
      
      // Input file with
      
      // ********************** GENERAL INPUT PARAMETERS **********************//

      // Enter parameter values from input file into paramValues[] array
      // and prints out on screen for error checking purposes
      string* paramNames = new string[numParameters];
      double* paramValues = new double[numParameters];

      if (myRank == 0)         
         cout << "Parameter Name" << "\t\t" << "Value" << endl;

      for (int i=0; i<numParameters; i++)
      {
         inputFile >> paramNames[i];
         while (paramNames[i].substr(0,2) == "//")
         {
            inputFile >> paramNames[i];
         }
         inputFile >> paramValues[i];
         if (myRank == 0)
         {
            cout << paramNames[i] << "\t\t" << paramValues[i] << endl;
         }
      }

      inputFileParams(&params, paramValues, mpiSize);

      // return memory to heap
      delete [] paramValues;
      delete [] paramNames;
      paramValues=0;
      paramNames=0;
           
      
      // ************************ INPUT SUBSTRATE MAP *************************//
      if (myRank == 0)
         cout << "Inputting Substrate Map " << endl;
      
      int arraySize = params.engParams_.width_ * params.engParams_.height_;
      substrateLocations = new int[arraySize];
      if (argc == 3)
      {
         if (mpiSize == 1)
            inputSubstrateMap(substrateLocations, &inputFile2, &params);
         else
            inputSubstrateMapMPI(substrateLocations, &inputFile2, &params, myRank, mpiSize);
      }
      else
      {
         for (int i=0; i<arraySize; i++)
            substrateLocations[i] = 0;
      }
      
      // ********************* ELEVATIONS/AGENT LOCATIONS *********************//
      
      // Input maps of elevation and agent location data from input files 
      xy_elevations = new int[arraySize];
      agentLocations = new int[arraySize];
      
      if (myRank == 0)
         cout << "Inputting Elevations " << endl;
      
      if (argc == 4)
      {
         if (mpiSize == 1)
         {
            inputMapData(xy_elevations, agentLocations, &inputFile2, &inputFile5, &params);
            inputAgentLocations(agentLocations, &inputFile6, &params);
         }
         else
         {
            inputMapDataMPI(xy_elevations, &inputFile2, &params, myRank, mpiSize);
            inputAgentLocationsMPI(agentLocations, &inputFile6, &params, myRank, mpiSize);
         }
      }
      else
      {
         for (int i=0; i<arraySize; i++)
         {
            xy_elevations[i] = 150;
            agentLocations[i] = 0;
         }
      }
         
      
      // Output size of agents for review purposes
      if (myRank == 0)
      {
         cout << "\nSize of Organism = " << sizeof(Organism) << " bytes\n";
         cout << "Size of Patch    = " << sizeof(Patch) << " bytes\n";
      }
    
      if (myRank == 0)
         cout << "\nInitializing... " << endl;
      //************************************************************************//
      
      //******************* 3. Initialize the System ***************************//

      if (myRank == 0)
      {
         simulation.init(params, &outputFile, &outputFile2, 
                                 outputFilesLB, numOutFilesLB, 
                                 outputFilesLB2, numOutFilesLB2, 
                                 outputFilesMG, numOutFilesMG, 
                                 xy_elevations, agentLocations,
                                 substrateLocations);
      }
      else
      {
         simulation.init(params, NULL, NULL, 
                           outputFilesLB, numOutFilesLB, 
                           outputFilesLB2, numOutFilesLB2, 
                           outputFilesMG, numOutFilesMG, 
                           xy_elevations, agentLocations,
                           substrateLocations);
      }


      
      //************************************************************************//
  
      //******************** 4. Run the simulation *****************************//
      if (myRank == 0)
         cout << "Go... " << endl;
      simulation.run();
     
           
      //************************************************************************//
      
      //*********************** 5. Close Down **********************************//
      // Clean up object pointers
      delete [] xy_elevations; delete [] agentLocations;
      xy_elevations=0; agentLocations = 0;
   
      if (myRank == 0)
      {
         outputFile.close();
         outputFile2.close();     
         for (int i=0; i<numOutFilesLB; i++)
         {
             outputFilesLB[i]->close();
             delete outputFilesLB[i];
         }
         for (int i=0; i<numOutFilesLB2; i++)
         {
             outputFilesLB2[i]->close();
             delete outputFilesLB2[i];
          }
         for (int i=0; i<numOutFilesMG; i++)
         {
             outputFilesMG[i]->close();
             delete outputFilesMG[i];
          }
          
          delete [] outputFilesLB; delete [] outputFilesLB2; delete [] outputFilesMG;
      }
     

      endTimekeeping(&simulation, start, myRank);
      MPI_Finalize();
   }
   catch( exception &error )
   {
      cerr << "Exception occurred: " << error.what() << endl ;
      return 1 ;
   }

   return 0 ;
}



//
// JM - Function inputs default parameters for program.
//

static void inputDefaultParams(Parameters *params)
{
   // Engine parameters
   params->engParams_.width_                  =    200 ;
   params->engParams_.height_                 =    200 ;
   params->engParams_.cellSize                =    1.0;
   params->engParams_.nbLoops_                =    8000 ;
   params->engParams_.displayPeriod_          =    100 ;
   params->engParams_.graphPeriod             =    730 ;
   params->engParams_.sampleLoop              =    365 ;
   params->engParams_.worldType_              =    OPEN   ;
   params->engParams_.nbRuns                  =    1;
   params->engParams_.mobility                =    0;  
   params->baitParams_.CUDALoops              =    2;

   // Model parameters
   params->baitParams_.initialOrganisms_      =    10  ;
   params->baitParams_.maxOrganisms_          =    10000 ;
   params->baitParams_.startStock             =    7.3E-04;
   params->baitParams_.probFertilise          =    2.0E-04;
   params->baitParams_.macroSize              =    0.2;
   params->baitParams_.meanDepth              =    1.0;
   params->baitParams_.photo_alpha            =    1.0;
   
   params->baitParams_.sporeHalfLife          =    24.0;
   params->baitParams_.sporeProd              =    0.25E+09;
   params->baitParams_.totalSpore             =    1.0E+10;
   params->baitParams_.probGerminate          =    1.0E-12;
   params->baitParams_.diffusionRate          =    0.0;
   params->baitParams_.diffusLoops            =    1;
   
   params->baitParams_.cos_units              =    365;
   params->baitParams_.tempAmp                =    4.14;
   params->baitParams_.temp_c                 =    229.99;
   params->baitParams_.temp_d                 =    13.33;
   params->baitParams_.solarAmp               =    9.05;
   params->baitParams_.solar_c                =    169.59;
   params->baitParams_.solar_d                =    11.36;
   params->baitParams_.DL_Amp                 =    3.77;
   params->baitParams_.DL_c                   =    171.135;
   params->baitParams_.DL_d                   =    12.25;
   params->baitParams_.kdPAR                  =    0.4;

   // Lattice Boltzmann parameters
   params->baitParams_.initRho                =    1.0;
   params->baitParams_.initUvelocity          =    0.05;
   params->baitParams_.initVvelocity          =    0.0;
   params->baitParams_.omega                  =    1.8;
   params->baitParams_.uForce                 =    0.9;
   params->baitParams_.obstXPos               =    50;
   params->baitParams_.obstYPos               =    100;
   params->baitParams_.obstLength             =    10;
   params->baitParams_.obstWidth              =    1;
   params->baitParams_.loopAddObst            =    -1;       // if -1, no obstacle
   params->baitParams_.loopRemoveObst         =    -1;       // if -1, don't remove
   params->baitParams_.rankAddObst            =    -1;
   params->baitParams_.ratio_hH               =    1.0;

   params->baitParams_.displayOrganismsNumber_=    true ;
   params->baitParams_.initDistribute       =    0;

   // Display parameters
   /// Probleme avec ca...
   params->xParams_.nbPatchColors_            =    32 ;
   params->xParams_.nbOrganismColors_         =    32 ;
   params->xParams_.scale_                    =    1 ;
}



//
// JM - Function for inputting program parameters from array paramValues
//

static void inputFileParams(Parameters *params, double *paramValues, int mpiSize)
{
   int i=0;
   // Engine parameters
   if (mpiSize > 1) 
   {  
      unsigned mpiDim = sqrt(mpiSize);
      unsigned totWidth = static_cast<unsigned>(paramValues[i++]);
      unsigned totHeight = static_cast<unsigned>(paramValues[i++]);
      unsigned bufferZone = ((mpiDim-2)*2)+2;
      
      // Error check: to make sure that environment is evnly dividable by mpiDim
      float width = (float)(totWidth+bufferZone)/mpiDim;
      float height = (float)(totHeight+bufferZone)/mpiDim;
      int testWidth = (int)(width*10);
      int testHeight = (int)(height*10);
      if (testWidth % 10 > 0 || testHeight % 10 > 0)
         cout << "Warning: Width/Height Cannot be divided equally by " << mpiDim << endl;
      
      params->engParams_.width_          =     (totWidth+bufferZone)/mpiDim;
      params->engParams_.height_         =     (totHeight+bufferZone)/mpiDim;
   }
   else
   {
      params->engParams_.width_          =    static_cast<unsigned>(paramValues[i++]);
      params->engParams_.height_         =    static_cast<unsigned>(paramValues[i++]);
   }
   params->engParams_.cellSize           =    paramValues[i++];
   params->engParams_.nbLoops_           =    static_cast<unsigned>(paramValues[i++]);
   params->engParams_.displayPeriod_     =    static_cast<unsigned>(paramValues[i++]);
   params->engParams_.graphPeriod        =    static_cast<unsigned>(paramValues[i++]);
   params->engParams_.sampleLoop         =    static_cast<unsigned>(paramValues[i++]);
   params->engParams_.worldType_         =    (WorldType_t)(static_cast<unsigned>(paramValues[i++])-1);
   params->engParams_.nbRuns             =    static_cast<unsigned>(paramValues[i++]);
   params->engParams_.mobility           =    static_cast<bool>(paramValues[i++]); 
   params->baitParams_.CUDALoops         =    static_cast<unsigned>(paramValues[i++]);

   // Bait parameters
   params->baitParams_.initialOrganisms_  =   static_cast<unsigned>(paramValues[i++]);
   params->baitParams_.maxOrganisms_      =   static_cast<unsigned>(paramValues[i++]);
   params->baitParams_.startStock         =   paramValues[i++];
   params->baitParams_.probFertilise      =   paramValues[i++];
   params->baitParams_.macroSize          =   paramValues[i++];
   params->baitParams_.meanDepth          =   paramValues[i++];
   params->baitParams_.photo_alpha        =   paramValues[i++];
   
   params->baitParams_.sporeHalfLife      =   paramValues[i++];
   params->baitParams_.sporeProd          =   paramValues[i++];
   params->baitParams_.totalSpore         =   paramValues[i++];
   params->baitParams_.probGerminate      =   paramValues[i++];
   params->baitParams_.diffusionRate      =   paramValues[i++];
   params->baitParams_.diffusLoops        =   paramValues[i++];
   
   params->baitParams_.cos_units          =   paramValues[i++];
   params->baitParams_.tempAmp            =   paramValues[i++];
   params->baitParams_.temp_c             =   paramValues[i++];
   params->baitParams_.temp_d             =   paramValues[i++];
   params->baitParams_.solarAmp           =   paramValues[i++];
   params->baitParams_.solar_c            =   paramValues[i++];
   params->baitParams_.solar_d            =   paramValues[i++];
   params->baitParams_.DL_Amp             =   paramValues[i++];
   params->baitParams_.DL_c               =   paramValues[i++];
   params->baitParams_.DL_d               =   paramValues[i++];
   params->baitParams_.kdPAR              =   paramValues[i++];

   // Lattice Boltzmann parameters
   params->baitParams_.initRho            =   paramValues[i++];
   params->baitParams_.initUvelocity      =   paramValues[i++];
   params->baitParams_.initVvelocity      =   paramValues[i++];
   params->baitParams_.omega              =   paramValues[i++];
   params->baitParams_.uForce             =   paramValues[i++];
   params->baitParams_.obstXPos           =   static_cast<unsigned>(paramValues[i++]);
   params->baitParams_.obstYPos           =   static_cast<unsigned>(paramValues[i++]);
   params->baitParams_.obstLength         =   static_cast<unsigned>(paramValues[i++]);
   params->baitParams_.obstWidth          =   static_cast<unsigned>(paramValues[i++]);
   params->baitParams_.loopAddObst        =   static_cast<int>(paramValues[i++]);
   params->baitParams_.loopRemoveObst     =   static_cast<int>(paramValues[i++]);
   params->baitParams_.rankAddObst        =   static_cast<int>(paramValues[i++]);
   params->baitParams_.ratio_hH           =   paramValues[i++];
   
   params->baitParams_.displayOrganismsNumber_ = true ;
   params->baitParams_.initDistribute     =  static_cast<int>(paramValues[i++]);

   // Display parameters
   /// Probleme avec ca...
   params->xParams_.nbPatchColors_        =   16 ;
   params->xParams_.nbOrganismColors_     =   16 ;
   params->xParams_.scale_                =   1.5 ;
}


//
// JM 18Dec2006 - Function for creating output file name
//
static char* nameOutputFile(time_t start, int myRank, string inFile, string identifier)
{
   ostringstream outputString;

   string inFilename = inFile.erase(0, 8);
   inFilename = inFilename.erase(inFilename.length()-3, 3);
   inFilename = inFilename + identifier;
   outputString << "./results/" << inFilename << ".out";
   string outputFilename = outputString.str();

   if (myRank == 0)
      cout << outputFilename << endl;

   int stringLength = outputFilename.length();
   char *filenamePtr = new char[stringLength+1];
   outputFilename.copy(filenamePtr, stringLength, 0);
   filenamePtr[stringLength] = '\0';

   return filenamePtr;
}

//
// JM 26APR11 - Function for creating output file names for Lattice Boltzmann files
//
static char* nameOutputFileLB(time_t start, int myRank, string inFile, string identifier)
{
   ostringstream outputString;

   string inFilename = inFile.erase(0, 8);
   inFilename = inFilename.erase(inFilename.length()-3, 3);
   inFilename = inFilename + "_" + identifier;
   outputString << "./results/lbe/" << inFilename << "_" << myRank << ".lbe";
   string outputFilename = outputString.str();

   if (myRank == 0)
      cout << outputFilename << endl;

   int stringLength = outputFilename.length();
   char *filenamePtr = new char[stringLength+1];
   outputFilename.copy(filenamePtr, stringLength, 0);
   filenamePtr[stringLength] = '\0';

   return filenamePtr;
}



static void endTimekeeping(Bait *simulation, time_t start, int myRank)
{
   time_t stop;
   unsigned duration=0;//, refTime=770;
   
   // JM - timekeeping code
   if (myRank == 0)
   {
      stop = time(NULL);
      duration = (stop-start);

      simulation->printDuration(duration);      
      cout << "Total program duration (time_t) = " << (duration/60) << "m"
            << (duration%60) << "s " << "\n" 
            /*<< "Rel to P4 530 (plate400,nut400,mol20) = " 
            << ((double)duration/refTime)*/ << endl;
   }

}


//
// JM 04Nov2011 - Code to input pixel values from input files 2 and 3 to elevation array
//                and agent locations array 
//                (Converted raster image file to pixel values with ImageJ software)
//
static void inputMapData(int *xy_elevationsPtr, int *agentLocationsPtr, ifstream *inFileElevPtr, ifstream *inFileAgentsPtr, Parameters * paramsPtr)
{
   // Input x-y elevation co-ordinates from text file
   // (Convert Raster Image with ImageJ software:
   // Plugins -> Tools -> Display Pixel Values)
   int width = paramsPtr->engParams_.width_ ;
   int height = paramsPtr->engParams_.height_;
   string junk="";
   
   for (int i=-1; i<height; i++)
   { 
      if (i < 0)     // ignore first line
      {
         for (int j=-1; j<width; j++)
         {
            (*inFileElevPtr) >> junk;
            //(*inFileAgentsPtr) >> junk;
         }
      }
      else
      {
         for (int j=-2; j<width; j++)
         {
            if (j < 0)
            {
               // ignore first two columns
               (*inFileElevPtr) >> junk;
            }
            else
            {
               (*inFileElevPtr) >> xy_elevationsPtr[(i*width)+j];
            }
         }
       }
   }  // end for
}



//
// JM 04Nov2011 - Code to input pixel values from input files 2 and 3 to elevation array 
//                and agent locations array (parallel code) 
//                (Converted raster image file to pixel values with ImageJ software:
//                Plugins -> Tools -> Display Pixel Values)
//
static void inputMapDataMPI(int *xy_elevationsPtr, ifstream *inFileElevPtr, Parameters * paramsPtr, int mpiRank, int mpiSize)
{
   // Input x-y elevation co-ordinates from text file
   int mpiDim = sqrt(mpiSize);
   
   // Total width/height includes all plates, minus buffer zones
   int width = paramsPtr->engParams_.width_;
   int height = paramsPtr->engParams_.height_;
   int totWidth = ((width-2) * mpiDim)+2;
   int totHeight = ((height-2) * mpiDim)+2;
   string junk="";
   int buffer = 2;
   
   int *all_xy_elevations = new int[totWidth*totHeight];
    
   // Temporarily store all pixel values from raster images in arrays
   for (int i=-1; i<totHeight; i++)
   { 
      // ignore first line
      if (i < 0)
      {
         for (int j=-1; j<totWidth; j++)
         {
            (*inFileElevPtr) >> junk;
         }
      }
      else
      {
         for (int j=-2; j<totWidth; j++)
         {
            if (j < 0)
            {
               // ignore first two columns
               (*inFileElevPtr) >> junk;
            }
            else
            {
               // input actual data to arrays
               (*inFileElevPtr) >> all_xy_elevations[(i*totWidth)+j];
            }
         }
         if (mpiRank == 0 && i<2)
            cout << endl;
      }
   }  // end for 
    
   
   int plateRow = mpiRank/mpiDim;
   int plateRowSize = totWidth*height;
   int plateCol = mpiRank % mpiDim;
   int index = 0, all_xy_index = 0;
   
   // Divide up raster image values into different overlapping plates:
   // First row of plates
   // Cycle through each patch row in plate
   index=0;
   for (int i=plateRowSize*plateRow; i<plateRowSize*(plateRow + 1); i+=totWidth)
   { 
      // Cycle through each patch colomn in plate
      for (int j=width*plateCol; j<width*(plateCol + 1); j++)
      {
         all_xy_index = (i-plateRow*totWidth*buffer)+(j-plateCol*buffer);

         xy_elevationsPtr[index] = all_xy_elevations[all_xy_index];
         index++;
      }     
   }  // end for

   delete [] all_xy_elevations;
}


//
// 16Jan12 - Input map of agent (Spartina) locations
// (Converted raster image file to xy coords with ImageJ software:
//                Analyze -> Tools -> Save XY Coordinates)
//
static void inputAgentLocations(int *agentLocationsPtr, ifstream *inFileAgentsPtr, Parameters * paramsPtr)
{
   int nbCells = paramsPtr->engParams_.width_ * paramsPtr->engParams_.height_;
   int x, y, r, g, b;
   
   // initialise agentLocations[] with no agents in any cell
   for (int i=0; i<nbCells; i++)
      agentLocationsPtr[i] = 0;
      
   while ( inFileAgentsPtr->good() )
   {
      // input pixel values (x, y, red, green, blue) from file (ImageJ software)
      (*inFileAgentsPtr) >> x;
      (*inFileAgentsPtr) >> y;
      (*inFileAgentsPtr) >> r;
      (*inFileAgentsPtr) >> g;
      (*inFileAgentsPtr) >> b;
   
      // if pixel red in colour then agent present
      if ((double)r/b > 5.0 && r > 200) //(b < 100 && r > b*2 && r > 50)
         agentLocationsPtr[y*paramsPtr->engParams_.width_+x] = 1;
   }
   
}


//
// 16Jan12 - Input map of agent (Spartina) locations (Parallelised code)
// (Converted raster image file to xy coords with ImageJ software:
//                Analyze -> Tools -> Save XY Coordinates)
//
static void inputAgentLocationsMPI(int *agentLocationsPtr, ifstream *inFileAgentsPtr, Parameters * paramsPtr, int mpiRank, int mpiSize)
{
   int mpiDim = sqrt(mpiSize);
   
   // Total width/height includes all plates, minus buffer zones
   int width = paramsPtr->engParams_.width_;
   int height = paramsPtr->engParams_.height_;
   int totWidth = ((width-2) * mpiDim)+2;
   int totHeight = ((height-2) * mpiDim)+2;
   int buffer = 2;
   
   int *allAgentLocations = new int[totWidth*totHeight];

   int nbCells = totWidth*totHeight;
   int x, y, r, g, b;
   
   
   // 1. initialise allAgentLocations[] with no agents in any cell            //
   for (int i=0; i<nbCells; i++)
      allAgentLocations[i] = 0;

   // 2. Input from file to allAgentLocations[]                               //
   while ( inFileAgentsPtr->good() )
   {
      // input pixel values (x, y, red, green, blue) from file (ImageJ software)
      (*inFileAgentsPtr) >> x;
      (*inFileAgentsPtr) >> y;
      (*inFileAgentsPtr) >> r;
      (*inFileAgentsPtr) >> g;
      (*inFileAgentsPtr) >> b;
   
      // if pixel reddish in colour then agent present
      if ((double)r/b > 1.5 && r > 50) //(b < 100 && r > b*2 && r > 50)
         allAgentLocations[y*totWidth+x] = 1;
   }
   
   int plateRow = mpiRank/mpiDim;
   int plateRowSize = totWidth*height;
   int plateCol = mpiRank % mpiDim;
   int index = 0, all_xy_index = 0;
             
   // 3. Divide up raster image values into different overlapping plates:     //
   // First row of plates
   // Cycle through each patch row in plate
   for (int i=plateRowSize*plateRow; i<plateRowSize*(plateRow + 1); i+=totWidth)
   { 
      // Cycle through each patch colomn in plate
      for (int j=width*plateCol; j<width*(plateCol + 1); j++)
      {
         all_xy_index = (i-plateRow*totWidth*buffer)+(j-plateCol*buffer);
         agentLocationsPtr[index] = allAgentLocations[all_xy_index];

         index++;
      }    
   }  // end for
   delete [] allAgentLocations;
}


//
// 30Mar12 - Input map of niche(Spartina) locations
// (Converted raster image file to xy coords with ImageJ software:
//                Analyze -> Tools -> Save XY Coordinates)
//
static void inputSubstrateMap(int *substrateLocationsPtr, ifstream *inFileNichePtr, Parameters * paramsPtr)
{
   int nbCells = paramsPtr->engParams_.width_ * paramsPtr->engParams_.height_;
   int x, y, r, g, b;
   
   // initialise substrateLocations[] with all cells negative
   for (int i=0; i<nbCells; i++)
      substrateLocationsPtr[i] = 0;
      
   while ( inFileNichePtr->good() )
   {
      // input pixel values (x, y, red, green, blue) from file (ImageJ software)
      (*inFileNichePtr) >> x;
      (*inFileNichePtr) >> y;
      (*inFileNichePtr) >> r;
      (*inFileNichePtr) >> g;
      (*inFileNichePtr) >> b;
   
      // if red pixel -> concrete (high affinity substrate)
      if ((float )r/b > 4.0 && r > 100 && g < (float)r/2) //(b < 100 && r > b*2 && r > 50)
         substrateLocationsPtr[y*paramsPtr->engParams_.width_+x] = 1;
      // else if blue pixel -> plastic (low affinity substrate)
      else if ((float)b/r > 4.0 && b > 100 && g < (float)b*3/4 )
         substrateLocationsPtr[y*paramsPtr->engParams_.width_+x] = 2;
      // else if black pixel -> obstacle
      else if (r < 50 && g < 50 && b < 50)
         substrateLocationsPtr[y*paramsPtr->engParams_.width_+x] = -1;
   }
   
}



//
// 30Mar12 - Input map of niche (Spartina) locations (Parallelised code)
// (Converted raster image file to xy coords with ImageJ software:
//                Analyze -> Tools -> Save XY Coordinates)
//
static void inputSubstrateMapMPI(int *substrateLocationsPtr, ifstream *inFileNichePtr, Parameters *paramsPtr, int mpiRank, int mpiSize)
{
   int mpiDim = sqrt(mpiSize);
   
   // Total width/height includes all plates, minus buffer zones
   int width = paramsPtr->engParams_.width_;
   int height = paramsPtr->engParams_.height_;
   int totWidth = ((width-2) * mpiDim)+2;
   int totHeight = ((height-2) * mpiDim)+2;
   int buffer = 2;
   
   int *allSubstrateLocations = new int[totWidth*totHeight];

   int nbCells = totWidth*totHeight;
   int x, y, r, g, b;
   
   
   // 1. initialise allAgentLocations[] with no agents in any cell            //
   for (int i=0; i<nbCells; i++)
      allSubstrateLocations[i] = 0;

   // 2. Input from file to allAgentLocations[]                               //
   while ( inFileNichePtr->good() )
   {
      // input pixel values (x, y, red, green, blue) from file (ImageJ software)
      (*inFileNichePtr) >> x;
      (*inFileNichePtr) >> y;
      (*inFileNichePtr) >> r;
      (*inFileNichePtr) >> g;
      (*inFileNichePtr) >> b;
      
      // if red pixel -> concrete (high affinity substrate)
      if ((float )r/b > 10.0 && r > 100 && g < (float)r/2) //(b < 100 && r > b*2 && r > 50)
         allSubstrateLocations[y*totWidth+x] = 1;
      // else if blue pixel -> plastic (low affinity substrate)
      else if ((float)b/r > 10.0 && b > 100 && g < (float)b*3/4 )
         allSubstrateLocations[y*totWidth+x] = 2;
      // else if black pixel -> obstacle
      else if (r < 50 && g < 50 && b < 50)
         allSubstrateLocations[y*totWidth+x] = -1;
   }
   
   int plateRow = mpiRank/mpiDim;
   int plateRowSize = totWidth*height;
   int plateCol = mpiRank % mpiDim;
   int index = 0, all_xy_index = 0;
             
   // 3. Divide up raster image values into different overlapping plates:     //
   // First row of plates
   // Cycle through each patch row in plate
   for (int i=plateRowSize*plateRow; i<plateRowSize*(plateRow + 1); i+=totWidth)
   { 
      // Cycle through each patch colomn in plate
      for (int j=width*plateCol; j<width*(plateCol + 1); j++)
      {
         all_xy_index = (i-plateRow*totWidth*buffer)+(j-plateCol*buffer);
         substrateLocationsPtr[index] = allSubstrateLocations[all_xy_index];

         index++;
      }    
   }  // end for 
   
   delete [] allSubstrateLocations;
}




// ------------------------------ End Of File ------------------------------- //
