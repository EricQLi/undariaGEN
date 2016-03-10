// -------------------------------------------------------------------------- //
// plate.cpp                                                                  //
// -------------------------------------------------------------------------- //

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

#include "plate.hpp"


// -------------------------------------------------------------------------- //
// Libraries                                                                  //
// -------------------------------------------------------------------------- //

#include "organism.hpp"
//#include "antibio.hpp"
#include "baitParameters.hpp"
#include <cmath>
using std::sqrt;

#include <mpi.h>

#include <iomanip>
using std::setprecision;

#include "bait.hpp"
#include <gsl/gsl_randist.h>

//#include <cuda.h>
//#include <cuda_runtime.h>
#include "math.h"


/*
// -------------------------------------------------------------------------- //
// CUDA Global Kernel Functions & Variables                                   //
// -------------------------------------------------------------------------- //
// CUDA constant variables declaration
__constant__ int widthDev;
__constant__ float halfLifeDev;
__constant__ float diffusionDev;

//
// CUDA Kernel for Diffusion algorithm to be launched on GPU device
//
__global__ void diffusion_kernel(float *sporeDev, float *newSporeDev)
{
   // Only process non-border patches: transform CUDA indices (i.e. all non-border cells: width-2, height-2)
   // to plate co-ordinates
   int idx = (widthDev+1) + (2*blockIdx.x) + (threadIdx.x + blockIdx.x * blockDim.x);
   
   newSporeDev[idx] = sporeDev[idx] 
               - diffusionDev * (sporeDev[idx]-sporeDev[idx-1])
               - diffusionDev * (sporeDev[idx]-sporeDev[idx+1])
               - diffusionDev * (sporeDev[idx]-sporeDev[idx-widthDev])
               - diffusionDev * (sporeDev[idx]-sporeDev[idx+widthDev])
               - (diffusionDev/1.414213562) * (sporeDev[idx]-sporeDev[idx-widthDev-1])
               - (diffusionDev/1.414213562) * (sporeDev[idx]-sporeDev[idx-widthDev+1])
               - (diffusionDev/1.414213562) * (sporeDev[idx]-sporeDev[idx+widthDev-1])
               - (diffusionDev/1.414213562) * (sporeDev[idx]-sporeDev[idx+widthDev+1]);

      
   // Apply half-life algorithm:
   newSporeDev[idx] = newSporeDev[idx] * pow((float)0.5, (float)halfLifeDev);
}

__global__ void update_kernel(float *sporeDev, float *newSporeDev)
{
   int idx = (widthDev+1) + (2*blockIdx.x) + (threadIdx.x + blockIdx.x * blockDim.x);
   
   sporeDev[idx] = newSporeDev[idx];
}
*/

// -------------------------------------------------------------------------- //
// Static member creation                                                     //
// -------------------------------------------------------------------------- //

BaitParameters* Plate::baitParamsPtr_ = 0 ;
const unsigned Plate::NUM_ADJ_PLATES=9;        // no. adjacent plates: 9 if 2D, 27 if 3D
const double Plate::CS2 = 1.0/3.0;
// Hour at which each calendar month starts (yr = 8766 h/365.25 days)
const unsigned Plate::stMonthsHr[] = {0,744,1422,2166,2886,3630,4350,5094,5838,6558,7302,8022};
// length of each calendar month in hours (feb = +6 hours to account for leap years)
const unsigned Plate::monthLength[] = {744,678,744,720,744,720,744,744,720,744,720,744};
// days per month
const unsigned Plate::daysPerMonth[12] = {31,28,31,30,31,30,31,31,30,31,30,31};


// -------------------------------------------------------------------------- //
// Constructor                                                                //
// -------------------------------------------------------------------------- //

Plate::Plate()
:World<Patch>()
{
   iStart = 0;
   iStop = 0;
   jStart = 0;
   jStop = 0;
   for (unsigned i=0; i<NUM_ADJ_PLATES;i++)
   {
      activeEdgeLength[i] = 0;
      adjacentPlate[i] = false;
      sndBufSize[i] = 0;
      rcvBufSize[i] = 0;
   }
}


// -------------------------------------------------------------------------- //
// Destructor                                                                 //
// -------------------------------------------------------------------------- //

Plate::~Plate()
{
   delete [] patchValues; delete [] biomassValues;
   delete [] elevationValues; delete [] sporeValues;
   delete [] substrateValues;
   delete [] allPatchValues; delete [] allBiomassValues;
   delete [] allElevationValues; delete [] allSporeValues;
   delete [] allSubstrateValues;
   delete [] dailyWaterTemp;
   delete [] dailySolarRad;
}


// -------------------------------------------------------------------------- //
// Initialiaze the static parameters pointer                                  //
// -------------------------------------------------------------------------- //

void Plate::initBaitParamsPtr( BaitParameters *inPtr)
{
   baitParamsPtr_ = inPtr ;
}



// -------------------------------------------------------------------------- //
// Specific initialisation                                                    //
// -------------------------------------------------------------------------- //

void Plate::init(int rank, unsigned mpiSize, int *xy_elevs, int *inAgentLocations, int *inSubstrateTypes, float *waterTempIn, float *solarRadIn, float *dayLengthIn)
{
   // Patch static pointer intialisation
   Patch::initBaitParamsPtr( baitParamsPtr_, engParamsPtr_ ) ;
   Patch::initWorldPtr( (World<Patch> *)this ) ;

   // World initialisation
   World<Patch>::init(mpiSize) ;
   unsigned plateSize = engParamsPtr_->width_ * engParamsPtr_->height_;
   patchValues = new float[plateSize];
   allPatchValues = new float[plateSize * mpiSize];
   biomassValues = new float[plateSize];
   allBiomassValues = new float[plateSize * mpiSize];
   elevationValues = new float[plateSize];
   allElevationValues = new float[plateSize * mpiSize];
   substrateValues = new int[plateSize];
   allSubstrateValues = new int[plateSize * mpiSize];
   sporeValues = new float[plateSize];
   allSporeValues = new float[plateSize * mpiSize];
   dailyWaterTemp = new float[365];
   dailySolarRad = new float[365];
   dayLengths = new float[365];
   waterLevel = 0.0;
   currWaterTemp = 0.0; 
   currSolarRad = 0.0;
   currDayLength = 0.0; 

   // initialise monthly sea water temp from inputted averages     
   for (unsigned i=0; i<365; i++)
   {
      dailyWaterTemp[i] = waterTempIn[i]; 
      dailySolarRad[i] = solarRadIn[i];
      dayLengths[i] = dayLengthIn[i];
   }
   //for (unsigned i=0; i<365*4; i++)
   //  dailyWaterTemp[i] = waterTempIn[i];

   myRank = rank;
   xy_elevations = xy_elevs;
   agentLocations = inAgentLocations;
   substrateTypes = inSubstrateTypes;
   numPlatesMPI = mpiSize;

   // Determine active patches (not in MPI buffer zones)
   setActiveRegion(mpiSize);  
   setActivePatches();
   
   // Determine upper/lower limits of S anglica growth (Gray's Model, 1989, 1995)
   //setNicheLimits();
   initElevNicheSubstrate(); 

   initUndariaNiche();

   int mpiLength = (int)sqrt((float)mpiSize);
   // identify adjacent plates, 1=E,2=N,3=W,4=S,5=NE,6=NW,7=SW,8=SE
   if ( (myRank+1)%mpiLength == 0 )         // if at right edge then move to left edge
   {                                      // periodic boundary condition
      dest[1] = myRank - (mpiLength-1);
      dest[5] = myRank - (mpiLength*2-1);
      dest[8] = myRank + 1;
   }       
   else
   {
      dest[1] = myRank+1;
      dest[5] = myRank - (mpiLength-1);
      dest[8] = myRank + (mpiLength+1);
   }
   if ( myRank%mpiLength == 0 )                 // if at left edge then move to right edge
   {                                      // periodic boundary condition
      dest[3] = myRank + (mpiLength-1);
      dest[6] = myRank - 1;
      dest[7] = myRank + (mpiLength*2-1);        
   }
   else
   {
      dest[3] = myRank - 1;
      dest[6] = myRank - (mpiLength+1);
      dest[7] = myRank + (mpiLength-1);
   }
   // Closed boundary at North/South edges (not periodic, yet)
   dest[2] = myRank-mpiLength;
   dest[4] = myRank+mpiLength;
   dest[0] = myRank;
   
   // Record total size of environment (if in parallel, all plates)
   if (mpiSize >= 4)
   {
      totEnvWidth = (engParamsPtr_->width_-2) * mpiLength + 2;
      totEnvHeight = (engParamsPtr_->height_-2) * mpiLength + 2;
   }
   else
   {
      totEnvWidth = engParamsPtr_->width_;
      totEnvHeight = engParamsPtr_->height_;
   }
   if (myRank == 0)
      cout << "Environment Size (x*y): " << totEnvWidth << " * " << totEnvHeight << endl; 
   
   // Meaning of x values: 0=NONE,1=E,2=N,3=W,4=S,5=NE,6=NW,7=SW,8=SE
   // Manually specify co-ords at which to start filling send/receive buffers 
   // for sendAndReceive() and parallelMovement()
   // For diffusion buffers (parDiff) filled from second-innermost row/col
   // For parallel movement of bac, filled from outermost row/col
   for (unsigned plate_dir=1; plate_dir<NUM_ADJ_PLATES; plate_dir++)
   {
      switch(plate_dir)
      {
         case 1:
            jPos_parDiff[plate_dir]=engParamsPtr_->width_-2;
            iPos_parDiff[plate_dir]=iStart;
            jPos_parMove[plate_dir]=engParamsPtr_->width_-1;
            iPos_parMove[plate_dir]=0;
            break;
         case 2:
            jPos_parDiff[plate_dir]=jStart;
            iPos_parDiff[plate_dir]=1;          
            jPos_parMove[plate_dir]=0;
            iPos_parMove[plate_dir]=0;
            break;
         case 3:
            jPos_parDiff[plate_dir]=1;
            iPos_parDiff[plate_dir]=iStart;
            jPos_parMove[plate_dir]=0;
            iPos_parMove[plate_dir]=0;
            break;
         case 4:
            jPos_parDiff[plate_dir]=jStart;
            iPos_parDiff[plate_dir]=engParamsPtr_->height_-2;
            jPos_parMove[plate_dir]=0;
            iPos_parMove[plate_dir]=engParamsPtr_->height_-1;
            break;
         case 5:
            jPos_parDiff[plate_dir]=engParamsPtr_->width_-2;
            iPos_parDiff[plate_dir]=1;
            jPos_parMove[plate_dir]=engParamsPtr_->width_-1;
            iPos_parMove[plate_dir]=0;
            break;
         case 6:
            jPos_parDiff[plate_dir]=1;
            iPos_parDiff[plate_dir]=1;
            jPos_parMove[plate_dir]=0;
            iPos_parMove[plate_dir]=0;
            break;
         case 7:
            jPos_parDiff[plate_dir]=1;
            iPos_parDiff[plate_dir]=engParamsPtr_->height_-2;
            jPos_parMove[plate_dir]=0;
            iPos_parMove[plate_dir]=engParamsPtr_->height_-1;
            break;
         case 8:
            jPos_parDiff[plate_dir]=engParamsPtr_->width_-2;
            iPos_parDiff[plate_dir]=engParamsPtr_->height_-2;
            jPos_parMove[plate_dir]=engParamsPtr_->width_-1;
            iPos_parMove[plate_dir]=engParamsPtr_->height_-1;
            break;
         default:
            jPos_parDiff[plate_dir]=0;
            iPos_parDiff[plate_dir]=0;
            jPos_parMove[plate_dir]=0;
            iPos_parMove[plate_dir]=0;
            break;
      }
   } // end for loop
   
   //Animp moive file initialization
   //char *file = (char *)"./results/movie_bac.an";
   //m = new movie (file);
   
   
   // Initialise Lattice Bolztmann (LBE) variables
   omega = baitParamsPtr_->omega;
   umoy = 0.0;
   vmoy = 0.0;
   avgDensity = 0.0;
   
   // Calculation of the viscosity (LBE)
   visc = (1.0/omega - 0.5) * CS2;
   if (myRank == 0)
      cout << "Viscosity (myRank 0): " << visc << endl;
   if (visc < 0.0)
      cout << "STOP! Omega out of (0,2) interval." << endl;
   
   // Calculation of constant applied force (Stokes problem) for Poisseulle flow (LBE)
   force = 8.0 * visc * baitParamsPtr_->uForce / (double)(totEnvHeight) / (double)(totEnvHeight);
   force = baitParamsPtr_->initRho * force/6; 
   if (myRank == 0)
      cout << "Intensity of the applied force (myRank 0): " << force << endl; 
}


// -------------------------------------------------------------------------- //
// Solve conflicts in each cell of the world                                  //
// Not needed - used to be used in Micro-Gen for Bac-Ab kinetic interactions  //
//                                                                            //
// -------------------------------------------------------------------------- //

void Plate::solveConflicts(unsigned nbLoops)
{
   // Solve all conflicts locally
   /*unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;
   for( unsigned i=0; i<nbCells; i++ )
   {
      if (this->cells_[i].getActiveRegion() == true)
         this->cells_[i].solveConflicts(nbLoops) ;
   }*/
}


//
// 31Jan14 - update patches (based on World::allCellsDevelop())
//
void Plate::allPatchesDevelop(unsigned nbLoops)
{
   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;
   
   for( unsigned i=0; i<nbCells; i++ )
   {
      if (this->cells_[i].getActiveRegion() && this->cells_[i].getSubstrateType() > -1)
         this->cells_[i].development(nbLoops);
   }
}


//
// 23Jan2014 - Cycle through patches - if substrateType > 0
//             probability of spore germination proportional to num spores.
//
void Plate::germinateSpores(unsigned nbLoops)
{
   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;
   
   // Spores only germinate under appropriate temperature conditions
   // see (D. Grulois thesis, 2010) and (Morita et al. 2003)
   for( unsigned i=0; i<nbCells; i++ )
      this->cells_[i].germinateSpores(i, currWaterTemp, myRank);
}


// -------------------------------------------------------------------------- //
// Update Animp 3D visualization display                                      //
// JM 12May2008                                                               //
//                                                                            //
// -------------------------------------------------------------------------- //

/*void Plate::animpDisplay( unsigned inLoopNumber, Fabric<Organism>& inOrganismsFabric, unsigned totalNumBac)
{	
   int size =  1; 
		
	// Adding Colors to the Movie
   m->add ( ap_color ( 0.9, 0.7, 0.3, 0.5 ) ); // Yellow

   // Adding frames and spheres to 3D visualisation movie
   ap_frame *f = new ap_frame(totalNumBac);
      
   for( int i=0; i<engParamsPtr_->height_; i++ )
   {
      for( int j=0; j<engParamsPtr_->width_; ++j )
      {
         unsigned numBac = (*this)(i,j).getNumOrganism() ;        
         if (numBac >= 1 )
            for (unsigned z=0; z<numBac; z+=4)
               f->add ( sphere_s(i,j,((unsigned)(z/4)),size,0) );
      }
   }	
   m->add ( f );
   delete f;
}

//
// JM 12May2008 - Write to Animp 3d visualization file.
//
void Plate::animpWrite()
{
   m->write();
}*/

//
// JM 26Nov2007 - Diffusion algorithm for serial program execution.
//
void Plate::diffusion(unsigned nbLoops)
{  
   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ; 

   // 1. Apply diffusion algorithm to each patch
   for( unsigned i=0; i<nbCells; i++ )
      if (this->cells_[i].getActiveRegion() == true) {
         this->cells_[i].sporeDiffusion(nbLoops);
      }
   // 2. Update spore level in each patch.
   for( unsigned i=0; i<nbCells; i++ )
   {
      if (this->cells_[i].getActiveRegion() == true) {
         this->cells_[i].updateSporeLevel(); 
      }
   }
   
}

//
// CUDA Diffusion algorithm for nvidia GPU (Feb2014)
// *Only works on non-border patches*
//
void Plate::diffusionCUDA(unsigned nbLoops)
{/*
   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;
   const int width = engParamsPtr_->width_;
   const int height = engParamsPtr_->height_;
   const float halfLife = 1.0/(float)baitParamsPtr_->sporeHalfLife;
   const float diffusionRate = (float)baitParamsPtr_->diffusionRate;
   
   // Don't include border patches
   unsigned nbCUDACells = nbCells - (width*2 + (height*2-4));   
   unsigned nb = nbCells * sizeof(float);
   
   float *sporeHost = new float[nbCells];
   float *newSporeHost = new float[nbCells];
   float *newSporeDev = NULL, 
         *sporeDev = NULL;
             

   // 1. Copy sporeLevel values for all patches to array for passing to CUDA GPU
   for( unsigned i=0; i<nbCells; i++ )
   {
      sporeHost[i] = (float)this->cells_[i].getSporeLevel();
      newSporeHost[i] = 0.0;
   }
   
   // 2. Allocates memory on GPU
   cudaError_t cuerr = cudaMalloc ( (void**)&newSporeDev, nb );
   cuerr = cudaMalloc ( (void**)&sporeDev, nb );
   cuerr = cudaMalloc( (void**)&widthDev, sizeof(int) );
   cuerr = cudaMalloc( (void**)&halfLifeDev, sizeof(float) );
   cuerr = cudaMalloc( (void**)&diffusionDev, sizeof(float) );
   if (cuerr != cudaSuccess)
   {
      cout << "Plate: Cannot allocate GPU memory: " << stderr 
         << cudaGetErrorString(cuerr) << endl;
      //return 1;
   }
   
   // copy integer value for plate width to GPU
   cuerr = cudaMemcpyToSymbol(widthDev, (const void *)&width, sizeof(int));//, 0, cudaMemcpyHostToDevice);
   cuerr = cudaMemcpyToSymbol(halfLifeDev, (const void *)&halfLife, sizeof(float));
   cuerr = cudaMemcpyToSymbol(diffusionDev, (const void *)&diffusionRate, sizeof(float));
   if (cuerr != cudaSuccess)
   {
      cout << "Plate: cudaMemcpyToSymbol failure: " << stderr << " " << cudaGetErrorString(cuerr) << endl;
      // return 1;
   }
   
   // 3. Set up the kernel launch configuration for n threads
   // CUDA is applied to all non-border patches, i.e. dimensions = (width-2) x (height-2) 
   const int BLOCK_SIZE = engParamsPtr_->width_-2;
   dim3 threads = dim3(BLOCK_SIZE, 1);
   unsigned n = nbCUDACells;
   dim3 blocks  = dim3(n / BLOCK_SIZE, 1);
   
   // 4. Copy sporeLevel data from host to GPU global memory
   cuerr = cudaMemcpy(sporeDev, sporeHost, nb, cudaMemcpyHostToDevice);
   cuerr = cudaMemcpy(newSporeDev, newSporeHost, nb, cudaMemcpyHostToDevice);
   cuerr = cudaGetLastError();
   if (cuerr != cudaSuccess)
   {
      cout << "Plate: cudaMemcpy failure: " << stderr << " " << cudaGetErrorString(cuerr) << endl;
      // return 1;
   }
   
   // 5. Execute diffusion kernel with the specified config and args
   for (unsigned i=0; i<baitParamsPtr_->CUDALoops;i++)
   {
      // a. Apply Fick's first law of diffusion
      diffusion_kernel<<<blocks, threads>>>(sporeDev, newSporeDev);
      cuerr = cudaGetLastError();
      if (cuerr != cudaSuccess)
      {
         cout << "Plate: Cannot launch CUDA kernel: " << stderr << " " << cudaGetErrorString(cuerr) << endl;
         // return 1;
      }
      cuerr = cudaDeviceSynchronize();
      
      // b. update spore levels
      update_kernel<<<blocks, threads>>>(sporeDev, newSporeDev);
      cuerr = cudaDeviceSynchronize();
   }
   
   // 6. Wait for the kernel to finish
   cuerr = cudaDeviceSynchronize();
   if (cuerr != cudaSuccess)
   {
       cout << "Plate: Cannot synchronize CUDA kernel: " << stderr << " " << cudaGetErrorString(cuerr) << endl;
      // return 1;
   }
   
   // 7. Copy the results back to CPU memory
   cuerr = cudaMemcpy (newSporeHost, newSporeDev, nb, cudaMemcpyDeviceToHost);
   if (cuerr != cudaSuccess)
   {
      cout << "Plate: Cannot copy data from dev to host: " << stderr << " " << cudaGetErrorString(cuerr) << endl;
      // return 1;
   }
   
   // Free GPU memory
   cudaFree ( sporeDev );   
   cudaFree ( newSporeDev );
   //cudaFree ( widthDev );
   
   // 8. Copy updated spore levels to patches
   for( unsigned i=0; i<nbCells; i++ )
      if (newSporeHost[i] != sporeHost[i])
      {
         if (newSporeHost[i] >= 0.001)
            this->cells_[i].setSporeLevel(newSporeHost[i]);
         else  // if sporeLevel<0.5 then change to 0.0
            this->cells_[i].setSporeLevel(0.0);
      }
         
   
   delete [] sporeHost;
   delete [] newSporeHost;*/
}






//
// JM 05Dec2006 - Parallelisation
//                - Adjacent plates overlap by two rows/columns
//                  - Diffusion not calculated for adjoining row/column.
//
void Plate::parallelDiffusion(int size, unsigned nbLoops)
{
   int width = (int)sqrt((float)size);

   if (nbLoops % baitParamsPtr_->diffusLoops == 0)
   {
      plateDiffusionNew(nbLoops);
      //updateLevels(nbLoops);

      // Communicate with adjacent plates.
      sendAndReceive(width, nbLoops);
   }
}

void Plate::parallelDiffusion_CUDA(int size, unsigned nbLoops)
{
   /*int width = (int)sqrt((float)size);

   diffusionCUDA(nbLoops);

   // Communicate with adjacent plates.
   sendAndReceive(width, nbLoops);*/
}

//
// JM 12MAY11 - No diffusion, just call MPI routine
//
void Plate::noDiffusionMPI(int size, unsigned nbLoops)
{
   int width = (int)sqrt((float)size);
   
   // Communicate with adjacent plates.
   sendAndReceive(width, nbLoops);
}


//
// JM 05Mar2007 - Determine the patches in which diffusion/movement calculations
//                  are carried out (overlapping borders between adjacent plates).
//
void Plate::setActiveRegion(unsigned size)
{
   unsigned width = (unsigned)sqrt((float)size);

   if (size >= 4)
   {  
         // Top left corner
         if ( myRank == 0 )
         {
            iStart = 0;
            iStop = engParamsPtr_->height_-1;
            jStart = 0;
            jStop = engParamsPtr_->width_-1;

            activeEdgeLength[1] = iStop - iStart;
            activeEdgeLength[4] = jStop - jStart;
            activeEdgeLength[8] = 1;
         }
         // Top right corner
         else if (myRank == width-1)
         {
            iStart = 0;
            iStop = engParamsPtr_->height_-1;
            jStart = 1;
            jStop = engParamsPtr_->width_;

            activeEdgeLength[3] = iStop - iStart;
            activeEdgeLength[4] = jStop - jStart;
            activeEdgeLength[7] = 1;
         }
         // Bottom left corner
         else if (myRank == size-width)
         {
            iStart = 1;
            iStop = engParamsPtr_->height_;
            jStart = 0;
            jStop = engParamsPtr_->width_-1;

            activeEdgeLength[1] = iStop - iStart;
            activeEdgeLength[2] = jStop - jStart;
            activeEdgeLength[5] = 1;
         }
         // Bottom right corner
         else if ( myRank == (size-1) )
         {
            iStart = 1;
            iStop = engParamsPtr_->height_;
            jStart = 1;
            jStop = engParamsPtr_->width_;

            activeEdgeLength[3] = iStop - iStart;
            activeEdgeLength[2] = jStop - jStart;  
            activeEdgeLength[6] = 1;
         }
         // Right edge
         else if ( (myRank+1)%width == 0)
         {
            iStart = 1;
            iStop = engParamsPtr_->height_-1;
            jStart = 1;
            jStop = engParamsPtr_->width_;

            activeEdgeLength[3] = iStop - iStart;
            activeEdgeLength[2] = jStop - jStart;
            activeEdgeLength[4] = jStop - jStart;
            activeEdgeLength[7] = 1;
            activeEdgeLength[6] = 1;
         }
         // Left edge
         else if ( myRank%width == 0)
         {
            iStart = 1;
            iStop = engParamsPtr_->height_-1;
            jStart = 0;
            jStop = engParamsPtr_->width_-1;

            activeEdgeLength[1] = iStop - iStart;
            activeEdgeLength[2] = jStop - jStart;
            activeEdgeLength[4] = jStop - jStart;
            activeEdgeLength[5] = 1;
            activeEdgeLength[8] = 1;
         }
         // Top edge
         else if ( myRank < width-1 )
         {
            iStart = 0;
            iStop = engParamsPtr_->height_-1;
            jStart = 1;
            jStop = engParamsPtr_->width_-1;

            activeEdgeLength[3] = iStop - iStart;
            activeEdgeLength[1] = iStop - iStart; 
            activeEdgeLength[4] = jStop - jStart;
            activeEdgeLength[8] = 1;
            activeEdgeLength[7] = 1;
         }
         // Bottom edge
         else if (size-myRank < width)
         {
            iStart = 1;
            iStop = engParamsPtr_->height_;
            jStart = 1;
            jStop = engParamsPtr_->width_-1;

            activeEdgeLength[3] = iStop - iStart;
            activeEdgeLength[1] = iStop - iStart;
            activeEdgeLength[2] = jStop - jStart;
            activeEdgeLength[5] = 1;
            activeEdgeLength[6] = 1;
         }
         // Interior plates
         else
         {
            iStart = 1;
            iStop = engParamsPtr_->height_-1;
            jStart = 1;
            jStop = engParamsPtr_->width_-1;

            activeEdgeLength[3] = iStop - iStart;
            activeEdgeLength[1] = iStop - iStart;
            activeEdgeLength[2] = jStop - jStart;
            activeEdgeLength[4] = jStop - jStart;
            for (unsigned i=5; i<NUM_ADJ_PLATES; i++)
               activeEdgeLength[i] = 1;
         }
   }
   // Some special (non-square) cases
   else if (size == 3)
   {
      // Left plate
      if ( myRank == 0 )
      {
         iStart = 0;
         iStop = engParamsPtr_->height_;
         jStart = 0;
         jStop = engParamsPtr_->width_-1;

         activeEdgeLength[1] = iStop - iStart;
      }
      // middle plate
      else if (myRank == 1)
      {
         iStart = 0;
         iStop = engParamsPtr_->height_;
         jStart = 1;
         jStop = engParamsPtr_->width_-1;

         activeEdgeLength[3] = iStop - iStart;
         activeEdgeLength[1] = iStop - iStart;
      }
      // Right plate
      else if (myRank == 2)
      {
         iStart = 0;
         iStop = engParamsPtr_->height_;
         jStart = 1;
         jStop = engParamsPtr_->width_;

         activeEdgeLength[3] = iStop - iStart;
      }
   }
   else if (size == 2)
   {
      // Left plate
      if ( myRank == 0 )
      {
         iStart = 0;
         iStop = engParamsPtr_->height_;
         jStart = 0;
         jStop = engParamsPtr_->width_-1;

         activeEdgeLength[1] = iStop - iStart;
      }
      // Right plate
      else if (myRank == 1)
      {
         iStart = 0;
         iStop = engParamsPtr_->height_;
         jStart = 1;
         jStop = engParamsPtr_->width_;

         activeEdgeLength[3] = iStop - iStart;
      }
   }
   // If not parallelised - all patches included in active region.
   else
   {
      iStart = 0;
      iStop = engParamsPtr_->height_;
      jStart = 0;
      jStop = engParamsPtr_->width_;
   }
   
   // LBE - allow East/West flow across plates for MPI
   if (size > 1 && engParamsPtr_->worldType_ == OPEN_EW)
   {
      jStart = 1;
      jStop = engParamsPtr_->width_-1;
      
      if (myRank >= width)                   // if NOT top row of plates, can send up
      {
         activeEdgeLength[2] = jStop - jStart;
         activeEdgeLength[5] = 1;
         activeEdgeLength[6] = 1;
      }
      if (size-myRank > width)              // if NOT bottom row of plates, can send down
      {
         activeEdgeLength[4] = jStop - jStart;
         activeEdgeLength[7] = 1;
         activeEdgeLength[8] = 1;
      }
      
      activeEdgeLength[1] = iStop - iStart;  // always send East/West
      activeEdgeLength[3] = iStop - iStart;
   }
   
   for (unsigned plate_dir=1; plate_dir<NUM_ADJ_PLATES; plate_dir++)
      if (activeEdgeLength[plate_dir] > 0)
         adjacentPlate[plate_dir] = true;
}

//
// JM 05Dec2006 - Carry out diffusion in patches (called by parallelDiffusion()).
//
void Plate::plateDiffusion(unsigned nbLoops)
{   
   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_;
     
   for( unsigned i=0; i<nbCells; i++ )
   {
      if (this->cells_[i].getActiveRegion() == true)
         this->cells_[i].sporeDiffusion(nbLoops);
   }
}

//
// JM 07Mar14 - Faster implementation of Fick's First Law Diffusion algorithm
//              first copy all spore values to local array for processing
//              can also run for a number of loops at a time (diffusLoops) like CUDA
//
void Plate::plateDiffusionNew(unsigned nbLoops)
{
   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;
   const unsigned width = engParamsPtr_->width_;
   const unsigned height = engParamsPtr_->height_;
   const float halfLife = 1/baitParamsPtr_->sporeHalfLife;
   const float diffusionRate = baitParamsPtr_->diffusionRate;
   
   float *spore = new float[nbCells];
   float *newSpore = new float[nbCells];
             
   // 1. Copy sporeLevel values for all patches to local array
   for( unsigned i=0; i<nbCells; i++ )
   {
      spore[i] = this->cells_[i].getSporeLevel();
      newSpore[i] = spore[i];
   }
   
   // 2. Run diffusion algorithm on local array and apply half-life
   for (unsigned   n=0; n<baitParamsPtr_->diffusLoops; n++)
   {
      for( unsigned i=0; i<nbCells; i++ )
      {   
         // Apply half-life algorithm:
         if (newSpore[i] > 1.0E-06)
            newSpore[i] = newSpore[i] * pow((float)0.5, (float)halfLife);
            
         if (newSpore[i] != spore[i])
            spore[i] = newSpore[i];
      }
      
      // Apply diffusion algorithm to all non-border patches
      unsigned index=0;
      for( unsigned i=1; i<height-1; i++)
      {
         for( unsigned j=1; j<width-1; j++)
         {
            index = (i*width) + j;
            
            // Fick's First Law of Diffusion:
            newSpore[index] = spore[index]
                  - diffusionRate * (spore[index]-spore[index-1])
                  - diffusionRate * (spore[index]-spore[index+1])
                  - diffusionRate * (spore[index]-spore[index-width])
                  - diffusionRate * (spore[index]-spore[index+width])
                  - (diffusionRate/1.414213562) * (spore[index]-spore[index-width-1])
                  - (diffusionRate/1.414213562) * (spore[index]-spore[index-width+1])
                  - (diffusionRate/1.414213562) * (spore[index]-spore[index+width-1])
                  - (diffusionRate/1.414213562) * (spore[index]-spore[index+width+1]);      
         }
      }
      
   }
   
   // 3. Copy updated spore levels back to patches
   for( unsigned i=width+1; i<nbCells-(width+1); i++ )
   {
      if (newSpore[i] != spore[i] && this->cells_[i].getSubstrateType() > -1 && this->cells_[i].getWaterDepth() > 0.0)
      {
         if (newSpore[i] >= 0.1)
            this->cells_[i].setSporeLevel(newSpore[i]);
         else  // if sporeLevel<0.5 then change to 0.0
            this->cells_[i].setSporeLevel(0.0);
      }
      else
         this->cells_[i].setSporeLevel(0.0);
   }
   delete [] spore;
   delete [] newSpore;
}


//
// JM 05Dec2006 - Update spore levels in patches.
//
void Plate::updateLevels(unsigned nbLoops)
{
   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;
   
   // Update spore and charge levels in each patch
   for( unsigned i=0; i<nbCells; i++ )
   {
      if (this->cells_[i].getActiveRegion() == true)
      {
         this->cells_[i].updateSporeLevel();
         //this->cells_[i].updateNutLevel();
         //this->cells_[i].updateCharge();
      }
   }
}

//
//   JM 05Dec2006 - Send & receive nutrient, Ab, beta-lactamase and bacterial biomass levels 
//                  in buffer zones between plates
//
void Plate::sendAndReceive(int width, unsigned nbLoops)
{   
   MPI_Status status[NUM_ADJ_PLATES], statusBac[NUM_ADJ_PLATES];
   MPI_Request request[NUM_ADJ_PLATES], requestBac[NUM_ADJ_PLATES];
  
   double *sndBuf[NUM_ADJ_PLATES], *rcvBuf[NUM_ADJ_PLATES];
   float *sndBacBuf[NUM_ADJ_PLATES], *rcvBacBuf[NUM_ADJ_PLATES]; 
   
   unsigned numParams = 2;
   //cout << myRank << ": Sending" << endl;
   //                                                                                 //
   // ------------------------ Sending Information --------------------------------   //
   //
   for (unsigned plate_dir=1; plate_dir<NUM_ADJ_PLATES; plate_dir++)
   {
      if (adjacentPlate[plate_dir])
      {
         sndBuf[plate_dir] = new double[activeEdgeLength[plate_dir]*numParams];
         sndBacBuf[plate_dir] = new float[activeEdgeLength[plate_dir]];
         
         // Fill send buffers
         if (plate_dir == 1 || plate_dir == 3)
         {
            // Fill vertical buffers
            for (unsigned i=0; i<activeEdgeLength[plate_dir]; i++)
            {
               //sndBuf[plate_dir][(i*numParams)] = (*this)(i+iPos_parDiff[plate_dir], jPos_parDiff[plate_dir]).getNutrientLevel();
               sndBuf[plate_dir][(i*numParams)+1] = (*this)(i+iPos_parDiff[plate_dir], jPos_parDiff[plate_dir]).getSporeLevel();
               sndBacBuf[plate_dir][i] = (*this)(i+iPos_parDiff[plate_dir], jPos_parDiff[plate_dir]).getSporoBiomass();
            }
         }
         else
         {
            // Fill horizontal buffers
            for (unsigned j=0; j<activeEdgeLength[plate_dir]; j++)
            {
               //sndBuf[plate_dir][(j*numParams)] = (*this)(iPos_parDiff[plate_dir], j+jPos_parDiff[plate_dir]).getNutrientLevel();
               sndBuf[plate_dir][(j*numParams)+1] = (*this)(iPos_parDiff[plate_dir], j+jPos_parDiff[plate_dir]).getSporeLevel();
               sndBacBuf[plate_dir][j] = (*this)(iPos_parDiff[plate_dir], j+jPos_parDiff[plate_dir]).getSporoBiomass();
            }
         }

         // Non-blocking Send
         MPI_Issend(sndBuf[plate_dir], activeEdgeLength[plate_dir]*numParams, MPI_DOUBLE, dest[plate_dir], (10+plate_dir), MPI_COMM_WORLD, &request[plate_dir]);
         MPI_Issend(sndBacBuf[plate_dir], activeEdgeLength[plate_dir], MPI_FLOAT, dest[plate_dir], (20+plate_dir), MPI_COMM_WORLD, &requestBac[plate_dir]);
      }
   }
   
                                                                                    //
   MPI_Barrier(MPI_COMM_WORLD);

   //cout << myRank << ": Receiving" << endl;
   //                                                                                 //
   // ----------------------- Receiving information ------------------------------- //
   //
   for (unsigned plate_dir=1; plate_dir<NUM_ADJ_PLATES; plate_dir++)
   {
      int messTag=0;
      if (plate_dir == 1 || plate_dir == 2 || plate_dir == 5 || plate_dir == 6)
         messTag = plate_dir+2;
      else
         messTag = plate_dir-2;
   
      if (adjacentPlate[plate_dir])
      {
         rcvBuf[plate_dir] = new double[activeEdgeLength[plate_dir]*numParams];
         rcvBacBuf[plate_dir] = new float[activeEdgeLength[plate_dir]];

         // Receive the data via MPI
         MPI_Recv(rcvBuf[plate_dir], activeEdgeLength[plate_dir]*numParams, MPI_DOUBLE, dest[plate_dir], (10+messTag), MPI_COMM_WORLD, &status[plate_dir]);
         MPI_Recv(rcvBacBuf[plate_dir], activeEdgeLength[plate_dir], MPI_FLOAT, dest[plate_dir], (20+messTag), MPI_COMM_WORLD, &statusBac[plate_dir]);

          
         // Transfer receive buffer to patches
         
         // Fill vertical buffers
         if (plate_dir == 1 || plate_dir == 3)
         {
            unsigned jPositionDiff;
            if (plate_dir==1)
               jPositionDiff = engParamsPtr_->width_-1;
            else
               jPositionDiff = 0;
            
            for (unsigned i=0; i<activeEdgeLength[plate_dir]; i++)
            {                
               //(*this)(i+iPos_parDiff[plate_dir], jPositionDiff).setNutrientLevel(rcvBuf[plate_dir][(i*numParams)]);
               (*this)(i+iPos_parDiff[plate_dir], jPositionDiff).setSporeLevel(rcvBuf[plate_dir][(i*numParams)+1]);
               (*this)(i+iPos_parDiff[plate_dir], jPositionDiff).setSporoBiomass(rcvBacBuf[plate_dir][i]);
            }
         }
         // Fill horizontal buffers
         else
         {
            unsigned iPositionDiff;
            unsigned tempJPos;
            if (plate_dir==2 || plate_dir==5 || plate_dir==6)
               iPositionDiff = 0;
            else
               iPositionDiff = engParamsPtr_->height_-1;

            // Special case for corner communication
            if (plate_dir == 5 || plate_dir==8)
               tempJPos = jPos_parDiff[plate_dir]+1;
            else if (plate_dir==7 || plate_dir==6)
               tempJPos = jPos_parDiff[plate_dir]-1;
            else
               tempJPos = jPos_parDiff[plate_dir];
               
            for (unsigned j=0; j<activeEdgeLength[plate_dir]; j++)
            {                  
               //(*this)(iPositionDiff, j+tempJPos).setNutrientLevel(rcvBuf[plate_dir][(j*numParams)]);
               (*this)(iPositionDiff, j+tempJPos).setSporeLevel(rcvBuf[plate_dir][(j*numParams)+1]);
               (*this)(iPositionDiff, j+tempJPos).setSporoBiomass(rcvBacBuf[plate_dir][j]);
            }
         }
      }  
   }
      
   //cout << myRank << ": waiting" << endl;                                                                           

   // --------------------------- Waiting... -------------------------------------//
   for (unsigned plate_dir=1; plate_dir<NUM_ADJ_PLATES; plate_dir++)
   {
      if (adjacentPlate[plate_dir])
      {
         int messTag;
         if (plate_dir == 1 || plate_dir == 2 || plate_dir == 5 || plate_dir == 6)
            messTag = plate_dir+2;
         else
            messTag = plate_dir-2;
         MPI_Wait(&request[plate_dir], &status[messTag]);
         MPI_Wait(&requestBac[plate_dir], &statusBac[messTag]);
      }
   }

   //cout << myRank << ": deleting" << endl;  
   for (unsigned plate_dir=1; plate_dir<NUM_ADJ_PLATES; plate_dir++)
   {
      if (adjacentPlate[plate_dir])
      {
         delete [] sndBuf[plate_dir]; 
         delete [] sndBacBuf[plate_dir]; 
         delete [] rcvBuf[plate_dir]; 
         delete [] rcvBacBuf[plate_dir];
      }
   }
   //delete [] sndBuf; delete [] sndBacBuf; delete [] rcvBuf; delete [] rcvBacBuf;
}




//
// 24Apr2007 - Retrieve total bacterial biomass level in all active patches on plate.
//
double Plate::calcTotalBiomass()
{
   double totBiomass=0.0;
   
   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;
   
   for( unsigned i=0; i<nbCells; i++)
   {         
      if (this->cells_[i].getActiveRegion())
         totBiomass += this->cells_[i].calcSporoBiomass();
   }

   return totBiomass;
}


//
// 27Feb2007 - Retrieve total seed level in all active patches on plate.
//
double Plate::getTotalSporeLevel(unsigned mpiSize)
{
   double totSpores=0.0;

   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;

   if (mpiSize > 1)
   {
      for( unsigned i=0; i<nbCells; i++)
      {         
         if (this->cells_[i].getSubstrateType() >= 0 && this->cells_[i].getActiveRegion())
            totSpores += this->cells_[i].getSporeLevel();
      }
   }
   else
      for( unsigned i=0; i<nbCells; i++)
         totSpores += this->cells_[i].getSporeLevel();
   
   return totSpores;
}

//
// 06Feb2014 - Check to see if any spores present on plate
//
bool Plate::checkSporesPresent(unsigned mpiSize)
{
   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;
   for( unsigned i=0; i<nbCells; i++)
   {
      if (this->cells_[i].getSporeLevel() > 0.0)
         return true;
   }
   
   return false;
}




//
// JM 12MAr2007 - Check if patch is in active region of plate.
//
void Plate::setActivePatches()
{
   for( unsigned i=iStart; i<iStop; i++ )
      for( unsigned j=jStart; j<jStop; j++ )
      {
         if ((*this)(i,j).getSubstrateType() > -1)
            (*this)(i,j).setActive(true);
      }
}

//
//   JM 05Dec2006 - Send & receive bacteria from buffer zone between plates.
//
int Plate::parallelMovement(int size)
{
   int width = (int)sqrt((float)size);

   MPI_Status  statusFloats[NUM_ADJ_PLATES],
               statusBools[NUM_ADJ_PLATES], statusInts[NUM_ADJ_PLATES], 
               statusPos[NUM_ADJ_PLATES];

   MPI_Request requestFloats[NUM_ADJ_PLATES],
               requestBools[NUM_ADJ_PLATES], requestInts[NUM_ADJ_PLATES], 
               requestPos[NUM_ADJ_PLATES];

   // Organism Float traits - stock, seedStock.
   float *sndFloatsBuf[NUM_ADJ_PLATES], *rcvFloatsBuf[NUM_ADJ_PLATES];  
   // Organism boolean traits - expressSeed.
   // Stored as short integers for sending via MPI (no MPI_BOOL).
   short *sndBoolsBuf[NUM_ADJ_PLATES], *rcvBoolsBuf[NUM_ADJ_PLATES];
   // Organism Integer traits - seedlingID.
   int *sndIntsBuf[NUM_ADJ_PLATES], *rcvIntsBuf[NUM_ADJ_PLATES];
   // Indexes of patches to transfer bacteria.
   int *sndPos[NUM_ADJ_PLATES], *rcvPos[NUM_ADJ_PLATES];   
   unsigned patchSndIndex=0, patchRcvIndex=0, bufIndex=0;

   // Calculate size of communication buffers (i.e. no. of bacteria to transfer).
   calcSndBufferSize();
   calcRcvBufferSize(width);
   
   // Keep track of how many agents are sent/received from/to node
   int numSent = 0;
   int numReceived = 0;
   int netChange = 0;

   MPI_Barrier(MPI_COMM_WORLD);
   
   for (unsigned plate_dir=1; plate_dir<NUM_ADJ_PLATES; plate_dir++)
   {  
      if (adjacentPlate[plate_dir] && sndBufSize[plate_dir] > 0)
      {
         //cout << myRank << ": Sending agents..." << endl;
         sndPos[plate_dir] = new int[sndBufSize[plate_dir]];
         sndFloatsBuf[plate_dir] = new float[(sndBufSize[plate_dir]*2)];
         sndBoolsBuf[plate_dir] = new short[sndBufSize[plate_dir]];
         sndIntsBuf[plate_dir] = new int[sndBufSize[plate_dir]];
         unsigned numBacAdded;
         
         // Fill vertical send buffers
         if (plate_dir == 1 || plate_dir == 3)
         {
            // Start to fill buffer at co-ords (iPos_parMove, jPos_parMove)
            patchSndIndex = iPos_parMove[plate_dir]; 
            bufIndex = 0;
            while (bufIndex < sndBufSize[plate_dir] && patchSndIndex < (unsigned)engParamsPtr_->height_)
            {
               numBacAdded = 0;

               Agent *ptrOrganism=0;   
               ptrOrganism = (*this)(patchSndIndex, jPos_parMove[plate_dir]).getAgents();

               // Add bacteria to send buffer
               if (ptrOrganism)
                  numBacAdded = addToBuffer(sndFloatsBuf[plate_dir], sndBoolsBuf[plate_dir], sndIntsBuf[plate_dir], ptrOrganism, bufIndex, sndBufSize[plate_dir] );
               
               // Get patch index of bacteria.
               if (numBacAdded)
               {
                  for (unsigned i=0; i<numBacAdded; i++)
                     sndPos[plate_dir][bufIndex+i] = patchSndIndex;
    
                  bufIndex += numBacAdded;
                  numSent += numBacAdded;
               }
               patchSndIndex++;
            }  // end while
         }
         // Fill horizontal send buffers
         else
         {
            // Start to fill buffer at co-ords (iPos_parMove, jPos_parMove) 
            patchSndIndex = jPos_parMove[plate_dir];

            // if not sending to diagonally adjacent plates
            if ( !(plate_dir==5 || plate_dir==6 || plate_dir==7 || plate_dir==8) )
               patchRcvIndex = patchSndIndex;
            else if (plate_dir==5 || plate_dir==8)
               patchRcvIndex = 1;
            else
               patchRcvIndex = engParamsPtr_->width_-2;
            
            bufIndex = 0;
            while (bufIndex < sndBufSize[plate_dir] && patchSndIndex < (unsigned)engParamsPtr_->width_)
            {
               numBacAdded = 0;

               Agent *ptrOrganism=0;   
               ptrOrganism = (*this)(iPos_parMove[plate_dir], patchSndIndex).getAgents();

               // Add bacteria to send buffer
               if (ptrOrganism)
                  numBacAdded = addToBuffer(sndFloatsBuf[plate_dir], sndBoolsBuf[plate_dir], sndIntsBuf[plate_dir], ptrOrganism, bufIndex, sndBufSize[plate_dir] );
               
               // Get patch index of bacteria.
               if (numBacAdded)
               {  
                  for (unsigned i=0; i<numBacAdded; i++)
                     sndPos[plate_dir][bufIndex+i] = patchRcvIndex;                        

                  bufIndex += numBacAdded;
                  numSent += numBacAdded;
               }
               patchSndIndex++;
               patchRcvIndex++;
            }  // end while
         }
         
         // Non-blocking Send
      MPI_Issend(sndPos[plate_dir], sndBufSize[plate_dir], MPI_INT, dest[plate_dir], 1000+plate_dir, MPI_COMM_WORLD, &requestPos[plate_dir]);
      MPI_Issend(sndFloatsBuf[plate_dir], (sndBufSize[plate_dir]*2), MPI_FLOAT, dest[plate_dir], 3000+plate_dir, MPI_COMM_WORLD, &requestFloats[plate_dir]);
      MPI_Issend(sndBoolsBuf[plate_dir], sndBufSize[plate_dir], MPI_SHORT, dest[plate_dir], 4000+plate_dir, MPI_COMM_WORLD, &requestBools[plate_dir]);
      MPI_Issend(sndIntsBuf[plate_dir], (sndBufSize[plate_dir]), MPI_INT, dest[plate_dir], 5000+plate_dir, MPI_COMM_WORLD, &requestInts[plate_dir]);
         
      }  // end if (adjacentPlate[plate_dir]...
      
   }  // end for loop
    
   
   //                                                                                    //
   // --------------------- Receiving information -------------------------------- //
   //                                                                                       //  
   for (unsigned plate_dir=1; plate_dir<NUM_ADJ_PLATES; plate_dir++)
   {
      if (adjacentPlate[plate_dir] && rcvBufSize[plate_dir] > 0)
      {
         //cout << myRank << ": Receiving agents..." << endl;
         int messTag;
         if (plate_dir == 1 || plate_dir == 2 || plate_dir == 5 || plate_dir == 6)
            messTag = plate_dir+2;
         else
            messTag = plate_dir-2;
            
         rcvPos[plate_dir] = new int[rcvBufSize[plate_dir]];
         rcvFloatsBuf[plate_dir] = new float[rcvBufSize[plate_dir]*2];
         rcvBoolsBuf[plate_dir] = new short[rcvBufSize[plate_dir]];
         rcvIntsBuf[plate_dir] = new int[rcvBufSize[plate_dir]];    
      
         MPI_Recv(rcvPos[plate_dir], rcvBufSize[plate_dir], MPI_INT, dest[plate_dir], 1000+messTag, MPI_COMM_WORLD, &statusPos[plate_dir]);
         MPI_Recv(rcvFloatsBuf[plate_dir], (rcvBufSize[plate_dir]*2), MPI_FLOAT, dest[plate_dir], 3000+messTag, MPI_COMM_WORLD, &statusFloats[plate_dir]);
         MPI_Recv(rcvBoolsBuf[plate_dir], rcvBufSize[plate_dir], MPI_SHORT, dest[plate_dir], 4000+messTag, MPI_COMM_WORLD, &statusBools[plate_dir]);
         MPI_Recv(rcvIntsBuf[plate_dir], (rcvBufSize[plate_dir]), MPI_INT, dest[plate_dir], 5000+messTag, MPI_COMM_WORLD, &statusInts[plate_dir]);
         

         // Fill from vertical receive buffers
         if (plate_dir == 1 || plate_dir == 3)
         {
            unsigned jPosition=1;
            if (plate_dir == 1)
               jPosition = engParamsPtr_->width_-2;

            // Transfer receive buffer to bacteria
            for (unsigned i=0; i<rcvBufSize[plate_dir]; i++) 
            {
               (*this)(rcvPos[plate_dir][i], jPosition).newOrganismMPI( 
                     rcvFloatsBuf[plate_dir][i*2],
                     rcvBoolsBuf[plate_dir][i]/*, rcvBoolsBuf[plate_dir][(i*3)+1], 
                     rcvBoolsBuf[plate_dir][(i*3)+2]*/, rcvIntsBuf[plate_dir][i] );
               numReceived++;
            }
         }
         // Fill from horizontal receive buffers
         else
         {
            unsigned iPosition=1;
            unsigned jPosition;
            if (plate_dir == 8 || plate_dir == 4 || plate_dir == 7)
               iPosition = engParamsPtr_->height_-2; 
            
            // Transfer receive buffer to bacteria
            for (unsigned i=0; i<rcvBufSize[plate_dir]; i++) 
            {   
               // Special case for corner movement
               if ( plate_dir == 5 || plate_dir == 8 )                     
                  jPosition = engParamsPtr_->width_-2;
               else if ( plate_dir == 6 || plate_dir == 7 )
                  jPosition = 1;
               else
                  jPosition = rcvPos[plate_dir][i];
                        
               (*this)(iPosition, jPosition).newOrganismMPI(
                     rcvFloatsBuf[plate_dir][i*2],
                     rcvBoolsBuf[plate_dir][i]/*, rcvBoolsBuf[plate_dir][(i*3)+1], 
                     rcvBoolsBuf[plate_dir][(i*3)+2]*/, rcvIntsBuf[plate_dir][i] );
               numReceived++;
            }
         }
  
      }
   }  // end for loop
   
   
   // --------------------------- Waiting... -------------------------------------//
   for (unsigned plate_dir=1; plate_dir<NUM_ADJ_PLATES; plate_dir++)
   {
      if (adjacentPlate[plate_dir] && sndBufSize[plate_dir] > 0) 
      {
         int messTag;
         if (plate_dir == 1 || plate_dir == 2 || plate_dir == 5 || plate_dir == 6)
            messTag = plate_dir+2;
         else
            messTag = plate_dir-2;
            
         MPI_Wait(&requestPos[plate_dir], &statusPos[messTag]);
         MPI_Wait(&requestFloats[plate_dir], &statusFloats[messTag]);
         MPI_Wait(&requestBools[plate_dir], &statusBools[messTag]);
         MPI_Wait(&requestInts[plate_dir], &statusInts[messTag]);
         //cout << myRank << ": Done" << endl;
      }
   }  // end for loop

   MPI_Barrier(MPI_COMM_WORLD);

   // Remove all organisms from buffer zones
   for (unsigned plate_dir=1; plate_dir<(NUM_ADJ_PLATES+1)/2; plate_dir++)
   {
      if (adjacentPlate[plate_dir])// && sndBufSize[plate_dir] > 0)
      {
         // Delete vertical send buffers
         if (plate_dir == 1 || plate_dir == 3)
         {
            unsigned jPosition=0;
            if (plate_dir == 1)
               jPosition = engParamsPtr_->width_-1;
               
            for (unsigned i=0; i<(unsigned)engParamsPtr_->height_; i++) 
            {
               Agent *ptrOrganism = (*this)(i, jPosition).getAgents();
               if (ptrOrganism)
                  (*this)(i, jPosition).removeAllOrganisms();
            }
         }
         else  // delete horizontal send buffers
         {
            unsigned iPosition=0;
            if (plate_dir == 4) //|| plate_dir == 3 || plate_dir == 5)
               iPosition = engParamsPtr_->height_-1;
               
            for (unsigned j=0; j<(unsigned)engParamsPtr_->width_; j++) 
            {
               Agent *ptrOrganism = (*this)(iPosition, j).getAgents();
               if (ptrOrganism)
                  (*this)(iPosition, j).removeAllOrganisms();
            }
         }
         //if (sndBufSize[plate_dir] > 0)
            //cout << myRank << ": Deleted bufs" << endl;
      }
   }
   
   MPI_Barrier(MPI_COMM_WORLD);

   // Clear memory...
   for (unsigned plate_dir=1; plate_dir<NUM_ADJ_PLATES; plate_dir++)
   {
      if (adjacentPlate[plate_dir] && sndBufSize[plate_dir] > 0) 
      {
         delete [] sndPos[plate_dir];
         delete [] sndFloatsBuf[plate_dir];
         delete [] sndBoolsBuf[plate_dir];
         delete [] sndIntsBuf[plate_dir];
      }
      else if (adjacentPlate[plate_dir] && rcvBufSize[plate_dir] > 0)
      {
         delete [] rcvPos[plate_dir];
         delete [] rcvFloatsBuf[plate_dir];
         delete [] rcvBoolsBuf[plate_dir];
         delete [] rcvIntsBuf[plate_dir];
      }
   }  // end for loop
   
   netChange = (numReceived - numSent); 
   //if (netChange != 0)
   //  cout << "netChange = " << netChange << endl;
   
   // return net change in number of agents as a result of sending/receiving
   return netChange;
}


//
// 12Mar2007 - Add organisms, from specified patch, to buffer for sending across nodes (MPI)
//
unsigned Plate::addToBuffer(float* sndFloatsBuf, short* sndBoolsBuf, int* sndIntsBuf, Agent* ptrOrganism, unsigned bufIndex, unsigned bufferSize)
{
   unsigned numBacAdded=0;
   // Go through list of organisms in patch
   while (ptrOrganism)
   {
      if (bufIndex >= bufferSize)
         cout << "Error: Plate::addToBuffer() - " << "send buffer(" << bufferSize 
               << ") too small for bufIndex(" << bufIndex << ")" << endl;

      if( ptrOrganism->getType() == ORGANISM  && ptrOrganism->isAlive())
      {
         // Input into buffers:
         // Get organism floating point traits.
         sndFloatsBuf[bufIndex*2] = ((Organism*)ptrOrganism)->getStock();
         //sndFloatsBuf[(bufIndex*2)+1] = ((Organism*)ptrOrganism)->getSeedStock();

         // Get organism boolean traits.
         sndBoolsBuf[bufIndex] = (short)((Organism*)ptrOrganism)->getReleaseSpore();
         //sndBoolsBuf[(bufIndex*3)+1] = (short)((Organism*)ptrOrganism)->getSOSResponse();
         //sndBoolsBuf[(bufIndex*3)+2] = (short)((Organism*)ptrOrganism)->getInitAbDeath();
         
         // Get organism integer traits.
         sndIntsBuf[bufIndex] = ((Organism*)ptrOrganism)->getSeedlingID();

         numBacAdded++;
         bufIndex++;
      }
      ptrOrganism = ptrOrganism->getNext();
   }
   return numBacAdded;
}

//
// JM 13Mar2007 - Determine size of send buffers (i.e. number of organisms).
//
void Plate::calcSndBufferSize()
{
   // Reset size of communication buffers to zero first
   //lBufSize=0, rBufSize=0, uBufSize=0, dBufSize=0;
   for (unsigned i=0; i<NUM_ADJ_PLATES; i++)
   {
      sndBufSize[i] = 0; rcvBufSize[i]=0;
   }
     
   for (unsigned plate_dir=1; plate_dir<NUM_ADJ_PLATES; plate_dir++)
   {
      switch(plate_dir)
      {
         case 1:
            for(unsigned i=0; i<engParamsPtr_->height_; i++)
               sndBufSize[plate_dir] += (*this)(i, engParamsPtr_->width_-1).countNumOrganism();
            break;
         case 2:
            for(unsigned j=0; j<engParamsPtr_->width_; j++)
               sndBufSize[plate_dir] += (*this)(0, j).countNumOrganism();
            break;
         case 3:
            for(unsigned i=0; i<engParamsPtr_->height_; i++)
               sndBufSize[plate_dir] += (*this)(i, 0).countNumOrganism();
            break;
         case 4:
            for(unsigned j=0; j<engParamsPtr_->width_; j++)
               sndBufSize[plate_dir] += (*this)(engParamsPtr_->height_-1, j).countNumOrganism();
            break;
         case 5:
            sndBufSize[plate_dir] += (*this)(0, engParamsPtr_->width_-1).countNumOrganism();
            break;
         case 6:
            sndBufSize[plate_dir] += (*this)(0, 0).countNumOrganism();
            break;
         case 7:
            sndBufSize[plate_dir] += (*this)(engParamsPtr_->height_-1, 0).countNumOrganism();
            break;
         case 8:
            sndBufSize[plate_dir] += (*this)(engParamsPtr_->height_-1, engParamsPtr_->width_-1).countNumOrganism();
            break;
         default:
            cout << "*Warning* Plate::calcSndBufferSize(): plate_dir out of range" << endl;
            sndBufSize[plate_dir] = 0;
            break;
      }
   }
}



//
// JM 13Mar2007 - Determine size of receive buffers.
//
void Plate::calcRcvBufferSize(int width)
{
   MPI_Status statusBufSize[NUM_ADJ_PLATES];
   MPI_Request requestBufSize[NUM_ADJ_PLATES];
         
   // Send buffer size to all adjacent plates    
   for (unsigned plate_dir=1; plate_dir<NUM_ADJ_PLATES; plate_dir++)
   {
      if(adjacentPlate[plate_dir])
         MPI_Issend(&sndBufSize[plate_dir], 1, MPI_INT, dest[plate_dir], 100+plate_dir, MPI_COMM_WORLD, &requestBufSize[plate_dir]);
   }
   MPI_Barrier(MPI_COMM_WORLD);
   
   // Receive buffer sizes from all adjacent plates
   int messTag;
   for (unsigned plate_dir=1; plate_dir<NUM_ADJ_PLATES; plate_dir++)
   {
      if (plate_dir == 1 || plate_dir == 2 || plate_dir == 5 || plate_dir == 6)
         messTag = plate_dir+2;
      else
         messTag = plate_dir-2;
         
      if(adjacentPlate[plate_dir])
         MPI_Recv(&rcvBufSize[plate_dir], 1, MPI_INT, dest[plate_dir], 100+messTag, MPI_COMM_WORLD, &statusBufSize[plate_dir]);
   }
   
   // Waiting
   // Send buffer size to all adjacent plates    
   for (unsigned plate_dir=1; plate_dir<NUM_ADJ_PLATES; plate_dir++)
   {
      if(adjacentPlate[plate_dir])
         MPI_Wait(&requestBufSize[plate_dir], &statusBufSize[plate_dir]);
   }                
}


//
// 12Aug14 - Count the total number of sporophytes (approx. >5cm)
//
unsigned Plate::countSporophytes(float minStock, float maxStock)
{
   unsigned count=0;
   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;
   
   for( unsigned i=0; i<nbCells; i++ )
   {
      // count organisms with stock_ > 0.1
      count += this->cells_[i].countNumSporophytes(minStock, maxStock);
   }
   return count;
}

//
// 12Aug14 - Count the total number of new sporophyte recruits
//          i.e. >~5cm and only discovered in last month
//
unsigned Plate::countNewRecruits(float minStock, unsigned age)
{
   unsigned count=0;
   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;
   
   for( unsigned i=0; i<nbCells; i++ )
   {
      // count organisms with stock_ > 0.1
      count += this->cells_[i].countNewRecruits(minStock, age);
   }
   return count;
}


unsigned Plate::countNewSporophytes(float minStock, unsigned age)
{
   unsigned count=0;
   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;
   
   for( unsigned i=0; i<nbCells; i++ )
   {
      // count organisms with stock_ > 0.1
      count += this->cells_[i].countNewSporophytes(minStock, age);
   }
   return count;
}

//
// 12Aug14 - Count the total number of gametophytes (or just mature ones)
//
unsigned Plate::countGametophytes(bool onlyMature)
{
   unsigned count=0;
   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;
   
   for( unsigned i=0; i<nbCells; i++ )
   {
      if (onlyMature)
         count += this->cells_[i].countMatureGametophytes(); 
      else  
         count += this->cells_[i].countNumGametophytes();
   }
   return count;
}




//
// JM 19Nov13 - Update water tide level 
//            - uses Sine wave as rough *approximation* of tidal wave
//
void Plate::updateWaterLevel(unsigned nbLoops)
{
   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;
   
   // Assuming 1 loop = 1 hour
   // Assuming semi-diurnal tidal cycle is 24h 50m
   unsigned cycleLength = 24.83;
   unsigned currentStep = nbLoops % cycleLength;
   // water level varies periodically according to sine wave (amplitude = 1.85 metres)
   waterLevel = 1.85 * sin( ( (2*M_PI * currentStep) / cycleLength ) * 4 );
   
   // Add 4.44m to match data from Aber Wrac'h (maree.info website, 09Jan2014)
   waterLevel += 4.44;
   
   // update current water level in all patches
   for( unsigned i=0; i<nbCells; i++ )
   {
      if (this->cells_[i].getSubstrateType() >= 0)
         this->cells_[i].updateWaterDepth(waterLevel);
   }
}

//
// JM 13Jan14 - Update water temperature and solar radiation (global values for all patches)
//
void Plate::updateWaterSolar(unsigned nbLoops)
{
   unsigned nbMonths = 12;
   unsigned yearLength = 8760; // compete orbit of sun = 8766.15 h (To Do: add leap years)
   unsigned currHour = nbLoops % yearLength;
   //unsigned currYear = nbLoops / yearLength;
   
   // temperature input data (Voisin thesis) is over four year cycle
  // unsigned tempCycle = 4;
   //if (currYear > tempCycle-1)
   //   currYear = currYear % tempCycle;

   currMonth = 0;
     
   // Find out what month it is so can look up temperature/solar data for months
   for (int i=nbMonths-1; i>=0; i--)
   {
      if (currHour >= stMonthsHr[i])
      {
         currMonth = i;
         // dayMod - modify based on day in month        
         //dayMod = (float)(currHour - stMonthsHr[i])/monthLength[i];
         break;
      }
   }
   
   // Current calendar day [0...364]
   currDay = (unsigned)(currHour/24.0);
   if (currDay >= 365)
   {
      cout << nbLoops << ": currDay = " << currDay << ", currHour = " << currHour << endl;
      currDay = 364;
   }

   // update current water temperature
   float oldTemp = currWaterTemp;
   currWaterTemp = dailyWaterTemp[currDay];
   currWaterTemp += sampleGaussian(currWaterTemp*0.005, 1);
   currDayLength = dayLengths[currDay];
   
   // convert solar radiation from MJ m^-2 day^-1 to MJ m^-2 h^-1    
   currSolarRad = dailySolarRad[currDay]/currDayLength;
   
   if (currWaterTemp < oldTemp*0.098)
      cout << currDay << ", " << currMonth << ": " << dailyWaterTemp[currDay] << ", " << daysPerMonth[currMonth] << ", " << ", " << currWaterTemp << endl;
}


//
// JM 20Jan14 - Update the annual sea water temperature curve based on inputted averages
//
void Plate::updateDailyWaterSolar(unsigned nbLoops, float *waterTempIn, float *solarRadIn)
{
   unsigned yearLength = 8766; // complete orbit of sun = 8766.15 h
   
   if (nbLoops % yearLength == 0)
   {
      for(unsigned i=0; i<365; i++)
      {
         dailyWaterTemp[i] = waterTempIn[i] + sampleGaussian(waterTempIn[i]*0.01, 2);
         dailySolarRad[i] = solarRadIn[i] + sampleGaussian(solarRadIn[i]*0.01, 2);
      }
   }
}


//
// JM 18APR11 - LBE: Run Lattice Boltzmann methods
//
/*void Plate::moveLBE(int sizeMPI, unsigned nbLoops)
{
   unsigned width = (int)sqrt((float)sizeMPI);
   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;
   //Patch *ptrPatTest1=NULL; Patch *ptrPatTest2=NULL;
   //ptrPatTest1 = &(this->smartAccess(11,31));
   //ptrPatTest2 = &(this->smartAccess(20,31));

   // Movement
   for( unsigned i=0; i<nbCells; i++ )
   {
      if (this->cells_[i].getActiveRegion() == true || this->cells_[i].getSubstrateType() == -1)
         this->cells_[i].moveLBE(engParamsPtr_, nbLoops);
   }
   
   for( unsigned i=0; i<nbCells; i++ )
   {
      if (this->cells_[i].getActiveRegion() == true || this->cells_[i].getSubstrateType() == -1)
         this->cells_[i].updateF(nbLoops);
   }


   // Communicate with adjacent plates.
   if (sizeMPI > 1)
      sendAndReceiveLBE(width, nbLoops);
}

void Plate::obstacleLBE(int sizeMPI, unsigned nbLoops)
{
   //int width = (int)sqrt((float)sizeMPI);
   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;

   if ( (int)nbLoops >= baitParamsPtr_->loopAddObst && 
        ((int)nbLoops < baitParamsPtr_->loopRemoveObst || baitParamsPtr_->loopRemoveObst == -1) )
   {
      for( unsigned i=0; i<nbCells; i++ )
      {
         if (this->cells_[i].getSubstrateType() == -1)
            this->cells_[i].obstacle(engParamsPtr_, nbLoops);
      }
   
      for( unsigned i=0; i<nbCells; i++ )
      {
         if (this->cells_[i].getSubstrateType() == -1)
            this->cells_[i].updateObst(nbLoops);
      }
   }
}

//
// JM 28APR11 - LBE: Code for inserting obstacle into flow
//
void Plate::addObstacleLBE(int sizeMPI)
{
   unsigned stPosI, stPosJ;

   stPosI = (engParamsPtr_->height_ - baitParamsPtr_->obstYPos);
   stPosJ = baitParamsPtr_->obstXPos-1;

   cout << "stPosI = " << stPosI << " to " << (stPosI+baitParamsPtr_->obstLength) <<  endl;
   Patch *ptrPatch=NULL;
        
    
   if ((int)myRank == baitParamsPtr_->rankAddObst || baitParamsPtr_->rankAddObst == -1)
   {
      cout << "(" << myRank << ") Obst (i,j): ";
      for (unsigned i=stPosI; i<=(stPosI+baitParamsPtr_->obstLength); i++)
      {
         for(unsigned j=stPosJ; j<(stPosJ+baitParamsPtr_->obstWidth); j++)
         {
            ptrPatch = &(this->smartAccess(i,j));
            ptrPatch->setIsObstacle(true);
            ptrPatch->setActive(false);
            if ( i == stPosI )
               ptrPatch->setIsTopObst(true); 
            else if ( i == (stPosI+baitParamsPtr_->obstLength) )
               ptrPatch->setIsBottObst(true);
            cout << "(" << i << "," << j << ") ";     
         }  
      }
      cout << endl;
   }  
}

//
// JM 28APR11 - LBE: Code for removing obstacle from flow
//
void Plate::removeObstacleLBE()
{
   if ((int)myRank == baitParamsPtr_->rankAddObst || baitParamsPtr_->rankAddObst == -1)
   {
      cout << myRank << ": Removing obstacle" << endl;
   
      unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;
      for( unsigned i=0; i<nbCells; i++ )
         if (this->cells_[i].getSubstrateType() == -1)
         {
            this->cells_[i].setIsObstacle(false);
            this->cells_[i].setIsTopObst(false);
            this->cells_[i].setIsBottObst(false);
            this->cells_[i].setActive(true);
         }
   } 
}



//
// JM 18APR11 - LBE: re-calculate hydrodynamic parameters based on updated values after movement
//              call hydrovar, equili, collis and force subroutines
//
void Plate::methodsLBE(unsigned nbLoops)
{
   //Patch *ptrPatTest1=NULL; Patch *ptrPatTest2=NULL;
   //ptrPatTest1 = &(this->smartAccess(11,31));
   //ptrPatTest2 = &(this->smartAccess(20,31));

   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;
   for( unsigned i=0; i<nbCells; i++ )
   {
      this->cells_[i].hydrovar(nbLoops);
      this->cells_[i].equili();
      this->cells_[i].collis(omega, nbLoops);
      this->cells_[i].forceLB(force);
   }
}



//
// JM 20APR11 - LBE: MPI commuication of lattice boltzmann frequencies (f[]) after movement
//
void Plate::sendAndReceiveLBE(int width, unsigned nbLoops)
{   
   MPI_Status status[NUM_ADJ_PLATES];
   MPI_Request request[NUM_ADJ_PLATES];
  
   double *sndBuf[NUM_ADJ_PLATES], *rcvBuf[NUM_ADJ_PLATES];   
   int numParams = 9;

   //                                                                                 //
   // ------------------------ Sending Information --------------------------------   //
   //
   for (unsigned plate_dir=1; plate_dir<NUM_ADJ_PLATES; plate_dir++)
   {
      if (adjacentPlate[plate_dir])
      {
         sndBuf[plate_dir] = new double[activeEdgeLength[plate_dir]*numParams];
         
         // Fill send buffers:
         // Fill vertical buffers
         if (plate_dir == 1 || plate_dir == 3)
         {
            for (unsigned i=0; i<activeEdgeLength[plate_dir]; i++)
               for (int k=0; k<numParams; k++)
                  sndBuf[plate_dir][(i*numParams)+k] = (*this)(i+iPos_parDiff[plate_dir], jPos_parDiff[plate_dir]).getF(k);
         }
         // Fill horizontal buffers
         else
         {
            for (unsigned j=0; j<activeEdgeLength[plate_dir]; j++)
               for (int k=0; k<numParams; k++)
                  sndBuf[plate_dir][(j*numParams)+k] = (*this)(iPos_parDiff[plate_dir], j+jPos_parDiff[plate_dir]).getF(k);
         }

         // Non-blocking Send
         MPI_Issend(sndBuf[plate_dir], activeEdgeLength[plate_dir]*numParams, MPI_DOUBLE, dest[plate_dir], (100+plate_dir), MPI_COMM_WORLD, &request[plate_dir]);
      }
   }
                                                                                       //
   //MPI_Barrier(MPI_COMM_WORLD);

   //                                                                                 //
   // ----------------------- Receiving information ------------------------------- //
   //
   for (unsigned plate_dir=1; plate_dir<NUM_ADJ_PLATES; plate_dir++)
   {
      int messTag=0;
      if (plate_dir == 1 || plate_dir == 2 || plate_dir == 5 || plate_dir == 6)
         messTag = plate_dir+2;
      else
         messTag = plate_dir-2;
   
      if (adjacentPlate[plate_dir])
      {
         rcvBuf[plate_dir] = new double[activeEdgeLength[plate_dir]*numParams];

         // Receive the data via MPI
         MPI_Recv(rcvBuf[plate_dir], activeEdgeLength[plate_dir]*numParams, MPI_DOUBLE, dest[plate_dir], (100+messTag), MPI_COMM_WORLD, &status[plate_dir]);
          
         // Transfer receive buffer to patches       
         // Fill vertical buffers
         if (plate_dir == 1 || plate_dir == 3)
         {
            unsigned jPositionDiff;
            if (plate_dir==1)
               jPositionDiff = engParamsPtr_->width_-1;
            else
               jPositionDiff = 0;
            
            for (unsigned i=0; i<activeEdgeLength[plate_dir]; i++)
               for (int k=0; k<numParams; k++)              
                  (*this)(i+iPos_parDiff[plate_dir], jPositionDiff).setF(rcvBuf[plate_dir][(i*numParams)+k], k);
         }
         // Fill horizontal buffers
         else
         {
            unsigned iPositionDiff;
            unsigned tempJPos;
            if (plate_dir==2 || plate_dir==5 || plate_dir==6)
               iPositionDiff = 0;
            else
               iPositionDiff = engParamsPtr_->height_-1;

            // Special case for corner communication
            if (plate_dir == 5 || plate_dir==8)
               tempJPos = jPos_parDiff[plate_dir]+1;
            else if (plate_dir==7 || plate_dir==6)
               tempJPos = jPos_parDiff[plate_dir]-1;
            else
               tempJPos = jPos_parDiff[plate_dir];
               
            for (unsigned j=0; j<activeEdgeLength[plate_dir]; j++)
               for (int k=0; k<numParams; k++)        
                  (*this)(iPositionDiff, j+tempJPos).setF(rcvBuf[plate_dir][(j*numParams)+k], k);
         }
      }  
   }
                                                                                 
   //cout << "waiting...";
   // --------------------------- Waiting... -------------------------------------//
   for (unsigned plate_dir=1; plate_dir<NUM_ADJ_PLATES; plate_dir++)
   {
      if (adjacentPlate[plate_dir])
      {
         int messTag;
         if (plate_dir == 1 || plate_dir == 2 || plate_dir == 5 || plate_dir == 6)
            messTag = plate_dir+2;
         else
            messTag = plate_dir-2;
         MPI_Wait(&request[plate_dir], &status[messTag]);
      }
   }

   //cout << "deleting...";
   for (unsigned plate_dir=1; plate_dir<NUM_ADJ_PLATES; plate_dir++)
   {
      if (adjacentPlate[plate_dir])
      {
         delete [] sndBuf[plate_dir]; 
         delete [] rcvBuf[plate_dir]; 
      }
   }
}


//
// JM 22APR11 - LBE - Diagnostics subroutine: get global average density and velocities (u and v)
//
void Plate::diag0D(unsigned nbLoops)//, ofstream *outputFile2)
{
   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;
   unsigned numActiveCells = 0;
   umoy=0.0, vmoy=0.0, avgDensity=0.0;
   int numDir = Patch::getNumDir();
   for( unsigned i=0; i<nbCells; i++ )
   {
      if (this->cells_[i].getActiveRegion() == true || this->cells_[i].getSubstrateType() == -1)
      {
         numActiveCells++;
         avgDensity += this->cells_[i].getRho_density();
         umoy += this->cells_[i].getU_vel();
         vmoy += this->cells_[i].getV_vel();
      }   
   }

   avgDensity = avgDensity / numActiveCells / numDir;
   umoy = umoy / numActiveCells;
   vmoy = vmoy / numActiveCells;
   
}

//
// JM 25APR11 - LBE - Config, print map of fluid velocities (u, v) for GNUPlot
//
void Plate::configLBE(unsigned nbLoops, ofstream *outputFiles[], int numOutFiles, unsigned index)
{
   double vor=0.0;
   int height = engParamsPtr_->height_;
   int width = engParamsPtr_->width_;
   Patch *ptrPatch=NULL, *ptrPatchNorth=NULL, *ptrPatchSouth=NULL,
         *ptrPatchWest=NULL, *ptrPatchEast=NULL;

   for (int i=height-1; i>=0; i--)
   {
      for(int j=0; j<width; j++)
      {
         ptrPatch = &(this->smartAccess(i,j));
         ptrPatchNorth = ptrPatch->getNeighbour(N, false);
         ptrPatchSouth = ptrPatch->getNeighbour(S, false);
         ptrPatchWest = ptrPatch->getNeighbour(W, false);
         ptrPatchEast = ptrPatch->getNeighbour(E, false);
      
         if ( ptrPatch->getActiveRegion() || (i >0 && i<height-1 && j>0 && j<width-1) 
               || ptrPatch->getSubstrateType() == -1 )
         {
            vor = ptrPatchNorth->getU_vel() - ptrPatchSouth->getU_vel() 
                  + ptrPatchEast->getV_vel() - ptrPatchWest->getV_vel();
            
            (*outputFiles[index]) << (j+1) << " " << (height-i) << " " << ptrPatch->getU_vel() << " " 
            << ptrPatch->getV_vel() << " " << vor << endl;
         }
      }
      (*outputFiles[index]) << endl;
   }
}

//
// JM 26APR11 - LBE - Profile: print various output files
//
void Plate::profileLBE(unsigned nbLoops, ofstream *outputFiles2[], int numOutFiles)
{
   unsigned height = engParamsPtr_->height_;
   unsigned width = engParamsPtr_->width_;
   Patch *ptrPatch=NULL, *ptrPatchL=NULL, *ptrPatchMid=NULL, *ptrPatchR=NULL, *ptrPatchProbe;

   for (int i=height-1; i>=0; i--)
   {
      ptrPatchL = &( this->smartAccess(i, (width/4)-1) );
      ptrPatchMid = &( this->smartAccess(i, (width/2)-1) );
      ptrPatchR = &( this->smartAccess(i, (3*width/4)-1) );
      
      if (ptrPatchL->getActiveRegion() == true || ptrPatchL->getSubstrateType() == -1)
      {
         // file .1 -> U velocities (.uy)
         (*outputFiles2[0]) << (height-i) << " " 
                     << ptrPatchL->getU_vel() << " " 
                     << ptrPatchMid->getU_vel() << " " 
                     << ptrPatchR->getU_vel() << endl;
         // file .2 -> V velocities (.vy)          
         (*outputFiles2[1]) << i  << " " << " " 
                     << ptrPatchL->getV_vel() << " " 
                     << ptrPatchMid->getV_vel() << " " 
                     << ptrPatchR->getV_vel() << endl;
         // file .3 -> Rho (.rhoy)            
         (*outputFiles2[2]) << (height-i) << " " 
                     << ptrPatchMid->getRho_density() << endl;
      }   
   }
   // need gap of two lines for producing GNUPlot movie with 'movie.gnu'
   (*outputFiles2[0]) << "\n" << endl;
   (*outputFiles2[1]) << "\n" << endl;
   (*outputFiles2[2]) << "\n" << endl;
   
   for (unsigned j=0; j<width; j++)
   {
      ptrPatch = &( this->smartAccess( (height/2)-1, j) );
      if (ptrPatch->getActiveRegion() == true || ptrPatch->getSubstrateType() == -1)
      {
         // file .4 -> u and v velocities at midway (.uvx)
         (*outputFiles2[3]) << (j+1) << " " << ptrPatch->getU_vel() << " " << ptrPatch->getV_vel() << endl;
         // file .5 -> f[1] and f[2] (.pop)
         (*outputFiles2[4]) << (j+1) << " " << ptrPatch->getF(1) << " " << ptrPatch->getF(3) << endl;
      }
   }
   (*outputFiles2[3]) << "\n" << endl;
   (*outputFiles2[4]) << "\n" << endl;
   
   // Sinusoidal Wave Decay test:
   // when initialise u_vel values with Sine function waves expected to decay exponentially over time:
   // e^(-visc*k^2*t)
   // visc (nu) = 1/3*(1/omega-0.5)
   // wave number k = 2*PI/height
   // t = loop number
   double expDecay = exp(-visc * pow((2*M_PI/height),2) * nbLoops );  
   double u_0 = sin( ( (2*M_PI * (height/4)) / height ) * 1 );
   
   ptrPatchProbe = &( this->smartAccess( (height/4)-1, (width/2)-1 ) );
   (*outputFiles2[5]) << nbLoops << " " << ptrPatchProbe->getU_vel()/u_0  << " " << expDecay << endl;  
   
}


//
// JM 10MAY11 - LBE: Update all lattice patches according to Lattice Boltzmann dynamics
//
void Plate::allCellsLBE(unsigned nbLoops)
{
   //double nutBefore = getTotalNutLevel();
   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_;

   // For all cells...
   for( unsigned i=0; i<nbCells; i++ )
      if (this->cells_[i].getActiveRegion() == true)
         this->cells_[i].fluidFlowLBE(nbLoops);
   
   for( unsigned i=0; i<nbCells; i++ )
      if (this->cells_[i].getActiveRegion() == true)
         this->cells_[i].updateSporeLevel();
}


//
// JM 19May11 - LBE: Code for setting up "Backward Facing Step Flow",
//
void Plate::backFacingStep(int mpiSize)
{
   unsigned h_up = (unsigned)((baitParamsPtr_->ratio_hH)*engParamsPtr_->height_);   // upstream channel
   //unsigned H_down = engParamsPtr_->height_;                         // downstream channel

   Patch *ptrPatch=NULL;
   
   if (myRank % mpiSize == 0)
   {        
      for (unsigned i=h_up+1; i<engParamsPtr_->height_; i++)
      {
         ptrPatch = &(this->smartAccess(i,0));
         ptrPatch->setIsStep(true);
         if (ptrPatch->getActiveRegion() == true)
            ptrPatch->setActive(false);
         //ptrPatch->setU_vel(0.0);
         //ptrPatch->setV_vel(0.0);    
      }
   }

}*/


//
// JM 25APR11 - Prints map of nutrients, organisms and spores for GNUPlot
//
void Plate::printMap(unsigned nbLoops, ofstream *outputFiles[], int numOutFiles, unsigned index)
{
   int height = engParamsPtr_->height_;
   unsigned width = engParamsPtr_->width_;
   Patch *ptrPatch=NULL;
   
   for (int i=height-1; i>=0; i--)
   {
      for(unsigned j=0; j<width; j++)
      {
         ptrPatch = &(this->smartAccess(i,j));
         
         if ( ptrPatch->getActiveRegion() || (i >0 && i<height-1 && j>0 && j<width-1) 
            || ptrPatch->getSubstrateType() == -1 )
         {
            (*outputFiles[index]) << (j+1) << " " << (height-i) << " " 
               << ptrPatch->getIsOrganism() << " " 
               << ptrPatch->getSporeLevel() << endl;
         }
      }
      (*outputFiles[index]) << endl;
   }
 
}


//
// JM 12Dec12 - Prints horizontal cross-section of nutrients, organisms and spores plotting
//
void Plate::printHorizCross(unsigned nbLoops, ofstream *ptrOutputFile2) //ofstream *outputFiles[], int numOutFiles, unsigned index)
{
   unsigned height = engParamsPtr_->height_;
   unsigned width = engParamsPtr_->width_;
   Patch *ptrPatch=NULL;

   unsigned i = (engParamsPtr_->height_-2)/2;
   
   for(unsigned j=0; j<width; j++)
   {
      ptrPatch = &(this->smartAccess(i,j));
      
      if ( ptrPatch->getActiveRegion() || (i >0 && i<height-1 && j>0 && j<width-1) 
         || ptrPatch->getSubstrateType() == -1 )
      {
         //(*ptrOutputFile2)  << ptrPatch->getPollenRcv() << " ";
      }
   }
   (*ptrOutputFile2) << endl; 
}

//
// 25Jan2013 - print stats (e.g. seedStock) for each individual plant colony (genet)
//
void Plate::printPlantStats(unsigned nbLoops, ofstream *ptrOutputFile2)
{
   unsigned numSeedlings = Organism::getNumSeedlings();
   cout << "Num Genets = " << numSeedlings << endl;
   int *seedlingIDs = new int[numSeedlings];
   unsigned *iPosLow = new unsigned[numSeedlings];
   unsigned *iPosHigh = new unsigned[numSeedlings];
   unsigned *jPosLow = new unsigned[numSeedlings];
   unsigned *jPosHigh = new unsigned[numSeedlings];
   unsigned *iPos = new unsigned[numSeedlings];
   unsigned *jPos = new unsigned[numSeedlings];
   unsigned *area = new unsigned[numSeedlings];
   //double *seedStock = new double[numSeedlings];
   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;
   
   // initialise the variables
   for(unsigned x=0; x<numSeedlings; x++)
   {
      seedlingIDs[x] = -1;
      iPosHigh[x] = 0;
      jPosHigh[x] = 0;
      iPosLow[x] = engParamsPtr_->height_+10;
      jPosLow[x] = engParamsPtr_->width_+10;
      area[x] = 0;
      //seedStock[x] = 0;
   }

   for( unsigned i=0; i<nbCells; i++ )
   {
      if (this->cells_[i].getIsOrganism() == true && this->cells_[i].getActiveRegion())
      {
         for(unsigned x=0; x<numSeedlings; x++)
         {
            // when find index with matching ID or else first empty index
            if ( seedlingIDs[x] == this->cells_[i].getSeedlingID() || seedlingIDs[x] == -1 )
            {
               seedlingIDs[x] = this->cells_[i].getSeedlingID();
               
               if (this->cells_[i].getIPos() > iPosHigh[x])
                  iPosHigh[x] = this->cells_[i].getIPos();
               if (this->cells_[i].getIPos() < iPosLow[x])
                  iPosLow[x] = this->cells_[i].getIPos();
               if (this->cells_[i].getJPos() > jPosHigh[x])
                  jPosHigh[x] = this->cells_[i].getJPos();
               if (this->cells_[i].getJPos() < jPosLow[x])
                  jPosLow[x] = this->cells_[i].getJPos();
                  
               area[x]++;
               break; 
            }
         }  
      }
   }
   
   // find mid-point for each colony
   for(unsigned x=0; x<numSeedlings; x++)
   {
      iPos[x] = (iPosHigh[x] + iPosLow[x]) / 2;
      jPos[x] = (jPosHigh[x] + jPosLow[x]) / 2;
   }
   
   
   // Print stats to output file
   for(unsigned x=0; x<numSeedlings; x++)
   {
      (*ptrOutputFile2) << nbLoops << " "
         << seedlingIDs[x] << " " 
         << iPos[x] << " " << jPos[x] << " " 
         << area[x] << " ";
   }
}


//
// JM 08Sep11 - Set elevation above sea level for each patch and calculate whether in S anglica niche
//
void Plate::initElevNicheSubstrate()
{
   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;

   for( unsigned i=0; i<nbCells; i++ )
   {
      if (this->cells_[i].getSubstrateType() > -1)
      {
         this->cells_[i].setElevation( xy_elevations[i] );
         // substrateTypes[i] == 0 means no substrate for Undaria to bind to
         if (substrateTypes[i] != 0)
            this->cells_[i].setSubstrateType(substrateTypes[i]);

         //this->cells_[i].setIsElevationNiche(upperLimit, lowerLimit);
      }
   }
}


//
// 28Jan2014 - true if appropriate elevation/niche type for Undaria pinnatifida
//
void Plate::initUndariaNiche()
{
   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;

   for( unsigned i=0; i<nbCells; i++ )
   {
      this->cells_[i].setIsUndariaNiche();
   }
}



//
// JM 12Dec13 - input upper/lower elevation niche limits from input file
//
/*void Plate::setNicheLimits()
{
   upperLimit = baitParamsPtr_->upperLimit;
   lowerLimit = baitParamsPtr_->lowerLimit;
}*/


//
// JM 02Nov2011 - calculate lower/upper limit for S anglica growth
//                Gray et al.'s (1989,1995) S. anglica Simple Niche Model:
//
/*void Plate::setNicheLimitsSpartina()
{
   upperLimit=0.0, lowerLimit=0.0;
   int numStdevs=1;

   // Gray's Simple Niche Model (don't use fetch/latitude/estuary area) (Gray et al., 1989, 1995)
   if (baitParamsPtr_->fetch == 0.0 || baitParamsPtr_->estuaryArea == 0.0
         || baitParamsPtr_->latitude == 0.0 )
   {
   
      while (upperLimit <= lowerLimit)
      {
         // upper elevation limit for growth (in metres OD Newlyn) 
         // Mean ~1.277 metres
         upperLimit = (-0.571 + sampleGaussian(0.132, numStdevs)) 
                     + (0.528+ sampleGaussian(0.018, numStdevs)) * baitParamsPtr_->springTidalRange;
        
         // lower elevation limit for growth (in metres OD Newlyn)
         // Mean ~0.7525
         lowerLimit = (-0.826 + sampleGaussian(0.127, numStdevs)) 
                     + (0.451+ sampleGaussian(0.018, numStdevs)) * baitParamsPtr_->springTidalRange;
      }
      
      if(myRank == 0)
      {
         cout << "\nGray's Simple Niche Model for S. anglica (Gray et al., 1989, 1995):\n"
         << "Upper Limit = " << upperLimit << " m (OD Newlyn)"
         << "\nLower Limit = " << lowerLimit << " m (OD Newlyn)\n" << endl;
      }
      
   }  // end if
   // Gray's Full Niche Model (fetch, latitude, estuary area included) (Gray et al., 1989, 1995)
   else
   {
   
      while (upperLimit <= lowerLimit)
      {
         // upper elevation limit for growth (in metres OD Newlyn) 
         // Mean ~1.277 metres
         upperLimit =  ( 4.74 + sampleGaussian(2.29, numStdevs) ) 
                     + ( 0.483 + sampleGaussian(0.028, numStdevs) ) * baitParamsPtr_->springTidalRange
                     + ( 0.068 + sampleGaussian(0.020, numStdevs) ) * baitParamsPtr_->fetch
                     - ( 0.099 + sampleGaussian(0.045, numStdevs) ) * baitParamsPtr_->latitude;
        
         // lower elevation limit for growth (in metres OD Newlyn)
         // Mean ~0.7525
         lowerLimit = (-0.805 + sampleGaussian(0.102, numStdevs) ) 
                     + ( 0.366 + sampleGaussian(0.018, numStdevs) ) * baitParamsPtr_->springTidalRange
                     + ( 0.053 + sampleGaussian(0.016, numStdevs) ) * baitParamsPtr_->fetch
                     + ( 0.135 + sampleGaussian(0.025, numStdevs) ) * log(baitParamsPtr_->estuaryArea);
      }  // end while
      cout << "\nGray's Full Niche Model for S. anglica (Gray et al., 1989, 1995):\n"
      << "Upper Limit = " << upperLimit << " m (OD Newlyn)"
      << "\nLower Limit = " << lowerLimit << " m (OD Newlyn)\n" << endl;
     
   }  // end else
   
}*/


//
// To gather and store all patch values from all plates on root node for display
//
void Plate::gatherPlatesMPI()
{
   int sendcnt = engParamsPtr_->width_ * engParamsPtr_->height_;
   int recvcount = sendcnt;  
   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;
   
   // Fill patchValues with the values for all patches in current plate
   for( unsigned i=0; i<nbCells; i++ )
   {
      if ( this->cells_[i].getActiveRegion() )
      {
         //patchValues[i] = this->cells_[i].getSeedLevel();
         patchValues[i] = this->cells_[i].getSporeLevel();
         elevationValues[i] = this->cells_[i].getElevation();
         substrateValues[i] = this->cells_[i].getSubstrateType();
         biomassValues[i] = this->cells_[i].calcSporoBiomass();
         sporeValues[i] = this->cells_[i].getSporeLevel();
      }
      else
      {
         patchValues[i] = 0.0;
         elevationValues[i] = 0.0;
         substrateValues[i] = 0;
         biomassValues[i] = 0.0;
         sporeValues[i] = 0;
      }
   }
   
   // Gather all patch values from all plates and store on Rank 0 (MPI)
   MPI_Gather(patchValues, sendcnt, MPI_FLOAT, allPatchValues, recvcount, MPI_FLOAT, 0, MPI_COMM_WORLD);
   MPI_Gather(elevationValues, sendcnt, MPI_FLOAT, allElevationValues, recvcount, MPI_FLOAT, 0, MPI_COMM_WORLD);
   MPI_Gather(substrateValues, sendcnt, MPI_INT, allSubstrateValues, recvcount, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Gather(biomassValues, sendcnt, MPI_FLOAT, allBiomassValues, recvcount, MPI_FLOAT, 0, MPI_COMM_WORLD);
   MPI_Gather(sporeValues, sendcnt, MPI_INT, allSporeValues, recvcount, MPI_INT, 0, MPI_COMM_WORLD);
}


//
// Count the number of patches that are occupied by an agent
//
int Plate::calcTotalCoverage()
{
   int totalCoverage = 0;
   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;
   float macroSize = baitParamsPtr_->macroSize;
      
   for( unsigned i=0; i<nbCells; i++ )
   {
      if ( this->cells_[i].getActiveRegion()
               && this->cells_[i].getIsOrganism() == true
               && this->cells_[i].countNumSporophytes(macroSize, 100.0) > 0 )
         totalCoverage++;
   }
   return totalCoverage;
}


//
// Count the number of patches on boundaries of colonies
//
int Plate::calcTotalPerimeter()
{
   int totalPerim = 0;
   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;
      
   for( unsigned i=0; i<nbCells; i++ )
   {
      if (this->cells_[i].getIsOrganism() == true && this->cells_[i].getActiveRegion() )
      {
         // Check all neighbouring patches to see if there are any empty ones (no agent)
         if ( this->cells_[i].getAnyFreeNeighbour() )
            totalPerim++;
      }
   }
   return totalPerim;
}



//
// JM - 20Feb2007 Method to sample from Gaussian (Normal) distribution using GSL
//        (GNU Scientific Library).
//
double Plate::sampleGaussian(double stdDev, int numDevs)
{
   double sample=0.0;
   sample = gsl_ran_gaussian(Bait::rng, stdDev);

   // Correct for samples outside biologically realistic range.
   if ( sample > stdDev*numDevs )
      sample = stdDev*numDevs;
   else if ( sample < -(stdDev*numDevs) )
      sample = -(stdDev*numDevs);

   return sample;
}


//
// 02Oct2012 - Find the MPI rank for the plate based on global i,j co-ords
//
int Plate::findPlate(int iGlobal, int jGlobal)
{
   int mpiDim = (int)sqrt(numPlatesMPI);
   int buffer = 2;
   int netHeight = engParamsPtr_->height_ - buffer;
   int netWidth = engParamsPtr_->width_ - buffer;
     
   int plateCol = (int) (jGlobal-1)/netWidth;
   // Correct for unequal buffer zone size in edge plates (where buffer = 1)
   if (plateCol < 0)
      plateCol = 0;
   if (plateCol >= mpiDim)
      plateCol = mpiDim -1;

   int plateRow = (int) (iGlobal-1)/netHeight;
   // Correct for unequal buffer zone size in edge plates (where buffer = 1)
   if (plateRow < 0)
      plateRow = 0;
   if (plateRow >= mpiDim)
      plateRow = mpiDim -1;

   int plateRank = (plateRow * mpiDim) + plateCol;
   
   return plateRank;
}

//
// 03Oct12 - Convert global i co-ord (all plates in MPI, not counting buffers) to local i co-ord
//
int Plate::convGlobalToLocali(int iGlob)
{
   int mpiDim = (int)sqrt(numPlatesMPI);
   int buffer = 2;
   int netHeight = engParamsPtr_->height_ - buffer;
   
   // Find what row of MPI plates it is in
   int plateRow = (int) (iGlob-1)/netHeight;
   // Correct for unequal buffer zone size in edge plates (where buffer = 1)
   if (plateRow < 0)
      plateRow = 0;
   if (plateRow >= mpiDim)
      plateRow = mpiDim -1;
   
   int iLocal = iGlob - (netHeight*plateRow);
   return iLocal;
}

//
// 03Oct12 - Convert global j co-ord (all plates in MPI, not counting buffers) to local j co-ord
//
int Plate::convGlobalToLocalj(int jGlob)
{
   int mpiDim = (int)sqrt(numPlatesMPI);
   int buffer = 2;
   int netWidth = engParamsPtr_->width_ - buffer;
   
   // Find what column of MPI plates it is in
   int plateCol = (int) (jGlob-1)/netWidth;
   // Correct for unequal buffer zone size in edge plates (where buffer = 1)
   if (plateCol < 0)
      plateCol = 0;
   if (plateCol >= mpiDim)
      plateCol = mpiDim -1;
   
   int jLocal = jGlob - (netWidth*plateCol);
   return jLocal;
}

//
// 03Oct12 - Convert local i co-ord to global i co-ord (all plates in MPI, not counting buffers)
//
int Plate::convLocalToGlobali(int iLoc, int rank)
{
   int height = engParamsPtr_->height_;
   int mpiDim = sqrt(numPlatesMPI);
   int plateRow = rank / mpiDim;
   int buffer = 2;
   
   // Check that jLoc not in buffer zones (don't count in global co-ords)
   if (iLoc == 0 && plateRow > 0)
      return -1;
   else if (iLoc == height-1 && plateRow < mpiDim-1)
      return -1;

   int iGlobal = iLoc + (height * plateRow);
   // don't count buffer zones in global co-ord
   if (plateRow > 0)
      iGlobal -= (buffer * plateRow);
      
   return iGlobal;
}


//
// 03Oct12 - Convert local j co-ord to global j co-ord (all plates in MPI, not counting buffers)
//
int Plate::convLocalToGlobalj(int jLoc, int rank)
{
   int width = engParamsPtr_->width_;
   int mpiDim = sqrt(numPlatesMPI);
   int plateCol = rank % mpiDim;
   int buffer = 2;
   
   // Check that jLoc not in buffer zones (don't count in global co-ords)
   if (jLoc == 0 && plateCol > 0)
      return -1;
   else if (jLoc == width-1 && plateCol < mpiDim-1)
      return -1;

   int jGlobal = jLoc + (width * plateCol);
   // don't count buffer zones in global co-ord
   if (plateCol > 0)
      jGlobal -= (buffer * plateCol);
      
   return jGlobal;
}



// Calculate the row num (i) of receiving patch based on distance/angle from source
int Plate::calcIPosDistAng(unsigned sourceIPos, float dist, float angle)
{
   int iPos = 0, y_int = 0;
   double y=0.0;

   y = sin(angle) * dist;
   //rounding
   y_int = round(y);
   
   // Apply y-offset to calculate iPos receiving pollen  
   iPos = sourceIPos - y_int;
   
   // Check that pollen lands within bounds of plate
   if (iPos >=0 && iPos < (int)engParamsPtr_->height_)
      return iPos;
   else
      return -1;
}

// Calculate the column num (j) of receiving patch based on distance/angle from source
int Plate::calcJPosDistAng(unsigned sourceJPos, float dist, float angle)
{
   int jPos = 0, x_int = 0;
   double x=0.0;

   x = cos(angle) * dist;
   //rounding
   x_int = round(x);
   
   // Apply y-offset to calculate iPos receiving pollen  
   jPos = sourceJPos - x_int;
   
   // Check that pollen lands within bounds of plate
   if (jPos >=0 && jPos < (int)engParamsPtr_->height_)
      return jPos;
   else
      return -1;
}


// ------------------------------ End Of File ------------------------------- //
