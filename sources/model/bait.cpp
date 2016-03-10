// -------------------------------------------------------------------------- //
// bait.cpp                                                                   //
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


#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR

// -------------------------------------------------------------------------- //
// Libraries                                                                  //
// -------------------------------------------------------------------------- //

#include "mpi.h"
#include "bait.hpp"
#include "parameters.hpp"
#include "patch.hpp"
#include "../engine/myException.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
using namespace std ;


// -------------------------------------------------------------------------- //
// Static member creation                                                     //
// -------------------------------------------------------------------------- //

const BaitParameters* Bait::baitParamsPtr_ = 0 ;
ofstream* Bait::ptrOutputFile=0;
ofstream* Bait::ptrOutputFile2=0;
gsl_rng* Bait::rng=0;
// days per month
const unsigned Bait::daysPerMonth[12] = {31,28,31,30,31,30,31,31,30,31,30,31};


// -------------------------------------------------------------------------- //
// Constructor                                                                //
// -------------------------------------------------------------------------- //

Bait::Bait()
: System()
{
}


// -------------------------------------------------------------------------- //
// Destructor                                                                 //
// -------------------------------------------------------------------------- //

Bait::~Bait()
{
   // Free memory associated with random number generator.
   gsl_rng_free (rng);
}


// -------------------------------------------------------------------------- //
// Initialiaze the static parameters pointer                                  //
// -------------------------------------------------------------------------- //

void Bait::initBaitParamsPtr( const BaitParameters *inPtr )
{
   baitParamsPtr_ = inPtr ;
}


// -------------------------------------------------------------------------- //
// All components initialisation                                              //
// -------------------------------------------------------------------------- //

void Bait::init(Parameters &inParams, ofstream *outputFile, ofstream *outputFile2, ofstream *outputFilesLB[], int numOutFilesLB_, ofstream *outputFilesLB2[], int numOutFilesLB2_, ofstream *outputFilesMG[], int numOutFilesMG, int *xy_elevs, int *agentLocations, int *substrateLocations)
{
   mpiSize=0;

   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
   
   // Static pointer intialisation
   initBaitParamsPtr( & inParams.baitParams_ ) ;

   // System initialisation
   System::initEngParamsPtr ( & inParams.engParams_  ) ;
   
   highestNumBac = 0;
   lastLoopNumber_ = 0;
   gettimeofday( &lastTime_, NULL );
   sporesPresent = false;
   plateSporeLevel = 0.0;

   // JM - initialize timer variables
   displayDuration=0;
   worldEvolveDuration=0;
   staticActivityDuration=0;
   startMoveDuration=0;
   moveDuration=0;
   solveConflictsDuration=0;
   interactionsDuration=0;
   reproductionDuration=0;

   if (outputFile != NULL)
      ptrOutputFile = outputFile;      // JM
   if (outputFile2 != NULL)
      ptrOutputFile2 = outputFile2;
    
   outFilesLB = outputFilesLB;
   numOutFilesLB = numOutFilesLB_;
   outFilesLB2 = outputFilesLB2;
   numOutFilesLB2 = numOutFilesLB2_; 
   outFilesMG = outputFilesMG;
   numOutFilesMG = numOutFilesMG;
          
   // Create instance of Mersenne Twister random number generator (RNG)
   // from GNU Scientific Library (GSL) (Matsumoto & Nishimura, 2002).
   gsl_rng_env_setup();
   rng = gsl_rng_alloc (gsl_rng_mt19937);
   if (rank == 0)
   {
      gsl_rng_set(rng, time(NULL));               // "Seed" the RNG.
      if (mpiSize > 1) {
         for (int i=1; i<mpiSize; i++)
         {
            int rngSeed = gsl_rng_get(rng);
            MPI_Send(&rngSeed, 1, MPI_INT, i, 17, MPI_COMM_WORLD);
         }
      }
   }
   else
   {
      MPI_Status status;
      int rngSeed;
      MPI_Recv(&rngSeed, 1, MPI_INT, 0, 17, MPI_COMM_WORLD, &status);
      gsl_rng_set(rng, rngSeed);
   }

   // Initiate water temperature, solar radiation & day length parameters
   initEnvParams();

   // Plate initialisation
   inParams.baitParams_.maxTotalNutrientLevel_ = 0 ;
   World<Patch>::initEngParamsPtr ( & inParams.engParams_  ) ;
   Plate::initBaitParamsPtr       ( & inParams.baitParams_ ) ;
   XPlate::initXParamsPtr         ( & inParams.xParams_    ) ;

   plate_.init(rank, mpiSize, xy_elevs, agentLocations, substrateLocations, waterTemp, solarRadiation, dayLengths) ;

   // Fabrics initialisation
   bacteriasFabric_.init( baitParamsPtr_->maxOrganisms_ ) ;

   // Organism Agents initilisation
   Organism::initStaticPtr( & inParams.baitParams_, & inParams.engParams_, &bacteriasFabric_ ) ;
   Organism::initWorldPtr(&plate_);

   // First bacterias birth
   if (baitParamsPtr_->initialOrganisms_ > (unsigned)mpiSize)
      currentNumBac = baitParamsPtr_->initialOrganisms_/mpiSize;
   else
      currentNumBac = baitParamsPtr_->initialOrganisms_;
   
   if (rank == 0)
   {
      if (baitParamsPtr_->initialOrganisms_ % mpiSize > 0)
      {
         cout << "Bait::init() - Warning: initialOrganisms not divisible by " << mpiSize << endl;
         cout << currentNumBac << " agents initialised in each node" << endl;
      }
   }
   currentNumAdult = 0;
   numNewRecruits = 0;
   
   
   if (baitParamsPtr_->initialOrganisms_ < 1)
   {
      cout << "Bait::init() - Error: initialOrganisms == 0 " << endl;
   }
   else
   {
      switch (baitParamsPtr_->initDistribute) 
      {
         case 0:
            randDistributeAgents();
            break;
         case 1:
            evenlyDistributeAgents(); 
            break;
         default:
            cout << "Bait::init(): Value of initDistrib out of bounds";
      }
   }
   MPI_Barrier(MPI_COMM_WORLD);
}


//
// Code to evenly distribute agents across plate at beginning of simulation.
//
void Bait::evenlyDistributeAgents()
{
   Organism *ptrOrganism;
   int iPos, jPos;
   int agentNum;

   for( unsigned i=0; i<currentNumBac; i++ )
   {
      // Create and initialise the agent
      ptrOrganism = bacteriasFabric_.newAgent() ;
      ptrOrganism->init(i, rank);
      
      if (baitParamsPtr_->initialOrganisms_%16 == 0)
      {
         agentNum = i % 16;

         iPos = gsl_rng_uniform_int(rng, (getEngParamsPtr()->height_-2)/4);
         jPos = gsl_rng_uniform_int(rng, (getEngParamsPtr()->width_-2)/4 );
         
         if (agentNum < 4)
            jPos = jPos + ((getEngParamsPtr()->width_-2)/4)*agentNum;
         else if (agentNum < 8)
         {
            jPos = jPos + ((getEngParamsPtr()->width_-2)/4)*(agentNum-4);
            iPos = iPos + ((getEngParamsPtr()->height_-2)/4)*1;
         }
         else if (agentNum < 12)
         {
            jPos = jPos + ((getEngParamsPtr()->width_-2)/4)*(agentNum-8);
            iPos = iPos + ((getEngParamsPtr()->height_-2)/4)*2;
         }
         else if (agentNum < 16)
         {
            jPos = jPos + ((getEngParamsPtr()->width_-2)/4)*(agentNum-12);
            iPos = iPos + ((getEngParamsPtr()->height_-2)/4)*3;
         }
         else
            cout << "More than 16 bacteria specified." << endl;
      }
      else if (baitParamsPtr_->initialOrganisms_ == 4)
      { 
         iPos = gsl_rng_uniform_int( rng, (getEngParamsPtr()->height_-2)/2 );
         jPos = gsl_rng_uniform_int( rng, (getEngParamsPtr()->width_-2)/2 );

         if (i < 2)
            jPos = jPos + ( (getEngParamsPtr()->width_-2)/2 )*i;
         else if (i < 4)
         {
            iPos = iPos + ( (getEngParamsPtr()->height_-2)/2 );
            jPos = jPos + ( (getEngParamsPtr()->width_-2)/2 )*(i-2);
         }            
      }
      else if (baitParamsPtr_->initialOrganisms_ == 1)
      {
         iPos = getEngParamsPtr()->height_/2 -1;
         jPos = getEngParamsPtr()->width_/2 -1;
      }
      else
      {
         // Insert organisms randomly
         iPos = gsl_rng_uniform_int(rng, (getEngParamsPtr()->height_-3) )+1;
         jPos = gsl_rng_uniform_int(rng, (getEngParamsPtr()->width_-3) )+1;

      }
      
      // Make sure agents not initialized in buffer zones (edges)
      if (iPos < 1)
         iPos = 1;
      if (jPos < 1)
         jPos = 1;
      
      // Insert the organisms to the plate
      plate_( iPos, jPos ).addOrganism( ptrOrganism ) ;
   
   }  // end for

}


//
// Code to randomly distribute agents across plate at beginning of simulation.
//
void Bait::randDistributeAgents()
{
   Organism *ptrOrganism;
   int iPos, jPos;

   for( unsigned i=0; i<currentNumBac; i++ )
   {
      // Create and initialise the new organism
      ptrOrganism = bacteriasFabric_.newAgent() ;
      ptrOrganism->init(i, rank);
      
      // Insert organisms randomly
      iPos = gsl_rng_uniform_int(rng, (getEngParamsPtr()->height_-3) )+1;
      jPos = gsl_rng_uniform_int(rng, (getEngParamsPtr()->width_-3) )+1;
      
      // only initialise agents in patches within Undaria niche
      int numAttempts =0;
      while ( !(plate_( iPos, jPos ).getIsUndariaNiche()) && numAttempts < 10000)
      {
         iPos = gsl_rng_uniform_int(rng, (getEngParamsPtr()->height_-3) )+1;
         jPos = gsl_rng_uniform_int(rng, (getEngParamsPtr()->width_-3) )+1;
         numAttempts++;
      }

      if (numAttempts >= 5000)
         cout << rank << ": Bait::randDistributeAgents(): Failed to find a niche site to initialise agent " 
               << i << "." << endl;
                          
      // Insert the organisms to the plate
      plate_( iPos, jPos ).addOrganism( ptrOrganism ) ;
      
   }  // end for
   ptrOrganism = NULL;
}





// -------------------------------------------------------------------------- //
// Run plate in system                                                        //
// -------------------------------------------------------------------------- //

void Bait::run()
{
   System::run( plate_ ) ;
}


// -------------------------------------------------------------------------- //
// 1. Display plate                                                              //
// -------------------------------------------------------------------------- //

void Bait::displayTime( unsigned nbLoops)
{
   //cout << "1. Display Time" << endl;
   start = clock();         // JM
   
   if (mpiSize > 1)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      plate_.gatherPlatesMPI();
   }
   
   // Update the bacterial biomass levels for all patches.
   if (rank == 0)
      plate_.display( nbLoops, bacteriasFabric_, mpiSize, rank) ;
   
   stop = clock();         // JM
   displayDuration += ((double)(stop-start)/CLOCKS_PER_SEC);
}


//
// 2. JM 20APR11 - Apply lattice boltzmann algorithm to update velocities at each patch
//
void Bait::latticeBoltzmannTime(unsigned nbLoops)
{
   /*MPI_Barrier(MPI_COMM_WORLD);
   start = clock();
   int outFileIndex=0;
    
   //cout << "\n" << rank << ": a.Latt" << endl;
   
   if ((int)nbLoops == baitParamsPtr_->loopAddObst)
      plate_.addObstacleLBE(mpiSize);      
   else if ((int)nbLoops == baitParamsPtr_->loopRemoveObst)
      plate_.removeObstacleLBE();

   if (baitParamsPtr_->initUvelocity > 0.0 || baitParamsPtr_->initVvelocity > 0.0)
   {
      // Optional: Code for initialising "backward facing step flow"
      if (baitParamsPtr_->ratio_hH < 1.0 && nbLoops == 1)
         plate_.backFacingStep(mpiSize);
   
      // 1. LB streaming function
      plate_.moveLBE(mpiSize, nbLoops);
      
      // 2. call hydrovar, equili, collis and force subroutines
      plate_.methodsLBE(nbLoops);

      // 3. LB obstacle function
      if (baitParamsPtr_->loopAddObst != -1)
         plate_.obstacleLBE(mpiSize, nbLoops);
   
      // 4. Print out LBE details for diagnostics   
      if (nbLoops % ((unsigned)getEngParamsPtr()->nbLoops_/10) == 0)
      {
         if (mpiSize < 2)  // configLBE() & printMap() not yet compatible with MPI (To Do)
         {
            outFileIndex = (nbLoops / (getEngParamsPtr()->nbLoops_/10))-1;    
            plate_.configLBE(nbLoops, outFilesLB, numOutFilesLB, outFileIndex);
            plate_.printMap(nbLoops, outFilesMG, numOutFilesMG, outFileIndex); 
         }        
      }
      if (nbLoops == 1 || nbLoops % ((unsigned)getEngParamsPtr()->nbLoops_/10) == 0)
         plate_.profileLBE(nbLoops, outFilesLB2, numOutFilesLB2);
      if(nbLoops % 1000 == 0)   
         plate_.diag0D(nbLoops);//, ptrOutputFile2);
   }
   else if ( getEngParamsPtr()->nbLoops_ > 5000 && nbLoops % ((unsigned)getEngParamsPtr()->nbLoops_/10) == 0)
   {
      outFileIndex = (nbLoops / (getEngParamsPtr()->nbLoops_/10))-1;
      //if (mpiSize < 2)  // printMap() not yet compatible with MPI (To Do)
         //plate_.printMap(nbLoops, outFilesMG, numOutFilesMG, outFileIndex);         
   }
   

   stop = clock();
   boltzmannDuration += ((double)(stop-start)/CLOCKS_PER_SEC);*/
}



// -------------------------------------------------------------------------- //
// 3. Invoke world evolution                                                     //
// -------------------------------------------------------------------------- //

void Bait::worldEvolutionTime(unsigned nbLoops)
{
   MPI_Barrier(MPI_COMM_WORLD);
   //cout << "2. World Evolution Time" << endl;
   start = clock();         // JM
   
   sporesPresent = plate_.checkSporesPresent(mpiSize);
   
   if (mpiSize > 1)
   {
      // Check to see if spores present in any other rank
      int sndNbSpore = (int)sporesPresent, rcvNbSpore = 0;
      MPI_Reduce(&sndNbSpore, &rcvNbSpore, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Bcast(&rcvNbSpore, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if (rcvNbSpore > 0)
         sporesPresent = true;
   }
   // Degradation of spores according to half-life (non-CUDA, non-MPI version)
   else if (sporesPresent && baitParamsPtr_->CUDALoops == 0)
      plate_.allPatchesDevelop(nbLoops);
      
   // Update the current water temperature (same for all patches)
   plate_.updateWaterSolar(nbLoops); 

   // JM - Diffusion of nutrients by Fick's Law of diffusion
   if (sporesPresent)
   {
      // Non-CUDA algorithm
      if (baitParamsPtr_->CUDALoops == 0)
      {
         if (baitParamsPtr_->diffusionRate > 1.0E-06)
         {      
            if (mpiSize > 1)
            {
               if (nbLoops == 1 && rank == 0)
                  cout << "Implementing parallel diffusion model." << endl;
               plate_.parallelDiffusion(mpiSize, nbLoops);
            }
            else if (plateSporeLevel > 0.0)
               plate_.diffusion(nbLoops);
         }
         else if (mpiSize > 1) // just call MPI routine for communicating between plates
            plate_.noDiffusionMPI(mpiSize, nbLoops);
      } 
      // Else use CUDA 
      else if (nbLoops % baitParamsPtr_->CUDALoops == 0)
      {
         if (mpiSize == 1)
            plate_.diffusionCUDA(nbLoops);
         else
            plate_.parallelDiffusion_CUDA(mpiSize, nbLoops);
      }
   }
   
   // Apply Lattice Boltzmann Equation (LBE) for fluid dynamics
   //if (baitParamsPtr_->initUvelocity > 0.0 || baitParamsPtr_->initVvelocity > 0.0)
      //plate_.allCellsLBE(nbLoops);

   stop = clock();         // JM
   worldEvolveDuration += ((double)(stop-start)/CLOCKS_PER_SEC);
}


// -------------------------------------------------------------------------- //
// 4. All agents do their static activity (survive cost)                      //
// -------------------------------------------------------------------------- //
void Bait::staticActivityTime(unsigned nbLoops)
{
   MPI_Barrier(MPI_COMM_WORLD);
   //cout << "3. Static Time" << endl;
   start = clock(); 
   unsigned loopOutResults = getEngParamsPtr()->graphPeriod;  // Write results to output file
   unsigned sampleLoop = getEngParamsPtr()->sampleLoop;    // Loop to print results
   float lps=0.0;
   
   // Release spores into environment.
   if (currentNumBac && baitParamsPtr_->sporeProd > 0.0)
      currentNumBac = bacteriasFabric_.aliveAgentsAction(&Organism::produceSpore);
     
   // Calculate loops per second
   if (((nbLoops-sampleLoop)%loopOutResults)==0)
      lps = calcLoopsPerSec(nbLoops);

   // Display the number of agents/other stats to console, if needed // JM
   if( baitParamsPtr_->displayOrganismsNumber_ && ( ((nbLoops-sampleLoop)%loopOutResults)==0 ) )
      outputResults(nbLoops, currentNumBac, lps);     
   
   // Keep track of when population increasing/decreasing (in case need to defragment)
   if (currentNumBac > highestNumBac)
      highestNumBac = currentNumBac;
         
   if (currentNumBac > 1000)   
   {
      // Defragment bacterial fabric if population size decreases by 5%.
      // This has very large performance benefit when population decreasing.
      if (currentNumBac < highestNumBac*0.8)
      {  	
         bacteriasFabric_.defragmentFabric();
         highestNumBac = currentNumBac;
      }
   }

   stop = clock(); 
   staticActivityDuration += ((double)(stop-start)/CLOCKS_PER_SEC);

}


// -------------------------------------------------------------------------- //
// 5. All alive agents can reproduce                                          //
// -------------------------------------------------------------------------- //

void Bait::reproductionTime(unsigned nbLoops)
{
   MPI_Barrier(MPI_COMM_WORLD);
   //cout << "4. Reprod Time" << endl;
   start = clock();
   
   // If spores present, may germinate into new agents
   if (sporesPresent)
      plate_.germinateSpores(nbLoops);

   stop = clock(); 
   reproductionDuration += ((double)(stop-start)/CLOCKS_PER_SEC);
}


// -------------------------------------------------------------------------- //
// 6. All agents movement/shunting                                            //
// -------------------------------------------------------------------------- //


void Bait::startMoveTime()
{
   start = clock();
   
   plate_.startMoveTime();
   
   stop = clock();
   startMoveDuration += ((double)(stop-start)/CLOCKS_PER_SEC);
}


void Bait::moveTime()
{
   MPI_Barrier(MPI_COMM_WORLD);
   //cout << "5. Move Time" << endl;
   start = clock(); 
      
   currentNumBac = bacteriasFabric_.aliveAgentsAction( &Organism::move ) ;
   
   if (mpiSize > 1)
      currentNumBac += plate_.parallelMovement(mpiSize);
        
   stop = clock();
   moveDuration += ((double)(stop-start)/CLOCKS_PER_SEC);
}



// -------------------------------------------------------------------------- //
// 7. Invoke plate conflict resolution method                                 //
// -------------------------------------------------------------------------- //

void Bait::solveConflictsTime(unsigned nbLoops)
{
   MPI_Barrier(MPI_COMM_WORLD);
   //cout << "6. Conflicts Time" << endl;
   start = clock();
   
   // Currently do nothing
   
   stop = clock();  
   solveConflictsDuration += ((double)(stop-start)/CLOCKS_PER_SEC);
}


// -------------------------------------------------------------------------- //
// 8. All alive agents interact with world (grow)                             //
// -------------------------------------------------------------------------- //
void Bait::interactionsTime(unsigned nbLoops)
{
   MPI_Barrier(MPI_COMM_WORLD);
   //cout << "7. Interactions Time" << endl;
   start = clock(); 
   
   float waterTemp = plate_.getWaterTemp();
   float dayLength = plate_.getCurrDayLength();
   float solarRad = plate_.getSolarRad();
   // Update gametophyte mean growth rate based on temp, day length, solar irradiance
   Organism::updGametGrowth(waterTemp, dayLength, solarRad);
   // Update sporophyte mean growth rate based on temp, day length, relative solar radiation
   Organism::updSporoGrowthMod(waterTemp, dayLength, solarRad);
   // Prob. of gametophyte maturity based on temp, day length
   Organism::updProbGamMature(waterTemp, dayLength);
   
   // All organisms grow
   if (currentNumBac > 0)
      currentNumBac = bacteriasFabric_.aliveAgentsAction( &Organism::growth, nbLoops ) ;
   
   // Update the total biomass levels in each patch
   plate_.calcTotalBiomass();
       
   stop = clock();  
   interactionsDuration += ((double)(stop-start)/CLOCKS_PER_SEC);
}


int Bait::calcPollenMonth(unsigned nbLoops)
{
   unsigned calDay = nbLoops % 365;
   int pollenMonth = -1;
   
   if (calDay > 90 && calDay <= 120)      // April
      pollenMonth = 0;
   else if (calDay > 120 && calDay <=151) // May
      pollenMonth = 1;
   else if (calDay > 151 && calDay <=181) // June
      pollenMonth = 2;
   else if (calDay > 181 && calDay <=212) // July
      pollenMonth = 3;
   else if (calDay > 212 && calDay <=243) // August
      pollenMonth = 4;
   else if (calDay > 243 && calDay <=273) // September
      pollenMonth = 5;
      
   return pollenMonth;
}


// JM - Method to sample from Gaussian (Normal) distribution using GSL
//        (GNU Scientific Library).
//
double Bait::sampleGaussian(double stdDev)
{
   double sample=0.0;
   sample = gsl_ran_gaussian(Bait::rng, stdDev);

   while (fabs(sample) > stdDev*3)
      sample = gsl_ran_gaussian(Bait::rng, stdDev);

   return sample;
}


//
// JM - method for printing out time taken for each of the main program functions.
//
void Bait::printDuration(double total)
{
   cout << "\nBreakdown of time for each step:\n";
   cout << "1. Duration of Display = " << displayDuration << setprecision(2)
         << " ("<< (((double)displayDuration/total)*100) <<"%)\n";
   cout << "2. Duration of World evolution = " << worldEvolveDuration << setprecision(2)
         << " ("<< (((double)worldEvolveDuration/total)*100) <<"%)\n";
   cout << "3. Duration of L Bolzmann = " << boltzmannDuration << setprecision(2)
         << " ("<< (((double)boltzmannDuration/total)*100) <<"%)\n";
   cout << "4. Duration of Static activity = " << staticActivityDuration << setprecision(2)
         << " ("<< (((double)staticActivityDuration/total)*100) <<"%)\n";
   cout << "5a. Duration of Start Move = " << startMoveDuration << setprecision(2)
         << " ("<< (((double)startMoveDuration/total)*100) <<"%)\n";
   cout << "5b. Duration of Movement = " << moveDuration << setprecision(2)
         << " ("<< (((double)moveDuration/total)*100) <<"%)\n";
   cout << "6. Duration of Solve conflicts = " << solveConflictsDuration << setprecision(2)
         << " ("<< (((double)solveConflictsDuration/total)*100) <<"%)\n";
   cout << "7. Duration of Interactions = " << interactionsDuration << setprecision(2)
         << " ("<< (((double)interactionsDuration/total)*100) <<"%)\n";
   cout << "8. Duration of Reproduction = " << reproductionDuration << setprecision(2)
         << " ("<< (((double)reproductionDuration/total)*100) <<"%)\n" << endl;
}



//
// Collate results from all nodes for output to file by node 0.
//
int Bait::mpiOutputResults(unsigned nbLoops, double totBio, int totCover, int totPerim, double totSpores, float lps, double avgDensity, double umoy, double vmoy)
{
   unsigned sampleLoop = getEngParamsPtr()->sampleLoop;    // Loop to print results

   int      sndNbAdult = (int)currentNumAdult, rcvNbAdult = 0,
            sndNbRecruit = (int)numNewRecruits, rcvNbRecruit = 0,
            sndNbBac = (int)currentNumBac, rcvNbBac = 0,
            sndNbSporo = (int)numSporophytes, rcvNbSporo = 0,
            sndNbNewSporo = (int)numNewSporophytes, rcvNbNewSporo = 0,
            sndNbGameto = (int)numGametophytes, rcvNbGameto = 0,
            sndSumAgeMature = sumAgeMaturity, rcvSumAgeMature = 0,
            sndCountMature = countMature, rcvCountMature = 0,
            sndSumAgeDeath = sumAgeDeath, rcvSumAgeDeath = 0,
            sndCountDead = countDead, rcvCountDead = 0,
            sndTotCover = totCover, rcvTotCover = 0;

   double   sndTotSpores = totSpores, rcvTotSpores = 0,
            sndAvgDensity = avgDensity/mpiSize, rcvAvgDensity = 0,
            sndUmoy = umoy/mpiSize, rcvUmoy = 0,
            sndVmoy = vmoy/mpiSize, rcvVmoy = 0;

   float    sndLps = lps, rcvLps = 0;


   // Collate results from all nodes
   MPI_Reduce(&sndNbAdult, &rcvNbAdult, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&sndNbRecruit, &rcvNbRecruit, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&sndNbBac, &rcvNbBac, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&sndNbSporo, &rcvNbSporo, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&sndNbNewSporo, &rcvNbNewSporo, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&sndNbGameto, &rcvNbGameto, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&sndTotSpores, &rcvTotSpores, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&sndLps, &rcvLps, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&sndSumAgeMature, &rcvSumAgeMature, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&sndCountMature, &rcvCountMature, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&sndSumAgeDeath, &rcvSumAgeDeath, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&sndCountDead, &rcvCountDead, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&sndTotCover, &rcvTotCover, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&sndAvgDensity, &rcvAvgDensity, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&sndUmoy, &rcvUmoy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&sndVmoy, &rcvVmoy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

   // Write collated results to output file
   if (rank == 0)
   {
      // Calculate the mean age at maturity/death
      float meanAgeMature = 0, meanAgeDeath = 0;
      if ( rcvCountMature > 0)
         meanAgeMature = (float)rcvSumAgeMature/rcvCountMature;
      if ( rcvCountDead > 0)
         meanAgeDeath = (float)rcvSumAgeDeath/rcvCountDead;
     
      if (nbLoops == sampleLoop) {
               (*ptrOutputFile)  << "\nNumLoops" << "\t"
               << "NumSporo" << "\t"
               << "NewSporo" << "\t"
               << "Recruit" << "\t"
               << "NumGameto" << "\t"
               << "Spores" << "\t"
               << "AgeMature" << "\t"
               << "AgeDeath" << "\t"
               << "Temp" << "\t"
               << "Solar" << "\t"
               << "DayL" << "\t"
               << "GamGrow" << "\t"
               << "SporGrowMod" << "\t"
               << "GamMatur" << "\t";
               (*ptrOutputFile) << "lps" << endl;
      }


      (*ptrOutputFile)  << nbLoops << "\t" 
               << rcvNbSporo << "\t"
               << rcvNbNewSporo << "\t"
               << rcvNbRecruit << "\t"
               << rcvNbGameto << "\t"
               << rcvTotSpores << "\t"
               << meanAgeMature << "\t"
               << meanAgeDeath << "\t"
               << plate_.getWaterTemp() << "\t"
               << plate_.getSolarRad() << "\t"
               << plate_.getCurrDayLength() << "\t"
               << Organism::getAvgGametGrowthMod() << "\t"
               << Organism::getAvgSporoGrowthMod() << "\t"
               << Organism::getProbGamMature() << "\t";
               (*ptrOutputFile) << (rcvLps/mpiSize) << endl;

      // Display in terminal
      if ( (nbLoops-sampleLoop) % getEngParamsPtr()->graphPeriod == 0 )
      {
         cout << "Diagnostics 0D (Lp, Density, Umoy, Vmoy): " << rcvAvgDensity << " " << rcvUmoy << " " << rcvVmoy << " " << endl; 
         cout << "Lp" << nbLoops << ": " 
            << " Ad=" << rcvNbSporo
            << " Rec=" << rcvNbRecruit
            << " Gam=" << rcvNbGameto
            << " Spore=" << rcvTotSpores
            << " Temp=" << plate_.getWaterTemp()
            << " Sol=" << plate_.getSolarRad()
            << " DL=" << plate_.getCurrDayLength()
            << " (lps= " << (rcvLps/mpiSize) << ")"
            << endl;
      }
      
      //
      // Do rough comparison with Brest harbour seasonal growth curve data (05/06) for sporophytes (Voisin, 2007) 
      // - see "UndariaGEN0.6 vs Brest_08Jul15.xlsx" 
      //   
      double   brestRecruit[12] = {0.0936842105, 0.0914979757, 0.2043724696, 0.2071429705, 0.2137334286, 0.192923712, 0.3448932285, 0.8011764706, 0.8153846154, 0.3550226244, 	0.1657918552, 0.022},
               brestSporo[12] ={0.06439046, 0.068187135, 0.158956365, 0.26070475, 0.203137, 0.27578348, 0.42689459, 0.80675214, 0.94222222, 0.517033685, 0.33395676, 0.0500654},
               brestSporo04_06[12] = {0.03961704,0.1021846,0.1827848,0.23469678,0.18934211,0.24480596,0.35464077,0.9323209,1.0,0.79544117,0.389625,0.05872643},
               brestRecruit04_06[12]={0.04796163,0.1237079,0.2212850,0.17001571,0.1591003,0.20573886,0.34697759,1.0,0.98401279,0.72818267,0.2289647,0.0393568};
      double brestRecruit_log[12],
             brestSporo_log[12],
             brestSporo04_06_log[12],
             brestRecruit04_06_log[12];
      
      // Get log values       
      for (int i=0; i<12; i++)
      {
         brestRecruit_log[i] = log(brestRecruit[i]);
         brestSporo_log[i] = log(brestSporo[i]);
         brestSporo04_06_log[i] = log(brestSporo04_06[i]);
         brestRecruit04_06_log[i] = log(brestRecruit04_06[i]);
      }


      // Initialise the counters
      if (nbLoops < 1000)
      {
         index=0;
         yrIndex=0;
         for (int i=0; i<4; i++)
         { 
            sporoMax[i]=0; 
            recruitMax[i]=0;
            maxSporoIndex[i]=0; 
            maxRecruitIndex[i]=0;
         }
         startLoopTest = 14235;
      } 
      
      // Start of each season August-July
      if (nbLoops == 22995 || nbLoops == 31755 || nbLoops == 40515)
      {
         startLoopTest = nbLoops;
         yrIndex++;
      }
            
      unsigned endLoop = startLoopTest+(730*11);
      
      // Recrod the monthly abundance/recruitment data        
      if ( nbLoops >= startLoopTest && nbLoops <= endLoop) 
      {
         if ( (nbLoops-sampleLoop) % getEngParamsPtr()->graphPeriod == 0 )
         {
            index = (int)((nbLoops - startLoopTest) / 730);
            numSporo[index][yrIndex] = rcvNbSporo;
            numRecruit[index][yrIndex] = rcvNbRecruit;
            
            if (numSporo[index][yrIndex] > sporoMax[yrIndex])
            {
               sporoMax[yrIndex] = numSporo[index][yrIndex];
               maxSporoIndex[yrIndex] = index;
            }
            if (numRecruit[index][yrIndex] > recruitMax[yrIndex])
            {
               recruitMax[yrIndex] = numRecruit[index][yrIndex];
               maxRecruitIndex[yrIndex] = index;
            }
         }
      } 
          
          
      // Calculate the Pearson correlation coefficient (vs. Brest 05/06)
      if (nbLoops >= endLoop && nbLoops < endLoop+1)
      {
         // Normalise abundance/recruitment values
         for (int i=0; i<12; i++)
         {
            for (int j=0; j<4; j++)
            {
               sporoNorm[i][j] = (float)numSporo[i][j] / sporoMax[j];
               recruitNorm[i][j] = (float)numRecruit[i][j] / recruitMax[j];
            }
         }
         
         // Get mean values (over 4 years) for each month
         for (int i=0; i<12; i++)
         {
            sporoAvg[i] = gsl_stats_mean(sporoNorm[i], 1, 4);
            recruitAvg[i] = gsl_stats_mean(recruitNorm[i], 1, 4);
            sporoAvg_log[i] = log(sporoAvg[i]);
            recruitAvg_log[i] = log(recruitAvg[i]);
         }
         
      }
      
      // Print R^2 values to output file
      if ( nbLoops > (getEngParamsPtr()->nbLoops_-730) )
      {
         // Pearson correlation coefficient
         const size_t stride = 1;
         r2_sporo = pow( gsl_stats_correlation(sporoAvg, stride, brestSporo, stride, 12), 2);
         r2_sporo04_06 = pow( gsl_stats_correlation(sporoAvg, stride, brestSporo04_06, stride, 12), 2);         
         r2_recruit = pow( gsl_stats_correlation(recruitAvg, stride, brestRecruit, stride, 12), 2); 
         r2_recruit04_06 = pow( gsl_stats_correlation(recruitAvg, stride, brestRecruit04_06, stride, 12), 2);
         
         r2_sporo_log = pow( gsl_stats_correlation(sporoAvg_log, stride, brestSporo_log, stride, 12), 2);
         r2_sporo04_06_log = pow( gsl_stats_correlation(sporoAvg_log, stride, brestSporo04_06_log, stride, 12), 2);         
         r2_recruit_log = pow( gsl_stats_correlation(recruitAvg_log, stride, brestRecruit_log, stride, 12), 2); 
         r2_recruit04_06_log = pow( gsl_stats_correlation(recruitAvg_log, stride, brestRecruit04_06_log, stride, 12), 2);  
         
         cout << "Model vs. Brest Harbour data (Voisin, 2007):" << endl;
         cout << "R^2 (Abundance 05/06) = " << r2_sporo << endl; 
         cout << "R^2 (Recruit 05/06) = " << r2_recruit << endl;
         cout << "\nR^2 (Abund 04-06) = " << r2_sporo04_06 << endl;
         cout << "R^2 (Recruit 04-06) = " << r2_recruit04_06 << endl;
         cout << endl;
         cout << "Model vs. Brest Harbour data (Voisin, 2007) (Log values):" << endl;
         cout << "R^2 (Log Abundance 05/06) = " << r2_sporo_log << endl; 
         cout << "R^2 (Log Recruit 05/06) = " << r2_recruit_log << endl;
         cout << "\nR^2 (Log Abund 04-06) = " << r2_sporo04_06_log << endl;
         cout << "R^2 (Log Recruit 04-06) = " << r2_recruit04_06_log << endl;
         cout << endl;  
      
         (*ptrOutputFile) << "R^2 (Abundance 05/06) = " << r2_sporo << "\n"
                          << "R^2 (Recruit 05/06) = " << r2_recruit << "\n"
                          << "R^2 (Abund 04-06) = " << r2_sporo04_06 << "\n"
                          << "R^2 (Recruit 04-06) = " << r2_recruit04_06 << endl;
         (*ptrOutputFile) << "R^2 (Log Abundance 05/06) = " << r2_sporo_log << "\n"
                          << "R^2 (Log Recruit 05/06) = " << r2_recruit_log << "\n"
                          << "R^2 (Log Abund 04-06) = " << r2_sporo04_06_log << "\n"
                          << "R^2 (Log Recruit 04-06) = " << r2_recruit04_06_log << endl;
      }
      
      //
      // End of correlation test
      //
   }  // endif (rank == 0)
   
   
   return rcvNbBac;
}


//
// Simple calculation of loops per second performance of program
//
float Bait::calcLoopsPerSec(unsigned nbLoops)
{
   struct timeval newTime;
   float lps=0.0;

   gettimeofday( &newTime, NULL );
   lps = ( nbLoops - lastLoopNumber_ ) /
      ( ( newTime.tv_sec - lastTime_.tv_sec )
      + (double) ( newTime.tv_usec - lastTime_.tv_usec ) / 1000000 ) ;
   lastTime_        = newTime     ;
   lastLoopNumber_  = nbLoops  ;
   return lps;
}


//
// Display various stats to the console
//
void Bait::outputResults(unsigned nbLoops, int totNumOrg, float lps)
{
   //unsigned loopsPerDay = 24;
   unsigned sampleLoop = getEngParamsPtr()->sampleLoop;    // Loop to print results
   unsigned maxAge = getEngParamsPtr()->graphPeriod;
   float macroSize = baitParamsPtr_->macroSize;
   
   // Retrieve total levels on plate
   currentNumAdult = plate_.countSporophytes(macroSize,100.0);
   numNewRecruits = plate_.countNewRecruits(macroSize, maxAge);
   numNewSporophytes = plate_.countNewSporophytes(baitParamsPtr_->startStock, maxAge);
   numSporophytes = plate_.countSporophytes(macroSize, 100.0);
   numGametophytes = plate_.countGametophytes(false);
   sumAgeMaturity = Organism::getSumAgeMaturity();
   countMature = Organism::getCountMature();
   sumAgeDeath = Organism::getSumAgeDeath();
   countDead = Organism::getCountDead();
   float plateWaterTemp = plate_.getWaterTemp();
   float plateSolarRad = plate_.getSolarRad();
   float plateDayLength = plate_.getCurrDayLength();
   plateBiomass = plate_.calcTotalBiomass();
   plateSporeLevel = plate_.getTotalSporeLevel(mpiSize);
   plateAvgDensity = plate_.getAvgDensity();
   plateUmoy = plate_.getUmoy();
   plateVmoy = plate_.getVmoy();
   plateCoverage = plate_.calcTotalCoverage();

   // Collate results from all nodes for output to file by node 0.
   if (mpiSize > 1)
   {
      totNumOrg = mpiOutputResults(nbLoops, plateBiomass, plateCoverage, colonyPerimeter, plateSporeLevel, lps, plateAvgDensity, plateUmoy, plateVmoy);
      
      // End of program test
      if( rank ==0 && !totNumOrg )
      {
         cout << "No Organism (time = " << time(NULL) << ")" << endl;

         throw myException( "No bacteria : end of program", __FILE__, __LINE__ ) ;
      }
   }
   else
   {      
      // JM - record results in output file
      if (ptrOutputFile != NULL)
      {
         if (nbLoops == 1) 
         {
            (*ptrOutputFile)  << "NumLoops" << "\t" 
            << "NumSporo" << "\t"
            << "NewSporo" << "\t"
            << "Recruit" << "\t"
            << "NumGameto" << "\t"
            << "SporeLevel" << "\t"
            << "AgeMature" << "\t"
            << "AgeDeath" << "\t"
            << "Temp" << "\t"
            << "Solar" << "\t"
            << "DayL" << "\t"
            << "GamGrow" << "\t"
            << "SporGrowMod" << "\t"
            << "GamMatur" << "\t";
            (*ptrOutputFile) << "lps" << endl;
         }
         else
         {
            (*ptrOutputFile)  << nbLoops << "\t" 
               << numSporophytes << "\t"
               << numNewSporophytes << "\t"
               << numNewRecruits << "\t"
               << numGametophytes << "\t"
               << plateSporeLevel << "\t"
               << (float)sumAgeMaturity/countMature << "\t"
               << (float)sumAgeDeath/countDead << "\t"
               << plateWaterTemp << "\t"
               << plateSolarRad << "\t"
               << plateDayLength << "\t"
               << Organism::getAvgGametGrowthMod() << "\t"
               << Organism::getAvgSporoGrowthMod() << "\t"
               << Organism::getProbGamMature() << "\t";
               (*ptrOutputFile) << lps << endl;
         }
      }

      if ( (nbLoops-sampleLoop) % getEngParamsPtr()->graphPeriod == 0 )
      {
         cout << "Diagnostics 0D (Lp, Density, Umoy, Vmoy): " << plateAvgDensity << " " << plateUmoy << " " << plateVmoy << " " << endl;
         
         cout << rank << ": Lp" << nbLoops << ": " 
            << "Adult=" << numSporophytes
            << " Rec=" << numNewRecruits
            << " Gam=" << numGametophytes
            << " Spor=" << plateSporeLevel
            << " Tmp=" << plateWaterTemp
            << " Sol=" << plateSolarRad
            << " DL=" << plateDayLength;
            cout << " (lps=" << lps << ")"
            << endl;
      }
      
      // Print horizontal cross-section of pollen counts
      //if (nbLoops%2000 == 0)
         //plate_.printHorizCross(nbLoops, ptrOutputFile2);
      //   plate_.printPlantStats(nbLoops, ptrOutputFile2);
      
   }  // end else
}

//
// Initiate water temperature, solar radiation & day length parameters
//
void Bait::initEnvParams()
{
   float tempAmp = baitParamsPtr_->tempAmp,
         temp_c = baitParamsPtr_->temp_c,
         temp_d = baitParamsPtr_->temp_d,
         solarAmp = baitParamsPtr_->solarAmp,
         solar_c = baitParamsPtr_->solar_c,
         solar_d = baitParamsPtr_->solar_d,
         DL_Amp = baitParamsPtr_->DL_Amp,
         DL_c = baitParamsPtr_->DL_c,
         DL_d = baitParamsPtr_->DL_d,
         cos_units = baitParamsPtr_->cos_units;
   unsigned day=1;
         
   // Add some noise to the temp/solar amplitude data
   //float ampNoise = sampleGaussian(0.07);
   //tempAmp += (ampNoise * tempAmp);
   //solarAmp += (ampNoise * solarAmp);

   // Cosine curves used to represent seasonal changes in environmental parameters.
   for (int i=0; i<365; i++)
   {
      day = i+1;
      waterTemp[i] = tempAmp * cos( (day-temp_c) * 2*M_PI/cos_units ) + temp_d;
      solarRadiation[i] = solarAmp * cos( (day-solar_c) * 2*M_PI/cos_units ) + solar_d;
      dayLengths[i] = DL_Amp * cos( (day-DL_c) * 2*M_PI/cos_units ) + DL_d;
   }
   
}






// ------------------------------ End Of File ------------------------------- //
