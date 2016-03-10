// -------------------------------------------------------------------------- //
// bait.hpp                                                                   //
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

#ifndef BAIT_HPP
#define BAIT_HPP


// -------------------------------------------------------------------------- //
// Libraries                                                                  //
// -------------------------------------------------------------------------- //

#include "organism.hpp"
#include "../interface/xPlate.hpp"
#include "../engine/system.hpp"
#include "../engine/fabric.hpp"
#include <time.h> 
#include <fstream> 
using std::ofstream;
#include <gsl/gsl_rng.h>

// -------------------------------------------------------------------------- //
// Forward declarations                                                       //
// -------------------------------------------------------------------------- //

class Patch ;
class Parameters ;
class BaitParameters ;


// -------------------------------------------------------------------------- //
// class Bait                                                                 //
// -------------------------------------------------------------------------- //

class Bait : public System
{
   public :

      Bait() ;
      virtual ~Bait() ;

      static void initBaitParamsPtr( const BaitParameters * ) ;
      void init(Parameters &, ofstream *, ofstream *, ofstream *outputFilesLB[], int, ofstream *outputFilesLB2[], int, ofstream *outputFilesMG[], int, int *, int *, int *) ;
      void initEnvParams();

      void run() ;

      void displayTime( unsigned ) ;
      void latticeBoltzmannTime (unsigned);
      void worldEvolutionTime(unsigned)    ;
      void staticActivityTime(unsigned)    ;
      void startMoveTime();
      void moveTime()              ;
      void solveConflictsTime(unsigned)    ;
      void interactionsTime(unsigned)      ;
      void reproductionTime(unsigned)      ;

      void printDuration(double);
      static gsl_rng* rng;
      void evenlyDistributeAgents();
      void randDistributeAgents();
      
      float calcLoopsPerSec(unsigned);
      void outputResults(unsigned, int, float);
      int calcPollenMonth(unsigned);
      double sampleGaussian(double);
      void convToMJm2h1(float *);


   protected :

      int mpiOutputResults(unsigned, double, int, int, double, float, double, double, double);

      static const BaitParameters* baitParamsPtr_ ;
      static ofstream* ptrOutputFile;
      static ofstream* ptrOutputFile2;
      static const unsigned daysPerMonth[]; 
               
      ofstream **outFilesLB, **outFilesLB2, **outFilesMG;
      int numOutFilesLB, numOutFilesLB2, numOutFilesMG;

      XPlate             plate_           ;
      Fabric<Organism>   bacteriasFabric_ ;
      clock_t            start, 
                         stop;
                 
      float waterTemp[365];
      float solarRadiation[365];    // in MJ m^-2 month^-1  
      float dayLengths[365];
      
      // Arrays for doing correlation test with Brest harbour data in Bait::mpiOutputResults()
      unsigned numSporo[12][4],
               numRecruit[12][4];         
      double   sporoNorm[12][4],
               recruitNorm[12][4];
      double   sporoAvg[12], 
               recruitAvg[12],
               sporoAvg_log[12], 
               recruitAvg_log[12];
      unsigned sporoMax[4], recruitMax[4], startLoopTest;
      int      index, yrIndex, maxSporoIndex[4], maxRecruitIndex[4];
      float    r2_sporo,
               r2_sporo04_06,
               r2_recruit,
               r2_recruit04_06,
               r2_sporo_log,
               r2_sporo04_06_log,
               r2_recruit_log,
               r2_recruit04_06_log;

      // JM - variables for tracking time taken for each function
      double   displayDuration, worldEvolveDuration, 
               boltzmannDuration, staticActivityDuration,
               startMoveDuration, 
               moveDuration, solveConflictsDuration, 
               interactionsDuration, reproductionDuration;
      
      // JM - variables to track totals of various constituent entities.
      double   plateBiomass,
               plateSporeLevel,
               plateAvgDensity,              // LBE variables
               plateUmoy,
               plateVmoy,
               platePollenRel,
               platePollenRcv,
               plateSeedStock;
      int      plateCoverage,
               colonyPerimeter,
               sumAgeMaturity,
               countMature,
               sumAgeDeath,
               countDead;
               
      double   *plateInactAbLevel,
               *plateAbActLevel,
               *plateMeanAbConc;     // micrograms/ml

      // JM 22Aug2006 - Variables for calculating loops per second
      unsigned lastLoopNumber_,
               highestNumBac,
               abInhibitedNum,
               currentNumBac,         // keep record of no. bacteria on plate
               currentNumAdult,       // num adult organisms
               numNewRecruits,        // num new recruits (Voisin definition)
               numGametophytes,
               numNewSporophytes,
               numSporophytes;
               
      
      bool     sporesPresent;       
                      
      // Variables for calculating colony growth rate
      float    cover_y1, 
               cover_y2;

      struct timeval  lastTime_;

      // MPICH variables
      int      err,
               rank,
               mpiSize;

};


// -------------------------------------------------------------------------- //

#endif                        // SYSTEM_HPP


// ------------------------------ End Of File ------------------------------- //
