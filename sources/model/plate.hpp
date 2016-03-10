// -------------------------------------------------------------------------- //
// plate.hpp                                                                  //
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

#ifndef PLATE_HPP
#define PLATE_HPP

// -------------------------------------------------------------------------- //
// Libraries                                                                  //
// -------------------------------------------------------------------------- //

#include "../engine/world.hpp"
#include "../engine/fabric.hpp"
#include "patch.hpp"
//#include <movie.h>
//#include <sphere.h>
//#include <frame.h>
#include <math.h>

// -------------------------------------------------------------------------- //
// Forward declarations                                                       //
// -------------------------------------------------------------------------- //

class Organism ;
class Antibio ;
class BaitParameters ;

// -------------------------------------------------------------------------- //
// class Plate                                                                //
// -------------------------------------------------------------------------- //

class Plate : public World<Patch>
{
   public :

      Plate() ;
      virtual ~Plate() ;

      static void initBaitParamsPtr( BaitParameters * ) ;
      void init(int, unsigned, int *, int *, int *, float *, float *, float *) ;

      void solveConflicts(unsigned) ;
      void allPatchesDevelop(unsigned);
      void germinateSpores(unsigned);

      virtual void display( unsigned, Fabric<Organism>&, int, int/*, Fabric<Antibio>&*/ ) {};
      //void animpDisplay(unsigned, Fabric<Organism>&, unsigned);
      //void animpWrite();

      void bacOvercrowding_calcCharge();

      void updateWaterLevel(unsigned);
      void updateWaterSolar(unsigned);
      void updateDailyWaterSolar(unsigned, float *, float *);
      void parallelDiffusion(int, unsigned);      // MPI Parallelisation
      void parallelDiffusion_CUDA(int, unsigned); // MPI && CUDA
      void diffusion(unsigned);
      int parallelMovement(int);
      void noDiffusionMPI(int, unsigned);
      void diffusionCUDA(unsigned);
      //unsigned countOrganisms(float);
      //unsigned countOrganisms(float, float);
      //unsigned countOrgsAge(unsigned, unsigned);
      //unsigned countNewRecruits(float, float);
      unsigned countSporophytes(float, float);
      unsigned countNewRecruits(float, unsigned);
      unsigned countNewSporophytes(float, unsigned);
      unsigned countGametophytes(bool);
           
      int calcTotalCoverage();
      int calcTotalPerimeter();
      
      unsigned randOrgEradication();
      unsigned perimEradication();
      unsigned perimRandEradication();
      unsigned rowByRowEradication(unsigned, unsigned);
      unsigned removeMeadowMPI(unsigned);
      unsigned removeOutliersMPI(unsigned);
       
      void printMap(unsigned, ofstream **, int, unsigned);
      void printHorizCross(unsigned, ofstream *);
      void printPlantStats(unsigned, ofstream *);
      void gatherPlatesMPI();
      double sampleGaussian(double, int);
      
      int calcIPosDistAng(unsigned, float, float);
      int calcJPosDistAng(unsigned, float, float);      
      
      int findPlate(int, int);
      int convGlobalToLocali(int);
      int convGlobalToLocalj(int);
      int convLocalToGlobali(int, int);
      int convLocalToGlobalj(int, int);
      
      // Lattice Boltzmann Methods
      /*void moveLBE(int, unsigned);
      void obstacleLBE(int, unsigned);
      void methodsLBE(unsigned);
      void sendAndReceiveLBE(int, unsigned);
      void diag0D(unsigned);//, ofstream *);
      void configLBE(unsigned, ofstream **, int, unsigned);
      void profileLBE(unsigned, ofstream **, int);
      void addObstacleLBE(int);
      void removeObstacleLBE();
      void allCellsLBE(unsigned);
      void backFacingStep(int);*/
      
      // Get functions
      double calcTotalBiomass();
      double getTotalSporeLevel(unsigned);
      bool checkSporesPresent(unsigned);
      unsigned getMyRank()
         { return myRank; }
      float getWaterLevel()
         { return waterLevel; }
      float getWaterTemp()
         { return currWaterTemp; }
      float getSolarRad()
         { return currSolarRad; }
      float getSolarRadiation()
         { return currSolarRad; } 
      float getCurrDayLength()
         { return currDayLength; }
      float getCurrDay()
         { return currDay; }
      float getDailySolarRad(unsigned i)
         {  return dailySolarRad[i-1]; }
      
      // LBE
      double getOmega()
         { return omega; }
      double getVisc()
         { return visc; }
      double getForce()
         { return force; }
      double getUmoy()
         { return umoy; }
      double getVmoy()
         { return vmoy; }
      double getAvgDensity()
         { return avgDensity; }
         
      void initElevNicheSubstrate();          // Initialise elevation, niche and substrate types
      void initUndariaNiche();
      //void setNicheLimitsSpartina();
      //void setNicheLimits();

   protected :

      void plateDiffusion(unsigned);
      void plateDiffusionNew(unsigned);
      void updateLevels(unsigned);
      void sendAndReceive(int, unsigned);

      unsigned addToBuffer(float*, short*, int*, Agent*, unsigned, unsigned);

      static BaitParameters* baitParamsPtr_ ;
      static const unsigned NUM_ADJ_PLATES; // no. adjacen plates: 8 if 2D, 26 if 3D
      // Starting days of each calendar month
      static const unsigned stMonthsHr[];
      static const unsigned monthLength[]; 
      static const unsigned daysPerMonth[];   

      void setActiveRegion(unsigned);
      void setActivePatches();

      void calcSndBufferSize();
      void calcRcvBufferSize(int);
      void calcSndBufSizeCorners();
      void calcRcvBufSizeCorners(int);

      // Variables required for parallelisation:
      // - start and end indexes (rows(j)/cols(i)) of active region.
      unsigned iStart, iStop, jStart, jStop;
      unsigned totEnvWidth, totEnvHeight;
      unsigned myRank, numPlatesMPI;

      // sendAndReceive() parameters - size of arrays being sent to adjacent plates
      // in parallelDiffusion. If zero, means no adjacent plate in that direction.
      unsigned activeEdgeLength[9];
      bool adjacentPlate[9];
      int dest[9];
      
      // Co-ords for parallel diffusion and parallel movement algorithm
      unsigned iPos_parDiff[9], jPos_parDiff[9],
               iPos_parMove[9], jPos_parMove[9];

      // Size of buffers for sending bacteria between adjacent nodes
      unsigned sndBufSize[9], rcvBufSize[9];
      
      // Variables for creating movie file for "Animp" 3D visualisation
      //movie *m;
      
      // Coastal variables
      int *xy_elevations; int *agentLocations; int *substrateTypes;
      float *allPatchValues, *allBiomassValues, *allElevationValues, *allSporeValues;
      int *substrateValues, *allSubstrateValues;
      float *patchValues, *biomassValues, *elevationValues, *sporeValues;
      //float upperLimit, lowerLimit;    // Elevation limits of S. anglica (in metres OD Newlyn)
      float waterLevel;             // current level of the tide
      
      // Water temperature (degrees celcius), solar radiation (MJ m^-2 h^-1) 
      // and day length (hours, sunrise to sunset) for current loop
      float currWaterTemp,
            currSolarRad,
            currDayLength;
      
      // Inputted daily values for water temperature (degrees celcius), 
      // solar radiation (MJ m^-2 month^-1), and day length (h, sunlight)
      float *dailyWaterTemp,
            *dailySolarRad,
            *dayLengths;

      unsigned currMonth,
               currDay;

      // Lattice Boltzmann Equation (LBE) variables
      // relaxation frequency omega, viscosity and uforce
      double omega, visc, force;
      static const double CS2;
      double umoy, vmoy, avgDensity;
};


// -------------------------------------------------------------------------- //

#endif                        // PLATE_HPP


// ------------------------------ End Of File ------------------------------- //
