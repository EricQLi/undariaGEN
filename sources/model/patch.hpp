// -------------------------------------------------------------------------- //
// patch.hpp                                                                  //
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

#ifndef PATCH_HPP
#define PATCH_HPP


// -------------------------------------------------------------------------- //
// Libraries                                                                  //
// -------------------------------------------------------------------------- //

#include "../engine/cell.hpp"
#include "../engine/world.hpp"
#include "../engine/fabric.hpp"
#include <fstream>      // JM
using std::ofstream;


// -------------------------------------------------------------------------- //
// Forward declarations                                                       //
// -------------------------------------------------------------------------- //

class Organism ;


// -------------------------------------------------------------------------- //
// class Patch                                                                //
// -------------------------------------------------------------------------- //

class Patch : public Cell
{
   public :

      Patch() ;
      virtual ~Patch() ;

      static void initBaitParamsPtr( BaitParameters * , EngineParameters *) ;
      static void initWorldPtr( World<Patch> * );
      void init(World<Patch> *, int) ;
      void init(World<Patch> *, int, int);
      //void initLBE();
      
      void development(unsigned);

      void addOrganism( Organism* ) ;
      void moveOrganism( Organism*, Direction_t) ;
      void addSpore(float);
      void diffusion(unsigned);
      void nutDiffusion(unsigned);
      void seedDiffusion(unsigned);
      void sporeDiffusion(unsigned);
      void updateNutLevel();
      void updateSporeLevel();
      void removeOrgFromPatch(float);
      void solveConflicts(unsigned) ;
      Patch* getNeighbour(Direction_t, bool);
      Patch* getPatch(unsigned, unsigned);
      void germinateSpores(unsigned, float, unsigned);

      void newOrganismMPI(float, short, int);
      void removeAllOrganisms();

      unsigned countNumOrganism();
      unsigned countNumSporophytes(float, float);
      unsigned countNumSporophytes();
      unsigned countNewRecruits(float, unsigned);
      unsigned countNewSporophytes(float, unsigned);
      unsigned countNumGametophytes();
      unsigned countMatureGametophytes();
      double sampleGaussian(double);     
      void updateAnyFreeNeighbour();
      

      // Get functions
      bool getIsOrganism() const
         { return isOrganism; }
      unsigned getNumOrg() const
         { return numOrg; }
      float getSporeLevel()
         { return sporeLevel; }
      float getSporoBiomass()
         { return sporoBiomass; }
      float calcSporoBiomass();
      bool getActiveRegion() const
         { return activeRegion; }
      bool getAnyFreeNeighbour() const
         { return anyFreeNeighbour; }
      bool getIsUndariaNiche() const
         { return isUndariaNiche; }
      float getElevation() const
         { return elevation;}
      int getSeedlingID();
      float getWaterDepth()
         { return waterDepth; }
      int getSubstrateType()
         { return substrateType; }


      // Set functions
      void setNumOrg(unsigned inNum)
         { numOrg = inNum; }
      void setSporeLevel(float inSporeLevel)
         { sporeLevel = inSporeLevel; }
      void setSporoBiomass(float inBiomass)
         { sporoBiomass = inBiomass; }
      void setActive(bool activeState)
         { activeRegion = activeState; }
      //void setIsObstacle(bool);
      //void setIsTopObst(bool inIsTopObst)
      //   { isTopObst = inIsTopObst; }
      //void setIsBottObst(bool inIsBottObst)
      //   { isBottObst = inIsBottObst; }         
      void setElevation(float);
      void updateWaterDepth(float);
      void setSubstrateType(int inSubstrate)
         { substrateType = inSubstrate; }
      void setIsUndariaNiche();
              
      // Lattice Boltzmann methods (adapted from S. Succi, lbe.f)
      //void equili();
      //void initPop();
      //void moveLBE(EngineParameters *, unsigned);
      //void obstacle(EngineParameters *, unsigned);
      //void hydrovar(unsigned);
      //void collis(double, unsigned);
      //void updateF(unsigned);
      //void updateObst(unsigned);
      //void forceLB(double);
      //void fluidFlowLBE(unsigned);
      //void printLBEDetails(ofstream *);
      
      // LBE Get functions   
      /*double getF(int index) const
         { return f[index]; }    
      double getRho_density() const
         {  return rho_density; }
      double getU_vel() const
         { return u_vel; }
      double getV_vel() const
         { return v_vel; }
      static int getNumDir()
         { return numDir; }
      bool getIsTopObst() const
         { return isTopObst; }
      bool getIsBottObst() const
         { return isBottObst; }
      bool getIsStep() const
         { return isStep; }
      Patch* getNeighbourLBE(Direction_t, EngineParameters *);*/
      
      // LBE Get functions
      /*void setF(double in_f, int index)
         { f[index] = in_f; }
      void setU_vel(double inU)
         { u_vel = inU; }
      void setV_vel(double inV)
         { v_vel = inV; }
      void setIsStep(bool inIsStep)
         { isStep = inIsStep; }*/


   private :
      static const int numDir = 9;
      static const int NUM_SPORE_AGES = 1;
      static BaitParameters* baitParamsPtr_ ;
      static EngineParameters* engParamsPtr_ ;
      static World<Patch> *ptrWorld;
           
      float  waterDepth;
             
      float sporeLevel,
            newSporeLevel;    // num spores

      float  temp,
             sporoBiomass;
                   
      float    elevation;       // elevation above sea level in metres                   
      short int      substrateType;
      short unsigned numOrg;
      bool     activeRegion,
               anyFreeNeighbour,
               isUndariaNiche,
               isOrganism;

      
      //
      // ***Lattice Boltzmann variables***
      //
      /*float u_vel,              // velocity in X direction
             v_vel,              // velocity in Y direction 
             rho_density;        // density
             
      double *f_eq,               // equilibrium values for each direction (9)
             *f,                 // freq of particles for each direction (9)
             *new_f;             // 0=0, 1=E, 2=N, 3=W, 4=S, 5=NE, 6=NW, 7=SW, 8=SE
                 
      static const double cs2,   // sound speed constants      
                          cs22,
                          cssq;
      
      static const double wt0,   // lattice weights
                          wt1,
                          wt2;
                          
      bool isTopObst, isBottObst, isStep; // isObstacle*/
                                                  
};


// -------------------------------------------------------------------------- //

#endif                        // PATCH_HPP


// ------------------------------ End Of File ------------------------------- //
