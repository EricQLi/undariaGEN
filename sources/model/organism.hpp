// -------------------------------------------------------------------------- //
// organism.hpp                                                               //
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

#ifndef ORGANISM_HPP
#define ORGANISM_HPP


// -------------------------------------------------------------------------- //
// Libraries                                                                  //
// -------------------------------------------------------------------------- //

#include "../engine/agent.hpp"
#include "../engine/fabric.hpp"
#include "baitParameters.hpp"
#include "patch.hpp"
#include <bitset>
using std::bitset;

#define PI 3.14159265

// -------------------------------------------------------------------------- //
//   Enumeration
// -------------------------------------------------------------------------- //

// JM 13Jun2006
enum Species { U_PINNATIFIDA, U_UNDARIOIDES };
enum TempType { NONE_T, PSYCHRO, MESO, THERMO };

// -------------------------------------------------------------------------- //
// Forward declarations                                                       //
// -------------------------------------------------------------------------- //

//class Molecule ;
class Patch;


// -------------------------------------------------------------------------- //
// class Organism                                                             //
// -------------------------------------------------------------------------- //

class Organism : public Agent
{
   public :

      Organism() ;
      virtual ~Organism() ;
      
      const Organism& operator=(const Organism &organismToCopy);

      static void initStaticPtr( const BaitParameters *, const EngineParameters *, Fabric<Organism>* ) ;
      static void initWorldPtr( World<Patch> * );

      static double sampleGaussian(double); 
      static void updGametGrowth(float, float, float);
      static void updSporoGrowthMod(float, float, float);
      static void updProbGamMature(float, float);
      static double heightOfWeibull(double, double, double);
      
      void init(int, int) ;
      void init(int, int, double) ;

      void move() ;
      void growth(unsigned) ;
      void gametophyteGrowth();
      void produceSpore();
      bool prematureDeath();

      bool isAlive()
         { if( stock_ > 1.0E-7 ) return true ; return false ; }

      void tryToReproduce();
      void deleteOrganism();
      void setNewTraits(double, bool, int);     
      void initAttributes();
      
      Direction_t nextDirection() ;
      void newRandomDirection() ;
      Direction_t getNextTowardDir() ;
      
      double phTempEffect(double);
      double weibullDistribution(float, TempType);
      int round(double);
      int round(float);

      // Get functions
      float getStock() const
         { return stock_ ; }
      float getSporeStock() const
         { return sporeStock; }
      bool getReleaseSpore() const
         { return releaseSpore; }
      unsigned getType() const
         { return agentType_ ; }
      unsigned getSeedlingID() const
         { return seedlingID; } 
      unsigned getAge() const
         { return age; }
      unsigned getAgeRecruit() const
         { return ageRecruit; } 
      unsigned getAgeGameto() const
         { return ageGameto; } 
      bool getIsGametophyte() const
         { return gametophyte; }
      bool getGamMature() const
         { return gamMature; }
                 
      static Fabric<Organism>* getFabricPtr()
         { return fabricPtr_; }
      static float getGametGrowth()
         { return gametGrowth; }
      static float getSporoGrowthMod()
         { return sporoGrowthMod; }
      static double getProbGamMature()
         { return probGamMature; }
      static int getNumSeedlings()
         { return numSeedlings; }
      static int* getColonySizes()
         { return colonySizes; }
      static int* getColonyIDs()
         { return colonyIDs; }
      static int getSumAgeMaturity()
         { return sumAgeMaturity; }
      static int getCountMature()
         { return countMature; }
      static int getSumAgeDeath()
         { return sumAgeDeath; }
      static int getCountDead()
         { return countDead; }
      static float getAvgGametGrowthMod()
         { return avgGametGrowthMod; }
      static float getAvgSporoGrowthMod()
         { return avgSporoGrowthMod; }

      
      // Set functions
      void setStock(float stk)
         { stock_ = stk; }
      void setSporeStock(float stk)
         { sporeStock = stk; }
      void setReleaseSpore(bool relSpore)
         { releaseSpore = relSpore; }
      void setShuntDir(Direction_t dir)
         { shuntDir = dir; }  
      void setSeedlingID(unsigned id)
         { seedlingID = id; }
      void resetAge()
      { 
         age = 0; 
         ageRecruit = 0;
         ageMaturity = 0;
         ageGameto = 0;
      }
      void setGametophyte(bool inGamet)
         { gametophyte = inGamet; }

   private :
   
      static const BaitParameters *baitParamsPtr_ ;
      static const EngineParameters *engParamsPtr_ ;
      static Fabric<Organism>     *fabricPtr_     ;
      static World<Patch> *ptrWorld;
      
      static AgentType_t          agentType_      ;      
      // To calculate mean age at maturity:
      static int sumAgeMaturity;
      static int countMature;
      static int sumAgeDeath;
      static int countDead;
      static int numSeedlings;
      static int *colonySizes;
      static int *colonyIDs;
      static float gametGrowth;
      static float sporoGrowthMod;
      static double probGamMature;
      static float sumSporoGrowthMod;
      static unsigned countSporoGrowthMod;
      static float avgSporoGrowthMod;
      static float avgGametGrowthMod;
     
      float stock_,
            sporeStock;
            
      short unsigned age, // number of loops as macroscopic sporophyte
              ageRecruit,
              ageGameto,
              ageMaturity; // age at maturity
      unsigned seedlingID; // unique identifier for agents that come from same original seedling
   
      // store traits of organism using bit fields to save memory.
      bool gametophyte:1;    
      bool releaseSpore:1;
      bool gamMature:1;
      TempType tempType: 2;           // thermophile, mesophile or psychrophile
      Species species: 1;             // U. pinnatifida, U. undarioides etc.
      Direction_t shuntDir:4;
    
};


// -------------------------------------------------------------------------- //

#endif                        // ORGANISM_HPP


// ------------------------------ End Of File ------------------------------- //
