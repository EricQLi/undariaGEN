// -------------------------------------------------------------------------- //
// patch.cpp                                                                  //
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

#include "patch.hpp"


// -------------------------------------------------------------------------- //
// Libraries                                                                  //
// -------------------------------------------------------------------------- //

#include "organism.hpp"
#include "../engine/myException.hpp"
#include <cstdlib>
#include <cmath>
using std::fabs;
#include <iostream>
using std::cout;
using std::cin;
using std::endl;
using std::scientific;
#include <iomanip>
using std::setprecision;
#include "bait.hpp"
#include <gsl/gsl_randist.h>


// -------------------------------------------------------------------------- //
// Static member creation                                                     //
// -------------------------------------------------------------------------- //

BaitParameters* Patch::baitParamsPtr_ = 0 ;
EngineParameters* Patch::engParamsPtr_ = 0 ;
World<Patch> * Patch::ptrWorld = 0 ;

// LBE variables
/*double const Patch::cs2 = 1.0/3.0;
double const Patch::cs22 = 2.0 * (1.0/3.0);
double const Patch::cssq = 2.0/9.0;
double const Patch::wt0 = 4.0/9.0;
double const Patch::wt1 = 1.0/9.0;
double const Patch::wt2 = 1.0/36.0;*/


// -------------------------------------------------------------------------- //
// Constructor                                                                //
// -------------------------------------------------------------------------- //

Patch::Patch()
: Cell()
{
   elevation = 0.0;
   sporoBiomass = 0.0;
   isOrganism = false;
   numOrg = 0;
   activeRegion = false;
   
   // LBE:
   //u_vel = 0.0;
   //v_vel = 0.0;
   //rho_density = 0.0;
   //isTopObst = false;
   //isBottObst = false;
}


// -------------------------------------------------------------------------- //
// Destructor                                                                 //
// -------------------------------------------------------------------------- //

Patch::~Patch()
{
}


// -------------------------------------------------------------------------- //
// Initialiaze the static parameters pointer                                  //
// -------------------------------------------------------------------------- //

void Patch::initBaitParamsPtr( BaitParameters *inPtr, EngineParameters *inEngPtr)
{
   baitParamsPtr_ = inPtr ;
   engParamsPtr_ = inEngPtr;   
}


// -------------------------------------------------------------------------- //
// Initialiaze the static parameters pointer                                  //
// -------------------------------------------------------------------------- //

void Patch::initWorldPtr( World<Patch> * inWorld )
{
   ptrWorld = inWorld;
}


// -------------------------------------------------------------------------- //
// Init the patch                                                             //
// -------------------------------------------------------------------------- //

void Patch::init(World<Patch> * inWorld, int mpiSize)
{
   anyFreeNeighbour = true;
   waterDepth = 15.0;
   substrateType = 0;
   isUndariaNiche = false;
   
   sporeLevel = 0.0;
   newSporeLevel = 0.0;

   // initialise Lattice Boltzmann variables
   //if (baitParamsPtr_->initUvelocity > 0.0 || baitParamsPtr_->initVvelocity > 0.0)
   //   initLBE();
}


//
// Apply half-life to spores in patches
//
void Patch::development(unsigned nbLoops)
{
   float currSpore = getSporeLevel();
   
   if (currSpore > 0.0 && getActiveRegion() && baitParamsPtr_->sporeHalfLife > 0.0)
   {
      if (currSpore > 10.0)
      {
         float sporeDegradeFactor = pow( 0.5, (1/baitParamsPtr_->sporeHalfLife) );
         setSporeLevel(currSpore * sporeDegradeFactor);
        
         if (getSporeLevel() < 1.0)
            setSporeLevel(0.0);
      }
      else
         setSporeLevel(0.0);
   }

}


//
// 23Jan2014 - prob of "new" agent proportional to num spores in patch & type of substrate
//
void Patch::germinateSpores(unsigned index, float temp, unsigned myRank)
{
   // if appropriate substrate for attachment/establishment present then can germinate
   if ( getIsUndariaNiche() )
   {
      float maxOrgInPatch = 10000;
      if (getSporeLevel() > 10.0 && numOrg < maxOrgInPatch && getActiveRegion() == true)
      {
         unsigned numNewAgents = 0;
         // Prob proportional to sporeLevel
         double probGerminate = baitParamsPtr_->probGerminate;
         
         // Density modifier
         double densityMod = 1.0;
         if(numOrg > 1)
            densityMod -= log(numOrg)/log((double)maxOrgInPatch);
         if(densityMod < 1.0E-06)
            densityMod = 0.0;
         probGerminate *= densityMod;
        
         numNewAgents = (unsigned)(probGerminate * sporeLevel);
         
         // If not enough spores to calculate directly
         if (numNewAgents < 1)
         {
            probGerminate = 1.0 - pow( (1.0-probGerminate), sporeLevel);
            if (gsl_rng_uniform(Bait::rng) < probGerminate)
               numNewAgents = 1;
         }

         // Check that do not add too many agents to patch
         if (numOrg >= maxOrgInPatch)
            numNewAgents = 0;
         else if (numOrg + numNewAgents > maxOrgInPatch)
            numNewAgents = maxOrgInPatch - numOrg;

         // Modify probability of establishment based on substrate type
         switch(getSubstrateType())
         {
            case 1:
               break;
            case 2:
               probGerminate *= 0.5;
               break;   
         }
         
         // Create new agent
         for (unsigned i=0; i<numNewAgents;i++)
         {
            Fabric<Organism> *fabricPtr;
            Organism   *ptrOrganism ;

            // Create and initialise new organism (from fabric array).
            fabricPtr = ((Organism*)getAgents())->getFabricPtr();
            ptrOrganism = fabricPtr->newAgent() ;

            // Initialise parameters for new agent
            if (ptrOrganism)
            {
               ptrOrganism->init(index, myRank);
               
               // Insert the organisms to the patch
               (*this).addOrganism(ptrOrganism);     
            }
            else
               cout << "Patch::germinateSpores(): Problem initialising new agent." << endl;

            ptrOrganism = NULL;
            fabricPtr= NULL;
         }
         
         setSporeLevel(sporeLevel - numNewAgents);
         
         if (numOrg > maxOrgInPatch)
            cout << numOrg << " agents in patch (>" << maxOrgInPatch << ") - " << numNewAgents << endl;
      }
   }
}

// -------------------------------------------------------------------------- //
// Add the organism to the patch. This method provides                        //
// additional instructions when an organism comes to a new patch.             //
// -------------------------------------------------------------------------- //

void Patch::addOrganism( Organism *inOrganism )
{
   // Add the bacteria to the cell
   Cell::addAgent( inOrganism ) ;

   isOrganism = true;
   setNumOrg(getNumOrg()+1);
}


// -------------------------------------------------------------------------- //
// Move organism and invoke addOrganism if organism actually moves (not NONE) //
// -------------------------------------------------------------------------- //

void Patch::moveOrganism( Organism *inOrganism, Direction_t inDir)
{
   if (inDir == NONE)
      Cell::addAgent( inOrganism ) ;   // Just add bacteria to the local list
   else
   {
      Patch *ptrPatch=NULL;
      
      ptrPatch = getNeighbour(inDir, true);
      ptrPatch->addOrganism(inOrganism);

      removeOrgFromPatch(inOrganism->getStock());
      ptrPatch = NULL;
   }

   countNumOrganism();
}




//
// JM 14Nov2013 - add to level of spores (of age 0) in patch
//
void Patch::addSpore(float spore)
{
   sporeLevel += spore;
}


//
// To Do: Solve conflict function
//
void Patch::solveConflicts(unsigned nbLoops)
{
   // [Empty]
}


//
// JM 13-17Dec2013 - Fick's law of diffusion for spore transport (all age categories)
//
void Patch::sporeDiffusion(unsigned nbLoops)
{
   // when tide is out there is no diffusion
   if (waterDepth > 0.0 && getSubstrateType() > -1)
   {
      Patch *ptrPatch=NULL;
      // Difference in particle concentrations between patches 
      float particlesDiff[9];//[2];  
      float diffusionRate = baitParamsPtr_->diffusionRate;
      float diffusionRateDiag = diffusionRate * (sqrt(2.0)/2);

      newSporeLevel = sporeLevel;

      // 1. Calculate differences in spore levels
      // Loop from index 1 (index 0 is current patch)
      // Directions: 0=NONE, 1=E, 2=N, 3=W, 4=S, 5=NE, 6=NW, 7=SW, 8=SE 
      for (int j=1; j<numDir; j++)
      {
         ptrPatch = getNeighbour((Direction_t) j, true);     

         particlesDiff[j] = ptrPatch->getSporeLevel() - getSporeLevel();
      }

      // 2. Apply diffusion algorithm according to Fick's first law.
      for (int j=1; j<(numDir+1)/2; j++)
         newSporeLevel += diffusionRate*particlesDiff[j];

      for (int j=5; j<numDir; j++)
         newSporeLevel += diffusionRateDiag * particlesDiff[j];
   }    
}



//
// JM 13Dec2013 - Update spore levels for each age category in patch.
//
void Patch::updateSporeLevel()
{
   setSporeLevel(newSporeLevel);
}


//
//   JM 06Jun2006 - Return pointer to neighbour patch.
//
Patch* Patch::getNeighbour(Direction_t dir, bool avoidObst)
{
      Patch *ptrPatch=NULL;

      switch (dir)
      {
         case N:
            ptrPatch = &(*ptrWorld).smartAccess(getIPos()-1,getJPos());
            break;
         case NE:
            ptrPatch = &(*ptrWorld).smartAccess(getIPos()-1,getJPos()+1);
            break;
         case E:
            ptrPatch = &(*ptrWorld).smartAccess(getIPos(),getJPos()+1);
            break;
         case SE:
            ptrPatch = &(*ptrWorld).smartAccess(getIPos()+1,getJPos()+1);
            break;
         case S:
            ptrPatch = &(*ptrWorld).smartAccess(getIPos()+1,getJPos());
            break;
         case SW:
            ptrPatch = &(*ptrWorld).smartAccess(getIPos()+1,getJPos()-1);
            break;
         case W:
            ptrPatch = &(*ptrWorld).smartAccess(getIPos(),getJPos()-1);
            break;
         case NW:
            ptrPatch = &(*ptrWorld).smartAccess(getIPos()-1,getJPos()-1);
            break;
         case NONE:
            ptrPatch = &(*ptrWorld).smartAccess(getIPos(),getJPos());
            break;
      }
      
      // JM 15May11: Not possible to move to "obstacle" patches
      if (ptrPatch->getSubstrateType() == -1 && avoidObst == true)
         ptrPatch = &(*ptrWorld).smartAccess(getIPos(),getJPos());
         
      return ptrPatch;
}

//
//   JM 23APR11: LBE - Special case for Lattice Boltzmann Equation
//
/*Patch* Patch::getNeighbourLBE(Direction_t dir, EngineParameters *engParamsPtr_)
{
      Patch *ptrPatch=NULL;

		if( engParamsPtr_->worldType_ == CLOSED || engParamsPtr_->worldType_ == OPEN_EW )
		{
			if (getIPos() > 0 && getIPos() < engParamsPtr_->height_-1)
				ptrPatch = getNeighbour(dir, false);
		   // if at the border (top/bottom) then special case: can't move
			else if (getIPos() == 0)
			{
				if (dir==N || dir ==NE || dir == NW)
					ptrPatch = this;
				else
					ptrPatch = getNeighbour(dir, false);
			}
			else if (getIPos() == engParamsPtr_->height_-1)
			{
				if (dir==S || dir ==SE || dir == SW)
					ptrPatch = this;
				else
					ptrPatch = getNeighbour(dir, false);
			}
      }

		
      return ptrPatch;
}*/


//
// JM 03Aug2012 - Return pointer to patch with coords i,j
//
Patch* Patch::getPatch(unsigned i, unsigned j)
{
      Patch *ptrPatch=NULL;

      ptrPatch = &(*ptrWorld).smartAccess(i,j);
      
      // JM 15May11: Not possible to move to "obstacle" patches
      //if (ptrPatch->getIsObstacle() == true && avoidObst == true)
         //ptrPatch = &(*ptrWorld).smartAccess(getIPos(),getJPos());
         
      return ptrPatch;
}



float Patch::calcSporoBiomass()
{
   Agent *ptrOrganism = getAgents();
   sporoBiomass = 0.0;
   
   // Go through list of bacteria in patch
   while (ptrOrganism)
   {                   
      if( ptrOrganism->getType() == ORGANISM && ptrOrganism->isAlive() )
      {
         if ( ((Organism*)ptrOrganism)->getIsGametophyte() == false )
            sporoBiomass += ((Organism*)ptrOrganism)->getStock();
      }
      
      ptrOrganism = ptrOrganism->getNext();
   }
   return sporoBiomass;
}

//
// JM 07Mar2007 - Remove organism from patch.
//
void Patch::removeOrgFromPatch(float inStock)
{
   Fabric<Organism> *fabricPtr;
   fabricPtr = ((Organism*)getAgents())->getFabricPtr();

   if (fabricPtr->getFirstRemoved() == false)
      fabricPtr->setFirstRemoved(true);

   setNumOrg(getNumOrg()-1);

   if (getNumOrg() < 1)
      isOrganism = false;
}


//
// 12Mar2007 - MPI - create new organism with supplied parameters
//               **Note - Only partially completed**
//
void Patch::newOrganismMPI(float inStock, short expSeed/*, short sos, short abDeath*/, int seedlingID)
{
   Fabric<Organism> *fabricPtr;
   Organism   *ptrOrganism ;

   // New organism created.
   fabricPtr = ((Organism*)getAgents())->getFabricPtr();
   ptrOrganism = fabricPtr->newAgent() ;

   // Set up parameters passed from other plate.
   ptrOrganism->setNewTraits(inStock, (bool)expSeed, seedlingID);     
   (*this).addOrganism(ptrOrganism);

   ptrOrganism = NULL;
   fabricPtr= NULL;
}

//
// 12Mar2007 - Remove all bacteria in patch
//
void Patch::removeAllOrganisms()
{
   Agent *ptrOrganism;   
   ptrOrganism = getAgents();

   // Go through list of bacteria in patch
   while (ptrOrganism)
   {
      if( ptrOrganism->getType() == ORGANISM )
         ((Organism*)ptrOrganism)->deleteOrganism();  

      ptrOrganism = ptrOrganism->getNext();
   }
   isOrganism = false;
   setNumOrg(0);
   setSporoBiomass(0.0);
}





//
// JM 06Feb10 - count and return number of bacteria present
//
unsigned Patch::countNumOrganism()
{
   Agent *ptrOrganism = getAgents();

   unsigned count =0;
   
   // Go through list of bacteria in patch
   while (ptrOrganism)
   {
      if( ptrOrganism->getType() == ORGANISM && ptrOrganism->isAlive() )
         count++;
      ptrOrganism = ptrOrganism->getNext();
   }
   if (count > 0)
      isOrganism = true;
   setNumOrg(count);

   return count;
}


//
// 12Aug14 - Count the number of sporophytes
//
unsigned Patch::countNumSporophytes(float minStock, float maxStock)
{
   Agent *ptrOrganism = getAgents();
   unsigned count =0;
   
   // Go through list of bacteria in patch
   while (ptrOrganism)
   {
      if( ptrOrganism->getType() == ORGANISM && ptrOrganism->isAlive() )
      {
         if (((Organism*)ptrOrganism)->getIsGametophyte() == false
            && ((Organism*)ptrOrganism)->getStock() > minStock
            && ((Organism*)ptrOrganism)->getStock() < maxStock )
            count++;
      }
      ptrOrganism = ptrOrganism->getNext();
   }
   if (count > 0)
      isOrganism = true;

   return count;
}

unsigned Patch::countNumSporophytes()
{
   Agent *ptrOrganism = getAgents();
   unsigned count=0;
   
   // Go through list of bacteria in patch
   while (ptrOrganism)
   {
      if( ptrOrganism->getType() == ORGANISM
            && ptrOrganism->isAlive()
            && ((Organism*)ptrOrganism)->getIsGametophyte() == false )
         count++;
      ptrOrganism = ptrOrganism->getNext();
   }
   if (count > 0)
      isOrganism = true;

   return count;
}


//
// 12Aug14 - Count the number of sporophytes
//
unsigned Patch::countNewRecruits(float minStock, unsigned maxAge)
{
   Agent *ptrOrganism = getAgents();
   unsigned count =0;
   
   // Go through list of bacteria in patch
   while (ptrOrganism)
   {
      if( ptrOrganism->getType() == ORGANISM && ptrOrganism->isAlive() )
      {
         if (((Organism*)ptrOrganism)->getIsGametophyte() == false
            && ((Organism*)ptrOrganism)->getStock() > minStock
            && ((Organism*)ptrOrganism)->getAgeRecruit() < maxAge )
            count++;
      }
      ptrOrganism = ptrOrganism->getNext();
   }
   if (count > 0)
      isOrganism = true;

   return count;
}


unsigned Patch::countNewSporophytes(float minStock, unsigned maxAge)
{
   Agent *ptrOrganism = getAgents();
   unsigned count =0;
   
   // Go through list of bacteria in patch
   while (ptrOrganism)
   {
      if( ptrOrganism->getType() == ORGANISM && ptrOrganism->isAlive() )
      {
         if (((Organism*)ptrOrganism)->getIsGametophyte() == false
            && ((Organism*)ptrOrganism)->getStock() > minStock
            && ((Organism*)ptrOrganism)->getAge() < maxAge )
            count++;
      }
      ptrOrganism = ptrOrganism->getNext();
   }
   if (count > 0)
      isOrganism = true;

   return count;
}


//
// 12Aug14 - Count the number of gametophytes
//
unsigned Patch::countNumGametophytes()
{
   Agent *ptrOrganism = getAgents();
   unsigned count =0;
   
   // Go through list of bacteria in patch
   while (ptrOrganism)
   {
      if( ptrOrganism->getType() == ORGANISM && ptrOrganism->isAlive() )
      {
         if (((Organism*)ptrOrganism)->getIsGametophyte() == true)
            count++;
      }
      ptrOrganism = ptrOrganism->getNext();
   }
   if (count > 0)
      isOrganism = true;

   return count;
}


unsigned Patch::countMatureGametophytes()
{
   Agent *ptrOrganism = getAgents();
   unsigned count=0;
   
   // Go through list of bacteria in patch
   while (ptrOrganism)
   {
      if( ptrOrganism->getType() == ORGANISM && ptrOrganism->isAlive() )
      {
         if (((Organism*)ptrOrganism)->getIsGametophyte() == true
               && ((Organism*)ptrOrganism)->getGamMature() )
            count++;
      }
      ptrOrganism = ptrOrganism->getNext();
   }
   if (count > 0)
      isOrganism = true;

   return count;
}



//
// JM 02Nov11: calculate elevation (in metres OD)
//
void Patch::setElevation(float inElev)
{   
   // For testing purposes just create simple gradient along y-axis from -20.0 to +10 m elevation
   //float totalRange = 20.0;
   //float underWaterRange = 15.0;   
   //elevation = (getIPos()/((float)engParamsPtr_->height_/totalRange)) - underWaterRange;
  
   //float horizMod = (float) getJPos() / engParamsPtr_->width_;   
   //elevation = elevation + (fabs(underWaterRange-elevation) * horizMod/2);
   
   elevation = -5.0;
}

//
// JM 12Dec13: Calculate the current water depth (in m) based on waterLevel and elevation (in OD)
//
void Patch::updateWaterDepth(float waterLevel)
{
   waterDepth = waterLevel - elevation;
}


//
// Check all neighbouring patches to see if any are empty, if so, return true
//
void Patch::updateAnyFreeNeighbour()
{
   Patch *ptrPatch=NULL;
   anyFreeNeighbour = false;
   
   for (int i=1; i<numDir; i++) 
   {
      ptrPatch = getNeighbour((Direction_t)i, true);
      if ( ptrPatch != this && ptrPatch->getIsOrganism()==false)
      {
         anyFreeNeighbour = true;
         break;
      }
   }
}



//
// Get seedling ID of first agent in patch
//
int Patch::getSeedlingID()
{
   int seedlingID = 0;
   
   Agent *ptrOrganism = getAgents();
   while (ptrOrganism)
   {
      if( ptrOrganism->getType() == ORGANISM && ptrOrganism->isAlive() )
      {
         seedlingID = ((Organism*)ptrOrganism)->getSeedlingID();
         break;
      }
      ptrOrganism = ptrOrganism->getNext();
   }

   return seedlingID;
}



//
// JM - 20Feb2007 Method to sample from Gaussian (Normal) distribution using GSL
//        (GNU Scientific Library).
//
double Patch::sampleGaussian(double stdDev)
{
   double sample=0.0;
   sample = gsl_ran_gaussian(Bait::rng, stdDev);

   // Correct for samples outside biologically realistic range.
   while (fabs(sample) > stdDev*3)
      sample = gsl_ran_gaussian(Bait::rng, stdDev);

   return sample;
}


//
// JM 28Jan2014 - Check to see if patch is within Undaria niche (eleveation, substrate type)
//
void Patch::setIsUndariaNiche()
{
   Patch *ptrPatch=NULL;
   isUndariaNiche = false;

   // if correct elevation & substrate, check for "empty" (no substrate) adjacent patch
   // - undaria only grows on "edge" of substrate  
   if (getSubstrateType() > 0 /*&& getIsElevationNiche()*/)
   {
      for (int i=0; i<numDir; i++)
      {
         ptrPatch = getNeighbour((Direction_t)i, true);
         if ( ptrPatch && ptrPatch != this && ptrPatch->getSubstrateType() == 0)
         {
            isUndariaNiche = true;
            break;
         }
      }
   }
}


//***********************************************************************************//
//                                                                                   //
//   Lattice Boltzmann Methods - adapted from lbe.f by S. Succi (Rome, April 2011)   //
//                                                                                   //
//***********************************************************************************//


//
// JM 21APR11 - LBE: initialise Lattice Boltzmann variables (inithydro, S. Succi)
//
/*void Patch::initLBE()
{
   rho_density = baitParamsPtr_->initRho;
   u_vel = baitParamsPtr_->initUvelocity;
     
   // Do the Sine test if no force is specified
   if (baitParamsPtr_->uForce < 1.0E-06)
      u_vel = sin( ((2*M_PI * getIPos()) / engParamsPtr_->height_)*1 );

   v_vel = baitParamsPtr_->initVvelocity;
   f_eq = new double[numDir];
   f = new double[numDir];
   new_f = new double[numDir];
   
   for (int i=0; i<numDir; i++)
   {
      f_eq[i] = 0.0;
      f[i] = 0.0;
      new_f[i] = 0.0;
   }
   
   // initialise LBE values for f_eq and f
   equili();
   initPop();
}


//
// JM 16APR11 - LBE: Initialise frequency values for each direction
//
void Patch::initPop()
{
   for (int i=0; i<numDir; i++)
   {
      f[i] = f_eq[i];
   }
}

//
// JM 18APR11 - LBE: Movement
//
void Patch::moveLBE(EngineParameters *engParamsPtr_, unsigned nbLoops)
{
   Patch *ptrPatch=NULL;
   int rcvDir = 0;
   
   for (int i=0; i<numDir; i++)
      new_f[i] = 0.0;
      
   // Particles with direction NONE don't move
   new_f[0] = f[0];

   for (int i=1; i<numDir; i++)
   {
      ptrPatch = getNeighbourLBE( (Direction_t) i, engParamsPtr_);

      if (i==1 || i==2 || i==5 || i==6)
         rcvDir = i+2;
      else
         rcvDir = i-2;

      if (ptrPatch != this)

         new_f[rcvDir] = ptrPatch->getF(rcvDir);
      // if at boundary of environment, no neighbour, just reverse direction
      else
         new_f[rcvDir] = f[i];
   }

}


//
// JM 03MAY11 - Obstacle:
//
void Patch::obstacle(EngineParameters *engParamsPtr_, unsigned nbLoops)
{
   Patch *ptrPatch=NULL;
   int rcvDir = 0;
      
   // Particles with direction NONE don't move
   //new_f[0] = f[0];
   
      
   for (int i=0; i<numDir; i++)
      new_f[i] = f[i];
   
   if (getSubstrateType() == -1)
   {     
      for (int i=1; i<numDir; i++)
      {
         ptrPatch = getNeighbourLBE( (Direction_t) i, engParamsPtr_);
         
         if (i==1 || i==2 || i==5 || i==6)
            rcvDir = i+2;
         else
            rcvDir = i-2;         
         
         // if is obstacle, direction of flow from neighbour is reversed
         if (ptrPatch != this)
         {                 
            if (i==1 || i==5 || i==8)
               new_f[i] = ptrPatch->getF(rcvDir);
            else if (i == 3)
            {
               new_f[3] = ptrPatch->getF(1);
               new_f[7] = ptrPatch->getF(5);
               new_f[6] = ptrPatch->getF(8);
            }
            
            // Special case for top and bottom of obstacle
            if (isTopObst == true)
            {
               if (i == 2)
               {
                  new_f[2] = ptrPatch->getF(4);
               }
               else if (i == 6)

                  new_f[6] = ptrPatch->getF(8);  
            }
            if (isBottObst == true)
            {
               if (i == 4)
               {
                  new_f[4] = ptrPatch->getF(2);
               }
               else if (i == 7)
                  new_f[7] = ptrPatch->getF(5);  
            }          
         }
         // if at boundary of environment, no neighbour, just reverse direction
         else
         {
            cout << "Obstacle at boundary!!" << endl;
            new_f[rcvDir] = f[i];
         }
      }
   }
}



//
// JM 18APR11 - LBE: Update f levels after movement completed
//
void Patch::updateF(unsigned nbLoops)
{
   for (int i=0; i<numDir; i++)
      f[i] = new_f[i];
      
   //if ( (getIPos() == 13 || getIPos() == 18) && getJPos() == 30 )
      //cout << setprecision(20) << "Upd: " << getIPos() << " " << ( (f[2]+f[5]+f[6]) - (f[4]+f[7]+f[8]) ) << endl; 
}

//
// JM 18APR11 - LBE: Update f levels after movement completed
//
void Patch::updateObst(unsigned nbLoops)
{
   f[1] = new_f[1];
   f[5] = new_f[5];
   f[8] = new_f[8];
   f[3] = new_f[3];
   f[7] = new_f[7];
   f[6] = new_f[6];
   
   if (isTopObst == true)
      f[2] = new_f[2];
   if (isBottObst == true)
      f[4] = new_f[4];
}


//
// JM 16APR11 - LBE: calculate equilibrium values for each direction (9) 
//              using density and velocity parameters
//
void Patch::equili()
{
   double inRho = baitParamsPtr_->initRho;
   double rl = rho_density/inRho;

   double usq = u_vel * u_vel;
   double vsq = v_vel * v_vel;
   double sumsq = (usq + vsq) / cs22;
   double sumsq2 = sumsq * (1.0 - cs2) / cs2;
   
   double u2 = usq / cssq; 
   double v2 = vsq / cssq;
   
   // **Removed use of local variables ui, vi because causes error with compiler (not bug in this code)**
   //double ui = (u_vel / cs2);
   //double vi = (v_vel / cs2);
	
	double uv = (u_vel / cs2) * (v_vel / cs2);
	
	//if ( (getIPos() == 13 || getIPos() == 18) && getJPos() == 30 )
      //cout << setprecision(24) << "Bef: " << getIPos() << " " << " " << vi << endl;
	
	f_eq[0] = rl*inRho*wt0 * (1.0 - sumsq);
	f_eq[1] = rl*inRho*wt1 * (1.0 - sumsq + u2 + (u_vel / cs2) );
   f_eq[2] = rl*inRho*wt1 * (1.0 - sumsq + v2 + (v_vel / cs2) );
   f_eq[3] = rl*inRho*wt1 * (1.0 - sumsq + u2 - (u_vel / cs2) );
   f_eq[4] = rl*inRho*wt1 * (1.0 - sumsq + v2 - (v_vel / cs2) );
   f_eq[5] = rl*inRho*wt2 * (1.0 + sumsq2 + (u_vel / cs2) + (v_vel / cs2) + uv);
   f_eq[6] = rl*inRho*wt2 * (1.0 + sumsq2 - (u_vel / cs2) + (v_vel / cs2) - uv);
   f_eq[7] = rl*inRho*wt2 * (1.0 + sumsq2 - (u_vel / cs2) - (v_vel / cs2) + uv);
   f_eq[8] = rl*inRho*wt2 * (1.0 + sumsq2 + (u_vel / cs2) - (v_vel / cs2) - uv);
   
   //if ( (getIPos() == 11 || getIPos() == 20) && getJPos() == 31 )
   //   cout << endl;
      
   //if ( (getIPos() == 13 || getIPos() == 18) && getJPos() == 30 )
     // cout << setprecision(24) << "Eq: " << getIPos() << " " << ( (f_eq[2]) - (f_eq[4]) ) << endl;

   // Check equilibrium
   double eps = 1.0E-03;
   double znorm = 1.0/rho_density;
   double u_eq = (f_eq[1] - f_eq[3] + f_eq[5] - f_eq[6] - f_eq[7] + f_eq[8])*znorm;
   //double v_eq = (f_eq[5] + f_eq[2] + f_eq[6] - f_eq[7] - f_eq[4] - f_eq[8])*znorm;
   
   // checking u_vel:
   double eu=0.0;
   if (u_vel > 1.0E-05)
   {   
      eu=fabs(u_vel / u_eq - 1.0);
      if (eu > eps)
      {
         cout << "Patch::equili() problem!" << endl;
         cout << "u_vel=" << u_vel << ", ueq=" << u_eq << ", eu=" << eu << endl;
      }    
   }
}

//
// JM 16APR11 - LBE: calculate density (rho) and velocities (u = x axis, v = y axis)
//
void Patch::hydrovar(unsigned nbLoops)
{
   u_vel = 0.0, v_vel = 0.0;
   //total density is the sum of particles frequencies in all 9 directions
   rho_density=0.0;
   for (int i=0; i<numDir; i++)
   {
      rho_density += f[i];
   }
   double rhoi = 1.0/rho_density;
   
   if (isStep != true)
   {
      u_vel = ( (f[1]+f[5]+f[8]) - (f[3]+f[7]+f[6]) )*rhoi ; //diff between east and west movement
      v_vel = ( (f[2]+f[5]+f[6]) - (f[4]+f[7]+f[8]) )*rhoi;  //diff between north and south movement
   }
   else
   {
      u_vel = 0.0;
      v_vel = 0.0;
   }  
   /*if (u_vel > 1.0)
   { 
      cout << "Patch::hydrover(): u_vel = " << u_vel  << " ( " << getIPos() << ", " << getJPos() << ")" << endl;
      cout << (f[1]+f[5]+f[8]) - (f[3]+f[7]+f[6]) << endl;
      cout << f[1] << " " << f[5] << " " << f[8] << " " << f[3] << " " << f[7] << " " << f[6]  << " " << rhoi << endl;
      cin.ignore();
   }*//*
}


//
// JM 16APR - LBE: collision subroutine
//
void Patch::collis(double omega, unsigned nbLoops)
{

   for (int i=0; i<numDir; i++)
   {
      double old_f = f[i];
      f[i] = f[i] * (1.0-omega) + (omega * f_eq[i]);     // 0.0214 * -0.8 + 0.016578
      if (f[i] < 0.0 || f_eq[i] < 0.0)
      {
         cout << "Collis ( " << i << ", " << getIPos() << ", " << getJPos() << "): " << old_f << " " << f[i] << " " << f_eq[i] << " " << omega << endl;
         cin.ignore(); 
      }
   }
   
}


//
// JM 20APR11 - LBE: Apply directional force to flow
//
void Patch::forceLB(double force)
{
   //double force = fpois;
   f[1] += force;
   f[5] += force;
   f[8] += force;
   f[3] -= force;
   f[6] -= force;
   f[7] -= force;
}


//
// JM 16MAY11 - LBE: Set patch as obstacle, reset all values to 0.0
//
void Patch::setIsObstacle(bool inIsObst)
{
   if (inIsObst)
      setSubstrateType(-1);
   else
      setSubstrateType(0);
   
   /*if (isObstacle == true)
   {
      nutrientLevel_ = 0.0;
      newNutrientLevel = 0.0;
   }*/ /*
}


//
// JM 19APR11 - For Testing purpose, print out LBE details
//
void Patch::printLBEDetails(ofstream *outputFile2)
{
   //for (int i=0; i<numDir; i++)     
     // cout << "f[" << i << "] = " << f[i] << ", f_eq[" << i << "] = " << f_eq[i] << endl;
      
   if (outputFile2 != NULL)
   {
      (*outputFile2) << u_vel << " ";
      //(*outputFile2) << "rho = " << rho_density << endl;
   }
   else
      cout << "ERROR" << endl;
   //cout << "u_vel = " << u_vel;
   //cout << "v_vel = " << v_vel;

}


//
// JM 20MAY111 - Algortithm moving molecules based on Lattice Boltzmann velocity field (u_vel and v_Vel)
//
void Patch::fluidFlowLBE(unsigned nbLoops)
{
   Patch *ptrPatchN=NULL, *ptrPatchS=NULL, *ptrPatchE=NULL, *ptrPatchW=NULL;

   newSporeLevel = sporeLevel;

   ptrPatchN = getNeighbourLBE(N, engParamsPtr_);
   ptrPatchS = getNeighbourLBE(S, engParamsPtr_);
   ptrPatchE = getNeighbourLBE(E, engParamsPtr_);
   ptrPatchW = getNeighbourLBE(W, engParamsPtr_);

    
   if (getSubstrateType() > -1)
   {
      // inward movement of fluid molecules
      if (ptrPatchW->getU_vel() > 0.0 )
      {
         if(ptrPatchW)
            newSporeLevel += ptrPatchW->getSporeLevel() * ptrPatchW->getU_vel();       // u*C(x-1,t)
      }
      if (ptrPatchE->getU_vel() < 0.0 )
      {
         if(ptrPatchE)
            newSporeLevel += ptrPatchE->getSporeLevel() * fabs(ptrPatchE->getU_vel()); // -u*C(x+1,t)
      }
      if (ptrPatchS->getV_vel() > 0.0)
      {
         if(ptrPatchS)
            newSporeLevel += ptrPatchS->getSporeLevel() * ptrPatchS->getV_vel();       // v*C(y-1,t)
      }
      if (ptrPatchN->getV_vel() < 0.0)
      {
         if(ptrPatchN)
            newSporeLevel += ptrPatchN->getSporeLevel() * fabs(ptrPatchN->getV_vel()); // -v*C(y+1,t)
      }
   
      // outward movement of fluid molecules
      if (u_vel > 0.0 && ptrPatchE->getSubstrateType() > -1)
         newSporeLevel -= sporeLevel * u_vel;
      else if (u_vel < 0.0 && ptrPatchW->getSubstrateType() > -1)
         newSporeLevel -= sporeLevel * fabs(u_vel);
         
      if (v_vel > 0.0 && ptrPatchN->getSubstrateType() > -1)
         newSporeLevel -= sporeLevel * v_vel;
      else if (v_vel < 0.0 && ptrPatchS->getSubstrateType() > -1)
         newSporeLevel -= sporeLevel * fabs(v_vel);
         
      if (u_vel + v_vel > 1.0)
         cout << "Patch::fluidFlowLBE() Problem: " << "u = " << u_vel << " v = " << v_vel << endl;
   }
}*/





// ------------------------------ End Of File ------------------------------- //
