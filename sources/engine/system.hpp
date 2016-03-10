// -------------------------------------------------------------------------- //
// system.hpp                                                                 //
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

#ifndef SYSTEM_HPP
#define SYSTEM_HPP


// -------------------------------------------------------------------------- //
// Librairies                                                                 //
// -------------------------------------------------------------------------- //

#include "world.hpp"


// -------------------------------------------------------------------------- //
// Forward declarations                                                       //
// -------------------------------------------------------------------------- //

class EngineParameters ;


// -------------------------------------------------------------------------- //
// class System                                                               //
// -------------------------------------------------------------------------- //

class System
{
   public :

      System() ;
      virtual ~System() ;

     // JM - intialize engine parameters
      static void initEngParamsPtr( const EngineParameters * ) ;

     // JM - simulation main loop
      template <typename _CellType_>
      void run( World<_CellType_> & ) ;

     // virtual functions implemented in class Bait
      virtual void displayTime( unsigned ) {}
      virtual void latticeBoltzmannTime(unsigned){}
      virtual void worldEvolutionTime(unsigned)    {}
      virtual void staticActivityTime(unsigned)=0;
      virtual void startMoveTime() {}
      virtual void moveTime()              {}
      virtual void solveConflictsTime(unsigned)    {}
      virtual void interactionsTime(unsigned)      {}
      virtual void reproductionTime(unsigned)      {}
      
      static const EngineParameters* getEngParamsPtr()
         { return engParamsPtr_; }


   private :

      static const EngineParameters*    engParamsPtr_ ;
};


// -------------------------------------------------------------------------- //
// void System::run()                                                         //
//                                                                            //
// Simulation main loop                                                       //
// -------------------------------------------------------------------------- //

template <typename _CellType_>
void System::run( World<_CellType_> &inWorld )
{
   unsigned   displayCounter = 0;
   //int stopRun = 0;

   // Main loop
   // 'nbLoops' limits how many loops the program goes through
   for( unsigned nbLoops=1 ; ! engParamsPtr_->nbLoops_ ||
      nbLoops<=engParamsPtr_->nbLoops_ ; ++nbLoops )
   {
      // 1 - Display the world
      // displays world every 'displayPeriod' number of loops
      if( engParamsPtr_->displayPeriod_ && !displayCounter )
      {
         displayTime( nbLoops ) ;

         displayCounter = engParamsPtr_->displayPeriod_ ;
      }
      --displayCounter ;
      
      // 2 - JM 20APR11  Lattice Boltzmann Methods (LBE)
      // for calculating hydrodynamic properties (velocities) 
      latticeBoltzmannTime(nbLoops) ;

      // 3 - World evolution
      // for evolving the characteristics of the world (diffusion)
      worldEvolutionTime(nbLoops) ;

      // 4 - Static activity of agents
      // subtract survival cost, display number of bacteria
      staticActivityTime(nbLoops);
      
      // 5 - Reproduction time
      reproductionTime(nbLoops) ;

      // Start the move time (see cell.hpp and cell.cpp)
      // Deletes agent list from each cell before new movement time 
      //     which will rebuild it
      startMoveTime() ;

      // 6 - Agents movment
      // moves bacteria and antibiotic agents
      moveTime() ;

      // 7 - Conflicts resolutions
      // e.g. antibiotic and bacteria in same patch/cell (not used)
      solveConflictsTime(nbLoops) ;

      // 8 - Lunch time
      // nbLoops passed because feeding during lag phase
      interactionsTime(nbLoops) ;
   
      //if (stopRun == 1)
         //break;
   }
}


// -------------------------------------------------------------------------- //

#endif                        // SYSTEM_HPP


// ------------------------------ End Of File ------------------------------- //
