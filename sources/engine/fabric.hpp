// -------------------------------------------------------------------------- //
// frabric.hpp                                                                //
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

#ifndef FABRIC_HPP
#define FABRIC_HPP


// -------------------------------------------------------------------------- //
// Libraries                                                                  //
// -------------------------------------------------------------------------- //

#include "myException.hpp"
#include "fabricAliveIterator.hpp"
#include <stdlib.h>


// -------------------------------------------------------------------------- //
// class Fabric                                                               //
// -------------------------------------------------------------------------- //

template <typename _AgentType_>
class Fabric
{
   public :

      Fabric() ;
      virtual ~Fabric() ;

      void init( unsigned ) ;
      void freeMemory() ;

      void setHighestIndexAgent();
      void setHighestIndexAgent(unsigned index)
      	{ highestIndexAgent = index; }     
      void defragmentFabric();
      void addToAgentBuf(double, float [], short [], int);            	
      void debugFabric();

      _AgentType_ *newAgent() ;

      FabricAliveIterator<_AgentType_> begin() ;

      unsigned aliveAgentsAction( void(_AgentType_::*)(unsigned), unsigned );
      unsigned aliveAgentsAction( void(_AgentType_::*)() );
      void aliveAgentsAction2( void(_AgentType_::*)(unsigned, float, float), unsigned, float, float );
      void aliveAgentsAction2( void(_AgentType_::*)(unsigned, double, unsigned, float, float, int), unsigned, double, unsigned, float, float, int );

      bool getFirstRemoved() const
         { return firstRemoved; }
      unsigned getHighestIndexAgent() const
      	{ return highestIndexAgent; }
      	
      void setFirstRemoved(bool removed)
         { firstRemoved = removed; }

   private :

      unsigned         size_    ,
                       last_    ,
                       highestIndexAgent;   // Position of highest alive agent.
      _AgentType_      *storage_ ;
      bool firstRemoved;                  // First agent removed (e.g. death)
};


// -------------------------------------------------------------------------- //
// Constructor                                                                //
// -------------------------------------------------------------------------- //

template <typename _AgentType_>
Fabric<_AgentType_>::
Fabric()
{
   storage_ = NULL ;

   // Start allocation from the beginning of the array
   last_ = 0 ;
}


// -------------------------------------------------------------------------- //
// Destructor                                                                 //
// -------------------------------------------------------------------------- //

template <typename _AgentType_>
Fabric<_AgentType_>::
~Fabric()
{
   freeMemory() ;
}


// -------------------------------------------------------------------------- //
// Allocate the storage array                                                 //
// -------------------------------------------------------------------------- //

template <typename _AgentType_>
void
Fabric<_AgentType_>::
init( unsigned inSize )
{
   // Save size
   size_ = inSize ;
   firstRemoved = false;
   highestIndexAgent = 0;

   // Allocate the storage array
   storage_ = new _AgentType_ [ size_ ] ;
   if( !storage_ )
      throw myException( "Cannot allocate fabric", __FILE__, __LINE__ ) ;
}


// -------------------------------------------------------------------------- //
// Free memory used by storage array                                          //
// -------------------------------------------------------------------------- //

template <typename _AgentType_>
void
Fabric<_AgentType_>::
freeMemory()
{
   if( storage_ )
      delete [] storage_ ;

   storage_ = NULL ;
}


//
// JM 23Mar2007 - Look for alive agent with highest index in fabric.
//                - Only iterate up to this index (Performance improvement).
//
template <typename _AgentType_>
void
Fabric<_AgentType_>::
setHighestIndexAgent()
{
   // DEBUG : Avoid SIGSEV with unallocated fabric
   if( !storage_ )
      throw myException( "Unallocate fabric", __FILE__, __LINE__ ) ;

   highestIndexAgent = 0;

   // Look for the highest index agent
   for( unsigned i=0; i<size_; ++i )
   {
      if( storage_[i].isAlive() )
         highestIndexAgent = i;
   }
}


// -------------------------------------------------------------------------- //
// Assign new agent                                                           //
// -------------------------------------------------------------------------- //

template <typename _AgentType_>
_AgentType_ *
Fabric<_AgentType_>::
newAgent()
{
   // DEBUG : Avoid SIGSEV with unallocated fabric
   if( !storage_ )
      throw myException( "Unallocate fabric", __FILE__, __LINE__ ) ;

   // Look for a dead agent from the last to the top
   for( unsigned i=highestIndexAgent; i<size_; ++i )
   {
      if( ! storage_[i].isAlive() )
      {
         highestIndexAgent = i ;
         return storage_+i ;
      }
   }

   // Look for a dead agent from the bottom to the last
   for( unsigned i=0; i<highestIndexAgent; ++i )
   {
      if( ! storage_[i].isAlive() )
         return storage_+i ;
   }

   // All agents are alive: out of memory
   throw myException( "Fabric is out of memory", __FILE__, __LINE__ ) ;

   //return NULL ;
}


// -------------------------------------------------------------------------- //
// Return an iterator on the first alive agent                                //
// JM - finds the first alive agent in the list
// -------------------------------------------------------------------------- //

template <typename _AgentType_>
FabricAliveIterator<_AgentType_>
Fabric<_AgentType_>::
begin()
{
   // Look for the first alive agent
   // Only search up to highestIndexAgent (instead of size_)
   // Fixes slowdown due to searching entire fabric every loop when no agents alive (MPI)
   unsigned searchLimit = 10;     // search at least first 10 spaces
   if (highestIndexAgent >= searchLimit)
      searchLimit = highestIndexAgent;
      
   for( unsigned i=0; i<searchLimit+2; ++i )
   {
      // Only need to iterate through fabric up to highestIndexAgent.
      if( storage_[i].isAlive() )
         return FabricAliveIterator<_AgentType_>( storage_, (highestIndexAgent+1), i ) ;
   }

   // No alive agent found - place iterator at end of fabric
   return FabricAliveIterator<_AgentType_>( storage_, size_, size_ ) ;
}


// -------------------------------------------------------------------------- //
// Call given method for each alive agent.                                    //
//                                                                            //
// Return the number of alive agent                                           //
// JM - iterates through list of agents, counts the alive ones           //
//      and makes them perform some action (*method) e.g. move/eat           //
// -------------------------------------------------------------------------- //

template< typename _AgentType_ >
unsigned
Fabric<_AgentType_>::aliveAgentsAction
( void (_AgentType_::*method) (unsigned), unsigned nbLoops )
{
   unsigned nbAliveAgents = 0 ;
   FabricAliveIterator<_AgentType_> it = begin() ;

   // For each alive agent in factory
   while( *it )
   {
      // Move it
      ((*it)->*method)(nbLoops)  ;

      // Count it
      ++nbAliveAgents ;

      // Go to next alive agent
      ++it ;
   }

   return nbAliveAgents ;
}


template< typename _AgentType_ >
unsigned
Fabric<_AgentType_>::aliveAgentsAction
( void (_AgentType_::*method) () )
{
   unsigned nbAliveAgents = 0 ;
   FabricAliveIterator<_AgentType_> it = begin() ;

   // For each alive agent in factory
   while( *it )
   {
      // Move it
      ((*it)->*method)()  ;

      // Count it
      ++nbAliveAgents ;

      // Go to next alive agent
      ++it ;
   }

   return nbAliveAgents ;
}


template< typename _AgentType_ >
void
Fabric<_AgentType_>::aliveAgentsAction2
( void (_AgentType_::*method) (unsigned, float, float), unsigned param1, float param2, float param3 )
{
   FabricAliveIterator<_AgentType_> it = begin() ;

   // For each alive agent in factory
   while( *it )
   {
      // Move it
      ((*it)->*method)(param1, param2, param3)  ;

      // Go to next alive agent
      ++it ;
   }
}


//
// JM 01Nov12 - for calling Organism::releasePollenMPI()
//
template< typename _AgentType_ >
void
Fabric<_AgentType_>::aliveAgentsAction2
( void (_AgentType_::*method) (unsigned, double, unsigned, float, float, int), unsigned param1, double param2, unsigned param3, float param4, float param5, int param6)
{
   FabricAliveIterator<_AgentType_> it = begin() ;

   // For each alive agent in factory
   while( *it )
   {
      // Move it
      ((*it)->*method)(param1, param2, param3, param4, param5, param6)  ;

      // Go to next alive agent
      ++it ;
   }
}


//
// JM 05Dec2007
//
template< typename _AgentType_ >
void
Fabric<_AgentType_>::defragmentFabric()
{
   // DEBUG : Avoid SIGSEV with unallocated fabric
   if( !storage_ )
      throw myException( "Unallocate fabric", __FILE__, __LINE__ ) ;

   if (firstRemoved == true)
   {
      // Look for a dead agent within the block of alive agents.     
      unsigned i=0; //firstDead;
      while (i < highestIndexAgent)
      {
         //cout << i << " ";
         if( !storage_[i].isAlive() )
         {
            while (!(storage_[highestIndexAgent].isAlive()) && highestIndexAgent > 0)
            {
               highestIndexAgent--;
            }

            if ( storage_[highestIndexAgent].isAlive()  && i < highestIndexAgent)
            {
               // Use overloaded assignment operator to assign object.
               storage_[i] = storage_[highestIndexAgent];
               storage_[highestIndexAgent].setStock(0.0);
               highestIndexAgent--;
            }
         }
         i++;
      }
   }
}


//
// 12Mar2007 - Add parameters for last agent in fabric to temporary buffers.
//
template< typename _AgentType_ >
void
Fabric<_AgentType_>::addToAgentBuf(double stockBuf, float floatsBuf [], short boolsBuf [], int integerBuf)
{
   // Input into buffers
   // Get agent's stock_ level.
   stockBuf = storage_[highestIndexAgent].getStock();

   // Get agent's floating point traits.
   floatsBuf[0] = storage_[highestIndexAgent].getXInCell();
   floatsBuf[1] = storage_[highestIndexAgent].getYInCell();
   floatsBuf[2] = storage_[highestIndexAgent].getXDelta();
   floatsBuf[3] = storage_[highestIndexAgent].getYDelta();
   floatsBuf[4] = storage_[highestIndexAgent].getStationaryLevel();
   floatsBuf[5] = storage_[highestIndexAgent].getPbp2aAcylated();

   // Get agent's boolean traits.
   boolsBuf[0] = (short)storage_[highestIndexAgent].getMutator();
   boolsBuf[1] = (short)storage_[highestIndexAgent].getAbResistant();
   boolsBuf[2] = (short)storage_[highestIndexAgent].getExpressLac();
   boolsBuf[3] = (short)storage_[highestIndexAgent].getSOSResponse();
   boolsBuf[4] = (short)storage_[highestIndexAgent].getInitAbDeath();
         
   // Get agent's integer traits.
   integerBuf = storage_[highestIndexAgent].getAbTimer();
}


template< typename _AgentType_ >
void
Fabric<_AgentType_>::debugFabric()
{
   cout << "highestIndexAgent = " << highestIndexAgent << endl;
}

// -------------------------------------------------------------------------- //

#endif                        // FABRIC_HPP


// ------------------------------ End Of File ------------------------------- //
