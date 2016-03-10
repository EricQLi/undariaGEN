// -------------------------------------------------------------------------- //
// fabricAliveIterator.hpp                                                    //
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

#ifndef FABRICALIVEITERATOR_HPP
#define FABRICALIVEITERATOR_HPP


// -------------------------------------------------------------------------- //
// Libraries                                                                  //
// -------------------------------------------------------------------------- //

#include <stdlib.h>
#include <iostream>
using std::cout;
using std::endl;


// -------------------------------------------------------------------------- //
// class FabricAliveIterator                                                  //
// -------------------------------------------------------------------------- //

template <typename _AgentType_>
class FabricAliveIterator
{
   public :

      FabricAliveIterator( _AgentType_*, unsigned, unsigned ) ;

      _AgentType_* operator* () ;

      FabricAliveIterator<_AgentType_>& operator++ () ;

   protected :

      _AgentType_   *ptrBase_  ;   // pointer to first object in list
      unsigned       size_     ,   // size of list
                     current_  ;   // current position in list
};


// -------------------------------------------------------------------------- //
// Constructor                                                                //                                //
// -------------------------------------------------------------------------- //

template <typename _AgentType_>
FabricAliveIterator<_AgentType_>::
FabricAliveIterator( _AgentType_ *inPtrBase, unsigned inSize, unsigned inPos )
: ptrBase_( inPtrBase ), size_( inSize ), current_( inPos )
{
}


// -------------------------------------------------------------------------- //
// Access operator                                                            //
// checks if current pointer is within bounds of list, else returns null      //
// -------------------------------------------------------------------------- //

template <typename _AgentType_>
inline
_AgentType_*
FabricAliveIterator<_AgentType_>::
operator* ()
{
   if( current_ < size_ )
      return ptrBase_ + current_ ;
   
   return NULL ;
}


// -------------------------------------------------------------------------- //
// Pre-incrementation operator                                                //
// used for moving to next alive agent in list                                //
// -------------------------------------------------------------------------- //

template <typename _AgentType_>
inline
FabricAliveIterator<_AgentType_>&
FabricAliveIterator<_AgentType_>::
operator++ ()
{
   // Looking for the next alive agent
   while( ++current_ < size_ )
   {
      if( ptrBase_[ current_ ].isAlive() )
         break ;
   }

   return *this ;
}


// -------------------------------------------------------------------------- //

#endif                        // FABRICALIVEITERATOR_HPP


// ------------------------------ End Of File ------------------------------- //
