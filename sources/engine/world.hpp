// -------------------------------------------------------------------------- //
// world.hpp                                                                  //
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

#ifndef WORLD_HPP
#define WORLD_HPP


// -------------------------------------------------------------------------- //
// Libraries                                                                  //
// -------------------------------------------------------------------------- //

#include "agent.hpp"
#include "fabric.hpp"
#include "fabricAliveIterator.hpp"
#include "engineParameters.hpp"
#include "../model/baitParameters.hpp"
#include <time.h>


// -------------------------------------------------------------------------- //
// Forward declarations                                                       //
// -------------------------------------------------------------------------- //

class Cell ;


// -------------------------------------------------------------------------- //
// class World                                                                //
// -------------------------------------------------------------------------- //

template< typename _CellType_ >
class World
{
   public :

      World() ;
      virtual ~World() ;

      static void initEngParamsPtr( EngineParameters * ) ;

      void init(int) ;
      void freeMemory() ;
      
      _CellType_& operator() ( int, int ) const ;
      _CellType_& smartAccess( int, int ) const ;

      void startMoveTime() ;

      void allCellsDevelop(unsigned) ;
      void allCellsDiffuse(unsigned) ;

      static EngineParameters* engParamsPtr_ ;
      static BaitParameters* baitParamsPtr_ ;

   protected :
      _CellType_    *cells_     ;
};


// -------------------------------------------------------------------------- //
// Post declaration librairies                                                //
// -------------------------------------------------------------------------- //

#include "cell.hpp"


// -------------------------------------------------------------------------- //
// Static member creation                                                     //
// -------------------------------------------------------------------------- //

template< typename _CellType_ >
EngineParameters* World<_CellType_>::engParamsPtr_ = 0 ;


// -------------------------------------------------------------------------- //
// Constructor                                                                //
// -------------------------------------------------------------------------- //

template< typename _CellType_ >
World<_CellType_>::World
()
{
   cells_ = NULL ;
}


// -------------------------------------------------------------------------- //
// Destructor                                                                 //
// -------------------------------------------------------------------------- //

template< typename _CellType_ >
World<_CellType_>::~World
()
{
}


// -------------------------------------------------------------------------- //
// Initialiaze the static parameters pointer                                  //
// -------------------------------------------------------------------------- //

template< typename _CellType_ >
void
World<_CellType_>::initEngParamsPtr
( EngineParameters *inPtr )
{
   engParamsPtr_ = inPtr ;
}


// -------------------------------------------------------------------------- //
// All components initialisation                                              //
// -------------------------------------------------------------------------- //

template< typename _CellType_ >
void
World<_CellType_>::init
(int mpiSize)
{
   // Allocations
   cells_ = new _CellType_[ engParamsPtr_->width_ * engParamsPtr_->height_ ] ;
   if( !cells_ )
      throw myException( "Cannot allocate moving grid", __FILE__, __LINE__ ) ;

   // Outloop (vertical)
   for( int i=0; i<engParamsPtr_->height_; ++i )
   {
      // Inloop (horizontal)
      for( int j=0; j<engParamsPtr_->width_; ++j )
      {
         // Coordinates initialisation
         (*this)(i,j).setCoordinates( i, j );
         // Cells initialisation
         (*this)(i,j).init(this, mpiSize);

         
      }
   }
}


// -------------------------------------------------------------------------- //
// Free memory used by storage array                                          //
// -------------------------------------------------------------------------- //

template< typename _CellType_ >
void
World<_CellType_>::freeMemory
()
{
   if( cells_ )
      delete [] cells_ ;

   cells_ = NULL ;
}


//----------------------------------------------------------------------------//
// Access method to the cells                                                 //
//----------------------------------------------------------------------------//

template< typename _CellType_ >
_CellType_&
World<_CellType_>::operator()
( int iPos, int jPos )
const
{
   /// DEBUG : Avoid out of range
   if( iPos >= engParamsPtr_->height_ )
   {
      cout << "iPos = " << iPos << ", jPos = " << jPos<< endl;
      throw myException( "iPos out of range (movingGrids)", __FILE__, __LINE__ ) ;
   }
   else if (jPos >= engParamsPtr_->width_)
   {
      cout << "iPos = " << iPos << ", jPos = " << jPos<< endl;
      throw myException( "jPos out of range (movingGrids)", __FILE__, __LINE__ ) ;
   }

   return cells_[ iPos * engParamsPtr_->width_ + jPos ] ;
}


//----------------------------------------------------------------------------//
// Access method to the cells but with possible extra range                   //
//----------------------------------------------------------------------------//

template< typename _CellType_ >
_CellType_&
World<_CellType_>::smartAccess
( int iPos, int jPos )
const
{
   // Determine iPos value
   if( iPos < 0 )
   {
       if( engParamsPtr_->worldType_ == CLOSED || engParamsPtr_->worldType_ == OPEN_EW )
          iPos = 0 ;
       else
          iPos = engParamsPtr_->height_ + iPos;
   }
   else if( iPos >= engParamsPtr_->height_ )
   {
       if( engParamsPtr_->worldType_ == CLOSED || engParamsPtr_->worldType_ == OPEN_EW )
          iPos = engParamsPtr_->height_ - 1 ;
       else
          iPos = iPos % engParamsPtr_->height_ ;
   }

   // Determine jPos value
   if( jPos < 0 )
   {
       if( engParamsPtr_->worldType_ == CLOSED || engParamsPtr_->worldType_ == OPEN_NS )
          jPos = 0 ;
       else
          jPos = engParamsPtr_->width_ + jPos;
   }
   else if( jPos >= engParamsPtr_->width_ )
   {
       if( engParamsPtr_->worldType_ == CLOSED || engParamsPtr_->worldType_ == OPEN_NS )
          jPos = engParamsPtr_->width_ - 1 ;
       else
          jPos = jPos % engParamsPtr_->width_ ;
   }

   // Return the correct cell
   return cells_[ iPos * engParamsPtr_->width_ + jPos ] ;
}


// -------------------------------------------------------------------------- //
// Delete all cells agents lists                                              //
// -------------------------------------------------------------------------- //

template< typename _CellType_ >
void
World<_CellType_>::startMoveTime
()
{        
   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;
   for( unsigned i=0; i<nbCells; i++ )
      this->cells_[i].resetAgentsPtr();
}


// -------------------------------------------------------------------------- //
// Each cell develop itself                                                   //
// -------------------------------------------------------------------------- //

template< typename _CellType_ >
void
World<_CellType_>::allCellsDevelop
(unsigned nbLoops)
{
   unsigned nbCells = engParamsPtr_->width_ * engParamsPtr_->height_ ;
   for( unsigned i=0; i<nbCells; i++ )
      this->cells_[i].development(nbLoops);
}


// -------------------------------------------------------------------------- //
// Each cell diffuse info to other cells                                      //
// -------------------------------------------------------------------------- //

template< typename _CellType_ >
void
World<_CellType_>::allCellsDiffuse
(unsigned nbLoops)
{
   //[Empty]
}


// -------------------------------------------------------------------------- //

#endif                        // WORLD_HPP


// ------------------------------ End Of File ------------------------------- //
