// -------------------------------------------------------------------------- //
// cell.hpp                                                                   //
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

#ifndef CELL_HPP
#define CELL_HPP


// -------------------------------------------------------------------------- //
// Libraries                                                                  //
// -------------------------------------------------------------------------- //

#include "fabric.hpp"
#include "world.hpp"


// -------------------------------------------------------------------------- //
// Enumeration                                                                //
// -------------------------------------------------------------------------- //

enum Direction_t { NONE, E, N, W, S, NE, NW, SW, SE };


// -------------------------------------------------------------------------- //
// Forward declarations                                                       //
// -------------------------------------------------------------------------- //

class Agent ;


// -------------------------------------------------------------------------- //
// class Cell                                                                 //
// -------------------------------------------------------------------------- //

class Cell
{
   public :

      Cell() ;
      virtual ~Cell() ;

      virtual void addAgent( Agent * ) ;
      void resetAgentsPtr(); 

      virtual void development(unsigned) {}
      template< typename _CellType_ >
      void diffusion( World<_CellType_> &, int, int ) {}

      void setCoordinates( unsigned inIPos, unsigned inJPos )
         { iPos_=inIPos; jPos_=inJPos; }
      unsigned short getIPos() const { return iPos_ ; }
      unsigned short getJPos() const { return jPos_ ; }

      virtual Cell* getNeighbour(Direction_t, bool)=0;
      Agent* getAgents()
         { return agents_; }


   private :
   
      Agent        *agents_    ;
      unsigned short iPos_, 
                     jPos_;
};


// -------------------------------------------------------------------------- //

#endif                        // CELL_HPP


// ------------------------------ End Of File ------------------------------- //
