// -------------------------------------------------------------------------- //
// agent.hpp                                                                  //
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

#ifndef AGENT_HPP
#define AGENT_HPP


// -------------------------------------------------------------------------- //
// Libraries                                                                  //
// -------------------------------------------------------------------------- //

#include "fabric.hpp"


// -------------------------------------------------------------------------- //
// Forward declarations                                                       //
// -------------------------------------------------------------------------- //

class Cell ;


// -------------------------------------------------------------------------- //
// class Agent                                                                //
// -------------------------------------------------------------------------- //

class Agent
{
   public :

      Agent() ;
      virtual ~Agent() ;

      virtual bool isAlive() = 0 ;
      virtual void init(int,int)    = 0 ;
      virtual void move()    = 0 ;
      virtual void antibioDamage(double){};   // JM 09Jan2006

      Agent* getNext() const { return next_ ; }
      void setNext( Agent *newPtr ) { next_ = newPtr ; }
      Agent** getNextAddr() { return &next_ ; }

      Cell* getCell() const { return cell_ ; }
      void setCell( Cell *newPtr ) { cell_ = newPtr ; }

      virtual unsigned getType() const = 0 ;


   private :

      Agent   *next_ ;
      Cell    *cell_ ;
};


// -------------------------------------------------------------------------- //

#endif                        // AGENT_HPP


// ------------------------------ End Of File ------------------------------- //
