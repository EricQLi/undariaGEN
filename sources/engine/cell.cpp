// -------------------------------------------------------------------------- //
// cell.cpp                                                                   //
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

#include "cell.hpp"


// -------------------------------------------------------------------------- //
// Libraries                                                                  //
// -------------------------------------------------------------------------- //

#include "agent.hpp"
#include "myException.hpp"
#include <stdlib.h>
#include <iostream>
using std::cout;
using std::endl;


// -------------------------------------------------------------------------- //
// Constructor                                                                //
// -------------------------------------------------------------------------- //

Cell::Cell()
{
   agents_    = NULL ;
}


// -------------------------------------------------------------------------- //
// Destructor                                                                 //
// -------------------------------------------------------------------------- //

Cell::~Cell()
{
}


// -------------------------------------------------------------------------- //
// Insert a agent to the agents list                                          //
// -------------------------------------------------------------------------- //

void Cell::addAgent( Agent *inPtrAgent )
{
   // Add the agent to the head of the list
   inPtrAgent->setNext( agents_ ) ;
   agents_ = inPtrAgent ;

   // Make the agent refer to this cell
   inPtrAgent->setCell( this ) ;
}


// -------------------------------------------------------------------------- //
// Delete agent list before new movement time which will rebuild it.          //
// -------------------------------------------------------------------------- //

void Cell::resetAgentsPtr()
{
   // Initialize agent list
   if (agents_)
      agents_ = NULL ;
}


// ------------------------------ End Of File ------------------------------- //
