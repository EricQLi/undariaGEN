// -------------------------------------------------------------------------- //
// system.cpp                                                                 //
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

#include "system.hpp"


// -------------------------------------------------------------------------- //
// Libraries                                                                  //
// -------------------------------------------------------------------------- //

#include "engineParameters.hpp"


// -------------------------------------------------------------------------- //
// Static member creation                                                     //
// -------------------------------------------------------------------------- //

const EngineParameters* System::engParamsPtr_ = 0 ;


// -------------------------------------------------------------------------- //
// Constructor                                                                //
// -------------------------------------------------------------------------- //

System::System()
{
}


// -------------------------------------------------------------------------- //
// Destructor                                                                 //
// -------------------------------------------------------------------------- //

System::~System()
{
}


// -------------------------------------------------------------------------- //
// Initialiaze the static parameters pointer                                  //
// -------------------------------------------------------------------------- //

void System::initEngParamsPtr( const EngineParameters *inPtr )
{
   engParamsPtr_ = inPtr ;
}


// ------------------------------ End Of File ------------------------------- //
