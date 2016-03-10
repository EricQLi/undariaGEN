// -------------------------------------------------------------------------- //
// engineParameters.hpp                                                       //
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

#ifndef ENGINEPARAMETERS_HPP
#define ENGINEPARAMETERS_HPP


// -------------------------------------------------------------------------- //
// Enumeration                                                                //
// -------------------------------------------------------------------------- //

enum WorldType_t { OPEN, CLOSED, OPEN_EW, OPEN_NS, OPEN_W, OPEN_E } ;


// -------------------------------------------------------------------------- //
// class EngineParameters                                                     //
// -------------------------------------------------------------------------- //

class EngineParameters
{
   public :

      EngineParameters() :
      width_             (  256 ),
      height_            (  256 ),
      nbLoops_           (    0 ),
      displayPeriod_     (    0 ),
      worldType_         ( OPEN )
      {
      } ;


   public :

      unsigned       width_         ,
                     height_        ,
                     nbLoops_       ,
                     displayPeriod_ ,
                     graphPeriod,
                     sampleLoop,
							nbRuns;
		float          cellSize;
		bool           mobility;
      WorldType_t    worldType_     ;
} ;


// -------------------------------------------------------------------------- //

#endif                        //  ENGINEPARAMETERS_HPP


// ------------------------------ End Of File ------------------------------- //
