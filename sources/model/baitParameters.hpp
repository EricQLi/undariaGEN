// -------------------------------------------------------------------------- //
// baitParameters.hpp                                                         //
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

#ifndef BAITPARAMETERS_HPP
#define BAITPARAMETERS_HPP

#include <string>
using std::string;

// -------------------------------------------------------------------------- //
// Enumeration                                                                //
// -------------------------------------------------------------------------- //

enum AgentType_t { ORGANISM, ANTIBIO } ;


// -------------------------------------------------------------------------- //
// class BaitParameters                                                       //
// -------------------------------------------------------------------------- //

class BaitParameters
{
   public :
   
      int        rankAddObst,
                 loopAddObst,
                 loopRemoveObst;

      unsigned   initialOrganisms_        ,
                 maxOrganisms_           ,
                 minMolecules_, maxMolecules_,
                 sndBufferSize,               
                 // Lattic Boltzmann:
                 obstXPos, obstYPos,
                 obstLength, obstWidth,
                 initDistribute,
                 CUDALoops,
                 diffusLoops;


      double     maxTotalNutrientLevel_  ,
                 diffusionRate,
                 probFertilise,

                 sporeProd,
                 kCatLactamase[2],
                 kMLactamase[2],
                 pbp2a_k2[2],
                 pbp2a_Kd[2],
                 
                 // Lattice Boltzmann method
                 initRho,
                 initUvelocity, initVvelocity,
                 omega,
                 uForce,
                 ratio_hH;

      float      sporeHalfLife,
                 probGerminate,
                 totalSpore,
                 macroSize,
                 startStock,
                 kdPAR,                // irradiance attenuation coefficient
                 meanDepth,
                 photo_alpha,
                 
                 // temperature, solar & day length input parameters (Cosine function)
                 cos_units,
                 tempAmp,
                 temp_c,
                 temp_d,
                 solarAmp,
                 solar_c,
                 solar_d,
                 DL_Amp,
                 DL_c,
                 DL_d;

      bool       displayOrganismsNumber_,
                 isAmp[2];

} ;


// -------------------------------------------------------------------------- //

#endif                        //  BAITPARAMETERS_HPP


// ------------------------------ End Of File ------------------------------- //
