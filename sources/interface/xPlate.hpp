// -------------------------------------------------------------------------- //
// xPlate.hpp                                                                 //
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

#ifndef XPLATE_HPP
#define XPLATE_HPP


// -------------------------------------------------------------------------- //
// Libraries                                                                  //
// -------------------------------------------------------------------------- //

#include "../model/plate.hpp"
#include "../engine/fabric.hpp"
#include <X11/Xlib.h>
#include <sys/time.h>


// -------------------------------------------------------------------------- //
// Forward declarations                                                       //
// -------------------------------------------------------------------------- //

class Organism ;
class Antibio ;
class InterfaceParameters ;


// -------------------------------------------------------------------------- //
// class XPlate                                                               //
// -------------------------------------------------------------------------- //

class XPlate : public Plate
{
   public :

      XPlate();
      virtual ~XPlate();

      static void initXParamsPtr( const InterfaceParameters * ) ;
      void init(int, int, int *, int *, int *, float *, float *, float *) ;

      void display( unsigned, Fabric<Organism>&, int, int) ;


   protected :

      void initColors() ;
      void keyHandler( XKeyEvent* ) ;
      void checkWindowEvent() ;
      //void displayPatches(unsigned) ;
      //void displayPatchesWin3();
      //void displayAntibios(unsigned, int);
      //void displayActAntibios(unsigned);
      void displaySpore(unsigned);
      template< typename _AgentType_ >
      void displayAgents( Fabric<_AgentType_>& , unsigned long ) ;
      void displayOrganisms(unsigned);
      //void displayPollen(unsigned);
      //void displayPollenMPI(unsigned);
      
      void displayAllPatchesMPI(unsigned);
      void displayAllBiomassMPI(unsigned);
      void displayAllElevationsMPI(unsigned);
      void displayElevations(unsigned);
      void displayWaterDepth(unsigned);
      void displayWaterTemp(unsigned);
      void displayWaterTempMPI(unsigned);
      void displaySubstrateTypes(unsigned);
      void displaySubstrateTypesMPI(unsigned);
      void displaySporeMPI(unsigned);


   protected :

      static const InterfaceParameters* xParamsPtr_ ;

      Window          win[3];
      Window          win_              ,
                      win2,
                      win3;
      Display        *dpy_              ;
      GC                gcontext[3];
      GC              gcontext_         ,
                      gcontext2,
                      gcontext3;
      XGCValues       gcv_              ;
      Colormap        cmap_             ;
      Pixmap          dbBuf[3];
      Pixmap          dbBuf_,
                      dbBuf2,
                      dbBuf3;
      unsigned long   gcmask_           ,
                     *patchesColors1_,
                     *patchesColors2_,
                     *patchesMedColors_,
                     *patchesLowColors_,
                     //*antibiosColours,
                     //*actAntibiosColours,
                     *sporeColours,
                     *bacteriasColours_  ,
                     *elevColours,
                     //*posChargeColours,
                     *namedColors_      ;
      unsigned        windowWidth_      ,
                      windowHeight_     ;
      double          xCelSize_         ,
                      yCelSize_         ;
      unsigned        lastLoopNumber_   ;
      struct timeval  lastTime_         ;
      bool            pause_            ;
      const unsigned short numWindows;
};


// -------------------------------------------------------------------------- //
// void XPlate::displayAgents( Fabric<_AgentType_> &, unsigned long )         //
//                                                                            //
// Display all alive agents in the specified color                            //
// -------------------------------------------------------------------------- //

template< typename _AgentType_ >
void
XPlate::displayAgents
( Fabric<_AgentType_> &inFabric, unsigned long inColor )
{
   unsigned                          iPos                  ,
                                     jPos                  ;
   FabricAliveIterator<_AgentType_>  it = inFabric.begin() ;

   while( *it )
   {
      // Retrieve agent coordinates
      iPos = (*it)->getCell()->getIPos() ;
      jPos = (*it)->getCell()->getJPos() ;

      // Display bacterias in Yellow
      XSetForeground( dpy_, gcontext_, inColor );

      // Display the rectangle
      XFillRectangle( dpy_                   ,
                      dbBuf_                 ,
                      gcontext_              ,
                      (int) (jPos*xCelSize_) ,
                      (int) (iPos*yCelSize_) ,
                      (int) xCelSize_        ,
                      (int) yCelSize_
                    );

      // Go to the next Agent
      ++it ;
   }
}


// -------------------------------------------------------------------------- //

#endif                        // XPLATE_HPP


// ------------------------------ End Of File ------------------------------- //
