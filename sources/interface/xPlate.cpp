// -------------------------------------------------------------------------- //
// xPlate.cpp                                                                 //
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

#include "xPlate.hpp"


// -------------------------------------------------------------------------- //
// Libraries                                                                 //
// -------------------------------------------------------------------------- //

#include "../engine/myException.hpp"
#include "../engine/engineParameters.hpp"
#include "../model/organism.hpp"
#include "../model/baitParameters.hpp"
#include "interfaceParameters.hpp"
#include <stdio.h>
#include <unistd.h>
#include <cmath>
using std::fabs;


// -------------------------------------------------------------------------- //
// Define key codes                                                           //
// -------------------------------------------------------------------------- //

#define END_KEY    0xff57
#define SPACE_KEY  0x20
#define UP_KEY     0xFF55
#define DOWN_KEY   0xFF56


// -------------------------------------------------------------------------- //
// Enumeration                                                                //
// -------------------------------------------------------------------------- //

enum { BLACK, RED_X, GREEN_X, YELLOW, BLUE_X, MAGENTA, CYAN, WHITE };


// -------------------------------------------------------------------------- //
// Static member creation                                                     //
// -------------------------------------------------------------------------- //

const InterfaceParameters* XPlate::xParamsPtr_ = 0 ;


// -------------------------------------------------------------------------- //
// Constructor                                                                //
// -------------------------------------------------------------------------- //

XPlate::XPlate()
:Plate(), 
 numWindows(1)
{
   lastLoopNumber_   = 0 ;
   lastTime_.tv_usec = 0 ;
   lastTime_.tv_sec  = 0 ;
   pause_ = false ;
   dpy_ = NULL ;
}


// -------------------------------------------------------------------------- //
// Destructor                                                                 //
// -------------------------------------------------------------------------- //

XPlate::~XPlate()
{
   // Closing the window
   if( dpy_ )
      XCloseDisplay( dpy_ );
}


// -------------------------------------------------------------------------- //
// Initialize the static parameters pointer                                   //
// -------------------------------------------------------------------------- //

void XPlate::initXParamsPtr( const InterfaceParameters *inPtr )
{
   xParamsPtr_ = inPtr ;
}


// -------------------------------------------------------------------------- //
// Initialise the graphical window                                            //
// -------------------------------------------------------------------------- //

void XPlate::init(int rank, int mpiSize, int *xy_elevs, int *agentLocations, int *substrateLocations, float *waterTempIn, float *solarRadIn, float *dayLengthIn)
{
   Window            root   ;
   int               ecran  ;
   unsigned long     bpx    ;

   // Initialise the plate
   Plate::init(rank, mpiSize, xy_elevs, agentLocations, substrateLocations, waterTempIn, solarRadIn, dayLengthIn) ;
   
   int mpiDimensions = sqrt(numPlatesMPI);

   // If a graphical window is needed
   if( engParamsPtr_->displayPeriod_ )
   {
      // Initial display is 1:1
      windowWidth_  = (unsigned) (engParamsPtr_->width_*mpiDimensions * xParamsPtr_->scale_) ;
      windowHeight_ = (unsigned) (engParamsPtr_->height_*mpiDimensions * xParamsPtr_->scale_) ;

      // Determine size of graphical cells
      xCelSize_ = (double) windowWidth_  / (engParamsPtr_->width_*mpiDimensions)  ;
      yCelSize_ = (double) windowHeight_ / (engParamsPtr_->height_*mpiDimensions) ;

      // Open the display
      dpy_ = XOpenDisplay( 0 ) ;
      if( !dpy_ )
         throw myException( "Cannot open display", __FILE__, __LINE__ ) ;

      // Retrieve graphic parameters
      ecran    = DefaultScreen     ( dpy_ );
      root     = DefaultRootWindow ( dpy_ );
      bpx      = BlackPixel        ( dpy_, ecran );

      // To Do: specify where to put windows using xPos and yPos:
      //int xPos = new int[numWindows];
      //int yPos = new int[numWindows];

      for (unsigned i=0; i<numWindows; i++)
      {
         // Create the main window
         win[i] = XCreateSimpleWindow( dpy_             ,
                                     root             ,
                                     0                ,
                                     0                ,
                                     windowWidth_     ,
                                     windowHeight_    ,
                                     0                ,
                                     bpx              ,
                                     bpx
                                    );
      }

      // Create a colormap
      cmap_ = DefaultColormap( dpy_, ecran );

      for (unsigned i=0; i<numWindows; i++)
         XSetWindowColormap( dpy_, win[i], cmap_ );

      initColors();

      // Select event inputs
      for (unsigned i=0; i<numWindows; i++)
         XSelectInput( dpy_, win[i], KeyPressMask | StructureNotifyMask );

      // Create a black background
      gcv_.background = namedColors_[BLACK];
      gcmask_         = GCForeground | GCBackground ;
      for (unsigned i=0; i<numWindows; i++)
      {
         gcontext[i] = XCreateGC( dpy_, win[i], gcmask_, &gcv_ );

         // Display the three windows
         XMapWindow( dpy_, win[i] );

         // Create buffer for double buffering display 1
         dbBuf[i] = XCreatePixmap( dpy_                    ,
                                 win[i]                  ,
                                 engParamsPtr_->width_*mpiDimensions   ,
                                 engParamsPtr_->height_*mpiDimensions  ,
                                 DefaultDepth( dpy_, DefaultScreen( dpy_ ) )
                                );
      }
   }
}


// -------------------------------------------------------------------------- //
// Refresh the window with a new display                                      //
// -------------------------------------------------------------------------- //

void XPlate::display( unsigned loopNumber, Fabric<Organism>& inOrganismsFabric, int mpiSize, int rank)
{
   char              title[256] ;
   struct timeval    newTime    ;
   double            lps        ;

   // Check window event
   checkWindowEvent() ;

   // Calculate loops per second
   gettimeofday( &newTime, NULL );
   lps = ( loopNumber - lastLoopNumber_ ) /
         ( ( newTime.tv_sec - lastTime_.tv_sec )
         + (double) ( newTime.tv_usec - lastTime_.tv_usec ) / 1000000 ) ;
   lastTime_        = newTime     ;
   lastLoopNumber_  = loopNumber  ;

   // Display the title
   //sprintf( title, "MG(R%d)-L:%d-LPS:%6.2f", rank, loopNumber,lps);

   for (unsigned i=0; i<numWindows; i++)
   {
      unsigned yearLength = 8766; // compete orbit of sun = 8766.15 h
      unsigned monthLength = 730;
      int currYear = (loopNumber/yearLength);
      int currHour = loopNumber % yearLength;
      int hourInMonth = currHour % (monthLength);
      int currDay = hourInMonth/24;
      
      //cout << currYear << ", " << currMonth << ", " << currHour << ", " << hourInMonth << ", " << currDay << endl;
   
      // Display the title
      //sprintf( title, "Mg(r%d)-L:%d Month:%d (lps:%6.2f) Win%d", rank, loopNumber, currMonth, lps, (i+1));

      sprintf( title, "UndariaGEN\tDate: %d/%d/%d\t\tTemp:%6.2f\t\tSolar:%6.2f\t\tLoop: %d (lps:%6.2f)", (currDay+1), (currMonth+1), (2014+currYear), currWaterTemp, currSolarRad, loopNumber, lps );
      
      XStoreName( dpy_, win[i], title );

      // Display the background
      XSetForeground( dpy_, gcontext[i], namedColors_[BLACK] );
      XFillRectangle( dpy_, dbBuf[i], gcontext[i], 0, 0, windowWidth_, windowHeight_ );
      
      if (i == 0)
      {
         // Display nutrient levels and bacteria biomass in window 1.
         /*if (loopNumber > baitParamsPtr_->loopAddAntibio)
         {
            displayAntibios(i, 0);
            displaySpore(i);
         }
         else*/
         if (mpiSize>1)
         {
            //displayAllPatchesMPI(i);
            displayAllElevationsMPI(i);
            displayWaterTempMPI(i);
            displaySporeMPI(i);
            displaySubstrateTypesMPI(i);
            //displayOrganisms(i);
            displayAllBiomassMPI(i);
         }
         else
         {
            displayElevations(i);
            displayWaterTemp(i);
            //displayWaterDepth(i);
            displaySpore(i);
            displaySubstrateTypes(i);
            displayOrganisms(i);
         }
      }
      else
      {
         // Display in window 2.
         //displayPatches(i);
         //displayAntibios(i, 0);
         //displayCharge(i);
      }


      // Refresh display
      XCopyArea( dpy_           ,
                 dbBuf[i]       ,
                 win[i]         ,
                 gcontext[i]    ,
                 0              ,
                 0              ,
                 windowWidth_   ,
                 windowHeight_  ,
                 0              ,
                 0
                );
   }

   XFlush( dpy_ );

}


// -------------------------------------------------------------------------- //
// Create color array                                                         //
// -------------------------------------------------------------------------- //

void XPlate::initColors()
{
   int     nbColor = 8 ;
   XColor  xc;

   // Pixel array allocation
   namedColors_ = (unsigned long *) malloc ( nbColor * sizeof( unsigned long) );

   // Allocate 8 colors
   XAllocNamedColor( dpy_, cmap_, "black",   &xc, &xc );
   namedColors_[ BLACK   ] = xc.pixel;
   XAllocNamedColor( dpy_, cmap_, "red",     &xc, &xc );
   namedColors_[ RED_X     ] = xc.pixel;
   XAllocNamedColor( dpy_, cmap_, "green",   &xc, &xc );
   namedColors_[ GREEN_X   ] = xc.pixel;
   XAllocNamedColor( dpy_, cmap_, "yellow",  &xc, &xc );
   namedColors_[ YELLOW  ] = xc.pixel;
   XAllocNamedColor( dpy_, cmap_, "blue",    &xc, &xc );
   namedColors_[ BLUE_X    ] = xc.pixel;
   XAllocNamedColor( dpy_, cmap_, "magenta", &xc, &xc );
   namedColors_[ MAGENTA ] = xc.pixel;
   XAllocNamedColor( dpy_, cmap_, "cyan",    &xc, &xc );
   namedColors_[ CYAN    ] = xc.pixel;
   XAllocNamedColor( dpy_, cmap_, "white",   &xc, &xc );
   namedColors_[ WHITE   ] = xc.pixel;

   unsigned darkest = xParamsPtr_->nbPatchColors_/8;

   // Pixel array allocation for patches
   patchesColors1_ = (unsigned long *) malloc( xParamsPtr_->nbPatchColors_
                                       * sizeof(unsigned long) );
   for( unsigned i=0; i<xParamsPtr_->nbPatchColors_; ++i )
   {
      xc.red   = ((int)(i*0.2)*256*256)/xParamsPtr_->nbPatchColors_ ;
      xc.green = ((int)(i*0.8)*256*256)/xParamsPtr_->nbPatchColors_ ;
      xc.blue  = (i*256*256)/xParamsPtr_->nbPatchColors_ ;
      XAllocColor( dpy_, cmap_, &xc );
      patchesColors1_[i] = xc.pixel ;
   }
   patchesColors2_ = (unsigned long *) malloc( xParamsPtr_->nbPatchColors_
                                       * sizeof(unsigned long) );
   for( unsigned i=0; i<xParamsPtr_->nbPatchColors_; ++i )
   {
      xc.red   = (i*256*256)/xParamsPtr_->nbPatchColors_ ;
      xc.green = ((int)(i*0.7)*256*256)/xParamsPtr_->nbPatchColors_ ;
      xc.blue  = ((int)(i*0.5)*256*256)/xParamsPtr_->nbPatchColors_ ;
      XAllocColor( dpy_, cmap_, &xc );
      patchesColors2_[i] = xc.pixel ;
   }
   
   // Pixel array allocation for medium nutrient patches
   patchesMedColors_ = (unsigned long *) malloc( xParamsPtr_->nbPatchColors_
                                       * sizeof(unsigned long) );
   for( unsigned i=0; i<xParamsPtr_->nbPatchColors_; ++i )
   {
      xc.red   = (0*256*256)/xParamsPtr_->nbPatchColors_ ;
      xc.green = (i*256*256)/xParamsPtr_->nbPatchColors_ ;
      xc.blue  = (0*256*256)/xParamsPtr_->nbPatchColors_ ;
      XAllocColor( dpy_, cmap_, &xc );
      patchesMedColors_[i] = xc.pixel ;
   }

   // Pixel array allocation for low nutrient patches
   patchesLowColors_ = (unsigned long *) malloc( xParamsPtr_->nbPatchColors_
                                       * sizeof(unsigned long) );
   for( unsigned i=0; i<xParamsPtr_->nbPatchColors_; ++i )
   {
      xc.red   = (i*256*256)/xParamsPtr_->nbPatchColors_ ;
      xc.green = (0*256*256)/xParamsPtr_->nbPatchColors_ ;
      xc.blue  = (0*256*256)/xParamsPtr_->nbPatchColors_ ;
      XAllocColor( dpy_, cmap_, &xc );
      patchesLowColors_[i] = xc.pixel ;
   }
   

   // Pixel array allocation for antibios
   /*antibiosColours = (unsigned long *) malloc( xParamsPtr_->nbPatchColors_
                                       * sizeof(unsigned long) );
   for( unsigned i=0; i<xParamsPtr_->nbPatchColors_; ++i )
   {
      xc.red   = (0*256*256)/xParamsPtr_->nbPatchColors_ ;
      xc.green = (0*256*256)/xParamsPtr_->nbPatchColors_ ;
      xc.blue  = (i*256*256)/xParamsPtr_->nbPatchColors_ ;
      XAllocColor( dpy_, cmap_, &xc );
      antibiosColours[i] = xc.pixel ;
   }

   // Pixel array allocation for antibios
   actAntibiosColours = (unsigned long *) malloc( xParamsPtr_->nbPatchColors_
                                       * sizeof(unsigned long) );
   for( unsigned i=0; i<xParamsPtr_->nbPatchColors_; ++i )
   {
      xc.red   = (0*256*256)/xParamsPtr_->nbPatchColors_ ;
      xc.green = (0*256*256)/xParamsPtr_->nbPatchColors_ ;
      xc.blue  = ((darkest+i)*256*256)/(xParamsPtr_->nbPatchColors_+darkest) ;
      XAllocColor( dpy_, cmap_, &xc );
      actAntibiosColours[i] = xc.pixel ;
   }*/

   // Pixel array allocation for spores.
   sporeColours = (unsigned long *) malloc( xParamsPtr_->nbPatchColors_
                                       * sizeof(unsigned long) );
   for( unsigned i=0; i<xParamsPtr_->nbPatchColors_; ++i )
   {
      xc.red   = (0*256*256)/xParamsPtr_->nbPatchColors_ ;
      xc.green = ((darkest+i)*256*256)/(xParamsPtr_->nbPatchColors_+darkest) ;
      xc.blue  = (0*256*256)/xParamsPtr_->nbPatchColors_ ;
      XAllocColor( dpy_, cmap_, &xc );
      sporeColours[i] = xc.pixel ;
   }

   // Pixel array allocation for bacterias
   bacteriasColours_ = (unsigned long *) malloc( xParamsPtr_->nbPatchColors_
                                       * sizeof(unsigned long) );
   for( unsigned i=0; i<xParamsPtr_->nbOrganismColors_; ++i )
   {
      //if (i > 0)
      //{
         xc.red   = (i*256*256)/xParamsPtr_->nbPatchColors_ ;
         xc.green = ((int)(i*0.54)*256*256)/xParamsPtr_->nbPatchColors_ ;
         xc.blue  = ((int)(i*0.1)*256*256)/xParamsPtr_->nbPatchColors_ ;
      /*}
      else
      {
         xc.red   = ((xParamsPtr_->nbOrganismColors_-1)*256*256)/xParamsPtr_->nbPatchColors_ ;
         xc.green = (0*256*256)/xParamsPtr_->nbPatchColors_ ;
         xc.blue  = (0*256*256)/xParamsPtr_->nbPatchColors_ ;
      }*/
      XAllocColor( dpy_, cmap_, &xc );
      bacteriasColours_[i] = xc.pixel ;
   }

   // Pixel array allocation for electrostatic charges.
   elevColours = (unsigned long *) malloc( xParamsPtr_->nbPatchColors_
                                       * sizeof(unsigned long) );
   for( unsigned i=0; i<xParamsPtr_->nbPatchColors_; ++i )
   {
      xc.red   = ((darkest+i)*256*256)/(xParamsPtr_->nbPatchColors_+darkest) ;
      xc.green = (0*256*256)/xParamsPtr_->nbPatchColors_ ;
      xc.blue  = (0*256*256)/xParamsPtr_->nbPatchColors_ ;
      XAllocColor( dpy_, cmap_, &xc );
      elevColours[i] = xc.pixel ;
   }

   // Pixel array allocation for electrostatic charges.
   /*posChargeColours = (unsigned long *) malloc( xParamsPtr_->nbPatchColors_
                                       * sizeof(unsigned long) );
   for( unsigned i=0; i<xParamsPtr_->nbPatchColors_; ++i )
   {
      xc.red   = ((darkest+i)*256*256)/(xParamsPtr_->nbPatchColors_+darkest) ;
      xc.green = ((darkest+i)*256*256)/(xParamsPtr_->nbPatchColors_+darkest) ;
      xc.blue  = (0*256*256)/(xParamsPtr_->nbPatchColors_) ;
      XAllocColor( dpy_, cmap_, &xc );
      posChargeColours[i] = xc.pixel ;
   }*/
}


// -------------------------------------------------------------------------- //
// Retrieve last pressed key and modify display parameter if it is needed.    //
// Launch an exception if the program must stops.                             //
// -------------------------------------------------------------------------- //

void XPlate::keyHandler( XKeyEvent* ev )
{
   KeySym     ks ;

   // Retrieve last key pressed
   ks = XLookupKeysym( ev, 0 );

   // Apply modification
   switch( ks )
   {
      case END_KEY :    // "End" Key
         throw myException( "End of program", __FILE__, __LINE__ ) ;

      case SPACE_KEY :  // Space Bar
         pause_ = !pause_ ;
         break ;

      case UP_KEY :     // "Page Up" Key
         if( engParamsPtr_->displayPeriod_ > 1)
            --engParamsPtr_->displayPeriod_ ;
         break ;

      case DOWN_KEY :   // "Page Down" Key
         ++engParamsPtr_->displayPeriod_ ;
         break ;
   }
}


// -------------------------------------------------------------------------- //
// Check the event list and then treat it. Return when no event left.         //
// -------------------------------------------------------------------------- //

void XPlate::checkWindowEvent()
{
   XEvent  ev ;
   int mpiDimensions = sqrt(numPlatesMPI);

   do
   {
      for (unsigned i=0; i<numWindows; i++)
      {
         // Check window event
         while( XCheckWindowEvent( dpy_, win[i], KeyPressMask|StructureNotifyMask ,
            &ev ) )
         {
            switch( ev.type )
            {
               case ConfigureNotify :

                  // New window size
                  windowWidth_  = ev.xconfigure.width  ;
                  windowHeight_ = ev.xconfigure.height ;
                  // New display buffer
                  XFreePixmap( dpy_, dbBuf[i] ) ;
                  dbBuf[i] = XCreatePixmap( dpy_, win[i], windowWidth_, windowHeight_,
                     DefaultDepth( dpy_, DefaultScreen( dpy_ ) ) );

                  // New cells size
                  xCelSize_ = (double) windowWidth_  / (engParamsPtr_->width_*mpiDimensions)  ;
                  yCelSize_ = (double) windowHeight_ / (engParamsPtr_->height_*mpiDimensions) ;

                  break ;

               case KeyPress :
                  // Look at the pressed key
                  keyHandler( (XKeyEvent*)&ev ) ;
                  break ;
            }
         }
      }

      // Avoid fast looping while pause
      if( pause_ )
         usleep( 100000 ) ;

   } while( pause_ ) ;
}


// -------------------------------------------------------------------------- //
// Display patches with color depending of the food level of each             //
//   Window 1
// -------------------------------------------------------------------------- //

/*void XPlate::displayPatches(unsigned winNum)
{
   unsigned     nbColors      ;
   double       nutrientLevel,
                max;
   bool isObstacle;

   // Preloop initialisation
   max = ( baitParamsPtr_->maxTotalNutrientLevel_ )*2;
   nbColors = xParamsPtr_->nbPatchColors_  ;

   // For each patch
   for( int i=0; i<engParamsPtr_->height_; ++i )
   {
      for( int j=0; j<engParamsPtr_->width_; ++j )
      {
         if ( (*this)(i,j).getActiveRegion() )      // Only display patches in active zone.
         {
            // Determine the color
            nutrientLevel = (*this)(i,j).getNutrientLevel() ;
            isObstacle = (*this)(i,j).getIsObstacle();

            if (isObstacle)
            {
               XSetForeground( dpy_, gcontext[winNum], bacteriasColours_[0] );
            }
            else
            {
               if (nutrientLevel < max)
                  XSetForeground( dpy_, gcontext[winNum],
                     patchesColors1_[ (unsigned) ( nutrientLevel*(nbColors-1)/max) ] );
               else if (nutrientLevel >= max)
                  XSetForeground( dpy_, gcontext[winNum],
                     patchesColors1_[ (unsigned)nbColors-1 ] );
               else
                  XSetForeground( dpy_, gcontext[winNum], namedColors_[WHITE] );
            }

            // Display the patch
            XFillRectangle( dpy_                 ,
                             dbBuf[winNum]        ,
                            gcontext[winNum]     ,
                            (int) (j*xCelSize_)  ,
                            (int) (i*yCelSize_)  ,
                            (int) xCelSize_ + 1  ,
                            (int) yCelSize_ + 1
                           );
         }
      }
   }
}*/

//
// Code to display patches from all plates in single window (MPI)
//
void XPlate::displayAllPatchesMPI(unsigned winNum)
{
   unsigned     nbColors      ;
   double       nutrientLevel,
                max;
   //int elevation;
   //bool isObstacle;
   int xPos=0, yPos=0;

   // Preloop initialisation
   max = ( baitParamsPtr_->maxTotalNutrientLevel_ )*2;
   nbColors = xParamsPtr_->nbPatchColors_  ;

   int totArraySize = engParamsPtr_->height_* engParamsPtr_->width_ * numPlatesMPI;
   int plateSize = engParamsPtr_->height_* engParamsPtr_->width_;
   int rankOfPlate = 0;
   int plateDim = sqrt(numPlatesMPI);
   int plateRow = rankOfPlate/plateDim;
   //int plateCol = rankOfPlate % plateDim;
   int width = engParamsPtr_->width_;
   int height = engParamsPtr_->height_;

   for (int i=0; i<totArraySize; i++)
   {
      nutrientLevel = allPatchValues[i];
      
      if (i > 1 && i % plateSize == 0)
         rankOfPlate++; 
              
      if (nutrientLevel < max)
         XSetForeground( dpy_, gcontext[winNum],
            patchesColors1_[ (unsigned) ( nutrientLevel*(nbColors-1)/max) ] );
      else if (nutrientLevel >= max)
         XSetForeground( dpy_, gcontext[winNum],
            patchesColors1_[ (unsigned)nbColors-1 ] );
      else
         XSetForeground( dpy_, gcontext[winNum], namedColors_[WHITE] );
      
      plateRow = rankOfPlate/plateDim;
      xPos = (width * (rankOfPlate - plateRow*plateDim)) + (i % width);
      yPos = (height * plateRow) + (i - (plateSize*rankOfPlate)) / width; 
         
      // Display the patch
      XFillRectangle( dpy_                 ,
                       dbBuf[winNum]        ,
                      gcontext[winNum]     ,
                      (int) (xPos*xCelSize_)  ,
                      (int) (yPos*yCelSize_)  ,
                      (int) xCelSize_ + 1  ,
                      (int) yCelSize_ + 1
                     );
   }

}

void XPlate::displayElevations(unsigned winNum)
{
   unsigned     nbColors      ;
   float       elevation,
                max;
   //bool         isSanglicaNiche;

   // Preloop initialisation
   max = 30;
   nbColors = xParamsPtr_->nbPatchColors_  ;

   // For each patch
   for( unsigned i=0; i<engParamsPtr_->height_; ++i )
   {
      for( unsigned j=0; j<engParamsPtr_->width_; ++j )
      {
         if ( (*this)(i,j).getActiveRegion() )      // Only display patches in active zone.
         {
            // Determine the color
            elevation = max - (*this)(i,j).getElevation();
            //isSanglicaNiche = (*this)(i,j).getIsSanglicaNiche();
            
            unsigned colorNum = (elevation*(nbColors-1)/(max*2)); 
            colorNum = nbColors - colorNum;
            
            //if ((*this)(i,j).getElevation() > baitParamsPtr_->lowerLimit && (*this)(i,j).getElevation() < baitParamsPtr_->upperLimit)
               XSetForeground( dpy_, gcontext[winNum],
                  elevColours[ colorNum] );
            //else
            //   XSetForeground( dpy_, gcontext[winNum],
            //      patchesColors2_[ colorNum] );

            // Display the patch
            XFillRectangle( dpy_                 ,
                             dbBuf[winNum]        ,
                            gcontext[winNum]     ,
                            (int) (j*xCelSize_)  ,
                            (int) (i*yCelSize_)  ,
                            (int) xCelSize_ + 1  ,
                            (int) yCelSize_ + 1
                           );
         }
      }
   }
}


//
// Code to display patches from all plates in single window (MPI)
//
void XPlate::displayAllElevationsMPI(unsigned winNum)
{
   unsigned     nbColors      ;
   double       elevation,
                max;
   //int elevation;
   //bool isObstacle;
   int xPos=0, yPos=0;

   // Preloop initialisation
   max = 20;
   nbColors = xParamsPtr_->nbPatchColors_  ;

   unsigned totArraySize = engParamsPtr_->height_* engParamsPtr_->width_ * numPlatesMPI;
   int plateSize = engParamsPtr_->height_* engParamsPtr_->width_;
   int rankOfPlate = 0;
   int plateDim = sqrt(numPlatesMPI);
   int plateRow = rankOfPlate/plateDim;
   //int plateCol = rankOfPlate % plateDim;
   int width = engParamsPtr_->width_;
   int height = engParamsPtr_->height_;

   for (unsigned i=0; i<totArraySize; i++)
   {
      //elevation = (double)allElevationValues[i];          
      elevation = max - (double)allElevationValues[i];
      
      if (i > 1 && i % plateSize == 0)
         rankOfPlate++;      
                
      unsigned colorNum = (elevation*(nbColors-1)/(max*2)); 
      colorNum = nbColors - colorNum;       
      XSetForeground( dpy_, gcontext[winNum],
            patchesColors2_[ colorNum ] );
      
      plateRow = rankOfPlate/plateDim;
      xPos = (width * (rankOfPlate - plateRow*plateDim)) + (i % width);
      yPos = (height * plateRow) + (i - (plateSize*rankOfPlate)) / width; 
         
      // Display the patch
      XFillRectangle( dpy_                 ,
                       dbBuf[winNum]        ,
                      gcontext[winNum]     ,
                      (int) (xPos*xCelSize_)  ,
                      (int) (yPos*yCelSize_)  ,
                      (int) xCelSize_ + 1  ,
                      (int) yCelSize_ + 1
                     );
   }
}

//
// Display the water depth (in metres)
//
void XPlate::displayWaterDepth(unsigned winNum)
{
   unsigned     nbColors      ;
   float       waterDepth,
                max;
   //bool         isSanglicaNiche;

   // Preloop initialisation
   max = 30;
   nbColors = xParamsPtr_->nbPatchColors_  ;

   // For each patch
   for( unsigned i=0; i<engParamsPtr_->height_; ++i )
   {
      for( unsigned j=0; j<engParamsPtr_->width_; ++j )
      {
         waterDepth = (*this)(i,j).getWaterDepth();
         if ( (*this)(i,j).getActiveRegion() && waterDepth > 0.0)      // Only display patches in active zone.
         {
            //isSanglicaNiche = (*this)(i,j).getIsSanglicaNiche();
            
            XSetForeground( dpy_, gcontext[winNum],
                     patchesColors1_[ (unsigned)(waterDepth*(nbColors-1)/max)] );

            // Display the patch
            XFillRectangle( dpy_                 ,
                             dbBuf[winNum]        ,
                            gcontext[winNum]     ,
                            (int) (j*xCelSize_)  ,
                            (int) (i*yCelSize_)  ,
                            (int) xCelSize_ + 1  ,
                            (int) yCelSize_ + 1
                           );
         }
      }
   }
}

//
// 20Jan14 - Display the sea water temperature 
//
void XPlate::displayWaterTemp(unsigned winNum)
{
   unsigned     nbColors      ;
   float       waterTempX,
                max;
   //bool         isSanglicaNiche;

   // Preloop initialisation
   max = 25;
   nbColors = xParamsPtr_->nbPatchColors_  ;
   waterTempX = (*this).getWaterTemp();

   // For each patch
   for( unsigned i=0; i<engParamsPtr_->height_; ++i )
   {
      for( unsigned j=0; j<engParamsPtr_->width_; ++j )
      {
         if ( (*this)(i,j).getActiveRegion() && waterTempX > 0.0)      // Only display patches in active zone.
         {
            //isSanglicaNiche = (*this)(i,j).getIsSanglicaNiche();
            
            XSetForeground( dpy_, gcontext[winNum],
                     patchesColors1_[ (unsigned)(waterTempX*(nbColors-1)/max)] );

            // Display the patch
            XFillRectangle( dpy_                 ,
                             dbBuf[winNum]        ,
                            gcontext[winNum]     ,
                            (int) (j*xCelSize_)  ,
                            (int) (i*yCelSize_)  ,
                            (int) xCelSize_ + 1  ,
                            (int) yCelSize_ + 1
                           );
         }
      }
   }
}

void XPlate::displayWaterTempMPI(unsigned winNum)
{
   unsigned nbColors;
   float waterTempX;
   int xPos=0, yPos=0;
   int totArraySize = engParamsPtr_->height_* engParamsPtr_->width_ * numPlatesMPI;
   int plateSize = engParamsPtr_->height_* engParamsPtr_->width_;
   int rankOfPlate = 0;
   int plateDim = sqrt(numPlatesMPI);
   int plateRow = 1;//rankOfPlate/plateDim;
   int width = engParamsPtr_->width_;
   int height = engParamsPtr_->height_;
   
   float max = 25;
   nbColors = xParamsPtr_->nbPatchColors_  ;
   waterTempX = (*this).getWaterTemp();

   for (int i=0; i<totArraySize; i++)
   {      
      if (i > 1 && i % plateSize == 0)
         rankOfPlate++; 
                 
      
      if (waterTempX > 0.0)  
      {       
         XSetForeground( dpy_, gcontext[winNum],
                  patchesColors1_[ (unsigned)(waterTempX*(nbColors-1)/max)] );
            
               
         for(int j=1; j<=sqrt(numPlatesMPI); j++)   
         {
            if (i < plateSize * (j*sqrt(numPlatesMPI)) )
            {
               plateRow = j;
               break;
            }
         }
      
         xPos = (width * (rankOfPlate - (plateRow-1)*plateDim)) + (i % width);
         yPos = (height * (plateRow-1)) + (i - (plateSize*rankOfPlate)) / width;

         // Display the patch
         XFillRectangle( dpy_                 ,
                          dbBuf[winNum]        ,
                         gcontext[winNum]     ,
                         (int) (xPos*xCelSize_)  ,
                         (int) (yPos*yCelSize_)  ,
                         (int) xCelSize_ + 1  ,
                         (int) yCelSize_ + 1
                        );
      }
   }
}



// -------------------------------------------------------------------------- //
// Display spore levels of each patch.                        //
// -------------------------------------------------------------------------- //

void XPlate::displaySpore(unsigned winNum)
{
   unsigned      nbColors;
   float        sporeLevelX, 
                 max, min,
                 logRatio;

   // Preloop initialisation
   max = baitParamsPtr_->sporeProd * 1.0E+04;
   min = 10.0;
   nbColors = xParamsPtr_->nbPatchColors_  ;

   // For each patch
   for( unsigned i=0; i<engParamsPtr_->height_; ++i )
   {
      for( unsigned j=0; j<engParamsPtr_->width_; ++j )
      {
         if ( (*this)(i,j).getActiveRegion() && (*this)(i,j).getSubstrateType() >= 0)      // Only display patches in active zone.
         {
            // Determine the color
            sporeLevelX = (*this)(i,j).getSporeLevel() ;

            //cout << sporeLevel << endl;
            // Only display when significant amount present (e.g. >10)
            if (sporeLevelX > (min+1.0))
            {
               //cout << "Here " << sporeLevelX << endl;
               logRatio = log(sporeLevelX-min)/log(max-min);
               //cout << "Log ratio = " << log(sporeLevelX-min) << " / " << log(max-min) << " = " << logRatio << endl;
               if (logRatio > 1.0)
                  XSetForeground( dpy_, gcontext[winNum], sporeColours[nbColors-1] );
               else
                  //XSetForeground( dpy_, gcontext[winNum],
                     //sporeColours[ (unsigned) ( sporeLevelX*(nbColors-1)/max) ] );
                  XSetForeground( dpy_, gcontext[winNum],
                     sporeColours[ (unsigned) (logRatio * (nbColors-1)) ] );

               // Display the patch
               XFillRectangle( dpy_                 ,
                            dbBuf[winNum]       ,
                            gcontext[winNum]     ,
                            (int) (j*xCelSize_)  ,
                            (int) (i*yCelSize_)  ,
                            (int) xCelSize_ + 1  ,
                            (int) yCelSize_ + 1
                              );
            }
         }
      }
   }
}

//
// Code to display patches from all plates in single window (MPI)
//
void XPlate::displaySporeMPI(unsigned winNum)
{
   unsigned nbColors;
   float spore;
   int xPos=0, yPos=0;
   int totArraySize = engParamsPtr_->height_* engParamsPtr_->width_ * numPlatesMPI;
   int plateSize = engParamsPtr_->height_* engParamsPtr_->width_;
   int rankOfPlate = 0;
   int plateDim = sqrt(numPlatesMPI);
   int plateRow = 1;//rankOfPlate/plateDim;
   int width = engParamsPtr_->width_;
   int height = engParamsPtr_->height_;
   
   float max, min, logRatio;
   max = baitParamsPtr_->sporeProd * 1.0E+04;
   min = 10.0;
   nbColors = xParamsPtr_->nbPatchColors_  ;

   for (int i=0; i<totArraySize; i++)
   {
      spore = allSporeValues[i];

      if (i > 1 && i % plateSize == 0)
         rankOfPlate++; 
       
      // Only display if organism present (>0)  
      if (spore > min + 1.0)
      {     
         logRatio = log(spore-min)/log(max-min);
         
         if (logRatio > 1.0)
            XSetForeground( dpy_, gcontext[winNum], sporeColours[nbColors-1] );
         else
            XSetForeground( dpy_, gcontext[winNum],
               sporeColours[ (unsigned) (logRatio * (nbColors-1)) ] );
               
               
         for(int j=1; j<=sqrt(numPlatesMPI); j++)   
         {
            if (i < plateSize * (j*sqrt(numPlatesMPI)) )
            {
               plateRow = j;
               break;
            }
         }
      
         xPos = (width * (rankOfPlate - (plateRow-1)*plateDim)) + (i % width);
         yPos = (height * (plateRow-1)) + (i - (plateSize*rankOfPlate)) / width;

         // Display the patch
         XFillRectangle( dpy_                 ,
                          dbBuf[winNum]        ,
                         gcontext[winNum]     ,
                         (int) (xPos*xCelSize_)  ,
                         (int) (yPos*yCelSize_)  ,
                         (int) xCelSize_ + 1  ,
                         (int) yCelSize_ + 1
                        );
      }
   }
}



// -------------------------------------------------------------------------- //
// Display all alive bacterias                                                //
// -------------------------------------------------------------------------- //

void XPlate::displayOrganisms( unsigned winNum )
{
   unsigned   nbColors;
   float      biomass,
              maxBiomass;
   
   // Preloop initialisation
   maxBiomass = 1.0;//baitParamsPtr_->maxBiomassInPatch;
   nbColors = xParamsPtr_->nbOrganismColors_  ;

   // Cycle through each patch
   for( unsigned i=0; i<engParamsPtr_->height_; ++i )
   {
      for( unsigned j=0; j<engParamsPtr_->width_; ++j )
      {
         if ((*this)(i,j).getIsOrganism() == true)
         {
            // Determine the colour based on total biomass per patch
            biomass = (*this)(i,j).calcSporoBiomass();
            unsigned bacColour = (unsigned) (2+biomass*(nbColors-1-2)/maxBiomass);
            if (bacColour > (nbColors-1) )
               bacColour = (nbColors-1);
            
            //if ((*this)(i,j).getNutrientLevel() > 500 && bacColour == 0)
            //   XSetForeground( dpy_, gcontext[winNum], GREEN_X );
            //else
               //XSetForeground( dpy_, gcontext[winNum], bacteriasColours_[bacColour] );
            if (biomass > 1.0)
               XSetForeground( dpy_, gcontext[winNum], namedColors_[RED_X] );
            else if (biomass > baitParamsPtr_->macroSize)
               XSetForeground( dpy_, gcontext[winNum], namedColors_[BLUE_X] );
            
            if (biomass > 0.1)
               // Display the bacteria
               XFillRectangle( dpy_                 ,
                               dbBuf[winNum]               ,
                               gcontext[winNum]           ,
                               (int) (j*xCelSize_)  ,
                               (int) (i*yCelSize_)  ,
                               (int) xCelSize_ + 1  ,
                               (int) yCelSize_ + 1
                               );
         }
      }
   }
}


//
// Code to display organism biomass levels from all plates in single window (MPI)
//
void XPlate::displayAllBiomassMPI(unsigned winNum)
{
   float biomass;
   int xPos=0, yPos=0;
   int totArraySize = engParamsPtr_->height_* engParamsPtr_->width_ * numPlatesMPI;
   int plateSize = engParamsPtr_->height_* engParamsPtr_->width_;
   int rankOfPlate = 0;
   int plateDim = sqrt(numPlatesMPI);
   int plateRow = 1;//rankOfPlate/plateDim;
   int width = engParamsPtr_->width_;
   int height = engParamsPtr_->height_;

   for (int i=0; i<totArraySize; i++)
   {
      biomass = allBiomassValues[i];

      if (i > 1 && i % plateSize == 0)
         rankOfPlate++; 
       
      // Only display if organism present (>0)  
      if (biomass > 1.0E-06)
      {     
         if (biomass > 1.0)
            XSetForeground( dpy_, gcontext[winNum], namedColors_[RED_X] );
         else if (biomass > baitParamsPtr_->macroSize)
            XSetForeground( dpy_, gcontext[winNum], namedColors_[BLUE_X] );
         
         for(int j=1; j<=sqrt(numPlatesMPI); j++)   
         {
            if (i < plateSize * (j*sqrt(numPlatesMPI)) )
            {
               plateRow = j;
               break;
            }
         }
      
         xPos = (width * (rankOfPlate - (plateRow-1)*plateDim)) + (i % width);
         yPos = (height * (plateRow-1)) + (i - (plateSize*rankOfPlate)) / width;

         // Display the patch
         XFillRectangle( dpy_                 ,
                          dbBuf[winNum]        ,
                         gcontext[winNum]     ,
                         (int) (xPos*xCelSize_)  ,
                         (int) (yPos*yCelSize_)  ,
                         (int) xCelSize_ + 1  ,
                         (int) yCelSize_ + 1
                        );
      }
   }
}







//
// Code to display pollen levels from all plates in single window (MPI)
//
/*void XPlate::displayPollenMPI(unsigned winNum)
{
   unsigned     nbColors      ;
   unsigned     pollen,
                maxPollen,
                minPollen;
   //int elevation;
   //bool isObstacle;
   int xPos=0, yPos=0;

   // Preloop initialisation
   minPollen = 100;
   maxPollen = baitParamsPtr_->pollPacketSize*20 + minPollen;
   nbColors = xParamsPtr_->nbOrganismColors_  ;

   unsigned totArraySize = engParamsPtr_->height_* engParamsPtr_->width_ * numPlatesMPI;
   unsigned plateSize = engParamsPtr_->height_* engParamsPtr_->width_;
   int rankOfPlate = 0;

   for (unsigned i=0; i<totArraySize; i++)
   {
      pollen = allPollenValues[i];
      if (i > 1 && i % plateSize == 0)
         rankOfPlate++; 
      
      if (pollen > minPollen)
      { 
         unsigned netPollen = pollen - minPollen;      
         unsigned polColour = (unsigned) (2+netPollen*(nbColors-1-2)/maxPollen);
            if (polColour > (nbColors-1) )
               polColour = (nbColors-1);
         
         XSetForeground( dpy_, gcontext[winNum], pollenColours[polColour] );
         
         // If in the first row of plates
         if ( i < plateSize * sqrt(numPlatesMPI))
         {
            xPos = (engParamsPtr_->width_*rankOfPlate) + (i % engParamsPtr_->width_);
            yPos = (i - (plateSize*rankOfPlate)) / engParamsPtr_->width_;  // yes, width_
         }
         // Row 2 of plates
         else if ( i < plateSize * (2*sqrt(numPlatesMPI)) )
         {
            xPos = ( engParamsPtr_->width_*(rankOfPlate-sqrt(numPlatesMPI)) ) + (i % engParamsPtr_->width_);
            yPos = engParamsPtr_->height_ + (i - (plateSize*rankOfPlate)) / engParamsPtr_->width_;   
         }
         // Row 3 of plates
         else if ( i < plateSize * (3*sqrt(numPlatesMPI)) )
         {
            xPos = ( engParamsPtr_->width_*(rankOfPlate-(3-1)*sqrt(numPlatesMPI)) ) + (i % engParamsPtr_->width_);
            yPos = engParamsPtr_->height_*(3-1) + (i - (plateSize*rankOfPlate)) / engParamsPtr_->width_;
         }
         // Row 4 of plates
         else if ( i < plateSize * (4*sqrt(numPlatesMPI)) )
         {
            xPos = ( engParamsPtr_->width_*(rankOfPlate-(4-1)*sqrt(numPlatesMPI)) ) + (i % engParamsPtr_->width_);
            yPos = engParamsPtr_->height_*(4-1) + (i - (plateSize*rankOfPlate)) / engParamsPtr_->width_;
         }
            
         // Display the patch
         XFillRectangle( dpy_                 ,
                          dbBuf[winNum]        ,
                         gcontext[winNum]     ,
                         (int) (xPos*xCelSize_)  ,
                         (int) (yPos*yCelSize_)  ,
                         (int) xCelSize_ + 1  ,
                         (int) yCelSize_ + 1
                        );
      }
   }

}*/

void XPlate::displaySubstrateTypes(unsigned winNum)
{
   int         substrate;

   // For each patch
   for( unsigned i=0; i<engParamsPtr_->height_; ++i )
   {
      for( unsigned j=0; j<engParamsPtr_->width_; ++j )
      {
         // Only display patches with substrate
         if ( (*this)(i,j).getActiveRegion() && (*this)(i,j).getSubstrateType() != 0)
         {
            substrate = (*this)(i,j).getSubstrateType();
            switch( substrate )
            {
               case -1 :     // Obstacle
                  XSetForeground( dpy_, gcontext[winNum], namedColors_[BLACK]);
                  break ;
               case 0 :       // no substate
                  break;
               case 1 :  // High Affinity substrate (e.g. concrete)
                  XSetForeground( dpy_, gcontext[winNum], namedColors_[GREEN_X] );
                  break ;
               case 2 :     // Medium Affinity Substrate (e.g. plastic)
                  XSetForeground( dpy_, gcontext[winNum], namedColors_[YELLOW] );
                  break ;
            }  

            // Display the patch
            XFillRectangle( dpy_                 ,
                             dbBuf[winNum]        ,
                            gcontext[winNum]     ,
                            (int) (j*xCelSize_)  ,
                            (int) (i*yCelSize_)  ,
                            (int) xCelSize_ + 1  ,
                            (int) yCelSize_ + 1
                           );
         }
      }
   }
}


//
// Code to display patches from all plates in single window (MPI)
//
void XPlate::displaySubstrateTypesMPI(unsigned winNum)
{
   int         substrate;
   int xPos=0, yPos=0;
   int totArraySize = engParamsPtr_->height_* engParamsPtr_->width_ * numPlatesMPI;
   int plateSize = engParamsPtr_->height_* engParamsPtr_->width_;
   int rankOfPlate = 0;
   int plateDim = sqrt(numPlatesMPI);
   int plateRow = 1;//rankOfPlate/plateDim;
   int width = engParamsPtr_->width_;
   int height = engParamsPtr_->height_;

   for (int i=0; i<totArraySize; i++)
   {
      substrate = allSubstrateValues[i];

      if (i > 1 && i % plateSize == 0)
         rankOfPlate++; 
       
      // Only display if substrate present (>0)  
      if (substrate != 0)
      {     
         switch( substrate )
         {
            case -1 :     // Obstacle
               XSetForeground( dpy_, gcontext[winNum], namedColors_[BLACK]);
               break ;
            case 0 :       // no substate
               break;
            case 1 :  // High Affinity substrate (e.g. concrete)
               XSetForeground( dpy_, gcontext[winNum], namedColors_[GREEN_X] );
               break ;
            case 2 :     // Medium Affinity Substrate (e.g. plastic)
               XSetForeground( dpy_, gcontext[winNum], namedColors_[YELLOW] );
               break ;
         }
         
         for(int j=1; j<=sqrt(numPlatesMPI); j++)   
         {
            if (i < plateSize * (j*sqrt(numPlatesMPI)) )
            {
               plateRow = j;
               break;
            }
         }
      
         xPos = (width * (rankOfPlate - (plateRow-1)*plateDim)) + (i % width);
         yPos = (height * (plateRow-1)) + (i - (plateSize*rankOfPlate)) / width;

         // Display the patch
         XFillRectangle( dpy_                 ,
                          dbBuf[winNum]        ,
                         gcontext[winNum]     ,
                         (int) (xPos*xCelSize_)  ,
                         (int) (yPos*yCelSize_)  ,
                         (int) xCelSize_ + 1  ,
                         (int) yCelSize_ + 1
                        );
      }
   }
}


// ------------------------------ End Of File ------------------------------- //
