//----------------------------------------------------------------------------//
// myException.hpp                                                            //
//----------------------------------------------------------------------------//

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

#ifndef MYEXCEPTION_HPP
#define MYEXCEPTION_HPP


//----------------------------------------------------------------------------//
// Libraries                                                                  //
//----------------------------------------------------------------------------//

#include <exception>
using std::exception ;

#include <string>
using std::string ;

#include <sstream>
using std::ostringstream ;


//----------------------------------------------------------------------------//
// myException                                                                //
//----------------------------------------------------------------------------//

class myException : public exception
{
   public :

      myException( const string       &inMessage    ,
                   const string       &inFileName   ,
                   const unsigned      inLineNumber
                 )
      : exception(), message_(inMessage), fileName_(inFileName),
      lineNumber_(inLineNumber)
      { }

      virtual ~myException() throw()
      { }

      const char * what () const throw()
      {
         ostringstream whatMessage ;

         whatMessage << "in " << fileName_ << " (" << lineNumber_
         << "): " << message_ ;

         return whatMessage.str().c_str() ;
      }


   protected :

      string      message_    ,
                  fileName_   ;
      unsigned    lineNumber_ ;
};


//----------------------------------------------------------------------------//

#endif                        // MYEXCEPTION_HPP


//------------------------------- End Of File --------------------------------//
