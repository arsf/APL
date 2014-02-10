//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef CARTESIANVECTOR_H
#define CARTESIANVECTOR_H

#include <blitz/array.h>
#include "conversions.h"

//#define CARTDEBUG

//Class to hold cartesian vectors to make passing between functions easier
//No safety checks are applied when inserting/retrieving data from X,Y,Z
//so as not to increase the processing overheads too much.
//If its decided to be made safer then make X,Y,Z private and have functions
//to Get/Set the data that also check numberofvectors and element number

class CartesianVector
{
   public:
      //Constructor with number of vectors
      CartesianVector(const unsigned int nvecs);

      //Constructor with number of vectors and origin location
      CartesianVector(const unsigned int nvecs,const double oX,const double oY,const double oZ);
      //Default destructor
      ~CartesianVector();

      //arrays that hold the xyz data
      double* X;
      double* Y;
      double* Z;

      //Return the index of the vector which is closest to nadir pointing
      unsigned int GetNadirIndex(blitz::TinyMatrix<double,3,1> nadir);

      //Return the origin 
      double const OriginX(){return originX;}
      double const OriginY(){return originY;}
      double const OriginZ(){return originZ;}

   private:
      //Number of vectors for X,Y,Z arrays
      unsigned int numberofvectors;

      //origin of vectors
      double originX,originY,originZ; 
};


#endif
