//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef LEVERBORE_H
#define LEVERBORE_H

#include "conversions.h"
#include "transformations.h"

//#define LEVERDEBUG

//----------------------------------------------------------------
// Borsight class to hold boresight values
//----------------------------------------------------------------
class Boresight
{
   public:
      //Constructor 
      Boresight(const double r, const double p, const double h);

      //Functions to return private member data
      double Roll() const {return roll;}
      double Pitch() const {return pitch;}
      double Heading() const {return heading;}

      //Function to apply the boresight values to the given attitude values
      void ApplyBoresight(double* const roll,double* const pitch,double* const heading);

   private:
      double roll,pitch,heading;
};



//----------------------------------------------------------------
// Lever arm class to hold lever arm values
//----------------------------------------------------------------
class Leverarm
{
   public:
      //Constructor
      Leverarm(const double lx, const double ly, const double lz);

      //Function to transform the lever arm and add it onto the GPS position
      void ApplyLeverArm(const double roll, const double pitch, const double heading,
                         double* const GPSlat, double* const GPSlon, double* const GPShei);

      //Functions to return private member data
      double X() const {return x;}
      double Y() const {return y;}
      double Z() const {return z;}

   private:
      //Function to add the ECEF XYZ lever arm coords onto the GPS
      void AddToGPS(double* const GPSlat, double* const GPSlon, double* const GPShei,Ellipsoid* ellipsoid);

      //Function to rotate the aircraft leverarm (x,y,z) to ECEF (tempX,tempY.tempZ)
      void LeverToECEF(const double roll, const double pitch, const double heading,double* const GPSlat, double* const GPSlon);
      
      //These hold the original lever arm values - they should not be changed
      double x,y,z; 

      //These hold ECEF converted lever arm values - considered temporary and private
      double tempX,tempY,tempZ; 

      //This is a blitz vector of the lever arm (for transformation functions)
      blitz::TinyMatrix<double,3,1> arm;
   
};


#endif
