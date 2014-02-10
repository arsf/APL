//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef VIEWVECTORS_H
#define VIEWVECTORS_H

#include <string>
#include <cmath>
#include "binfile.h"
#include "commonfunctions.h"

//-------------------------------------------------------------------------
// Class to handle view vectors in sensor/aircraft coordinate system
//-------------------------------------------------------------------------
class ViewVectors
{
   public:
      ViewVectors(std::string fname);//constructor
      ViewVectors(std::string fname,std::string lev1fname);//constructor
      ViewVectors(ViewVectors &ref); //copy constructor
      ~ViewVectors();//destructor

      //Dynamic vars to contain rotation angles about XYZ axes
      //Rotations about X=roll, Y=pitch, Z=yaw
      double* rotX;
      double* rotY;
      double* rotZ; 

      unsigned int NumberItems(){return ccdrows*ccdcols;}

      //Function to apply angular offsets to view vectors and overwrite rotX,Y,Z
      int ApplyAngleRotations(const double rx, const double ry, const double rz);

      //Function to return absolute maximum of the rotX end points
      double AbsMaxX()const{return std::max(fabs(rotX[0]),fabs(rotX[ccdrows*ccdcols-1]));}

      //This function has only been implemented for the python bindings to get around 
      //dereferencing the rotX pointer in python. Can write GetY, GetZ functions as 
      //and when required
      double GetX(const unsigned int i)const
      {
         if(i >= ccdrows*ccdcols)
            throw std::string("Index out of bounds in viewvectors GetX(): ")+ToString(i);
         else
            return rotX[i];
      }

   private:
      //Size of the sensor CCD - i.e. will define the number of view vectors
      unsigned int ccdrows,ccdcols; 
      //BILReader object to handle reading in view vectors
      BinFile* binf; 
      //filename 
      std::string filename;
      //Function to read in the given view vector file
      int ReadVVFile();
};

#endif
