//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef CONVERSIONS_H
#define CONVERSIONS_H

#include <cmath>
#include <string>

//-------------------------------------------------------------------------
//This is a file of the simple coordinate conversions used within 
//the applcorr geocorrection software
//-------------------------------------------------------------------------

#define GEODETIC 0
#define GEOGRAPHIC 1

//-------------------------------------------------------------------------
//Ellipsoid models supported in geolocation software / ellipsoid class
//-------------------------------------------------------------------------
enum ELIPMODEL{WGS84};

//-------------------------------------------------------------------------
//Ellipsoid class
//-------------------------------------------------------------------------
class Ellipsoid{
   public:
      Ellipsoid(const double a, const double b,const double f);
      Ellipsoid(ELIPMODEL ellmodel);
      const double a(){return A;}
      const double b(){return B;}
      const double f(){return F;}
      const double ee(){return EE;}
      const std::string Name(){return name;}
      //Approximation for a 1 degree meridional distance (lat-0.5 to lat+0.5)
      const double meridional_degree(const double lat){return md1 - md2*(cos(2*lat)) + md3*(cos(4*lat));} 
   private:
      double A,B,F,EE;
      double md1,md2,md3;
      std::string name;
};

//-------------------------------------------------------------------------
//Global Constants
//-------------------------------------------------------------------------
#ifndef PI_
#define PI_
const double PI=atan(1.0)*4;
#endif

//-------------------------------------------------------------------------
//Convert Geodetic or Geocentric Lat/Lon/Hei to ECEF XYZ
//-------------------------------------------------------------------------
int ConvertLLH2XYZ(const double* const Lat, const double* const Lon, const double* const Hei, double* const X, double* const Y, double* const Z, const unsigned int npoints, unsigned int type,Ellipsoid* const ell,const int baddatavalue=-9999);

//-------------------------------------------------------------------------
//Convert ECEF XYZ to Geodetic or Geocentric Lat/Lon/Hei
//-------------------------------------------------------------------------
int ConvertXYZ2LLH(const double* const X, const double* const Y, const double* const Z, double* const Lat, double* const Lon, double* const Hei, const unsigned int npoints, unsigned int type,Ellipsoid* const ell,const int baddatavalue=-9999999);


#endif
