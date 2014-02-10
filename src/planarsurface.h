//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef PLANARSURFACE
#define PLANARSURFACE

//#define PLANARDEBUG

#include <cmath>
#include <iostream>

#ifndef PI_
#define PI_
const double PI=4*atan(1.0);
#endif

//-------------------------------------------------------------------------
// Class to describe a planar surface - essentially a plane equation
//-------------------------------------------------------------------------
class PlanarSurface
{
public:
   PlanarSurface(const double* const p1,const double* const p2,const double* const p3);
   //Function to assign the local up vector (array of length 3)
   void AssignLocalUp(const double* const up);
   //Function to calculate the slope angle of the plane
   double CalculateSlope();
   //Function to calculate the azimuth angle of the plane
   double CalculateAzimuth(const double* const north);

protected:
   //Normal vector
   double nx,ny,nz;
   //Point in plane
   double px,py,pz;
   //Local up vector
   double ux,uy,uz;
};


//-------------------------------------------------------------------------
// A type of planar surface defined by, and whose extent is, a triangle
//-------------------------------------------------------------------------
class TriangularPlane:public PlanarSurface
{
public:
   TriangularPlane(const double* const p1,const double* const p2,const double* const p3):PlanarSurface(p1, p2, p3)
   {
      //Set the 3 points to the given point data
      P1[0]=p1[0];
      P1[1]=p1[1];
      P1[2]=p1[2];
      P2[0]=p2[0];
      P2[1]=p2[1];
      P2[2]=p2[2];
      P3[0]=p3[0];
      P3[1]=p3[1];
      P3[2]=p3[2];
   }

   //Calculate where the given vector intersects the plane and return true if it is within the triangular bounds
   //vX,vY,vZ are arrays of 2 elements, containing data of 2 points which the vector passes through
   //px,py,pz is the returned intersect point
   bool Intersect(const double* const vX, const double* const vY, const double* const vZ, double* const px, double* const py, double* const pz);

   const double* Point1()const{return P1;}
   const double* Point2()const{return P2;}
   const double* Point3()const{return P3;}


private:
   double P1[3],P2[3],P3[3];
   //Are points x and y on the same side of the line defined by a and b [these are 2D points]
   bool SameSide(const double* const a, const double* const b,const double* const x,const double* const y);
   //Check if a point x is within a triangle defined by a, b and c [these are 2D points]
   bool InTriangle(const double* const a, const double* const b,const double* const c,const double* const x);
   //Check if a point x is in triangle a,b,c 
   bool Barycentric(const double* const a, const double* const b,const double* const c,const double* const x);
};

#endif
