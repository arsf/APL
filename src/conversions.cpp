//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#include "conversions.h"

//-------------------------------------------------------------------------
//Constructor for ellipsoid class. Pass a and b radius values to initiate.
//-------------------------------------------------------------------------
Ellipsoid::Ellipsoid(const double a, const double b,const double f)
{
   this->A=a;  //major radius
   this->B=b;  //minor radius
   this->EE=(1-(this->B*this->B)/(this->A*this->A));  //eccentricty squared
   this->F=f; //flattening (or 1 over flattening 1/f)
   this->md1=0;
   this->md2=0;
   this->md3=0;
}

//-------------------------------------------------------------------------
//Constructor for ellipsoid class. Pass ellipsoid model keyword to initiate.
//-------------------------------------------------------------------------
Ellipsoid::Ellipsoid(ELIPMODEL ellmodel)
{
   switch(ellmodel)
   {
   case WGS84:
      this->A=6378137; //major radius
      this->B=6356752.3142; //minor radius
      this->EE=(1-(this->B*this->B)/(this->A*this->A));  //eccentricty squared
      this->F=(this->A - this->B)/this->A; //flattening
      this->name="WGS-84";
      //These are the first 3 terms in the 1 degree meridional distance approximations from 
      //the delambre series expansion of the meridional distance integral for a distance of 1 degree 
      //at latitude in the centre of this degree. eg the distance between lat-0.5 and lat+0.5 
      this->md1=111132.954;
      this->md2=559.822;
      this->md3=1.175;
      break;
   default:
      throw "Ellipsoid model not recognised";
      break;
   }

}

//-------------------------------------------------------------------------
//Convert Geodetic or Geographic Lat/Lon/Hei to ECEF XYZ
//Pass pointers to arrays containing the Lat/Lon/hei data
//and also pointers to the storage arrays of XYZ data
//and the number of points to be converted (length of arrays).
//Type should be either GEOGRAPHIC or GEODETIC
//-------------------------------------------------------------------------
int ConvertLLH2XYZ(const double* const Lat, const double* const Lon, const double* const Hei, double* const X, double* const Y, double* const Z, const unsigned int npoints, unsigned int type,Ellipsoid* const ell)
{
   double a=ell->a();
   //double b=ell->b();
   double ee=ell->ee(); //eccentricty squared

   //Test if geographic or geodetic conversion is required
   switch(type) 
   {
   case GEODETIC:
      {//Need braces because of declarations of vars in the specific case
         //Convert the geodetic latitude, longitude and height to ECEF cartesians
         double C=0; //denominator for the normal i.e. distance from Z to surface along ellipsoid normal
         double lorad=0, larad=0;
         for(unsigned int i=0;i<npoints;i++)
         {
            //Need Lat/Lon/Hei in radians here so convert from degrees
            larad=Lat[i]*PI/180.0;
            lorad=Lon[i]*PI/180.0;
            C=sqrt(1-ee*sin(larad)*sin(larad));
            X[i]=(a/C + Hei[i])*cos(larad)*cos(lorad);
            Y[i]=(a/C + Hei[i])*cos(larad)*sin(lorad);
            Z[i]=((a*(1-ee))/C + Hei[i])*sin(larad);
         }
      }
      break;
   case GEOGRAPHIC:
      //Convert the geographic latitude, longitude and height to ECEF cartesians
      throw "This GEOGRAPHIC conversion method has not been implemented";
      break;
   default:
      // "Unrecognised conversion type in ConvertLLH2XYZ"
      return -1; //error flag
   }

   return 1;
}

//-------------------------------------------------------------------------
//Convert ECEF XYZ to Geodetic or Geographic Lat/Lon/Hei
//-------------------------------------------------------------------------
int ConvertXYZ2LLH(const double* const X, const double* const Y, const double* const Z, double* const Lat, double* const Lon, double* const Hei, const unsigned int npoints, unsigned int type,Ellipsoid* const ell)
{
   double a=ell->a();
   double b=ell->b();
   double ee=ell->ee(); //eccentricty squared

   //Test if geocentric or geodetic conversion is required
   switch(type) 
   {
   case GEODETIC:
      {//Need braces because of declarations of vars in the specific case
         //Convert the ECEF cartesians into geodetic latitude, longitude and height
         double bb=b*b; //b squared
         double aa=a*a; //a squared
         double eeprime=((aa)/(bb)-1); //second eccentricty squared
         double EE=aa-bb;

         double rr=0,F=0,G=0,C=0,S=0,P=0,Q=0,r0=0,U=0,V=0,Z0=0;
         
         for(unsigned int i=0;i<npoints;i++)
         {
            rr=(X[i]*X[i]+Y[i]*Y[i]);
            F=54*bb*Z[i]*Z[i];
            G=rr+(1-ee)*Z[i]*Z[i]-ee*EE;
            C=(ee*ee*F*rr)/(G*G*G);
            S=pow((1+C+sqrt(C*C+2*C)),(1/3.0));
            P=(F)/(3*G*G*(S+1+1/S)*(S+1+1/S));
            Q=sqrt(1+2*ee*ee*P);
            r0=((-P*ee*sqrt(rr))/(1+Q)) + sqrt(0.5*aa*(1+1/Q) - (P*(1-ee)*Z[i]*Z[i])/(Q*(1+Q)) -0.5*P*rr );
            U=sqrt( (sqrt(rr)-ee*r0)*(sqrt(rr)-ee*r0) + Z[i]*Z[i]);
            V=sqrt( (sqrt(rr)-ee*r0)*(sqrt(rr)-ee*r0) +(1-ee)*Z[i]*Z[i] );
            Z0=(bb*Z[i])/(a*V);
            Hei[i]=U*(1-(bb/(a*V)));
            Lat[i]=atan((Z[i] + eeprime*Z0)/(sqrt(rr)));
            Lon[i]=atan2(Y[i],X[i]);
         }
      }
      break;
   case GEOGRAPHIC:
      //Convert the ECEF cartesians into geographic latitude, longitude and height 
      throw "This GEOGRAPHIC conversion method has not bee implemented";
      break;
   default:
      // "Unrecognised conversion type in ConvertXYZ2LLH"
      return -1; //error flag
   }

   return 1;
}



