//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

//Geodesic distance, azimuth and zenith calculations
//Bowring and Vincentys methods

#include "geodesics.h"

//-------------------------------------------------------------------------
// Function to get the geodesic distance between two points on the ellipsoid
// surface (if hei1=hei2=0) or if heights given pythogoras of geodesic & heights
// using the Bowring inverse solution (non-iterative). Also returns the
// azimuth and zenith angles from P2 to P1.
// Bowring is suitable (accurate to mm) for points upto ~150km apart
// Assumes input points lon1,lat1 lon2,lat2 are in RADIANS
//-------------------------------------------------------------------------
void GetGeodesicDistance_Bowring(const double lon1,const double lat1,const double hei1,const double lon2,const double lat2,const double hei2,double& distance,double& azimuth,double& zenith,Ellipsoid* ell)
{
   //Use the bowring method for getting azimuth and geodesic distance - this should
   //be ok for distances upto ~150km which is easily good enough for ARSF data.
   //Look at Thomas method or iterative vincenty if need more accurate over longer distances

   //Define a constant ee' based on ellipsoid flattening
   double eep=((2/ell->f())-1) / ((1/ell->f() -1)*(1/ell->f() -1));
   double dphi=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0;
   double w=0,hs=0,S=0;
   //Difference in latitudes (between aircraft and point i)
   dphi=lat1-lat2;
   A=sqrt(1+eep*pow(cos(lat2),4));
   B=sqrt(1+eep*pow(cos(lat2),2));
   C=sqrt(1+eep);
   w=0.5*A*(lon1-lon2);
   D=(dphi/(2*B))*(1+((3*eep*dphi*sin(2*lat2+2*dphi/3))/(4*B*B)));
   E=sin(D)*cos(w);
   F=(1/A)*sin(w)*(B*cos(lat2)*cos(D) - sin(lat2)*sin(D));
   G=atan2(F,E);
   hs=asin(sqrt(E*E+F*F));
   H=atan((1/A)*(sin(lat2) + B*cos(lat2)*tan(D))*tan(w));
   azimuth=(G-H)*180/PI;
   if(azimuth<0)
      azimuth += 360; //to get between 0-360 instead of +/-180
   S=ell->a()*C*2*hs/(B*B);
   zenith=(PI-atan(S/(hei2-hei1)))*180/PI;
   distance=sqrt(S*S+(hei2-hei1)*(hei2-hei1));
}

//-------------------------------------------------------------------------
// Function to get the geodesic distance between two points on the ellipsoid
// surface (if hei1=hei2=0) or if heights given pythogoras of geodesic & heights
// using the Vincenty inverse solution (iterative). Also returns the
// azimuth and zenith angles from P2 to P1.
// Vincenty is suitable (accurate to mm) for all points except ones opposite
// sides of the ellipsoid (as solution may not converge)
// Assumes input points lon1,lat1 lon2,lat2 are in RADIANS
//-------------------------------------------------------------------------
void GetGeodesicDistance_Vincenty(const double lon1,const double lat1,const double hei1,const double lon2,const double lat2,const double hei2,double& distance,double& azimuth,double& zenith,Ellipsoid* ell)
{
   //Use the Vincenty method for calculating geodesics
   //An iterative method but is accurate for all distances
   //may not converge for points at opposite ends of the globe

   double u1=atan((1-ell->f())*tan(lat1));
   double u2=atan((1-ell->f())*tan(lat2));
   double dlat=lat1-lat2;
   double lambda=dlat;
   double sinsigma=0,cossigma=0,sigma=0;
   double sinalpha=0,cossqalpha=0,cos2sigmam=0,C=0;
   double usq=0,A=0,B=0,dsigma=0,S=0;
   double difflambda=lambda,prevlambda=lambda;
   while((difflambda > 0.000000000001)||(difflambda < -0.000000000001))
   {
      sinsigma=sqrt(pow((cos(u2)*sin(lambda)),2) + pow(cos(u1)*sin(u2)-sin(u1)*cos(u2)*cos(lambda),2));
      cossigma=sin(u1)*sin(u2)+cos(u1)*cos(u2)*cos(lambda);
      sigma=atan(sinsigma / cossigma);

      sinalpha=(cos(u1)*cos(u2)*sin(lambda))/ (sinsigma);
      cossqalpha=1-sinalpha*sinalpha;
      cos2sigmam=cossigma - (2*sin(u1)*sin(u2)/cossqalpha);
      C=ell->f()*cossqalpha*(4+ell->f()*(4-3*cossqalpha))/16.0;
      prevlambda=lambda;
      lambda=dlat + (1-C)*ell->f()*sinalpha*(sigma+C*sinsigma*(cos2sigmam + C*cossigma*(-1+2*cos2sigmam*cos2sigmam)));
      difflambda=lambda-prevlambda;
   }

   usq=cossqalpha*(ell->a()*ell->a() - ell->b()*ell->b())/(ell->b()*ell->b());
   A=1+usq*(4096+usq*(-768+usq*(320-175*usq)))/16384.0;
   B=usq*(256+usq*(-128+usq*(74-47*usq)))/1024.0;
   dsigma=B*sinsigma*(cos2sigmam+0.25*B*(cossigma*(-1+2*cos2sigmam*cos2sigmam)) - B*cos2sigmam*(-3+4*sinsigma*sinsigma)*(-3+4*cos2sigmam*cos2sigmam)/6.0 );
   S=ell->b()*A*(sigma-dsigma);
   zenith=(PI-atan(S/(hei2-hei1)))*180/PI;
   distance=sqrt(S*S+(hei2-hei1)*(hei2-hei1));
   azimuth=atan((cos(u2)*sin(lambda))/(cos(u1)*cos(u2)-sin(u1)*cos(u2)*cos(lambda)));
}

//-------------------------------------------------------------------------
// Function to find P2 (on ellipse) given P1 (on ellipse), a distance (in m) 
// and azimuth. Uses the Bowring equations (direct solution). Assumes 
// the input lon1, lat1 are in RADIANS - returns points in DEGREES
//-------------------------------------------------------------------------
void GetDestinationPoint_Bowring(const double lon1,const double lat1,const double distance,const double azimuth,double &lon2,double &lat2,Ellipsoid* ell)
{
   //Use the Bowring method for getting the final destination given a start point, azimuth and distance
   //Define a constant ee' based on ellipsoid flattening
   double eep=((2/ell->f())-1) / ((1/ell->f() -1)*(1/ell->f() -1));
   double A=0,B=0,C=0,D=0;
   double w=0,sigma=0;
   A=sqrt(1+eep*pow(cos(lat1),4));
   B=sqrt(1+eep*pow(cos(lat1),2));
   C=sqrt(1+eep);

   sigma=(distance*B*B)/(ell->a()*C);
   lon2=(lon1+atan((A*tan(sigma)*sin(azimuth))/(B*cos(lat1)-tan(sigma)*sin(lat1)*cos(azimuth)))/A);
   w=0.5*A*(lon2-lon1);
   D=0.5*asin(sin(sigma)*(cos(azimuth)-(sin(lat1)*sin(azimuth)*tan(w))/A));
   lat2=(lat1+2*D*(B-1.5*eep*D*sin(2*lat1+4*B*D/3.0)));

   //now convert to degrees
   lon2=lon2*180/PI;
   lat2=lat2*180/PI;

}

//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
void GetDestinationPoint_Vincenty()
{
   //Use the Vincenty method for getting the final destination given a start point, azimuth and distance
   
}
