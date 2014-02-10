//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#include "leverbore.h"

#ifndef LEVERDEBUG
   #define DEBUGPRINT(X)  
#else
   #define DEBUGPRINT(X) std::cout<< X <<std::endl;
#endif

//-------------------------------------------------------------------------
//Boresight constructor
//-------------------------------------------------------------------------
Boresight::Boresight(const double r, const double p, const double h)
{
   this->roll=r;
   this->pitch=p;
   this->heading=h;
}

//-------------------------------------------------------------------------
// Apply the boresight values in-place onto the given roll, pitch, heading
//-------------------------------------------------------------------------
void Boresight::ApplyBoresight(double* const roll,double* const pitch,double* const heading)
{
   *roll+=this->roll;
   *pitch+=this->pitch;
   *heading+=this->heading;
}


//-------------------------------------------------------------------------
//Leverarm constructor
//-------------------------------------------------------------------------
Leverarm::Leverarm(const double lx, const double ly, const double lz)
{
   this->x=lx;
   this->y=ly;
   this->z=lz;

   //Create a blitz vector of the lever arm values (for use in conversions)
   arm=x,y,z;
}

//-------------------------------------------------------------------------
// Function to rotate the aircraft leverarm (x,y,z) to ECEF (tempX,tempY.tempZ)
//-------------------------------------------------------------------------
void Leverarm::LeverToECEF(const double roll, const double pitch, const double heading,double* const GPSlat, double* const GPSlon)
{
   //Apply a 2-part transformation to the lever arms
   //1 ... Rotate by aircraft attitude to get into a local level reference frame
   //2 ... Rotate by lat/lon to get into ECEF reference frame

   //Variabe to hold the result of the transformation
   double ECEFXYZ[3]={0,0,0};

   //Call the transformation function - this does both parts of the 2-part transformation
   GetVVinECEFXYZ(&arm,ECEFXYZ,*GPSlat,*GPSlon,roll,pitch,heading);

   //Store the ECEFXYZ in the temp variables
   tempX=ECEFXYZ[0];
   tempY=ECEFXYZ[1];
   tempZ=ECEFXYZ[2];

   DEBUGPRINT("ECEF lever arm values:"<<tempX<<" "<<tempY<<" "<<tempZ<<" Magnitude: "<<sqrt(tempX*tempX+tempY*tempY+tempZ*tempZ))
}

//-------------------------------------------------------------------------
// Function to add the ECEF lever arms onto the GPS
//-------------------------------------------------------------------------
void Leverarm::AddToGPS(double* const GPSlat, double* const GPSlon, double* const GPShei,Ellipsoid* ellipsoid)
{
   //Will need to convert the GPS from lat/lon to ECEF,
   //add on the ECEF leverarm values and then convert back
   double ECEFGPS[3]={0,0,0};
   double newgps[3]={0,0,0};
   
   //Convert GPS to ECEF cartesians
   ConvertLLH2XYZ(GPSlat,GPSlon,GPShei,&ECEFGPS[0],&ECEFGPS[1],&ECEFGPS[2],1,GEODETIC,ellipsoid);

   //Add on the ECEF lever arm
   ECEFGPS[0]+=this->tempX;
   ECEFGPS[1]+=this->tempY;
   ECEFGPS[2]+=this->tempZ;

   //Convert ECEF coordinates back to lat/lon/hei
   ConvertXYZ2LLH(&ECEFGPS[0],&ECEFGPS[1],&ECEFGPS[2],&newgps[0],&newgps[1],&newgps[2],1,GEODETIC,ellipsoid);

   //Set the new point to overwrite the old one
   *GPSlat=newgps[0]*180/PI;
   *GPSlon=newgps[1]*180/PI;
   *GPShei=newgps[2];

}

//-------------------------------------------------------------------------
// Function to transform the lever arm and add it onto the GPS position
//-------------------------------------------------------------------------
void Leverarm::ApplyLeverArm(const double roll, const double pitch, const double heading,
                   double* const GPSlat, double* const GPSlon, double* const GPShei)
{
   //Need to convert the lever arm from aircraft coords and add it onto
   //the GPS position, returning the values in place
   
   //Convert the lever arm to ECEF
   LeverToECEF(roll,pitch,heading,GPSlat,GPSlon);

   //Create an ellipsoid model to use in the coordinate transformations
   Ellipsoid ellipsoid(WGS84);

   //Add on to the GPS - this currently OVERWRITES the GPSlat/lon/hei values
   AddToGPS(GPSlat,GPSlon,GPShei,&ellipsoid);
}
