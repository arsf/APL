//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------


#include "transformations.h"

#define TRANSFORMDEBUG

#ifndef TRANSFORMDEBUG
   #define DEBUGPRINT(X)  
#else
   #define DEBUGPRINT(X) std::cout<< X <<std::endl;
#endif


//----------------------------------------------------------------------------
// Function to transform the vector V from aircraft coords into ECEF XYZ reference frame
//----------------------------------------------------------------------------
void GetVVinECEFXYZ(blitz::TinyMatrix<double,3,1>* const V,double ECEFXYZ[3],const double lat, const double lon,
                     const double roll,const double pitch,const double heading)
{
   blitz::TinyMatrix<double,3,1> look; //sensor look vectors
   blitz::TinyMatrix<double,3,3> rotmat; //rotation matrix
   blitz::TinyMatrix<double,3,1> retmat; //return matrix

   look=0;
   retmat=0;
   rotmat=0;

   //ref frame rotations
   double rotx=0;
   double roty=-(90+lat);
   double rotz=lon;

   //Get rotation matrix to convert V vector into the actual look direction
   rotmat=Create3DRotMatrix(roll,pitch,heading,RZXY); //RZXY

   //Convert V to look vector
   look=blitz::product(rotmat,*V);

   //get matrix to transform aircraft to ECEF XYZ
   rotmat=Create3DRotMatrix(rotx,roty,rotz,RXZY); //RXYZ

   //Now apply the rot mats to convert the sensor vv to ECEF vv  
   retmat=blitz::product(rotmat,look);

   //Now update the returned ECEF viewvectors
   ECEFXYZ[0]=retmat(0,0);
   ECEFXYZ[1]=retmat(1,0);      
   ECEFXYZ[2]=retmat(2,0);

}


//----------------------------------------------------------------------------
// Function to transform the vector V from aircraft coords into ECEF XYZ reference frame
// Note that that variable names have swapped compared to those that call this function :(
//----------------------------------------------------------------------------
void GetVVinECEFXYZ(blitz::TinyMatrix<double,3,1>* const V,double ECEFXYZ[3],const double lat, const double lon,
                     const double theta,const double phi,const double kappa,const double roll,const double pitch,const double heading)
{
   blitz::TinyMatrix<double,3,1> look; //sensor look vectors
   blitz::TinyMatrix<double,3,3> rotmat; //rotation matrix
   blitz::TinyMatrix<double,3,1> retmat; //return matrix
   blitz::TinyMatrix<double,3,1> sensor; //sensor vector

   look=0;
   retmat=0;
   rotmat=0;
   sensor=0;

   //ref frame rotations - these get from local-level to ECEF
   double rotx=0;
   double roty=-(90+lat);
   double rotz=lon;

   //Sensor look vector
   rotmat=Create3DRotMatrix(theta,phi,kappa,RZXY);//RYZX
   sensor=blitz::product(rotmat,*V);   

   //Get rotation matrix to convert V vector into the actual look direction
   rotmat=Create3DRotMatrix(roll,pitch,heading,RZXY); //RZXY
   //Convert sensor to look vector
   look=blitz::product(rotmat,sensor);

   //get matrix to transform local level to ECEF XYZ
   rotmat=Create3DRotMatrix(rotx,roty,rotz,RXZY); //RXYZ

   //Now apply the rot mats to convert the sensor vv to ECEF vv  
   retmat=blitz::product(rotmat,look);

   //Now update the returned ECEF viewvectors
   ECEFXYZ[0]=retmat(0,0);
   ECEFXYZ[1]=retmat(1,0);      
   ECEFXYZ[2]=retmat(2,0);

}


//----------------------------------------------------------------------------
//Function to return a 3d rotation matrix - expects rx,ry,rz to be in degrees
//This uses CLOCKWISE rotations: 
//    - a positive rx rotates Z towards Y
//    - a positive ry rotates X towards Z
//    - a positive rz rotates Y towards X
//----------------------------------------------------------------------------
blitz::TinyMatrix<double,3,3> Create3DRotMatrix(const double rx, const double ry, const double rz,const rotation_order orderflag)
{
   //Matrix to be returned
   blitz::TinyMatrix<double,3,3> retmat;
   //Temporary matrix to be used for first part of calculation
   blitz::TinyMatrix<double,3,3> tmpmat;
   // Get the sine and cosine of the angles
   //Convert the angles from degrees to radians first
   double degradscale=PI/180.0;
   double cx=cos(rx*degradscale);
   double sx=sin(rx*degradscale);
   double cy=cos(ry*degradscale);
   double sy=sin(ry*degradscale);
   double cz=cos(rz*degradscale);
   double sz=sin(rz*degradscale);

   //Create the 3 rotation matrices
   blitz::TinyMatrix<double,3,3> rotX,rotY,rotZ;

   //Clockwise rotations
   rotX=1,0,0,
         0,cx,-sx,
         0,sx,cx ;

   rotY=cy,0,sy,
         0,1,0,
         -sy,0,cy ;

   rotZ=cz,-sz,0,
         sz,cz,0,
         0,0,1 ;

   //Depending on the order of rotations, generate a 3d rotation matrix
   switch(orderflag)
   {  
   case RXYZ://RxRyRz(V)
      tmpmat=blitz::product(rotX,rotY); //RxRy
      retmat=blitz::product(tmpmat,rotZ); //(RxRy)Rz
      break;
   case RXZY:
      tmpmat=blitz::product(rotX,rotZ);
      retmat=blitz::product(tmpmat,rotY);//(RxRz)Ry
      break;
   case RYXZ:
      tmpmat=blitz::product(rotY,rotX);
      retmat=blitz::product(tmpmat,rotZ);//(RyRx)Rz
      break;
   case RYZX:
      tmpmat=blitz::product(rotY,rotZ);
      retmat=blitz::product(tmpmat,rotX);//(RyRz)Rx
      break;
   case RZXY:
      tmpmat=blitz::product(rotZ,rotX);
      retmat=blitz::product(tmpmat,rotY);//(RzRx)Ry
      break;
   case RZYX:
      tmpmat=blitz::product(rotZ,rotY);
      retmat=blitz::product(tmpmat,rotX);//(RzRy)Rx
      break;
   default:
      std::cout<<"Unknown case in Create3DRotMatrix switch"<<std::endl;
      retmat=0;//set the returned matrix to be all zero
      break;
   }

   return retmat;
}
