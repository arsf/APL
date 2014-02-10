//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#include "planarsurface.h"

#ifndef PLANARDEBUG
   #define DEBUGPRINT(X)  
#else
   #define DEBUGPRINT(X) std::cout<< X <<std::endl;
#endif

//----------------------------------------------------------
//
// PlanarSurface base class methods
//
//----------------------------------------------------------


//-------------------------------------------------------------------------
//Create a plane equation from the 3 given points (all to be in the plane)
//Input points are in ECEF XYZ (expect pointers to arrays of size 3 elements)
//-------------------------------------------------------------------------
PlanarSurface::PlanarSurface(const double* const p1,const double* const p2,const double* const p3)
{

   //Set one of the points to p (p is the point in the plane eqn)
   this->px=p1[0];
   this->py=p1[1];
   this->pz=p1[2];

   //Create 2 vectors in the plane from the 3 given points
   //P2 - P1
   double v1x=p2[0]-p1[0];
   double v1y=p2[1]-p1[1];
   double v1z=p2[2]-p1[2];
   //P3 - P1
   double v2x=p3[0]-p1[0];
   double v2y=p3[1]-p1[1];
   double v2z=p3[2]-p1[2];

   //Calculate the cross product to get normal vector to the plane
   this->nx=(v1y*v2z)-(v1z*v2y);
   this->ny=-(v2z*v1x)+(v2x*v1z);
   this->nz=(v1x*v2y)-(v1y*v2x);

   //Normailise it to get unit normal for easier calculations
   double magn=sqrt(nx*nx+ny*ny+nz*nz);
   this->nx=nx/magn;
   this->ny=ny/magn;
   this->nz=nz/magn;

   //Assign up to 0
   this->ux=0;
   this->uy=0;   
   this->uz=0;
}

//-------------------------------------------------------------------------
// Function to assign the local up vector (ECEF XYZ)
//-------------------------------------------------------------------------
void PlanarSurface::AssignLocalUp(const double* const up)
{
   //Are there any checks I can do on the up vector?
   ux=up[0];
   uy=up[1];
   uz=up[2];
   //Force it to be a normal vector
   double magn=sqrt(ux*ux+uy*uy+uz*uz);
   if(magn!=1)
   {
      ux=ux/magn;
      uy=uy/magn;
      uz=uz/magn;
   }
   //We can check the normal vector against up here and change the sign
   //if the direction of the normal is down. If n.u < 0 then n is below horizon
   double dp=ux*nx + uy*ny + uz*nz;
   if(dp<0)
   {
      nx=-nx;
      ny=-ny;
      nz=-nz;
   }
   DEBUGPRINT(ux<<" "<<uy<<" "<<uz<<" "<<dp)
}

//-------------------------------------------------------------------------
// Function to calculate the azimuth angle of the plane: 0 = North facing
//-------------------------------------------------------------------------
double PlanarSurface::CalculateAzimuth(const double* const North)
{
   //Check that the north vector is normal
   double magn=sqrt(North[0]*North[0] + North[1]*North[1] +North[2]*North[2]);
   double north[3]={0};
   if(magn!=1)
   {
      north[0]=North[0]/magn;
      north[1]=North[1]/magn;
      north[2]=North[2]/magn;
   }

   //Get the vector perpendicular to n and up, call it P
   double Px=ny*uz - nz*uy;
   double Py=-(nx*uz - nz*ux);
   double Pz=nx*uy - ny*ux;

   //Now get the normal vector projected in the lat/lon plane (plane perpendicular to up)
   //This is the solution to 3 eqns:
   //np . P =0
   //np . u =0
   //np . n =cos(theta)

   //Get the slope of the plane and convert it to the angle from horizon up
   double theta=PI/2.0 - CalculateSlope();

   //Define this constant to make eqns "nicer"
   double A=-(uy - (ux*Py/Px)) / (uz-(ux*Pz/Px));

   //np is the normal vector projected into the lat/lon plane
   double npy=cos(theta) / (A*nz + ny - A*nx*Pz/Px - Py*nx/Px);
   double npz=A*npy;
   double npx=(-Pz*npz - Py*npy)/Px;

   //Now we need the dot product between this and North vector to get azimuth angle
   double cos_az=north[0]*npx + north[1]*npy + north[2]*npz;
   //And the cross product
   double sin_az=sqrt((north[1]*npz - north[2]*npy)*(north[1]*npz - north[2]*npy) 
               + (north[0]*npz - north[2]*npx)*(north[0]*npz - north[2]*npx)
                + (north[0]*npy - north[1]*npx)*(north[0]*npy - north[1]*npx));
   DEBUGPRINT("Up: "<<ux<<" "<<uy<<" "<<uz)   
   DEBUGPRINT("North: "<<north[0]<<" "<<north[1]<<" "<<north[2]) 
   DEBUGPRINT("normal: "<<nx<<" "<<ny<<" "<<nz)     
   DEBUGPRINT("normal_projected: "<<npx<<" "<<npy<<" "<<npz)   
   DEBUGPRINT("Cosine: "<<cos_az<<" Sine: "<<sin_az)
   //To return the azimuth (between 0 and 2PI)
   //since atan2 is in range -Pi to Pi (with +ve counterclockwise)
   if(atan2(sin_az,cos_az)>0)
      return -(atan2(sin_az,cos_az)-2*PI);
   else
      return -(atan2(sin_az,cos_az));
}

//-------------------------------------------------------------------------
// Function to calculate the slope of the plane: 0 = flat | 90 = vertical
//-------------------------------------------------------------------------
double PlanarSurface::CalculateSlope()
{
   //Cos angle is the dot product between the normal and up vector (both are unit vectors)
   //This returns the slope in radians
   double theta=acos(ux*nx+uy*ny+uz*nz);
   //We want to return the slope between 0 and Pi/2 (the accute angle between vectors)
   if(theta < PI/2.0)
      return theta;
   else
      return PI - theta;
}

//----------------------------------------------------------
//
// TriangularPlane Methods
//
//----------------------------------------------------------

//-------------------------------------------------------------------------
//Function to find if 2 points x and y are on the same side of a line defined by a and b 
//-------------------------------------------------------------------------
bool TriangularPlane::SameSide(const double* const a, const double* const b,const double* const x,const double* const y)
{
   //Need to cross product (b-a)^(b-x) and (b-a)^(b-y)
   //Each vector has 3 elements
   double cp1=((b[1]-a[1])*(b[2]-x[2]) - (b[2]-a[2])*(b[1]-x[1])
              -(b[0]-a[0])*(b[2]-x[2]) + (b[2]-a[2])*(b[0]-x[0])
              +(b[0]-a[0])*(b[1]-x[1]) - (b[1]-a[1])*(b[0]-x[0]));

   double cp2=((b[1]-a[1])*(b[2]-y[2]) - (b[2]-a[2])*(b[1]-y[1])
              -(b[0]-a[0])*(b[2]-y[2]) + (b[2]-a[2])*(b[0]-y[0])
              +(b[0]-a[0])*(b[1]-y[1]) - (b[1]-a[1])*(b[0]-y[0]));

   DEBUGPRINT(cp1<<" "<<cp2<<" "<<cp1*cp2<<" "<<floor(1000*(cp1 * cp2)))

   //If both are +ve or both -ve then x,y are on the same side
   //14/01/2011 - added the ceil(()-0.9999) to try and account for when intersect is "on the line" of plane boundary (and very small -ve) so we'll accept it.
   if(ceil((cp1 * cp2)-0.9999) >= 0)
      return true;
   else
      return false;
}

//-------------------------------------------------------------------------
//Function to check if a point x is within a triangle defined by a, b and c 
//-------------------------------------------------------------------------
bool TriangularPlane::InTriangle(const double* const a, const double* const b,const double* const c,const double* const x)
{
   //If the point is on the same side of each triangle edge as the other point then it is within the triangle
   if(SameSide(a,b,c,x) && SameSide(b,c,a,x) && SameSide(c,a,b,x))
      return true;
   else
      return false;
}

//-------------------------------------------------------------------------
//Function to check if a point is within a triangle using barycentric coords
//-------------------------------------------------------------------------
bool TriangularPlane::Barycentric(const double* const a, const double* const b,const double* const c,const double* const x)
{
   //We can describe a point as being A+u*(B-A)+v*(C-A)
   //We need 3 vectors: 2 sides of the triangle (b-a, c-a) and vector to the point (x-a)
   double ab[3]={b[0]-a[0],b[1]-a[1],b[2]-a[2]};
   double ac[3]={c[0]-a[0],c[1]-a[1],c[2]-a[2]};
   double ax[3]={x[0]-a[0],x[1]-a[1],x[2]-a[2]};
   //Now calculate dot products and magnitude square
   double magsq_ab=ab[0]*ab[0]+ab[1]*ab[1]+ab[2]*ab[2];
   double magsq_ac=ac[0]*ac[0]+ac[1]*ac[1]+ac[2]*ac[2];
   double abdotac=ab[0]*ac[0]+ab[1]*ac[1]+ab[2]*ac[2];
   double abdotax=ab[0]*ax[0]+ab[1]*ax[1]+ab[2]*ax[2];
   double acdotax=ac[0]*ax[0]+ac[1]*ax[1]+ac[2]*ax[2];
   //Now calculate u and v
   double denom=magsq_ab*magsq_ac-abdotac*abdotac;
   double u=(magsq_ab*acdotax - abdotac*abdotax)/denom;
   double v=(magsq_ac*abdotax-abdotac*acdotax)/denom;
   //We will allow points that fall on the boundary line to be classed as intersecting
   //introduced zero variable which is -1x10^(-8) to try and solve
   //floating point error when intersect falls just short of 0 and therefore is not counted as intersecting
   double zero=-0.00000001;
   return ((u>=zero)&&(v>=zero)&&(u+v<=1));
}

//-------------------------------------------------------------------------
//Function to calculate if a vector and triangular plane intersect within 
//the bounds defined by the triangle. 
//Inputs are:vx,vy,vz - 2 points which vector passes through
//             px,py,pz - returned point (intersect) 
//-------------------------------------------------------------------------
bool TriangularPlane::Intersect(const double* const vX, const double* const vY, const double* const vZ, double* const px, double* const py, double* const pz)
{
   //Calculate the line scalar to find intersect point
   double t=0;
   double numer=this->nx*(this->px-vX[0])+this->ny*(this->py-vY[0])+this->nz*(this->pz-vZ[0]); //n . (np - lp1)
   double denom=this->nx*(vX[1]-vX[0])+this->ny*(vY[1]-vY[0])+this->nz*(vZ[1]-vZ[0]);//n . (lp2 - lp1)

   if(denom !=0) //if line is not parallel to plane
      t=(numer)/(denom);
   else
      return false;

   //Find intersect point
   *px=vX[0]+t*(vX[1]-vX[0]);
   *py=vY[0]+t*(vY[1]-vY[0]);
   *pz=vZ[0]+t*(vZ[1]-vZ[0]);

   double PXYZ[3]={*px,*py,*pz};
   DEBUGPRINT("Intersect point XYZ: "<<px[0]<<" "<<py[0]<<" "<<pz[0]);

   //Is point within given bounds 
   //test if the point falls within the given bounds
   //if(this->InTriangle(this->P1,this->P2,this->P3,PXYZ))
   if(this->Barycentric(this->P1,this->P2,this->P3,PXYZ))
   {
      //Yes it does - return true
      return true;
   }
   else
   {
      //No it doesn't - return false
      return false;
   }

}
