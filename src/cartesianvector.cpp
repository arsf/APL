//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#include "cartesianvector.h"

#ifndef CARTDEBUG
   #define DEBUGPRINT(X)  
#else
   #define DEBUGPRINT(X) std::cout<< X <<std::endl;
#endif

//-------------------------------------------------------------------------
// Constructor for vectors, using the number of vectors as input
//-------------------------------------------------------------------------
CartesianVector::CartesianVector(const unsigned int nvecs)
{
   if(nvecs != 0)
   {
      //Assign the number of vectors and create arrays of this size
      numberofvectors=nvecs;
      X=new double[nvecs];
      Y=new double[nvecs];
      Z=new double[nvecs];

      //Initialise the vectors to 0
      for(unsigned int i=0;i<nvecs;i++)
      {
         X[i]=0;
         Y[i]=0;
         Z[i]=0;
      }
   }
   else
      throw "There must be at least 1 element in the CartesianVector.";

   //Set the origin of the system to 0
   this->originX=0;
   this->originY=0;
   this->originZ=0;
}

//-------------------------------------------------------------------------
// Constructor for vectors, using the number of vectors as input and the origin
//-------------------------------------------------------------------------
CartesianVector::CartesianVector(const unsigned int nvecs,const double oX,const double oY,const double oZ)
{
   if(nvecs != 0)
   {
      //Assign the number of vectors and create arrays of this size
      numberofvectors=nvecs;
      X=new double[nvecs];
      Y=new double[nvecs];
      Z=new double[nvecs];

      //Initialise the vectors to 0
      for(unsigned int i=0;i<nvecs;i++)
      {
         X[i]=0;
         Y[i]=0;
         Z[i]=0;
      }
   }
   else
      throw "There must be at least 1 element in the CartesianVector.";

   //Set the origin of the system to that given
   this->originX=oX;
   this->originY=oY;
   this->originZ=oZ;
}

//-------------------------------------------------------------------------
// Destructor to destroy the allocated memory
//-------------------------------------------------------------------------
CartesianVector::~CartesianVector()
{
   //If the arrays still exist then they must be destroyed
   if(X!=NULL)
      delete[] X;
   if(Y!=NULL)
      delete[] Y;
   if(Z!=NULL)
      delete[] Z;
}

//-------------------------------------------------------------------------
// Function to return the index of the view vector closest to nadir pointing
//-------------------------------------------------------------------------
unsigned int CartesianVector::GetNadirIndex(blitz::TinyMatrix<double,3,1> nadir)
{
   //This function is called to find the view vector with the 
   //closest angle to nadir pointing and return its index

   //For each view vector calculate dot product
   //between it and nadir and find one with smallest angle
   //vectors are in sequential order so will be (semi)"monotonic"
   double dotp=0,grad=1;
   double dotpm1=this->X[0]*nadir(0,0)+this->Y[0]*nadir(1,0)+this->Z[0]*nadir(2,0);
   unsigned int pixel=0; //counter
   
   //Grad should always start +ve unless pixel 0 is closest to nadir
   //when it turns -ve we know we have reached the maximum (at previous pixel)
   while((pixel<this->numberofvectors)&&(grad>0))
   {
      //increase counter at start since pixel 0 done above
      pixel++; 

      //Get dot product of vv and nadir
      dotp=this->X[pixel]*nadir(0,0)+this->Y[pixel]*nadir(1,0)+this->Z[pixel]*nadir(2,0);

      //calculate gradient 
      grad=dotp - dotpm1;

      //Set up for next iteration
      dotpm1=dotp;
   }

   return (pixel-1);
}

