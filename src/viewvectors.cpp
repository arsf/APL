//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#include "viewvectors.h"

//-------------------------------------------------------------------------
// Constructor from view vector filename
//-------------------------------------------------------------------------
ViewVectors::ViewVectors(std::string fname)
{
   //Need to create the bil reader and assign ccd size
   //And read in the view vectors from the file
   this->filename=fname;
   this->binf=new BinFile(fname);

   //Get the size of the ccd (the part that we're interested in) from
   //the bil header file
   this->ccdrows=StringToUINT(binf->FromHeader("samples"));
   this->ccdcols=StringToUINT(binf->FromHeader("lines"));

   //Create the dynamic variables
   this->rotX=new double[ccdrows*ccdcols];
   this->rotY=new double[ccdrows*ccdcols];
   this->rotZ=new double[ccdrows*ccdcols];

   //Populate the rotX,rotY,rotZ variables - read in the view vector file
   this->ReadVVFile();
}

//-------------------------------------------------------------------------
// Constructor from view vector filename and level 1 data
//-------------------------------------------------------------------------
ViewVectors::ViewVectors(std::string fname,std::string lev1fname)
{
   //Need to create the bil reader and assign ccd size
   //And read in the view vectors from the file
   this->filename=fname;
   this->binf=new BinFile(fname);

   //Get the size of the ccd from the bil header file
   this->ccdrows=StringToUINT(binf->FromHeader("samples"));
   this->ccdcols=StringToUINT(binf->FromHeader("lines"));

   //Create the dynamic variables
   this->rotX=new double[ccdrows*ccdcols];
   this->rotY=new double[ccdrows*ccdcols];
   this->rotZ=new double[ccdrows*ccdcols];

   //Populate the rotX,rotY,rotZ variables - read in the view vector file
   this->ReadVVFile();

   //Now we need to trim off any unwanted pixels (e.g. fodis) - get this info from lev1 file
   //And apply spatial binning to the data
   BinFile lev1(lev1fname);
   unsigned int l1samps=StringToUINT(lev1.FromHeader("samples"));
   unsigned int xstart=StringToUINT(lev1.FromHeader("x start"));
   unsigned int spatbin=StringToUINT(lev1.FromHeader("binning",1));
   if(spatbin==0)
      throw "Spatial binning is missing from the level-1 file header. Please add a line in the header file containing: binning = VALUE where VALUE is the correct spatial and spectral binning of the data in the form, e.g., {1,1}";

   lev1.Close();
   //The binned vector row size
   unsigned int newsizerows=((ccdrows)/spatbin);

   //Temporary arrays to fill with binned vvs
   double* tmpX=new double[newsizerows*ccdcols];
   double* tmpY=new double[newsizerows*ccdcols];
   double* tmpZ=new double[newsizerows*ccdcols];

   for(unsigned int c=0;c<ccdcols;c++)
   {
      for(unsigned int p=0;p<newsizerows;p++)
      {
         //Assign zero value
         tmpX[c*newsizerows + p]=0;
         tmpY[c*newsizerows + p]=0;
         tmpZ[c*newsizerows + p]=0;
         for(unsigned int b=0;b<spatbin;b++)
         {
            //Add on the pixels for this bin
            tmpX[c*newsizerows + p]+=rotX[c*ccdrows+spatbin*p+b];
            tmpY[c*newsizerows + p]+=rotY[c*ccdrows+spatbin*p+b];
            tmpZ[c*newsizerows + p]+=rotZ[c*ccdrows+spatbin*p+b];
         }
         //Average the value in the bin
         tmpX[c*newsizerows + p]=tmpX[c*newsizerows + p]/spatbin;
         tmpY[c*newsizerows + p]=tmpY[c*newsizerows + p]/spatbin;
         tmpZ[c*newsizerows + p]=tmpZ[c*newsizerows + p]/spatbin;
      }
   }
   //Now we need to trim these and assign to the rotX,rotY,rotZ vars

   //Clear the current values
   delete[] rotX;
   delete[] rotY;
   delete[] rotZ;
   //Create the dynamic variables
   this->ccdrows=(newsizerows-xstart);

   //Check that this size agrees with lev1 size
   if(ccdrows!=l1samps)
   {
      std::cout<<"Number of samples: binned FOV file , level 1 file: "<<ccdrows<<" , "<<l1samps<<std::endl;
      throw "Binned view vector file does not have the same number of samples as the level 1 file. Have you given the correct -vvfile filename?\n";
   }

   this->rotX=new double[ccdrows*ccdcols];
   this->rotY=new double[ccdrows*ccdcols];
   this->rotZ=new double[ccdrows*ccdcols];  

   for(unsigned int c=0;c<ccdcols;c++)
   { 
      for(unsigned int i=0;i<ccdrows;i++)
      {
         rotX[c*ccdrows + i]=tmpX[xstart+i];
         rotY[c*ccdrows + i]=tmpY[xstart+i];
         rotZ[c*ccdrows + i]=tmpZ[xstart+i];
      }
   }
   //Remove tmp variables
   delete[] tmpX;
   delete[] tmpY;
   delete[] tmpZ;
   
}


//-------------------------------------------------------------------------
// ViewVector copy constructor
//-------------------------------------------------------------------------
ViewVectors::ViewVectors(ViewVectors &ref) 
{
   //DO NOT allow access to the BIL for a copy
   this->binf=NULL;

   //Get the size of the ccd (the part that we're interested in) from
   //the bil header file
   this->ccdrows=ref.ccdrows;
   this->ccdcols=ref.ccdcols;

   //Create the dynamic variables
   this->rotX=new double[ccdrows*ccdcols];
   this->rotY=new double[ccdrows*ccdcols];
   this->rotZ=new double[ccdrows*ccdcols];

   //Populate the rotX,rotY,rotZ variables - read in the view vector file
   for(unsigned int i=0;i<this->ccdrows*this->ccdcols;i++)
   {
      this->rotX[i]=ref.rotX[i];
      this->rotY[i]=ref.rotY[i];
      this->rotZ[i]=ref.rotZ[i];
   }
}

//-------------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------------
ViewVectors::~ViewVectors()
{
   //Destroy the bil reader if it exists
   if(this->binf!=NULL)
      delete this->binf;

   //Destroy the dynamic vars if they exist
   if(this->rotX!=NULL)
      delete[] this->rotX;

   if(this->rotY!=NULL)
      delete[] this->rotY;

   if(this->rotZ!=NULL)
      delete[] this->rotZ;
}

//-------------------------------------------------------------------------
//This function is called to read in the view vectors from the bil file
//and store the values in the dynamic variables.
//-------------------------------------------------------------------------
int ViewVectors::ReadVVFile()
{
   //The bil file should contain 3 bands (X,Y,Z angles)
   if(StringToUINT(this->binf->FromHeader("bands"))!=3)
      throw "View vector file should contain 3 bands: rotations about X, Y and Z";

   //Read in the information into the arrays if they have been defined
   if((this->rotX==NULL)||(this->rotY==NULL)||(this->rotZ==NULL))
      throw "Cannot read view vector file - rotX,Y,Z arrays not defined";

   //read in every line for each band in turn
   this->binf->Readband((char*)this->rotX,0);
   this->binf->Readband((char*)this->rotY,1);
   this->binf->Readband((char*)this->rotZ,2);

   return 1;
}

//-------------------------------------------------------------------------
//This function is called to apply angular rotations to the view vectors
//-------------------------------------------------------------------------
int ViewVectors::ApplyAngleRotations(const double rx, const double ry, const double rz)
{
   unsigned int arrsize=(this->ccdrows * this->ccdcols);
   for(unsigned int i=0;i<arrsize;i++)
   {
      this->rotX[i]=this->rotX[i] + rx;
      this->rotY[i]=this->rotY[i] + ry;
      this->rotZ[i]=this->rotZ[i] + rz;
   }
   return 1;
}


