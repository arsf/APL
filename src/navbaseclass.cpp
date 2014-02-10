//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
//Member function definitions for the navbaseclass
//-------------------------------------------------------------------------

#include "navbaseclass.h"


//-------------------------------------------------------------------------
//Constructor using nav filename
//-------------------------------------------------------------------------
NavBaseClass::NavBaseClass(std::string fname)
{
   //Create the BILreader object
   this->binf=new BinFile(fname);
      
   //Test BIL file has expected dimensions
   //That is: 1 samples per scan line
   if(StringToUINT(this->binf->FromHeader("samples"))!=1)
      throw "This BIL file does not have expected dimensions. Processed Navigation should have 1 samples per scan.";
   //And number of bands should be 7 (la,lo,hei,p,r,hea)
   if(StringToUINT(this->binf->FromHeader("bands"))!=7)
      throw "This BIL file does not have expected dimensions. Processed Navigation should have 7 bands.";

   //Read in the first scan line
   ReadScan(0);

   //Initialise min/maxs
   this->minlat=this->minlon=this->minhei=this->minroll=9999;
   this->maxlat=this->maxlon=this->maxhei=this->maxroll=-9999;
}

//-------------------------------------------------------------------------
//destructor
//-------------------------------------------------------------------------
NavBaseClass::~NavBaseClass()
{
   if(this->binf!=NULL)
      delete this->binf;
}

//-------------------------------------------------------------------------
//Function to read in a scan's data from the per-scan navigation BIL file
//-------------------------------------------------------------------------
int NavBaseClass::ReadScan(const unsigned int scannumber)
{
   //Check that scannumber is less than total number of scans
   if(scannumber >= this->TotalScans())
   {
      //Throw an exception
      throw "Cannot read scan from file, scan number is larger than total number of scans in file.";
   }

   //navfile should be 1 sample, 7 bands, N scans
   if( (StringToUINT(this->binf->FromHeader("bands"))!=7 )||(StringToUINT(this->binf->FromHeader("samples"))!=1))
      throw "Processed Navigation BIL file should have 7 bands and 1 sample. Failed to read scan.";

   //Read in the 7 bands of data for the given scan number
   double tmp[7]={0};
   this->binf->Readline((char*)tmp,scannumber);
   //Store the data in the 'proper' variables
   this->time=tmp[0];
   this->lat=tmp[1];
   this->lon=tmp[2];
   this->hei=tmp[3];
   this->roll=tmp[4];
   this->pitch=tmp[5];
   this->heading=tmp[6];
   //set the id of the scan that has been read in
   this->scanid=scannumber;
   //return the error status
   return 1; //success
}

//-------------------------------------------------------------------------
//Find the min and max of elements of the navigation file
//-------------------------------------------------------------------------
const void NavBaseClass::FindLimits()
{
   FindLimits(0,this->TotalScans());
}

//-------------------------------------------------------------------------
//Find the min and max of elements of a section of the navigation file
//-------------------------------------------------------------------------
const void NavBaseClass::FindLimits(unsigned int lowerscan,unsigned int upperscan)
{
   //Check upperscan is > lowerscan
   if(upperscan <= lowerscan)
   {
      throw "Upperscan cannot be less than lowerscan in NavBaseClass::FindLimits.";
   }
   //Initialise min/maxs - this is incase run multiple times using different scan line bounds (e.g. 0->1000, 1001->5000)
   this->minlat=this->minlon=this->minhei=this->minroll=9999;
   this->maxlat=this->maxlon=this->maxhei=this->maxroll=-9999;

   //Scan through the nav file and update min/max accordingly
   for(unsigned int i=lowerscan;i<upperscan;i++)
   {
      //Read in the scan line
      this->ReadScan(i);
      //Check values for latitude
      if(this->lat > this->maxlat)
         this->maxlat=this->lat;
      if(this->lat < this->minlat)
         this->minlat=this->lat;
      //Check values for longitude
      if(this->lon > this->maxlon)
         this->maxlon=this->lon;
      if(this->lon < this->minlon)
         this->minlon=this->lon;
      //Check values for height
      if(this->hei > this->maxhei)
         this->maxhei=this->hei;
      if(this->hei < this->minhei)
         this->minhei=this->hei;
      //Check values for roll
      if(this->roll > this->maxroll)
         this->maxroll=this->roll;
      if(this->roll < this->minroll)
         this->minroll=this->roll;
   }
}
