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
//Class definition for the navigation base class.
//This will contain information on the per-scan line navigation
//written for use in geocorrection software and will be inherited
//from in the navigation processing software.
//-------------------------------------------------------------------------

#ifndef NAVBASECLASS_H
#define NAVBASECLASS_H

#include <string>
#include "binfile.h"
#include "commonfunctions.h"

class NavBaseClass
{
   public:
      NavBaseClass(std::string fname); //constructor
      virtual ~NavBaseClass(); //destructor

      int ReadScan(const unsigned int scannumber);//function to read info for given scan

      //Following are functions to access private/protected member data
      double Time()const{return time;}
      double Lat()const{return lat;}
      double Lon()const{return lon;}
      double Hei()const{return hei;}
      double Roll()const{return roll;}
      double Pitch()const{return pitch;}
      double Heading()const{return heading;}    
      unsigned int ScanID()const{return scanid;}  

      double MinLat()const{return minlat;}
      double MinLon()const{return minlon;}
      double MinHei()const{return minhei;}
      double MinRoll()const{return minroll;}
      double MaxLat()const{return maxlat;}
      double MaxLon()const{return maxlon;}
      double MaxHei()const{return maxhei;}
      double MaxRoll()const{return maxroll;}

      //Scan through the nav file and find min/max of lat/lon/height
      const void FindLimits();
      const void FindLimits(unsigned int lowerscan,unsigned int upperscan);

      //Functions to access member data from BIL - if too many then just make bil public
      //This function returns the total number of scans in the file
      unsigned int TotalScans() const {return StringToUINT(binf->FromHeader("lines"));}

   protected:
      //Data variables
      double lat,lon,hei,roll,pitch,heading,time; //one scans worth of information
      unsigned int scanid; //scan id - runs from 0 to number_of_scans in file

      //Min/Max values of the nav file
      double maxlat,maxlon,maxhei,maxroll,minlat,minlon,minhei,minroll;

      BinFile* binf; //To read in the per-scan line navigation from Binary (BIL) file

};

#endif
