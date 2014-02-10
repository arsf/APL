//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef DATAHANDLER_H
#define DATAHANDLER_H

#include <cstdio>
#include <string>
#include <sstream>
#include <iostream>
#include <cmath>
#include <typeinfo>
#include "logger.h"
#include "commonfunctions.h"


//Constants for checking plausibilty of the data - currently these values are used
//for both SBET and .nav file and also for interpolated scan data. May be better to have different sets.
const double plausibleTimeDifference=0.05;
const double plausibleHeightDifference=0.5;
const double plausibleLatDifference=0.0001;
const double plausibleLonDifference=0.0001;
const double plausibleRollDifference=0.12;
const double plausiblePitchDifference=1.0;
const double plausibleHeadingDifference=1.0;


//-------------------------------------------------------------------------
// NavDataLine class - holds data from an epoch
//-------------------------------------------------------------------------
class NavDataLine
{
public:
   NavDataLine();
   double lat, lon, hei, roll, pitch, heading, time;
   char quality; //quality flag 0 is good data, bad data has below values
   enum {BADLAT=1,BADLON=2,BADHEI=4,BADROLL=8,BADPITCH=16,BADHEADING=32,BADTIME=64};
};


//-------------------------------------------------------------------------
// NavDataCollection class - an array of NavDataLines
//-------------------------------------------------------------------------
class NavDataCollection
{
public:
   NavDataCollection(const unsigned long len);
   ~NavDataCollection();
   const unsigned long SizeOfArray(){return size_of_array;}
   NavDataLine* GetLine(const unsigned long l)
   {
      if((l>=0)&&(l<size_of_array))
         return &navarray[l];
      else
         return NULL;
   };
   enum NavDataItem {LAT,LON,TIME,ROLL,PITCH,HEADING,HEI,QUALITY};
   void SetValues(const unsigned long item,const NavDataLine* line);
   void SetValues(const unsigned long item,const NavDataItem key,const double value);
   const double GetValue(const unsigned long item,const NavDataItem key);
   const char GetFlag(const unsigned long item){return navarray[item].quality;}
   double* GetReferenceValue(const unsigned long item,const NavDataItem key);
   char* GetFlagReference(const unsigned long item){return &navarray[item].quality;}
   //Function to check the nav data looks plausible
   void CheckPlausibility();

private:   
   unsigned long size_of_array;
   NavDataLine* navarray;

};

//-------------------------------------------------------------------------
// DataHandler class - base class for the SBETData and SpecimNavData classes
//-------------------------------------------------------------------------
class DataHandler
{
public:
   DataHandler();
   DataHandler(unsigned long n);
   virtual ~DataHandler();

   //Function that reads in the navigation file 
   virtual void Reader()=0;

   //Function that returns a ptr to NavDataLine object relating to data from line l.
   NavDataLine* GetLine(const unsigned long l){return navcollection->GetLine(l);}

   //Return the numentries value
   unsigned long GetNumEntries(){return navcollection->SizeOfArray();}
   //Get the sync time i
   virtual double GetSyncTime(const unsigned long i){return 0;}
   //Get the number of sync messages found
   virtual unsigned long GetNumSyncs(){return 0;}
   //Function to smooth the data (once it has been read in)
   void Smooth(void (*f)(const unsigned long ,DataHandler* ,NavDataLine* ,const int),const unsigned int smoothkernelsize);
   void Smooth(void (*f)(const unsigned long ,DataHandler* ,NavDataLine* ,const int),const int element, NavDataLine* store,const unsigned int smoothkernelsize);
   void CheckPlausibility(){navcollection->CheckPlausibility();}

   virtual std::string GetInformation()
   {
      std::stringstream stream(std::stringstream::in | std::stringstream::out);
      stream<<"Start and end times of file: Type "<<typeid(*this).name()<<" "<<navcollection->GetValue(0,NavDataCollection::TIME)
            <<" "<<navcollection->GetValue(GetNumEntries()-1,NavDataCollection::TIME)<<std::endl;
      return stream.str();
   }

protected:
   std::string filename;
   NavDataCollection* navcollection;

};


#endif
