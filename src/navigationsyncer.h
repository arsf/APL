//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef NAVIGATIONSYNCER_H
#define NAVIGATIONSYNCER_H

#include <string>
#include <cerrno>
#include <map>
#include <vector>
#include <limits>
#include "logger.h"
#include "datahandler.h"
#include "navfileclasses.h"
#include "binfile.h"
#include "binaryreader.h"

#define DEBUGNAVSYNCER

//-------------------------------------------------------------------------
// Class to hold leap seconds since January 2006
//-------------------------------------------------------------------------
class LeapSecond
{
public:
   //-------------------------------------------------------------------------
   // Leap second class constructor
   //-------------------------------------------------------------------------
   LeapSecond()
   {
      data["01-01-2006"]=14;
      data["01-01-2009"]=15;
      data["01-07-2012"]=16;
   };

   int GetLeapSeconds(std::string mydate);

private:
   std::map<std::string,int> data;
   bool CheckDateFormat(std::string testdate);
};


//-------------------------------------------------------------------------
// NavigationSyncer - class that does the scan syncing and interpolation
//-------------------------------------------------------------------------
class NavigationSyncer
{
public:
   //default constructor
   NavigationSyncer();
   //constructor taking input nav data and lev1 filenames
   NavigationSyncer(std::string navfilename,std::string lev1filename);
   //destructor
   ~NavigationSyncer();

   //Function to find a per scan line time
   void FindScanTimes();
   //Function to return a pointer to the times array
   double* PtrToTimes(){return time;};
   //Function to apply a time shift to the data
   void ApplyTimeShift(const double shift);
   //Function to apply Lead seconds
   void ApplyLeapSeconds(){for(unsigned int i=0;i<nscans;i++){time[i]+=leapseconds;}};

   double GetCropTimeOffset()const{return croptimeoffset;}

private:

   //Handle the input nav data
   SpecimFileChooser* navfile;
   float NOSYNCINHDR;
   //scan line times
   double* time;

   unsigned long nscans; // number of scan lines from level 1 hdr file
   double hdrsync; //sync time from level 1 hdr file
   double framerate; //frame rate from level 1 hdr file
   double croptimeoffset; //time offset to apply if level 1 data has been cropped

   std::string acquisitiondate; //string of the acquisition date from the hdr
   std::string gpsstarttime; //Strings of the start and stop times from the hdr
   std::string gpsstoptime;

   unsigned int leapseconds;
   int lev1firstscanmaxexpectedsize;

   //Function to compare the level1 and nav times in gps sec of week
   void CompareLev1ToNavTimes(double first_scan_time);
};

#endif
