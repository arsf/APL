//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef NAVINTERPOLATOR_H
#define NAVINTERPOLATOR_H

#include <string>
#include <cerrno>
#include <fstream>
#include "datahandler.h"
#include "navfileclasses.h"
#include "binfile.h"
#include "binaryreader.h"
#include "bilwriter.h"
#include "leverbore.h"

#include "interpolationfunctions.h"
enum FILETYPE{SBET,SPECIMNAV,SOL,BADFILE};
//-------------------------------------------------------------------------
// NavigationInterpolator - class that does the scan syncing and interpolation
//-------------------------------------------------------------------------
class NavigationInterpolator
{
public:
   //default constructor
   NavigationInterpolator();
   //constructor taking input nav data and lev1 filenames
   NavigationInterpolator(std::string navfilename,std::string lev1filename);
   //destructor
   ~NavigationInterpolator();

   //Function to write the data out into a BIL file for applcorr
   void Writer(std::string outfilename,std::string extrainfo="");
   void WriteFlags(std::string outfilename,std::string extrainfo="");
   //Function to interpolate the navigation data to the level 1 scan times
   //void Interpolate(void (*f)(double,DataHandler*,NavDataLine*,std::string,std::string));
   void Interpolate(void (*f)(const double*,const int,DataHandler*,NavDataCollection*,std::string,std::string));
   //Function to interpolate the navigation data to time t
   //void Interpolate(const double t, void (*f)(double,DataHandler*,NavDataLine*,std::string,std::string),NavDataLine* store);
   //Function to assign scan line times
   void SetTimes(const double* const times);
   //Function to apply boresight offsets to scanline data
   void ApplyBoresight(Boresight* boresight);
   //Function to apply leverarm offsets to scanline data
   void ApplyLeverarm(Leverarm* leverarm);
   //Function to smooth the navigation data (use for the raw data)
   //void Smooth(void (*f)(const unsigned long ,DataHandler* ,NavDataLine* ,const int),const int element, NavDataLine* store);
   //NavDataLine* Smooth(void (*f)(const unsigned long ,DataHandler* ,NavDataLine* ,const int));
   void SmoothNavData(void (*f)(const unsigned long ,DataHandler* ,NavDataLine* ,const int),const unsigned int smoothkernelsize){dhandle->Smooth(f,smoothkernelsize);};
   //Apply a shift between the positions and attitude data
   void PosAttShift(void (*f)(const double*,const int,DataHandler*,NavDataCollection*,std::string,std::string),const double toffset);
   //Check the plausibilty of the interpolated data
   void CheckPlausibility();
private:
   //Function to set up the navdataline and scanid arrays
   void SetupArrays();
   //Handle the input nav data
   DataHandler* dhandle;
   //Store the interpolated per scan line data
//   NavDataLine* navdata;
   NavDataCollection* navcollection;

   double* scanid;
   unsigned long nscans; // number of scan lines from level 1 hdr file
   std::string gpsstarttime; //Strings of the start and stop times from the hdr
   std::string gpsstoptime;

   FILETYPE DetectFileType(std::string filename);
};

#endif
