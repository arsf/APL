//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef SPECIMSENSORS_H
#define SPECIMSENSORS_H

#include <iostream>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cerrno>
#include <cmath>
#include <vector>

#include "sensor.h"
#include "binfile.h"
#include "bilwriter.h"
#include "commonfunctions.h"
#include "logger.h"


//-------------------------------------------------------------------------
//Specialist BIL reader for Specim files. Includes all normal BILreader 
//functions but with extra functions to read in the dark lines and
//calculate the total number of missing frames.
//-------------------------------------------------------------------------
class SpecimBILReader: public BILReader
{
public:
   //Constructor
   SpecimBILReader(std::string fname,bool force=false):BILReader(fname)
   {
      DARKFORCE=force;
      sensorid_string=TestFile();
      darkscalar=1.0;
   }
 
   //-------------------------------------------------------------------------
   //Specim Specific functions but probably dont belong in general BILReader class
   //-------------------------------------------------------------------------

   std::string TestFile(); //Test the file for possible errors etc

   void AverageAllDarkFrames(double* const data,std::string externalFileName="");
   void DarkFramesStdDeviation(double* const stdev,const double* const mean, std::string externalFileName="");
   void AverageRefinedDarkFrames(double* const data,const double* const stdev,const double* const mean, std::string externalFileName="");
   void AverageAllDarkFramesFromFile(double* const darkframes,std::string strDarkFileName);

   short int GetTotalMissingFrames() const {return totalmissing;}
   unsigned short GetMissingFramesBetweenLimits(unsigned int start,unsigned int end){return TotalMissingFrames(start,end);}
   unsigned int GetNumDarkFrames() const {return ndarklines;}
   unsigned int GetNumImageFrames() {return StringToUINT(FromHeader("lines")) - ndarklines;}
   double DarkScalar() const {return darkscalar;}

private:
   bool DARKFORCE;//force the use of the autodarkstartline given in hdr file
   short int totalmissing;//total number of missing frames
   unsigned int ndarklines;
   unsigned int darklinestart;
   const static unsigned int MAXFRAMECOUNT=65535;
   double darkscalar; //scalar if separate dark file used of different integration time
   std::string sensorid_string;

   //Reads the first/last frame IDs and compare to number of lines to get total number of missing frames
   void TotalMissingFrames(); 
   //Reads the start/end frame IDs and compare to number of lines within these limits to get number of missing frames
   unsigned short TotalMissingFrames(const unsigned int start, const unsigned int end);
};



//-------------------------------------------------------------------------
// Specim BinFile class that inherits from BinFile to allow a SpecimBILReader
//-------------------------------------------------------------------------
class SpecimBinFile: public BinFile
{
public:
   SpecimBinFile(std::string strFilename,bool force=false)
   {
      br=new BinaryReader(strFilename);
      BinaryReader::interleavetype ft=br->GetFileStyle();
      delete br;

      if(ft==BinaryReader::BSQ)
         throw "SpecimBSQReader not implemented yet.";
      else if(ft==BinaryReader::BIL)
         br=new SpecimBILReader(strFilename,force);      
      else
         throw "Error. Iterleave type in file: "+strFilename+" is not bsq or bil.";
   }

private:
};


//-------------------------------------------------------------------------
// Silly little class made for handling the qcfailures
//-------------------------------------------------------------------------
class Pair
{
public:
   Pair(unsigned int b,unsigned int s){band=b;sample=s;}
   ~Pair(){};
   unsigned int band;
   unsigned int sample;
};

//-------------------------------------------------------------------------
// Specim object
//-------------------------------------------------------------------------
class Specim: public Sensor
{
public:
   //constructor with sensor id and calibration file prefix,and raw filename
   Specim(std::string strFilename,bool force=false);  

   //destructor
   virtual ~Specim(); 

   //Enum the various mask pixel values
   enum MaskType {Good=0, UnderFlow=1, OverFlow=2, Badpixel=4, SmearAffected=8,DroppedScan=16,CorruptData=32,QCFailure=64};

   SpecimBILReader* br;

   unsigned int SpectralBinning()const{return spectralbinning;}
   unsigned int SpatialBinning()const{return spatialbinning;}

   unsigned int NumBands()const{return numbands;}
   unsigned int NumLines()const{return numlines;}
   unsigned int NumSamples()const{return numsamps;}
   double IntegrationTime() const {return tint;}
   unsigned int RadianceScalar() const {return radscalar;}

   int SensorID() const {return SENSOR_ID;}
   unsigned short int RawMax() const {return rawmax;}
   unsigned short int CalibratedMax() const {return calibratedmax;}

   unsigned int LowerScanlineLimit() const{return scanlinelowerlimit;}
   unsigned int UpperScanlineLimit() const{return scanlineupperlimit;}

   std::string CalibratedUnits()const{return calibratedunits;}
   std::string RawFilename() const {return strRawFilename;}

   //List of band,sample that failed a quality control check and should be added to the mask
   std::vector<Pair> qcfailures;

protected:

   //-------------------------------------------------------------------------
   //Constants relating to the Specim sensors
   //-------------------------------------------------------------------------
   unsigned short int rawmax; //maximum value
   unsigned short int calibratedmax; //maximum value of a uint16
   unsigned int radscalar; //used to multiply the radiance so can store a more accurate value as uint16
   int SENSOR_ID; //id string for the Eagle sensor (serial no.)
   std::string calibratedunits;

   //Raw data file name  - used to derive other filenames
   std::string strRawFilename;

   //-------------------------------------------------------------------------
   //General Specim sensor variables
   //-------------------------------------------------------------------------
   //char sensorid; //'e'agle or 'h'awk
   unsigned int prevfnum; //previous frame number - unsigned int (not short) so can have a flag value
   double tint; //integration time of sensor
   unsigned int spatialbinning,spectralbinning; // spatial and spectral binning values
   unsigned int numbands, numsamps, numlines;   //number of bands, samples and lines relating to the image data array
   unsigned int scanlineupperlimit, scanlinelowerlimit; //himg {} values from specim header
};


class Eagle: public Specim
{
public:
   Eagle(std::string strFilename,bool force=false);
   ~Eagle();

   double FrameTransferTime() const {return trant;}
   unsigned long LowerFodis()const{return lowerfodis;}
   unsigned long UpperFodis()const{return upperfodis;}
   std::string FodisUnits()const{return fodisunits;}
protected:
   double trant;//frame transfer time for Eagle in ms - from CCD document (recommended operating rate - not necessarily true rate!)
   unsigned int lowerfodis,upperfodis; // hold the lower and upper fodis bounds
   std::string fodisunits;
};


class Hawk: public Specim
{
public:
   Hawk(std::string strFilename,bool force=false);
   ~Hawk();

   void DecodeBadPixel(std::map<int,int> revbandmap); //function to fill in the bad pixel array from specim files

   unsigned int NumBadPixels()const{return nbadpixels;}
   const unsigned int* BadPixels(){return badpixels;}
   const unsigned char* BadPixelMethod(){return badpixelmethod;}
   unsigned int BandNotInUse()const{return bandnotinuse;}
   const std::string* MethodDescriptor(){return bpmethod_descriptor;}
   std::stringstream badpixelstream;//stringstream to hold bad pixels

   bool arsfbadpixelfiletype;//true if ARSF file false if specim file
   enum BadPixelMethodName {None=0, A=1, B=2, C=4, D=8,E=16};

protected:
   //-------------------------------------------------------------------------
   //Bad pixel variables
   //-------------------------------------------------------------------------
   unsigned int bandnotinuse;//flag value used for hawk bad pixels - to show when a band is not in use(e.g spectral binning of 2, band numbers 127->256 dont exist)
   unsigned int* badpixels;  //array to hold the bad pixels in (sample,band)
   unsigned char* badpixelmethod; //array to hold method used to detect bad pixel for ARSF calibrated bad pixel file
   std::string* bpmethod_descriptor;//array of method descriptors
   unsigned int nbadpixels; //number of bad pixels

   void DecodeARSFBadPixels(std::map<int,int> revbandmap);
   void DecodeSpecimBadPixels(std::map<int,int> revbandmap); 
};


#endif
