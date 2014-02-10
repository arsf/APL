//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef CALIBRATION_H
#define CALIBRATION_H

#include <typeinfo>

#include "logger.h"
#include "binfile.h"
#include "bilwriter.h"
#include "specimsensors.h"

 

//-------------------------------------------------------------------------
// Class to hold the various data used in the calibration
//-------------------------------------------------------------------------
class Data
{
public:
   Data(unsigned long size);
   ~Data();

   enum TransformArray{BAND,SAMPLE};

   double* Fodis()const{return fodis;}
   unsigned char* Mask()const{return mask;}
   unsigned char* BadPixMethod()const{return badpixmethod;}
   double* Image()const{return image;}
   double* AverageDark()const{return avdark;}
   double* Gains()const{return gains;};

   unsigned long ArraySize()const{return arraysize;}
   void InitialiseFodis();
   void InitialiseMask();
   void InitialiseBadPixMethod();
   void InitialiseDarkFrames();
   void InitialiseGains();

   void TransformArrays(const unsigned int bands, const unsigned int samples,TransformArray order);

   template<class T>
   void FlipBandData(T* array, const unsigned int bands, const unsigned int samples)
   {
      T* datacopy=new T[2*samples]; //create an array the size of two lines
      unsigned int icount=0;
      unsigned int max=bands-1;  
      unsigned int rb=max;//set rb (reverse band) variable to be max of array position

      for(unsigned int b=0;b<rb;b++)// for each band (upto rb) of data array
      {
         for(icount=0;icount<samples;icount++)
         {
            //copy the frame of the band of data
            datacopy[icount]=array[b*samples + icount]; //enter a value into copy array of band x
            datacopy[icount+samples]=array[rb*samples + icount]; //enter a value into copy array of band y 
         }

         //Now need to swap this data over
         for(icount=0;icount<samples;icount++)
         {
            array[b*samples + icount]=datacopy[icount+samples];
            array[rb*samples + icount]=datacopy[icount];
         }
         
         //Now decrease rb
         rb--;  
      }
      delete[] datacopy;
   }

   template<class T>
   void FlipSampleData(T* array, const unsigned int bands, const unsigned int samples)
   {
      T* datacopy=new T[samples]; //create an array the size of one line
      unsigned int icount=0,jcount=0;
      unsigned int max=samples-1;
      
      jcount=max;//set jcount variable to be max of array position
 
      for(unsigned int b=0;b<bands;b++)// for each band of data array
      {
         jcount=max;//set jcount variable to be max of array position
         for(icount=0;icount<samples;icount++)
         {
            datacopy[icount]=array[b*samples + jcount]; //enter a value into copy array
            jcount--; //decrease original data array position
         }
         //Now swap the data over from copy back to orig
         for(icount=0;icount<samples;icount++)
         {
            array[b*samples + icount]=datacopy[icount];
         }
      }
      delete[] datacopy;
   }

   template<class T>
   void ClearArray(T* array)
   {
      for(unsigned int i=0;i<arraysize;i++)
      {
         array[i]=0;
      }
   }
 

private:
   double* fodis;
   unsigned char* mask;
   unsigned char* badpixmethod;
   double* image;
   double* avdark;
   double* gains;

   unsigned long arraysize;

};

//-------------------------------------------------------------------------
// Class to perform calibration algorithms
//-------------------------------------------------------------------------
class Calibration
{
public:
   
   Calibration(Specim* sensor,std::string calFile="");
   ~Calibration();

   void RemoveDarkFrames();
   bool SmearCorrect();
   void ApplyGains();
   void FlagPixels();
   bool AverageFodis();

   void TestCalfile();
   void InitialiseDarkFrames(std::string darkfile="");
   void ReadLineOfRaw(unsigned int line);
   int CheckFrameCounter(unsigned int start,unsigned int end);
   void ClearPerlineData();

   const Data* pData()const{return data;}
   std::string CalibrationFile()const{return calibrationFilenamePrefix;}

   void InitialiseMask();
   void InitialiseBadPixMethod();
   void InitialiseFodis();
   bool ReadBadPixelFile();
   void ReadQCFailureFile(std::string qcfailurefile);

private:
   std::string calibrationFilenamePrefix;
   Data* data;
   Specim* sensor;
   
   void CheckCalWavelengths(float* const wl_cal, const unsigned int numwl_cal);
   void ReadBinAndTrimGains(double* const trimmedcal);
   void AssignMaskValue(const unsigned int ele,const Specim::MaskType type);

   //Maps to hold the mapping from raw band -> cal file bands and vice versa
   std::map<int,int>bandmap; //raw to cal
   std::map<int,int>revbandmap; //reverse of bandmap (cal to raw)
};



#endif
