//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef BILWRITER_H
#define BILWRITER_H

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <ctime>
#include "commonfunctions.h"

#ifdef EXE_NAME
const std::string bildescriptionstring="BIL file created by "+std::string(EXE_NAME);
#else
const std::string bildescriptionstring="BIL file created by ARSF BILWriter";
#endif


class BILWriter
{
public:
   BILWriter(std::string filename,unsigned int datasize); //constructor to open a new file filename and write data of type datasize 
                                                          //(e.g. 1 byte (char),2 byte (uint) etc)
                                                          //if file exists it will throw an error
   BILWriter(std::string filename,unsigned int datasize,char cmethod); //constructor to open file using method cmethod ('a' - append, 'w' -overwrite)

   //constructor to assign BIL dimensions at "start up" as well as above methods
   BILWriter(std::string filename, unsigned int datatype, unsigned int nrows,unsigned int nsamps,unsigned int nbands,char cmethod);

   virtual ~BILWriter(); //destructor

   //int WriteBytes(char* data,unsigned long int numbytes); //function to write numbytes worth of data array to file
   int WriteBandLine(char* const data);//write a line of data fro 1 band
   int WriteLine(char* const data);

   //Function to write a number of lines of data (for all bands)
   void WriteLines(char* const data,const unsigned int nl)
   {
      int ret=0;
      for(unsigned int i=0;i<nl;i++)
      {
         ret=WriteLineSection(&data[i*this->numsamples*this->numbands*this->datasize],this->numsamples,0,numsamples-1);
         if(ret==-1){std::cout<<bilinfo.str()<<std::endl;}
      }
   }

   //Function to write a line from given start point to end point.
   //i.e. write a section of a line (for all bands) from an array of a larger size - also need the number_of_samples of the data array
   //since these may be different (WILL BE) to output number of samples
   //NOTE numsamples_array IS ONLY THE NUMBER OF SAMPLES (NOT SAMPLES*BANDS*LINES etc)
   int WriteLineSection(char* const data,const unsigned int numsamples_array,const unsigned int start, const unsigned int end);
   //int WriteLineSection(const char* data,unsigned int start, unsigned int end); 

   int WriteLineWithValue(const char xval); //write out data with constant value over the line

   int Close(); //Close the file, output the hdr information, check nrows/samples/bands agree with written filesize and call the destructor

   void SetBILDimensions(unsigned int nrows,unsigned int nsamps, unsigned int nbands) //Set the private BIL size variables
      {numrows=nrows;numsamples=nsamps;numbands=nbands;}       //Should be no need to use this function if coded correctly.

   void AddToHdr(std::string item){hdrtext<<item<<std::endl;} //Adds the string item to the stringstream of header text
   bool IsGood()const{return isgood;}//Get the status of the BILWriter (1=good 0=bad)
   std::string GetBilInfo() const{return bilinfo.str();}
   virtual void PrepareHeader(); //Prepare the header file stringstream to contain the relevant info
   unsigned int GetDataSize() const {return datasize;}

   //Exception class
   class BILexception
   {
   public:
      std::string info; //pointer to BILwriter stringstream bilinfo
      BILexception();
      BILexception(std::string ss){info=ss;};  

      const char* what() const throw()
      {
         return "A BIL Exception has occurred.";
      }
   
   };

   enum DataType {uchar8, char8, uint16, int16, uint32, int32, float32, float64 };

   unsigned int GetDataType() const {return datatype;}

   //numsamples_array = in terms of if array = [s*b*l] then numsamples_array=s NOT s*b*l
   template<class T>
   void WriteDataToLineSection(T const data,const unsigned int numsamples_array,const unsigned int start, const unsigned int end)
   {
      unsigned char* cp=NULL;
      unsigned short* usp=NULL;
      unsigned int* uip=NULL;
      short* sp=NULL;
      int* ip=NULL;
      float* fp=NULL;
      double* dp=NULL;
      switch(datatype)
      {
      case 1:
         cp=new unsigned char[numsamples_array*numbands];
         for(unsigned int i=0;i<numsamples_array*numbands;i++)
         {
            cp[i]=static_cast<unsigned char>(data[i]+0.5);// add on 0.5 to round to nearest
         }
         WriteLineSection((char*)cp,numsamples_array,start,end);
         delete[] cp;
         break;
      case 2:
         sp=new short[numsamples_array*numbands];
         for(unsigned int i=0;i<numsamples_array*numbands;i++)
         {
            sp[i]=static_cast<short>(data[i]+0.5);// add on 0.5 to round to nearest
         }
         WriteLineSection((char*)sp,numsamples_array,start,end);
         delete[] sp;
         break;
      case 3:
         ip=new int[numsamples_array*numbands];
         for(unsigned int i=0;i<numsamples_array*numbands;i++)
         {
            ip[i]=static_cast<int>(data[i]+0.5);// add on 0.5 to round to nearest
         }
         WriteLineSection((char*)ip,numsamples_array,start,end);
         delete[] ip;
         break;
      case 4:
         fp=new float[numsamples_array*numbands];
         for(unsigned int i=0;i<numsamples_array*numbands;i++)
         {
            fp[i]=static_cast<float>(data[i]);
         }
         WriteLineSection((char*)fp,numsamples_array,start,end);
         delete[] fp;
         break;
      case 5:
         dp=new double[numsamples_array*numbands];
         for(unsigned int i=0;i<numsamples_array*numbands;i++)
         {
            dp[i]=static_cast<double>(data[i]);
         }
         WriteLineSection((char*)dp,numsamples_array,start,end);
         delete[] dp;
         break;
      case 12:
         usp=new unsigned short[numsamples_array*numbands];
         for(unsigned int i=0;i<numsamples_array*numbands;i++)
         {
            usp[i]=static_cast<unsigned short>(data[i]+0.5); // add on 0.5 to round to nearest
         }
         WriteLineSection((char*)usp,numsamples_array,start,end);
         delete[] usp;
         break;
      case 13:
         uip=new unsigned int[numsamples_array*numbands];
         for(unsigned int i=0;i<numsamples_array*numbands;i++)
         {
            uip[i]=static_cast<unsigned int>(data[i]+0.5);// add on 0.5 to round to nearest
         }
         WriteLineSection((char*)uip,numsamples_array,start,end);
         delete[] uip;
         break;
      default:
         break;
      }
   }

private:
   unsigned int numrows,numsamples,numbands,datasize,datatype;
   std::string filename; //name of output bil file (without an extension, so will output to filename.bil, filename.hdr)
   std::ofstream fileout; //stream object

   std::stringstream bilinfo;//string stream object to hold status information/error messages etc of BIL

   std::stringstream hdrtext; //string stream to hold additional text to output to header file

   bool isgood;//status flag
   void SetSizeFromType(const unsigned int type); //set the datatype/datasize values
   int WriteHeader(); //Creates and outputs the associated header file
};

#endif
