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
#include "filewriter.h"

#ifdef EXE_NAME
const std::string bildescriptionstring="BIL file created by "+std::string(EXE_NAME);
#else
const std::string bildescriptionstring="BIL file created by ARSF BILWriter";
#endif


class BILWriter : public FileWriter
{
public:
   //constructor to assign BIL dimensions at "start up" as well as above methods
   BILWriter(std::string filename, unsigned int datatype, unsigned int nrows,unsigned int nsamps,unsigned int nbands,char cmethod);
   virtual ~BILWriter(); //destructor

   int WriteBandLine(char* const data);//write a line of data for 1 band
   int WriteLine(char* const data);
   int WriteBandLineSection(char* const data,const unsigned int numsamples_array,const unsigned int start, const unsigned int end);
   int WriteBandLineWithValue(const char xval){throw "This function is no longer used - stub remains for base class compatibility.";}; //write out data with constant value over the line

   template<class T>
   void WriteBandLineWithValue(T xval) //write out data with constant value over the line
   {
      T* data=new T[this->numsamples]; //create array of data of size of 1 line for all bands
      for(unsigned int i=0;i<this->numsamples;i++)
         data[i]=static_cast<T>(xval);
      WriteDataToBandLineSection(data,this->numsamples,0,this->numsamples-1);
      delete data;
      data=NULL;
   }

   int Close(); //Close the file, output the hdr information, check nrows/samples/bands agree with written filesize and call the destructor

   void AddToHdr(std::string item){hdrtext<<item<<std::endl;} //Adds the string item to the stringstream of header text
   void AddMetadata(std::string name,std::string value);

   bool IsGood()const{return isgood;}//Get the status of the BILWriter (1=good 0=bad)
   unsigned int GetDataSize() const {return datasize;}

   //Exception class
   class BILexception : public FileWriterException
   {
   public:
      BILexception() : FileWriterException(){}
      BILexception(std::string ss) : FileWriterException(ss){}
      BILexception(std::string ss,const char* extra): FileWriterException(ss,extra){}
      const char* what() const throw()
      {
         return "A BIL Exception has occurred.";
      }
   };

   unsigned int GetDataType() const {return datatype;}

   //numsamples_array = in terms of if array = [s*b*l] then numsamples_array=s NOT s*b*l
   //Data is an array of only 1 line for 1 band length
   template<class T>
   void WriteDataToBandLineSection(T const data,const unsigned int numsamples_array,const unsigned int start, const unsigned int end)
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
         cp=new unsigned char[numsamples_array];
         for(unsigned int i=0;i<numsamples_array;i++)
         {
            cp[i]=static_cast<unsigned char>(data[i]+0.5);// add on 0.5 to round to nearest
         }
         WriteBandLineSection((char*)cp,numsamples_array,start,end);
         delete[] cp;
         break;
      case 2:
         sp=new short[numsamples_array];
         for(unsigned int i=0;i<numsamples_array;i++)
         {
            sp[i]=static_cast<short>(data[i]+0.5);// add on 0.5 to round to nearest
         }
         WriteBandLineSection((char*)sp,numsamples_array,start,end);
         delete[] sp;
         break;
      case 3:
         ip=new int[numsamples_array];
         for(unsigned int i=0;i<numsamples_array;i++)
         {
            ip[i]=static_cast<int>(data[i]+0.5);// add on 0.5 to round to nearest
         }
         WriteBandLineSection((char*)ip,numsamples_array,start,end);
         delete[] ip;
         break;
      case 4:
         fp=new float[numsamples_array];
         for(unsigned int i=0;i<numsamples_array;i++)
         {
            fp[i]=static_cast<float>(data[i]);
         }
         WriteBandLineSection((char*)fp,numsamples_array,start,end);
         delete[] fp;
         break;
      case 5:
         dp=new double[numsamples_array];
         for(unsigned int i=0;i<numsamples_array;i++)
         {
            dp[i]=static_cast<double>(data[i]);
         }
         WriteBandLineSection((char*)dp,numsamples_array,start,end);
         delete[] dp;
         break;
      case 12:
         usp=new unsigned short[numsamples_array];
         for(unsigned int i=0;i<numsamples_array;i++)
         {
            usp[i]=static_cast<unsigned short>(data[i]+0.5); // add on 0.5 to round to nearest
         }
         WriteBandLineSection((char*)usp,numsamples_array,start,end);
         delete[] usp;
         break;
      case 13:
         uip=new unsigned int[numsamples_array];
         for(unsigned int i=0;i<numsamples_array;i++)
         {
            uip[i]=static_cast<unsigned int>(data[i]+0.5);// add on 0.5 to round to nearest
         }
         WriteBandLineSection((char*)uip,numsamples_array,start,end);
         delete[] uip;
         break;
      default:
         break;
      }
   }

   //Data is an array of only 1 line for all bands length
   template<class T>
   void WriteDataToLineSection(T const data,const unsigned int numsamples_array,const unsigned int start, const unsigned int end)
   {
      for(unsigned int b=0;b<numbands;b++)
         WriteDataToBandLineSection(&data[numsamples_array*b],numsamples_array,start,end);
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
   virtual void PrepareHeader(); //Prepare the header file stringstream to contain the relevant info
};

#endif
