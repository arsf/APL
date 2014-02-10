//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef BINARYREADER_H 
#define BINARYREADER_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <string>
#include <exception>
#include <stdint.h>

#include "commonfunctions.h"

//Class that holds the common code for reading in BSQ and BIL files including all the header information.
//The reading functions are overloaded in the bilreader and bsqreader classes.
//Note this is not actually abstract - can create a instance of it but only to tell you the filestyle (bsq/bil) - no reading methods and header destroyed.

#ifdef _W32
   #define FileSeek(X,Y,Z) fseeko64(X,Y,Z);
   #define FileTell(X) ftello64(X);
#else
   #define FileSeek(X,Y,Z) fseeko(X,Y,Z);
   #define FileTell(X) ftello(X);
#endif


class BinaryReader
{
public:
   BinaryReader(); //constructor   
   BinaryReader(std::string fname); //constructor with filename - only to return the filetype (BIL or BSQ)
   
   virtual ~BinaryReader();//destructor

   std::string FromHeader(std::string key) //Search the map for the string key and return  the string item. (Return a string from the map)
   {
      std::map<std::string,std::string>::iterator it;
      it=this->Header.find(ToLowerCase(key));
      if(it!=this->Header.end())
         return it->second;
      else
         return "";
   }       

   std::string FromHeader(std::string key,int itemnum); //Search the map for the string key and return the itemnum variable on the string. (for multiple objects)
   unsigned int GetDataSize() const {return datasize;}
   unsigned int GetDataType() const {return datatype;}

   virtual void Readline(char* const chdata){throw BRexception("Undefined function call.");} //Read 1 line for all bands [equivalent number of bytes from current position]
   virtual void Readline(char* const chdata, unsigned int line) //Read specified line for all bands
         {Readlines(chdata,line,1);} //This is essentially a special case of readlines where numlines=1
   virtual void Readlines(char* const chdata, unsigned int startline, unsigned int numlines){throw BRexception("Undefined function call.");}; //Read numlines lines of data from startline
   virtual void Readbytes(char* const chdata, unsigned long int bytes){throw BRexception("Undefined function call.");}; //Read the specified number of bytes from current position
   virtual int Readband(char* const chdata, unsigned int band){throw BRexception("Undefined function call.");}; //Reads the specified band
   virtual int Readbandline(char* const chdata, unsigned int band, unsigned int line){throw BRexception("Undefined function call.");}; //reads the given line for the given band

   virtual double ReadCell(const unsigned int band,const unsigned int line, const unsigned int col){throw BRexception("Undefined function call.");};
   virtual void ReadlineToDoubles(double* const ddata,unsigned int line){throw BRexception("Undefined function call.");};

   virtual bool IsOpen()const{if(filein==NULL){return false;} else {return true;}} //test if the binary BIL file is open
   virtual bool IsGood()const{if(ferror(filein)!=0){return false;} else if(feof(filein)!=0){return false;}else{return true;} }//Test if eof/fail/bad bits set
   virtual int Ferror()const{return ferror(filein);}
   virtual int Feof()const{return feof(filein);}
   void Close(); //Close the file, free any allocated memory, call the destructor

   //Class to handle exceptions
   class BRexception
   {
   public:
      std::string info; //pointer to BILReader stringstream bilinfo
      BRexception(std::string ss){info=ss;};  

      const char* what() const throw()
      {
         return "A BinaryReader Exception has occurred.";
      }
   
   };

   enum interleavetype {BSQ,BIL,UNKNOWN};
   interleavetype GetFileStyle(){return FILESTYLE;};

   std::string HeaderDump(bool ret);

   //Including this as it is required by DEMBinFile - maybe better to inheri this class as DEMBinaryReader?
   virtual int ReadRect(char* const chdata, const int minrow, const int maxrow,const int mincol,const int maxcol) {throw "Undefined function call.";};
   //Tidy up strings passed to header file
   virtual std::string TidyForHeader(std::string totidy);

   uint64_t GetFileSize() const {return filesize;}
   std::map<std::string, std::string, cmpstr> CopyHeader()const{return Header;}
   //Function to return the filename
   std::string GetHeaderFilename(){return hdrfilename;}


protected:

   interleavetype FILESTYLE;

   FILE* filein;
   std::string filename; //Name of the opened file.
   std::string hdrfilename;
   uint64_t filesize; //Size of the bil (binary) file
   uint64_t numrows,numsamples, numbands; //description of the binary BIL file - made into uint64_t so as to avoid casting everywhere they're used
   unsigned long int datasize; //description of the binary BIL file
   unsigned int datatype; //data type value from envi hdr (related to datasize variable)  
   std::stringstream brinfo; // a stringstream object to store information/warnings/errors in. This is to keep from having to put couts in the functions

   std::map<std::string, std::string, cmpstr> Header; //Associative array (map) for containing the BIL header. Keys and mapped values are strings
   int ReadHeader(); //Reads in the hdr file and stores the parameters as required. Perform check on BIL filesize.
   void DatasizeFromDatatype(); //uses the ENVI datatype to get the data size
   void CannotFind(std::string key); //function to throw an exception because an element from Header does not exist

   bool CheckCapacity(uint64_t numbytes); //Checks there are enough bytes after file pointer to read in numbytes of data
   void GotoPreviousPosition(); //move the file pointer to the previous position (i.e. the position in prevpoi nterloc)
   off_t prevpointerloc; //the previous location of the file pointer. Useful in case you want to read a chunk of data "out of order"
                                     //but then be able to return to the previous position to carry on reading "in order"
   //Function to get the header filename from BIL filename
   std::string GetHeaderFileName(const short type);

   double DerefToDouble(char* cbuffer);

};



#endif


