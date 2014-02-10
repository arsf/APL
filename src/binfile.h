//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef BINFILE_H
#define BINFILE_H

#include "binaryreader.h"
#include "bil.h"
#include "bsq.h"

//-------------------------------------------------------------------------
// Class to use for BIL/BSQ reader - make an instance of this one to use
// for either a BIL or BSQ file
//-------------------------------------------------------------------------

class BinFile
{
public:
   BinFile(){};//default constructor for inheritence
   BinFile(std::string fname);
   ~BinFile();

   void Readline(char* const chdata){br->Readline(chdata);}
   void Readline(char* const chdata, unsigned int line){br->Readline(chdata,line);} 
   void Readlines(char* const chdata, unsigned int startline, unsigned int numlines){br->Readlines(chdata,startline,numlines);} 
   void Readbytes(char* const chdata, unsigned long int bytes){br->Readbytes(chdata,bytes);}
   int Readband(char* const chdata, unsigned int band){return br->Readband(chdata,band);} 
   int Readbandline(char* const chdata, unsigned int band, unsigned int line){return br->Readbandline(chdata,band,line);}

   double ReadCell(const unsigned int band,const unsigned int line, const unsigned int col){return br->ReadCell(band,line,col);}
   void ReadlineToDoubles(double* const ddata,unsigned int line){return br->ReadlineToDoubles(ddata,line);}

   std::string FromHeader(std::string key){return br->FromHeader(key);}
   std::string FromHeader(std::string key,int itemnum){return br->FromHeader(key,itemnum);} 

   unsigned int GetDataSize() const {return br->GetDataSize();}
   unsigned int GetDataType() const {return br->GetDataType();}
   void Close(){br->Close();}
   std::string HeaderDump(bool ret){return br->HeaderDump(ret);}
   std::string TidyForHeader(std::string totidy){return br->TidyForHeader(totidy);}
   uint64_t GetFileSize() const {return br->GetFileSize();}
   std::map<std::string, std::string, cmpstr> CopyHeader()const{return br->CopyHeader();}

   std::map<std::string, std::string, cmpstr> CopyHeaderExcluding();

   std::string GetHeaderFilename(){return br->GetHeaderFilename();}

protected:
   BinaryReader* br;
};


#endif
