//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#include "binfile.h"

//-------------------------------------------------------------------------
// BinFile Constructor - takes filename as argument
//-------------------------------------------------------------------------


BinFile::BinFile(std::string filename)
{
   br=NULL;
   br=new BinaryReader(filename);
   BinaryReader::interleavetype ft=br->GetFileStyle();

   delete br;
   br=NULL;

   if(ft==BinaryReader::BSQ)
      br=new BSQReader(filename);
   else if(ft==BinaryReader::BIL)
      br=new BILReader(filename);      
   else
      throw "Error. Iterleave type in file: "+filename+" is not bsq or bil.";
}

//-------------------------------------------------------------------------
// BinFile dectructor
//-------------------------------------------------------------------------

BinFile::~BinFile()
{
   if(br!=NULL)
      delete br;
}

//-------------------------------------------------------------------------
// Function to return a copy of the header excluding parts which define
// the binary data - for use when copying parts of files or applying
// masks etc.
//-------------------------------------------------------------------------
std::map<std::string, std::string, cmpstr> BinFile::CopyHeaderExcluding()
{
   std::map<std::string, std::string, cmpstr> header=this->CopyHeader();
   for(std::map<std::string, std::string, cmpstr>::iterator iter=header.begin();iter!=header.end();iter++)
   {
      if(((*iter).first.compare("bands")==0) ||
         ((*iter).first.compare("samples")==0) ||
         ((*iter).first.compare("lines")==0) ||
         ((*iter).first.compare("data type")==0) ||
         ((*iter).first.compare("ENVI")==0) ||
         ((*iter).first.compare("interleave")==0) ||
         ((*iter).first.compare("description")==0)||
         ((*iter).first.compare("")==0) )
      {
         header.erase(iter);
      }        
   }
   return header;
}

