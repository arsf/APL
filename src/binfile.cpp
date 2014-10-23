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
   //Need to identify if a single bil/bsq has been passed or if it is a
   //ascii file with a list of bil/bsq filenames within
   bool SINGLEFILE=true;
   br=NULL;
   try
   {
      //If it opens correctly we assume it is a normal bil/bsq
      br=new BinaryReader(filename);
   }
   catch(BinaryReader::BRexception bre)
   {
      //something other than no hdr file happened ... throw it
      if(bre.info.find("Unable to open a hdr file, does one exist?")==std::string::npos)
         throw bre;

      //A no hdr found error - could be because it is a multifile
      SINGLEFILE=false;
      br=NULL;
   }

   if(SINGLEFILE)
   {
      //Get the interleave type and open an instance of that type of file reader
      BinaryReader::interleavetype ft=br->GetFileStyle();

      delete br;
      br=NULL;

      if(ft==BinaryReader::BSQ)
         br=new BSQReader(filename);
      else if(ft==BinaryReader::BIL)
         br=new BILReader(filename);      
      else
         throw "Error. Interleave type in file: "+filename+" is not bsq or bil.";
   }
   else
   {
      //Try and open a multi file type reader
      br=new MultiFile(filename);
   }
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

