//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

//Methods for the BSQReader class

#include "bsq.h"

#ifndef BILREADERDEBUG
   #define DEBUGPRINT(X)  
#else
   #define DEBUGPRINT(X) std::cout<< X <<std::endl;
#endif


BSQReader::BSQReader(std::string fname)
{
   //Asssign filename
   this->filename=fname;
   //Want to open bsq file, read in hdr file, and assign all the class variables
   filein=fopen(this->filename.c_str(),"rb");
   if(!this->IsOpen())
   {
      //File has failed to open
      brinfo<<"An error has occured opening the file: "<<this->filename<<std::endl;
      this->filesize=-1;
      this->numrows=-1;
      this->numsamples=-1;
      this->numbands=-1;
      this->datasize=-1;
      throw BRexception(brinfo.str()); //Throw an exception
   }
   else
   {
      //File exists and has been opened successfully
      //Get size of file
      FileSeek(filein,0,SEEK_END); //seek to end of file
      this->filesize=FileTell(filein); //get the position (in bytes)
      FileSeek(filein,0,SEEK_SET); //seek back to start of file
      //Now try and read in hdr file
      int retval=this->ReadHeader();
      //An error occured, store a general error message into the brinfo stream
      if(retval!=1)
      {
         brinfo<<"An error has occurred whilst reading from the hdr file."<<std::endl;      
         throw BRexception(brinfo.str()); //Throw an exception
      }
      std::map<std::string,std::string>::iterator it; //use an iterator to access the Header map
      //set the numsamps, number of lines etc and also perform other checks against file
      it=this->Header.find("interleave"); //Search for interleave type
      if(it==this->Header.end())//check if present in header
      {
         this->CannotFind("interleave");
      }
      if( ToLowerCase(it->second) != "bsq")
      {
         brinfo<<"Header file states file "<<this->filename<<" is not a bsq file but a "<<it->second<<std::endl;
         throw BRexception(brinfo.str()); //Throw an exception
      }
      //get the number of lines
      it=this->Header.find("lines");
      if(it==this->Header.end()) //doesnt exist in Header map/file
         this->CannotFind("lines");
      else
         this->numrows=StringToUINT(it->second);

      it=this->Header.find("samples");
      if(it==this->Header.end()) //doesnt exist in Header map/file
         this->CannotFind("samples");
      else
         this->numsamples=StringToUINT(it->second);

      it=this->Header.find("bands");
      if(it==this->Header.end()) //doesnt exist in Header map/file
         this->CannotFind("bands");
      else
         this->numbands=StringToUINT(it->second);

      it=this->Header.find("data type");
      if(it==this->Header.end()) //doesnt exist in Header map/file
         this->CannotFind("data type");
      else
         this->datatype=StringToUINT(it->second);

      DatasizeFromDatatype(); //assign data size

      //Test filesize on disk against filesize from header file
      uint64_t fsize=static_cast<uint64_t>(this->numrows)*this->numsamples*this->numbands*this->datasize;

      if(fsize != this->filesize)
      {
         brinfo<<"Filesize from header: "<<fsize<<" does not agree with filesize from disk (for "<<this->filename<<"): "<<this->filesize<<std::endl;
         throw BRexception(brinfo.str()); //Throw an exception
      }

      //Test for datatypes 6 or 9 - these are not capable of being read
      if((this->datatype==6)||(this->datatype==9))
      {
         brinfo<<"Data type is complex (real/imaginary). This file type is unsupported"<<std::endl;
         throw BRexception(brinfo.str()); //Throw an exception
      }
   }
}

//destructor
BSQReader::~BSQReader()
{
   Close();
}

//Read 1 line for all bands [equivalent number of bytes from current position]
void BSQReader::Readline(char* const chdata) 
{
   throw BRexception("This method is deprecated and should not be used (It doesn't even exist - maybe it should).");
}

//Read numlines lines of data from startline
//Data is stored in chdata as band0line0;band0line1;band0linex;band1line0;band1line1;...;bandx;linex
void BSQReader::Readlines(char* const chdata, unsigned int startline, unsigned int numlines) 
{
   //check error bit is not set
   if(this->IsGood() == false)
   {
      DEBUGPRINT("Read band line has failed...stream not good.") 
      brinfo<<"Read band line has failed...stream not good."<<std::endl;   
      throw BRexception(brinfo.str());
   }
   
   //DEBUG statement
   DEBUGPRINT("Reading Lines from BSQ...from "<<startline<<" reading "<<numlines<<" lines.");
   //number of bytes that should be read from the file
   uint64_t nbytestoread=(this->numsamples * this->datasize * this->numbands) * numlines;
   uint64_t bytesperline=(this->numsamples * this->datasize); 
   //position at which to start reading from
   //unsigned long int pos=bytesperline * startline;
   //Move the file pointer to this position
   //this->filein.seekg(pos,std::ios::beg);

   //Cannot do the check on capacity (using the binary reader function) for this function as 
   //it is not reading a continuous block of data
   for(unsigned int band=0;band<numbands;band++)
   {
      for(unsigned int line=startline;line<startline+numlines;line++)
      {
         Readbandline(&chdata[bytesperline*(line-startline)+band*numlines*bytesperline], band, line);
      }
   }
}

//Read the specified number of bytes from current position
void BSQReader::Readbytes(char* const chdata, uint64_t nbytes) 
{
   if(this->IsGood() == false)
   {
      DEBUGPRINT("Read bytes has failed...stream not good.") 
      brinfo<<"Read bytes has failed...stream not good."<<std::endl;   
      throw BRexception(brinfo.str());
   }

   //DEBUG statement
   DEBUGPRINT("Reading "<<nbytes<<" bytes from BSQ...");
   if(this->CheckCapacity(nbytes))
   {
      this->prevpointerloc=FileTell(filein);
      fread(chdata,sizeof(char),nbytes,filein);
   }
   else
   {
      brinfo<<"Attempted to read "<<nbytes<<" but file reports there are less bytes remaining to be read"<<std::endl;   
      throw BRexception(brinfo.str());
   }
}

//Reads the specified band
int BSQReader::Readband(char* const chdata, unsigned int band) 
{
   if(this->IsGood() == false)
   {
      DEBUGPRINT("Read band line has failed...stream not good.") 
      brinfo<<"Read band line has failed...stream not good."<<std::endl;   
      throw BRexception(brinfo.str());
   }

   //DEBUG statement
   DEBUGPRINT("Reading Line from BSQ...");

   //number of bytes that should be read from the file
   uint64_t nbytestoread=(this->numsamples * this->datasize * this->numrows);
   if(this->CheckCapacity(nbytestoread))
   {
      this->prevpointerloc=FileTell(filein);
      fread((chdata),sizeof(char),nbytestoread,filein);
   }
   else
   {
      brinfo<<"Attempted to read "<<nbytestoread<<" but file reports there are less bytes remaining to be read"<<std::endl;
      throw BRexception(brinfo.str());
   }
   return 1;
}

//reads the given line for the given band
int BSQReader::Readbandline(char* const chdata, unsigned int band, unsigned int line) 
{
   if(this->IsGood() == false)
   {
      DEBUGPRINT("Read band line has failed...stream not good.") 
      brinfo<<"Read band line has failed...stream not good."<<std::endl;   
      throw BRexception(brinfo.str());
   }

   //DEBUG statement
   DEBUGPRINT("Reading line "<<line<<" from band "<<band);
   //number of bytes that should be read from the file
   uint64_t nbytestoread=this->numsamples * this->datasize;
   //position at which to start reading from
   uint64_t pos=(this->numsamples * this->datasize * this->numrows)*band + (this->numsamples * this->datasize * line);
   //Move the file pointer to this position
   //this->filein.seekg(pos,std::ios::beg);
   FileSeek(filein,pos,SEEK_SET);

   if(this->CheckCapacity(nbytestoread))
      fread(chdata,sizeof(char),nbytestoread,filein);
   else
   {
      brinfo<<"Attempted to read "<<nbytestoread<<" but file reports there are less bytes remaining to be read"<<std::endl;   
      throw BRexception(brinfo.str());
   }
   return 1;
}

//Read in a cell from the file and return it as a double
double BSQReader::ReadCell(const unsigned int band,const unsigned int line, const unsigned int col)
{
   if(this->IsGood() == false)
   {
      DEBUGPRINT("ReadCell has failed...stream not good.") 
      brinfo<<"ReadCell has failed...stream not good."<<std::endl;   
      throw BRexception(brinfo.str());
   }   
   //Get the number of bytes to skip to get to point to read
   uint64_t nbytestoskip=((this->numsamples * this->numrows)*band + (this->numsamples * line) + col) * this->datasize;
   //Create a char array of the size of one cell
   char* cbuffer=new char[this->datasize];
   //Seek to and read that cell of data
   //this->filein.seekg(nbytestoskip,std::ios::beg);;
   FileSeek(filein,nbytestoskip,SEEK_SET);

   fread(cbuffer,this->datasize,1,filein);
   //Need to dereference this, convert to double and return
   double val=DerefToDouble(cbuffer);
   delete[] cbuffer;
   return val;  
}

void BSQReader::ReadlineToDoubles(double* const ddata,unsigned int line)
{
   char* chtmp=new char[this->numsamples*this->numbands* this->datasize];
   Readlines(chtmp,line,1);

   for(unsigned int sample=0;sample<this->numsamples*this->numbands;sample++)
   {
      //The following switch statement allows the data in chtmp to be derefenced
      //to many types, allowing different data formats in the BSQ to be used 
      switch(this->datatype)
      {
      case 1: //8-bit
         {
            char* cderef=(char*)(chtmp);      
            ddata[sample]=static_cast<double>(cderef[sample]);
            break;   
         }
      case 2: //16 bit signed short int
         {
            short int* sideref=(short int*)(chtmp);      
            ddata[sample]=static_cast<double>(sideref[sample]);
            break;
         }
      case 3: //32 bit signed short int
         {
            int* ideref=(int*)(chtmp);      
            ddata[sample]=static_cast<double>(ideref[sample]);
            break;
         }
      case 4: //float
         {
            float* fderef=(float*)(chtmp);
            ddata[sample]=static_cast<double>(fderef[sample]);
            break;
         }
      case 5: //double
         {
            double* dderef=(double*)(chtmp);
            ddata[sample]=static_cast<double>(dderef[sample]);
            break;
         }
      case 12: //16 bit unsigned short int
         {
            unsigned short int* usideref=(unsigned short int*)(chtmp);
            ddata[sample]=static_cast<double>(usideref[sample]);
            break;
         }
      case 13: //32 bit unsigned int
         {
            unsigned int* uideref=(unsigned int*)(chtmp);
            ddata[sample]=static_cast<double>(uideref[sample]);
            break;
         }
      default:
         throw "Unrecognised data type for BSQ file. Currently supports 8-bit, signed and unsigned 16-bit and 32-bit integer, 32-bit and 64-bit float";
         break;
      }
   }   
   delete[] chtmp;
}

//-------------------------------------------------------------------------
// Functon to read a rectangle from a BSQ (DEM) file 
// - This is just a direct copy of the BILReader function since number
// of bands has to be 1 for the DEM file. Done this way to allow future
// update to multiple band files if ever required
//-------------------------------------------------------------------------
int BSQReader::ReadRect(char* const chdata,const int minrow, const int maxrow,const int mincol,const int maxcol)
{
   if(numbands!=1)
      throw "ReadRect should currently only be used for files with 1 band.";

   //Check input parameters and throw exception if they are not in correct order
   if ((maxrow < minrow) || (maxcol < mincol))
      throw "Order of elements in BSQReader.ReadRect should be llx,lly urx,ury. Min row/col is greater than max row/col.";

   //Get the extent of the area in rows/cols
   unsigned int numlines=(unsigned int)(maxrow-minrow+1);
   unsigned int numsamps=(unsigned int)(maxcol-mincol+1);

   //Convert these to byte positions within the BSQ file and the
   //number of bytes for each line to read. Then use bil->Readbytes
   //methodology.
   uint64_t numofbytesperline=numsamps*this->GetDataSize();
   DEBUGPRINT("Will read "<<numsamps<<" samples of data with size "<<this->GetDataSize()<<" bytes from "<<numlines<<" lines.")

   //Loop through each of the lines reading in the required part, storing in chdata
   for(unsigned int l=0;l<numlines;l++)
   {
      //For each line to read from - read in some data
      ReadPartOfLine(&chdata[l*numofbytesperline],l+(unsigned int)minrow,(unsigned int)mincol,numsamps);
   }
   return 0;

}

//-------------------------------------------------------------------------
// Function to read in part of a line from a BSQ (DEM) file
//-------------------------------------------------------------------------
void BSQReader::ReadPartOfLine(char* const chdata,const unsigned int lineno,const unsigned int sampleno,const unsigned int nsamps)
{
   if(this->IsGood() == false)
   {
      DEBUGPRINT("ReadPartOfLine has failed...stream not good.") 
      this->brinfo<<"ReadPartOfLine has failed...stream not good."<<std::endl;   
      throw BRexception(brinfo.str());
   }

   //DEBUG statement
   DEBUGPRINT("Reading Part Of Line from BinaryFile..."<<nsamps<<" samples starting at sample "<<sampleno<<" of line "<<lineno);

   //number of bytes that should be read from the file
   uint64_t nbytestoread=(nsamps * this->GetDataSize() );

   //position to move file get pointer to
   uint64_t pos=0;
   
   //Check that there are nbytestoread left in the file, if not then throw exception
   if(this->CheckCapacity(nbytestoread))
   {
      //Store the previous position
      this->prevpointerloc=FileTell(filein);
      //Skip pointer to new location
      pos=lineno*this->numsamples*this->GetDataSize() + sampleno*this->GetDataSize();
      //this->filein.seekg(pos,std::ios::beg);
      FileSeek(filein,pos,SEEK_SET);
      //Read in the length of data
      fread(chdata,sizeof(char),nbytestoread,filein);
   }
   else
   {
      brinfo<<"Attempted to read "<<nbytestoread<<" but file reports there are less bytes remaining to be read"<<std::endl;
      throw BRexception(brinfo.str());
   }
}

