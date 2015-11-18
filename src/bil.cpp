//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

//Methods for the BILReader class

#include "bil.h"
//#define BILREADERDEBUG
#ifndef BILREADERDEBUG
   #define DEBUGPRINT(X)  
#else
   #define DEBUGPRINT(X) std::cout<< X <<std::endl;
#endif

BILReader::BILReader(std::string fname)
{

   //Asssign filename
   this->filename=fname;
   //Want to open bil file, read in hdr file, and assign all the class variables
   this->filein=fopen(this->filename.c_str(),"rb");
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
      this->filesize=FileTell(filein);//get the position (in bytes)
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
      if( ToLowerCase(it->second) != "bil")
      {
         brinfo<<"Header file states file "<<this->filename<<" is not a bil file but a "<<it->second<<std::endl;
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

      //Test filesize on disk against filesize from header file - casting to return as a uint64 type
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

BILReader::~BILReader()
{
   DEBUGPRINT("Entering BILReader destructor...")
   Close();
}

//Read in a lines worth of data for each band from the current file pointer position
void BILReader::Readline(char* const chdata)
{
   if(this->IsGood() == false)
   {
      DEBUGPRINT("Read line has failed...stream not good.") 
      brinfo<<"Read line function has failed because 'file stream is not good.'"<<std::endl;  
      brinfo<<"Further information for debugging: FERROR value: "<<this->Ferror()<<" FEOF: "<<this->Feof()<<std::endl;
      brinfo<<"File in question: "<<this->filename<<std::endl;
      brinfo<<strerror(errno)<<std::endl;    
      throw BRexception(brinfo.str());
   }

   //DEBUG statement
   DEBUGPRINT("Reading Line from BIL...");

   //number of bytes that should be read from the file
   int64_t nbytestoread=(this->numsamples * this->datasize * this->numbands);
   if(this->CheckCapacity(nbytestoread))
   {
      this->prevpointerloc=FileTell(this->filein);
      fread(chdata,sizeof(char),nbytestoread,filein);
   }
   else
   {
      brinfo<<"Attempted to read "<<nbytestoread<<" but file reports there are less bytes remaining to be read"<<std::endl;
      throw BRexception(brinfo.str());
   }
}

//Read in numlines lines of data for each band from startline
void BILReader::Readlines(char* const chdata, unsigned int startline, unsigned int numlines)
{
   //check error bit is not set
   if(this->IsGood() == false)
   {
      DEBUGPRINT("Readlines has failed...stream not good.") 
      brinfo<<"Readlines function has failed because 'file stream is not good.'"<<std::endl; 
      brinfo<<"Further information for debugging: FERROR value: "<<this->Ferror()<<" FEOF: "<<this->Feof()<<std::endl;     
      brinfo<<"File in question: "<<this->filename<<std::endl;    
      brinfo<<strerror(errno)<<std::endl;
      throw BRexception(brinfo.str());
   }
   
   //DEBUG statement
   DEBUGPRINT("Reading Lines from BIL...from "<<startline<<" reading "<<numlines<<" lines.");
   DEBUGPRINT("Position in file (before seeking):"<<ftello(filein));
   //number of bytes that should be read from the file
   uint64_t nbytestoread=(this->numsamples * this->datasize * this->numbands) * numlines;
   //position at which to start reading from
   uint64_t pos=(this->numsamples * this->datasize * this->numbands) * startline;

   //Move the file pointer to this position
   FileSeek(filein,pos,SEEK_SET);

   if(this->CheckCapacity(nbytestoread))
   {
      fread(chdata,sizeof(char),nbytestoread,filein);
   }
   else
   {
      brinfo<<"Attempted to read "<<nbytestoread<<" but file reports there are less bytes remaining to be read"<<std::endl;
      throw BRexception(brinfo.str());
   }
}

//Read in a certain number of bytes from the current file pointer position
void BILReader::Readbytes(char* const chdata, unsigned long int nbytes)
{
   if(this->IsGood() == false)
   {
      DEBUGPRINT("Readbytes has failed...stream not good.") 
      brinfo<<"Readbytes function has failed because 'file stream is not good.'"<<std::endl; 
      brinfo<<"Further information for debugging: FERROR value: "<<this->Ferror()<<" FEOF: "<<this->Feof()<<std::endl;     
      brinfo<<"File in question: "<<this->filename<<std::endl;    
      brinfo<<strerror(errno)<<std::endl;
      throw BRexception(brinfo.str());
   }

   //DEBUG statement
   DEBUGPRINT("Reading "<<nbytes<<" bytes from BIL...");
   if(this->CheckCapacity(nbytes))
   {
      this->prevpointerloc=FileTell(this->filein);
      fread(chdata,sizeof(char),nbytes,filein);
   }
   else
   {
      brinfo<<"Attempted to read "<<nbytes<<" but file reports there are less bytes remaining to be read"<<std::endl;   
      throw BRexception(brinfo.str());
   }
}

//Reads the specified band
int BILReader::Readband(char* const chdata, unsigned int band)
{
   //DEBUG statement
   DEBUGPRINT("Reading band "<<band<<" from BIL...");
   for(unsigned int line=0;line<this->numrows;line++)
   {
      Readbandline(&chdata[this->numsamples*line*this->datasize],band,line);
   }
   return 1;
}

//reads the given line for the given band
int BILReader::Readbandline(char* const chdata, unsigned int band, unsigned int line)
{
   if(this->IsGood() == false)
   {
      DEBUGPRINT("Readbandline has failed...stream not good.") 
      brinfo<<"Readbandline function has failed because 'file stream is not good.'"<<std::endl;   
      brinfo<<"Further information for debugging: FERROR value: "<<this->Ferror()<<" FEOF: "<<this->Feof()<<std::endl;   
      brinfo<<"File in question: "<<this->filename<<std::endl;    
      brinfo<<"Errno says: "<<strerror(errno)<<std::endl;
      throw BRexception(brinfo.str());
   }

   //DEBUG statement
   DEBUGPRINT("Reading line "<<line<<" from band "<<band);
   //number of bytes that should be read from the file
   uint64_t nbytestoread=this->numsamples * this->datasize;
   //position at which to start reading from
   uint64_t pos=(this->numsamples * this->datasize * this->numbands) * line  + (this->numsamples * this->datasize * band);

   //Move the file pointer to this position
   FileSeek(filein,pos,SEEK_SET);

   if(this->CheckCapacity(nbytestoread))
      fread(chdata,sizeof(char),nbytestoread,filein);
   else
   {
      brinfo<<"Attempted to read "<<nbytestoread<<" but file reports there are less bytes remaining to be read."<<std::endl;   
      throw BRexception(brinfo.str());
   }

   return 1;
}

//Read in a cell from the file and return it as a double
double BILReader::ReadCell(const unsigned int band,const unsigned int line, const unsigned int col)
{
   if(this->IsGood() == false)
   {
      DEBUGPRINT("ReadCell has failed...stream not good.") 
      brinfo<<"ReadCell function has failed because 'file stream is not good.'"<<std::endl;   
      brinfo<<"Further information for debugging: FERROR value: "<<this->Ferror()<<" FEOF: "<<this->Feof()<<std::endl;   
      brinfo<<"File in question: "<<this->filename<<std::endl;    
      brinfo<<strerror(errno)<<std::endl;
      throw BRexception(brinfo.str());
   }   

   if((line>numrows)||(col>numsamples))
   {
      brinfo<<"ReadCell has failed...row/col out of bounds: row: "<<line<<" col: "<<col<<std::endl;   
      throw BRexception(brinfo.str());      
   }

   //Get the number of bytes to skip to get to point to read
   uint64_t nbytestoskip=((this->numsamples * this->numbands)*line + (this->numsamples * band) + col) * this->datasize;
   //Create a char array of the size of one cell
   char* cbuffer=new char[this->datasize];
   //Seek to and read that cell of data
   FileSeek(filein,nbytestoskip,SEEK_SET);

   fread(cbuffer,this->datasize,1,filein);
   //Need to dereference this, convert to double and return
   double val=DerefToDouble(cbuffer);
   delete[] cbuffer;
   return val;
}


void BILReader::ReadlineToDoubles(double* const ddata,unsigned int line)
{
   char* chtmp=new char[this->numsamples*this->numbands* this->datasize];
   Readlines(chtmp,line,1);
   for(unsigned int sample=0;sample<this->numsamples*this->numbands;sample++)
   {

      ddata[sample]=DerefToDouble(&chtmp[sample*datasize]);
   }   
   delete[] chtmp;
}

//-------------------------------------------------------------------------
//Read the BIL data that falls within a rectangular area defined by the lower left (ll)
//and upper right (ur) coordinates (Usually here BIL file is a DEM)
//-------------------------------------------------------------------------
int BILReader::ReadRect(char* const chdata,const int minrow, const int maxrow,const int mincol,const int maxcol)
{
   if(numbands!=1)
      throw "ReadRect should currently only be used for files with 1 band.";

   //Check input parameters and throw exception if they are not in correct order
   if ((maxrow < minrow) || (maxcol < mincol))
      throw "Order of elements in BILReader.ReadRect should be llx,lly urx,ury. Min row/col is greater than max row/col.";

   //Get the extent of the area in rows/cols
   unsigned int numlines=(maxrow+1)-minrow;
   unsigned int numsamps=(maxcol+1)-mincol;

   //Convert these to byte positions within the DEM file and the
   //number of bytes for each line to read. Then use bil->Readbytes
   //methodology.
   unsigned long int numofbytesperline=numsamps*GetDataSize();
   DEBUGPRINT("Will read "<<numsamps<<" samples of data with size "<<GetDataSize()<<" bytes from "<<numlines<<" lines.")

   //Loop through each of the lines reading in the required part, storing in chdata
   for(unsigned int l=0;l<numlines;l++)
   {
      //For each line to read from - read in some data
      ReadPartOfLine(&chdata[l*numofbytesperline],l+(unsigned int)minrow,(unsigned int)mincol,numsamps);
   }

   return 0;
}

//-------------------------------------------------------------------------
//Read in a section nsamps long starting from sampleno of a given line from 
//the dem and store in chdata array. chdata should be the size 
//of nsamps * size of data (Usually here BIL file is a DEM)
//-------------------------------------------------------------------------
void BILReader::ReadPartOfLine(char* const chdata,const unsigned int lineno,const unsigned int sampleno,const unsigned int nsamps)
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
   uint64_t nbytestoread=(nsamps * GetDataSize() );

   //position to move file get pointer to
   uint64_t pos=0;
   
   //Check that there are nbytestoread left in the file, if not then throw exception
   if(this->CheckCapacity(nbytestoread))
   {
      //Store the previous position
      this->prevpointerloc=FileTell(filein);
      //Skip pointer to new location
      pos=lineno*this->numsamples*GetDataSize() + sampleno*GetDataSize();
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
