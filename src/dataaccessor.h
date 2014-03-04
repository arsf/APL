//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef DATAACCESSOR_H
#define  DATAACCESSOR_H

#include "binfile.h"
#include "logger.h"
//#include "linesegment.h"

//-------------------------------------------------------------------------
// Template class that defines a block of data of type T
// - A block is essentially a section of continuous data that can
//   be 'safely' accessed by a pointer
//-------------------------------------------------------------------------
template <class T>
class Block
{
public:
   Block(T* bd,unsigned int nl,unsigned int ns,unsigned int nb,unsigned int fr,unsigned int er)
   {
      Logger::Verbose("Constructing Block: number of lines="+ToString(nl)+" number of samples="+ToString(ns)+" number of bands="+ToString(nb)
                        +" first and last rows:"+ToString(fr)+" "+ToString(er));
      blockdata=bd;
      nlines=nl;
      nsamples=ns;
      nbands=nb;
      firstrow=fr;
      endrow=er;
   }
   ~Block(){Logger::Verbose("Destructing block...");}

   T* const Data() const{return blockdata;}
   const unsigned int FirstRow()const{return firstrow;}
   const unsigned int EndRow()const{return endrow;}
   const unsigned int Samples()const{return nsamples;}
   const unsigned int Bands()const{return nbands;}
   const unsigned int Lines()const{return nlines;}

private:
   T* blockdata;
   unsigned int nlines,nsamples,nbands;
   unsigned int firstrow,endrow;

};


//-------------------------------------------------------------------------
// Class to wrap the reading of data from a file or a block if it is in RAM
//-------------------------------------------------------------------------
template <class T>
class DataAccessor
{
public:
   DataAccessor(Block<T>* b,std::string level1filename,unsigned int* blist);
   ~DataAccessor();

   T GetData(unsigned int band, unsigned int row, unsigned int col);

private:
   Block<T>* block;
   BinFile* file;
   unsigned int* bandlist;
};

//-------------------------------------------------------------------------
// Constructor - takes pointer to block, level1 filename and list of bands to map
//-------------------------------------------------------------------------
template <class T>
DataAccessor<T>::DataAccessor(Block<T>* b,std::string filename,unsigned int* blist)
{
   block=NULL;
   file=NULL;
   bandlist=NULL;

   block=b;
   if(filename != "")
      file=new BinFile(filename);
   bandlist=blist;
}

//-------------------------------------------------------------------------
// Destructor - deletes BinFile object
//-------------------------------------------------------------------------
template <class T>
DataAccessor<T>::~DataAccessor()
{
   if (file != NULL)
      delete file;
}

//-------------------------------------------------------------------------
// Function that does the work and returns the requested data pixel
//-------------------------------------------------------------------------
template <class T>
T DataAccessor<T>::GetData(unsigned int band, unsigned int row, unsigned int col)
{
   if((block==NULL) || (row < block->FirstRow()) ||(col<0) || (row >= block->EndRow()) || (col >=block->Samples()))
   {
      //Need to access from file
      return this->file->ReadCell(band,row,col);       
   }
   else
   {
      //Can access from block
      //Need to convert real band number to block band number
      //made b static - this will be quicker if calls to this function are consequtive band
      //requests - which is the most likely procedure.
      static unsigned int b=0;
      while((b<block->Bands())&&(bandlist[b]!=band))
      {
         b++;
      }
      if(b==block->Bands())
      {
         //Repeat the while loop starting from b=0 since b is static and the band we want may be before the initial b value
         b=0;
         while((b<block->Bands())&&(bandlist[b]!=band))
         {
            b++;
         }
         if(b==block->Bands())
         {
            //Logger::Error("Trying to read from a band not in memory in DataAccessor - defaulting to read from file.");
            return this->file->ReadCell(band,row,col);       
         }
      }
      
      return (block->Data()[block->Bands()*block->Samples()*(row-block->FirstRow()) + b*block->Samples() + col]);
   }
}

#endif
