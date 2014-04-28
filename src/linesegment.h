//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef LINESEGMENT_H
#define LINESEGMENT_H

#include "level3grid.h"
#include "dataaccessor.h"
#include <limits>


//-------------------------------------------------------------------------
// Class that defines a continuous section of a flight line
// Is now a template class so that it can support different level1 types
//-------------------------------------------------------------------------
template <class T>
class LineSegment
{
public:
   LineSegment(unsigned int fr,unsigned int er,unsigned int ob,double psx,double psy,std::string bandlist,std::string igmfilename,std::string level1filename,const Area* const region);
   ~LineSegment();

   Level3GridInfo* segmentinfo;
   Level3Outline* outline;

   Block<double>* igm;
   Block<T>* level1;

   void OffsetToGrid(Level3GridInfo* ginfo);

private:
   double* igmblock;
   T* level1block;
};


template <class T>
LineSegment<T>::LineSegment(unsigned int fr,unsigned int er,unsigned int overlap_lines,double psx,double psy,std::string bandlist,std::string igmfilename,std::string level1filename,const Area* const region=NULL)
{
   //May as well set all to NULL although shouldn't matter.
   segmentinfo=NULL;
   outline=NULL;
   igm=NULL;
   level1=NULL;
   igmblock=NULL;
   level1block=NULL;

   //firstrow=fr;
   //endrow=er; // note this is one after the last row
   uint64_t nlines=er-fr;
   unsigned int nbandsl1=static_cast<unsigned int>(GetNumberOfItemsFromString(bandlist," "));

   Logger::Verbose("Constructing LineSegment starting and ending at rows: "+ToString(fr)+" "+ToString(er)+" using filename: "+igmfilename);

   //Create an IGM worker and read in the required section of the IGM file
   Basic_IGM_Worker igmr(igmfilename);
   unsigned int nsamples=igmr.Samples();

   //The overlap region is to allow a buffer of IGM data (and level 1 data) to be kept in RAM
   //to speed up searches near the start/end of the line segments. The Block data cover this
   //area but the level3grid and outline use only the data without the overlap buffer. 

   //variables for use with overlapped region
   uint64_t first_line=0;
   uint64_t end_line=0;
   //unsigned int igmblock_start_of_segment=0;

   //If we can apply the overlap lines to the start of the segment then do so
   //having to cast as signed else wraps to very large number - could maybe improve test instead
   if((signed)fr - (signed)overlap_lines >= 0)
   {
      first_line=fr-overlap_lines;
      //igmblock_start_of_segment=2*nsamples*overlap_lines;      
   }
   else //else just use the original bounds
   {
      first_line=fr;
      //igmblock_start_of_segment=0;    
   }

   //If we can apply the overlap lines to the end of the segment then do so
   if(er + overlap_lines < igmr.Lines())
   {
      end_line=er+overlap_lines;
   }
   else
   {
      end_line=er;
   }

   Logger::Debug("Using first_line and end_line of: "+ToString(first_line)+" "+ToString(end_line));

   //Update the number of lines to include the overlap region
   uint64_t nlines_with_overlap=end_line - first_line;

   Logger::Verbose("Creating IGM block of size (bytes): "+ToString(nsamples*nlines_with_overlap*2*sizeof(double)));
   igmblock=new double[nsamples*nlines_with_overlap*2]; //2 bands
   //Initialise min values to large numbers and max values to small numbers
   double maxx=-std::numeric_limits<double>::max(),minx=std::numeric_limits<double>::max();
   double maxy=-std::numeric_limits<double>::max(),miny=std::numeric_limits<double>::max();
   double tmaxx=0,tminx=0,tmaxy=0,tminy=0;
   
   //Firstly read in all the data including the overlaps
   for(uint64_t i=0;i<nlines_with_overlap;i++)
   {
      //Read in X for line i
      igmr.fin->Readbandline((char*)&(igmblock[2*i*nsamples]),0,first_line+i);
      //Read in Y for line i
      igmr.fin->Readbandline((char*)&(igmblock[(2*i+1)*nsamples]),1,first_line+i);

   }

   //Then calculate min/max only using the non-overlap data
   //non-overlap data either starts at 0 or 2*overlap_lines*nsamples into the igmblock
   for(uint64_t i=(fr-first_line);i<nlines+(fr-first_line);i++)
   {
      //Calculate the min/max X values of all the data being read in
      GetArrayLimits(&(igmblock[2*i*nsamples]),nsamples,tminx,tmaxx,igmr.IgnoreValue());
      if(tmaxx > maxx)
         maxx=tmaxx;
      if(tminx < minx)
         minx=tminx;

      //Calculate the min/max Y values of all the data being read in
      GetArrayLimits(&(igmblock[(2*i+1)*nsamples]),nsamples,tminy,tmaxy,igmr.IgnoreValue());
      if(tmaxy > maxy)
         maxy=tmaxy;
      if(tminy < miny)
         miny=tminy;
   }

   //Check if any part of segment possibly lies within the supplied area
   //If not then return NULL instead of creating a linesegment object
   // ...as we cannot return a value from a constructor we throw an exception instead 
   // (please do not change the content of exception without changing matching catch)
   if((region!=NULL)&&((miny > region->MaxY()) || (maxy < region->MinY()) 
         || (minx > region->MaxX()) || (maxx < region->MinX())))
      throw std::string("LineSegment not created as it is outside given region");

   //Set the IGM block up
   igm=new Block<double>(igmblock,nlines_with_overlap,nsamples,2,first_line,end_line);

   //Create a grid info object describing this section
   Logger::Verbose("Creating grid info using min x,y max x,y: "+ToString(minx)+" "+ToString(miny)+" "+ToString(maxx)+" "+ToString(maxy));
   segmentinfo=new Level3GridInfo(minx,miny,maxx,maxy,psx,psy,bandlist,false);

   //Need to somewhere offset the segmentinfo bounds Tx,Ty such that they are
   //a divisor of the pixel size and the level3GridInfo Tx,Ty
   //We shall do this as a call in the Map::MapLineSegment function

   //Create an outline describing this section
   //outline=new Level3Outline(segmentinfo,&igmblock[igmblock_start_of_segment],nlines,nsamples,igmr.IgnoreValue());
   outline=new Level3Outline(segmentinfo,igm,nlines,abs(first_line-fr),igmr.IgnoreValue());

   //Read in the level1 data
   Logger::Verbose("Creating line segment level1 block of size (bytes): "+ToString(sizeof(T)*(nsamples*nlines_with_overlap*nbandsl1)));
   try
   {
      level1block=new T[nsamples*nlines_with_overlap*nbandsl1];
   }
   catch(std::bad_alloc) //only catch std::bad_alloc here - leave it for further up try/catch chain if other errors occur
   {
      throw "Trying to allocate more RAM than system will allow - please use -buffersize to set a smaller amount (in MB) to use. Default is 1024MB. Note that 32-bit applications can not address more than 3GB of RAM even if your system has more avaiable.";
   }

   //----------------------------------------------------------------------
   // Check the level 1 and IGM file have the same number of lines
   //----------------------------------------------------------------------
   BinFile lev1(level1filename);
   if(igmr.Lines() != StringToUINT(lev1.FromHeader("lines")))
      throw "Number of lines in level 1 file does not agree with number of lines in IGM file. Are you sure these are for the same flight line?";      

   for(uint64_t i=0;i<nlines_with_overlap;i++)
   {
      for(unsigned int b=0;b<nbandsl1;b++)
      {
         lev1.Readbandline((char*)&(level1block[i*nsamples*nbandsl1 + b*nsamples]),segmentinfo->Bands()[b],first_line+i);
      }
   }

   //Set the level1 block up
   level1=new Block<T>(level1block,nlines_with_overlap,nsamples,nbandsl1,first_line,end_line);
   lev1.Close();
}


//-------------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------------
template <class T>
LineSegment<T>::~LineSegment()
{
   Logger::Verbose("Destructing line segment.");
   if(igmblock!=NULL)
      delete[] igmblock;
   if(level1block!=NULL)
      delete[] level1block;
   if(segmentinfo!=NULL)
      delete segmentinfo;
   if(outline!=NULL)
      delete outline;
   if(igm!=NULL)
      delete igm;
   if(level1!=NULL)
      delete level1;
}

//-------------------------------------------------------------------------
// Function to adjust the segmentgrid such that it exactly
// overlaps the ginfo grid lines
//-------------------------------------------------------------------------
template <class T>
void LineSegment<T>::OffsetToGrid(Level3GridInfo* ginfo)
{
   const double EPSILON=0.000000001;
   Logger::Verbose("Offsetting LineSegment to Grid.");

   if((ginfo->PixelSizeX()!=this->segmentinfo-> PixelSizeX()) || (ginfo->PixelSizeY()!=this->segmentinfo-> PixelSizeY()))
   {
      throw "Level3Grid and segmentinfo do not have the same pixel size - they can not be overlaid.";
   }

   //Need to offset the segmentinfo area bounds TX and TY such that they are
   //divisors of the ginfo bounds TX and TY and multiples of pixelsize
   //The bounds should ALWAYS be extended, rather than reduced, upto a max
   //of the level3 grid TX and TY

   double modTLX=fmod(segmentinfo->TopLeftX()-ginfo->TopLeftX(),segmentinfo->PixelSizeX());
   double modTLY=fmod(ginfo->TopLeftY()-segmentinfo->TopLeftY(),segmentinfo->PixelSizeY());

   //This is required (esp.) for mapping in degrees as fmod returns the decimal part not the 
   //modulous of the number (in mathematical terms).
   if((fabs(modTLX)<EPSILON)||(fabs(fabs(modTLX)-fabs(segmentinfo-> PixelSizeX()))<EPSILON))
   {
      modTLX=0;
   }

   if((fabs(modTLY)<EPSILON)||(fabs(fabs(modTLY)-fabs(segmentinfo-> PixelSizeY()))<EPSILON))
   {
      modTLY=0;
   }


   if(modTLX==0) 
   {
      //Nothing needs to be done for X - already overlays
      if(modTLY==0)
      {      
         //Nothing needs to be done for Y - already overlays
         //return;
      }
      else
      {
         //Need to offset Y
         segmentinfo->UpdateTopLeftY(segmentinfo->TopLeftY()+modTLY);
      }
   }
   else
   {
      //Need to offset X
      segmentinfo->UpdateTopLeftX(segmentinfo->TopLeftX()-modTLX);

      if(modTLY==0)
      {      
         //Nothing needs to be done for Y - already overlays
         //return;
      }
      else
      {
         //Need to offset Y
         segmentinfo->UpdateTopLeftY(segmentinfo->TopLeftY()+modTLY);
      }
   }

   Logger::Verbose("Top Left grid overlay X:"+ToString(ginfo->TopLeftX())+" "+ToString(segmentinfo->TopLeftX())+" "+ToString(modTLX));
   Logger::Verbose("Top Left grid overlay Y:"+ToString(ginfo->TopLeftY())+" "+ToString(segmentinfo->TopLeftY())+" "+ToString(modTLY));
}

#endif

