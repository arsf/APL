//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#include "basic_igm_worker.h"

//-------------------------------------------------------------------------
//Constructor taking IGM file name
//-------------------------------------------------------------------------
Basic_IGM_Worker::Basic_IGM_Worker(std::string fname)
{
   filename=fname;   
   fin=new BinFile(filename);
   //Check that data are in double precision - currently only supported type for igm files
   if(fin->GetDataType() != 5)
   {
      throw "IGM files are currently only supported for double precision (float64) data.";
   }
   minx=1000000000;
   maxx=-100000000;
   miny=1000000000;
   maxy=-1000000000;
   data=NULL;

   //Get the ignore value for bad lat/lon points
   std::string nodata=fin->FromHeader("data ignore value");
   if(nodata.compare("")==0)
   {
      //No "data ignore value" in hdr file - set to max of double
      Logger::Debug("No 'data ignore value' in igm hdr file - setting to maximum double value.");
      nodatavalue=std::numeric_limits<double>::max();
   }
   else
   {
      //convert the no data value string to a double
      Logger::Verbose("Assigning data ignore value as given in igm hdr file: "+nodata);
      nodatavalue=StringToDouble(nodata);
   }

   //Check if min/max is in hdr
   std::string mix=fin->FromHeader(";min x");
   std::string max=fin->FromHeader(";max x");
   std::string miy=fin->FromHeader(";min y");
   std::string may=fin->FromHeader(";max y");
   if((mix.empty())||(max.empty())||(miy.empty())||(may.empty()))
   {
      //If not in hdr, calculate them here
      GetMinMax();
   }
   else
   {
      try
      {
         //Convert min/max to doubles
         minx=(StringToDouble(mix));
         maxx=(StringToDouble(max));
         miny=(StringToDouble(miy));
         maxy=(StringToDouble(may));       
      }
      catch(std::string e)
      {
         throw e+"\n"+"Min / Max values in igm header file do not appear to be 'good' numeric values.";
      }
   }

   //Get projection info
   proj=fin->FromHeader("projection");
   ell=fin->FromHeader("datum ellipsoid");

   //Get file info for speed increases
   nsamples=StringToUINT(fin->FromHeader("samples"));
   nlines=StringToUINT(fin->FromHeader("lines"));
   nbands=StringToUINT(fin->FromHeader("bands"));

   if(nlines==1)
      ISARSF=false;
   else
      ISARSF=true;

   Logger::Debug("Basic igm worker opened. "+ToString(nsamples)+" "+ToString(nlines)+" "+ToString(nbands));
}

//-------------------------------------------------------------------------
// Copy Constructor
//-------------------------------------------------------------------------
Basic_IGM_Worker::Basic_IGM_Worker(Basic_IGM_Worker* b)
{
   filename=b->FileName();
   fin=new BinFile(filename);
   data=NULL;

   minx=b->MinX();
   maxx=b->MaxX();
   miny=b->MinY();
   maxy=b->MaxY();

   //Get projection info
   proj=b->Projection();
   ell=b->Ellipse();

   //Get file info for speed increases
   nsamples=b->Samples();
   nlines=b->Lines();
   nbands=b->Bands();

   //Get the ignore value for bad lat/lon points
   nodatavalue=b->IgnoreValue();

   if(nlines==1)
      ISARSF=false;
   else
      ISARSF=true;
}

//-------------------------------------------------------------------------
//Destructor 
//-------------------------------------------------------------------------
Basic_IGM_Worker::~Basic_IGM_Worker()
{
   if(data!=NULL)
      delete[] data;
   if(fin!=NULL)
      fin->Close();
   Logger::Debug("Basic igm worker destructed.");
}

//-------------------------------------------------------------------------
//Function to get the min/max of the data
//-------------------------------------------------------------------------
void Basic_IGM_Worker::GetMinMax()
{
   //The first band of the IGM file should be x, second band y, third band z (we ignore this one)
   unsigned int nsamples=StringToUINT(fin->FromHeader("samples"));
   unsigned int nbands=StringToUINT(fin->FromHeader("bands"));
   unsigned int nlines=StringToUINT(fin->FromHeader("lines"));

   Logger::Log("Calculating min/max bounds from file data.");

   double* databuffer=new double[nsamples*nbands];
   double tmaxy=0,tmaxx=0;
   double tminy=0,tminx=0;
   //Read in per line of all bands
   for(unsigned int line=0;line<nlines;line++)
   {
      //Read in a line for all bands
      fin->Readline((char*)databuffer,line);
      //Get min/max x for this line
      //Update variables if a new min or max is found
      GetArrayLimits(databuffer,nsamples,tminx,tmaxx,nodatavalue);
      if(tmaxx > maxx)
         maxx=tmaxx;
      if(tminx < minx)
         minx=tminx;

      //Get min/max y for this line
      //Update variables if a new min or max is found
      GetArrayLimits(databuffer+nsamples,nsamples,tminy,tmaxy,nodatavalue);
      if(tmaxy > maxy)
         maxy=tmaxy;
      if(tminy < miny)
         miny=tminy;
   }

   Logger::Debug("GetMinMax: MinX:"+ToString(minx)+" MaxX: "+ToString(maxx)+" MinY: "+ToString(miny)+" MaxY: "+ToString(maxy));
   delete[] databuffer;
}

//-------------------------------------------------------------------------
// Function to get a line of data from the IGM file and store in data array
// Also takes care of along track oversampling and returns a pointer to the 
// data array
//-------------------------------------------------------------------------
double* const Basic_IGM_Worker::GetLine(const unsigned int line)
{
   if(data==NULL)
   {
      data=new double[Samples()*Bands()];
   }
   fin->Readline((char*)data,line);
   return data;
}

//-------------------------------------------------------------------------
// Function to estimate the pixel size from the data in the IGM file
// This assumes a regular ARSF IGM file - e.g. each scan line is a complete 
// across track swath...minimum tests against this so beware if not using
// ARSF igm files
//-------------------------------------------------------------------------
bool Basic_IGM_Worker::GetPixelSize(unsigned int pixelid,double* pixsize)
{
   double emean=0,emin=9999,emax=-9999;
   double nmean=0,nmin=9999,nmax=-9999;
   double diffEal=0,diffNal=0;
   double diffEac=0,diffNac=0;
   double mean_al=0,mean_ac=0;
   
   double* tempdata1=new double[Samples()*Bands()];
   double* tempdata2=new double[Samples()*Bands()];   

   if(ISARSF==false)
   {
      //Only 1 line - this is not a regular ARSF IGM file so we cannot assume
      //what we need to in the below algorithms
      Logger::Log("IGM file is not in 'regular ARSF swath' format so cannot calculate pixel sizes.");
      return BadPixelSizeCalculation(pixsize);
   }   

   for(unsigned int line=0;line<Lines()-1;line++)
   {
      //Read in two lines of data and estimate the values from them
      fin->Readline((char*)tempdata1,line);
      fin->Readline((char*)tempdata2,line+1);
      //E is band 1
      if((tempdata2[pixelid]==IgnoreValue())||(tempdata1[pixelid]==IgnoreValue())
            ||(tempdata2[pixelid + Samples()]==IgnoreValue())||(tempdata1[pixelid+Samples()]==IgnoreValue()))
      {
         //Skip using this line as it contains nodata values in the cells we're using
         continue;
      }

      //Difference between neighbours of consecutive scans (along track)
      diffEal=tempdata2[pixelid]-tempdata1[pixelid];
      //Difference between neighbours within a scan (across track)
      if(pixelid>0)
         diffEac=tempdata2[pixelid]-tempdata2[pixelid-1];
      else
         diffEac=tempdata2[pixelid]-tempdata2[pixelid+1];
      
      //N is band 2
      //Difference between neighbours of consecutive scans (along track)
      diffNal=tempdata2[pixelid + Samples()]-tempdata1[pixelid + Samples()];
      //Difference between neighbours within a scan (across track)
      if(pixelid>0)
         diffNac=tempdata2[pixelid + Samples()]-tempdata2[pixelid-1 + Samples()];
      else
         diffNac=tempdata2[pixelid + Samples()]-tempdata2[pixelid+1 + Samples()];
	 
      //Track the statistics
      mean_al=mean_al+sqrt(diffEal*diffEal+diffNal*diffNal);
      mean_ac=mean_ac+sqrt(diffEac*diffEac+diffNac*diffNac);      
      emean=emean+sqrt(diffEal*diffEal+diffEac*diffEac);
      nmean=nmean+sqrt(diffNal*diffNal+diffNac*diffNac); 

      emin=std::min(emin,sqrt(diffEal*diffEal+diffEac*diffEac));
      nmin=std::min(nmin,sqrt(diffNal*diffNal+diffNac*diffNac));
      emax=std::max(emax,sqrt(diffEal*diffEal+diffEac*diffEac));
      nmax=std::max(nmax,sqrt(diffNal*diffNal+diffNac*diffNac));

   }

   if((mean_al==0)||(mean_ac==0))
      return BadPixelSizeCalculation(pixsize);

   //Fill in the array that is returned
   pixsize[0]=mean_al/(Lines()-1);
   pixsize[1]=mean_ac/(Lines()-1);
   pixsize[2]=emin;
   pixsize[3]=emean/(Lines()-1);
   pixsize[4]=emax;
   pixsize[5]=nmin;
   pixsize[6]=nmean/(Lines()-1);
   pixsize[7]=nmax;

   Logger::Debug("Pixel sizes from IGM at pixel "+ToString(pixelid)+": "+ToString(pixsize[0])+" "+ToString(pixsize[1])+" "+
         ToString(pixsize[2])+" "+ToString(pixsize[3])+" "+ToString(pixsize[4])+" "+ToString(pixsize[5])+" "+ToString(pixsize[6])+" "+
         ToString(pixsize[7]));

   //Tidy up
   delete[] tempdata1;
   delete[] tempdata2;

   return true;
}

bool Basic_IGM_Worker::BadPixelSizeCalculation(double* pixsize)
{
      Logger::Verbose("Unable to calculate pixel separation from IGM data.");
      pixsize[0]=pixsize[1]=pixsize[2]=-1;
      pixsize[3]=pixsize[4]=pixsize[5]=-1;
      pixsize[6]=pixsize[7]=-1;
      return false;
}

//-------------------------------------------------------------------------
// Function to read a specific cell from the IGM and return it as a double
//-------------------------------------------------------------------------
double Basic_IGM_Worker::ReadCell(const unsigned int band,const unsigned int line, const unsigned int col)
{
   return fin->ReadCell(band,line,col);
}

