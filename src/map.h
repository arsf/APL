//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef MAP_H
#define MAP_H

#include <limits>
#include "binfile.h"
#include "bilwriter.h"
#include "TreeGrid.h"
#include "level3grid.h"
#include "interpolator.h"
#include "linesegment.h"
#include "dataaccessor.h"

//-------------------------------------------------------------------------
//Template function used to copy data from a double buffer into a new buffer
//of any type (e.g. int / float). This is used for the output buffer.
//-------------------------------------------------------------------------
template <class OUTTYPE>
void TransferData(OUTTYPE* newdata,double* olddata,unsigned long int start, unsigned long int end)
{
   for(unsigned long int i=start;i<end;i++)
   {
      newdata[i]=static_cast<OUTTYPE>(olddata[i]);
   }
}

//-------------------------------------------------------------------------
// Abstract base class for Map - for polymorphism
//-------------------------------------------------------------------------
class AbstractMap
{
public:
   virtual ~AbstractMap(){};
   virtual void MapLineSegments(TreeGrid* tg,std::string igmfilename,std::string level1filename)=0;
   //virtual void SetMaxInterpolationDistance(double mid){};
   virtual void AssignProjection(std::string proj)=0;
   virtual unsigned int GetOutputDataSize()=0;

   Level3Grid* grid;
};


//-------------------------------------------------------------------------
// Map class for doing all the main processing - now a template class
// so that it can map data of different types from level 1 files 
//-------------------------------------------------------------------------
template<class T>
class Map : public AbstractMap
{
public:
   Map();
   Map(std::string outfname,double Xpixelsize,double Ypixelsize, std::string strBandList,Area* output_area,
            std::string lev1fname,Interpolators::InterpolatorType itype,int npoints,uint64_t buffsize,
            BILWriter::DataType datatype,std::string strRowColMappingFilename,T ndv);
   ~Map();

   virtual void AssignProjection(std::string proj);
   virtual void MapLineSegments(TreeGrid* tg,std::string igmfilename,std::string level1filename);
   virtual void SetMaxInterpolationDistance(double mid)const{interpolator->SetMaxInterpDistance(mid);}

   virtual unsigned int GetOutputDataSize(){return writer->GetDataSize();}

   virtual void SetInterpolatorIgnoreValue(double ig){interpolator->SetIgnoreValue(ig);}
   virtual void SetInterpolatorIgnoreFlag(bool f){interpolator->SetIgnoreFlag(f);}
   virtual void SetInterpolatorNoDataValue(T ndv){interpolator->SetNoDataValue(ndv);}
   virtual void SetNoDataValue(T val);

private:
   double* buffer;
   uint64_t length_of_buffer;

   //Number of nearest points to find from tree grid
   int numpoints;
   FILE* boundedwriter;
   uint64_t maximum_segment_memory;

   char* outputbuffer;

   BILWriter* writer;
   Interpolator<T>* interpolator;

   void FillPixel(long int row, long int col, std::vector<Item>* points,DataAccessor<T>* lev1data);
   unsigned long WriteBuffer();
   void WriteBuffer(int* bounds,Level3GridInfo* segment,int row);
   void TransformDataType(unsigned int datatype,unsigned long int start,unsigned long int end);
   void SetSegmentSize(uint64_t bs){maximum_segment_memory=bs;}
   void AssignWavelengths(std::string level1fname);

   //Writers for row/col level1 mapping
   FILE* boundedl1mappingwriter;
   BILWriter* l1mappingwriter;
   int* l1mapping_rows;
   int* l1mapping_cols;
   void UpdateL1MappingHeader();
   //value to insert into map as nodata
   T nodatavalue;

};


//----------------------------------------------------------------------
// Default constructor
//----------------------------------------------------------------------
template<class T>
Map<T>::Map()
{
   Logger::Verbose("Constructing Map.");
   grid=NULL;
   writer=NULL;
   interpolator=NULL;
   buffer=NULL;
   length_of_buffer=0;
   outputbuffer=NULL;
   l1mappingwriter=NULL;
   boundedl1mappingwriter=NULL;
   l1mapping_rows=NULL;
   l1mapping_cols=NULL;
   nodatavalue=0;
}

//----------------------------------------------------------------------
// Destructor
//----------------------------------------------------------------------
template<class T>
Map<T>::~Map()
{
   if(grid!=NULL)
      delete grid;
   if(writer!=NULL)
   {
      writer->Close();
      delete writer;
   }
   if(interpolator!=NULL)
      delete interpolator;
   if(buffer!=NULL)
      delete[] buffer;
   if(outputbuffer!=NULL)
      delete[] outputbuffer;
   if(boundedwriter!=NULL)
      fclose(boundedwriter);

   if(l1mappingwriter!=NULL)
   {
      l1mappingwriter->Close();
      delete l1mappingwriter;
   }
   if(boundedl1mappingwriter!=NULL)
      fclose(boundedl1mappingwriter);
   if(l1mapping_rows!=NULL)
      delete[] l1mapping_rows;
   if(l1mapping_cols!=NULL)
      delete[] l1mapping_cols;


}

//----------------------------------------------------------------------
// Constructor taking output filename and level 3 grid specifications
//----------------------------------------------------------------------
template<class T>
Map<T>::Map(std::string outfname,double Xpixelsize,double Ypixelsize, std::string strBandList,Area* output_area,std::string lev1fname,
            Interpolators::InterpolatorType itype,int npoints,uint64_t buffsize,BILWriter::DataType datatype,std::string strRowColMappingFilename,T ndv=0)
{
   Logger::Verbose("Constructing Map.");
   Logger::Verbose("Building Map Level 3 grid ... ");
   grid=new Level3Grid(Xpixelsize,Ypixelsize,strBandList,output_area);
   Logger::Verbose("Building Map data writer ... ");
   writer=new BILWriter(outfname,datatype,grid->NumRows(),grid->NumCols(),grid->NumBands(),'w');
   boundedwriter=fopen(outfname.c_str(), "wb" );
   if(boundedwriter==NULL)
   {
      throw "Failed to open boundedwriter in Map."+outfname;
   }

   Logger::Verbose("Building Map interpolator ... ");
   if(itype==Interpolators::NEARESTNEIGHBOUR)
   {
      interpolator=new NearestNeighbour<T>(grid->NumBands());
      if(npoints !=1)
         Logger::Log("Nearest neighbour only uses 1 point - overriding given number of points.");
      numpoints=1;
   }
   else if(itype==Interpolators::IDW)
   {
      interpolator=new IDW<T>(grid->NumBands());
      numpoints=npoints;
   }
   else if(itype==Interpolators::BILINEAR)
   {
      interpolator=new Bilinear<T>(grid->NumBands());
      if(npoints !=10)
         Logger::Log("Bilinear uses a fixed number of points - overriding given number of points.");
      //numpoints=1;
      numpoints=npoints;
   }
   else if(itype==Interpolators::BILINEARLEVEL3)
   {
      interpolator=new BilinearLevel3<T>(grid->NumBands());
      if(npoints !=1)
         Logger::Log("BilinearLevel3 uses a fixed number of points - overriding given number of points.");
      numpoints=1;
   }
   else if(itype==Interpolators::CUBIC)
   {
      interpolator=new Cubic<T>(grid->NumBands());
      if(npoints !=1)
         Logger::Log("Cubic uses a fixed number of points - overriding given number of points.");
      numpoints=4;
   }
   else
      throw "Unknown interpolator type in Map construction.";

   interpolator->SetNumPoints(numpoints);

   //Get the wavelengths of the data
   AssignWavelengths(lev1fname);

   //Set the segment buffer size
   SetSegmentSize(buffsize);

   //Set the no data value
   SetNoDataValue(ndv);

   try
   {
      length_of_buffer=grid->NumCols()*grid->NumBands();
      buffer=new double[length_of_buffer];
      outputbuffer=new char[length_of_buffer*writer->GetDataSize()];
      //Assign buffer values of nodatavalue
      for(unsigned int i=0;i<length_of_buffer;i++)
         buffer[i]=nodatavalue;

      for(unsigned int i=0;i<length_of_buffer*writer->GetDataSize();i++)
         outputbuffer[i]=0;

      Logger::Verbose("Created a map buffer of size: "+ToString(length_of_buffer));
   }
   catch(std::bad_alloc& ba)
   {
      //Write out error and exit
      Logger::Log(ba.what());
      Logger::Error("Failed to allocate RAM for a single line of output grid. Will have to write per pixel - this will be slower. Not YET IMPLEMENTED");
      exit(1);
   }

   //If level-1 row/col mapping is to be written out - open the writers here
   //and create the arrays to hold the data
   if(strRowColMappingFilename.compare("")!=0)
   {
      l1mappingwriter=new BILWriter(strRowColMappingFilename,BILWriter::int32,grid->NumRows(),grid->NumCols(),2,'w');
      boundedl1mappingwriter=fopen(strRowColMappingFilename.c_str(), "wb" );
      if(boundedl1mappingwriter==NULL)
      {
         throw "Failed to open boundedl1mappingwriter in Map."+strRowColMappingFilename;
      }

      l1mapping_rows=new int[grid->NumCols()];
      l1mapping_cols=new int[grid->NumCols()];

      for(unsigned int i=0;i<grid->NumCols();i++)
      {
         l1mapping_rows[i]=-1;
         l1mapping_cols[i]=-1;
      }

      UpdateL1MappingHeader();
   }
   else
   {
      l1mappingwriter=NULL;
      boundedl1mappingwriter=NULL;
      l1mapping_rows=NULL;
      l1mapping_cols=NULL;  
   }
}

//----------------------------------------------------------------------
// Assign the no data value
//----------------------------------------------------------------------
template<class T>
void Map<T>::SetNoDataValue(T val)
{
   nodatavalue=val;
   SetInterpolatorNoDataValue(nodatavalue);
   writer->AddToHdr("data ignore value = "+ToString(nodatavalue));
}

//----------------------------------------------------------------------
// Assign the band wavelengths to the level 3 data writer
//----------------------------------------------------------------------
template<class T>
void Map<T>::AssignWavelengths(std::string lev1fname)
{
   unsigned int ba=0;
   BinFile lev1(lev1fname);
   for(unsigned int b=0;b<grid->NumBands();b++)
   {
      ba=grid->Bands()[b];
      //get the wavelength of the data for this band
      grid->AddWavelength(lev1.FromHeader("wavelength",ba));
   }
   writer->AddToHdr("band names = {"+grid->Wavelengths()+"}");
   //Also copy over the wavelengths too incase further envi analysis
   writer->AddToHdr("wavelength = {"+grid->Wavelengths()+"}");
}

//----------------------------------------------------------------------
// Assign the projection information to the level 3 data writer
//----------------------------------------------------------------------
template<class T>
void Map<T>::AssignProjection(std::string proj)
{
   writer->AddToHdr("map info = "+grid->MapInfo(proj));   
   //Also add it to the l1mapping file if created
   if(l1mappingwriter!=NULL)
      l1mappingwriter->AddToHdr("map info = "+grid->MapInfo(proj));   
}

//----------------------------------------------------------------------
// Update the l1mapping row/col file to have some more header information
//----------------------------------------------------------------------
template<class T>
void Map<T>::UpdateL1MappingHeader()
{
   l1mappingwriter->AddToHdr(";This file contains the level-1 row/column identifiers that were used to fill in the mapped file.");
   l1mappingwriter->AddToHdr("band names = {row value, column value}");   
}

//----------------------------------------------------------------------
// Function to add an interpolated value to the map buffer
//----------------------------------------------------------------------
template<class T>
void Map<T>::FillPixel(long int row, long int col, std::vector<Item>* points,DataAccessor<T>* lev1data=NULL)
{
   if((row)!=0)
      throw "Row should be NULL until a multiple row buffer writing method is implemented";

   //If column is -ve it is not in required region so no point mapping it
   if((col < 0))
      return;

   //return if no points found or if closest point is further than max interpolation distance, unless BilinearLevel3 used.
   if( (points==NULL)|| //OR
       (points->size()==0)|| //OR
       (((*points)[0].distance>interpolator->GetMaxInterpDistanceSq())&&(interpolator->interpolator_type!=Interpolators::BILINEARLEVEL3)) )
   {
      return;
   }
   else
   {
      //Interpolate the value
      interpolator->Interpolate(points,grid->Bands(),lev1data);
      //Get the interpolated value
      const double* const dp=interpolator->Data();
      //Assign the interpolated value (for each band)
      for(unsigned int b=0;b<grid->NumBands();b++)
      {
         buffer[b*grid->NumCols()+(col)]=dp[b];
      }
   }
}

//----------------------------------------------------------------------
// Function to write out the contents of the buffer to the output file
//----------------------------------------------------------------------
template<class T>
unsigned long Map<T>::WriteBuffer()
{
   //Added static integer to keep track of number of lines written (i.e. number 
   //of times that this function has been called). This is then used to check
   //that the correct number of lines are written out during the processing of
   //segment 0
   static unsigned long lineswritten=0;

   TransformDataType(writer->GetDataType(),0,length_of_buffer);

   writer->WriteLine(outputbuffer);
   //Clear buffer here - then assume buffer has values of NODATA for speed in FillPixel
   for(unsigned int i=0;i<length_of_buffer;i++)
      buffer[i]=nodatavalue;

   for(unsigned int i=0;i<length_of_buffer*writer->GetDataSize();i++)
      outputbuffer[i]=0;

   //If the level1 mapping is requested then write it out also
   if(l1mappingwriter!=NULL)
   {
      l1mappingwriter->WriteBandLine((char*)l1mapping_rows);
      l1mappingwriter->WriteBandLine((char*)l1mapping_cols);
      //Reset arrays for next line
      for(unsigned int i=0;i<grid->NumCols();i++)
      {
         l1mapping_rows[i]=-1;
         l1mapping_cols[i]=-1;
      }
   }

   lineswritten++;
   return lineswritten;
}

//----------------------------------------------------------------------
// Function to transform the double* data buffer and store the data
// in the char* outputbuffer as a different data type as requested
//----------------------------------------------------------------------
template<class T>
void Map<T>::TransformDataType(unsigned int datatype,unsigned long int start,unsigned long int end)
{  
   //Using the datatypes as defined by ENVI and used in the bilwriter
   switch(datatype)
   {
      case 1://char
         TransferData<char>((char*)outputbuffer,buffer,start,end);
         break;
      case 2://signed 16-bit int
         TransferData<int16_t>((int16_t*)outputbuffer,buffer,start,end);
         break;
      case 3://signed 32 bit int
         TransferData<int32_t>((int32_t*)outputbuffer,buffer,start,end);
         break;
      case 4:
         TransferData<float>((float*)outputbuffer,buffer,start,end);
         break;
      case 5:
         TransferData<double>((double*)outputbuffer,buffer,start,end);
         break;
      case 12:
         TransferData<uint16_t>((uint16_t*)outputbuffer,buffer,start,end);
         break;
      case 13:
         TransferData<uint32_t>((uint32_t*)outputbuffer,buffer,start,end);
         break;
      default:
         break;
   }
}

//----------------------------------------------------------------------
// Function to perform the mapping by splitting the line up into segments
//----------------------------------------------------------------------
template<class T>
void Map<T>::MapLineSegments(TreeGrid* tg,std::string igmfilename,std::string level1filename)
{
   //hold the bounds of columns to process for each row
   int bounds[2]={0,0};
   //Pointer to returned vector of points 
   std::vector<Item>* dp=NULL;
   L3Point rc(0,0);
   IGMPoint xy(0,0);
   //Percent counter
   float perccount=0;
   int upperbound=0,lowerbound=0;

   if(numpoints <=0)
      throw "Number of points for TreeGrid search cannot be less than 1.";
   if(numpoints > 10)
      Logger::Log("Requested number of near points to find is large - will continue but may be slow.");

   //*******************************************
   //** Set up interpolator for mask interpolation

   //Set the tree grid to the interpolator so that it can search to fill in mask holes
   interpolator->tg=tg;
   double searchradius=0;
   searchradius=sqrt(interpolator->GetMaxInterpDistanceSq());

   interpolator->SetSearchRadius(searchradius);
   //**
   //*******************************************


   unsigned long number_of_lines_written=0;
   uint64_t segmentsize=std::numeric_limits<uint64_t>::max();
   int nsegments=0;
   int* segbounds=NULL;
   bool lowerboundfound=false;
   std::vector<int> colbounds;

   //Need to calculate the segments here that are to be processed   
   Logger::Verbose("Calculating the number of segments to split the processing up into - based on maximum allowed memory buffer size.");
   Basic_IGM_Worker igmw(igmfilename);
 
   //Calculate the amount of RAM to store the entire required data in memory and store in fullsize (cast to enforce uint64_t)
   uint64_t fullsize=(static_cast<uint64_t>(grid->NumBands())*igmw.Samples()*igmw.Lines()*sizeof(T)) + (static_cast<uint64_t>(2*igmw.Samples()*igmw.Lines()*sizeof(double)));
   Logger::Debug("Amount of RAM required for processing in a single segment (bytes): "+ToString(fullsize));
   if(fullsize <  maximum_segment_memory)
   {
      //We can process in one segment
      nsegments=1;
      segbounds=new int[nsegments+1];
   }
   else
   {
      //We need to process in multiple segments
      int ncount=1;
      while( segmentsize > maximum_segment_memory)
      {
         ncount++;
         segmentsize=fullsize/ncount;
      }
      nsegments=ncount;
      segbounds=new int[nsegments+1];
   }
   Logger::Verbose("Splitting up into "+ToString(nsegments)+" segments.");
   Logger::Debug("Segment bounds: ");
   for(int i=0;i<=nsegments;i++)
   {
      segbounds[i]=i*igmw.Lines()/nsegments;
      Logger::Debug("   "+ToString(segbounds[i]));
   }
   
   LineSegment<T>* linesegment=NULL;
   //Loop through for each segment mapping it in turn
   for(int seg=0;seg<nsegments;seg++)
   {
      //Try and create a line segment
      try
      {
         unsigned int segment_buffer=10;
         if(seg==0)
            linesegment=new LineSegment<T>(segbounds[seg],segbounds[seg+1],segment_buffer,grid->PixelSizeX(),grid->PixelSizeY(),grid->GetGridInfo()->BandList(),igmfilename,level1filename,grid->GetGridInfo()->GetBounds());
         else //start from the previous line so as to avoid discontinuity or missing data between segments
            linesegment=new LineSegment<T>(segbounds[seg]-1,segbounds[seg+1],segment_buffer,grid->PixelSizeX(),grid->PixelSizeY(),grid->GetGridInfo()->BandList(),igmfilename,level1filename,grid->GetGridInfo()->GetBounds());
      }
      catch(std::string e)
      {
         //Something went wrong - check if it is the OK exception of being outside given region
         if(e.compare("LineSegment not created as it is outside given region")==0)
         {
            //No part of segment rectangle intersects with level3grid rectangle
            if(seg==0)
            {
               //We need to write out the full file of 0's before skipping
               //as following segments expect the file to exist
               for(unsigned int i=0;i<grid->NumRows();i++)
               {
                  number_of_lines_written=WriteBuffer();
               }
            }
            Logger::Verbose("Skipping this segment - it is not required for the specified output region.");
            //Add a percentage counter message to the terminal to update the counter for this segment
            perccount+=100.0/nsegments;
            Logger::Log("Approximate percent complete: "+ToString((int)perccount));
            //Skip onto next segment
            continue;
         }
         else
            throw e;
      }

      //Offset the line segment such that it overlays the level3Grid lines exactly
      linesegment->OffsetToGrid(grid->GetGridInfo());

      //Apply the line segment IGM data to the treegrid
      tg->itemdata.Set(linesegment->igm->Data(),linesegment->igm->FirstRow(),0,linesegment->igm->Samples(),linesegment->igm->Lines(),igmfilename);

      //Create a data accessor object - this is for the interpolation and tree grid searches when it needs to know the data
      //values and compare against ignore values. This object is an interface to the binary file and linesegment block data.
      DataAccessor<T>* da=new DataAccessor<T>(linesegment->level1,level1filename,const_cast<unsigned int*>(grid->Bands()));
      DataAccessor<T>* dummy=NULL; //Annoying NULL pointer of type T needed as cannot just pass "NULL" in function call later on

      //Output a chunk of 0 data if this first segment is not at the very top of the level3grid and is within the level3grid
      if((seg==0)&&(linesegment->segmentinfo->TopLeftY()<grid->TopLeftY()))
      {
         Logger::Log("Outputting zero'ed buffer data until reach first image data ...");
         if(linesegment->segmentinfo->TopLeftY()>grid->BottomRightY())
         {
            //Output some lines of 0's upto the first line of the segment
            for(int i=0;i<(int)((grid->TopLeftY() - linesegment->segmentinfo->TopLeftY())/grid->PixelSizeY());i++)
            {
               number_of_lines_written=WriteBuffer();
            }
         }
         else
         {
            //Output some lines of 0's upto the end of the level3grid
            for(int i=0;i<(int)((grid->TopLeftY()-grid->BottomRightY())/grid->PixelSizeY());i++)
            {
               number_of_lines_written=WriteBuffer();
            }
         }
      }

      //Buffer offset for outputting segment 0 to shift the buffer from segment position to level3grid position in fillgrid
      int buffoffset=(int)((linesegment->segmentinfo->TopLeftX()-grid->TopLeftX())/grid->PixelSizeX());

      //Need a row offset to add onto segment to get into level3 grid row values
      int rowoffset=(int)((grid->TopLeftY()-linesegment->segmentinfo->TopLeftY())/grid->PixelSizeY());
      Logger::Debug("Using a row offset of "+ToString(rowoffset)+" and a buffoffset of "+ToString(buffoffset)); 

      //linesegment->outline->WriteEdge("edges.txt");

      //For each row of the segment grid
      for(unsigned int row=0;row<linesegment->segmentinfo->NumRows();row++)
      {
         //Bool to test if the lower bound in colbound has been set or not in the lowerbound variable
         lowerboundfound=false;
         //Assign the row value to rc such that it is the level 3 grid row - not segment row
         rc.row=(int)row+rowoffset;
         //Update the percent complete counter here - this is incase rows are being skipped below
         //means the percentage counter is still updated.
         if(row % (unsigned int)(linesegment->segmentinfo->NumRows()/10.0) == (unsigned int)(linesegment->segmentinfo->NumRows()/10.0)-1)
         {
            perccount+=10.0/nsegments;
            Logger::Log("Approximate percent complete: "+ToString((int)perccount));
         }
         //Skip if out of bounds of the level3 grid
         if((rc.row < 0)||(rc.row >= (long)grid->NumRows()))
            continue;

         //Fill in the array with the edge points
         linesegment->outline->GetEdgeIntersectsOfRow(row,colbounds);

         //Check that there is an even number of bounds
         if(colbounds.size()%2!=0)
         {
            Logger::Warning("Problem using internal column bounds calculation - falling back to slower mapping of each column for this row: "+ToString(rc.row));
            colbounds.clear();
            colbounds[0]=0;
            colbounds[1]=grid->NumCols()-1;
         }

//FIXME
//Commented out 3 lines below and added bounds[0]/[1] lines to revert back to old method of choosing which pixels
//to actually process. This is because of problems when running through data that overlaps intself, c.f. spiral flight line.
         //Loop through the edge points to process only between these bounds
//         for(std::vector<int>::iterator it=colbounds.begin();it<colbounds.end();it=it+2)

if(colbounds.size()!=0)
         {
            //Set the bounds to be the next 2 in the colbounds for this row
//            bounds[0]=(*it);//lowerbound
//            bounds[1]=(*(it+1));//upperbound
bounds[0]=*(colbounds.begin());
bounds[1]=*(colbounds.end()-1);

            //Add on the buffoffset to offset them to the level3 grid col values
            bounds[0]+=buffoffset;
            bounds[1]+=buffoffset;

            //Check if they are in limits and update appropriately
            if(bounds[0] < 0)
               bounds[0]=0;
            if(bounds[1] >= (long)grid->NumCols())
               bounds[1]=grid->NumCols()-1;

            //Set the values to -1 for "impossible" cases - this signifies we don't want to map the row, but still output a row of 0's
            if ((bounds[0] > bounds[1]) || (bounds[1] < 0))
            {
               //Dont think we need the -1 case anymore so continue-ing instead
               //bounds[0]=bounds[1]=-1;
               continue;
            }

            Logger::Debug("Mapping grid row: "+ToString(rc.row)+" between columns: "+ToString(bounds[0])+" "+ToString(bounds[1]));

            //Process the data for these columns and row
            for(int col=bounds[0];col<=bounds[1];col++)
            {
               //If -ve then skip this column value
               if(col < 0)
                  continue;

               //Assign the column value to rc
               rc.col=col;

               //Convert the rc to xy 
               grid->ConvertRC2XY(&rc,&xy); 

               if((interpolator->interpolator_type==Interpolators::BILINEARLEVEL3)
               ||(interpolator->interpolator_type==Interpolators::CUBIC))
               {
                  //Get the nearest 4*numpoints items that create a quad surrounding xy searching
                  //upto a distance of searchradius
                  dp=tg->GetQuadItems(numpoints,&xy,searchradius,dummy,0,0); 
               }
               else
               {
                  //Get the nearest numpoints items to xy searching upto searchradius distance
                  //the last 2 values can be anything as the level1 file is set to NULL
                  dp=tg->GetNearestXItems(numpoints,&xy,searchradius,dummy,0,0); 
               }

               //If level1 rowcol mapping is to be output and we are using nearest neighbour interpolation
               if((l1mappingwriter!=NULL)&&(interpolator->interpolator_type==Interpolators::NEARESTNEIGHBOUR)&&(dp!=NULL))
               {
                  //Store the values for output at end of row processing
                  l1mapping_rows[col]=(*dp)[0].igmrow; //append the row
                  l1mapping_cols[col]=(*dp)[0].igmcol; //append the column
               }

               //Only insert a value into the buffer if col < number of columns of grid
               if(col < (long)grid->NumCols())
               {
                  interpolator->SetL3Pos(&xy); //the interpolator may need to know the lev3 pixel position in X,Y
                  FillPixel(0,col,dp,da);
               }
            }

            //We need to keep track of lowest bound and upper most bound for the writing out
            //If this is the first iteration of loop - we can get the lower bound
            if(lowerboundfound==false)
            {
               lowerbound=bounds[0];
               lowerboundfound=true;
            }

            //Upperbound is simply the bound[1] after the last loop
            //However - since the last loop may be before the last colbound, we'll just update this everytime in the loop
            upperbound=bounds[1];   
         }

         Logger::Debug("Lower and upper bounds used for writing data: "+ToString(lowerbound)+" "+ToString(upperbound));

         //Write out the line of data currently in the buffer
         if(seg==0)//just write out the full buffer
         {
            if((linesegment->segmentinfo->TopLeftY() - row*grid->PixelSizeY() <= grid->TopLeftY())
               &&(linesegment->segmentinfo->TopLeftY() - row*grid->PixelSizeY() >= grid->BottomRightY()))
            {
               number_of_lines_written=WriteBuffer();
            }
         }
         else//write out sections of the buffer - beware bounds are changed by this function
         {
            //need to set bounds to min/max here - else implement a loop to write data for each section of bounds.
            //assuming buffer is always reset to 0s then it is fine to output between min/max bounds.
            bounds[0]=lowerbound;
            bounds[1]=upperbound;//bounds[1] is already equal to this but put here to make it clearer whats going on

            WriteBuffer(bounds,linesegment->segmentinfo,rc.row);
         }
         //Clear the vector to ensure all used bounds are removed
         colbounds.clear();
         //reset bounds
         lowerbound=upperbound=-1;
      }

      //Output a chunk of 0 data if this first segment is not at the very bottom of the level3grid but within the grid
      if((seg==0)&&(linesegment->segmentinfo->BottomRightY()>grid->BottomRightY())&&(linesegment->segmentinfo->BottomRightY()<grid->TopLeftY()))
      {
         //Output some lines of 0's upto the last line of the level3grid
         //Commented out for loop and replaced with the while loop - this is because
         //sometimes the for loop gives 1 line too few and so the output file is missing
         //the final line. The 'number_of_lines_written' should keep better track of the amount
         //of lines to write out here.

//         for(int i=0;i<(int)((linesegment->segmentinfo->BottomRightY()-grid->BottomRightY())/grid->PixelSizeY());i++)
//         {
//            WriteBuffer();
//         }

         while(number_of_lines_written < grid->NumRows())
            number_of_lines_written=WriteBuffer();

      }
      delete da;
      delete linesegment;      
   }
   delete[] segbounds;
}


//----------------------------------------------------------------------
// Function to write out the contents of the buffer to the output file
// but only in the section described by bounds
//----------------------------------------------------------------------
template<class T>
void Map<T>::WriteBuffer(int* bounds,Level3GridInfo* segment,int thisrow)
{
   //Suppose we don't know if the segment is contained within the level3grid
   //for example, if we are mapping a segment which is outside of the area of interest

   //We need to find any overlapping part and only write that bit.
   //Find the overlap between map->grid and segment->grid.

   if((segment->BottomRightY() > grid->TopLeftY()) || (segment->TopLeftY() < grid->BottomRightY()) 
            || (segment->TopLeftX() > grid->BottomRightX()) || (segment->BottomRightX() < grid->TopLeftX()))
   {
      return; //line segment does not intersect with level 3 grid
   }

   //Number of rows to skip through file to reach the start of this row to write
   uint64_t numrowstoskip=thisrow; 

   if((numrowstoskip < 0) || (numrowstoskip > grid->NumRows()))
   {
      return; //row of segment is outside of level3 grid
   }
   
   double firstX=grid->TopLeftX() + bounds[0]*grid->PixelSizeX();
   double lastX=grid->TopLeftX() + bounds[1]*grid->PixelSizeX();
   Logger::Debug("First and last X for output: "+ToString(firstX)+" "+ToString(lastX));

   if((lastX < grid->TopLeftX()) || (firstX > grid->BottomRightX()))
   {
      return; //column bounds for this row are outside the level 3 grid
   }

   if(firstX < grid->TopLeftX())
      bounds[0]=0;

   if(lastX > grid->BottomRightX())
      bounds[1]=grid->NumCols()-1;

   uint64_t startofboundedbuffer=numrowstoskip*grid->NumCols()*grid->NumBands() + bounds[0];
   int boundedlength=bounds[1]-bounds[0] +1; //+1 because bounds are inclusive
   Logger::Debug("Writing data of length: "+ToString(boundedlength)+" at point "+ToString(startofboundedbuffer));

   if(boundedwriter==NULL)
   {
      throw "Bounded writer is not open!?!";
   }

   unsigned long datastartcell=0;
   fseeko64(boundedwriter,startofboundedbuffer*writer->GetDataSize(),SEEK_SET);
   for(unsigned int b=0;b<grid->NumBands();b++)
   {
      datastartcell=(grid->NumCols()*b + bounds[0]);
      TransformDataType(writer->GetDataType(),datastartcell,(datastartcell+boundedlength));
      fwrite(&(outputbuffer[datastartcell*writer->GetDataSize()]),writer->GetDataSize(),boundedlength,boundedwriter);
      fseeko64(boundedwriter,(grid->NumCols()-boundedlength)*writer->GetDataSize(),SEEK_CUR);
   }

   if((ferror(boundedwriter))||(feof(boundedwriter)))
      Logger::Log("Error writing to bounded file.");

   //If the level1 mapping is requested then write it out also
   if(boundedl1mappingwriter!=NULL)
   {
      int mappingdatasize=sizeof(int);
      fseeko64(boundedl1mappingwriter,startofboundedbuffer*mappingdatasize,SEEK_SET);
      fwrite(&(l1mapping_rows[bounds[0]*mappingdatasize]),mappingdatasize,boundedlength,boundedl1mappingwriter);
      fseeko64(boundedl1mappingwriter,(grid->NumCols()-boundedlength)*mappingdatasize,SEEK_CUR);
      fwrite(&(l1mapping_cols[bounds[0]*mappingdatasize]),mappingdatasize,boundedlength,boundedl1mappingwriter);

      //Clear the row/cols buffer arrays for next loop
      for(unsigned int i=0;i<grid->NumCols();i++)
      {
         l1mapping_rows[i]=-1;
         l1mapping_cols[i]=-1;
      }    
   }

   //Clear buffer here - then assume buffer has values of 0 for speed in FillPixel
   for(unsigned int i=0;i<length_of_buffer;i++)
   {
      buffer[i]=nodatavalue;
      //outputbuffer[i]=0;
   }

   for(unsigned int i=0;i<length_of_buffer*writer->GetDataSize();i++)
      outputbuffer[i]=0;
   
}

#endif


