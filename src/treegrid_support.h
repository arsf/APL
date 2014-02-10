//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef TREEGRID_SUPPORT_H
#define TREEGRID_SUPPORT_H

#include <fstream>
#include "logger.h"
#include "basic_igm_worker.h"

//-------------------------------------------------------------------------
// Class which holds the information about where to get the item XY data from
//-------------------------------------------------------------------------
class ItemData
{
public:
   ItemData(){first=NULL;nsamples=nrows=row=col=0;igm=NULL;};
   ItemData(double* data,unsigned long int r, unsigned long int c, unsigned long int ns,unsigned long int nr,Basic_IGM_Worker* i);
   ~ItemData(){first=NULL;if(igm!=NULL){delete igm;}}

   void Set(double* data,unsigned long int r, unsigned long int c, unsigned long int ns,unsigned long int nr,Basic_IGM_Worker* i);
   void Set(double* data,unsigned long int r, unsigned long int c, unsigned long int ns, unsigned long int nr,std::string igmfilename);

   inline const double* const First() const{return first;}
   inline const unsigned int NSamples() const {return nsamples;}
   inline const unsigned int NRows() const {return nrows;}
   inline const unsigned int FirstRow() const {return row;}
   inline const unsigned int FirstCol() const {return col;}//nsamples must be the full length of a row of igm data - so this will be 0
   Basic_IGM_Worker* igm;

   bool IsInRAM(long r,long c=0) const //only check rows since a full row of data should be contained i.e. col 0 -> nsamples
   {
      if((r<row) || (r>=row+nrows)||(c<0)||(c>=nsamples)) 
         return false;
      else
         return true;
   };

   inline const double GetX(const long r,const long c)const{if(IsInRAM(r,c)) { return *(first+(r-row)*(2*nsamples) + (c-col));} else {return igm->ReadCell(0,r,c);}}
   inline const double GetY(const long r,const long c)const{if(IsInRAM(r,c)) { return *(first+(r-row)*(2*nsamples) +nsamples + (c-col));} else {return igm->ReadCell(1,r,c);}}

private:
   double* first;
   unsigned int row,col,nsamples,nrows; //nsamples must be the full length of a row of igm data
};

//-------------------------------------------------------------------------
// Class which gets inserted into the TreeGrid 
//-------------------------------------------------------------------------
class Item
{
public:
   Item();
   Item(ItemData* id,unsigned int r,unsigned int c);
   ~Item(){data=NULL;}

   bool operator< (const Item rhs) const {return distance < rhs.distance;}
   //void SaveToDisk(std::ofstream* fout);
   //void LoadFromDisk(std::ifstream* fin);

   //Distance (squared) between this and the last item tested
   float distance;

   //igm and level1 file row/col
   unsigned int igmrow,igmcol; 

   inline double X()
   {
      if((igmrow<data->FirstRow()) || (igmcol<data->FirstCol()) || (igmrow>=data->FirstRow()+data->NRows()) || (igmcol>=data->FirstCol()+data->NSamples()) )
      {
         //Data is not currently in the ItemData data array
         //Need to do something here like read in the data
         Logger::Debug("Reading IGM x value: "+ToString(igmrow)+" "+ToString(igmcol));
         return data->igm->ReadCell(0,igmrow,igmcol);
      }
      else
      {
         //return the X data from the data array
         if(data->First()!=NULL)
            return *(data->First()+(igmrow-data->FirstRow())*(2*data->NSamples()) + (igmcol-data->FirstCol()));
         else
         {
            std::cout<<data->FirstRow()<<" "<<data->FirstCol()<<" "<<data->NRows()<<" "<<data->NSamples()<<" "<<igmrow<<" "<<igmcol<<std::endl;
            throw "Trying to access NULL pointer in Item::X().";
         }
      }
   }

   inline double Y()
   {
      if((igmrow<data->FirstRow()) || (igmcol<data->FirstCol()) || (igmrow>=data->FirstRow()+data->NRows()) || (igmcol>=data->FirstCol()+data->NSamples()) )
      {
         //Data is not currently in the ItemData data array
         //Need to do something here like read in the data
         Logger::Debug("Reading IGM y value: "+ToString(igmrow)+" "+ToString(igmcol));
         return data->igm->ReadCell(1,igmrow,igmcol);
      }
      else
      {
         //return the Y data from the data array
         if(data->First()!=NULL)
            return *(data->First()+(igmrow-data->FirstRow())*(2*data->NSamples()) + data->NSamples() + (igmcol-data->FirstCol()));
         else
         {
            std::cout<<data->FirstRow()<<" "<<data->FirstCol()<<" "<<data->NRows()<<" "<<data->NSamples()<<" "<<igmrow<<" "<<igmcol<<std::endl;
            throw "Trying to access NULL pointer in Item::X().";
         }
      }
   }

   void SetData(ItemData* d){data=d;}
   const ItemData* const GetData()const{return data;}

private:
   ItemData* data;
};

//-------------------------------------------------------------------------
// Defines a bounded area
//-------------------------------------------------------------------------
class Area
{
public:
   Area();
   Area(Area* a);
   Area(double mix,double max,double miy,double may);
   ~Area();

   //Function to test is a point x,y is within the area
   bool Inside(double x,double y);
   
   //Functions to return bounds of area
   const double MaxX()const{return maxx;}
   const double MinX()const{return minx;}
   const double MaxY()const{return maxy;}
   const double MinY()const{return miny;}

private:
   double minx,maxx,miny,maxy;
};

#endif
