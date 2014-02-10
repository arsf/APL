//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#include "treegrid_support.h"


//------------------------------------------------------------------------
//Default constructor - assign area to 0
//------------------------------------------------------------------------
Area::Area()
{
   minx=0;
   maxx=0;
   miny=0;
   maxy=0;
}

//------------------------------------------------------------------------
// Copy constructor
//------------------------------------------------------------------------
Area::Area(Area* a)
{
   minx=a->minx;
   maxx=a->maxx;
   miny=a->miny;
   maxy=a->maxy;
}

//------------------------------------------------------------------------
// Constructor taking top left, bottom right coordinates 
//------------------------------------------------------------------------
Area::Area(double mix,double max,double miy,double may)
{
   if((mix>=max)||(miy>=may))
   {
      throw "Cannot create an area whose minimum x/y is greater than or equal to its maximum x/y."+ToString(mix)+" "+
      ToString(max)+" "+ToString(miy)+" "+ToString(may);
   }

   minx=mix;
   maxx=max;
   miny=miy;
   maxy=may;
}

//------------------------------------------------------------------------
// Area Destructor
//------------------------------------------------------------------------
Area::~Area()
{
}


//------------------------------------------------------------------------
// Check if a point (x,y) is within the rectangular area
//------------------------------------------------------------------------
bool Area::Inside(double x,double y)
{
   if((minx<x)&&(x<maxx)&&(miny<y)&&(y<maxy))
      return true;
   else
      return false;
}

//------------------------------------------------------------------------
// Constructor for ItemData when all information is known
//------------------------------------------------------------------------
ItemData::ItemData(double* data,unsigned long int r, unsigned long int c, unsigned long int ns, unsigned long int nr,Basic_IGM_Worker* i)
{
   first=data;
   row=r;
   if(c!=0)
      throw "This cannot be handled properly - first col should always be 0.";
   col=c;
   nsamples=ns;
   nrows=nr;
   if(igm!=NULL)
      delete igm;
   igm=new Basic_IGM_Worker(i);
}

//------------------------------------------------------------------------
// Set the information for ItemData after construction
//------------------------------------------------------------------------
void ItemData::Set(double* data,unsigned long int r, unsigned long int c, unsigned long int ns, unsigned long int nr,Basic_IGM_Worker* i)
{
   first=data;
   row=r;
   if(c!=0)
      throw "This cannot be handled properly - first col should always be 0.";
   col=c;
   nsamples=ns;
   nrows=nr;
   if(igm!=NULL)
      delete igm;
   igm=new Basic_IGM_Worker(i);
}

//------------------------------------------------------------------------
// Set the information for ItemData after construction
//------------------------------------------------------------------------
void ItemData::Set(double* data,unsigned long int r, unsigned long int c, unsigned long int ns, unsigned long int nr,std::string igmfilename)
{
   first=data;
   row=r;
   if(c!=0)
      throw "This cannot be handled properly - first col should always be 0.";
   col=c;
   nsamples=ns;
   nrows=nr;
   if(igm!=NULL)
      delete igm;
   igm=new Basic_IGM_Worker(igmfilename);
}

//------------------------------------------------------------------------
// Construct an empty Item
//------------------------------------------------------------------------
Item::Item()
{
   data=NULL;
}

//------------------------------------------------------------------------
// Construct an Item when all info is known prior
//------------------------------------------------------------------------
Item::Item(ItemData* id,unsigned int r,unsigned int c)
{
   data=id;
   igmrow=r;
   igmcol=c;
}


////------------------------------------------------------------------------
//// Function to write a DatoPoint to a stream
////------------------------------------------------------------------------
//void Item::SaveToDisk(std::ofstream* fout)
//{
//   if(!fout->is_open())
//   {
//      Logger::Log("Error writing Item - stream not open");
//   }
//   else
//   {
//      fout->write((char*)&igmrow,sizeof(unsigned int));
//      fout->write((char*)&igmcol,sizeof(unsigned int));
//   }
//}

////------------------------------------------------------------------------
//// Function to load a data point object from a stream
////------------------------------------------------------------------------
//void Item::LoadFromDisk(std::ifstream* fin)
//{
//   if(!fin->is_open())
//   {
//      Logger::Log("Error reading Item - stream not open");
//   }
//   else
//   {
//      fin->read((char*)&igmrow,sizeof(unsigned int));
//      fin->read((char*)&igmcol,sizeof(unsigned int));
//   }
//}

