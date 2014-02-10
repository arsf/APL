//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#include "TreeGrid.h"

//-------------------------------------------------------------------------
// Default TreeGrid Constructor for an empty treegrid
//-------------------------------------------------------------------------
TreeGrid::TreeGrid()
{
   Logger::Verbose("Constructing TreeGrid.");
   rows=0;
   cols=0;
   sizeX=sizeY=0;
   topLeftX=topLeftY=0;
   bottomRightX=0;
   bottomRightY=0;
   collection=NULL;
   NullCollection=NULL;
}

//-------------------------------------------------------------------------
// TreeGrid Constructor for when number of rows/cols known
//-------------------------------------------------------------------------
TreeGrid::TreeGrid(unsigned int r, unsigned int c,Ellipsoid* ell)
{
   Logger::Verbose("Constructing TreeGrid.");
   rows=r;
   cols=c;
   sizeX=sizeY=0;
   topLeftX=topLeftY=0;
   bottomRightX=0;
   bottomRightY=0;
   ellipse=ell;
   //Create an empty collection
   NullCollection=new Collection();

   //Create the grid
   collection=new Collection**[rows];
   for(unsigned int i=0;i<rows;i++)
   {
      collection[i]=new Collection*[cols];
   }
}

//-------------------------------------------------------------------------
// TreeGrid Destructor
//-------------------------------------------------------------------------
TreeGrid::~TreeGrid()
{
   Logger::Verbose("Destructing TreeGrid.");
   for(unsigned int i=0;i<rows;i++)
   {
      for(unsigned int j=0;j<cols;j++)
      {
         if((((collection[i][j])) != NULL)   &&(collection[i][j]!=NullCollection))
            delete ((collection[i][j]));
      }
      delete[] collection[i];
   }

   if(collection!=NULL)
      delete[] collection;   

   if(NullCollection!=NULL)
      delete NullCollection;
}

//-------------------------------------------------------------------------
// Function to set the treegrid up and initialise with size and top left
//-------------------------------------------------------------------------
void TreeGrid::SetUpGrid(double sX,double sY,double tlX,double tlY)
{
   //Set up the regular grid
   topLeftX=tlX;
   topLeftY=tlY;
   sizeX=sX;
   sizeY=sY;

   if(NullCollection==NULL)
   {
      //Create an empty collection
      NullCollection=new Collection();
   }

   for(unsigned long int r=0;r<rows;r++)
   {
      for(unsigned long int c=0;c<cols;c++)
      {
         collection[r][c]=NullCollection;        
      }
   }
}

//-------------------------------------------------------------------------
// Function to set the treegrid up and initialise fully
//-------------------------------------------------------------------------
void TreeGrid::SetUpGrid(unsigned int r, unsigned int c,double sX,double sY,double tlX,double tlY,double brX,double brY)
{
   //Set up the regular grid if no collection already
   if(collection != NULL)
   {
      throw "Cannot set up collection more than once.";
   }
   else
   {
      rows=r;
      cols=c;
      //Create the grid
      collection=new Collection**[rows];
      for(unsigned int i=0;i<rows;i++)
         collection[i]=new Collection*[cols];
   }

   topLeftX=tlX;
   topLeftY=tlY;
   bottomRightX=brX;
   bottomRightY=brY;
   sizeX=sX;
   sizeY=sY;

   if(NullCollection==NULL)
   {
      //Create an empty collection
      NullCollection=new Collection();
   }

   for(unsigned long int r=0;r<rows;r++)
   {
      for(unsigned long int c=0;c<cols;c++)
      {
         collection[r][c]=NullCollection;        
      }
   }
}

//-------------------------------------------------------------------------
// Function to insert item data into the treegrid
//-------------------------------------------------------------------------
void TreeGrid::InsertData(Item* data,unsigned int nitems)
{
   unsigned int to_col=0,to_row=0;

   double posX=0,posY=0;
   for(unsigned int i=0;i<nitems;i++)
   {
      posX=(data[i].X()-topLeftX)/sizeX;
      posY=(topLeftY-data[i].Y())/sizeY;
      to_col=static_cast<unsigned int>(floor(posX));
      to_row=static_cast<unsigned int>(floor(posY));

      if(collection[to_row][to_col]==NullCollection)
         collection[to_row][to_col]=InitialiseNewCollection(to_row,to_col);
      
      collection[to_row][to_col]->Insert(&data[i]);
   }
}

//-------------------------------------------------------------------------
// Function to return a pointer to a specific collection of the treegrid
//-------------------------------------------------------------------------
Collection* TreeGrid::GetCollection(unsigned int r,unsigned int c)
{
   return (collection[r][c]);
}



////-------------------------------------------------------------------------
//// Function to save the treegrid to disk
////-------------------------------------------------------------------------
//void TreeGrid::SaveToDisk(std::string fname)
//{
//   //Function to save tree grid to disk
//   std::ofstream fout;
//   fout.open(fname.c_str(),std::ios::binary);
//   fout.write((char*)&rows,sizeof(unsigned long int));
//   fout.write((char*)&cols,sizeof(unsigned long int));
//   fout.write((char*)&topLeftX,sizeof(double));
//   fout.write((char*)&topLeftY,sizeof(double));
//   fout.write((char*)&sizeX,sizeof(double));
//   fout.write((char*)&sizeY,sizeof(double));
//   for(unsigned int r=0;r<rows;r++)
//   {
//      for(unsigned int c=0;c<cols;c++)
//      {
//         collection[r][c]->SaveToDisk(&fout);
//      }
//   }
//   fout.close();
//   fout.clear();
//}

////-------------------------------------------------------------------------
//// Function to load treegrid to disk
////-------------------------------------------------------------------------
//void TreeGrid::LoadFromDisk(std::string fname)
//{
//   //Function to save tree grid to disk
//   std::ifstream fin;
//   fin.open(fname.c_str(),std::ios::binary);
//   fin.read((char*)&rows,sizeof(unsigned long int));
//   fin.read((char*)&cols,sizeof(unsigned long int));
//   fin.read((char*)&topLeftX,sizeof(double));
//   fin.read((char*)&topLeftY,sizeof(double));
//   fin.read((char*)&sizeX,sizeof(double));
//   fin.read((char*)&sizeY,sizeof(double));
//   for(unsigned int r=0;r<rows;r++)
//   {
//      for(unsigned int c=0;c<cols;c++)
//      {
//         collection[r][c]->LoadFromDisk(&fin);
//      }
//   }
//   fin.close();
//   fin.clear();
//}

//-------------------------------------------------------------------------
// Function to return a pointer to a collection after creating it
//-------------------------------------------------------------------------
Collection* TreeGrid::InitialiseNewCollection(unsigned int r, unsigned int c)
{
   double cX=topLeftX+c*sizeX+0.5*sizeX;
   double cY=topLeftY-r*sizeY-0.5*sizeY;
   //Create a new collection
   Collection* coll=NULL;
   if(ellipse==NULL)
   {
      coll=new Collection(cX,cY,sizeX,sizeY);   
   }
   else
   {
      coll=new EllipsoidCollection(cX,cY,sizeX,sizeY,this->ellipse);      
   }
   return coll;
}


//-------------------------------------------------------------------------
// Default Collection constructor
//-------------------------------------------------------------------------
Collection::Collection()
{
   centreX=centreY=0;
   sizeX=sizeY=0;
}

//-------------------------------------------------------------------------
// Collection constructor when size and centre known
//-------------------------------------------------------------------------
Collection::Collection(double cX,double cY,double sX,double sY)
{
   centreX=cX;
   centreY=cY;
   sizeX=sX;
   sizeY=sY;
}

//-------------------------------------------------------------------------
// Collection destructor
//-------------------------------------------------------------------------
Collection::~Collection()
{
   items.clear();
}

//-------------------------------------------------------------------------
// Ellipsoid Collection destructor
//-------------------------------------------------------------------------
EllipsoidCollection::~EllipsoidCollection()
{
   items.clear();
}

////-------------------------------------------------------------------------
//// Function to save a collection to disk
////-------------------------------------------------------------------------
//int Collection::SaveToDisk(std::ofstream* fout)
//{
//   if(!fout->is_open())
//   {
//      std::cout<<"Error - cannot write collection, stream closed."<<std::endl;
//      return -1;
//   }
//   else
//   {
//      fout->write((char*)&centreX,sizeof(double));
//      fout->write((char*)&centreY,sizeof(double));
//      fout->write((char*)&sizeX,sizeof(double));
//      fout->write((char*)&sizeY,sizeof(double));
//      fout->write((char*)items.size(),sizeof(size_t));
//      for(std::vector<Item>::iterator it=items.begin();it<items.end();it++)
//      {
//         (*it).SaveToDisk(fout);
//      }
//   }
//   return 0;
//}

////-------------------------------------------------------------------------
//// Function to load a collection from disk
////-------------------------------------------------------------------------
//void Collection::LoadFromDisk(std::ifstream* fin)
//{
//   if(!fin->is_open())
//   {
//      std::cout<<"Error - cannot read collection, stream closed."<<std::endl;
//   }
//   else
//   {
//      fin->read((char*)&centreX,sizeof(double));
//      fin->read((char*)&centreY,sizeof(double));
//      fin->read((char*)&sizeX,sizeof(double));
//      fin->read((char*)&sizeY,sizeof(double));
//      fin->read((char*)items.size(),sizeof(size_t));
//      for(std::vector<Item>::iterator it=items.begin();it<items.end();it++)
//      {
//         (*it).LoadFromDisk(fin);
//      }
//   }
//}

//-------------------------------------------------------------------------
// Function to set the collection centre and size
//-------------------------------------------------------------------------
void Collection::SetInfo(double cX,double cY,double sX,double sY)
{
   centreX=cX;
   centreY=cY;
   sizeX=sX;
   sizeY=sY;
}

//-------------------------------------------------------------------------
// Function to insert an item into the collection
//-------------------------------------------------------------------------
inline void Collection::Insert(Item* item)
{
   //Function to insert the data into items vector
   items.push_back(*item);
}

//-------------------------------------------------------------------------
// Function to return all items without making distance calculation
//-------------------------------------------------------------------------
std::vector<Item>* Collection::GetAllItems()
{
   if(retItems.empty()!=true)
      retItems.clear();

   for(std::vector<Item>::iterator it=items.begin();it<items.end();it++)
   {
      retItems.push_back(*it);
   }
   return &retItems;
}



//-------------------------------------------------------------------------
// Function to check whether an area intersects with a collection
//-------------------------------------------------------------------------
bool Collection::Intersect(Area* area)
{
   if((area->MinX() > (centreX+0.5*sizeX))||
      (area->MaxX() < (centreX-0.5*sizeX))||
      (area->MinY() > (centreY+0.5*sizeY))||
      (area->MaxY() < (centreY-0.5*sizeY)))
      return false;
   else
      return true;
}

//------------------------------------------------------------------------
// Constructor to create and fill an IGMTreeGrid from an IGM file
//------------------------------------------------------------------------
IGMTreeGrid::IGMTreeGrid(std::string fname,std::vector<unsigned int> dropscanvector,Area* region=NULL)
{
   Logger::Verbose("... using IGM file: "+fname);
   //Create an igm worker to read in the data
   //and define the grid size
   igm=new Basic_IGM_Worker(fname);
   NullCollection=NULL;

   //Get projection from IGM file 
   std::string p=igm->Projection();
   std::string e=igm->Ellipse();
   islatlon=false;

   if(p.compare("Geographic Lat/Lon")==0)
   {
      Logger::Log("Using a geographic data set - that is in latitude/longitude.");
      //Projection is geographic lat/lon
      //Get the ellipsoid used
      if(e.compare("WGS-84")==0)
      {
         //WGS84
         ellipse=new Ellipsoid(WGS84);
         islatlon=true;
      }
      else
      {
         //Unknown as yet
         throw "Unknown ellipsoid detected in IGM hdr file.";
      }
   }
   else
   {
      Logger::Log("Assuming a projected data set - not in latitude/longitude. Uses projection keyword in igm hdr and tests vs 'Geographic Lat/Lon'.");
      ellipse=NULL;
   }

   //Calculate pixel sizes for nadir-ish pixel
   double pixsize[8]={0};
   if(!igm->GetPixelSize(igm->Samples()/2,pixsize))
   {
      //Failed to get pixel size from IGM file - maybe not an ARSF IGM file?
      Logger::Log("Pixel size calculation failed, therefore cannot use pixel sizes for TreeGrid collection size calculation.");
      //Read it from the IGM header file is next best option
      if(igm->fin->FromHeader("TreeGridSize").compare("")!=0)
      {
         pixsize[3]=pixsize[6]=StringToDouble(igm->fin->FromHeader("TreeGridSize"));
      }
      else
      {
         throw "Cannot calculate approximate pixel size from IGM. Please enter 'TreeGridSize=X' in IGM hdr file and replace X with a value to use for the "\
               "TreeGrid collection size. A value close to average pixel spacing is fine. Note that this value only affects time efficiency of the algorithms.";
      }
   }
   else
   {
      Logger::Log("Average nadir pixel size East: "+ToString(pixsize[3])+" North: "+ToString(pixsize[6]));
   }
   //Store these for use by external functions
   pixsepX=pixsize[3];
   pixsepY=pixsize[6];   
      
   //Multiply by 5 to get container size have approx 25-30 items in 
   double sX=(5*pixsize[3]);
   double sY=(5*pixsize[6]);  

   //Calculate number of rows and cols for the treegrid
   unsigned int c=0,r=0;

   //If using the IGM file (i.e. region = NULL) then we need to use
   //the min/max of the IGM with an extra buffer so that the min/max points
   //are included inside this region. Defining the buffer as a small fraction 
   //of the X/Y TreeGrid cell size
   double ttlx=igm->MinX();//-sX/20.0;
   double ttly=igm->MaxY();//+sY/20.0;
   double tbrx=igm->MaxX();//+sX/20.0;
   double tbry=igm->MinY();//-sY/20.0;

   if(region==NULL)
   {
      //Use IGM bounds (actually these are the bounds as described with the additional buffer as above) 
      c=static_cast<unsigned int>((tbrx-ttlx)/sX +1);   
      r=static_cast<unsigned int>((ttly-tbry)/sY +1);   
      Logger::Debug("Number of columns and rows of Tree Grid set up from IGM region: "+ToString(c)+" "+ToString(r));
   }
   else
   {
      //Use area region bounds
      c=static_cast<unsigned int>((region->MaxX()-region->MinX())/sX +1);   
      r=static_cast<unsigned int>((region->MaxY()-region->MinY())/sY +1); 
      Logger::Debug("Number of columns and rows of Tree Grid set up from user region: "+ToString(c)+" "+ToString(r));
   }

   Logger::Verbose("Container size: "+ToString(sX)+" "+ToString(sY));

   //Set up the grid
   if(region==NULL)
   {
      //Use the bounds of the IGM file data plus a small offset to ensure it includes the data on the bounds
      SetUpGrid(r,c,sX,sY,ttlx,ttly,tbrx,tbry);
   }
   else
   {
      //Use the bounds of the region
      SetUpGrid(r,c,sX,sY,region->MinX(),region->MaxY(),region->MaxX(),region->MinY());
   }

   //If the ellipsoid is not NULL then we should set the dateline - where the longitude values wrap
   if(ellipse!=NULL)
   {
      double radttlx=floor(ttlx)*PI/180;
      double radtbrx=ceil(tbrx)*PI/180;
      if((sin(radttlx)==sin(radtbrx))&&(cos(radttlx)==cos(radtbrx)))
      {
         upperdateline=floor(ttlx);
         lowerdateline=ceil(tbrx);
      }
      else
      {
         //Just assign to -9999
         upperdateline=-9999;
         lowerdateline=-9999;
      }
   }
   else
   {
      //Just assign to -9999
      upperdateline=-9999;
      lowerdateline=-9999;
   }

   //Insert the IGM data into the grid
   InsertData(dropscanvector,region);

   //Set the item data to the IGM file (this can be overridden at a later date)
   itemdata.Set(NULL,0,0,0,0,fname);
}

//------------------------------------------------------------------------
// IGMTreeGrid destructor
//------------------------------------------------------------------------
IGMTreeGrid::~IGMTreeGrid()
{
   if(igm!=NULL)
      delete igm;
   if(ellipse!=NULL)
      delete ellipse;
}

void IGMTreeGrid::GetAveragePixelSeparation(double &x,double &y)
{
   x=pixsepX;
   y=pixsepY;   
}

//------------------------------------------------------------------------
// Function to insert data into the IGMTreeGrid using the IGM file data
//------------------------------------------------------------------------
void IGMTreeGrid::InsertData(std::vector<unsigned int> dropscanvector,Area* region=NULL)
{
   Item* item=new Item();
   double* igmdata=NULL;
   unsigned int to_col=0,to_row=0;

   for(unsigned int myrow=0;myrow<igm->Lines();myrow++)
   {
      //Do some filtering of dropped scans here: - dropped scans affect all bands and so easiest way
      //to "get rid of them" or interpolate over them, is to not include them in the TreeGrid
      if((dropscanvector.empty()==false)&&(myrow==dropscanvector.front()))
      {
         //remove this element from the vector
         dropscanvector.erase(dropscanvector.begin());
         //skip loading in these points and go onto next line of data
         continue;
      }

      igmdata=igm->GetLine(myrow);
      for(unsigned int mycol=0;mycol<igm->Samples();mycol++)
      {
         if((igmdata[mycol]==igm->IgnoreValue())||(igmdata[mycol+igm->Samples()]==igm->IgnoreValue()))
         {
            //We don't want to insert this point into the tree as it is to be ignored
            continue;
         }

         if((region==NULL)||(region->Inside(igmdata[mycol],igmdata[mycol+igm->Samples()])))
         {
            item->igmrow=myrow;
            item->igmcol=mycol;
            item->SetData(&(this->itemdata));

            to_col=static_cast<unsigned int>(floor((igmdata[mycol]-topLeftX)/sizeX));
            to_row=static_cast<unsigned int>(floor((topLeftY-igmdata[mycol+igm->Samples()])/sizeY));

            //This needs to be here because of floating point error in above calculations
            //topleftX/Y differ slightly from min/max igmdata[] value - even though come from same value
            if(to_col==std::numeric_limits<unsigned int>::max())//((igmdata[mycol]-topLeftX) < 0)
               to_col=0;

            if(to_row==std::numeric_limits<unsigned int>::max())//((topLeftY-igmdata[mycol+igm->Samples()])<0)
               to_row=0;

            //Check that the values are within bounds - they should be but check anyway
            if((to_col >= cols)||(to_row >= rows))
            {
               Logger::Log("Error inserting IGM data into TreeGrid: out of bounds: (col,row)="+ToString(to_col)+" "+ToString(to_row));
               continue;
            }

            //Check if collection exists - it is not the null collection
            if(collection[to_row][to_col]==NullCollection)
            {
               collection[to_row][to_col]=InitialiseNewCollection(to_row,to_col);
            }

            collection[to_row][to_col]->Insert(item); 
         }        
      }
   }
   igmdata=NULL;
   delete item;
}

//-------------------------------------------------------------------------
// Function to get the map info string
//-------------------------------------------------------------------------
std::string IGMTreeGrid::GetMapInfo()
{
   std::string mapinfo;
   mapinfo=igm->Projection()+std::string(" ")+igm->Ellipse();
   return mapinfo;
}

//FIXME Needs to take into account dateline and poles if geographic lat/lon
bool TreeGrid::GetAllCollectionsWithinRadius(std::vector<Collection*>* colls,const IGMPoint* const searchpoint,double searchradius)
{
   //Identify the collection to search (containing searchX,searchY)
   unsigned long int r=0,c=0;
   r=static_cast<unsigned long int>(floor((topLeftY-searchpoint->Y)/sizeY));
   c=static_cast<unsigned long int>(floor((searchpoint->X-topLeftX)/sizeX));

   if((r>=rows)||(c>=cols))
   {
      //This collection does not exist
      Logger::Debug("Collection does not exist at row,col: "+ToString(r)+" "+ToString(c));
      return false;
   }

   //Clear the vector of collections
   colls->clear();

   //Add this collection to the list
   colls->push_back(collection[r][c]); 

   //Create an area the size of the search radius centred on searchX,searchY
   Area* search_box=NULL;  
   if(ellipse==NULL)
   {
      search_box=new Area(searchpoint->X-searchradius,searchpoint->X+searchradius,searchpoint->Y-searchradius,searchpoint->Y+searchradius);  
   }
   else
   {
      //Searchradius is in metres but here we want it in degrees to create a box to intersect with the treegrid
      double search_xdegrees=0,search_ydegrees=0;
      double destlat=0,destlon=0;

      //FIXME Need some check here to see if searchpoint + radius wraps around the dateline/pole
      //get a destination lat,lon in a pure x direction (remember azimuth is in radians) - is this azimuth correct?
      GetDestinationPoint_Bowring(searchpoint->X*PI/180,searchpoint->Y*PI/180,searchradius,90.0*PI/180,destlon,destlat,ellipse);
      search_xdegrees=fabs(destlon-searchpoint->X);
      //get a destination lat,lon in a pure y direction (remember azimuth is in radians)
      GetDestinationPoint_Bowring(searchpoint->X*PI/180,searchpoint->Y*PI/180,searchradius,0,destlon,destlat,ellipse);
      search_ydegrees=fabs(destlat-searchpoint->Y);
      search_box=new Area(searchpoint->X-search_xdegrees,searchpoint->X+search_xdegrees,searchpoint->Y-search_ydegrees,searchpoint->Y+search_ydegrees);  
   }

   //Find which collections intersect with the search box and test those 
   bool intersection=true;
   int offset=0;

   while(intersection == true)
   {
      offset++; //increase the offset to search the greater rectangle around the previous one
      intersection=false;

      //Search rows above and below
      for(long i=(long)(c-offset);i<=(long)(c+offset);i++)
      {
         if((i<0)||(i>=(long)cols))
            continue;

         //row above
         if((r-offset>=0)&&(r-offset<rows)&&(collection[r-offset][i]->Intersect(search_box)))
         {
            colls->push_back((collection[r-offset][i]));   
            intersection=true;
         }
         //row below
         if((r+offset>=0)&&(r+offset<rows)&&(collection[r+offset][i]->Intersect(search_box)))
         {
            colls->push_back((collection[r+offset][i]));  
            intersection=true;
         }
      }
      
      //Search columns to the left and right
      for(long i=(long)(r-(offset-1));i<=(long)(r+(offset-1));i++)
      {
         if((i<0)||(i>=(long)rows))
            continue;

         //col left
         if((c-offset>=0)&&(c-offset<cols)&&(collection[i][c-offset]->Intersect(search_box)))
         {
            colls->push_back((collection[i][c-offset]));   
            intersection=true;
         }
         //col right
         if((c+offset>=0)&&(c+offset<cols)&&(collection[i][c+offset]->Intersect(search_box)))
         {
            colls->push_back((collection[i][c+offset]));   
            intersection=true;
         }
      }
   }
   //TODO FIXME If geographic we also need to check against dateline and poles for wrapping and check the relevant collections

   delete search_box;
   return true;
}
