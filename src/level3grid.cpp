//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#include "level3grid.h"

//Small number for truncating doubles to ints in Level3Offset
const double EPSILON=0.0001;

//------------------------------------------------------------------------
//Constructor taking min/max X and Y and pixel size
//------------------------------------------------------------------------
Level3GridInfo::Level3GridInfo(double minX, double minY, double maxX, double maxY,double psx,double psy,std::string bandlist)
{
   double tlX=minX-psx; //add on a pixel offset
   double tlY=maxY+psy; //add on a pixel offset
   double brX=maxX+psx;
   double brY=minY-psy;

   pixelsizeX=psx;
   pixelsizeY=psy;

   //Keep the band list as a string for future use
   strBandList=bandlist;

   //We want to round the top left coords to a value such that pixelsize divides it
   //but also such that nrows and ncols will still be correct
   //remember that pix size could be a fraction so rounding to integer is not advisable
   //eg lat/lon grid with 0.00001 pix size, dont want to round 12.34567 to 12.0
   double rounded=(int)(tlY/pixelsizeY)*pixelsizeY;
   if(tlY-rounded<0.5*pixelsizeY)
      tlY=rounded;    
   else
      tlY=rounded+pixelsizeY;

   rounded=(int)(tlX/pixelsizeX)*pixelsizeX;
   if(tlX-rounded<0.5*pixelsizeX)
      tlX=rounded;    
   else
      tlX=rounded+pixelsizeX;

   bounds=new Area(tlX,brX,brY,tlY);

   nrows=(unsigned int)ceil(fabs(tlY-minY)/pixelsizeY)+1;
   ncols=(unsigned int)ceil(fabs(maxX-tlX)/pixelsizeX)+1;

   nbands=GetNumberOfItemsFromString(strBandList," ");
   bands=new unsigned int[nbands];   
   std::string b;
   unsigned int bi=0;
   while(bi<nbands)
   {
      b=GetItemFromString(strBandList,bi,' ');
      bands[bi]=StringToUINT(b);
      bi++;
   }

}

//------------------------------------------------------------------------
// Grid Info Copy Constructor 
//------------------------------------------------------------------------
Level3GridInfo::Level3GridInfo(Level3GridInfo* ref)
{
   this->nrows=ref->nrows;
   this->ncols=ref->ncols;
   this->pixelsizeX=ref->pixelsizeX;
   this->pixelsizeY=ref->pixelsizeY;
   this->bounds=new Area(ref->bounds);
   this->nbands=ref->nbands;
   this->strBandList=ref->strBandList;
   this->wavelengths=ref->wavelengths;
   this->bands=new unsigned int[this->nbands];
   for(unsigned int bi=0;bi<nbands;bi++)
   {
      bands[bi]=ref->bands[bi];
   }

}

//------------------------------------------------------------------------
// Grid Info Destructor
//------------------------------------------------------------------------
Level3GridInfo::~Level3GridInfo()
{
   if(bands!=NULL)
      delete[] bands;
   if(bounds!=NULL)
      delete bounds;
}

//------------------------------------------------------------------------
// Function to add the correct wavelengths to the string
//------------------------------------------------------------------------
void Level3GridInfo::AddWavelength(const std::string w)
{
   if(wavelengths.length()==0)
      wavelengths=wavelengths+w;
   else
      wavelengths=wavelengths+",\n"+w; 
}

//-------------------------------------------------------------------------
// Function to get the map info string
//-------------------------------------------------------------------------
std::string Level3GridInfo::GetMapInfo(std::string mapkey)
{
   std::string mapinfo;
   if((mapkey.find("utm_wgs84N")==0)||(mapkey.find("utm_wgs84S")==0))
   {
      //UTM projection detected
      mapinfo="{UTM,1,1,"+ToString(TopLeftX())+","+ToString(TopLeftY())+","+ToString(PixelSizeX())
               +","+ToString(PixelSizeY())+","+GetItemFromString(mapkey,1)+","+GetItemFromString(mapkey,2)+","+GetItemFromString(mapkey,3)+",units=Meters}";
   }
   else if(mapkey.find("Geographic Lat/Lon")!=std::string::npos)
   {
      //Lat/Lon detected
      mapinfo="{Geographic Lat/Lon, 1, 1, "+ToString(TopLeftX())+","+ToString(TopLeftY())
               +","+ToString(PixelSizeX())+","+ToString(PixelSizeY())+", "+GetItemFromString(mapkey,2)+"}";
   }
   else if(mapkey.find("osng")!=std::string::npos)
   {
      //OS national grid
      mapinfo = "{OSNG, 1, 1," +ToString(TopLeftX())+","+ToString(TopLeftY())+","+ToString(PixelSizeX())
         +","+ToString(PixelSizeY())+", Ordnance Survey of Great Britain '36, units=Meters}\n projection info = {3, 6377563.396, 6356256.910, 49.000000, -2.000000, 400000.0, -100000.0, 0.9996012717, Ordnance Survey of Great Britain '36, OSNG, units=Meters}";
   }
   else
   {
      //Unknown projection
      Logger::Log("\nUnknown map projection: You will have to fill in the projection name and datum in the map info in the .hdr file yourself.\n");
      mapinfo="{Arbitrary,1,1,"+ToString(TopLeftX())+","+ToString(TopLeftY())+","+ToString(PixelSizeX())
         +","+ToString(PixelSizeY())+",}";
   }

   return mapinfo;
}


//------------------------------------------------------------------------
// Convert the Level3Point (row,col) to an IGMpoint (X,Y)
//------------------------------------------------------------------------
int Level3GridInfo::ConvertRC2XY(L3Point* l3p, IGMPoint* igmp)
{
   //Find X - add 0.5 a pixel size to get the centre of the cell - this assumes TopLeftX is lower bound not centre of cell
   igmp->X=(l3p->col*PixelSizeX()) + TopLeftX() + 0.5*PixelSizeX();
   //Find Y - subtract 0.5 a pixel size to get the centre of the cell - this assumes TopLeftY is upper bound not centre of cell
   igmp->Y=TopLeftY() - (l3p->row*PixelSizeY()) - 0.5*PixelSizeY();
   //Test if in bounds
   if((igmp->X < TopLeftX())||(igmp->X > BottomRightX())
      ||(igmp->Y > TopLeftY())||(igmp->Y < BottomRightY()))
   {
      //xy point out of bounds
      return -1;
   }      

   //xy point in bounds
   return 0; 
}



//------------------------------------------------------------------------
// Construct an empty Level3Grid - this shouldn't be used explicitly
//------------------------------------------------------------------------
Level3Grid::Level3Grid()
{
   info=NULL;
}

//------------------------------------------------------------------------
// Construct a level 3 grid from a level3gridinfo object
//------------------------------------------------------------------------
Level3Grid::Level3Grid(Level3GridInfo* gi)
{
   info=new Level3GridInfo(gi);
}


//------------------------------------------------------------------------
// Construct a level 3 grid from a defined output area
//------------------------------------------------------------------------
Level3Grid::Level3Grid(double output_pixel_sizeX,double output_pixel_sizeY,std::string bandlist,Area* outputrect)
{
   if(outputrect != NULL)
   {
      //Create the grid info object
      info=new Level3GridInfo(outputrect->MinX(),outputrect->MinY(),outputrect->MaxX(),outputrect->MaxY(),output_pixel_sizeX,output_pixel_sizeY,bandlist);
   }
   else
   {
      throw "NULL pointer - cannot create level 3 grid from a NULL output rectangle.";      
   }

}

//------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------
Level3Grid::~Level3Grid()
{
   if(info!=NULL)
      delete info;
}

//------------------------------------------------------------------------
// Constructor for Level 3 Outline from igm filename
// This assumes an ARSF style IGM file with each scan line being a full
// across track swath...outline will fail if not.
// i.e. it assumes that the edge pixels of the IGM file describe the edge
// of the swath
//------------------------------------------------------------------------
Level3Outline::Level3Outline(Level3GridInfo* const gi,std::string igmfilename)
{
   //Set up a bool for keeping track of calls to GetIntersectEdgeRow
   firstcalltogetintersects=true;
   //Copy the level 3 info
   info=new Level3GridInfo(gi);
   //IGM reader
   Basic_IGM_Worker igm(igmfilename); 

   if(!igm.IsARSFStyle())
   {
      Logger::Log("Outline cannot be created from an IGM file not in ARSF style - faking outline to be the full size of the grid.");
      FakeEdges();
   }
   else
   {
      std::vector<IGMPoint> edges;
      //Read in the edge points (X,Y) from the IGM file and store in the edges vector
      if(!ReadEdges(&igm,&edges))
      {
         //Something has gone wrong - just use every pixel by faking the edges
         FakeEdges();
      }
      else
      {
         //Iterate through this edge vector and convert points from X,Y into Row,Col of
         //the level 3 grid. If required, fill in any holes in the outline edges. Also
         //add the lines in that reflect the edge at the start and end of IGM (i.e. first
         //and last rows of the IGM)
         //Initialise(edges);
         InitialiseForScanlineFill(edges);
      }
   }
}

//------------------------------------------------------------------------
// Constructor for Level 3 Outline from igm data array
//------------------------------------------------------------------------
Level3Outline::Level3Outline(Level3GridInfo* const gi,Block<double>* igmblock,unsigned int nlines, unsigned int start_line_offset,double ignoreval)
{
   //Set up a bool for keeping track of calls to GetIntersectEdgeRow
   firstcalltogetintersects=true;
   //Copy the level 3 info
   info=new Level3GridInfo(gi);

   //FIXME Need a better test here for non-ARSF style flightline
   if(igmblock->Lines()==1)
   {
      Logger::Log("Outline cannot be created from an IGM file not in ARSF style - faking outline to be the full size of the grid.");
      FakeEdges();
   }
   else
   {
      //Make a dataaccessor object - then we can use same functions 
      //for both outlines created by (a) igm files and (b) data blocks
      std::vector<IGMPoint> edges;
      unsigned int bandlist[2]={0,1}; //we only need read bands 0 and 1 of the igm file
      //Create a data accessor object based on the block of data - igm filename is null
      DataAccessor<double> data(igmblock,"",bandlist);
      //Get the IGM index for the first line of data in the block discarding any overlap region if present 
      //(e.g. block line (0+overlap) is line index X in IGM)
      unsigned int firstlineofdata=igmblock->FirstRow()+start_line_offset;
      //Increment counter for loop
      unsigned int lineinc=0;
      //For each line of the block - add the swath edge points to the 'edge' array
      while(lineinc<nlines)
      {
         if(!AddPairToEdgeArray(&data,ignoreval,edges,firstlineofdata+lineinc,igmblock->Samples()))
         {
            break;
         }
         lineinc++;
      }
      if(lineinc!=nlines)
      {
         //Something has gone wrong - just use every pixel by faking the edges
         FakeEdges();
      }
      else
      {

         //Now add the top/bottom edges
         if(!AddTopBottomToEdgeArray(&data,ignoreval,edges,firstlineofdata,nlines,igmblock->Samples()))
         {
            //Something has gone wrong - just use every pixel by faking the edges
            FakeEdges();
         }
         else
         {
            //Iterate through this edge vector and convert points from X,Y into Row,Col of
            //the level 3 grid. If required, fill in any holes in the outline edges. Also
            //add the lines in that reflect the edge at the start and end of IGM (i.e. first
            //and last rows of the IGM)
            //Initialise(edges);
            InitialiseForScanlineFill(edges);
         }
      }
   }
}

//------------------------------------------------------------------------
// Convert the IGM point (X,Y) to a level3 point (row,col)
//------------------------------------------------------------------------
inline int Level3Outline::ConvertXYtoRC(IGMPoint* igmp,L3Point* l3p)
{
   l3p->row=static_cast<long>(floor((info->TopLeftY()-igmp->Y)/info->PixelSizeY()));
   l3p->col=static_cast<long>(floor((igmp->X-info->TopLeftX())/info->PixelSizeX()));

   //Check if in bounds or not
   if((l3p->row > info->NumRows()) || (l3p->col > info->NumCols())||
      (l3p->row < 0) || (l3p->col <0))
   {
      return -1;
   }
   return 0;
}

//------------------------------------------------------------------------
// Add a point to the edge array
//------------------------------------------------------------------------
bool Level3Outline::AddPointToEdgeArray(DataAccessor<double>* data, double ignoreval,std::vector<IGMPoint> &edges,unsigned int line,unsigned int sample)
{
   IGMPoint p(0,0);

   //Read the near edge values from the IGM and assign to point
   p.X=(data->GetData(0,line,sample));
   p.Y=(data->GetData(1,line,sample));

   //Add the point if not to be ignored
   if((p.X != ignoreval)&&(p.Y != ignoreval))
   {
      edges.push_back(p); 
   }
   else
   {
      return false;
   }
   return true;
}

//------------------------------------------------------------------------
// Add a pair of points (nearside and farside) to the edge array
//------------------------------------------------------------------------
bool Level3Outline::AddPairToEdgeArray(DataAccessor<double>* data, double ignoreval,std::vector<IGMPoint> &edges,unsigned int line,unsigned int nsamples)
{
   //Add the near-side point
   if(!AddPointToEdgeArray(data, ignoreval,edges,line,0))
      return false;
   //Add the far-side point
   if(!AddPointToEdgeArray(data, ignoreval,edges,line,nsamples-1))
      return false;
   return true;
}

//------------------------------------------------------------------------
// Add the top/bottom (i.e. across swath direction) edges
//------------------------------------------------------------------------
bool Level3Outline::AddTopBottomToEdgeArray(DataAccessor<double>* data, double ignoreval,std::vector<IGMPoint> &edges,unsigned int firstlineofdata,unsigned int nlines,unsigned int nsamples)
{
   for(unsigned int sample=0;sample<nsamples;sample++)
   {
      //Add the top-side point
      if(!AddPointToEdgeArray(data, ignoreval,edges,firstlineofdata,sample))
         return false;
      //Add the bottom-side point
      if(!AddPointToEdgeArray(data, ignoreval,edges,firstlineofdata+nlines-1,sample))
         return false;      
   }
   return true;
}

//------------------------------------------------------------------------
// Read in the edges from the IGM file and store in the edge array
//------------------------------------------------------------------------
bool Level3Outline::ReadEdges(Basic_IGM_Worker* igm,std::vector<IGMPoint>* edges)
{
   unsigned int bandlist[2]={0,1}; //we only need read bands 0 and 1 of the igm file
   DataAccessor<double> data(NULL,igm->FileName(),bandlist);
   //Add the port/starboard edges
   for(unsigned int line=0;line<igm->Lines();line++)
   {
      if(!AddPairToEdgeArray(&data,igm->IgnoreValue(),*edges,line,igm->Samples()))
         return false;
   }
   //Now add the top/bottom edges
   if(!AddTopBottomToEdgeArray(&data,igm->IgnoreValue(),*edges,0,igm->Lines(),igm->Samples()))
      return false;

   return true;
}

//------------------------------------------------------------------------
// Insert edge into vector
//------------------------------------------------------------------------
inline void Level3Outline::InsertEdge(Edge* p1)
{
   edgetable.push_back((*p1));
}


//------------------------------------------------------------------------
// Create and add the edge to the edge table 
//------------------------------------------------------------------------
inline void Level3Outline::CreateEdge(L3Point* p1,L3Point* p2)
{
   //Need to create the Edge object and populate
   //An Edge always has: the X of the point with minimum Y value
   //                    the Y of both points 
   //                    the gradient (in terms of X)
   Edge EDGE;

   if(p1->row > p2->row)
   {
      EDGE.minX=p2->col;
      EDGE.maxY=p1->row;
      EDGE.minY=p2->row;
      //gradient from lowest (min y) to highest vertex
      //note we want gradient for x not y (e.g. x=(1/m)*y +c)
      EDGE.grad=(double)(p1->col-p2->col)/(p1->row-p2->row);
   }
   else if(p1->row < p2->row)
   {
      EDGE.minX=p1->col;
      EDGE.maxY=p2->row;
      EDGE.minY=p1->row;
      //gradient from lowest (min y) to highest vertex
      //note we want gradient for x not y (e.g. x=(1/m)*y +c)
      EDGE.grad=(double)(p2->col-p1->col)/(p2->row-p1->row);
   }
   else
   {
      //We dont want horizontal edges in the list
      return;
   }
   //Add the edge to the table
   InsertEdge(&EDGE);
}

//------------------------------------------------------------------------
// Function to fake the egde array, creating a full-size edge container
// - i.e. such that the whole level3 grid is contained
//------------------------------------------------------------------------
void Level3Outline::FakeEdges()
{
   Logger::Debug("Faking edges...");
   for(unsigned int row=0;row<info->NumRows()-1;row++)
   {
      L3Point p1(row,0);
      L3Point p2(row+1,0);
      CreateEdge(&p1,&p2);
      L3Point p3(row, info->NumCols()-1);
      L3Point p4(row+1, info->NumCols()-1);
      CreateEdge(&p3,&p4);
   }
}

//------------------------------------------------------------------------
// Initialisation code similar to both constructors
//   For generating a table of edges to perform a scanline fill style 
//   algorithm
//------------------------------------------------------------------------
void Level3Outline::InitialiseForScanlineFill(std::vector<IGMPoint> edges)
{
   std::vector<IGMPoint>::iterator it;
   std::vector<IGMPoint>::iterator startpoint;
   L3Point p1(0,0);
   L3Point p2(0,0);
   bool badpoint=false;

   //Check the number of edges are even - this should always be the case
   if(edges.size()%2!=0)
   {
      Logger::Log("There are "+ToString(edges.size())+" outline points");
      throw "Error in Level3Outline::InitialiseForScanlineFill. Uneven number of edge points - this should be impossible!";
   }
   else
   {
      Logger::Verbose("There are "+ToString(edges.size())+" outline points");
   }


   //2 iterations: 1 each for port and starboard sides
   for(int i=0;i<2;i++)
   {
      if(i==0)
      {  
         //Initialise for one side of swath
         startpoint=edges.begin();
      }
      else if(i==1)
      {
         //Initialise for other side of swath
         startpoint=edges.begin()+1;
      }
      else
      {
         throw "InitialiseForScanlineFill(): More than 2 iterations of loop - This error should not happen!";
      }

      //Iterate through stepping every other point 
      for(it=startpoint;it<edges.end()-2;it=it+2)
      {
         //Is first point in area
         if(ConvertXYtoRC(&(*it),&p1) == -1)
         {
            badpoint=true;
            Logger::Debug("Bad point detected: "+ToString(p1.row)+" "+ToString(p1.col)+" "+ToString((*(it)).X)+" "+ToString((*(it)).Y)+" "+ToString(std::distance(edges.begin(),it)));
            //Maybe skip until we get a good point
         }

         //Is second point in area
         if(ConvertXYtoRC(&(*(it+2)),&p2) == -1)
         {
            badpoint=true;
            Logger::Debug("Bad point detected: "+ToString(p2.row)+" "+ToString(p2.col)+" "+ToString((*(it+2)).X)+" "+ToString((*(it+2)).Y)+" "+ToString(std::distance(edges.begin(),it+2)));
            //Maybe skip until we get a good point
         }

         //Calculate gradient of line connecting vertices and store maxY (maximum y vertex value), X of vertex with min Y and gradient (in X not Y)
         if(badpoint==false)
         {
            CreateEdge(&p1,&p2);
         }
         //reset badpoint bool
         badpoint=false;
      }
   }

   //Now finalise by adding the start/end edges 
   it=edges.begin();
   if(ConvertXYtoRC(&(*it),&p1) == -1)
   {
      badpoint=true;
      Logger::Debug("Bad point detected: "+ToString(p1.row)+" "+ToString(p1.col)+" "+ToString((*(it)).X)+" "+ToString((*(it)).Y)+" "+ToString(std::distance(edges.begin(),it)));
   }

   if(ConvertXYtoRC(&(*(it+1)),&p2) == -1)
   {
      badpoint=true;
      Logger::Debug("Bad point detected: "+ToString(p2.row)+" "+ToString(p2.col)+" "+ToString((*(it+1)).X)+" "+ToString((*(it+1)).Y)+" "+ToString(std::distance(edges.begin(),it+1)));
   }

   if(badpoint == false)
   {
      CreateEdge(&p1,&p2);
   }

   //Reset badpoint and go to end edge
   badpoint=false;
   it=edges.end();
   if(ConvertXYtoRC(&(*(it-2)),&p1) == -1)
   {
      badpoint=true;
      Logger::Debug("Bad point detected: "+ToString(p1.row)+" "+ToString(p1.col)+" "+ToString((*(it-2)).X)+" "+ToString((*(it-2)).Y)+" "+ToString(std::distance(edges.begin(),it-2)));
   }
   if(ConvertXYtoRC(&(*(it-1)),&p2) == -1)
   {
      badpoint=true;
      Logger::Debug("Bad point detected: "+ToString(p2.row)+" "+ToString(p2.col)+" "+ToString((*(it-1)).X)+" "+ToString((*(it-1)).Y)+" "+ToString(std::distance(edges.begin(),it-1)));
   }

   if(badpoint == false)
   {
      CreateEdge(&p1,&p2);
   }

   //Now destroy the grid info
   delete info;

   //We now have a vector containing row/col points of the 
   //edges of the flight line.
   //Sort it by row then by column
   std::sort(edgetable.begin(),edgetable.end(),Edge::compare);
   Logger::Verbose("There are "+ToString(edgetable.size())+" edges.");
}

//------------------------------------------------------------------------
// Function to find the intersects between rows and outline edge to use
// as the processing bounds in the map
// NOTE: row is NOT necessarily relative to the final output level3 grid
//       but only to the level3 grid used to construct level3outline
//       i.e. the segment grid not the output grid
//------------------------------------------------------------------------
void Level3Outline::GetEdgeIntersectsOfRow(int row,std::vector<int> &intersects)
{
   std::vector<Edge>::iterator it;
   //Check if given row is less than max row in outline - else return
   //This assumes outline has been ordered by min -> max row
   it=edgetable.end();
   if((*(it-1)).maxY < row)
   {
      Logger::Debug("Call to get edge intersects of row past last row in edgetable: "+ToString(row));
      return;
   }
   //Check if row is before first row in edgetable
   it=edgetable.begin();
   if((*it).minY > row)
   {
      Logger::Debug("Call to get edge intersects of row before first row in edgetable: "+ToString(row));
      return;
   }   

   //Check if first run and if so we may need to add previous rows that should be activated already
   //but arent because we may be starting from row != 0. This is the case for when using an area.
   if((firstcalltogetintersects==true)&&(row!=0))
   {
      Logger::Debug("Setting up active edge table."+ToString(row)+" "+ToString((*it).minY));
      while(((*it).minY<row)&&(it<edgetable.end()))
      {
         //If the maxY is > than row then it is still active and should be added
         //As this is the first loop we know it does not already exist in the array
         if((*it).maxY > row)
            activeedges.push_back(*(it));

         it++;
      }
   }

   //Update by removing old 'deactivated' edges - these are edges which no longer affect the output
   //for the row in question (or any other rows IF they are done in an ordered fashion)
   for(it=activeedges.begin();it<activeedges.end();it++)
   {
      if(row >= (*it).maxY)
      {
         activeedges.erase(it);
         it--;//send loop variable back one to take into account the fact one has been deleted from the array
      }
   }

   //Scan until we get to first record on given row
   it=edgetable.begin();
   while(((*it).minY<row)&&(it<edgetable.end()))
      it++;

   //Update by adding new 'activated' edges - these are edges which have come into effect on this row
   while(((*it).minY==row)&&(it<edgetable.end()))
   {
      activeedges.push_back(*(it));
      it++;
   }

   //Get row intersects with edges and return
   for(it=activeedges.begin();it<activeedges.end();it++)
   {
      //If the vertex is on this row then use minX for the col intersect
      //as this is the X value for the vertex
      if((*it).minY==row)
      {
         (intersects).push_back(((*it).minX));
      }
      //Else we need to calculate the col intersect value - do this from the increase in row value * gradient + X1
      else
      {
         (intersects).push_back(ceil((*it).minX+(row-(*it).minY)*(*it).grad));
      }
   }   

   //Sort from lower col to higher col
   std::sort((intersects).begin(),(intersects).end());
   //Set this to false - only needs to be done once but this does it every call
   firstcalltogetintersects=false;
}

//------------------------------------------------------------------------
// Write out the outline to an ascii file
//------------------------------------------------------------------------
void Level3Outline::WriteEdge(std::string filename)
{
   std::ofstream fout;
   fout.open(filename.c_str());
   std::vector<Edge>::iterator it;
   //Give the data a 3rd column of value 255 - used for importing into grass 
   //as a raster with r.in.xyz
   int x=0;
   for(it=edgetable.begin();it<edgetable.end();it++)
   {
      for(int y=(*it).minY;y<(*it).maxY;y++)
      {
         //For each edge write out all the outline points
         x=ceil((*it).minX+(y-(*it).minY)*(*it).grad);
         fout<<x<<" "<<y<<" "<<"255"<<std::endl;
      }
   }
   fout.close();
   fout.clear();
}

