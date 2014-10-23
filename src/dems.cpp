//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#include "dems.h"

#ifndef DEMDEBUG
   #define DEBUGPRINT(X)  
#else
   #define DEBUGPRINT(X) std::cout.precision(20);std::cout<< X <<std::endl;
#endif

#ifndef PI_
#define PI_
const double PI=4*atan(1.0);
#endif


/******************************
   DEM_AOI Methods
*******************************/

//-------------------------------------------------------------------------
//Constructor with no inputs
//-------------------------------------------------------------------------
DEM_AOI::DEM_AOI()
{
   //Default to zeros
   this->Set(0,0,0,0);
}

//-------------------------------------------------------------------------
//Constructor with AOI defined
//-------------------------------------------------------------------------
DEM_AOI::DEM_AOI(const double illx, const double illy,const double iurx,const double iury)
{
   //Set the bounds to those given
   this->Set(illx,illy,iurx,iury);
}

//-------------------------------------------------------------------------
//Retrieve the component of the vertex as specified
//-------------------------------------------------------------------------
double DEM_AOI::Get(const vertex v)
{
   switch(v)
   {
   case LLX: //lower left x
      return this->llx;
   case LLY: //lower left y
      return this->lly;
   case URX: //uper right x
      return this->urx;
   case URY: //upper right y
      return this->ury;
   default:
      return 0;
   }
}

//-------------------------------------------------------------------------
//Set the AOI to the lower left / upper right extents as given
//-------------------------------------------------------------------------
bool DEM_AOI::Set(const double illx, const double illy,const double iurx,const double iury)
{
   //If lower left x and y are less than (or equal to) upper 
   //right x and y, set the bounds and return true
   if((illx <= iurx) && (illy <= iury))
   {
      llx=illx; 
      lly=illy; 
      urx=iurx; 
      ury=iury;
      return true;
   }
   else //else set to zeros and return false
   {
      llx=0; 
      lly=0; 
      urx=0; 
      ury=0;
      return false;
   }
}

/***********************************
   DEM Class Methods 
************************************/

//-------------------------------------------------------------------------
//Constructor for dem object. pass it the filename of the DEM file
//-------------------------------------------------------------------------
DEM::DEM(std::string strFilename)
{
   DEBUGPRINT("Entering DEM constructor...")

   //Set pointer to NULL
   this->data=NULL;

   //create the bil object
   this->file=new DEMBinFile(strFilename);

   //check the dem file is currently supported
   this->CheckSupport();

   //Get the map info into the associated variables
   //0: projection, 
   //1 and 2: ref pixel c,r 
   //3 and 4: ref pixel x,y
   //5 and 6: x,y spacing
   //7: datum

   //Get the reference pixel column value
   this->refc=StringToDouble(this->file->FromHeader("map info",1));  
   if((errno==ERANGE) || (refc <= 0))
   {
      throw "Error trying to convert reference point column to double. Value should be greater than 0 (top left pixel is 1,1)";
   }
   else
   {
      DEBUGPRINT("refc: "<<refc)
   }

   //Get the reference pixel row value
   this->refr=StringToDouble(this->file->FromHeader("map info",2));  
   if((errno==ERANGE)|| (refr <= 0))
   {
      throw "Error trying to convert reference point row to double. Value should be greater than 0 (top left pixel is 1,1)";
   }
   else
   {
      DEBUGPRINT("refr: "<<refr)
   }

   //Get the reference pixel x value
   this->refx=StringToDouble(this->file->FromHeader("map info",3));  
   if(errno==ERANGE)
   {
      throw "Error trying to convert reference point x to double.";
   }
   else
   {
      DEBUGPRINT("refx: "<<refx)
   }

   //Get the reference pixel y value
   this->refy=StringToDouble(this->file->FromHeader("map info",4));  
   if(errno==ERANGE)
   {
      throw "Error trying to convert reference point y to double.";
   }
   else
   {
      DEBUGPRINT("refy: "<<refy)
   }

   //Get the x spacing value
   this->xspace=StringToDouble(this->file->FromHeader("map info",5));  
   if(errno==ERANGE)
   {
      throw "Error trying to convert x spacing to double.";
   }
   else
   {
      DEBUGPRINT("xspace: "<<xspace)
   }

   //Get the y spacing value
   this->yspace=StringToDouble(this->file->FromHeader("map info",6));  
   if(errno==ERANGE)
   {
      throw "Error trying to convert y spacing to double.";
   }
   else
   {
      DEBUGPRINT("yspace: "<<yspace)
   }

   //Now calculate the bounds of the DEM file
   //get file size rows/cols
   nrows=StringToUINT(this->file->FromHeader("lines"));
   ncols=StringToUINT(this->file->FromHeader("samples"));

   //We use refc-1 becuase the value of top left pixel is 1,1 (not 0,0)
   this->minx=refx - (refc - 1) * xspace;
   this->maxx=refx + (ncols - (refc - 1)) * xspace;

   this->miny=refy - (nrows - (refr - 1)) * yspace;
   this->maxy=refy + (refr - 1) * yspace;   

   DEBUGPRINT("DEM file bounds are; minx:="<<minx<<" maxx:="<<maxx<<" miny:="<<miny<<" maxy:="<<maxy)

}

//-------------------------------------------------------------------------
// DEM destructor
//-------------------------------------------------------------------------
DEM::~DEM()
{
   DEBUGPRINT("Entering DEM destructor...")

   //Close and destroy the bilreader
   this->file->Close();

   //Free up data array
   if(this->data != NULL)
      delete[] this->data;
}

//-------------------------------------------------------------------------
//Check the DEM file to see if it is supported
//-------------------------------------------------------------------------
void DEM::CheckSupport()
{
   //Check there is only 1 band in the BIL file
   if(this->file->FromHeader("bands") != "1")
      throw "DEM Unsupported: Expected DEM to be a 1-band BIL/BSQ file. Number of bands reported in hdr is not 1";

   //Check if "map info" is present in bil header
   std::string mapinfo=this->file->FromHeader("map info");
   DEBUGPRINT("map info string from hdr: "+mapinfo)
   if(mapinfo.compare("") == 0)
   {
      //Empty string returned - no map info in hdr file
      throw "DEM Unsupported: No map information in hdr file for DEM. Expected 'map info' entry.";
   }
   else
   {
      //Map info found - should probably do some further tests
      //Check that data is Geographic latitude and longitude
      if(mapinfo.find("Geographic Lat/Lon") == std::string::npos)
         throw "DEM Unsupported: Projection is not Geographic Lat/Lon. Expected first part of map info in hdr file to be 'Geographic Lat/Lon'.";

      //Check that the DEM is in WGS84 datum
      if(mapinfo.find("WGS-84") == std::string::npos)
         throw "DEM Unsupported: Projection is not in WGS84. Expected map info in hdr file to contain 'WGS-84'.";
   }
}

//-------------------------------------------------------------------------
//Convert X coordinate to file column
//-------------------------------------------------------------------------
double DEM::X2C(const double x)
{
   //Get the number of columns in the file as a sanity check
   //unsigned int ncols=StringToUINT(this->file->FromHeader("samples"));
   //Calculate the column
   double c=((x-this->minx)/this->xspace);
   //test and return
   if ((c >= 0)&&(c <= ncols+0.5)) //need to fix this floating point error problem
      return c;
   else
      return -1;
}

//-------------------------------------------------------------------------
//Convert the Y coordinate to file row
//-------------------------------------------------------------------------
double DEM::Y2R(const double y)
{
   //Get the number of rows as a sanity check
   //unsigned int nrows=StringToUINT(this->file->FromHeader("lines"));
   //Calculate the row
   double r=((this->maxy-y)/this->yspace);
   //test and return
   if ((r >= 0)&&(r <= nrows+0.5)) //need to fix this floating point error problem
      return r;
   else
      return -1;
}

//-------------------------------------------------------------------------
//Convert the file column to X coordinate
//-------------------------------------------------------------------------
double DEM::C2X(const double c)
{
   double x=this->minx + c * this->xspace;
   //Not testing against minx/maxx
   return x;   
}

//-------------------------------------------------------------------------
//Convert the file row to Y coordinate
//-------------------------------------------------------------------------
double DEM::R2Y(const double r)
{
   double y=this->maxy - r * this->yspace;
   //Not testing against miny/maxy
   return y;   
}

//-------------------------------------------------------------------------
//Read in a rectangles worth of data from the DEM file, described by AOI.
//-------------------------------------------------------------------------
int DEM::ReadRect(char* const chdata)
{
   return ReadRect(chdata, this->AOI.Get(LLX), this->AOI.Get(LLY),this->AOI.Get(URX),this->AOI.Get(URY));
}

//-------------------------------------------------------------------------
//Read in a rectangles worth of data from the DEM file.
//-------------------------------------------------------------------------
int DEM::ReadRect(char* const chdata, const double llx, const double lly,const double urx,const double ury)
{
   //Convert the x,y coordinates to file row,col
   int minrow=rounded(Y2R(ury));
   int maxrow=rounded(Y2R(lly));
   int mincol=rounded(X2C(llx));
   int maxcol=rounded(X2C(urx));

   DEBUGPRINT("Area defined by lat/lon is: lats:"<<ury<<" "<<lly<<" lons:"<<urx<<" "<<llx)
   DEBUGPRINT("Area defined by row/col is: rows:"<<minrow<<" "<<maxrow<<" cols:"<<mincol<<" "<<maxcol)

   //Check the mins are less than the maxs
   if ((minrow > maxrow) || (mincol > maxcol))
      throw "Order of elements in DEM.ReadRect should be llx,lly urx,ury. Min row/col is greater than max row/col.";
   
   //Check all the columns and rows are positive
   if ((minrow < 0)||(maxrow < 0)||(mincol < 0)||(maxcol < 0))
      throw "Negative row/col of file is not allowed. Check X,Y corners of rectangle in DEM.ReadRect.";   

   //Check that the passed upper right coordinates are below the maxima of the DEM coverage
   if((ury > this->maxy)||(urx > this->maxx))
      throw "Upper right coordinate elements passed to DEM.ReadRect are greater than DEM bounds.";

   //Read in the rectangle of data and store in chdata array
   file->ReadRect(chdata,minrow,maxrow,mincol,maxcol);

   return 0;
}

//-------------------------------------------------------------------------
//Return the size of the AOI in bytes (in case want to pass data arrays to read data)
//-------------------------------------------------------------------------
unsigned long int DEM::SizeOf()
{
   return SizeOf(this->AOI.Get(LLX), this->AOI.Get(LLY),this->AOI.Get(URX),this->AOI.Get(URY));
}

//-------------------------------------------------------------------------
//Calculates the size of a rectangle, in bytes, to determine the size of an array for ReadRect
//-------------------------------------------------------------------------
unsigned long int DEM::SizeOf(const double llx, const double lly,const double urx,const double ury)
{
   //Convert the x,y coordinates to file row,col
   int minrow=rounded(Y2R(ury));
   int maxrow=rounded(Y2R(lly));
   int mincol=rounded(X2C(llx));
   int maxcol=rounded(X2C(urx));

   //Check the mins are less than the maxes
   if ((minrow > maxrow) || (mincol > maxcol))
      throw "Order of elements in DEM.SizeOf should be llx,lly urx,ury. Min row/col is greater than max row/col.";

   if ((minrow < 0)||(maxrow < 0)||(mincol < 0)||(maxcol < 0))
      throw "Negative row/col of file is not allowed. Check X,Y corners of rectangle in DEM.SizeOf.";  

   //Check that the passed upper right coordinates are below the maxima of the DEM coverage
   if((ury > this->maxy)||(urx > this->maxx))
      throw "Upper right coordinate elements passed to DEM.SizeOf are greater than DEM bounds.";

   //Calculate the required size in bytes - cast to int here else nastiness may occur 
   //should probably round up for maxes and round down for mins?
   unsigned long int bytesize=(maxrow-minrow +1)*(maxcol-mincol+1)*file->GetDataSize();

   DEBUGPRINT("min/max row/col: "<<minrow<<" "<<maxrow<<" "<<mincol<<" "<<maxcol)
   DEBUGPRINT("Size of area in bytes is:"<<bytesize)

   //Return the size in bytes
   return bytesize;   
}

//-------------------------------------------------------------------------
// Set the Area of Interest to the given rectangle
//-------------------------------------------------------------------------
bool DEM::SetAOI(const double llx, const double lly,const double urx,const double ury)
{
   //Check the area is within the DEM bounds
   if((llx >= this->minx)&&(lly >= this->miny)&&(urx<=this->maxx)&&(ury<=this->maxy))
   {
      //Set the AOI
      this->AOI.Set(llx,lly,urx,ury);
      return FitAOIToGrid();
   }
   else //not fuly contained in DEM
   {
      return false; //return false
   }
}

//-------------------------------------------------------------------------
// Extend the AOI such that it fits to the DEM grid cell boundaries so
// that we can work with integer amounts of cells
//-------------------------------------------------------------------------
bool DEM::FitAOIToGrid()
{
   //Need to find the nearest grid line n/e/s/w of the AOI
   //Truncate to get the closest grid line west of the AOI
   double nc=floor((AOI.Get(LLX) - minx) / xspace);
   double newllx=minx+nc*xspace;
   //Get the closest grid line east of the AOI
   nc=ceil((AOI.Get(URX) - minx) / xspace);
   double newurx=minx+nc*xspace;
   //Get the closest grid line north of the AOI
   nc=floor((maxy - AOI.Get(URY))/yspace);
   double newury=maxy - nc*yspace;
   //Get the closest grid line south of the AOI
   nc=ceil((maxy - AOI.Get(LLY))/yspace);
   double newlly=maxy - nc*yspace;
   
   //Extend the AOI
   return AOI.Set(newllx,newlly,newurx,newury);
}


//-------------------------------------------------------------------------
//Return the value of the given AOI vertex
//-------------------------------------------------------------------------
double DEM::GetAOI(const vertex v)
{
   return this->AOI.Get(v);
}

//-------------------------------------------------------------------------
// Fill the chdata array with DEM data
//-------------------------------------------------------------------------
void DEM::FillArray()
{
   //Create the char* data array and fill it
   //with data from the DEM file relating to the AOI
   if(data != NULL)
   {
      //free up the array
      delete[] data;
      data=NULL;
   }
   
   //Create array the size of the AOI
   data=new char[this->SizeOf()];
   
   //Read in the data
   this->ReadRect(this->data);
}

//-------------------------------------------------------------------------
//Return the cell index that relates to the given longitude/latitude
//-------------------------------------------------------------------------
unsigned int DEM::GetAOICell(const double lon,const double lat)
{
   
   //Check if given point is within the AOI. If not, return -1
   double toonorth=this->AOI.Get(URY)-lat;
   double toosouth=lat-this->AOI.Get(LLY);
   double toowest=lon-this->AOI.Get(LLX);
   double tooeast=this->AOI.Get(URX)-lon;

   if((toowest<0)||(tooeast<0)||(toonorth<0)||(toosouth<0))
   {
      //throw "Error! Requested position is outside of imported DEM data array.";
      return UINT_MAX; //This is UINT_MAX since previously returned -1 and is an unsigned int
   }

   //Point is within bounds so get lat/lon pos wrt boundary top left
   double latdiff=this->AOI.Get(URY)-lat;
   double londiff=lon-this->AOI.Get(LLX);

   //Cell that the data point falls in (if array was a 2d array)
   double ycell=(latdiff / this->yspace);
   double xcell=(londiff / this->xspace);

   //Number of cells in a line of the AOI
   //Round each one to nearest integer and plus 1
   unsigned int ncells=(floor(X2C(this->AOI.Get(URX))+0.5) - floor(X2C(this->AOI.Get(LLX))+0.5)+1);
   DEBUGPRINT("Number of cells in AOI per line is: "<<ncells)

   //Cell number we want data from (since array is actually a 1d array)
   //Adding 0.01 on below before truncating due to double rounding errors
   unsigned int cell=floor(ycell+0.01)*ncells+floor(xcell+0.01);
   DEBUGPRINT("Lat,Lon is: "<<lat<<" "<<lon<<" xcell,yxell is: "<<xcell<<" "<<ycell<<" AOI Cell is: "<<cell)
   return cell;
}

//-------------------------------------------------------------------------
//Function to retrieve height at a given lon,lat, from the data array
//-------------------------------------------------------------------------
double DEM::GetHeight(const double lon,const double lat)
{
   //Get the AOI cell that relates to the given latitude/longitude
   unsigned int cell=GetAOICell(lon,lat);

   //Check that the lat/lon is within the AOI
   if(cell==UINT_MAX) //largest value since cell is unsigned
   {
      return DEMOutOfBounds; //error flag
   }
   char* cp=NULL;
   short int* sip=NULL;
   float* fp=NULL;
   unsigned short int* usip=NULL;
   int* ip=NULL;
   unsigned int* uip=NULL;
   double* dp=NULL;

   double retval=0;

   if(data==NULL)
   {
      throw "Attempt to read from DEM.data when it is still NULL.";      
   }

   switch(file->GetDataType())
   {
   case 1: //8-bit
      cp=(char*)(data);
      retval=(double)cp[cell];
      break;
   case 2: //16 bit signed int
      sip=(short int*)(data);      
      retval=(double)sip[cell];
      break;
   case 3:
      ip=(int*)(data);
      retval= (double)ip[cell];
      break;
   case 4: //float
      fp=(float*)(data);
      retval= (double)fp[cell];
      break;
   case 5: //double
      dp=(double*)(data);
      retval= dp[cell];
      break;
   case 12: //16 bit unsigned short int
      usip=(unsigned short int*)(data);
      retval= (double)usip[cell];
      break;
   case 13: //32 bit unsigned int
      uip=(unsigned int*)(data);
      retval= (double)uip[cell];
      break;
   default:
      throw "Unrecognised data type for DEM. Currently supports 8-bit, both signed and unsigned 16 & 32-bit integer, and 32 & 64-bit float";
      break;
   }

   //Check data is not NULL
   if(retval==file->GetDataIgnoreValue())
   {
      throw "Null value encountered: "+ToString(retval)+". DEMs with a null data value ('data ignore value') cannot yet be used within aplcorr. Please ensure that your DEM has been interpolated to remove any null values and try running again.";
   }

   DEBUGPRINT("Cell and Height from lat,lon: "<<lat<<" "<<lon<<" "<<cell<<" "<<retval)
   return retval;
}


//-------------------------------------------------------------------------
//Return a formatted string containing information about the DEM file
//-------------------------------------------------------------------------
std::string DEM::Info()
{
   //Get the file dimensions
   std::string dims=this->file->FromHeader("samples")+" "+this->file->FromHeader("lines");

   //The format of map_info in bil hdr file is not standard - that is, sometimes there may
   //be more information than others. Will assume only Geographic Lat/Lon to be used.

   //But I'll test it anyway
   std::string mapinfo=this->file->FromHeader("map info");
   if(mapinfo.find("Geographic Lat/Lon")!= std::string::npos)
   {
      //Form the returned string
      std::string deminfo="DEM Information:\n   ";
      deminfo=deminfo+"Number of rows and columns: "+dims+"\n   ";
      deminfo=deminfo+"Projection: "+this->file->FromHeader("map info",0)+"\n   ";
      deminfo=deminfo+"Reference pixel (c,r): "+this->file->FromHeader("map info",1)+", "+this->file->FromHeader("map info",2)+"\n   ";
      deminfo=deminfo+"Reference pixel location (x,y): "+this->file->FromHeader("map info",3)+", "+this->file->FromHeader("map info",4)+"\n   ";
      deminfo=deminfo+"Grid cell spacing (x,y): "+this->file->FromHeader("map info",5)+", "+this->file->FromHeader("map info",6)+"\n   ";
      deminfo=deminfo+"Other information: "+this->file->FromHeader("map info",7)+"\n   ";
      return deminfo;
   }
   else
   {
      std::string deminfo="DEM Information: Unable to format DEM imformation.\n   ";
      deminfo=deminfo+"Number of rows and columns: "+dims+"\n   ";
      deminfo=deminfo+mapinfo;
      return deminfo;
   }
}

//-------------------------------------------------------------------------
//Find the 3 nearest DEM cell centre points to a given (lon,lat) position
//This will be used for creating a plane for testing DEM intersection
//-------------------------------------------------------------------------
bool DEM::GetNearest3Points(const double lon,const double lat,double* const pLat,double* const pLon, double* const pHei)
{
   //Need to find the 3 nearest cell centres to the given point
   //First find the file row,col of the given point - only need integer part here
   double c=X2C(lon);
   double r=Y2R(lat);
   if((r == -1)||(c == -1))
   {
      DEBUGPRINT("Row/Col does not fall within DEM.")
      return false;
   }
   
   //Round the r,c to nearest point (not just truncate)
   c=(int)(c+0.5);
   r=(int)(r+0.5);
   DEBUGPRINT("lat,lon "<<lat<<" "<<lon<<" corresponds to cell row,col: "<<r<<" "<<c)

   //Now convert this back to lat,lon (gives the centre of cell position)
   double clat=R2Y(r);
   double clon=C2X(c);
   DEBUGPRINT("centre of cell r,c "<<r<<" "<<c<<" corresponds to lat,lon: "<<clat<<" "<<clon)
   //Set the cell centre as position 1
   pLat[0]=clat;
   pLon[0]=clon;

   //Now need to find the other 2 nearest cell centres
   //This is based on the sign of the differences between lats and lons
   double latdiff=lat-clat;
   double londiff=lon-clon;

   if(londiff > 0)
   {
      pLon[1]=clon+this->xspace;
      pLat[1]=clat;
   }
   else
   {
      pLon[1]=clon-this->xspace;
      pLat[1]=clat;
   }
  
   if(latdiff > 0)
   {
      pLon[2]=clon;
      pLat[2]=clat+this->yspace;
   }
   else
   {
      pLon[2]=clon;
      pLat[2]=clat-this->yspace;
   }

   //Get the height at these 3 points
   pHei[0]=this->GetHeight(pLon[0],pLat[0]);
   pHei[1]=this->GetHeight(pLon[1],pLat[1]);
   pHei[2]=this->GetHeight(pLon[2],pLat[2]);

   //Test the points were within the area
   if((pHei[0]==DEMOutOfBounds)||(pHei[1]==DEMOutOfBounds)||(pHei[2]==DEMOutOfBounds))
   {
      //Not all points are within the DEM AOI bounds
//      std::cout<<pLon[0]<<" "<<pLat[0]<<std::endl;
//      std::cout<<pLon[1]<<" "<<pLat[1]<<std::endl;
//      std::cout<<pLon[2]<<" "<<pLat[2]<<std::endl;
//      std::cout<<pHei[0]<<" "<<pHei[1]<<" "<<pHei[2]<<std::endl;
      return false;
   }

   return true;
}


//-------------------------------------------------------------------------
// To get DEM slope using the 'neighbourhood method' 
// see [Horn 1981, Srinivasan and Engel 1991, Warren 2004]
//-------------------------------------------------------------------------

//From intersect point
   //Get row,col of nearest cell centre
   //Get height values for 8 neighbours [r-1,r+1,c-1,c+1]

//Cell neighbourhood is described as below for slope calculations:
//    0 1 2 
//    3 4 5
//    6 7 8

//Use neighbourhood method to calculate the DEM slope gradients   
//dz/dx: ((z2 + 2*z5 +z8) - (z0 + 2*z3 +z6)) / (8*dX)  
//dz/dy: ((z0+2*z1+z2) - (z6+2*z7+z8)) / (8*dY)


//-------------------------------------------------------------------------
// Function to calculate the DEM slope and aspect values for a given array
// of latitude and longitude values with 'length' elements.
//-------------------------------------------------------------------------
void DEM::CalculateSlopeAndAzimuth(const double* const lat, const double* const lon, double* const slope, double* const aspect,const int length)
{
   //The correct functional procedure to calculate the DEM slope and aspect is:
   // GetNeighbourhood()
   // CalculateGradient()
   // Slope() , Aspect()
   
   double neighbourhood[9]={0};
   double gradient[2]={0};

   //Because APL DEMs are in lat/lon and DEM heights are in metres we need to scale 
   //the x,y such that they are in metres. This will depend on latitude of the data
   //and the fact the spheroid is flattened not spherical. We can use the same scalar
   //for all points in the array - the difference will be negligable - but will calculate
   //it for each new array (call to this function)
   Ellipsoid ell(WGS84);
   //for xscalar use the reduced latitude to get an estimate
   double beta=atan((ell.b()/ell.a())*tan(lat[0]));
   double xscalar=(PI/180.0)*(ell.a()*cos(beta));
   //For yscalar, use approximation of difference of integrals of meridional distance
   double yscalar=ell.meridional_degree(lat[0]);

   //Loop through each lon/lat pair and get the slope/aspect value
   for(int item=0;item<length;item++)
   {
      GetNeighbourhood(lat[item],lon[item],neighbourhood);
      CalculateGradient(neighbourhood,gradient,xscalar,yscalar);
      slope[item]=Slope(gradient[0],gradient[1])*180/PI;
      //Check if slope is 0, if so then set aspect to 0.
      if(slope[item]==0)
         aspect[item]=0;
      else //calculate the aspect from the gradients
      {
         aspect[item]=Aspect(gradient[0],gradient[1])*180/PI;
         aspect[item]=(90-aspect[item]); //now it goes from North to east to south to west 0-270, and North to west 0 - -90
         if(aspect[item]<0)
            aspect[item] += 360;

      }
   }
}

//-------------------------------------------------------------------------
// Function to return the heights for a neighbourhood of cells from the DEM 
// for a given latitude and longitude (in radians).
//-------------------------------------------------------------------------
void DEM::GetNeighbourhood(const double lat,const double lon,double* const neighbourhood)
{
   //i is a counting variable for the nested for loops
   int i=0;
   //Get the cell column of the DEM relating to longitude value
   double a=X2C(lon*180/PI);
   //Get the cell row of the DEM relating to the latitude value
   double b=Y2R(lat*180/PI);

   //Test if this is an edge cell (this only tests 2 edges - need to test all 4 really)
   if((a<1) || (b<1))
   {
      //We can't calculate slope here due to being on the very edge
      //so just return a zero array
      for(int i=0;i<9;i++)
         neighbourhood[i]=0;
      return;
   }

   //Enter the loop to get the height for each neighbouring cell
   for(int r=b-1;r<=b+1;r++)
   {
      for(int c=a-1;c<=a+1;c++)
      {
         //For each cell in the neighbourhood get the height from the DEM AOI
         neighbourhood[i]=GetHeight(C2X(c),R2Y(r));
         //test if we got the height alright
         if(neighbourhood[i]==DEMOutOfBounds)
         {
            //We didn't get the height value
            std::cout<<"Dem row/col, lat/lon: "<<r<<" "<<c<<" "<<R2Y(r)<<" "<<C2X(c)<<std::endl;
            throw "DEM out of bounds error in DEM::GetNeighbourhood - inspecting a point outside of DEM AOI.";
         }
         //Increase the loop counter for the array index
         i++;
      }
   }
}

//-------------------------------------------------------------------------
// Function to calculate the gradients from the neighbourhood of cells
//-------------------------------------------------------------------------
void DEM::CalculateGradient(const double* const neighbourhood,double* const gradients, const double xscalar, const double yscalar)
{
   //Use neighbourhood method to calculate the DEM slope gradients   
   //dz/dx: ((z2 + 2*z5 +z8) - (z0 + 2*z3 +z6)) / (8*dX)  
   //dz/dy: ((z0+2*z1+z2) - (z6+2*z7+z8)) / (8*dY)

   //Gradient dz/dx
   gradients[0]=((neighbourhood[2] + 2*neighbourhood[5] + neighbourhood[8])
         -(neighbourhood[0] + 2*neighbourhood[3] + neighbourhood[6])) / (8*xspace*xscalar);

   //Gradient dz/dy
   gradients[1]=((neighbourhood[0] + 2*neighbourhood[1] + neighbourhood[2])
         -(neighbourhood[6] + 2*neighbourhood[7] + neighbourhood[8])) / (8*yspace*yscalar);
}

//-------------------------------------------------------------------------
// To return DEM slope from neighbourhood gradient (in radians)
//-------------------------------------------------------------------------
double DEM::Slope(const double dzdx, const double dzdy)
{
   //Use neighbourhood method to calculate the DEM slope   
   // sqrt[ (dz/dx)^2 + (dz/dy)^2 ]  gives the % gradient, so atan to get the angle
   // (in radians)
   return atan(sqrt(dzdx*dzdx + dzdy*dzdy));
}

//-------------------------------------------------------------------------
// To return DEM aspect/azimuth from neighbourhood gradient (in radians)
//-------------------------------------------------------------------------
double DEM::Aspect(const double dzdx, const double dzdy)
{
   //Use the partial results of the above neighbourhood method to get the
   //tangent tan A = (dz/dy / dz/dx)
   //Note if slope gradients are 0 then this is undefined and could return anything.
   //Note2 these are -ve since I think gradients are calculated wrong way around (i.e. A-B instead of B-A)
   return atan2(-dzdy,-dzdx);
}

//-------------------------------------------------------------------------
// Return true if the lat/lon is within epsilon of a DEM cell boundary
//-------------------------------------------------------------------------
bool DEM::OnCellBound(const double lat,const double lon,short* xory=NULL)
{
   //Define a small value for this DEM - use as a bound test
   const double epsilon=std::min(xspace,yspace) / 100.0;
   //Get x,y position normalised by cell size
   double xpos= (lon - refx) / (xspace);
   xpos = xpos - (int)xpos;
   double ypos= (refy - lat) / (yspace);
   ypos = ypos - (int)ypos;

   //Test if these positions are on bounds or not
   if((xpos > -epsilon) && (xpos < epsilon))
   {
      if((ypos > -epsilon) && (ypos < epsilon))
      {
         //Both X and Y are on bounds
         if(xory!=NULL)
            *xory=3;
      }
      else //only x is on bounds
      {
         if(xory!=NULL)
            *xory=1;
      }
      return true;

   }
   else if((ypos > -epsilon) && (ypos < epsilon))
   { 
      // only y is on bounds
      if(xory!=NULL)
         *xory=2;

      return true;
   }

   //Nothing is on the bounds
   if(xory!=NULL)
      *xory=0;

   return false;   
}

