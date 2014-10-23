//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
//Defines the DEM and associated classes used to handle Digital Elevation 
//Model data. Currently the only supported format for the DEM files is 
//1-band BIL data.
//-------------------------------------------------------------------------

//Also assume DEM values are for centre of cell (I think this is consistent with Grass)
//but this may not be true.

#ifndef DEMS_H
#define DEMS_H

#include <cerrno>
#include <string>
#include <climits>
#include <cmath>
#include "viewvectors.h"
#include "binfile.h"
#include "binaryreader.h"
#include "bil.h"
#include "bsq.h"
#include "commonfunctions.h"
#include "conversions.h"

//#define DEMDEBUG

//-------------------------------------------------------------------------
//Error flag for if DEM is accessed out of bounds
//-------------------------------------------------------------------------
const double DEMOutOfBounds=-99900999;

//-------------------------------------------------------------------------
//Function to round doubles to ints
//-------------------------------------------------------------------------
inline int rounded(double x){return int(x+0.5);}

//-------------------------------------------------------------------------
//Used with the AOI to get the value of Lower Left and Upper Right vertices
//-------------------------------------------------------------------------
enum vertex {LLX,LLY,URX,URY};

//-------------------------------------------------------------------------
/*
   Area of Interest Class
   Stores lower left x,y and upper right x,y coordinates
   Can be used directly with DEM to read in a specific part of the file
*/
//-------------------------------------------------------------------------

class DEM_AOI {
public:

   //Constructors
   DEM_AOI();
   DEM_AOI(const double illx, const double illy,const double iurx,const double iury);

   //Set the values of the lower left and upper right vertices of the AOI
   bool Set(const double illx, const double illy,const double iurx,const double iury);

   //Get the value of a component of a vertex
   double Get(const vertex v);

private:
   //The components of the lower left and upper right vertices
   double llx,lly,urx,ury;
};


//-------------------------------------------------------------------------
// DEM Binary File class that derives from BinFile to allow BSQ or BIL DEMs
//-------------------------------------------------------------------------

class DEMBinFile: public BinFile{
public:
   //New constructor to check that only one band exists in the file
   DEMBinFile(std::string strFilename) : BinFile(strFilename)
   {
      //Check there is only 1 band in the BIL file
      if(this->FromHeader("bands") != "1")
         throw "Expected DEM to be a 1-band BIL file. Number of bands reported in hdr is not 1";

      //Get the data ignore value if it exists
      if(FromHeader("data ignore value")!="")
         dataignore=StringToDouble(FromHeader("data ignore value"));
      else
         dataignore=-99999999; //need to give it a number - this is unlikely to ever be used as real DEM values 
   } 

   int ReadRect(char* const chdata, const double minrow, const double maxrow,const double mincol,const double maxcol)
   {
      return br->ReadRect(chdata,minrow,maxrow,mincol,maxcol);
   }

   double GetDataIgnoreValue() const {return dataignore;}

private:
   double dataignore;
};



//-------------------------------------------------------------------------
// Digital elevation model handling class
//-------------------------------------------------------------------------

class DEM {
public:

   //constructor: pass it the DEM filename
   DEM(std::string strFilename);
   ~DEM(); //destructor

   //Read in a rectangle worth of data from the DEM into chdata
   //using the DEM AOI object
   int ReadRect(char* const chdata); 

   //Read in a rectangle worth of data from the DEM into chdata where
   //area is defined by lower left, upper right corners (in lat/lon).
   int ReadRect(char* const chdata,const double llx, const double lly,const double urx,const double ury); 

   //Get the size (in bytes) required to read in a rectangular area
   //using the DEM AOI object
   unsigned long int SizeOf();

   //Get the size (in bytes) required to read in a rectangular area where
   //area is defined by lower left, upper right corners (in lat/lon).
   unsigned long int SizeOf(const double llx, const double lly,const double urx,const double ury);

   //Functions for interaction with the Area Of Interest (AOI) object
   bool SetAOI(const double llx, const double lly,const double urx,const double ury);
   double GetAOI(const vertex v);

   //Read some data into the char* data using the AOI to define the area to be imported 
   void FillArray();

   //Functions to get the height from an AOI cell at lon/lat
   unsigned int GetAOICell(const double lon,const double lat);
   double GetHeight(const double lon,const double lat);

   //Functions to return the grid spacing of the DEM
   double GetXSpace()const {return xspace;}
   double GetYSpace()const {return yspace;}

   //Function to return a formatted string containing information about the DEM
   std::string Info();

   //Find the 3 nearest DEM cell centre points to a given (lon,lat) position
   //This will be used for creating a plane for testing DEM intersection
   bool GetNearest3Points(const double lon,const double lat,double* const pLat,double* const pLon, double* const pHei);

   double X2C(const double x); //r,c are doubles not int for future interpolation to sub-pixel
   double Y2R(const double y); //r,c are doubles not int for future interpolation to sub-pixel
   double C2X(const double c);
   double R2Y(const double r);

   //Function to calculate the slope and aspect values for given lat/lon arrays
   void CalculateSlopeAndAzimuth(const double* const lat, const double* const lon, double* const slope, double* const aspect,const int length);

   bool OnCellBound(const double lat,const double lon,short* xory);

private:

   //To read in the DEM file include a DEM BIL Reader
   DEMBinFile* file;

   //Area of Interest to describe part of the DEM for working with
   DEM_AOI AOI;

   //Edge bounds of the DEM
   double maxy,maxx,miny,minx; 

   //Cell spacing of the DEM (e.g. DEM pixel resolution)
   double xspace,yspace; 

   //Reference pixel row,col and x,y values 
   double refx,refy,refr,refc;

   //Check that the DEM supplied is of a suitable format
   void CheckSupport(); 

   //Array to store DEM data in
   char* data;

   //To speed up functions store number of rows/cols of dem instead of accessing from header
   unsigned int ncols,nrows;

   //To extend the AOI to fit the grid lines of the DEM
   bool FitAOIToGrid();

   //Functions involved with generating the slope and aspect values
   void GetNeighbourhood(const double lat,const double lon,double* const neighbourhood);
   void CalculateGradient(const double* const neighbourhood,double* const gradients,const double xscalar,const double yscalar);
   double Slope(const double dzdx, const double dzdy);
   double Aspect(const double dzdx, const double dzdy);

};

#endif
