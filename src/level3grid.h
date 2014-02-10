//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef LEVEL3GRID_H
#define LEVEL3GRID_H

#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include "commonfunctions.h"
#include "treegrid_support.h"
#include "basic_igm_worker.h"
#include "dataaccessor.h"

//-------------------------------------------------------------------------
// Class to describe a point in a level 3 grid
//-------------------------------------------------------------------------
class L3Point
{
public:
   L3Point(long r,long c){row=r;col=c;}
   long int row, col;

   static bool compare(L3Point a, L3Point b)
   {
      if(a.row == b.row)
         return a.col < b.col;
      else
         return a.row < b.row;
   }
};

//-------------------------------------------------------------------------
// Class to describe an IGM file point
//-------------------------------------------------------------------------
class IGMPoint
{
public:
   IGMPoint(double x,double y){X=x;Y=y;}
   double X,Y;
};

//-------------------------------------------------------------------------
// Class to describe an outline edge
//-------------------------------------------------------------------------
class Edge
{
public:
   Edge(){minX=0;maxY=0;minY=0;grad=0;}
   int minX,maxY,minY;
   double grad;

   static bool compare(Edge a, Edge b)
   {
      if(a.minY == b.minY)
         return a.minX < b.minX;
      else
         return a.minY < b.minY;
   }
};

//-------------------------------------------------------------------------
// Defines the output grid information e.g. size
//-------------------------------------------------------------------------
class Level3GridInfo
{
public:
   Level3GridInfo(double minX, double minY, double maxX, double maxY,double psx,double psy,std::string bandlist);
   Level3GridInfo(Level3GridInfo* ref);
   ~Level3GridInfo();

   const unsigned long int NumRows()const{return nrows;}
   const unsigned long int NumCols()const{return ncols;}
   const unsigned long int NumBands()const{return nbands;}
   const double PixelSizeX()const{return pixelsizeX;}
   const double PixelSizeY()const{return pixelsizeY;}
   const double TopLeftX()const{return bounds->MinX();}
   const double TopLeftY()const{return bounds->MaxY();}
   const double BottomRightX()const{return bounds->MaxX();}
   const double BottomRightY()const{return bounds->MinY();}
   const unsigned int* Bands()const{return bands;}

   const std::string Wavelengths()const{return wavelengths;}
   const std::string BandList()const{return strBandList;}
   void AddWavelength(const std::string w);

   const bool Inside(double x,double y)const{return bounds->Inside(x,y);}

   std::string GetMapInfo(std::string mapkey);
   int ConvertRC2XY(L3Point* l3p, IGMPoint* igmp);

   void UpdateTopLeftX(double tx){Area copy(bounds); delete bounds; bounds=new Area(tx,copy.MaxX(),copy.MinY(),copy.MaxY());}
   void UpdateTopLeftY(double ty){Area copy(bounds); delete bounds; bounds=new Area(copy.MinX(),copy.MaxX(),copy.MinY(),ty);}

   const Area* const GetBounds()const{return bounds;}

private:
   //Size of cell in terms of X,Y
   double pixelsizeY,pixelsizeX;
   //Size of grid in terms of cells
   unsigned long int ncols,nrows,nbands;
   //Bounds of grid in terms of X,Y
   Area* bounds;

   unsigned int* bands;
   std::string strBandList;   
   std::string wavelengths;

};

//-------------------------------------------------------------------------
// Class that defines a level 3 Grid
//-------------------------------------------------------------------------
class Level3Grid
{
public:
   Level3Grid();
   Level3Grid(Level3GridInfo* gi);
   Level3Grid(double output_pixel_sizeX,double output_pixel_sizeY,std::string bandlist,Area* outputrect);
   ~Level3Grid();

   const unsigned long int NumRows()const{return info->NumRows();}
   const unsigned long int NumCols()const{return info->NumCols();}
   const unsigned long int NumBands()const{return info->NumBands();}
   const double PixelSizeX()const{return info->PixelSizeX();}
   const double PixelSizeY()const{return info->PixelSizeY();}
   const double TopLeftX()const{return info->TopLeftX();}
   const double TopLeftY()const{return info->TopLeftY();}
   const double BottomRightX()const{return info->BottomRightX();}
   const double BottomRightY()const{return info->BottomRightY();}
   const unsigned int* Bands(){return info->Bands();}

   const bool InArea(double a,double b)const{return info->Inside(a,b);}

   const std::string Wavelengths(){return info->Wavelengths();}
   void AddWavelength(const std::string w){info->AddWavelength(w);}
   std::string MapInfo(std::string mapkey) const{return info->GetMapInfo(mapkey);}

   Level3GridInfo* const GetGridInfo()const{return info;}

   int ConvertRC2XY(L3Point* l3p, IGMPoint* igmp){return info->ConvertRC2XY(l3p,igmp);}

private:
   Level3GridInfo* info;

};

//-------------------------------------------------------------------------
// Class to contain an outline of the data in a level 3 grid
//-------------------------------------------------------------------------
class Level3Outline
{
public:
   Level3Outline(Level3GridInfo* const gi,std::string igmfilename);
   //Level3Outline(Level3GridInfo* const gi,double* igmdata,unsigned int nlines,unsigned int nsamples,double ignoreval);
   Level3Outline(Level3GridInfo* const gi,Block<double>* igmblock,unsigned int nlines, unsigned int start_line,double ignoreval);
   ~Level3Outline(){};

   //Write out the outline as an ASCII file
   void WriteEdge(std::string filename);
   void GetEdgeIntersectsOfRow(int row,std::vector<int> &intersects);

private:
   Level3GridInfo* info;
   //vector containing edges for scanlinefill
   std::vector<Edge> edgetable;
   //vector containing list of the active edges
   std::vector<Edge> activeedges;
   inline void InsertEdge(Edge* p1);
   void CreateEdge(L3Point* p1,L3Point* p2);
   void InitialiseForScanlineFill(std::vector<IGMPoint> edges);
   //function to convert the XY to row col
   inline int ConvertXYtoRC(IGMPoint*,L3Point*);
   void ReadEdges(Basic_IGM_Worker* igm,std::vector<IGMPoint>* edges);
   //void Initialise(std::vector<IGMPoint> edges);
   bool firstcalltogetintersects;
   void FakeEdges();
   void AddPointToEdgeArray(DataAccessor<double>* data, double ignoreval,std::vector<IGMPoint> &edges,unsigned int line,unsigned int numsamples);
};

#endif
