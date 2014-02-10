//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef BASIC_IGM_WORKER_H
#define BASIC_IGM_WORKER_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <limits>
#include "binfile.h"
#include "bilwriter.h"
#include "logger.h"


class Basic_IGM_Worker
{
public:
   Basic_IGM_Worker(std::string fname);
   Basic_IGM_Worker(Basic_IGM_Worker* b);
   ~Basic_IGM_Worker();

   //Object to read the IGM BIL file
   BinFile* fin;

   //Functions to retrieve the min and max data values
   inline const double MaxX(){return maxx;}
   inline const double MinX(){return minx;}
   inline const double MaxY(){return maxy;}
   inline const double MinY(){return miny;}

   //Functions to retrieve the file dimensions
   inline const unsigned int Samples(){return nsamples;}
   inline const unsigned int Lines(){return nlines;}
   inline const unsigned int Bands(){return nbands;}

   //Function to return projection info from IGM file
   inline const std::string Projection(){return proj;}
   inline const std::string Ellipse(){return ell;}

   std::string FileName(){return filename;}
   uint64_t GetFileSize(){return fin->GetFileSize();}

   bool GetPixelSize(unsigned int pixelid,double* pixsize);

   //Functions to get data from IGM
   double* const GetLine(const unsigned int line);
   double ReadCell(const unsigned int band,const unsigned int line, const unsigned int col);

   bool IsARSFStyle(){return ISARSF;}
   double IgnoreValue(){return nodatavalue;}

private:
   std::string filename;

   double minx,maxx,miny,maxy;
   void GetMinMax();

   unsigned int nsamples,nlines,nbands;
   bool BadPixelSizeCalculation(double* pixsize);
   double* data;
   std::string proj;
   std::string ell;

   bool ISARSF;
   double nodatavalue;
};

#endif
