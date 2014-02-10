//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef INTERPOLATIONFUNCTIONS_H
#define INTERPOLATIONFUNCTIONS_H

#include "datahandler.h"
#include "commonfunctions.h"
#include <iostream>
#include <string>

//#define DEBUGINTERPFUNC

//Straight forward linear interpolation using one data point either side of desired
void Linear(const double* const times,const int len, DataHandler* const dhandle,NavDataCollection* store,std::string start,std::string stop);

//Smoothes the navigation (raw) data to try and remove any jumps in the data
void Triangle(const unsigned long element,DataHandler* const dhandle,NavDataLine* store,const int kernelsize);

void CubicSpline(const double* const times,const int len,DataHandler* const dhandle,NavDataCollection* store,std::string start,std::string stop);


#endif

