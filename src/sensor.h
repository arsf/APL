//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef SENSOR_H
#define SENSOR_H

#include <iostream>
#include <sstream>
#include <map>
#include <string>
#include <limits>
#include <vector>

#include "binfile.h"

//-------------------------------------------------------------------------
//This is an abstract data class so don't try and make an instance of it!
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Sensor object class
//-------------------------------------------------------------------------
class Sensor
{
public:
   Sensor(); //constructor
   virtual ~Sensor(); //destructor

   //string stream object to hold information on status
   std::stringstream ssinfo;

   //-------------------------------------------------------------------------
   //Adding a BIL reader to the class to make reading the raw data and 
   //giving the sensor access to the file dimensions and hdr values rather
   //than havng to pass them all over as variables between the two classes
   //-------------------------------------------------------------------------
   
protected:

};

const unsigned short EAGLE_RAW_MAX=4095;
const double FRAME_TRANSFER_TIME=0.002;

const unsigned short HAWK_RAW_MAX=16383;

const unsigned int RADIANCE_DATA_SCALAR=1000;
const unsigned short CALIBRATED_DATA_MAX=std::numeric_limits<unsigned short>::max();


const int EAGLE_SENSOR_IDS_INTARR[2]={100022,110001};
const int HAWK_SENSOR_IDS_INTARR[2]={300011,310018};

const std::vector<int> EAGLE_SENSOR_IDS(EAGLE_SENSOR_IDS_INTARR,EAGLE_SENSOR_IDS_INTARR + sizeof(EAGLE_SENSOR_IDS_INTARR) / sizeof(int));
const std::vector<int> HAWK_SENSOR_IDS(HAWK_SENSOR_IDS_INTARR,HAWK_SENSOR_IDS_INTARR + sizeof(HAWK_SENSOR_IDS_INTARR) / sizeof(int));

enum SENSORTYPE {EAGLE,HAWK};
bool CheckSensorID(SENSORTYPE s,int id);

#endif
