//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#include "sensor.h"

//-------------------------------------------------------------------------
//Constructor with no input 
//-------------------------------------------------------------------------
Sensor::Sensor()
{
}
//-------------------------------------------------------------------------
//Destructor
//-------------------------------------------------------------------------
Sensor::~Sensor()
{

}

//-------------------------------------------------------------------------
// Function to check if sensor id relates to the sensortype s
//-------------------------------------------------------------------------
bool CheckSensorID(SENSORTYPE s,int id)
{
   switch(s)
   {
   case EAGLE:
      for(std::vector<int>::const_iterator it=EAGLE_SENSOR_IDS.begin();it != EAGLE_SENSOR_IDS.end();it++)
      {
         if(id == *it)
         {
            return true; //id is a member of this sensortype
         }
      }
      break;
   case HAWK:
      for(std::vector<int>::const_iterator it=HAWK_SENSOR_IDS.begin();it != HAWK_SENSOR_IDS.end();it++)
      {
         if(id == *it)
         {
            return true; //id is a member of this sensortype
         }
      }
      break;
   default:
      throw "Unknown SENSORTYPE in CheckSensorID";
      break;
   }
   //If it reaches here then id is not from sensortype s
   return false;
}

