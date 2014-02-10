//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#include "datahandler.h"

#ifdef DEBUGDATAHANDLER
   #define DEBUGPRINT(x) std::cout<<x<<std::endl;
#else
   #define DEBUGPRINT(x)
#endif

//FIXME TODO Remove the GLOBAL_FORCE variable and do this properly
extern bool GLOBAL_FORCE;

NavDataLine::NavDataLine()
{
   lat=0;
   lon=0;
   time=0;
   hei=0;
   roll=0;
   pitch=0;
   heading=0;
   quality=0;
}


//-------------------------------------------------------------------------
// Constructor for NavDataCollection
//-------------------------------------------------------------------------
NavDataCollection::NavDataCollection(unsigned long length)
{
   navarray=NULL;
   size_of_array=length;
   navarray=new NavDataLine[size_of_array];
}

//-------------------------------------------------------------------------
// Destructor for NavDataCollection
//-------------------------------------------------------------------------
NavDataCollection::~NavDataCollection()
{
   if(navarray!=NULL)
      delete[] navarray;
}

//-------------------------------------------------------------------------
// Set the value for an item in the NavDataCollection
//-------------------------------------------------------------------------
void NavDataCollection::SetValues(const unsigned long item,const NavDataLine* line)
{
   //Check that the item can be set
   if((item<size_of_array)&&(navarray!=NULL))
   {
      navarray[item].lat=line->lat;
      navarray[item].lon=line->lon;
      navarray[item].time=line->time;
      navarray[item].roll=line->roll;
      navarray[item].pitch=line->pitch;
      navarray[item].heading=line->heading;
      navarray[item].hei=line->hei;
      navarray[item].quality=line->quality;
   }
}

//-------------------------------------------------------------------------
// Set the single value for an item in the NavDataCollection
//-------------------------------------------------------------------------
void NavDataCollection::SetValues(const unsigned long item,const NavDataItem key,const double value)
{
   //Check that the item can be set
   if((item<size_of_array)&&(navarray!=NULL))
   {
      switch(key)
      {
      case LAT:
         navarray[item].lat=value;
         break;
      case LON:
         navarray[item].lon=value;
         break;
      case TIME:
         navarray[item].time=value;
         break;
      case ROLL:
         navarray[item].roll=value;
         break;
      case PITCH:
         navarray[item].pitch=value;
         break;
      case HEADING:
         navarray[item].heading=value;
         break;
      case HEI:
         navarray[item].hei=value;
         break;
      case QUALITY:
         navarray[item].quality=static_cast<char>(value);
         break;
      default:
         throw "Unknown NavDataItem in SetValues()";
      }
   }
}

//-------------------------------------------------------------------------
// Return a colection items value based on given key word
//-------------------------------------------------------------------------
const double NavDataCollection::GetValue(const unsigned long item,const NavDataItem key)
{
   //Check that the item can be got
   if((item<size_of_array)&&(navarray!=NULL))
   {
      switch(key)
      {
      case LAT:
         return navarray[item].lat;
      case LON:
         return navarray[item].lon;
      case TIME:
         return navarray[item].time;
      case ROLL:
         return navarray[item].roll;
      case PITCH:
         return navarray[item].pitch;
      case HEADING:
         return navarray[item].heading;
      case HEI:
         return navarray[item].hei;
      default:
         throw "Unknown NavDataItem in GetValue()";
      }
   }
   return 0;
}

//-------------------------------------------------------------------------
// Return a reference to the specified value
//-------------------------------------------------------------------------
double* NavDataCollection::GetReferenceValue(const unsigned long item,const NavDataItem key)
{
   //Check that the item can be got
   if((item<size_of_array)&&(navarray!=NULL))
   {
      switch(key)
      {
      case LAT:
         return &navarray[item].lat;
      case LON:
         return &navarray[item].lon;
      case TIME:
         return &navarray[item].time;
      case ROLL:
         return &navarray[item].roll;
      case PITCH:
         return &navarray[item].pitch;
      case HEADING:
         return &navarray[item].heading;
      case HEI:
         return &navarray[item].hei;
      default:
         throw "Unknown NavDataItem in GetValue()";
      }
   }
   return NULL;
}

//-------------------------------------------------------------------------
// Function to check the plausibilty of the data in the Collection
// This is really only useful if the data is consecutive .e.g. as if 
// read in from the a file and stored as in the file
//-------------------------------------------------------------------------
void NavDataCollection::CheckPlausibility()
{
   //The first check to make is to see if each consecutive entry
   //is within a plausible distance from previous

   bool badt=false,badr=false,badp=false,badhead=false,badlo=false,badla=false,badhei=false;

   //Assume first flag is good since we have no previous point to test against
   navarray[0].quality=0;
   for(unsigned int epoch=1;epoch<size_of_array;epoch++)
   {
      //Set flag to zero before running
      navarray[epoch].quality=0;
      if ((navarray[epoch].time - navarray[epoch-1].time) > plausibleTimeDifference)
      {
         if(badt==false)
         {
            badt=true;
            Logger::Log("Time difference between consecutive epochs larger than acceptable threshold. Further warnings of this type supressed.");
            DEBUGPRINT("Time difference between consecutive epochs larger than acceptable threshold. Expected < "
                     <<plausibleTimeDifference<<", got  "<<(navarray[epoch].time - navarray[epoch-1].time)
                     <<" for epoch "<<epoch)
         }
         navarray[epoch].quality+=NavDataLine::BADTIME;
      }

      //Throw an exception if time goes backwards - this is not allowed
      if((navarray[epoch].time - navarray[epoch-1].time) < 0)
      {
         if(GLOBAL_FORCE == false)
         {
            throw "Time goes backwards in navigation file at epoch "+ToString(epoch)+"\nPrevious time: "
                  +ToString(navarray[epoch-1].time)+" Current time: "+ToString(navarray[epoch].time) +
                  "\nIf this occurs in a .nav file then it is probably OK to continue if you are not outputting the real-time navigation."
                  "\nTo try and force APL to continue use the -force command line option.";
         }
         else
         {
            Logger::Warning("Time goes backwards in navigation file at epoch "+ToString(epoch)+"\nPrevious time: "
                  +ToString(navarray[epoch-1].time)+" Current time: "+ToString(navarray[epoch].time) +
                  "\nThis is being ignored as user has specified the -force command line option ...");            
         }
      }
      
      if(fabs(navarray[epoch].hei - navarray[epoch-1].hei) > plausibleHeightDifference)
      {

         if(badhei==false)
         {
            badhei=true;
            Logger::Log("Height difference between consecutive epochs larger than acceptable threshold. Further warnings of this type supressed.");
            DEBUGPRINT("Height difference between consecutive epochs larger than acceptable threshold. Expected < "
                     <<(plausibleHeightDifference)<<", got  "<<(navarray[epoch].hei - navarray[epoch-1].hei)
                     <<" for epoch "<<epoch);
            DEBUGPRINT("Epoch: "<<epoch-1<<" Time: "<<navarray[epoch-1].time<<" Height: "<<navarray[epoch-1].hei
            <<" Lat: "<<navarray[epoch-1].lat<<" Lon: "<<navarray[epoch-1].lon<<" Roll: "<<navarray[epoch-1].roll
            <<" Pitch: "<<navarray[epoch-1].pitch<<" Heading: "<<navarray[epoch-1].heading);
            DEBUGPRINT("Epoch: "<<epoch<<" Time: "<<navarray[epoch].time<<" Height: "<<navarray[epoch].hei
            <<" Lat: "<<navarray[epoch].lat<<" Lon: "<<navarray[epoch].lon<<" Roll: "<<navarray[epoch].roll
            <<" Pitch: "<<navarray[epoch].pitch<<" Heading: "<<navarray[epoch].heading);
         }
         navarray[epoch].quality+=NavDataLine::BADHEI;
      }
      
      if(fabs(navarray[epoch].lat - navarray[epoch-1].lat)> plausibleLatDifference)
      {
         if(badla==false)
         {
            badla=true;
            Logger::Log("Latitude difference between consecutive epochs larger than acceptable threshold. Further warnings of this type supressed.");
            DEBUGPRINT("Latitude difference between consecutive epochs larger than acceptable threshold. Expected < "
                     <<(plausibleLatDifference)<<", got  "<<(navarray[epoch].lat - navarray[epoch-1].lat)
                     <<" for epoch "<<epoch);
         }
         navarray[epoch].quality+=NavDataLine::BADLAT;
      }

      if(fabs(navarray[epoch].lon - navarray[epoch-1].lon) > plausibleLonDifference)
      {
         if(badlo==false)
         {
            badlo=true;
            Logger::Log("Longitude difference between consecutive epochs larger than acceptable threshold. Further warnings of this type supressed.");
            DEBUGPRINT("Longitude difference between consecutive epochs larger than acceptable threshold. Expected < "
                     <<(plausibleLonDifference)<<", got  "<<(navarray[epoch].lon - navarray[epoch-1].lon)
                     <<" for epoch "<<epoch);
         }
         navarray[epoch].quality+=NavDataLine::BADLON;
      }

      if(fabs(navarray[epoch].roll - navarray[epoch-1].roll) > plausibleRollDifference)
      {
         if(badr==false)
         {
            badr=true;
            Logger::Log("Roll difference between consecutive epochs larger than acceptable threshold. Further warnings of this type supressed.");
            DEBUGPRINT("Roll difference between consecutive epochs larger than acceptable threshold. Expected < "
                     <<(plausibleRollDifference)<<", got  "<<(navarray[epoch].roll - navarray[epoch-1].roll)
                     <<" for epoch "<<epoch);
         }
         navarray[epoch].quality+=NavDataLine::BADROLL;
      }

      if(fabs(navarray[epoch].pitch - navarray[epoch-1].pitch) > plausiblePitchDifference)
      {
         if(badp==false)
         {
            badp=true;
            Logger::Log("Pitch difference between consecutive epochs larger than acceptable threshold. Further warnings of this type supressed.");
            DEBUGPRINT("Pitch difference between consecutive epochs larger than acceptable threshold. Expected < "
                     <<(plausiblePitchDifference)<<", got  "<<(navarray[epoch].pitch - navarray[epoch-1].pitch)
                     <<" for epoch "<<epoch);
         }
         navarray[epoch].quality+=NavDataLine::BADPITCH;
      }

      //Two tests for heading since it wraps around 0->360 degrees
      if((fabs(navarray[epoch].heading - navarray[epoch-1].heading) > plausibleHeadingDifference)
      &&( (fabs(navarray[epoch].heading - navarray[epoch-1].heading))-360 > plausibleHeadingDifference ))
      {
         if(badhead==false)
         {
            badhead=true;
            Logger::Log("Heading difference between consecutive epochs larger than acceptable threshold. Further warnings of this type suppressed.");
            DEBUGPRINT("Heading difference between consecutive epochs larger than acceptable threshold. Expected < "
                     <<(plausibleHeadingDifference)<<", got  "<<(navarray[epoch].heading - navarray[epoch-1].heading)
                     <<" for epoch "<<epoch);
         }
         navarray[epoch].quality+=NavDataLine::BADHEADING;
      }
   }
}


//-------------------------------------------------------------------------
// Default constructor for DataHandler
//-------------------------------------------------------------------------
DataHandler::DataHandler()
{
   navcollection=NULL;
}

//-------------------------------------------------------------------------
// Constructor if number of entries in data file is known
//-------------------------------------------------------------------------
DataHandler::DataHandler(unsigned long n)
{
   navcollection=new NavDataCollection(n);
}

DataHandler::~DataHandler()
{
   if(navcollection!=NULL)
      delete navcollection;
}


//-------------------------------------------------------------------------
// Function to smooth input data to remove high frequency jitters
//-------------------------------------------------------------------------
void DataHandler::Smooth(void (*f)(const unsigned long ,DataHandler* ,NavDataLine* ,const int),const unsigned int smoothkernelsize)
{
   //DEBUGPRINT("Entering DataHandler::Smooth(x)...")
   NavDataLine* smoothednav=new NavDataLine[GetNumEntries()];
   unsigned int morethanhalf=(smoothkernelsize+1)/2.0;
   //FIXME fix the morethanhalf below so that everything is smoothed - this is tricky
   //The morethanhalf come from the length of the kernel used in the smoothing function
   for(unsigned long i=morethanhalf;i<GetNumEntries()-morethanhalf;i++)
   {
      Smooth(f,i,&smoothednav[i],smoothkernelsize);
   }

   //For the elements that couldn't be smoothed we shall just insert the original data since zero-ing causes extra problems
   //Insert from the front and back on the same loop
   for(unsigned long i=0,j=GetNumEntries()-1;i<morethanhalf;i++,j--)
   {
      smoothednav[i].time=navcollection->GetValue(i,NavDataCollection::TIME);
      smoothednav[i].lat=navcollection->GetValue(i,NavDataCollection::LAT);
      smoothednav[i].lon=navcollection->GetValue(i,NavDataCollection::LON);
      smoothednav[i].hei=navcollection->GetValue(i,NavDataCollection::HEI);
      smoothednav[i].roll=navcollection->GetValue(i,NavDataCollection::ROLL);
      smoothednav[i].pitch=navcollection->GetValue(i,NavDataCollection::PITCH);
      smoothednav[i].heading=navcollection->GetValue(i,NavDataCollection::HEADING);

      smoothednav[j].time=navcollection->GetValue(j,NavDataCollection::TIME);
      smoothednav[j].lat=navcollection->GetValue(j,NavDataCollection::LAT);
      smoothednav[j].lon=navcollection->GetValue(j,NavDataCollection::LON);
      smoothednav[j].hei=navcollection->GetValue(j,NavDataCollection::HEI);
      smoothednav[j].roll=navcollection->GetValue(j,NavDataCollection::ROLL);
      smoothednav[j].pitch=navcollection->GetValue(j,NavDataCollection::PITCH);
      smoothednav[j].heading=navcollection->GetValue(j,NavDataCollection::HEADING);
   }
   
   //Now copy the smoothed navigation in place of the original (do a hard copy)
   for(unsigned int i=0;i<GetNumEntries();i++)
   {
      navcollection->SetValues(i,&smoothednav[i]);
   }
   delete[] smoothednav;
}

//-------------------------------------------------------------------------
// Function to apply smoothing function to data element 
//-------------------------------------------------------------------------
void DataHandler::Smooth(void (*f)(const unsigned long ,DataHandler* ,NavDataLine* ,const int),const int element, NavDataLine* store,const unsigned int smoothkernelsize)
{
   //FIXME change the kernel size to me a parameter
   f(element,this,store,smoothkernelsize);
}

