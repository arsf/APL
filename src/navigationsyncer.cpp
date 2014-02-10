//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#include "navigationsyncer.h"

#ifdef DEBUGNAVSYNCER
   #define DEBUGPRINT(x) std::cout<<x<<std::endl;
#else
   #define DEBUGPRINT(x)
#endif


//-------------------------------------------------------------------------
// Function to return the leap seconds for the given date
// Date should be in format dd-mm-yyyy
//-------------------------------------------------------------------------
int LeapSecond::GetLeapSeconds(std::string mydate)
{
   //Check date format
   if(!CheckDateFormat(mydate))
   {
      throw "Acquisition date in unexpected format. Require 'dd-mm-yyyy' but received:"+mydate;
   }

   //Set up 2 empty time structs
   struct tm mytime={0};
   struct tm testday={0};
   //Fill in the time struct for the date of the data
   FillTimeStruct(&mytime,mydate);
   //Now get the seconds since epoch for data collection time (to the nearest day)
   time_t collectiontime = mktime(&mytime);
   //Check while myday is greater than date
   std::map<std::string,int>::iterator it=data.begin();

   for(it=data.begin();it!=data.end();it++)
   {
      DEBUGPRINT("Leap second data value "<<(*it).first<<" "<<(*it).second)  
      FillTimeStruct(&testday,(*it).first);
      //get the seconds since epoch for leap second increase day
      time_t leaptime=mktime(&testday);

      //Test date is not before start of data 
      if((it==data.begin())&&(collectiontime < leaptime))
         throw "Given date is before the first date in leap second class, received:"+mydate+" and first date is:"+(*it).first;

      if(collectiontime < leaptime)
      {
         break;
      } 
   }


   //Return the previous entry's second value
   return (*(--it)).second;
}

//-------------------------------------------------------------------------
// Function to check the format of the date string matches the dd-mm-yyyy format
//-------------------------------------------------------------------------
bool LeapSecond::CheckDateFormat(std::string testdate)
{
   std::string nicechars="0123456789-";
   //Test if any chars not of nicechars in date string
   if(testdate.find_first_not_of(nicechars)!=std::string::npos)
      return false;

   //Test there are only 2 -'s
   if(TotalOccurence(testdate,'-')!=2)
      return false;   

   //Test if - are in correct place   
   if((testdate.find_first_of('-')!=2)&&(testdate.find_last_of('-')!=5))
      return false;

   return true;
}

//-------------------------------------------------------------------------
// NavigationSyncer constructor using nav file and lev1 file to set up
//-------------------------------------------------------------------------
NavigationSyncer::NavigationSyncer(std::string navfilename, std::string lev1filename)
{
   //Set to NULL here anyway just to be safe - they should all be none-null by the end of this function
   nscans=0;
   time=NULL;
   navfile=NULL;
   hdrsync=0;
   framerate=0;
   lev1firstscanmaxexpectedsize=30; //30 seconds 

   if(navfilename.compare("NULL")!=0)
   {
      //Create a new Specim nav file object
      navfile=new SpecimFileChooser(navfilename);
      //Read in the nav file to get the time syncs
      navfile->Reader();
   }

   //Need to read in the level1 file hdr info to get nscans
   BinFile bilin(lev1filename);
   //Get the number of scan lines
   nscans=StringToUINT(bilin.FromHeader("lines"));
   DEBUGPRINT("Number of scans: "<<nscans)
   //Get the NavSync timing from the hdr
   hdrsync=(StringToUINT(TrimWhitespace(bilin.FromHeader("NavSync Timing",1)))/1000.0);
   DEBUGPRINT("Sync from header:"<<hdrsync)
   //Get the acquisition date of the data
   acquisitiondate=bilin.FromHeader("acquisition date");
   //Remove the start of the date string 
   acquisitiondate=TrimWhitespace(acquisitiondate.substr(acquisitiondate.find_first_of(':')+1));
   DEBUGPRINT("Date: "<<acquisitiondate)
   //Get the start and stop times from the header to use in case no sync messages 
   //found in the specim nav file
   gpsstarttime=bilin.FromHeader("GPS Start Time");
   gpsstarttime=RemoveAllBut(gpsstarttime,"1234567890.:");
   gpsstarttime=TrimWhitespace(ReplaceAllWith(&gpsstarttime,':',' '));
   gpsstoptime=bilin.FromHeader("GPS Stop Time");
   gpsstoptime=RemoveAllBut(gpsstoptime,"1234567890.:");
   gpsstoptime=TrimWhitespace(ReplaceAllWith(&gpsstoptime,':',' '));
   
   //Get the frame rate from the hdr
   framerate=StringToDouble(bilin.FromHeader("fps"));
   DEBUGPRINT("Frame rate from header:"<<framerate)

   //Get the Processed_crop_start from the header file - this tells us if 
   //nav for the full line or crop of the line is required
   std::string cropstart=bilin.FromHeader("y start");
   if(cropstart.compare("")==0)
   {
      Logger::Warning("No y start found in level 1 header, if data was cropped in previous stages navigation may be wrongly synced.");
      //Set the time offset to 0
      croptimeoffset=0;
   }
   else
   {
      //Convert to a double
      croptimeoffset=StringToDouble(cropstart);
      //Get the number of dropped scans that occurred in the crop prior to y start (if y start = 0 so will this)
      std::string prevdropscans=bilin.FromHeader("dropped scans before y start");
      if(prevdropscans.compare("")==0)
      {
         Logger::Warning("No 'dropped scans before y start' found in level 1 header, if y start is non-zero navigation may be wrongly synced.");
         //Set the time offset to 0
         prevdropscans="0";
      }
      double previousdroppedscans=StringToDouble(prevdropscans);
      //Convert the sum of the frames (cropstart and prevdropscans) to a time offset to add onto start time
      croptimeoffset=(croptimeoffset + previousdroppedscans)/framerate;
      Logger::Log("Using cropped level-1 data - will add a time offset relating to number of lines cropped (y start + dropped scans values in hdr): "+ToString(croptimeoffset));
   }

   //Close the level 1 file
   bilin.Close();

   //Create the scan time array
   time=new double[nscans];

   //Get the number of leap seconds for the data
   LeapSecond leap;
   leapseconds=leap.GetLeapSeconds(acquisitiondate);
   DEBUGPRINT("Using leap seconds of:"<<leapseconds);
}

//-------------------------------------------------------------------------
// NavigationSyncer destructor
//-------------------------------------------------------------------------
NavigationSyncer::~NavigationSyncer()
{
   if(navfile!=NULL)
      delete navfile;
   if(time!=NULL)
      delete[] time;
}


//-------------------------------------------------------------------------
// Function to get times for each scan line
//-------------------------------------------------------------------------
void NavigationSyncer::FindScanTimes()
{
   //Get the time interval, in seconds, between scan lines
   //nscans should include ALL scans i.e. dropped scans too

   //A negative usersynctime means that it will use the hdr/nav files
   //to find the scan times, +ve means the first scan line will have time usersynctime

   double scanseparation= 1.0 /framerate;

   //Now need to find the time for the first scan line
   //To do this we need the time from the nav file when the sync occurred
   int syncnumber=0;
   std::vector<int> syncrecord; //use a vector to hold incase of multiple entries
   bool syncfound=false;
   
  // if((usersynctime<0)&&(navfile!=NULL)) // Use the nav/hdr files to find the sync message
   if(navfile!=NULL) // Use the nav/hdr files to find the sync message
   {

      //Find which sync relates to this line if multiple sync messages
      for(syncnumber=0;syncnumber<navfile->GetNumSyncs();syncnumber++)
      {  
         if(navfile->GetSyncTime(syncnumber) == hdrsync)
         {
            syncfound=true;
            //add the syncnumber to the vector of sync records
            syncrecord.push_back(syncnumber);
            //Now check to see if there are other identical sync messages
         }
      }

      //Did we escape the loop 
      if(syncfound == false)
      {
         throw "Error in NavigationSyncer::FindScanTimes() - No matching sync message found in .nav file that matches value in hdr file. Re-run with -nonav and manually find the sync.";
      }
      else //we found a matching sync message
      {
         int syncrecord_index=0;
         //If there are more than 1 records found then syncrecord_index may be updated in the if statement below
         //else just use the value of 0
         if(syncrecord.size() > 1)
         {
            //We need to test which of the sync messages is for this file - can look at GPS time of sync message
            //and select the one that is closest to the header start time
            double lev1starttime=GetSecOfWeek(acquisitiondate,ReplaceAllWith(&gpsstarttime,':',' ')); 
            int mindiff=std::numeric_limits<int>::max();
            //Need to loop through all syncs to find the one closest to start time
            for(unsigned int vecindex=0;vecindex<syncrecord.size();vecindex++)
            {
               if( abs(lev1starttime - navfile->GetGPSSync(syncrecord[vecindex])) < mindiff)
               {
                  mindiff=abs(lev1starttime - navfile->GetGPSSync(syncrecord[vecindex]));
                  syncrecord_index=vecindex;
               }
            }    
            //Now have a second test here to make sure that the selected index is within X seconds 
            //of the level 1 start time. This is hard coded below!
            if(mindiff > lev1firstscanmaxexpectedsize)
            {
               //Closest is still greater than "lev1firstscanmaxexpectedsize" seconds away - just list choices and let user select.
               Logger::Log("Multiple identical sync delay values in Specim nav file. None fall within the "+ToString(lev1firstscanmaxexpectedsize)+" seconds window of start time of level 1 file.");
               Logger::Log("Header start time: "+ToString(lev1starttime)+" sync delay value: "+ToString(hdrsync));
               for(unsigned int vecindex=0;vecindex<syncrecord.size();vecindex++)   
               {
                  //The value to add onto the hdr start time to get the "synced" time assuming the hdr start time is used as scan line 0 time
                  double scntoff=(navfile->GetGPSSync(syncrecord[vecindex]) -hdrsync) - lev1starttime;
                  Logger::Log("Sync message time: "+ToString(navfile->GetGPSSync(syncrecord[vecindex]))+" Scantimeoffset value to use: "+ToString(scntoff));
               }
               throw "Multiple possibilities for sync time. Try using one of the suggested values above as a -scantimeoffset and processing with -nonav";
            }
         }
         //output a message for the sync index used
         Logger::Log("Using the sync message from index "+ToString(syncrecord[syncrecord_index])+" which probably means it is for flight line "+ToString(syncrecord[syncrecord_index])+" in this specim nav file (note these are referenced from 0 not 1).");
         //Get the GPS time of the sync message
         double gpssync=navfile->GetGPSSync(syncrecord[syncrecord_index]);
         //Convert this into seconds of day
         //gpssync=((int)gpssync % (24*3600));

         //Subtract the sync message time
         double firstscantime=gpssync-hdrsync;

         //Add on GPS leap seconds to get into UTC if from Specim Nav file
         //firstscantime+=leapseconds;

         //Add on the cropped start time offset - this is 0 if no cropping (at the start of the line) has taken place
         Logger::Log("Applying crop time offset of: "+ToString(croptimeoffset));
         firstscantime=firstscantime + croptimeoffset;         

         //Compare to start time of level-1 header - in GPS sec of week - to prevent data from a different weekday being used 
         CompareLev1ToNavTimes(firstscantime);

         //Set the per scan times
         for(unsigned int i=0;i<nscans;i++)
            time[i]=firstscantime + i*scanseparation;

         Logger::Log("First scan will have time: "+ToString(time[0])+".  Last scan will have time: "+ToString(time[nscans-1]));
      }
   }
   else //Use a user-specified value as the first scan time
   {
      double firstscantime=0;
      //if(usersynctime>=0) 
      //{
      //   //Use the given sync time in conjunction with the start time from the header
      //   firstscantime=GetSecOfWeek(acquisitiondate,ReplaceAllWith(&gpsstarttime,':',' '));
      //   firstscantime=firstscantime+usersynctime;
      //}
      //else
      {
         //Use the start time from the header file as the first scan time
         firstscantime=GetSecOfWeek(acquisitiondate,ReplaceAllWith(&gpsstarttime,':',' '));
         Logger::Log("Using a first scan time of "+ToString(firstscantime));

         //Add on the cropped start time offset - this is 0 if no cropping (at the start of the line) has taken place
         Logger::Log("Applying crop time offset of: "+ToString(croptimeoffset));
         firstscantime=firstscantime + croptimeoffset;   
      }

      //Set the per scan times
      for(unsigned int i=0;i<nscans;i++)
         time[i]=firstscantime + i*scanseparation;

      DEBUGPRINT("First scan time: "<<time[0]<<" Last scan time: "<<time[nscans-1]);

   }
}

//-------------------------------------------------------------------------
// Function to apply a given time shift to the scan times
//-------------------------------------------------------------------------
void NavigationSyncer::ApplyTimeShift(const double shift)
{
   for(unsigned int i=0;i<nscans;i++)
   {
      time[i]+=shift;
   }
}

//-------------------------------------------------------------------------
// Function to compare the level 1 start time (GPS SecOfWeek) with nav data times
//-------------------------------------------------------------------------
void NavigationSyncer::CompareLev1ToNavTimes(double first_scan_time)
{
   const int numsecsperday=3600*24;
   //Get the level 1 start time in GPS sec of week
   double lev1starttime=GetSecOfWeek(acquisitiondate,ReplaceAllWith(&gpsstarttime,':',' ')); 
   DEBUGPRINT("GPS Second of week for level1 start time: "<<lev1starttime)

   Logger::Log("Difference between start time from level-1 header file and start time from navigation: "+ToString(lev1starttime-first_scan_time)+" seconds.");

   if(fabs(lev1starttime-first_scan_time) > lev1firstscanmaxexpectedsize)
   {
      Logger::Warning("This time difference appears excessive: "+ToString(lev1starttime-first_scan_time));
   }

   //Do some comparisons
   // if the end time of the nav data is greater than the start time from the lev1 header
   if(navfile->GetLine(navfile->GetNumEntries()-1)->time < lev1starttime)
   {
         throw "Error: The level-1 start time is after the end of the navigation data: nav end time: "+ToString(navfile->GetLine(navfile->GetNumEntries()-1)->time)+" lev 1 time: "+ToString(lev1starttime);
   }
   //Test if data is more than a day (of week) out from both start and end times of nav data
   if((fabs(navfile->GetLine(0)->time - lev1starttime) > numsecsperday )&&(fabs(navfile->GetLine(navfile->GetNumEntries()-1)->time - lev1starttime) > numsecsperday))
   {
      throw "Error: The level-1 start time is on a different week day to both the navigation start time and end times";
   }
   else if (fabs(navfile->GetLine(0)->time - lev1starttime) > numsecsperday)
   {
      Logger::Warning("The level-1 start time is on a different week day to the navigation start time.");
   }
   else if(fabs(navfile->GetLine(navfile->GetNumEntries()-1)->time - lev1starttime) > numsecsperday)
   {
      Logger::Warning("The level-1 start time is on a different week day to the navigation end time.");
   }

}
