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
   NOSYNCINHDR=-999;
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
   try
   {
      hdrsync=StringToUINT(TrimWhitespace(bilin.FromHeader("NavSync Timing",1,"true")))/1000.0; 
   }
   catch(std::string e)
   {
      if(e.compare(0,bilin.MissingHeaderItemError().length(),bilin.MissingHeaderItemError())==0)
      {
         //Set to a value that means "no value in header"
         hdrsync=NOSYNCINHDR;
      }
      else
      {
         throw e;         
      }
   }

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
   if((framerate <= 0)||(framerate>100000))
   {
      throw "Frame rate (fps) in hdr file seems erroneous - will only process for frame rates >0 and <100000.";
   }

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
   
   if(navfile!=NULL) // Use the nav/hdr files to find the sync message
   {
      //Check for single versus multiple time stamps
      bool persec=navfile->UsePerSecondForSync();
      int syncrecord_index=0;

      //We need to test which of the sync messages is for this file - can look at GPS time of sync message
      //and select the one that is closest to the header start time
      double lev1starttime=GetSecOfWeek(acquisitiondate,ReplaceAllWith(&gpsstarttime,':',' ')); 
      int mindiff=std::numeric_limits<int>::max();
      //Need to loop through all syncs to find the one closest to start time
      for(unsigned int vecindex=0;vecindex<navfile->GetNumSyncs();vecindex++)
      {
         if( abs(lev1starttime - navfile->GetGPSSync(vecindex)) < mindiff)
         {
            mindiff=abs(lev1starttime - navfile->GetGPSSync(vecindex));
            syncrecord_index=vecindex;
         }
      }    
      //Now have a second test here to make sure that the selected index is within X seconds 
      //of the level 1 start time. This is hard coded below!
      if(mindiff > lev1firstscanmaxexpectedsize)
      {
         //Are there multiple sync messages?
         if((navfile->GetNumSyncs()>1)&&(persec==false))
         {
            //Closest is still greater than "lev1firstscanmaxexpectedsize" seconds away - just list choices and let user select.
            Logger::Log("Multiple identical sync delay values in Specim nav file. None fall within the "+ToString(lev1firstscanmaxexpectedsize)+" seconds window of start time of level 1 file.");
            Logger::Log("Header start time: "+ToString(lev1starttime)+" sync delay value: "+ToString(hdrsync)+" ["+ToString(NOSYNCINHDR)+" means no header sync value]");
            for(unsigned int vecindex=0;vecindex<navfile->GetNumSyncs();vecindex++)   
            {
               //The value to add onto the hdr start time to get the "synced" time assuming the hdr start time is used as scan line 0 time
               double scntoff=(navfile->GetGPSSync(vecindex) -navfile->GetSyncDelay(vecindex)) - lev1starttime;
               Logger::Log("Sync message time: "+ToString(navfile->GetGPSSync(vecindex))+" Scantimeoffset value to use: "+ToString(scntoff));
            }
            throw "Multiple possibilities for sync time. Try using one of the suggested values above as a -scantimeoffset and processing with -nonav";
         }
         else if((navfile->GetNumSyncs()>1)&&(persec==true))
         {
            //This is to be expected (more than 1 sync in file)
            //Should be closer to start time than lev1firstscanmaxexpectedsize though?
            Logger::Warning("1 Sync message found but greater than 'maximum expected size' of "+ToString(lev1firstscanmaxexpectedsize)+" seconds away from level1 start time.");
         }
         else 
         {
            //Only 1 sync message but > lev1firstscanmaxexpectedsize seconds away - use with a warning
            Logger::Warning("1 Sync message found but greater than 'maximum expected size' of "+ToString(lev1firstscanmaxexpectedsize)+" seconds away from level1 start time.");
         }
      }

      //output a message for the sync index used
      if(persec==false)
      {
         Logger::Log("Using the sync message from index "+ToString(syncrecord_index)+" which probably means it is for flight line "+ToString(syncrecord_index)+" in this specim nav file (note these are referenced from 0 not 1).");
         Logger::Log("Sync value: "+ToString(navfile->GetSyncDelay(syncrecord_index)));
      }
      else
      {
         Logger::Log("Using the per second sync messages in this specim nav file.");
         Logger::Log("First sync delay value: "+ToString(navfile->GetSyncDelay(0)));
      }

      //Get the GPS time of the (first?) sync message
      double gps_integer_time=navfile->GetGPSSync(syncrecord_index);

      //Subtract the sync message time
      double firstscantime=gps_integer_time - navfile->GetSyncDelay(syncrecord_index);

      //Add on the cropped start time offset - this is 0 if no cropping (at the start of the line) has taken place
      Logger::Log("Applying crop time offset of: "+ToString(croptimeoffset));
      firstscantime=firstscantime + croptimeoffset;         

      //Compare to start time of level-1 header - in GPS sec of week - to prevent data from a different weekday being used 
      CompareLev1ToNavTimes(firstscantime);

      if(persec==false)
      {
         //Set the per scan times
         for(unsigned int i=0;i<nscans;i++)  
            time[i]=firstscantime + i*scanseparation;
      }
      else
      {
         unsigned int thissyncframe=navfile->GetFrame(syncrecord_index);
         unsigned int nextsyncframe=navfile->GetFrame(syncrecord_index+1);
         //Set the per scan times
         for(unsigned int i=0;i<nscans;i++)  
         {
            if((i<nextsyncframe)&&(i>=thissyncframe))
               time[i]=firstscantime + (i-thissyncframe)*scanseparation;
            else if(i<thissyncframe)
            {
               //This case shouldn't occur unless there is a problem with the nav file or our understanding of it
               //So for the moment we will not handle it other than to exit
               throw "(Software Bug): Frame number is less than current sync message frame - will need to write code for this eventuallity.";
            }
            else
            {
               //need to update to new sync message
               try
               {
                  syncrecord_index++;
                  if(syncrecord_index+1 < navfile->GetNumSyncs())
                  {
                     thissyncframe=navfile->GetFrame(syncrecord_index);
                     nextsyncframe=navfile->GetFrame(syncrecord_index+1);
                  }
                  else
                  {
                     //Here we don't actually have a next sync message - so all following scans shall use previous PPS to sync to
                     //adding on 1000 as an arbitrarily large number - this could essentially be anything > 1
                     nextsyncframe+=1000;
                     //Decrease the syncrecord index back to what it was
                     syncrecord_index--;
                  }
                  gps_integer_time=navfile->GetGPSSync(syncrecord_index);
                  firstscantime=gps_integer_time - navfile->GetSyncDelay(syncrecord_index);
                  firstscantime=firstscantime + croptimeoffset;         
               }
               catch(...)
               {
                  throw "Error in persecond syncing for each scantime - probably accessing past a bounded array. This is a coding bug - report it.";
               }
               time[i]=firstscantime + (i-thissyncframe)*scanseparation;              
            }               
         }
      }
      Logger::Log("First scan will have time: "+ToString(time[0])+".  Last scan will have time: "+ToString(time[nscans-1]));
      
   }
   else //Use a user-specified value as the first scan time
   {
      double firstscantime=0;

      //Use the start time from the header file as the first scan time
      firstscantime=GetSecOfWeek(acquisitiondate,ReplaceAllWith(&gpsstarttime,':',' '));
      Logger::Log("Using a first scan time of "+ToString(firstscantime));

      //Add on the cropped start time offset - this is 0 if no cropping (at the start of the line) has taken place
      Logger::Log("Applying crop time offset of: "+ToString(croptimeoffset));
      firstscantime=firstscantime + croptimeoffset;   

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
