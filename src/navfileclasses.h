//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef NAVFILECLASSES_H
#define NAVFILECLASSES_H

#include "datahandler.h"
#include "commonfunctions.h"
#include "logger.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <list>

#ifndef PI_
#define PI_
const double PI=4.0*atan(1.0);
#endif

//SBET file consists of a sequence of records each containing 17 doubles.

//-------------------------------------------------------------------------
// SBET handling class to read / handle the navigation from an SBET file
//-------------------------------------------------------------------------
class SBETData : public DataHandler
{
public:
   SBETData(const std::string filename);
   void Reader();

   static unsigned int GetRecordSize(){return sizeofrecord;}
private:
   static const unsigned int sizeofrecord=136; //17*8 bytes
};


//-------------------------------------------------------------------------
// Base class for SpecimNavData objects
//-------------------------------------------------------------------------
class SpecimNavData : public DataHandler
{
public:
   SpecimNavData();
   SpecimNavData(const std::string filename);
   virtual ~SpecimNavData();
   virtual void Reader(){throw "Unassigned reader function in SpecimNavData.";};
   double GetSyncDelay(const unsigned long i)
   {
      if(i<numsyncs)
         return syncdelay[i];
      else
         throw "Requested sync time index is out of bounds in GetSyncDelay().";
   };

   double GetGPSSync(const unsigned long i)
   {
      if(i<numsyncs)
         return syncgps[i];
      else
         throw "Requested GPS sync index is out of bounds in GetGPSSync(). Requested "+ToString(i)+" of "+ToString(numsyncs);
   };

   int GetFrame(const unsigned long i)
   {
      if(i<numsyncs)
         return syncframe[i];
      else
         throw "Requested frame index is out of bounds in GetFrame(). Requested "+ToString(i)+" of "+ToString(numsyncs);
   };

   unsigned long GetNumSyncs(){return numsyncs;}
   bool UsePerSecondForSync();
protected:
   
   unsigned long numsyncs; //number of sync messages
   double* syncdelay; //sync times from message
   int * syncgps; //integer GPS time of sync message
   int* syncframe; //frame number that delay corresponds too

   std::vector<int> persecond_frame;
   std::vector<int> persecond_syncgps;
   std::vector<double> persecond_syncdelay;
   bool use_persecond;

};


//-------------------------------------------------------------------------
// Class to hold SOL record headers and read in from sol file stream
// This is only used for easy reading of the records from the sol file
//  - data is actually stored in navdatarecords 
//-------------------------------------------------------------------------
class SOLRecordHeader
{
public:
   SOLRecordHeader(){};

   void Read(std::ifstream &fin);
   unsigned int GetSize(){return sizeofrecordheader;}

   unsigned char TimeType1(){return timetype1;}
   unsigned char TimeType2(){return timetype2;}

   double Time(int type)
   {
      if(type==1)
         {return time1;}
      else if(type==2)
         {return time2;}
      else {throw "Unrecognised time type ("+ToString(type)+") requested in SOLRecordHeader.";}
   }
private:
   unsigned char version,dataversion,sourceid,destinationid,status,reserved,timetype1,timetype2,checksum;
   unsigned short int preamble,messagelength,transactionid,messageid,gpsweek;
   double time1,time2;
   static const unsigned int sizeofrecordheader=35; //bytes 
};

//-------------------------------------------------------------------------
// Class to hold SOL data records and read in from sol file stream
// This is only used for easy reading of the records from the sol file
//  - data is actually stored in navdatarecords 
//-------------------------------------------------------------------------
class SOLRecordData
{
public:
   SOLRecordData(){};

   void Read(std::ifstream &fin);
   unsigned int GetSize(){return sizeofrecord;}
   double Height(){return hei;}
   double Lat(){return lat;}
   double Lon(){return lon;}
   double Roll(){return roll;}
   double Pitch(){return pitch;}
   double Heading(){return truehead;}

private:
   unsigned char solutionstatus,geoidmodel,dummy[512];
   unsigned short int datum,solutionorigin,solutionlevel;
   float stdev_lat,stdev_lon,stdev_hei,stdev_roll,stdev_pitch,stdev_head,stdev_north,stdev_east,stdev_up,geoidundulation;
   double lat,lon,hei,roll,pitch,truehead,northvel,eastvel,upvel,rollrate,pitchrate,headrate;
   static const unsigned int sizeofrecord=149;  //bytes
};

//-------------------------------------------------------------------------
// Class to hold a full SOL record [header + data] 
// This is only used for easy reading of the records from the sol file
//  - data is actually stored in navdatarecords 
//-------------------------------------------------------------------------
class SOLRecord
{
public:
   SOLRecord(){};
   SOLRecord(std::ifstream &fin){header.Read(fin); data.Read(fin);}

   SOLRecordHeader header;
   SOLRecordData data;

   unsigned int GetSize(){return (header.GetSize()+data.GetSize());}
   double Height(){return data.Height();}
   double Time()
   {
      //We want to return the GPS time from the SOL file - this has timetype value=1
      //check which of times 1 or 2 has this time type
      if(header.TimeType1() == 1) //GPS time
         return header.Time(1);
      else if(header.TimeType2() == 1)//GPS time
         return header.Time(2);
      else
         throw "Neither timetype in SOL Record header is GPS time (1).";
   }
   double Latitude(){return data.Lat();}
   double Longitude(){return data.Lon();}
   double Roll(){return data.Roll();}
   double Pitch(){return data.Pitch();}
   double Heading(){return data.Heading();}

private:
};

//-------------------------------------------------------------------------
// SOL handling class to read / handle the navigation from a SOL file
// For each record in the sol file there is an accompanying header
//-------------------------------------------------------------------------
class SOLData : public DataHandler
{
public:
   SOLData(const std::string filename);
   void Reader();
   
private:
   SOLRecord* record;
   unsigned int completerecordsize;
};


//-------------------------------------------------------------------------
//Abstract class to inherit from for different message types
//-------------------------------------------------------------------------
class Message
{
public:
   unsigned short wordcount, flags, headerchecksum, delayvalue, datachecksum,framecounter;
private:
};

//-------------------------------------------------------------------------
//describes the 999 message type - reads it in through a passed file stream
//-------------------------------------------------------------------------
class M999 : public Message
{
public:
   M999(std::ifstream &fin)
   {
      if(fin.is_open())
      {
         //Read in the #999 message
         fin.read((char*)&wordcount,sizeof(wordcount));
         fin.read((char*)&flags,sizeof(wordcount));
         fin.read((char*)&headerchecksum,sizeof(wordcount));
         fin.read((char*)&delayvalue,sizeof(wordcount));
         fin.read((char*)&datachecksum,sizeof(wordcount));
         framecounter=0;//this does not exist in 999 message types
      }
      else
      {
         throw "Cannot read #999 message as file stream is not open.";
      }
   }

};

//-------------------------------------------------------------------------
//describes the 998 message type - reads it in through a passed file stream
//-------------------------------------------------------------------------
class M998 : public Message
{
public:
   M998(std::ifstream &fin)
   {
      if(fin.is_open())
      {
         //Read in the #998 message
         fin.read((char*)&wordcount,sizeof(wordcount));
         fin.read((char*)&flags,sizeof(wordcount));
         fin.read((char*)&headerchecksum,sizeof(wordcount));
         fin.read((char*)&delayvalue,sizeof(wordcount));
         fin.read((char*)&framecounter,sizeof(wordcount));
         fin.read((char*)&datachecksum,sizeof(wordcount));
      }
      else
      {
         throw "Cannot read #998 message as file stream is not open.";
      }
   }
};

//-------------------------------------------------------------------------
// Sync message class
//-------------------------------------------------------------------------
class SyncMessage
{
public:
   SyncMessage(std::ifstream &fin)
   {
      message=NULL;
      //Read in next 2 bytes to check what kind of message this is
      if(fin.is_open())
      {
         fin.read((char*)&id,sizeof(id));
         if(id==999)
         {
            //Create a new #999 message
            message=new M999(fin);
         }
         else if(id==998)
         {
            //Create a new #998 message
            message=new M998(fin);
         }
         else
         {
            message=NULL;
            //std::cout<< "Unrecognised message type in SyncMessage. Got id: "<<id<<std::endl;
         }
      }
      else
      {
         message=NULL;
         throw "Cannot create SyncMessage as file stream is not open.";
      }
   };

   ~SyncMessage(){if(message != NULL) delete message;}

   //Return private members
   const Message* GetMessage() {return (const Message*)message;}
   unsigned short int ID() const {return id;}

//This function can not be here because it only applies for eagle data
//and the hawk data changes the flag too but is apparently meaningless.
//this will need to go into main code if/when we deal with 998 messages

//   //Function to get the frame counter value to apply the delay value too
//   unsigned short int GetFrameCounter()
//   {
//      if(id != 998)
//         return 0; //there is no frame counter for this message
//      else
//      {
//         if(message->flags == 32768)//just return the frame counter
//            return (message->framecounter);
//         else if(message->flags == 32769)//need to add one onto frame counter
//            return (message->framecounter+1);
//         else//something somewhere is wrong with this
//            return 0;
//      }
//   }

private:
   unsigned short id;
   Message* message;
};



//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//
// Specim NMEA style message stuff here
//
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Base class for NMEA style Specim nav files
//-------------------------------------------------------------------------
class MessageNMEA
{
public:
   MessageNMEA(){bad=false;}
   virtual ~MessageNMEA(){}

   std::string Id()const{return id;}
   bool Bad() const {return bad;}

protected:   
   virtual std::string StripChecksum(std::string str)
   {
      size_t index=str.find('*');
      if(index!=std::string::npos)
      {
         str=str.substr(0,index);
      }
      return str;
   }
   std::string id;  
   bool bad; 
};

//-------------------------------------------------------------------------
// Class to parse the GPZDA Date message 
// Sample message: $GPZDA,120003.150,19,01,2012,,,*76
// ID,Time,day,month,year,,,
//
// Note that the time is to be ignored as irrelevant in this message
//-------------------------------------------------------------------------
class GPZDA : public MessageNMEA
{
public:
   GPZDA(std::string message,char delim)
   {
      std::string delimstr=std::string(1,delim);
      id="$GPZDA";
      if(id.compare(GetItemFromString(message,0,delim))!=0)
         throw "Given message does not contain the GPZDA id tag in position 0.";

      //Test there are 8 objects separated by delim character in string
      if(GetNumberOfItemsFromString(message,delimstr)!=8)
      {
         bad=true;
         return;
      }

      day=StringToUINT(GetItemFromString(message,2,delim));
      month=StringToUINT(GetItemFromString(message,3,delim));
      year=StringToUINT(GetItemFromString(message,4,delim));     

      std::string date=GetItemFromString(message,2,delim)+"-"+GetItemFromString(message,3,delim)+"-"+GetItemFromString(message,4,delim);
      dayofweek=GetDayOfWeek(date);
      secofweek_to_startofday=dayofweek*3600*24;
   }
   unsigned int day,month,year;
   unsigned int dayofweek,secofweek_to_startofday;
};

//-------------------------------------------------------------------------
// Class to parse the PRDID Attitude message 
// Sample message: $PRDID,-1.20,+1.20,181.20*75
// ID,Pitch,Roll,Heading
//-------------------------------------------------------------------------
class PRDID : public MessageNMEA
{
public:
   PRDID(std::string message,char delim)
   {      
      std::string delimstr=std::string(1,delim);
      id="$PRDID";
      if(id.compare(GetItemFromString(message,0,delim))!=0)
         throw "Given message does not contain the PRDID id tag in position 0.";

      //Test there are 4 objects separated by delim character in string
      if(GetNumberOfItemsFromString(message,delimstr)!=4)
      {
         bad=true;
         return;
      }

      pitch=StringToDouble(GetItemFromString(message,1,delim),false);
      roll=StringToDouble(GetItemFromString(message,2,delim),false);
      heading=StringToDouble(StripChecksum(GetItemFromString(message,3,delim)),false);
   }
   double pitch,roll,heading;
};

//-------------------------------------------------------------------------
// Class to parse GPGGA Time and Position messages
// Sample message:
// $GPGGA,120003.000,6500.30000,N,2499.70000,E,1,13,3.5,1000.3,M,17.8,M,0.0,1023*7b
// ID,time (hhmmss),lat(*100),lon(*100),Fix quality,Num satellites,HDOP,ALT,Metres,HeightofGeoid,Metres,--,--
//-------------------------------------------------------------------------
class GPGGA : public MessageNMEA
{
public:
   GPGGA(std::string message,char delim)
   {      
      std::string delimstr=std::string(1,delim);
      id="$GPGGA";
      if(id.compare(GetItemFromString(message,0,delim))!=0)
         throw "Given message does not contain the GPGGA id tag in position 0.";

      //Test there are 15 objects separated by delim character in string
      if(GetNumberOfItemsFromString(message,delimstr)!=15)
      {
         bad=true;
         return;
      }

      //Get the TIME and convert into seconds of day
      std::string time=GetItemFromString(message,1,delim);
      int hh=StringToINT(time.substr(0,2));
      int mm=StringToINT(time.substr(2,2));
      int ss=StringToINT(time.substr(4,2));
      //double ss=StringToDouble(time.substr(4)); //As advised by Specim - just use the rounded down version of seconds - not full decimal
      secofday=hh*3600 + mm*60 +ss;
      //Get the Latitude and convert to degrees
      //Lat and Lon are stored in GPGGA as (d)ddmm.mmmm
      int latdeg=int(StringToDouble(GetItemFromString(message,2,delim),false)/100.0);
      double latmin=StringToDouble(GetItemFromString(message,2,delim),false) - latdeg*100 ;
      lat= latdeg + (latmin / 60.0);

      if(GetItemFromString(message,3,delim).compare("S")==0)
         lat=-lat;

      //Get the Longitude and convert to degrees
      int londeg=int(StringToDouble(GetItemFromString(message,4,delim),false)/100.0);
      double lonmin=StringToDouble(GetItemFromString(message,4,delim),false) - londeg*100 ;
      lon=londeg + (lonmin / 60.0);

      if(GetItemFromString(message,5,delim).compare("W")==0)
         lon=-lon;
      //Get the altitude and convert to above ellipsoid
      double height_mean_sea=StringToDouble(GetItemFromString(message,7,delim),false);
      double geoid_ellipsoid_sep=StringToDouble(GetItemFromString(message,9,delim),false);
      alt=height_mean_sea + geoid_ellipsoid_sep;

   }

   double secofday,lat,lon,alt;
};

//-------------------------------------------------------------------------
// Class to define start of recording (sync message)
// This is the equivalent of #999 messages
// Sample message: $SPTSMP2,500*32
//-------------------------------------------------------------------------
class SPTSMP2 : public MessageNMEA
{
public:
   SPTSMP2(std::string message,char delim)
   {      
      std::string delimstr=std::string(1,delim);
      id="$SPTSMP2";
      if(id.compare(GetItemFromString(message,0,delim))!=0)
         throw "Given message does not contain the SPTSMP2 id tag in position 0.";

      //Test there are 2 objects separated by delim character in string
      if(GetNumberOfItemsFromString(message,delimstr)!=2)
      {
         bad=true;
         return;
      }
      
      delayvalue=StringToDouble(StripChecksum(GetItemFromString(message,1,delim)),false)/1000.0;
   }
   double delayvalue;
};

//-------------------------------------------------------------------------
// Class to define periodic time stamps
// This is the equivalent of #998 messages
// Sample message: $SPTSMP,65,300,0*35
// ID,delay from latest 1pps pulse, current frame number, frame trigger pulse active flag
//-------------------------------------------------------------------------
class SPTSMP : public MessageNMEA
{
public:
   SPTSMP(std::string message,char delim)
   {      
      std::string delimstr=std::string(1,delim);
      id="$SPTSMP";
      if(id.compare(GetItemFromString(message,0,delim))!=0)
         throw "Given message does not contain the SPTSMP id tag in position 0.";

      //Test there are 4 objects separated by delim character in string
      if(GetNumberOfItemsFromString(message,delimstr)!=4)
      {
         bad=true;
         return;
      }
      //Test that there are only numbers in the string message components
      //If it throws an exception then skip these messages
      try
      {
         CheckNumbersOnly(GetItemFromString(message,1,delim));
         CheckNumbersOnly(GetItemFromString(message,2,delim));
      }
      catch(...)
      {
         throw "There appears to be a non-numeric value in a specim time stamp SPTSMP message in the raw .nav file. Please correct this and re-run.";
         //bad=true;
         //return;
      }

      delay=StringToDouble(GetItemFromString(message,1,delim),false)/10000.0;
      framenumber=StringToUINT(GetItemFromString(message,2,delim));
      triggerflag=StringToUINT(StripChecksum(GetItemFromString(message,3,delim)));
   }

   unsigned int framenumber,triggerflag;
   double delay;
};

//-------------------------------------------------------------------------
// Class to read in ascii NMEA-style Specim nav files
//-------------------------------------------------------------------------
class NMEASpecimNavData : public SpecimNavData
{
public:
   NMEASpecimNavData(std::string filename);//set dayofweeksecs to first zda record in file here
   virtual void Reader();

   unsigned int dayofweeksecs;
protected:
   unsigned long ScanNumSyncs();
   unsigned long GetNumRecords();
   unsigned int GetFirstDayOfWeekSecs();
   char delim;

};

//-------------------------------------------------------------------------
// Class to read in binary Specim nav files
//-------------------------------------------------------------------------
class BinSpecimNavData : public SpecimNavData
{
public:
   BinSpecimNavData(std::string filename);
   virtual void Reader();
protected:
   static const unsigned int sizeofrecord=28; //size of normal record
   static const unsigned int sizeofsync=14; //size of sync message

};

//-------------------------------------------------------------------------
// Class to decide what kind of nav file is being used
//-------------------------------------------------------------------------
class SpecimFileChooser 
{
public:
   SpecimFileChooser();
   SpecimFileChooser(std::string filename);
   ~SpecimFileChooser();

   void Reader(){spnav->Reader();}
   double GetSyncDelay(const unsigned long i){return spnav->GetSyncDelay(i);}
   double GetGPSSync(const unsigned long i){return spnav->GetGPSSync(i);}
   int GetFrame(const unsigned long i){return spnav->GetFrame(i);}
   unsigned long GetNumSyncs(){return spnav->GetNumSyncs();}
   unsigned long GetNumEntries(){return spnav->GetNumEntries();}
   NavDataLine* GetLine(const unsigned long l){return spnav->GetLine(l);}
   bool UsePerSecondForSync(){return spnav->UsePerSecondForSync();}
   bool IsASCII(){return asciifile;}

protected:
   SpecimNavData* spnav;
   bool asciifile;
};

#endif
