//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#include "navfileclasses.h"


//-------------------------------------------------------------------------
// Constructor for SBETData class
//-------------------------------------------------------------------------
SBETData::SBETData(const std::string filename)
{
   //Get the number of entries from the file somehow,
   //and set up the navdata array
   //Use file size and fact record = 17*8 =136 bytes
   unsigned long numentries=0;

   std::ifstream fin;
   fin.open(filename.c_str(),std::ios::binary);
   if(!fin.is_open())
   {
      throw "Cannot open SBET with filename "+filename;
   }
   else
   {
      fin.seekg(0,std::ios::end);
      unsigned long length=fin.tellg();
      fin.close();

      if(length % (sizeofrecord)!=0)
      {
         throw "SBET file may be corrupt - incomplete records assuming record size of 17 doubles";  
      }
      else
      {
         numentries=length / sizeofrecord;
      }
   }
   //Create array of nav data lines
   navcollection=new NavDataCollection(numentries);

   //Copy the filename over
   this->filename=filename;
}

//-------------------------------------------------------------------------
// Function to read in the SBET file and store it 
//-------------------------------------------------------------------------
void SBETData::Reader()
{
   //Read in the SBET file data and store in the navdataline array
   //SBET has 17 items in each record:
   // Time, Lat, Long, Alt, X vel, Y vel, Z vel, roll, pitch, platform heading,
   // wander angle, x body accel, y body accel, z body accel, x body angular rate,
   // y body angular rate, z body angular rate

   //We will only keep the Time, Lat, Lon, Alt, roll, pitch, heading for now

   //Test that numentries is not 0
   if(this->GetNumEntries() == 0)
      throw "Trying to read data into 0 sized arrays in SBETData::Reader()";

   std::ifstream fin;
   fin.open(filename.c_str(),std::ios::binary);
   if(!fin.is_open())
   {
      throw "Cannot open SBET "+filename+" in SBETData::Reader()";
   }
   else
   {
      //Create a buffer to read data into
      char buffer[sizeofrecord];
      //pointer to derefernce buffer as doubles
      double* dbuffer=(double*)buffer;
      for(unsigned long record=0;record<GetNumEntries();record++)
      {
         //read in a record into the buffer
         fin.read((char*)buffer,this->sizeofrecord);
         //Fill in the navdata array
         navcollection->SetValues(record,NavDataCollection::TIME,dbuffer[0]);
         navcollection->SetValues(record,NavDataCollection::LAT,dbuffer[1]*180/PI);
         navcollection->SetValues(record,NavDataCollection::LON,dbuffer[2]*180/PI);
         navcollection->SetValues(record,NavDataCollection::HEI,dbuffer[3]);
         navcollection->SetValues(record,NavDataCollection::ROLL,dbuffer[7]*180/PI);
         navcollection->SetValues(record,NavDataCollection::PITCH,dbuffer[8]*180/PI);
         if(dbuffer[9]<0)
            navcollection->SetValues(record,NavDataCollection::HEADING,(dbuffer[9]*180/PI+360-dbuffer[10]*180/PI)); //Heading plus 360 minus wander angle
         else
            navcollection->SetValues(record,NavDataCollection::HEADING,(dbuffer[9]*180/PI-dbuffer[10]*180/PI));//Heading minus wander angle


      }
      //Test file position
      unsigned long fileat=fin.tellg(); //get current pos
      fin.seekg(0,std::ios::end);
      unsigned long fileend=fin.tellg(); //get end of file pos
      if(fileat!=fileend)
         throw "SBET Reader has finished reading before the end of the file. Suggests numentries is wrong: "+ToString(this->GetNumEntries());
      
      //Close the file
      fin.close();
   }

   //Check the plausibility of the file
   CheckPlausibility();

   //Output some information on the file
   Logger::Log(GetInformation());
}

//-------------------------------------------------------------------------
// Specim NAV file methods
//-------------------------------------------------------------------------

SpecimNavData::SpecimNavData()
{
   this->synctime=NULL;
   this->syncgps=NULL;
   use_persecond=false;
}

//-------------------------------------------------------------------------
// SpecimNAVData Constructor
//-------------------------------------------------------------------------
SpecimNavData::SpecimNavData(const std::string filename)
{   
   //Copy the filename over
   this->filename=filename;
   this->synctime=NULL;
   this->syncgps=NULL;
   use_persecond=false;
}

//-------------------------------------------------------------------------
// SpecimNavData Destructor
//-------------------------------------------------------------------------
SpecimNavData::~SpecimNavData()
{
   if(synctime!=NULL)
      delete[] synctime;
   if(syncgps!=NULL)
      delete[] syncgps;
}

// Function to use the per second sync messages (SPTSMP) instead of the SPTSMP2 one
// Primarily for use when there are no SPTSMP messages in the file 
void SpecimNavData::UsePerSecondForSync()
{
   if(use_persecond==true)
   {
      //Delete the syncgps and times from the SPTSMP2 messages if there are any
      if(syncgps!=NULL)
         delete[] syncgps;
      if(synctime!=NULL)
         delete[] synctime;
      //Now create arrays of size 1 for the 1 sync message we will create
      numsyncs=1;
      syncgps=new int[1];
      synctime=new double[1];
      //Get the first SPTSMP message delay value - to make it like an SPTSMP2 one we want 1-delay
      synctime[0]=1-persecond_synctime.front();
      //Get the corresponding GPS time stamp for message - to make it same code for SPTSMP2 we want to plus 1
      syncgps[0]=persecond_syncgps.front() + 1;
      std::cout<<"sync time: "<<synctime[0]<<" sync_gps: "<<syncgps[0]<<std::endl;
   }
}


//-------------------------------------------------------------------------
// Constructor for SOL class
//-------------------------------------------------------------------------
SOLData::SOLData(const std::string filename)
{
   //Get the number of entries from the file somehow,
   //and set up the navdata array

   //Use file size and fact record = 149 + 35 bytes
   SOLRecord sol_record;
   completerecordsize=sol_record.GetSize();
   //int completerecordsize=sizeofheader+sizeofrecord;
   uint64_t numentries=0;

   std::ifstream fin;
   fin.open(filename.c_str(),std::ios::binary);
   if(!fin.is_open())
   {
      throw "Cannot open SOL file with filename "+filename;
   }
   else
   {
      fin.seekg(0,std::ios::end);
      uint64_t length=static_cast<uint64_t>(fin.tellg());
      fin.close();

      if(length % (completerecordsize)!=0)
      {
         throw "SOL file may be corrupt - incomplete records assuming record size of (header+record): "+ToString(completerecordsize);  
      }
      else
      {
         numentries=length / completerecordsize;
      }
   }
   //Create array of nav data lines
   navcollection=new NavDataCollection(numentries);

   //Copy the filename over
   this->filename=filename;
}


//-------------------------------------------------------------------------
// Function to read in the SOL file and store it 
//-------------------------------------------------------------------------
void SOLData::Reader()
{
   //Test that numentries is not 0
   if(this->GetNumEntries() == 0)
      throw "Trying to read data into 0 sized arrays in SOLData::Reader()";

   //int completerecordsize=sizeofheader+sizeofrecord;
   std::ifstream fin;
   fin.open(filename.c_str(),std::ios::binary);
   if(!fin.is_open())
   {
      throw "Cannot open SOL "+filename+" in SOLData::Reader()";
   }
   else
   {
   
      for(unsigned long recordid=0;recordid<GetNumEntries();recordid++)
      {
         //read in a record into the buffer
         this->record=new SOLRecord(fin);
         //Now only extract the parts we are interested in and throw away the rest
         
         navcollection->SetValues(recordid,NavDataCollection::TIME,this->record->Time());
         navcollection->SetValues(recordid,NavDataCollection::LAT,this->record->Latitude()*180/PI);
         navcollection->SetValues(recordid,NavDataCollection::LON,this->record->Longitude()*180/PI);
         navcollection->SetValues(recordid,NavDataCollection::HEI,this->record->Height());
         navcollection->SetValues(recordid,NavDataCollection::ROLL,this->record->Roll()*180/PI);
         navcollection->SetValues(recordid,NavDataCollection::PITCH,this->record->Pitch()*180/PI);
         if(this->record->Heading()<0)
            navcollection->SetValues(recordid,NavDataCollection::HEADING,this->record->Heading()*180/PI+360); //Heading plus 360 
         else
            navcollection->SetValues(recordid,NavDataCollection::HEADING,this->record->Heading()*180/PI);//Heading

         //std::cout<<recordid<<" "<<fin.tellg()<<" "<<this->record->Latitude()*180/PI<<" "<<this->record->Longitude()*180/PI<<" "<<this->record->Height()<<" "<<this->record->Roll()<<" "<<this->record->Pitch()<<" "<<this->record->Heading()<<std::endl;

         //destroy this record - no longer required
         delete this->record;
         this->record=NULL;      

      }

      //Test file position
      std::streampos fileat=fin.tellg(); //get current pos
      fin.seekg(0,std::ios::end);
      std::streampos fileend=fin.tellg(); //get end of file pos
      if(fileat!=fileend)
         throw "SOL Reader has finished reading before the end of the file. Suggests numentries is wrong: "+ToString(this->GetNumEntries());
      
      //Close the file
      fin.close();
   }

   //Check the plausibility of the file
   CheckPlausibility();

   //Output some information on the file
   Logger::Log(GetInformation());
}

//-------------------------------------------------------------------------
// Function to read in a SOL file record data from the current position 
// in the stream - no checking on whether it is a valid record
//-------------------------------------------------------------------------
void SOLRecordData::Read(std::ifstream &fin)
{
   //Check the stream is open
   if(!fin.is_open())
   {
      throw "Cannot read header from SOL file in SOLRecordData::Read()";
   }
   else
   {
      //Check there are enough bytes to read in
      std::streampos here=fin.tellg();
      fin.seekg(0,std::ios::end);
      std::streampos lengthtoendoffile=fin.tellg()-here;
      fin.seekg(here,std::ios::beg);
      if(lengthtoendoffile < sizeofrecord)
      {
         throw "Length of header is more than length to end of file in SOLRecordData::Read()";
      }
      //read in the record
      fin.read((char*)&datum,sizeof(datum));
      fin.read((char*)&solutionorigin,sizeof(solutionorigin));
      fin.read((char*)&solutionlevel,sizeof(solutionlevel));
      fin.read((char*)&solutionstatus,sizeof(solutionstatus));
      fin.read((char*)&lat,sizeof(lat));
      fin.read((char*)&lon,sizeof(lon));
      fin.read((char*)&hei,sizeof(hei));
      fin.read((char*)&stdev_lat,sizeof(stdev_lat));
      fin.read((char*)&stdev_lon,sizeof(stdev_lon));
      fin.read((char*)&stdev_hei,sizeof(stdev_hei));
      fin.read((char*)&roll,sizeof(roll));
      fin.read((char*)&pitch,sizeof(pitch));
      fin.read((char*)&truehead,sizeof(truehead));
      fin.read((char*)&stdev_roll,sizeof(stdev_roll));
      fin.read((char*)&stdev_pitch,sizeof(stdev_pitch));
      fin.read((char*)&stdev_head,sizeof(stdev_head));
      fin.read((char*)&northvel,sizeof(northvel));
      fin.read((char*)&eastvel,sizeof(eastvel));
      fin.read((char*)&upvel,sizeof(upvel));
      fin.read((char*)&stdev_north,sizeof(stdev_north));
      fin.read((char*)&stdev_east,sizeof(stdev_east));
      fin.read((char*)&stdev_up,sizeof(stdev_up));
      fin.read((char*)&rollrate,sizeof(rollrate));
      fin.read((char*)&pitchrate,sizeof(pitchrate));
      fin.read((char*)&headrate,sizeof(headrate));
      fin.read((char*)&geoidundulation,sizeof(geoidundulation));
      fin.read((char*)&geoidmodel,sizeof(geoidmodel));
      //5 further bytes here - at the moment these are ignored (checksum bytes?)
      fin.read((char*)dummy,5*sizeof(char));
   }
}

//-------------------------------------------------------------------------
// Function to read in a SOL file record header from the current position 
// in the stream - no checking on whether it is a valid record
//-------------------------------------------------------------------------
void SOLRecordHeader::Read(std::ifstream &fin)
{
   //Check the stream is open
   if(!fin.is_open())
   {
      throw "Cannot read header from SOL in SOLRecordHeader::Read()";
   }
   else
   {
      //Check there are enough bytes to read in
      std::streampos here=fin.tellg();
      fin.seekg(0,std::ios::end);
      std::streampos lengthtoendoffile=fin.tellg()-here;
      fin.seekg(here,std::ios::beg);
      if(lengthtoendoffile < sizeofrecordheader)
      {
         throw "Length of header is more than length to end of file in SOLRecordHeader::Read()";
      }
      //read in the header record
      fin.read((char*)&preamble,sizeof(preamble));
      fin.read((char*)&messagelength,sizeof(messagelength));
      fin.read((char*)&version,sizeof(version));
      fin.read((char*)&dataversion,sizeof(dataversion));
      fin.read((char*)&sourceid,sizeof(sourceid));
      fin.read((char*)&destinationid,sizeof(destinationid));
      fin.read((char*)&status,sizeof(status));
      fin.read((char*)&reserved,sizeof(reserved));
      fin.read((char*)&transactionid,sizeof(transactionid));
      fin.read((char*)&messageid,sizeof(messageid));
      fin.read((char*)&timetype1,sizeof(timetype1));
      fin.read((char*)&timetype2,sizeof(timetype2));
      fin.read((char*)&gpsweek,sizeof(gpsweek));
      fin.read((char*)&time1,sizeof(time1));
      fin.read((char*)&time2,sizeof(time2));
      fin.read((char*)&checksum,sizeof(checksum));
   }
}

//-------------------------------------------------------------------------
// Constructor for NMEA style specim nav file class
//-------------------------------------------------------------------------
NMEASpecimNavData::NMEASpecimNavData(std::string filename)
{
   delim=','; //hardcoded specim nav file delimiter here - should be easy enough to change later if needed
   //Copy the filename over
   this->filename=filename;
   numsyncs=ScanNumSyncs();
   unsigned long numentries=GetNumRecords();
   dayofweeksecs=GetFirstDayOfWeekSecs();

   navcollection=new NavDataCollection(numentries);
   if(numsyncs!=0)
   {
      synctime=new double[numsyncs];
      syncgps=new int[numsyncs];
   }
   else
   {
      //We must use the per-second messages if there are any - construct a sync from (one of) these
      Logger::Warning("There are no SPTSMP2 sync messages in this .nav file - will try to use SPTSMP instead.");
      use_persecond=true;
   }
}

//-------------------------------------------------------------------------
// Function to get the number of sync messages in file
//-------------------------------------------------------------------------
unsigned long NMEASpecimNavData::ScanNumSyncs()
{
   std::string line;
   unsigned long syncs=0;
   std::ifstream fin;
   fin.open(filename.c_str());
   if(!fin.is_open())
   {
      throw "Cannot open .nav "+filename+" in NMEASpecimNavData::Reader()";
   }
   else
   {   
      while(!fin.eof())
      {
         getline(fin,line);
         if(GetItemFromString(line,0,delim).compare("$SPTSMP2")==0)
         {
            SPTSMP2 s(line,delim);
            //Only use good messages
            if(s.Bad() == false)
               syncs++;
         }
      }
      fin.close();
      fin.clear();
   }
   return syncs;
}

//-------------------------------------------------------------------------
// Function to get the number of records in nav file (based on number of
// GPGGA records in the file)
//-------------------------------------------------------------------------
unsigned long NMEASpecimNavData::GetNumRecords()
{
   std::string line;
   unsigned long records=0;
   std::ifstream fin;
   fin.open(filename.c_str());
   if(!fin.is_open())
   {
      throw "Cannot open .nav "+filename+" in NMEASpecimNavData::Reader()";
   }
   else
   {   
      while(!fin.eof())
      {
         getline(fin,line);
         if(GetItemFromString(line,0,delim).compare("$GPGGA")==0)
         {
            GPGGA g(line,delim);
            //Only use good messages
            if(g.Bad() == false)
               records++;
         }
      }
      fin.close();
      fin.clear();
   }
   return records;
}

//-------------------------------------------------------------------------
// Function to return the number of seconds in the week to the start of
// this day (day of data in nav file)
//-------------------------------------------------------------------------
unsigned int NMEASpecimNavData::GetFirstDayOfWeekSecs()
{
   std::string line;
   unsigned int secs=0;
   std::ifstream fin;
   fin.open(filename.c_str());
   if(!fin.is_open())
   {
      throw "Cannot open .nav "+filename+" in NMEASpecimNavData::Reader()";
   }
   else
   {   
      while(!fin.eof())
      {
         getline(fin,line);
         if(GetItemFromString(line,0,delim).compare("$GPZDA")==0)
         {
            GPZDA g(line,delim);
            //If the data is bad then do nothing and get next record
            if(g.Bad() == false)
            {
               secs=g.secofweek_to_startofday;
               break; //we only want first occurence in file
            }
         }
      }
      fin.close();
      fin.clear();
   }
   return secs;
}


//-------------------------------------------------------------------------
// Function to read in records and sync messages 
//-------------------------------------------------------------------------
void NMEASpecimNavData::Reader()
{
   std::string message_key;
   std::string line;
   std::ifstream fin;

   bool record_complete=true;
   unsigned long record=0;
   unsigned long sync=0;

   GPGGA* gpgga=NULL;
   PRDID* prdid=NULL;
   GPZDA* gpzda=NULL;
   SPTSMP* sptsmp=NULL;
   SPTSMP2* sptsmp2=NULL;

   //Test that numentries is not 0
   if(this->GetNumEntries() == 0)
      throw "Trying to read data into 0 sized arrays in SpecimNavData::Reader()";

   fin.open(filename.c_str());
   if(!fin.is_open())
   {
      throw "Cannot open .nav "+filename+" in NMEASpecimNavData::Reader()";
   }
   else
   {
      while(!fin.eof())
      {
         //Get a line from the file - maybe getline is better in case of rogue spaces?
         fin>>line;
         //Get the key to the mssage type
         message_key=GetItemFromString(line,0,delim);

         if(message_key.compare("$GPGGA")==0)
         {
            //This is a position message - use this as a nav record
            if(record_complete==true) //ie not mid way through a record
            {
               gpgga=new GPGGA(line,delim);
               //Test if object is bad - most likely because it is incomplete
               if(gpgga->Bad()==false)
               {
                  navcollection->SetValues(record,NavDataCollection::TIME,gpgga->secofday+this->dayofweeksecs);
                  navcollection->SetValues(record,NavDataCollection::LAT,gpgga->lat);
                  navcollection->SetValues(record,NavDataCollection::LON,gpgga->lon);   
                  navcollection->SetValues(record,NavDataCollection::HEI,gpgga->alt);
                  record_complete=false;
               }

               delete gpgga;        
            }    
         }        
         else if(message_key.compare("$PRDID")==0)
         {
            //This is an attitude message - ONLY use this if it follows a GPGGA record
            //because there appears to be ~10x the amount of prdids than gpggas
            if(record_complete==false)   
            {
               prdid=new PRDID(line,delim);
               //Test if object is bad - most likely because it is incomplete
               if(prdid->Bad()==false)
               {
                  navcollection->SetValues(record,NavDataCollection::ROLL,prdid->roll);
                  navcollection->SetValues(record,NavDataCollection::PITCH,prdid->pitch);
                  navcollection->SetValues(record,NavDataCollection::HEADING,prdid->heading);    
                  record_complete=true;
                  record++;
               }        

               delete prdid;
            }
         }
         else if(message_key.compare("$GPZDA")==0)
         {
            //This is a date message - use this to get secofday into secofweek
            gpzda=new GPZDA(line,delim);
            //Test if object is bad - most likely because it is incomplete
            if(gpzda->Bad()==false)
            {
               if(dayofweeksecs!=gpzda->secofweek_to_startofday)
               {
                  //Day of week has changed mid way through file?
                  if(dayofweeksecs == gpzda->secofweek_to_startofday - 3600*24)
                  {
                     //Gone over to next day - this is plausible
                     //Update the dayofweeksecs
                     dayofweeksecs=gpzda->secofweek_to_startofday;
                  }
                  else
                     throw "Day of week has changed by more than 1 day (or gone backwards?) in nav file.";
               }
            }

            delete gpzda;
         }
         else if(message_key.compare("$SPTSMP")==0)
         {
            //This is a specim time stamp message (equivalent to #998)
            sptsmp=new SPTSMP(line,delim);
            //Test if object is bad - most likely because it is incomplete
            if(sptsmp->Bad()==false)
            {
               persecond_synctime.push_back(sptsmp->delay);
               if(record!=0)
                  persecond_syncgps.push_back((int)(navcollection->GetValue(record-1,NavDataCollection::TIME)));
               else
                  persecond_syncgps.push_back(-1); //Get at a later point when we have the next record - store as -1 (assume time is never -ve)
            }
            delete sptsmp;            
         }
         else if(message_key.compare("$SPTSMP2")==0)
         {
            //This is a specim sync message (equivalent to #999)
            //Get the sync time
            sptsmp2=new SPTSMP2(line,delim);
            //Test if object is bad - most likely because it is incomplete
            if(sptsmp2->Bad()==false)
            {
               synctime[sync]=sptsmp2->delayvalue;
               //Also store GPS time (integer part) 
               if(record!=0)
                  syncgps[sync]=(int)(navcollection->GetValue(record-1,NavDataCollection::TIME) +1);
               else
                  syncgps[sync]=-1; //Get at a later point when we have the next record - store as -1 (assume time is never -ve)

               sync++;
            }
            delete sptsmp2;
         }
         else
         {
            //Unrecoginsed message id
            Logger::Log("Unrecognised specim nav message ID: "+message_key+". Assuming corrupt record and trying again ...");
         }
      }
   }

   //Now check that -ve records do not exist in syncgps. If they do we need to replace them with the record value of TIME 
   //as it did not exist in the array at the time the sync message was found. Note that the record number will always be 0 as this
   //only occurrs if the sync message comes before the first record
   for(unsigned int ii=0; ii < sync;ii++)
   {  
      if(syncgps[ii] == -1)
      {
         syncgps[ii]=(int)(navcollection->GetValue(0,NavDataCollection::TIME));
      }
   }
   //And do the same for the persecond message arrays too
   for(std::list<int>::iterator it=persecond_syncgps.begin();it!=persecond_syncgps.end();it++)
   {
      if((*it)==-1)
      {
         (*it)=(int)(navcollection->GetValue(0,NavDataCollection::TIME));
      }
   }
   
   //Check the plausibility of the file
   CheckPlausibility();

   //Output some information on the file
   Logger::Log(GetInformation());
}

BinSpecimNavData::BinSpecimNavData(std::string filename)
{
   //Get the number of entries from the file somehow,
   //and set up the navdata array
   //Use file size and fact that normal  record = 28 bytes
   // sync record = 14 bytes
   char buffer[sizeofrecord];
   short* sbuff=(short*)buffer;
   unsigned long numentries=0;
   unsigned long nrecords=0;
   unsigned long nsyncs=0;
   bool corrupt=false;

   unsigned long corruptbytes=0;
   size_t flen=0;

   std::ifstream fin;
   fin.open(filename.c_str(),std::ios::binary);
   if(!fin.is_open())
   {
      throw "Cannot open Specim .nav file with filename "+filename;
   }
   else
   {
      //Get size of file
      fin.seekg(0,std::ios::end);
      flen=fin.tellg();

      fin.seekg(0,std::ios::beg);
      //Read in the first flag
      fin.read((char*)buffer,sizeof(short));

      //Need to read in record flags to find how many sync messages are in file
      while(!fin.eof())
      {
         if(sbuff[0] == -28160)
         {
            unsigned long pos=fin.tellg();
            fin.seekg(pos+(sizeofrecord-sizeof(short)),std::ios::beg); //move file pointer to start of next record
            if((!fin.good())||(fin.eof()))//failed to move the pointer along - premature end of file? break before increasing nrecords
            {
               std::cout<<"failed to move file pointer."<<std::endl;
               break;
            }
            fin.read((char*)buffer,sizeof(short));
            nrecords++;
         }
         else if (sbuff[0] == -32257)
         {
            unsigned long pos=fin.tellg();
            fin.seekg(pos+(sizeofsync-sizeof(short)),std::ios::beg); //move file pointer to start of next record
            if(fin.fail())//failed to move the pointer along - premature end of file? break before increasing nrecords
               break;
            fin.read((char*)buffer,sizeof(short));
            nsyncs++;
         }
         else
         {
            //corrupt file?
            unsigned long pos=fin.tellg();
            if(corrupt==false)
            {
               //We only want to output this message once
               corrupt=true;
               std::cout<< "Unrecognised flag in Specim nav file at "+ToString(pos)+" bytes, .nav file may be corrupt."<<std::endl;
            }

            //Added a while loop here and only search for next record (ignore any sync message since they will be no good anyway without the previous record)
            //This safe guards against getting a fluke -32257 message but could still get fluke -28160 which will cause a bad record
            while((sbuff[0]!=-28160)&&(fin.good()))
            {
               //Advance by 1 byte at a time until a recognisable flag? - how else can I deal with this.
               pos=fin.tellg();
               if(pos != 0)
                  fin.seekg(pos-1,std::ios::beg); //move file pointer position back 1 byte
               else
                  fin.seekg(1,std::ios::beg); //move file pointer position forward 1 byte
               fin.read((char*)buffer,sizeof(short));//read in 2 bytes

               //Count how many corrupt bytes there are
               corruptbytes+=1;
            }
         }
      }

      //std::cout<<"Size of file:"<<flen<<" Size from records etc:"<<nrecords*sizeofrecord+nsyncs*sizeofsync+corruptbytes<<std::endl;
      while(nrecords*sizeofrecord+nsyncs*sizeofsync+corruptbytes > flen)
      {
         nrecords=nrecords-1;
         //std::cout<<"Size of file:"<<flen<<" Size from records etc:"<<nrecords*sizeofrecord+nsyncs*sizeofsync+corruptbytes<<std::endl;
      }

      fin.close();
      numentries=nrecords;
      numsyncs=nsyncs;

   }
   //Copy the filename over
   this->filename=filename;

   //Create array of nav data lines
   navcollection=new NavDataCollection(numentries);
   synctime=new double[nsyncs];
   syncgps=new int[nsyncs];
}


void BinSpecimNavData::Reader()
{
   //Read in the nav file data and store in the navdataline array
   //.nav has 9 items in each 'normal' record:
   // flag, time, roll, pitch, heading, lat, lon, alt, vel
   // short double short short (unsigned) short, int, int, short, short 

   //We will only keep the Time, Lat, Lon, Alt, roll, pitch, heading 

   //Test that numentries is not 0
   if(this->GetNumEntries() == 0)
      throw "Trying to read data into 0 sized arrays in SpecimNavData::Reader()";

   std::ifstream fin;
   fin.open(filename.c_str(),std::ios::binary);
   if(!fin.is_open())
   {
      throw "Cannot open .nav "+filename+" in SpecimNavData::Reader()";
   }
   else
   {
      //Create a buffer to read data into
      char buffer[sizeofrecord];
      //create pointers to dereference buffer
      short* sbuff=(short*)buffer;
      unsigned short* usbuff=(unsigned short*)buffer;
      double* dbuff=(double*)buffer;
      int* ibuff=(int*)buffer;      
      //nuber of records read in
      unsigned long record=0;
      unsigned long sync=0;
      //Keep track of whether the file is suspected corrupt or not
      bool corrupt=false;

      fin.seekg(0,std::ios::beg);
      //Read in the first flag
      fin.read((char*)buffer,sizeof(short));
      //Need to read in record flags to find how many sync messages are in file
      while((!fin.eof())&&(record<navcollection->SizeOfArray()))
      {
         if(sbuff[0] == -28160)
         {
            fin.read(buffer,sizeof(double));
            if(fin.eof())  
            {
               std::cout<<"End of file reached prematurely (reading Time) - mid way through record "<<record<<" of "<<this->GetNumEntries()<<" records."<<std::endl;
               break;
            }
            navcollection->SetValues(record,NavDataCollection::TIME,dbuff[0]);
            fin.read(buffer,sizeof(short)*3);
            if(fin.eof())  
            {
               std::cout<<"End of file reached prematurely (reading Roll/Pitch/Heading) - mid way through record "<<record<<" of "<<this->GetNumEntries()<<" records."<<std::endl;
               break;
            }
            navcollection->SetValues(record,NavDataCollection::ROLL,(double)(sbuff[0]/100.0));
            navcollection->SetValues(record,NavDataCollection::PITCH,(double)(sbuff[1]/100.0));
            navcollection->SetValues(record,NavDataCollection::HEADING,(double)(usbuff[2]/100.0));
            fin.read(buffer,sizeof(int)*2);
            if(fin.eof())  
            {
               std::cout<<"End of file reached prematurely (reading Lat/Lon) - mid way through record "<<record<<" of "<<this->GetNumEntries()<<" records."<<std::endl;
               break;
            }
            navcollection->SetValues(record,NavDataCollection::LAT,(double)(ibuff[0]/3600000.0));
            navcollection->SetValues(record,NavDataCollection::LON,(double)(ibuff[1]/3600000.0));            
            fin.read(buffer,sizeof(short)*2); //read in 2 shorts - 2nd is not needed
            if(fin.eof())  
            {
               std::cout<<"End of file reached prematurely (reading Height) - mid way through record "<<record<<" of "<<this->GetNumEntries()<<" records."<<std::endl;
               break;
            }

            navcollection->SetValues(record,NavDataCollection::HEI,(double)(usbuff[0]/10.0));

            //increase record counter
            record++;
            //read in next record flag
            fin.read(buffer,sizeof(short)); 

         }
         else if (sbuff[0] == -32257)
         {
            //A sync message has been found - read it in
            SyncMessage syncm(fin);
            //message 999 occurs at start of the flight line
            if(syncm.ID()==999)
            {
               //Store the delay value (converted to seconds) 
               synctime[sync]=(syncm.GetMessage()->delayvalue/1000.0);
               //Also store GPS time (integer part) 
               syncgps[sync]=(int)(navcollection->GetValue(record-1,NavDataCollection::TIME) +1);
               //increase sync counter
               sync++;
               //read in next record flag
               fin.read(buffer,sizeof(short));
            }   
            else if(syncm.ID()==998)
            {
               //message 998 occurs every second of flight line - not implemented here yet
               //readin the next flag and continue
               fin.read(buffer,sizeof(short));

// need to read in the 998 message here and store the synctime and syngps


            }
            else
            {
               //Unrecognised message
               std::cout<< "Unrecognised message type in SyncMessage. Got id: "<<syncm.ID()<<std::endl;
               //readin the next flag and continue
               fin.read(buffer,sizeof(short));
            }
         }
         else
         {
            //corrupt file?
            unsigned long pos=fin.tellg();
            if(corrupt==false)
            {
               //We only want to output this message once
               corrupt=true;
               std::cout<< "Unrecognised flag in Specim nav file at "+ToString(pos)+" bytes, .nav file may be corrupt."<<std::endl;
            }

            //Added a while loop here and only search for next record (ignore any sync message since they will be no good anyway without the previous record)
            //This safe guards against getting a fluke -32257 message but could still get fluke -28160 which will cause a bad record
            while(sbuff[0]!=-28160)
            {
               //Advance by 1 byte at a time until a recognisable flag? - how else can I deal with this.
               pos=fin.tellg();
               if(pos != 0)
                  fin.seekg(pos-1,std::ios::beg); //move file pointer position back 1 byte
               else
                  fin.seekg(1,std::ios::beg); //move file pointer position forward 1 byte
               fin.read((char*)buffer,sizeof(short));//read in 2 bytes
               //std::cout<<sbuff[0]<<std::endl;
               if(sbuff[0]==-32257)
               {
                  //Found a -32257 message BUT this occurs after a load of corrupt bytes.
                  //Do we trust it? - Say no for now and output a message to log the fact a "possible" sync message found
                  //If we decide to trust it future need to change "999" section above so that record-1 is not used to get times
                  //as the actual record - 1 was corrupt and ignored 
                  SyncMessage syncm(fin);
                  if(syncm.ID()==999)
                  {
                     double cordelay=(syncm.GetMessage()->delayvalue/1000.0);
                     int corsec=(int)(navcollection->GetValue(record-1,NavDataCollection::TIME) +1);
                     Logger::Warning("Possible sync message found but in stream of corrupt bytes - not trusted." 
                     "\nSync delay value: = "+ToString(cordelay)+
                     "\nGPS Second value: = "+ToString(corsec)+
                     "\nStart time of instrument data capture: = "+ToString(corsec-cordelay)+
                     "\nOffset to apply to scantimeoffset (excluding the usual timing error) if not using nav file: = "+ToString(corsec-cordelay)+" - lev1 data start time");
                  }   
               }
            }
         }
      }

      //Test file position
      unsigned long fileat=fin.tellg(); //get current pos
      fin.seekg(0,std::ios::end);
      unsigned long fileend=fin.tellg(); //get end of file pos
      if(fileat!=fileend)
         Logger::Log("SpecimNav Reader has finished reading before the end of the file. Suggests numentries is wrong or nav file contains corrupt records: "+ToString(this->GetNumEntries()));
      
      //Close the file
      fin.close();
   }

   //Check the plausibility of the file
   CheckPlausibility();

   //Output some information on the file
   Logger::Log(GetInformation());
}


SpecimFileChooser::SpecimFileChooser(std::string filename)
{
   spnav=NULL;
   asciifile=false;
   std::string line;
   std::ifstream fin;
   fin.open(filename.c_str());
   if(!fin.is_open())
   {
      throw "Specim nav file failed to open: "+filename;
   }
   else
   {
      while(!fin.eof())
      {
         std::getline(fin,line);
         if(GetItemFromString(line,0,',').compare("$GPGGA")==0)
         {
            //File is most likely ascii
            asciifile=true;
            break;
         }
      }
      fin.close();
      fin.clear();
   }

   if(asciifile==true)
   {
      //Ascii nav file
      spnav=new NMEASpecimNavData(filename);
   }
   else
   {
      //Binary nav file
      spnav=new BinSpecimNavData(filename);
   }
}

SpecimFileChooser::~SpecimFileChooser()
{
   if(spnav!=NULL)
      delete spnav;
}

