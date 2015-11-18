//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#include "navigationinterpolator.h"

#ifdef DEBUGNAVINTERP
   #define DEBUGPRINT(x) std::cout<<x<<std::endl;
#else
   #define DEBUGPRINT(x)
#endif

//-------------------------------------------------------------------------
// Default NavigationInterpolator constructor
//-------------------------------------------------------------------------
NavigationInterpolator::NavigationInterpolator()
{
   scanid=NULL;
   navcollection=NULL;
   nscans=0;
   dhandle=NULL;
}

//-------------------------------------------------------------------------
// NavigationInterpolator constructor using nav file and lev1 file to set up
//-------------------------------------------------------------------------
NavigationInterpolator::NavigationInterpolator(std::string navfilename, std::string lev1filename)
{
   //Set to NULL here anyway just to be safe - they should all be none-null by the end of this function
   scanid=NULL;
   navcollection=NULL;
   nscans=0;
   dhandle=NULL;
   //hdrsync=0;

   //Call function to detect file type
   FILETYPE ftype=DetectFileType(navfilename);
 
   //Use the navfilename to auto-detect the nav data file type
   //std::string::reverse_iterator riter=navfilename.rbegin();
   //std::string postfix;
   //Get the last 4 chars of filename
   //for(unsigned int i=0;i<4;i++)
   //   postfix=*(riter++)+postfix;
   
   //if(postfix.compare(".nav")==0)
   if(ftype==SPECIMNAV)
   {
      //Nav file is a specim .nav file
      SpecimFileChooser spf(navfilename);
      if(spf.IsASCII()==true)
         dhandle=new NMEASpecimNavData(navfilename);
      else
         dhandle=new BinSpecimNavData(navfilename);
   }
   //else if(postfix.compare(".out")==0)
   else if(ftype==SBET)
   {
      //Nav file is a SBET file
      dhandle=new SBETData(navfilename);
   }
   else if(ftype==SOL)
   {
      dhandle=new SOLData(navfilename);
   }
   else if(ftype==BADFILE)
   {
      throw "Problem finding Navigation file - are you sure it exists?: "+navfilename;
   }
   else
   {
      //Unrecognised file format
      throw "Navigation file is of an unrecognised format (in NavigationInterpolator()): "+navfilename;
   }

   //Need to read in the level1 file hdr info to get nscans
   BinFile bilin(lev1filename);
   //Get the number of scan lines
   nscans=StringToUINT(bilin.FromHeader("lines"));
   DEBUGPRINT("Number of scans: "<<nscans)
   //Get the start and stop times from the header to use in case no sync messages 
   //found in the specim nav file
   gpsstarttime=bilin.FromHeader("GPS Start Time");
   gpsstarttime=RemoveAllBut(gpsstarttime,"1234567890.:");
   gpsstarttime=TrimWhitespace(ReplaceAllWith(&gpsstarttime,':',' '));
   gpsstoptime=bilin.FromHeader("GPS Stop Time");
   gpsstoptime=RemoveAllBut(gpsstoptime,"1234567890.:");
   gpsstoptime=TrimWhitespace(ReplaceAllWith(&gpsstoptime,':',' '));
   DEBUGPRINT("GPS Start/Stop times:"<<gpsstarttime<<" "<<gpsstoptime)
   bilin.Close();

   //Setup the data arrays
   SetupArrays();

   //read in the data
   dhandle->Reader();

}

//-------------------------------------------------------------------------
// NavigationInterpolator destructor
//-------------------------------------------------------------------------
NavigationInterpolator::~NavigationInterpolator()
{
   if(scanid!=NULL)
      delete[] scanid;
   if(navcollection!=NULL);
      delete navcollection;
   if(dhandle!=NULL)
      delete dhandle;
}

//-------------------------------------------------------------------------
// Function to set up the data arrays 
//-------------------------------------------------------------------------
void NavigationInterpolator::SetupArrays()
{
   //Check if the number of scans has been initialised
   if(nscans==0)
   {
      throw "Error - Cannot set up data arrays in NavigationInterpolator::SetupArrays(), number of scan lines is zero.";
   }
   
   //Check that arrays are not already initialised
   if((scanid !=NULL)||(navcollection!=NULL))
   {
      throw "Error - The scanid array or navcollection object are already initialised in NavigationInterpolator::SetupArrays().";
   }

   //Set up the arrays
   scanid=new double[nscans];
   navcollection=new NavDataCollection(nscans);
}

//-------------------------------------------------------------------------
// Function to set the navdata scan times to the given array values
//-------------------------------------------------------------------------
void NavigationInterpolator::SetTimes(const double* const times)
{
   //Assuming the length of times and the navdata array are the same - they should be
   for(unsigned int i=0;i<nscans;i++)
   {
      navcollection->SetValues(i,NavDataCollection::TIME,times[i]);
   }
}


//-------------------------------------------------------------------------
// Function to interpolate the navigation data to the scan times using the given function
//-------------------------------------------------------------------------

void NavigationInterpolator::Interpolate(void (*f)(const double*,const int,DataHandler*,NavDataCollection*,std::string,std::string))
{
   //Create an array of times
   double* times=new double[this->nscans];
   for(unsigned int i=0;i<this->nscans;i++)
   {
      times[i]=this->navcollection->GetValue(i,NavDataCollection::TIME);
      //Added a test to see if all interpolated times are within nav data 
      if((times[i] > dhandle->GetLine(dhandle->GetNumEntries()-1)->time) || (times[i] < dhandle->GetLine(0)->time))
      {
         throw "Interpolated time is outside the range of the navigation data for scan line: "+ToString(i)+" and time:"+ToString(times[i]);
      }      
   }

   //Call the function
   f(times,this->nscans,dhandle,navcollection,this->gpsstarttime,this->gpsstoptime);

   //clean up after ourselves.
   delete[] times;
}




//-------------------------------------------------------------------------
// Function to write out the data to a BIL file
//-------------------------------------------------------------------------
void NavigationInterpolator::Writer(std::string outfilename,std::string extrainfo)
{
   //Write out the data as doubles to a 7 band, 1 sample, N lines BIL file
   BILWriter bilout(outfilename,FileWriter::float64,nscans,1,7,'w');
   if (extrainfo.compare("")!=0)
      bilout.AddToHdr(extrainfo);

   bilout.AddToHdr("band names = {Time, Latitude, Longitude, Altitude, Roll, Pitch, Heading}");

   for(unsigned int i=0;i<nscans;i++)
   {
      //Set the height values to 2 dp only
      navcollection->SetValues(i,NavDataCollection::HEI,((int)(navcollection->GetValue(i,NavDataCollection::HEI)*100))/100.0);

      //Convert times to GPS second of day from second of week
      navcollection->SetValues(i,NavDataCollection::TIME,((int)navcollection->GetValue(i,NavDataCollection::TIME) % (24*3600)) 
                                                       + navcollection->GetValue(i,NavDataCollection::TIME)
                                                       - (int)navcollection->GetValue(i,NavDataCollection::TIME) );

      //Write out per band to ensure correct order maintained
      bilout.WriteBandLine((char*)navcollection->GetReferenceValue(i,NavDataCollection::TIME));
      bilout.WriteBandLine((char*)navcollection->GetReferenceValue(i,NavDataCollection::LAT));
      bilout.WriteBandLine((char*)navcollection->GetReferenceValue(i,NavDataCollection::LON));
      bilout.WriteBandLine((char*)navcollection->GetReferenceValue(i,NavDataCollection::HEI));
      bilout.WriteBandLine((char*)navcollection->GetReferenceValue(i,NavDataCollection::ROLL));
      bilout.WriteBandLine((char*)navcollection->GetReferenceValue(i,NavDataCollection::PITCH));
      bilout.WriteBandLine((char*)navcollection->GetReferenceValue(i,NavDataCollection::HEADING));
   } 
   //Close the bil writer
   bilout.Close();
}

//-------------------------------------------------------------------------
// Function to write out the navigation quality flag data to a BIL file
//-------------------------------------------------------------------------
void NavigationInterpolator::WriteFlags(std::string outfilename,std::string extrainfo)
{
   BILWriter bilout(outfilename,FileWriter::uchar8,nscans,1,1,'w');
   if (extrainfo.compare("")!=0)
      bilout.AddToHdr(extrainfo);

   bilout.AddToHdr("band names = {Quality Flag}");
   //Add the bad flag information to the header file
   std::stringstream stream(std::stringstream::in | std::stringstream::out);
   stream<<"flags = {Good data = 0, Bad Latitude = "<<NavDataLine::BADLAT<<", Bad Longitude = "<<NavDataLine::BADLON<<", Bad Height = "<<NavDataLine::BADHEI
         <<", Bad Roll = "<<NavDataLine::BADROLL<<", Bad Pitch = "<<NavDataLine::BADPITCH<<", Bad Heading = "<<NavDataLine::BADHEADING
         <<", Bad Time = "<<NavDataLine::BADTIME<<"}";
   bilout.AddToHdr(stream.str());
   //Write out the data
   for(unsigned int i=0;i<nscans;i++)
   {
      bilout.WriteBandLine((char*)navcollection->GetFlagReference(i));      
   }
}

//-------------------------------------------------------------------------
// Function to add on the boresight offset onto the per scan data
//-------------------------------------------------------------------------
void NavigationInterpolator::ApplyBoresight(Boresight* boresight)
{

   for(unsigned int i=0;i<nscans;i++)
   {
      boresight->ApplyBoresight(navcollection->GetReferenceValue(i,NavDataCollection::ROLL),
                                navcollection->GetReferenceValue(i,NavDataCollection::PITCH),
                                navcollection->GetReferenceValue(i,NavDataCollection::HEADING));
   }
}

//-------------------------------------------------------------------------
// Function to add on the lever arm offsets to the per scan position data
//-------------------------------------------------------------------------
void NavigationInterpolator::ApplyLeverarm(Leverarm* leverarm)
{
   //For each epoch (a scan line of data) we need to apply the lever arm
   //onto the GPS position, the lever arm will be different for each epoch
   //depending on the roll, pitch, heading of the aircraft

   for(unsigned int i=0;i<nscans;i++)
   {
      leverarm->ApplyLeverArm(navcollection->GetValue(i,NavDataCollection::ROLL),navcollection->GetValue(i,NavDataCollection::PITCH),
                              navcollection->GetValue(i,NavDataCollection::HEADING),
                              navcollection->GetReferenceValue(i,NavDataCollection::LAT),
                              navcollection->GetReferenceValue(i,NavDataCollection::LON),
                              navcollection->GetReferenceValue(i,NavDataCollection::HEI));
   }
}

//-------------------------------------------------------------------------
// Function to apply a shift between the position and attitude data
// that keeps the position the same but moves the attitude data
//-------------------------------------------------------------------------
void NavigationInterpolator::PosAttShift(void (*f)(const double*,const int,DataHandler*,NavDataCollection*,std::string,std::string),const double toffset)
{
   //Create a temporary navdataline object to hold the shifted data
   NavDataCollection tmpnavdata(this->nscans);

   //Add on the position-attitude time shift to the scan times
   double* times=new double[this->nscans];
   for(unsigned int i=0;i<this->nscans;i++)
      times[i]=this->navcollection->GetValue(i,NavDataCollection::TIME)+toffset;

   //Call the interpolation function, storing the result in tmpnavdata
   f(times,this->nscans,this->dhandle,&tmpnavdata,this->gpsstarttime,this->gpsstoptime);

   //Now copy the attitude data from the tmp to "proper" storage cells
   for(unsigned int i=0;i<this->nscans;i++)
   {
      this->navcollection->SetValues(i,NavDataCollection::ROLL,tmpnavdata.GetValue(i,NavDataCollection::ROLL));
      this->navcollection->SetValues(i,NavDataCollection::PITCH,tmpnavdata.GetValue(i,NavDataCollection::PITCH));
      this->navcollection->SetValues(i,NavDataCollection::HEADING,tmpnavdata.GetValue(i,NavDataCollection::HEADING));
   }

   //Free up the memory
   delete[] times;
}

//-------------------------------------------------------------------------
// Function to check the navcollection for possible bad nav data
//-------------------------------------------------------------------------
void NavigationInterpolator::CheckPlausibility()
{
   navcollection->CheckPlausibility();
}

//-------------------------------------------------------------------------
// Function to detect the file type of filename [SBET or SPECIM NAV]
//-------------------------------------------------------------------------
FILETYPE NavigationInterpolator::DetectFileType(std::string filename)
{
   FILETYPE ftype;
   std::ifstream fin;
   //Check file exists and if so open it
   if(DoesPathExist(filename))
      fin.open(filename.c_str());
   else 
      return BADFILE; //return bad file flag
   
   //SBET file should have x records of size 17*double
   //get the size of the file
   fin.seekg(0,std::ios::end);
   unsigned long int filesize=fin.tellg();
   //Check against sbet record size if the remainder of division is 0
   if(fmod(filesize,SBETData::GetRecordSize())==0)
   {
      //Likely an SBET file (not definite)

      //Also check against SOL here as detection method is not infalible
      if(fmod(filesize,SOLRecord().GetSize())==0)
      {
         //Could be either SOL or SBET - spew out warning and default to file extension detection
         Logger::Warning("Unable to detect whether navigation file is sol or sbet: "+filename);
         Logger::Log("Using file extension instead to detect instead...");
         if(filename.find(".sol")!=std::string::npos)
         {
            Logger::Log("Filename suggests sol file");
            ftype=SOL;
         }
         else if(filename.find(".out")!=std::string::npos)
         {
            Logger::Log("Filename suggests SBET file");
            ftype=SBET;
         }
         else
         {
            Logger::Log("Unable to detect file type from filename.");
            return BADFILE; //return bad file flag
         }   
      }
      else
      {
         Logger::Log("Detected this file as an SBET as it is divisible by sbet record size: "+filename);
         ftype=SBET;
      }

   }
   else if(fmod(filesize,SOLRecord().GetSize())==0)
   {
      //Likely a SOL file (not definite)
      Logger::Log("Detected this file as a SOL as it is divisible by sol record size: "+filename);
      ftype=SOL;
   }
   else
   {
      //Likely a SPECIM NAV file (not definite)
      Logger::Log("Detected this file as Specim nav as it failed sbet and sol tests: "+filename);
      ftype=SPECIMNAV;      
   }

   //Close the file
   fin.close();
   fin.clear();
   //Return the file type
   return ftype;   
}
