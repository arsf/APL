//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cerrno>
#include <cmath>
#include <vector>

#include "specimsensors.h"

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
   case FENIX:
      for(std::vector<int>::const_iterator it=FENIX_SENSOR_IDS.begin();it != FENIX_SENSOR_IDS.end();it++)
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

//-------------------------------------------------------------------------
// Specim Binary File (BinFile) constructor
//-------------------------------------------------------------------------
SpecimBinFile::SpecimBinFile(std::string strFilename) : BinFile(strFilename)
{
   sensorid=TestFile();
}

//-------------------------------------------------------------------------
// Specim Binary File (BinFile) destructor
//-------------------------------------------------------------------------
SpecimBinFile::~SpecimBinFile()
{
   if(br!=NULL)
      delete br;
   br=NULL;
}


//-------------------------------------------------------------------------
// Function to test the suitability of the file
//-------------------------------------------------------------------------
std::string SpecimBinFile::TestFile()
{
   //Test if unsigned 16-bit integer
   if(br->GetDataSize()!=2)  
      throw BinaryReader::BRexception("Cannot read from Specim File as data size was expected to be 2 bytes");   

   //Test frame rate for sanity
   double fps=StringToDouble(br->FromHeader("fps","true"));
   //Test it is greater than 0 and less than a suitably large number
   if((fps <= 0)||(fps > 100000))
      Logger::Warning("Frame rate (fps) in hdr file seems to be an erroneous value. Continuing as fps is not required for calibration.");

   //Insert other checks for existing header items here too

   //return the sensorid to identify if eagle, hawk or fenix file
   return br->FromHeader("sensorid");
}

//-------------------------------------------------------------------------
// EagleHawkBinFile Constructor
//-------------------------------------------------------------------------
EagleHawkBinFile::EagleHawkBinFile(std::string strFilename) : SpecimBinFile(strFilename)
{
   //Test 'tint' exists
   br->FromHeader("tint","true");
}


//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
//void EagleHawkBinFile::SensorSpecificTestFile()
//{
//   //Insert other checks for existing header items here too  
//}

//-------------------------------------------------------------------------
// Function to get particular items from the file such that the calling
// function does not need to use the FromHeader function. Main reason
// for this function existing is to overload it for Fenix files
//-------------------------------------------------------------------------
std::string EagleHawkBinFile::GetFromFile(std::string keyword)
{
   if(keyword.compare("spatialbinning")==0)
   {
      return FromHeader("binning",1,"true");
   }
   else if(keyword.compare("spectralbinning")==0)
   {
      return FromHeader("binning",0,"true");
   }
   else if(keyword.compare("binningForHeader")==0)
   {
      std::string strspec=FromHeader("binning",0,"true");
      std::string strspat=FromHeader("binning",1,"true");
      return "binning = {"+strspec+","+strspat+"}";
   }
   else if(keyword.compare("lowerhimg")==0)
   {
      return FromHeader("himg",0,"true");
   }
   else if(keyword.compare("upperhimg")==0)
   {
      return FromHeader("himg",1,"true");
   }
   else if(keyword.compare("lowervimg")==0)
   {
      return FromHeader("vimg",0,"true");
   }
   else if(keyword.compare("uppervimg")==0)
   {
      return FromHeader("vimg",1,"true");
   }
   else if(keyword.compare("integrationtime")==0)
   {
      return FromHeader("tint","true");
   }
   else if(keyword.compare("tintForHeader")==0)
   {
      return "tint = "+FromHeader("tint","true");
   }
   else if(keyword.compare("subsensorBandsForHeader")==0)
   {
      return "bandrange = {"+GetFromFile("lowervimg")+","+GetFromFile("uppervimg")+"}";
   }
   else if(keyword.compare("Wavelength")==0)
   {
      std::string wlv="";
      for(unsigned int b=0;b<NumBands();b++)
      {
         wlv+=FromHeader("Wavelength",b) +";";
      }
      return TrimPunctuation(wlv);
   }
   else if(keyword.compare("fwhm")==0)
   {
      std::string fwhm="";
      for(unsigned int b=0;b<NumBands();b++)
      {
         fwhm+=FromHeader("fwhm",b) +";";
      }
      return TrimPunctuation(fwhm);
   }
   else
      return FromHeader(keyword);
}


//-------------------------------------------------------------------------
// FenixBinFile Constructor - defaults to the VNIR section of the file
// to change it to SWIR call the SetSubSensor function
//-------------------------------------------------------------------------
FenixBinFile::FenixBinFile(std::string strFilename) : SpecimBinFile(strFilename)
{
   //Default to VNIR
   subsensor=VNIR;
   subsenlowerband=StringToUINT(GetFromFile("lowervimg"))-1;//minus 1 to normalize to 0
   subsenupperband=StringToUINT(GetFromFile("uppervimg"))-1;//minus 1 to normalize to 0
   subsennumberofbands=subsenupperband-subsenlowerband+1;
   //Test the two 'tint' values exist
   br->FromHeader("tint1","true");
   br->FromHeader("tint2","true");
}

//-------------------------------------------------------------------------
// FenixBinFile destructor
//-------------------------------------------------------------------------
FenixBinFile::~FenixBinFile()
{
}

//-------------------------------------------------------------------------
// 
//-------------------------------------------------------------------------
//void FenixBinFile::SensorSpecificTestFile()
//{
//   //Insert other checks for existing header items here too  
//}

//-------------------------------------------------------------------------
// SetSubSensor function used to switch between VNIR and SWIR reading
//-------------------------------------------------------------------------
void FenixBinFile::SetSubSensor(const Subsensor sub)
{
   subsensor=sub;
   subsenlowerband=StringToUINT(GetFromFile("lowervimg"))-1;//minus 1 to normalize to 0
   subsenupperband=StringToUINT(GetFromFile("uppervimg"))-1;//minus 1 to normalize to 0
   subsennumberofbands=subsenupperband-subsenlowerband+1;
}

//-------------------------------------------------------------------------
// SetSubSensor function used to switch between VNIR and SWIR reading
//-------------------------------------------------------------------------
void FenixBinFile::SetSubSensor(const unsigned int sub)
{
   subsensor=static_cast<Subsensor>(sub);
   subsenlowerband=StringToUINT(GetFromFile("lowervimg"))-1;//minus 1 to normalize to 0
   subsenupperband=StringToUINT(GetFromFile("uppervimg"))-1;//minus 1 to normalize to 0
   subsennumberofbands=subsenupperband-subsenlowerband+1;
}

//-------------------------------------------------------------------------
// Overloaded function for use by Fenix file reader to get the 
// correct item from the header for the bands that are being processed
//-------------------------------------------------------------------------
std::string FenixBinFile::GetFromFile(std::string keyword)
{
   if(keyword.compare("spatialbinning")==0)
   {
      if(subsensor==VNIR)
         return br->FromHeader("binning",1,"true");
      else if(subsensor==SWIR)
         return br->FromHeader("binning2",1,"true");
   }
   else if(keyword.compare("spectralbinning")==0)
   {
      if(subsensor==VNIR)
         return br->FromHeader("binning",0,"true");
      else if(subsensor==SWIR)
         return br->FromHeader("binning2",0,"true");
   }
   else if(keyword.compare("binningForHeader")==0)
   {
      if(subsensor==VNIR)
      {
         std::string strspecvnir=br->FromHeader("binning",0,"true");
         std::string strspatvnir=br->FromHeader("binning",1,"true");
         return "binning_VNIR = {"+strspecvnir+","+strspatvnir+"}";
      }
      else if(subsensor==SWIR)
      {
         std::string strspecswir=br->FromHeader("binning2",0,"true");
         std::string strspatswir=br->FromHeader("binning2",1,"true");
         return "binning_SWIR = {"+strspecswir+","+strspatswir+"}";
      }
   }
   else if(keyword.compare("lowerhimg")==0)
   {
      if(subsensor==VNIR)
         return br->FromHeader("himg1",0,"true");
      else if(subsensor==SWIR)
         return br->FromHeader("himg2",0,"true");
   }
   else if(keyword.compare("upperhimg")==0)
   {
      if(subsensor==VNIR)
         return br->FromHeader("himg1",1,"true");
      else if(subsensor==SWIR)
         return br->FromHeader("himg2",1,"true");
   }
   else if(keyword.compare("lowervimg")==0)
   {
      if(subsensor==VNIR)
         return br->FromHeader("vimg1",0,"true");
      else if(subsensor==SWIR)
         return br->FromHeader("vimg2",0,"true");
   }
   else if(keyword.compare("uppervimg")==0)
   {
      if(subsensor==VNIR)
         return br->FromHeader("vimg1",1,"true");
      else if(subsensor==SWIR)
         return br->FromHeader("vimg2",1,"true");
   }
   else if(keyword.compare("integrationtime")==0)
   {
      if(subsensor==VNIR)
         return br->FromHeader("tint1","true");
      else if(subsensor==SWIR)
         return br->FromHeader("tint2","true");
   }
   else if(keyword.compare("tintForHeader")==0)
   {
      if(subsensor==VNIR)
         return "tint_VNIR = "+br->FromHeader("tint1","true");
      else if(subsensor==SWIR)
         return "tint_SWIR = "+br->FromHeader("tint2","true");
   }
   else if(keyword.compare("subsensorBandsForHeader")==0)
   {
      if(subsensor==VNIR)
         return "bandrange_VNIR = {"+GetFromFile("lowervimg")+","+GetFromFile("uppervimg")+"}";
      else if(subsensor==SWIR)
         return "bandrange_SWIR = {"+GetFromFile("lowervimg")+","+GetFromFile("uppervimg")+"}";
   }
   else if(keyword.compare("Wavelength")==0)
   {
      std::string wlv="";
      for(unsigned int b=0;b<subsennumberofbands;b++)
      {
         wlv+=br->FromHeader("Wavelength",b+subsenlowerband) +";";
      }
      return TrimPunctuation(wlv);
   }
   else if(keyword.compare("fhwm")==0)
   {
      std::string fwhm="";
      for(unsigned int b=0;b<subsennumberofbands;b++)
      {
         fwhm+=br->FromHeader("fhwm",b+subsenlowerband) +";";
      }
      return TrimPunctuation(fwhm);
   }
   else
      return br->FromHeader(keyword);
   
   //It should never reach here but just in case
   return "";
}

//-------------------------------------------------------------------------
// Overloaded BinFile function to take into account number of bands of 
// subsensor is not the same as that given in hdr file
//-------------------------------------------------------------------------
double FenixBinFile::ReadCell(const unsigned int band,const unsigned int line, const unsigned int col)
{
   unsigned int bandoffset=0;
   if((band == 0) && (col == 0))
   {
      //We will force the function to read from band 0, sample 0 of the file (not the band) - this is the frame counter
      bandoffset=0;
   }
   else
   {
      //We will read in the requested cell from the file using the subsensor lower band as offset
      bandoffset=subsenlowerband;
   }

   unsigned int updatedband=band+bandoffset;
   return br->ReadCell(updatedband,line,col);;
}

//-------------------------------------------------------------------------
// Overloaded BinFile function to take into account number of bands of 
// subsensor is not the same as that given in hdr file
//-------------------------------------------------------------------------
void FenixBinFile::Readlines(char* const chdata, unsigned int startline, unsigned int numlines)
{
   unsigned int startband=subsenlowerband;
   uint64_t nbytestoread_onebandline=NumSamples()*GetDataSize();
   for(unsigned int line=0;line<numlines;line++)
   {
      for(unsigned int band=0;band<SubsensorNumOfBands();band++)
      {
         br->Readbandline(&(chdata[line*band*nbytestoread_onebandline + band*nbytestoread_onebandline]),startband+band,startline+line);
      }
   }     
}

//-------------------------------------------------------------------------
// Overloaded BinFile function to take into account number of bands of 
// subsensor is not the same as that given in hdr file
//-------------------------------------------------------------------------
void FenixBinFile::ReadlineToDoubles(double* const ddata,unsigned int line)
{
   char* chtmp=new char[NumSamples()*SubsensorNumOfBands()*GetDataSize()];
   this->Readlines(chtmp,line,1);
   for(unsigned int sample=0;sample<NumSamples()*SubsensorNumOfBands();sample++)
   {
      ddata[sample]=br->DerefToDouble(&chtmp[sample*GetDataSize()]);
   }   
   delete[] chtmp;
}


//-------------------------------------------------------------------------
// Specim Sensor constructor
//-------------------------------------------------------------------------
Specim::Specim(std::string strFilename,bool force)
{
   strRawFilename=strFilename;
   bin=new SpecimBinFile(strRawFilename);

   //Get the number of bands and samples of the raw image 
   //The data array MUST conform to these sizes
   numbands=StringToUINT(bin->FromHeader("bands","true").c_str()); 
   numsamps=StringToUINT(bin->FromHeader("samples","true").c_str());
   numlines=StringToUINT(bin->FromHeader("lines","true").c_str());

   //Create a fodis and then destroy it if the region is not defined
   fodis=new Fodis(bin);
   if(fodis->RegionSize()<=0)
   {
      Logger::Log("Sensor has no FODIS region defined.");
      delete fodis;
      fodis=NULL;
   }

   DARKFORCE=force;
   darkscalar=1.0;
   //Test if dark frames are present
   if(bin->FromHeader("autodarkstartline").compare("")==0)
   {
      Logger::WarnOnce("No autodarklinestart in hdr file. Assuming no dark frames in file.");
      ndarklines=0;
      darklinestart=bin->NumLines(); // set to maximum row number for use in totalmissingframes
   }
   else
   {
      darklinestart=StringToUINT(bin->FromHeader("autodarkstartline"));         
      ndarklines=StringToUINT(bin->FromHeader("lines")) - darklinestart; 
      //Note that Eagle/Hawk will do extra tests here and may update ndarklines and darkstartline
   }

   //Test if dark frames number seems reasonable
   if(ndarklines > StringToUINT(bin->FromHeader("lines"))/2.0 )
   {
      Logger::Warning("More than half the number of scan lines are dark frames - seems a bit odd.");      
   }
  
   //Set these to zero for specim sensor - only the specific sensor type can read these (e.g. Hawk, Eagle, Fenix)
   //Maybe they should be removed from the specim sensor class - along with access functions
   tint=0;
   spectralbinning=0;
   spatialbinning=0; 
   scanlinelowerlimit=0;
   scanlineupperlimit=0;
   rawmax=0;
   calibratedmax=0;
   radscalar=0;
   SENSOR_ID=StringToINT(bin->FromHeader("sensorid").c_str());;
   calibratedunits="nW/(cm)^2/(sr)/(nm)";
}

//-------------------------------------------------------------------------
// Specim destructor
//-------------------------------------------------------------------------
Specim::~Specim()
{
   if(bin!=NULL)
      bin->Close();
   if(fodis!=NULL)
      delete fodis;
}

//-------------------------------------------------------------------------
// Function to be called from a derived class (EAGLE, HAWK, FENIX)
// to read in the tint,binnings and other info from the header file
//-------------------------------------------------------------------------
void Specim::GetExtraInfoFromHeader()
{
   //Get the integration time of the data and store as a double for easy access
   tint=StringToDouble(TrimWhitespace(bin->GetFromFile("integrationtime").c_str()));
   if(tint == 0)
   {
      throw "Integration time from raw file "+strRawFilename+" is 0. This is not good.";
   }

   //Get the spectral and spatial binning value - store as uint for easy future access
   spectralbinning=StringToUINT(TrimWhitespace(bin->GetFromFile("spectralbinning").c_str()));
   spatialbinning=StringToUINT(TrimWhitespace(bin->GetFromFile("spatialbinning").c_str()));

   if((spectralbinning==0)||(spatialbinning==0))
   {
      Logger::Log("Either spectral or spatial binning is 0");
      throw "Error: Either spectral or spatial binning is 0. Spectral: "+ToString(spectralbinning)+" Spatial: "+ToString(spatialbinning);
   }

   //Get the upper and lower scan line values - store as uint for easy future access
   scanlinelowerlimit=StringToUINT(TrimWhitespace(bin->GetFromFile("lowerhimg").c_str())) -1; //minus 1 to normalize to 0
   scanlineupperlimit=StringToUINT(TrimWhitespace(bin->GetFromFile("upperhimg").c_str())) -1; //minus 1 to normalize to 0
}

//-------------------------------------------------------------------------
// TotalMissingFrames function to use frame counter to identify missing
// frames between the passed start and end frame variables
//-------------------------------------------------------------------------
unsigned short Specim::TotalMissingFrames(const unsigned int start, const unsigned int end)
{
   unsigned short int total=0;
   //Start should be less than end - just return 0 as it is true there are 0. Could maybe warn user also?
   if(start >= end)
      return 0;

   //Updating end to be end-1 as we only process in the main loop from start < end (not start <= end)
   unsigned int newend=end-1;

   double first=bin->ReadCell(0,start,0);
   double last=bin->ReadCell(0,newend,0);

   if((last > first)&&((newend-start)<MAXFRAMECOUNT))
   {
      //last is greater than first (this is good!)
      //and end-start is less than MAXFRAMECOUNT ... therefore at most only one wrap around (due to uint16)
      //so since last > first there can't be a wrap around
      total=last-first +1;//total number of frames (inc. missing ones)
   }
   else if((last < first)&&((newend-start)<MAXFRAMECOUNT))
   {
      //last is less than first but still less than MAXFRAMECOUNT lines
      //therefore only one wrap around
      Logger::Log("Frame count has wrapped around due to short integer. This should be ok.");
      total=(last + (MAXFRAMECOUNT+1)) - first +1;
   }
   else
   {
      //more than one wrap around possibly occurred - estimate the number of wraps
      int numofwraps=(newend-start)/MAXFRAMECOUNT;
      Logger::Log("Frame count has wrapped around more than once - estimated total number of frames according to counter.");
      total=last+(numofwraps*(MAXFRAMECOUNT+1)) - first +1;
   }

   unsigned short tmp=(total-(newend - start+1));//total number of missing frames
   Logger::Log("First frame: "+ToString(first)+" Last frame: "+ToString(last)+" Total frames: "+ToString(total));
   Logger::Log("Total number of missing frames are: "+ToString(tmp));

   return tmp;   
}

//-------------------------------------------------------------------------
// Function to read in the dark frames and return in dlstorage
//-------------------------------------------------------------------------
void Specim::ReadInAllDarkFrames(double* dlstorage,std::string externalFileName,uint64_t linecellsize,unsigned int subsensor)
{
   //Check the linecellsize - this is the number of cells of 1 line of the dlstorage array (i.e. numsamples*nbands)
   //has a default value of 0 meaning to use file dimensions   

   if(externalFileName.compare("")==0)
   {
      if(linecellsize==0)
         linecellsize=NumSamples()*NumBands();

      if(linecellsize!=NumSamples()*NumBands())
         Logger::Log("Reading in dark frames using an array of different size to number of samples * number of bands. This should only happen for Fenix sensors.");

      //Read in numdarklines of data
      for(unsigned int dl=0;dl<this->ndarklines;dl++)
      {
         bin->ReadlineToDoubles(&dlstorage[dl*linecellsize],darklinestart+dl);    
      }
   }
   else
   {
      //Read dark frames from another file
      //First try to open the BIL/BSQ file
      SpecimBinFile* dark=NULL;
      //If we're looking at a fenix raw file then we expect fenix dark frames - we need a fenix bil reader to read them
      if(CheckSensorID(FENIX,SENSOR_ID))
      {
         dark=new FenixBinFile(externalFileName);
      }
      else
      {
         dark=new SpecimBinFile(externalFileName);
      }

      //Set up data arrays
      //Get the size in bytes of 1 lines worth of data (for all bands)
      unsigned int nsamples=StringToUINT(dark->FromHeader("samples"));
      unsigned int nbands=StringToUINT(dark->FromHeader("bands"));
      unsigned int nlines=StringToUINT(dark->FromHeader("lines"));

      if(linecellsize==0)
         linecellsize=nsamples*nbands;

      if(linecellsize!=nsamples*nbands)
      {
         Logger::Log("Reading in dark frames using an array of different size to number of samples * number of bands. This should only happen for Fenix sensors.");
         //Open up for reading specific subsensor as defined by the subsensor parameter
         dark->SetSubSensor(subsensor);
      }

      //Read in darklines from file 
      Logger::Log("Reading in dark lines...");

      //For each line of data
      for(unsigned int line=0;line<nlines;line++)
      {
         //Read in the line of data for all bands
         dark->ReadlineToDoubles(&dlstorage[line*linecellsize],line); 
      }
      delete dark;
   }
}

//-------------------------------------------------------------------------
// Function to average all the dark frames into data array
//-------------------------------------------------------------------------
void Specim::AverageAllDarkFrames(double* const data,std::string externalFileName,uint64_t linecellsize,unsigned int subsensor)
{
   double* dlstorage=NULL;
   unsigned int nlines=0;
   unsigned int nsamples=0;
   if(externalFileName.compare("")==0)
   {
      //Do not use an external file and read from this one
      if(linecellsize==0)
         linecellsize=NumSamples()*NumBands();
      nsamples=NumSamples();
      nlines=this->ndarklines;      
   }
   else
   {
      //Read dark frames from another file
      //First try to open the BIL/BSQ file
      SpecimBinFile* dark=NULL;
      if(CheckSensorID(FENIX,SENSOR_ID))
      {
         dark=new FenixBinFile(externalFileName);
      }
      else
      {
         dark=new SpecimBinFile(externalFileName);
      }
      //This does nothing if not a FenixBinFile
      dark->SetSubSensor(subsensor);

      //Set up data arrays
      //Get the size in bytes of 1 lines worth of data (for all bands)
      nsamples=StringToUINT(dark->FromHeader("samples"));
      unsigned int nbands=StringToUINT(dark->FromHeader("bands"));
      nlines=StringToUINT(dark->FromHeader("lines"));

      if(linecellsize==0)
         linecellsize=nsamples*nbands;

      //Update darkscalar value
      double raw_integration_time=StringToDouble(bin->GetFromFile("integrationtime"));
      double dark_integration_time=StringToDouble(dark->GetFromFile("integrationtime"));
      this->darkscalar=raw_integration_time/dark_integration_time;
      Logger::Log("Updated dark scalar based on integration times:"+ToString(this->darkscalar));
      dark->Close();
      delete dark;
   }

   //Read in the dark frames casting as doubles
   dlstorage=new double[nlines*linecellsize];
   ReadInAllDarkFrames(dlstorage, externalFileName,linecellsize,subsensor);
   //Sum up the darkframes for each sample of each band
   for(unsigned int line=0;line<nlines;line++)
   {
      data[0]=0; //first element is 1st sample of 1st band which is frame number, so set to 0
      for(unsigned int s=1;s<linecellsize;s++)//loop from 1 onwards
      {
         data[s]=data[s]+dlstorage[line*linecellsize + s];
      }  
   }

   //Average them up
   for(unsigned int s=0;s<linecellsize;s++)
   {
      data[s]=data[s]/static_cast<double>(nlines);
   }

   delete[] dlstorage;

}

//-------------------------------------------------------------------------
// Function to calculate the standard deviation of the dark frames
//-------------------------------------------------------------------------
void Specim::DarkFramesStdDeviation(double* const stdev,const double* const mean, std::string externalFileName,uint64_t linecellsize,unsigned int subsensor)
{
   double* dlstorage=NULL;
   unsigned int nlines=0;
   unsigned int nsamples=0;
   if(externalFileName.compare("")==0)
   {
      //Do not use an external file and read from this one
      if(linecellsize==0)
         linecellsize=NumSamples()*NumBands();
      nsamples=NumSamples();
      nlines=this->ndarklines;
   }
   else
   {
      //Read dark frames from another file
      //First try to open the BIL/BSQ file
      SpecimBinFile* dark=NULL;
      if(CheckSensorID(FENIX,SENSOR_ID))
      {
         dark=new FenixBinFile(externalFileName);
      }
      else
      {
         dark=new SpecimBinFile(externalFileName);
      }
      //This does nothing if not a FenixBinFile
      dark->SetSubSensor(subsensor);

      //Set up data arrays
      //Get the size in bytes of 1 lines worth of data (for all bands)
      nsamples=StringToUINT(dark->FromHeader("samples"));
      unsigned int nbands=StringToUINT(dark->FromHeader("bands"));
      nlines=StringToUINT(dark->FromHeader("lines"));

      if(linecellsize==0)
         linecellsize=nsamples*nbands;

      //Update darkscalar value
      double raw_integration_time=StringToDouble(bin->GetFromFile("integrationtime"));
      double dark_integration_time=StringToDouble(dark->GetFromFile("integrationtime"));
      this->darkscalar=raw_integration_time/dark_integration_time;
      Logger::Log("Updated dark scalar based on integration times:"+ToString(this->darkscalar));
      dark->Close();
      delete dark;
   }

   //Read in the dark frames casting as doubles
   dlstorage=new double[nlines*linecellsize];
   ReadInAllDarkFrames(dlstorage, externalFileName,linecellsize,subsensor);

   for(unsigned int line=0;line<nlines;line++)
   {  
      stdev[0]=0; //first element is 1st sample of 1st band which is frame number, so set to 0
      for(unsigned int s=1;s<linecellsize;s++)//loop from 1 onwards
      {
         stdev[s]=stdev[s]+pow((dlstorage[line*linecellsize + s]-mean[s]),2);
      }  
   }   

   //Now finalise the standard devs
   for(unsigned int s=0;s<linecellsize;s++)
   {
      if(nlines > 1)
         stdev[s]=sqrt(stdev[s]/(nlines-1.0));
      else
         stdev[s]=sqrt(stdev[s]/(nlines));
   }

   delete[] dlstorage;
}

//-------------------------------------------------------------------------
// Function to check dark values against standard deviation and mean, then
// remove ones it considers outliers before calculating the average
//-------------------------------------------------------------------------
void Specim::AverageRefinedDarkFrames(double* const data,const double* const stdev,const double* const mean, std::string externalFileName,uint64_t linecellsize,unsigned int subsensor)
{
   unsigned long int* numitems=NULL;
   double* dlstorage=NULL;
   unsigned int nlines=0;
   unsigned int nsamples=0;
   if(externalFileName.compare("")==0)
   {
      //Do not use an external file and read from this one
      if(linecellsize==0)
         linecellsize=NumSamples()*NumBands();
      nsamples=NumSamples();
      nlines=this->ndarklines;      
   }
   else
   {
      //Read dark frames from another file
      //First try to open the BIL/BSQ file
      SpecimBinFile* dark=NULL;
      if(CheckSensorID(FENIX,SENSOR_ID))
      {
         dark=new FenixBinFile(externalFileName);
      }
      else
      {
         dark=new SpecimBinFile(externalFileName);
      }
      //This does nothing if not a FenixBinFile
      dark->SetSubSensor(subsensor);

      //Set up data arrays
      //Get the size in bytes of 1 lines worth of data (for all bands)
      nsamples=StringToUINT(dark->FromHeader("samples"));
      unsigned int nbands=StringToUINT(dark->FromHeader("bands"));
      nlines=StringToUINT(dark->FromHeader("lines"));

      if(linecellsize==0)
         linecellsize=nsamples*nbands;

      //Update darkscalar value
      double raw_integration_time=StringToDouble(bin->GetFromFile("integrationtime"));
      double dark_integration_time=StringToDouble(dark->GetFromFile("integrationtime"));
      this->darkscalar=raw_integration_time/dark_integration_time;
      Logger::Log("Updated dark scalar based on integration times, using a scalar of:"+ToString(this->darkscalar));
      dark->Close();
      delete dark;
   }
   //Read in the dark frames casting as doubles
   dlstorage=new double[nlines*linecellsize];
   ReadInAllDarkFrames(dlstorage, externalFileName,linecellsize,subsensor);

   //Set up the numitems array
   numitems=new unsigned long int[linecellsize];
   for(unsigned int s=0;s<linecellsize;s++)
   {
      numitems[s]=0;
   }

   for(unsigned int line=0;line<nlines;line++)
   {
      data[0]=0; //first element is 1st sample of 1st band which is frame number, so set to 0
      for(unsigned int s=1;s<linecellsize;s++)//loop from 1 onwards
      {
         //If the dark value is within 3 times the stdev from the mean, use the data
         if( (dlstorage[line*linecellsize + s] <= (mean[s] + 3*stdev[s])) && (dlstorage[line*linecellsize + s] >= (mean[s] - 3*stdev[s])))
         {
            data[s]=data[s]+dlstorage[line*linecellsize + s];
            numitems[s]++;//individual counter for each pixel
         }
      }  
   }

   //average them up - skip the first entry
   for(unsigned int s=1;s<linecellsize;s++)
   {
      if(numitems[s] < (nlines/2.0))
      {
         Logger::Warning("Less than half the dark values for this pixel (of the ccd - i.e. 0 to samples*bands) have been used to calculate the average: "+ToString(s));
      }
      data[s]=data[s]/static_cast<double>(numitems[s]);
   }

   delete[] dlstorage;
   delete[] numitems;
}

//-------------------------------------------------------------------------
// Function to check and update darklinestart variable if required
//-------------------------------------------------------------------------
void Specim::DarkFrameSanityCheck()
{
   if(ndarklines != 0)
   {
      //Updated test - check frame counter between (darkstart-1,darkstart) and (darkstart,darkstart+1)
      //when dark frames start there is a jump in frame counter, can use this as a check to make
      //sure that the darklinestart relates to first dark rather than last light
      Logger::Log("Checking frame counter just before auto dark start...");
      unsigned int prevframejump=this->GetMissingFramesBetweenLimits(darklinestart-1,darklinestart+1);
      Logger::Log("Checking frame counter just after auto dark start...");
      unsigned int nextframejump=this->GetMissingFramesBetweenLimits(darklinestart,darklinestart+2);
      if(prevframejump > nextframejump)
      {
         //Do not need to do anything as darklinestart points to first dark line
         Logger::Log("auto dark start appears to point at first dark line.");
      }
      else if(nextframejump > prevframejump)
      {
         //Need to increase darklinestart by 1 as points to last light line
         //This should be the case for most Eagle data but it appears to be not all.
         Logger::Log("auto dark start appears to point at last light line - adding one on to point at first dark line.");
         darklinestart++;          
      }
      else
      {
         //There is the same frame counter jump in each case - darkstart is probably incorrect
         //However - Eagle 110001 appears to behave differently and have no frame jump at dark start!
         if((bin->FromHeader("sensorid").compare("110001")!=0)&&(this->DARKFORCE==false))
         {
            throw "Autodarkstartline may not be correct - there is no frame counter jump between light and dark frames. "
                   "Please check the autodarkstartline in the raw hdr file is correct (check vs the raw data in ENVI where the dark frames start). "
                   "It should be the first line of dark data. "
                   "If you are sure you wish to process the data using this autodarkstartline then use -darkforce on command line.";
         }
         else if((bin->FromHeader("sensorid").compare("110001")!=0)&&(this->DARKFORCE==true))
         {
            //If the user has opted to process using the autodarkstart value when there is no jump in frame number and not sensor 110001
            Logger::Log("User has requested to force the use of the autodarkstartline in the hdr file even though there was "
                        "no frame counter jump between dark and light frames. Will now assume autodarkstartline is pointing to first dark frame.");
            
         }
         else
         {
            //Now do a different test - this is required because 2nd Eagle (SN110001) has no frame counter jump
            //Can't think of a suitable test here - looking at some data it appears to point to previous light line
            //as with the other Eagle (SN100022) so will increase by 1 until new test can be found.
            Logger::Log("auto dark start appears to point at last light line - adding one on to point at first dark line.");
            darklinestart++;
         }
      }
         
      ndarklines=StringToUINT(bin->FromHeader("lines")) - darklinestart; 
   }
}

//-------------------------------------------------------------------------
// Eagle constructor
//-------------------------------------------------------------------------
Eagle::Eagle(std::string strFilename,bool force)  : Specim(strFilename,force)
{
   rawmax=EAGLE_RAW_MAX;
   calibratedmax=CALIBRATED_DATA_MAX;
   radscalar=RADIANCE_DATA_SCALAR;
   trant=FRAME_TRANSFER_TIME;

   //Check the sensor id
   bool iseagle=CheckSensorID(EAGLE,SENSOR_ID);
   if(!iseagle)
   {
      Logger::Warning("Sensor ID from Header file is not an Eagle ID - Constructing Eagle object from a non-eagle data file.");
   }

   //Create an EagleBiNFile
   if(bin!=NULL)
      delete bin;
   bin=new EagleHawkBinFile(strFilename);

   GetExtraInfoFromHeader();
   TotalMissingFrames();
   //Get the fodis lower and upper bounds - store as uints for future access  
   //Subtract 1 from them because they appear to be indexed from 1 - n in hdr file
   //lowerfodis=StringToUINT(TrimWhitespace(br->FromHeader("fodis",0).c_str())) -1;
   //upperfodis=StringToUINT(TrimWhitespace(br->FromHeader("fodis",1).c_str())) -1;
   //fodisunits="nW/(cm)^2/(nm)";
}

//-------------------------------------------------------------------------
// Eagle destructor
//-------------------------------------------------------------------------
Eagle::~Eagle()
{
}

//-------------------------------------------------------------------------
//Reads the first/last frame IDs and compare to number of lines to get total number of missing frames
//-------------------------------------------------------------------------
void Eagle::TotalMissingFrames()
{
   //Frame numbers are the first 2 bytes of each frame on band 1.
   //we need to read the first and last frame (excluding dark frames)
   DarkFrameSanityCheck();
   //Use end limits to do whole flight line
   totalmissing=Specim::TotalMissingFrames(0,darklinestart);
} 


//-------------------------------------------------------------------------
// Hawk constructor
//-------------------------------------------------------------------------
Hawk::Hawk(std::string strFilename,bool force) : Specim(strFilename,force)
{
   rawmax=HAWK_RAW_MAX;
   calibratedmax=CALIBRATED_DATA_MAX;
   radscalar=RADIANCE_DATA_SCALAR;

   //Check the sensor id
   bool ishawk=CheckSensorID(HAWK,SENSOR_ID);
   if(!ishawk)
   {
      Logger::Warning("Sensor ID from Header file is not a Hawk ID - Constructing Hawk object from a non-hawk data file.");
   }

   //Create an EagleHawkBilReader
   if(bin!=NULL)
      delete bin;
   bin=new EagleHawkBinFile(strFilename);

   GetExtraInfoFromHeader();
   TotalMissingFrames();
}

//-------------------------------------------------------------------------
// Hawk destructor
//-------------------------------------------------------------------------
Hawk::~Hawk()
{
}

//-------------------------------------------------------------------------
//Reads the first/last frame IDs and compare to number of lines to get total number of missing frames
//-------------------------------------------------------------------------
void Hawk::TotalMissingFrames()
{
   //Frame numbers are the first 2 bytes of each frame on band 1.
   //we need to read the first and last frame (excluding dark frames)
   
   DarkFrameSanityCheck();
   //Use end limits to do whole flight line
   this->totalmissing=Specim::TotalMissingFrames(0,darklinestart);
} 

//-------------------------------------------------------------------------
// Fenix sensor constructor
//-------------------------------------------------------------------------
Fenix::Fenix(std::string strFilename): Specim(strFilename)
{
   calibratedmax=CALIBRATED_DATA_MAX;
   radscalar=RADIANCE_DATA_SCALAR;

   //Check the sensor id
   if(!CheckSensorID(FENIX,SENSOR_ID))
   {
      Logger::Warning("Sensor ID from Header file is not a Fenix ID - Constructing Fenix object from a non-fenix data file.");
   }

   //Create a FenixBinFile for both parts of file
   vnbr=new FenixBinFile(strFilename);
   swbr=new FenixBinFile(strFilename);
   swbr->SetSubSensor(SWIR);

   //Clean up br if not already (it is probably not)
   if(bin!=NULL)
      delete bin;

   //Default to VNIR section of data 
   SetUpFenixFor(VNIR);

   TotalMissingFrames();
}

//-------------------------------------------------------------------------
// Fenix sensor destructor
//-------------------------------------------------------------------------
Fenix::~Fenix()
{
   if(vnbr!=NULL)
      delete vnbr;
   if(swbr!=NULL)
      delete swbr;
   bin=NULL;
}

//-------------------------------------------------------------------------
// Function to set up the Fenix class to read one of the subsensor sections
// of data (VNIR or SWIR)
//-------------------------------------------------------------------------
void Fenix::SetUpFenixFor(Subsensor sen)
{
   //Set up the specimbilreader to be a fenixbilreader
   if(sen==VNIR)
   {
      bin=vnbr;
      rawmax=FENIX_VNIR_RAW_MAX;
   }
   else if(sen==SWIR)
   {
      bin=swbr;
      rawmax=FENIX_SWIR_RAW_MAX;
   }
   else
      throw "Unknown subsensor type given in SetUpFenixFor().";

   GetExtraInfoFromHeader();

   //Fenix spatial binning is software binning - this means that raw maximum value can be 4095*(spat_binning)
   if(static_cast<uint64_t>(rawmax)*spatialbinning < std::numeric_limits<uint16_t>::max())
   {
      rawmax=rawmax*spatialbinning;
   }
   else
   {
      Logger::WarnOnce("The spatial binning of this data means that the maximum raw value would be > 16-bit maximum. Beware that data may not be correct if the software binning wraps rather than caps at the maximum.");
      rawmax=std::numeric_limits<uint16_t>::max();
   }

}

//-------------------------------------------------------------------------
//Reads the first/last frame IDs and compare to number of lines to get total number of missing frames
//-------------------------------------------------------------------------
void Fenix::TotalMissingFrames()
{
   //Frame numbers are the first 2 bytes of each frame on band 1.
   //we need to read the first and last frame (excluding dark frames)
   
   //We start from after first frame as the fenix drops frames until the whole second pps
   //But it collects one frame before dropping loads - we ignore missing frames between scanline 1 and scanline 2
   this->totalmissing=Specim::TotalMissingFrames(1,darklinestart);
} 


//-------------------------------------------------------------------------
// Fodis constructor
//-------------------------------------------------------------------------
Fodis::Fodis(SpecimBinFile* bin)
{
   lowerfodis=StringToUINT(TrimWhitespace(bin->FromHeader("fodis",0).c_str())) -1;
   upperfodis=StringToUINT(TrimWhitespace(bin->FromHeader("fodis",1).c_str())) -1;
   fodisunits="nW/(cm)^2/(nm)";   
}

