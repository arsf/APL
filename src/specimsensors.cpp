//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------


#include "specimsensors.h"

#ifndef SPECIMDEBUG
   #define DEBUGPRINT(X)  
#else
   #define DEBUGPRINT(X) std::cout<< X <<std::endl;
#endif


//-------------------------------------------------------------------------
// Function to test the suitability of the file
//-------------------------------------------------------------------------
std::string SpecimBILReader::TestFile()
{
   //Test if dark frames are present
   if(this->FromHeader("autodarkstartline").compare("")==0)
   {
      Logger::Log("WARNING: No autodarklinestart in hdr file. Assuming no dark frames in file.");
      ndarklines=0;
      darklinestart=this->numrows; // set to maximum row number for use in totalmissingframes
   }
   else
   {
      darklinestart=StringToUINT(this->Header["autodarkstartline"]);
      //Annoying thing - in Eagle files it appears that autodarkstartline is the (0-based) 
      //line number for the last image line, whereas for Hawk it appears that it is the
      //line number for the first dark frame
      //if(StringToINT(this->FromHeader("sensorid"))==EAGLE_SENSOR_ID)
      //   darklinestart++; 

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
         if((this->FromHeader("sensorid").compare("110001")!=0)&&(this->DARKFORCE==false))
         {
            throw "Autodarkstartline may not be correct - there is no frame counter jump between light and dark frames. "
                   "Please check the autodarkstartline in the raw hdr file is correct (check vs the raw data in ENVI where the dark frames start). "
                   "It should be the first line of dark data. "
                   "If you are sure you wish to process the data using this autodarkstartline then use -darkforce on command line.";
         }
         else if((this->FromHeader("sensorid").compare("110001")!=0)&&(this->DARKFORCE==true))
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
         
      ndarklines=StringToUINT(FromHeader("lines")) - darklinestart; 
   }

   //Test if dark frames number seems reasonable
   if(ndarklines > StringToUINT(FromHeader("lines"))/2.0 )
   {
      Logger::Log("WARNING: More than half the number of scan lines are dark frames - seems a bit odd.");      
   }

   //Test if unsigned 16-bit integer
   if(this->datasize!=2)  
      throw BRexception("Cannot read from Specim File as data size was expected to be 2 bytes");   

   //Calculate the number of missing frames
   TotalMissingFrames();

   //return the sensorid to identify if eagle or hawk file
   return this->FromHeader("sensorid");
}


//-------------------------------------------------------------------------
// Function to average all the dark frames into data array
//-------------------------------------------------------------------------
void SpecimBILReader::AverageAllDarkFrames(double* const data,std::string externalFileName)
{

   if(externalFileName.compare("")==0)
   {
      //Do not use an external file and read from this one
      //unsigned int dlinesstart=StringToUINT(this->Header["autodarkstartline"]);
      //dlinesstart++; //increase to get a more sensible value (i.e the start of dark lines if counting from 0)      
      //unsigned int numdarklines=this->ndarklines;

      unsigned short int* dlstorage=NULL;  
      dlstorage=new unsigned short int[this->numsamples*this->numbands];

      //Skip to start of dark lines
      uint64_t bytestoskip=(darklinestart)*this->numsamples*this->numbands*this->datasize;
      FileSeek(filein,bytestoskip,SEEK_SET);
      //Read in numdarklines of data
      for(unsigned int dl=0;dl<this->ndarklines;dl++)
      {
         //Check its ok to read from the file (this doesnt check anything to do with size of data to be read)
         if(this->IsGood())
            fread((char*)dlstorage,sizeof(char),this->numsamples*this->numbands*this->datasize,filein);
         else
            throw BRexception("Trying to read from a bad stream in GetAverageDarkFrames");         

         data[0]=0; //first element is 1st sample of 1st band which is frame number, so set to 0
         for(unsigned int s=1;s<this->numsamples*this->numbands;s++)//loop from 1 onwards
         {
            data[s]=data[s]+static_cast<double>(dlstorage[s]);
         }  
      }

      //Average them up
      for(unsigned int s=0;s<this->numsamples*this->numbands;s++)
      {
         data[s]=data[s]/static_cast<double>(this->ndarklines);
      }

      delete[] dlstorage;
   }
   else
   {
      //Read dark frames from another file
      //First try to open the BIL/BSQ file
      BinFile dark(externalFileName);

      //Set up data arrays
      //Get the size in bytes of 1 lines worth of data (for all bands)
      unsigned int nsamples=StringToUINT(dark.FromHeader("samples"));
      unsigned int nbands=StringToUINT(dark.FromHeader("bands"));
      unsigned int nlines=StringToUINT(dark.FromHeader("lines"));

      //Update darkscalar value
      double raw_integration_time=StringToDouble(this->FromHeader("tint"));
      double dark_integration_time=StringToDouble(dark.FromHeader("tint"));
      this->darkscalar=raw_integration_time/dark_integration_time;
      Logger::Log("Updated dark scalar based on integration times:"+ToString(this->darkscalar));

      //Read in darklines from file 
      Logger::Log("Reading in dark lines...");

      double* filedata=new double[nsamples*nbands];
      //For each line of data
      for(unsigned int line=0;line<nlines;line++)
      {
         //Read in the line of data for all bands
         dark.ReadlineToDoubles(filedata,line);
         //sum up a running total per pixel per band
         for(unsigned int sample=0;sample<nsamples*nbands;sample++)
         {
            data[sample]=data[sample] + filedata[sample];
         }       
      }

      Logger::Log("Averaging dark lines...");
      //Now convert the sum to the average
      for(unsigned int sample=0;sample<nsamples*nbands;sample++)
      {
         data[sample]=data[sample]/static_cast<double>(nlines);
      }   

      Logger::Log("Have averaged "+ToString(nlines)+" lines of dark data from file "+externalFileName);

      //Tidy up
      delete[] filedata;
      dark.Close();
   }

}

//-------------------------------------------------------------------------
// Function to calculate the standard deviation of the dark frames
//-------------------------------------------------------------------------
void SpecimBILReader::DarkFramesStdDeviation(double* const stdev,const double* const mean, std::string externalFileName)
{
   if(externalFileName.compare("")==0)
   {
      //Do not use an external file and read from this one
      //unsigned int dlinesstart=StringToUINT(this->Header["autodarkstartline"]);
      //dlinesstart++; //increase to get a more sensible value (i.e the start of dark lines if counting from 0)      
      //unsigned int numdarklines=this->ndarklines;

      unsigned short int* dlstorage=NULL;  
      dlstorage=new unsigned short int[this->numsamples*this->numbands];

      //Skip to start of dark lines
      uint64_t bytestoskip=(darklinestart)*this->numsamples*this->numbands*this->datasize;
      FileSeek(filein,bytestoskip,SEEK_SET);
      //Read in numdarklines of data
      for(unsigned int dl=0;dl<this->ndarklines;dl++)
      {
         //Check its ok to read from the file (this doesnt check anything to do with size of data to be read)
         if(this->IsGood())
            fread((char*)dlstorage,sizeof(char),this->numsamples*this->numbands*this->datasize,filein);
         else
            throw BRexception("Trying to read from a bad stream in GetAverageDarkFrames");         

         stdev[0]=0; //first element is 1st sample of 1st band which is frame number, so set to 0
         for(unsigned int s=1;s<this->numsamples*this->numbands;s++)//loop from 1 onwards
         {
            stdev[s]=stdev[s]+pow((static_cast<double>(dlstorage[s])-mean[s]),2);
         }  
      }   

      //Now finalise the standard devs
      for(unsigned int s=0;s<this->numsamples*this->numbands;s++)
      {
         if(this->ndarklines > 1)
            stdev[s]=sqrt(stdev[s]/(this->ndarklines-1.0));
         else
            stdev[s]=sqrt(stdev[s]/(this->ndarklines));
      }

      delete[] dlstorage;

   }
   else
   {
      //Read dark frames from another file
      //First try to open the BIL/BSQ file
      BinFile dark(externalFileName);

      //Set up data arrays
      //Get the size in bytes of 1 lines worth of data (for all bands)
      unsigned int nsamples=StringToUINT(dark.FromHeader("samples"));
      unsigned int nbands=StringToUINT(dark.FromHeader("bands"));
      unsigned int nlines=StringToUINT(dark.FromHeader("lines"));

      //Read in darklines from file 
      Logger::Log("Reading in dark lines for std deviation calculation...");

      double* filedata=new double[nsamples*nbands];
      //For each line of data
      for(unsigned int line=0;line<nlines;line++)
      {
         //Read in the line of data for all bands
         dark.ReadlineToDoubles(filedata,line);

         for(unsigned int s=0;s<this->numsamples*this->numbands;s++)
         {
            stdev[s]=stdev[s]+pow((filedata[s])-mean[s],2);
         }        
      }

      //Now calculate the (unbiased) stdev
      for(unsigned int sample=0;sample<nsamples*nbands;sample++)
      {
         if(nlines > 1)
            stdev[sample]=sqrt(stdev[sample]/(nlines-1.0));
         else
            stdev[sample]=sqrt(stdev[sample]/(nlines));
      }   

      //Tidy up
      delete[] filedata;
      dark.Close();
   }
}

//-------------------------------------------------------------------------
// Function to check dark values against standard deviation and mean, then
// remove ones it considers outliers before calculating the average
//-------------------------------------------------------------------------
void SpecimBILReader::AverageRefinedDarkFrames(double* const data,const double* const stdev,const double* const mean, std::string externalFileName)
{
   if(externalFileName.compare("")==0)
   {
      //Do not use an external file and read from this one
      //unsigned int dlinesstart=StringToUINT(this->Header["autodarkstartline"]);
      //dlinesstart++; //increase to get a more sensible value (i.e the start of dark lines if counting from 0)      
      //unsigned int numdarklines=this->ndarklines;

      unsigned short int* dlstorage=NULL;  
      unsigned long int* numitems=NULL;
      dlstorage=new unsigned short int[this->numsamples*this->numbands];
      numitems=new unsigned long int[this->numsamples*this->numbands];

      for(unsigned int s=0;s<this->numsamples*this->numbands;s++)
      {
         numitems[s]=0;
      }

      //Skip to start of dark lines
      uint64_t bytestoskip=(darklinestart)*this->numsamples*this->numbands*this->datasize;
      FileSeek(filein,bytestoskip,SEEK_SET);
      //Read in numdarklines of data
      for(unsigned int dl=0;dl<this->ndarklines;dl++)
      {
         //Check its ok to read from the file (this doesnt check anything to do with size of data to be read)
         if(this->IsGood())
            fread((char*)dlstorage,sizeof(char),this->numsamples*this->numbands*this->datasize,filein);
         else
            throw BRexception("Trying to read from a bad stream in GetAverageDarkFrames");         

         data[0]=0; //first element is 1st sample of 1st band which is frame number, so set to 0
         for(unsigned int s=1;s<this->numsamples*this->numbands;s++)//loop from 1 onwards
         {
            //If the dark value is within 3 times the stdev from the mean, use the data
            if( (dlstorage[s] <= (mean[s] + 3*stdev[s])) && (dlstorage[s] >= (mean[s] - 3*stdev[s])))
            {
               data[s]=data[s]+static_cast<double>(dlstorage[s]);
               numitems[s]++;//individual counter for each pixel
            }
         }  
      }

      //average them up - skip the first entry
      for(unsigned int s=1;s<this->numsamples*this->numbands;s++)
      {
         if(numitems[s] < (this->ndarklines/2.0))
         {
            Logger::Log("WARNING: Less than half the dark values for this pixel (of the ccd - i.e. 0 to samples*bands) have been used to calculate the average: "+ToString(s));
         }
         data[s]=data[s]/static_cast<double>(numitems[s]);
      }

      delete[] dlstorage;
      delete[] numitems;
   }
   else
   {
      //Read dark frames from another file
      //First try to open the BIL/BSQ file
      BinFile dark(externalFileName);

      //Set up data arrays
      //Get the size in bytes of 1 lines worth of data (for all bands)
      unsigned int nsamples=StringToUINT(dark.FromHeader("samples"));
      unsigned int nbands=StringToUINT(dark.FromHeader("bands"));
      unsigned int nlines=StringToUINT(dark.FromHeader("lines"));

      //Read in darklines from file 
      Logger::Log("Reading in dark lines for refined average calculation...");

      double* filedata=new double[nsamples*nbands];
      unsigned long int* numitems=new unsigned long int[nsamples*nbands];
      for(unsigned int s=0;s<this->numsamples*this->numbands;s++)
      {
         numitems[s]=0;
      }
      //For each line of data
      for(unsigned int line=0;line<nlines;line++)
      {
         //Read in the line of data for all bands
         dark.ReadlineToDoubles(filedata,line);

         for(unsigned int s=0;s<nsamples*nbands;s++)
         {
            //If the dark value is within 3 times the stdev from the mean, use the data
            if( (filedata[s] <= (mean[s] + 3*stdev[s])) && (filedata[s] >= (mean[s] - 3*stdev[s])))
            {
               data[s]=data[s]+filedata[s];
               numitems[s]++;//individual counter for each pixel
            }
         }       
      }

      //average them up 
      for(unsigned int s=0;s<nsamples*nbands;s++)
      {
         if(numitems[s] < (nlines/2.0))
         {
            Logger::Log("WARNING: Less than half the dark values for this pixel (of the ccd - i.e. 0 to samples*bands) have been used to calculate the average: "+ToString(s));
         }
         data[s]=data[s]/static_cast<double>(numitems[s]);
      }  

      //Tidy up
      delete[] filedata;
      delete[] numitems;
      dark.Close();
   }
}



//-------------------------------------------------------------------------
//Reads the first/last frame IDs and compare to number of lines to get total number of missing frames
//-------------------------------------------------------------------------
void SpecimBILReader::TotalMissingFrames()
{
   //Frame numbers are the first 2 bytes of each frame on band 1.
   //we need to read the first and last frame (excluding dark frames)
   
   //unsigned int dlinesstart=0;
   //Check if autodarkstartline is in the hdr file
   if(this->FromHeader("autodarkstartline").compare("")==0)
   {
      Logger::Log("WARNING: No autodarkstartline found in header file. Assuming no dark lines and using max line number in TotalMissingFrames.");
      darklinestart=this->numrows;
   }
   else
   {
      //dlinesstart=StringToUINT(this->Header["autodarkstartline"]);
      if(darklinestart==0)
      {
         //Obvious error has occurred (maybe missing darklines)
         throw BRexception("Error in TotalMissingFrames: Dark frames suggested to start at line 0.");
      }
      //Subtract one off this value since we want the last image line rather than the first darkframe line
      //else it may add dropped scans after the last image line (but before the first dark frame)
      //dlinesstart--;

   }
   //Use end limits to do whole flight line
   this->totalmissing=TotalMissingFrames(0,darklinestart);
} 

//-------------------------------------------------------------------------
//Reads the start/end frame IDs and compare to number of lines within these limits 
//to get total number of missing frames within these limits
//-------------------------------------------------------------------------
unsigned short SpecimBILReader::TotalMissingFrames(const unsigned int start, const unsigned int end)
{
   DEBUGPRINT("In totalmissingframes() using start: "<<start<<" end: "<<end)   
   unsigned short int first=0,last=0,total=0;
   //Get the current file position
   int64_t pos=FileTell(this->filein);

   //Updating end to be end-1 as we only process in the main loop from start < end (not start <= end)
   unsigned int newend=end-1;

   //Find file position of start line
   uint64_t startbytes=start*this->numsamples*this->numbands*this->datasize;
   //move to beginning of start line
   FileSeek(filein,startbytes,SEEK_SET);
   //read 2 bytes
   fread((char*)&first,sizeof(short),1,filein);
   //now find begining position of end line
   uint64_t bytestoskip=(newend)*this->numsamples*this->numbands*this->datasize;
   //move to start of last line band 1
   FileSeek(filein,bytestoskip,SEEK_SET);
   //Check if bytes to skip is greater than the size of the file and exit if true
   if(this->filesize <= bytestoskip)
      throw BRexception("Error in TotalMissingFrames: Trying to read past end of file. Suggests dark frame start line is wrong.");

   //read 2 bytes
   fread((char*)&last,sizeof(short),1,filein);
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

   //now reset the file position
   FileSeek(filein,pos,SEEK_SET);

   return tmp;
} 


//-------------------------------------------------------------------------
// Specim constructor
//-------------------------------------------------------------------------
Specim::Specim(std::string strFilename,bool force)
{
   strRawFilename=strFilename;
   br=new SpecimBILReader(strRawFilename,force);

   //Get the integration time of the data and store as a double for easy access
   tint=StringToDouble(br->FromHeader("tint"));
   if(tint == 0)
   {
      throw "Integration time from raw file "+strFilename+" is 0. This is not good.";
   }

   //Get the spectral and spatial binning value - store as uint for easy future access
   spectralbinning=StringToUINT(TrimWhitespace(br->FromHeader("binning",0).c_str()));
   spatialbinning=StringToUINT(TrimWhitespace(br->FromHeader("binning",1).c_str()));

   if((spectralbinning==0)||(spatialbinning==0))
   {
      Logger::Log("Either spectral or spatial binning is 0");
      throw "Error: Either spectral or spatial binning is 0. Spectral: "+ToString(spectralbinning)+" Spatial: "+ToString(spatialbinning);
   }

   //Get the number of bands and samples of the raw image 
   //The data array MUST conform to these sizes
   numbands=StringToUINT(br->FromHeader("bands").c_str()); 
   numsamps=StringToUINT(br->FromHeader("samples").c_str());
   numlines=StringToUINT(br->FromHeader("lines").c_str());

   //Get the spectral and spatial binning value - store as uint for easy future access
   scanlinelowerlimit=StringToUINT(TrimWhitespace(br->FromHeader("himg",0).c_str())) -1; //minus 1 to normalize to 0
   scanlineupperlimit=StringToUINT(TrimWhitespace(br->FromHeader("himg",1).c_str())) -1; //minus 1 to normalize to 0

   rawmax=0;
   calibratedmax=0;
   radscalar=0;
   SENSOR_ID=StringToINT(br->FromHeader("sensorid").c_str());;

   calibratedunits="nW/(cm)^2/(sr)/(nm)";
}

//-------------------------------------------------------------------------
// Specim destructor
//-------------------------------------------------------------------------
Specim::~Specim()
{
   if(br!=NULL)
      br->Close();
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

   //Get the fodis lower and upper bounds - store as uints for future access  
   //Subtract 1 from them because they appear to be indexed from 1 - n in hdr file
   lowerfodis=StringToUINT(TrimWhitespace(br->FromHeader("fodis",0).c_str())) -1;
   upperfodis=StringToUINT(TrimWhitespace(br->FromHeader("fodis",1).c_str())) -1;

   fodisunits="nW/(cm)^2/(nm)";
}

//-------------------------------------------------------------------------
// Eagle destructor
//-------------------------------------------------------------------------
Eagle::~Eagle()
{
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
   bandnotinuse = 9999;
   nbadpixels=0;
   badpixels=NULL;
   badpixelmethod=NULL;
   bpmethod_descriptor=NULL;
   arsfbadpixelfiletype=false;
}

//-------------------------------------------------------------------------
// Hawk destructor
//-------------------------------------------------------------------------
Hawk::~Hawk()
{
   if(badpixels!=NULL)
      delete[] badpixels;
   if(badpixelmethod!=NULL)
      delete[] badpixelmethod;
   if(bpmethod_descriptor!=NULL)
      delete[] bpmethod_descriptor;
}

//-------------------------------------------------------------------------
// Function to decode the bad pixel stringstream into an unsigned int array
// requires the reverse band mapping from the calibration to identify 
// non-used bands (i.e. map raw bands to calibration bands).
// This version is for Specim Calibrated Bad Pixel files
//-------------------------------------------------------------------------
void Hawk::DecodeSpecimBadPixels(std::map<int,int> revbandmap)
{
   //Bad pixel file is assumed to be in the following format
   //1st line description , ccdsamples, ccdbands
   //all other lines have 6 elements
   //id bsample bband rsample rband GON

   //id increase by 1 each time
   //bsample,bband,rsample,rband describe the bad pixel and the pixel to replace it with
   //we will just mask (set=0) rather than replace it
   //GON ....who knows?

   unsigned int id=0, bsample=0, bband=0, rsample=0, rband=0;
   unsigned int previd=0;
   std::string isgon="";
   int bnew=0;
   std::string tempstr;

   //Get the first line
   std::getline(badpixelstream,tempstr);

   //read in the first id
   badpixelstream>>id; 

   //Loop through to check that the stringstream is as expected
   //and get the number of bad pixels
   while(badpixelstream.good())
   {
      //get info from stringstream
      if(id!=previd+1)
      {
         //An error has occurred.
         throw std::string("An error occurred decoding bad pixel file ... id does not increase by 1 in file at ID: ")+ToString(id);     
      }

      badpixelstream>>bsample;
      badpixelstream>>bband;
      badpixelstream>>rsample;
      badpixelstream>>rband;
      badpixelstream>>isgon;
      if(isgon.compare("GON")!=0)
      {
         //An error has occurred...expected GON   
         throw std::string("An error occurred decoding bad pixel file ... 6th word is not 'GON' at ID: ")+ToString(id);     
      }
      previd=id;
      //get next id from next line
      badpixelstream>>id;
   }

   //Number of bad pixels
   nbadpixels=previd;

   //Create an array of required size to store bad pixel ccd row/col
   badpixels=new unsigned int[nbadpixels*2];
   previd=0;
   badpixelstream.clear();
   badpixelstream.seekg(0,std::ios::beg);
   std::getline(badpixelstream,tempstr);

   badpixelstream>>id;//get first id

   //Loop through the stream again, this time storing the bad pixels into the array
   while(badpixelstream.good())
   {
      //get info from stringstream
      if(id!=previd+1)
      {
         //An error has occurred. 
         throw std::string("An error occurred decoding bad pixel file ... id does not increase by 1 in file at ID: ")+ToString(id);     
      }
      badpixelstream>>bsample;
      badpixelstream>>bband;
      badpixelstream>>rsample;
      badpixelstream>>rband;
      badpixelstream>>isgon;
      if(isgon.compare("GON")!=0)
      {
         //An error has occurred...expected GON   
         throw std::string("An error occurred decoding bad pixel file ... 6th word is not 'GON' at ID: ")+ToString(id);     
      }
      previd=id;

      std::map<int,int>::iterator it;
      it=revbandmap.find(bband-1);

      if(it!=revbandmap.end())
      {
         //Convert bad pixel band number to raw band number
         bnew=it->second;
      }
      else
      {
         //flag the band with a not in use flag, used later on when masking bad pixel raw data
         bnew=bandnotinuse; 
      }

      badpixels[(id-1)*2+0]=bsample-1; //minus 1 to normalise between 0 and max-1
      badpixels[(id-1)*2+1]=bnew; //no minus 1 because this is taken care of in the band mapping   

      //read in next id   
      badpixelstream>>id; 
   }
}

//-------------------------------------------------------------------------
// Function to decode the bad pixel stringstream into an unsigned int array
// requires the reverse band mapping from the calibration to identify 
// non-used bands (i.e. map raw bands to calibration bands).
// This version is for ARSF Calibrated Bad Pixel files
//-------------------------------------------------------------------------
void Hawk::DecodeARSFBadPixels(std::map<int,int> revbandmap)
{
   //It is assumed that the ARSF calibrated bad pixel files are
   //referenced to 0 - i.e. the band and sample numbers run from 0->n-1
   int id=0, bsample=0, bband=0;
   int previd=-1;
   int bnew=0;
   std::string method;

   std::string tempstr;
   int nheaderlines=0;

   //Get the first line and the number of headerlines to skip
   std::getline(badpixelstream,tempstr);
   tempstr=TrimLeadingChars(tempstr,"headerlines=");
   nheaderlines=StringToINT(tempstr); 
   if(nheaderlines==0)
   {
      //This means something has gone wrong somewhere
       throw std::string("An error occurred decoding bad pixel file ... cannot get number of headerlines: ")+tempstr;          
   }  

   int method_counter=0;
   //We have already read 1 line so read nheaderlines - 1 more
   // and test if they are method descriptors
   for(int i=0;i<nheaderlines-1;i++)
   {
      std::getline(badpixelstream,tempstr);
      if(tempstr.find_first_of("method",0)==0)
      {
         method_counter++;
      }
   }

   //We have found method_counter number of descriptors
   Logger::Log("Will create method descriptors: "+ToString(method_counter));
   bpmethod_descriptor=new std::string[method_counter];   

   //read in the first id
   badpixelstream>>id;

   //Loop through to check that the stringstream is as expected
   //and get the number of bad pixels
   while(badpixelstream.good())
   {
      //get info from stringstream
      if(id!=previd+1)
      {
         //An error has occurred.
         throw std::string("An error occurred decoding bad pixel file ... id does not increase by 1 in file at ID: ")+ToString(id);     
      }

      badpixelstream>>bband;
      badpixelstream>>bsample;
      badpixelstream>>method;

      previd=id;

      //get next id from next line
      badpixelstream>>id;
   }

   //Number of bad pixels
   nbadpixels=previd+1; //+1 since previd referenced to 0 not 1

   //Create an array of required size to store bad pixel ccd row/col
   badpixels=new unsigned int[nbadpixels*2];
   badpixelmethod=new unsigned char[nbadpixels];

   //reset variables for reading in again
   previd=-1;
   badpixelstream.clear();
   badpixelstream.seekg(0,std::ios::beg);
   //We need to skip nheaderlines 
   //Actually - we want to read in and store the method descriptors to output later
   int meth=0;
   for(int i=0;i<nheaderlines;i++)
   {
      std::getline(badpixelstream,tempstr);
      //Test if line starts with method - these should be method descriptors
      if(tempstr.find("method")==0)
      {
         if(meth >= method_counter)
         {
            throw "More bad pixel detection method descriptions have been detected than initially found.";
         }
         Logger::Log(" "+tempstr);
         bpmethod_descriptor[meth]=tempstr;
         meth++;
      }
   }

   badpixelstream>>id;//get first id
   std::map<int,int>::iterator it;

   if(!badpixelstream.good())
      throw "Bad pixel stream is not good - but it should be so fix the code.";

   //Loop through the stream again, this time storing the bad pixels into the array
   while(badpixelstream.good())
   {
      //get info from stringstream
      if(id!=previd+1)
      {
         //An error has occurred. 
         throw std::string("An error occurred decoding bad pixel file ... id does not increase by 1 in file at ID: ")+ToString(id);     
      }
      badpixelstream>>bband;
      badpixelstream>>bsample;
      badpixelstream>>method;

      previd=id;
      it=revbandmap.find(bband);

      if(it!=revbandmap.end())
      {
         //Convert bad pixel band number to raw band number
         bnew=it->second;
      }
      else
      {
         //flag the band with a not in use flag, used later on when masking bad pixel raw data
         bnew=bandnotinuse; 
      }

      badpixels[(id)*2+0]=bsample; 
      badpixels[(id)*2+1]=bnew;   

      //Store the method as a char bit flag rather than actual method value
      //Method is now a comma delimited string of methods e.g. A,B,E
      //so loop through each method in string "method" and set all that apply
      char meth='z';
      badpixelmethod[id]=0;
      for(size_t m=0;m<GetNumberOfItemsFromString(method,",");m++)
      {
         std::string strItem=GetItemFromString(method,m,',').c_str();
         if(strItem.length() > 1)
            throw "Expected bad pixel method to be 1 char length, I got: "+strItem;
         else
            strItem.copy(&meth,1);
         switch(meth)
         {
         case 'A':   
            if((badpixelmethod[id]&A) == 0) //check if A is set, if NOT then set it.
               badpixelmethod[id]+=A;
            break;
         case 'B':
            if((badpixelmethod[id]&B) == 0)
               badpixelmethod[id]+=B;
            break;
         case 'C':
            if((badpixelmethod[id]&C) == 0)
               badpixelmethod[id]+=C;
            break;
         case 'D':
            if((badpixelmethod[id]&D) == 0)
               badpixelmethod[id]+=D;
            break;
         case 'E':
            if((badpixelmethod[id]&E) == 0)
               badpixelmethod[id]+=E;
            break;
         default:
            throw std::string("Unrecognised bad pixel detection method in bad pixel file. Expected one of A,B,C,D,E but got: ")+method;
         }
      }

      //read in next id   
      badpixelstream>>id; 
   }
}

//-------------------------------------------------------------------------
// Function wrapper to decode the bad pixel stringstream into an unsigned 
// int array. Calls the function method based on the stringstream info
//-------------------------------------------------------------------------   
void Hawk::DecodeBadPixel(std::map<int,int> revbandmap)
{
   //The first line of ARSF and Specim bad pixel files differ
   //ARSF: headerlines = x
   //Specim: NERC_Hawk_BPR_NUC2_GON0.txt  320 256 (or something similar)
   
   //Get the first line
   std::string tempstr;
   std::getline(badpixelstream,tempstr);
   
   //Seek back to the start of the stream
   badpixelstream.seekg(0,std::ios::beg);

   //Test for suspected file type
   if(tempstr.find("headerlines")!=std::string::npos)
   {
      //ARSF Calibrated file detected
      Logger::Log("Detected ARSF calibrated bad pixel file.");
      arsfbadpixelfiletype=true;
      DecodeARSFBadPixels(revbandmap);
   }
   else if(tempstr.find("320 256")!=std::string::npos)
   {
      //Specim calibrated file detected
      Logger::Log("Detected Specim calibrated bad pixel file.");
      arsfbadpixelfiletype=false;
      DecodeSpecimBadPixels(revbandmap);
   }
   else
   {
      //Unable to detect
      Logger::Log(tempstr);
      throw "Unable to detect bad pixel type - failed ARSF and Specim file type test.";
   }
}

