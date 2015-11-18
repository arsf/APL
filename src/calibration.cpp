//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#include <cmath>
#include <typeinfo>

#include "calibration.h"
//-------------------------------------------------------------------------
// Constructor for data object taking size of array 
// which would usually be (number of bands * number of samples)
//-------------------------------------------------------------------------
Data::Data(unsigned long size)
{
   arraysize=size;
   fodis=NULL;
   mask=NULL;
   badpixmethod=NULL;
   image=new double[size];
   avdark=NULL;
   gains=NULL;

   for(unsigned long i=0;i<size;i++)
   {
      image[i]=0;
   }
}

//-------------------------------------------------------------------------
// destructor for data object
//-------------------------------------------------------------------------
Data::~Data()
{
   if(fodis!=NULL)
      delete[] fodis;
   if(mask!=NULL)
      delete[] mask;
   if(badpixmethod!=NULL)
      delete[] badpixmethod;
   if(image!=NULL)
      delete[] image;
   if(avdark!=NULL)
      delete[] avdark;
   if(gains!=NULL)
      delete[] gains;
}

//-------------------------------------------------------------------------
// Initialise the dark frames array
//-------------------------------------------------------------------------
void Data::InitialiseDarkFrames()
{
   if(avdark==NULL) 
   {
      avdark=new double[arraysize];
      for(unsigned long i=0;i<arraysize;i++)
         avdark[i]=0;
   }
}

//-------------------------------------------------------------------------
// Initialise the mask array
//-------------------------------------------------------------------------
void Data::InitialiseMask()
{
   if(mask==NULL) 
   {
      mask=new unsigned char[arraysize];
      for(unsigned long i=0;i<arraysize;i++)
         mask[i]=0;
   }
}

//-------------------------------------------------------------------------
// Initialise the bad pixel method array
//-------------------------------------------------------------------------
void Data::InitialiseBadPixMethod()
{
   if(badpixmethod==NULL) 
   {
      badpixmethod=new unsigned char[arraysize];
      for(unsigned long i=0;i<arraysize;i++)
         badpixmethod[i]=0;
   }
}

//-------------------------------------------------------------------------
// Initialise the fodis array
//-------------------------------------------------------------------------
void Data::InitialiseFodis()
{
   if(fodis==NULL) 
   {
      fodis=new double[arraysize];
      for(unsigned long i=0;i<arraysize;i++)
         fodis[i]=0;
   }
}

//-------------------------------------------------------------------------
// Initialise the gains array
//-------------------------------------------------------------------------
void Data::InitialiseGains()
{
   if(gains==NULL) 
   {
      gains=new double[arraysize];
      for(unsigned long i=0;i<arraysize;i++)
         gains[i]=0;
   }
}


//-------------------------------------------------------------------------
//Function to transform the ordering of data in arrays
//-------------------------------------------------------------------------
void Data::TransformArrays(const unsigned int bands, const unsigned int samples,TransformArray order)
{
   if(samples*bands != arraysize)
      throw "Cannot transform data arrays as given number of bands/samples does not agree with array size.";

   if(fodis!=NULL)
   {
      if(order == BAND) 
         FlipBandData(fodis,bands,samples);
      else if(order == SAMPLE)
         FlipSampleData(fodis,bands,samples);
   }
   if(mask!=NULL)
   {
      if(order == BAND) 
         FlipBandData(mask,bands,samples);
      else if(order == SAMPLE)
         FlipSampleData(mask,bands,samples);
   }
   if(badpixmethod!=NULL)
   {
      if(order == BAND) 
         FlipBandData(badpixmethod,bands,samples);
      else if(order == SAMPLE)
         FlipSampleData(badpixmethod,bands,samples);
   }
   if(image!=NULL)
   {
      if(order == BAND) 
         FlipBandData(image,bands,samples);
      else if(order == SAMPLE)
         FlipSampleData(image,bands,samples);
   }
}

//-------------------------------------------------------------------------
// Calibration constructor 
//-------------------------------------------------------------------------
Calibration::Calibration(Specim* sen,std::string calFile)
{
   //Assign variables
   sensor=sen;
   thissubsensor=0;
   calibrationFilenamePrefix=calFile;
   badpixels=NULL;

   if(CheckSensorID(FENIX,this->sensor->SensorID()))
   {      
      //Get size of 1 line (all bands) of data
      unsigned int vnb=dynamic_cast<Fenix*>(sensor)->NumBandsVnir();
      unsigned int snb=dynamic_cast<Fenix*>(sensor)->NumBandsSwir();
      unsigned int ns=sensor->NumSamples();
      unsigned long sizevnir=vnb*ns;
      unsigned long sizeswir=snb*ns;
      //Create new data objects of correct size
      numofsensordata=2;
      sensordata=new Data*[numofsensordata];
      sensordata[0]=new Data(sizevnir);
      sensordata[1]=new Data(sizeswir);
      sensorbandmap=new std::map<int,int>[2];
      sensorrevbandmap=new std::map<int,int>[2];
      //Set to subsensor 0 by default
      data=sensordata[0];
      bandmap=&(sensorbandmap[0]);
      revbandmap=&(sensorrevbandmap[0]);
   }
   else
   {
      //Get size of 1 line (all bands) of data
      unsigned int nb=sensor->NumBands();
      unsigned int ns=sensor->NumSamples();
      unsigned long size=nb*ns;
      //Create new data object of correct size
      numofsensordata=1;
      sensordata=new Data*[numofsensordata];
      sensordata[0]=new Data(size);
      sensorbandmap=new std::map<int,int>[1];
      sensorrevbandmap=new std::map<int,int>[1];
      //Set to subsensor 0 (as it is the only one)
      data=sensordata[0];
      bandmap=sensorbandmap;
      revbandmap=sensorrevbandmap;
   }

   //If mask is requested - initialise it - always do this?
   InitialiseMask();
}

void Calibration::ChangeSubSensor(unsigned int sensorindex)
{
   //If only one subsensor then return
   if(numofsensordata==1)
      return;
   //Else update to the alternative
   if(sensorindex < numofsensordata)
   {
      if(sensorindex==0)
         dynamic_cast<Fenix*>(sensor)->SetUpFenixFor(VNIR);
      else if (sensorindex==1)
         dynamic_cast<Fenix*>(sensor)->SetUpFenixFor(SWIR);
      else
         throw "Calibration is only set up for Fenix with 2 subsensors";

      data=sensordata[sensorindex];
      subsensorlowerband=dynamic_cast<Fenix*>(sensor)->SubSenLowerBand();
      bandmap=&(sensorbandmap[sensorindex]);
      revbandmap=&(sensorrevbandmap[sensorindex]);
      thissubsensor=sensorindex;
   }
   else
      throw "Sensorindex in changesubsensor is greater than number of sensors: "+ToString(sensorindex);

}

//-------------------------------------------------------------------------
// Mask initialiser
//-------------------------------------------------------------------------
void Calibration::InitialiseMask()
{
   //If mask is requested - initialise it
   for(unsigned int i=0;i<numofsensordata;i++)
      sensordata[i]->InitialiseMask();
}

void Calibration::InitialiseBadPixMethod()
{
   //If mask is requested - initialise it
   for(unsigned int i=0;i<numofsensordata;i++)
      sensordata[i]->InitialiseBadPixMethod();
}

//-------------------------------------------------------------------------
// Fodis initialiser
//-------------------------------------------------------------------------
void Calibration::InitialiseFodis()
{
   //If fodis is requested - initialise it
   for(unsigned int i=0;i<numofsensordata;i++)
      sensordata[i]->InitialiseFodis();
}

//-------------------------------------------------------------------------
// Calibration destructor
//-------------------------------------------------------------------------
Calibration::~Calibration()
{
   if(sensordata!=NULL)
   {
      for(unsigned int s=0;s<numofsensordata;s++)
         delete sensordata[s];
   }
   delete sensordata;
   delete[] sensorbandmap;
   delete[] sensorrevbandmap;

   if(badpixels!=NULL)
      delete[] badpixels;

   if(data!=NULL)
      data=NULL;
   if(sensor!=NULL)
      sensor=NULL;
}

//-------------------------------------------------------------------------
// Function to test the calibration file is suitable for this sensor data
//-------------------------------------------------------------------------
void Calibration::TestCalfile()
{
   if(calibrationFilenamePrefix.compare("")==0)
      return;

   std::string calibrationFilename=calibrationFilenamePrefix + ".cal";
   Logger::Log("Reading calibration file: "+calibrationFilename);
   Logger::Warning("Note that integration time of calibration file is ignored - assume it relates to a value of 1.0.");
   BinFile calin(calibrationFilename);

   //Get dimensions of the file - these are likely eagle: 1024,1,1024 or hawk: 320,1,256
   unsigned int nsamps=StringToUINT(calin.FromHeader("samples"));
   unsigned int nlines=StringToUINT(calin.FromHeader("lines"));
   unsigned int nbands=StringToUINT(calin.FromHeader("bands"));

   Logger::Log("Calibration file contains:  "+ToString(nsamps)+" samples, "+ToString(nlines)+" lines and "+ToString(nbands)+" bands.");      

   if(nlines!=1)
   {
      //Error expects cal file to be 1 line of data
      throw "Error. Calibration file must have 1 line of data only.";
   }   

   if(TrimWhitespace(calin.FromHeader("data type"))!="4")//float data
   {
      //Error can only handle floating point cal files
      throw "Error. Calibration file specified is not floating point data. Cannot handle this data type.";
   }

   if(StringToINT(TrimWhitespace(calin.FromHeader("sensorid"))) != sensor->SensorID())
   {
      throw "Sensor ID for calibration file disagrees with Sensor id for raw file.\nCalibration: "+calin.FromHeader("sensorid")+" Raw: "+ToString(sensor->SensorID());
   }

   unsigned int icount=0;
   //Need to get the number of wavelength entries
   while(calin.FromHeader("Wavelength",icount).compare("")!=0)
   {
      icount++;
   }
   if(icount==0)
   {
      //No entries....error
   }

   //set the number of wavelength entries
   unsigned int numwl_cal=icount;

   //Check if this equals the number of bands in the file
   if(numwl_cal != nbands)
   {
      throw "Number of bands in the calibration file does not agree with the number of wavelengths in the calibration file.";
   }

   //create array for calibration wavelengths (per band)
   float* wl_cal=new float[numwl_cal];

   for(icount=0;icount<numwl_cal;icount++)
   {
      wl_cal[icount]=static_cast<float>(StringToDouble(calin.FromHeader("Wavelength",icount))); //convert to double (cast as float)
   }

   for(unsigned int subsensor=0;subsensor<numofsensordata;subsensor++)
   {
      Logger::Log("Checking calibration wavelengths for subsensor: "+ToString(subsensor));
      ChangeSubSensor(subsensor);
      CheckCalWavelengths(wl_cal,numwl_cal);
   }
   
   //Close the Bil reader
   calin.Close();
   //delete the wavelength array
   delete[] wl_cal;
}

//-------------------------------------------------------------------------
// Function to get the ratio of binning values raw : cal
//-------------------------------------------------------------------------
unsigned int Calibration::GetBinningRatio(std::string bintype)
{
   std::string calibrationFilename=calibrationFilenamePrefix+".cal";
   BinFile calfile(calibrationFilename);

   //Need to also take into account if calibration file has been binned (such as with FENIX)
   unsigned int calbinning=0;
   float calratio=0;
   if(bintype.compare("spectral")==0)
   {
      if(data==sensordata[0])
         calbinning=StringToUINT(TrimWhitespace(calfile.FromHeader("binning",0).c_str()));
      else if(data==sensordata[1])
         calbinning=StringToUINT(TrimWhitespace(calfile.FromHeader("binning2",0).c_str()));
      //Instead of just using the raw data binning values - we use a ratio of raw/cal
      calratio=(float)sensor->SpectralBinning() / calbinning;
   }
   else if(bintype.compare("spatial")==0)
   {
      if(data==sensordata[0])
         calbinning=StringToUINT(TrimWhitespace(calfile.FromHeader("binning",1).c_str()));
      else if(data==sensordata[1])
         calbinning=StringToUINT(TrimWhitespace(calfile.FromHeader("binning2",1).c_str()));
      //Instead of just using the raw data binning values - we use a ratio of raw/cal
      calratio=(float)sensor->SpatialBinning() / calbinning;
   }
   else
   {
      throw "Unrecognised binning type in GetBinningRatio: "+bintype;
   }

   //The calratio should always be an integer if the two files are compatible. If less than 1 then it means
   //that the cal file has been binned too far to use for this raw data

   //Check there are no less than 1 ratios
   if(calratio < 1)
   {
      throw "Calibration file has a binning higher than the raw data file - the calibration file is therefore not suitable for use.";
   }
   //Check that the ratios are integers
   else if(floor(calratio)!=calratio)
   {
      throw "Ratio of raw:calibration binning is not an integer.\nThis suggests some odd binnings and that the calibration file may be unsuitable for the raw data.";
   }
   //As these are integer values we'll store them as ints 
   return static_cast<unsigned int>(calratio);  
}

//-------------------------------------------------------------------------
//Function to check the calibration file and see if it is suitable for the raw data
//-------------------------------------------------------------------------
void Calibration::CheckCalWavelengths(float* const wl_cal, const unsigned int numwl_cal)
{
   const double searchbound=0.006; // search +/- this on wavelength when looking for matches (since exact match is very near impossible)

   //Check the centre wavelengths to see if they agree for the required number of bands
   //The returned value is the number of agreements between raw and cal wavelengths
   unsigned int numwl_raw=0;
   //Need to get the number of wavelength entries
   std::string rawwaves=TrimWhitespace(sensor->bin->GetFromFile("Wavelength"));
   numwl_raw=GetNumberOfItemsFromString(rawwaves,";");
   if(numwl_raw==0)
   {
      //No entries....error
      throw "An error has occurred...there doesn't appear to be any wavelengths in the raw hdr file.";
   }

   unsigned int numagree=0;//counter for number of agreeing wavelengths

   //The below should give an integer number so test it - we use the ratio of spectral binnings incase the cal file is binned too
   unsigned int specbinratio=GetBinningRatio("spectral");
   float numbinnedcal=numwl_cal / specbinratio;

   if(numbinnedcal - int(numbinnedcal) != 0)
      throw "Binned calibration data is not an integer value. This is not good. Num wl: "+ToString(numwl_cal)+" spectral binning: "+ToString(sensor->SpectralBinning());

   if(numwl_raw > numbinnedcal)
      throw "Number of raw file wavelengths is larger than number of bands in binned calibration file: "+ToString(numwl_raw)+" vs "+ToString(numbinnedcal);

   float* wl_raw=NULL;
   wl_raw=new float[numwl_raw];//create array for raw wavelengths (per band)
   for(unsigned int icount=0;icount<numwl_raw;icount++)
   {
      wl_raw[icount]=static_cast<float>(StringToDouble(GetItemFromString(rawwaves,icount,';'))); //convert to double (cast as float)
   }

   //Create a new array to hold the binned wavelengths
   float* binnedcal=new float[int(numbinnedcal)];

   //Spectral binning needs to adjust wl_cal (the wavelengths from the calibration file hdr)
   for(unsigned int j=0;j<numbinnedcal;j++)
   {
      binnedcal[j]=0;
      for(unsigned int i=0;i<specbinratio;i++)
      {
         binnedcal[j]=binnedcal[j]+wl_cal[specbinratio*j + i];
      }
      binnedcal[j]=binnedcal[j] / specbinratio;
   }

   //We now have the binned wavelength values...assuming identical fwhm (bandwidths actually)
   //Logger::Log("WARNING: This assumes equal fwhms. If they vary then this calibration is wrong."); 
  
   unsigned int iterator=0;//iterating variable for loop below
   //We now need to go through all the wavelengths of the raw image data and check
   //that each one appears in the list of (binned) calibration wavelengths
   for(unsigned int w=0;w<numwl_raw;w++)
   {
      while(iterator<numbinnedcal)
      {
         //Check if raw wavelength is within +/- searchbound of the binned calibration wavelengths since exact floating point is dangerous
         if((wl_raw[w]<=binnedcal[iterator]+ searchbound)&&(wl_raw[w]>=binnedcal[iterator] - searchbound))
         {
            (*bandmap)[w]=iterator; //maps raw band 'w' to binned calibration band 'iterator'
            (*revbandmap)[iterator]=w;//reverse mapping from binned calibration band to raw band

            numagree++;//increase the number of raw/cal agreements
            break; //break out of the while loop
         }
         else
         {
            iterator++;//go onto the next cal wavelength to see if that agrees with the raw wavelength
         }
      }
      iterator=0; //reset iterator
   }   

   //Numagree should be the same as the number of bands in the raw file if the correct calibration file is used.
   if(numagree != sensor->NumBands())
   {
      throw "A number of bands in the raw file disagree with the calibration file wavelengths. Number that agree: "+ToString(numagree)+" Number of bands: "+ToString(sensor->NumBands())
               +"\nThis probably means that the calibration file and raw file are not compatible. Maybe the wrong bandset / config file has been used at data collection.";
   }
   else
   {
      //Output this for infomation
      Logger::Log("There are "+ToString(numagree)+" bands in (binned) calibration file whose centre wavelengths agree with the raw image.");
   }
   
   //Remove the temporary arrays
   delete[] binnedcal;
   delete[] wl_raw;

}

//-------------------------------------------------------------------------
//Function to read line of raw
//-------------------------------------------------------------------------
void Calibration::ReadLineOfRaw(unsigned int line)
{
   sensor->bin->ReadlineToDoubles(data->Image(),line);
}

//-------------------------------------------------------------------------
//Function to initialise the dark frame array for subtraction
//-------------------------------------------------------------------------
void Calibration::InitialiseDarkFrames(std::string darkfile)
{
   if(data->AverageDark()!=NULL)
      throw "Average dark frame array already initialised - cannot do it twice.";

   if(darkfile.compare("")!=0)
   {
      BinFile dark(darkfile);

      unsigned int nsamples=StringToUINT(dark.FromHeader("samples"));
      unsigned int nbands=StringToUINT(dark.FromHeader("bands"));
      //We should sum all subsensor array sizes before checking against dark file array sizes
      uint64_t sumarraysize=0;
      for(unsigned int i=0;i<numofsensordata;i++)
      {
         ChangeSubSensor(i);
         sumarraysize += data->ArraySize();
      }

      if(sumarraysize != nsamples*nbands)
      {
         throw "Number of samples or bands of dark file disagree's with raw file. Cannot use this dark file and raw file together.";
      }
   }
   else if(sensor->GetNumDarkFrames()==0)
   {
      //No dark frames in file and no external file given - force exit
      throw "No dark frames found in file and no external dark file given - will not proceed. Use -NODARK on commandline to override.";
   }

   Logger::Log("\nStarting dark frame analysis...");
   //We want to do this for all subsensors here
   for(unsigned int i=0;i<numofsensordata;i++)
   {
      ChangeSubSensor(i);
      data->InitialiseDarkFrames();

      double* mean=NULL;
      double* stdev=NULL;

      mean=new double[data->ArraySize()];
      stdev=new double[data->ArraySize()];

      for(unsigned int looper=0;looper<data->ArraySize();looper++)
      {
         mean[looper]=stdev[looper]=0;
      }

      sensor->AverageAllDarkFrames(mean,darkfile,data->ArraySize(),i);
      sensor->DarkFramesStdDeviation(stdev,mean,darkfile,data->ArraySize(),i);
      sensor->AverageRefinedDarkFrames(data->AverageDark(),stdev,mean,darkfile,data->ArraySize(),i);

      //Scale the darkframes if scalar is not = 1 - actually will exit now
      if(sensor->DarkScalar() != 1)
      {
         throw "Dark frames have a different integration time to 'light' frames. Please choose appropriate dark frames or scale your dark frames. Scalar: "+ToString(sensor->DarkScalar());
      }

      Logger::Log("Finished dark frame analysis for sensor index "+ToString(i)+" ...\n");

      delete[] mean;
      mean=NULL;
      delete[] stdev;
      stdev=NULL;
   }
}

//-------------------------------------------------------------------------
//Function to remove dark frames
//-------------------------------------------------------------------------
void Calibration::RemoveDarkFrames()
{

   //Test if data avdark and image arrays exist
   if(data->AverageDark()==NULL)
      throw "The data -> average dark value array has not been initialised.";

   if(data->Image()==NULL)
      throw "The data -> image array has not been initialised.";

   //For each element of the data buffer
   for(unsigned int ele=0;ele<(sensor->NumBands()*sensor->NumSamples());ele++)
   {  
      if((data->Image()[ele]!=0)&&(data->Image()[ele]<sensor->CalibratedMax()))//Data is not under flown or above max calibrated value
      {
         //If the average dark value is larger than the maximum allowed this is bad
         if(data->AverageDark()[ele]>sensor->RawMax())
         {
            throw "ERROR in RemoveDarkFrames. Average dark value is greater than raw maximum: "+ToString(data->AverageDark()[ele])+" > "+ToString(sensor->RawMax())+" element number "+ToString(ele);
         }
         else if(data->AverageDark()[ele]==sensor->RawMax())
         {
            //Send a message that the average dark for this pixel is the maximum value - only need to output it once though
            Logger::WarnOnce("Average dark value for pixel "+ToString(ele)+" is the maximum raw value: "+ToString(sensor->RawMax()));
         }

         //If the data after subtraction is less than (or equals) 0 then set as underflow
         if( (data->Image()[ele]-data->AverageDark()[ele]) <= 0) 
         {
            //DEBUGPRINT("Dark value: "<<data->AverageDark()[p]<<" Image value: "<<data->Image()[ele]<<" element: "<<ele);
            data->Image()[ele]=0; // assign as an underflow
            AssignMaskValue(ele,sensor->UnderFlow);   
         }
         else 
            data->Image()[ele]=((data->Image()[ele]- data->AverageDark()[ele])); //subtract dark values to normalise to 0
      }
   }
}



//-------------------------------------------------------------------------
//Function to apply a frame smear correction to (eagle raw data - dark frames)
//-------------------------------------------------------------------------
bool Calibration::SmearCorrect()
{
   //Frame smear correction scalar - frame transfer time / integration time
   //Check we have Eagle data
   if(!CheckSensorID(EAGLE,this->sensor->SensorID()))
   {
      Logger::Warning("Cannot apply smear correction to this sensor data - it claims not to be from an Eagle sensor. Skipping in future.");
      return false;
   }
   else
   {
      //Logger::Log("Eagle data has been identified in Smear Correction.");
   }

   //Its OK to cast the sensor type to Eagle from here on - dynamic casting allows further safety
   Eagle* eagle=dynamic_cast<Eagle*>(this->sensor);

   //Calculate smear scalar
   double fsc=(eagle->FrameTransferTime()/eagle->IntegrationTime())*(eagle->SpectralBinning());
   
   //calculate per pixel correction and apply
   //The fsc is multiplied by sum of (spectrally) previous bands. This is removed from current band.
   //Use the corrected bands not the uncorrected ones to correct the current band.

   double bandsum=0; //value to hold the sum of the previous bands
   unsigned long current_band_sample=0 ,previous_band_sample=0;

   //For each sample of the frame
   for(unsigned int s=0;s<eagle->NumSamples();s++)
   {         
      bandsum=0;//reset the bandsum variable to zero before entering band loop
      //First band has no correction applied
      for(unsigned int b=1;b<eagle->NumBands();b++)
      {
         //PROBLEM...DOES THIS NEED TO RUN FROM THE 1ST BAND OF CCD (IE 1024 BANDS)
         //OR FROM THE FIRST BAND OF THE DATA. WILL ASSUME DATA FOR THE MOMENT SINCE NOTHING WE CAN DO IF REQUIRED FOR WHOLE CCD 

         //Assign sample indices to variables to make the following clearer code
         current_band_sample=b*eagle->NumSamples()+s;
         previous_band_sample=(b-1)*eagle->NumSamples()+s;

         //for each previous band, sum up the total - use a rolling total adding previous band each iteration
         bandsum=bandsum+data->Image()[previous_band_sample];

         //Now correct the current band 
         data->Image()[current_band_sample]=((data->Image()[current_band_sample] - fsc*(bandsum)));
         if(data->Image()[current_band_sample]<0)
         {
            //Assign as underflow
            data->Image()[current_band_sample]=0;  
            AssignMaskValue(current_band_sample,sensor->UnderFlow);   
         }
      }
   }

   return true;
}

//-------------------------------------------------------------------------
// Function to apply the gains from the calibration file to the image data
// Also handles some checking of the gains, binning etc
//-------------------------------------------------------------------------
void Calibration::ApplyGains()
{
   //Do some checks that everything is in place
   if(this->calibrationFilenamePrefix.compare("")==0)
   {
      throw "Cannot apply gains if calibration file is not set.";
   }
   
   //Read in the gains from file
   //Bin the gains to the correct size (checking cal vs Raw)
   if(data->Gains()==NULL)
   {
      data->InitialiseGains();
      ReadBinAndTrimGains(data->Gains());
   }

   //Apply the gains
   if((sensor->IntegrationTime() == 0))
   {
      throw "Error integration time is 0 in ApplyGains()."; 
   }

   // we need to scale the radiances by this value
   double radmultiplier= sensor->RadianceScalar() / sensor->IntegrationTime(); 

   //temporary double to test value
   double tempval=0; 

   //For each element of the data buffer we want to apply the calibration
   for(unsigned int ele=0;ele<sensor->NumBands()*sensor->NumSamples();ele++)
   {  
      if((data->Image()[ele]!=0)&&(data->Image()[ele]!=sensor->CalibratedMax()))
      {
         tempval=(data->Image()[ele]*data->Gains()[ele]*radmultiplier);//scale by the calibration multipliers
         if((tempval) >= sensor->CalibratedMax())
         {
            data->Image()[ele]=sensor->CalibratedMax(); // assign as an overflow
            AssignMaskValue(ele,sensor->OverFlow);   
         }
         else
         {
            data->Image()[ele]=tempval;
         }
      }  
   }
}

//-------------------------------------------------------------------------
// Function to read in the calibration gains, bin them and trim unused bands
//-------------------------------------------------------------------------
void Calibration::ReadBinAndTrimGains(double* const trimmedcal)
{
   std::string calibrationFilename=calibrationFilenamePrefix+".cal";
   BinFile calfile(calibrationFilename);
   unsigned int nsamps=StringToUINT(calfile.FromHeader("samples"));
   unsigned int nlines=StringToUINT(calfile.FromHeader("lines"));
   unsigned int nbands=StringToUINT(calfile.FromHeader("bands"));

   //Use the ratio of binnings rather than just the raw binning value
   unsigned int specbinratio=GetBinningRatio("spectral");
   unsigned int spatbinratio=GetBinningRatio("spatial");

   //Lets spectrally and spatially bin the data here.
   unsigned int numbinnedband=nbands / specbinratio; //the number of binned calibration bands
   unsigned int numbinnedsamps=nsamps / spatbinratio; //the number of binned calibration samples

   //Check that the number of samples of the raw image agree with the binned cal samples. If not - exit
   if(numbinnedsamps!=sensor->NumSamples())
      throw "Number of binned calibration samples is not equal to number of raw image samples";   

   if(nlines!=1)
      throw std::string("Calibration file should only have one line of data - got: ")+ToString(nlines);

   //Read in the calibration gains into a temporary array - note only a 1 line bil file
   double* tmpgains=new double[nsamps*nbands];
   calfile.ReadlineToDoubles(tmpgains,0);

   //Sample count is a counting variable to keep track of the sample number of the binned data
   int samplecount=0;
   unsigned long binned_index=0;
   double* binnedgains=new double[numbinnedband*numbinnedsamps];

   Logger::Log("Will bin calibration file so that it has "+ToString(numbinnedband)+" bands and "+ToString(numbinnedsamps)+" samples.");

   //For each band of the binned data
   for(unsigned int j=0;j<numbinnedband;j++)
   {
      samplecount=0;// counter for sample number of new array
      for(unsigned int s=0;s<nsamps;s=s+spatbinratio)//For each sample
      {
         //Assign an index variable for clearer code
         binned_index=j*numbinnedsamps + samplecount;
         binnedgains[binned_index]=0;//initialise to zero
         for(unsigned int i=0;i<specbinratio;i++) // loop through the bands that will be binned for this cell
         {
            for(unsigned int p=0;p<spatbinratio;p++) //loop throuhg the samples that will be binned for this cell
            {
               //Sum up the values that will be binned into one of the new array cells
               binnedgains[binned_index]=binnedgains[binned_index]+tmpgains[(j*specbinratio+i)*nsamps + (s+p)];
            }
         }

         if(CheckSensorID(FENIX,this->sensor->SensorID()))
         {
            //Fenix handles binning differently to Eagle + hawk
            //SPECTRAL is averaged - only divide by 1 factor of specbinratio (for the binned gains)
            //SPATIAL is summed - divide by 2 factors of specbinratio (1 for the binned gains and 1 for the binned raw data)
            binnedgains[binned_index]=binnedgains[binned_index] / ( (specbinratio * spatbinratio * spatbinratio) );
         }
         else if((CheckSensorID(EAGLE,this->sensor->SensorID()))||(CheckSensorID(EAGLE,this->sensor->SensorID())))
         {
            //Divide by spectralbin ^ 2 since the actual binned raw data is the sum of the binned bands not the mean of the binned bands.
            //Similarly for spatial binning too. i.e. both SPETCRAL and SPATIAL are summed.
            binnedgains[binned_index]=binnedgains[binned_index] / ( (specbinratio*specbinratio) * (spatbinratio*spatbinratio) );
         }
         else
         {
            throw std::string("Unrecognised sensor in calibration gains binning. Sensor id: ")+ToString(this->sensor->SensorID());
         }

         samplecount++;//increase onto next sample of array
      }
   }

   //Now that we have a calibration file that is binned to the same binning as the raw data
   //We need to subset the calibration file to get the same dimensions as the raw image
   //For each band in the raw image
   Logger::Log("Will trim calibration file so that it has "+ToString(sensor->NumBands())+" bands and "+ToString(sensor->NumSamples())+" samples.");
   for(unsigned int b=0;b<sensor->NumBands();b++)
   {
      for(unsigned int s=0;s<sensor->NumSamples();s++) //for each sample 
      {
         trimmedcal[b*sensor->NumSamples() + s] = binnedgains[(*bandmap)[b]*sensor->NumSamples() + s];
      }
   }   


   delete[] tmpgains;
   delete[] binnedgains;
   calfile.Close();
}

//-------------------------------------------------------------------------
// Function to apply pixel flags to the data
//-------------------------------------------------------------------------
void Calibration::FlagPixels()
{
   unsigned int p=0;

   //For each element of this scan line of data we need to flag erroneous pixels
   for(unsigned int band=0;band<sensor->NumBands();band++)
   {
      for(unsigned int sample=0;sample<sensor->NumSamples();sample++)
      {
         p=band*sensor->NumSamples() + sample;
         //This is a pre-calibration test of the raw data - will not work for post calibrated data
//         if((data->Image()[p] > sensor->RawMax())&&(data->Image()[p]!=sensor->CalibratedMax())&&((band!=0)&&(sample!=0)))
//         {
//            std::string error=std::string("Image contains data greater than the raw maximum. This should not happen...but has...Check the raw data at band,sample:")
//                              +ToString(band)+" "+ToString(sample)+" value: "+ToString(data->Image()[p]);
//            //throw error;
//            Logger::Warning(error);
//         }
         if((band==0)&&((sample==0)||(sample==1)))
         {
            //this is the frame counting ccd pixel or the 0 pixel next to it
            //  ... unless it is a fenix subsensor 1(SWIR). So check if the lower band 
            // limit is 0 - meaning that band=0 truely means band 0 of the raw file
            if(sensor->LowerBandLimit()==0)
            {
               data->Image()[p]=0; 
               AssignMaskValue(p,sensor->Badpixel);   
            }
         }
         else if(data->Image()[p] == sensor->RawMax())
         {
            //Overflowed data
            AssignMaskValue(p,sensor->OverFlow);   
            //also need to flag each corresponding sample for each lower band for the Eagle sensor
            if(CheckSensorID(EAGLE,sensor->SensorID()))
            {
               for(unsigned int b=band+1;b<sensor->NumBands();b++)
               {
                  AssignMaskValue(b*sensor->NumSamples() + sample,sensor->SmearAffected);
               }
            }        
         }
         else if(data->AverageDark()!=NULL) //dark frames have been initialised
         {
            if(data->Image()[p] <= data->AverageDark()[p])
            {
               //Data have underflown
               data->Image()[p]=0;
               AssignMaskValue(p,sensor->UnderFlow);
            }
         }
         else if((data->Image()[p] != sensor->CalibratedMax())&&(data->Image()[p] != sensor->RawMax()))
         {
            AssignMaskValue(p,sensor->Good);   
         }
      }
   }

   //Flag "bad pixels" from the bad pixel file - if no bad pixels this loop is essentially skipped
   if((badpixels==NULL)&&(CheckSensorID(EAGLE,sensor->SensorID())))
   {
      //This is OK - Eagle doesn't have bad pixels
   }
   else if((badpixels==NULL)&&(calibrationFilenamePrefix.compare("")==0))
   {
      //No calibration file has been given - hence no bad pixel file
      Logger::WarnOnce("As no calibration file has been given - will not be able to mask bad pixels.");
   }
   else if((badpixels==NULL)&&(!CheckSensorID(EAGLE,sensor->SensorID())))
   {
      //Other sensors use bad pixel files - throw exception
      throw "Bad pixels array has not been declared. All sensors except Eagle should use a bad pixel file.";
   }
   else //badpixels is not null
   {
      for(unsigned int bp=0;bp<badpixels[thissubsensor].NumBadPixels();bp++) //loop through bad pixel array
      {
         //Test to make sure band is in image range - since fenix band range is synthetic (i.e 0 -> x, 0 -> y and not 0 -> x+y)
         //Only mask if the band is in use in the raw data
         if(badpixels[thissubsensor].GetBadPixels()[2*bp+1]!=badpixels[thissubsensor].BandNotInUse())
         {
            AssignMaskValue((badpixels[thissubsensor].GetBadPixels()[2*bp+1])*sensor->NumSamples()
                            + badpixels[thissubsensor].GetBadPixels()[2*bp],sensor->Badpixel);   
           //Assign the detection method if ARSF bad pixel file
            if((badpixels[thissubsensor].arsfbadpixelfiletype == true)&&(data->BadPixMethod()!=NULL))
            {
               data->BadPixMethod()[badpixels[thissubsensor].GetBadPixels()[2*bp+1]*sensor->NumSamples() 
                                    + badpixels[thissubsensor].GetBadPixels()[2*bp]] = badpixels[thissubsensor].BadPixelMethod()[bp];
            }
         }
      } 
   }

   //Flag any quality control failing pixels here - these have been externally identified and supplied via a file
   if(sensor->qcfailures.size()!=0)
   {
      //If there are some qc failures - loop through them all and flag
      for(std::vector<Pair>::iterator it=sensor->qcfailures.begin();it!=sensor->qcfailures.end();it++)
      {
         AssignMaskValue((*it).band*sensor->NumSamples() + (*it).sample,sensor->QCFailure);
      }
   }
}

//-------------------------------------------------------------------------
// Function to read in the bad pixel file for sensor data
//-------------------------------------------------------------------------
bool Calibration::ReadBadPixelFile()
{
   //Currently we only do this for hawk but could potentially want to run
   //on other sensors too?
   if(CheckSensorID(EAGLE,sensor->SensorID()))
   {
      Logger::Warning("Currently don't apply bad pixels to Eagle data - skipping reading bad pixel file.");
      return false;
   }

   if(calibrationFilenamePrefix.compare("")==0)
   {
      //No calibration file has been given - also means we can't flag bad pixels
      Logger::Warning("As no calibration file was given - cannot identify or flag badpixels - skipping reading bad pixel file.");
      return false;
   }

   std::string badfile=calibrationFilenamePrefix+".bad";
   badpixels=new BadPixels[numofsensordata];
   for(unsigned int i=0;i<numofsensordata;i++)
   {
      ChangeSubSensor(i);
      badpixels[i].SetUpBadPixels(badfile,*revbandmap);
      Logger::Log("Number of bad pixels decoded from file: "+ToString(badpixels[i].NumBadPixels()));
   }

   return true;
}

//-------------------------------------------------------------------------
// Function to check if a bit has been set in the mask, if not (and it should be) it sets it
//-------------------------------------------------------------------------
void Calibration::AssignMaskValue(const unsigned int ele,const Specim::MaskType type)
{
   //Check the mask has been initialised
   if(data->Mask()==NULL)
   {
      throw "Error assigning mask value prior to mask being initialised.";
   }

   //Simple check to see if the mask has any value yet
   if(data->Mask()[ele]==0)
   {
      data->Mask()[ele]+=type;
      return;
   }
   //the mask has a value - check if the bit we want to set is already set
   if(data->Mask()[ele]&type)
   {
      //This bit has already been set - do nothing
      return;
   }
   else
   {
      //Set the bit
      data->Mask()[ele]+=type;   
      return;
   }
}

//-------------------------------------------------------------------------
// Function to average the FODIS region (over the samples) for this line
// Get an average value per band of the data
//-------------------------------------------------------------------------
bool Calibration::AverageFodis()
{
   //Only run this for Eagle sensors
   if(sensor->fodis==NULL)
   {
      Logger::Warning("No FODIS defined in sensor object - skipping in future.");
      return false;
   }
   //Check that the fodis array exists
   if(data->Fodis()==NULL)
   {
      throw "Error averaging the FODIS - fodis array has not been initialised.";
   }

   //Check that the sensor fodis values are useable
   if((sensor->fodis->LowerFodis()==sensor->fodis->UpperFodis())||(sensor->fodis->LowerFodis()>sensor->fodis->UpperFodis()))
   {
      Logger::Warning("Attempt to average FODIS for sensor without FODIS or with incorrect FODIS values in raw hdr file.");
      return false;      
   }

   double sum=0;
   unsigned int numbertoaverageover=0;
   for(unsigned int band=0;band<sensor->NumBands();band++)   
   {
      numbertoaverageover=0;
      //for each pixel of the fodis region for this line
      for(unsigned int p=sensor->fodis->LowerFodis();p<sensor->fodis->UpperFodis();p++)
      {
         if(data->Image()[band*sensor->NumSamples() + p] != 0)
         {
            //Sum up the fodis pixels per band per line
            sum=sum + data->Image()[band*sensor->NumSamples() + p]; 
            //Keep track of number of non-zero values  -only average over these
            numbertoaverageover++;
         }
      }

      //Average this value - note that FODIS array is the same size as
      //other arrays (nbands*nsamples) even though we're only interested
      //in the 1 sample per scan line.
      if(numbertoaverageover!=0)
      {
         sum=(sum / numbertoaverageover);
      }
      else
      {
         Logger::Warning("There are no valid FODIS pixels to average - all have value 0 in band:"+ToString(band));
         sum=0;
      }

      if(sum<sensor->CalibratedMax())
         data->Fodis()[band*sensor->NumSamples()]=(sum);
      else
         data->Fodis()[band*sensor->NumSamples()]=sensor->CalibratedMax();

      //reset sum for next band iteration
      sum=0;  
   } 
   return true;
}

//-------------------------------------------------------------------------
// Function to check how many frames have occured between start and end frame
// Note that start and end can be anything, not necessarily first and last frames
//-------------------------------------------------------------------------
int Calibration::CheckFrameCounter(unsigned int start,unsigned int end)
{
   //If they're the same then there are no frames between them
   if(start==end)
      return 0;

   //Start frame must occur before end frame
   if(start > end)
      throw "Error in checkframecounter - start frame should be less than end frame.";

   //Get frame counter for frame start
   double startcounter=sensor->bin->ReadCell(0,start,0);

   //Get frame counter for frame end
   double endcounter=sensor->bin->ReadCell(0,end,0);

   return static_cast<int>(endcounter-startcounter);
}

//-------------------------------------------------------------------------
// Function to clear the data arrays that are updated per line of the
// calibration procedure
//-------------------------------------------------------------------------
void Calibration::ClearPerlineData()
{
   if(data->Image()!=NULL)
   {
      data->ClearArray(data->Image());
   }
   if(data->Fodis()!=NULL)
   {
      data->ClearArray(data->Fodis());
   }
   if(data->Mask()!=NULL)
   {
      data->ClearArray(data->Mask());
   }
   if(data->BadPixMethod()!=NULL)
   {
      data->ClearArray(data->BadPixMethod());
   }
}

//-------------------------------------------------------------------------
// Function to read in the file that contains bad pixels to mask that
// have been detected at QC time and assign to the sensor
//-------------------------------------------------------------------------
void Calibration::ReadQCFailureFile(std::string qcfailurefile)
{
   std::ifstream filein;
   size_t index;
   std::string bstr;
   std::string sstr;
   unsigned int b=0,s=0;

   filein.open(qcfailurefile.c_str());
   if(!filein.is_open())
   {
      //Error opening file
      throw "An error occured whilst opening the qc failure pixel file - are you sure it exists?: " + qcfailurefile;
   }
   else
   {
      std::string tempstr="";
      //read a line from the file
      std::getline(filein,tempstr);
      while(filein.good())
      {
         //Trim whitespace to remove any errors due to this
         tempstr=TrimWhitespace(tempstr);
         //check that are 2 items on the line (band,sample)
         if(GetNumberOfItemsFromString(tempstr," ")!=2)
         {
            throw "An error occured whilst reading the qc failure pixel file - format of file should be ASCII: space separated band sample per line, I got: " + tempstr;           
         }
         //Get the band and sample and do lots of tests
         bstr=GetItemFromString(tempstr,0);
         //Test if safe to convert to unsigned integer
         index=bstr.find_first_not_of("0123456789");
         if(index == std::string::npos)
            b=StringToUINT(bstr);
         else
            throw "An error occured whilst reading the qc failure pixel file - non integer exists in band number: " + bstr;

         //Check number is less than total number of bands
         if(b >= sensor->NumBands())
         {
            throw "An error occured whilst reading the qc failure pixel file. Given band is greater than number of bands in file (0 indexed): " + ToString(b);
         }
         sstr=GetItemFromString(tempstr,1);
         //Test if safe to convert to unsigned integer
         index=sstr.find_first_not_of("0123456789");
         if(index == std::string::npos)
            s=StringToUINT(sstr);
         else
            throw "An error occured whilst reading the qc failure pixel file - non integer exists in sample number: " + sstr;

         if(s >= sensor->NumSamples())
         {
            throw "An error occured whilst reading the qc failure pixel file. Given sample is greater than number of samples in file (0 indexed): " + ToString(s);
         }
         //add the band sample pair to the vector in the sensor object
         sensor->qcfailures.push_back(Pair(StringToUINT(GetItemFromString(tempstr,0)),StringToUINT(GetItemFromString(tempstr,1))));

         //read a line from the file
         std::getline(filein,tempstr);
      }
      filein.close();
      filein.clear();
   } 
   //Log the information to the output
   Logger::Log("Will apply QC Failure flags to the mask for the following band, sample pairs: ");
   for(std::vector<Pair>::iterator it=sensor->qcfailures.begin();it!=sensor->qcfailures.end();it++)
   {
      Logger::Log(" "+ToString((*it).band)+" "+ToString((*it).sample));
   }
}


//-------------------------------------------------------------------------
// BadPixels constructor- SetUpBadPixels needs to be run independantly 
//-------------------------------------------------------------------------
BadPixels::BadPixels()
{
   badpixels=NULL;
   badpixelmethod=NULL;
   bpmethod_descriptor=NULL;
   bandnotinuse=9999;
}

//-------------------------------------------------------------------------
// BadPixels destructor
//-------------------------------------------------------------------------
BadPixels::~BadPixels()
{
   if(badpixels!=NULL)
      delete[] badpixels;
   if(badpixelmethod!=NULL)
      delete[] badpixelmethod;
   if(bpmethod_descriptor!=NULL)
      delete[] bpmethod_descriptor;
}

//-------------------------------------------------------------------------
// Set up the bad pixel object from a file and a band mapping map
//-------------------------------------------------------------------------
void BadPixels::SetUpBadPixels(std::string badpixelfilename,std::map<int,int> &revbandmap)
{
   std::ifstream badin;
   badin.open(badpixelfilename.c_str());
   if(!badin.is_open())
   {
      //Error opening file
      throw "An error occured whilst opening the bad pixel calibration file - are you sure it exists?" + badpixelfilename;
   }
   else
   {
      std::string tempstr="";
      while(badin.good())
      {
         //store the bad pixels in a string stream object
         std::getline(badin,tempstr);
         badpixelstream<<tempstr<<std::endl; 
      }
      badin.close();
      badin.clear();

      //Now call a function to decode the file - i.e. interpret the data in the stringstream
      DecodeBadPixel(revbandmap);
   }  
}


//-------------------------------------------------------------------------
// Function to decode the bad pixel stringstream into an unsigned int array
// requires the reverse band mapping from the calibration to identify 
// non-used bands (i.e. map raw bands to calibration bands).
// This version is for Specim Calibrated Bad Pixel files
//-------------------------------------------------------------------------
void BadPixels::DecodeSpecimBadPixels(std::map<int,int> revbandmap)
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
void BadPixels::DecodeARSFBadPixels(std::map<int,int> revbandmap)
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

   unsigned int method_counter=0;
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
   num_method_descriptors=method_counter;
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
   unsigned int meth=0;
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

   if(!badpixelstream.good() && (nbadpixels!=0))
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
         case 'F':
            if((badpixelmethod[id]&F) == 0)
               badpixelmethod[id]+=F;
            break;
         default:
            throw std::string("Unrecognised bad pixel detection method in bad pixel file. Expected one of A,B,C,D,E,F but got: ")+method;
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
void BadPixels::DecodeBadPixel(std::map<int,int> revbandmap)
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
