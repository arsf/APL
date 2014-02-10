//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#include "mainworker.h"

//-------------------------------------------------------------------------
// Constructor for mainworker class
//-------------------------------------------------------------------------
MainWorker::MainWorker(std::string rawfile,std::string outfile,std::string cl,bool force)
{
   cal=NULL;
   bwmask=NULL;
   bwmaskmethod=NULL;
   bwfodis=NULL;
   bwimage=NULL;
   eagle=NULL;
   hawk=NULL;
   sensor=NULL;

   sensor=new Specim(rawfile,force);
   numbercalibratedframes=sensor->br->GetNumImageFrames() + sensor->br->GetTotalMissingFrames();
   outputfileprefix=outfile;
   lowersample=sensor->LowerScanlineLimit();
   uppersample=sensor->UpperScanlineLimit();
   startline=0;
   nummissingscanspriortostartline=0;
   endline=sensor->br->GetNumImageFrames();
   commandlinetext=cl;

   if(CheckSensorID(EAGLE,this->sensor->SensorID()))
   {
      Logger::Log("Eagle sensor detected - created eagle object.");
      delete sensor;
      sensor=NULL;
      eagle=new Eagle(rawfile);
      sensor=eagle;
   }
   else if(CheckSensorID(HAWK,this->sensor->SensorID()))
   {
      Logger::Log("Hawk sensor detected - created hawk object.");
      delete sensor;
      sensor=NULL;
      hawk=new Hawk(rawfile);
      sensor=hawk;
   }
   else
      throw "Unrecognised sensor type in raw file, got ID: "+ToString(sensor->SensorID());

   //Initialise task map to be false for everything
   tasks[remove_dark_frames]=false;
   tasks[apply_gains]=false;
   tasks[calibrate_fodis]=false;
   tasks[insert_missing_scans]=false;
   tasks[smear_correct]=false;
   tasks[output_mask]=false;
   tasks[flip_bands]=false;
   tasks[flip_samples]=false;
   tasks[apply_qcfailures]=false;

}

//-------------------------------------------------------------------------
// Constructor to try and force using a particular given sensor type
//-------------------------------------------------------------------------
MainWorker::MainWorker(std::string rawfile,std::string outfile,char csensor,std::string cl,bool force)
{
   cal=NULL;
   bwmask=NULL;
   bwmaskmethod=NULL;
   bwfodis=NULL;
   bwimage=NULL;
   eagle=NULL;
   hawk=NULL;
   sensor=NULL;

   if(csensor=='e')
   {
      eagle=new Eagle(rawfile,force);
      sensor=eagle;
   }
   else if(csensor=='h')
   {
      hawk=new Hawk(rawfile,force);
      sensor=hawk;
   }
   else
      throw "Unrecognised sensor type in raw file, got sensor char: "+ToString(csensor);

   numbercalibratedframes=sensor->br->GetNumImageFrames() + sensor->br->GetTotalMissingFrames();
   outputfileprefix=outfile;
   lowersample=sensor->LowerScanlineLimit();
   uppersample=sensor->UpperScanlineLimit();
   startline=0;
   nummissingscanspriortostartline=0;
   endline=sensor->br->GetNumImageFrames();
   commandlinetext=cl;

   //Initialise task map to be false for everything
   tasks[remove_dark_frames]=false;
   tasks[apply_gains]=false;
   tasks[calibrate_fodis]=false;
   tasks[insert_missing_scans]=false;
   tasks[smear_correct]=false;
   tasks[output_mask]=false;
   tasks[output_mask_method]=false;
   tasks[flip_bands]=false;
   tasks[flip_samples]=false;
   tasks[apply_qcfailures]=false;
}

//-------------------------------------------------------------------------
// Destructor for mainworker
//-------------------------------------------------------------------------
MainWorker::~MainWorker()
{
   if(eagle!=NULL)
      delete eagle;
   if(hawk!=NULL)
      delete hawk;   

   if(bwmask!=NULL)
      bwmask->Close();
   if(bwmaskmethod!=NULL)
      bwmaskmethod->Close();
   if(bwfodis!=NULL)
      bwfodis->Close();
   if(bwimage!=NULL)
      bwimage->Close();

   if(cal!=NULL)
      delete cal;
}

//-------------------------------------------------------------------------
// Set up the calibration object with an optional calibration file
//-------------------------------------------------------------------------
void MainWorker::InitialiseCalibration(std::string calfile,std::string externaldarkframes,std::string qcfailurefile)
{
   //Set up the calibration object
   cal=new Calibration(sensor,calfile);
   //Test the calibration file is suitable
   cal->TestCalfile();
   //calculate the dark frame averages, if requested, from the optional
   //externaldarkframes filename (uses the raw file if no external file given)
   if(tasks[remove_dark_frames])
      cal->InitialiseDarkFrames(externaldarkframes);

   //Now set up the fodis if required
   if(tasks[calibrate_fodis])
      cal->InitialiseFodis();

   if((tasks[calibrate_fodis]==false)&&(CheckSensorID(EAGLE,this->sensor->SensorID())))
   {
      //We shall leave the fodis region on the image - no way to get this
      //info externally so have to hard code it here
      lowersample=0;
   }

   //If hawk - read in the bad pixel file
   if(CheckSensorID(HAWK,this->sensor->SensorID()))
   {
      cal->ReadBadPixelFile();
      //If the bad pixel file is a specim one then do not output the method file
      if(dynamic_cast<Hawk*>(sensor)->arsfbadpixelfiletype == false)
      {
         tasks[output_mask_method]=false;
      }
      else
      {
         cal->InitialiseBadPixMethod();
      }
      //Force the FODIS and Smear Correction to be false here
      tasks[calibrate_fodis]=false;
      tasks[smear_correct]=false;
   }

   //If a qcfailure file of bad pixels has been given, read them in and set up sensor
   if(tasks[apply_qcfailures])
   {
      cal->ReadQCFailureFile(qcfailurefile);
   }

}

//-------------------------------------------------------------------------
// Return the number of calibrated image lines the final image will have
//-------------------------------------------------------------------------
unsigned int MainWorker::GetNumCalibratedImageLines()
{
   //Check if full line is being processed or a user defined chunk
   if((startline==0)&&(endline==sensor->br->GetNumImageFrames()))
   {
      //Full line is being processed - check if missing scans are to be inserted
      if(tasks[insert_missing_scans])
         return sensor->br->GetNumImageFrames()+sensor->br->GetTotalMissingFrames();
      else
         return sensor->br->GetNumImageFrames();
   }
   else
   {
      //User defined chunk of a line being processed
      //check if missing scans are to be inserted
      if(tasks[insert_missing_scans])
         return (endline-startline)+sensor->br->GetMissingFramesBetweenLimits(startline,endline);
      else
         //return sensor->br->GetNumImageFrames();      
         return (endline-startline);
   }
}

//-------------------------------------------------------------------------
// Return the number of calibrated image samples the final image will have
//-------------------------------------------------------------------------
unsigned int MainWorker::GetNumCalibratedImageSamples()
{
   //if(tasks[define_sample_bounds]) //user defined sample limits
      return uppersample-lowersample+1;
   //else
   //   return sensor->UpperScanlineLimit() - sensor->LowerScanlineLimit() +1;
}

//-------------------------------------------------------------------------
// Apply all the calibration routines that have been requested to the
// specific line of data passed to the function
//-------------------------------------------------------------------------
void MainWorker::DoCalibrationForLine(unsigned int line)
{
   //Read in the line of data from the raw file
   cal->ReadLineOfRaw(line);

   //Flag the pixels for under/over flows etc
   //Note that average dark frames have to be initialised before flagging pixels
   cal->FlagPixels();

   //If required then remove the average dark frame values from the image data
   if(tasks[remove_dark_frames])
      cal->RemoveDarkFrames();

   //If required then smear correct the image data - returns false if not
   //eagle data and updates the tasks so that it is not run next time
   if(tasks[smear_correct])
      tasks[smear_correct]=cal->SmearCorrect();

   //If required apply the gain values to the image data
   if(tasks[apply_gains])
      cal->ApplyGains();

   //If required average the fodis data up - returns false if not
   //eagle data and updates the tasks so that it is not run next time
   if(tasks[calibrate_fodis])
      tasks[calibrate_fodis]=cal->AverageFodis();

   //If required flip the data spectraly from red to blue
   if(tasks[flip_bands])
      const_cast<Data*>(cal->pData())->TransformArrays(sensor->NumBands(),sensor->NumSamples(),Data::BAND);

   //If required flip the data spatially from left to right
   if(tasks[flip_samples])
      const_cast<Data*>(cal->pData())->TransformArrays(sensor->NumBands(),sensor->NumSamples(),Data::SAMPLE);

   //Write out the required datasets and clear the data in memory
   WriteOutData(Normal);
   cal->ClearPerlineData();
}

//-------------------------------------------------------------------------
// Function to set up the bil writers for writing out the image, fodis, mask
//-------------------------------------------------------------------------
void MainWorker::InitialiseWriters()
{
   std::string strOutputFilename="";
   strOutputFilename=outputfileprefix;

   if(bwimage!=NULL)
   {
      //Already been initialised
      return;
   }

   //----------------------------------------------------------------------
   //Create a bil writer for the main image data
   //----------------------------------------------------------------------
   Logger::Log("Will write calibrated image data to: "+strOutputFilename);
   bwimage=new BILWriter(strOutputFilename,bwimage->uint16,GetNumCalibratedImageLines(),GetNumCalibratedImageSamples(),sensor->NumBands(),'w');
   //Copy over some information from the raw hdr file to the calibrated hdr file
   TransferHeaderInfo(bwimage);
   //Add x and y start values to the level1 bil file - both zero based
   bwimage->AddToHdr("x start = "+ToString(lowersample));
   bwimage->AddToHdr("y start = "+ToString(startline));
   bwimage->AddToHdr("dropped scans before y start = "+ToString(nummissingscanspriortostartline));
   //Add command line
   bwimage->AddToHdr(";The command line used to process the data: "+commandlinetext);   
   //Add raw file name
   bwimage->AddToHdr(";Raw data file: "+sensor->RawFilename());   
   //Also add the name of the calibration file that is to be used
   bwimage->AddToHdr(";The data has been calibrated using the file: "+cal->CalibrationFile());   
   //Add the wavelength units - this will need to be hard coded and assumed to be nm
   bwimage->AddToHdr("Wavelength units = nm");   
   //Add the calibrated data units 
   bwimage->AddToHdr("Radiance data units = "+sensor->CalibratedUnits());

   //----------------------------------------------------------------------
   //Create a bil writer for the fodis data if it has been requested
   //----------------------------------------------------------------------
   if(tasks[calibrate_fodis])
   {
      strOutputFilename=outputfileprefix+"_FODIS.bil";
      Logger::Log("Will write calibrated FODIS data to: "+strOutputFilename);
      bwfodis=new BILWriter(strOutputFilename,bwfodis->uint16,GetNumCalibratedImageLines(),1,sensor->NumBands(),'w');
      //Add some hdr information about what the data is
      bwfodis->AddToHdr(";File containing averaged per-scan radiometrically calibrated data from the fibre optic downwelling irradiance sensor.");
      //Units the data is in
      bwfodis->AddToHdr("Units = "+(dynamic_cast<Eagle*>(sensor))->FodisUnits());     
   }

   //----------------------------------------------------------------------
   //Create a bil writer for the mask data if it has been requested
   //----------------------------------------------------------------------
   if(tasks[output_mask])
   {
      strOutputFilename=outputfileprefix+"_mask.bil";
      Logger::Log("Will write calibrated image mask data to: "+strOutputFilename);
      bwmask=new BILWriter(strOutputFilename,bwmask->uchar8,GetNumCalibratedImageLines(),GetNumCalibratedImageSamples(),sensor->NumBands(),'w');

      //Add some information to the header file
      //Add x and y start values to the level1 bil file - both zero based
      bwmask->AddToHdr("x start = "+ToString(lowersample));
      bwmask->AddToHdr("y start = "+ToString(startline));
      bwmask->AddToHdr("dropped scans before y start = "+ToString(nummissingscanspriortostartline));
      std::string waves=sensor->br->FromHeader("Wavelength");      
      if(tasks[flip_bands])
      {
         waves=ReverseWavelengthOrder(waves);
      }
      //Tidy up wavelength array to remove {, and ,} and replace with { and }
      waves=sensor->br->TidyForHeader(waves);
      bwmask->AddToHdr("Wavelength = "+waves);
      bwmask->AddToHdr("Wavelength units = nm");
      bwmask->AddToHdr(";Mask file for "+outputfileprefix);
      bwmask->AddToHdr(";Values of: \n; "+
                      ToString(Specim::Good)+" = Good data.\n; "+
                      ToString(Specim::UnderFlow)+" = Underflows.\n; "+
                      ToString(Specim::OverFlow)+" = Overflows.\n; "+
                      ToString(Specim::Badpixel)+" = Hawk CCD bad pixels.\n; "+
                      ToString(Specim::SmearAffected)+" = Pixel affected by uncorrected smear.\n; "+
                      ToString(Specim::DroppedScan)+" = Dropped scans.\n; "+
                      ToString(Specim::CorruptData)+" = Corrupt raw data.\n;"+
                      ToString(Specim::QCFailure)+" = Quality control failure.\n");
   }
   //----------------------------------------------------------------------
   //Create a bil writer for the bad pixel method data if it has been requested
   //----------------------------------------------------------------------
   if(tasks[output_mask_method])
   {
      strOutputFilename=outputfileprefix+"_mask-badpixelmethod.bil";
      Logger::Log("Will write bad pixel method data to: "+strOutputFilename);
      hawk=dynamic_cast<Hawk*>(sensor);
      bwmaskmethod=new BILWriter(strOutputFilename,bwmaskmethod->uchar8,GetNumCalibratedImageLines(),GetNumCalibratedImageSamples(),sensor->NumBands(),'w');

      //Add some information to the header file
      //Add x and y start values to the level1 bil file - both zero based
      bwmaskmethod->AddToHdr("x start = "+ToString(lowersample));
      bwmaskmethod->AddToHdr("y start = "+ToString(startline));
      bwmaskmethod->AddToHdr("dropped scans before y start = "+ToString(nummissingscanspriortostartline));
      std::string waves=sensor->br->FromHeader("Wavelength");      
      if(tasks[flip_bands])
      {
         waves=ReverseWavelengthOrder(waves);
      }
      //Tidy up wavelength array to remove {, and ,} and replace with { and }
      waves=sensor->br->TidyForHeader(waves);
      bwmaskmethod->AddToHdr("Wavelength = "+waves);
      bwmaskmethod->AddToHdr("Wavelength units = nm");
      bwmaskmethod->AddToHdr(";Bad CCD pixel detection methods. Values of: \n; "+
                      ToString(Hawk::None)+" = Not flagged as a bad ccd pixel.\n; "+
                      ToString(Hawk::A)+" = "+hawk->MethodDescriptor()[0]+"\n; "+
                      ToString(Hawk::B)+" = "+hawk->MethodDescriptor()[1]+"\n; "+
                      ToString(Hawk::C)+" = "+hawk->MethodDescriptor()[2]+"\n; "+
                      ToString(Hawk::D)+" = "+hawk->MethodDescriptor()[3]+"\n; "+
                      ToString(Hawk::E)+" = "+hawk->MethodDescriptor()[4]+"\n");
   }
}

//-------------------------------------------------------------------------
// Function to write out data using the given file prefix
//-------------------------------------------------------------------------
void MainWorker::WriteOutData(OutputDataFlag flag)
{
   //This only needs to be called once - probably better to move it somewhere else
   InitialiseWriters();

   //For each initialised element of the data object - write it out
   if(cal->pData()->Image()!=NULL)
   {
      if(flag==Normal)
         bwimage->WriteDataToLineSection(cal->pData()->Image(),sensor->NumSamples(),lowersample,uppersample);
      else if((flag==MissingScan)||(flag==CorruptData))
         bwimage->WriteLineWithValue(0);
      else
         throw "Unrecognised flag in WriteOutData - writing Image().";

   }
   if((cal->pData()->Mask()!=NULL)&&(tasks[output_mask]))
   {
      if(flag==Normal)
         bwmask->WriteDataToLineSection(cal->pData()->Mask(),sensor->NumSamples(),lowersample,uppersample);
      else if(flag==MissingScan)
         bwmask->WriteLineWithValue(sensor->DroppedScan);
      else if(flag==CorruptData)
         bwmask->WriteLineWithValue(sensor->CorruptData);
      else
         throw "Unrecognised flag in WriteOutData - writing Mask().";
   }
   if((cal->pData()->BadPixMethod()!=NULL)&&(tasks[output_mask_method]))
   {
      if(flag==Normal)
         bwmaskmethod->WriteDataToLineSection(cal->pData()->BadPixMethod(),sensor->NumSamples(),lowersample,uppersample);
      else if((flag==MissingScan)||(flag==CorruptData))
         bwmaskmethod->WriteLineWithValue(0);
      else
         throw "Unrecognised flag in WriteOutData - writing BadPixelMethod().";
   }
   if((cal->pData()->Fodis()!=NULL)&&(tasks[calibrate_fodis]))
   {
      if(flag==Normal)
         bwfodis->WriteDataToLineSection(cal->pData()->Fodis(),sensor->NumSamples(),0,0);
      else if((flag==MissingScan)||(flag==CorruptData))
         bwfodis->WriteLineWithValue(0);
      else
         throw "Unrecognised flag in WriteOutData - writing Fodis().";
   }
}

//-------------------------------------------------------------------------
// Function to reverse the order of the wavelengths in the bil hdr file
// Use this if you have used the flip_bands task
//-------------------------------------------------------------------------
std::string MainWorker::ReverseWavelengthOrder(std::string wavelengths)
{
   std::string revwavelengths;//String to store the reversed wavelengths in
   std::string tmpstr;//temporary buffer string
   size_t index1=0,index2=0;

   index1=wavelengths.find(';',index2);
   index2=wavelengths.find(';',index1+1);

   //wavelengths string is delimitted by ;
   while((index1!=std::string::npos)&&(index2!=std::string::npos))
   {
      tmpstr=wavelengths.substr(index1,index2-index1);
      revwavelengths=revwavelengths.insert(0,tmpstr);
      index1=wavelengths.find(';',index2);
      index2=wavelengths.find(';',index1+1);
   }

   //Append the {} onto the string
   revwavelengths=revwavelengths.insert(0,"{");
   revwavelengths=revwavelengths.append(";}");
   //Also need to add a comma to first and remove a comma from last wavelength
   index1=revwavelengths.find(';',0);
   index2=revwavelengths.find(';',index1+1);
   revwavelengths=revwavelengths.insert(index2,",");
   index1=revwavelengths.find_last_of(',');
   revwavelengths.erase(index1,1);

   return revwavelengths; //return the reversed wavelengths
}

//-------------------------------------------------------------------------
// Function to transfer the important parts of the raw hdr file to the 
// calibrated hdr file.
//-------------------------------------------------------------------------
void MainWorker::TransferHeaderInfo(BILWriter* const bw)
{
   std::string strtransfer("");

   strtransfer=sensor->br->FromHeader("sensor type"); //get the sensor type
   if(!strtransfer.empty()) //if the string is not empty
   {
      strtransfer.insert(0,"sensor type = ");
      //Need some way of removing the ; from the string and replacing with \n
      strtransfer=sensor->br->TidyForHeader(strtransfer);
      bw->AddToHdr(strtransfer); //add the string to the header
   }

   strtransfer=sensor->br->FromHeader("acquisition date"); //get the acquisition date
   if(!strtransfer.empty()) //if the string is not empty
   {
      strtransfer.insert(0,"acquisition date = ");
      //Need some way of removing the ; from the string and replacing with \n
      strtransfer=sensor->br->TidyForHeader(strtransfer);
      bw->AddToHdr(strtransfer); //add the string to the header
   }

   strtransfer=sensor->br->FromHeader("fps"); //get the frames per second type
   if(!strtransfer.empty()) //if the string is not empty
   {
      strtransfer.insert(0,"fps = ");
      //Need some way of removing the ; from the string and replacing with \n
      strtransfer=sensor->br->TidyForHeader(strtransfer);
      bw->AddToHdr(strtransfer); //add the string to the header
   }

   strtransfer=sensor->br->FromHeader("binning"); //get the spatial/spectral binning
   if(!strtransfer.empty()) //if the string is not empty
   {
      strtransfer.insert(0,"binning = ");
      //Need some way of removing the ; from the string and replacing with \n
      strtransfer=sensor->br->TidyForHeader(strtransfer);
      bw->AddToHdr(strtransfer); //add the string to the header
   }

   strtransfer=sensor->br->FromHeader("tint"); //get the intergration time
   if(!strtransfer.empty()) //if the string is not empty
   {
      strtransfer.insert(0,"tint = ");
      //Need some way of removing the ; from the string and replacing with \n
      strtransfer=sensor->br->TidyForHeader(strtransfer);
      bw->AddToHdr(strtransfer); //add the string to the header
   }

   strtransfer=sensor->br->FromHeader("sensorid"); //get the sensor id
   if(!strtransfer.empty()) //if the string is not empty
   {
      strtransfer.insert(0,"sensorid = ");
      //Need some way of removing the ; from the string and replacing with \n
      strtransfer=sensor->br->TidyForHeader(strtransfer);
      bw->AddToHdr(strtransfer); //add the string to the header
   }

   strtransfer=sensor->br->FromHeader("GPS Start Time"); //get the GPS start time
   if(!strtransfer.empty()) //if the string is not empty
   {
      strtransfer.insert(0,"GPS Start Time = ");
      //Need some way of removing the ; from the string and replacing with \n
      strtransfer=sensor->br->TidyForHeader(strtransfer);
      bw->AddToHdr(strtransfer); //add the string to the header
   }

   strtransfer=sensor->br->FromHeader("GPS Stop Time"); //get the GPS stop time
   if(!strtransfer.empty()) //if the string is not empty
   {
      strtransfer.insert(0,"GPS Stop Time = ");
      //Need some way of removing the ; from the string and replacing with \n
      strtransfer=sensor->br->TidyForHeader(strtransfer);
      bw->AddToHdr(strtransfer); //add the string to the header
   }

   strtransfer=sensor->br->FromHeader("NavSync Timing"); //get the nav sync timing
   if(!strtransfer.empty()) //if the string is not empty
   {
      strtransfer.insert(0,"NavSync Timing = ");
      //Need some way of removing the ; from the string and replacing with \n
      strtransfer=sensor->br->TidyForHeader(strtransfer);
      bw->AddToHdr(strtransfer); //add the string to the header
   }

   strtransfer=sensor->br->FromHeader("Wavelength"); //get the wavelength array
   if(!strtransfer.empty()) //if the string is not empty
   {
      if(tasks[flip_bands])
      {
         strtransfer=ReverseWavelengthOrder(strtransfer);
      }
      strtransfer.insert(0,"Wavelength = ");
      //Need some way of removing the ; from the string and replacing with \n
      strtransfer=sensor->br->TidyForHeader(strtransfer);
      bw->AddToHdr(strtransfer); //add the string to the header
   }

   strtransfer=sensor->br->FromHeader("fwhm"); //get the fwhm array
   if(!strtransfer.empty()) //if the string is not empty
   {
      if(tasks[flip_bands])
      {
         strtransfer=ReverseWavelengthOrder(strtransfer);
      }
      strtransfer.insert(0,"fwhm = ");
      //Need some way of removing the ; from the string and replacing with \n
      strtransfer=sensor->br->TidyForHeader(strtransfer);
      bw->AddToHdr(strtransfer); //add the string to the header
   }   
}

//-------------------------------------------------------------------------
// Function to return the tasks as a formatted string
//-------------------------------------------------------------------------
std::string MainWorker::TasksAsString()
{
   std::stringstream formattedtasks;

   formattedtasks<<" Remove dark frames: "<<tasks[remove_dark_frames]<<std::endl;
   formattedtasks<<" Calibrate FODIS (if exists): "<<tasks[calibrate_fodis]<<std::endl;
   formattedtasks<<" Insert missing scans as line of 0's: "<<tasks[insert_missing_scans]<<std::endl;
   formattedtasks<<" Apply radiometric gains: "<<tasks[apply_gains]<<std::endl;
   formattedtasks<<" Smear correct the data (if Eagle): "<<tasks[smear_correct]<<std::endl;
   formattedtasks<<" Output the mask file: "<<tasks[output_mask]<<std::endl;
   formattedtasks<<" Output the mask method file: "<<tasks[output_mask_method]<<std::endl;
   formattedtasks<<" Flip the raw data spectrally (red to blue): "<<tasks[flip_bands]<<std::endl;
   formattedtasks<<" Flip the raw data spatially (left to right): "<<tasks[flip_samples]<<std::endl;
   formattedtasks<<" Apply QC failure bad pixels: "<<tasks[apply_qcfailures]<<std::endl;

   return formattedtasks.str();
}

