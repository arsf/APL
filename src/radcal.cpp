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
#include <cerrno>
#include <cstdlib>
#include <iostream>

#include "os_dependant.h"
#include "logger.h"
#include "commandline.h"
#include "calibration.h"
#include "mainworker.h"
#include "specimsensors.h"


//-------------------------------------------------------------------------
// Software description
//-------------------------------------------------------------------------
const std::string DESCRIPTION="Radiometric Calibration Software";

//-------------------------------------------------------------------------
//Number of options that can be on command line
//-------------------------------------------------------------------------
const int number_of_possible_options = 22;

//-------------------------------------------------------------------------
//Option names that can be on command line
//-------------------------------------------------------------------------
const std::string availableopts[number_of_possible_options]={
"-input",
"-calfile",
"-output",
"-darkfile",
"-sensor",
"-lines",
"-FLIPSAMPLES",
"-FLIPBANDS",
"-NOFLIP",
"-NOFODIS",
"-NOMASK",
"-NOMISSSCAN",
"-NODARK",
"-NORAD",
"-NOSMEAR",
"-avdark",
"-gains",
"-corruptscans",
"-qcfailures",
"-darkforce",
"-help"
};  

//-------------------------------------------------------------------------
//Descriptions of option that can be on command line
//-------------------------------------------------------------------------
const std::string optsdescription[number_of_possible_options]={
"Raw Eagle/Hawk filename",
"Calibration filename (excluding .cal extension)",
"Name of output radiometrically calibrated BIL file.",
"Name of file containing dark frames. Default is to use dark frames from the input raw image.",
"Sensor ID 'e' for Eagle, 'h' for Hawk. Default is to auto detect from file.",
"Define a section of the image to process using a start and end scan line number. Default is for full image.",
"Flips the data spatially (left to right). This is default for Hawk.",
"Flips the data spectrally (red to blue). This is default for Eagle.",
"Do no flipping of data at all - i.e. Eagle will not flip bands, Hawk will not flip samples.",
"Do not output the FODIS data.",
"Do not output a mask file.",
"Do not insert missing frames into the output data.",
"Do not subtract the dark frames from the raw data.",
"Do not convert the DN to radiance.",
"Do not apply the smear correction to Eagle data (no effect for Hawk data)",
"Outputs the average dark values to a BIL file of given name",
"Outputs the binned calibration gain values to a BIL file of given name.",
"A space separated list of scan lines to ignore and mark as corrupt.",
"Give a filename containing a list of space separated band sample pairs (one pair per line) of pixels to mask as QCFailure. NOTE these pixels are in the raw data geometry with band/sample starting from 0.",
"Force the use of the autodarkstartline number given in the hdr file. If using this please check first that it is correct - this should be used as a last resort.",
"Show this help text."
}; 



//-------------------------------------------------------------------------
// Main program function for the radiometric calibration
//-------------------------------------------------------------------------
int main(int argc,char* argv[])
{
   //----------------------------------------------------------------------
   //Set up variables and initialise default values
   //----------------------------------------------------------------------

   //Set the precision of cout streaming
   std::cout.precision(10);

   unsigned int startline=0,endline=0; //these are used to define the lines of data to be processed
   std::string strRawFileName; //raw file filename and path to it eg /users/rsg/arsf/eagle/VNIRse1.raw
   std::string strCalibFileName; //calibration file prefix and path (eg .cal file but without the .cal extension)
   std::string strDarkFileName; //dark file to use if dont want to use raw image dark lines  
   std::string strMaskFileName; //mask file name to use to output the mask to
   std::string strOutputDataFileName; //output filename to use
   std::string strAverageDarkFramesFile; //filename to output averaged dark frames to
   std::string strBinnedGainsFile; //filename to output binned gains to
   std::string strQCFailureFileName; //filename that contains qcfailure bad pixels to be masked detected at QC time 
   char sensor='0';

   unsigned int number_of_corrupt_lines=0;
   unsigned int* corrupt_lines=NULL;

   bool OUTPUT_AVERAGED_DARKFRAMES=false;
   bool OUTPUT_BINNED_GAINS=false;

   //Do we want to flip the bands from red to blue
   bool DO_FLIP_BANDS=false;
   //Do we want to flip the samples from right to left
   bool DO_FLIP_SAMPLES=false;
   bool DO_INSERT_MISSINGSCANS=true;
   bool DO_REMOVE_DARK=true;
   bool DO_CALIBRATE_FODIS=true;
   bool DO_SMEAR_CORRECT=true;
   bool DO_APPLY_GAINS=true;
   bool DO_OUTPUT_MASK=true;
   bool DO_APPLY_QCFAILURES=false;

   //Do we want to force the use of autodarkstart frame number
   bool DARKFORCE=false;

   Logger log; //create logger to terminal only
   std::stringstream strout; //string to hold text messages in

   //Get exe name without the path
   std::string niceexename=std::string(argv[0]);
   niceexename=niceexename.substr(niceexename.find_last_of("/\\")+1);

   //Output a "nice" header
   Logger::FormattedInformation(niceexename,VERSION,DESCRIPTION);

   //----------------------------------------------------------------------
   // Enter the command line handling loop and sort out the user defined data
   //----------------------------------------------------------------------

   //Create a pointer to a command line object
   CommandLine* cl=NULL;
   try
   {
      Logger::Log("Reading command line and sorting options.");
      //Create the command line object 
      cl=new CommandLine(argv,argc);
      //Check if anything went wrong (most likely an exception will have been thrown but lets be safe)
      if(cl->IsGood()==false)
      {
         Logger::Log("An error has occurred with the command line"); //Log that an error has occurred
      }      

      //----------------------------------------------------------------------
      // Get some information about the system running the program and log it
      //----------------------------------------------------------------------
      ComputerInfo cinfo;
      Logger::Log(cinfo.GetOutput());

      //Get additional info and extra checks if in debug mode
      #ifdef DEBUG
         Logger::Log(cl->ReviewCL(availableopts,number_of_possible_options,true)); 
      #endif

      //check the command line options given against the list of available options in the software
      std::string badoptions;
      int retval=cl->CheckAvailableOptions(availableopts,number_of_possible_options,&badoptions); 
      if(retval<0)
      {
         //Options are on commnd line that are not available in this software
         strout<<"There are "<<(-retval)<<" unrecognised options on command line: "<<badoptions;
         throw CommandLine::CommandLineException(strout.str());
      }

      //get the arguments from the command line
      //----------------------------------------------------------------------
      //Check for help first since if there we dont do anything other than output help text
      //----------------------------------------------------------------------
      if(cl->OnCommandLine("-help"))
      {
         //CommandLineUsage();
         Logger::Log(cl->ProgramUsage(number_of_possible_options,availableopts,optsdescription));
         delete cl;
         log.Close();//close the logger
         exit(0);                  
      } 

      //----------------------------------------------------------------------
      //These 2 (input, output) must be on command line - check this
      //----------------------------------------------------------------------
      if((cl->OnCommandLine("-input"))&&(cl->OnCommandLine("-output")))
      {
         if(cl->GetArg("-input").compare(optiononly)!=0) //check that there is an argument following this
         {
            strRawFileName=cl->GetArg("-input");   //get the raw file filename (+path)        
            strRawFileName=GetExistingFilePath(strRawFileName,true);
         }
         else
         {
            //Throw an exception because there is no raw filename given
            throw CommandLine::CommandLineException("Argument -input must precede the sensor raw filename.");
         }

         if(cl->GetArg("-output").compare(optiononly)!=0)
            strOutputDataFileName=cl->GetArg("-output"); //get the output filename.
         else
         {
            //Throw an exception because there is no output filename given
            throw CommandLine::CommandLineException("Argument -output must precede the calibrated output BIL filename.");
         }
      }
      else
      {
         //Throw an exception because there is not both -input and -calfile given
         throw CommandLine::CommandLineException("Arguments -input and -output must be present on the command line.");         
      }


      //----------------------------------------------------------------------
      //The following are optional
      //----------------------------------------------------------------------

      if(cl->OnCommandLine("-calfile"))
      {
         if(cl->GetArg("-calfile").compare(optiononly)!=0)
         {
            strCalibFileName=cl->GetArg("-calfile"); //get the calibration file prefix (+path)
         }
         else
         {
            //Throw an exception because there is no calibration filename given
            throw CommandLine::CommandLineException("Argument -calfile must precede the sensor calibration prefix.");
         }
      }

      //-------------------------------------------------------------------
      // Specify (force) a particular sensor type
      //-------------------------------------------------------------------
      if(cl->OnCommandLine("-sensor"))
      {
         std::string strSensor=cl->GetArg("-sensor"); //get the sensor id ('e'agle or 'h'awk)
         if((strSensor.compare("h")!=0)&&(strSensor.compare("e")!=0))
         {
            //unrecognised sensor type
            throw CommandLine::CommandLineException("Sensor type "+strSensor+" is unrecognised.");         
         }
         else
            sensor=*strSensor.c_str(); //convert to char
      }

      //----------------------------------------------------------------------
      // Which lines to calibrate from the raw image, give it start and end lines
      //----------------------------------------------------------------------
      if(cl->OnCommandLine("-lines"))
      {
         if(cl->NumArgsOfOpt("-lines")!=2)
         {
            //lines option should have 2 arguments
            throw CommandLine::CommandLineException("Lines option should have exactly 2 parameters. Got: "+cl->GetArg("-lines"));      
         }
         //Else there were 2 parameters
         startline=StringToUINT(TrimWhitespace(cl->GetArg("-lines",0))); 
         endline=StringToUINT(TrimWhitespace(cl->GetArg("-lines",1)));

         if(startline>=endline)
         {
            //startline is greater than or same as endline. This cannot be
            throw CommandLine::CommandLineException("Lines start must be less than end. Got: "+ToString(startline)+", "+ToString(endline));         
         }
         Logger::Log("Using line limits of "+ToString(startline)+" "+ToString(endline));
      }
      else
      {
         //-lines option was not on command line
         startline=0;
         endline=0;
      }

      //----------------------------------------------------------------------
      // Get the dark frame file if on command line
      //----------------------------------------------------------------------
      if(cl->OnCommandLine("-darkfile"))
      {
         if(cl->NumArgsOfOpt("-darkfile")!=1)
         {
            throw CommandLine::CommandLineException("-darkfile should immediately preceed the filename of the dark file to use. Got: "+cl->GetArg("-darkfile"));      
         }
         strDarkFileName=cl->GetArg("-darkfile"); //get the dark filename that will be used to get the dark lines from instead of the raw image
         strDarkFileName=GetExistingFilePath(strDarkFileName,true);
      }

      //----------------------------------------------------------------------
      // Get the filename to write average dark frames to
      //----------------------------------------------------------------------
      if(cl->OnCommandLine("-avdark"))
      {
         if(cl->NumArgsOfOpt("-avdark")!=1)
         {
            throw CommandLine::CommandLineException("-avdark should immediately preceed the filename of the file to create. Got: "+cl->GetArg("-avdark"));      
         }     
         strAverageDarkFramesFile=cl->GetArg("-avdark"); //get the filename that will be used to write the average dark frames to 
         OUTPUT_AVERAGED_DARKFRAMES=true;
      }

      //----------------------------------------------------------------------
      // Get the filename to write binned gains to
      //----------------------------------------------------------------------
      if(cl->OnCommandLine("-gains"))
      {
         if(cl->NumArgsOfOpt("-gains")!=1)
         {
            throw CommandLine::CommandLineException("-gains should immediately preceed the filename of the file to create. Got: "+cl->GetArg("-gains"));      
         }     
         strBinnedGainsFile=cl->GetArg("-gains"); //get the filename that will be used to write the average dark frames to 
         OUTPUT_BINNED_GAINS=true;
      }

      //----------------------------------------------------------------------
      // Get any lines to be marked as corrupt data
      //----------------------------------------------------------------------
      if(cl->OnCommandLine("-corruptscans"))
      {
         if(cl->NumArgsOfOpt("-corruptscans")==0)
         {
            throw CommandLine::CommandLineException("-corruptscans should immediately preceed a list of scan lines to mark as corrupt.");      
         }           
         else
            number_of_corrupt_lines=cl->NumArgsOfOpt("-corruptscans");

         //Create array to hold corrupt line numbers
         corrupt_lines=new unsigned int[number_of_corrupt_lines];
         for(unsigned int i=0;i<number_of_corrupt_lines;i++)
         {
            corrupt_lines[i]=StringToUINT(cl->GetArg("-corruptscans",i));
         }
         Logger::Warning("These "+ToString(number_of_corrupt_lines)+"lines will be marked as corrupt and set to 0 in output file: "+cl->GetArg("-corruptscans"));
      }

      //----------------------------------------------------------------------
      // Get a filename containing pixels to be masked with QCFailure flag
      //----------------------------------------------------------------------
      if(cl->OnCommandLine("-qcfailures"))
      {
         if(cl->NumArgsOfOpt("-qcfailures")!=1)
         {
            throw CommandLine::CommandLineException("-qcfailures should immediately preceed the filename of the bad pixels to use. Got: "+cl->GetArg("-qcfailures"));      
         }
         strQCFailureFileName=cl->GetArg("-qcfailures"); 
         strQCFailureFileName=GetExistingFilePath(strQCFailureFileName,true);
         DO_APPLY_QCFAILURES=true;
      }
      else
      {
         DO_APPLY_QCFAILURES=false;
      }

      //----------------------------------------------------------------------
      // Force the use of the autodarkstartline value from the hdr file
      //----------------------------------------------------------------------
      if(cl->OnCommandLine("-darkforce"))
      {
         if(cl->NumArgsOfOpt("-darkforce")!=0)
         {
            throw CommandLine::CommandLineException("-darkforce should not have a parameter. Got: "+cl->GetArg("-darkforce"));      
         }
         DARKFORCE=true;
      }
      else
      {
         DARKFORCE=false;
      }



      //----------------------------------------------------------------------
      // Check for flags on the command line to say what processes should be run
      //----------------------------------------------------------------------
      if(cl->OnCommandLine("-NOFODIS"))
         DO_CALIBRATE_FODIS=false; //dont process the FODIS
      if(cl->OnCommandLine("-NOMASK"))
         DO_OUTPUT_MASK=false;  //dont output the mask
      if(cl->OnCommandLine("-NOMISSSCAN"))
         DO_INSERT_MISSINGSCANS=false; //dont insert missing scans
      if(cl->OnCommandLine("-NODARK"))
         DO_REMOVE_DARK=false; //dont remove average dark frames from raw data
      if(cl->OnCommandLine("-NORAD"))
         DO_APPLY_GAINS=false; //dont do the radiometric calibration

      if((((cl->OnCommandLine("-FLIPBANDS"))&&(cl->OnCommandLine("-NOFLIP")))
         || ((cl->OnCommandLine("-FLIPSAMPLES"))&&(cl->OnCommandLine("-NOFLIP")))))
      {
         throw "Cannot use -NOFLIP together with FLIPBANDS or FLIPSAMPLES";
      }
      else
      {
         //Get the values from the command line at a later point after job creation
      }

      if(cl->OnCommandLine("-NOSMEAR"))
         DO_SMEAR_CORRECT=false; //Do not do smear correction on Eagle data
   }
   catch(CommandLine::CommandLineException e)
   {
      Logger::Error(std::string(e.what())+"\n"+e.info);
      log.Flush();
      log.Close();
      delete cl;
      exit(1);
   }
   catch(std::string e)
   {
      Logger::Error(e);
      log.Flush();
      log.Close();
      delete cl;
      exit(1);
   }
   catch(std::exception& e)
   {
      Logger::Error(e.what());
      log.Flush();
      log.Close();
      delete cl;
      exit(1);
   }

   Logger::Log("Command line used to run: "+cl->ReturnCLAsString());

   //This is for the for loop later on - creating it here so that it is available
   //for the catch statements to output extra info
   unsigned int line=0;

   try
   {
      //----------------------------------------------------------------------
      //Set up a main worker object which deals with the task of giving jobs
      //to the calibration object and essentially is the "main" function
      //Give it the raw file name and the filename to save image data to
      //----------------------------------------------------------------------
      MainWorker* job=NULL;
      if(sensor == '0')
      {
         //Use the default auto detecting sensor from raw hdr file
         job=new MainWorker(strRawFileName,strOutputDataFileName,cl->ReturnCLAsString(),DARKFORCE);
      }
      else
      {
         //Force the use of a particular sensor - not guaranteed to calibrate
         //if for example eagle cal file used with hawk data etc
         job=new MainWorker(strRawFileName,strOutputDataFileName,sensor,cl->ReturnCLAsString(),DARKFORCE);         
      }

      //----------------------------------------------------------------------
      // Set the job processing item flags
      //----------------------------------------------------------------------
      job->SetTask(MainWorker::insert_missing_scans,DO_INSERT_MISSINGSCANS);
      job->SetTask(MainWorker::remove_dark_frames,DO_REMOVE_DARK);
      job->SetTask(MainWorker::smear_correct,DO_SMEAR_CORRECT);
      job->SetTask(MainWorker::apply_gains,DO_APPLY_GAINS);
      job->SetTask(MainWorker::calibrate_fodis,DO_CALIBRATE_FODIS);
      job->SetTask(MainWorker::output_mask,DO_OUTPUT_MASK);
      job->SetTask(MainWorker::apply_qcfailures,DO_APPLY_QCFAILURES);

      //If hawk and outputting mask also output the bad pixel method file
      if((DO_OUTPUT_MASK==true)&&(CheckSensorID(HAWK,job->sensor->SensorID())))
      {
         job->SetTask(MainWorker::output_mask_method,true);
      }

      if(cl->OnCommandLine("-NOFLIP"))
      {
         DO_FLIP_BANDS=false;
         DO_FLIP_SAMPLES=false;
      }
      else
      {
         if((cl->OnCommandLine("-FLIPBANDS"))||(CheckSensorID(EAGLE,job->sensor->SensorID()) ))
            DO_FLIP_BANDS=true; 

         if((cl->OnCommandLine("-FLIPSAMPLES"))||(CheckSensorID(HAWK,job->sensor->SensorID()) ))
            DO_FLIP_SAMPLES=true; 
      }
  
      job->SetTask(MainWorker::flip_bands,DO_FLIP_BANDS);
      job->SetTask(MainWorker::flip_samples,DO_FLIP_SAMPLES);
         
      //----------------------------------------------------------------------
      //Set the limits of the number of lines to process
      //if the startline and endline are 0 then process the full file
      //if not 0 then use the values in the variables
      //----------------------------------------------------------------------
      if((startline == 0)&&(endline == 0))
      {
         //We want to set endline to be the last line in the raw image to calibrate
         //i.e. do not include missing scans etc
         endline = job->sensor->br->GetNumImageFrames();
      }
      else
      {
         job->SetLineLimits(startline,endline);
         //If startline is not 0 check earlier lines for missing frames
         //this is important for syncing the data in aplnav
         if(startline!=0)
            job->SetDroppedScansPriorToStartLine(job->sensor->br->GetMissingFramesBetweenLimits(0,startline));
      }

      //----------------------------------------------------------------------
      //Initialise the calibration by passing the (optional) filename to
      //use as the calibration gains file, and (optional) dark frame file
      //this also sets up all the arrays required based on the job tasks
      //----------------------------------------------------------------------
      job->InitialiseCalibration(strCalibFileName,strDarkFileName,strQCFailureFileName);

      Logger::Log("Number of frames of image (minus dark frames) should be: "+ToString(job->sensor->br->GetNumImageFrames()));
      Logger::Log("Number of dark frames is: "+ToString(job->sensor->br->GetNumDarkFrames()));      
      Logger::Log("Number of missing frames is: "+ToString(job->sensor->br->GetTotalMissingFrames()));
      Logger::Log("\nNumber of frames in final calibrated image will be: "+ToString(job->GetNumCalibratedImageLines()));
      Logger::Log("Number of samples in final calibrated image will be: "+ToString(job->GetNumCalibratedImageSamples()));

      //----------------------------------------------------------------------
      // Do we want to output the averaged dark values 
      //----------------------------------------------------------------------      
      if(OUTPUT_AVERAGED_DARKFRAMES)
      {
         Logger::Log("Outputting average dark frames to file: "+strAverageDarkFramesFile);
         BILWriter avdark(strAverageDarkFramesFile,avdark.float32,1,job->sensor->NumSamples(),job->sensor->NumBands(),'w');
         avdark.WriteDataToLineSection(job->cal->pData()->AverageDark(),job->sensor->NumSamples(),0,job->sensor->NumSamples()-1);
         avdark.Close();
      }

      //----------------------------------------------------------------------
      // Output to the terminal the list of tasks that are to be performed
      //----------------------------------------------------------------------     
      Logger::Log("\nThe job will complete the following tasks:");
      Logger::Log(job->TasksAsString());

      //----------------------------------------------------------------------
      // Start the main processing loop - calibrate each line of data
      //----------------------------------------------------------------------
      int frameincrease=0;
      bool linecorrupt=false;//pesky bool to keep track of whether we have marked this line as corrupt
      for(line=startline;line<endline;line++)
      {
         //----------------------------------------------------------------------
         //For every line after the first one check if the frame counter is
         //increasing by 1 or more and insert missing scans if appropriate
         //only do this if missing scans are requested to be inserted
         //----------------------------------------------------------------------
         if(job->GetTask(MainWorker::insert_missing_scans))
         {
            if(line > startline)
            {
               frameincrease=job->cal->CheckFrameCounter(line-1,line);   
               if(frameincrease > 1)
               {
                  //Missing frames - need to insert (frameincrease - 1) frames
                  Logger::Log("Missing scan detected: "+ToString(frameincrease-1)+" line(s) at raw line "+ToString(line));
                  for(int i=0;i<frameincrease-1;i++)
                     job->WriteOutData(MainWorker::MissingScan);
               }
               //Allow the case where it wraps (0 - 65535) and also allow for 2 dropped frame here, else throw an error
               else if((frameincrease <= -65533) &&  (frameincrease > -65535))
               {
                  //Missing frames - need to insert frames
                  frameincrease = frameincrease + 65536;
                  Logger::Log("Missing scan detected: "+ToString(frameincrease-1)+" line(s) at raw line "+ToString(line));
                  for(int i=0;i<frameincrease-1;i++)
                     job->WriteOutData(MainWorker::MissingScan);                  
               }
               else if((frameincrease <= 0) &&  (frameincrease > -65533)) 
               {
                  //This doesn't sound too good does it.
                  throw "Frame counter has changed by: "+ToString(frameincrease)+" This seems odd?";
               }
            }
         }

         for(unsigned int c=0;c<number_of_corrupt_lines;c++)
         {
            if(line == corrupt_lines[c])
            {
               //----------------------------------------------------------------------
               //If this line is to be marked as corrupt then output 0s and skip calibration
               //----------------------------------------------------------------------
               Logger::Log("Skipping line: "+ToString(line)+" and marking it as corrupt");
               job->WriteOutData(MainWorker::CorruptData);
               linecorrupt=true;
               break;
            }
         }
         
         //----------------------------------------------------------------------
         //Calibrate this line of data and write it out
         //----------------------------------------------------------------------
         if(linecorrupt==false)
            job->DoCalibrationForLine(line);
         else
            linecorrupt=false;
      }

      //----------------------------------------------------------------------
      // Do we want to output the binned gain values
      //----------------------------------------------------------------------
      if(OUTPUT_BINNED_GAINS)
      {
         Logger::Log("Outputting binned gains to file: "+strBinnedGainsFile);
         BILWriter gains(strBinnedGainsFile,gains.float32,1,job->sensor->NumSamples(),job->sensor->NumBands(),'w');
         gains.WriteDataToLineSection(job->cal->pData()->Gains(),job->sensor->NumSamples(),0,job->sensor->NumSamples()-1);
         gains.Close();
      }

      //----------------------------------------------------------------------     
      // Close the job and finalise the output bil hdr files
      //----------------------------------------------------------------------
      delete job;
      Logger::Log("Calibration processing completed.\n");

   }
   catch(std::string e)
   {
      Logger::Error("Error on line: "+ToString(line)+"\n"+e);
   }
   catch(const char* e)
   {
      Logger::Error("Error on line: "+ToString(line)+"\n"+e);
   }
   catch(std::exception& e)
   {
      Logger::Error("Error on line: "+ToString(line)+"\n"+e.what());
   }
   catch(BinaryReader::BRexception e)
   {
      Logger::Error(std::string(e.what())+"\n"+e.info);
   }
   catch(BILWriter::BILexception e)
   {
      Logger::Error(std::string(e.what())+"\n"+e.info);
   }
   delete cl;
   if(corrupt_lines!=NULL)
      delete[] corrupt_lines;
}
