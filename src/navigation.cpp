//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Navigation.cpp - the main function for the navigation processing software
//    The main tasks of this software are to match up level 1 scan lines
//    with the post-processed navigation data. This is then fed into the 
//    geocorrection software to create level 3 imagery.
//-------------------------------------------------------------------------

#include <cstdlib>
#include <cmath>
#include <iostream>

#include "navigationsyncer.h"
#include "navigationinterpolator.h"
#include "commandline.h"
#include "logger.h"
#include "leverbore.h"
#include "os_dependant.h"


//-------------------------------------------------------------------------
// Software description
//-------------------------------------------------------------------------
const std::string DESCRIPTION="Navigation Interpolation Software";

//----------------------------------------------------------------
//Number of options that can be on command line
//----------------------------------------------------------------
const int number_of_possible_options = 14;

//----------------------------------------------------------------
//Option names that can be on command line
//----------------------------------------------------------------
const std::string availableopts[number_of_possible_options]={
"-procnav",
"-nav",
"-output",
"-lev1",
"-scantimeoffset",
"-leverarm",
"-boresight",
"-smooth",
"-interp",
"-posattoff",
"-nonav",
"-qualityfile",
"-force",
"-help"
}; 

//----------------------------------------------------------------
//One line description of each option - in same order as availableopts
//----------------------------------------------------------------
const std::string optsdescription[number_of_possible_options]={
"Post-processed navigation data file either in SBET or SOL format",
"Real time Specim navigation file",
"Output BIL filename to save navigation to",
"Level-1 data BIL file",
"Offset to apply to scan times, shifts the navigation data w.r.t the scan lines (default 0 seconds)",
"Lever arm corrections to apply: X Y Z",
"Boresight corrections to apply: Pitch Roll Heading",
"Smooth the input data using a Triangular filter of kernel length x (default is no smoothing)",
"Method of interpolation to use, either 'Linear' or 'Spline' (default is Linear)",
"Offset to apply (in seconds) to offset the Position and Attitude data by",
"If no Specim navigation file exists for this line.",
"An optional BIL filename to output the quality flags to for the navigation.",
"Force the processing when 'time goes backwards' in a navigation file (only use without processed nav when the data is not used for further processing). DO NOT USE FOR TYPICAL DATA PROCESSING.",
"Display this help"
}; 

//----------------------------------------------------------------
//Nasty global define to force processing in various situations. Currently:
// *  Time goes backwards in navigation file (may want to continue processing if specim nav file corrupt)
//----------------------------------------------------------------
bool GLOBAL_FORCE=false;

int main(int argc,char* argv[])
{
   //Set the precision of the data output to the terminal
   std::cout.precision(10);

   //A command line handling object
   CommandLine* cl=NULL;
   //Create a logger to output the data to terminal
   Logger log;

   //A stringstream object to store output information in
   std::stringstream strout;

   //Specim navigation filename
   std::string strSpecimNavFile="";
   //Output filename to write results to
   std::string strOutputFile="";
   //Level-1 data file name to read data properties from
   std::string strLevel1File="";
   //Post-processed nav file name 
   std::string strPostProcNavFile="";
   //Pointer to a NavigationInterpolator
   NavigationInterpolator* interpolator=NULL;
   //Pointer to a boresight onject
   Boresight* boresight=NULL;
   //Pointer to a lever arm object
   Leverarm* leverarm=NULL;
   //String to hold interpolation choice
   std::string strinterpmethod="";
   //Value to hold user specified time to use as first scan line time
   //double usersynctime=-1;

   //Scan offset time
   double scantimeoffset=0;
   //Position-Attitude time offset
   double posattoffset=0;
   //Smoothing filter kernel size
   unsigned int smoothkernelsize=0;

   //Shall we write out the navigation quality flags
   bool WRITE_QUALITY=false;
   std::string strOutputFlagFile;

   //String to contain info for output hdr file
   std::string info="";

   //Get exe name without the path
   std::string niceexename=std::string(argv[0]);
   niceexename=niceexename.substr(niceexename.find_last_of("/\\")+1);

   //Output a "nice" header
   Logger::FormattedInformation(niceexename,VERSION,DESCRIPTION);

   try
   {
      //----------------------------------------------------------------------------
      //Create the command line object
      //----------------------------------------------------------------------------
      cl=new CommandLine(argv,argc);

      //----------------------------------------------------------------------------
      //Check if anything went wrong (most likely an exception will have been thrown but lets be safe)
      //----------------------------------------------------------------------------
      if(cl->IsGood()==false)
      {
         throw "An error has occurred with the command line\n"; //throw exception that an error has occurred
      }    

      //----------------------------------------------------------------------
      // Get some information about the system running the program and log it
      //----------------------------------------------------------------------
      ComputerInfo cinfo;
      Logger::Log(cinfo.GetOutput());

      //----------------------------------------------------------------------------
      //check the command line options given against the list of available options in the software
      //----------------------------------------------------------------------------
      std::string badoptions;
      int retval=cl->CheckAvailableOptions(availableopts,number_of_possible_options,&badoptions); 
      if(retval<0)
      {
         //Options are on commnd line that are not available in this software
         strout<<"There are "<<(-retval)<<" unrecognised options on command line: "<<badoptions<<std::endl;
         throw CommandLine::CommandLineException(strout.str());
      }

      Logger::Log("Command line used to run: "+cl->ReturnCLAsString());

      //----------------------------------------------------------------------------
      // Go through each possible option in turn and set up the required data / response
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-help"))
      {
         log.Add(cl->ProgramUsage(number_of_possible_options,availableopts,optsdescription));
         log.Flush();
         throw "";
      }
   
      //-------------------------------------------------------------------
      // Is the Specim .nav file on the command line (IT MUST BE PRESENT)
      //-------------------------------------------------------------------
      if((cl->OnCommandLine("-nav"))&&(!cl->OnCommandLine("-nonav")))
      {
         //Check that an argument follows the nav option - and get it if it exists
         if(cl->GetArg("-nav").compare(optiononly)!=0)
         {
            strSpecimNavFile=cl->GetArg("-nav");
            log.Add("Will use Specim Navigation file: "+strSpecimNavFile);
         }
         else
            throw CommandLine::CommandLineException("Argument -nav must immediately precede the Specim .nav filename.\n");
      }
      else if((!cl->OnCommandLine("-nav"))&&(cl->OnCommandLine("-nonav")))
      {
         //Throw an exception
         //throw CommandLine::CommandLineException("Argument -nav [Specim navigation file] must be present on the command line.\n");  
         if(cl->GetArg("-nonav").compare(optiononly)!=0)
         {
            throw CommandLine::CommandLineException("Option -nonav does not take any arguments.\n");
         }
         else
            strSpecimNavFile="NULL";
      }
      else
      {
         //Throw an exception
         throw CommandLine::CommandLineException("Argument -nav [Specim navigation file] or -nonav must be present on the command line.\n");  
      }

      //-------------------------------------------------------------------
      // Is the output filename on the command line (IT MUST BE PRESENT)
      //-------------------------------------------------------------------
      if(cl->OnCommandLine("-output"))
      {
         //Check that an argument follows the output option - and get it if it exists
         if(cl->GetArg("-output").compare(optiononly)!=0)
         {
            strOutputFile=cl->GetArg("-output");
            log.Add("Will write to output BIL file: "+strOutputFile);
         }
         else
            throw CommandLine::CommandLineException("Argument -output must immediately precede the output filename.\n");
      }
      else
      {
         //Throw an exception
         throw CommandLine::CommandLineException("Argument -output [the file to write to] must be present on the command line.\n");  
      }

      //-------------------------------------------------------------------
      // Is the level 1 filename on the command line (IT MUST BE PRESENT)
      //-------------------------------------------------------------------
      if(cl->OnCommandLine("-lev1"))
      {
         //Check that an argument follows the lev1 option - and get it if it exists
         if(cl->GetArg("-lev1").compare(optiononly)!=0)
         {
            strLevel1File=cl->GetArg("-lev1");
            log.Add("Will read the sensor data properties from the level-1 file: "+strLevel1File);
         }
         else
            throw CommandLine::CommandLineException("Argument -lev1 must immediately precede the level-1 filename.\n");
      }
      else
      {
         //Throw an exception
         throw CommandLine::CommandLineException("Argument -lev1 [the level-1 hyperspectral file] must be present on the command line.\n");  
      }

      //-------------------------------------------------------------------
      // Is the lever arm on the command line (IT MUST BE PRESENT)
      //-------------------------------------------------------------------
      if(cl->OnCommandLine("-leverarm"))
      {
         //Check that 3 arguments follow the leverarm option
         if(cl->NumArgsOfOpt("-leverarm") != 3 )
            throw CommandLine::CommandLineException("Error: There should be 3 arguments following the -leverarm option.\n");
         
         //Get the 3 arguments into the leverarm parameters
         double lax=StringToDouble(cl->GetArg("-leverarm",0));
         double lay=StringToDouble(cl->GetArg("-leverarm",1));
         double laz=StringToDouble(cl->GetArg("-leverarm",2));

         log.Add("Will apply lever arm corrections of (X,Y,Z): "+ToString(lax)+" "+ToString(lay)+" "+ToString(laz));
         leverarm=new Leverarm(lax,lay,laz);

      }
      else
      {
         throw CommandLine::CommandLineException("Argument -leverarm [the sensor lever arm values] must be present on the command line.\n");  
      }

      //-------------------------------------------------------------------
      // Is the boresight on the command line (IT MUST BE PRESENT)
      //-------------------------------------------------------------------
      if(cl->OnCommandLine("-boresight"))
      {
         //Check that 3 arguments follow the boresight option
         if(cl->NumArgsOfOpt("-boresight") != 3 )
            throw CommandLine::CommandLineException("Error: There should be 3 arguments following the -boresight option.\n");
         
         //Get the 3 arguments into the boresight parameters
         double p=StringToDouble(cl->GetArg("-boresight",0));
         double r=StringToDouble(cl->GetArg("-boresight",1));
         double h=StringToDouble(cl->GetArg("-boresight",2));
         log.Add("Will apply boresight corrections of (R,P,H): "+ToString(r)+" "+ToString(p)+" "+ToString(h));
         boresight=new Boresight(r,p,h);
      }
      else
      {
         throw CommandLine::CommandLineException("Argument -boresight [the sensor boresight values] must be present on the command line.\n");  
      }

      //-------------------------------------------------------------------
      // Is the post-processed nav file on the command line (ITS OPTIONAL)
      //-------------------------------------------------------------------
      if(cl->OnCommandLine("-procnav"))
      {
         //Check that an argument follows the procnav option - and get it if it exists
         if(cl->GetArg("-procnav").compare(optiononly)!=0)
         {
            strPostProcNavFile=cl->GetArg("-procnav");
            log.Add("Will read navigation data from SBET/SOL file: "+strPostProcNavFile);
            info=info+";Navigation from post-processed SBET/SOL file. \n";
         }
         else
            throw CommandLine::CommandLineException("Argument -procnav must immediately precede the SBET/SOL filename.\n");
      }
      else
      {            
         log.Add("No SBET/SOL file has been given, therefore will read real-time navigation data from Specim .nav file: "+strSpecimNavFile);         
         info=info+";Navigation from real-time Specim .nav file. \n";
      }

      //-------------------------------------------------------------------
      // Is the ScanTime Offset on the command line (ITS OPTIONAL)
      //-------------------------------------------------------------------
      if(cl->OnCommandLine("-scantimeoffset"))
      {
         //Check that an argument follows the scantimeoffset option - and get it if it exists
         if(cl->GetArg("-scantimeoffset").compare(optiononly)!=0)
         {
            std::string cl_tmp=cl->GetArg("-scantimeoffset");
            scantimeoffset=StringToDouble(cl_tmp.c_str());
            if(errno==ERANGE)
            {
               throw "An error has occurred in the conversion of scan time offset in CommandLine of main().";
            }
            log.Add("Will apply a user-specified scan offset of: "+cl_tmp);
         }
         else
            throw CommandLine::CommandLineException("Argument -scantimeoffset must immediately precede the scan time offset value.");
      }
      else
      {
         log.Add("No user-supplied scan time offset to be applied.");
      }
      
      //-------------------------------------------------------------------
      // Is smoothing requested on the command line (ITS OPTIONAL)
      //-------------------------------------------------------------------
      if(cl->OnCommandLine("-smooth"))
      {
         //Check that an argument follows the smooth option - and get it if it exists
         if(cl->GetArg("-smooth").compare(optiononly)!=0)
         {
            std::string cl_tmp=cl->GetArg("-smooth");
            smoothkernelsize=StringToUINT(cl_tmp);
            if(errno==ERANGE)
            {
               throw "An error has occurred in the conversion of smooth kernel size in CommandLine of main().";
            }
            if(smoothkernelsize%2==0)
               throw "Smoothing kernel size must be an odd number.";
               
            log.Add("Will apply a smoothing of kernel size: "+cl_tmp);
         }
         else
            throw CommandLine::CommandLineException("Argument -smooth must immediately precede the smoothing kernel size.");
      }
      else
      {
         log.Add("No smoothing of navigation data to be applied.");
      }

      //-------------------------------------------------------------------
      // What method of interpolation is to be used
      //-------------------------------------------------------------------
      if(cl->OnCommandLine("-interp"))
      {
         //Check that an argument follows the interp option - and get it if it exists
         if(cl->GetArg("-interp").compare(optiononly)!=0)
         {
            std::string cl_tmp=cl->GetArg("-interp");
            strinterpmethod=cl_tmp;

            log.Add("Will use the interpolation method: "+cl_tmp);
         }
         else
            throw CommandLine::CommandLineException("Argument -interp must immediately precede the interpolation method keyword.");
      }
      else
      {
         strinterpmethod="Linear";
         log.Add("No interpolation method of navigation data supplied, will use Linear.");
      }

      //-------------------------------------------------------------------
      // Add a Position-Attitude time offset 
      //-------------------------------------------------------------------
      if(cl->OnCommandLine("-posattoff"))
      {
         //Check that an argument follows the posattoff option - and get it if it exists
         if(cl->GetArg("-posattoff").compare(optiononly)!=0)
         {
            std::string cl_tmp=cl->GetArg("-posattoff");
            posattoffset=StringToDouble(cl_tmp.c_str());
            if(errno==ERANGE)
            {
               throw "An error has occurred in the conversion of position attitude time offset in CommandLine of main().";
            }
            log.Add("Will use the given position attitude offset: "+cl_tmp);
         }
         else
            throw CommandLine::CommandLineException("Argument -posattoff must immediately precede the position/attitude shift value.");
      }
      else
      {
         posattoffset=0; 
         log.Add("Will not use a position-attitude offset.");
      }


      //-------------------------------------------------------------------
      // Output the quality flags
      //-------------------------------------------------------------------
      if(cl->OnCommandLine("-qualityfile"))
      {
         //Check that an argument follows the qualityfile option - and get it if it exists
         if(cl->GetArg("-qualityfile").compare(optiononly)!=0)
         {
            strOutputFlagFile=cl->GetArg("-qualityfile");
            log.Add("Will write quality flags to: "+strOutputFlagFile);
         }
         else
            throw CommandLine::CommandLineException("Argument -qualityfile must immediately precede the output quality flag filename.\n");

         //Set the bool to true
         WRITE_QUALITY=true;
      }
      else
      {
         WRITE_QUALITY=false;
      }

      //-------------------------------------------------------------------
      // Force processing of aplnav
      //-------------------------------------------------------------------
      if(cl->OnCommandLine("-force"))
      {
         if(cl->GetArg("-force").compare(optiononly)!=0)
         {
            throw CommandLine::CommandLineException("Option -force does not take any arguments.\n");
         }
         else
            GLOBAL_FORCE=true;
      }
      else
      {
         GLOBAL_FORCE=false;
      }

      //-------------------------------------------------------------------
      // ENTER NEW COMMAND LINE OPTIONS HERE
      //-------------------------------------------------------------------

   }
   catch(CommandLine::CommandLineException e)
   {
      Logger::Error(std::string(e.what())+"\n"+e.info);
      delete cl;
      exit(1);
   }
   catch(std::string e)
   {
      Logger::Error(e);
      delete cl;
      exit(1);
   }
   catch(const char* e)
   {
      Logger::Error(e);
      delete cl;
      exit(1);
   }
   catch(BinaryReader::BRexception e)
   {
      Logger::Error(std::string(e.what())+"\n"+e.info);
      delete cl;
      exit(1); //exit the program
   }
   catch(std::exception &e)
   {
      Logger::Error(e.what());
      delete cl;
      exit(1);
   }

   //Flush the log
   log.Flush();
   //Output a blank line
   Logger::Log("");

   //Add some information to the header files
   info=info+";Command line used to process data: "+cl->ReturnCLAsString()+"\n";
   info=info+";boresight (P,R,H) = "+ToString(boresight->Pitch())+" "+ToString(boresight->Roll())+" "+ToString(boresight->Heading())+"\n";
   info=info+";leverarm (X,Y,Z) = "+ToString(leverarm->X())+" "+ToString(leverarm->Y())+" "+ToString(leverarm->Z())+"\n";


   try
   { 
      //-------------------------------------------------------------------
      // In this section we deal with getting per scan times
      //-------------------------------------------------------------------

      //Set up navigation syncer to get per scan line times
      NavigationSyncer syncer(strSpecimNavFile,strLevel1File);

      //Add the y start value to the output header file.
      info=info+"y start = "+ToString(syncer.GetCropTimeOffset())+"\n";

      //Get the perscan line times
      log.Add("Finding per-scan times...");
      log.Flush();
      syncer.FindScanTimes();

      //Correct the times for GPS leapseconds if from SBET/SOL file
      //This is because the times from the Specim nav file are in GPS time
      //and SBET/SOL times are in UTC. So this adds leapseconds onto the times 
      //which have been derived from the Specim data, making them relevant
      //for querying the data from the SBET/SOL file.
      //Ignore (for the moment) if from Specim Nav file - adds time on
      //to that data before writing out (see below)
      if(strPostProcNavFile!="")
         syncer.ApplyLeapSeconds();

      //Apply a user-defined time offset to the scans
      if(scantimeoffset!=0)
      {
         log.Add("\nApplying user defined timing offset...");
         log.Flush();
         syncer.ApplyTimeShift(scantimeoffset);
         info=info+";User defined scan timing offset added onto data: "+ToString(scantimeoffset)+"\n";
      }

      //-------------------------------------------------------------------
      // In this section we deal with getting per scan navigation data
      //-------------------------------------------------------------------

      //Set up an interpolator to interpolate the navigation to the scan lines
      log.Add("Creating Navigation Interpolation object...");
      log.Flush();
      if(strPostProcNavFile!="")
         interpolator=new NavigationInterpolator(strPostProcNavFile,strLevel1File);
      else
         interpolator=new NavigationInterpolator(strSpecimNavFile,strLevel1File);         

      //Assign the scan times to the interpolator - these are just 
      //pointing to the data so dont delete the syncer
      log.Add("\nSetting times to interpolation object...");
      log.Flush();
      double* scantimes=NULL;
      scantimes=syncer.PtrToTimes();
      interpolator->SetTimes(scantimes);

      //Smooth the nav data
      if(smoothkernelsize!=0)
      {
         log.Add("Smoothing the data using a triangular low-pass filter...");
         log.Flush();
         info=info+";Smoothed input navigation data using a triangular low-pass filter with kernel size: "+ToString(smoothkernelsize)+"\n";
         interpolator->SmoothNavData(Triangle,smoothkernelsize);
      }

      //Interpolate the data to the scan times
      log.Add("\nInterpolating navigation data to scan times...");
      log.Flush();
      if(strinterpmethod.compare("Linear")==0)
      {
         interpolator->Interpolate(Linear);
         if(posattoffset!=0)
            interpolator->PosAttShift(Linear,posattoffset);
      }
      else if(strinterpmethod.compare("Spline")==0)
      {
         interpolator->Interpolate(CubicSpline);
         if(posattoffset!=0)
            interpolator->PosAttShift(Linear,posattoffset);
      }
      else
         throw "Unknown interpolation method. Expected 'Linear' or 'Spline'";

      //Apply the lever arm offsets to the interpolated position data
      log.Add("Adding leverarm correction...");
      log.Flush();
      interpolator->ApplyLeverarm(leverarm);    

      //Apply the boresight offsets to the interpolated attitude data
      log.Add("Adding boresight correction...");
      log.Flush();
      interpolator->ApplyBoresight(boresight);

      //If the times are from the Specim nav file then they must be converted
      //from GPS time to UTC by adding on the leap seconds.
      if(strPostProcNavFile=="")
         syncer.ApplyLeapSeconds();
      scantimes=syncer.PtrToTimes();
      interpolator->SetTimes(scantimes);

      //Check the plausibilty of the interpolated data
      interpolator->CheckPlausibility();

      //Write out the data
      log.Add("\nWriting data out...");
      log.Flush();
      interpolator->Writer(strOutputFile,info);

      //Write out the quality flags
      if(WRITE_QUALITY==true)
      {
         interpolator->WriteFlags(strOutputFlagFile);
      }

      //delete the interpolator
      delete interpolator;

   }
   catch(std::string e)
   {
      Logger::Error(e);
      exit(1);
   }
   catch(const char* e)
   {
      Logger::Error(e);
      exit(1);
   }
   catch(BinaryReader::BRexception e)
   {
      Logger::Error(std::string(e.what())+"\n"+e.info);
      exit(1); 
   }
   catch(BILWriter::BILexception e)
   {
      Logger::Error(std::string(e.what())+"\n"+e.info);
      exit(1); 
   }
   catch(std::exception &e)
   {
      Logger::Error(e.what());
      exit(1);
   }

   Logger::Log("Navigation processing completed. \n \n");

   //Delete the command line object
   if(cl!=NULL)
      delete cl;
   if(boresight!=NULL)
      delete boresight;
   if(leverarm!=NULL)
      delete leverarm;
}

