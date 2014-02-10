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
#include <cmath>

#include "commandline.h"
#include "navbaseclass.h"
#include "datahandler.h"
#include "interpolationfunctions.h"
#include "bilwriter.h"
#include "logger.h"

//-------------------------------------------------------------------------
// Software description
//-------------------------------------------------------------------------
const std::string DESCRIPTION="Apply time offset to APL navigation file";


//----------------------------------------------------------------
//Number of options that can be on command line
//----------------------------------------------------------------
const int number_of_possible_options = 6;

//----------------------------------------------------------------
//Option names that can be on command line
//----------------------------------------------------------------
const std::string availableopts[number_of_possible_options]={
"-nav",
"-output",
"-fps",
"-timeoffset",
"-lev1file",
"-help"
}; 

//----------------------------------------------------------------
//One line description of each option - in same order as availableopts
//----------------------------------------------------------------
const std::string optsdescription[number_of_possible_options]={
"7-Band APL Navigation BIL input data filename",
"Output data filename.",
"Input data frame rate in frames per second",
"The time offset to apply to the navigation data.",
"The level-1 filename. If given, a new level-1 file is created that contains the data trimmed to fit the new navigation data coverage.",
"Display this help."
}; 

//----------------------------------------------------------------
// This is not used it is a horrible "fix" due to having a global
// variable in aplnav which shares cpp files with aplshift
// FIXME TODO REMOVE THIS AND FIX THE GLOBAL IN APLNAV TO WORK!!
//----------------------------------------------------------------
bool GLOBAL_FORCE=false;

//----------------------------------------------------------------
// Class for reading and handling the navigation data that is
// already in APL format
//----------------------------------------------------------------
class APLNavFile : public DataHandler
{
public:
   APLNavFile(std::string strNavigationFile);
   void Reader();
private:
};

//----------------------------------------------------------------
// Constructor for APLNavFile
//----------------------------------------------------------------
APLNavFile::APLNavFile(std::string strNavigationFile)
{
   //Get the number of entries from the file,
   //and set up the navdata array
   navcollection=NULL;
   //Open an APL navigation file for getting number of scans
   NavBaseClass inNav(strNavigationFile);

   //Create array of nav data lines
   navcollection=new NavDataCollection(inNav.TotalScans());

   //Copy the filename over
   this->filename=strNavigationFile;
}

//----------------------------------------------------------------
// Reader for APLNavFile - essentially a wrapper for NavBaseClass
// reading function
//----------------------------------------------------------------
void APLNavFile::Reader()
{
   //Read in an APL navigation file
   NavBaseClass inNav(this->filename);

   for(unsigned int record=0;record<inNav.TotalScans();record++)
   {
      //Read the scan of data
      inNav.ReadScan(record);
      //Fill in the navdata array
      navcollection->SetValues(record,NavDataCollection::TIME,inNav.Time());
      navcollection->SetValues(record,NavDataCollection::LAT,inNav.Lat());
      navcollection->SetValues(record,NavDataCollection::LON,inNav.Lon());
      navcollection->SetValues(record,NavDataCollection::HEI,inNav.Hei());
      navcollection->SetValues(record,NavDataCollection::ROLL,inNav.Roll());
      navcollection->SetValues(record,NavDataCollection::PITCH,inNav.Pitch());
      navcollection->SetValues(record,NavDataCollection::HEADING,inNav.Heading()); 
   }
}

//-------------------------------------------------------------------------
// Function that will read the level-1 data and write out to a new file.
// It only reads/writes the lines that fall within the inclusive bounds
// given by start and end
//-------------------------------------------------------------------------
template<class T>
void ReadWriteData(BinFile* bilIn,BILWriter* bilOut,unsigned int start,unsigned int end)
{
   unsigned int inBands=StringToUINT(bilIn->FromHeader("bands"));
   unsigned int inSamples=StringToUINT(bilIn->FromHeader("samples"));
   //Read in from the required start line until the end line, writing the data out
   T* databuffer=new T[inBands*inSamples];
   for(unsigned int line=start;line<=end;line++)
   {
      bilIn->Readline((char*)databuffer,line);
      bilOut->WriteLine((char*)databuffer);
   }
   delete[] databuffer;
}

//-------------------------------------------------------------------------
// Function called to trim the level-1 data to match the coverage of the
// new navigation file. This is an optional procedure.
//-------------------------------------------------------------------------
void TrimLevel1Data(std::string strLevel1File,std::string strTrimmedLevel1File,unsigned int start,unsigned int end)
{
   //Open the level-1 input file
   BinFile bilIn(strLevel1File);
   //Create an output file
   unsigned int inBands=StringToUINT(bilIn.FromHeader("bands"));
   unsigned int inSamples=StringToUINT(bilIn.FromHeader("samples"));
   BILWriter* bilOut=NULL;
   switch(bilIn.GetDataType())
   {
   case 1:
      bilOut=new BILWriter(strTrimmedLevel1File,BILWriter::char8,(end-start+1),inSamples,inBands,'w');
      ReadWriteData<char>(&bilIn,bilOut,start,end);
      break;
   case 2:
      bilOut=new BILWriter(strTrimmedLevel1File,BILWriter::int16,(end-start+1),inSamples,inBands,'w');
      ReadWriteData<short int>(&bilIn,bilOut,start,end);
      break;
   case 3:
      bilOut=new BILWriter(strTrimmedLevel1File,BILWriter::int32,(end-start+1),inSamples,inBands,'w');
      ReadWriteData<int>(&bilIn,bilOut,start,end);
      break;
   case 4:
      bilOut=new BILWriter(strTrimmedLevel1File,BILWriter::float32,(end-start+1),inSamples,inBands,'w');
      ReadWriteData<float>(&bilIn,bilOut,start,end);
      break;
   case 5:
      bilOut=new BILWriter(strTrimmedLevel1File,BILWriter::float64,(end-start+1),inSamples,inBands,'w');
      ReadWriteData<double>(&bilIn,bilOut,start,end);
      break;
   case 12:
      bilOut=new BILWriter(strTrimmedLevel1File,BILWriter::uint16,(end-start+1),inSamples,inBands,'w');
      ReadWriteData<unsigned short int>(&bilIn,bilOut,start,end);
      break;
   case 13:
      bilOut=new BILWriter(strTrimmedLevel1File,BILWriter::uint32,(end-start+1),inSamples,inBands,'w');
      ReadWriteData<unsigned int>(&bilIn,bilOut,start,end);
      break;
   default:
      throw "Unrecognised data type in Level-1 file header.";
      break;
   }

   //Also need to copy the hdr items over - excluding certain ones
   std::map<std::string, std::string, cmpstr> header=bilIn.CopyHeaderExcluding();
   for(std::map<std::string, std::string, cmpstr>::iterator iter=header.begin();iter!=header.end();iter++)
   {
      std::string tmp="";
      if((*iter).first.at(0)!=';')
      {
         if((start != 0)&&((*iter).first.compare("y start")==0))
         {
            //Then we need to update the y start value in the hdr file
            unsigned int original=StringToUINT((*iter).second);
            (*iter).second=ToString(original+start);
         }
         tmp=(*iter).first+" = "+(*iter).second;
      }
      else
         tmp=(*iter).first; //a comment - does not have a second

      bilOut->AddToHdr(bilIn.TidyForHeader(tmp));
   }

   //Close the files
   bilIn.Close();
   bilOut->Close();
}

//----------------------------------------------------------------
// Main function
//----------------------------------------------------------------
int main (int argc, char* argv[])
{
   std::cout.precision(10);
   //A command line handling object
   CommandLine* cl=NULL;
   std::stringstream strout;

   std::string strNavigationFile="";
   std::string strOutputNavFile="";
   std::string strTrimmedLevel1File="";
   std::string strLevel1File="";

   double time_offset=0;
   double fps=0;
   double scan_separation=0;

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
         throw "An error has occurred with the command line\n";
      }    

      //----------------------------------------------------------------------------
      //check the command line options given against the list of available options in the software
      //----------------------------------------------------------------------------
      std::string badoptions;
      int retval=cl->CheckAvailableOptions(availableopts,number_of_possible_options,&badoptions); 
      if(retval<0)
      {
         //Options are on commnd line that are not available in this software
         strout<<"There are "<<(-retval)<<" unrecognised options on command line: "<<badoptions;
         throw CommandLine::CommandLineException(strout.str());
      }

      //----------------------------------------------------------------------------
      // Go through each possible option in turn and set up the required data / response
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-help"))
      {
         Logger::Log(cl->ProgramUsage(number_of_possible_options,availableopts,optsdescription));
         throw "";
      }

      Logger::Log("Command line used to run: "+cl->ReturnCLAsString());

      //----------------------------------------------------------------------------
      // Get the input APL navigation file name
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-nav"))
      {
         //Check that an argument follows the nav option - and get it if it exists
         if(cl->GetArg("-nav").compare(optiononly)!=0)
         {
            strNavigationFile=cl->GetArg("-nav");
            //This is an existing file - use GetExistingFilePath to check file status
            strNavigationFile=GetExistingFilePath(strNavigationFile,true);
            Logger::Log("Will use input nav BIL file: "+strNavigationFile);
         }
         else
            throw CommandLine::CommandLineException("Argument -nav must immediately precede the APL navigation filename.\n");         
      }
      else
      {
         //Throw an exception
         throw CommandLine::CommandLineException("Argument -nav [the Navigation BIL file] must be present on the command line.\n");  
      }

      //----------------------------------------------------------------------------
      // Get the output navigation file name
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-output"))
      {
         //Check that an argument follows the output option - and get it if it exists
         if(cl->GetArg("-output").compare(optiononly)!=0)
         {
            strOutputNavFile=cl->GetArg("-output");
            Logger::Log("Will write to BIL file: "+strOutputNavFile);
         }
         else
            throw CommandLine::CommandLineException("Argument -output must immediately precede the output navigation filename.\n");         
      }
      else
      {
         //Throw an exception
         throw CommandLine::CommandLineException("Argument -output [the output Navigation BIL file] must be present on the command line.\n");  
      }
      //----------------------------------------------------------------------------
      // Get the time offset to apply
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-timeoffset"))
      {
         //Check that an argument follows the timeoffset option - and get it if it exists
         if(cl->GetArg("-timeoffset").compare(optiononly)!=0)
         {
            time_offset=StringToDouble(cl->GetArg("-timeoffset"));
            Logger::Log("Will apply time offset: "+ToString(time_offset));
         }
         else
            throw CommandLine::CommandLineException("Argument -timeoffset must immediately precede the time offst value.\n");         
      }
      else
      {
         //Throw an exception
         throw CommandLine::CommandLineException("Argument -timeoffset must be present on the command line.\n");  
      }
      //----------------------------------------------------------------------------
      // Get the frame rate of the data
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-fps"))
      {
         //Check that an argument follows the fps option - and get it if it exists
         if(cl->GetArg("-fps").compare(optiononly)!=0)
         {
            fps=StringToDouble(cl->GetArg("-fps"));
            Logger::Log("Will use a frame rate of: "+ToString(fps));
         }
         else
            throw CommandLine::CommandLineException("Argument -fps must immediately precede the frame rate value.\n");         
      }
      else
      {
         fps=0; //default - means read it from file or calculate from data
      }
      //----------------------------------------------------------------------------
      // Get the level-1 file name
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-lev1file"))
      {
         //Check that an argument follows the lev1file option - and get it if it exists
         if(cl->GetArg("-lev1file").compare(optiononly)!=0)
         {
            strLevel1File=cl->GetArg("-lev1file");
            strLevel1File=GetExistingFilePath(strLevel1File,true);
            //Create a new level-1 output filename for the trimmed data
            strTrimmedLevel1File=strLevel1File+"_trimmed_"+ToString(time_offset)+".bil";
            Logger::Log("Will create a new level-1 file using the image data from: "+strLevel1File+" and write it to "+strTrimmedLevel1File);
         }
         else
            throw CommandLine::CommandLineException("Argument -lev1file must immediately precede the level-1 data filename.\n");         
      }
      else
      {
         strLevel1File=""; //default - do not output a new level 1 file
      }
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
   catch(char const* e)
   {
      Logger::Error(e);
      delete cl;
      exit(1);      
   }
   catch(BinaryReader::BRexception e)
   {
      Logger::Error(std::string(e.what())+"\n"+e.info);
      delete cl;
      exit(1);
   }
   catch(std::exception& e)
   {
      Logger::Error(e.what());
      delete cl;
      exit(1);
   }


   //Read in an APL navigation file
   //Create an object to hold navdata
   APLNavFile *inNav=NULL;
   try
   {
      inNav= new APLNavFile(strNavigationFile);
      inNav->Reader();
   }
   catch(std::string e)
   {
      Logger::Error(e);
      exit(1);
   }
   catch(char const* e)
   {
      Logger::Error(e);
      exit(1);      
   }
   catch(BinaryReader::BRexception e)
   {
      Logger::Error(std::string(e.what())+"\n"+e.info);
      exit(1); 
   }
   catch(std::exception& e)
   {
      Logger::Error(e.what());
      exit(1);
   }

   //Get the total number of scans in the nav data
   unsigned int number_of_scans=inNav->GetNumEntries();

   //get the start and end times of the nav data
   NavDataLine* navdata=NULL;
   navdata=inNav->GetLine(0);
   double start_time=navdata->time;

   navdata=inNav->GetLine(number_of_scans-1);
   double end_time=navdata->time;

   //Get the scan separation (= 1 / fps)
   if(fps != 0)
   {
      scan_separation=1/fps;
   }
   else
   {
      //We need to get the information from somewhere
      //Calculate an estimate from the data
      //use number_of_scans-1 as start and end times are inclusive (i.e. there are n-1 increments between them)
      scan_separation=(end_time - start_time) / (number_of_scans-1);
   }

   Logger::Log("\nStart and end times of original (input) navigation file: "+ToString(start_time)+" "+ToString(end_time));
   Logger::Log("Number of scans of original navigation: "+ToString(number_of_scans));
   Logger::Log("Separation (in seconds) of scans of original navigation: "+ToString(scan_separation) +" Equivalent to (frames per second): "+ToString(1/scan_separation));

   //Take a time offset and apply it to the data to shift it
   //Get the new start time (with the offset applied). Note this is the minimum start time
   //due to offset being either +/- and there may not be the data to cover this period
   double offset_start_time=start_time + time_offset;

   //Now find the time for the first scan of the new nav data. This is either the
   //offset_start_time if this is > start_time else it is the first time that is
   //greater than start_time that is the offset_start_time + multiple of scan separation
   double first_scan_time=offset_start_time;
   unsigned int i=0;
   while(first_scan_time < start_time)
   {
      i++;
      first_scan_time=offset_start_time+i*scan_separation;
   }

   //If time_offset is -ve then we lose the data at the start of the image
   //else if it is +ve we lose data at the end
   unsigned int number_of_lost_scans=ceil(fabs(time_offset) / scan_separation);
   unsigned int number_of_offset_scans=number_of_scans - number_of_lost_scans;

   Logger::Log("\nStart and end times of offset (output) navigation file: "+ToString(first_scan_time)+" "+ToString(first_scan_time+(number_of_offset_scans-1)*scan_separation));
   Logger::Log("Number of scans of output navigation file: "+ToString(number_of_offset_scans));
   Logger::Log("Number of lost image scans (as there will be no nav data for these): "+ToString(number_of_lost_scans));

   //Output which level-1 image data scan lines the start and end times relate to
   unsigned int level1_startscan=0,level1_endscan=0;
   if(time_offset > 0)
   {
      level1_startscan=0;
      level1_endscan=number_of_offset_scans-1;
   }
   else
   {
      level1_startscan=number_of_lost_scans;
      level1_endscan=number_of_scans-1;
   }
   Logger::Log("\nStart time is for level-1 image scan line (0 based): "+ToString(level1_startscan)+" \nEnd time is for level-1 image scan line (0 based): "+ToString(level1_endscan));

   //Create an array to store the interpolated navigation in
   NavDataCollection interpolated_nav(number_of_offset_scans);
   //Create an array to store the interpoalted times in
   double* offset_times=new double[number_of_offset_scans];

   //Fill the times array with the new times for each scan using the scan separation and start time
   for(unsigned int scan=0;scan<number_of_offset_scans;scan++)
   {
      offset_times[scan]=first_scan_time+scan*scan_separation;
   }

   //Interpolate data to new scan times
   Logger::Log("Interpolating data to new times...");
   try
   {
      Linear(offset_times,number_of_offset_scans,inNav,&interpolated_nav,"","");
   }
   catch(std::string e)
   {
      Logger::Error(e);
      exit(1);
   }
   catch(char const* e)
   {
      Logger::Error(e);
      exit(1);      
   }
   catch(BinaryReader::BRexception e)
   {
      Logger::Error(std::string(e.what())+"\n"+e.info);
      exit(1); 
   }
   catch(std::exception& e)
   {
      Logger::Error(e.what());
      exit(1);
   }

   //Write out new navigation file
   try
   {
      //try and open the output file
      BILWriter navOut(strOutputNavFile,BILWriter::float64,number_of_offset_scans,1,7,'w');
      //Create a buffer of size 1 line * 1 sample * 7 bands
      double* buffer=new double[7];
      Logger::Log("Writing out new navigation file...");
      for(unsigned int scan=0;scan<number_of_offset_scans;scan++)
      {
         //Fill the buffer with the interpolated data
         buffer[0]=offset_times[scan];
         buffer[1]=interpolated_nav.GetValue(scan,NavDataCollection::LAT);
         buffer[2]=interpolated_nav.GetValue(scan,NavDataCollection::LON);
         buffer[3]=interpolated_nav.GetValue(scan,NavDataCollection::HEI);
         buffer[4]=interpolated_nav.GetValue(scan,NavDataCollection::ROLL);
         buffer[5]=interpolated_nav.GetValue(scan,NavDataCollection::PITCH);
         buffer[6]=interpolated_nav.GetValue(scan,NavDataCollection::HEADING);
         //write out the buffer
         navOut.WriteLine((char*)buffer);
      }

      //Close the BIL writer
      navOut.Close();
      delete[] buffer;
   }
   catch(std::string e)
   {
      Logger::Error(e);
      exit(1);
   }
   catch(char const* e)
   {
      Logger::Error(e);
      exit(1);      
   }
   catch(BILWriter::BILexception e)
   {
      Logger::Error(std::string(e.what())+"\n"+e.info);
      exit(1); 
   }
   catch(std::exception& e)
   {
      Logger::Error(e.what());
      exit(1);
   }


   //If requested we can trim the level-1 image data here
   if(strLevel1File.compare("")!=0)
   {
      Logger::Log("Writing out new level 1 data file...");
      try
      {
         TrimLevel1Data(strLevel1File,strTrimmedLevel1File,level1_startscan,level1_endscan);
      }
      catch(std::string e)
      {
         Logger::Error(e);
         exit(1);
      }
      catch(char const* e)
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
      catch(std::exception e)
      {
         Logger::Error(e.what());
         exit(1);
      }
   }

   //Tidy up
   delete inNav;
   delete[] offset_times;

   return 0;
}
