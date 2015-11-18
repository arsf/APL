//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#include <fstream>
#include <typeinfo>
#include <sstream>
#include "TreeGrid.h"
#include "level3grid.h"
#include "binfile.h"
#include "commandline.h"
#include "map.h"
#include "os_dependant.h"
#include "filewriter.h"

//-------------------------------------------------------------------------
// Software description
//-------------------------------------------------------------------------
const std::string DESCRIPTION="Mapping / Gridding Software";

//----------------------------------------------------------------
//Number of options that can be on command line
//----------------------------------------------------------------
const int number_of_possible_options = 17;

//----------------------------------------------------------------
//Option names that can be on command line
//----------------------------------------------------------------
const std::string availableopts[number_of_possible_options]={
"-igm",
"-lev1",
"-mapname",
"-bandlist",
"-pixelsize",
"-area",
"-interpolation",
"-outputlevel",
"-buffersize",
"-maxinterpdistance",
"-outputdatatype",
"-ignorescan",
"-ignorevalue",
"-nodata",
"-rowcolmap",
"-ignorediskspace",
"-help"
}; 

//----------------------------------------------------------------
//One line description of each option - in same order as availableopts
//----------------------------------------------------------------
const std::string optsdescription[number_of_possible_options]={
"3-Band IGM BIL data file (Band 1: X, Band 2: Y, Band 3: Z) created by aplcorr / apltran.",
"Level 1 data file containing the data to map. This data must correspond to the IGM file created.",
"The filename of the mapped (output) BIL file.",
"The list of bands that you want to map from the level 1 file.",
"The size of the output mapped pixels in X and Y.",
"Defines a rectangular region to map a specific part of a flight line. ",
"Specify the interpolation algorithm to use. ",
"Set the level of the text output to the terminal. ",
"Set a buffersize in MB to use for storing level1 data while processing - note that the total required RAM will be larger than this. Default is 1024.",
"Set the maximum ground distance over which to interpolate. Defaults to three times the average nadir pixel spacing.",
"Select the output data format type. ",
"A space separated string of scans to ignore. These should be those identified as dropped scans.",
"A value to ignore in the level1 data and interpolate over. The default value is 0. Use 'NONE' to not ignore any value.",
"A value to set as the nodata value inserted into the mapped image. Default is 0.",
"Specify this, followed by an output filename, to output an additional BIL file that contains 2 bands: row and col values of the level-1 image in the mapped grid. This will only run with interpolation method 'nearest'",
"Process even if insufficient disk space is reported. Only use if the disk space reported is incorrect.",
"Display this help."
};

//--------------------------------------------------------------------------
//Function to return a more detailed help string about the available options
//--------------------------------------------------------------------------
std::string GetHelpFor(const std::string str)
{
   std::map <std::string,std::string> helpdoc;

   helpdoc["igm"]="\nThis is used to specify the IGM file to read the geolocation information from, as created by aplcorr or apltran. "
                  "The IGM file should contain 3 bands ordered such that band 1 = X, band 2 = Y and band 3 = Z. "
                  "The projection of the IGM data determines the projection of the mapped image data.\n";

   helpdoc["lev1"]="\nThis is used to specify the level-1 data file to map. This should correspond to the given IGM file.\n";

   helpdoc["maxinterpdistance"]="\nThis optional argument can be used to specify the maximum distance over which to interpolate. "
                  "If there are not enough points for the selected interpolation method which are within this distance from the centre of the grid cell, "
                  "a value of 0 will be entered.\n"
                  "\nThe default value for this is three times the average spacing of the central IGM pixels.\n";

   helpdoc["outputdatatype"]="\nThis is used to select the output format of the mapped data. Note however that data will be truncated and not scaled. "
                  "If a lower level is selected than the input level-1 data loss of information will occur.\nOptions:\n"
                  " uchar8 - for 8-bit unsigned data.\n"                  
                  " int16 - for 16-bit signed integer data.\n"
                  " uint16 - for 16-bit unsigned integer data.\n"
                  " int32 - for 32-bit signed integer data.\n"
                  " uint32 - for 32-bit unsigned integer data.\n"
                  " float32 - for 32-bit floating point data.\n"
                  " float64 - for 64-bit floating point data.\n";

   helpdoc["bandlist"]="\nThis is used to request which bands of the level-1 data will be mapped. This can be entered in 3 different formats:\n"
                  " ALL - using the keyword ALL will map every band of the level-1 data file.\n"
                  " x-y - using the '-' will map all bands between band x and band y, e.g. 5-15 will map bands 5 through 15 inclusive.\n"
                  " a b c - using a space separated list will map these individual bands, e.g. 1 7 17 107 99 will map only bands 1, 7, 17, 107 and 99.\n";

   helpdoc["interpolation"]="\nThis optional argument can be used to specify the algorithm to use for interpolating / resampling the gridded level-1 data.\n"
                  "\nOptions:\n"
                  " nearest - will use the level-1 pixel that is nearest to the centre of the mapped grid cell.\n"
                  " bilinear - will use a bilinear interpolation of the 4 nearest points - one from each quadrant surrounding the grid cell. \n"
                  " idw <max> - will use an inverse distance weighted algorithm with maximum number of points equal to <max>. \n"
                  " cubic - will use hermite cubic splines with the 16 nearest points, 4 in each quadrant surrounding grid cell.\n"
                  "\nThe default interpolation method that will be used is nearest. Further information can be found in the software documentation.\n";

   helpdoc["pixelsize"]="\nThis is used to specify the size of the output mapped pixels in X and Y. The units will depend on the projection of the input IGM file. "
                  "If the IGM file is in metres then the pixelsize should be given in metres. If the IGM file is in degrees then the pixel size should "
                  "be given in degrees also.\n"
                  "\nFor example, to map to square pixels of 2 metres, use: -pixelsize 2 2.\n"
                  "\nNote that -area has a option to specify the number of rows and columns of the output image which can be defined instead of a pixel size.";
                  

   helpdoc["outputlevel"]="\nThis optional argument is used to set the amount of detail that is output to the terminal during a run.\n"
                  "Options:\n"
                  " standard - Output the standard information to the terminal.\n"
                  " verbose - Outputs more information than standard - can be useful for tracking down where crashes occur.\n"
                  " debug - Can be very verbose and would not be recommended for general running.\n"
                  "\nThe default output level is standard.\n";

   helpdoc["area"]="\nThis optional argument can be used if only a portion of the flight line is required to be mapped."
                  "The area can be defined by a rectangle of coordinates in the same projection as the IGM file.\n"
                  "\nThe order the coordinates are given is: minimumX maximumX minimumY maximumY\n"
                  "Note that the full flight line is still read through by the software to account for repeat area coverage.\n"
                  "A 6 parameter option is also available where the 5th and 6th parameter describes the number of rows and columns that the final image "
                  "should have. This can be specified instead of the -pixelsize option if the dimensions of the output image are to be fixed rather than resolution.";

   helpdoc["ignorescan"]="\nThis allows a space separated list of scan lines to be entered that will be ignored in the mapper. This is useful if "
                         "there are suspicious scan lines in the level-1 data that you wish to be excluded from the mapped data.\n"
                         "For example -ignorescan 100 500 800  will ignore scan lines 100, 500 and 800 (referenced to 0) from the level-1 data.\n";

   helpdoc["ignorevalue"]="\nThis allows a value to be entered as the value that will be ignored by the mapper and interpolated over. For example for masked pixels."
                           "The value entered is converted to the data type of the level-1 file - e.g. -1 will be converted to 65535 in unsigned 16-bit integer data.\n"
                           "If you do not wish to ignore any data (i.e. do not fill in masked areas) then you can set this value to 'NONE'."
                           "\nThe default value for this is 0.\n";

   helpdoc["nodata"]="\nThis allows a value to be specified into the mapped data product to describe a cell with no valid data. "
                     "The default value used is 0.\n";

   helpdoc["rowcolmap"]="\nThis outputs an additional BIL file that contains the Level-1 data row/col that corresponds to each map grid pixel. "
                        "A negative value means that there is no image data for that map grid pixel. This only accounts for when data is not "
                        "interpolated over, e.g. it assumes that all data values (included masked ones) will be mapped. Therefore it gives the mapping "
                        "as if ignorevalue is set to NONE regardless of it's true value. Currently this method only runs with interpolation method 'nearest'."
                        "\nThis option must be followed by the name of the BIL file you wish to create to store the data in.";

   helpdoc["ignorediskspace"]="\nThis flag will continue processing even if there is insufficient space reported. "
                              "Required for some shared file systems where the amount of free space is not correctly reported. "
                              "Use with caution as if there isn't sufficient space for the output file, data will be written until the disk is filled.";


   //If the special keyword FULL is given then concatenate all the help strings and return
   if(str.compare("FULL")==0)
   {
      std::string ret="";
      for(std::map<std::string,std::string>::iterator it=helpdoc.begin();it!=helpdoc.end();it++)
         ret += "\n---------------------------\n-"+(*it).first+"\n---------------------------\n"+(*it).second+"\n";

      return ret;
   }
   else //just return the help string requested
   {
      if(helpdoc[str].compare("")==0)
         helpdoc[str]="No extra help for this topic yet.";

      return helpdoc[str];
   }
}

//----------------------------------------------------------------
// Program main loop starts here
//----------------------------------------------------------------
int main(int argc, char* argv[])
{
   std::cout.precision(10);
   int main_return_value=0; //return value of the program - 0 is good
   //A command line handling object
   CommandLine* cl=NULL;
   //Create a logger to output the data to terminal
   Logger log(0);
   int logoutputlevel=0;
   //Input IGM filename
   std::string strInputIGMFilename;
   //Level 1 filename
   std::string strLev1File;
   //Output mapped image name
   std::string strMapName;
   //List of bands to map
   std::string strBandList;
   unsigned int nbtomap=0;
   //Pixel size to output to
   double Xpixelsize=0,Ypixelsize=0;
   //User defined area to map - defaults to full IGM file
   Area* user_area=NULL;   
   //hold the choice of the interpolation method to use
   Interpolators::InterpolatorType interpolation_method;
   int numpoints=0;
   //Buffer size in MB to use for processing
   uint64_t process_buffer_sizeMB=0;
   //Maximum distance of which points are used in interpolation
   double maxinterpdist=0;

   //value to ignore in the level-1 data
   double ignore_lev1_value=0;

   //by default we will ignore data (that has value 0)
   bool IGNOREDATA=true;

   //value to set as no data
   double nodata_value=0;

   //should we round the top left corner of level 3 grid to be a multiple of pixel size
   //this will be true when using the IGM to define the treegrid (i.e. when using all data)
   //and false when using the -area flag to specify an output region
   bool TOROUND=false;

   //by default will exit if there is insufficient space for output files
   //can ignore this for systems where disk space reporting is incorrect.
   bool IGNORE_DISKSPACE=false;

   //Filename for rowcol mapping file if requested
   std::string strRowColMapFilename="";

   //vector to hold dropped scans in passed to IGMTreeGrid
   //these rows will then be ignored when loading into Tree
   std::vector<unsigned int> dropscanvector;

   //Get exe name without the path
   std::string niceexename=std::string(argv[0]);
   niceexename=niceexename.substr(niceexename.find_last_of("/\\")+1);

   //data type to output (e.g. float32)
   FileWriter::DataType output_data_type=FileWriter::float32;

   //Output a "nice" header
   Logger::FormattedInformation(niceexename,VERSION,DESCRIPTION);

   try
   {      
      //----------------------------------------------------------------------------
      // Create the command line object
      //----------------------------------------------------------------------------
      cl=new CommandLine(argv,argc);

      //----------------------------------------------------------------------------
      // Check if anything went wrong (most likely an exception will have 
      // been thrown but lets be safe)
      //----------------------------------------------------------------------------
      if(cl->IsGood()==false)
      {
         throw "An unknown error has occurred with the command line\n"; 
      }   

      //----------------------------------------------------------------------------
      // check the command line options given against the list of available 
      // options in the software
      //----------------------------------------------------------------------------
      std::string badoptions;
      int retval=cl->CheckAvailableOptions(availableopts,number_of_possible_options,&badoptions); 
      if(retval<0)
      {
         //Options are on command line that are not available in this software
         std::stringstream strout(std::stringstream::out);
         strout<<"There are "<<(-retval)<<" unrecognised options on command line: "<<badoptions;
         throw CommandLine::CommandLineException(strout.str());
      }

      Logger::Log("Command line used to run: "+cl->ReturnCLAsString());

      //----------------------------------------------------------------------------
      // Go through each possible option in turn and set up the required data / response
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-help"))
      {

         if(cl->NumArgsOfOpt("-help")==0)
         {
            log.Add(cl->ProgramUsage(number_of_possible_options,availableopts,optsdescription));
            log.Add("\nUse -help <argname> for further information about argument.\nUse -help FULL for a full help listing.");
            log.Flush();
         }
         else
         {
            if(cl->GetArg("-help").compare("FULL")==0)
               log.Add(cl->ProgramUsage(number_of_possible_options,availableopts,optsdescription));

            Logger::Log(GetHelpFor(cl->GetArg("-help")));
         }
         throw "";
      }

      //----------------------------------------------------------------------------
      // IGM file name
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-igm"))
      {
         //Check that an argument follows the igm option - and get it if it exists
         if(cl->GetArg("-igm").compare(optiononly)!=0)
         {
            strInputIGMFilename=cl->GetArg("-igm");
            //This is an existing file - use GetExistingFilePath to check file status
            strInputIGMFilename=GetExistingFilePath(strInputIGMFilename,true);
            log.Add("Will use input IGM BIL file: "+strInputIGMFilename);
         }
         else
            throw CommandLine::CommandLineException("Argument -igm must immediately precede the igm filename.\n");         
      }
      else
      {
         //Throw an exception
         throw CommandLine::CommandLineException("Argument -igm [the IGM BIL file to get geocorrection information from] must be present on the command line.\n");  
      }

      //----------------------------------------------------------------------------
      // Get the level-1 file name
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-lev1"))
      {
         //Check that an argument follows the lev1 option - and get it if it exists
         if(cl->GetArg("-lev1").compare(optiononly)!=0)
         {
            strLev1File=GetExistingFilePath(cl->GetArg("-lev1"),true);
            Logger::Log("Will use input Level-1 BIL file: "+strLev1File);
         }
         else
            throw CommandLine::CommandLineException("Argument -lev1 must immediately precede the level-1 filename.\n");         
      }
      else
      {
         //Throw an exception
         throw CommandLine::CommandLineException("Argument -lev1 [the Level-1 data BIL file to geocorrect] must be present on the command line.\n");  
      }

      //----------------------------------------------------------------------------
      // Get the pixel size to map too
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-pixelsize"))
      {
         //Check that two arguments follows the pixelsize option - and get them if they exist
         if(cl->NumArgsOfOpt("-pixelsize")!=2)
            throw CommandLine::CommandLineException("Argument -pixelsize must immediately precede the X and Y pixelsize.\n"); 
         else
         {
            //Use another try statement here to allow more information on where this error occured for users benefit.
            try
            {
               Xpixelsize=StringToDouble(cl->GetArg("-pixelsize",0));
               Ypixelsize=StringToDouble(cl->GetArg("-pixelsize",1));
            }
            catch(std::string e)
            {
               throw CommandLine::CommandLineException(std::string("Error with -pixelsize:\n")+e);
            }

            Logger::Log("Will use pixelsize: "+ToString(Xpixelsize)+" "+ToString(Ypixelsize));
         }
      }
      //If the pixel size is not on the command line the -area with 6 options is not present then exit
      else if(!((cl->GetArg("-area").compare(optiononly)!=0)&&(cl->NumArgsOfOpt("-area")==6)))
      {
         throw CommandLine::CommandLineException("Argument -pixelsize must be present on the command line.\n");
      }

      //----------------------------------------------------------------------------
      // Get the list of bands to map
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-bandlist"))
      {
         //Check that an argument follows the bandlist option - and get it if it exists
         if(cl->GetArg("-bandlist").compare(optiononly)!=0)
         {
            nbtomap=cl->NumArgsOfOpt("-bandlist");
            strBandList=cl->GetArg("-bandlist");
            Logger::Log("Will map bands: "+strBandList+" which sums to "+ToString(nbtomap)+" group(s) of bands.");
         }
         else
            throw CommandLine::CommandLineException("Argument -bandlist must immediately precede the list of bands to process.\n");         
      }
      else
      {
         throw CommandLine::CommandLineException("Argument -bandlist must be present on the command line.\n");  
      }

      //----------------------------------------------------------------------------
      // Get the filename of the output map
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-mapname"))
      {
         //Check that an argument follows the mapname option - and get it if it exists
         if(cl->GetArg("-mapname").compare(optiononly)!=0)
         {
            strMapName=CreatePath(cl->GetArg("-mapname"));
            Logger::Log("Will write map to: "+strMapName);
         }
         else
            throw CommandLine::CommandLineException("Argument -mapname must immediately precede the name of the output map file.\n");         
      }
      else
      {
         throw CommandLine::CommandLineException("Argument -mapname must be present on the command line.\n");  
      }

      //----------------------------------------------------------------------------
      // Get the area to map
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-area"))
      {
         //Check that 4 numbers follow the area option - and get them if they exist
         if((cl->GetArg("-area").compare(optiononly)!=0)&&(cl->NumArgsOfOpt("-area")==4))
         {
            double mix=StringToDouble(cl->GetArg("-area",0));
            double max=StringToDouble(cl->GetArg("-area",1));
            double miy=StringToDouble(cl->GetArg("-area",2));
            double may=StringToDouble(cl->GetArg("-area",3));
            user_area=new Area(mix,max,miy,may);

            Logger::Log("Will only map inside coordinates defined by area: Min X: "+ToString(mix)+" Max X: "+ToString(max)+" Min Y: "+ToString(miy)+" Max Y: "+ToString(may));
         }
         //Or check that 6 numbers follow the area option - and get them if they exist
         else if((cl->GetArg("-area").compare(optiononly)!=0)&&(cl->NumArgsOfOpt("-area")==6))
         {            
            double mix=StringToDouble(cl->GetArg("-area",0));
            double max=StringToDouble(cl->GetArg("-area",1));
            double miy=StringToDouble(cl->GetArg("-area",2));
            double may=StringToDouble(cl->GetArg("-area",3));
            user_area=new Area(mix,max,miy,may);
            //User has specified number of rows and columns they want in the output image
            //Use these to calculate the pixel size - throw an exception if pixel size has been specified also 
            unsigned int nrows=StringToUINT(cl->GetArg("-area",4));
            unsigned int ncols=StringToUINT(cl->GetArg("-area",5));
            if(cl->OnCommandLine("-pixelsize"))
            {
               throw "Argument -area cannot be specified with 6 parameters if the -pixelsize has been specified also.";
            }
            Xpixelsize=(max-mix)/ncols;
            Ypixelsize=(may-miy)/nrows;
            Logger::Log("Derived pixel sizes (X,Y): "+ToString(Xpixelsize)+" "+ToString(Ypixelsize));
            //Check for possible floating point errors - I'm not sure if this is required or if it prevents problems further down the line.
            //We round the result, as if there is a floating point error (e) it is assumed to be small, giving a result of
            //(ncols - e) or (ncols +e)  - HOWEVER, maybe rounding actually makes this test worthless as it hides the error we're looking for?
            //Hence I have changed the test to ceil as this will show when an error is present, and is the method used in the grid dimension calculator.
            if((ceil((max-mix)/Xpixelsize) !=ncols) || (ceil((may-miy)/Ypixelsize) !=nrows))
            {
               throw "Rounding error is preventing the correct number of rows/columns to be generated from the derived pixel sizes. "
                     "Derived pixel size X: "+ToString(Xpixelsize)+"  pixel size Y: "+ToString(Ypixelsize)+"\nIt may be possible to define pixel size "
                     "instead to get desired dimensions, report this as an error.";
            }            
         }
         else
            throw CommandLine::CommandLineException("Argument -area must immediately precede the bounds to define the area rectangle.\n");         
      }
      else
      {
         //Default to full flight line - i.e. do nothing here
      }

      //----------------------------------------------------------------------------
      // Specify the interpolation algorithm to use
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-interpolation"))
      {
         //Check that 2 args follow the -interpolation option - and get them if they exist
         if((cl->GetArg("-interpolation").compare(optiononly)!=0))
         {
            std::string strInterpolate=cl->GetArg("-interpolation",0);
            if((strInterpolate.compare("idw")==0)&&(cl->NumArgsOfOpt("-interpolation")!=2))
            {
               throw CommandLine::CommandLineException("Argument -interpolation must immediately precede the algorithm name, and also the number of points to use if algorithm is idw.\n");         
            }
            //get the number of points to use if idw
            if(strInterpolate.compare("idw")==0)
               numpoints=StringToINT(cl->GetArg("-interpolation",1));
            else if(strInterpolate.compare("bilinearlev1")==0)
               numpoints=10; // search nearest 10 points to give us a selection of seeds
            else
               numpoints=1;

            Logger::Log("Will map using the interpolation algorithm defined by: "+strInterpolate);//+" and using "+ToString(numpoints)+" points.");
            if(strInterpolate.compare("nearest")==0)
               interpolation_method=Interpolators::NEARESTNEIGHBOUR;
            else if(strInterpolate.compare("idw")==0)
               interpolation_method=Interpolators::IDW;
            else if(strInterpolate.compare("bilinear")==0)
               interpolation_method=Interpolators::BILINEARLEVEL3;
            else if(strInterpolate.compare("bilinearlev1")==0)
               interpolation_method=Interpolators::BILINEAR;
            else if(strInterpolate.compare("cubic")==0)
               interpolation_method=Interpolators::CUBIC;
            else
               throw CommandLine::CommandLineException("Unknown interpolation algorithm requested.\n");  
         }     
         else
            throw CommandLine::CommandLineException("Argument -interpolation must immediately precede the algorithm name.\n");         
      }
      else
      {
         //Default to nearest neighbour
         interpolation_method=Interpolators::NEARESTNEIGHBOUR;
         Logger::Log("Using default interpolator of: nearest neighbour");
      }

      //----------------------------------------------------------------------------
      // Specify the maximum interpolation distance
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-maxinterpdistance"))
      {
         //Check that an arg follows the -maxinterpdistance option 
         if((cl->GetArg("-maxinterpdistance").compare(optiononly)!=0)&&(cl->NumArgsOfOpt("-maxinterpdistance")==1))
         {
            try
            {
               maxinterpdist=StringToDouble(cl->GetArg("-maxinterpdistance",0));
            }
            catch(std::string e)
            {
               throw CommandLine::CommandLineException("Error with -maxinterpdistance:\n"+e);         
            }
            Logger::Log("Will use a maximum interpolator distance of: "+ToString(maxinterpdist));
         }     
         else
            throw CommandLine::CommandLineException("Argument -maxinterpdistance must immediately precede the maximum distance value, with no extra parameters.\n");         
      }
      else
      {
         //Will use a default value
         maxinterpdist=-1;
//         Logger::Log("Using default maximum interpolator distance of: "+ToString(maxinterpdist));
      }
      

      //----------------------------------------------------------------------------
      // Specify the output level to use
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-outputlevel"))
      {
         //Check that a code word follows the option
         if((cl->GetArg("-outputlevel").compare(optiononly)!=0)&&(cl->NumArgsOfOpt("-outputlevel")==1))
         {
            std::string outputlevel=cl->GetArg("-outputlevel",0);
            if(outputlevel.compare("standard")==0)
               logoutputlevel=0;
            else if(outputlevel.compare("verbose")==0)
               logoutputlevel=1;
            else if(outputlevel.compare("debug")==0)
               logoutputlevel=2;
            else
               throw CommandLine::CommandLineException("Unknown output level.\n");  
         }     
         else
            throw CommandLine::CommandLineException("Argument -outputlevel must immediately precede the level to use: standard, verbose or debug. Nothing else should follow.\n");         
      }
      else
      {
         //Default to standard
         logoutputlevel=0;
      }

      //----------------------------------------------------------------------------
      // Specify the processing buffer size
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-buffersize"))
      {
         //Check that a value follows the option
         if((cl->GetArg("-buffersize").compare(optiononly)!=0)&&(cl->NumArgsOfOpt("-buffersize")==1))
         {
            process_buffer_sizeMB=StringToUINT(cl->GetArg("-buffersize",0));
            Logger::Log("Will assign a buffer size for input data of (MB): "+ToString(process_buffer_sizeMB));
         }     
         else
            throw CommandLine::CommandLineException("Argument -buffersize must immediately precede the value of the buffer size in MB, with no extra parameters.\n");         
      }
      else
      {
         //Default to 1GB
         process_buffer_sizeMB=1024;
      }


      //----------------------------------------------------------------------------
      // Specify the output data type
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-outputdatatype"))
      {
         //Check that a value follows the option
         if((cl->GetArg("-outputdatatype").compare(optiononly)!=0)&&(cl->NumArgsOfOpt("-outputdatatype")==1))
         {
            std::string strOutputdatatype=(cl->GetArg("-outputdatatype",0));
            //Is there a better way to do this?
            if(strOutputdatatype.compare("char")==0)
            {
               output_data_type=FileWriter::uchar8;
            }            
            else if(strOutputdatatype.compare("int16")==0)
            {
               output_data_type=FileWriter::int16;
            }
            else if(strOutputdatatype.compare("uint16")==0)
            {
               output_data_type=FileWriter::uint16;
            }
            else if(strOutputdatatype.compare("int32")==0)
            {
               output_data_type=FileWriter::int32;
            }
            else if(strOutputdatatype.compare("uint32")==0)
            {
               output_data_type=FileWriter::uint32;
            }
            else if(strOutputdatatype.compare("float32")==0)
            {
               output_data_type=FileWriter::float32;
            }
            else if(strOutputdatatype.compare("float64")==0)
            {
               output_data_type=FileWriter::float64;
            }
            else
            {
               throw "Unrecognised data type given in -outputdatatype. Refer to -help outputdatatype for accepted keywords.";
            }
            Logger::Log("Will write out data as type: "+strOutputdatatype);
         }     
         else
            throw CommandLine::CommandLineException("Argument -outputdatatype must immediately precede the data type of the output data, with no extra parameters.\n");         
      }
      else
      {
         //Default to float32
         output_data_type=FileWriter::float32;
         Logger::Log("Will write out data as default type: float32");
      }

      //----------------------------------------------------------------------------
      // Get the list of dropped scans to ignore
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-ignorescan"))
      {
         //Check that an argument follows the ignorescan option - and get it if it exists
         if(cl->GetArg("-ignorescan").compare(optiononly)!=0)
         {
            std::string strtmp=cl->GetArg("-ignorescan");
            int value=0,position=0;
            Logger::Log("Will ignore level 1 scans (from command line): "+strtmp);

            //Get number of lines of level1 file to check against
            BinFile tmp(strLev1File); 
            int l1lines=StringToINT(tmp.FromHeader("lines"));
            tmp.Close();

            //Loop through the list to get the scans that will be ignored
            while(GetItemFromString(strtmp,position,';').compare("")!=0)
            {
               std::string ascan=GetItemFromString(strtmp,position,';');
               value=StringToINT(CheckNumbersOnly(ascan));
               if((value < 0) || (value > l1lines))
                  throw "Scanline is outside range of number of scans in file: "+ascan;

               dropscanvector.push_back(value);
               position++;
            }
            strtmp="";
            for(unsigned int i=0;i<dropscanvector.size();i++)
               strtmp += ToString(dropscanvector[i])+" ";
            Logger::Log("Will ignore level 1 scans (after parsing string): "+strtmp);

         }
         else
            throw CommandLine::CommandLineException("Argument -ignorescan must immediately precede the list of scan lines to ignore - these should be dropped scans.\n");         
      }

      //----------------------------------------------------------------------------
      // Specify the data ignore value in the level1 file
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-ignorevalue"))
      {
         //Check that a value follows the option
         if((cl->GetArg("-ignorevalue").compare(optiononly)!=0)&&(cl->NumArgsOfOpt("-ignorevalue")==1))
         {
            std::string tmpstr=cl->GetArg("-ignorevalue",0);
            if(tmpstr.compare("NONE")==0)
            {
               //Do not ignore 0's - map what is there
               IGNOREDATA=false;
            }
            else
            {
               ignore_lev1_value=StringToDouble(tmpstr);
               Logger::Log("Will ignore data values (in the level1 file) of: "+ToString(ignore_lev1_value));
            }
         }     
         else
            throw CommandLine::CommandLineException("Argument -ignorevalue must immediately precede the value to be ignored in the level-1 data.\n");         
      }
      else
      {
         //Default to 0
         ignore_lev1_value=0;
      }

      //----------------------------------------------------------------------------
      // Specify the no data value in the mapped file
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-nodata"))
      {
         //Check that a value follows the option
         if((cl->GetArg("-nodata").compare(optiononly)!=0)&&(cl->NumArgsOfOpt("-nodata")==1))
         {
            std::string tmpstr=cl->GetArg("-nodata",0);
            nodata_value=StringToDouble(tmpstr);
            Logger::Log("Will set nodata value in mapped image to be: "+ToString(nodata_value));
         }     
         else
            throw CommandLine::CommandLineException("Argument -nodata must immediately precede the value to be set as nodata.\n");         
      }
      else
      {
         //Default to 0
         nodata_value=0;
      }

      //----------------------------------------------------------------------------
      // Specify output of a row/col mapping file
      // This shows which level-1 pixel was used to fill each map grid pixel
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-rowcolmap"))
      {
         //Check that a value follows the option
         if((cl->GetArg("-rowcolmap").compare(optiononly)!=0)&&(cl->NumArgsOfOpt("-rowcolmap")==1))
         {
            strRowColMapFilename=CreatePath(cl->GetArg("-rowcolmap"));
            Logger::Log("Will output Level-1 row/col mapping to file: "+strRowColMapFilename);
         }     
         else
            throw CommandLine::CommandLineException("Argument -rowcolmap must immediately precede the name of the file to write row/col data to.\n");         
      
      }
      else
      {
         //Nothing - we won't write row/col data to a file
      }
      //-------------------------------------------------------------------
      // Ignore if insufficient disk space is reported
      //-------------------------------------------------------------------
      if(cl->OnCommandLine("-ignorediskspace"))
      {
         if(cl->GetArg("-ignorediskspace").compare(optiononly)!=0)
         {
            throw CommandLine::CommandLineException("Option -ignorediskspace does not take any arguments.\n");
         }
         else
         {
            IGNORE_DISKSPACE=true;
         }
      }
      else
      {
         IGNORE_DISKSPACE=false;
      }

 
   }
   catch(CommandLine::CommandLineException e)
   {
      Logger::Error(std::string(e.what())+"\n"+e.info);
      delete cl;
      exit(1);
   }
   catch(BinaryReader::BRexception e)
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
   catch(...)
   {
      PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString());
      delete cl;
      exit(1);
   } 

   //----------------------------------------------------------------------
   // END OF COMMAND LINE HANDLING - START OF REAL CODE
   //----------------------------------------------------------------------

   //Set the output level of the logger to the requested value
   log.SetLevel(logoutputlevel);


   //----------------------------------------------------------------------
   // Check the requested bands to map exist and get the wavelengths
   //----------------------------------------------------------------------
   try
   {
      if(strBandList.compare("ALL")==0)
      {
         strBandList.clear(); //remove the all keyword
         //Need to make a list of all bands
         BinFile tmp(strLev1File); 
         for(unsigned int b=0;b<StringToUINT(tmp.FromHeader("bands"));b++)
         {
            strBandList=strBandList+ToString(b)+" ";
         }
         tmp.Close();
      }
      else
      {
         //Need to make the band list space separated and zero based
         //assuming the bands from command line start referenced from 1
         strBandList=ReplaceAllWith(&strBandList,';',' ');   
         std::string tmpBandList; //tmp copy of list

         //Go through each item of string and check if a single band or group
         int position=0;
         while(GetItemFromString(strBandList,position,' ').compare("")!=0)
         {
            std::string aband=GetItemFromString(strBandList,position,' ');
            if(aband.find("-")==std::string::npos)
            {
               //Single band - add it on to tmp string
               tmpBandList=tmpBandList+aband+" ";
            }
            else
            {
               //group of bands - add them all onto tmp string individually
               std::string firstband=GetItemFromString(aband,0,'-');
               std::string lastband=GetItemFromString(aband,1,'-');
               //Test the values to see if they make sense
               if((firstband.compare("")==0)||(lastband.compare("")==0))
               {
                  Logger::Log("There was a problem getting the bands to map from group a-b form. I got: "+firstband+"-"+lastband);
                  exit(1);
               }

               if((StringToUINT(firstband)==0)||(StringToUINT(lastband)==0))
               {
                  Logger::Log("At least one of the requested bands to map does not exist: "+firstband+" "+lastband);
                  exit(1);               
               }

               if(StringToUINT(firstband)>=StringToUINT(lastband))
               {
                  Logger::Error("In band groupings with a-b form, a must be less than b.");
                  exit(1);
               }
               
               for(unsigned int b=StringToUINT(firstband);b<=StringToUINT(lastband);b++)
               {
                  tmpBandList=tmpBandList+ToString(b)+" ";
               }
            }
            position++;
         }

         strBandList.clear();//clear original band list
         //update the number of bands to map
         nbtomap=GetNumberOfItemsFromString(TrimWhitespace(tmpBandList)," ");

         unsigned int ba=0;
         BinFile tmp(strLev1File); 
         for(unsigned int b=0;b<nbtomap;b++)
         {
            //Get the band and subtract one from it 
            ba=StringToUINT(GetItemFromString(tmpBandList,b))-1;
            if(ba >= StringToUINT(tmp.FromHeader("bands")))
            {
               Logger::Log("Requested band to map does not exist: "+GetItemFromString(tmpBandList,b));
               exit(1);
            }

            //convert back to the band list string
            strBandList=strBandList+ToString(ba)+" ";
         }
         tmp.Close();
      }
      //Remove any whitespace at the end of the string
      strBandList=TrimWhitespace(strBandList);

      Logger::Debug("Bands to map (zero based): "+strBandList);

   }
   catch(BinaryReader::BRexception e)
   {
      Logger::Error(std::string(e.what())+"\n"+e.info);
      delete cl;
      exit(1); 
   }
   catch(...)
   {
      PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString());
      delete cl;
      exit(1);
   } 

   //----------------------------------------------------------------------
   // Set up the TreeGrid
   //----------------------------------------------------------------------
   IGMTreeGrid* tg=NULL;
   Area* fulltree=NULL;
   try
   {
      tg=new IGMTreeGrid(strInputIGMFilename,dropscanvector,user_area);
      //Create an area based on the full IGM gridtree
      //fulltree = new Area(tg->TopLeftX(),tg->TopLeftX()+(tg->SizeX()*tg->NumCols()),tg->TopLeftY()-(tg->SizeY()*tg->NumRows()),tg->TopLeftY());
      fulltree = new Area(tg->TopLeftX(),tg->BottomRightX(),tg->BottomRightY(),tg->TopLeftY());
      Logger::Log("\nArea defined by IGM file fits within rectangle: "+ToString(fulltree->MinX())+" < X < "+ToString(fulltree->MaxX())+" "+
                  ToString(fulltree->MinY())+" < Y < "+ToString(fulltree->MaxY()));

      //If user area is undefined - define it as the full tree area to pass to the mapper
      if(user_area==NULL)
      {
         user_area=fulltree;
         TOROUND=true;
      }
   }
   catch(BinaryReader::BRexception e)
   {
      Logger::Error(std::string(e.what())+"\n"+e.info);
      delete cl;
      delete tg;
      exit(1);
   }
   catch(std::string e)
   {
      Logger::Error(e);
      delete cl;
      delete tg;
      exit(1);
   }
   catch(char const* e)
   {
      Logger::Error(e);
      delete cl;
      delete tg;
      exit(1);      
   }
   catch(std::exception &e)
   {
      PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString(),&e);
      delete cl;
      delete tg;
      exit(1);       
   } 
   catch(...)
   {
      PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString());
      delete cl;
      delete tg;
      exit(1);
   } 

   //----------------------------------------------------------------------
   // Set up the level 3 (mapped) grid object and data writer
   // This is done using the map class
   //----------------------------------------------------------------------
   //We need to know the data type of the level1 data to map
   BinFile tmp(strLev1File); 
   unsigned int intype=tmp.GetDataType();
   tmp.Close();

   //Use an abstract base map to construct the map we want
   AbstractMap* map=NULL;
   try
   {
      switch(intype)
      {     
      case 1:// byte data
         map=new Map<char>(strMapName,Xpixelsize,Ypixelsize,strBandList,user_area,strLev1File,interpolation_method,numpoints,process_buffer_sizeMB*1024*1024,output_data_type,strRowColMapFilename,nodata_value,TOROUND);
         break;
      case 2://signed 16bit integer data
         map=new Map<int16_t>(strMapName,Xpixelsize,Ypixelsize,strBandList,user_area,strLev1File,interpolation_method,numpoints,process_buffer_sizeMB*1024*1024,output_data_type,strRowColMapFilename,nodata_value,TOROUND);
         break;
      case 3://signed 32bit integer byte data
         map=new Map<int32_t>(strMapName,Xpixelsize,Ypixelsize,strBandList,user_area,strLev1File,interpolation_method,numpoints,process_buffer_sizeMB*1024*1024,output_data_type,strRowColMapFilename,nodata_value,TOROUND);
         break;
      case 4://float32 data
         map=new Map<float>(strMapName,Xpixelsize,Ypixelsize,strBandList,user_area,strLev1File,interpolation_method,numpoints,process_buffer_sizeMB*1024*1024,output_data_type,strRowColMapFilename,nodata_value,TOROUND);
         break;
      case 5:// double data
         map=new Map<double>(strMapName,Xpixelsize,Ypixelsize,strBandList,user_area,strLev1File,interpolation_method,numpoints,process_buffer_sizeMB*1024*1024,output_data_type,strRowColMapFilename,nodata_value,TOROUND);
         break;
      case 12://16-bit unsigned short int data
         map=new Map<uint16_t>(strMapName,Xpixelsize,Ypixelsize,strBandList,user_area,strLev1File,interpolation_method,numpoints,process_buffer_sizeMB*1024*1024,output_data_type,strRowColMapFilename,nodata_value,TOROUND);
         break; 
      case 13://32-bit unsigned int data
         map=new Map<uint32_t>(strMapName,Xpixelsize,Ypixelsize,strBandList,user_area,strLev1File,interpolation_method,numpoints,process_buffer_sizeMB*1024*1024,output_data_type,strRowColMapFilename,nodata_value,TOROUND);
         break; 
      default:
         throw "Unrecognised data type in level 1 file. Cannot create a map of this data type. Got: "+ToString(intype);
         break;
      }
      map->AssignProjection(tg->GetMapInfo());
   }   
   catch(BinaryReader::BRexception e)
   {
      Logger::Error(std::string(e.what())+"\n"+e.info);
      delete cl;
      delete tg;
      exit(1);
   }
   catch(FileWriter::FileWriterException e)
   {
      Logger::Error(std::string(e.what())+"\n"+e.info);
      delete cl;
      delete tg;
      exit(1);      
   }
   catch(std::string e)
   {
      Logger::Error(e);
      delete cl;
      delete tg;
      exit(1);
   }
   catch(char const* e)
   {
      Logger::Error(e);
      delete cl;
      delete tg;
      exit(1);      
   }
   catch(std::exception &e)
   {
      PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString(),&e);
      delete cl;
      delete tg;
      exit(1);       
   }   
   catch(...)
   {
      PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString());
      delete cl;
      delete tg;
      exit(1);
   } 

   //----------------------------------------------------------------------
   // Output some information at this point
   //----------------------------------------------------------------------

   Logger::Log("Initialisation is now complete. Starting mapping ...");
   Logger::Log("Number of rows, columns and bands of output mapped image file: "+ToString(map->grid->NumRows())+" "+ToString(map->grid->NumCols())+" "+ToString(map->grid->NumBands()));
   Logger::Log("Final mapped image will require "+ToString(ceil((uint64_t)map->grid->NumRows()*map->grid->NumCols()*map->grid->NumBands()*map->GetOutputDataSize()/1024.0/1024.0))
                  +" megabytes of disk space.");

   //----------------------------------------------------------------------
   //Test disk space
   //----------------------------------------------------------------------
   DiskSpace ds;
   uint64_t avds=ds.GetAvailableSpace(strMapName);
   uint64_t tds=ds.GetTotalSpace(strMapName);
   //Convert to MB and print to terminal
   Logger::Log("Total amount of disk space: "+ToString(tds/(1024*1024)));   
   Logger::Log("Amount of free disk space available: "+ToString(avds/(1024*1024)));

   if(avds < (uint64_t)map->grid->NumRows()*map->grid->NumCols()*map->grid->NumBands()*map->GetOutputDataSize())
   {
      //There is not enough disk space
      if(IGNORE_DISKSPACE)
      {
         Logger::Warning("Insufficient disk space was found for processing. Ignoring and carrying on as '--ignorediskspace' option used.");
      }
      else
      {
         Logger::Error("There is not enough available disk space for this data file.");
         exit(1);
      }
   }

   //----------------------------------------------------------------------
   // Start the mapping loop
   //----------------------------------------------------------------------

   std::vector<Item>::iterator iter;

   try
   {
      //Get nadir pixel separation (in metres)
      double xsep=0,ysep=0;
      tg->GetAveragePixelSeparationMetres(xsep,ysep);
      //Get a default value for the maxinterpdist if required (if maxinterpdist is currently -ve)
      if(maxinterpdist <=0)
      {
         //Use a value of three times nadir pixel separation
         maxinterpdist=3*sqrt(xsep*xsep+ysep*ysep);
         //maxinterpdist=0.6*sqrt(tg->SizeX()*tg->SizeX() + tg->SizeY()*tg->SizeY());
         Logger::Log("Will use a default value for maximum interpolation of three times the average separation of a nadir pixel. This is: "+ToString(maxinterpdist)+" metres.");
      }

      //Send a warning if the maximum interpolation distance is rather large - only do for no geographic tree grids as maxinterpdist is always in metres now
      if(maxinterpdist > 50*sqrt(xsep*xsep + ysep*ysep))
      {
         Logger::Warning("Maximum interpolation distance is large compared to pixel spacing [approx nadir pixel spacing: "+ToString(sqrt(xsep*xsep + ysep*ysep))+"], maybe this is a mistake?");
      }

      //Warn if pixel size looks vastly incorrect - use a (randomly picked) value of 5*TreeGrid cell size
      if((Xpixelsize > 5*tg->SizeX()) || (Ypixelsize > 5*tg->SizeY()))
      {
         Logger::Warning("Given pixel size appears very large compared to data resolution of IGM file. Are you sure the pixel size is in the correct units for the IGM file?");
      }

      //Have to go back to specific Map type rather than using base AbstractMap type
      //Would be better if we didn't have to ...
      Map<char>* cmap=NULL;
      Map<int16_t>* smap=NULL;
      Map<uint16_t>* usmap=NULL;
      Map<int32_t>* imap=NULL;
      Map<uint32_t>* uimap=NULL;
      Map<float>* fmap=NULL;
      Map<double>* dmap=NULL;

      switch(intype)
      {
      case 1://8-bit data
         //Now cast back to the actual data-type Map object
         cmap=dynamic_cast<Map<char>*>(map);
         //Apply the maximum interpolation distance to the map interpolator
         if(maxinterpdist > 0)
         {
            cmap->SetMaxInterpolationDistance(maxinterpdist);
            cmap->SetInterpolatorIgnoreFlag(IGNOREDATA);
            if(IGNOREDATA==true)
               cmap->SetInterpolatorIgnoreValue(ignore_lev1_value);
         }
         else
            throw "Maximum interpolation distance is not greater than 0 - select a suitable size on command line.";

         //Call the mapping function
         Logger::Verbose("Calling MapLineSegments to create the mapped image.");
         cmap->MapLineSegments(tg,strInputIGMFilename,strLev1File);
         break;

      case 2://16-bit data
         smap=dynamic_cast<Map<int16_t>*>(map);
         //Apply the maximum interpolation distance to the map interpolator
         if(maxinterpdist > 0)
         {
            smap->SetMaxInterpolationDistance(maxinterpdist);
            smap->SetInterpolatorIgnoreFlag(IGNOREDATA);
            if(IGNOREDATA==true)
               smap->SetInterpolatorIgnoreValue(ignore_lev1_value);
         }
         else
            throw "Maximum interpolation distance is not greater than 0 - select a suitable size on command line.";

         //Call the mapping function
         Logger::Verbose("Calling MapLineSegments to create the mapped image.");
         smap->MapLineSegments(tg,strInputIGMFilename,strLev1File);
         break;

      case 3://32-bit integer data
         imap=dynamic_cast<Map<int32_t>*>(map);
         //Apply the maximum interpolation distance to the map interpolator
         if(maxinterpdist > 0)
         {
            imap->SetMaxInterpolationDistance(maxinterpdist);
            imap->SetInterpolatorIgnoreFlag(IGNOREDATA);
            if(IGNOREDATA==true)
               imap->SetInterpolatorIgnoreValue(ignore_lev1_value);
         }
         else
            throw "Maximum interpolation distance is not greater than 0 - select a suitable size on command line.";

         //Call the mapping function
         Logger::Verbose("Calling MapLineSegments to create the mapped image.");
         imap->MapLineSegments(tg,strInputIGMFilename,strLev1File);
         break;
      
      case 4://float32 data
         //Now cast back to the actual data-type Map object
         fmap=dynamic_cast<Map<float>*>(map);
         //Apply the maximum interpolation distance to the map interpolator
         if(maxinterpdist > 0)
         {
            fmap->SetMaxInterpolationDistance(maxinterpdist);
            fmap->SetInterpolatorIgnoreFlag(IGNOREDATA);
            if(IGNOREDATA==true)
               fmap->SetInterpolatorIgnoreValue(ignore_lev1_value);
         }
         else
            throw "Maximum interpolation distance is not greater than 0 - select a suitable size on command line.";

         //Call the mapping function
         Logger::Verbose("Calling MapLineSegments to create the mapped image.");
         fmap->MapLineSegments(tg,strInputIGMFilename,strLev1File);
         break;

      case 5://float64 data
         //Now cast back to the actual data-type Map object
         dmap=dynamic_cast<Map<double>*>(map);
         //Apply the maximum interpolation distance to the map interpolator
         if(maxinterpdist > 0)
         {
            dmap->SetMaxInterpolationDistance(maxinterpdist);
            dmap->SetInterpolatorIgnoreFlag(IGNOREDATA);
            if(IGNOREDATA==true)
               dmap->SetInterpolatorIgnoreValue(ignore_lev1_value);
         }
         else
            throw "Maximum interpolation distance is not greater than 0 - select a suitable size on command line.";

         //Call the mapping function
         Logger::Verbose("Calling MapLineSegments to create the mapped image.");
         dmap->MapLineSegments(tg,strInputIGMFilename,strLev1File);
         break;

      case 12://16-bit unsigned short int data
         //Now cast back to the actual data-type Map object
         usmap=dynamic_cast<Map<uint16_t>*>(map);
         //Apply the maximum interpolation distance to the map interpolator
         if(maxinterpdist > 0)
         {
            usmap->SetMaxInterpolationDistance(maxinterpdist);
            usmap->SetInterpolatorIgnoreFlag(IGNOREDATA);
            if(IGNOREDATA==true)
               usmap->SetInterpolatorIgnoreValue(ignore_lev1_value);
         }
         else
            throw "Maximum interpolation distance is not greater than 0 - select a suitable size on command line.";

         //Call the mapping function
         Logger::Verbose("Calling MapLineSegments to create the mapped image.");
         usmap->MapLineSegments(tg,strInputIGMFilename,strLev1File);
         break; 

      case 13://32-bit unsigned short int data
         //Now cast back to the actual data-type Map object
         uimap=dynamic_cast<Map<uint32_t>*>(map);
         //Apply the maximum interpolation distance to the map interpolator
         if(maxinterpdist > 0)
         {
            uimap->SetMaxInterpolationDistance(maxinterpdist);
            uimap->SetInterpolatorIgnoreFlag(IGNOREDATA);
            if(IGNOREDATA==true)
               uimap->SetInterpolatorIgnoreValue(ignore_lev1_value);
         }
         else
            throw "Maximum interpolation distance is not greater than 0 - select a suitable size on command line.";

         //Call the mapping function
         Logger::Verbose("Calling MapLineSegments to create the mapped image.");
         uimap->MapLineSegments(tg,strInputIGMFilename,strLev1File);
         break; 

      default:
         throw "Unrecognised data type in level 1 file. Cannot create a map of this data type. Got: "+ToString(intype);
         break;
      }

   }
   catch(BinaryReader::BRexception e)
   {
      Logger::Error(std::string(e.what())+"\n"+e.info);
      main_return_value=1;
   }
   catch(std::string e)
   {
      Logger::Error(e);
      main_return_value=1;
   }
   catch(char const* e)
   {
      Logger::Error(e);
      main_return_value=1;
   }
   catch(std::exception &e)
   {
      PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString(),&e);
      main_return_value=1;   
   }
   catch(...)
   {
      PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString());
      main_return_value=1;   
   }

   if(map!=NULL)
      delete map;

   //Clear up the user area
   if(user_area!=NULL)
      delete user_area;
   //clear up the tree grid
   delete tg;

   return main_return_value;
}
