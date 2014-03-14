//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

//This is the main code file for the geolocation software

#include "commandline.h"
#include "viewvectors.h"
#include "navbaseclass.h"
#include "conversions.h"
#include "bilwriter.h"
#include "binfile.h"
#include "cartesianvector.h"
#include "dems.h"
#include "logger.h"
#include "planarsurface.h"
#include "leverbore.h"
#include "transformations.h"
#include "geodesics.h"

#include <string>
#include <cerrno>
#include <cstdlib>
#include <blitz/array.h>
#include <algorithm>


enum conversion_type {ECEF};
enum vvmethods {NONE,SPLIT, COMBINED};

//-------------------------------------------------------------------------
// Value to assign to pixels that cannot be located properly. For example
// their view vector may be above the horizon
//-------------------------------------------------------------------------
const int BADDATAVALUE=-9999;

//-------------------------------------------------------------------------
// Maximum allowed viewvector angle (degrees) - all others above this will be given 
// the bad data value and not geocorrected
//-------------------------------------------------------------------------
const float defaultmaxallowedvvangle=80.0;

//-------------------------------------------------------------------------
// Software description
//-------------------------------------------------------------------------
const std::string DESCRIPTION="Geocorrection Software";

//----------------------------------------------------------------
//Number of options that can be on command line
//----------------------------------------------------------------
const int number_of_possible_options = 13;

//----------------------------------------------------------------
//Option names that can be on command line
//----------------------------------------------------------------
const std::string availableopts[number_of_possible_options]={
"-vvfile",
"-navfile",
"-boresight",
"-igmfile",
"-ellipsoid",
"-vvSPLIT",
"-vvCOMBINE",
"-heightoffset",
"-dem",
"-lev1file",
"-atmosfile",
"-maxvvangle",
"-help"
}; 

//----------------------------------------------------------------
//One line description of each option - in same order as availableopts
//----------------------------------------------------------------
const std::string optsdescription[number_of_possible_options]={
"Instrument View Vector file",
"Navigation BIL file",
"Boresight corrections to apply: Pitch Roll Heading (defaults 0 0 0)",
"Output pixel position IGM file",
"Ellipsoid model to use: WGS84 (default is WGS84)",
"Use the \"split\" view vector method (this is the default method)",
"Use the \"combined\" view vector method",
"Add a height offset to the ellipsoid surface when mapping to the ellipsoid.",
"Digital Elevation Model to use for geocorrection. A 1 band BSQ/BIL file with heights in WGS84 Latitude/Longitude referenced to the WGS84 ellipsoid.)",
"Level 1 data filename - uses this to bin and trim view vector file to fit level 1 data set",
"Filename to output extra parameters to which are useful for atmospheric correction. These are: view azimuth and zenith, dem slope and dem aspect at intersect dem cell.",
"Maximum allowed view vector look angle in degrees. Sometimes if mapping on a tight bank of the aircraft view vectors can reach above the horizon. To prevent this cap the viewvectors to this maximum value. Default is "+ToString(defaultmaxallowedvvangle),
"Display this help"
}; 

//----------------------------------------------------------------
// Function definitions
//----------------------------------------------------------------

//Delete Arrays or Objects if they have been created, to free up memory
void TidyArrays(double* Plat,double* Plon,double* Pheight,double* Px,double* Py,double* Pz,double* hdist);
void TidyObjects(CommandLine* cl, Boresight* b,ViewVectors* v,NavBaseClass* n,ViewVectors* v2,Ellipsoid* e, DEM* d,BILWriter* bl);

//For a scan line gives the distance along view vectors to ellipsoid intersection
void GetDistanceToEllipsoid(const double X, const double Y, const double Z, const double height, CartesianVector* ECEF_vector,Ellipsoid* const ellipsoid,
                     const double hoff,double* const hdistance,const unsigned int npixels);

//Return a vector in the nadir direction
blitz::TinyMatrix<double,3,1> GetNadirVector(const double lat,const double lon);

//Function to get the scan line view vectors into cartesian ECEF XYZ
void GetScanLineViewVectorsInECEFXYZ(CartesianVector* ECEF_XYZ, ViewVectors* const vv,const double lat, const double lon,vvmethods method,
                                       const float maxallowedvvangle,uint64_t &numofbadpixels,const double theta=0,const double phi=0,const double kappa=0);

//Function to get the intersect point between the DEM and view vector
void FindIntersect(double* const px,double* const py,double* const pz,double* const seedlat,double* const seedlon,
                  Ellipsoid* ellipsoid, DEM* dem,CartesianVector* ECEF_vectors,const unsigned int pixel);

//Function to return a triangular plane created using the seed position
TriangularPlane* CreatePlaneFromNearestDEMPoints(const double* const seedlat,const double* const seedlon,Ellipsoid* ellipsoid,DEM* dem);

//Function to return a triangular plane that 'completes the square' with the plane created from the seed position
TriangularPlane* CompleteTheSquare(double* const seedlat,double* const seedlon,Ellipsoid* ellipsoid,DEM* dem);

void ShuffleSeed(double* const seedlat,double* const seedlon, DEM* dem);

//Function to set the dem aoi for reading
bool SetDEMAreaToReadIn(NavBaseClass* nav, ViewVectors* vv, DEM* dem,Ellipsoid* ellipsoid,bool quiet);

//----------------------------------------------------------------------------
// Delete objects if they exist
//----------------------------------------------------------------------------
void TidyObjects(CommandLine* cl, Boresight* b,ViewVectors* v,NavBaseClass* n,ViewVectors* v2,Ellipsoid* e, DEM* d,BILWriter* bl)
{
   if(cl!=NULL)
      delete cl;
   if(b!=NULL)
      delete b;
   if(v!=NULL)
      delete v;
   if(n!=NULL)
      delete n;
   if(v2!=NULL)
      delete v2;
   if(e!=NULL)
      delete e;
   if(d!=NULL)
      delete d;
   if(bl!=NULL)
      delete bl;
}

//----------------------------------------------------------------------------
// Delete arrays if they exist
//----------------------------------------------------------------------------
void TidyArrays(double* Plat,double* Plon,double* Pheight,double* Px,double* Py,double* Pz,double* hdist)
{
   if(Plat!=NULL)
      delete[] Plat;
   if(Plon!=NULL)
      delete[] Plon;
   if(Pheight!=NULL)
      delete[] Pheight;
   if(Px!=NULL)
      delete[] Px;
   if(Py!=NULL)
      delete[] Py;
   if(Pz!=NULL)
      delete[] Pz;
   if(hdist!=NULL)
      delete[] hdist;
}


int main(int argc,char* argv[])
{
   //----------------------------------------------------------------------------
   //The first task is to read in the command line and sort through the arguments
   //setting up the processes and performing checks on the data.
   //----------------------------------------------------------------------------

   //Set the output precision
   std::cout.precision(10);
   std::cout<<std::endl;

   //Create a pointer to the command line object
   CommandLine* cl=NULL;
   //Pointer to a boresight object
   Boresight* boresight=NULL;
   //Pointer to a leverarm object
   //Leverarm* leverarm=NULL;
   //Pointer to a view vector table object
   ViewVectors* viewvectors=NULL;
   //Pointer to a second view vector object
   ViewVectors* viewvectorsscanline=NULL;
   //Pointer to a navigation object
   NavBaseClass* navigation=NULL;
   //Pointer to an ellipsoid model
   Ellipsoid* ellipsoid=NULL;
   //Pointer to a DEM object
   DEM* dem=NULL;
   //Logging object - one for the entire program
   Logger log;
   //Output bil writer
   BILWriter* bilout=NULL;
   //the level 1 filename
   std::string strLevel1FileName;
   //Filename if atmospheric software parameters are required to be output 
   std::string strAtmosOutFilename;

   // Maximum allowed viewvector angle (degrees) - all others above this will be given the bad data value and not geocorrected
   float maxallowedvvangle=defaultmaxallowedvvangle;
   uint64_t numofbadpixels=0;

   std::stringstream strout;  //string to hold text messages in
   int retval=0; //return values stored here

   double height_offset=0; //value to use as a height offset for ellipsoid surface
   vvmethods vvmethod=NONE;//set value to largest possible - this is a flag for which view vector method to use

   std::string strViewVectorFileName; //view vector filename string
   std::string strNavFileName; //processed per-scan navigation filename string
   std::string strppoutFileName; //per-pixel position (lat/lon) filename to output to
   std::string strDEMFileName; //DEM filename

   double minlat=180,minlon=500,maxlon=-500,maxlat=-500;
   double tminlat=0,tminlon=0,tmaxlat=0,tmaxlon=0;

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
      //struct utsname information;
      //uname(&information);
      //strout<<cl->ExeName()<<" compiled at: "<<__TIME__<<" on "<<__DATE__<<std::endl<<std::endl;
      //strout<<"Operating system: "<<information.sysname<<" Release: "<<information.release<< " Version: "<<information.version<<std::endl;
      //strout<<"Machine information: "<<information.machine<<" Host name: "<<information.nodename<<" Domain name: "<<information.domainname<<std::endl;
      //strout<<cl->ExeName()<<" compiled on: "<<__MYDATE__<<std::endl<<std::endl;
      Logger::Log(strout.str());
      strout.str("");
      log.Flush();

      //----------------------------------------------------------------------------
      //check the command line options given against the list of available options in the software
      //----------------------------------------------------------------------------
      std::string badoptions;
      retval=cl->CheckAvailableOptions(availableopts,number_of_possible_options,&badoptions); 
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
         //ProgramUsage(cl->ExeName(),log);
         Logger::Log(cl->ProgramUsage(number_of_possible_options,availableopts,optsdescription));
         log.Flush();
         throw ""; //hack - change this to exit more gracefully - but then why bother?
      }

      //----------------------------------------------------------------------------
      //Check for view vector file - THIS MUST BE ON COMMAND LINE
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-vvfile"))
      {
         //Check that an argument follows the vvfile option - and get it if it exists
         if(cl->GetArg("-vvfile").compare(optiononly)!=0)
         {
            strViewVectorFileName=GetExistingFilePath(cl->GetArg("-vvfile"),true);

            Logger::Log("Will use view vector file: "+strViewVectorFileName);
         }
         else
            throw CommandLine::CommandLineException("Argument -vvfile must immediately precede the view vector filename.\n");
      }
      else
      {
         //Throw an exception
         throw CommandLine::CommandLineException("Argument -vvfile [view vector file] must be present on the command line.\n");  
      }

      //----------------------------------------------------------------------------
      //Check for processed navigation file - THIS MUST BE ON COMMAND LINE
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-navfile"))
      {
         //Check that an argument follows the navfile option - and get it if it exists
         if(cl->GetArg("-navfile").compare(optiononly)!=0)
         {
            strNavFileName=GetExistingFilePath(cl->GetArg("-navfile"),true);

            Logger::Log("Will use navigation data from file: "+strNavFileName);
         }
         else
            throw CommandLine::CommandLineException("Argument -navfile must immediately precede the processed navigation filename.\n");
      }
      else
      {
         //Throw an exception
         throw CommandLine::CommandLineException("Argument -navfile [processed navigation file] must be present on the command line.\n");  
      }

      //----------------------------------------------------------------------------
      //Check for boresight values - if not on command line then use 0,0,0
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-boresight"))
      {
         //Check that 3 arguments follows the boresight option
         if(cl->NumArgsOfOpt("-boresight") != 3 )
            throw CommandLine::CommandLineException("Error: There should be 3 arguments following the -boresight option.\n");
         
         //Get the 3 arguments into the boresight parameters
         double p=StringToDouble(cl->GetArg("-boresight",0));
         double r=StringToDouble(cl->GetArg("-boresight",1));
         double h=StringToDouble(cl->GetArg("-boresight",2));
         boresight=new Boresight(r,p,h);
         Logger::Log("Will apply boresight corrections of (R,P,H): "+ToString(r)+" "+ToString(p)+" "+ToString(h));
      }
      else
      {
         //Use default values for the boresight parameters
         boresight=new Boresight(0,0,0);
         Logger::Log("Will apply default boresight corrections of (X,Y,Z): 0 0 0");
      }

      //----------------------------------------------------------------------------
      //Check for per-pixel position output file name 
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-igmfile"))
      {
         //Check that an argument follows the igmfile option - and get it if it exists
         if(cl->GetArg("-igmfile").compare(optiononly)!=0)
         {
            //Get the filename
            strppoutFileName=CreatePath(cl->GetArg("-igmfile"));
            //Check if the file already exists - if so, then exit.
            std::ifstream testopen(strppoutFileName.c_str());
            if(testopen.is_open())
            {
               //Exit gracefully ... or not
               throw "Output file already exists. Please delete it or choose a new output file and rerun.\nFile Name: "+cl->GetArg("-igmfile");
            }
         }
         else
            throw CommandLine::CommandLineException("Argument -igmfile must immediately precede the output per-pixel position filename.\n");
      }
      else
      {
         //Throw an exception
         throw CommandLine::CommandLineException("Argument -igmfile [output per-pixel position filename] must be present on the command line.\n");  
      }

      Logger::Log("Will write per-pixel positions (longitude,latitude,height) to: "+strppoutFileName);
   
      //----------------------------------------------------------------------------
      //Check for ellipsoid model to use 
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-ellipsoid"))
      {
         //Check that an argument follows the ellipsoid option - and get it if it exists
         if(cl->GetArg("-ellipsoid").compare(optiononly)!=0)
         {   
            std::string keyword=cl->GetArg("-ellipsoid");
            if(keyword.compare("WGS84") ==0)
            {
               ellipsoid=new Ellipsoid(WGS84);
               Logger::Log("Using WGS-84 Ellipsoid.");
            }
            else
               throw CommandLine::CommandLineException("Unrecognised ellipsoid model.\n");
         }
         else
            throw CommandLine::CommandLineException("Argument -ellipsoid must immediately precede the ellipsoid model keyword.\n");
      }
      else
      {
         //Default to WGS84 ellipsoid model
         ellipsoid=new Ellipsoid(WGS84);
         Logger::Log("Using default Ellipsoid of WGS-84.");
      }

      //----------------------------------------------------------------------------
      // Is a height offset to the ellipsoid required
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-heightoffset"))
      {
         //Check that an argument follows the heightoffset option - and get it if it exists
         if(cl->GetArg("-heightoffset").compare(optiononly)!=0)
         {   
            std::string keyword=cl->GetArg("-heightoffset");
            //Check that it is a number
            if(keyword.find_first_not_of("0123456789.") == std::string::npos)
            {
               height_offset=StringToDouble(keyword);
               Logger::Log("Will use a height above the ellipsoid of "+ToString(height_offset)+"m to map to.");
            }
            else
               throw CommandLine::CommandLineException("Unrecognised heightoffset value - should be a number.\n");
         }
         else
            throw CommandLine::CommandLineException("Argument -heightoffset must immediately precede the height offset value.\n");
      }
      else
      {
         //Do not use a height offset - map to surface of ellipsoid
         height_offset=0;
         Logger::Log("Will not add a height correction to the ellipsoid surface.");
      }

      //----------------------------------------------------------------------------
      // Which view vector method is to be used - give slightly different results
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-vvSPLIT") && !cl->OnCommandLine("-vvCOMBINE"))
      {
         //use the split method
         vvmethod=SPLIT;
         Logger::Log("Will use split view vector method");
      }
      else if(cl->OnCommandLine("-vvCOMBINE") && !cl->OnCommandLine("-vvSPLIT"))
      {
         //Use the combined vv method
         vvmethod=COMBINED;
         Logger::Log("Will use combined view vector method");
      }
      else if(cl->OnCommandLine("-vvCOMBINE") && cl->OnCommandLine("-vvSPLIT"))
      {
         //Error should only specify one
         throw CommandLine::CommandLineException("Error: There should only be one vvmethod given on the command line.\n");
      }
      else
      {
         //Use default of split
         Logger::Log("Will use default view vector method: split");
         vvmethod=SPLIT;
      }

      //----------------------------------------------------------------------------
      // Is a DEM to be used - get the filename if so
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-dem"))
      {
         //Check that an argument follows the dem option - and get it if it exists
         if(cl->GetArg("-dem").compare(optiononly)!=0)  
         {
            strDEMFileName=GetExistingFilePath(cl->GetArg("-dem"),true);

            dem=new DEM(strDEMFileName);
            Logger::Log("Will use heights from Digital Elevation Model: "+strDEMFileName);
            Logger::Log("\n"+dem->Info());

         }
         else
            throw CommandLine::CommandLineException("Argument -dem must immediately precede the Digital Elevation Model filename.\n");
      }
      else
      {
         //Do not use a dem - map to surface of ellipsoid
         strDEMFileName="";
         Logger::Log("Warning: No Digital Elevation Model given. Will map to Ellipsoid surface (+offset if given).");
      }

      //-------------------------------------------------------------------------
      // Level 1 filename - needed for setting up view vectors
      //-------------------------------------------------------------------------  
      if(cl->OnCommandLine("-lev1file"))
      {
         //Check that an argument follows the lev1file option - and get it if it exists
         if(cl->GetArg("-lev1file").compare(optiononly)!=0)  
         {
            strLevel1FileName=GetExistingFilePath(cl->GetArg("-lev1file"),true);

            Logger::Log("Will trim view vectors to fit the level 1 filename: "+strLevel1FileName);
         }
         else
            throw CommandLine::CommandLineException("Argument -lev1file must immediately precede the level 1 data filename.\n");
      }
      else
      {
         throw "Level 1 filename (-lev1file) must be given on the command line.";
      }

      //-------------------------------------------------------------------------
      // Filename for writing atmospheric parameters for external software
      //-------------------------------------------------------------------------  
      if(cl->OnCommandLine("-atmosfile"))
      {
         //Check that an argument follows the atmosfile option - and get it if it exists
         if(cl->GetArg("-atmosfile").compare(optiononly)!=0)  
         {
            strAtmosOutFilename=CreatePath(cl->GetArg("-atmosfile"));
            //Check if the file already exists - if so, then exit.
            std::ifstream testopen(strAtmosOutFilename.c_str());
            if(testopen.is_open())
            {
               //Exit gracefully ... or not
               testopen.close();
               throw "Extra atmospheric parameter output file already exists. Please delete it or choose a new output file and rerun.\n File name: "+cl->GetArg("-atmosfile");
            }

            Logger::Log("Will write extra parameters out to: "+strAtmosOutFilename);
         }
         else
            throw CommandLine::CommandLineException("Argument -atmosfile must immediately precede the atmospheric parameters output data filename.\n");
      }
      else
      {
         //Assign it nothing
         strAtmosOutFilename="";
      }

      if(cl->OnCommandLine("-maxvvangle"))
      {
         //Check that an argument follows the maxvvangle option - and get it if it exists
         if(cl->GetArg("-maxvvangle").compare(optiononly)!=0)
         {   
            std::string keyword=cl->GetArg("-maxvvangle");
            //Check that it is a number
            if(keyword.find_first_not_of("0123456789.") == std::string::npos)
            {
               maxallowedvvangle=static_cast<float>(StringToDouble(keyword));
               Logger::Log("Will use a maximum view vector of "+ToString(maxallowedvvangle)+"degrees.");
            }
            else
               throw CommandLine::CommandLineException("Unrecognised maxvvangle value - should be a number <90 in degrees.\n");
         }
         else
            throw CommandLine::CommandLineException("Argument -maxvvangle must immediately precede the maximum angle in degrees value.\n");
      }
      //*****************************************************************
      // ENTER NEW COMMAND LINE OPTION CODE HERE
      //*****************************************************************

   }
   catch(CommandLine::CommandLineException e)
   {
      Logger::Error(std::string(e.what())+"\n"+e.info);
      TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
      exit(1);
   }
   catch(char const* e)
   {
      Logger::Error(e);
      TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
      exit(1);      
   }
   catch(std::string e)
   {
      Logger::Error(e);
      TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
      exit(1);      
   }
   catch(BinaryReader::BRexception e)
   {
      Logger::Error(std::string(e.what())+"\n"+e.info);
      TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
      exit(1); //exit the program
   }
   catch(std::exception& e)
   {
      PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString(),&e);
      TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
      exit(1);
   }
   catch(...)
   {
      PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString());
      TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
      exit(1);
   } 

   //Flush the log here
   log.Flush();   

   //Convert the maxallowedvvangle to radians
   maxallowedvvangle=maxallowedvvangle*PI/180;

   //----------------------------------------------------------------------------
   //Set up the view vectors 
   //----------------------------------------------------------------------------
   try
   {
      //Create a view vector object from the view vector file
      Logger::Log("Creating view vector object.");
      viewvectors=new ViewVectors(strViewVectorFileName,strLevel1FileName);
      //Update the view vectors with the boresight values
      viewvectors->ApplyAngleRotations(boresight->Roll(),boresight->Pitch(),boresight->Heading());
   }
   catch(BinaryReader::BRexception e)
   {
      Logger::Error(std::string(e.what())+"\n"+e.info);
      TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
      exit(1); //exit the program
   }
   catch(char const* e)
   {
      Logger::Error(e);
      TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
      exit(1);      
   }
   catch(...)
   {
      PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString());
      TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
      exit(1);
   } 
   
   //----------------------------------------------------------------------------
   //Create a navigation object and open the nav file for reading
   //----------------------------------------------------------------------------
   try
   {
      navigation=new NavBaseClass(strNavFileName);
   }
   catch(BinaryReader::BRexception e)
   {
      Logger::Error(std::string(e.what())+"\n"+e.info);
      TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
      exit(1); //exit the program
   }
   catch(char const* e)
   {
      Logger::Error(e);
      TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
      exit(1);
   }
   catch(...)
   {
      PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString());
      TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
      exit(1);
   } 

   //----------------------------------------------------------------------------
   // Create a copy of the view vectors to use for per scan line offsets
   //----------------------------------------------------------------------------
   viewvectorsscanline=new ViewVectors(*viewvectors);

   //----------------------------------------------------------------------------
   // Set up some variables for the main loop
   //----------------------------------------------------------------------------

   double X=0,Y=0,Z=0; //Aircraft ECEF XYZ position
   double lat=0,lon=0,hei=0;//Aircraft LLH position

   //ECEF view vectors
   CartesianVector* ECEF_vectors=NULL;

   //Position on surface in X,Y,Z
   double* Px=NULL;
   double* Py=NULL;
   double* Pz=NULL;

   //Geodetic lat/lon/hei of the new pixel position
   double* Plat=NULL;
   double* Plon=NULL;
   double* Pheight=NULL;
   //distance from sensor to surface
   double* hdist=NULL; 

   //Surface intersect points in ECEF
   Px=new double[viewvectorsscanline->NumberItems()];
   Py=new double[viewvectorsscanline->NumberItems()];
   Pz=new double[viewvectorsscanline->NumberItems()];
   //Surface intersect points in LLH
   Plat=new double[viewvectorsscanline->NumberItems()];
   Plon=new double[viewvectorsscanline->NumberItems()];
   Pheight=new double[viewvectorsscanline->NumberItems()];
   //Distance from aircraft to ellipsoid surface for each scan line pixel
   hdist=new double[viewvectorsscanline->NumberItems()];

   //Variables for DEM processing loops
   double seedlat=0; //seed position latitude
   double seedlon=0; //seed position longitude
   //Nadir vector
   blitz::TinyMatrix<double,3,1> myNadir;
   //index for vv closest to nadir
   unsigned int nadirindex=0; 
   try
   {
      //----------------------------------------------------------------------------
      //Output some information about the navigation data before going into the loop
      //----------------------------------------------------------------------------
      navigation->ReadScan(0);
      lat=navigation->Lat();
      lon=navigation->Lon();
      hei=navigation->Hei();
      Logger::Log("\nAircraft start position (Lon,Lat,Hei): "+ToString(lon)+" "+ToString(lat)+" "+ToString(hei));
      Logger::Log("Aircraft navigation start time: "+ToString(navigation->Time()));
      navigation->ReadScan(navigation->TotalScans()-1);
      lat=navigation->Lat();
      lon=navigation->Lon();
      hei=navigation->Hei();
      Logger::Log("\nAircraft end position (Lon,Lat,Hei): "+ToString(lon)+" "+ToString(lat)+" "+ToString(hei));
      Logger::Log("Aircraft navigation end time: "+ToString(navigation->Time()));
      Logger::Log("Total number of navigation scan lines to map: "+ToString(navigation->TotalScans()));

      //Test extent of nav file here and log the min/max
      navigation->FindLimits();
      Logger::Log("\nNavigation Min/Max Latitude: "+ToString(navigation->MinLat())+" "+ToString(navigation->MaxLat()));
      Logger::Log("Navigation Min/Max Longitude: "+ToString(navigation->MinLon())+" "+ToString(navigation->MaxLon()));
      Logger::Log("Navigation Min/Max Height: "+ToString(navigation->MinHei())+" "+ToString(navigation->MaxHei()));
      Logger::Log("Navigation Min/Max Roll: "+ToString(navigation->MinRoll())+" "+ToString(navigation->MaxRoll()));
      //Flush log to output before entering loop
      log.Flush();

   }
   catch(BinaryReader::BRexception e)
   {
      Logger::Error(std::string(e.what())+"\n"+e.info);
      TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
      exit(1); //exit the program
   }
   catch(char const* e)
   {
      Logger::Error(e);
      TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
      exit(1);
   }
   catch(...)
   {
      PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString());
      TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
      exit(1);
   } 

   //----------------------------------------------------------------------------
   //If DEM is to be used then do the following section of code
   //----------------------------------------------------------------------------
   std::vector<unsigned int> sectionscanlimits;
   unsigned int lowerscan=0,upperscan=0;
   if(strDEMFileName!="")
   {
      try
      {
         //Get DEM bounds for area to read in
         bool DEMAREAOK=SetDEMAreaToReadIn(navigation,viewvectorsscanline,dem,ellipsoid,false);

         //If the DEM is not big enough for some reason
         if(!DEMAREAOK)
         {
            Logger::Log("WARNING: It appears that the DEM does not cover the area of the navigation file.");
            Logger::Log("Exiting...");
            TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
            TidyArrays(Plat,Plon,Pheight,Px,Py,Pz,hdist);
            exit(1);      
         }
         else
         {
            //Now read in the required portion from the DEM
            Logger::Log("Reading in DEM will require approx. memory (MB) of: "+ToString(dem->SizeOf()/(1024.0*1024.0)));
         }

         //Check size to be read in is not ridiculously large to cause memory issues
         //If it is too large then we may need to reduce the size of it.

         //This AOI is too large so we need to reduce the size of it
         //Loop until we have a list of limits that describe each subsection of the DEM AOI to be read in
         lowerscan=0;
         upperscan=navigation->TotalScans();
         //Logger::Log("\nSplitting up processing into smaller sections to keep DEM RAM usage down.\n");            

         //Enter loop to possibly split processing up if DEM AOI is too large for RAM (32-bit machines)
         while(lowerscan != navigation->TotalScans())
         {
            //Break for while loop is within loop - when sizeof the DEM AOI is less than 2GB
            while(true)
            {
               //Calculate min/max limits for this section of navigation
               navigation->FindLimits(lowerscan,upperscan);
               //Now try DEM AOI with these limits
               DEMAREAOK=SetDEMAreaToReadIn(navigation,viewvectorsscanline,dem,ellipsoid,true);     
               if(!DEMAREAOK)
                  Logger::Error("DEM AOI is not OK - This should never happen and is a bug - please notify ARSF.");

               //If less than 2GB we're happy so break out of this loop
               if(dem->SizeOf()/(1024.0*1024.0) < 2048)
                  break;

               Logger::WarnOnce("Original DEM AOI is too large (2GB enforced limit for 32-bit support). Will split up and do processing in chunks.");

               //Update the upper scan limit to be smaller 
               //but we need to ensure upperscan will be > lowerscan, and that upperscan values depend on lowerscan values
               //e.g. if 250 lines is maximum size, and initial lowerscan=0 upperscan=1000:
               // then 0,250 ok,  
               //      250,500 ok, 
               //      500 500 not ok.
               upperscan=upperscan - 0.5*(upperscan-lowerscan);

            }
            //Limits passed the test - AOI is less than 2GB for this subsection
            //Store the upper/lower limits
            Logger::Log("Section to be processed using scan bounds: "+ToString(lowerscan)+" : "+ToString(upperscan));
            Logger::Log("will require approx. memory (MB) of: "+ToString(dem->SizeOf()/(1024.0*1024.0)));
            sectionscanlimits.push_back(lowerscan);
            sectionscanlimits.push_back(upperscan);

            //Now repeat for next AOI subsection
            lowerscan=upperscan;
            upperscan=navigation->TotalScans();
         }
      }
      catch(const char* e)
      {
         Logger::Error(e);
         TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
         exit(1);
      }
      catch(std::string e)
      {
         Logger::Error(e);
         TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
         exit(1);
      }
      catch(std::bad_alloc& e)
      {
         Logger::Error("Exception: trying to allocate more RAM than is available. Current work around - use a lower resolution (in lat/lon) DEM.");
         TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
         exit(1);
      }
      catch(std::exception& e)
      {
         PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString(),&e);
         TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
         exit(1);
      }
      catch(...)
      {
         PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString());
         TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
         exit(1);
      }      
   }
   else
   {
      Logger::Log("Warning - no Digital Elevation Model was given on command line. Will map to ellipsoid surface.");
      sectionscanlimits.push_back(0);
      sectionscanlimits.push_back(navigation->TotalScans());  
   }


   //----------------------------------------------------------------------
   // Create a new BILwriter to output the results with
   //----------------------------------------------------------------------
   try
   {
      bilout=new BILWriter(strppoutFileName,BILWriter::float64,navigation->TotalScans(),viewvectorsscanline->NumberItems(),3,'a');
      bilout->AddToHdr("projection = Geographic Lat/Lon");
      bilout->AddToHdr("datum ellipsoid = "+ellipsoid->Name());
      //Add band names
      bilout->AddToHdr("band names = {Longitude, Latitude, Height}");
      //Add the x start and y start in case further mapping in ENVI is required
      //as the IGM should match the level 1 (without x/y start ENVI misunderstands)
      BinFile lev1(strLevel1FileName);
      std::string xstart=lev1.FromHeader("x start");
      std::string ystart=lev1.FromHeader("y start");
      lev1.Close();
      bilout->AddToHdr(";These describe which pixels from the original raw image the IGM file positions relate to.");
      bilout->AddToHdr("x start = "+xstart);
      bilout->AddToHdr("y start = "+ystart);
      bilout->AddToHdr("data ignore value = "+ToString(BADDATAVALUE));

   }
   catch(BILWriter::BILexception e)
   {
      Logger::Error(std::string(e.what())+"\n"+e.info);
      TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
      exit(1);
   }
   catch(std::exception& e)
   {
      PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString(),&e);
      TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
      exit(1);
   }
   catch(...)
   {
      PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString());
      TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
      exit(1);
   } 

   //----------------------------------------------------------------------------
   //Loop through the navigation file reading per scan
   //Note this now takes into account possibility of multiple DEM readings
   //----------------------------------------------------------------------------

   for(std::vector<unsigned int>::iterator it=sectionscanlimits.begin();it!=sectionscanlimits.end();it=it+2)
   {
      //Get the lower/upper scan limits for this subsection of DEM
      lowerscan=*it;
      upperscan=*(it+1);
      Logger::Log("Processing section with scan bounds: "+ToString(lowerscan)+" : "+ToString(upperscan));
      //If we are using a DEM then we need to set the AOI to match this region defined by the lower/upper scans
      if(strDEMFileName!="")
      {
         try
         {
            //Set DEM AOI to this subsection
            navigation->FindLimits(lowerscan,upperscan);
            SetDEMAreaToReadIn(navigation,viewvectorsscanline,dem,ellipsoid,true);
            //Read in DEM data
            dem->FillArray();
         }
         catch(char const* e)
         {
            Logger::Error(e);
            TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
            TidyArrays(Plat,Plon,Pheight,Px,Py,Pz,hdist);
            exit(1);
         }
         catch(std::string e)
         {
            Logger::Error(e);
            TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
            TidyArrays(Plat,Plon,Pheight,Px,Py,Pz,hdist);
            exit(1);
         }
      }
      for(unsigned int scan=lowerscan;scan<upperscan;scan++)
      {
         Logger::Verbose("Starting scan: "+ToString(scan));
         //Read in the current scan
         navigation->ReadScan(scan);
         //update the view vectors accordingly
         if(vvmethod==COMBINED)
            viewvectorsscanline->ApplyAngleRotations(navigation->Roll(),navigation->Pitch(),navigation->Heading());

         //get the lat/lon/hei of aircraft 
         lat=navigation->Lat();
         lon=navigation->Lon();
         hei=navigation->Hei();
         Logger::Verbose("Aircraft position: latitude: "+ToString(lat)+" longitude: "+ToString(lon)+" height: "+ToString(hei));

         //Now convert the aircraft pos to ECEF XYZ
         ConvertLLH2XYZ(&lat, &lon, &hei, &X, &Y, &Z, 1, GEODETIC,ellipsoid);

         //Create a cartesian vector object to store the ecef vectors in, with an origin at the aircraft
         ECEF_vectors=new CartesianVector(viewvectorsscanline->NumberItems(),X,Y,Z);

         //Convert the view vectors into earth centred earth fixed cartesians
         if(vvmethod==COMBINED)
         {
            GetScanLineViewVectorsInECEFXYZ(ECEF_vectors,viewvectorsscanline,lat,lon,vvmethod,maxallowedvvangle,numofbadpixels);
         }
         else if(vvmethod==SPLIT)
         {
            GetScanLineViewVectorsInECEFXYZ(ECEF_vectors,viewvectorsscanline,lat,lon,vvmethod,maxallowedvvangle,numofbadpixels,navigation->Roll(),navigation->Pitch(),navigation->Heading());
         }

         //If Ellipsoid mapping has been requested on command line then map to the ellipsoid
         //Currently this depends on whether a DEM was given on the command line
         if(strDEMFileName=="") //no dem given on command line
         {
            GetDistanceToEllipsoid(X,Y,Z,hei,ECEF_vectors,ellipsoid,height_offset,hdist,viewvectorsscanline->NumberItems());
            //For each of the pixels per scan line, project onto the ellipsoid surface
            for(unsigned int p=0;p<viewvectorsscanline->NumberItems();p++)
            {
               if(hdist[p]==BADDATAVALUE)
               {
                  Px[p]=Py[p]=Pz[p]=BADDATAVALUE;
               }
               else
               {
                  Px[p]=X+ECEF_vectors->X[p]*hdist[p];
                  Py[p]=Y+ECEF_vectors->Y[p]*hdist[p];
                  Pz[p]=Z+ECEF_vectors->Z[p]*hdist[p];
               }
            }
         }
         else
         //Else if DEM has been supplied and requested for use in mapping
         {
            try
            {
               //For this scan line - get the nadir vector and find the vv closest to nadir            
               myNadir=GetNadirVector(lat,lon);
               nadirindex=ECEF_vectors->GetNadirIndex(myNadir);
               //Logger::Debug("Nadir Vector: X:"+ToString(myNadir(0,0))+" Y:"+ToString(myNadir(1,0))+" Z:"+ToString(myNadir(2,0)));
               //Logger::Debug("ViewVector index closest to nadir: "+ToString(nadirindex));
               //Logger::Debug("Vector ("+ToString(nadirindex)+"): "+ToString(ECEF_vectors->X[nadirindex])+" "+ToString(ECEF_vectors->Y[nadirindex])
               //            +" "+ToString(ECEF_vectors->Z[nadirindex]));

               //Set up seed position for DEM intersection search - use aircraft lat/lon   
               seedlat=lat;
               seedlon=lon;
               //"Shuffle" these to add on an offset if required - only if this falls on a boundary line of the dem,
               //we want to shift it so that it is definitely on one side of the boundary (either of the sides is ok)
               ShuffleSeed(&seedlat,&seedlon,dem);

               //For each pixel from nadir index to end of ccd
               for(unsigned int pixel=nadirindex;pixel<viewvectorsscanline->NumberItems();pixel++)
               {
                  //FindIntersect(dem,&Px[pixel],&Py[pixel],&Pz[pixel],&seedlat,&seedlon,ellipsoid,ECEF_vectors,pixel);           
                  FindIntersect(&Px[pixel],&Py[pixel],&Pz[pixel],&seedlat,&seedlon,ellipsoid,dem,ECEF_vectors,pixel);
               }
               //Reset seed lat/lon to aircraft/nadir position
               seedlat=lat;
               seedlon=lon;
               //"Shuffle" these to add on an offset if required - only if this falls on a boundary line of the dem,
               //we want to shift it so that it is definitely on one side of the boundary (either of the sides is ok)
               ShuffleSeed(&seedlat,&seedlon,dem);

               //Now go through the other viewvectors i.e. nadir index to start of ccd
               for(int pixel=nadirindex-1;pixel>=0;pixel--)
               {
                  //FindIntersect(dem,&Px[pixel],&Py[pixel],&Pz[pixel],&seedlat,&seedlon,ellipsoid,ECEF_vectors,pixel);             
                  FindIntersect(&Px[pixel],&Py[pixel],&Pz[pixel],&seedlat,&seedlon,ellipsoid,dem,ECEF_vectors,pixel);
               }
            }
            catch(char const* e)
            {
               Logger::Error(e);
               TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
               TidyArrays(Plat,Plon,Pheight,Px,Py,Pz,hdist);
               exit(1);
            }
            catch(std::string e)
            {
               Logger::Error(e);
               TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
               TidyArrays(Plat,Plon,Pheight,Px,Py,Pz,hdist);
               exit(1);
            }
            catch(std::exception& e)
            {
               PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString(),&e);
               TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
               exit(1);
            }
            catch(...)
            {
               PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString());
               TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
               exit(1);
            } 
         }

         //Now convert the new pixel position to LLH for each pixel of scanline - THESE ARE RETURNED AS RADIANS
         ConvertXYZ2LLH(Px, Py, Pz,Plat, Plon, Pheight,viewvectorsscanline->NumberItems(),GEODETIC,ellipsoid,BADDATAVALUE);

         //if atmospheric correction software parameters are to be output - do this here
         if(strAtmosOutFilename.compare("")!=0)
         {
            const int nbandsatmosfile=5;
            //Create arrays to hold the dem slope and aspect for this scan line worth of data
            double* atmosout=new double[viewvectorsscanline->NumberItems()*nbandsatmosfile];
            //constant pointers for simplicty pointing into atmosout parameter
            double* const azimuth=&atmosout[0];
            double* const zenith=&atmosout[viewvectorsscanline->NumberItems()];
            double* const distance=&atmosout[viewvectorsscanline->NumberItems()*2];
            double* const demslope=&atmosout[viewvectorsscanline->NumberItems()*3];
            double* const demaspect=&atmosout[viewvectorsscanline->NumberItems()*4];

            //Get the view vectors in azimuth,zenith 
            //For each of the pixels in the arrays, calculate the geodesic distance, azimuith and zenith
            for(unsigned int i=0;i<viewvectorsscanline->NumberItems();i++)
            {
               GetGeodesicDistance_Bowring(Plon[i],Plat[i],Pheight[i],lon*PI/180,lat*PI/180,hei,distance[i],azimuth[i],zenith[i],ellipsoid);
            }

            //Get the DEM slope and aspect values for cells in which the intersect is contained
            if(strDEMFileName.compare("")!=0)
            {
               dem->CalculateSlopeAndAzimuth(Plat,Plon,demslope,demaspect,viewvectorsscanline->NumberItems());
            }
            else
            {
               //No dem is being used - set to zero for outputs
               for(unsigned int iv=0;iv<viewvectorsscanline->NumberItems();iv++)
               {
                  demslope[iv]=0;
                  demaspect[iv]=0;
               }
            }

            //Open up a bil file to write to
            try
            {
               BILWriter* dataout=new BILWriter(strAtmosOutFilename,BILWriter::float64,navigation->TotalScans(),viewvectorsscanline->NumberItems(),nbandsatmosfile,'a');
               //Add band names
               dataout->AddToHdr("band names = {View azimuth, View zenith, Distance, DEM slope, DEM aspect}");
               dataout->AddToHdr(";View azimuth and DEM aspect (azimuth) are measured clockwise from North in degrees.");
               dataout->AddToHdr(";View zenith is measured in degrees from the vertical to the nadir.");
               dataout->AddToHdr(";DEM slope is measured in degrees from the horizontal.");
               dataout->AddToHdr(";Distance is the distance from sensor to ground intersect and measured in metres.");
               //Write out the data
               dataout->WriteLine((char*)atmosout);
               //Tidy up 
               dataout->Close();
               delete dataout;
            }
            catch(BILWriter::BILexception e)
            {
               //At the moment just say what went wrong
               Logger::Error(std::string(e.what())+"\n"+e.info);
            }
            catch(std::exception& e)
            {
               PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString(),&e);
               TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
               exit(1);
            }
            catch(...)
            {
               PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString());
               TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
               exit(1);
            } 

            delete[] atmosout;
         }

         //Output the pixel lat/lon/hei arrays
         try
         {
            //Need to convert from radians to degrees
            for(unsigned int j=0;j<viewvectorsscanline->NumberItems();j++)
            {
               if((Plon[j]!=BADDATAVALUE)||(Plat[j]!=BADDATAVALUE))
               {
                  Plon[j]=Plon[j]*180/PI;
                  Plat[j]=Plat[j]*180/PI;
               }
            }

            //Keep track of min/max lat and longs here
            tmaxlat=*std::max_element(Plat,Plat+viewvectorsscanline->NumberItems());
            tminlat=*std::min_element(Plat,Plat+viewvectorsscanline->NumberItems());
            tmaxlon=*std::max_element(Plon,Plon+viewvectorsscanline->NumberItems());
            tminlon=*std::min_element(Plon,Plon+viewvectorsscanline->NumberItems());
            //Update variables if a new min or max is found
            if(tmaxlat > maxlat)
               maxlat=tmaxlat;
            if(tminlat < minlat)
               minlat=tminlat;
            if(tmaxlon > maxlon)
               maxlon=tmaxlon;
            if(tminlon < minlon)
               minlon=tminlon;

            //Create a BIL writer to handle the output
            //BILWriter bil(strppoutFileName,BILWriter::float64,navigation->TotalScans(),viewvectorsscanline->NumberItems(),3,'a');
            bilout->WriteBandLine((char*)Plon);
            bilout->WriteBandLine((char*)Plat);
            bilout->WriteBandLine((char*)Pheight);
            //bil.Close();
         }
         catch(BILWriter::BILexception e)
         {
            //At the moment just say what went wrong
            Logger::Error(std::string(e.what())+"\n"+e.info);
         }
         catch(std::exception& e)
         {
            PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString(),&e);
            TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
            exit(1);
         }
         catch(...)
         {
            PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString());
            TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);
            exit(1);
         } 

         //Delete the current copy of viewvectorsscanline and replace with a copy of viewvectors 
         delete viewvectorsscanline;
         viewvectorsscanline=new ViewVectors(*viewvectors);
         delete ECEF_vectors;

         //Percent done counter
         PercentProgress(scan,navigation->TotalScans());

      }
   }

   Logger::Log("Geocorrection processing completed. \n\n");

   //Output the min/max to the bil header
   bilout->AddToHdr(";Min X = "+ToString(minlon));
   bilout->AddToHdr(";Max X = "+ToString(maxlon));
   bilout->AddToHdr(";Min Y = "+ToString(minlat));
   bilout->AddToHdr(";Max Y = "+ToString(maxlat));

   bilout->Close();

   //Output the number of bad pixels
   if(numofbadpixels>0)
   {
      Logger::Warning("There were some pixels which were not mapped because their view vector angle was greater"
                      " than the maximum allowed (set by -maxvvangle. Total number: "+ToString(numofbadpixels));
   }
   TidyArrays(Plat,Plon,Pheight,Px,Py,Pz,hdist);
   TidyObjects(cl,boresight,viewvectors,navigation,viewvectorsscanline,ellipsoid,dem,bilout);

}


//-------------------------------------------------------------------------
// Derive the view vectors in ECEF cartesian coordinates by applying rotational
// transformations. THETA, PHI, KAPPA are the navigation attitude values
// These can be ignored if method is COMBINED
//-------------------------------------------------------------------------
void GetScanLineViewVectorsInECEFXYZ(CartesianVector* ECEF_XYZ, ViewVectors* const vv,const double lat, const double lon,
                                       vvmethods method,const float maxallowedvvangle,uint64_t &numofbadpixels,const double theta,const double phi,const double kappa)
{
   //Store the returned XYZ from the GetVVinECEFXYZ function
   double ECEFXYZ[3]={0,0,0};
   double roll=0,pitch=0,heading=0;

   blitz::TinyMatrix<double,3,1> nadir;
   nadir=0,0,1;
  
   //We need a nadir unit vector in ECEF for the platform
   //Use the aircraft position vector - not perfect as assumes a spheroidal Earth but may be good enough
   //The aircraft position is the origin in the ECEF vector array
   double down_vector[3]={-ECEF_XYZ->OriginX(),-ECEF_XYZ->OriginY(),-ECEF_XYZ->OriginZ()};
   double down_vector_magnitude=sqrt(down_vector[0]*down_vector[0] + down_vector[1]*down_vector[1] + down_vector[2]*down_vector[2]);
   double unit_down_vector[3]={down_vector[0]/down_vector_magnitude,down_vector[1]/down_vector_magnitude,down_vector[2]/down_vector_magnitude};   
   double cos_al=0;
   //For each pixel in the scan line, convert the view vector to ECEF XYZ
   for(unsigned int pixel=0;pixel<vv->NumberItems();pixel++)
   {
      //These are angles that rotate nadir to local level: sensor look ( + aircraft attitude if COMBINED)
      roll=vv->rotX[pixel];
      pitch=vv->rotY[pixel];
      heading=vv->rotZ[pixel];

      //Get the nadir view vector in ECEF XYZ
      if(method==COMBINED)
      {
         GetVVinECEFXYZ(&nadir,ECEFXYZ,lat,lon,roll,pitch,heading);
      }
      else if(method==SPLIT)
      {
         GetVVinECEFXYZ(&nadir,ECEFXYZ,lat,lon,roll,pitch,heading,theta,phi,kappa);         
      }

      //Test if view vector is above the horizon (eeek) or not (phew)
      //dot nadir vs vv and test the angle is less than maxallowedvvangle
      cos_al=unit_down_vector[0]*ECEFXYZ[0] + unit_down_vector[1]*ECEFXYZ[1] + unit_down_vector[2]*ECEFXYZ[2];
      if((cos_al < 0) || (acos(cos_al) > maxallowedvvangle))
      {
         //BAD?
         ECEF_XYZ->X[pixel]=BADDATAVALUE;
         ECEF_XYZ->Y[pixel]=BADDATAVALUE;      
         ECEF_XYZ->Z[pixel]=BADDATAVALUE; 
         numofbadpixels++;  
      }
      else
      {
         //Now update the returned ECEF viewvectors
         ECEF_XYZ->X[pixel]=ECEFXYZ[0];
         ECEF_XYZ->Y[pixel]=ECEFXYZ[1];      
         ECEF_XYZ->Z[pixel]=ECEFXYZ[2];  
      }    
   }
}


//----------------------------------------------------------------------------
//Function to calculate the distance from aircraft to ellipsoid surface and return it
//Need to find vector magnitude to interesect with the ellipsoid surface - solve the quadratic eqn to get h
//----------------------------------------------------------------------------
void GetDistanceToEllipsoid(const double X, const double Y, const double Z, const double height,CartesianVector* ECEF_vector,
            Ellipsoid* const ellipsoid,const double hoff,double* const hdistance,const unsigned int npixels)
{
   //Get the semi-major/minor axes squared and take into account the height offset for surface
   double aa=(ellipsoid->a()*ellipsoid->a()) + 2*hoff*ellipsoid->a() + hoff*hoff;
   double bb=(ellipsoid->b()*ellipsoid->b()) + 2*hoff*ellipsoid->b() + hoff*hoff; 

   //Quadratic equation variables
   double A=0,B=0,C=0;
   //Possible height distance
   double h1=0,h2=0;
   //Test variables
   double h1test=0,h2test=0;

   //For each pixel
   for(unsigned int pixel=0;pixel<npixels;pixel++)
   {
      if((ECEF_vector->X[pixel]==BADDATAVALUE)||(ECEF_vector->Y[pixel]==BADDATAVALUE)||(ECEF_vector->Z[pixel]==BADDATAVALUE))
      {
         hdistance[pixel]=BADDATAVALUE;   
         continue;
      }

      //Set up the quadratic equation components
      A=(ECEF_vector->X[pixel]*ECEF_vector->X[pixel] + ECEF_vector->Y[pixel]*ECEF_vector->Y[pixel]) / aa + (ECEF_vector->Z[pixel]*ECEF_vector->Z[pixel])/bb;
      B=2*(X*ECEF_vector->X[pixel] + Y*ECEF_vector->Y[pixel])/aa + (2*Z*ECEF_vector->Z[pixel])/bb;
      C=(X*X+Y*Y)/aa+ (Z*Z)/bb -1;

      //The two possible distances from aircraft to surface
      h1=(-B+sqrt(B*B -4*A*C))/(2*A);
      h2=(-B-sqrt(B*B -4*A*C))/(2*A);

      //Test which is correct solution - we want closest one to the height
      //i.e. smallest of h1test/h2test gives the correct height
      h1test=(h1-height)*(h1-height);
      h2test=(h2-height)*(h2-height);

      if(h1test < h2test)
         hdistance[pixel]=h1;
      else
         hdistance[pixel]=h2;

      //We should check for negative distances here
      if(hdistance[pixel]<0)
         hdistance[pixel]=BADDATAVALUE;
   }
}

//----------------------------------------------------------------------------
//Calculate the nadir vector at aircraft position
//----------------------------------------------------------------------------
blitz::TinyMatrix<double,3,1> GetNadirVector(const double lat,const double lon) 
{
   //Reframe vector - nadir is defined as +ve z axis in aircrft coordinates 
   blitz::TinyMatrix<double,3,1> nadir;
   nadir=0,0,1;

   //Store the result in an array
   double ECEFXYZ[3]={0,0,0};
   //Transform the nadir in aircraft coords to ECEF coords
   GetVVinECEFXYZ(&nadir,ECEFXYZ,lat,lon,0,0,0);
   //update the nadir vactor and return this
   nadir=ECEFXYZ[0],ECEFXYZ[1],ECEFXYZ[2];
   return nadir;
}


//-------------------------------------------------------------------------
// Find an intersect between dem and a view vector
// The result in ECEF cartesians is returned in px,py,pz, seedlat and seelon 
// are updated to be the intersect in ellipsiodal lat/lon
// ECEF_vectors[pixel] is the view vector
//-------------------------------------------------------------------------
void FindIntersect(double* const px,double* const py,double* const pz,double* const seedlat,double* const seedlon,
                  Ellipsoid* ellipsoid, DEM* dem,CartesianVector* ECEF_vectors,const unsigned int pixel)
{
   //Get the 3 nearest points from the DEM to the seed position
   //Convert these into ECEF XYZ
   //Create a planar surface from these 3 points
   TriangularPlane* triplane=CreatePlaneFromNearestDEMPoints(seedlat,seedlon,ellipsoid,dem);
   if(triplane==NULL)
      throw "DEM does not cover the entire flight line (Actually - could not find an intersect with the DEM, so could be due to other issues too).\n";      

   //Get 2 points that the view vector passes through
   double pvX[2],pvY[2],pvZ[2];
   double plat=0,plon=0,pheight=0;
   //Point 1 is the Origin itself
   pvX[0]=ECEF_vectors->OriginX();
   pvY[0]=ECEF_vectors->OriginY();
   pvZ[0]=ECEF_vectors->OriginZ();
   //Point 2 is the vector + Origin
   pvX[1]=ECEF_vectors->OriginX() + ECEF_vectors->X[pixel];
   pvY[1]=ECEF_vectors->OriginY() + ECEF_vectors->Y[pixel];
   pvZ[1]=ECEF_vectors->OriginZ() + ECEF_vectors->Z[pixel];

   //Set up loop parameters including a counter and initial x,y
   //values and directions (dx and dy)
   unsigned int loopcounter=0;
   int64_t x=0,y=0;
   int dx=0,dy=-1,tmpdx=0;

   double origseedlon= *seedlon;
   double origseedlat= *seedlat;

//   static int debugcalls=0;
//   debugcalls++;

   //Test if vector and plane intersect within the bounds of the plane
   while(!(triplane->Intersect(pvX,pvY,pvZ,px,py,pz)))
   {
//      if(debugcalls > 3127*1024)// 2837*512)
//      {
//         std::cout<<"Pixel: "<<pixel<<" lat,lon: "<<*seedlat<<" "<<*seedlon<<std::endl;
//         std::cout<<"Point1: "<<triplane->Point1()[0]<<" "<<triplane->Point1()[1]<<" "<<triplane->Point1()[2]<<std::endl;
//         std::cout<<"Point2: "<<triplane->Point2()[0]<<" "<<triplane->Point2()[1]<<" "<<triplane->Point2()[2]<<std::endl;
//         std::cout<<"Point3: "<<triplane->Point3()[0]<<" "<<triplane->Point3()[1]<<" "<<triplane->Point3()[2]<<std::endl;
//         std::cout<<"Intersect: "<<px[0]<<" "<<py[0]<<" "<<pz[0]<<std::endl;
//      }

      //No - Update the planar surface and try again
      delete triplane;
      triplane=NULL;

      if((loopcounter % 2)==0)
      {
         //If this is an even numbered iteration of the loop then we should
         //"complete the square" - this means use the triangle plane that would
         //form a square together with the previous triangle plane
         triplane=CompleteTheSquare(seedlat,seedlon,ellipsoid,dem);

      }
      else
      {
         //This is an odd numbered iteration of the loop - so we should get
         //the triangular plane that is further along the travel path on the dem.
         //This is in a while loop that repeats until a triangular plane is created
         //that is fully contained within the DEM AOI.
         while(triplane==NULL)
         {
            //If any of these conditions are met we need to chnage direction (anticlockwise)
            if((x==y)||((x<0)&&(x==-y))||((x>0)&&(x==1-y)))
            {
               //Update the movement direction through array
               tmpdx=dx;
               dx=-dy;
               dy=tmpdx;
            }
            //Update cell position of array query
            x=x+dx;
            y=y+dy;

            //Update the seed position - this will accumulate floating point error
            //may be better to use, e.g., newlon=origseedlon + x*dem->GetXSpace();
            //*seedlon=(*seedlon)+dx*dem->GetXSpace();
            //*seedlat=(*seedlat)+dy*dem->GetYSpace();

            //Multiplying by 0.99 to cater for rare cases where seed point is close to 
            //dem bound and adding 100% of X/YSpace takes seed point more than 1 dem cell away
            //due to floating point rounding. Should not have much of an affect unless there are
            //many loops in which case a grid cell maybe checked twice.
            *seedlon=origseedlon+x*dem->GetXSpace()*0.99;
            *seedlat=origseedlat+y*dem->GetYSpace()*0.99;

            //Create a new plane from this seed position
            triplane=CreatePlaneFromNearestDEMPoints(seedlat,seedlon,ellipsoid,dem);

//            if(debugcalls > 3127*1024) // 2837*512)
//            {
//               std::cout<<"Pixel: "<<pixel<<" X,Y: "<<x<<" "<<y<<" lat,lon: "<<*seedlat<<" "<<*seedlon<<" "<<triplane<<std::endl;
//            }
         }  

      }
      //The loop counter is only used to keep track of even/odd iterations
      //it could also be used in the while() definition to prevent infinite loops
      loopcounter++;
   }

   //Yes - return x,y,z position of intersect [implicit in function definition]
   //Update seed position to be lat/lon of this intersect
   ConvertXYZ2LLH(px, py, pz,&plat, &plon, &pheight,1,GEODETIC,ellipsoid);
   *seedlat=plat*180/PI;
   *seedlon=plon*180/PI;

   //delete the triplane object before leaving this function
   delete triplane;
   triplane=NULL;
}


//-------------------------------------------------------------------------
// Function to return a triangular plane that is the opposite half of the
// square of the current triangle plane
//
// e.g. If we currently have plane made of vertices ACD, then we return
//      the plane constructed of vertices ABC
//       A                       A     B
//
//       D     C                       C
//
//-------------------------------------------------------------------------
TriangularPlane* CompleteTheSquare(double* const seedlat,double* const seedlon,Ellipsoid* ellipsoid,DEM* dem)
{
   //3 Points form the plane - use the seed point to start
   double pplon[3]={0};
   double pplat[3]={0};
   double pphei[3]={0};
   bool pointsindem=dem->GetNearest3Points(*seedlon,*seedlat,pplat,pplon,pphei);

   //Check the dem heights do not contain the flag value for reading outside the limits of the DEM
   if(!pointsindem)
   {
      //Querying points outside the DEM AOI so return NULL.
      return NULL;
   }

   double newp[3]={0},oldp[3]={0};
   //New point is made from the least occuring lat and lon:
   //Old point is from the most occuring lat and lon
   if(pplon[0] == pplon[1])
   {
      newp[0]=pplon[2];
      oldp[0]=pplon[0];
   }
   else if (pplon[0] == pplon[2])
   {
      newp[0]=pplon[1];
      oldp[0]=pplon[0];
   }
   else
   {
      newp[0]=pplon[0];
      oldp[0]=pplon[1];
   }

   if(pplat[0] == pplat[1])
   {
      newp[1]=pplat[2];
      oldp[1]=pplat[0];
   }
   else if (pplat[0] == pplat[2])
   {
      newp[1]=pplat[1];
      oldp[1]=pplat[0];
   }
   else
   {          
      newp[1]=pplat[0];
      oldp[1]=pplat[1];
   }

   //Replace old point with new point
   for(int i=0;i<3;i++)
   {
      if((pplon[i] == oldp[0])&&(pplat[i] == oldp[1]))
      {
         pplon[i]=newp[0];
         pplat[i]=newp[1];
         pphei[i]=dem->GetHeight(pplon[i],pplat[i]);
         if((pphei[i]==DEMOutOfBounds))
         {
            //Point is not within the DEM AOI bounds - return NULL
            //Logger::Log("Trying to access DEM out of bounds in CompleteTheSquare");
            return NULL;
         }
         break;
      }
   }

   //Create a plane in ECEF XYZ for these 3 points - so convert the LLH to XYZ
   double tmpX[3],tmpY[3],tmpZ[3];
   double P1_XYZ[3], P2_XYZ[3], P3_XYZ[3];

   ConvertLLH2XYZ(pplat, pplon, pphei, tmpX, tmpY, tmpZ, 3, GEODETIC,ellipsoid);
   //Point 1 in XYZ
   P1_XYZ[0]=tmpX[0];
   P1_XYZ[1]=tmpY[0];
   P1_XYZ[2]=tmpZ[0];
   //Point 2 in XYZ
   P2_XYZ[0]=tmpX[1];
   P2_XYZ[1]=tmpY[1];
   P2_XYZ[2]=tmpZ[1];
   //Point 3 in XYZ
   P3_XYZ[0]=tmpX[2];
   P3_XYZ[1]=tmpY[2];
   P3_XYZ[2]=tmpZ[2];

   //Create a triangular plane from these points and return it
   TriangularPlane* triplane=new TriangularPlane(P1_XYZ,P2_XYZ,P3_XYZ);
   return triplane;
}

//-------------------------------------------------------------------------
// Get the 3 nearest points from the DEM to the seed position
// Convert these into ECEF XYZ
// Create a planar surface from these 3 points
//-------------------------------------------------------------------------
TriangularPlane* CreatePlaneFromNearestDEMPoints(const double* const seedlat,const double* const seedlon, 
                        Ellipsoid* ellipsoid,DEM* dem)
{
   //3 Points form the plane - use the seed point to start
   double planepointslon[3]={0};
   double planepointslat[3]={0};
   double planepointshei[3]={0};
   bool pointsindem=dem->GetNearest3Points(*seedlon,*seedlat,planepointslat,planepointslon,planepointshei);

   //Check the dem heights do not contain the flag value for reading outside the limits of the DEM - retun NULL if they do
   if(!pointsindem)
   {
      //std::cout<<"View Vector item: "<<pixel<<std::endl;
//      std::cout<<"Point 1: "<< planepointslon[0]<<" "<< planepointslat[0]<<" "<<planepointshei[0]<<std::endl;
//      std::cout<<"Point 2: "<< planepointslon[1]<<" "<< planepointslat[1]<<" "<<planepointshei[1]<<std::endl;
//      std::cout<<"Point 3: "<< planepointslon[2]<<" "<< planepointslat[2]<<" "<<planepointshei[2]<<std::endl;
//      throw "DEM does not cover the entire flight line (Actually - could not find an intersect with the DEM, so could be due to other issues too).\n";
      return NULL;
   }

   //Create a plane in ECEF XYZ for these 3 points - so convert the LLH to XYZ
   double tmpX[3],tmpY[3],tmpZ[3];
   double P1_XYZ[3], P2_XYZ[3], P3_XYZ[3];

   ConvertLLH2XYZ(planepointslat, planepointslon, planepointshei, tmpX, tmpY, tmpZ, 3, GEODETIC,ellipsoid);
   //Point 1 in XYZ
   P1_XYZ[0]=tmpX[0];
   P1_XYZ[1]=tmpY[0];
   P1_XYZ[2]=tmpZ[0];
   //Point 2 in XYZ
   P2_XYZ[0]=tmpX[1];
   P2_XYZ[1]=tmpY[1];
   P2_XYZ[2]=tmpZ[1];
   //Point 3 in XYZ
   P3_XYZ[0]=tmpX[2];
   P3_XYZ[1]=tmpY[2];
   P3_XYZ[2]=tmpZ[2];

   //Create a triangular plane from these points and return
   TriangularPlane* triplane=new TriangularPlane(P1_XYZ,P2_XYZ,P3_XYZ);
   return triplane;
}


//-------------------------------------------------------------------------
// Function to add on a small offset (if required) to shift the seed position
// such that it does not sit on a boundary line (or centre?) of the DEM cells
// i.e. it is defined within a triangle of 3 dem cells
//-------------------------------------------------------------------------
void ShuffleSeed(double* const seedlat,double* const seedlon, DEM* dem)
{
   short XorY=-1;
   //Check if the seed position falls on the bounds of a dem cell
   //and if so offset either the x or the y or both
   if(dem->OnCellBound(*seedlat,*seedlon,&XorY))
   {
      if(XorY==1)
         *seedlon += dem->GetXSpace()/100.0;
      else if(XorY==2)
         *seedlat += dem->GetYSpace()/100.0;
      else if(XorY==3)
      {
         *seedlat += dem->GetYSpace()/100.0;
         *seedlon += dem->GetXSpace()/100.0;
      }
      else
         throw "This should never happen in Shuffle Seed - OnCellBound returned a value other than 1, 2 or 3.";
   }
}

//-------------------------------------------------------------------------
// Function to calculate and set the area of interest for a DEM based
// on the given navigation extents and view vectors.
// Note that the navigation may just be a subsection of a flight line
// and not the entire navigation of the line
//-------------------------------------------------------------------------
bool SetDEMAreaToReadIn(NavBaseClass* nav, ViewVectors* vv, DEM* dem, Ellipsoid* ellipsoid,bool quiet)
{
   //maximum absolute vector for ccd pixels - assumed to be either pixel 0 or pixel n in X vector
   double maxview=vv->AbsMaxX();

   //swath buffer is a predicted maximum half swath width based on maxheight,maxroll + maxview degree offset
   double posswathbuffer=nav->MaxHei() * tan(fabs((nav->MaxRoll()) + maxview)*(PI/180));
   double negswathbuffer=nav->MaxHei() * tan((fabs(nav->MinRoll()) + maxview)*(PI/180));

   //Add the positive and negative buffers together to get a total swath width buffer
   double swathbuffer=posswathbuffer+negswathbuffer;

   //Now convert this distance into degrees - will use the mid latitude / longitude / height of the aircraft to get a better estimate of this conversion
   double tmpx[2],tmpy[2],tmpz[2],tmplat[2],tmplon[2],tmphei[2];
   tmplat[0]=tmplat[1]=nav->MinLat() + 0.5*(nav->MaxLat()-nav->MinLat());
   tmplon[0]=nav->MinLon() + 0.5*(nav->MaxLon()-nav->MinLon());
   tmplon[1]=tmplon[0]+0.1;
   tmphei[0]=tmphei[1]=nav->MinHei() + 0.5*(nav->MaxHei()-nav->MinHei());
   ConvertLLH2XYZ(tmplat, tmplon, tmphei, tmpx, tmpy, tmpz, 2, GEODETIC,ellipsoid);
   double tmpdist=sqrt((tmpx[0]-tmpx[1])*(tmpx[0]-tmpx[1]) + (tmpy[0]-tmpy[1])*(tmpy[0]-tmpy[1]) + (tmpz[0]-tmpz[1])*(tmpz[0]-tmpz[1]) );

   //Swath buffer in degrees is therefore
   double sbdegrees=(swathbuffer / tmpdist)*(tmplon[1]-tmplon[0]);

   if(!quiet)
   {
      Logger::Log("Maximum view angle of sensor (assuming level flying): "+ToString(maxview));
      Logger::Log("At this latitude a distance of: "+ToString(swathbuffer)+" metres is equivalent to "+ToString(sbdegrees)+" degrees.");
      Logger::Log("This will be used as a buffer added onto the Navigation min/max for DEM reading.");
      Logger::Log("Setting DEM area bounds to: min long: "+ToString(nav->MinLon()-sbdegrees) + " min lat: "+ToString(nav->MinLat()-sbdegrees) +
               " max long: " +ToString(nav->MaxLon()+sbdegrees) + " max lat: "+ ToString(nav->MaxLat()+sbdegrees));
   }
   //If the dem does not fully contain these limits then exit now - since the DEM must contain these
   //plus the actual swath on the ground (which will be outside of these limits)
   //we should add on 1 dem cell width extra in each direction to ensure 
   //that there is always a 1 pixel buffer for interpolation incase the view vector hits the final dem cell in the region
   if(!dem->SetAOI(nav->MinLon()-sbdegrees-dem->GetXSpace(),nav->MinLat()-sbdegrees-dem->GetYSpace(),
         nav->MaxLon()+sbdegrees+dem->GetXSpace(),nav->MaxLat()+sbdegrees+dem->GetYSpace()))
   {
      return false; 
   }

   return true;
}
