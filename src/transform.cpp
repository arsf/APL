//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

//Command line to get input and output projection systems
//Then set up the projection systems
//We need to read in the IGM file (assume lat/lon for the moment - read it from header)
//Convert line at a time
//Write out to new IGM file

#include <proj_api.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <algorithm>
#include "logger.h"
#include "commandline.h"
#include "commonfunctions.h"
#include "binfile.h"
#include "bilwriter.h"
#include "basic_igm_worker.h"

const double PI=4*atan(1.0);

//-------------------------------------------------------------------------
// Software description
//-------------------------------------------------------------------------
const std::string DESCRIPTION="Coordinate Transformation Software";

//----------------------------------------------------------------
//Number of options that can be on command line
//----------------------------------------------------------------
const int number_of_possible_options = 7;

//----------------------------------------------------------------
//Option names that can be on command line
//----------------------------------------------------------------
const std::string availableopts[number_of_possible_options]={
"-igm",
"-output",
"-inproj",
"-outproj",
"-inprojstr",
"-outprojstr",
"-help"
}; 

//----------------------------------------------------------------
//One line description of each option - in same order as availableopts
//----------------------------------------------------------------
const std::string optsdescription[number_of_possible_options]={
"3-Band IGM BIL input data filename (Band 1: X, Band 2: Y, Band 3: Z).",
"3-Band IGM BIL output data filename.",
"The projection of the input IGM file (default is WGS84 Geographic Lat/Lon).",
"The projection of the output IGM file.",
"The input projection in the format of a PROJ string.",
"The output projection in the format of a PROJ string.",
"Display this help."
}; 

std::string GetHelpFor(const std::string str);

int main(int argc,char* argv[]) 
{
   std::cout.precision(10);
   //A command line handling object
   CommandLine* cl=NULL;
   //Create a logger to output the data to terminal
   Logger log;

   //PROJ projection objects
   projPJ proj_out, proj_in;
   //Further PROJ projection objects for osng projection (since it needs two projections)
   projPJ proj_out_2nd=NULL, proj_in_2nd=NULL;
   
   std::stringstream strout;
   //Input IGM filename
   std::string strInputIGMFilename;
   //Output IGM filename
   std::string strOutputIGMFilename;

   //Output projection keyword
   std::string strOutProjection;
   //Output ellipsoid to use
   std::string strOutEllipsoid;
   //Input projection keyword
   std::string strInProjection;
   //Input ellipsoid to use
   std::string strInEllipsoid;
   //UTM zone to use for UTM transformations
   std::string strUTMZone;
   //Which hemisphere (north or south)
   std::string strHemisphere;

   //Input projection strings
   std::string projin;
   std::string projin_2nd;
   //Output projection strings
   std::string projout;
   std::string projout_2nd="";

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

      #ifdef DEBUG
      Logger::Log(cl->ReviewCL(availableopts,number_of_possible_options,true)); 
      #endif
      //----------------------------------------------------------------------------
      // Go through each possible option in turn and set up the required data / response
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-help"))
      {
         if(cl->NumArgsOfOpt("-help")==0)
         {
            log.Add(cl->ProgramUsage(number_of_possible_options,availableopts,optsdescription));
            log.Add("\nUse -help <argname> for further information about argument.");
            log.Flush();
         }
         else
            Logger::Log(GetHelpFor(cl->GetArg("-help")));

         throw "";
      }

      Logger::Log("Command line used to run: "+cl->ReturnCLAsString());

      //----------------------------------------------------------------------------
      // Get the IGM file name
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
      // Get the OUTPUT IGM file name
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-output"))
      {
         //Check that an argument follows the output option - and get it if it exists
         if(cl->GetArg("-output").compare(optiononly)!=0)
         {
            strOutputIGMFilename=CreatePath(cl->GetArg("-output"));
            log.Add("Will use output IGM BIL file: "+strOutputIGMFilename);
         }
         else
            throw CommandLine::CommandLineException("Argument -output must immediately precede the output filename.\n");         
      }
      else
      {
         //Throw an exception
         throw CommandLine::CommandLineException("Argument -output [the output IGM BIL file] must be present on the command line.\n");  
      }

      //----------------------------------------------------------------------------
      // Get the input data projection - defaults to WGS84 lat/lon if not specified
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-inproj"))
      {
         //Check that an argument follows the inproj option - and get it if it exists
         if(cl->GetArg("-inproj").compare(optiononly)!=0)
         {
            strInProjection=cl->GetArg("-inproj",0); 
            if(strInProjection.compare("latlong")==0)
            {
               //Lat/lon coordinate system
               strInEllipsoid=cl->GetArg("-inproj",1);
               if((strInEllipsoid.compare("wgs84")==0)||(strInEllipsoid.compare("WGS84")==0))
                  strInEllipsoid="WGS84";                     
               else
                  throw CommandLine::CommandLineException("-inproj currently only supports wgs84 ellipsoid.\nMaybe you should use -inprojstr?");  

               projin="+proj="+strInProjection+" +ellps="+strInEllipsoid;
            }
            
         }
         else
            throw CommandLine::CommandLineException("Argument -inproj must immediately precede the input projection system.\n");         
      }
      else if(cl->OnCommandLine("-inprojstr"))
      {
         //Check that an argument follows the inprojstr option - and get it if it exists
         if(cl->GetArg("-inprojstr").compare(optiononly)!=0)
         {
            projin=cl->GetArg("-inprojstr");
            projin=ReplaceAllWith(&projin,';',' ');
            Logger::Log("Using PROJ formatted projection string: "+projin);
         }
         else
            throw CommandLine::CommandLineException("Argument -inprojstr must immediately precede the input PROJ projection string.\n");         
      }
      else
      {
         //Assume wgs84 lat/lon 
         strInProjection="latlong";
         strInEllipsoid="WGS84";
         projin="+proj="+strInProjection+" +ellps="+strInEllipsoid+" +datum=WGS84 +towgs84=0,0,0";      
      }

      //----------------------------------------------------------------------------
      // Get the output data projection
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-outproj"))
      {
         //Check that an argument follows the outproj option - and get it if it exists
         if(cl->GetArg("-outproj").compare(optiononly)!=0)
         {
            strOutProjection=cl->GetArg("-outproj",0);
            if(strOutProjection.compare("utm_wgs84N")==0)
            {
               //UTM North requested
               strUTMZone=cl->GetArg("-outproj",1);
               if((StringToUINT(strUTMZone)<1)||(StringToUINT(strUTMZone)>60))
                  throw CommandLine::CommandLineException("UTM zone number should be between 1 and 60.\n");         
               strHemisphere="North";
               strOutEllipsoid="WGS84";      
               Logger::Log("Will reproject data into UTM coordinate system: using Zone "+strUTMZone+" North");
               projout="+proj=utm +ellps="+strOutEllipsoid+" +zone="+strUTMZone;   

            }
            else if(strOutProjection.compare("utm_wgs84S")==0)
            {
               //UTM South requested
               strUTMZone=cl->GetArg("-outproj",1);
               if((StringToUINT(strUTMZone)<1)||(StringToUINT(strUTMZone)>60))
                  throw CommandLine::CommandLineException("UTM zone number should be between 1 and 60.\n");         
               strHemisphere="South";
               strOutEllipsoid="WGS84";      
               Logger::Log("Will reproject data into UTM coordinate system: using Zone "+strUTMZone+" South");
               projout="+proj=utm +ellps="+strOutEllipsoid+" +zone="+strUTMZone+" +south";   
            }
            else if(strOutProjection.compare("osng")==0)
            {
//               if((cl->OnCommandLine("-inproj"))||(cl->OnCommandLine("-inprojstr")))
//               {
//                  throw CommandLine::CommandLineException("No input projection allowed for OS National Grid projection. Only accepts WGS84 Geographic Lat/Lon.");
//               }        
//               std::string strGridFile=cl->GetArg("-outproj",1);
//               if(strGridFile.compare("")==0)
//                  throw CommandLine::CommandLineException("Projection Grid filename must follow the 'osng' keyword.");
//               std::ifstream fin;
//               fin.open(strGridFile.c_str());
//               if(!fin.is_open())
//               {
//                 throw CommandLine::CommandLineException("Projection Grid file does not exist or will not open.");
//               }
//               else
//                  fin.close();

//               //Ordnance survey national grid projection
//               //This consists of two projections - first is WGS84 latlon to OSGB36 latlon using a grid
//               // - second is OSGB36 latlon to OSGB36 EN
//               strInEllipsoid="WGS84";
//               strOutEllipsoid="airy";
//               //First projection
//               projin="+proj=latlong +ellps="+strInEllipsoid+" +towgs84=0,0,0";
//               projout="+proj=latlong +nadgrids="+strGridFile;
//               //Second projection
//               projin_2nd="+proj=latlong +datum=OSGB36 +ellps=airy";
//               projout_2nd="+proj=tmerc +datum=OSGB36 +ellps=airy +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000";

               if((cl->OnCommandLine("-inproj"))||(cl->OnCommandLine("-inprojstr")))
               {
                  //Test if it is = latlong WGS84 , if not then exit
                  //std::cout<<cl->GetArg("-inproj")<<std::endl;
                  if((cl->GetArg("-inproj",0).compare("latlong"))&&(cl->GetArg("-inproj",0).compare("WGS84")!=0))
                     throw CommandLine::CommandLineException("No input projection allowed for OS National Grid projection. Only accepts WGS84 Geographic Lat/Lon.");
                  
               } 

               //Test if a gridshift filename was given on command line
               if(cl->GetArg("-outproj",1).compare("")==0)
                  throw CommandLine::CommandLineException("Projection Grid filename must follow the 'osng' keyword.");
       
               //Test if it exists - to deal with spaces we need to make sure we get everything for -outproj (not just 2nd element)
               std::string strGridFileWithProjection=cl->GetArg("-outproj");
               //Find the first ';' - this marks the end of the projection keyword
               size_t pos=strGridFileWithProjection.find_first_of(';');
               //The grid file should be the path given after this ';' - create a substring containing this
               std::string strFullGridFile=GetExistingFilePath(strGridFileWithProjection.substr(pos+1),true);

               //Convert this to absolute path for PROJ
               strFullGridFile=AbsolutePath(strFullGridFile);

               //Check the file exists and can be opened
               std::ifstream fin;
               fin.open(strFullGridFile.c_str());
               if(!fin.is_open())
               {
                 throw CommandLine::CommandLineException("Projection Grid file does not exist or will not open: "
                                                +strFullGridFile+"\nHave you specified the file using an absolute path?");
               }
               else
                  fin.close();

               //Now split up into dirname and filename
               std::string strGridFilePath=DirName(strFullGridFile);
               std::string strGridFile=strFullGridFile.substr(strGridFilePath.length()+1);
               //Assign the grid file path to the search path for proj
               //this allows spaces in path names whilst the 'nadgrids=' keyword does not
               const char** searchpath=new const char*[1];
               searchpath[0]=strGridFilePath.c_str();
               pj_set_searchpath(1, searchpath); 

               //Output this information
               Logger::Log("\nHave set grid shift file search path to: "+strGridFilePath);
               Logger::Log("And using grid shift file name: "+strGridFile);

               //Ordnance survey national grid projection
               //This consists of two projections - first is WGS84 latlon to OSGB36 latlon using a grid
               // - second is OSGB36 latlon to OSGB36 EN
               strInEllipsoid="WGS84";
               strOutEllipsoid="airy";
               //First projection
               projin="+proj=latlong +ellps="+strInEllipsoid+" +towgs84=0,0,0";
               projout="+proj=latlong +ellps="+strOutEllipsoid+" +nadgrids="+strGridFile;
               //Second projection
               projin_2nd="+proj=latlong +ellps="+strOutEllipsoid;
               projout_2nd="+proj=tmerc +ellps="+strOutEllipsoid+" +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000";              

            }
            else
            {
               throw CommandLine::CommandLineException("Unknown output projection. Currently supported: utm_wgs84N, utm_wgs84S, osng.\nMaybe you should use -outprojstr instead?\n");                        
            }
         }
         else
            throw CommandLine::CommandLineException("Argument -outproj must immediately precede the output projection system.\n");         
      }
      else
      {
         //----------------------------------------------------------------------------
         // Get the output data projection PROJ string
         //----------------------------------------------------------------------------
         if(cl->OnCommandLine("-outprojstr"))
         {
            //Check that an argument follows the outprojstr option - and get it if it exists
            if(cl->GetArg("-outprojstr").compare(optiononly)!=0)
            {
               projout=cl->GetArg("-outprojstr");
               projout=ReplaceAllWith(&projout,';',' ');
               Logger::Log("Using PROJ formatted projection string: "+projout);
            }
            else
               throw CommandLine::CommandLineException("Argument -outprojstr must immediately precede the output PROJ projection string.\n");         
         }
         else
         {
            //Throw an exception
            throw CommandLine::CommandLineException("Argument -outproj [the output IGM coordinate projection] must be present on the command line.\n");   
         }
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
      exit(1); //exit the program
   }
   catch(std::exception& e)
   {
      PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString(),&e);
      delete cl;
      exit(1);
   }
   //Flush the log
   log.Flush();

   delete cl;


   Logger::Log("\nPROJ format input projection string: "+projin);
   if (!(proj_in = pj_init_plus(projin.c_str())) )
   {
      Logger::Error("There is a problem with the input projection string:\n"+projin +
                    "\nThe problem was:\n"+std::string(pj_strerrno(*(pj_get_errno_ref()))));
      exit(1);
   }

   //Now do a test by retrieving the string from proj and comparing to the entered string
   Logger::Log("Input projection test returned from proj: "+std::string(pj_get_def(proj_in,0)));   
   if(TrimWhitespace(std::string(pj_get_def(proj_in,0))).compare(projin)!=0)
   {
      Logger::Warning("Input projection entered into and retrieved from PROJ appear to be different.");
   }


   Logger::Log("PROJ format output projection string: "+projout);
   if (!(proj_out = pj_init_plus(projout.c_str())) )
   {
      Logger::Error("There is a problem with the output projection string:\n"+projout +
                    "\nThe problem was:\n"+std::string(pj_strerrno(*(pj_get_errno_ref()))));
      exit(1);
   }

   //Now do a test by retrieving the string from proj and comparing to the entered string
   Logger::Log("Output projection test returned from proj: "+std::string(pj_get_def(proj_out,0)));
   if(TrimWhitespace(std::string(pj_get_def(proj_out,0))).compare(projout)!=0)
   {
      Logger::Warning("Output projection entered into and retrieved from PROJ appear to be different.");
   }

   if(strOutProjection.compare("osng")==0)
   {
      Logger::Log("PROJ format input projection string (part 2): "+projin_2nd);
      if (!(proj_in_2nd = pj_init_plus(projin_2nd.c_str())) )
      {
         Logger::Error("There is a problem with the 2nd input projection string:\n"+projin_2nd +
                    "\nThe problem was:\n"+std::string(pj_strerrno(*(pj_get_errno_ref()))));
         exit(1);
      }

      //Now do a test by retrieving the string from proj and comparing to the entered string
      Logger::Log("Input2 projection test returned from proj: "+std::string(pj_get_def(proj_in_2nd,0)));   
      if(TrimWhitespace(std::string(pj_get_def(proj_in_2nd,0))).compare(projin_2nd)!=0)
      {
         Logger::Warning("Input2 projection entered into and retrieved from PROJ appear to be different.");
      }

      Logger::Log("PROJ format output projection string (part 2): "+projout_2nd);
      if (!(proj_out_2nd = pj_init_plus(projout_2nd.c_str())) )
      {
         Logger::Error("There is a problem with the 2nd output projection string:\n"+projout_2nd +
                    "\nThe problem was:\n"+std::string(pj_strerrno(*(pj_get_errno_ref()))));
         exit(1);
      }
      //Now do a test by retrieving the string from proj and comparing to the entered string
      Logger::Log("Output2 projection test returned from proj: "+std::string(pj_get_def(proj_out_2nd,0)));   
      if(TrimWhitespace(std::string(pj_get_def(proj_out_2nd,0))).compare(projout_2nd)!=0)
      {
         Logger::Warning("Output2 projection entered into and retrieved from PROJ appear to be different.");
      }
   }                     
   else
   {
   }
  
   BinFile* br=NULL;
   try
   {
      //Open a reader to read the data
      br=new BinFile(strInputIGMFilename);
   }
   catch(BinaryReader::BRexception e)
   {
      Logger::Error(std::string(e.what())+"\n"+e.info);
      exit(1); //exit the program
   }

   //Get the number of lines from the header
   unsigned long nlines=StringToUINT(br->FromHeader("lines"));
   //Get the number of samples from the header
   unsigned long nsamps=StringToUINT(br->FromHeader("samples"));
   //Get the number of bands from the header
   unsigned long nbands=StringToUINT(br->FromHeader("bands"));
   if(nbands!=3)
   {
      Logger::Log("Input file has more than 3 bands.");
      exit(1);
   }

   //Get the no data value if in the IGM header file - if there is not one then assign a temporary value
   std::string ignorestr=br->FromHeader("data ignore value");
   double nodatavalue=0;
   if(ignorestr.compare("")==0)
   {
      //There is no "data ignore value" in the header so set our own to be the proj HUGE_VAL
      //And then change to a more suitable value when we know the actual contents of the file.
      nodatavalue=HUGE_VAL;
   }
   else
   {
      //Use the no data value given in file (try/catch in case its dodgy)
      try
      {
         nodatavalue=StringToDouble(ignorestr);
         Logger::Log("Will use a 'no data value' of: "+ToString(nodatavalue));
         //Test for the edge case where the ignore value is the HUGE_VAL value - we use this value when there is no 
         //ignore value in the header file. So better to change the ignore value to something else in the data set and re-run?
         if(nodatavalue==HUGE_VAL)
         {
            Logger::Warning("The data ignore value can be any value except "+ToString(HUGE_VAL)+". Please use a different ignore value in your data and re-run.");
            exit(1);
         }
      }
      catch(std::string e)
      {
         Logger::Error(e);
         exit(1);
      }
   }

   //Create arrays to hold the data
   double* X=new double[nsamps];
   double* Y=new double[nsamps];
   double* Z=new double[nsamps];

   double tmaxx=0,tminx=0,tmaxy=0,tminy=0;
   double maxx=-99999999, minx=99999999,maxy=-999999999,miny=99999999;

   //Create an output writer
   BILWriter* bw=NULL;
   try
   {
      bw=new BILWriter(strOutputIGMFilename,FileWriter::float64,nlines,nsamps,3,'w');
   }
   catch(BILWriter::BILexception e)
   {
      Logger::Error(std::string(e.what())+"\n"+e.info);
      exit(1); //exit the program
   }

   Logger::Log("\nPlease note that Z values are not transformed and will remain in the input reference.");

   for(unsigned int line=0;line<nlines;line++)
   {
      //std::cout<<"Converting line:"<<line<<std::endl;
      
      br->Readbandline((char*)X,0,line);
      br->Readbandline((char*)Y,1,line);
      br->Readbandline((char*)Z,2,line);

      //Check for no data value and set to HUGE_VAL if there are any
      for(unsigned int s=0;s<nsamps;s++)
      {
         if((X[s]==nodatavalue)||(Y[s]==nodatavalue))
         {
            X[s]=HUGE_VAL;
            Y[s]=HUGE_VAL;
         }
      }     

      //std::cout<<"Changing to radians:"<<line<<std::endl;
      if(pj_is_latlong( proj_in))
      {
         for(unsigned int s=0;s<nsamps;s++)
         {
            if((X[s]==nodatavalue)||(Y[s]==nodatavalue))
               continue; //skip converting these as they're no data value
            else
            {
               X[s]=X[s]*PI/180.0;
               Y[s]=Y[s]*PI/180.0;
               Z[s]=Z[s];
            }
         }
      }

      int ret=pj_transform(proj_in,proj_out,nsamps,1,X,Y,NULL);
      if(ret!=0)
      {
         Logger::Error("Error in transformation: " + std::string(pj_strerrno(ret)));
         exit(1);  // FIXME: lame
      }

      //If OSNG then do second transformation
      if(strOutProjection.compare("osng")==0)
      {       
         ret=pj_transform(proj_in_2nd,proj_out_2nd,nsamps,1,X,Y,NULL);
         if(ret!=0)
         {
            Logger::Error("Error in osng 2nd transformation: " + std::string(pj_strerrno(ret)));
            exit(1);  // FIXME: lame
         }
      }

      //Test projection is working / suitable - e.g. if a usable utm zone is given
      for(unsigned int s=0;s<nsamps;s++)
      {
         if((X[s]==HUGE_VAL)||(Y[s]==HUGE_VAL))
         {
            //If nodatavalue is HUGE_VAL then there was no ignore data in hdr file - something gone wrong with transformation
            if(nodatavalue==HUGE_VAL)
            {
               Logger::Error("Error in transformation - maybe selected projection is unsuitable for data. ");               
               exit(1);
            }
            else
            {
               //This is probably just the ignore value - but warn just to make user aware something could be incorrect
               Logger::WarnOnce("Possible error in transformation - probably due to NO DATA VALUE existing in IGM file - but could be incorrect projection for data.");
            }
         }
      }

      //If projection is still in lat/long we need to convert from radians to degrees before output
      //NEED A BETTER TEST THAN THIS
      if((pj_is_latlong(proj_out)==true)&&(proj_out_2nd==NULL))
      {
         //std::cout<<"pj_is_latlong:"<<pj_is_latlong(proj_out)<<"proj-out2:"<<((proj_out_2nd))<<" projout_2nd:"<<projout_2nd<<"."<<std::endl;
         for(unsigned int s=0;s<nsamps;s++)
         {
            X[s]=X[s]*180.0/PI;
            Y[s]=Y[s]*180.0/PI;
            Z[s]=Z[s];
         }   
      }

      //Get min/max x for this line - ignoring HUGE_VAL
      GetArrayLimits(X,nsamps,tminx,tmaxx,HUGE_VAL);
      //Update variables if a new min or max is found
      if(tmaxx > maxx)
         maxx=tmaxx;
      if(tminx < minx)
         minx=tminx;

      //Get min/max y for this line - ignoring HUGE_VAL
      GetArrayLimits(Y,nsamps,tminy,tmaxy,HUGE_VAL);
      //Update variables if a new min or max is found
      if(tmaxy > maxy)
         maxy=tmaxy;
      if(tminy < miny)
         miny=tminy;

      //Convert HUGE_VAL back into NODATAVALUE
      //Note if we assigned nodatavlue=HUGE_VAL then this SHOULD not matter as
      //there was no ignore value - so there should be no HUGE_VAL being written out
      //unless something went wrong with the transformation.
      for(unsigned int s=0;s<nsamps;s++)
      {
         if((X[s]==HUGE_VAL)||(Y[s]==HUGE_VAL))
         {
            X[s]=Y[s]=nodatavalue;
         }
      } 

      bw->WriteBandLine((char*)X);
      bw->WriteBandLine((char*)Y);
      bw->WriteBandLine((char*)Z);

      //Percent done counter
      PercentProgress(line,nlines);

   }

   //Add the projection information to the IGM file
   if((pj_is_latlong(proj_out))&&(proj_out_2nd==NULL))
   {
      bw->AddToHdr("projection = Geographic Lat/Lon");
      std::string projdef=std::string(pj_get_def(proj_out,0));
      std::string ellipsesearchterm="+ellps=";
      std::string projdefafterellps=projdef.substr(projdef.find(ellipsesearchterm)+ellipsesearchterm.length());
      std::string strEllps=GetItemFromString(projdefafterellps,0);
      bw->AddToHdr("datum ellipsoid = "+strEllps);
   }
   else
   {
      bw->AddToHdr("projection = "+strOutProjection+" "+strUTMZone+" "+strHemisphere);
      bw->AddToHdr("datum ellipsoid = "+strOutEllipsoid);
   }

   if(proj_out_2nd!=NULL)
   {
      bw->AddToHdr("proj4 projection string 1 = "+std::string(pj_get_def(proj_out,0)));
      bw->AddToHdr("proj4 projection string 2 = "+std::string(pj_get_def(proj_out_2nd,0)));
      //Deallocate grid 
      pj_deallocate_grids();
   }
   else
      bw->AddToHdr("proj4 projection string = "+std::string(pj_get_def(proj_out,0)));

   //Copy the xtsart and ystart from the original IGM file
   std::string xstart=br->FromHeader("x start");
   std::string ystart=br->FromHeader("y start");
   bw->AddToHdr(";These describe which pixels from the original raw image the IGM file positions relate to.");
   bw->AddToHdr("x start = "+xstart);   
   bw->AddToHdr("y start = "+ystart);   

   //Need to add min/max bounds to the IGM
   bw->AddToHdr(";Min X = "+ToString(minx));
   bw->AddToHdr(";Max X = "+ToString(maxx));
   bw->AddToHdr(";Min Y = "+ToString(miny));
   bw->AddToHdr(";Max Y = "+ToString(maxy));
   //Label the bands (vague)
   bw->AddToHdr("band names = {X,Y,Height}");

   //Copy the no data value over to the new IGM (or add it if it was missing from original)
   //If the value is HUGE_VAL then we assume that there was not a value in the hdr file [possibly could have been HUGE_VAL].
   if(nodatavalue==HUGE_VAL)
   {
      //We want to select a value that is not in the file and can be represented easily in
      //the header file such that it can be converted from string to double with no loss.
      nodatavalue= floor( -fabs(std::min(minx,miny)) -1);
   }
   bw->AddToHdr(";In most cases there are no data with the 'data ignore' value. However it is always included in the hdr for consistency.");
   bw->AddToHdr("data ignore value = "+ToString(nodatavalue));
   br->Close();
   bw->Close();

   delete br;
   delete bw;

   delete[] X;
   delete[] Y;
   delete[] Z;

   pj_free(proj_in);
   pj_free(proj_out);

   //Could output some pixel size information here that could be used for the later resampling stages
   //Open the file that has been created
   Basic_IGM_Worker igm(strOutputIGMFilename);
   //Get the pixel size information for the centre pixel
   double pixsize[7]={0};
   igm.GetPixelSize(igm.Samples()/2,pixsize);
   Logger::Log("\nAverage nadir pixel sizes in along track, across track are: "+ToString(pixsize[0])+" "+ToString(pixsize[1]));
   Logger::Log("Average nadir pixel sizes in projected X,Y are: "+ToString(pixsize[3])+" "+ToString(pixsize[6]));
   Logger::Log("Coordinate Transformation Complete");
}

//--------------------------------------------------------------------------
//Function to return a more detailed help string about the available options
//--------------------------------------------------------------------------
std::string GetHelpFor(const std::string str)
{
   std::map <std::string,std::string> helpdoc;

   helpdoc["igm"]="The input IGM data to apply the coordinate transform to. This should be a 3-band BIL (binary) data file.\n"
                  "Band 1 is the X data (e.g. longitude)\n"
                  "Band 2 is the Y data (e.g. latitude)\n"
                  "Band 3 is the Z data (e.g. height). IMPORTANT This band is disregarded in the transformation.\n";

   helpdoc["output"]="The output IGM data in the new coordinate system. This will be a 3-band BIL (binary) data file.\n"
                  "Band 1 is the X data (e.g. easting)\n"
                  "Band 2 is the Y data (e.g. northing)\n"
                  "Band 3 is the Z data (e.g. height). IMPORTANT This band has been disregarded in the transformation.\n";

   helpdoc["inproj"]="To easily select a common projection to describe the input IGM file. \nOptions:\n"
                     "   latlong <ellipsoid> - Data is in latitude and longitude referenced to ellipsoid <ellipsoid>.\n"
                     "\nDefault if the -inproj option is missing from command line is to use latlong WGS84.\n";

   helpdoc["outproj"]="To easily select a common projection for the output IGM file. Uses PROJ to reproject the data. \nOptions:\n"
                      "   utm_wgs84N <zone> - Output to UTM North projection using the WGS84 ellipsoid, for zone <zone>.\n"
                      "   utm_wgs84S <zone> - Output to UTM South projection using the WGS84 ellipsoid, for zone <zone>.\n"
                      "   osng <gridfile> - Output to Ordnance Survey National Grid (OSGB36/OSTN02) projection, using the gridfile to apply the transformation.\n";


   if(helpdoc[str].compare("")==0)
      helpdoc[str]="No extra help for this topic yet.";

   return helpdoc[str];
}

