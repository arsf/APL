//-------------------------------------------------------------------------
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK
//
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0
//
//You should have received a copy of the Licence along with the APL source
//If not, please contact arsf-processing@pml.ac.uk
//-------------------------------------------------------------------------

#include "binfile.h"
#include "bilwriter.h"
#include "commandline.h"
#include "logger.h"
#include <vector>
#include <string>
#include <bitset>
#include <cmath>

//-------------------------------------------------------------------------
// Software description
//-------------------------------------------------------------------------
const std::string DESCRIPTION="Apply an APL mask to Level-1 data";

const std::string strProgramBlurb="\nThis utility will scan through the given mask file and apply the masked data value (0) "
                                  "to the given level-1 data - outputting to a new file. The pixels that will be masked are "
                                  "determined by the flag values given. By default, all pixels that are not marked good in the mask file "
                                  "(have a value other than 0) will be masked.\n";


//----------------------------------------------------------------
//Number of options that can be on command line
//----------------------------------------------------------------
const int number_of_possible_options = 6;

//----------------------------------------------------------------
//Option names that can be on command line
//----------------------------------------------------------------
const std::string availableopts[number_of_possible_options]={
"-lev1",
"-mask",
"-output",
"-flags",
"-onlymaskmethods",
"-help"
}; 

//----------------------------------------------------------------
//One line description of each option - in same order as availableopts
//----------------------------------------------------------------
const std::string optsdescription[number_of_possible_options]={
"Level-1 data file containing the data to map.",
"Mask data file relating to level 1 image (this currently must be an 8-bit char file).",
"The filename of the masked (output) level-1 BIL file.",
"List of space separated flag values to use as mask. Select from: 1 2 4 8 16 32 64. Default is to mask all types.",
"Only mask bad CCD pixels (i.e. with value 4 in mask file) detected by certain methods. Select from: A=1 B=2 C=4 D=8 E=16",
"Display this help."
};

//-------------------------------------------------------------------------
// Template function to apply mask to data file
//-------------------------------------------------------------------------
template <class L,class M>
void MaskData(L* lev1data,M* maskdata,BinFile* lev1,BinFile* mask, BILWriter* output, unsigned int flagsum,L MASKVALUE,BinFile* maskmethod,std::vector<unsigned char> methods)
{
   unsigned int lines=StringToUINT(lev1->FromHeader("lines"));
   unsigned int bands=StringToUINT(lev1->FromHeader("bands"));
   unsigned int samples=StringToUINT(lev1->FromHeader("samples"));

   //Create arrays to hold 1 scan line of data
   lev1data=new L[bands*samples];
   maskdata=new M[bands*samples];
   char* methoddata=NULL;

   //Check the maskmethod stuff
   if(((maskmethod!=NULL)&&(methods.size()==0))||((maskmethod==NULL)&&(methods.size()!=0)))
   {
      throw "MaskData Function: Mask method file given but no methods (or vice versa). Either specify both or neither.";
   }  
   else if((maskmethod!=NULL)&&(methods.size()!=0))
   {
      methoddata=new char[bands*samples];
   }

   //Set the bits to test against (sizeof in bytes not bits so multiply by 8)
   std::bitset<sizeof(M)*8> flagtester(flagsum);
   int* bit=new int[flagtester.count()];
   int position=0;
   //std::cout<<"flagsum "<<flagsum<<" count "<<flagtester.count()<<std::endl;
   //Get the positions of 'true' bits and store in bits array
   for (unsigned int i=0; i<flagtester.size(); i++)
   {
      if(flagtester.test(i)==true)
      {
         bit[position]=i;
         position++;
         Logger::Log("Will test against mask bit position: "+ToString(i)+" equivalent to value: "+ToString(pow(2,i)));
      }
   }

   std::vector<unsigned int> droppedscan;
   bool SKIPREMAININGFLAGS=false;
   //Loop through each line of file
   for(unsigned int l=0;l<lines;l++)
   {
      //Read in a line from level1 
      lev1->Readline((char*)lev1data,l);
      //Read in a line from mask
      mask->Readline((char*)maskdata,l);
      //Read a line from the mask method file if given
      if(maskmethod!=NULL)
         maskmethod->Readline((char*)methoddata,l); 

      //Loop through each band and sample and set level1 data to flag value
      for(unsigned int b=0;b<bands;b++)
      {
         for(unsigned int s=0;s<samples;s++)
         {
            std::bitset<sizeof(M)*8> maskbits(maskdata[b*samples + s]);
            for(unsigned int pos=0;pos<flagtester.count();pos++)
            {
               //Test if any of the true bits in bits array are true in maskbits
               //if so then mask the level 1 data
               if(maskbits[bit[pos]]==true)
               {
                  if((methoddata!=NULL)&&(bit[pos]==2))
                  {
                     //This means we need to use the method data in the mask method file
                     //Loop through the methods vector and test until we find a match
                     for(std::vector<unsigned char>::iterator it=methods.begin();it<methods.end();it++)
                     {
                        //bitwise AND should detect if methoddata and method are both set to true
                        if((methoddata[b*samples + s]&(*it)) > 0)
                        {
                           //std::cout<<"Masking: "<<b<<" "<<s<<" "<<(unsigned int)methoddata[b*samples + s]<<std::endl;
                           lev1data[b*samples + s] = MASKVALUE;
                           SKIPREMAININGFLAGS=true;
                           break;
                        }
                     }
                  }
                  else
                  {
                     lev1data[b*samples + s] = MASKVALUE;
                     SKIPREMAININGFLAGS=true;
                  }

                  //As the data has been masked for this flag we can skip the rest of the flags
                  if(SKIPREMAININGFLAGS==true)
                  {
                     SKIPREMAININGFLAGS=false;
                     break;
                  }
               }
            }
         }
      }

      //Add a test for dropped scans - this uses the 1st sample of the 1st band to check
      std::bitset<sizeof(M)*8> maskbits(maskdata[0]);
      if(maskbits[4]==true)
      {
         droppedscan.push_back(l);
      }

      //Write out this line of data
      output->WriteLine((char*)lev1data);

      //Percent done counter
      PercentProgress(l,lines);

   }

   //Output the dropped scans that have been identified (referenced to 0)
   if(droppedscan.size()!=0)
   {
      Logger::Log("Dropped scans: ");
      for(std::vector<unsigned int>::iterator it=droppedscan.begin();it<droppedscan.end();it++)
         std::cout<<(*it)<<" ";
      Logger::Log("");
   }
   if(lev1data != NULL)
   {
      delete[] lev1data;
      lev1data=NULL;
   }

   if(maskdata != NULL)
   {
      delete[] maskdata;   
      maskdata=NULL;
   }
   
   if(methoddata != NULL)
   {
      delete[] methoddata;
      methoddata=NULL;
   }

   if(bit!=NULL)
      delete[] bit;
}



int main (int argc, char* argv[])
{
   //A command line handling object
   CommandLine* cl=NULL;
   //Create a logger to output the data to terminal
   Logger log(0);

   int status=0;
   std::string strLev1File="";
   std::string strOutputName="";
   std::string strMaskFile="";
   std::string strMaskMethodFile="";
   std::vector<unsigned char>methods;

   unsigned int flagsum=0;

   //Get exe name without the path
   std::string niceexename=std::string(argv[0]);
   niceexename=niceexename.substr(niceexename.find_last_of("/\\")+1);

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

            Logger::Log(strProgramBlurb);

            Logger::Log(cl->ProgramUsage(number_of_possible_options,availableopts,optsdescription));
            //log.Flush();
         }
         throw "";
      }

      //----------------------------------------------------------------------------
      // Get the level 1 file to apply mask to
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
         throw CommandLine::CommandLineException("Argument -lev1 [the Level-1 data BIL file] must be present on the command line.\n");  
      }

      //----------------------------------------------------------------------------
      // Get the name of the output masked level1 file
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-output"))
      {
         if(cl->GetArg("-output").compare(optiononly)!=0)
         {
            strOutputName=CreatePath(cl->GetArg("-output"));
            Logger::Log("Will write masked data to: "+strOutputName);
         }
         else
            throw CommandLine::CommandLineException("Argument -output must immediately precede the name of the output masked file.\n");         
      }
      else
      {
         throw CommandLine::CommandLineException("Argument -output must be present on the command line.\n");  
      }

      //----------------------------------------------------------------------------
      // Get the name of the mask file
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-mask"))
      {
         if(cl->GetArg("-mask").compare(optiononly)!=0)
         {
            strMaskFile=GetExistingFilePath(cl->GetArg("-mask"),true);
            Logger::Log("Will use input mask BIL file: "+strMaskFile);
         }
         else
            throw CommandLine::CommandLineException("Argument -mask must immediately precede the mask filename.\n");         
      }
      else
      {
         throw CommandLine::CommandLineException("Argument -mask [the Level-1 mask BIL file] must be present on the command line.\n");  
      }
      //----------------------------------------------------------------------------
      // Get the list of flags to mask out
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-flags"))
      {
         if(cl->GetArg("-flags").compare(optiononly)!=0)
         {
            std::string strFlags=cl->GetArg("-flags");
            Logger::Log("Will mask only flag values of: "+strFlags);
            int value=0,position=0;
            //Loop through the list to get the flag values
            while(GetItemFromString(strFlags,position,';').compare("")!=0)
            {
               std::string aflag=GetItemFromString(strFlags,position,';');
               value=StringToINT(aflag);
               //Check flag values
               if(value==0)
                  Logger::Warning("Will mask good values (value=0) are you sure you want to do this?");
               switch(value)
               {
               case 1:
               case 2:
               case 4:
               case 8:
               case 16:
               case 32:
               case 64:
                  break;
               default:
                  throw "Unknown flag value: "+ToString(value);
               }
               flagsum = flagsum + value;
               position++;
            }
         }
         else
            throw CommandLine::CommandLineException("Argument -flags must immediately precede the list of flag values to use.\n");         
      }
      else
      {
         Logger::Log("Will flag all pixels which are not good (have non-zero value in mask file).");  
         flagsum=1+2+4+8+16+32+64;
      }
      //----------------------------------------------------------------------------
      // Get the name of the mask method file if given and which methods to mask
      //----------------------------------------------------------------------------
      if(cl->OnCommandLine("-onlymaskmethods"))
      {
         if(cl->GetArg("-onlymaskmethods").compare(optiononly)!=0)
         {
            strMaskMethodFile=GetExistingFilePath(cl->GetArg("-onlymaskmethods",0),true);
            Logger::Log("Will use input mask method BIL file: "+strMaskMethodFile);
            //Now get methods to match against and put into vector 'methods'
            for(int i=1;i<cl->NumArgsOfOpt("-onlymaskmethods");i++)
            {
               std::string m=cl->GetArg("-onlymaskmethods",i);
               if(m.size() > 1)
                  throw "Unknown mask bad CCD pixel detection method: "+m;
               else    
                  methods.push_back(*(m.c_str()));
            }

            //Now check these are sensible values and convert from the ascii char value to the "bit" value
            //eg A=1, B=2, C=4, D=8 ...
            Logger::Log("Only masking bad CCD pixels that were detected by method:");
            for(std::vector<unsigned char>::iterator it=methods.begin();it<methods.end();it++)
            {
               //Output the value we're checking against
               std::cout<<" "<<(*it)<<std::endl;
               //Check for appropriate values and then convert to 'bit' value
               switch((*it))
               {
               case 'A':
                  (*it)=static_cast<unsigned char>(1); //00000001
                  break;
               case 'B':
                  (*it)=static_cast<unsigned char>(2); //00000010
                  break;
               case 'C':
                  (*it)=static_cast<unsigned char>(4); //00000100
                  break;
               case 'D':
                  (*it)=static_cast<unsigned char>(8); //00001000
                  break;
               case 'E':
                  (*it)=static_cast<unsigned char>(16);//00010000
                  break;
               default:
                  throw "Unknown mask bad CCD detection method value. Please review the -onlymaskmethods option.";
               }
            }
         }
         else
            throw CommandLine::CommandLineException("Argument -onlymaskmethods must immediately precede the mask method filename and list of methods to match against.\n");         
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
   catch(...)
   {
      PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString());
      delete cl;
      exit(1);
   } 

   char* clev1data=NULL;
   short int* silev1data=NULL;
   unsigned short int* usilev1data=NULL;
   int* ilev1data=NULL;
   unsigned int* uilev1data=NULL;
   float* flev1data=NULL;
   double* dlev1data=NULL;

   //THIS NEEDS TO BE CHANGED TO BE WHATEVER THE MASK FILE IS
   //CURRENTLY WILL ONLY SUPPORT CHAR MASK FILES
   char* maskdata=NULL;

   try
   {
      //Open level 1 BIL file and get dimensions
      BinFile lev1(strLev1File);
      unsigned int lines=StringToUINT(lev1.FromHeader("lines"));
      unsigned int bands=StringToUINT(lev1.FromHeader("bands"));
      unsigned int samples=StringToUINT(lev1.FromHeader("samples"));

      //Open Mask BIL file and get dimensions
      BinFile mask(strMaskFile);
      unsigned int masklines=StringToUINT(mask.FromHeader("lines"));
      unsigned int maskbands=StringToUINT(mask.FromHeader("bands"));
      unsigned int masksamples=StringToUINT(mask.FromHeader("samples"));

      BinFile* maskmethod=NULL;
      //If given, open the Mask Method BIL file
      if(strMaskMethodFile!="")
      {
         maskmethod=new BinFile(strMaskMethodFile);
         unsigned int maskmethodlines=StringToUINT(mask.FromHeader("lines"));
         unsigned int maskmethodbands=StringToUINT(mask.FromHeader("bands"));
         unsigned int maskmethodsamples=StringToUINT(mask.FromHeader("samples"));

         if((masklines != maskmethodlines) || (maskbands!=maskmethodbands) || (masksamples!=maskmethodsamples))
            throw "Mask dimensions do not match mask methods file dimensions. Are you sure these mask and method files are a pair?";
      }

      if((masklines != lines) || (maskbands!=bands) || (masksamples!=samples))
         throw "Mask dimensions do not match level 1 file dimensions. Are you sure this mask is for this level1 file?";

      if(mask.GetDataType()!=1)
         throw "Can only support mask files containing 1-byte data type at the moment.";

      //Get the output data type - this will match the input level 1 data type
      BILWriter::DataType output_data_type;
      switch(lev1.GetDataType())
      {
      case 1: //8-bit
         {
            output_data_type=BILWriter::uchar8;
            break;   
         }
      case 2: //16 bit signed short int
         {
            output_data_type=BILWriter::int16;
            break;
         }
      case 3: //32 bit signed short int
         {
            output_data_type=BILWriter::int32;
            break;
         }
      case 4: //float
         {
            output_data_type=BILWriter::float32;   
            break;
         }
      case 5: //double
         {
            output_data_type=BILWriter::float64;
            break;
         }
      case 12: //16 bit unsigned short int
         {
            output_data_type=BILWriter::uint16;
            break;
         }
      case 13: //32 bit unsigned int
         {
            output_data_type=BILWriter::uint32;
            break;
         }
      default:
         throw "Unrecognised data type for level-1 BIL file. Currently supports 8-bit, signed and unsigned 16-bit and 32-bit integer, 32-bit and 64-bit float";
         break;
      }

      //Open the file to output to
      BILWriter output(strOutputName,output_data_type,lines,samples,bands,'w');

      //Copy over the level-1 header contents to the new level-1 file
      std::map<std::string, std::string, cmpstr> header=lev1.CopyHeaderExcluding();
      for(std::map<std::string, std::string, cmpstr>::iterator iter=header.begin();iter!=header.end();iter++)
      {
         std::string tmp="";
         if((*iter).first.at(0)!=';')
            tmp=(*iter).first+" = "+(*iter).second;
         else
            tmp=(*iter).first; //a comment - does not have a second

         output.AddToHdr(lev1.TidyForHeader(tmp));
      }
      
      switch(lev1.GetDataType())
      {
      case 1: //8-bit
         {
            char MASKVALUE=0;
            MaskData(clev1data,maskdata,&lev1,&mask,&output,flagsum,MASKVALUE,maskmethod,methods);
            break;   
         }
      case 2: //16 bit signed short int
         {
            short int MASKVALUE=0;
            MaskData(silev1data,maskdata,&lev1,&mask,&output,flagsum,MASKVALUE,maskmethod,methods);
            break;
         }
      case 3: //32 bit signed int
         {
            int MASKVALUE=0;
            MaskData(ilev1data,maskdata,&lev1,&mask,&output,flagsum,MASKVALUE,maskmethod,methods);
            break;
         }
      case 4: //float
         {
            float MASKVALUE=0;
            MaskData(flev1data,maskdata,&lev1,&mask,&output,flagsum,MASKVALUE,maskmethod,methods);   
            break;
         }
      case 5: //double
         {
            double  MASKVALUE=0;
            MaskData(dlev1data,maskdata,&lev1,&mask,&output,flagsum,MASKVALUE,maskmethod,methods);
            break;
         }
      case 12: //16 bit unsigned short int
         {
            unsigned short int MASKVALUE=0;
            MaskData(usilev1data,maskdata,&lev1,&mask,&output,flagsum,MASKVALUE,maskmethod,methods);
            break;
         }
      case 13: //32 bit unsigned int
         {
            unsigned int MASKVALUE=0;
            MaskData(uilev1data,maskdata,&lev1,&mask,&output,flagsum,MASKVALUE,maskmethod,methods);
            break;
         }
      default:
         throw "Unrecognised data type for level-1 BIL file. Currently supports 8-bit, signed and unsigned 16-bit and 32-bit integer, 32-bit and 64-bit float";
         break;
      }

      if(maskmethod!=NULL)
         delete maskmethod;
   }
   catch(std::string e)
   {
      Logger::Error(e);
      status=1;
   }
   catch(const char* e)
   {
      Logger::Error(e);
      status=1;
   }
   catch(BILWriter::BILexception e)
   {
      Logger::Error(e.info);
      status=1;
   }
   catch(BinaryReader::BRexception e)
   {
      Logger::Error(e.info);
      status=1;
   }
   catch(std::exception& e)
   {
      PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString(),&e);
      status=1;
   }
   catch(...)
   {
      PrintAbnormalExitMessage(__FILE__,__LINE__,niceexename,VERSION,CONTACTEMAIL,cl->ReturnCLAsString());
      status=1;
   } 

   if(maskdata!=NULL)
      delete[] maskdata;

   return status;
}
