//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#include "binaryreader.h"

#ifndef BILREADERDEBUG
   #define DEBUGPRINT(X)  
#else
   #define DEBUGPRINT(X) std::cout<< X <<std::endl;
#endif


BinaryReader::BinaryReader() //constructor
{
   //Assign error string and set items to NULL that are undefined
   //may need to do a more full setting of items to 0 here.
   missingheaderitemerror="Trying to read a key from hdr which does not exist: ";
   this->filein=NULL;

}

BinaryReader::BinaryReader(std::string fname)//constructor with filename - only to return the filetype (BIL or BSQ)
{
   //Assign error string
   missingheaderitemerror="Trying to read a key from hdr which does not exist: ";
   //Asssign filename
   this->filename=fname;
   //Want to open bsq file, read in hdr file, and assign all the class variables
   this->filein = fopen(this->filename.c_str(),"rb");
   if(!this->IsOpen())
   {
      //File has failed to open
      brinfo<<"An error has occured opening the file: "<<this->filename<<std::endl;
      throw BRexception(brinfo.str()); //Throw an exception
   }
   else
   {
      //File exists and has been opened successfully
      //Now try and read in hdr file
      int retval=this->ReadHeader();
      //An error occured, store a general error message into the brinfo stream
      if(retval!=1)
      {
         brinfo<<"An error has occurred whilst trying to read from the hdr file for data file: "<<this->filename<<std::endl;      
         throw BRexception(brinfo.str()); //Throw an exception
      }
      std::map<std::string,std::string>::iterator it; //use an iterator to access the Header map
      //set the numsamps, number of lines etc and also perform other checks against file
      it=this->Header.find("interleave"); //Search for interleave type
      if(it==this->Header.end())//check if present in header
      {
         this->CannotFind("interleave");
      }
      if(it->second == "bsq")
      {
         FILESTYLE=BSQ;
      }
      else if(it->second == "bil")
      {
         FILESTYLE=BIL;
      }
      else
      {
         FILESTYLE=UNKNOWN;
      }
   }

   //Destroy the header here so as not to tempt anyone to use this to read data
   //It is purely to get the file stye (bil/bsq) only
   Header.clear();

}

BinaryReader::~BinaryReader()//destructor
{
   Close();
}

std::string BinaryReader::GetHeaderFileName(const short type)
{
   std::string hdrfilename=this->filename;
   switch(type)
   {
   case 0:
      //Try trimming off the file extension and replacing with hdr
      hdrfilename.replace(this->filename.length() - 3, 3, "hdr");
      break;
   case 1:
      //Try appending on .hdr to filename
      hdrfilename.append(".hdr");
      break;
   default:
      throw "Unknown type given to GetHeaderFileName.";
      break;
   }
   return hdrfilename;
}

//-------------------------------------------------------------------------
// Function to read in the information from the header file and store in map
// Keys are now converted to lowercase , Values are not converted
//-------------------------------------------------------------------------
int BinaryReader::ReadHeader()
{
   //DEBUG statement
   DEBUGPRINT("Reading BIL Header file...");
   //Open the header file 
   //It is expected to be the same as filename(-'raw' + 'hdr')
   //std::string hdrfile=this->filename;
   //hdrfile.replace(this->filename.length() - 3, 3, "hdr");
   std::string hdrfile;
   //Now try and open the file
   std::ifstream hdrin;
   int failcount=0,maxtypes=2;

   //Loop through the types of header file name to try and open until one works
   for(short int type=0;type<maxtypes;type++)
   {
      this->hdrfilename=GetHeaderFileName(type);
      hdrin.open(this->hdrfilename.c_str());
      if(!hdrin.is_open())
      {
         failcount++;
         //An error has occured
         //brinfo<<"An error has occurred opening the hdr file: "<<this->hdrfilename<<std::endl;
         continue;
         //return -1;
      }
      else
         break;
   }

   if(failcount==maxtypes)
   {
      //An error has occured
      brinfo<<"Unable to open a hdr file, does one exist?"<<std::endl;
      return -1;
   }
   else
   {
      //the hdr file has been opened successfully
      //Now search for the number of lines/samples/bands first
      std::string strbuffer;//buffer to hold lines from hdr file
      std::string strkey; //string for storing map key temporarally
      std::string strval; //string for storing map value temporarly
      bool multilineval=false; // bool to test if values are over multiple file lines
      while(!hdrin.eof())
      {
         //get a line from the hdr file 
         std::getline(hdrin,strbuffer);
         //std::cout<<strbuffer<<std::endl;
         //search for an equals
         if(strbuffer.find('=') != std::string::npos)
         {
            //An = has been found in this string
            if(strbuffer.find('{') != std::string::npos)
            {
               //A { has been found. This could mean that the value is on more than one line 
               strkey=strbuffer.substr(0,strbuffer.find('=')); //can set key
               strval.assign(strbuffer.substr(strbuffer.find('=')+1)); //this is part of value
               //trim whitespace
               strkey=TrimWhitespace(strkey);
               strval=TrimWhitespace(strval);

               if(strbuffer.find('}')!=std::string::npos) // } also found on same line
               {
                  //The value is only on one line
                  this->Header[ToLowerCase(strkey)]=strval;
                  strkey.clear();
                  strval.clear();
                  continue;
                  
               }
               else
                  multilineval=true;// value is spread over multiple lines of header file
            }
            else
            {
               //The value is only on this line
               //Split up string at equals
               strkey=strbuffer.substr(0,strbuffer.find('='));
               strval=strbuffer.substr(strbuffer.find('=')+1);
               //trim whitespace
               strkey=TrimWhitespace(strkey);
               strval=TrimWhitespace(strval);
               this->Header[ToLowerCase(strkey)]=strval;
               strkey.clear();
               strval.clear();//empty the string
               continue;
            }
         }         
         else if(strbuffer.find('}') != std::string::npos)
         {
            //there is a }. This could mean that the end of a value has been reached
            //Check if found at the end of the line (whitespace removed)
            strbuffer=TrimWhitespace(strbuffer);
            if(strbuffer.find('}')==strbuffer.length()-1)
            {
               //At the end - get everything upto there as values - comma separates values
               strbuffer=ReplaceAllWith(&strbuffer,',',';');
               strval=strval+";"+strbuffer; // add the string to the value string using ; as delimiter 
               this->Header[ToLowerCase(strkey)]=strval; //add to the header map
               multilineval=false; //reset the bool
               strval.clear();     
            }
            else
            {
               //Not at the end - this is odd
               throw "Error. Problem with hdr file - contains data after '}' value: "+strbuffer;
            }

         }
         if((strbuffer.find('=') == std::string::npos) && (strbuffer.find('}') == std::string::npos)) // no = and no }
         {
            strbuffer=TrimWhitespace(strbuffer);
            //No = in this string
            if(multilineval==false)//For the moment will just pass this as the key AND value into the map            
            {   
               this->Header[strbuffer]=strbuffer;
            }
            else if(!strbuffer.empty())
            {   
               // 14/02/2011 TrimPunctuation and then remove commas and replace with ';' - this could be dangerous
               strbuffer=TrimPunctuation(strbuffer);
               strbuffer=ReplaceAllWith(&strbuffer,',',';');
               strval=strval+";"+strbuffer; // add the string to the value string using ; as delimiter
            }
         }
      }

      //close the file
      hdrin.close();
      hdrin.clear();
   }
   return 1;
}

void BinaryReader::DatasizeFromDatatype()
{
   //Function that assigns the datasize from the datatype
   //datatype is an ENVI variable that describes the type of data in the BIL
   //datasize in the size in bytes
   switch(this->datatype)
   {
      case 1:
         //8 bit
         this->datasize=1; 
         break;
      case 2:
         //16 bit signed int
         this->datasize=2;
         break;
      case 3:
         //32 bit signed long int
         this->datasize=4;
         break;
      case 4:
         //32 bit floating point
         this->datasize=4;
         break;
      case 5:
         //64 bit double
         this->datasize=8;
         break;
      case 6:
         //2*32 bit real/imag
         this->datasize=4; //although 2 of size 4
         break;
      case 9:
         //2*64 bit real/imag
         this->datasize=8; //although 2 of size 8
         break;
      case 12:
         //16 bit unsigned
         this->datasize=2; //2 bytes
         break;
      case 13:
         //32 bit unsigned long int
         this->datasize=4;
         break;
      case 14:
         //64 bit signed long int
         this->datasize=8;
         break;
      case 15:
         //64 bit unsigned long int
         this->datasize=8;
         break;
      default:
         brinfo<<"Unknown data type "<<datatype<<std::endl;
         throw BRexception(brinfo.str());
         break;
      }
}


// Function now converts the passed key to lower case to search the header
std::string BinaryReader::FromHeader(std::string key,int itemnum,std::string THROW)
{
   //DEBUG statement
   DEBUGPRINT("Retrieving item "<<itemnum<<" from value relating to key "<<key);
   //Function to return the value string mapped by the key string
   //itemnum starts from 0 and items are separated by ;
   std::map<std::string,std::string>::iterator it;
   std::string str="";
   //Find the value from the map
   it=this->Header.find(ToLowerCase(key));
   if(it==this->Header.end()) //if the key doesnt exist return an empty string
   {
      if(THROW.compare("true")==0)
      {
         throw missingheaderitemerror+key;
      }
      return "";
   }
   else //else continue and get the related value to the key
      str=it->second;

   std::string retstr="";
   if(!str.compare(""))//This should never get called now that above check is in place, unless a key has an empty string as value (shouldnt happen)
   {
      //Empty string - probably means key is unknown
      return retstr;
   }
   else
   {
      if(str.find(';')!=std::string::npos)
      {         
         //This occurs when multiple objects where on the multiple lines
         //eg {;1
         //    ;56
         //    ;76}
         //Not an empty string - search for itemnums ; delimiter 
         unsigned int pos=0; //set position to start of string
         for(int i=0;i<=itemnum;i++)
         {
            pos=str.find(';',pos+1);
            if(pos == std::string::npos)
            {
               //there are less than itemnum ; in string
               return retstr;
            }
         }
         //found itemnum ;
         int index=str.find(';',pos+1);
         retstr=str.substr(pos+1,index-(pos+1));
         retstr=TrimPunctuation(retstr);
         return retstr;
      }
      else if(str.find(',')!=std::string::npos)
      {
         //This occurs when multiple objects were on the same line
         //eg {1, 56, 76} - only returns numbers
      
         //Commented out below 2 lines 10/09/2010 and replaced with new fnction to use comma delims and return anything 
         //str=RemoveAllBut(str,"0123456789.-+ ");
         //return GetItemFromString(str,itemnum);

         str=RemoveAllOf(str,"{}");
         return GetItemFromString(str,itemnum,',');
      }
      else
      {
         //unknown delimitng type.
         return "";
      }
   }
}

void BinaryReader::CannotFind(std::string key)
{
   brinfo<<"Cannot find \""<<key<<"\" in header file."<<std::endl;
   throw BRexception(brinfo.str()); //Throw an exception  

}

//Function to manually close the file and tidy up anything else
void BinaryReader::Close()
{
   //Close the binary file if open
   if(this->IsOpen())
   {
      //DEBUG statement
      DEBUGPRINT("Closing BinaryReader...");
      clearerr(filein);
      if(fclose(filein)==0)
         filein=NULL;   
   }   
}

bool BinaryReader::CheckCapacity(uint64_t numbytes)
{
   //Function to check number of bytes left to be read from file against number of bytes user wants to read
   //returns false if not enough bytes left in file
   uint64_t by=this->filesize - FileTell(this->filein);
   //DEBUGPRINT(by<<" "<<numbytes);
   if (by >= numbytes)
      return true;
   else
      return false;
}

//-------------------------------------------------------------------------
// Function to move file pointer to previous position
//-------------------------------------------------------------------------
void BinaryReader::GotoPreviousPosition()
{
   FileSeek(this->filein,this->prevpointerloc,SEEK_SET);
}

//-------------------------------------------------------------------------
// Function to dereference a char* buffer containing a single item read 
// from the data file, and return as a double
//-------------------------------------------------------------------------
double BinaryReader::DerefToDouble(char* cbuffer)
{
   //Pointers to various data types
   char* cp=NULL;
   short int* sip=NULL;
   float* fp=NULL;
   unsigned short int* usip=NULL;
   int* ip=NULL;
   unsigned int* uip=NULL;
   double* dp=NULL;

   //Check buffer is not null
   if(cbuffer==NULL)
   {
      throw "Attempt to convert NULL to double in DerefToDouble().";      
   }

   //The following switch statement allows the data in chtmp to be derefenced
   //to many types, allowing different data formats in the BIL to be used 
   switch(this->GetDataType())
   {
   case 1: //8-bit
      cp=(char*)(cbuffer);      
      return static_cast<double>(*cp);
   case 2: //16 bit signed int
      sip=reinterpret_cast<short int*>(cbuffer);      
      return static_cast<double>(*sip);
   case 3:
      ip=reinterpret_cast<int*>(cbuffer);
      return static_cast<double>(*ip);
   case 4: //float
      fp=reinterpret_cast<float*>(cbuffer);
      return static_cast<double>(*fp);
   case 5: //double
      dp=reinterpret_cast<double*>(cbuffer);
      return (*dp);
   case 12: //16 bit unsigned short int
      usip=reinterpret_cast<unsigned short int*>(cbuffer);
      return static_cast<double>(*usip);
   case 13: //32 bit unsigned int
      uip=reinterpret_cast<unsigned int*>(cbuffer);
      return static_cast<double>(*uip);
   default:
      throw "Unrecognised data type. Currently supports 8-bit, both signed and unsigned 16 & 32-bit integer, and 32 & 64-bit float";
      break;
   }

}

//-------------------------------------------------------------------------
// Function to dump the header to terminal or return in a string
//-------------------------------------------------------------------------
std::string BinaryReader::HeaderDump(bool ret)
{
   std::string dump="";
   std::map<std::string,std::string>::iterator it;

   for(it=Header.begin();it!=Header.end();it++)
   {
      if(((*it).first=="")&&((*it).second==""))
         continue;//ignore empty items 
      dump+=((*it).first)+" = "+((*it).second)+"\n";
   }
   if(ret==true)
      return dump;
   else
      std::cout<<dump<<std::endl;
   
   return "";
}

//-------------------------------------------------------------------------
// Function to tidy up strings when added to the header file
// - makes {, -> { and ,} -> } 
//-------------------------------------------------------------------------
std::string BinaryReader::TidyForHeader(std::string totidy,bool wrapinbraces)
{
   std::string ret=totidy;
   //replace any ; with , unless string starts with ; - may be a comment
   if(ret.at(0)!=';')
      ret=ReplaceAllWith(&ret,";",",\n");
   ret=ReplaceAllWith(&ret,"{,","{");
   ret=ReplaceAllWith(&ret,"{\n,","{\n");
   ret=ReplaceAllWith(&ret,",}","}");  
   ret=ReplaceAllWith(&ret,",\n}","\n}"); 
   if(wrapinbraces)
   {
      if(ret.find("{")!=0)
         ret.insert(0,"{");
      if(ret.rfind("}")!=ret.length()-1)
         ret.push_back('}');
   }
   return ret;
}

