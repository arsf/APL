//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

/**************************************************************************
*     Collection of functions that are required for various other 
*     cpp files (such as BILReader.cpp, commandline.cpp) 
***************************************************************************/

//#include <string>
#include "commonfunctions.h"
//#include "logger.h"


//String function to strip leading and trailing whitespace 
std::string TrimWhitespace(std::string str)
{
   //This will trim both leading and trailing whitespace characters
   std::string whitespace(" \t\n\r\f\v");
   size_t index=str.find_last_not_of(whitespace); //find last occurrance of non-whitespace
   if(index!=std::string::npos)
      str.erase(index+1); //erase everything after this character
   else
      str.clear();   //entire string is whitespace

   //Now trim leading characters
   index=str.find_first_not_of(whitespace);
   if(index!=std::string::npos)
      str.erase(0,index); //erase everything after this character

   return str;
}

//This function will trim the trailing punctuation of a string
std::string TrimPunctuation(std::string str)
{
   //This will trim trailing punctuation characters
   std::string punct(",;.:{}");
   size_t index=str.find_last_not_of(punct); //find last occurrance of non-whitespace
   if(index!=std::string::npos)
      str.erase(index+1); //erase everything after this character
   else
      str.clear();   //entire string is punctuation

   //Now trim leading characters
   //index=str.find_first_not_of(punct);
   //if(index!=std::string::npos)
   //   str.erase(0,index); //erase everything after this character

   return str;
}

//This function will trim the trailing punctuation of a string
std::string TrimTrailingChars(std::string str,const char* chars)
{
   //This will trim trailing characters
   std::string char_string(chars);
   size_t index=str.find_last_not_of(char_string); //find last occurrance of non-whitespace
   if(index!=std::string::npos)
      str.erase(index+1); //erase everything after this character
   else
      str.clear();   //entire string is made up of only characters from chars array

   return str;
}

//This function will trim the trailing punctuation of a string
std::string TrimLeadingChars(std::string str,const char* chars)
{
   //This will trim trailing characters
   std::string char_string(chars);
   size_t index=str.find_first_not_of(char_string); //find last occurrance of non-whitespace
   if(index!=std::string::npos)
      str.erase(0,index); //erase everything upto this character
   else
      str.clear();   //entire string is made up of only characters from chars array

   return str;
}


//This function will remove all occurrences of chars in str
std::string RemoveAllBut(std::string str,const char* chars)
{
   std::string char_string(chars);
   size_t index=str.find_first_not_of(char_string); //find first non-chars char
   while(index!=std::string::npos)
   {
      str.erase(index,1);
      index=str.find_first_not_of(char_string);
   }
   return str;
}

//This function will remove all occurrences of chars in str
std::string RemoveAllOf(std::string str,const char* chars)
{
   std::string char_string(chars);
   size_t index=str.find_first_of(char_string); //find first chars char
   while(index!=std::string::npos)
   {
      str.erase(index,1);
      index=str.find_first_of(char_string);
   }
   return str;
}

//This function will get the itemnum'th element from the string delimited by space
std::string GetItemFromString(std::string str,const unsigned int itemnum)
{
   size_t i=0,index=0,pos=0,endindex=0;
   if(itemnum==0)
   {
      index=str.find(' ',pos); //find first space char
      return str.substr(pos,index-pos);
   }
   else
   {
      while(i<itemnum)
      {
         index=str.find(' ',pos); //find first space char
         pos=index+1;
         i++;
      }   
      if(index!=std::string::npos)
      {
         endindex=str.find(' ',pos);
         if(endindex!=std::string::npos)
            return str.substr(pos,endindex-pos);
         else
            return str.substr(pos);
      }
      else
         return "";
   } 
   return "";
}

//This function will get the itemnum'th element from the string delimited by delim
std::string GetItemFromString(std::string str,const unsigned int itemnum,const char delim)
{
   size_t i=0,index=0,pos=0,endindex=0;
   if(itemnum==0)
   {
      index=str.find(delim,pos); //find first delim char
      return str.substr(pos,index-pos);
   }
   else
   {
      while(i<itemnum)
      {
         index=str.find(delim,pos); //find first delim char
         pos=index+1;
         i++;
      }   
      if(index!=std::string::npos)
      {
         endindex=str.find(delim,pos);
         if(endindex!=std::string::npos)
            return str.substr(pos,endindex-pos);
         else
            return str.substr(pos);
      }
      else
         return "";
   } 
   return "";
}

//Function to check that the string only contains numbers - throw exception otherwise
//returns the string given if passes.
//if NEGATIVE is TRUE (default) then it will allow negative numbers
std::string CheckNumbersOnly(std::string str,bool NEGATIVE)
{
   std::string integers("0123456789"); //create a string containing integers   
   size_t index=str.find_first_not_of(integers);//find first non-occurance
   std::string exception_string="Test failed - not an integer number: "+str;
   if(index == std::string::npos)
   {
      //Only integers
      return str;
   }
   else if((index == 0) && (str[0]=='-') && (NEGATIVE==true)) 
   {
      //Minus sign at start - test from reverse end
      index=str.find_last_not_of(integers);
      if(index == 0)
      {
         //Only integers and -ve sign
         return str;
      }
      else
      {
         throw exception_string;
      }
   }
   else
   {
      //Non integers found therefore unsafe to convert
      throw exception_string;
   }
   
}

//Function to convert a string to unsigned int
unsigned int StringToUINT(std::string str)
{
   std::string integers("0123456789"); //create a string containing integers
   //Ignore whitespace
   std::string nospace=TrimWhitespace(str);
   size_t index=nospace.find_first_not_of(integers);//find first non-occurance
   unsigned int retval=0; //return value
   if(index == std::string::npos)
   {
      //Only integers found...therefore safe(-ish) to convert
      retval=static_cast<unsigned int>(atoi(nospace.c_str()));
   }
   else
   {
      //Non integers found therefore unsafe to convert
      retval=0; //will have to return 0 as no other suitable value can be returned
   }
   return retval;    
}


//Function to convert a string to a signed int
int StringToINT(std::string str)
{
   std::string integers("0123456789"); //create a string containing integers
   //Ignore whitespace
   std::string nospace=TrimWhitespace(str);

   size_t index=nospace.find_first_not_of(integers);//find first non-occurance
   int retval=0; //return value
   if(index == std::string::npos)
   {
      //Only integers found...therefore safe(-ish) to convert
      retval=atoi(nospace.c_str());
   }
   else if((index == 0) && (nospace[0]=='-' )) //minus sign at start of string
   {
      index=nospace.find_first_not_of(integers,1);//check for other non-integers
      if(index == std::string::npos) //none found so convert to number
         retval=atoi(nospace.c_str());
   }
   else
   {
      //Non integers found therefore unsafe to convert
      retval=0; //will have to return 0 as no other suitable value can be returned
   }
   return retval;    
}

// function to count the number of occurences of c in str
int TotalOccurence(std::string str, const char c)
{
   //unsigned int l=0;
   size_t l=0;
   int count=0;

   while(l!=std::string::npos)
   {
      l=str.find(c,++l);
      ++count;
   }
   return count-1;
}

//Function to replace all occurences of "replace" with "sub" in string strtosearch
std::string ReplaceAllWith(const std::string* strtosearch, const std::string replace, const std::string sub)
{
   std::string newstr("");
   newstr.assign(*strtosearch);//copy strtosearch to newstr
   size_t index=strtosearch->find_first_of(replace);
   if(index==std::string::npos)
   {
      // replace is not in strtosearch
      return newstr;
   }
   else
   {
      size_t c_count=0;
      while(c_count<newstr.length())
      {
         if(newstr.substr(c_count,replace.length()).compare(replace)==0)
         {
            newstr.replace(c_count,replace.length(),sub);         
            c_count+=sub.length(); //add on length of inserted string
         }
         else
            c_count++;
      }
   }
   return newstr;
}

//Function to replace all occurences of "replace" with "sub" in string strtosearch
std::string ReplaceAllWith(const std::string* strtosearch, const char replace, const char sub)
{
   std::string newstr("");
   newstr.assign(*strtosearch);//copy strtosearch to newstr
   size_t index=strtosearch->find_first_of(replace);
   if(index==std::string::npos)
   {
      // replace is not in strtosearch
      return newstr;
   }
   else
   {
      size_t c_count=0;
      while(c_count<strtosearch->length())
      {
         if(strtosearch->at(c_count)==replace)
            newstr.replace(c_count,1,&sub,1);
         c_count++;
      }
   }
   return newstr;
}


//-------------------------------------------------------------------------
// Convert string to seconds of day - cannot use strptime() since not available
// under windows Time should be in format "hh mm ss.ssss"
//-------------------------------------------------------------------------
double ConvertTime(std::string timestr)
{
   double secofday=0;
   try
   {
      std::string hours=timestr.substr(0,timestr.find_first_of(' '));      
      std::string minutes=timestr.substr(timestr.find_first_of(' ')+1,(timestr.find_last_of(' ')-timestr.find_first_of(' ')-1));
      std::string seconds=timestr.substr(timestr.find_last_of(' ')+1);

      secofday=StringToUINT(hours)*3600 + StringToUINT(minutes)*60 + StringToDouble(seconds);
   }
   catch(std::exception& e)
   {
      std::cout<<"Failed to convert String time into a double: "<<timestr;
      exit(1);
   }

   return secofday;
}



//-------------------------------------------------------------------------
// Convert a date and time into GPS seconds of week
// date format "dd-mm-yyyy" time format "hh mm ss.ssss"
//-------------------------------------------------------------------------
double GetSecOfWeek(std::string datestr,std::string timestr)
{
   //First get the day of the week
   int dayofweek=GetDayOfWeek(datestr);
   //std::cout<<"Day of week: "<<datestr<<" "<<dayofweek<<std::endl;

   if((dayofweek<0)||(dayofweek>6))
      throw "Error converting date string to day of week";

   //Now get the seconds of day
   double secofday=ConvertTime(timestr);
   //std::cout<<"Sec of Day: "<<timestr<<" "<<secofday<<std::endl;
   //Now return the second of week   
   return dayofweek*24*3600 + secofday;
}

//-------------------------------------------------------------------------
// Return a numeric value for day of week (0-sunday, 6-saturday)
//-------------------------------------------------------------------------
//Datestr of format "dd-mm-yyyy"
int GetDayOfWeek(std::string datestr)
{
   //Create a time struct
   struct tm mytime={0};

   //Fill the time struct with the date
   FillTimeStruct(&mytime,datestr);

   //Create a time_t object - fills in the week day in the mytime tm struct
   //dont actually bother getting the time_t object, just run mktime
   mktime(&mytime);

   //std::cout<<"The day of the week is: "<<mytime.tm_wday<<std::endl;
   return mytime.tm_wday;
}

//-------------------------------------------------------------------------
// Fill in the time struct with the given date string
//-------------------------------------------------------------------------
void FillTimeStruct(tm* this_time,std::string datestr)
{
   //Deconstruct the date string
   std::string daystr=datestr.substr(0,datestr.find_first_of('-'));
   std::string monthstr=datestr.substr(datestr.find_first_of('-')+1,datestr.find_last_of('-')-datestr.find_first_of('-')-1);
   std::string yearstr=datestr.substr(datestr.find_last_of('-')+1);   

   //Fill in the day/month/year
   this_time->tm_mday=StringToINT(daystr); //1-31
   this_time->tm_mon=StringToINT(monthstr)-1; //0-11
   this_time->tm_year=StringToINT(yearstr)-1900;//years since 1900
}



size_t GetNumberOfItemsFromString(std::string str,const char delim)
{
   //Remove everything but the delimiter
   std::string temp=RemoveAllBut(str,&delim);
   //Count the delimiters and add 1 since e.g. "a" has no spaces in,"a b" has 1 space etc
   return temp.length()+1;
}

//-------------------------------------------------------------------------
// Function for use with command line to get path names with spaces in
// for windows compatibiltiy
//-------------------------------------------------------------------------
std::string GetExistingFilePath(std::string p,bool THROW)
{
   std::string path=CreatePath(p);
   //Now check if this is a viable path
   if(DoesPathExist(path))
   {
      //path is a valid existing file
      return path;
   }
   else
   {
      //path is not an existing file
      if(THROW==false)
         return "";
      else
         throw "File does not exist: "+p;
   }
}

//-------------------------------------------------------------------------
// Function to convert a ; separated string into space separated
// THIS IS NOT SAFE BUT BEST AT MOMENT
//-------------------------------------------------------------------------
std::string CreatePath(std::string p)
{
   std::string path;
   //Replace all ; with spaces
   path=ReplaceAllWith(&p,';',' ');
   return path;   
}

//std::string GetWritableFilePath(std::string p)
//{
//   std::string path;
//   //Replace all ; with spaces
//   path=ReplaceAllWith(&p,';',' ');
//   //Now check if this is a viable path   
//   if(IsFileWritable(path))
//   {
//      //File path exists and can be written to
//      return path;
//   }
//   else
//   {
//      return "";
//   }
//}

//-------------------------------------------------------------------------
// Function to check if a path exists and can be opened for reading
//-------------------------------------------------------------------------
bool DoesPathExist(std::string filename)
{
   std::ifstream fin;
   fin.open(filename.c_str(),std::ios::in);
   if(fin.is_open())
   {
      fin.close();
      fin.clear();
      return true;
   }
   else
      return false;
}

//-------------------------------------------------------------------------
// Function to check if a file or path is writable
//-------------------------------------------------------------------------
//bool IsFileWritable(std::string filename)
//{
//   std::ofstream fout;
//   fout.open(filename.c_str(),std::ios::out);
//   if(fout.is_open())
//   {
//      fout.close();
//      fout.clear();
//      return true;
//   }
//   else
//      return false;
//}

//-------------------------------------------------------------------------
// Function to convert string to double safely 
// If except is false it will not throw exceptions (true by default)
//-------------------------------------------------------------------------
double StringToDouble(std::string str,bool except)
{
   //Input string stream object
   std::istringstream is(str);
   double val=0;
   char ch;
   //Check if str fits into a double and there are no chars left over
   if ((!(is >> val)) || (is.get(ch)))
   {
      if(except==true)
         throw "Error in StringToDouble() conversion. Converting: "+str+" and got "+ToString(val);
   }
   return val;   
}

//-------------------------------------------------------------------------
// Function to make a string all lowercase characters
//-------------------------------------------------------------------------
std::string ToLowerCase(std::string tolower)
{
   //Need to scan string for all upper case chars and replace with lower case versions
   std::string UPPERCASE="ABCDEFGHIJKLMNOPQRSTUVWXYZ";
   std::string LOWERCASE="abcdefghijklmnopqrstuvwxyz";
   std::string::iterator it_u;
   std::string::iterator it_l;
   std::string inlower=tolower;

   for(it_u=UPPERCASE.begin(), it_l=LOWERCASE.begin();it_u<UPPERCASE.end(),it_l<LOWERCASE.end() ;it_u++,it_l++)
   {
       inlower=ReplaceAllWith(&inlower, (*it_u), (*it_l));
   }
   return inlower;
}

//-------------------------------------------------------------------------
// Function to create an absolute path from a relative path
//-------------------------------------------------------------------------
std::string AbsolutePath(std::string relpath)
{
   //Check if a relative path or absolute
   if((relpath[0]=='/')||(relpath[0]=='\\')||((relpath[1]==':')&&(relpath[2]=='\\')))
   {
      //Starts with a / so is an absolute path
      return relpath;
   }

   //Nasty limit to 1024 chars length
   char currentdir[1024];
   //Get the current (working) directory
   if(!GetCurrentDir(currentdir,1023))
   {
      throw "Unable to find current working path in AbsolutePath.";
   }

   //Join relative path to current working path
   #ifdef _W32
   std::string abspath=std::string(currentdir)+"\\"+relpath;
   #else
   std::string abspath=std::string(currentdir)+"/"+relpath;
   #endif

   return abspath;
}

//-------------------------------------------------------------------------
// Function to calculate progress (i.e. a progress counter)
// returns 0 if not at a 10% interval
//-------------------------------------------------------------------------
unsigned int PercentProgress(unsigned int line,unsigned int nlines,bool quiet)
{
   static unsigned int perccount=0;
   //Test for divide by 0 in modulo statement - if <1 then return 0
   if((nlines/10.0) < 1)
      return 0;      

   if( line % (unsigned int)(nlines/10.0) == (unsigned int)(nlines/10.0) -1 )
   {
      perccount+=(10.0);
      if(quiet == false)
      {
         #ifdef LOGGER_H
         Logger::Log("Approximate percent complete: "+ToString(perccount)+" at line "+ToString(line));      
         #else
         std::cout<<"Approximate percent complete: "<<perccount<<" at line "<<line<<std::endl;   
         #endif
      }
   
      return perccount;
   }
   else
      return 0;
}

//-------------------------------------------------------------------------
// Function to trim off filename and return directory name
//-------------------------------------------------------------------------
std::string DirName(std::string filename)
{
   //Find last of '/' or '\'
   size_t index=filename.find_last_of("\\/");
   std::string path=""; 
   if(index!=std::string::npos)
      path=filename.substr(0,index);
   else
      path="."; //no directory in filename
   return path;
}

//-------------------------------------------------------------------------
// Function to print out the error message for abnormal exit 
// with an unknown reason - includes __FILE__ and __LINE__ of function call
//-------------------------------------------------------------------------
void PrintAbnormalExitMessage(const char* file, int line, std::string& exename, const char* vers, const char* contactdetails,std::string cl,std::exception* exc)
{
   std::stringstream errorstream;
   errorstream<<"Exception: "<<exename<<", Version: "<<vers<<" has encountered an unexpected error and is aborting abnormaly. "
              <<"This IS a bug and we would appreciate being told about it so that we can fix it. Please contact "<<contactdetails
              <<" sending the command line used: \n "<<cl
              <<"\nwith the following additional information: "<<exc->what()<<" in file: "<<file<<" at line "<<line<<std::endl;
   #ifdef LOGGER_H
   Logger::Error(errorstream.str());
   #else
   std::cout<<errorstream.str()<<std::endl;
   #endif
}


//-------------------------------------------------------------------------
// Function to determine limits of an array of data when you want to ignore
// particular data values. minval and maxval will be updated with the 
// derived min and max values 
//-------------------------------------------------------------------------
void GetArrayLimits(double* array,unsigned int length, double &minval, double &maxval, double ignorevalue)
{
   //Set the min and max values to the highest/lowest possible
   minval=std::numeric_limits<double>::max();
   maxval=-std::numeric_limits<double>::max();
   //Loop through each element and find min/max ignoring data points with ignore value
   for(unsigned int i=0;i<length;i++)
   {
      if(array[i]==ignorevalue)
         continue;
      else
      {
         if(array[i]>maxval)
            maxval=array[i];
         if(array[i]<minval)
            minval=array[i];
      }
   }
}
