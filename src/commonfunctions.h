//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef COMMONFUNCTIONS_H
#define COMMONFUNCTIONS_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <stdint.h>
#include <limits>

//This is included to get absolute paths
#ifdef _W32
    #include <direct.h>
    #define GetCurrentDir _getcwd
#else
    #include <unistd.h>
    #define GetCurrentDir getcwd
#endif

//Header file to contain common functions that may be useful elsewhere

std::string TrimPunctuation(std::string str);
std::string TrimWhitespace(std::string str);
unsigned int StringToUINT(std::string str);
int StringToINT(std::string str);
double StringToDouble(std::string str,bool except=true);
std::string TrimTrailingChars(std::string str,const char* chars);
std::string TrimLeadingChars(std::string str,const char* chars);
int TotalOccurence(std::string str, const char c);

std::string ToLowerCase(std::string tolower);

//Functions to convert the int/float/etc ival into a string
//a template function is "nicer" but would need to be included in all files
inline std::string ToString(short int ival) {std::stringstream ss;  ss<<ival;   return ss.str();}
inline std::string ToString(unsigned short int ival) {std::stringstream ss;  ss<<ival;   return ss.str();}
inline std::string ToString(unsigned int ival) {std::stringstream ss;  ss<<ival;   return ss.str();}
inline std::string ToString(long ival) {std::stringstream ss;  ss<<ival;   return ss.str();}
inline std::string ToString(unsigned long ival) {std::stringstream ss;  ss<<ival;   return ss.str();}
#if __WORDSIZE != 64
   //Only define this for 32 bit since it is the same as unsigned long on 64 bit
   inline std::string ToString(uint64_t ival) {std::stringstream ss;  ss<<ival;   return ss.str();}
   inline std::string ToString(int64_t ival) {std::stringstream ss;  ss<<ival;   return ss.str();}
#endif
inline std::string ToString(int ival) {std::stringstream ss;  ss<<ival;   return ss.str();}
inline std::string ToString(float ival) {std::stringstream ss;  ss<<ival;   return ss.str();}
inline std::string ToString(double ival) {std::stringstream ss;  ss<<std::setprecision(10)<<ival;   return ss.str();}


std::string RemoveAllOf(std::string str,const char* chars);
std::string ReplaceAllWith(const std::string* strtosearch, const char replace, const char sub);
std::string ReplaceAllWith(const std::string* strtosearch, const std::string replace, const std::string sub);
std::string RemoveAllBut(std::string str,std::string char_string);
std::string GetItemFromString(std::string str,const unsigned int itemnum);
std::string GetItemFromString(std::string str,const unsigned int itemnum,const char delim);
size_t GetNumberOfItemsFromString(std::string str,std::string delim);

//Function to convert a Time string HH mm ss.ss into a double second of day
double ConvertTime(std::string timestr);
//Function to convert date string "dd-mm-yyyy" to day of week 
int GetDayOfWeek(std::string datestr,std::string format="dd-mm-yyyy");
//Function to get the GPS second of day from a date and time string
double GetSecOfWeek(std::string datestr,std::string timestr,std::string dateformat);
std::string pad(std::string str,char p,size_t len,bool atend=false);

//Function to check for file exist and readabilty
bool DoesPathExist(std::string filename);
//Function to check if file/path can be written to
bool IsFileWritable(std::string filename);
//Function to check status of files/paths
std::string CreatePath(std::string p);
std::string GetExistingFilePath(std::string p,bool THROW=false);
std::string GetWritableFilePath(std::string p);

std::string AbsolutePath(std::string relpath);

void FillTimeStruct(tm* this_time,std::string datestr,std::string format="");
void DesconstructDateString(std::string datestr,std::string format,int* date);

unsigned int PercentProgress(unsigned int line,unsigned int nlines,bool quiet=false);
std::string CheckNumbersOnly(std::string str,bool NEGATIVE=true);

std::string DirName(std::string filename);

void PrintAbnormalExitMessage(const char* file, int line, std::string& exename, const char* vers, const char* contactdetails,std::string cl,std::exception* exc=NULL);

void GetArrayLimits(double* array,unsigned int length, double &minval, double &maxval, double ignorevalue);

//struct to compare two strings using the standard string.h functions {for use in map methods in CommandLine}
#ifndef FUNCTION_CMPSTR
#define FUNCTION_CMPSTR
struct cmpstr
{
   bool operator()(const std::string c1, const std::string c2)const
   {
      return c1.compare(c2) < 0;
   }
};
#endif


#endif
