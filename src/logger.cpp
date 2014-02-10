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
*     Logging class. Adds strings or char arrays to a stringstream object
*     Option to output the text to terminal or terminal + file]
*
*     Converted to static class - 6th February 2009
*     Added the Verbose and Debug options August 2011
***************************************************************************/

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>

#include "logger.h"

//Define static variables here
bool Logger::tofile=false;
short Logger::VERBOSE=0;
std::string Logger::filename("");
std::stringstream Logger::logtext;
std::ofstream Logger::logfile;
std::vector<std::string> Logger::warnonce;

//Basic constructor for logging to the terminal
Logger::Logger()
{
   tofile=false;
   VERBOSE=0;
   std::cout.precision(10);
   std::cout.setf(std::ios::unitbuf); //for unbuffered output into the GUI
}

Logger::Logger(short verbose=0)
{
   tofile=false;
   VERBOSE=verbose;
   std::cout.precision(10);
   std::cout.setf(std::ios::unitbuf); //for unbuffered output into the GUI
}

//constructor for logging to the terminal and file
Logger::Logger(std::string fname,short verbose=0)
{
   try
   {
      filename=fname;
      tofile=true;
      VERBOSE=verbose;
      logfile.open(fname.c_str(),std::ios::out);//open the log file
      if(!logfile.is_open())
      {
         throw LogException("Log file: "+fname+" failed to open.");//need to throw an exception
      }
   }
   catch(Logger::LogException e)
   {
      std::cout<<e.what()<<std::endl;
      std::cout<<e.info<<std::endl;
   }
}

Logger::~Logger()
{
   //if the file is open lets close it.
   Close();
}

void Logger::Close()
{
   //if the file is open lets close it.
   if(logfile.is_open())
   {
      logfile.close();
      logfile.clear();
   } 
}


//function to add text string to the stringstream object without outputting
void Logger::Add(const std::string text)
{
   logtext<<text<<std::endl;
}

//function to add text char to the stringstream object without outputting
void Logger::Add(const char* text)
{
   logtext<<text<<std::endl;
}

//function to add text string to the stringstream object with outputting
void Logger::Log(const std::string text)
{
   logtext<<text<<std::endl;
   Flush();
}

//function to add text char to the stringstream object with outputting
void Logger::Log(const char* text)
{
   logtext<<text<<std::endl;
   Flush();
}

//Function that outputs information if either verbose or debug levels set
void Logger::Verbose(const std::string text)
{
   if(VERBOSE>0)
   {
      logtext<<text<<std::endl;  
      Flush();
   }
}

//Function that outputs information if either verbose or debug levels set
void Logger::Verbose(const char* text)
{
   if(VERBOSE>0)
   {
      logtext<<text<<std::endl;  
      Flush();
   }
}

//Function that outputs information if only debug level set
void Logger::Debug(const std::string text)
{
   if(VERBOSE>=2)
   {
      logtext<<text<<std::endl;  
      Flush();
   }
}

//Function that outputs information if only debug level set
void Logger::Debug(const char* text)
{
   if(VERBOSE>=2)
   {
      logtext<<text<<std::endl;  
      Flush();
   }
}

//function to output string as error text
void Logger::Error(const std::string text)
{
   //Only output if text is not empty
   if(text.compare("")!=0)
   {
      logtext<<"########## ERROR #########"<<std::endl;
      logtext<<text<<std::endl;
      logtext<<"########## ERROR #########"<<std::endl;
      Flush();
   }
}

void Logger::Error(const char* text)
{
   //Only output if text is not empty
   if(strcmp(text,"")!=0)
   {
      logtext<<"########## ERROR #########"<<std::endl;
      logtext<<text<<std::endl;
      logtext<<"########## ERROR #########"<<std::endl;
      Flush();
   }
}

//function to output string as warning text
void Logger::Warning(const std::string text)
{
   //Only output if text is not empty
   if(text.compare("")!=0)
   {
      logtext<<"---------- WARNING ----------"<<std::endl;
      logtext<<text<<std::endl;
      logtext<<"---------- ------- ----------"<<std::endl;
      Flush();
   }
}

void Logger::Warning(const char* text)
{
   //Only output if text is not empty
   if(strcmp(text,"")!=0)
   {
      logtext<<"---------- WARNING ----------"<<std::endl;
      logtext<<text<<std::endl;
      logtext<<"---------- ------- ----------"<<std::endl;
      Flush();
   }
}

//function to output a warning once only, no matter how many times called
void Logger::WarnOnce(const char* text)
{
   //Only output if text is not empty
   if(strcmp(text,"")!=0)
   {
      //Check if message has appeared in warnonce vector before
      for(std::vector<std::string>::iterator it=warnonce.begin();it!=warnonce.end();it++)
      {
         if((*it).compare(text)==0)
            return;
      }
      //add it to the vector for next calls
      warnonce.push_back(std::string(text));
      //Output the message
      logtext<<"---------- WARNING ----------"<<std::endl;
      logtext<<text<<std::endl;
      logtext<<"---------- ------- ----------"<<std::endl;
      Flush();
   }  
}

//function to output a warning once only, no matter how many times called
void Logger::WarnOnce(const std::string text)
{
   //Only output if text is not empty
   if(text.compare("")!=0)
   {
      //Check if message has appeared in warnonce vector before
      for(std::vector<std::string>::iterator it=warnonce.begin();it!=warnonce.end();it++)
      {
         if((*it).compare(text)==0)
            return;
      }
      //add it to the vector for next calls
      warnonce.push_back(text);
      //Output the message
      logtext<<"---------- WARNING ----------"<<std::endl;
      logtext<<text<<std::endl;
      logtext<<"---------- ------- ----------"<<std::endl;
      Flush();
   }  
}


//function to write the stringstream to the file
void Logger::WriteToLogfile()
{
   //Check file is open .. this should always be the case
   try
   {
      if(logfile.is_open())
         logfile<<logtext.str().c_str();
      else
         throw LogException("Trying to write to a closed log file");
   }
   catch(LogException e)
   {
      std::cout<<e.what()<<std::endl;
      std::cout<<e.info<<std::endl;
   }
} 

 //function to output the string stream object (if tofile true outputs to screen + file else just to screen)
void Logger::Flush()
{
   //If logging to file do this first
   if(tofile==true)
   {
      WriteToLogfile();
   }
   //else skip the file logging
   WriteToScreen();
   logtext.str("");//empty the text string stream
}

//Function to output a formatted header for the start of a program including
//software name, version and description.
void Logger::FormattedInformation(std::string name,std::string version,std::string desc)
{
   std::string formatted;
   formatted= "/***************************************\n";
   formatted+="* Software Name: "+name+"\n";
   formatted+="* Version: "+version+"\n";
   formatted+="* Description: "+desc+"\n";
   formatted+= "\\***************************************\n";

   Log(formatted);
}
