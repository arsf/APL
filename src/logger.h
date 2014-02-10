//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

class Logger
{
public:
   Logger();
   Logger(short verbose); //constructor for logging to the screen only
   Logger(std::string filename,short verbose); //constructor to logging to the screen and file

   ~Logger(); //destructor

   //Might be better to make Log function a template function to accept anything...?
   static void WarnOnce(const std::string text);
   static void WarnOnce(const char* text);
   static void Warning(const std::string text);
   static void Warning(const char* text);
   static void Error(const std::string text);//function to output string as error text
   static void Error(const char* text);//function to output string as error text
   static void Log(const std::string text); //function to add text string to the stringstream object and output
   static void Log(const char* text); //function to add text string to the stringstream object and output
   static void Add(const std::string text); //function to add text string to the stringstream object
   static void Add(const char* text); //function to add text string to the stringstream object
   static void Flush(); //function to output the string stream object (if tofile true outputs to screen + file else just to screen)
   static void Close();

   static void Verbose(const std::string text);
   static void Verbose(const char* text);
   static void Debug(const std::string text);
   static void Debug(const char* text);

   static void SetLevel(short v){VERBOSE=v;}

   static void FormattedInformation(std::string name,std::string version,std::string desc);

   //Exception for logger
   class LogException
   {
   public:
      LogException();
      LogException(std::string str):info(str){}; 
      std::string info; //contains optional information about exception
      const char* what() const throw()
      {
         return "A Logger Exception has occurred.";
      }
   };

private:

   static bool tofile; //true if output to file else false (set in constructor)
   static std::string filename; //log file name
   static std::stringstream logtext; //string stream to hold the text to output
   static void WriteToScreen(){std::cout<<logtext.str();} //function to write the stringstream to the screen
   static void WriteToLogfile(); //function to write the stringstream to the screen

   static std::ofstream logfile; //ofstream object for file output
   static short VERBOSE;
   static std::vector<std::string> warnonce;
};


#endif
