//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------


#ifndef COMMANDLINE_H
#define COMMANDLINE_H

#include <iostream>
#include <string>
#include <map>
#include "commonfunctions.h"

#define optiononly "_NO_ARGUMENT_"

class CommandLine
{
public:
   CommandLine(char** args,int numargs); //Takes the command line and fills the elements map with the options and arguements (by calling Handle function)
   ~CommandLine();

   bool OnCommandLine(const std::string opt);//check if string opt is present on command line
   
   int CheckAvailableOptions(const std::string* avopts,const int, std::string* const badopts=NULL); //function to check command line options against avopts (strings of available options eg "-f")
                                                   //returns -ve number if options on cl are not in avopts [value == to number of incorrect options ] 
                                                   //returns 0 otherwise
   bool IsGood()const{return clgood;} // Check if error flag is set (1=good, 0=bad)

   std::string GetArg(const std::string opt);
   std::string ExeName(){return strExecutableName;} // Get the executable program name

   std::string GetArg(const std::string opt,const int argn);// returns the nth argument of the option opt
   int NumArgsOfOpt(const std::string key);// returns the number of arguments for the given key

   std::string ReviewCL(const std::string* availopts,const int number_of_possible_options,bool printout=false);//function to output the command line args/opts

   std::string ProgramUsage(const int number_of_possible_options, const std::string* availableopts,const std::string* optsdescription);

   std::string ReturnCLAsString();

   //Command line exception
   class CommandLineException
   {
      public:
      std::string info; //string of information about exception
      CommandLineException(){info="";} 
      CommandLineException(std::string ss){info=ss;} 

      const char* what() const throw()
      {
         return "A CommandLine Exception has occurred. Please use -help to get usage information.";
      }
   };

private:
   std::map<std::string,std::string,cmpstr> elements; //map to contain the options (keys) and arguments (values)
   bool Handle(char** args,int numargs); //Function to iterate through the numargs of argv[] and sort into options and arguements
   bool clgood;  //true if constructor finsihed successfully
   bool IsOpt(const std::string str);//Check if string str is an option or argument

   int* argcount; //number of arguments for each option
   void CountArgs(); //function to count the number of arguments per option

   std::string strExecutableName; //name of the program running
   std::string strFormattedCL;//formatted string of command line
};

#endif
