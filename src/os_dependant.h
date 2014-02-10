//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

//Classes for cross-platform OS specific stuff
#ifndef OS_DEPENDANT_H
#define OS_DEPENDANT_H

#include <iostream>
#include <string>
#include <stdint.h>
#include "commonfunctions.h"
#include "logger.h"


#ifdef _W32
   #include <windows.h> //For windows functions
#else
   #include <sys/statvfs.h> //For DiskSpace class
   #include <sys/utsname.h> //For ComputerInfo class
#endif

//-------------------------------------------------------------------------
// Class to find available disk space for a given location
//-------------------------------------------------------------------------
class DiskSpace
{
public:
   DiskSpace();
   ~DiskSpace();

   uint64_t GetTotalSpace(std::string diskname);
   uint64_t GetAvailableSpace(std::string diskname);

private:
   #ifdef _W32
      ULARGE_INTEGER totalall;
      ULARGE_INTEGER freeavailable;
   #else
      struct statvfs* databuffer; //Linux buffer to hold return from statvfs
   #endif
};

//-------------------------------------------------------------------------
// Class to get information about the computer
//-------------------------------------------------------------------------
class ComputerInfo
{
public:
   ComputerInfo();
   ~ComputerInfo();
   
   std::string GetOutput();

private:
   std::string host,domain,machine,system,version,release;
};


#endif
