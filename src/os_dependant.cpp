//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#include "os_dependant.h"

//-------------------------------------------------------------------------
//Constructor - create new databuffer struct
//-------------------------------------------------------------------------
DiskSpace::DiskSpace()
{
   #ifdef _W32
//      totalall=0;
//      freeavailable=0;
   #else
      databuffer=new struct statvfs;
   #endif
   
}

//-------------------------------------------------------------------------
//Desctuctor - delete databuffer struct
//-------------------------------------------------------------------------
DiskSpace::~DiskSpace()
{
   #ifndef _W32
      if(databuffer != NULL)
         delete databuffer;
   #endif
}

//-------------------------------------------------------------------------
//Return the available space 
//-------------------------------------------------------------------------
uint64_t DiskSpace::GetAvailableSpace(std::string diskname)
{
   #ifdef _W32 
   {
      if(GetDiskFreeSpaceEx(DirName(diskname).c_str(),&freeavailable,NULL,NULL))
         return (freeavailable.QuadPart);
      else
      {
         Logger::Warning("Failed to find amount of free disk space - there may not be enough for the file.");
         return 0;
      }
   }
   #else
   {
      statvfs(diskname.c_str(),databuffer);
      //return the size in MB
      return (databuffer->f_bsize * databuffer->f_bfree); 
   }
   #endif
   
}

//-------------------------------------------------------------------------
//Return the total space of the disk (in truncated MB)
//-------------------------------------------------------------------------
uint64_t DiskSpace::GetTotalSpace(std::string diskname)
{
   #ifdef _W32 
   {
      if(GetDiskFreeSpaceEx(DirName(diskname).c_str(),NULL,&totalall,NULL))
         return (totalall.QuadPart);
      else
      {
         Logger::Warning("Failed to find amount of total disk space.");
         return 0;
      }
   }
   #else
   {
      statvfs(diskname.c_str(),databuffer);
      //return the size in MB
      return (databuffer->f_frsize * databuffer->f_blocks); 
   }
   #endif   
}


//-------------------------------------------------------------------------
//Constructor for ComputerInfo class
//-------------------------------------------------------------------------
ComputerInfo::ComputerInfo()
{
   #ifdef _W32
   {
      TCHAR buffer[256]=TEXT("");
      DWORD dwSize = sizeof(buffer);
      //Get host name
      GetComputerNameEx(ComputerNameDnsHostname,buffer,&dwSize);
      host=std::string(buffer);
      //Clean up the buffer
      dwSize = _countof(buffer);
      ZeroMemory(buffer, dwSize);
      //Get domain name
      GetComputerNameEx(ComputerNameDnsDomain,buffer,&dwSize);
      domain=std::string(buffer);
      //Clean up the buffer
      dwSize = _countof(buffer);
      ZeroMemory(buffer, dwSize);

      //Get the windows version information
      OSVERSIONINFO osvi;
      osvi.dwOSVersionInfoSize = sizeof(OSVERSIONINFO);
      GetVersionEx(&osvi);
      //Get system
      system="Microsoft Windows";
      //Get version
      release=ToString(osvi.dwMajorVersion)+"."+ToString(osvi.dwMinorVersion);
      version=ToString(osvi.dwBuildNumber);

      //Get the processor architecture (machine)
      SYSTEM_INFO si;
      GetSystemInfo(&si);
      switch(si.wProcessorArchitecture)
      {
      case PROCESSOR_ARCHITECTURE_INTEL:
         machine="x86";
         break;
      case PROCESSOR_ARCHITECTURE_AMD64:
         machine="x64";
         break;
      case PROCESSOR_ARCHITECTURE_ARM:
         machine="ARM";
         break;
      case PROCESSOR_ARCHITECTURE_IA64:
         machine="Intel Itanium-based";
         break;
      case PROCESSOR_ARCHITECTURE_UNKNOWN:
         machine="Unknown";
         break;
      default:
         machine="Unknown";
         break;
      }
   }
   #else
   {
      struct utsname information;
      uname(&information);
      system=information.sysname;
      release=information.release;
      version=information.version;
      machine=information.machine;
      host=information.nodename;
      domain=information.domainname;
   }
   #endif
}

//-------------------------------------------------------------------------
//Destructor for ComputerInfo class
//-------------------------------------------------------------------------
ComputerInfo::~ComputerInfo()
{
}

//-------------------------------------------------------------------------
//Return a string containing the formatted output
//-------------------------------------------------------------------------
std::string ComputerInfo::GetOutput()
{
   std::stringstream strout;
   strout<<"Operating system: "<<system<<" Release: "<<release<< " Version: "<<version<<std::endl;
   strout<<"Machine information: "<<machine<<" Host name: "<<host<<" Domain name: "<<domain<<std::endl;
   return strout.str();
}

