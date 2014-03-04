//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef BIL_H
#define BIL_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <string>
#include <cstring>
#include <exception>
#include <cerrno>
#include "binaryreader.h"

//Class for reading in BIL files. Supercedes the original BILReader class.

class BILReader:public BinaryReader
{
public:
   BILReader(std::string fname);
   ~BILReader();

   virtual void Readline(char* const chdata); //Read 1 line for all bands [equivalent number of bytes from current position]
   virtual void Readline(char* const chdata, unsigned int line) //Read specified line for all bands
         {Readlines(chdata,line,1);} //This is essentially a special case of readlines where numlines=1
   virtual void Readlines(char* const chdata, unsigned int startline, unsigned int numlines); //Read numlines lines of data from startline
   virtual void Readbytes(char* const chdata, unsigned long int bytes); //Read the specified number of bytes from current position
   virtual int Readband(char* const chdata, unsigned int band); //Reads the specified band
   virtual int Readbandline(char* const chdata, unsigned int band, unsigned int line); //reads the given line for the given band

   virtual double ReadCell(const unsigned int band,const unsigned int line, const unsigned int col);

   virtual void ReadlineToDoubles(double* const ddata,unsigned int line);

protected:

};

#endif
