//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef WRITER_H
#define WRITER_H

#include "filewriter.h"
#include "netcdfhelperclasses.h"

class NCWriter : public FileWriter
{
public:
   NCWriter(std::string filename,char openflag);
   NCWriter(std::string filename,unsigned int nrows,unsigned int ncols,unsigned int nbands,char openflag,DataType dt=float32);
   virtual ~NCWriter();

   //--------------------------------------
   //Implementation required for FileWriter
   //--------------------------------------
   int WriteLine(char* const data);
   int WriteBandLine(char* const data);
   int WriteBandLineSection(char* const data,const unsigned int numsamples_array,const unsigned int start, const unsigned int end);
   int WriteBandLineWithValue(const char xval); 
   int Close(){dataFile->~NcFile();return 0;}
   unsigned int GetDataSize()const{return datasize;}
   unsigned int GetDataType()const{return datatype;}

   template <class T>
   void SetNoDataValue(T value){dataFile->getVar(variable_to_write_to).setFill(true,value);}

   //--------------------------------------
   // Extra stuff for NCWriter
   //--------------------------------------
   void CreateNewVariable(std::string name,std::string type,std::vector<std::string> dimnames);
   void CreateNewVariable(std::string name,std::string type,std::string dimname);
   void ChangeVariableToWriteTo(std::string variablename);
   void AddMetadata(std::string name,std::string value){AddMetadata(name,value,"");}
   void AddMetadata(std::string name,std::string value,std::string variablename="");
   void AddMetadata(std::string name,void* vals,size_t len, NcType type, std::string variablename="");

   void WriteDataAt(char* const data,const unsigned band,const unsigned line,const unsigned int sample, 
                           const unsigned int samplelength,const unsigned int bandlength,const unsigned int linelength);
   //Exception class
   class NCexception : public FileWriterException
   {
   public:
      NCexception() : FileWriterException(){}
      NCexception(std::string ss) : FileWriterException(ss){}
      NCexception(std::string ss,const char* extra): FileWriterException(ss,extra){}
      const char* what() const throw()
      {
         return "A netCDF Writer Exception has occurred.";
      }
   };

private:
   NcFile* dataFile;
   std::string variable_to_write_to;
   std::vector<unsigned int> dims;
   size_t previouslinewritten,previousbandwritten;

   void WriteDataToVariable(char* const data, std::vector<size_t> startpoint, std::vector<size_t> numineachdir);
   void IncrementBandsWritten();
   NcVar GetVariableFromFile(std::string variablename) const;
   unsigned int NcTypeToEnviType(NcType type) const;
   void SetDataType();
   void SetDataSize();
};



#endif
