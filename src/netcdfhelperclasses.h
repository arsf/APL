//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef NETCDF_HELPER_CLASSES
#define NETCDF_HELPER_CLASSES

#include <netcdf>
#include <string>
#include <vector>

using namespace netCDF;

enum State {DEFINE, DATA};

//-------------------------------------------------------------------------
// NCAttribute class - a mini class for storing attribute data
//-------------------------------------------------------------------------
class NCAttribute
{   
public:
//   NCAttribute(std::string n,std::string &v);
   NCAttribute(std::string n,std::string v);
   NCAttribute(std::string n,NcType t,size_t l,void* v);
   ~NCAttribute();

   std::string Name()const{return name;}
   NcType* Type()const{return type;}
   size_t Length()const{return len;}
   char* Data()const{return values;}

private:
   std::string name;
   NcType* type;
   size_t len;
   char* values;

};

//-------------------------------------------------------------------------
// NCVariableData - a mini class for storing actual data to write into
// a variable
//-------------------------------------------------------------------------
class NCVariableData
{
public:
   NCVariableData(void* data,NcType* t,bool hard=false);   
   NCVariableData(void* data,NcType* t,std::vector<size_t>startpoint,std::vector<size_t>numineachdir,bool hard=false);
   ~NCVariableData();

   void* Data()const{return data;}
   std::vector<size_t>StartPoint()const{return startpoint;}
   std::vector<size_t>NumInEachDir()const{return numineachdir;}
private:
   std::vector<size_t>startpoint;
   std::vector<size_t>numineachdir;
   char* data;
   char* refdata;
   char* harddata;
   unsigned int datasize;

   void HardCopy();
};

//-------------------------------------------------------------------------
// NCVariable class - a class for storing variable information including 
//  functions for adding the variable to a netCDF file
//-------------------------------------------------------------------------
class NCVariable
{
public:
   NCVariable(std::string n,unsigned int ndims,unsigned int* dims,std::string* dimnames=NULL);   
   NCVariable(std::string n,std::vector<unsigned int> dims,std::string* dimnames=NULL);
   ~NCVariable();

   std::string Name()const{return name;}
   NcType* Type()const{return type;}

   void SetChunksize(std::vector<size_t> c){chunksizes=c;}
   void SetChunksizeForDim(unsigned int dim,unsigned int c){chunksizes[dim]=c;}
   void SetType(NcType t){type=new NcType(t);}
   void SetDimensionName(unsigned int d,std::string n){dimension_names[d]=n;}
   void SetFillValue(void* f){fillvalue=f;}

   void AddToFile(NcFile* nc);
   void AddAttribute(NCAttribute* att){attributes.push_back(att);}
   void AddDataRef(void* data,std::vector<size_t>startpoint,std::vector<size_t>numineachdir,NcType* ntype=NULL);
   void AddDataRef(void* data,NcType* ntype=NULL);
   void AddDataHard(void* data,std::vector<size_t>startpoint,std::vector<size_t>numineachdir,NcType* ntype=NULL);
   void AddDataHard(void* data,NcType* ntype=NULL);

   //void* pData()const{return Data->Data();}
   NCVariableData* Data;

private:
   void CommonConstructorCode(std::string n,std::string* dimnames);

   std::string name;
   std::vector<unsigned int> dimensions;
   std::vector<size_t> chunksizes;
   std::vector<std::string> dimension_names;
   std::vector<NcDim> ncdimensions;

   NcType* type;
   std::vector<NCAttribute*> attributes;
   void* fillvalue;

};

//-------------------------------------------------------------------------
// Various functions to help with netCDF files
//-------------------------------------------------------------------------
NcFile* OpenNetCDF(std::string netcdfFileName,char openflag);
void SetFileState(NcFile* nc,State s);
void WriteVariableDataToFile(NcFile* nc,NCVariable &var);

#endif


