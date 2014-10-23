//------------------------------------------------------------------------- 
//Copyright (c) 2014 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#include "netcdfhelperclasses.h"

//-------------------------------------------------------------------------
// Function to open or create a netCDF file
//   arguments: Filename and Opening flag (a: append, w:overwrite)
//-------------------------------------------------------------------------
NcFile* OpenNetCDF(std::string netcdfFileName,char openflag)
{
   NcFile* dataFile=NULL;

   switch(openflag)
   {
      case 'a':
         std::cout<<"Will keep existing file and append to it."<<std::endl;
         dataFile=new NcFile(netcdfFileName,NcFile::write,NcFile::nc4classic);
         break;
      case 'r':
         std::cout<<"Attempt to open existing file read only."<<std::endl;
         dataFile=new NcFile(netcdfFileName,NcFile::read,NcFile::nc4classic);
         break;
      case 'o':
         std::cout<<"Will create a new file - overwrite existing file."<<std::endl;
         dataFile=new NcFile(netcdfFileName,NcFile::replace,NcFile::nc4classic);
         break;
      case 'w':
         std::cout<<"Will create a new file - fail if it exists file."<<std::endl;
         dataFile=new NcFile(netcdfFileName,NcFile::newFile,NcFile::nc4classic);
         break;
      default:
         break;
   }

   return dataFile;
}

//-------------------------------------------------------------------------
// Function to set the netCDF file to be in define mode or data mode
//-------------------------------------------------------------------------
void SetFileState(NcFile* nc,State s)
{
   switch(s)
   {
   case DATA:
      try{ncCheck(nc_enddef(nc->getId()),(char*)__FILE__,__LINE__);}
      catch(netCDF::exceptions::NcNotInDefineMode){std::cout<<"Already in data mode"<<std::endl;}
      break;
   case DEFINE:
      try{ncCheck(nc_redef(nc->getId()),(char*)__FILE__,__LINE__);}
      catch(netCDF::exceptions::NcInDefineMode){std::cout<<"Already in define mode."<<std::endl;}
      break;
   }
}


//-------------------------------------------------------------------------
// Function to write variable data (values) in to the netCDF file
// If the variable data has a startpoint defined it will use it together
// with the number of points to write in each dimension.
// If not, then it will write all the data from point 0
//-------------------------------------------------------------------------
void WriteVariableDataToFile(NcFile* nc,NCVariable &variable)
{
   NcVar var=nc->getVar(variable.Name());
   if(variable.Data->StartPoint().empty())
   {
      //std::cout<<"Writing array"<<std::endl;
      var.putVar(variable.Data->Data());
   }
   else
   {
      //std::cout<<"Writing section"<<std::endl;
      var.putVar(variable.Data->StartPoint(),variable.Data->NumInEachDir(),variable.Data->Data());
   }
}


//-------------------------------------------------------------------------
// Attribute Constructor takes name, type, length and data.
// The data are stored as chars
//-------------------------------------------------------------------------
NCAttribute::NCAttribute(std::string n,NcType t,size_t l,void* v)
{
   name=n;
   type=new NcType(t);
   len=l;
   unsigned int nbytes=t.getSize()*len;
   values=new char[nbytes];

   for(unsigned int i=0;i<nbytes;i++)
   {
      values[i]=static_cast<char*>(v)[i];
   }

}

//NCAttribute::NCAttribute(std::string n,std::string &v)
//{
//   name=n;
//   type=NULL;
//   len=v.length();
//   values=new char[v.length()];

//   for(unsigned int i=0;i<v.length();i++)
//   {
//      values[i]=v[i];
//   }
//}

NCAttribute::NCAttribute(std::string n,std::string v)
{
   name=n;
   type=NULL;
   len=v.length();
   values=new char[v.length()];

   for(unsigned int i=0;i<v.length();i++)
   {
      values[i]=v[i];
   }
}

//-------------------------------------------------------------------------
// Attribute destructor
//-------------------------------------------------------------------------
NCAttribute::~NCAttribute()
{
   if(values!=NULL)
      delete[] values;
}

//-------------------------------------------------------------------------
// Constructor taking the name string and vector of uint dimensions
//-------------------------------------------------------------------------
NCVariable::NCVariable(std::string n,std::vector<unsigned int> dims,std::string* dimnames)
{
   dimensions=dims;
   CommonConstructorCode(n,dimnames);
}

//-------------------------------------------------------------------------
// Constructor taking name string, number of dimensions and uint array 
// of dimensions
//-------------------------------------------------------------------------
NCVariable::NCVariable(std::string n,unsigned int ndims,unsigned int* dims,std::string* dimnames)
{
   for(unsigned int i=0;i<ndims;i++)
   {
      dimensions.push_back(dims[i]);
   }
   CommonConstructorCode(n,dimnames);
}

//-------------------------------------------------------------------------
// Code common to all constructors - called after constructor specific code
//-------------------------------------------------------------------------
void NCVariable::CommonConstructorCode(std::string n,std::string* dimnames)
{
   name=n;
   chunksizes.resize(dimensions.size());
   //Default type to float - can be changed with SetType method
   type=new NcType(ncFloat);
   fillvalue=NULL;
   if(dimnames!=NULL)
   {
      for(unsigned int i=0;i<dimensions.size();i++)
      {
         dimension_names.push_back(dimnames[i]);
      }
   }
   else
   {
      dimension_names.resize(dimensions.size());
   }
   //Default chunksize to the dimension size for all dimensions 
   //this can be changed with a call to Setchunksize
   for(unsigned int i=0;i<dimensions.size();i++)
   {
      chunksizes[i]=dimensions[i];
   }
   Data=NULL;
}

//-------------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------------
NCVariable::~NCVariable()
{
   delete type;
   if(Data!=NULL)
      delete Data;
}

//-------------------------------------------------------------------------
// Function to add the variable object to the open netCDF file object nc
// No longer assumes the file is in define mode - will change it everytime
// This may cause slow down but allows the option of adding data to the 
// variable at the same time to make function calls 'nicer'
//-------------------------------------------------------------------------
void NCVariable::AddToFile(NcFile* nc)
{
   SetFileState(nc,DEFINE);
   //Add the dimensions if not already in file
   for(unsigned int i=0;i<dimensions.size();i++)
   {
      //If the dimension already exists don't try and re-add it
      if(nc->getDim(dimension_names[i]).isNull())
      {
         std::cout<<"Dimension does not exist: "<<dimension_names[i]<<" Will add it."<<std::endl;
         nc->addDim(dimension_names[i],dimensions[i]);
      }
      else
      {
         std::cout<<"Dimension already exists: "<<dimension_names[i]<<std::endl;
      }
      ncdimensions.push_back(nc->getDim(dimension_names[i]));
   }

   //Add the variable if not in the file - else silently fail
   if(nc->getVar(name).isNull())
   {
      NcVar var=nc->addVar(name,(*type),ncdimensions);
      var.setChunking(NcVar::nc_CHUNKED,chunksizes);
      //Add the associated attributes
      for(std::vector<NCAttribute*>::iterator it=attributes.begin();it!=attributes.end();it++)
      {
         if((*it)->Type()==NULL)
         {
            //Use the string adding function
            var.putAtt((*it)->Name(),std::string((*it)->Data(),(*it)->Length()));
         }
         else
         {
            var.putAtt((*it)->Name(),*((*it)->Type()),(*it)->Length(),(*it)->Data());
         }
      }
      //Set the fill value
      if(fillvalue!=NULL)
         var.setFill(true,(void*)fillvalue);
      //Set the compression - currently hard coded
      var.setCompression(false,true,6);
      //Add data to the file if associated with variable
      if(Data!=NULL)
      {
         SetFileState(nc,DATA);
         WriteVariableDataToFile(nc,*this);
      }
   }
   else
   {
      std::cout<<"Variable already exists"<<std::endl;
   }
   
}

//-------------------------------------------------------------------------
// Function to add data with only a reference (pointer) to the data array
//-------------------------------------------------------------------------
void NCVariable::AddDataRef(void* data,std::vector<size_t>startpoint,std::vector<size_t>numineachdir,NcType* ntype)
{
   //startpoint is a vector containing info on where to write the data to in each dimension
   //numineachdir is the number of data points to write out in each direction
   if(Data!=NULL)
      delete Data;
   if(ntype==NULL)
      Data=new NCVariableData(data,type,startpoint,numineachdir,false);
   else
      Data=new NCVariableData(data,ntype,startpoint,numineachdir,false);
}

void NCVariable::AddDataRef(void* data,NcType* ntype)
{
   //startpoint is a vector containing info on where to write the data to in each dimension
   //numineachdir is the number of data points to write out in each direction
   if(Data!=NULL)
      delete Data;
   if(ntype==NULL)
      Data=new NCVariableData(data,type,false);
   else
      Data=new NCVariableData(data,ntype,false);
}

//-------------------------------------------------------------------------
// Function to add data via a hard copy of the data array
//-------------------------------------------------------------------------
void NCVariable::AddDataHard(void* data,std::vector<size_t>startpoint,std::vector<size_t>numineachdir,NcType* ntype)
{
   //startpoint is a vector containing info on where to write the data to in each dimension
   //numineachdir is the number of data points to write out in each direction
   if(Data!=NULL)
      delete Data;
   if(ntype==NULL)
      Data=new NCVariableData(data,type,startpoint,numineachdir,true);
   else
      Data=new NCVariableData(data,ntype,startpoint,numineachdir,true);
}

void NCVariable::AddDataHard(void* data,NcType* ntype)
{
   //startpoint is a vector containing info on where to write the data to in each dimension
   //numineachdir is the number of data points to write out in each direction
   if(Data!=NULL)
      delete Data;
   if(ntype==NULL)
      Data=new NCVariableData(data,type,true);
   else
      Data=new NCVariableData(data,ntype,true);
}

//-------------------------------------------------------------------------
// Constructor for data object taking start point and num of points in
// each dimension vectors
//-------------------------------------------------------------------------
NCVariableData::NCVariableData(void* d,NcType* t,std::vector<size_t>sp,std::vector<size_t>ndir,bool hard)
{
   startpoint=sp;
   numineachdir=ndir;
   datasize=t->getSize();
   refdata=static_cast<char*>(d);
   harddata=NULL;
   data=refdata;
   if(hard)
      HardCopy();
}

//-------------------------------------------------------------------------
// Constructor for data object
//-------------------------------------------------------------------------
NCVariableData::NCVariableData(void* d,NcType* t,bool hard)
{
   //Set the vectors to be empty
   startpoint.clear();
   numineachdir.clear();
   datasize=t->getSize();
   refdata=static_cast<char*>(d);
   harddata=NULL;
   data=refdata;
   if(hard)
      HardCopy();
}

//-------------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------------
NCVariableData::~NCVariableData()
{
   if(harddata!=NULL)
      delete[] harddata;
   refdata=NULL;
   data=NULL;
}


//-------------------------------------------------------------------------
// Function to take a hard copy of the data passed rather than only pointer
//-------------------------------------------------------------------------
void NCVariableData::HardCopy()
{
   if(numineachdir.empty())
      throw "Cannot use hard copy without specifying number of data points to write out.";

   unsigned int nbytes=0;
   for(std::vector<size_t>::iterator it=numineachdir.begin();it!=numineachdir.end();it++)
   {
      nbytes+=(*it)*datasize;
   }
   harddata=new char[nbytes];
   for(unsigned int i=0;i<nbytes;i++)
   {
      harddata[i]=static_cast<char*>(refdata)[i];
   }
   data=harddata;
}




