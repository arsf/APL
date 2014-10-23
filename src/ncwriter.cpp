//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#include "ncwriter.h"

//----------------------------------------------------------------------
// Constructor for NetCDF writer
//----------------------------------------------------------------------
NCWriter::NCWriter(std::string filename,char openflag)
{
   //Open / create a netcdf file
   dataFile=OpenNetCDF(filename,openflag);
   previouslinewritten=0;
   previousbandwritten=0;
   variable_to_write_to="";
}

//----------------------------------------------------------------------
// Constructor for NetCDF writer and gridded variable dimensions
//----------------------------------------------------------------------
NCWriter::NCWriter(std::string filename,unsigned int nrows,unsigned int ncols,unsigned int nbands,char openflag,DataType dt)
{
   //Open / create a netcdf file
   dataFile=OpenNetCDF(filename,openflag);

   //Need to add dimensions and variables
   //These variables are: 
   //    gridded_data with dimensions X,Y,Bands 
   //    X with dimension number of elements in X direction
   //    Y with dimension number of elements in Y direction
   
   //Get the dimensions into a vector
   dims.push_back(nbands);
   dims.push_back(nrows);
   dims.push_back(ncols);

   //Get a matched array of dimension names
   std::string* dimnames=new std::string[3];
   dimnames[0]="Bands";
   dimnames[1]="Rows";
   dimnames[2]="Columns"; 

   variable_to_write_to="Gridded_data";
   previouslinewritten=0;
   previousbandwritten=0;

   //Pointers to general attributes
   NCAttribute* coordinates=NULL;

   //----------------------------------------------------------------------
   //Create the gridded data variable object
   //----------------------------------------------------------------------
   NCVariable mappedVar(variable_to_write_to,dims,dimnames);
   mappedVar.SetChunksizeForDim(0,1);
   mappedVar.SetChunksizeForDim(1,1);
   mappedVar.SetChunksizeForDim(2,ncols);

   coordinates=new NCAttribute("coordinates","bands_rows_columns");
   mappedVar.AddAttribute(coordinates);

   //Set the type of the variable
   switch(dt)
   {
   case char8:
      mappedVar.SetType(ncChar);
      break;
   case int16:
      mappedVar.SetType(ncShort);
      break;
   case int32:
      mappedVar.SetType(ncInt);
      break;
   case float32:
      mappedVar.SetType(ncFloat);
      break;
   case float64:
      mappedVar.SetType(ncDouble);
      break;
   case uint16:
      mappedVar.SetType(ncUshort);
      break;
   case uint32:
      mappedVar.SetType(ncUint);
      break;
   default:
      throw "Unrecognised datatype in NCWriter Constructor.";
   }

   //add the variable to the file
   try
   {
      mappedVar.AddToFile(dataFile);
   }
   catch(netCDF::exceptions::NcStrictNc3& e)
   {
      //Make a more sensible error message for when unsigned data are requested
      throw "Unable to open the netCDF file in strict mode. This is probably due to requesting an unsigned data variable, as these are unsupported in strict mode.";
   }

   //Clean up attributes
   delete coordinates;

   //Final set up of class parameters
   SetDataSize();
   SetDataType();

}

//----------------------------------------------------------------------
// Destructor
//----------------------------------------------------------------------
NCWriter::~NCWriter()
{
   Close();
}

//----------------------------------------------------------------------
// Function to increment the band number which has just been written
// When it reaches max num bands it flips to 0 and increments the 
// line number that has just been written
//----------------------------------------------------------------------
void NCWriter::IncrementBandsWritten()
{
   previousbandwritten++;
   if(previousbandwritten==dims[0])
   {
      previousbandwritten=0;
      previouslinewritten++;
   }
}

//----------------------------------------------------------------------
// Write a line of data to the file based on current position in file
//----------------------------------------------------------------------
int NCWriter::WriteLine(char* const data)
{
   //Set the start points for writing - these are based on previous writes
   //and complete lines (for all bands) being written

   size_t spoint[3]={0,previouslinewritten,0};
   std::vector<size_t>startpoint(spoint,spoint+3);
   //Set the number in each dimension to write
   size_t ndir[3]={dims[0],1,dims[2]};
   std::vector<size_t>numineachdir(ndir,ndir+3);

   WriteDataToVariable(data,startpoint,numineachdir);
   //Increment the number of lines written counter as we have just
   //written a new line at this point
   previouslinewritten++;

   return 0; //required due to original declaration as returning int for bilwriter
}

//----------------------------------------------------------------------
// Function that actually does the writing
//----------------------------------------------------------------------
void NCWriter::WriteDataToVariable(char* const data, std::vector<size_t> startpoint, std::vector<size_t> numineachdir)
{
   //Put in to data mode if not already
   try{ncCheck(nc_enddef(dataFile->getId()),(char*)__FILE__,__LINE__);}
   catch(netCDF::exceptions::NcNotInDefineMode){}//silence

   //Get the variable
   NcVar var=GetVariableFromFile(variable_to_write_to);
   if(var.getType().getName().compare("byte")==0)
   {
      var.putVar(startpoint,numineachdir,(NcByte*)data);
   } 
   else if(var.getType().getName().compare("char")==0)
   {
      var.putVar(startpoint,numineachdir,(NcChar*)data);
   }
   else if(var.getType().getName().compare("short")==0)
   {
      var.putVar(startpoint,numineachdir,(NcShort*)data);
   }
   else if(var.getType().getName().compare("ushort")==0)
   {
      var.putVar(startpoint,numineachdir,(NcUshort*)data);
   }
   else if(var.getType().getName().compare("int")==0)
   {
      var.putVar(startpoint,numineachdir,(NcUshort*)data);
   }
   else if(var.getType().getName().compare("uint")==0)
   {
      var.putVar(startpoint,numineachdir,(NcUshort*)data);
   }
   else if(var.getType().getName().compare("float")==0)
   {
      var.putVar(startpoint,numineachdir,(NcFloat*)data);
   }
   else if(var.getType().getName().compare("double")==0)
   {
      var.putVar(startpoint,numineachdir,(NcDouble*)data);
   }
   else
   {
      throw "Data type cannot be written yet - easy to fix - do it!";
   }
}

//----------------------------------------------------------------------
// Write a band worth of data for a line to the current position
//----------------------------------------------------------------------
int NCWriter::WriteBandLine(char* const data)
{
   //Set the start point for writing - band 0, current line and starting sample
   size_t spoint[3]={previousbandwritten,previouslinewritten,0};
   std::vector<size_t>startpoint(spoint,spoint+3);
   //Set the number in each dimension to write - all bands, 1 line , number of samples defined by end and start
   size_t ndir[3]={1,1,dims[2]};
   std::vector<size_t>numineachdir(ndir,ndir+3);      

   WriteDataToVariable(data,startpoint,numineachdir); 
   IncrementBandsWritten();
   return 0; //required due to original declaration as returning int for bilwriter
}

//----------------------------------------------------------------------
// Write a section of a band of a line at current position
//----------------------------------------------------------------------
int NCWriter::WriteBandLineSection(char* const data,const unsigned int numsamples_array,const unsigned int start, const unsigned int end)
{
   //Set the start point for writing - current band, current line and starting sample
   size_t spoint[3]={previousbandwritten,previouslinewritten,start};
   std::vector<size_t>startpoint(spoint,spoint+3);
   //Set the number in each dimension to write - 1 band, 1 line , number of samples defined by end and start
   size_t ndir[3]={1,1,end-start+1};
   std::vector<size_t>numineachdir(ndir,ndir+3);      

   WriteDataToVariable(data,startpoint,numineachdir); 
   IncrementBandsWritten();
   return 0; //required due to original declaration as returning int for bilwriter
}

//----------------------------------------------------------------------
// Write a select amount of data at given position
//----------------------------------------------------------------------
void NCWriter::WriteDataAt(char* const data,const unsigned band,const unsigned line,const unsigned int sample, 
                           const unsigned int samplelength,const unsigned int bandlength,const unsigned int linelength)
{
   //Set the start point for writing
   size_t spoint[3]={band,line,sample};
   std::vector<size_t>startpoint(spoint,spoint+3);
   //Set the number in each dimension to write
   size_t ndir[3]={bandlength,linelength,samplelength};
   std::vector<size_t>numineachdir(ndir,ndir+3);      

   WriteDataToVariable(data,startpoint,numineachdir); 
}



//----------------------------------------------------------------------
// Private function to get the variable from the file or throw exception
//----------------------------------------------------------------------
NcVar NCWriter::GetVariableFromFile(std::string variablename) const
{
   if(dataFile->getVar(variablename).isNull())
   {
      throw "Variable does not exist in netCDF file: "+variablename;
   }
   return dataFile->getVar(variablename);
}

//----------------------------------------------------------------------
// Write a band worth of data for a line of using a particular value
//----------------------------------------------------------------------
int NCWriter::WriteBandLineWithValue(const char xval)
{
   NcVar var=GetVariableFromFile(variable_to_write_to);
   size_t size=var.getType().getSize();

   if((xval=='0')||(size==1))
   {
      char* zdata=new char[dims[2]*size]; //create array of data of size of 1 line for all bands
      for(unsigned int i=0;i<dims[2]*size;i++)
         zdata[i]=xval;
      int retval=WriteBandLineSection(zdata,0,0,dims[2]-1);
      delete[] zdata;
      return retval;
   }
   else
      throw "I don't think this function works how you want for the set up of this BIL file. FIX IT!";
} 

//----------------------------------------------------------------------
// Function to add metadata to netcdf file
//----------------------------------------------------------------------
void NCWriter::AddMetadata(std::string name,std::string value,std::string variablename)
{
   //File needs to be in define mode
   try{ncCheck(nc_redef(dataFile->getId()),(char*)__FILE__,__LINE__);}
   catch(netCDF::exceptions::NcInDefineMode){}//silence
   //Get the variable from the netcdf file
   NcVar var;
   if(variablename.compare("")!=0)
      var=GetVariableFromFile(variablename);
   else
      var=GetVariableFromFile(variable_to_write_to);

   std::cout<<name<<" "<<value<<std::endl;
   
   var.putAtt(name,value);
}

//----------------------------------------------------------------------
// Function to add metadata to netcdf file
//----------------------------------------------------------------------
void NCWriter::AddMetadata(std::string name,void* vals,size_t len, NcType type, std::string variablename)
{
   //File needs to be in define mode
   try{ncCheck(nc_redef(dataFile->getId()),(char*)__FILE__,__LINE__);}
   catch(netCDF::exceptions::NcInDefineMode){}//silence
   //Get the variable from the netcdf file
   NcVar var;
   if(variablename.compare("")!=0)
   {
      var=GetVariableFromFile(variablename);
   }
   else
   {
      var=GetVariableFromFile(variable_to_write_to);
   }

   //Get the number of bytes that the meta data take up   
   unsigned int nbytes=type.getSize()*len;
   //Create a char array of this size
   char* values=new char[nbytes];
   //Fill the array with the data cast as char
   for(unsigned int i=0;i<nbytes;i++)
   {
      values[i]=static_cast<char*>(vals)[i];
   }   
   var.putAtt(name,type,len,values);
}

//----------------------------------------------------------------------
// Function to add a new variable to the netcdf file
//----------------------------------------------------------------------
void NCWriter::CreateNewVariable(std::string name,std::string type,std::vector<std::string> dimnames)
{   
   //File needs to be in define mode
   try{ncCheck(nc_redef(dataFile->getId()),(char*)__FILE__,__LINE__);}
   catch(netCDF::exceptions::NcInDefineMode){}//silence
   if(dataFile->getVar(name).isNull())
   {
      try{NcVar var=dataFile->addVar(name,type,dimnames);}
      catch(netCDF::exceptions::NcNullDim)
      {
         for(std::vector<std::string>::iterator it=dimnames.begin();it!=dimnames.end();it++)
            dataFile->addDim(*it);
         NcVar var=dataFile->addVar(name,type,dimnames);
      }
   }
   else
      throw "Variable already exists in netCDF file: "+name;
}

//----------------------------------------------------------------------
// Function to add a new variable to the netcdf file
//----------------------------------------------------------------------
void NCWriter::CreateNewVariable(std::string name,std::string type,std::string dimname)
{   
   //File needs to be in define mode
   try{ncCheck(nc_redef(dataFile->getId()),(char*)__FILE__,__LINE__);}
   catch(netCDF::exceptions::NcInDefineMode){}//silence
   if(dataFile->getVar(name).isNull())
   {
      try{NcVar var=dataFile->addVar(name,type,dimname);}
      catch(netCDF::exceptions::NcNullDim)
      {
         //Try adding the dimension first
         dataFile->addDim(dimname);
         NcVar var=dataFile->addVar(name,type,dimname);
      }
   }
   else
      throw "Variable already exists in netCDF file: "+name;
}

void NCWriter::ChangeVariableToWriteTo(std::string variablename)
{
   variable_to_write_to=variablename;
   SetDataSize();
   SetDataType();
}

void NCWriter::SetDataSize()
{
   NcVar var=GetVariableFromFile(variable_to_write_to);
   datasize=(unsigned int)var.getType().getSize();
}

void NCWriter::SetDataType()
{
   NcVar var=GetVariableFromFile(variable_to_write_to);
   NcType type=var.getType();
   datatype=NcTypeToEnviType(type);
}

unsigned int NCWriter::NcTypeToEnviType(NcType type)const
{
   switch(type.getId())
   {
   case NC_BYTE:
   case NC_UBYTE:
      return 1; //ENVI data type only specifies 8-bit (not signed/unsigned)
   case NC_CHAR:   
      return 1; //ENVI data type only specifies 8-bit (not signed/unsigned)      
   case NC_USHORT:
      return 12; //ENVI data type
   case NC_SHORT:
      return 2; //ENVI data type 
   case NC_UINT:
      return 13; //ENVI data type 
   case NC_INT:
      return 3; //ENVI data type 
   case NC_FLOAT:
      return 4; //ENVI data type 
   case NC_DOUBLE:
      return 5; //ENVI data type 
   default:
      throw NCexception(type.getName()," is an unknown datatype for netCDF to envi type conversion.");
      break;
   }
}
