//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#include "bilwriter.h"

#ifndef BILWRITERDEBUG
   #define DEBUGPRINT(X)  
#else
   #define DEBUGPRINT(X) std::cout<< X <<std::endl;
#endif

//constructor to open a new file filename and write data of type datasize 
//(e.g. 1 byte (char),2 byte (uint) etc)
//if file exists it will throw an error
BILWriter::BILWriter(std::string fname,unsigned int dsize)
{
   //DEBUG output
   DEBUGPRINT("Entering BILWriter Constructor...");
   this->isgood=true;
   //Assign name and data size
   this->filename=fname;
   this->datasize=dsize;
   //Assign other vars to 0
   this->numrows=0;
   this->numsamples=0;
   this->numbands=0;
   
   //Add the "bil" interleave format to the hdr object
   std::string forhdr;
   forhdr="ENVI"; //needed so that ENVI recognises the file
   AddToHdr(forhdr);
   forhdr="interleave = bil";
   AddToHdr(forhdr);

   //Chceck if the file already exist
   std::ifstream fin;
   fin.open(this->filename.c_str());
   if(fin.is_open())
   {
      //File exists and is open
      this->isgood=false;
      fin.close();
      fin.clear();
      //DEBUG output
      DEBUGPRINT("BIL output file already exists. Will not overwrite unless told to do so explicitly (use 'a'/'w' to force append/overwrite)");
      //Throw an exception
      throw BILexception("BIL output file already exists. Will not overwrite unless told to do so explicitly (use 'a'/'w' to force append/overwrite)");
   }
   else
   {
      //File hasnt opened (still might exist)
      //Assume it doesnt exist...?      
      //Open output file
      this->fileout.open(this->filename.c_str(),std::ios::binary);
      //Add a check here to see if it opened
      if(!fileout.is_open())
      {
         this->isgood=false;
         throw BILexception("BIL output file failed to open: "+filename);         
      }
   }

}

//constructor to open file using method cmethod ('a' - append, 'w' -overwrite)                                                      
BILWriter::BILWriter(std::string fname,unsigned int dsize,const char cmethod)
{
   //DEBUG output
   DEBUGPRINT("Entering BILWriter Constructor...trying to open file "<<fname<<" using method -"<<cmethod<<" to output data of size "<<dsize<<" bytes.");
   this->isgood=true;
   //Assign name and data size
   this->filename=fname;
   this->datasize=dsize;
   //Assign other vars to 0
   this->numrows=0;
   this->numsamples=0;
   this->numbands=0;
   
   //Add the "bil" interleave format to the hdr object
   std::string forhdr;
   forhdr="ENVI"; //needed so that ENVI recognises the file
   AddToHdr(forhdr);
   forhdr="interleave = bil";
   AddToHdr(forhdr);

   //Check cmethod
   switch(cmethod)
   {
   case 'a':
      this->fileout.open(this->filename.c_str(),std::ios::binary | std::ios::app);
      break;
   case 'w':
      this->fileout.open(this->filename.c_str(),std::ios::binary | std::ios::trunc);
      break;
   default:
      this->isgood=false; //unrecognized option
      throw BILexception("Unknown method of writing to file - use 'a' or 'w'.");
      break;
   }

   //Add a check here to see if the file opened
   if(!fileout.is_open())
   {
      this->isgood=false;
      throw BILexception("BIL output file failed to open: "+filename);         
   }

}

//constructor to assign BIL dimensions at "start up" as well as above methods
//This is the recommended way of using the BILWriter to avoid disaster
BILWriter::BILWriter(std::string fname, const unsigned int dtype,const unsigned int nrows,const unsigned int nsamps,const unsigned int nbands,const char cmethod)
{
   //DEBUG output
   DEBUGPRINT("Entering BILWriter Constructor...trying to open file "<<fname<<" using method -"<<cmethod<<" to output data of type "<<dtype<<" and"
               <<nsamps<<" samples, "<<nrows<<" lines and "<<nbands<<" bands.");
   this->isgood=true;
   //Assign name and data size
   this->filename=fname;
   //Set the datasize from datatype
   this->SetSizeFromType(dtype);//do this before assigning datatype since func includes checks
   //this->datatype=dtype;    //this is done in n above function now
   //Assign other vars
   this->numrows=nrows;
   this->numsamples=nsamps;
   this->numbands=nbands;

   //Add these to the hdr object too
   //These are the ENVI standard
   std::string forhdr;
   forhdr="ENVI"; //needed so that ENVI recognises the file
   AddToHdr(forhdr);
   time_t nowtime; //time obj
   time(&nowtime);//get the time now
   
   forhdr="description = { "+bildescriptionstring+" on "+ctime(&nowtime)+" }";
   AddToHdr(forhdr);
   forhdr="lines = "+ToString(this->numrows);
   AddToHdr(forhdr);
   forhdr="samples = "+ToString(this->numsamples);
   AddToHdr(forhdr);
   forhdr="bands = "+ToString(this->numbands);
   AddToHdr(forhdr);
   forhdr="interleave = bil";
   AddToHdr(forhdr);
   forhdr="data type = "+ToString(this->datatype);
   AddToHdr(forhdr);   

   //Check cmethod
   switch(cmethod)
   {
   case 'a':
      this->fileout.open(this->filename.c_str(),std::ios::binary | std::ios::app);
      break;
   case 'w':
      this->fileout.open(this->filename.c_str(),std::ios::binary | std::ios::trunc);
      break;
   default:
      this->isgood=false; //unrecognized option
      throw BILexception("Unknown method of writing to file - use 'a' or 'w'.");
      break;
   }

   //Add a check here to see if the file opened
   if(!fileout.is_open())
   {
      this->isgood=false;
      throw BILexception("BIL output file failed to open: "+filename);         
   }
}

//Destructor 
BILWriter::~BILWriter()
{
   //DEBUG output
   DEBUGPRINT("Entering BILWriter destructor...");

   this->Close();
}

//Close the file, output the hdr information, check nrows/samples/bands agree with written filesize and call the destructor
int BILWriter::Close()
{
   //DEBUG output
   DEBUGPRINT("Closing the BIL file...");
   int retval=1;
   //Close the file if still open
   if(this->fileout.is_open())
   {
      this->fileout.close();
      this->fileout.clear();
      //Output the hdr file
      PrepareHeader();
      retval=WriteHeader();
   }  
   return retval;
}
 
//Creates and outputs the associated header file
int BILWriter::WriteHeader()
{
   //DEBUG output
   DEBUGPRINT("Writing out the Headerfile...");
   DEBUGPRINT(this->hdrtext.str());
   //Open the hdr file 
   std::string hdrfilename="";
   //hdrfilename=this->filename.substr(0,this->filename.length()-4)+".hdr";
   hdrfilename=this->filename+".hdr";
   std::ofstream hdrout;
   hdrout.open(hdrfilename.c_str());
   if(!hdrout.is_open())
   {
      //An error has meant the file cannot be opened
      //Add some info to a stream info obj
      throw BILexception("Cannot open output header file: "+hdrfilename);
      return -1;
   }
   else
   {
      //File has opened ok
      hdrout<<this->hdrtext.str();
      hdrout.close();
      hdrout.clear();
   }
   return 1;
}

//Function to write out blank (zeroed) data
//int BILWriter::WriteLineZeroes()
//{  
//   char* zdata=new char[this->numbands*this->numsamples*this->datasize]; //create array of zero data of size of 1 line for all bands
//   for(unsigned int i=0;i<this->numbands*this->numsamples*this->datasize;i++)
//      zdata[i]=0;
//   int retval=WriteLineSection(zdata,this->numsamples,0,numsamples-1);
//   delete[] zdata; //delte array
//   zdata=NULL;
//   return retval;
//}

//write out data with constant value over the line
int BILWriter::WriteLineWithValue(const char xval)
{
   char* zdata=new char[this->numbands*this->numsamples*this->datasize]; //create array of data of size of 1 line for all bands
   for(unsigned int i=0;i<this->numbands*this->numsamples*this->datasize;i++)
      zdata[i]=xval;
   int retval=WriteLineSection(zdata,this->numsamples,0,numsamples-1);
   delete[] zdata; //delte array
   zdata=NULL;
   return retval;
}

//Function to write a line of data of 1 band to the BIL file
//return +ve means it has written the data, -ve a problem occurred
int BILWriter::WriteBandLine(char* const data)
{
   DEBUGPRINT("trying to write a line of data...");
   //Check nsamps is known
   if(this->numsamples==0)
   {
      //Cannot out a line of data if we dont know its size
      this->bilinfo<<"Number of samples is unknown so cannot write out a line of data."<<std::endl;
      return -1;      
   }     
   //Check datasize is known
   if(this->datasize==0)
   {
      //Cannot output a line if we dont know what bytesize to use (e.g. 1btye data, 4 byte data etc)
      this->bilinfo<<"Size of data to output is unknown so cannot output a line of data."<<std::endl;
      return -1;
   }   
   //Check file is open
   if(!this->fileout.is_open())
   {
      //File is not open for some reason
      this->bilinfo<<"The BIL file is closed. Cannot output a line of data."<<std::endl;
      return -1;
   }

   //Lets try and write the data
   this->fileout.write(data,numsamples*datasize);
   if(this->fileout.bad())
   {
      this->bilinfo<<"A problem has occurred writing the line of data to file: "<<this->filename<<std::endl;
      return -1;   
   }
   return 1;
}

//Function to write a line of data of all bands to the BIL file
//return +ve means it has written the data, -ve a problem occurred
//int BILWriter::WriteLine(const char* data)
//{
//   DEBUGPRINT("trying to write a line of data...");
//   //Check nsamps is known
//   if(this->numsamples==0)
//   {
//      //Cannot out a line of data if we dont know its size
//      this->bilinfo<<"Number of samples is unknown so cannot write out a line of data."<<std::endl;
//      return -1;      
//   }     
//   //Check nbands is known
//   if(this->numbands==0)
//   {
//      //Cannot out a line of data for all bands if we dont know its size
//      this->bilinfo<<"Number of bands is unknown so cannot write out a line of data."<<std::endl;
//      return -1;      
//   } 
//   //Check datasize is known
//   if(this->datasize==0)
//   {
//      //Cannot output a line if we dont know what bytesize to use (e.g. 1btye data, 4 byte data etc)
//      this->bilinfo<<"Size of data to output is unknown so cannot output a line of data."<<std::endl;
//      return -1;
//   }   
//   //Check file is open
//   if(!this->fileout.is_open())
//   {
//      //File is not open for some reason
//      this->bilinfo<<"The BIL file is closed. Cannot output a line of data."<<std::endl;
//      return -1;
//   }
//
//   //Lets try and write the data
//   this->fileout.write(const_cast<char*>(data),numsamples*datasize*numbands);
//   if(this->fileout.bad())
//   {
//      this->bilinfo<<"A problem has occurred writing the line of data to file: "<<this->filename<<std::endl;
//      return -1;   
//   }
//   return 1;
//}

//Function to write a section of line of data of all bands to the BIL file
//return +ve means it has written the data, -ve a problem occurred
//start and end are the cells to start outputting from and to end the output
//eg 0,n-1 outputs the full line
//  50,150 outputs 101 cells of data, from (and including) cell 50 to cell 150
int BILWriter::WriteLineSection(char* const data,const unsigned int numsamples_array,const unsigned int start, const unsigned int end)
{
   DEBUGPRINT("trying to write a section of a line of data...");

   //check start < end < numsamples of array
   if(start > end)
   {
      //Return failure      
      this->bilinfo<<"Start pixel number to output must be less than end pixel number."<<std::endl<<"Start: "<<start<<" End: "<<end<<std::endl;
      return -1;  
   }
   if(end >= numsamples_array)
   {      
      //Return failure      
      this->bilinfo<<"End pixel number must be less than the number of samples of the passed array"<<std::endl<<"End: "<<end<<" Num samples: "<<this->numsamples<<std::endl;
      return -1; 
   }
   //Check nsamps is known
   if(this->numsamples==0) //Jan 19th 2009 - This check is no longer required now that we pass start and end points to output
   {
      //Cannot out a line of data if we dont know its size
      this->bilinfo<<"Number of samples is unknown so cannot write out a line of data."<<std::endl;
      return -1;      
   }     
   //Check nbands is known
   if(this->numbands==0)
   {
      //Cannot out a line of data for all bands if we dont know its size
      this->bilinfo<<"Number of bands is unknown so cannot write out a line of data."<<std::endl;
      return -1;      
   } 
   //Check datasize is known
   if(this->datasize==0)
   {
      //Cannot output a line if we dont know what bytesize to use (e.g. 1btye data, 4 byte data etc)
      this->bilinfo<<"Size of data to output is unknown so cannot output a line of data."<<std::endl;
      return -1;
   }   
   //Check file is open
   if(!this->fileout.is_open())
   {
      //File is not open for some reason
      this->bilinfo<<"The BIL file is closed. Cannot output a line of data."<<std::endl;
      return -1;
   }

   DEBUGPRINT("Start: "<<start<<" End: "<<end<<" datasize "<<datasize<<" numsamplesarray "<<numsamples_array);
   //Lets try and write the data
   for(unsigned int b=0;b<numbands;b++)
   { 
      DEBUGPRINT("Writing from sample..."<<b*numsamples_array + start<<" for "<<(end - start +1)*datasize<<"bytes.");
      //Write out the x bytes for each line of each band using the same starting pixel
      //need cell (i*numsamples_array + start)*datasize because each iteration we need to skip 'i' bands, and then move to the first cell to output
      //which is 'start' cells along. But data is in chars so need to multiply by datasize to get to right point
      this->fileout.write((&data[(b*numsamples_array + start)*datasize]),(end - start +1)*datasize);
   }

   if(this->fileout.bad())
   {
      this->bilinfo<<"A problem has occurred writing the line of data to file: "<<this->filename<<std::endl;
      return -1;   
   }
   return 1;
}

int BILWriter::WriteLine(char* const data)
{
   //Check nsamps is known
   if(this->numsamples==0) //Jan 19th 2009 - This check is no longer required now that we pass start and end points to output
   {
      //Cannot out a line of data if we dont know its size
      this->bilinfo<<"Number of samples is unknown so cannot write out a line of data."<<std::endl;
      return -1;      
   }     
   //Check nbands is known
   if(this->numbands==0)
   {
      //Cannot out a line of data for all bands if we dont know its size
      this->bilinfo<<"Number of bands is unknown so cannot write out a line of data."<<std::endl;
      return -1;      
   } 
   //Check datasize is known
   if(this->datasize==0)
   {
      //Cannot output a line if we dont know what bytesize to use (e.g. 1btye data, 4 byte data etc)
      this->bilinfo<<"Size of data to output is unknown so cannot output a line of data."<<std::endl;
      return -1;
   }   
   //Check file is open
   if(!this->fileout.is_open())
   {
      //File is not open for some reason
      this->bilinfo<<"The BIL file is closed. Cannot output a line of data."<<std::endl;
      return -1;
   }

   //Lets try and write the data
   this->fileout.write(data,numbands*numsamples*datasize);

   if(this->fileout.bad())
   {
      this->bilinfo<<"A problem has occurred writing the line of data to file: "<<this->filename<<std::endl;
      return -1;   
   }
   return 1;
}

//Prepare the header file stringstream to contain the relevant info
void BILWriter::PrepareHeader()
{
   DEBUGPRINT("Preparing header for output...");
   //First check that the dimensions are in the stream
   std::string hdrstring=this->hdrtext.str();
   if(hdrstring.find("samples = ")==std::string::npos)
   {
      //number of samples is not present ... try to add it
      if(this->numsamples > 0)
         this->hdrtext<<"samples = "<<ToString(this->numsamples)<<std::endl;
      else
      {
         //Number of samples is set to 0 or -ve. this is not good
         this->bilinfo<<"The number of samples for this BIL is "<<ToString(this->numsamples)<<std::endl;
         throw BILexception(bilinfo.str());
      }
   }

   if(hdrstring.find("bands = ")==std::string::npos)
   {
      //number of bands is not present ... try to add it
      if(this->numbands > 0)
         this->hdrtext<<"bands = "<<ToString(this->numbands)<<std::endl;
      else
      {
         //Number of bands is set to 0 or -ve. this is not good
         this->bilinfo<<"The number of bands for this BIL is "<<ToString(this->numbands)<<std::endl;
         throw BILexception(bilinfo.str());
      }
   }   

   if(hdrstring.find("lines = ")==std::string::npos)
   {
      //number of lines is not present ... try to add it
      if(this->numrows > 0)
         this->hdrtext<<"lines = "<<ToString(this->numrows)<<std::endl;
      else
      {
         //Number of lines is set to 0 or -ve. this is not good
         this->bilinfo<<"The number of lines for this BIL is "<<ToString(this->numrows)<<std::endl;
         throw BILexception(bilinfo.str());
      }
   }

   // now check that interleave is set to bil
   if(hdrstring.find("interleave ")==std::string::npos)
      this->hdrtext<<"interleave = bil"<<std::endl;            //Set interleave to bil
   else if(hdrstring.find("interleave = bil")==std::string::npos)
   {
      //Interlave is set to something other than bil
      this->bilinfo<<"Header string has interleave not equal to bil "<<std::endl;
      throw BILexception(bilinfo.str());
   }
   
   //Datatype should be output for ENVI hdr compatibility
   if(hdrstring.find("data type = ")==std::string::npos)
      this->hdrtext<<"data type = "<<this->datatype<<std::endl;  
   
}

void BILWriter::SetSizeFromType(const unsigned int type)
{
   DEBUGPRINT("Setting data type "<<type);
   //Set the size of data given by the type, and also check if a known datatype
   switch(type)
   {
   case uchar8:
      this->datasize=1; //1byte data
      this->datatype=1; //ENVI data type only specifies 8-bit (not signed/unsigned)
      break;
   case char8:
      this->datasize=1; //1byte data
      this->datatype=1; //ENVI data type only specifies 8-bit (not signed/unsigned)      
      break;
   case uint16:
      this->datasize=2; //2 byte data
      this->datatype=12; //ENVI data type
      break;
   case int16:
      this->datasize=2; //2 byte data
      this->datatype=2; //ENVI data type 
      break;
   case uint32:
      this->datasize=4; //4 byte data
      this->datatype=13; //ENVI data type 
      break;
   case int32:
      this->datasize=4; //4 byte data
      this->datatype=3; //ENVI data type 
      break;
   case float32:
      this->datasize=4; //4 byte data
      this->datatype=4; //ENVI data type 
      break;
   case float64:
      this->datasize=8; //8 byte data
      this->datatype=5; //ENVI data type 
      break;
   default:
      this->bilinfo<<"Unknown datatype for BIL file. Type given: "<<type<<std::endl;
      throw BILexception(bilinfo.str());
      break;
   }
}
