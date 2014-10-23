#include "multifile.h"

//-------------------------------------------------------------------------
// MultiFileReader Constructor - takes a string filename that is an ascii
// file which contains a list of BIL/BSQ filenames. Note that all BIL/BSQ
// files should have same dimensions and data type etc
//-------------------------------------------------------------------------
MultiFile::MultiFile(std::string filenamelist) : BinaryReader()
{
   //Filename is an ascii list of bil/bsq files of single band (this we will determine)
   br=NULL;
   //tmp storage for filename
   std::string tmpfilename;
   //Need to read the list of files
   std::ifstream fin;
   fin.open(filenamelist.c_str(),std::ifstream::in);
   if(!fin.is_open())
   {
      throw "Unable to open file (which should be an ascii list of files if I'm correct.): "+filenamelist;
   }
   else
   {
      while(fin>>tmpfilename)
      {
         //add the filename to the list
         filenames.push_back(tmpfilename);
      }
      //close the file
      fin.close();
      fin.clear();
   }
   numoffiles=filenames.size();
   numbands=numoffiles;
   br=new BinaryReader*[numoffiles];
   //This method requires all files to have number of bands =1
   //therefore bil is equivalent to bsq also
   int i=0;
   for(std::list<std::string>::iterator it=filenames.begin();it!=filenames.end();it++)
   {
      //Get the interleave type
      br[i]=new BinaryReader((*it));
      BinaryReader::interleavetype ft=br[i]->GetFileStyle();
      delete br[i];

      br[i]=NULL;
      if(ft==BinaryReader::BSQ)
         br[i]=new BSQReader((*it));
      else if(ft==BinaryReader::BIL)
         br[i]=new BILReader((*it));      
      else
         throw "Error. Iterleave type in file: "+(*it)+" is not bsq or bil.";

      i++;
   }

   CheckAllFilesConsistent();
}

//-------------------------------------------------------------------------
// Multifile destructor
//-------------------------------------------------------------------------
MultiFile::~MultiFile()
{
   if(br!=NULL)
   {
      for(unsigned int i=0;i<numoffiles;i++)
      {
         if(br[i]!=NULL)
            delete br[i];
      }
      delete[] br;
   }
}

//-------------------------------------------------------------------------
// Function to check the consistency of all the files to be 'combined'
// Test dimensions and data type and size are the same.
//-------------------------------------------------------------------------
bool MultiFile::CheckAllFilesConsistent()
{
   unsigned int tnumsamples=0,tnumlines=0;
   unsigned int tdatasize=0,tdatatype=0;
   for(unsigned int i=0;i<numoffiles;i++)
   {
      if(tnumsamples==0)
         tnumsamples=StringToUINT(br[i]->FromHeader("samples"));
      else if(tnumsamples!=StringToUINT(br[i]->FromHeader("samples")))
         throw "Number of samples is different between previous file and following file in MultiFileReader: "+br[i]->GetFileName();

      if(StringToUINT(br[i]->FromHeader("bands"))!=1)
         throw "Number of bands must be 1 in MultiFileReader, check file: "+br[i]->GetFileName();

      if(tnumlines==0)
         tnumlines=StringToUINT(br[i]->FromHeader("lines"));
      else if(tnumlines!=StringToUINT(br[i]->FromHeader("lines")))
         throw "Number of lines is different between previous file and following file in MultiFileReader:"+br[i]->GetFileName();

      if(tdatasize==0)
         tdatasize=br[i]->GetDataSize();
      else if(tdatasize!=br[i]->GetDataSize())
         throw "Data size (e.g. 8-bit, 16-bit etc) is different between previous file and following file in MultiFileReader:"+br[i]->GetFileName();

      if(tdatatype==0)
         tdatatype=br[i]->GetDataType();
      else if(tdatatype!=br[i]->GetDataType())
         throw "Data type (e.g. 2, 12, etc) is different between previous file and following file in MultiFileReader:"+br[i]->GetFileName();

   }
   numsamples=tnumsamples;
   numrows=tnumlines;
   bytesinline=tnumsamples*tdatasize;
   datasize=tdatasize;
   datatype=tdatatype;
   filesize=tnumsamples*tnumlines*tdatasize*numoffiles;
   return true;
}

//-------------------------------------------------------------------------
// The following functions are overloaded data reading functions
//-------------------------------------------------------------------------
void MultiFile::Readline(char* const chdata)
{
   for(unsigned int i=0;i<numoffiles;i++)
   {
      br[i]->Readline(&(chdata[i*bytesinline]));
   }
}

void MultiFile::Readline(char* const chdata, unsigned int line)
{
   for(unsigned int i=0;i<numoffiles;i++)
   {
      br[i]->Readline(&(chdata[i*bytesinline]),line);
   }
} 

void MultiFile::Readlines(char* const chdata, unsigned int startline, unsigned int numlines)
{
   for(unsigned int line=0;line<numlines;line++)
   {
      for(unsigned int i=0;i<numoffiles;i++)
      {
         br[i]->Readline(&(chdata[line*numoffiles*bytesinline + i*bytesinline]),startline+line);
      }
   }
}
 
int MultiFile::Readband(char* const chdata, unsigned int band)
{
   return br[band]->Readband(chdata,0);
} 

int MultiFile::Readbandline(char* const chdata, unsigned int band, unsigned int line)
{
   return br[band]->Readbandline(chdata,0,line);
}

double MultiFile::ReadCell(const unsigned int band,const unsigned int line, const unsigned int col)
{
   return br[band]->ReadCell(0,line,col);
}

void MultiFile::ReadlineToDoubles(double* const ddata,unsigned int line)
{
   for(unsigned int i=0;i<numoffiles;i++)
   {
      br[i]->ReadlineToDoubles(&(ddata[i*numsamples]),line);
   }
}

void MultiFile::Close()
{
   for(unsigned int i=0;i<numoffiles;i++)
   {
      br[i]->Close();
      br[i]=NULL;
   }
}

std::string MultiFile::FromHeader(std::string key,std::string THROW)
{
   std::string lowerkey=ToLowerCase(key);
   std::string stuff=br[0]->FromHeader(lowerkey,THROW);
   std::string tmpstuff="";
   //---------------------------------------------------------------------
   //Special cases
   //---------------------------------------------------------------------
   if(lowerkey.compare("bands")==0)
   {
      return ToString(numoffiles);
   }

   if((lowerkey.compare("wavelength")==0))
   {
      stuff="";
      for(unsigned int i=0;i<numoffiles;i++)
      {
         tmpstuff=br[i]->FromHeader(lowerkey,THROW);
         //Assuming since num bands should be 1 - therefore should only have 1 wavelength in the file
         stuff += RemoveAllOf(tmpstuff,"{}")+" ";
      }   
      return stuff;
   }

   //---------------------------------------------------------------------
   // Check all headers give the same value for the remaining keywords
   //---------------------------------------------------------------------
   for(unsigned int i=1;i<numoffiles;i++)
   {
      tmpstuff=br[i]->FromHeader(lowerkey,THROW);
      if(stuff.compare(tmpstuff)!=0)
      {
         throw "Header item is different in different files - wont allow this for now in function Fromheader(key)"
               " . If this is a mistake update the code with the special case for this keyword: "+key;
      }
   }
   return stuff;
}

std::string MultiFile::FromHeader(std::string key,int itemnum,std::string THROW)
{
   //Tricky function to write in a good way - will do some special cases and 
   //add to them when problems arise / new cases are required
   std::string lowerkey=ToLowerCase(key);
   if((lowerkey.compare("wavelength")==0)||(lowerkey.compare("fwhm")==0))
   {
      if((unsigned)itemnum<NumBands())
         return RemoveAllOf(br[itemnum]->FromHeader(key,THROW),"{};");
      else
         return "";
   }

   //Now check for items that are in each header (and should be the same for each one) but not
   //spread across the different headers, e.g. binning

   //---------------------------------------------------------------------
   // Check all headers give the same value for the remaining keywords
   //---------------------------------------------------------------------
   std::string stuff=br[0]->FromHeader(lowerkey,itemnum,THROW);
   std::string tmpstuff="";
   for(unsigned int i=1;i<numoffiles;i++)
   {
      tmpstuff=br[i]->FromHeader(lowerkey,itemnum,THROW);
      if(stuff.compare(tmpstuff)!=0)
      {
         throw "Header item is different in different files - wont allow this for now in function FromHeader(key,num)"
               " . If this is a mistake update the code with the special case for this keyword: "+key;
      }
   }
   return stuff;
}
