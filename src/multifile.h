#ifndef MULTIFILE_H
#define MULTIFILE_H

#include <list>
#include <string>

#include "binaryreader.h"
#include "bil.h"
#include "bsq.h"

class MultiFile : public BinaryReader
{
public:
   MultiFile(std::string filenamelist);
   ~MultiFile();

   virtual void Readline(char* const chdata); 
   virtual void Readline(char* const chdata, unsigned int line);
   virtual void Readlines(char* const chdata, unsigned int startline, unsigned int numlines); 
   virtual void Readbytes(char* const chdata, unsigned long int bytes)
               {throw BinaryReader::BRexception("Undefined MultiFile function call: ",__PRETTY_FUNCTION__);} 
   virtual int Readband(char* const chdata, unsigned int band); 
   virtual int Readbandline(char* const chdata, unsigned int band, unsigned int line); 
   virtual double ReadCell(const unsigned int band,const unsigned int line, const unsigned int col);
   virtual void ReadlineToDoubles(double* const ddata,unsigned int line);
   virtual int ReadRect(char* const chdata, const int minrow, const int maxrow,const int mincol,const int maxcol)
               {throw BinaryReader::BRexception("Undefined MultiFile function call: ",__PRETTY_FUNCTION__);} 
   virtual void Close();
   virtual std::string FromHeader(std::string key,std::string THROW="false");
   virtual std::string FromHeader(std::string key,int itemnum,std::string THROW="false");
   virtual std::map<std::string, std::string, cmpstr> CopyHeader()const{throw "Cannot call CopyHeader for a MultiFile";}

private:   
   virtual void ReadPartOfLine(char* const chdata,const unsigned int lineno,const unsigned int sampleno,const unsigned int nsamps)
               {throw BinaryReader::BRexception("Undefined MultiFile function call: ",__PRETTY_FUNCTION__);}

   bool CheckAllFilesConsistent();
   BinaryReader** br;
   std::list<std::string> filenames;
   unsigned int numoffiles;
   uint64_t bytesinline;

};

#endif
