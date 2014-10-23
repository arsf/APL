#ifndef FILEWRITER_H
#define FILEWRITER_H

#include <string>

//Base class for writing to different file types
//Used as an interface to the classes that contain implementation
class FileWriter
{
public:
   FileWriter(std::string fname){filename=fname;datasize=0;datatype=0;}
   FileWriter(){}
   virtual ~FileWriter(){};
   virtual int WriteLine(char* const data)=0;
   virtual int WriteBandLine(char* const data)=0;
   virtual int WriteBandLineSection(char* const data,const unsigned int numsamples_array,const unsigned int start, const unsigned int end)=0;
   virtual int WriteBandLineWithValue(const char xval)=0; 

   virtual int Close()=0;
   virtual void AddToHdr(std::string item)
      {throw FileWriterException("Undefined function call: ",__PRETTY_FUNCTION__);}

   virtual unsigned int GetDataSize() const =0;
   virtual unsigned int GetDataType() const =0;

   enum DataType {uchar8, char8, uint16, int16, uint32, int32, float32, float64 };

   //Exception class
   class FileWriterException
   {
   public:
      std::string info;
      FileWriterException();
      FileWriterException(std::string ss){info=ss;}
      FileWriterException(std::string ss,const char* extra){info=ss+std::string(extra);}

      const char* what() const throw()
      {
         return "A FileWriterException has occurred.";
      }
   };

   virtual void AddMetadata(std::string name,std::string value)
      {throw FileWriterException("Undefined function call: ",__PRETTY_FUNCTION__);}

protected:
   std::string filename;
   unsigned int datasize,datatype;
};


#endif


