//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef SPECIMSENSORS_H
#define SPECIMSENSORS_H

#include "binfile.h"
#include "bilwriter.h"
#include "commonfunctions.h"
#include "logger.h"

//-------------------------------------------------------------------------
// Constants defined for sensors used in aplcal
// - maximum raw data values (Eagle, Hawk, Fenix based at software binning 1)
// - frame transfer time (Eagle CCD)
// - scalar to apply to radiance for storing as 16-bit
// - sensor id / serial numbers
//-------------------------------------------------------------------------
const unsigned short HAWK_RAW_MAX=16383;
const unsigned short EAGLE_RAW_MAX=4095;
const unsigned short FENIX_VNIR_RAW_MAX=4095;
const unsigned short FENIX_SWIR_RAW_MAX=65535;

const double FRAME_TRANSFER_TIME=0.002;

const unsigned int RADIANCE_DATA_SCALAR=1000;
const unsigned short CALIBRATED_DATA_MAX=std::numeric_limits<unsigned short>::max();

const int EAGLE_SENSOR_IDS_INTARR[2]={100022,110001};
const int HAWK_SENSOR_IDS_INTARR[2]={300011,310018};
const int FENIX_SENSOR_IDS_INTARR[2]={350005,350002};

//Constant vectors that contain the sensor id numbers for each sensor type
const std::vector<int> EAGLE_SENSOR_IDS(EAGLE_SENSOR_IDS_INTARR,EAGLE_SENSOR_IDS_INTARR + sizeof(EAGLE_SENSOR_IDS_INTARR) / sizeof(int));
const std::vector<int> HAWK_SENSOR_IDS(HAWK_SENSOR_IDS_INTARR,HAWK_SENSOR_IDS_INTARR + sizeof(HAWK_SENSOR_IDS_INTARR) / sizeof(int));
const std::vector<int> FENIX_SENSOR_IDS(FENIX_SENSOR_IDS_INTARR,FENIX_SENSOR_IDS_INTARR+sizeof(FENIX_SENSOR_IDS_INTARR)/sizeof(int));

//Defines for functions
enum SENSORTYPE {EAGLE,HAWK,FENIX};
bool CheckSensorID(SENSORTYPE s,int id);

//Enum a type to use in functions to specify which data to use
enum Subsensor {VNIR=0,SWIR=1};


//-------------------------------------------------------------------------
//This is an abstract data class so don't try and make an instance of it!
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Sensor object class
//-------------------------------------------------------------------------
class Sensor
{
public:
   Sensor(){}; //constructor
   virtual ~Sensor(){}; //destructor

   //string stream object to hold information on status
   std::stringstream ssinfo;

   //-------------------------------------------------------------------------
   //Adding a BIL reader to the class to make reading the raw data and 
   //giving the sensor access to the file dimensions and hdr values rather
   //than havng to pass them all over as variables between the two classes
   //-------------------------------------------------------------------------
   
protected:

};

//-------------------------------------------------------------------------
// Silly little class made for handling the qcfailures
//-------------------------------------------------------------------------
class Pair
{
public:
   Pair(unsigned int b,unsigned int s){band=b;sample=s;}
   ~Pair(){};
   unsigned int band;
   unsigned int sample;
};


//-------------------------------------------------------------------------
// SpecimBinFile class that inherits from BinFile to allow file reading
//-------------------------------------------------------------------------
class SpecimBinFile: public BinFile
{
public:
   SpecimBinFile(std::string strFilename);
   virtual ~SpecimBinFile();

   std::string TestFile(); //Test the file for possible errors etc

   virtual std::string GetFromFile(std::string keyword){Logger::Warning("No GetFromFile declared for this sensor - using FromHeader: "+keyword);return FromHeader(keyword);}
   virtual void SetSubSensor(const Subsensor sub){throw "Trying to set subsensor for a SpecimBinFile - did you mean to use a FenixBinFile?";}
   virtual void SetSubSensor(const unsigned int sub){throw "Trying to set subsensor for a SpecimBinFile - did you mean to use a FenixBinFile?";}

private:
   std::string sensorid;
};

//-------------------------------------------------------------------------
// EagleHawkBinFile class that inherits from SpecimBinFile to allow file reading
//-------------------------------------------------------------------------
class EagleHawkBinFile: public SpecimBinFile
{
public:
   EagleHawkBinFile(std::string strFilename);
   virtual ~EagleHawkBinFile(){};
   virtual std::string GetFromFile(std::string keyword);

private:
};

//-------------------------------------------------------------------------
// FenixBinFile class that inherits from SpecimBinFile to allow file reading
//-------------------------------------------------------------------------
class FenixBinFile: public SpecimBinFile
{
public:
   FenixBinFile(std::string strFilename);
   virtual ~FenixBinFile();

   virtual std::string GetFromFile(std::string keyword);
   virtual void SetSubSensor(const Subsensor sub);
   virtual void SetSubSensor(const unsigned int sub);

   unsigned int SubsensorLowerBand()const{return subsenlowerband;}
   unsigned int SubsensorUpperBand()const{return subsenupperband;}
   unsigned int SubsensorNumOfBands()const{return subsennumberofbands;}

   //Will need to overload some (all?) of the reading functions to take into account
   //the fact that the number of bands and first band are not the same as those in hdr file.
   virtual double ReadCell(const unsigned int band,const unsigned int line, const unsigned int col);
   virtual void ReadlineToDoubles(double* const ddata,unsigned int line);
   virtual void Readlines(char* const chdata, unsigned int startline, unsigned int numlines);

   //The following functions are not guaranteed to work with Fenix files and until they've been implemented throw an exception.
   virtual void Readline(char* const chdata){throw BinaryReader::BRexception("Undefined function call: ",__PRETTY_FUNCTION__);}
   virtual void Readline(char* const chdata, unsigned int line){throw BinaryReader::BRexception("Undefined function call: ",__PRETTY_FUNCTION__);}
   virtual void Readbytes(char* const chdata, unsigned long int bytes){throw BinaryReader::BRexception("Undefined function call: ",__PRETTY_FUNCTION__);}
   virtual int Readband(char* const chdata, unsigned int band){throw BinaryReader::BRexception("Undefined function call: ",__PRETTY_FUNCTION__);}
   virtual int Readbandline(char* const chdata, unsigned int band, unsigned int line){throw BinaryReader::BRexception("Undefined function call: ",__PRETTY_FUNCTION__);}

private:
   Subsensor subsensor;
   unsigned int subsenlowerband,subsenupperband,subsennumberofbands;
};

//-------------------------------------------------------------------------
// Class to describe the FODIS region
//-------------------------------------------------------------------------
class Fodis
{
public:
   Fodis(SpecimBinFile* bin);
   ~Fodis(){};
   unsigned long LowerFodis()const{return lowerfodis;}
   unsigned long UpperFodis()const{return upperfodis;}
   std::string FodisUnits()const{return fodisunits;}
   unsigned long RegionSize()const{return upperfodis-lowerfodis;}
private:
   unsigned int lowerfodis,upperfodis; // hold the lower and upper fodis bounds
   std::string fodisunits;
};

//-------------------------------------------------------------------------
// Specim object
//-------------------------------------------------------------------------
class Specim: public Sensor
{
public:
   //constructor with sensor id and calibration file prefix,and raw filename
   Specim(std::string strFilename,bool force=false);  
   //destructor
   virtual ~Specim(); 
   //Enum the various mask pixel values
   enum MaskType {Good=0, UnderFlow=1, OverFlow=2, Badpixel=4, SmearAffected=8,DroppedScan=16,CorruptData=32,QCFailure=64};

   SpecimBinFile* bin;
   Fodis* fodis;

   virtual unsigned int TotalNumBands()const{return NumBands();}
   virtual unsigned int NumBands()const{return numbands;}
   unsigned int NumLines()const{return numlines;}
   unsigned int NumSamples()const{return numsamps;}
   int SensorID() const {return SENSOR_ID;}
   std::string CalibratedUnits()const{return calibratedunits;}
   std::string RawFilename() const {return strRawFilename;}

   //List of band,sample that failed a quality control check and should be added to the mask
   std::vector<Pair> qcfailures;

   double IntegrationTime() const {return tint;}
   virtual unsigned int SpectralBinning()const{return spectralbinning;}
   virtual unsigned int SpatialBinning()const{return spatialbinning;}
   unsigned int LowerScanlineLimit() const{return scanlinelowerlimit;}
   unsigned int UpperScanlineLimit() const{return scanlineupperlimit;}

   unsigned int RadianceScalar() const {return radscalar;}
   unsigned short int RawMax() const {return rawmax;}
   unsigned short int CalibratedMax() const {return calibratedmax;}

   unsigned int GetNumImageFrames() {return StringToUINT(bin->FromHeader("lines")) - ndarklines;}
   short int GetTotalMissingFrames() const {return totalmissing;}
   unsigned int GetNumDarkFrames() const {return ndarklines;}
   double DarkScalar() const {return darkscalar;}

   void ReadInAllDarkFrames(double* dlstorage,std::string externalFileName,uint64_t linecellsize=0,unsigned int subsensor=0);
   void AverageAllDarkFrames(double* const data,std::string externalFileName="",uint64_t linecellsize=0,unsigned int subsensor=0);
   void DarkFramesStdDeviation(double* const stdev,const double* const mean, std::string externalFileName="",uint64_t linecellsize=0,unsigned int subsensor=0);
   void AverageRefinedDarkFrames(double* const data,const double* const stdev,const double* const mean, std::string externalFileName="",uint64_t linecellsize=0,unsigned int subsensor=0);
   unsigned short GetMissingFramesBetweenLimits(unsigned int start,unsigned int end){return TotalMissingFrames(start,end);}

   //Reads the first/last frame IDs and compare to number of lines to get total number of missing frames
   virtual void TotalMissingFrames(){throw "Undefined function: TotalMissingFrames";}
   //Reads the start/end frame IDs and compare to number of lines within these limits to get number of missing frames
   virtual unsigned short TotalMissingFrames(const unsigned int start, const unsigned int end);
   //Function to return the lower band in the band range - normalised to 0.
   unsigned int LowerBandLimit()const{return StringToUINT(bin->GetFromFile("lowervimg"))-1;}

protected:

   void DarkFrameSanityCheck();

   //-------------------------------------------------------------------------
   //Constants relating to the Specim sensors
   //-------------------------------------------------------------------------
   unsigned short int rawmax; //maximum value
   unsigned short int calibratedmax; //maximum value of a uint16
   unsigned int radscalar; //used to multiply the radiance so can store a more accurate value as uint16
   int SENSOR_ID; //id string for the Eagle sensor (serial no.)
   std::string calibratedunits;
   const static unsigned int MAXFRAMECOUNT=65535;
   //Raw data file name  - used to derive other filenames
   std::string strRawFilename;

   //-------------------------------------------------------------------------
   //General Specim sensor variables
   //-------------------------------------------------------------------------
   unsigned int prevfnum; //previous frame number - unsigned int (not short) so can have a flag value
   double tint; //integration time of sensor
   unsigned int spatialbinning,spectralbinning; // spatial and spectral binning values
   unsigned int numbands, numsamps, numlines;   //number of bands, samples and lines relating to the image data array
   unsigned int scanlineupperlimit, scanlinelowerlimit; //himg {} values from specim header

   short int totalmissing;//total number of missing frames
   unsigned int ndarklines;
   unsigned int darklinestart;
   double darkscalar; //scalar if separate dark file used of different integration time
   bool DARKFORCE;
   //-------------------------------------------------------------------------
   // Not that the below function should NOT be called from this class - only
   // from derived ones - but it is the same for all derived classes so it's
   // included here rather than duplicating code.
   //-------------------------------------------------------------------------
   void GetExtraInfoFromHeader();
};


//-------------------------------------------------------------------------
// Eagle object
//-------------------------------------------------------------------------
class Eagle: public Specim
{
public:
   Eagle(std::string strFilename,bool force=false);
   ~Eagle();

   double FrameTransferTime() const {return trant;}

   //Reads the first/last frame IDs and compare to number of lines to get total number of missing frames
   virtual void TotalMissingFrames(); 

protected:
   double trant;//frame transfer time for Eagle in ms - from CCD document (recommended operating rate - not necessarily true rate!)
};

//-------------------------------------------------------------------------
// Hawk object
//-------------------------------------------------------------------------
class Hawk: public Specim
{
public:
   Hawk(std::string strFilename,bool force=false);
   ~Hawk();

   //Reads the first/last frame IDs and compare to number of lines to get total number of missing frames
   virtual void TotalMissingFrames(); 
protected:

};

//-------------------------------------------------------------------------
// Fenix object
//-------------------------------------------------------------------------
class Fenix: public Specim
{
public:
   Fenix(std::string strFilename);
   ~Fenix();

   //Function to set the Fenix object up for reading / calibrating either swir or vnir
   void SetUpFenixFor(Subsensor sen);
   //Trick the NumBands into returning only the number for the subsensor we're working with
   virtual unsigned int NumBands()const{return dynamic_cast<FenixBinFile*>(bin)->SubsensorNumOfBands();}
   virtual unsigned int TotalNumBands()const{return NumBandsVnir()+NumBandsSwir();}
   unsigned int NumBandsVnir() const {return vnbr->SubsensorNumOfBands();}
   unsigned int NumBandsSwir() const {return swbr->SubsensorNumOfBands();}

   virtual unsigned int SpectralBinning()const{return StringToUINT(TrimWhitespace(bin->GetFromFile("spectralbinning").c_str()));}
   virtual unsigned int SpatialBinning()const{return StringToUINT(TrimWhitespace(bin->GetFromFile("spatialbinning").c_str()));}

   unsigned int SubSenLowerBand()const{return dynamic_cast<FenixBinFile*>(bin)->SubsensorLowerBand();}
   unsigned int SubSenUpperBand()const{return dynamic_cast<FenixBinFile*>(bin)->SubsensorLowerBand();}

   //Reads the first/last frame IDs and compare to number of lines to get total number of missing frames
   virtual void TotalMissingFrames(); 
protected:
   //Two readers - one for each CCD - only one of which is assigned to br at a time
   //based on the SetUpFenixFor function call
   FenixBinFile* vnbr;
   FenixBinFile* swbr;

};



#endif
