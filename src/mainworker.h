//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef MAINWORKER_H
#define MAINWORKER_H

#include "specimsensors.h"
#include "calibration.h"


class MainWorker
{
public:

   MainWorker(std::string rawfile,std::string outfile,std::string cl="",bool force=false);
   MainWorker(std::string rawfile,std::string outfile,char sensor,std::string cl="",bool force=false);
   ~MainWorker();

   Specim* sensor;
   Calibration* cal;

   enum OutputDataFlag {Normal,MissingScan,CorruptData};
   enum Task {remove_dark_frames,apply_gains,calibrate_fodis,insert_missing_scans,smear_correct,output_mask,flip_bands,flip_samples,output_mask_method,apply_qcfailures};

   void WriteOutData(OutputDataFlag flag);
   void InitialiseCalibration(std::string calfile,std::string externaldarkframes,std::string qcfailurefile);
   unsigned int GetNumCalibratedImageLines();
   unsigned int GetNumCalibratedImageSamples();
   void DoCalibrationForLine(unsigned int line);

   bool GetTask(Task o){return tasks[o];}

   void SetTask(Task t, bool b){tasks[t]=b;}
   void SetSampleLimits(unsigned int l, unsigned int u){lowersample=l;uppersample=u;}
   void SetLineLimits(unsigned int l, unsigned int u){startline=l;endline=u;}
   void SetDroppedScansPriorToStartLine(unsigned int d){nummissingscanspriortostartline=d;}
   std::string TasksAsString();

private:
   BILWriter* bwmaskmethod;   
   BILWriter* bwmask;
   BILWriter* bwimage;
   BILWriter* bwfodis;

   Eagle* eagle;
   Hawk* hawk;

   std::string outputfileprefix;
   void InitialiseWriters();
   std::string ReverseWavelengthOrder(std::string wavelengths);
   void TransferHeaderInfo(BILWriter* const bw);

   std::string commandlinetext;

   unsigned int numbercalibratedframes;
   unsigned int lowersample,uppersample;
   unsigned int startline,endline,nummissingscanspriortostartline;
   std::map<int,bool> tasks;
};

#endif
