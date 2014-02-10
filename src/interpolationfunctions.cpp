//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#include "interpolationfunctions.h"

#ifdef DEBUGINTERPFUNC
   #define DEBUGPRINT(x) std::cout<<x<<std::endl;
#else
   #define DEBUGPRINT(x)
#endif

enum {LATITUDE,LONGITUDE,HEIGHT,ROLL,PITCH,HEADING};
const double HEADING_DISCONTINUITY_CHECK=350;

//Functions only available to others within this cpp
//Function to generate the kernel for the Smooth function
void TriangleKernel(double* const kernel,const int length);
//Functions related to the Spline interpolation method
unsigned long GetSecondDerivatives(double start, double stop, DataHandler* const dhandle,const short dataflag,double* const derivatives);
double GetSplineResult(const double t,DataHandler* const dhandle,const short dataflag,const double* const derivatives,unsigned long startpoint);
unsigned long GetSecondDerivativesWrapper(std::string start, std::string stop, DataHandler* const dhandle,const short dataflag,double* const derivatives);
unsigned long GetSecondDerivativesWrapper(double start, double stop, DataHandler* const dhandle,const short dataflag,double* const derivatives);

//-------------------------------------------------------------------------
// Interpolation functions below followed by smoothing functions
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Function to do a linear interpolation to the navigation data using 
// only the entries directly either side of the interpolated point
//-------------------------------------------------------------------------
void Linear(const double* const times,const int len, DataHandler* const dhandle,NavDataCollection* store,std::string start=NULL,std::string stop=NULL)
{

   //Need to find the entries either side of the given time to use for interpolation
   //get the first time value (assume this is the minimum time)
   unsigned long l=0; //line counter
   NavDataLine* tmpLine=NULL;
   double t=dhandle->GetLine(l)->time;

   for(int ti=0;ti<len;ti++)
   {     
      //Check if time is after first navigation data point
      if(times[ti] < t)   
      {
         throw "Error - The given time falls before the navigation data in Linear(): "+ToString(times[ti]);
      }

      while((t < times[ti])&&(l < dhandle->GetNumEntries()))
      {
         l++;
         tmpLine=dhandle->GetLine(l);
         if(tmpLine==NULL)
         {
            throw "Navigation file does not contain enough data to cover the flight line, assuming start time is correct.";
         }
         t=tmpLine->time;
         //std::cout<<l<<" "<<t<<" "<<times[ti]<<" "<<ti<<" "<<dhandle->GetNumEntries()<<" extra:"<<tmpLine->hei<<std::endl;
      }

      //Check if time is before last navigation data point
      if (l >= dhandle->GetNumEntries())
      {
         //The given time does not fall within the nav data
         throw "Error - The given time falls after the navigation data in Linear(): Item "+ToString(ti)+" : Time "+ToString(times[ti]) 
               + "\n This suggests that the navigation data does not cover the entirety of the flight line.";

         //return NULL;
      }

      //Added a special case for if time to interpolate to is == time of nav data
      //and nav line is the first one (==0)
      NavDataLine* before=NULL;
      NavDataLine* after=NULL;
      if((times[ti] == t)&& (l==0))
      {
         before=dhandle->GetLine(l);
         after=dhandle->GetLine(l+1);
      }
      else
      {
         //We want data from line l and l-1
         before=dhandle->GetLine(l-1);
         after=dhandle->GetLine(l);
      }

      //use the "distance" weighting to interpolate the data
      double timespan=after->time - before->time;
      double interppoint=times[ti] - before->time;
      double scalar=interppoint / timespan;
     
      //store->SetValues(ti,NavDataCollection::TIME,times[ti]); //Do not need to store times - they are overwritten later anyway
      store->SetValues(ti,NavDataCollection::LAT,before->lat + (after->lat - before->lat)*scalar);
      store->SetValues(ti,NavDataCollection::LON,before->lon + (after->lon - before->lon)*scalar);
      store->SetValues(ti,NavDataCollection::HEI,before->hei + (after->hei - before->hei)*scalar);
      store->SetValues(ti,NavDataCollection::ROLL,before->roll + (after->roll - before->roll)*scalar);
      store->SetValues(ti,NavDataCollection::PITCH,before->pitch + (after->pitch - before->pitch)*scalar);

      //Need to do something different with heading since wrap around (0-360)
      if((after->heading - before->heading)<-HEADING_DISCONTINUITY_CHECK)
      {
         //The heading has wrapped around clockwise
         store->SetValues(ti,NavDataCollection::HEADING,before->heading + ((after->heading+360) - before->heading)*scalar);
         //Now correct the value if it is larger than 360
         if(store->GetValue(ti,NavDataCollection::HEADING)>360)
            store->SetValues(ti,NavDataCollection::HEADING,store->GetValue(ti,NavDataCollection::HEADING)-360);
      }
      else if((after->heading - before->heading)>HEADING_DISCONTINUITY_CHECK)
      {
         //The heading has wrapped around counter-clockwise
         store->SetValues(ti,NavDataCollection::HEADING,before->heading + ((after->heading-360) - before->heading)*scalar);
         //Now correct the value if it is less than 0
         if(store->GetValue(ti,NavDataCollection::HEADING)<0)
            store->SetValues(ti,NavDataCollection::HEADING,store->GetValue(ti,NavDataCollection::HEADING)+360);
      }
      else
      {
         //The heading has not wrapped around
         store->SetValues(ti,NavDataCollection::HEADING,before->heading + (after->heading - before->heading)*scalar);
      }
      //DEBUGPRINT("Line values: "<<l-1<<" "<<l<<" Time values: Before "<<before->time<<" To Get "<<times[ti]<<" After "<<after->time)

      //Restore the variables for a new search - dont assume time array is in time order
      l=0; //line counter
      t=dhandle->GetLine(l)->time;

   }
}

//-------------------------------------------------------------------------
// Cubic Spline Wrapper function: to call GetSecondDerivatives using string times
//-------------------------------------------------------------------------
unsigned long GetSecondDerivativesWrapper(std::string start, std::string stop, DataHandler* const dhandle,const short dataflag,double* const derivatives)
{
   //Convert time strings to second of day
   double starttime=0,stoptime=0;
   starttime=ConvertTime(start);
   stoptime=ConvertTime(stop);
   DEBUGPRINT("Start time: "<<start<<" as seconds: "<<starttime<<" Stop time: "<<stop<<" as seconds: "<<stoptime)
   return GetSecondDerivatives(starttime,stoptime,dhandle,dataflag,derivatives);
}

//-------------------------------------------------------------------------
// Cubic Spline Wrapper function: to call GetSecondDerivatives using double times
//-------------------------------------------------------------------------
unsigned long GetSecondDerivativesWrapper(double start, double stop, DataHandler* const dhandle,const short dataflag,double* const derivatives)
{
   //make sure in sec of day
   double starttime=(int)start%(24*3600)+start-(int)start;;
   double stoptime=(int)stop%(24*3600)+stop-(int)stop;;
   return GetSecondDerivatives(starttime,stoptime,dhandle,dataflag,derivatives);
}


//-------------------------------------------------------------------------
// Cubic Spline function: to calculate the second derivatives of the data function
//-------------------------------------------------------------------------
unsigned long GetSecondDerivatives(double starttime, double stoptime, DataHandler* const dhandle,const short dataflag,double* const derivatives)
{
   //We want to solve the matrix equation Ax=y where 
   //     [4 1 0 ..0 0 0                              [y1 - 2y2 +y3
   // A =  1 4 1 ..0 0 0   x=[second derivatives]  y=  y2 - 2y3 +y4
   //      0 1 4 ..0 0 0                                .......
   //      .............
   //      0 0 0 ..0 1 4]                              yn-2 - 2yn-1 + yn]

   //We shall use only the times [start:stop] rather than the full set of points

   //Get the number of data points in the entire file
   //const long nfull=dhandle->GetNumEntries();

   //Get the number of data points from [start:stop]
   double testtime=0;
   unsigned long l=0;
   while(testtime < starttime)
   {
      testtime=dhandle->GetLine(l)->time;
      testtime=(int)testtime%(24*3600)+testtime-(int)testtime; //make sure in sec of day
      l++;
   }
   double tmph=testtime;
   unsigned long startpoint=0; //line index of nav data before the data we want
   //Test if start is before first point and set it's value
   if(l==0)
   {
      //start time is prior to first data point of nav
      Logger::Warning("Start time of image data is before the start time of the navigation data.");      
      startpoint=l;
   }
   else
   {
      startpoint=l-1;
   }

   while(testtime < stoptime)
   {
      testtime=dhandle->GetLine(l)->time;
      testtime=(int)testtime%(24*3600)+testtime-(int)testtime; //make sure in sec of day
      l++;
   }
   unsigned long stoppoint=l; //line index of nav data after data we want
   const long n=stoppoint-startpoint+1; //number of data points we want to use
   double h=(testtime-tmph)/n;   //assume regularly spaced data

   DEBUGPRINT("First point at: "<<startpoint<<" Second point at: "<<stoppoint<<" Number of points "<<n<<" Time separation between consecutive points: "<<h)

   //Set the end derivatives to 0 (boundary condition for natural splines)
   derivatives[0]=0;
   derivatives[n-1]=0;

   //Set up the data arrays
   //Y:
   double* y=new double[n];
   double* a=new double[n];
   double* b=new double[n];
   double* c=new double[n];
   double* d=new double[n];
   NavDataLine* ndl1=NULL;
   NavDataLine* ndl2=NULL;
   NavDataLine* ndl3=NULL;
   for(int i=1;i<n-1;i++)
   {
      ndl1=dhandle->GetLine(startpoint+i-1);
      ndl2=dhandle->GetLine(startpoint+i);
      ndl3=dhandle->GetLine(startpoint+i+1);
      switch(dataflag)
      {
         case LATITUDE: 
            y[i]=(6/(h*h))*(ndl1->lat - 2*ndl2->lat + ndl3->lat);
            break;
         case LONGITUDE: 
            y[i]=(6/(h*h))*(ndl1->lon - 2*ndl2->lon + ndl3->lon);
            break;
         case HEIGHT: 
            y[i]=(6/(h*h))*(ndl1->hei - 2*ndl2->hei + ndl3->hei);
            break;
         case ROLL:
            y[i]=(6/(h*h))*(ndl1->roll - 2*ndl2->roll + ndl3->roll);
            break;
         case PITCH: 
            y[i]=(6/(h*h))*(ndl1->pitch - 2*ndl2->pitch + ndl3->pitch);
            break;
         case HEADING:
            //Heading needs to be handled differently due to wrap around 360->0 
            if((((ndl2->heading-ndl1->heading)>HEADING_DISCONTINUITY_CHECK)||((ndl2->heading-ndl1->heading)<-HEADING_DISCONTINUITY_CHECK))&&(((ndl3->heading-ndl2->heading)>HEADING_DISCONTINUITY_CHECK)||((ndl3->heading-ndl2->heading)<-HEADING_DISCONTINUITY_CHECK)))
            {
               //There is a double crossover of the 0/360 discontinuity
               //e.g. 359,0,359 or 0,359,0 
               if((ndl2->heading-ndl1->heading)>HEADING_DISCONTINUITY_CHECK)
               {
                  //This means it is 0,359,0 so subtract 360 from ndl2
                  y[i]=(6/(h*h))*(ndl1->heading - 2*(ndl2->heading-360) + ndl3->heading);                 
               }
               else
               {
                  //This means it is 359,0,359 so add 360 onto ndl2
                  y[i]=(6/(h*h))*(ndl1->heading - 2*(ndl2->heading+360) + ndl3->heading);                 
               }
            }
            else if(((ndl2->heading-ndl1->heading)>HEADING_DISCONTINUITY_CHECK)||((ndl2->heading-ndl1->heading)<-HEADING_DISCONTINUITY_CHECK))
            {
               //There is a single crossover and its between ndl2 and ndl1
               if(ndl2->heading>HEADING_DISCONTINUITY_CHECK) //Then ndl1 is around 0 or slightly higher so add 360 onto this (since ndl3 ~359 too) - we need all 3 to be of similar magn
                  y[i]=(6/(h*h))*(ndl1->heading+360 - 2*ndl2->heading + ndl3->heading); 
               else //ndl2 is around 0 so subtract 360 from ndl1 (since ndl3 will be ~0 too) - we need all 3 to be of similar magn
                  y[i]=(6/(h*h))*(ndl1->heading-360 - 2*ndl2->heading + ndl3->heading); 
            }
            else if(((ndl3->heading-ndl2->heading)>HEADING_DISCONTINUITY_CHECK)||((ndl3->heading-ndl2->heading)<-HEADING_DISCONTINUITY_CHECK))
            {
               //There is a single crossover and its between ndl3 and ndl2
               if(ndl3->heading>HEADING_DISCONTINUITY_CHECK) //Then ndl2 is around 0 or slightly higher (as will be ndl1) so subtract 360 from ndl3
                  y[i]=(6/(h*h))*(ndl1->heading - 2*ndl2->heading + ndl3->heading-360); 
               else //ndl3 is around 0 so add 360 onto that value since ndl2 (and ndl1) are around 359
                  y[i]=(6/(h*h))*(ndl1->heading - 2*ndl2->heading + ndl3->heading+360); 
            }
            else 
            {
               //No headings cross the 0/360 discontinuity - dont need to add anything
               y[i]=(6/(h*h))*(ndl1->heading - 2*ndl2->heading + ndl3->heading);
            }
            break;
         default:
            throw "Unrecognised dataflag in CubicSplines";
            break;
      }
   }

   //a,b:
   for(int i=0;i<n;i++)
   {
      a[i]=1;
      b[i]=4;
   }
   //set a[0] to be 0
   a[0]=0;

   //c:
   c[0]=1/b[0];
   d[0]=y[0]/b[0];
   for(int i=1;i<n;i++)
   {
      c[i]=1/(b[i]-c[i-1]*a[i]);
      d[i]=(y[i]-d[i-1]*a[1])/(b[i]-c[i-1]*a[i]);
   }

   derivatives[n-1]=y[n-1];
   for(int i=n-2;i>0;i--)
   {
      derivatives[i]=y[i]-c[i]*derivatives[i+1];
   }
  
   delete[] y;
   delete[] a;
   delete[] b;
   delete[] c;
   delete[] d;

   return startpoint;
}

//-------------------------------------------------------------------------
// Cubic Spline function: to calculate the interpolated value for time t
//-------------------------------------------------------------------------
void CubicSpline(const double* const times,const int len,DataHandler* const dhandle,NavDataCollection* store,std::string start,std::string stop)
{
   //Set the arrays for the 2nd derivatives to have the correct size
   double* lat_2der=new double[dhandle->GetNumEntries()];
   double* lon_2der=new double[dhandle->GetNumEntries()];   
   double* hei_2der=new double[dhandle->GetNumEntries()];
   double* roll_2der=new double[dhandle->GetNumEntries()];
   double* pitch_2der=new double[dhandle->GetNumEntries()];
   double* heading_2der=new double[dhandle->GetNumEntries()];
   
   //Get the second derivatives to use in the spline calculation
   //Call this function once only and store the result for each call for spline interpolation
   //and use it as a look up table for the 2nd derivatives   

   //Unsure as to why using the string start/stop times from header and not from the times array
   //assumming this was done for a reason have written a wrapper function that accepts either string  
   //or double times - but am using the double times now. If this goes wrong horribly then change 
   //mystart / mystop back to start and stop (although these were causing problems hence the change)
   double mystart=times[0];
   double mystop=times[len-1];

   unsigned long latfirstpoint=GetSecondDerivativesWrapper(mystart,mystop,dhandle,LATITUDE,lat_2der);
   unsigned long lonfirstpoint=GetSecondDerivativesWrapper(mystart,mystop,dhandle,LONGITUDE,lon_2der);
   unsigned long heifirstpoint=GetSecondDerivativesWrapper(mystart,mystop,dhandle,HEIGHT,hei_2der);
   unsigned long rollfirstpoint=GetSecondDerivativesWrapper(mystart,mystop,dhandle,ROLL,roll_2der);
   unsigned long pitchfirstpoint=GetSecondDerivativesWrapper(mystart,mystop,dhandle,PITCH,pitch_2der);
   unsigned long headingfirstpoint=GetSecondDerivativesWrapper(mystart,mystop,dhandle,HEADING,heading_2der);

   //Calculate the value of y at x
   for(int t=0;t<len;t++)
   {  
      store->SetValues(t,NavDataCollection::LAT,GetSplineResult(times[t],dhandle,LATITUDE,lat_2der,latfirstpoint));
      store->SetValues(t,NavDataCollection::LON,GetSplineResult(times[t],dhandle,LONGITUDE,lon_2der,lonfirstpoint));
      store->SetValues(t,NavDataCollection::HEI,GetSplineResult(times[t],dhandle,HEIGHT,hei_2der,heifirstpoint));
      store->SetValues(t,NavDataCollection::ROLL,GetSplineResult(times[t],dhandle,ROLL,roll_2der,rollfirstpoint));
      store->SetValues(t,NavDataCollection::PITCH,GetSplineResult(times[t],dhandle,PITCH,pitch_2der,pitchfirstpoint));
      store->SetValues(t,NavDataCollection::HEADING,GetSplineResult(times[t],dhandle,HEADING,heading_2der,headingfirstpoint));

      DEBUGPRINT(t<<" "<<times[t]<<" "<<store->GetValue(t,NavDataCollection::LAT)
                  <<" "<<store->GetValue(t,NavDataCollection::LON)<<" "<<store->GetValue(t,NavDataCollection::HEI)
                  <<" "<<store->GetValue(t,NavDataCollection::ROLL)<<" "<<store->GetValue(t,NavDataCollection::PITCH)
                  <<" "<<store->GetValue(t,NavDataCollection::HEADING))
   }

   //Free up the memory
   delete[] lat_2der;
   delete[] lon_2der;
   delete[] hei_2der;
   delete[] roll_2der;
   delete[] pitch_2der;
   delete[] heading_2der;
}


double GetSplineResult(const double t,DataHandler* const dhandle,const short dataflag,const double* const derivatives,unsigned long startpoint)
{
   //Y=A*yi + B*yi+1 +C*y''i +D*y''i+1
   static unsigned long l=startpoint; //made this static - can't think of any reason why l would need to go 'backwards'
   double testtime=dhandle->GetLine(l)->time;
   double testtime_prev=0;

   //Added these tests incase the static l above causes problems
   if(l!=0)
      testtime_prev=dhandle->GetLine(l-1)->time;
   else
      testtime_prev=testtime;

   //if the time from the nav file (testtime) is greater than one time point from the time we're looking for (t), l will equal the start of the nav file
   //actually compare against testtime_prev as this is the previous time point in the file
   if(testtime_prev > t)
   {
      l=startpoint;
      testtime=dhandle->GetLine(l)->time;
   }

   while((testtime < t)&&(l < dhandle->GetNumEntries()))
   {
      l++;
      testtime=dhandle->GetLine(l)->time;
      //testtime=(int)testtime%(24*3600)+testtime-(int)testtime; //make sure in sec of day
   }

   if(l == dhandle->GetNumEntries())
   {
      //Then we exited the above loop before we found the correct time
      //as we reached the end of the data in the nav array - this probably
      //means that the wrong nav data has been given for this file or it is not long enough etc
      throw "Time of scan line to interpolate to does not fall within navigation data: "+ToString(t);
   }

   unsigned long lhigh=l;
   unsigned long llow=l-1;

   double tlow=dhandle->GetLine(llow)->time;
   double thigh=dhandle->GetLine(lhigh)->time;

   //DEBUGPRINT("Line values: "<<llow<<" "<<lhigh<<" Time values: Before "<<tlow<<" To Get "<<t<<" After "<<thigh)

   double A=(thigh-t)/(thigh-tlow);
   double B=1-A;
   double C=(A*A*A-A)*(thigh-tlow)*(thigh-tlow)/6.0;
   double D=(B*B*B-B)*(thigh-tlow)*(thigh-tlow)/6.0;

   double ylow=0;
   double yhigh=0;

   switch(dataflag)
   {
      case LATITUDE: 
         ylow=dhandle->GetLine(llow)->lat;
         yhigh=dhandle->GetLine(lhigh)->lat;
         break;
      case LONGITUDE: 
         ylow=dhandle->GetLine(llow)->lon;
         yhigh=dhandle->GetLine(lhigh)->lon;
         break;
      case HEIGHT: 
         ylow=dhandle->GetLine(llow)->hei;
         yhigh=dhandle->GetLine(lhigh)->hei;
         break;
      case ROLL:
         ylow=dhandle->GetLine(llow)->roll;
         yhigh=dhandle->GetLine(lhigh)->roll;
         break;
      case PITCH: 
         ylow=dhandle->GetLine(llow)->pitch;
         yhigh=dhandle->GetLine(lhigh)->pitch;
         break;
      case HEADING: 
         ylow=dhandle->GetLine(llow)->heading;
         yhigh=dhandle->GetLine(lhigh)->heading;
         //Need to treat heading differently due to wrap around discontinuity at 0/360
         if(ylow-yhigh>HEADING_DISCONTINUITY_CHECK)
         {
            //ylow is ~359 and yhigh is ~0 - subtract 360 from ylow
            double tempvalue=A*(ylow-360)+B*yhigh+C*derivatives[l-startpoint]+D*derivatives[l-startpoint+1];
            if(tempvalue<0)
               tempvalue+=360;
            return tempvalue;//return here as a special case
         }
         else if(ylow-yhigh<-HEADING_DISCONTINUITY_CHECK)
         {
            //yhigh is ~359 and ylow ~0 - subtract 360 from yhigh
            double tempvalue=A*ylow+B*(yhigh-360)+C*derivatives[l-startpoint]+D*derivatives[l-startpoint+1];
            if(tempvalue<0)
               tempvalue+=360;
            return tempvalue;//return here as a special case  
         }
         //Else just break out of the loop and return like normal below
         break;
      default:
         throw "Unrecognised dataflag in GetSplineResult";
         break;
   }

   //DEBUGPRINT(dataflag<<" "<<A<<" "<<ylow<<" "<<B<<" "<<yhigh<<" "<<C<<" "<<derivatives[l-startpoint]<<" "<<D<<" "<<derivatives[l-startpoint+1]<<" "<<l<<" "<<startpoint)
   return A*ylow+B*yhigh+C*derivatives[l-startpoint]+D*derivatives[l-startpoint+1];
   
}



//-------------------------------------------------------------------------
// Smoothing functions below
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
// Function used to smooth the raw navigation to remove any jumps in the data
// of element, using a triangular kernel of size kernelsize
//-------------------------------------------------------------------------
void Triangle(const unsigned long element,DataHandler* const dhandle,NavDataLine* store,const int kernelsize)
{   
   //Do some checks on input data
   if(kernelsize%2 == 0)
      throw "Kernel size in Smooth function should be an odd number.";

   double* kernel=new double[kernelsize];
   int halflen=(kernelsize-1)/2;

   //Generate the kernel to smooth by
   TriangleKernel(kernel,kernelsize);

   NavDataLine** data;
   data=new NavDataLine*[kernelsize];

   //create a variable to hold heading wrap status
   // 0 - no wrap around
   // 1 - counter-clockwise wrap
   // 2 - clockwise wrap
   short int headingwrap=0;
   //Get the data
   for(int i=0;i<kernelsize;i++)
   {
      data[i]=dhandle->GetLine(element+i-halflen);
      //Check for heading wrap arounds
      if(i>=1)
      {     
         if(((data[i]->hei-data[i-1]->hei)>HEADING_DISCONTINUITY_CHECK)||((data[i]->hei-data[i-1]->hei)<-HEADING_DISCONTINUITY_CHECK))
         {
            //A wrap around has occured
            if(data[i]->hei > data[i-1]->hei)
            {
               //counterclockwise e.g. 1 -> 359
               headingwrap=1;
            }
            else
            {
               //clockwise e.g. 359 -> 1
               headingwrap=2;            
            }
         }
      }
   }

   //Do the smoothing
   for(int i=0;i<kernelsize;i++)
   {
      DEBUGPRINT(i<<" "<<store->time<<" "<<store->lat<<" "<<kernel[i])
      store->time=data[halflen]->time; //dont smooth time
      store->lat=store->lat + data[i]->lat*kernel[i];
      store->lon=store->lon + data[i]->lon*kernel[i];
      store->hei=store->hei + data[i]->hei*kernel[i];
      store->roll=store->roll + data[i]->roll*kernel[i];
      store->pitch=store->pitch + data[i]->pitch*kernel[i];
      //Heading needs to be handled differently due to discontinuity at 0/360
      if(headingwrap == 0)
         store->heading=store->heading + data[i]->heading*kernel[i];
      else if(headingwrap == 1)
      {
         //counter clockwise - so subtract 360 from the large headings
         if(data[i]->heading > HEADING_DISCONTINUITY_CHECK)
            store->heading=store->heading + (data[i]->heading-360)*kernel[i];

      }
      else if(headingwrap == 2)
      {
         //clockwise -
         if(data[i]->heading > HEADING_DISCONTINUITY_CHECK)
            store->heading=store->heading + (data[i]->heading-360)*kernel[i];
      }
   }   

   //Check that the heading is not -ve
   if(store->heading < 0)
      store->heading += 360;

   delete[] data;
   delete[] kernel;
}

//-------------------------------------------------------------------------
// Generates the kernel used to smooth the raw navigation in Triangle function
//-------------------------------------------------------------------------
void TriangleKernel(double* const kernel,const int length)
{
   int halflen=(length-1)/2.0;
   //set the points with weightings (end points are 0)
   for(int i=0;i<halflen;i++)
   {
      kernel[halflen+i]=kernel[halflen-i]=((halflen-i)/(double)(halflen*halflen));
   }
   //Set end points
   kernel[0]=0.0;
   kernel[length-1]=0.0;
}
