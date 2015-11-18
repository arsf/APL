//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include <string>
#include <typeinfo>
#include "TreeGrid.h"
#include "dataaccessor.h"

//-------------------------------------------------------------------------
//Keywords to describe Interpolator method
//-------------------------------------------------------------------------
namespace Interpolators
{
   enum InterpolatorType{NEARESTNEIGHBOUR,IDW,BILINEAR,BILINEARLEVEL3,CUBIC};
}

//-------------------------------------------------------------------------
// Base interpolator class
// Is now a template class so that it can support different level1 types
//-------------------------------------------------------------------------
template<class T>
class Interpolator
{
public:
   const double* const Data()const{return data;}

   Interpolator(int numbands)
   {
      length=numbands;
      data=new double[length];
      SetMaxInterpDistance(0.0); // default value to 0m squared
      this->tg=NULL;
      this->IGNOREVALUE=0;
      this->NODATAVALUE=0;
      this->l3pos=NULL;
      this->ignoredata=true;//default to interpolating over "ignore values"
   }

   virtual ~Interpolator()
   {
      if(data!=NULL)
         delete[] data;
   }

   virtual void Interpolate(std::vector<Item>* dp,const unsigned int* const bands,DataAccessor<T>* lev1data)=0;
   void SetL3Pos(IGMPoint* lp){this->l3pos=lp;}

   //set the maximum distance (squared) to interpolate over (units depends on input IGM data projection)
   void SetMaxInterpDistance(double mid) {sqmaxinterpdistance=mid*mid;}
   const double GetMaxInterpDistanceSq()const{return sqmaxinterpdistance;}

   Interpolators::InterpolatorType interpolator_type;

   virtual std::vector<Item>* UpdatePoints(unsigned int band,DataAccessor<T>* lev1data)=0;
   TreeGrid* tg;
   virtual void SetNumPoints(unsigned int np){numpoints=np;}
   virtual void SetSearchRadius(double sr){searchradius=sr;}
   virtual void SetIgnoreValue(double ig)
   {
      IGNOREVALUE=static_cast<T>(ig);
      Logger::Log("Interpolating level1 data ignoring values of: "+ToString(IGNOREVALUE)+" This may be different to entered value as it has been cast to level-1 data type.");
   }
   virtual void SetNoDataValue(double ndv)
   {
      NODATAVALUE=ndv;
   }
   virtual void SetIgnoreFlag(bool f){this->ignoredata=f;}

protected:

   T IGNOREVALUE;
   bool ignoredata; 
   double NODATAVALUE;  

   double* data;
   unsigned int length;  
   double sqmaxinterpdistance;
   IGMPoint* l3pos;  
   double searchradius;
   unsigned int numpoints;
};

//-------------------------------------------------------------------------
// Nearest Neighbour interpolator
//-------------------------------------------------------------------------
template<class T>
class NearestNeighbour : public Interpolator<T>
{
public:
   NearestNeighbour(int numbands) : Interpolator<T>(numbands) {this->interpolator_type=Interpolators::NEARESTNEIGHBOUR;}
   virtual void Interpolate(std::vector<Item>* dp,const unsigned int* const bands,DataAccessor<T>* lev1data=NULL);
protected:
   virtual std::vector<Item>* UpdatePoints(unsigned int band,DataAccessor<T>* lev1data);
};

//-------------------------------------------------------------------------
// Inverse Distance Weighted interpolator
//-------------------------------------------------------------------------
template<class T>
class IDW : public Interpolator<T>
{
public:
   IDW(int numbands) : Interpolator<T>(numbands) {this->interpolator_type=Interpolators::IDW;}
   virtual void Interpolate(std::vector<Item>* dp,const unsigned int* const bands,DataAccessor<T>* lev1data=NULL);

   //Function to calculate the sum of weights AND trim out unused points from the vector dp
   double CalculateSumWeights(std::vector<Item>* dp)
   {
      double sumweight=0;
      for(std::vector<Item>::iterator it=(*dp).begin();it<(*dp).end();it++)
      {
         if((*it).distance > this->sqmaxinterpdistance)
         {
            (*dp).erase(it,(*dp).end());
            Logger::Debug("Less than max number of points used closer than interpolation distance for this cell.");
            break;
         }
         sumweight += 1/((*it).distance);
      }
      return sumweight;
   }
protected:
   virtual std::vector<Item>* UpdatePoints(unsigned int band,DataAccessor<T>* lev1data);
};

//-------------------------------------------------------------------------
// Bilinear interpolator using 4 neighbouring level 1 points
//-------------------------------------------------------------------------
template<class T>
class Bilinear : public Interpolator<T>
{
public:
   Bilinear(int numbands) : Interpolator<T>(numbands) {this->interpolator_type=Interpolators::BILINEAR;}
   virtual void Interpolate(std::vector<Item>* dp,const unsigned int* const bands,DataAccessor<T>* lev1data);

protected:
   void FallBack(const unsigned int line,const unsigned int column,const unsigned int* const bands,DataAccessor<T>* lev1data=NULL);
   bool GetNeighbouringPoints(const IGMPoint* l3p,const Item* l1p, long* const igmrows, long* const igmcols);
   double GetInterpolatedValue(const double* const l3pos, const double* const datavalues);
   void GetUV(const double* const xpos,const double* const ypos, double* const interppos);
   virtual std::vector<Item>* UpdatePoints(unsigned int band,DataAccessor<T>* lev1data){return NULL;}
};

//-------------------------------------------------------------------------
// Bilinear interpolator using 4 nearest level 3 points
//-------------------------------------------------------------------------
template<class T>
class BilinearLevel3 : public Bilinear<T>
{
public:
   BilinearLevel3(int numbands) : Bilinear<T>(numbands) {this->interpolator_type=Interpolators::BILINEARLEVEL3;}
   virtual void Interpolate(std::vector<Item>* dp,const unsigned int* const bands,DataAccessor<T>* lev1data);

protected:
   virtual std::vector<Item>* UpdatePoints(unsigned int band,DataAccessor<T>* lev1data);
   virtual bool PreparePoints(std::vector<Item>* dp,double* interppos);
};

//-------------------------------------------------------------------------
// Cubic interpolator using 4*4 nearest level 3 points
//-------------------------------------------------------------------------
template<class T>
class Cubic : public Interpolator<T>
{
public:
   Cubic(int numbands) : Interpolator<T>(numbands) {this->interpolator_type=Interpolators::CUBIC;}
   virtual void Interpolate(std::vector<Item>* dp,const unsigned int* const bands,DataAccessor<T>* lev1data=NULL);

protected:
   virtual std::vector<Item>* UpdatePoints(unsigned int band,DataAccessor<T>* lev1data);
   void GetYInterpAtX(double* ypts,double x,double* xpts, double* fpts);
   double GetInterpValue(double x,double* xpts,double* fpts);
   double GetCardinal(int id,double* xpts, double* fpts);
   double GetFiniteDiff(int id,double* xpts, double* fpts);
   unsigned int GetH(double* h,double x,double* xpts);
   void OrderPoints(std::vector<Item>*dp,std::vector<Item> &ordered);
};



//----------------------------------------------------------------------
// Function to apply nearest neighbour interpolation
//----------------------------------------------------------------------
template <class T>
void NearestNeighbour<T>::Interpolate(std::vector<Item>* dp, const unsigned int* const bands,DataAccessor<T>* lev1data)
{
   //Get the first point in the vector - this assumes vector has been sorted by distance
   unsigned int line=(*dp)[0].igmrow;
   unsigned int column=(*dp)[0].igmcol;
   unsigned int new_line=0,new_column=0;

   //Vector to hold results of updated search (if required)
   std::vector<Item>* new_points=NULL;

   //Hold the value from the level1 data
   double value=0;

   for(unsigned int b=0;b<this->length;b++)
   {
      value=static_cast<double>(lev1data->GetData(bands[b],line,column));
      //Test if the value needs interpolating over - i.e. will not be used in the output map
      if((this->ignoredata==true)&&(value==this->IGNOREVALUE))
      {
         //Update the points
         new_points=UpdatePoints(bands[b],lev1data);
         //If no points were found set the value to NODATAVALUE
         if(new_points==NULL)
            value=this->NODATAVALUE;
         else
         {
            //update the line/column for this new point
            new_line=(*new_points)[0].igmrow;
            new_column=(*new_points)[0].igmcol;  
            //Get the value for this new point
            value=static_cast<double>(lev1data->GetData(bands[b],new_line,new_column));
         }
      }
      //Assign the value to the data array
      this->data[b]=value;
   }
}

//----------------------------------------------------------------------
// Function to update the points for the nearest neighbour interpolation
//----------------------------------------------------------------------
template <class T>
std::vector<Item>* NearestNeighbour<T>::UpdatePoints(unsigned int band,DataAccessor<T>* lev1data)
{
   //Find the nearest non-zero point 
   return this->tg->GetNearestXItems(this->numpoints,this->l3pos,this->searchradius,lev1data,band,this->IGNOREVALUE);
}

//----------------------------------------------------------------------
// Function to apply inverse distance weighted interpolation
//----------------------------------------------------------------------
template <class T>
void IDW<T>::Interpolate(std::vector<Item>* dp,const unsigned int* const bands,DataAccessor<T>* lev1data)
{
   //Clean data array
   for(unsigned int b=0;b<this->length;b++)
      this->data[b]=0;   

   //Remember that 'distance' member of Item is actually the distance squared
   //Calculate the sum of the inverse weights ( e.g.sum of  1 / d*d)
   double sumweight=CalculateSumWeights(dp);

   //Store this value incase we need to update points
   double orig_sumweight=sumweight;
   //Copy the original points in case we update (the vector will be overwritten if we do)
   std::vector<Item> orig_points(*dp);
   std::vector<Item>* points=&orig_points;
   //Vector to hold results of updated search (if required)
   std::vector<Item>* new_points=NULL;
   //To keep track on whether we break out of a loop or not
   bool loopbreaker=false;

   //To hold the level 1 data value in
   double value=0;

   for(unsigned int b=0;b<this->length;b++)
   {
      //If no points were found set the value to NODATAVALUE
      if(points==NULL)
      {
         this->data[b]=this->NODATAVALUE;
      }
      else
      {
         for(std::vector<Item>::iterator it=(*points).begin();it<(*points).end();it++)
         {
               //Get the level 1 value
               value=static_cast<double>(lev1data->GetData(bands[b],(*it).igmrow,(*it).igmcol));

               //Check if value = IGNOREVALUE - if so should re-search for new non-zero points
               if((this->ignoredata==true)&&(value==this->IGNOREVALUE))
               {
                  //Need to remove this point - re-search for nearest non-zero points
                  new_points=this->UpdatePoints(bands[b],lev1data);
                  loopbreaker=true;
                  break;
               }
               //Add onto the sum of weighted data values
               this->data[b] += ((1/(*it).distance) / sumweight)*(value);
         }
      }

      //Check if we broke out of the loop - if so then update to use new values for the band 2nd time around
      if(loopbreaker==true)
      {
         //skip if no new points were found
         if(new_points==NULL)
            continue;

         //recalculate sum of weights
         sumweight=CalculateSumWeights(new_points);
         //point to the new updated vector
         points=new_points;
         this->data[b]=0; //reset the interpolated sum value to 0
         b=b-1; //subtract one off b to "trick" loop into doing the last band again
         loopbreaker=false;
      }
      else
      {
         //use the original sum of weights
         sumweight=orig_sumweight;
         //point to the original set of points
         points=&orig_points;
      }
   }
}

//----------------------------------------------------------------------
// Function to update the points for the IDW interpolation
//----------------------------------------------------------------------
template<class T>
std::vector<Item>* IDW<T>::UpdatePoints(unsigned int band,DataAccessor<T>* lev1data)
{
   //Find the nearest non-zero point 
   return this->tg->GetNearestXItems(this->numpoints,this->l3pos,this->searchradius,lev1data,band,this->IGNOREVALUE);
}

//----------------------------------------------------------------------
// Function to apply Bilinear interpolation
//----------------------------------------------------------------------
template <class T>
void Bilinear<T>::Interpolate(std::vector<Item>* dp,const unsigned int* const bands,DataAccessor<T>* lev1data)
{
   //Clean data array
   for(unsigned int b=0;b<this->length;b++)
      this->data[b]=0; 

   //First find the 4 level 1 pixels that we need to interpolate between
   //These are found by taking the nearest point and then the 3 level 1 file neighbours

   //Store the 4 points row and cols in these arrays   
   long igmrows[2]={0};
   long igmcols[2]={0};   
   //Store the 4 points x,y positions in these arrays
   double xpos[4]={0};
   double ypos[4]={0};
   //store the interpolated position in this array - l3pos converted to a row/col offset
   double interppos[2]={0};
   //Store the 4 data values of the band being interpolated into this array
   double datavalues[4]={0};
   //Get the itemdata to read X and Y values from the row/col pairs
   const ItemData* const idata=(*dp)[0].GetData();

   Logger::Debug("Nearest point is: "+ToString((*dp)[0].igmrow)+" "+ToString((*dp)[0].igmcol));
   Logger::Debug("Level3pos: "+ToString(this->l3pos->X)+" "+ToString(this->l3pos->Y));

   //Get the 4 IGM points that we need to interpolate between and put 
   //into igmrows and igmcols variables, testing return value for edge points
   if((!GetNeighbouringPoints(this->l3pos,&((*dp)[0]),igmrows,igmcols)))
   {
      //This is an edge cell so no neighbours to use for interpolation
      //TODO: Update to allow interpolation for points in edge cells that lie in an internal quadrant and can be interpolated
      //Do nearest neighbour - actually just enter 0's for now
      //Logger::Debug("Using nearest neighbour for this point - suspected edge point.");
      for(unsigned int b=0;b<this->length;b++)
      {
         //FallBack((*dp)[0].igmrow,(*dp)[0].igmcol,bands,lev1data);
         this->data[0]=this->NODATAVALUE;
      }
      return;
   }
   else
   {
//      //Do bilinear interpolation
      int counter=0;
//      //Get the X,Y positions of the data values and store in xpos, ypos
//      for(unsigned int rr=igmrows[0];rr<=igmrows[1];rr++)
//      {
//         for(unsigned int cc=igmcols[0];cc<=igmcols[1];cc++)
//         {
//            xpos[counter]=idata->GetX(rr,cc);
//            ypos[counter]=idata->GetY(rr,cc);
//            Logger::Debug("X,Y of point "+ToString(counter)+": "+ToString(xpos[counter])+" "+ToString(ypos[counter]));
//            counter++; 
//         }
//      }

//      //Calculate U and V such that Xx=Ax+V*(Bx-Ax),  Xy=Ay+V*(By-Ay)
//      //and Ax=Px+U(Rx-Px),  Ay=Py+U(Ry-Py), Bx=Qx+U(Sx-Qx), By=Qy+U(Sy-Qy)
//      //U and V then describe the position of X in terms of fraction from P [U for PQ and V for AB]
//      GetUV(xpos,ypos,interppos);

//      //Check if the point returned was bad
////      if((interppos[0]==-9999)&&(interppos[1]==-9999))
////      {
////         std::cout<<"Bad point detected - Need to do something about it!"<<std::endl;
////         for(unsigned int b=0;b<length;b++)
////         {
////            data[b]=100000;
////         }
////         return;
////      }

//      Logger::Debug("First identified U,V: "+ToString(interppos[0])+" "+ToString(interppos[1]));
//      //Check if the values of U and V are between 0 and 1. If not, enter loop and try and correct by 
//      //updating the 4 neighbours.
      
//      if((interppos[0] < 0) || (interppos[1] < 0) || (interppos[0] >1) || (interppos[1] >1))
      if(true)
      {
         //int itercounter=0;

         for(unsigned int ptcount=0;ptcount<(*dp).size();ptcount++)
         {

            if(!GetNeighbouringPoints(this->l3pos,&((*dp)[ptcount]),igmrows,igmcols))
            {
               continue;
            }

            //Do bilinear interpolation
            int counter=0;
            //Get the X,Y positions of the data values and store in xpos, ypos
            for(unsigned int rr=igmrows[0];rr<=igmrows[1];rr++)
            {
               for(unsigned int cc=igmcols[0];cc<=igmcols[1];cc++)
               {      
                  //Skip if any item is an edge point - this is to give us "tidy edges"
                  //but will mean that data at the very edge pixel will not be mapped
                  if((cc==0)||(cc==idata->igm->Samples()-1))
                     return;

                  xpos[counter]=idata->GetX(rr,cc);
                  ypos[counter]=idata->GetY(rr,cc);
                  Logger::Debug("X,Y of point "+ToString(counter)+": "+ToString(xpos[counter])+" "+ToString(ypos[counter]));
                  counter++; 
               }
            }

            //Calculate U and V such that Xx=Ax+V*(Bx-Ax),  Xy=Ay+V*(By-Ay)
            //and Ax=Px+U(Rx-Px),  Ay=Py+U(Ry-Py), Bx=Qx+U(Sx-Qx), By=Qy+U(Sy-Qy)
            //U and V then describe the position of X in terms of fraction from P [U for PQ and V for AB]
            GetUV(xpos,ypos,interppos);

            //go onto next iteration
            int itercounter=0;
            while(((interppos[0] < 0) || (interppos[1] < 0) || (interppos[0] >1) || (interppos[1] >1))&&(itercounter<10))
            {
               //std::cout<<"Iteration number: "<<itercounter<<std::endl;
               //if((interppos[0]>50)||(interppos[1]>50))
               //   break; //these differences are too large to make sense really

               //Update the pixels that we want to interpolate
               igmcols[0]+=floor(interppos[0]);
               igmcols[1]+=floor(interppos[0]);
               igmrows[0]+=floor(interppos[1]);
               igmrows[1]+=floor(interppos[1]);

               //Check if out of bounds of IGM file
               if((igmrows[0] < 0)||(igmrows[1] < 0)||(igmcols[0] < 0)||(igmcols[1] < 0)||
                  (igmrows[0] >= idata->igm->Lines())||(igmrows[1] >= idata->igm->Lines())||(igmcols[0] >= idata->igm->Samples())||(igmcols[1] >= idata->igm->Samples()))
                  {break;}

               //Need to get new x,y positions and recalculate V,U if needs to be precise.
               counter=0;

               for(unsigned int rr=igmrows[0];rr<=igmrows[1];rr++)
               {
                  for(unsigned int cc=igmcols[0];cc<=igmcols[1];cc++)
                  {                  
                     xpos[counter]=idata->GetX(rr,cc);
                     ypos[counter]=idata->GetY(rr,cc);
                     Logger::Debug("New X,Y of point "+ToString(counter)+"("+ToString(rr)+", "+ToString(cc)+"): "+ToString(xpos[counter])+" "+ToString(ypos[counter]));
                     counter++;
                  }
               }
               //Get updated UV positions
               GetUV(xpos,ypos,interppos); 

               itercounter++;
            }

            if((interppos[0] < 0) || (interppos[1] < 0) || (interppos[0] >1) || (interppos[1] >1))
            {
            }
            else
            {
               //these are OK
               break;
            }
         }

         //UV are still bad after looking at all 4 pts
         if((interppos[0] < 0) || (interppos[1] < 0) || (interppos[0] >1) || (interppos[1] >1))
         {
            std::cout<<"Bilinear failed - falling back to nearest neighbour for this pixel: "<<this->l3pos->X<<" "<<this->l3pos->Y<<std::endl;
            FallBack((*dp)[0].igmrow,(*dp)[0].igmcol,bands,lev1data);
            return;
         }



//         while((interppos[0] < 0) || (interppos[1] < 0) || (interppos[0] >1) || (interppos[1] >1))
//         {
//            //Update the pixels that we want to interpolate
//            igmcols[0]+=floor(interppos[0]);
//            igmcols[1]+=floor(interppos[0]);
//            igmrows[0]+=floor(interppos[1]);
//            igmrows[1]+=floor(interppos[1]);

//            if((igmrows[0] < 0)||(igmrows[1] < 0)||(igmcols[0] < 0)||(igmcols[1] < 0)||
//               (igmrows[0] >= idata->igm->Lines())||(igmrows[1] >= idata->igm->Lines())||(igmcols[0] >= idata->igm->Samples())||(igmcols[1] >= idata->igm->Samples()))
//            {
//               //New points are outside limits so do not do bilinear interpolation 
//               Logger::Debug("Using nearest neighbour for this point - suspected edge point.");
//               for(unsigned int b=0;b<length;b++)
//               {
//                  FallBack((*dp)[0].igmrow,(*dp)[0].igmcol,bands,lev1data);
//               }
//               return;
//            }

//            //Need to get new x,y positions and recalculate V,U if needs to be precise.
//            counter=0;
//            for(unsigned int rr=igmrows[0];rr<=igmrows[1];rr++)
//            {
//               for(unsigned int cc=igmcols[0];cc<=igmcols[1];cc++)
//               {
//                  xpos[counter]=idata->GetX(rr,cc);
//                  ypos[counter]=idata->GetY(rr,cc);
//                  Logger::Debug("New X,Y of point "+ToString(counter)+": "+ToString(xpos[counter])+" "+ToString(ypos[counter]));
//                  counter++;
//               }
//            }

//            //Get updated UV positions
//            GetUV(xpos,ypos,interppos);            

//            if((itercounter>10)||((interppos[0]==-9999)&&(interppos[1]==-9999)))
//            {
////               //Stuck in an infinite loop or another bad pixel returned
////               std::cout<<"Bad point detected - Need to do something about it!"<<std::endl;
////               for(unsigned int b=0;b<length;b++)
////               {
////                  data[b]=100000;
////               }
//               std::cout<<"Bilinear failed - falling back to nearest neighbour for this pixel: "<<l3pos->X<<" "<<l3pos->Y<<std::endl;
//               FallBack((*dp)[0].igmrow,(*dp)[0].igmcol,bands,lev1data);
//               return;
//            }
//            itercounter++;
//         }
         
      }

      Logger::Debug("Interpolate position to be used in col / row offset is: "+ToString(interppos[0])+" "+ToString(interppos[1]));

      for(unsigned int b=0;b<this->length;b++)
      {
         //Get the 4 data values from RAM or disk and store in datavalues [ordered row0,col0 row0,col1 row1,col0 row1,col1]
         counter=0;
         //Iterate through each igmrow/col pair of the 4 points
         for(unsigned int rr=igmrows[0];rr<=igmrows[1];rr++)
         {
            for(unsigned int cc=igmcols[0];cc<=igmcols[1];cc++)
            {
               datavalues[counter]=static_cast<double>(lev1data->GetData(bands[b],rr,cc));
               //Increase the counter
               counter++;
            }
         }

         //Interpolate to get the data value we want
         this->data[b]=GetInterpolatedValue(interppos,datavalues);
      }
   }
}

//----------------------------------------------------------------------
// Function to perform the bilinear interpolation using given values
//----------------------------------------------------------------------
template <class T>
inline double Bilinear<T>::GetInterpolatedValue(const double* const pos, const double* const datavalues)
{
   //Need to interpolate linearly horizontally twice and then once vertically
   // a  m'  b     i.e between a and b, and c and d to get m' and m''
   //    p             and then between m' and m'' to get p
   // c  m'' d
   //This simplifies to the below assuming pixels are on a regular spaced grid
   //of equal distances separating the values with x,y positions of 0,0 0,1 1,0 1,1

   return (datavalues[0]*(1-pos[1])*(1-pos[0]) + datavalues[1]*(1-pos[1])*pos[0] 
                  + datavalues[2]*pos[1]*(1-pos[0]) + datavalues[3]*pos[1]*pos[0]);
}

//----------------------------------------------------------------------
// Function to get row/col indices of the points required for bilinear 
//----------------------------------------------------------------------
template <class T>
bool Bilinear<T>::GetNeighbouringPoints(const IGMPoint* l3p,const Item* l1p,long* const igmrows,long* const igmcols)
{
   //get the itemdata
   const ItemData* const idata=l1p->GetData();
   long int i=0,j=0,k=0;
   double xt=0,yt=0;
   double distancesq[2]={0};

   //If level1 pixel is not on an edge of the level1 file
   if((l1p->igmrow > 0) && (l1p->igmrow < idata->igm->Lines()-1) && (l1p->igmcol > 0) && (l1p->igmcol < idata->igm->Samples()-1))
   {
      //Get row and col of point above and below l1p
      j=l1p->igmcol;
      for(i=l1p->igmrow-1;i<=l1p->igmrow+1;i+=2)
      {
         //Get the x,y for this point
         xt=idata->GetX(i,j);
         yt=idata->GetY(i,j);
         //Calculate the squared distance to l3p
         distancesq[k]=(xt-l3p->X)*(xt-l3p->X) + (yt-l3p->Y)*(yt-l3p->Y);
         k++;
      }

      if(distancesq[0] < distancesq[1])
      {
         //Minimum is distancesq[0] so we want (row-1)  
         igmrows[1]=l1p->igmrow;
         igmrows[0]=igmrows[1] -1;
      }
      else
      {
         //Minimum is distancesq[1] so we want (row+1)  
         igmrows[0]=l1p->igmrow;
         igmrows[1]=igmrows[0] +1;         
      }

      //Get row and col of point left and right of l1p
      k=0;
      i=l1p->igmrow;
      for(j=l1p->igmcol-1;j<=l1p->igmcol+1;j+=2)
      {
         //Get the x,y for this point
         xt=idata->GetX(i,j);
         yt=idata->GetY(i,j);
         //Calculate the squared distance to l3p
         distancesq[k]=(xt-l3p->X)*(xt-l3p->X) + (yt-l3p->Y)*(yt-l3p->Y);
         k++;
      }
      if(distancesq[0] < distancesq[1])
      {
         //Minimum is distancesq[0] so we want (col-1)  
         igmcols[1]=l1p->igmcol;
         igmcols[0]=igmcols[1] -1;
      }
      else
      {
         //Minimum is distancesq[1] so we want (col+1)  
         igmcols[0]=l1p->igmcol;
         igmcols[1]=igmcols[0] +1;         
      }
   }
   else
   {
      //An edge pixel - lets just do nearest neighbour for this.
      return false;
   }
   Logger::Debug("Using the pixels at these rows and cols to interpolate: "+ToString(igmrows[0])+" "+ToString(igmrows[1])+" "+ToString(igmcols[0])+" "+ToString(igmcols[1]));
   return true;
}


//bool Bilinear::GetNeighbouringPoints(const IGMPoint* l3p,const Item* l1p,long* const igmrows,long* const igmcols)
//{
//   //get the itemdata
//   const ItemData* const idata=l1p->GetData();
//   double xta=0,yta=0;
//   double xtb=0,ytb=0;
//   double xtr=0,ytr=0;
//   double lev3X=0,lev3Y=0;
//   double cp_ref=0,cp_test=0;

//   //If level1 pixel is not on an edge of the level1 file
//   if((l1p->igmrow > 0) && (l1p->igmrow < idata->igm->Lines()-1) && (l1p->igmcol > 0) && (l1p->igmcol < idata->igm->Samples()-1))
//   {
//      //Get XY of ref point (level 1 point)
//      double l1pX=idata->GetX(l1p->igmrow,l1p->igmcol);
//      double l1pY=idata->GetY(l1p->igmrow,l1p->igmcol);
//      //Get distance to ref point from l3point
//      lev3X=l3p->X-l1pX;
//      lev3Y=l3p->Y-l1pY;

//      //Get distance to ref point from XY of point above ref point (i.e. point previous in IGM file)
//      xta=idata->GetX(l1p->igmrow-1,l1p->igmcol)-l1pX;
//      yta=idata->GetY(l1p->igmrow-1,l1p->igmcol)-l1pY;
//      //Test cross product and record "sign"
//      cp_ref=xta*lev3Y - yta*lev3X;

//      //Get XY of point 'to the right' (i.e. to the column + 1 in IGM)
//      xtr=idata->GetX(l1p->igmrow,l1p->igmcol+1)-l1pX;
//      ytr=idata->GetY(l1p->igmrow,l1p->igmcol+1)-l1pY;
//      //Test cross product - if sign change then we want points above, right and above-right
//      cp_test=xta*ytr - yta*xtr;

//      if(cp_test*cp_ref >= 0)
//      {
//         //Level 3 point is on the same side as 'point to the right'
//         igmcols[0]=l1p->igmcol;
//         igmcols[1]=l1p->igmcol+1;
//      }
//      else
//      {
//         //Level 3 point is on the opposite side as 'point to the right'
//         igmcols[0]=l1p->igmcol-1;
//         igmcols[1]=l1p->igmcol;         
//      }

//      //Get distance to ref point from XY of point below
//      xtb=idata->GetX(l1p->igmrow+1,l1p->igmcol)-l1pX;
//      ytb=idata->GetY(l1p->igmrow+1,l1p->igmcol)-l1pY;

//      //Test cross product and record "sign"
//      cp_ref=xtr*lev3Y - ytr*lev3X;
//      cp_test=xtr*ytb - ytr*xtb;

//      if(cp_test*cp_ref >= 0)
//      {
//         //Level 3 point is on the same side as 'point below'
//         igmrows[0]=l1p->igmrow;
//         igmrows[1]=l1p->igmrow+1;
//      }
//      else
//      {
//         //Level 3 point is on the opposite side as 'point to the right'
//         igmrows[0]=l1p->igmrow-1;
//         igmrows[1]=l1p->igmrow;         
//      }

//      //Now we should check using the opposite corner point to see if above remains true
//      //or if the point is not actually within this quad at all. TODO


//   }
//   else
//   {
//      //An edge pixel - lets just do nearest neighbour for this.
//      return false;
//   }
//   Logger::Debug("Using the pixels at these rows and cols to interpolate: "+ToString(igmrows[0])+" "+ToString(igmrows[1])+" "+ToString(igmcols[0])+" "+ToString(igmcols[1]));
//   return true;
//}


//----------------------------------------------------------------------
// Function to get row/col position of level3 pixel in level1 coordinates
//----------------------------------------------------------------------
template <class T>
inline void Bilinear<T>::GetUV(const double* const xpos,const double* const ypos,double* const interppos)
{
   //Convert the l3pos (X in diagram below) to a l1(c,r)
   //P   Q
   //
   //
   // A X B
   //
   //  R   S

   //Find the distances required for the later calculations
   double PQx=xpos[1] - xpos[0];
   double PQy=ypos[1] - ypos[0];
   double PRx=xpos[2] - xpos[0];
   double PRy=ypos[2] - ypos[0];
   double RSx=xpos[3] - xpos[2];
   double RSy=ypos[3] - ypos[2];      
   double PXx=this->l3pos->X - xpos[0];
   double PXy=this->l3pos->Y - ypos[0];  

   //We need to solve a quadratic eqn - split into a, b, c and use formula
   double qa=PQy*(RSx-PQx) - PQx*(RSy-PQy);
   double qb=PRx*PQy - PRy*PQx + PXx*(RSy-PQy) - PXy*(RSx-PQx);
   double qc=PRy*PXx - PRx*PXy;

   //Work through the possible solutions for U/V
   if(qa==0)
   {
      //Linear equation - u^2 part is 0
      if(qb!=0)
         interppos[0]=-qc/qb;
      else
         throw "Error calculating U,V position of level3 pixel in terms of level1 row,col. There appears to be no U term in equation.";
   }
   else
   {
      //Quadratic equation - check square root will not be -ve
      double sqrtpart=qb*qb - 4*qa*qc;
      if(sqrtpart>=0)
      {
         double U1=(-qb + sqrt(sqrtpart))/(2*qa); //+ve root
         double U2=(-qb - sqrt(sqrtpart))/(2*qa); //-ve root
         //Select the first root we find that is 0 <= root <= 1
         if((U1>=0)&&(U1<=1))
         {
            //Select this root
            interppos[0]=U1;
         }
         else if((U2>=0)&&(U2<=1))
         {
            //Select this root
            interppos[0]=U2;
         }
         else
         {
            //Neither root of U is within the bounds of 0<=U<=1
            //This probably means we have the incorrect 4 nearest neighbours - we should try and correct this
            //Select the root with smallest absolute value. This will be used as a start point to find new neighbours.
            Logger::Debug("Neither U root is within the 0,1 bounds - try and find new 4 neighbours.");
            if(fabs(U1)<fabs(U2))
            {
               interppos[0]=U1;  
            }
            else 
            {
               interppos[0]=U2;
            }
            Logger::Debug("Selecting U root: "+ToString(interppos[0]));
            //Now use it to find the 4 points that surround the level3 position according to the U,V values
            //This is in the main bilinear function
         }
      }
      else
      {
         //-ve square root
         Logger::Debug("Negative square root in bilinear UV calculation. Quadratic parameters a,b and c: "+ToString(qa)+" "+ToString(qb)+" "+ToString(qc));

//         if(qa < 0.05)
//         {
            //Make a small value assumption - u^2 is close to 0 so get an approx UV from using a linear eqn
            interppos[0]=(-qc/qb);
//         }
//         else
//         {
//            //qa is not that small - set as a bad pixel and deal with it back in the main bilinear function
//            interppos[0]=-9999;
//            interppos[1]=-9999;
//            return;
//         }
      }
   }

   //Now find V - check X denominator will not be 0, if it is then use the Y eqn instead.
   if((PRx + interppos[0]*(RSx-PQx))!=0)
   {
      interppos[1]=(PXx-interppos[0]*PQx)/(PRx + interppos[0]*(RSx-PQx));
   }
   else if((PRy + interppos[0]*(RSy-PQy))!=0)
   {
      interppos[1]=(PXy-interppos[0]*PQy)/(PRy + interppos[0]*(RSy-PQy));
   }
   else
   {
      //Both will give divide by 0 errors - v undefined
      throw "Error in Bilinear UV calculation - there is no solution for V - it is undefined.";
   }
}

//----------------------------------------------------------------------
// Function to fall back to if bilinear fails - nearest neighbour
// - Note that this is a hack fix and should figure out a nice way
//   to fix the bilinear algorithm
//----------------------------------------------------------------------
template <class T>
void Bilinear<T>::FallBack(const unsigned int line,const unsigned int column,const unsigned int* const bands,DataAccessor<T>* lev1data)
{
   //Now check if the data exists in RAM or if we need to read it from disk
   for(unsigned int b=0;b<this->length;b++)
   {
      this->data[b]=static_cast<double>(lev1data->GetData(bands[b],line,column));
   }
}



//----------------------------------------------------------------------
// Function to apply Bilinear interpolation using a quadrilateral 
// formed by the 4 nearest level 3 points that create a quad containing
// the point to be interpolated
//----------------------------------------------------------------------
template <class T>
void BilinearLevel3<T>::Interpolate(std::vector<Item>* dp,const unsigned int* const bands,DataAccessor<T>* lev1data)
{
   //We know the 4 points that we want to use
   //Get UV
   double orig_interppos[2]={0};
   double new_interppos[2]={0};
   double* interppos=NULL;

   //Clean data array
   for(unsigned int b=0;b<this->length;b++)
      this->data[b]=0; 

   if((*dp).size()!=4)
      throw "Bilinear expects only 4 points to be used. Got: "+ToString((*dp).size());

   //Copy the original points in case we update (the vector will be overwritten if we do)
   std::vector<Item> orig_points((*dp));

   std::vector<Item>* points=&orig_points;
   //Vector to hold results of updated search (if required)
   std::vector<Item>* new_points=NULL;
   bool loopbreaker=false;
   bool nonewpoints=false;

   //Store the 4 data values of the band being interpolated into this array
   double datavalues[4]={0};

   //Get interppos
   if(!this->PreparePoints(&orig_points,orig_interppos))
      return;

   //point the interppos pointer to the data for the loop below
   interppos=orig_interppos;

   for(unsigned int b=0;b<this->length;b++)
   {
      //Get the 4 data values from RAM or disk and store in datavalues [ordered row0,col0 row0,col1 row1,col0 row1,col1]
      if(points==NULL)
      {
         for(unsigned int i=0;i<4;i++)
            datavalues[i]=this->NODATAVALUE;
      }
      else
      {
         //Iterate through each igmrow/col pair of the 4 points
         for(unsigned int i=0;i<4;i++)
         {
            datavalues[i]=static_cast<double>(lev1data->GetData(bands[b],(*points)[i].igmrow,(*points)[i].igmcol));

            //Test if the data value is OK to use
            if((this->ignoredata==true)&&(datavalues[i]==this->IGNOREVALUE))
            {
               //Update the points
               new_points=this->UpdatePoints(bands[b],lev1data);
               if(new_points==NULL)
                  nonewpoints=true;

               loopbreaker=true;
               break;
            }
         }
      }

      if(loopbreaker==true)
      {
         //reset the bool
         loopbreaker=false;

         if(nonewpoints==true)
         {
            nonewpoints=false;
            this->data[b]=this->NODATAVALUE;
         }
         else
         {
            //point to the new collection of points
            points=new_points;

            //Prepare the points using the new data and point to new array
            if(!this->PreparePoints(new_points,new_interppos))
            {
               points=&orig_points;
               this->data[b]=this->NODATAVALUE;
            }
            else
            {
               //subtrcat one off bands to trick loop into redoing this band
               b=b-1;
               interppos=new_interppos;
            }
         }
      }
      else
      {
         //Interpolate to get the data value we want
         this->data[b]=this->GetInterpolatedValue(interppos,datavalues);
         //reset to original data for next band
         points=&orig_points;
         interppos=orig_interppos;
      }
   }
}

//-------------------------------------------------------------------------
//Function to get the interpolate position
//-------------------------------------------------------------------------
template<class T>
bool BilinearLevel3<T>::PreparePoints(std::vector<Item>* dp,double* interppos)
{
   //Get the itemdata to read X and Y values from the row/col pairs
   const ItemData* const idata=(*dp)[0].GetData();

   double xpos[4]={0};
   double ypos[4]={0};

   for(int i=0;i<4;i++)
   {
      //Skip if any item is an edge point - this is to give us "tidy edges"
      //but will mean that data at the very edge pixel will not be mapped
      if(((*dp)[i].igmcol==0)||((*dp)[i].igmcol==idata->igm->Samples()-1))
         return false;

      //Get the X and Y positions from the IGM data for this row/col
      xpos[i]=idata->GetX((*dp)[i].igmrow,(*dp)[i].igmcol);
      ypos[i]=idata->GetY((*dp)[i].igmrow,(*dp)[i].igmcol);
   }

   //Calculate the U,V position that the interpolate point has within the 
   //quadrilateral formed by the 4 xpos/ypos pairs
   this->GetUV(xpos,ypos,interppos);

   return true;
}


//-------------------------------------------------------------------------
// Function to update the points used
//-------------------------------------------------------------------------
template<class T>
std::vector<Item>* BilinearLevel3<T>::UpdatePoints(unsigned int band,DataAccessor<T>* lev1data)
{
   //Find the nearest non-zero points
   return this->tg->GetQuadItems(1,this->l3pos,this->searchradius,lev1data,band,this->IGNOREVALUE);
}


//-------------------------------------------------------------------------
// CUBIC STUFF
// Method assumes array data is monotonic increasing
//
// Not the most efficient algorithm - involves lots of sorting / reordering
// points and transfering from vectors / arrays / classes
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
//Hermite multipliers
//-------------------------------------------------------------------------
template<class T>
unsigned int Cubic<T>::GetH(double* h,double x,double* xpts)
{
   unsigned int i=0;
   double t=0;

   //Get the element of the xpts array that the interpolated
   //point lies between (i and i-1)
   while((x > xpts[i])&&(i < 4))
      i++;

   if(i >= 4)
   {
      //Use WarnOnce so as not to output for each band being mapped - will still warn for other points as includes the interpolate point value
      Logger::WarnOnce("Interpolated point is outside the maximum bound in cubic interpolation - this should not happen.\n"
                       "Will use maximum bound to calculate Hermite multiplier for this point (with x or y value): "+ToString(x));
      //Set the variables to be relative to the maximum bound to prevent crashing
      i=3;
      t=1;
   }
   else if(i==0)
   {
      //Use WarnOnce so as not to output for each band being mapped - will still warn for other points as includes the interpolate point value
      Logger::WarnOnce("Interpolated point is outside the minimum bound in cubic interpolation - this should not happen.\n"
                       "Will use minimum bound to calculate Hermite multiplier for this point (with x or y value): "+ToString(x));
      //Set the variables to be relative to the minimum bound to prevent crashing
      i=1;
      t=0;
   }
   else
   {
      //t is the position between 0-1 where the interpolated point lies between i and i-1
      t=(x-xpts[i-1])/(xpts[i]-xpts[i-1]);
   }

   //Caclulate the Hermite multipliers
   h[0]=2*t*t*t - 3*t*t + 1;
   h[1]=t*t*t - 2*t*t + t;
   h[2]=-2*t*t*t + 3*t*t;
   h[3]=t*t*t - t*t;

   return i;
}

//-------------------------------------------------------------------------
//Get the tangents using finite difference method
//-------------------------------------------------------------------------
template<class T>
double Cubic<T>::GetFiniteDiff(int id,double* xpts, double* fpts)
{
   double val=0;
   if(id==0)
      val=(fpts[1]-fpts[0])/(2*(xpts[1]-xpts[0]));
   else if(id==1)
      val=((fpts[2]-fpts[1]) / (2*(xpts[2]-xpts[1])) + (fpts[1]-fpts[0])/(2*(xpts[1]-xpts[0])));
   else if(id==2)
      val=((fpts[3]-fpts[2]) / (2*(xpts[3]-xpts[2])) + (fpts[2]-fpts[1])/(2*(xpts[2]-xpts[1])));
   else if(id==3)
      val=(fpts[3]-fpts[2]) / (2*(xpts[3]-xpts[2]));
   return val;
}

//-------------------------------------------------------------------------
//Get the tangents for using cardinal splines (actually catmull-rom as tension=1 here)
//-------------------------------------------------------------------------
template<class T>
double Cubic<T>::GetCardinal(int id,double* xpts, double* fpts)
{
   double val=0;
   if(id==0)
      val=0;
   else if(id==1)
      val=((fpts[2]-fpts[0]) / (xpts[2]-xpts[0])) ;
   else if(id==2)
      val=((fpts[3]-fpts[1]) / (xpts[3]-xpts[1]));
   else if(id==3)
      val=0;
   return val;
}

//-------------------------------------------------------------------------
//Get the interpolated value at point x from the xpts and fpts
//-------------------------------------------------------------------------
template<class T>
double Cubic<T>::GetInterpValue(double x,double* xpts,double* fpts)
{
   double h[4]={0};
   unsigned int i=GetH(h,x,xpts);
   //double m1=GetFiniteDiff(1,xpts,fpts);
   //double m2=GetFiniteDiff(2,xpts,fpts);
   double m1=GetCardinal(i-1,xpts,fpts);
   double m2=GetCardinal(i,xpts,fpts);
   double f=h[0]*fpts[i-1] + h[1]*m1 + h[2]*fpts[i] + h[3]*m2;
   return f;
}

//-------------------------------------------------------------------------
//Get the 4 interpolated values (one from each set of splines)
//-------------------------------------------------------------------------
template<class T>
void Cubic<T>::GetYInterpAtX(double* ypts,double x,double* xpts, double* fpts)
{
   //xpts and fpts are 16 in size, ypts 4 in size
   for(int i=0;i<4;i++)
   {
      //Interpolate at this x to get the y values
      ypts[i]=GetInterpValue(x,&xpts[4*i],&fpts[4*i]);
   }
}

//-------------------------------------------------------------------------
//Class to allow sorting of 2 arrays based on the elements of one of them
//In this case it is used to sort the interpolated Y and F values based on
//the Y data, such that Y will be monotonic increasing.
//-------------------------------------------------------------------------
class Pair
{
public:
   Pair(){y=0;f=0;}
   Pair(double a,double b)
   {
      y=a;
      f=b;
   }
   double y,f;
   static bool compare(Pair p1,Pair p2)
   {
      return (p1.y < p2.y);
   }
};

//-------------------------------------------------------------------------
// Function to order the points returned from the TreeGrid search into
// an order suitable to enter the cubic spline functions
// Need to order points such they are
//  1 2  3 4
//  5 6  7 8
//      x
//  9 10 11 12
// 13 14 15 16
//
// currently they are
//  1-4 5-8
//     x
// 9-12 13-16
//
// Each set of 4 (quadrant) needs to be ordered by N-S and then by E-W to get:
// 1 2
// 3 4
//-------------------------------------------------------------------------
template<class T>
void Cubic<T>::OrderPoints(std::vector<Item>*dp,std::vector<Item> &ordered)
{
   //Create a vector to hold the ordered points of each quadrant and set size to 4
   std::vector<Item> semi;
   semi.reserve(4);

   //Do not clear ordered vector before using as it has been assigned a size
   //and as long as no inserts are used it should be ok (i.e. only use [] to assign / overwrite)

   std::vector<Item>::iterator it=(*dp).begin();

   //Sort the elements in the UL quad 
   //elements 0-3 of original array 
   //elements 0,1,4,5 in ordered array
   for(it=(*dp).begin();it<(*dp).begin()+4;it++)
   {
      //This loop orders by N-S
      std::vector<Item>::iterator sit=semi.begin();
      while((sit<semi.end())&&((*it).Y() < (*sit).Y()))
        sit++;
      semi.insert(sit,(*it));
   }

   //Now insert elements 0,1 based on E-W
   if(semi[0].X() < semi[1].X())
   {
      ordered[0]=semi[0];
      ordered[1]=semi[1];
   }
   else
   {
      ordered[0]=semi[1];
      ordered[1]=semi[0];
   }

   //Now insert elements 4,5 based on E-W
   if(semi[2].X() < semi[3].X())
   {
      ordered[4]=semi[2];
      ordered[5]=semi[3];
   }
   else
   {
      ordered[4]=semi[3];
      ordered[5]=semi[2];
   }

   //UR quad
   semi.clear();
   for(it=(*dp).begin()+4;it<(*dp).begin()+8;it++)
   {
      std::vector<Item>::iterator sit=semi.begin();
      while((sit<semi.end())&&((*it).Y() < (*sit).Y()))
         sit++;
      semi.insert(sit,(*it));
   }

   if(semi[0].X() < semi[1].X())
   {
      ordered[2]=semi[0];
      ordered[3]=semi[1];
   }
   else
   {
      ordered[2]=semi[1];
      ordered[3]=semi[0];
   }

   if(semi[2].X() < semi[3].X())
   {
      ordered[6]=semi[2];
      ordered[7]=semi[3];
   }
   else
   {
      ordered[6]=semi[3];
      ordered[7]=semi[2];
   }

   //BL quad
   semi.clear();
   for(it=(*dp).begin()+8;it<(*dp).begin()+12;it++)
   {
      std::vector<Item>::iterator sit=semi.begin();
      while((sit<semi.end())&&((*it).Y() < (*sit).Y()))
         sit++;
      semi.insert(sit,(*it));
   }

   if(semi[0].X() < semi[1].X())
   {
      ordered[8]=semi[0];
      ordered[9]=semi[1];
   }
   else
   {
      ordered[8]=semi[1];
      ordered[9]=semi[0];
   }

   if(semi[2].X() < semi[3].X())
   {
      ordered[12]=semi[2];
      ordered[13]=semi[3];
   }
   else
   {
      ordered[12]=semi[3];
      ordered[13]=semi[2];
   }

   //BR quad
   semi.clear();
   for(it=(*dp).begin()+12;it<(*dp).begin()+16;it++)
   {
      std::vector<Item>::iterator sit=semi.begin();
      while((sit<semi.end())&&((*it).Y() < (*sit).Y()))
         sit++;
      semi.insert(sit,(*it));
   }

   if(semi[0].X() < semi[1].X())
   {
      ordered[10]=semi[0];
      ordered[11]=semi[1];
   }
   else
   {
      ordered[10]=semi[1];
      ordered[11]=semi[0];
   }

   if(semi[2].X() < semi[3].X())
   {
      ordered[14]=semi[2];
      ordered[15]=semi[3];
   }
   else
   {
      ordered[14]=semi[3];
      ordered[15]=semi[2];
   }
}

//-------------------------------------------------------------------------
// Function to handle the cubic interpolation 
//-------------------------------------------------------------------------
template<class T>
void Cubic<T>::Interpolate(std::vector<Item>* dp, const unsigned int* const bands,DataAccessor<T>* lev1data)
{
   double finterppts[4]={0};
   double yinterppts[4]={0};

   double xpts[16]={0};
   double ypts[16]={0};
   double fpts[16]={0};

   std::vector<Item> ordered;
   ordered.reserve(16);
   ordered.resize(16);

   std::vector<Item> ordered_new;
   ordered_new.reserve(16);
   ordered_new.resize(16);

   //Clean the data array
   for(unsigned int b=0;b<this->length;b++)
      this->data[b]=0; 

   //check the size of the vector
   if((*dp).size()!=16)
      throw "Cubic expects 16 points to be used. Got: "+ToString((*dp).size());

   //Order the points in dp
   this->OrderPoints(dp,ordered);

   //Point to the ordered points for loop below
   std::vector<Item>* points=&ordered;

   //Vector to hold results of updated search (if required)
   std::vector<Item>* new_points=NULL;
   bool loopbreaker=false;
   bool nonewpoints=false;
   unsigned int line=0, column=0;

   for(unsigned int b=0;b<this->length;b++)
   {
      //Get the data points into an array
      for(unsigned int i=0;i<points->size();i++)
      {
         line= (*points)[i].igmrow;
         column= (*points)[i].igmcol;
         //If data not in memory then read from file
         fpts[i]=static_cast<double>(lev1data->GetData(bands[b],line,column));

         //Test if the data value is OK to use
         if((this->ignoredata==true)&&(fpts[i]==this->IGNOREVALUE))
         {
            //Update the points
            new_points=this->UpdatePoints(bands[b],lev1data);
            if(new_points==NULL)
               nonewpoints=true;
               //this->data[b]=0;
               //return;

            loopbreaker=true;
            break;
         }
      }

      if(loopbreaker==true)
      {
         if(nonewpoints==true)
         {
            //Search did not find enough points with good data
            //set this band to NODATAVALUE and move onto next one
            loopbreaker=false;
            nonewpoints=false;
            this->data[b]=this->NODATAVALUE;
            continue;
         }

         //Order this new collection of points
         this->OrderPoints(new_points,ordered_new);
         //point to the new collection of points
         points=&ordered_new;
         //subtract one off bands to trick loop into redoing this band
         b=b-1;
         //reset the bool
         loopbreaker=false;
      }
      else
      {
         //Interpolate to get the data value we want

         //Now extract the x,y values from the vector
         //and insert into arrays
         for(unsigned int i=0;i<(*points).size();i++)
         {
            xpts[i]=(*points)[i].X();
            ypts[i]=(*points)[i].Y();
         }
    
         //Get 4 values of interpolated data by interpolating across 4 splines at value x
         //xpts are the 4 x values of the data with value fpts [1 set of horiz points]
         //this will give us the interpolated function value from each spline
         GetYInterpAtX(finterppts,this->l3pos->X,xpts,fpts);

         //Get 4 y values at these points by interpolation
         //xpts,ypts the x,y of the 4 points per spline
         //this will give the y value at the interpolated point
         GetYInterpAtX(yinterppts,this->l3pos->X,xpts,ypts);

         //Sort interpolated y,f into increasing order based on y
         //to satisfy the monotonic increasing assumption
         std::vector<Pair> mypair;
         for(int i=0;i<4;i++)
         {
            mypair.push_back(Pair(yinterppts[i],finterppts[i]));
         }
         //Sort the pairs based on y
         std::sort(mypair.begin(),mypair.end(),Pair::compare);
         //extract the data again from the pairs
         for(int i=0;i<4;i++)
         {
            yinterppts[i]=mypair[i].y;
            finterppts[i]=mypair[i].f;
         }

         //Get the interpolated value at y from spline made from y values
         //of the interpolated f values - this is the interpolated point (x,y)
         this->data[b]=GetInterpValue(this->l3pos->Y,yinterppts,finterppts);

         //reset to original data for next band
         points=&ordered;
      }
   }
}

//-------------------------------------------------------------------------
// Function to search for updated points that have no ignore value in the set
//-------------------------------------------------------------------------
template<class T>
std::vector<Item>* Cubic<T>::UpdatePoints(unsigned int band,DataAccessor<T>* lev1data)
{
   return this->tg->GetQuadItems(4,this->l3pos,this->searchradius,lev1data,band,this->IGNOREVALUE);
}

#endif
