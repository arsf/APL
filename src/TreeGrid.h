//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef TREEGRID_H
#define TREEGRID_H

#include <vector>
#include <list>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>

#include "basic_igm_worker.h"
#include "commonfunctions.h"
#include "treegrid_support.h"
#include "dataaccessor.h"
#include "geodesics.h"

//-------------------------------------------------------------------------
// Class to contain point data based on the X and Y locations of the point
//-------------------------------------------------------------------------
class Collection
{
public:
   Collection();
   Collection(double,double,double,double);
   virtual ~Collection();
   //function to insert an Item i into collection
   inline void Insert(Item* i); 
   void SetInfo(double cX,double cY,double sX,double sY);
   //int SaveToDisk(std::ofstream* fout);
   //void LoadFromDisk(std::ifstream* fin);

   //return the centre X and Y of the collection
   const double CentreX() const {return centreX;}
   const double CentreY() const {return centreY;}

   //Return the size of the items vector
   inline const unsigned long int GetSize()const{return items.size();}

   //Test if an area intersects with this collection
   bool Intersect(Area* area);

   std::vector<Item>* GetAllItems();
   void ClearRetItems(){retItems.clear();}

   template<class T>
   std::vector<Item>* GetNearestXItems(unsigned int num,IGMPoint* searchpoint,DataAccessor<T>* level1,unsigned int band,double IGNOREVALUE);

   virtual double GetDistance(const IGMPoint &searchpoint1,const IGMPoint &searchpoint2)
   {
      double distance=pow(searchpoint1.X-searchpoint2.X,2) + pow(searchpoint1.Y-searchpoint2.Y,2);    
      return distance;
   }

protected:
   //Container for the items
   std::vector<Item> items;
   //XY coordinates of centre point
   double centreX,centreY;
   //size of collection in X and Y (e.g. 5m x 5m)
   double sizeX,sizeY; 
   //A vector that contains nearby points to a test point which is returned
   //from GetNearestXItems
   std::vector<Item> retItems;
};

//-------------------------------------------------------------------------
// Class to contain Lat/Lon point data based on the X and Y locations of the point
// Derives from collection but calculates distances using geodesics
//-------------------------------------------------------------------------
class EllipsoidCollection : public Collection
{
public:
   EllipsoidCollection() : Collection() {ell=NULL;};
   EllipsoidCollection(double cx,double cy,double sx,double sy,Ellipsoid* ellipse) : Collection(cx,cy,sx,sy) {ell=ellipse;};
   ~EllipsoidCollection();

   double GetDistance(const IGMPoint &searchpoint1,const IGMPoint &searchpoint2)
   {
      //We assume points on ellipsoid surface so heights are = 0
      double distance=0,height1=0,height2=0;
      double azimuth=0,zenith=0;
      //Calculate geodesic distance
      GetGeodesicDistance_Bowring(searchpoint1.X*PI/180,searchpoint1.Y*PI/180,height1,searchpoint2.X*PI/180,searchpoint2.Y*PI/180,height2,distance,azimuth,zenith,ell);
      return distance*distance;//we return distance squared
   }   
private:
   Ellipsoid* ell;
};

//-------------------------------------------------------------------------
// Function to return a pointer to a vector containing the nearest N items
// in the collection to the given search X and Y.
//-------------------------------------------------------------------------
template<class T>
std::vector<Item>* Collection::GetNearestXItems(unsigned int num,IGMPoint* searchpoint,DataAccessor<T>* level1,unsigned int band,double IGNOREVALUE)
{
   if(retItems.empty()!=true)
      retItems.clear();

   std::vector<Item>::iterator it;
   std::vector<Item>::iterator ret_iter;
   double currdist=0,largestretdist=0,newlarge=0;
   unsigned int icount=0;

   if(num >= items.size())
   {
      //We want more points than are in here so return all points
      for(it=items.begin();it<items.end();it++)
      {
         if((level1==NULL)||(level1->GetData(band,it->igmrow,it->igmcol) != IGNOREVALUE))
         {
            currdist=GetDistance(*searchpoint,IGMPoint(it->X(),it->Y()));
            (*it).distance=static_cast<float>(currdist);
            retItems.push_back(*it);
            icount++;
         }
      }
      return &retItems;
   }
   else
   {
      //We need to search and return only num items
      for(it=items.begin();it<items.end();it++)
      {
         //If masking data we dont want to include points = 0
         if((level1!=NULL)&&(level1->GetData(band,it->igmrow,it->igmcol) == IGNOREVALUE))
         {
            continue; //we don't want this point as it is equal to IGNOREVALUE in the level 1
         }

         currdist=GetDistance(*searchpoint,IGMPoint(it->X(),it->Y()));
         if((currdist < largestretdist)||(icount<num))
         {
            if(icount != num)
            {
               (*(it)).distance=static_cast<float>(currdist);
               retItems.push_back(*it);
               icount++;

               if(largestretdist<currdist)
                  largestretdist=currdist;

               //Sort the vector here if it contains the required number of items
               if(icount==num)
                  std::sort(retItems.begin(),retItems.end());

               //break out here since we dont need do the next stuff
               continue;
            }
            else
            {
               //remove the last item as this is the largest
               retItems.pop_back();
               //If there are still items in the vector - get the largest distance
               if(retItems.size()>0)
                  newlarge=retItems.back().distance;
               //assign the current distance to the new item
               (*(it)).distance=static_cast<float>(currdist);
               //If the current distance is larger than the largest distance
               //insert the new item at the end of the array
               if(currdist > newlarge)
                  retItems.push_back(*it);
               else
               {
                  //else loop through the vector and insert at the correct place
                  for(ret_iter=retItems.begin();ret_iter<retItems.end();ret_iter++)
                  {
                     if((*ret_iter).distance > currdist)
                     {
                        retItems.insert(ret_iter,(*it));
                        break;
                     }
                  }
               }
               //Recalculate the largest return distance
               if(newlarge < currdist)
                  largestretdist=currdist;
               else
                  largestretdist=newlarge;
            }
         }
      }
   }

   return &retItems;
}

//-------------------------------------------------------------------------
// Class to contain the point data from the data* array in collections determined
// by X and Y position of point data
//-------------------------------------------------------------------------
class TreeGrid
{
public:
   TreeGrid();
   TreeGrid(unsigned int r,unsigned int c,Ellipsoid* ell=NULL);
   virtual ~TreeGrid();
   void SetUpGrid(double sX,double sY,double tlX,double tlY);
   void SetUpGrid(unsigned int r, unsigned int c,double sX,double sY,double tlX,double tlY,double brX,double brY);
   //Insert data into the treegrid
   void InsertData(Item* data,unsigned int nitems);
   //Return a pointer to a specific collection
   Collection* GetCollection(unsigned int r,unsigned int c);
   //Functions for import / output of the TreeGrid
   //void SaveToDisk(std::string fname);
   //void LoadFromDisk(std::string fname);

   const unsigned long int NumRows() const {return rows;}
   const unsigned long int NumCols() const {return cols;}
   const double TopLeftX() const {return topLeftX;}
   const double TopLeftY() const {return topLeftY;}
   const double BottomRightX() const {return bottomRightX;}
   const double BottomRightY() const {return bottomRightY;}
   const double SizeX() const {return sizeX;}
   const double SizeY() const {return sizeY;}
   bool GetAllCollectionsWithinRadius(std::vector<Collection*>* colls,const IGMPoint* const searchpoint,double searchradius);
   ItemData itemdata;

   //Return a pointer to a vector that contains the nearest num points to position of searchpoint
   template<class T>
   std::vector<Item>* GetNearestXItems(unsigned int num,IGMPoint* searchpoint,double searchradius,DataAccessor<T>* level1,unsigned int band,double IGNOREVALUE);
   template<class T>
   std::vector<Item>* GetQuadItems(unsigned int npoints,const IGMPoint* searchpoint,double searchradius,DataAccessor<T>* level1,unsigned int band,double IGNOREVALUE);
   
   bool IsGeographic(){return islatlon;}
protected:
   //size of grid
   unsigned long int rows,cols; 
   //location of top left corner
   double topLeftX,topLeftY; 
   //location of bottom right corner of data, that is, the minY,maxX (not necessarily tree grid coverage)
   //This has been included to calculate the number of rows of output level 3 image data correctly 
   double bottomRightX,bottomRightY;

   //size of grid cells
   double sizeX,sizeY;
   //collection array which holds the data
   Collection*** collection; 
   //pointer to a null collection i.e. one that exists but is empty
   Collection* NullCollection;

   std::vector<Item> retItems;

   Collection* InitialiseNewCollection(unsigned int r, unsigned int c);
   Ellipsoid* ellipse;
   double upperdateline,lowerdateline;
   bool islatlon;
   void CheckSearchBoxForWraps(std::list<Area> &search_areas);
};


//-------------------------------------------------------------------------
// Function to return a pointer to a vector containing the nearest N
// items to the given X,Y position
//-------------------------------------------------------------------------
template<class T>
std::vector<Item>* TreeGrid::GetNearestXItems(unsigned int num,IGMPoint* searchpoint,double searchradius,DataAccessor<T>* level1,unsigned int band,double IGNOREVALUE)
{
   //Identify the collection to search
   unsigned long int r=0,c=0,size=0;
   std::vector<Item>::iterator iter;
   std::vector<Item>* retvec=NULL;   
   r=static_cast<unsigned long int>(floor((topLeftY-searchpoint->Y)/sizeY));
   c=static_cast<unsigned long int>(floor((searchpoint->X-topLeftX)/sizeX));
   //Clear retItems - the vector that a pointer to is returned
   retItems.clear();

   //Check that the colection in which the search point lies exists (the search point may be outside the treegrid)
   if((r>=rows)||(c>=cols))
   {
      //This collection does not exist - but that does not mean we shouldn't search to see if any do exist within
      //the search radius of the search point
      std::vector<Item>* retvec_additional=NULL;
      std::vector<Collection*> colls;
      std::vector<Collection*>::iterator coll_it;
      //Get all collections within search radius
      GetAllCollectionsWithinRadius(&colls,searchpoint,searchradius);
      coll_it=colls.begin();

      while((coll_it<colls.end())) //This will search all collections - could be made better if you know you have searched all nearest ones
      {
         if((*coll_it)->GetSize()!=0)
         {
            retvec_additional=(*coll_it)->GetNearestXItems(num,searchpoint,level1,band,IGNOREVALUE);     
            //Now insert these new values into retItems
            iter=retItems.begin();
            retItems.insert(iter,retvec_additional->begin(),retvec_additional->end());
            //update size variable
            size=retItems.size();            
         }
         //Point to next collection
         coll_it++;
      }
   }
   else
   {
      //Search the current collection, as it does exist, if it is not empty
      if(collection[r][c]->GetSize()!=0)
      {
         //The collection contains some data - lets search it
         retvec=collection[r][c]->GetNearestXItems(num,searchpoint,level1,band,IGNOREVALUE);
         //Insert the found items into the return vector
         iter=retItems.begin();
         retItems.insert(iter,retvec->begin(),retvec->end());
         size=retItems.size();
      }
      else //it is empty - so we have no items as of yet
         size=0;

      //We need to check we have num items before we return - so if we have less
      //we need to check surrounding containers too
      if(size==num)
      {
         //We should now check that furthest point is closer to the search point
         //than the edges of the collection - if not we need to search these neighbouring
         //collections too
         //Get the furthest Item
         Item furthest=retItems.back();
         //Get all collections within furthest distance (note sqrt because distance is squared)
         std::vector<Collection*> colls;

         if((!GetAllCollectionsWithinRadius(&colls,searchpoint,sqrt(furthest.distance))) || (colls.size()==0)) //there are no intersects
         {
            std::sort(retItems.begin(),retItems.end());
            return &retItems;
         }
         else
         {  
            //Collections were found with potentially nearer points
            std::vector<Collection*>::iterator it;
            //Iterate through the collections and get the nearest n points from each one
            //Then before returning - order them and remove the last (total - n) points
            //The first element of colls is the collection we have just searched so skip it
            for(it=colls.begin()+1;it<colls.end();it++)
            {
               retvec=(*it)->GetNearestXItems(num,searchpoint,level1,band,IGNOREVALUE);
               iter=retItems.begin();
               retItems.insert(iter,retvec->begin(),retvec->end());
               size=retItems.size();
            }
         }
         //Clear the colls vector
         colls.clear();
      }
      else
      {
         //Haven't found enough nearby points yet so search
         //neighbouring collections 
         std::vector<Item>* retvec_additional=NULL;
         //Get all collections within search radius
         std::vector<Collection*> colls;
         GetAllCollectionsWithinRadius(&colls,searchpoint,searchradius);
         std::vector<Collection*>::iterator coll_it;
         //we've already searched first collection so point to second collection 
         coll_it=colls.begin()+1;

         while((coll_it<colls.end())) //This will search all collections - could be made better if you know you have searched all nearest ones
         {
            if((*coll_it)->GetSize()!=0)
            {
               retvec_additional=(*coll_it)->GetNearestXItems(num,searchpoint,level1,band,IGNOREVALUE);     
               //Now insert these new values into retItems
               iter=retItems.begin();
               retItems.insert(iter,retvec_additional->begin(),retvec_additional->end());
               //update size variable
               size=retItems.size();            
            }
            //Point to next collection
            coll_it++;
         }
      }  
   }

   //Sort the vector to get in order
   std::sort(retItems.begin(),retItems.end());

   if(retItems.size()>num)
   {
      //Remove any points past num - we only want num points
      retItems.erase(retItems.begin()+num,retItems.end());
   }

   //Now prune off any that are further away than the search radius
   for(std::vector<Item>::iterator it=retItems.begin();it<retItems.end();it++)
   {
      if((*it).distance > searchradius*searchradius)
      {
         retItems.erase(it,retItems.end());
         break;
      }
   }

   if(retItems.size()==0)
      return NULL;
   
   return &retItems;
}



//-------------------------------------------------------------------------
// Function to return a pointer to a vector containing the nearest 4
// items to the given X,Y position that form a quad around the XY position
//-------------------------------------------------------------------------
template<class T>
std::vector<Item>* TreeGrid::GetQuadItems(unsigned int npoints,const IGMPoint* searchpoint,const double searchradius,DataAccessor<T>* level1,unsigned int band,double IGNOREVALUE)
{
   //Identify the collection to search
   unsigned long int r=0,c=0;
   std::vector<Item>::iterator iter;
   r=static_cast<unsigned long int>(floor((topLeftY-searchpoint->Y)/sizeY));
   c=static_cast<unsigned long int>(floor((searchpoint->X-topLeftX)/sizeX));

   if((r>=rows)||(c>=cols))
   {
      //This collection does not exist
      Logger::Debug("Collection does not exist at row,col: "+ToString(r)+" "+ToString(c));
      return NULL;
   }

   //Clear retItems - the vector that a pointer to is returned
   retItems.clear();

   //MAY NEED TO ALLOCATE A VECTOR OF SIZE 4*npoints HERE
   retItems.resize(4*npoints);

   //Set position pointers
   int URquad=1;
   int BRquad=3;
   int BLquad=2;
   int ULquad=0;

   //bool UR=false, BR=false, UL=false, BL=false;
   unsigned int URsum=0,BRsum=0, ULsum=0,BLsum=0;

   std::vector<Item>* retvec=NULL;
   std::vector<Collection*> colls;

   double dx=0,dy=0;

   //Get all the collections within the search radius
   if(!GetAllCollectionsWithinRadius(&colls,searchpoint,searchradius))
   {
      //No collection exists containing the search point
      return NULL;
   }

   std::vector<Collection*>::iterator coll_it=colls.begin();
   std::vector<Collection*> checked;
   checked.reserve(colls.size());

   //Bool to check if we need to refine the search area
   bool refined_search=false;

   //Check if there are 4*npoints pts (npoints in each quadrant)
   while(coll_it<colls.end())
   {
      //Get the items from this collection
      if((*coll_it)->GetSize()!=0)
      {
         retvec=(*coll_it)->GetAllItems();
      }
      else
      {
         checked.push_back(*coll_it);
         coll_it++;
         continue;
      }

      //loop through all the items and check if they are suitable for returning
      for(std::vector<Item>::iterator it=retvec->begin();it<retvec->end();it++)
      {
         if((level1!=NULL)&&(level1->GetData(band,it->igmrow,it->igmcol)==IGNOREVALUE))
            continue; //we don't want this point as it has value IGNOREVALUE

         dx=it->X()-searchpoint->X;
         dy=it->Y()-searchpoint->Y;
         (*it).distance=static_cast<float>((*coll_it)->GetDistance(*searchpoint,IGMPoint(it->X(),it->Y())));
         if(dx>=0)
         {
            if(dy>=0)
            {
               if(URsum==npoints)
               {
                  //Check if this point is closer than the others in this quadrant
                  Item* furthest=std::max_element(&retItems[URquad*npoints],&retItems[URquad*npoints+npoints]);

                  if((*it).distance < furthest->distance)
                  {
                     //Update the current value                  
                     *furthest=(*it);
                  }
               }
               else
               {
                  retItems[URquad*npoints+URsum]=(*it);                     
                  URsum++;
               }
            }
            else
            {
               if(BRsum==npoints)
               {
                  //Check if this point is closer than the others in this quadrant
                  Item* furthest=std::max_element(&retItems[BRquad*npoints],&retItems[BRquad*npoints+npoints]);

                  if((*it).distance < furthest->distance)
                  {
                     //Update the current value                  
                     *furthest=(*it);
                  }
               }
               else
               {
                  retItems[BRquad*npoints+BRsum]=(*it);
                  BRsum++;
               }
            }
         }
         else
         {
            if(dy>=0)
            {
               if(ULsum==npoints)
               {
                  //Check if this point is closer than the others in this quadrant
                  Item* furthest=std::max_element(&retItems[ULquad*npoints],&retItems[ULquad*npoints+npoints]);

                  if((*it).distance < furthest->distance)
                  {
                     //Update the current value                  
                     *furthest=(*it);
                  }
               }
               else
               {
                  retItems[ULquad*npoints+ULsum]=(*it);     
                  ULsum++;
               }
            }
            else
            {
               if(BLsum==npoints)
               {
                  //Check if this point is closer than the others in this quadrant
                  Item* furthest=std::max_element(&retItems[BLquad*npoints],&retItems[BLquad*npoints+npoints]);

                  if((*it).distance < furthest->distance)
                  {
                     //Update the current value                  
                     *furthest=(*it);
                  }
               }
               else
               {
                  retItems[BLquad*npoints+BLsum]=(*it);  
                  BLsum++;
               }
            }            
         }
      }

      //Finished checking this collection - add it to the checked list
      checked.push_back(*coll_it);

      //Clear the retItems vector of this collection
      (*coll_it)->ClearRetItems();
      retvec=NULL;

      //If we have npoints in each quadrant then we only need search collections within
      //the furthest distance from the point XY. Else we need to carry on searching colections within search radius
      if((refined_search==false)&&((URsum==npoints)&&(ULsum==npoints)&&(BRsum==npoints)&&(BLsum==npoints)))
      {
         float furthest_distance=(*retItems.begin()).distance;
         for(iter=retItems.begin()+1;iter<retItems.end();iter++)
         {
            if((*iter).distance > furthest_distance)
               furthest_distance = (*iter).distance;
         }
         furthest_distance=sqrt(furthest_distance); //as distance is squared remember
         //Find all collections within this distance from the search point
         GetAllCollectionsWithinRadius(&colls,searchpoint,furthest_distance);
         //remove any collections already searched
         for(std::vector<Collection*>::iterator it=checked.begin();it<checked.end();it++)
         {
            for(std::vector<Collection*>::iterator it2=colls.begin();it2<colls.end();it2++)
            {
               //if this exists in the checked vector and the new vector
               if((*it)==(*it2))
               {
                  //this collection has already been checked so remove it from the new vector
                  colls.erase(it2);
                  break;
               }
            }
         }

         //reset the pointer to the start of the new collection
         coll_it=colls.begin();
         refined_search=true;
      }
      else
         coll_it++;

   }
   //Exit of while loop - we either have all 4*npoints or there aren't 4*npoints in the search region

   //Clean up
   if(!colls.empty())
      colls.clear();
   if(!checked.empty())
      checked.clear();

   if((URsum==npoints)&&(ULsum==npoints)&&(BRsum==npoints)&&(BLsum==npoints))
      return &retItems;
   else
      return NULL;
}






//-------------------------------------------------------------------------
// Class to contain the point data from an IGM file in collections determined
// by X and Y position of point data, derived from TreeGrid
//-------------------------------------------------------------------------
class IGMTreeGrid : public TreeGrid
{
public:
   IGMTreeGrid(std::string fname,std::vector<unsigned int> dropscanvector,Area* region);
   ~IGMTreeGrid();

   void InsertData(std::vector<unsigned int> dropscanvector,Area* region);
   std::string GetMapInfo();
   void GetAveragePixelSeparation(double &x,double &y);
   void GetAveragePixelSeparationMetres(double &x,double &y);
private:
   Basic_IGM_Worker* igm;
   double pixsepX,pixsepY;
};




#endif
