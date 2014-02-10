//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef GEODESICS_H
#define GEODESICS_H

#include<cmath>
#include "conversions.h"
#include "level3grid.h"

#ifndef PI_
#define PI_
const double PI=4*atan(1.0);
#endif

void GetGeodesicDistance_Bowring(const double lon1,const double lat1,const double hei1,const double lon2,const double lat2,const double hei2,double& distance,double& azimuth,double& zenith,Ellipsoid* ell);
void GetGeodesicDistance_Vincenty(const double lon1,const double lat1,const double hei1,const double lon2,const double lat2,const double hei2,double& distance,double& azimuth,double& zenith,Ellipsoid* ell);
void GetDestinationPoint_Bowring(const double lon1,const double lat1,const double distance,const double azimuth,double &lon2,double &lat2,Ellipsoid* ell);
#endif
