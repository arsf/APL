//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

#ifndef TRANSFORMATIONS_H
#define TRANSFORMATIONS_H

#include <blitz/array.h>
#include "cartesianvector.h"

enum rotation_order {RXYZ, RXZY, RYXZ, RYZX, RZXY, RZYX};

//-------------------------------------------------------------------------
// Get the vector V into ECEF Cartesian coordinates [COMBINED method]
//-------------------------------------------------------------------------
void GetVVinECEFXYZ(blitz::TinyMatrix<double,3,1>* const V,double* const ECEFXYZ,const double lat, const double lon,
                     const double roll,const double pitch,const double heading);

//-------------------------------------------------------------------------
// Get the vector V into ECEF Cartesian coordinates [SPLIT method]
//-------------------------------------------------------------------------
void GetVVinECEFXYZ(blitz::TinyMatrix<double,3,1>* const V,double* const ECEFXYZ,const double lat, const double lon,
                     const double theta,const double phi,const double kappa,const double roll,const double pitch,const double heading);

//-------------------------------------------------------------------------
// Return a 3D rotation matrix for the given angles 
//-------------------------------------------------------------------------
blitz::TinyMatrix<double,3,3> Create3DRotMatrix(const double rx, const double ry, const double rz,const rotation_order orderflag);

#endif 
