#pragma once

#include <map>
#include "_seg_common.h"

void FMM(bool *Seeds,
         SegPrecisionTYPE *SpeedI,
         SegPrecisionTYPE * GeoTime,
         SegPrecisionTYPE Max,
         int * Long_2_Short_Indices,
         int * Short_2_Long_Indices,
         ImageSize * CurrSizes);

float * DoubleEuclideanDistance_3D(bool *Lable, float * speedptr,
                                   ImageSize * CurrSizes);

SegPrecisionTYPE CalcGeoTime(int index,
                             SegPrecisionTYPE *GeoTime,
                             SegPrecisionTYPE * SpeedI,
                             int * neighbour,
                             int * Short_2_Long_Indices,
                             SegPrecisionTYPE Max);

SegPrecisionTYPE CalcGeoTime_long(int index,
                                  SegPrecisionTYPE *GeoTime,
                                  SegPrecisionTYPE SpeedI,
                                  int * neighbour,
                                  SegPrecisionTYPE Max);

void TransformGeoTime(SegPrecisionTYPE *GeoTime,
                      SegPrecisionTYPE MaxGeoTime,
                      int * L2S,
                      int * S2L,
                      ImageSize * CurrSizes);
