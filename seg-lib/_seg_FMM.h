#pragma once

#include <map>
#include "_seg_common.h"

void FMM(bool *Seeds,
         segPrecisionTYPE *SpeedI,
         segPrecisionTYPE * GeoTime,
         segPrecisionTYPE Max,
         int * Long_2_Short_Indices,
         int * Short_2_Long_Indices,
         ImageSize * CurrSizes);

float * DoubleEuclideanDistance_3D(bool *Lable, float * speedptr,
                                   ImageSize * CurrSizes);

segPrecisionTYPE CalcGeoTime(int index,
                             segPrecisionTYPE *GeoTime,
                             segPrecisionTYPE * SpeedI,
                             int * neighbour,
                             int * Short_2_Long_Indices,
                             segPrecisionTYPE Max);

segPrecisionTYPE CalcGeoTime_long(int index,
                                  segPrecisionTYPE *GeoTime,
                                  segPrecisionTYPE SpeedI,
                                  int * neighbour,
                                  segPrecisionTYPE Max);

void TransformGeoTime(segPrecisionTYPE *GeoTime,
                      segPrecisionTYPE MaxGeoTime,
                      int * L2S,
                      int * S2L,
                      ImageSize * CurrSizes);
