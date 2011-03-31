#pragma once

#include <map>
#include "_seg_common.h"

void FMM(bool *Seeds,
         PrecisionTYPE *SpeedI,
         PrecisionTYPE * GeoTime,
         PrecisionTYPE Max,
         int * Long_2_Short_Indices,
         int * Short_2_Long_Indices,
         ImageSize * CurrSizes);

PrecisionTYPE CalcGeoTime(int index,
                          PrecisionTYPE *GeoTime,
                          PrecisionTYPE * SpeedI,
                          int * neighbour,
                          int * Short_2_Long_Indices,
                          PrecisionTYPE Max);

void TransformGeoTime(PrecisionTYPE *GeoTime,
                      PrecisionTYPE MaxGeoTime,
                      int * L2S,
                      int * S2L,
                      ImageSize * CurrSizes);
