#pragma once

#include "NiftySegWinExportHeader.h"
#include <map>
#include "_seg_common.h"

NIFTYSEG_WINEXPORT void FMM(bool *Seeds,
         segPrecisionTYPE *SpeedI,
         segPrecisionTYPE * GeoTime,
         segPrecisionTYPE Max,
         int * Long_2_Short_Indices,
         int * Short_2_Long_Indices,
         ImageSize * CurrSizes);

NIFTYSEG_WINEXPORT float * DoubleEuclideanDistance_3D(bool *Lable, float * speedptr,
                                   ImageSize * CurrSizes);

NIFTYSEG_WINEXPORT segPrecisionTYPE CalcGeoTime(int index,
                             segPrecisionTYPE *GeoTime,
                             segPrecisionTYPE * SpeedI,
                             int * neighbour,
                             int * Short_2_Long_Indices,
                             segPrecisionTYPE Max);

NIFTYSEG_WINEXPORT segPrecisionTYPE CalcGeoTime_long(int index,
                                  segPrecisionTYPE *GeoTime,
                                  segPrecisionTYPE SpeedI,
                                  int * neighbour,
                                  segPrecisionTYPE Max);

NIFTYSEG_WINEXPORT void TransformGeoTime(segPrecisionTYPE *GeoTime,
                      segPrecisionTYPE MaxGeoTime,
                      int * L2S,
                      int * S2L,
                      ImageSize * CurrSizes);
