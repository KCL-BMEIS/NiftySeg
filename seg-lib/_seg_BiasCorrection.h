#ifndef _SEG_BIASCORRECTION_H
#define _SEG_BIASCORRECTION_H

#include <iostream>
#include <math.h>
#include "_seg_common.h"
#include "_seg_matrix.h"

#ifdef _OPENMP
#include "omp.h"
#endif



using namespace std;

inline SegPrecisionTYPE pow_int(const SegPrecisionTYPE x,
                                int exp);

void get_xyz_pow_int(SegPrecisionTYPE xpos,
                     SegPrecisionTYPE ypos,
                     SegPrecisionTYPE zpos,
                     SegPrecisionTYPE currxpower[10],
                     SegPrecisionTYPE currypower[10],
                     SegPrecisionTYPE currzpower[10],
                     int maxorder);

void BiasCorrection(SegPrecisionTYPE * BiasField,
                    SegPrecisionTYPE * BiasFieldCoefs,
                    nifti_image * T1,
                    SegPrecisionTYPE * Expec,
                    SegPrecisionTYPE * Outlierness,
                    SegPrecisionTYPE * M,
                    SegPrecisionTYPE * V,
                    int biasOrder,
                    ImageSize * CurrSizes,
                    bool flag_Bias,
                    int verbose_level);

void BiasCorrection_mask(SegPrecisionTYPE * BiasField,
                         SegPrecisionTYPE * BiasFieldCoefs,
                         nifti_image * T1,
                         int * Long_2_Short_Indices,
                         SegPrecisionTYPE * Expec,
                         SegPrecisionTYPE * Outlierness,
                         SegPrecisionTYPE * M,
                         SegPrecisionTYPE * V,
                         int biasOrder,
                         ImageSize * CurrSizes,
                         bool flag_Bias,
                         int verbose_level);

void get_xy_pow_int(SegPrecisionTYPE xpos,
                    SegPrecisionTYPE ypos,
                    SegPrecisionTYPE currxpower[10],
                    SegPrecisionTYPE currypower[10],
                    int maxorder);

void BiasCorrection2D(SegPrecisionTYPE * BiasField,
                      SegPrecisionTYPE * BiasFieldCoefs,
                      nifti_image * T1,
                      SegPrecisionTYPE * Expec,
                      SegPrecisionTYPE * Outlierness,
                      SegPrecisionTYPE * M,
                      SegPrecisionTYPE * V,
                      int biasOrder,
                      ImageSize * CurrSizes,
                      bool flag_Bias,
                      int verbose_level);

void BiasCorrection_mask2D(SegPrecisionTYPE * BiasField,
                           SegPrecisionTYPE * BiasFieldCoefs,
                           nifti_image * T1,
                           int * Long_2_Short_Indices,
                           SegPrecisionTYPE * Expec,
                           SegPrecisionTYPE * Outlierness,
                           SegPrecisionTYPE * M,
                           SegPrecisionTYPE * V,
                           int biasOrder,
                           ImageSize * CurrSizes,
                           bool flag_Bias,
                           int verbose_level);

#endif
