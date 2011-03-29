#ifndef _SEG_BIASCORRECTION_H
#define _SEG_BIASCORRECTION_H

#include <iostream>
#include <math.h>
#include "_seg_common.h"
#include "_seg_matrix.h"



using namespace std;

inline PrecisionTYPE pow_int(const PrecisionTYPE x,
                     int exp);

void get_xyz_pow_int(PrecisionTYPE xpos,
                     PrecisionTYPE ypos,
                     PrecisionTYPE zpos,
                     PrecisionTYPE currxpower[10],
                     PrecisionTYPE currypower[10],
                     PrecisionTYPE currzpower[10],
                     int maxorder);

void BiasCorrection(PrecisionTYPE * BiasField,
                    PrecisionTYPE * BiasFieldCoefs,
                    nifti_image * T1,
                    PrecisionTYPE * Expec,
                    PrecisionTYPE * M,
                    PrecisionTYPE * V,
                    int biasOrder,
                    ImageSize * CurrSizes,
                    bool flag_Bias,
                    int verbose_level);

void BiasCorrection_mask(PrecisionTYPE * BiasField,
                    PrecisionTYPE * BiasFieldCoefs,
                    nifti_image * T1,
                    int * Long_2_Short_Indices,
                    PrecisionTYPE * Expec,
                    PrecisionTYPE * M,
                    PrecisionTYPE * V,
                    int biasOrder,
                    ImageSize * CurrSizes,
                    bool flag_Bias,
                    int verbose_level);

void get_xy_pow_int(PrecisionTYPE xpos,
                     PrecisionTYPE ypos,
                     PrecisionTYPE currxpower[10],
                     PrecisionTYPE currypower[10],
                     int maxorder);

void BiasCorrection2D(PrecisionTYPE * BiasField,
                    PrecisionTYPE * BiasFieldCoefs,
                    nifti_image * T1,
                    PrecisionTYPE * Expec,
                    PrecisionTYPE * M,
                    PrecisionTYPE * V,
                    int biasOrder,
                    ImageSize * CurrSizes,
                    bool flag_Bias,
                    int verbose_level);

void BiasCorrection_mask2D(PrecisionTYPE * BiasField,
                    PrecisionTYPE * BiasFieldCoefs,
                    nifti_image * T1,
                    int * Long_2_Short_Indices,
                    PrecisionTYPE * Expec,
                    PrecisionTYPE * M,
                    PrecisionTYPE * V,
                    int biasOrder,
                    ImageSize * CurrSizes,
                    bool flag_Bias,
                    int verbose_level);

#endif
