#pragma once

#include "_seg_common.h"

void MRFregularization_mask(const SegPrecisionTYPE * Expec,
                       const SegPrecisionTYPE * G,
                       const SegPrecisionTYPE * H,
                       SegPrecisionTYPE * MRFbeta,
                       SegPrecisionTYPE * MRFprior,
                       SegPrecisionTYPE * AtlasPrior,
                       int * Long_2_Short_Indices,
                        int * Short_2_Long_Indices,
                       ImageSize * CurrSizes,
                       bool MRFflag,
                       int verbose_level);

void MRFregularization(const SegPrecisionTYPE * Expec,
                       const SegPrecisionTYPE * G,
                       const SegPrecisionTYPE * H,
                       SegPrecisionTYPE * MRFbeta,
                       SegPrecisionTYPE * MRFprior,
                       SegPrecisionTYPE * AtlasPrior,
                       ImageSize * CurrSizes,
                       bool MRFflag,
                       int verbose_level);

void MRFregularization_mask2D(const SegPrecisionTYPE * Expec,
                       const SegPrecisionTYPE * G,
                       const SegPrecisionTYPE * H,
                       SegPrecisionTYPE * MRFbeta,
                       SegPrecisionTYPE * MRFprior,
                       SegPrecisionTYPE * AtlasPrior,
                       int * Long_2_Short_Indices,
                        int * Short_2_Long_Indices,
                       ImageSize * CurrSizes,
                       bool MRFflag,
                       int verbose_level);

void MRFregularization2D(const SegPrecisionTYPE * Expec,
                       const SegPrecisionTYPE * G,
                       const SegPrecisionTYPE * H,
                       SegPrecisionTYPE * MRFbeta,
                       SegPrecisionTYPE * MRFprior,
                       SegPrecisionTYPE * AtlasPrior,
                       ImageSize * CurrSizes,
                       bool MRFflag,
                       int verbose_level);


