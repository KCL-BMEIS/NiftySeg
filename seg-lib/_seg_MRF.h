#pragma once

#include "_seg_common.h"

void MRFregularization_mask(const PrecisionTYPE * Expec,
                       const PrecisionTYPE * G,
                       const PrecisionTYPE * H,
                       PrecisionTYPE * MRFbeta,
                       PrecisionTYPE * MRFprior,
                       PrecisionTYPE * AtlasPrior,
                       int * Long_2_Short_Indices,
                        int * Short_2_Long_Indices,
                       ImageSize * CurrSizes,
                       bool MRFflag,
                       int verbose_level);

void MRFregularization(const PrecisionTYPE * Expec,
                       const PrecisionTYPE * G,
                       const PrecisionTYPE * H,
                       PrecisionTYPE * MRFbeta,
                       PrecisionTYPE * MRFprior,
                       PrecisionTYPE * AtlasPrior,
                       ImageSize * CurrSizes,
                       bool MRFflag,
                       int verbose_level);

void MRFregularization_mask2D(const PrecisionTYPE * Expec,
                       const PrecisionTYPE * G,
                       const PrecisionTYPE * H,
                       PrecisionTYPE * MRFbeta,
                       PrecisionTYPE * MRFprior,
                       PrecisionTYPE * AtlasPrior,
                       int * Long_2_Short_Indices,
                        int * Short_2_Long_Indices,
                       ImageSize * CurrSizes,
                       bool MRFflag,
                       int verbose_level);

void MRFregularization2D(const PrecisionTYPE * Expec,
                       const PrecisionTYPE * G,
                       const PrecisionTYPE * H,
                       PrecisionTYPE * MRFbeta,
                       PrecisionTYPE * MRFprior,
                       PrecisionTYPE * AtlasPrior,
                       ImageSize * CurrSizes,
                       bool MRFflag,
                       int verbose_level);


