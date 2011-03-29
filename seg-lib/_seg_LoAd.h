#pragma once

#include "_seg_common.h"
#include "_seg_BiasCorrection.h"
#include "_seg_tools.h"
#include "_seg_MRF.h"
#include "_seg_FMM.h"


nifti_image * LoAd_Segment(nifti_image * T1,
                           nifti_image * Mask,
                           nifti_image * Priors,
                           SEG_PARAM * segment_param);
