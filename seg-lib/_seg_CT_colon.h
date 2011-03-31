#pragma once

#include "_seg_common.h"
#include "_seg_BiasCorrection.h"
#include "_seg_tools.h"
#include "_seg_MRF.h"
#include "_seg_FMM.h"
#include "_seg_Topo.h"


nifti_image * CT_colon(nifti_image * CT, char* filename_out,SEG_PARAM * segment_param);
nifti_image * Copy_ShortExpec_to_Result_Lungs(PrecisionTYPE * Expec,
                                              PrecisionTYPE * BiasField,
                                              int * Short_2_Long_Indices,
                                              nifti_image * T1,
                                              char * filename,
                                              ImageSize * CurrSizes);



