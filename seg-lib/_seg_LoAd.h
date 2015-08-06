#pragma once

#include "NiftySegWinExportHeader.h"
#include <iostream>
#include <cmath>
#include "_seg_common.h"
#include "_seg_tools.h"
#include "_seg_matrix.h"
#include "_seg_FMM.h"

#ifdef _OPENMP
#include "omp.h"
#endif


NIFTYSEG_WINEXPORT int calcM_mask(nifti_image * T1,
               segPrecisionTYPE * Expec,
               segPrecisionTYPE * BiasField,
               segPrecisionTYPE * Outlierness,
               int * S2L,
               segPrecisionTYPE * M,
               segPrecisionTYPE * V,
               segPrecisionTYPE * M_MAP,
               segPrecisionTYPE * V_MAP,
               segPrecisionTYPE reg_factor,
               ImageSize * CurrSizes,
               int verbose);

NIFTYSEG_WINEXPORT int calcM_mask_LoAd(nifti_image * T1,
                    segPrecisionTYPE * Expec,
                    segPrecisionTYPE * BiasField,
                    int * S2L,
                    segPrecisionTYPE * M,
                    segPrecisionTYPE * V,
                    ImageSize * CurrSizes,
                    int verbose,
                    bool PVon);

NIFTYSEG_WINEXPORT int calcM(nifti_image * T1,
          segPrecisionTYPE * Expec,
          segPrecisionTYPE * BiasField,
          segPrecisionTYPE * Outlierness,
          segPrecisionTYPE * M,
          segPrecisionTYPE * V,
          segPrecisionTYPE * M_MAP,
          segPrecisionTYPE * V_MAP,
          segPrecisionTYPE reg_factor,
          ImageSize * CurrSizes,
          int verbose);

NIFTYSEG_WINEXPORT int calcE(nifti_image * T1,
          segPrecisionTYPE * MRF,
          segPrecisionTYPE * Expec,
          double * loglik,
          segPrecisionTYPE * BiasField,
          segPrecisionTYPE * Outlierness,
          segPrecisionTYPE OutliernessThreshold,
          segPrecisionTYPE * M,
          segPrecisionTYPE * V,
          ImageSize * CurrSizes,
          int verbose);

NIFTYSEG_WINEXPORT int calcE_mask(nifti_image * T1,
               segPrecisionTYPE * IterPrior,
               segPrecisionTYPE * Expec,
               double * loglik,
               segPrecisionTYPE * BiasField,
               segPrecisionTYPE * Outlierness,
               segPrecisionTYPE OutliernessThreshold,
               int * S2L,
               segPrecisionTYPE * M,
               segPrecisionTYPE * V,
               ImageSize * CurrSizes,
               int verbose);
/*
int calcE_aprox(nifti_image * T1,
                 SegPrecisionTYPE * MRF,
                 SegPrecisionTYPE * Expec,
                 SegPrecisionTYPE * loglik,
                 SegPrecisionTYPE * BiasField,
                 SegPrecisionTYPE * M,
                 SegPrecisionTYPE * V,
                 ImageSize * CurrSizes,
                 int verbose);

int calcE_mask_aprox(nifti_image * T1,
                      SegPrecisionTYPE * IterPrior,
                      SegPrecisionTYPE * Expec,
                      SegPrecisionTYPE * loglik,
                      SegPrecisionTYPE * BiasField,
                      int * S2L,
                      SegPrecisionTYPE * M,
                      SegPrecisionTYPE * V,
                      ImageSize * CurrSizes,
                      int verbose);
*/
NIFTYSEG_WINEXPORT int Relax_Priors(segPrecisionTYPE * Priors,
                 segPrecisionTYPE * Expec,
                 segPrecisionTYPE * MRF,
                 int * S2L,
                 int * L2S,
                 float RelaxFactor,
                 segPrecisionTYPE * G,
                 segPrecisionTYPE ba,
                 segPrecisionTYPE be,
                 ImageSize * CurrSizes,
                 seg_EM_Params * segment_param);

NIFTYSEG_WINEXPORT int Normalize_Image(nifti_image * T1,
                    ImageSize * CurrSizes,
                    bool verbose);

NIFTYSEG_WINEXPORT int Normalize_Image_mask(nifti_image * T1,
                         nifti_image * Mask,
                         ImageSize * CurrSizes,
                         bool verbose);

NIFTYSEG_WINEXPORT int RelaxPriors( segPrecisionTYPE * Expec,
                 segPrecisionTYPE * MRF,
                 int * S2L ,int * L2S ,
                 ImageSize CurrSizes,
                 int class_with_CSF);

NIFTYSEG_WINEXPORT int Convert_to_PV(nifti_image * T1,
                  segPrecisionTYPE * BiasField,
                  segPrecisionTYPE * ShortPrior,
                  segPrecisionTYPE * Expec,
                  segPrecisionTYPE * MRF,
                  segPrecisionTYPE * M,
                  segPrecisionTYPE * V,
                  int * S2L,
                  int * L2S,
                  ImageSize * CurrSizes,
                  seg_EM_Params * segment_param);

NIFTYSEG_WINEXPORT int Create_GH_5class(segPrecisionTYPE * G,
                     segPrecisionTYPE * H,
                     segPrecisionTYPE ba,
                     segPrecisionTYPE be,
                     segPrecisionTYPE ratio,
                     seg_EM_Params * segment_param);

NIFTYSEG_WINEXPORT int Create_GH_7class(segPrecisionTYPE * G,
                     segPrecisionTYPE * H,
                     segPrecisionTYPE ba,
                     segPrecisionTYPE be,
                     segPrecisionTYPE ratio,
                     seg_EM_Params * segment_param);

NIFTYSEG_WINEXPORT int Relax_Priors_Share(segPrecisionTYPE * Priors,
                       segPrecisionTYPE * Expec,
                       float RelaxFactor,
                       segPrecisionTYPE * G,
                       segPrecisionTYPE be,
                       ImageSize * CurrSizes);



NIFTYSEG_WINEXPORT int Normalize_T1_and_MV(nifti_image * T1,
                        nifti_image * Mask,
                        segPrecisionTYPE * M,
                        segPrecisionTYPE * V,
                        ImageSize * CurrSizes);

NIFTYSEG_WINEXPORT int Normalize_NaN_Priors(nifti_image * Priors,
                         bool verbose);

NIFTYSEG_WINEXPORT int Normalize_NaN_Priors_mask(nifti_image * Priors,
                              nifti_image * Mask,
                              bool verbose);

NIFTYSEG_WINEXPORT int Convert_WM_and_GM_to_PV(nifti_image * T1,
                            segPrecisionTYPE * BiasField,
                            segPrecisionTYPE * ShortPrior,
                            segPrecisionTYPE * Expec,
                            int * S2L,
                            segPrecisionTYPE * M,
                            segPrecisionTYPE * V,
                            ImageSize * CurrSize);

NIFTYSEG_WINEXPORT nifti_image * LoAd_Segment(nifti_image * T1,
                           nifti_image * Mask,
                           nifti_image * Priors,
                           seg_EM_Params * segment_param);


NIFTYSEG_WINEXPORT int Create_diagonal_GH_Nclass(segPrecisionTYPE * G,
                              segPrecisionTYPE * H,
                              segPrecisionTYPE ratio,
                              seg_EM_Params * segment_param);

NIFTYSEG_WINEXPORT int Sulci_and_gyri_correction(segPrecisionTYPE * MRF_Beta,
                              segPrecisionTYPE * ShortPrior,
                              segPrecisionTYPE * Expec,
                              segPrecisionTYPE *MRF,
                              int * S2L,
                              int * L2S,
                              ImageSize *CurrSizes);

NIFTYSEG_WINEXPORT nifti_image * Copy_ShortExpec_to_Result(nifti_image * T1,
                                        segPrecisionTYPE * Expec,
                                        segPrecisionTYPE * BiasField,
                                        segPrecisionTYPE * BiasFieldCoefs,
                                        int * S2L,
                                        nifti_image * Priors,
                                        seg_EM_Params * segment_param,
                                        segPrecisionTYPE * M,
                                        ImageSize * CurrSizes);

NIFTYSEG_WINEXPORT int * Create_Long_2_Short_Matrix_from_NII(nifti_image * Mask);

NIFTYSEG_WINEXPORT int * Create_Short_2_Long_Matrix_from_NII(nifti_image * Mask,
                                          long *shortsize);

NIFTYSEG_WINEXPORT segPrecisionTYPE * Create_cArray_from_Prior_mask(nifti_image * Mask,
                                                 nifti_image * Priors,
                                                 long numclass,
                                                 bool PV_ON);

NIFTYSEG_WINEXPORT segPrecisionTYPE * Create_cArray_from_Prior(nifti_image * Priors,
                                            long numclass,
                                            bool PV_ON);

NIFTYSEG_WINEXPORT segPrecisionTYPE * Create_cArray_from_3D_image(nifti_image * Mask,
                                               nifti_image * SourceImage);

NIFTYSEG_WINEXPORT int * Create_Short_2_Long_Matrix_from_Carray(bool * Mask,
                                             int * shortsize,
                                             int nvox);

NIFTYSEG_WINEXPORT int *  Create_Long_2_Short_Matrix_from_Carray(bool * Mask,
                                              int nvox);

NIFTYSEG_WINEXPORT nifti_image * Copy_Single_ShortImage_to_Result(segPrecisionTYPE * SingleImage,
                                               int * Short_2_Long_Indices,
                                               nifti_image * Priors,
                                               char * filename,
                                               ImageSize * CurrSizes);


NIFTYSEG_WINEXPORT nifti_image * Copy_Expec_and_BiasCorrected_to_Result_mask(segPrecisionTYPE * Expec,
                                                          segPrecisionTYPE * BiasField,
                                                          int * Short_2_Long_Indices,
                                                          nifti_image * T1,
                                                          char * filename,
                                                          ImageSize * CurrSizes);
NIFTYSEG_WINEXPORT nifti_image * Copy_Expec_and_BiasCorrected_to_Result(segPrecisionTYPE * Expec,
                                                     segPrecisionTYPE * BiasField,
                                                     nifti_image * T1,
                                                     char * filename,
                                                     ImageSize * CurrSizes);

NIFTYSEG_WINEXPORT nifti_image * Copy_Expec_to_Result_mask(segPrecisionTYPE * Expec,
                                        int * Short_2_Long_Indices,
                                        nifti_image * T1,
                                        char * filename,
                                        ImageSize * CurrSizes);

NIFTYSEG_WINEXPORT nifti_image * Copy_Expec_to_Result_Neonate_mask(segPrecisionTYPE * Expec,
                                                int * Short_2_Long_Indices,
                                                int * Long_2_Short_Indices,
                                                nifti_image * T1,
                                                float * Biasfield,
                                                float * M,
                                                char * filename,
                                                ImageSize * CurrSizes);


NIFTYSEG_WINEXPORT nifti_image * Copy_Expec_to_Result(segPrecisionTYPE * Expec,
                                   nifti_image * T1,
                                   char * filename,
                                   ImageSize * CurrSizes);

NIFTYSEG_WINEXPORT int PriorWeight_mask(float * ShortPrior,
                     nifti_image * Priors,
                     float * Expec,
                     float GaussKernelSize,
                     float RelaxFactor,
                     int * S2L,
                     int * L2S,
                     ImageSize * CurrSizes,
                     int verbose_level);

NIFTYSEG_WINEXPORT nifti_image * Get_Bias_Corrected(float * BiasField, nifti_image * T1, char * filename, ImageSize * CurrSizes);

NIFTYSEG_WINEXPORT nifti_image * Get_Bias_Corrected_mask(float * BiasFieldCoefs, nifti_image * T1, nifti_image *Mask, char * filename, ImageSize * CurrSizes, int biasOrder);

NIFTYSEG_WINEXPORT int printloglik(int iter,
                segPrecisionTYPE loglik,
                segPrecisionTYPE oldloglik);


NIFTYSEG_WINEXPORT bool * binarise_image(segPrecisionTYPE * SingleImage,
                      segPrecisionTYPE Threshold,
                      ImageSize * CurrSizes);


NIFTYSEG_WINEXPORT void MRFregularization_mask(const segPrecisionTYPE * Expec,
                            const segPrecisionTYPE * G,
                            const segPrecisionTYPE * H,
                            segPrecisionTYPE * MRFbeta,
                            segPrecisionTYPE * MRFprior,
                            segPrecisionTYPE * AtlasPrior,
                            int * Long_2_Short_Indices,
                            int * Short_2_Long_Indices,
                            ImageSize * CurrSizes,
                            bool MRFflag,
                            int verbose_level);

NIFTYSEG_WINEXPORT void MRFregularization(const segPrecisionTYPE * Expec,
                       const segPrecisionTYPE * G,
                       const segPrecisionTYPE * H,
                       segPrecisionTYPE * MRFbeta,
                       segPrecisionTYPE * MRFprior,
                       segPrecisionTYPE * AtlasPrior,
                       ImageSize * CurrSizes,
                       bool MRFflag,
                       int verbose_level);

NIFTYSEG_WINEXPORT void MRFregularization_mask2D(const segPrecisionTYPE * Expec,
                              const segPrecisionTYPE * G,
                              const segPrecisionTYPE * H,
                              segPrecisionTYPE * MRFbeta,
                              segPrecisionTYPE * MRFprior,
                              segPrecisionTYPE * AtlasPrior,
                              int * Long_2_Short_Indices,
                              int * Short_2_Long_Indices,
                              ImageSize * CurrSizes,
                              bool MRFflag,
                              int verbose_level);

NIFTYSEG_WINEXPORT void MRFregularization2D(const segPrecisionTYPE * Expec,
                         const segPrecisionTYPE * G,
                         const segPrecisionTYPE * H,
                         segPrecisionTYPE * MRFbeta,
                         segPrecisionTYPE * MRFprior,
                         segPrecisionTYPE * AtlasPrior,
                         ImageSize * CurrSizes,
                         bool MRFflag,
                         int verbose_level);

NIFTYSEG_WINEXPORT void BiasCorrection_SPARCS(float * BiasField,
                           float * T1,
                           float * Expec,
                           float * Mask,
                           float * M,
                           float * V,
                           int biasOrder,
                           int nrOfClasses,
                           int aceletation_factor,
                           int xyzsize[3]);

NIFTYSEG_WINEXPORT void get_xyz_pow_int(segPrecisionTYPE xpos,
                     segPrecisionTYPE ypos,
                     segPrecisionTYPE zpos,
                     segPrecisionTYPE currxpower[10],
                     segPrecisionTYPE currypower[10],
                     segPrecisionTYPE currzpower[10],
                     int maxorder);

NIFTYSEG_WINEXPORT void BiasCorrection(segPrecisionTYPE * BiasField,
                    segPrecisionTYPE * BiasFieldCoefs,
                    nifti_image * T1,
                    segPrecisionTYPE * Expec,
                    segPrecisionTYPE * Outlierness,
                    segPrecisionTYPE * M,
                    segPrecisionTYPE * V,
                    int biasOrder,
                    ImageSize * CurrSizes,
                    bool flag_Bias,
                    int verbose_level);

NIFTYSEG_WINEXPORT void BiasCorrection_mask(segPrecisionTYPE * BiasField,
                         segPrecisionTYPE * BiasFieldCoefs,
                         nifti_image * T1,
                         int * Long_2_Short_Indices,
                         segPrecisionTYPE * Expec,
                         segPrecisionTYPE * Outlierness,
                         segPrecisionTYPE * M,
                         segPrecisionTYPE * V,
                         int biasOrder,
                         ImageSize * CurrSizes,
                         bool flag_Bias,
                         int verbose_level);

NIFTYSEG_WINEXPORT void get_xy_pow_int(segPrecisionTYPE xpos,
                    segPrecisionTYPE ypos,
                    segPrecisionTYPE currxpower[10],
                    segPrecisionTYPE currypower[10],
                    int maxorder);

NIFTYSEG_WINEXPORT void BiasCorrection2D(segPrecisionTYPE * BiasField,
                      segPrecisionTYPE * BiasFieldCoefs,
                      nifti_image * T1,
                      segPrecisionTYPE * Expec,
                      segPrecisionTYPE * Outlierness,
                      segPrecisionTYPE * M,
                      segPrecisionTYPE * V,
                      int biasOrder,
                      ImageSize * CurrSizes,
                      bool flag_Bias,
                      int verbose_level);

NIFTYSEG_WINEXPORT void BiasCorrection_mask2D(segPrecisionTYPE * BiasField,
                           segPrecisionTYPE * BiasFieldCoefs,
                           nifti_image * T1,
                           int * Long_2_Short_Indices,
                           segPrecisionTYPE * Expec,
                           segPrecisionTYPE * Outlierness,
                           segPrecisionTYPE * M,
                           segPrecisionTYPE * V,
                           int biasOrder,
                           ImageSize * CurrSizes,
                           bool flag_Bias,
                           int verbose_level);

