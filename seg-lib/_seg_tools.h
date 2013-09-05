#pragma once

#include "_seg_common.h"
#include "_seg_matrix.h"
#include "_seg_FMM.h"
#include "_seg_Topo.h"
#include "_seg_BiasCorrection.h"

#if (defined(_WIN32) || defined(_WINDOWS)) && !defined(__CYGWIN__)
#define SEP "\\0"
#else
#define SEP "/\0"
#endif

int printloglik(int iter,
                SegPrecisionTYPE loglik,
                SegPrecisionTYPE oldloglik);
int printTrace(int iter,
                SegPrecisionTYPE loglik,
                SegPrecisionTYPE oldloglik);

int calcM_mask(nifti_image * T1,
               SegPrecisionTYPE * Expec,
               SegPrecisionTYPE * BiasField,
               SegPrecisionTYPE * Outlierness,
               int * S2L,
               SegPrecisionTYPE * M,
               SegPrecisionTYPE * V,
               SegPrecisionTYPE * M_MAP,
               SegPrecisionTYPE * V_MAP,
               SegPrecisionTYPE reg_factor,
               ImageSize * CurrSizes,
               int verbose);

int calcM_mask_LoAd(nifti_image * T1,
                    SegPrecisionTYPE * Expec,
                    SegPrecisionTYPE * BiasField,
                    int * S2L,
                    SegPrecisionTYPE * M,
                    SegPrecisionTYPE * V,
                    ImageSize * CurrSizes,
                    int verbose,
                    bool PVon);

int calcM(nifti_image * T1,
          SegPrecisionTYPE * Expec,
          SegPrecisionTYPE * BiasField,
          SegPrecisionTYPE * Outlierness,
          SegPrecisionTYPE * M,
          SegPrecisionTYPE * V,
          SegPrecisionTYPE * M_MAP,
          SegPrecisionTYPE * V_MAP,
          SegPrecisionTYPE reg_factor,
          ImageSize * CurrSizes,
          int verbose);

int calcE(nifti_image * T1,
          SegPrecisionTYPE * MRF,
          SegPrecisionTYPE * Expec,
          double * loglik,
          SegPrecisionTYPE * BiasField,
          SegPrecisionTYPE * Outlierness,
          SegPrecisionTYPE OutliernessThreshold,
          SegPrecisionTYPE * M,
          SegPrecisionTYPE * V,
          ImageSize * CurrSizes,
          int verbose);

int calcE_mask(nifti_image * T1,
               SegPrecisionTYPE * IterPrior,
               SegPrecisionTYPE * Expec,
               double * loglik,
               SegPrecisionTYPE * BiasField,
               SegPrecisionTYPE * Outlierness,
               SegPrecisionTYPE OutliernessThreshold,
               int * S2L,
               SegPrecisionTYPE * M,
               SegPrecisionTYPE * V,
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
int Relax_Priors(SegPrecisionTYPE * Priors,
                 SegPrecisionTYPE * Expec,
                 SegPrecisionTYPE * MRF,
                 int * S2L,
                 int * L2S,
                 float RelaxFactor,
                 SegPrecisionTYPE * G,
                 SegPrecisionTYPE ba,
                 SegPrecisionTYPE be,
                 ImageSize * CurrSizes,
                 SEG_PARAM * segment_param);

int Normalize_Image(nifti_image * T1,
                    ImageSize * CurrSizes,
                    bool verbose);

int Normalize_Image_mask(nifti_image * T1,
                         nifti_image * Mask,
                         ImageSize * CurrSizes,
                         bool verbose);

int RelaxPriors( SegPrecisionTYPE * Expec,
                 SegPrecisionTYPE * MRF,
                 int * S2L ,int * L2S ,
                 ImageSize CurrSizes,
                 int class_with_CSF);

int Convert_to_PV(nifti_image * T1,
                  SegPrecisionTYPE * BiasField,
                  SegPrecisionTYPE * ShortPrior,
                  SegPrecisionTYPE * Expec,
                  SegPrecisionTYPE * MRF,
                  SegPrecisionTYPE * M,
                  SegPrecisionTYPE * V,
                  int * S2L,
                  int * L2S,
                  ImageSize * CurrSizes,
                  SEG_PARAM * segment_param);

bool * binarise_image(SegPrecisionTYPE * SingleImage,
                      SegPrecisionTYPE Threshold,
                      ImageSize * CurrSizes);

int Create_GH_5class(SegPrecisionTYPE * G,
                     SegPrecisionTYPE * H,
                     SegPrecisionTYPE ba,
                     SegPrecisionTYPE be,
                     SegPrecisionTYPE ratio,
                     SEG_PARAM * segment_param);

int Create_GH_7class(SegPrecisionTYPE * G,
                     SegPrecisionTYPE * H,
                     SegPrecisionTYPE ba,
                     SegPrecisionTYPE be,
                     SegPrecisionTYPE ratio,
                     SEG_PARAM * segment_param);

int Relax_Priors_Share(SegPrecisionTYPE * Priors,
                       SegPrecisionTYPE * Expec,
                       float RelaxFactor,
                       SegPrecisionTYPE * G,
                       SegPrecisionTYPE be,
                       ImageSize * CurrSizes);

int seg_convert2binary(nifti_image *image,
                       float thresh);

int Normalize_T1_and_MV(nifti_image * T1,
                        nifti_image * Mask,
                        SegPrecisionTYPE * M,
                        SegPrecisionTYPE * V,
                        ImageSize * CurrSizes);

int Normalize_NaN_Priors(nifti_image * Priors,
                         bool verbose);

int Normalize_NaN_Priors_mask(nifti_image * Priors,
                              nifti_image * Mask,
                              bool verbose);

int Convert_WM_and_GM_to_PV(nifti_image * T1,
                            SegPrecisionTYPE * BiasField,
                            SegPrecisionTYPE * ShortPrior,
                            SegPrecisionTYPE * Expec,
                            int * S2L,
                            SegPrecisionTYPE * M,
                            SegPrecisionTYPE * V,
                            ImageSize * CurrSize);

int Gaussian_Filter_Short_4D(SegPrecisionTYPE * ShortData,
                             int * S2L,
                             int * L2S,
                             SegPrecisionTYPE gauss_std,
                             ImageSize * CurrSizes,
                             int class_with_CSF);

int Gaussian_Filter_4D(SegPrecisionTYPE * LongData,
                       SegPrecisionTYPE gauss_std,
                       ImageSize * CurrSizes);


SegPrecisionTYPE * Gaussian_Filter_4D_inside_mask(SegPrecisionTYPE * LongData,
                                                  bool * mask,
                                                  SegPrecisionTYPE gauss_std,
                                                  ImageSize * CurrSizes);

void GaussianSmoothing(nifti_image * Data,int * mask,float gauss_std);
void GaussianSmoothing_carray(float * DataPTR,int * mask,float gauss_std_in, ImageSize *Currentsize);

int Create_diagonal_GH_Nclass(SegPrecisionTYPE * G,
                              SegPrecisionTYPE * H,
                              SegPrecisionTYPE ratio,
                              SEG_PARAM * segment_param);

int Sulci_and_gyri_correction(SegPrecisionTYPE * MRF_Beta,
                              SegPrecisionTYPE * ShortPrior,
                              SegPrecisionTYPE * Expec,
                              SegPrecisionTYPE *MRF,
                              int * S2L,
                              int * L2S,
                              ImageSize *CurrSizes);

nifti_image * Copy_ShortExpec_to_Result(nifti_image * T1,
                                        SegPrecisionTYPE * Expec,
                                        SegPrecisionTYPE * BiasField,
                                        SegPrecisionTYPE * BiasFieldCoefs,
                                        int * S2L,
                                        nifti_image * Priors,
                                        SEG_PARAM * segment_param,
                                        SegPrecisionTYPE * M,
                                        ImageSize * CurrSizes);

int * Create_Long_2_Short_Matrix_from_NII(nifti_image * Mask);

int * Create_Short_2_Long_Matrix_from_NII(nifti_image * Mask,
                                          long *shortsize);

SegPrecisionTYPE * Create_cArray_from_Prior_mask(nifti_image * Mask,
                                                  nifti_image * Priors,
                                                  long numclass,
                                                  bool PV_ON);

SegPrecisionTYPE * Create_cArray_from_Prior(nifti_image * Priors,
                                            long numclass,
                                            bool PV_ON);

SegPrecisionTYPE * Create_cArray_from_3D_image(nifti_image * Mask,
                                               nifti_image * SourceImage);

int * Create_Short_2_Long_Matrix_from_Carray(bool * Mask,
                                             int * shortsize,
                                             int nvox);

int *  Create_Long_2_Short_Matrix_from_Carray(bool * Mask,
                                              int nvox);

nifti_image * Copy_Single_ShortImage_to_Result(SegPrecisionTYPE * SingleImage,
                                               int * Short_2_Long_Indices,
                                               nifti_image * Priors,
                                               char * filename,
                                               ImageSize * CurrSizes);

template <class DTYPE> int seg_convert2binary_data(nifti_image *image,
                                                   float thresh);

nifti_image * Copy_Expec_and_BiasCorrected_to_Result_mask(SegPrecisionTYPE * Expec,
                                                          SegPrecisionTYPE * BiasField,
                                                          int * Short_2_Long_Indices,
                                                          nifti_image * T1,
                                                          char * filename,
                                                          ImageSize * CurrSizes);
nifti_image * Copy_Expec_and_BiasCorrected_to_Result(SegPrecisionTYPE * Expec,
                                                     SegPrecisionTYPE * BiasField,
                                                     nifti_image * T1,
                                                     char * filename,
                                                     ImageSize * CurrSizes);

nifti_image * Copy_Expec_to_Result_mask(SegPrecisionTYPE * Expec,
                                        int * Short_2_Long_Indices,
                                        nifti_image * T1,
                                        char * filename,
                                        ImageSize * CurrSizes);

nifti_image * Copy_Expec_to_Result_Neonate_mask(SegPrecisionTYPE * Expec,
                                                int * Short_2_Long_Indices,
                                                int * Long_2_Short_Indices,
                                                nifti_image * T1,
                                                float * Biasfield,
                                                float * M,
                                                char * filename,
                                                ImageSize * CurrSizes);


nifti_image * Copy_Expec_to_Result(SegPrecisionTYPE * Expec,
                                   nifti_image * T1,
                                   char * filename,
                                   ImageSize * CurrSizes);

int PriorWeight_mask(float * ShortPrior,
                     nifti_image * Priors,
                     float * Expec,
                     float GaussKernelSize,
                     float RelaxFactor,
                     int * S2L,
                     int * L2S,
                     ImageSize * CurrSizes,
                     int verbose_level);

nifti_image * Copy_single_image_to_Result(int * Mask,
                                          nifti_image * Original,
                                          char * filename);

nifti_image * Copy_single_image_to_Result(bool * Mask,
                                          nifti_image * Original,
                                          char * filename);
nifti_image * Copy_single_image_to_Result(float * Mask,
                                          nifti_image * Original,
                                          char * filename);
nifti_image * Copy_single_image_to_Result(double * Mask,
                                          nifti_image * Original,
                                          char * filename);

int quickSort(int *arr, int elements);
int quickSort(float *arr, int elements);

int * quickSort_order(int *arr, int elements);
int * quickSort_order(float *arr, int elements);

void HeapSort(float * a,int n);

nifti_image * Get_Bias_Corrected(float * BiasField,
                                 nifti_image * T1,
                                 char * filename,
                                 ImageSize * CurrSizes);

nifti_image * Get_Bias_Corrected_mask(float * BiasFieldCoefs,
                                      nifti_image * T1,
                                      char * filename,
                                      ImageSize * CurrSizes,
                                      int biasOrder);

unsigned char * seg_norm_4D_GNCC(nifti_image * BaseImage,nifti_image * NCC,int numberordered,ImageSize * CurrSizes,int verbose);
float seg_norm3GNCC(nifti_image * BaseImage,nifti_image * Template,nifti_image * Mask,int verbose);
unsigned char * seg_norm4ROINCC(nifti_image * LableImage,nifti_image * BaseImage,nifti_image * LNCC,int numberordered,ImageSize * CurrSizes,int DilSize, int verbose);
unsigned char * seg_norm4LNCC(nifti_image * BaseImage,nifti_image * LNCC,float distance,int numberordered,ImageSize * CurrSizes,int verbose);
unsigned char * seg_norm4MLLNCC(nifti_image * BaseImage, nifti_image * LNCC,float distance,int labels, int numberordered,ImageSize * CurrSizes,int verbose);

extern "C++" template <class NewTYPE>
int seg_changeDatatype(nifti_image *image);







/* *************************************************************** */
template<class SourceTYPE, class FieldTYPE>
void Resample_NN_with_weights(  nifti_image *sourceImage,
                                nifti_image *deformationField,
                                nifti_image *resultImage,
                                nifti_image *resultImageWeights,
                                int *mask,
                                float bgValue);
extern "C++" template <class DTYPE>
void seg_mat44_mul(mat44 *mat,
                   DTYPE *in,
                   DTYPE *out);
int get_all_files_and_folders_in_dir (string dir, vector<string> &files , vector<string> &folders);
int get_all_files_that_match_string (string dir, vector<string> &files , string string_to_match);
int get_all_files_that_match_2_strings(string dir, vector<string> &files , string string_to_match, string string_to_match2);
int get_all_files_in_dir_without_extension(string dir, vector<string> &files);
float * getHeatWij(float * DistanceMatrix,int size_matrix, float temperature);
float * getNN(float * DistanceMatrix,int size_matrix,int sizeneig);

void LTS_Vecs(float * Y, float * X,int * mask, float percentOutliers,int maxNumbIter, float convergenceRatio, unsigned int size, float *a, float *b);
void LS_Vecs(float * Y, float * X,int * mask, unsigned int size, float *a, float *b);

void otsu(float * Image, int * mask, ImageSize *Currentsize );
void BiasCorrect(float * Image,ImageSize *Currentsize);
