#pragma once

#include "_seg_common.h"
#include "_seg_FMM.h"

#if (defined(_WIN32) || defined(_WINDOWS)) && !defined(__CYGWIN__)
#define SEP "\\0"
#else
#define SEP "/\0"
#endif

// Gaussian convolution functions
void GaussianFilter4D_cArray(segPrecisionTYPE * ShortData, int * S2L, int * L2S, segPrecisionTYPE gauss_std, ImageSize * CurrSizes);
void GaussianFilter4D_cArray(segPrecisionTYPE * LongData, segPrecisionTYPE gauss_std, ImageSize * CurrSizes);
void GaussianSmoothing4D_nifti(nifti_image * Data,int * mask,float gauss_std);
void BlockSmoothing(nifti_image * Data,int * mask,int side_size);
void SmoothLab(float * DataPTR, float factor, ImageSize * CurrSizes);

// Data conversion tools
int seg_convert2binary(nifti_image *image, float thresh);
template <class DTYPE> int seg_convert2binary_data(nifti_image *image, float thresh);
template <class NewTYPE>int seg_changeDatatype(nifti_image *image);

// Sorting algrithms (with and without providing the order), for different data types
int quickSort(int *arr, int elements);
int quickSort(float *arr, int elements);
int * quickSort_order(int *arr, int elements);
int * quickSort_order(float *arr, int elements);
void HeapSort(float * a,int n);

// Estimate the order of multiple types of NCC similarity (local, global, regional, multilevel) given the input target image and a 4D set of resampled images. To be used for label fusion.
unsigned char * estimateNCC4D(nifti_image * BaseImage,nifti_image * NCC,int numberordered,ImageSize * CurrSizes,int verbose);
unsigned char * estimateROINCC4D(nifti_image * LableImage,nifti_image * BaseImage,nifti_image * LNCC,int numberordered,ImageSize * CurrSizes,int DilSize, int verbose);
unsigned char * estimateLNCC4D(nifti_image * BaseImage,nifti_image * LNCC,float distance,int numberordered,ImageSize * CurrSizes,int verbose);
unsigned char * estimateMLNCC4D(nifti_image * BaseImage, nifti_image * LNCC,float distance,int labels, int numberordered,ImageSize * CurrSizes,int verbose);
float estimateNCC3D(nifti_image * BaseImage,nifti_image * Template,nifti_image * Mask,int verbose);

// Data scraping tools, used to get files/folders inside directories (with or without string match)
int get_all_files_and_folders_in_dir(string dir, vector<string> &files , vector<string> &folders);
int get_all_files_that_match_string(string dir, vector<string> &files , string string_to_match);
int get_all_files_that_match_2_strings(string dir, vector<string> &files , string string_to_match, string string_to_match2);
int get_all_files_in_dir_without_extension(string dir, vector<string> &files);

// Estimate the LTS or the LS vetween two images (<X>,<Y>), inside a <mask>
void LTS_Vecs(float * Y, float * X,int * mask, float percentOutliers,int maxNumbIter, float convergenceRatio, unsigned int size, float *a, float *b);
void LS_Vecs(float * Y, float * X,int * mask, unsigned int size, float *a, float *b);

// Estimate some mathematical morphology operators, such as the connected components, dilations, erosions, etc.
void ConnectComp(int * Old, int * New, int dimensions[3],int varin);
template <class OldType, class NewType> void Largest_ConnectComp(void * Old, void * New, ImageSize * Currentsize);
template <class OldType, class NewType> void Close_Forground_ConnectComp(void * Old, void * New, ImageSize * Currentsize);
void Dillate(float * Image,int kernel,ImageSize * Currentsize );
void Erosion(float * Image,int kernel,ImageSize * Currentsize );
void Dillate_const(bool * Image, bool * Const, int kernel, int dimensions[3], int direction);

// Otsu intensity thresholding.
void otsu(float * Image, int * mask, ImageSize *Currentsize );
