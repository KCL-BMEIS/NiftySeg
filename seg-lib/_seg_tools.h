#pragma once
#include "NiftySegWinExportHeader.h"

#include "_seg_common.h"
#include "_seg_FMM.h"

#if (defined(_WIN32) || defined(_WINDOWS)) && !defined(__CYGWIN__)
#define SEP "\\"
#else
#define SEP "/\0"
#endif


// Gaussian convolution functions
NIFTYSEG_WINEXPORT void GaussianFilter4D_cArray(segPrecisionTYPE * ShortData, int * S2L, int * L2S, segPrecisionTYPE gauss_std, ImageSize * CurrSizes);
NIFTYSEG_WINEXPORT void GaussianFilter4D_cArray(segPrecisionTYPE * LongData, segPrecisionTYPE gauss_std, ImageSize * CurrSizes);
NIFTYSEG_WINEXPORT void GaussianSmoothing4D_Nan_nifti(nifti_image * Data, nifti_image * mask);
NIFTYSEG_WINEXPORT void GaussianSmoothing5D_nifti(nifti_image * Data,int * mask,float gauss_std);
NIFTYSEG_WINEXPORT void BlockSmoothing(nifti_image * Data,int * mask,int side_size);
NIFTYSEG_WINEXPORT void SmoothLab(float * DataPTR, float factor, ImageSize * CurrSizes);

// Data conversion tools
NIFTYSEG_WINEXPORT int seg_convert2binary(nifti_image *image, float thresh);

template <class DTYPE> NIFTYSEG_WINEXPORT
	int seg_convert2binary_data(nifti_image *image, float thresh);

template <class NewTYPE> NIFTYSEG_WINEXPORT
	int seg_changeDatatype(nifti_image *image);

//template NIFTYSEG_WINEXPORT int seg_changeDatatype<unsigned char>(nifti_image *);
//template NIFTYSEG_WINEXPORT int seg_changeDatatype<float>(nifti_image *);
//template NIFTYSEG_WINEXPORT int seg_changeDatatype<double>(nifti_image *);

// Sorting algrithms (with and without providing the order), for different data types
NIFTYSEG_WINEXPORT int quickSort(int *arr, int elements);
NIFTYSEG_WINEXPORT int quickSort(float *arr, int elements);
NIFTYSEG_WINEXPORT int * quickSort_order(int *arr, int elements);
NIFTYSEG_WINEXPORT int * quickSort_order(float *arr, int elements);
NIFTYSEG_WINEXPORT void HeapSort(float * a,int n);

// Estimate the order of multiple types of NCC similarity (local, global, regional, multilevel) given the input target image and a 4D set of resampled images. To be used for label fusion.
NIFTYSEG_WINEXPORT unsigned char * estimateNCC4D(nifti_image * BaseImage,nifti_image * NCC,int numberordered,ImageSize * CurrSizes,int verbose);
NIFTYSEG_WINEXPORT unsigned char * estimateROINCC4D(nifti_image * LableImage,nifti_image * BaseImage,nifti_image * LNCC,int numberordered,ImageSize * CurrSizes,int DilSize, int verbose);
NIFTYSEG_WINEXPORT unsigned char * estimateLNCC5D(nifti_image * BaseImage,nifti_image * LNCC,float distance,int numberordered,ImageSize * CurrSizes,int verbose);
NIFTYSEG_WINEXPORT unsigned char * estimateLNCC4D(nifti_image * BaseImage,nifti_image * LNCC,float distance,int numberordered,ImageSize * CurrSizes,int verbose);
NIFTYSEG_WINEXPORT unsigned char * estimateMLNCC4D(nifti_image * BaseImage, nifti_image * LNCC,float distance,int labels, int numberordered,ImageSize * CurrSizes,int verbose);
NIFTYSEG_WINEXPORT float estimateNCC3D(nifti_image * BaseImage,nifti_image * Template,nifti_image * Mask,int verbose);
NIFTYSEG_WINEXPORT float seg_getNMIValue(nifti_image *referenceImage, nifti_image *warpedImage, unsigned char *referenceMask);

// Data scraping tools, used to get files/folders inside directories (with or without string match)
NIFTYSEG_WINEXPORT int get_all_files_and_folders_in_dir(string dir, vector<string> &files , vector<string> &folders);
NIFTYSEG_WINEXPORT int get_all_files_that_match_string(string dir, vector<string> &files , string string_to_match);
NIFTYSEG_WINEXPORT int get_all_files_that_match_2_strings(string dir, vector<string> &files , string string_to_match, string string_to_match2);
NIFTYSEG_WINEXPORT int get_all_files_in_dir_without_extension(string dir, vector<string> &files);

// Estimate the LTS or the LS vetween two images (<X>,<Y>), inside a <mask>
NIFTYSEG_WINEXPORT void LTS_Vecs(float * Y, float * X,int * mask, float percentOutliers,int maxNumbIter, float convergenceRatio, unsigned int size, float *a, float *b);
NIFTYSEG_WINEXPORT void LS_Vecs(float * Y, float * X,int * mask, unsigned int size, float *a, float *b);

// Estimate some mathematical morphology operators, such as the connected components, dilations, erosions, etc.
NIFTYSEG_WINEXPORT void ConnectComp(int * Old, int * New, int dimensions[3],int varin);

template <class OldType, class NewType> NIFTYSEG_WINEXPORT
  void Largest_ConnectComp(void * Old, void * New, ImageSize * Currentsize);

template <class OldType, class NewType> NIFTYSEG_WINEXPORT
	void ConnectComp26NN(void * Old, void * New, ImageSize * Currentsize);

template <class OldType, class NewType> NIFTYSEG_WINEXPORT
	void ConnectComp6NN(void * Old, void * New, ImageSize * Currentsize);

template <class OldType, class NewType> NIFTYSEG_WINEXPORT 
	void Close_Forground_ConnectComp(void * Old, void * New, ImageSize * Currentsize);

NIFTYSEG_WINEXPORT void Dillate(float * Image,int kernel,ImageSize * Currentsize );
NIFTYSEG_WINEXPORT void Erosion(float * Image,int kernel,ImageSize * Currentsize );
NIFTYSEG_WINEXPORT void TopologicalErosion(float * Image, int kernel, ImageSize *Currentsize );
NIFTYSEG_WINEXPORT void Dillate_const(bool * Image, bool * Const, int kernel, int dimensions[3], int direction);

// Otsu intensity thresholding.
NIFTYSEG_WINEXPORT void otsu(float * Image, int * mask, ImageSize *Currentsize );

template <class DTYPE> NIFTYSEG_WINEXPORT void seg_mat44_mul(mat44 const* mat, DTYPE const* in,DTYPE *out);

