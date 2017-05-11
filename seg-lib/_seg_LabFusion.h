#pragma once
#include "NiftySegWinExportHeader.h"

#include "_seg_common.h"
extern "C++" {
#include "_seg_tools.h"
}


class NIFTYSEG_WINEXPORT seg_LabFusion
{
protected:
    nifti_image*    inputCLASSIFIER; // pointer to external
    bool    inputImage_status;
    string  FilenameOut;
    int     verbose_level;

    // Size
    int     dimentions;
    int     nx;
    int     ny;
    int     nz;
    int     nt;
    int     nu;
    float     dx;
    float     dy;
    float     dz;
    long     numel;
    int     iter;
    ImageSize * CurrSizes;

    int TYPE_OF_FUSION; // 1 - STEPS/STAPLE ; 2 - MV ; 3 - SBA
    long NumberOfLabels;
    float Thresh_IMG_value;
    bool Thresh_IMG_DO;
    // SegParameters
    int LableCorrespondences_big_to_small[5000];
    int LableCorrespondences_small_to_big[5000];
    segPrecisionTYPE * ConfusionMatrix;
    segPrecisionTYPE * W;
    segPrecisionTYPE * FinalSeg;
    segPrecisionTYPE * Prop;
    long * maskAndUncertainIndeces;
    long  sizeAfterMaskingAndUncertainty;
    bool    uncertainflag;
    float uncertainthresh;
    int dilunc;
    bool    Fixed_Prop_status;
    bool    PropUpdate;
    int     numb_classif;
    int     numb_nummod;
    float   tracePQ;
    float   oldTracePQ;
    float   loglik;
    float   oldloglik;
    int     maxIteration;
    float   Conv;

    // LNCC
    unsigned char * LNCC;
    bool LNCC_status;
    unsigned char * NCC;
    bool NCC_status;
    int Numb_Neigh;

    // MRF
    bool    MRF_status;
    segPrecisionTYPE   MRF_strength;
    segPrecisionTYPE*  MRF;
    segPrecisionTYPE* MRF_matrix;

    // Private funcs
    int Create_CurrSizes();
    int EstimateInitialDensity();
    int UpdateDensity();
    int UpdateDensity_noTest();

    //    int Find_WMax();
    int STAPLE_STEPS_Multiclass_Maximization();
    int STAPLE_STEPS_Multiclass_Expectation();
    int STAPLE_STEPS_Multiclass_Expectation_Maximization();
    int MV_Estimate();
    int SBA_Estimate();
    int UpdateMRF();
    int Allocate_Stuff_MV();
    int Allocate_Stuff_SBA();
    int Allocate_Stuff_STAPLE();

public:
    seg_LabFusion(int _numb_classif,int _numb_labels,int _numb_neigh, int _numb_modalities);
    ~seg_LabFusion();
    int SetinputCLASSIFIER(nifti_image * LABELS, bool UNCERTAINflag);
    int SetMLLNCC(nifti_image * LNCC,nifti_image * BaseImage,float distance,int levels, int Numb_Neigh);
    int SetLNCC(nifti_image * LNCC,nifti_image * BaseImage,float distance,int Numb_Neigh);
    int SetLMETRIC(nifti_image * _METRIC,int Numb_Neigh);
    int SetGNCC(nifti_image * _GNCC,nifti_image * BaseImage,int Numb_Neigh);
    int SetROINCC(nifti_image * _ROINCC,nifti_image * BaseImage,int Numb_Neigh, int DilSize);
    int SetProp(float prior);
    int SetConv(float conv);
    int SetDilUnc(int _dilunc);
    int SetUncThresh(float _uncthresh);
    int SetImgThresh(float Thresh_IMG_value);
    int SetFilenameOut(char *);
    int SetPQ(float tmpP,float tmpQ);
    int SetMask(nifti_image * Mask);

    int Turn_MRF_ON(float MRF_strenght);
    int Turn_Prop_Update_ON();
    int SetMaximalIterationNumber(unsigned int numberiter);
    int SetVerbose(unsigned int verblevel);

    int CheckParameters_EM();
    int Initisalise_EM();
    int Run_STAPLE_or_STEPS();
    int Run_MV();
    int Run_SBA();
    nifti_image *GetResult_label();
    nifti_image *GetResult_probability();
};



