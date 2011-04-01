#pragma once

#include "_seg_common.h"
#include "_seg_BiasCorrection.h"
#include "_seg_tools.h"
#include "_seg_MRF.h"
#include "_seg_FMM.h"
#include "_seg_Topo.h"



class seg_EM
{
protected:

    nifti_image*    inputImage; // pointer to external
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
    int     numel;
    int     iter;
    ImageSize * CurrSizes;


    // SegParameters
    float*  M;
    float*  V;
    float*  Expec;
    float*  ShortPrior;
    int*    Short_2_Long_Indices;
    int*    Long_2_Short_Indices;
    int     numb_classes;
    float   loglik;
    float   oldloglik;
    int     maxIteration;
    bool    aprox;

    // Mask
    nifti_image*    Mask; // pointer to external
    bool    maskImage_status;
    int     numelmasked;

    // Priors Specific
    bool    Priors_status;
    nifti_image*  Priors;

    // MRF Specific
    bool    MRF_status;
    float   MRF_strength;
    float*  MRF;
    float*  MRF_beta;
    float*  MRF_transitionMatrix; // G matrix

    // BiasField Specific
    bool    BiasField_status;
    int     BiasField_order;
    float*  BiasField;
    float*  BiasField_coeficients;
    float   BiasField_ratio;
    int     numelbias;

    // LoAd Specific
    int LoAd_phase;
    bool PV_model_status;
    bool SG_deli_status;
    bool Relax_status;
    float  Relax_factor;
    float RelaxGaussKernelSize;
    bool MAP_status;
    float* MAP_M;
    float* MAP_V;

    // Private funcs
    int Create_diagonal_MRF_transitionMatrix();
    int Normalize_T1();
    int Create_CurrSizes();
    int Maximization();
    int Expectation();
    int UpdateMRF();
    int UpdatePriorWeight();
    int UpdateBiasField();
    int Normalize_Image_and_Priors();
    int Allocate_Expec_and_ShortPrior();

public:
    seg_EM(int numb_classes,int NumbMultiSpec,int NumbTimePoints);
    ~seg_EM();
    int SetInputImage(nifti_image *);
    int SetMaskImage(nifti_image *);
    int SetPriorImage(nifti_image *);
    int SetFilenameOut(char *);
    int SetMAP(float *M, float* V);
    int Turn_Relaxation_ON(float relax_factor,float relax_gauss_kernel);
    int Turn_MRF_ON(float MRF_strenght);
    int Turn_BiasField_ON(int BiasField_order,float ratiothresh);
    int SetLoAd(float RelaxFactor,bool Relax_status,bool PV_model_status,bool SG_deli_status);
    int SetMaximalIterationNumber(unsigned int numberiter);
    int SetAprox(bool aproxval);
    int SetVerbose(unsigned int verblevel);


    int CheckParameters_EM();
    int Initisalise_EM();
    int* Run_EM();


    virtual int CheckMemoryMB_EM(){return 0;};
    nifti_image *GetResult();
    nifti_image *GetResultNeonate();
    nifti_image *GetBiasCorrected(char * filename);
};

#include "_seg_EM.cpp"

