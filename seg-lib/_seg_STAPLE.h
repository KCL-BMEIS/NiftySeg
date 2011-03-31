#pragma once

#include "_seg_common.h"
#include "_seg_tools.h"


template <class T>
        class seg_STAPLE
{
protected:

    nifti_image*    inputLABLES; // pointer to external
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
    float  Q[MaxSTAPLElable];
    float  P[MaxSTAPLElable];
    float*  W;
    float   Prop;
    bool    Fixed_Prop_status;
    bool    PropUpdate;
    int     numb_lables;
    float   loglik;
    float   oldloglik;
    int     maxIteration;
    float   Conv;

    // LNCC
    float * LNCC;
    bool LNCC_status;
    float lnccthresh;

    // MRF
    bool    MRF_status;
    float   MRF_strength;
    float*  MRF;

    // Private funcs
    int Maximization();
    int Create_CurrSizes();
    int EstimateInitialDensity();
    int UpdateDensity();
    int Expectation();
    int UpdateMRF();
    int Allocate_Expec_and_MRF();

public:
    seg_STAPLE(int numb_lables);
    ~seg_STAPLE();
    int SetPrior(nifti_image * _Prior,float lnccthresh);
    int SetInputLables(nifti_image * LABLES);
    int SetLNCC(nifti_image * LNCC,nifti_image * BaseImage,int distance,float lnccthresh,bool saveLNCC);
    int SetProp(float prior);
    int SetConv(float conv);
    int SetFilenameOut(char *);
    int SetPQ(float tmpP,float tmpQ);
    int Turn_MRF_ON(float MRF_strenght);
    int Turn_Prop_Update_ON();
    int SetMaximalIterationNumber(unsigned int numberiter);
    int SetVerbose(unsigned int verblevel);

    int CheckParameters_EM();
    int Initisalise_EM();
    int* Run_EM();


    virtual int CheckMemoryMB_EM(){return 0;};
    nifti_image *GetResult();
};

#include "_seg_STAPLE.cpp"

