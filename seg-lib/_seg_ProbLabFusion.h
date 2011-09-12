#pragma once

#include "_seg_common.h"
#include "_seg_tools.h"


class seg_LabFusion
{
protected:
    nifti_image*    inputLABELS; // pointer to external
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

    float Thresh_IMG_value;
    bool Thresh_IMG_DO;
    // SegParameters
    LabFusion_datatype * Q;
    LabFusion_datatype * P;
    LabFusion_datatype*  W;
    LabFusion_datatype   Prop;
    bool*  uncertainarea;
    bool    Fixed_Prop_status;
    bool    PropUpdate;
    int     numb_lables;
    float   loglik;
    float   oldloglik;
    int     maxIteration;
    float   Conv;

    // LNCC
    char * LNCC;
    bool LNCC_status;
    char * NCC;
    bool NCC_status;
    int Numb_Neigh;

    // MRF
    bool    MRF_status;
    LabFusion_datatype   MRF_strength;
    LabFusion_datatype*  MRF;

    // Private funcs
    int STAPLE_Maximization();
    int Create_CurrSizes();
    int EstimateInitialDensity();
    int UpdateDensity();
    int STAPLE_Expectation();
    int MV_Estimate();
    int UpdateMRF();
    int Allocate_Stuff();

public:
    seg_LabFusion(int numb_lables);
    ~seg_LabFusion();
    int SetInputLabels(nifti_image * LABELS, bool UNCERTAINflag);
    int SetLNCC(nifti_image * LNCC,nifti_image * BaseImage,float distance,int Numb_Neigh);
    int SetGNCC(nifti_image * _GNCC,nifti_image * BaseImage,int Numb_Neigh);
    int SetROINCC(nifti_image * _ROINCC,nifti_image * BaseImage,int Numb_Neigh, int DilSize);
    int SetProp(float prior);
    int SetConv(float conv);
    int SetImgThresh(float Thresh_IMG_value);
    int SetFilenameOut(char *);
    int SetPQ(float tmpP,float tmpQ);
    int Turn_MRF_ON(float MRF_strenght);
    int Turn_Prop_Update_ON();
    int SetMaximalIterationNumber(unsigned int numberiter);
    int SetVerbose(unsigned int verblevel);

    int CheckParameters_EM();
    int Initisalise_EM();
    int Run_STAPLE();
    int Run_MV();
    int Run_SBA();
    nifti_image *GetResult();
};



