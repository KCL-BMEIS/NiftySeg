#pragma once

#include "_seg_common.h"
#include "_seg_tools.h"


class seg_LabFusion
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
  int     numel;
  int     iter;
  ImageSize * CurrSizes;

  int NUMBER_OF_CLASSES;
  float Thresh_IMG_value;
  bool Thresh_IMG_DO;
  // SegParameters
  int LableCorrespondences_big_to_small[5000];
  int LableCorrespondences_small_to_big[5000];
  LabFusion_datatype * ConfusionMatrix;
  LabFusion_datatype * W;
  LabFusion_datatype *  Prop;
  bool*  uncertainarea;
  bool    uncertainflag;
  int dilunc;
  bool    Fixed_Prop_status;
  bool    PropUpdate;
  int     numb_classif;
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
  LabFusion_datatype* MRF_matrix;

  // Private funcs
  int Create_CurrSizes();
  int EstimateInitialDensity();
  int UpdateDensity();
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
  seg_LabFusion(int _numb_classif,int _numb_labels,int _numb_neigh);
  ~seg_LabFusion();
  int SetinputCLASSIFIER(nifti_image * LABELS, bool UNCERTAINflag);
  int SetMLLNCC(nifti_image * LNCC,nifti_image * BaseImage,float distance,int levels, int Numb_Neigh);
  int SetLNCC(nifti_image * LNCC,nifti_image * BaseImage,float distance,int Numb_Neigh);
  int SetGNCC(nifti_image * _GNCC,nifti_image * BaseImage,int Numb_Neigh);
  int SetROINCC(nifti_image * _ROINCC,nifti_image * BaseImage,int Numb_Neigh, int DilSize);
  int SetProp(float prior);
  int SetConv(float conv);
  int SetDilUnc(int _dilunc);
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
  nifti_image *GetResult(int ProbOutput);
};



