/**
 * @file _seg_PatchMatch.h
 * @author Ferran Prados
 * @date 16/06/2016
 *
 * Copyright (c) 2016, University College London. All rights reserved.
 * Centre for Medical Image Computing (CMIC)
 * See the LICENSE.txt file in the nifty_seg root folder
 *
 */
 
#ifndef _SEG_PATCHMATCH_H
#define _SEG_PATCHMATCH_H
#include "_seg_tools.h"
#include "_seg_PatchMatchResult.h"

#define SegPrecisionTYPE float

template <class T>
class seg_PatchMatch
{
protected:
    nifti_image *inputImage;
    nifti_image *inputMask;
    nifti_image *ouputImageSample;
    vector<string> imageDatabaseFileList;
    vector<string> outputDatabaseFileList;
    vector<string> maskDatabaseFileList;
    float *result;
    std::vector<PatchMatchResult *> KNN;
    float *outputMask;
    int *matchingCount;
    int verbose;
    int debug;
    int better_match;
    int pmatch_executions;
    int pm_iter;
    int csa_size;
    int patchSize;
    int distance;
    int numTP;
    T *imgPtr;
    T *maskPtr;
    T *outputMaskPtr;
    ImageSize *currSize;
    ImageSize *outputSize;
    int count;
    bool filling;

    float *db_list_image;
    float *db_list_mask;
    float *db_list_output;

    void  normalizeImageIntesities(float,float,T*);
    void  normalizeImageIntesitiesOutliers(float,float,T*);
    int   getNumTP();
    void  setNumTP(int);
    long  getSingleVolumSize();
    long  getTotalVolumSize();
    long  getSliceSize();
    long  getXAxisSize();
    long  getYAxisSize();
    long  getZAxisSize();
    void  saveImagePtr(int *,nifti_image *,char *);
    void  saveImagePtr(float *,nifti_image *,char *);
    int   getDebug();
    int   getVerbose();
    float calculateSSDDistance(long,long,int,int,float);
    float calculateSSDDistanceSliced(long,long,int,int,int*,int*,int*);
    float calculateSSDistance1Sliced(long,long,int,long,long,int,float,float,int[6][2]);
    void  calculate(int, int,int,int,float &,int &,float &,float &,float &,float &);
    float calculateLNCCDistance1(long,long,int,float);
    float calculateSSDistance1(long,long,int,float);
    float calculateDistance(long,long,int,float);
    void  getNextRecordDatabase(int num,float *&,float *&);
    void  loadFile(nifti_image *&,string,float *&,ImageSize *&);
    void  fusePatch(long,long,int,long);
    void  saveResults();
    void  computePatchMatch(int);
    void  getRandomIndex(long,long,long &,int &s);
    void  propagationSTEP(int,long,long);
    void  propagationSTEP(int,long,long,int[6][2]);
    void  sortResults(long);
    void  constrainedRandomSearch(int,long,long);
    void  loadInputDatabase();
    void  loadOuputDatabase();
    void  patchMatch();
    void  sortKNNResults();
    void  loadingInputData();
    void  initKNNVectors();
    void  saveDebugResults();
    void  labelFusion();
    long  recomputeInputData(int);

public:
    seg_PatchMatch();
   ~seg_PatchMatch();

    void setInputImage(nifti_image *);
    void setInputMask(nifti_image *);
    void setInputImageDatabase(vector<string>);
    void setInputMaskDatabase(vector<string>);
    void setOutputFilesDatabase(vector<string>);
    void setPatchSize(int);
    int  getPatchSize();
    void setBetterMatch(int);
    int  getBetterMatch();
    void setPatchMatchIterations(int);
    int  getPatchMatchIterations();
    void setPatchMatchExecutions(int);
    int  getPatchMatchExecutions();
    void setConstrainedSearchAreaSize(int);
    int  getConstrainedSearchAreaSize();
    void setDistance(int);
    int  getDistance();
    void  setFilling(bool);
    bool  isFilling();
    float* getOutputResult();
    void setVerbose(int);
    void setDebug(int);
    void init();
    void runIt();
    void saveImage(nifti_image *,char *);

};

#include "_seg_PatchMatch.cpp"


#endif // _SEG_PATCHMATCH_H

