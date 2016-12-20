#ifndef _SEG_FILL_LESIONS_H
#define _SEG_FILL_LESIONS_H
#include "_seg_tools.h"

template <class T>
class seg_fill_lesions
{
protected:
    nifti_image *inputImage;
    nifti_image *normImage;
    nifti_image *meanImage;
    nifti_image *outputImage;
    nifti_image *inputLesionMask;
    nifti_image *inputMask;
    int verbose;
    int debug;
    float mult;
    int patchSearchAreaSize;
    float percentage;
    float expanding;
    int numTP;
    T *imgPtr;
    float *normImgPtr;
    float *meanImgPtr;
    T *lesMaskPtr;
    T *maskPtr;
    int *tmpLesMask;
    int *originalLesMask;
    float *tmpImg;
    float *tmpNorm;
    ImageSize * currSize;
    float k;
    bool patch2D;

    int   getPatchSearchAreaSize();
    float calculateDistance(int,int,int,int);
    int   getNumTP();
    float getK();
    void  setNumTP(int);
    long  isMaskInAnyTP(long,int *);
    bool  isValidVoxelInAnyTP(long);
    void  calculateEuclideanDistance(float *);
    void  expandPatches(float *);
    void  normalizeImageIntesities(float,float);
    long  getSingleVolumSize();
    long  getTotalVolumSize();
    long  getSliceSize();
    long  getXAxisSize();
    long  getYAxisSize();
    long  getZAxisSize();
    float getSmoothing();
    float getExpandingPercentage();
    void  saveImagePtr(int *,nifti_image *,char *);
    void  saveImagePtr(float *,nifti_image *,char *);
    int   countLesionVoxels(int,int);
    bool  isIt2D();
    int   getDebug();
    int   getVerbose();

public:
    seg_fill_lesions();
    seg_fill_lesions(const seg_fill_lesions &o);
    int operator=(const seg_fill_lesions &o);
    ~seg_fill_lesions();

    void setInputImage(nifti_image *);
    void setInputMask(nifti_image *);
    void setInputLesionMask(nifti_image *);
    void setPatchSearchAreaSize(int);
    void setPatchPercentage(float);
    void setExpandingPercentage(float);
    void setSmoothing(float);
    void setVerbose(int);
    void setDebug(int);
    void setK(float);
    void setDimensionality(bool);
    void init();
    void runIt();
    void saveImage(nifti_image *,char *);
};

#include "_seg_fill_lesions.cpp"

#endif // _SEG_FILL_LESIONS_H
