#ifndef _SEG_FILL_LESIONS_OTHER_H
#define _SEG_FILL_LESIONS_OTHER_H
#include "_seg_tools.h"

template <class T>
class seg_fill_lesions_other
{
protected:
    nifti_image *inputImage;
    nifti_image *normImage;
    nifti_image *outputImage;
    nifti_image *inputLesionMask;
    int verbose;
    int debug;
    float mult;
    int patchSearchAreaSize;
    int patchSize;
    int numTP;
    T *imgPtr;
    float *normImgPtr;
    T *lesMaskPtr;
    int *tmpLesMask;
    int *originalLesMask;
    float *tmpImg;
    float *tmpNorm;
    ImageSize * currSize;

    int   getSearchArea();
    int   getPatchSize();
    float calculateDistance(int,int);
    int   getNumTP();
    void  setNumTP(int);
    long  isMaskInAnyTP(long,int *);
    bool  isValidVoxelInAnyTP(long);
    void  normalizeImageIntesities(float,float);
    long  getSingleVolumSize();
    long  getTotalVolumSize();
    long  getSliceSize();
    long  getXAxisSize();
    long  getYAxisSize();
    long  getZAxisSize();
    void  saveImagePtr(int *,nifti_image *,char *);
    void  saveImagePtr(float *,nifti_image *,char *);
    int   countLesionVoxels(int,int);
    int   getDebug();
    int   getVerbose();
    void  calculateEuclideanDistance(float *);

public:
    seg_fill_lesions_other();
    seg_fill_lesions_other(const seg_fill_lesions_other &o);
    int operator=(const seg_fill_lesions_other &o);
    ~seg_fill_lesions_other();

    void setInputImage(nifti_image *);
    void setInputLesionMask(nifti_image *);
    void setPatchSize(int);
    void setSearchArea(int);
    void setVerbose(int);
    void setDebug(int);
    void init();
    void runIt();
    void saveImage(nifti_image *,char *);
};

#include "_seg_fill_lesions_other.cpp"

#endif // _SEG_FILL_LESIONS_OTHER_H
