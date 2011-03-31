#include "_seg_common.h"
#include "_seg_BiasCorrection.h"
#include "_seg_tools.h"
#include "_seg_MRF.h"
#include "_seg_FMM.h"
#include <iostream>
#include <time.h>
using namespace std;
#define PrecisionTYPE float



int main(int argc, char **argv)
{


    nifti_image * T1=nifti_image_read(argv[1],true);
    if(T1 == NULL){
        fprintf(stderr,"* Error when reading the T1 image\n");
        return 1;
    }

    nifti_image * Mask = nifti_image_read(argv[2],true);
    if(Mask == NULL){
        fprintf(stderr,"* Error when reading the mask image\n");
        return 1;
    }
    nifti_image * img1 = nifti_image_read(argv[3],true);
    if(img1 == NULL){
        fprintf(stderr,"* Error when reading the image1\n");
        return 1;
    }
    nifti_image * img2 = nifti_image_read(argv[4],true);
    if(img2 == NULL){
        fprintf(stderr,"* Error when reading the image2\n");
        return 1;
    }
    char * outputname = argv[5];


    seg_changeDatatype<PrecisionTYPE>(T1);
    seg_changeDatatype<PrecisionTYPE>(Mask);
    seg_changeDatatype<PrecisionTYPE>(img1);
    seg_changeDatatype<PrecisionTYPE>(img2);

    int NumClass=non_PV_numclass;
    ImageSize * CurrSizes = new ImageSize [1]();
    CurrSizes->numel=(int)(T1->nx*T1->ny*T1->nz);
    CurrSizes->xsize=T1->nx;
    CurrSizes->ysize=T1->ny;
    CurrSizes->zsize=T1->nz;
    CurrSizes->numclass=5;
    CurrSizes->numelmasked=0;
    CurrSizes->numelbias=0;

    if(Mask->datatype!=DT_BINARY){seg_convert2binary(Mask,0);}
    Normalize_Image_mask(T1,Mask,CurrSizes,1);
    int * Short_2_Long_Indices = Create_Short_2_Long_Matrix_from_NII(Mask,&(CurrSizes->numelmasked));
    int * Long_2_Short_Indices = Create_Long_2_Short_Matrix_from_NII(Mask);

    PrecisionTYPE * img1_short = Create_cArray_from_3D_image(Mask,img1);
    PrecisionTYPE * img2_short = Create_cArray_from_3D_image(Mask,img2);

    // TEST AREA

    bool * img1_short_thresholded=binarise_image(img1_short,0.8,CurrSizes);
    PrecisionTYPE * GeoTime = new PrecisionTYPE [CurrSizes->numelmasked]();

    for(int i=0; i<(CurrSizes->numelmasked); i++){img2_short[i]=0.0f;}
    FMM(img1_short_thresholded, img2_short, GeoTime,100, Long_2_Short_Indices, Short_2_Long_Indices, CurrSizes);

    string str_out="Test.nii";
    nifti_image * Test = Copy_Single_ShortImage_to_Result(GeoTime,Short_2_Long_Indices,T1,outputname,CurrSizes);

    delete [] GeoTime;

    // Finish Test Area

    nifti_image_write(Test);
    nifti_image_free(Test);
    nifti_image_free(img1);
    nifti_image_free(img2);
    nifti_image_free(T1);
    nifti_image_free(Mask);
    delete [] CurrSizes;
    delete [] img1_short_thresholded;
    delete [] img1_short;
    delete [] img2_short;
    delete [] Short_2_Long_Indices;
    delete [] Long_2_Short_Indices;

    return 0;
}

