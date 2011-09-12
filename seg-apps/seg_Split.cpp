#include "_seg_common.h"
#include <iostream>
#include <time.h>
using namespace std;
#define PrecisionTYPE float

void Usage(char *exec)
{
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    printf("Usage:\t%s -in <filename> -out <filename>.\n\n",exec);
    printf("\t* * Mandatory * *\n");
    printf("\t-in <filename>\t\tFilename of the input image segmentation\n\n");
    printf("\t-out <filename>\t\tFilename of the brainmask of the input image\n");
    printf("\t* * Options * *\n");
    printf("\t-mode <int>\t\tOutput Mode [0 = Bias Corrected Image, 1 = Cortical GM, 2 = WM,  3 = Brain] (default = 0)\n");
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    return;
}

int main(int argc, char **argv)
{
    char * filename_out=NULL;
    char * filename_in=NULL;
    int mode=0;
    for(int i=1;i<argc;i++){
        if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
           strcmp(argv[i], "-HELP")==0 || strcmp(argv[i], "-h")==0 ||
           strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0){
            Usage(argv[0]);
            return 0;
        }

        else if(strcmp(argv[i], "-in") == 0){
            filename_in=argv[++i];
        }
        else if(strcmp(argv[i], "-out") == 0){
            filename_out = argv[++i];
        }
        else if(strcmp(argv[i], "-mode") == 0){
            mode=(int)atoi(argv[++i]);
        }
    }



    if(filename_in == NULL){
        fprintf(stderr,"* Error: No input defined\n");
        return 1;
    }

    nifti_image * Segmentation=nifti_image_read(filename_in,true);
    if(Segmentation == NULL){
        fprintf(stderr,"* Error when reading the input Segmentation image\n");
        return 1;
    }


    if(filename_out == NULL){
        fprintf(stderr,"* Error: No output defined\n");
        return 1;
    }



    nifti_image * Result = nifti_copy_nim_info(Segmentation);
    Result->dim[0]=3;
    Result->dim[4]=0;
    Result->datatype=DT_FLOAT32;
    Result->cal_max=1;
    nifti_update_dims_from_array(Result);
    nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
    nifti_set_filenames(Result,filename_out,0,0);
    Result->data = (void *) calloc(Result->nvox, sizeof(PrecisionTYPE));
    PrecisionTYPE * Result_PTR = static_cast<PrecisionTYPE *>(Result->data);
    PrecisionTYPE * Segmentation_PTR = static_cast<PrecisionTYPE *>(Segmentation->data);


    if(mode==3){
        for(unsigned int i=0; i<Result->nvox; i++){
        Result_PTR[i]=Segmentation_PTR[i+(WMclass+1)*Result->nvox]+Segmentation_PTR[i+(GMclass+1)*Result->nvox]+Segmentation_PTR[i+(dGMclass+1)*Result->nvox]+Segmentation_PTR[i+(iCSFclass+1)*Result->nvox];
        }

    }
    else if(mode==2){
        for(unsigned int i=0; i<Result->nvox; i++){
        Result_PTR[i]=Segmentation_PTR[i+(WMclass+1)*Result->nvox];
        }

    }
    else if(mode==1){
        for(unsigned int i=0; i<Result->nvox; i++){
        Result_PTR[i]=Segmentation_PTR[i+(GMclass+1)*Result->nvox];
        }

    }
    else if(mode==0){
        float max=-1000000;
        float min=1000000;
        for(unsigned int i=0; i<Result->nvox; i++){
        Result_PTR[i]=Segmentation_PTR[i];
        if(Result_PTR[i]>max){
            max=Result_PTR[i];
        }
        if(Result_PTR[i]<min){
            min=Result_PTR[i];
        }
        Result->cal_max=max;
        Result->cal_min=min;
        }

    }
    else{
        fprintf(stderr,"* Error: Mode not implemented\n");
        return 1;
    }
    // Finish Test Area
    nifti_image_write(Result);
    nifti_image_free(Result);
    nifti_image_free(Segmentation);
    return 0;
}
