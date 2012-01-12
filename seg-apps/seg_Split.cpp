#include "_seg_common.h"
#include <iostream>
#include <time.h>
#include "_seg_common.h"
#include "_seg_tools.h"
#include "_seg_Topo.h"

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
    printf("\t-tp <int>\t\tTimepoint\n");
    printf("\t-bin\t\tInteger image with the highest probability value.\n");
    printf("\t-fill\t\tFill masked image.\n");
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    return;
}
int main(int argc, char **argv)
{
    char * filename_out=NULL;
    char * filename_in=NULL;
    int mode=0;
    int tp=-1;
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
        else if(strcmp(argv[i], "-tp") == 0){
            mode=0;
            tp=(int)atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-bin") == 0){
            mode=1;
        }
        else if(strcmp(argv[i], "-fill") == 0){
            mode=2;
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
    seg_changeDatatype<PrecisionTYPE>(Segmentation);


    if(mode==0){
        for(unsigned int i=0; i<Result->nvox; i++){
            Result_PTR[i]=Segmentation_PTR[i+(tp)*Result->nx*Result->ny*Result->nz];
        }
    }
    if(mode==1){
        for(int i=0; i<(Result->nx*Result->ny*Result->nz); i++){
            float lab_at_index=0;
            float val_at_index=0.01;
            for(int cl=0; cl<Segmentation->nt; cl++){
                if(Segmentation_PTR[i+(cl)*(Result->nx*Result->ny*Result->nz)]>val_at_index){
                    val_at_index=Segmentation_PTR[i+cl*(Result->nx*Result->ny*Result->nz)];
                    lab_at_index=(int)cl+1;
                    //cout<<val_at_index<<" "<<lab_at_index<<endl;
                }
                //cout<<val_at_index<<" "<<lab_at_index<<endl;
            }
            Result_PTR[i]=lab_at_index;
        }
    }
    if(mode==2){
        int dimensions[3];
        dimensions[0]=Segmentation->nx;
        dimensions[1]=Segmentation->ny;
        dimensions[2]=Segmentation->nz;

        seg_changeDatatype<int>(Segmentation);
        seg_changeDatatype<int>(Result);
        int * Result_PTR2 = static_cast<int *>(Result->data);
        int * Segmentation_PTR2 = static_cast<int *>(Segmentation->data);
        Close_Forground_ConnectComp(Segmentation_PTR2,Result_PTR2,dimensions);
    }



    // Finish Test Area
    nifti_image_write(Result);
    nifti_image_free(Result);
    nifti_image_free(Segmentation);
    return 0;
}
