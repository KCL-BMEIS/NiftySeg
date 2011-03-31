#include "_seg_CT_lungs.h"
#include <iostream>
#include <time.h>


void Usage(char *exec)
{
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    printf("Usage:\t%s -in <filename> -out <filename>.\n\n",exec);
    printf("\t* * Mandatory * *\n");
    printf("\t-in <filename>\t\tFilename of the input image segmentation\n\n");
    printf("\t-mask <filename>\t\tFilename of the Lung mask\n\n");
    printf("\t-out <filename>\t\tFilename of the brainmask of the input image\n");

    printf("\t* * Options * *\n");
    printf("\t-max_iter <int>\t\t\tMaximum number of iterations (default = 100)\n");
    printf("\t-v <int>\t\t\tVerbose level [0 = off, 1 = verbose, 2 = debug] (default = 0)\n");
    printf("\t-mrf_beta <float>\t\tMRF prior strength [off = 0] (default = 0.2) \n");
    printf("\t-bc_order <int>\t\t\tPolynomial order for the bias field [off = 0, max = 5] (default = 4) \n");
    //printf("\t* * Options * *\n");
    //printf("\t-mode <int>\t\tOutput Mode [0 = Split All (Default), 1 = Cortical GM, 2 = Brain] (default = 0)\n");
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    return;
}


int main(int argc, char **argv)
{


    SEG_PARAM * segment_param = new SEG_PARAM [1]();
    //SEG_PARAM segment_param;
    segment_param->maxIteration=100;
    segment_param->flag_T1=0;
    segment_param->flag_out=0;
    segment_param->flag_mask=0;
    segment_param->flag_MRF=1;
    segment_param->flag_Bias=1;
    segment_param->flag_SG_deli=1;
    segment_param->flag_bc_out=0;
    segment_param->relax_factor=1;
    segment_param->flag_PV_model=0;
    segment_param->verbose_level=0;
    segment_param->flag_manual_priors=0;
    segment_param->bias_order=4;
    segment_param->MRF_strength=0.2f;
    segment_param->numb_classes=2;



    char * filename_out=NULL;
    char * filename_in=NULL;
    char * filename_mask=NULL;
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
        else if(strcmp(argv[i], "-mask") == 0){
            filename_mask=argv[++i];
        }
        else if(strcmp(argv[i], "-out") == 0){
            filename_out = argv[++i];
        }
        else if(strcmp(argv[i], "-bc_order") == 0){
            segment_param->bias_order=(int)(atof(argv[++i]));
            if(segment_param->bias_order>maxallowedpowerorder){
                cout << "Maximum polynomial basis function order = 6. \n Assuming bc_order = 6" << endl;
                segment_param->bias_order=6;
            }
            if(segment_param->bias_order==0){
            segment_param->flag_Bias=0;
            }
        }
        else if(strcmp(argv[i], "-v") == 0){
            segment_param->verbose_level=(int)atoi(argv[++i]);
        }
       /* else if(strcmp(argv[i], "-bc_out") == 0){
            segment_param->filename_bc_T1=argv[++i];
            segment_param->flag_bc_out=1;
        }*/
        else if(strcmp(argv[i], "-mrf_beta") == 0){
            segment_param->MRF_strength=(PrecisionTYPE)atof(argv[++i]);
            if(segment_param->bias_order==0){
            segment_param->flag_MRF=0;
            }
        }
        else if(strcmp(argv[i], "-max_iter") == 0){
            segment_param->maxIteration=atoi(argv[++i]);
        }
        else{
            fprintf(stderr,"Err:\tParameter %s unknown->\n",argv[i]);
            Usage(argv[0]);
            return 1;
        }
    }


    if(filename_in == NULL){
        fprintf(stderr,"* Error: No input defined\n");
        return 1;
    }

    nifti_image * LungCT=nifti_image_read(filename_in,true);
    if(LungCT == NULL){
        fprintf(stderr,"* Error when reading the input Segmentation image\n");
        return 1;
    }
    seg_changeDatatype<PrecisionTYPE>(LungCT);

    nifti_image * Mask=nifti_image_read(filename_mask,true);
    if(Mask == NULL){
        fprintf(stderr,"* Error when reading the input Mask image\n");
        return 1;
    }
    seg_changeDatatype<PrecisionTYPE>(Mask);

    if(filename_out == NULL){
        fprintf(stderr,"* Error: No output defined\n");
        return 1;
    }

     nifti_image * Result = CT_Lung(LungCT, Mask, filename_out, segment_param);

    nifti_image_write(Result);
    nifti_image_free(Result);
    nifti_image_free(LungCT);
    return 0;
}


