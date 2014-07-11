/**
 * @file seg_EM.cpp
 * @author M. Jorge Cardoso
 * @date 01/01/2014
 * @brief The main EM segmentation tool of nifty_seg. This function replaces seg_LoAd (seg_LoAd is deprecated).
 *
 * Copyright (c) 2014, University College London. All rights reserved.
 * Centre for Medical Image Computing (CMIC)
 * See the LICENSE.txt file in the nifty_seg root folder
 *
 */

#include "_seg_EM.h"
#include "seg_EM_CLIxml.h"

#include <iostream>
#include <time.h>
#include <stdlib.h>

using namespace std;
#define SegPrecisionTYPE float

// Executable usage message
void Usage(char *exec)
{
    printf("\nEM Statistical Segmentation:\nUsage ->\t%s -in <filename> [OPTIONS]\n\n",exec);
    printf("\t* * * * * * * * * * * * * * * * * Mandatory * * * * * * * * * * * * * * * * * *\n\n");
    printf("\t-in <filename>\t\t| Filename of the input image\n");
    printf("\t-out <filename>\t\t| Filename of the segmented image\n");
    printf("\t\t\t\t| The input image should be 2D, 3D or 4D images. 2D images should be on the XY plane.\n");
    printf("\t\t\t\t| 4D images are segmented as if they were multimodal.\n\n");
    printf("  \t\t- Select one of the following (mutually exclusive) -\n\n");
    printf("\t-priors <n> <fnames>\t| The number of priors (n>0) and their filenames. Priors should be registerd to the input image\n");
    printf("\t-priors4D <fname>\t| 4D image with the piors stacked in the 4th dimension. Priors should be registerd to the input image\n");
    printf("\t-nopriors <n>\t\t| The number of classes (n>0) \n\n");
    printf("\t* * * * * * * * * * * * * * * * General Options * * * * * * * * * * * * * * * * *\n\n");
    printf("\t-mask <filename>\t| Filename of the brain-mask of the input image\n");
    printf("\t-max_iter <int>\t\t| Maximum number of iterations (default = 100)\n");
    printf("\t-min_iter <int>\t\t| Minimum number of iterations (default = 0)\n");
    printf("\t-v <int>\t\t| Verbose level [0 = off, 1 = on, 2 = debug] (default = 0)\n");
    printf("\t-mrf_beta <float>\t| MRF prior strength [off = 0, max = 1] (default = 0.4) \n");
    printf("\t-bc_order <int>\t\t| Polynomial order for the bias field [off = 0, max = 5] (default = 3) \n");
    printf("\t-bc_thresh <float>\t| Bias field correction will run only if the ratio of improvement is below bc_thresh (default=0 [OFF]) \n");
    printf("\t-bc_out <filename>\t| Output the bias corrected image\n");
    printf("\t-reg <float>\t\t| Amount of regularization over the diagonal of the covariance matrix [above 1]\n");
    printf("\t-outlier <fl1> <fl2>\t| Outlier detection as in (Van Leemput TMI 2003). <fl1> is the Mahalanobis threshold [recommended between 3 and 7] \n");
    printf("\t\t\t\t| <fl2> is a convergence ratio below which the outlier detection is going to be done [recommended 0.01].\n");
    printf("\t-out_outlier <filename>\t| Output outlierness image \n");
    printf("\t-rf <rel> <gstd>\t| Relax Priors [relaxation factor: 0<rf<1 (recommended=0.5), gaussian regularization: gstd>0 (recommended=2.0)] /only 3D/\n");
    printf("\t-MAP <M V M V ...> \t| MAP formulation: M and V are the parameters (mean & variance) of the semiconjugate prior over the class mean\n");
#ifdef _GIT_HASH
    printf("\t--version\t\t|Print current source code git hash key and exit\n\t\t\t\t(%s)\n",_GIT_HASH);
#endif
    printf("\n\t* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    return;
}

void Merge_Priors(nifti_image * Priors, nifti_image ** Priors_temp, seg_EM_Params * segment_param)
{
    long img_size= Priors->nx * Priors->ny * Priors->nz;
    SegPrecisionTYPE * Prior_ptr_start = static_cast<SegPrecisionTYPE *>(Priors->data);
    for(long cl=0; cl<(segment_param->numb_classes); cl++)
    {
        SegPrecisionTYPE * Prior_tmp_ptr = static_cast<SegPrecisionTYPE *>(Priors_temp[cl]->data);
        SegPrecisionTYPE * Prior_ptr = &Prior_ptr_start[cl*img_size];
        for(long i=0; i<img_size; i++)
        {
            (*Prior_ptr)=(*Prior_tmp_ptr);
            Prior_ptr++;
            Prior_tmp_ptr++;
        }
    }
    return;
}

void no_memory ()
{
    cout << "Failed to allocate memory!\n";
    exit (1);
}

bool GetCommaSeparatedFloats(char *str, float &f1, float &f2)
{
    f1 = f2 = 0.;
    if (sscanf(str, "%f%*1[x,]%f", &f1, &f2) != 2)
    {
        std::cerr << "Incorrectly formatted string: '" << str << "' - should be '<a>,<b>'."
                  << std::endl;
        return false;
    }
    return true;
}

// Allocate and initialise to the arrays for MAP (if they haven't been allocated previously)
void AllocAndInitialiseMAP( float *&MAP_M, float *&MAP_V, bool &flag_MAP, long &numb_classes,
                            char *str, int iClass )
{
    // Has the number of classes been set large enough?
    if ( (long)iClass >= numb_classes )
    {
        // Array isn't allocated so we can just do it now
        if ( ! MAP_M )
        {
            MAP_M = (float*) calloc( iClass + 1, sizeof(float) );
        }
        // We have to extend the array and copy the values across
        else
        {
            float *New_MAP_M = (float*) calloc( iClass + 1, sizeof(float) );
            for(long classnum=0; classnum<numb_classes; classnum++)
                New_MAP_M[classnum]=MAP_M[classnum];

            free(MAP_M);
            MAP_M = New_MAP_M;
        }
        // Array isn't allocated so we can just do it now
        if ( ! MAP_V )
        {
            MAP_V = (float*) calloc( iClass + 1, sizeof(float) );
        }
        // We have to extend the array and copy the values across
        else
        {
            float *New_MAP_V = (float*) calloc( iClass + 1, sizeof(float) );
            for(long classnum=0; classnum<numb_classes; classnum++)
                New_MAP_V[classnum]=MAP_V[classnum];
            free(MAP_V);
            MAP_V = New_MAP_V;
        }
        numb_classes = iClass + 1;
    }
    // Yes we have sufficient classes but do we have to allocate the arrays?
    else
    {
        if ( ! MAP_M )
            MAP_M = (float*) calloc( numb_classes, sizeof(float) );
        if ( ! MAP_V )
            MAP_V = (float*) calloc( numb_classes, sizeof(float) );
    }
    if ( GetCommaSeparatedFloats(str, MAP_M[iClass], MAP_V[iClass]) )
        flag_MAP = true;
    else
        flag_MAP = false;
}


int main(int argc, char **argv)
{
    try
    {

        set_new_handler(no_memory);

        if (argc <=1)
        {
            Usage(argv[0]);
            return 0;
        }

        seg_EM_Params * segment_param = new seg_EM_Params [1]();

        // Set defaults (SEG_PARAM constructor sets everything to zero)
        segment_param->MRF_strength=0.4f;
        segment_param->bias_order=3;
        segment_param->flag_Bias=1;
        segment_param->flag_MRF=1;
        segment_param->flag_PV_model=1;
        segment_param->flag_SG_deli=1;
        segment_param->maxIteration=100;
        segment_param->minIteration=4;
        segment_param->OutliernessRatio=0.01;

        int To_do_MAP_index_argv=0;
        float regularization_amount=1;


        /* read the input parameter */
        for(long i=1; i<(long)argc; i++)
        {
            if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
                    strcmp(argv[i], "-HELP")==0 || strcmp(argv[i], "-h")==0 ||
                    strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0)
            {
                Usage(argv[0]);
                return 0;
            }
            else if( (strcmp(argv[i], "-xml")==0) || (strcmp(argv[i], "--xml")==0) )
            {
                cout << xml_segEM;
                return 0;
            }
            else if((strcmp(argv[i], "-in") == 0 || strcmp(argv[i], "--in") == 0) && (i+1)<(long)argc)
            {
                segment_param->filename_T1 = argv[++i];
                segment_param->flag_T1=1;
            }
            else if((strcmp(argv[i], "-out") == 0 || strcmp(argv[i], "--out") == 0) && (i+1)<(long)argc)
            {
                segment_param->filename_out = argv[++i];
                segment_param->flag_out=1;
            }
            else if((strcmp(argv[i], "-mask") == 0 || strcmp(argv[i], "--mask") == 0) && (i+1)<(long)argc)
            {
                segment_param->filename_mask=argv[++i];
                segment_param->flag_mask=1;
            }
            else if((strcmp(argv[i], "-priors") == 0 || strcmp(argv[i], "--priors") == 0) && (i+1)<(long)argc)
            {
                segment_param->numb_classes=atoi(argv[++i]);
                if(segment_param->numb_classes<2)
                {
                    cout<<"Number of classes has to be bigger than 1";
                    return 0;
                }
                if((i+segment_param->numb_classes)<(long)argc)
                {
                    segment_param->filename_priors= (char **) calloc(segment_param->numb_classes,sizeof(char *));
                    for(long classnum=0; classnum<segment_param->numb_classes; classnum++)
                    {
                        segment_param->filename_priors[classnum]=argv[++i];
                    }
                }
                else
                {
                    fprintf(stderr,"Err:\tParameter -priors are incomplete\n");
                    Usage(argv[0]);
                    return 1;
                }
                segment_param->flag_manual_priors=1;
            }
            else if((strcmp(argv[i], "-priors4D") == 0 || strcmp(argv[i], "--priors4D") == 0) && (i+1)<(long)argc)
            {
                segment_param->filename_priors= (char **) calloc(1,sizeof(char *));
                segment_param->filename_priors[0]=argv[++i];
                segment_param->flag_manual_priors=1;
                nifti_image * tmpread = nifti_image_read(segment_param->filename_priors[0],false);
                segment_param->numb_classes=tmpread->nt;
                if(segment_param->numb_classes<2)
                {
                    cout<<"Number of classes has to be bigger than 1";
                    return 0;
                }
                nifti_image_free(tmpread);
                segment_param->flag_priors4D=true;
            }
            else if((strcmp(argv[i], "-nopriors") == 0 || strcmp(argv[i], "--nopriors") == 0) && (i+1)<(long)argc)
            {
                if(!segment_param->numb_classes)
                    segment_param->numb_classes=atoi(argv[++i]);
            }
            else if((strcmp(argv[i], "-outlier") == 0 || strcmp(argv[i], "--outlier") == 0) && (i+2)<(long)argc)
            {
                segment_param->flag_Outlierness=1;
                segment_param->OutliernessThreshold=atof(argv[++i]);
                segment_param->OutliernessRatio=atof(argv[++i]);
            }
            // Additional options "-outlier_thresh" and "-outlier_ratio"
            // added for operation with command line interface XML
            // interface
            else if((strcmp(argv[i], "-outlier_thresh") == 0 || strcmp(argv[i], "--outlier_thresh") == 0) && (i+1)<(long)argc)
            {
                segment_param->OutliernessThreshold=atof(argv[++i]);
                if (segment_param->OutliernessThreshold)
                    segment_param->flag_Outlierness=1;
            }
            else if((strcmp(argv[i], "-outlier_ratio") == 0 || strcmp(argv[i], "--outlier_ratio") == 0) && (i+1)<(long)argc)
            {
                segment_param->OutliernessRatio=atof(argv[++i]);
            }
            else if((strcmp(argv[i], "-out_outlier") == 0 || strcmp(argv[i], "--out_outlier") == 0) && (i+1)<(long)argc)
            {
                segment_param->filename_out_outlier = argv[++i];
                segment_param->flag_out_outlier=1;
            }
            // Additional options "-MAP_MV1/2/3/4/5" added for operation
            // with command line interface XML interface. Format: -MAP_MV1 <mean>,<var> etc.
            else if(strcmp(argv[i], "-MAP_MV1") == 0 || strcmp(argv[i], "--MAP_MV1") == 0)
            {
                if((i+1)<(long)argc)
                {
                    AllocAndInitialiseMAP(segment_param->MAP_M, segment_param->MAP_V,
                                          segment_param->flag_MAP, segment_param->numb_classes,
                                          argv[++i], 0);
                }
            }
            else if(strcmp(argv[i], "-MAP_MV2") == 0 || strcmp(argv[i], "--MAP_MV2") == 0)
            {
                if((i+1)<(long)argc)
                {
                    AllocAndInitialiseMAP(segment_param->MAP_M, segment_param->MAP_V,
                                          segment_param->flag_MAP, segment_param->numb_classes,
                                          argv[++i], 1);
                }
            }
            else if(strcmp(argv[i], "-MAP_MV3") == 0 || strcmp(argv[i], "--MAP_MV3") == 0)
            {
                if((i+1)<(long)argc)
                {
                    AllocAndInitialiseMAP(segment_param->MAP_M, segment_param->MAP_V,
                                          segment_param->flag_MAP, segment_param->numb_classes,
                                          argv[++i], 2);
                }
            }
            else if(strcmp(argv[i], "-MAP_MV4") == 0 || strcmp(argv[i], "--MAP_MV4") == 0)
            {
                if((i+1)<(long)argc)
                {
                    AllocAndInitialiseMAP(segment_param->MAP_M, segment_param->MAP_V,
                                          segment_param->flag_MAP, segment_param->numb_classes,
                                          argv[++i], 3);
                }
            }
            else if(strcmp(argv[i], "-MAP_MV5") == 0 || strcmp(argv[i], "--MAP_MV5") == 0)
            {
                if((i+1)<(long)argc)
                {
                    AllocAndInitialiseMAP(segment_param->MAP_M, segment_param->MAP_V,
                                          segment_param->flag_MAP, segment_param->numb_classes,
                                          argv[++i], 4);
                }
            }
            // This option overides "-MAP_MV1/2/3/4/5"
            else if(strcmp(argv[i], "-MAP") == 0 || strcmp(argv[i], "--MAP") == 0)
            {
                if(segment_param->numb_classes>0 && (i+2*segment_param->numb_classes)<(long)argc)
                {
                    segment_param->MAP_M= (float*) calloc(segment_param->numb_classes,sizeof(float));
                    segment_param->MAP_V= (float*) calloc(segment_param->numb_classes,sizeof(float));
                    for(long classnum=0; classnum<segment_param->numb_classes; classnum++)
                    {
                        segment_param->MAP_M[classnum]=atof(argv[++i]);
                        segment_param->MAP_V[classnum]=atof(argv[++i]);
                    }
                }
                else
                {
                    To_do_MAP_index_argv=i;
                }
            }
            else if((strcmp(argv[i], "-bc_order") == 0 || strcmp(argv[i], "--bc_order") == 0) && (i+1)<(long)argc)
            {
                segment_param->bias_order=(int)(atof(argv[++i]));
                if(segment_param->bias_order==0)
                {
                    segment_param->flag_Bias=0;
                }
            }
            else if((strcmp(argv[i], "-bc_thresh") == 0 || strcmp(argv[i], "--bc_thresh") == 0) && (i+1)<(long)argc)
            {
                segment_param->Bias_threshold=atof(argv[++i]);
            }
            else if((strcmp(argv[i], "-reg") == 0 || strcmp(argv[i], "--reg") == 0) && (i+1)<(long)argc)
            {
                regularization_amount=atof(argv[++i]);
            }
            else if((strcmp(argv[i], "-rf") == 0 || strcmp(argv[i], "--rf") == 0) && (i+1)<(long)argc)
            {
                segment_param->relax_factor=atof(argv[++i]);
                segment_param->relax_gauss_kernel=atof(argv[++i]);
            }
            else if((strcmp(argv[i], "-rf_factor") == 0 || strcmp(argv[i], "--rf_factor") == 0) && (i+1)<(long)argc)
            {
                segment_param->relax_factor=atof(argv[++i]);
            }
            else if((strcmp(argv[i], "-rf_gauss_kernel") == 0 || strcmp(argv[i], "--rf_gauss_kernel") == 0) && (i+1)<(long)argc)
            {
                segment_param->relax_gauss_kernel=atof(argv[++i]);
            }
            //else if(strcmp(argv[i], "-aprox_off") == 0 || strcmp(argv[i], "--aprox_off") == 0){
            //    segment_param->aprox=false;
            //}

            else if((strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--v") == 0) &&(i+1)<(long)argc)
            {
                segment_param->verbose_level=(int)atoi(argv[++i]);
            }
            else if((strcmp(argv[i], "-bc_out") == 0 || strcmp(argv[i], "--bc_out") == 0) && (i+1)<(long)argc)
            {
                segment_param->filename_bias = argv[++i];
                segment_param->flag_bc_out=1;
            }
            else if((strcmp(argv[i], "-mrf_beta") == 0 || strcmp(argv[i], "--mrf_beta") == 0) && (i+1)<(long)argc)
            {
                segment_param->MRF_strength=(SegPrecisionTYPE)atof(argv[++i]);

            }
            else if((strcmp(argv[i], "-max_iter") == 0 || strcmp(argv[i], "--max_iter") == 0) && (i+1)<(long)argc)
            {
                segment_param->maxIteration=atoi(argv[++i]);
            }
            else if((strcmp(argv[i], "-min_iter") == 0 || strcmp(argv[i], "--min_iter") == 0) && (i+1)<(long)argc)
            {
                segment_param->minIteration=atoi(argv[++i]);
            }
#ifdef _GIT_HASH
            else if( strcmp(argv[i], "--version")==0)
            {
                printf("%s\n",_GIT_HASH);
                return 0;
            }
#endif
            else
            {
                fprintf(stderr,"Err:\tParameter %s unknown or incomplete \n",argv[i]);
                Usage(argv[0]);
                return 1;
            }

        }
        if(segment_param->MRF_strength<=0)
        {
            segment_param->flag_MRF=0;
        }

        //string envVarName="NIFTYSEG";
        //char * currpath_char=NULL;
        //currpath_char=getenv(envVarName.c_str());


        if(To_do_MAP_index_argv>0)
        {
            int i=To_do_MAP_index_argv;
            if(segment_param->numb_classes>0 && (i+2*segment_param->numb_classes)<(long)argc)
            {
                segment_param->MAP_M= (float*) calloc(segment_param->numb_classes,sizeof(float));
                segment_param->MAP_V= (float*) calloc(segment_param->numb_classes,sizeof(float));
                for(long classnum=0; classnum<segment_param->numb_classes; classnum++)
                {
                    segment_param->MAP_M[classnum]=atof(argv[++i]);
                    segment_param->MAP_V[classnum]=atof(argv[++i]);
                }
                segment_param->flag_MAP=true;
            }
        }

        // Check that all the MAP value pairs are non-zero
        if (segment_param->MAP_M && segment_param->MAP_V)
        {
            for(long classnum=0; classnum<segment_param->numb_classes; classnum++)
            {
                if ((! segment_param->MAP_M[classnum]) && (! segment_param->MAP_V[classnum]))
                {
                    fprintf(stderr,"Err:\tMAP value pair (class %d) are zero.\n",(int)classnum+1);
                    Usage(argv[0]);
                    return 1;
                }
            }
        }

        if(!segment_param->flag_T1)
        {
            fprintf(stderr,"Err:\tThe T1 image name has to be defined.\n");
            Usage(argv[0]);
            return 1;
        }

        if(segment_param->flag_out==0)
        {
            fprintf(stderr,"Err:\tThe output image name has to be defined.\n");
            Usage(argv[0]);
            return 1;
        }

        // READING T1
        nifti_image * InputImage=nifti_image_read(segment_param->filename_T1,true);
        if(InputImage == NULL)
        {
            fprintf(stderr,"* Error when reading the T1 image: %s\n",segment_param->filename_T1);
            return 1;
        }
        if(InputImage->datatype!=NIFTI_TYPE_FLOAT32)
            seg_changeDatatype<SegPrecisionTYPE>(InputImage);

        InputImage->dim[4]=InputImage->nt=(InputImage->nt<1)?1:InputImage->nt;
        InputImage->dim[5]=InputImage->nu=(InputImage->nu<1)?1:InputImage->nu;
        if(InputImage->nu>1)
        {
            InputImage->dim[5]=InputImage->nu;
            InputImage->dim[4]=1;
        }
        else if(InputImage->nt>1)
        {
            InputImage->dim[5]=InputImage->nt;
            InputImage->dim[4]=1;
        }
        else if(InputImage->nt>1)
        {
            InputImage->dim[5]=1;
            InputImage->dim[4]=1;
        }
        InputImage->dim[0]=5;
        nifti_update_dims_from_array(InputImage);


        nifti_image * Mask=NULL;
        if(segment_param->flag_mask)
        {
            Mask = nifti_image_read(segment_param->filename_mask,true);
            if(Mask->datatype!=NIFTI_TYPE_FLOAT32)
                seg_changeDatatype<SegPrecisionTYPE>(Mask);

            if(Mask == NULL)
            {
                fprintf(stderr,"* Error when reading the mask image: %s\n",segment_param->filename_mask);
                return 1;
            }
            if( (Mask->nx != InputImage->nx) || (Mask->ny != InputImage->ny) || (Mask->nz != InputImage->nz) )
            {
                fprintf(stderr,"* Error: Mask image not the same size as input\n");
                return 1;
            }
            Mask->dim[4]=Mask->nt=(Mask->nt<1)?1:Mask->nt;
            Mask->dim[5]=Mask->nu=(Mask->nu<1)?1:Mask->nu;
            nifti_update_dims_from_array(Mask);
        }

        nifti_image ** Priors_temp=new nifti_image * [segment_param->numb_classes];
        nifti_image * Priors=NULL;

        if(segment_param->flag_manual_priors)
        {

            if(segment_param->flag_priors4D)
            {
                int i=0;
                Priors = nifti_image_read(segment_param->filename_priors[i],true);
                if(Priors == NULL)
                {
                    fprintf(stderr,"* Error when reading the Prior image: %s\n",segment_param->filename_priors[i]);
                    return 1;
                }
                if( (Priors->nx != InputImage->nx) || (Priors->ny != InputImage->ny) || (Priors->nz != InputImage->nz) )
                {
                    fprintf(stderr,"* Error: Prior image ( %s ) not the same size as input\n",segment_param->filename_priors[i]);
                    return 1;
                }
                seg_changeDatatype<SegPrecisionTYPE>(Priors);

            }
            else
            {

                for(long i=0; i<segment_param->numb_classes; i++)
                {
                    Priors_temp[i] = nifti_image_read(segment_param->filename_priors[i],true);

                    if(Priors_temp[i] == NULL)
                    {
                        fprintf(stderr,"* Error when reading the Prior image: %s\n",segment_param->filename_priors[i]);
                        return 1;
                    }
                    if( (Priors_temp[i]->nx != InputImage->nx) || (Priors_temp[i]->ny != InputImage->ny) || (Priors_temp[i]->nz != InputImage->nz) )
                    {
                        fprintf(stderr,"* Error: Prior image ( %s ) not the same size as input\n",segment_param->filename_priors[i]);
                        return 1;
                    }
                    seg_changeDatatype<SegPrecisionTYPE>(Priors_temp[i]);
                }

                Priors=nifti_copy_nim_info(InputImage);
                Priors->dim[0]=4;
                Priors->nt=Priors->dim[4]=segment_param->numb_classes;
                Priors->datatype=NIFTI_TYPE_FLOAT32;
                Priors->cal_max=1;

                nifti_update_dims_from_array(Priors);
                nifti_datatype_sizes(Priors->datatype,&Priors->nbyper,&Priors->swapsize);
                Priors->data = (void *) calloc(Priors->nvox, sizeof(SegPrecisionTYPE));

                Merge_Priors(Priors,Priors_temp,segment_param);
                for(long i=0; i<segment_param->numb_classes; i++)
                {
                    nifti_image_free(Priors_temp[i]);
                    Priors_temp[i]=NULL;
                }
            }
        }
        else
        {
            if(segment_param->numb_classes<2)
            {
                fprintf(stderr,"Please use either the -priors or -nopriors options with at least 2 classes\n");
                return 1;
            }
            else
            {
                if(segment_param->verbose_level>0)
                    printf("WARNING: No priors selected. Results might be inconsistent\n");
            }

        }

        delete [] Priors_temp;

        if ( segment_param->verbose_level>1 )
            segment_param->Print( cout );

        seg_EM SEG(segment_param->numb_classes,InputImage->dim[5],InputImage->dim[4]);
        SEG.SetInputImage(InputImage);
        if(segment_param->flag_mask)
            SEG.SetMaskImage(Mask);
        if(segment_param->flag_manual_priors)
            SEG.SetPriorImage(Priors);

        SEG.SetVerbose(segment_param->verbose_level);
        SEG.SetFilenameOut(segment_param->filename_out);
        SEG.SetMaximalIterationNumber(segment_param->maxIteration);
        SEG.SetMinIterationNumber(segment_param->minIteration);

        if(segment_param->flag_Outlierness)
            SEG.SetOutlierness(segment_param->OutliernessThreshold,segment_param->OutliernessRatio);
        if(segment_param->flag_Bias)
            SEG.SetBiasField(segment_param->bias_order,segment_param->Bias_threshold);
        if(segment_param->flag_MRF)
            SEG.SetMRF(segment_param->MRF_strength);
        if(segment_param->relax_factor>0)
            SEG.SetRelaxation(segment_param->relax_factor,segment_param->relax_gauss_kernel);
        if(segment_param->flag_MAP)
            SEG.SetMAP(segment_param->MAP_M,segment_param->MAP_V);
        if(regularization_amount>0)
            SEG.SetRegValue(regularization_amount);

        SEG.Run_EM();

        if(segment_param->verbose_level>0)
        {
            cout << "Saving Segmentation to file: " << segment_param->filename_out <<endl;
        }
        nifti_image * Result = SEG.GetResult();
        nifti_image_write(Result);
        nifti_image_free(Result);
        nifti_image_free(Priors);


        nifti_image * OutliernessImage=NULL;
        if(segment_param->flag_out_outlier && segment_param->flag_Outlierness)
        {
            if(segment_param->verbose_level>0)
            {
                cout << "Saving Outlierness Image to: " << segment_param->filename_out_outlier <<endl;
            }
            OutliernessImage=SEG.GetOutlierness(segment_param->filename_out_outlier);
            nifti_image_write(OutliernessImage);
            nifti_image_free(OutliernessImage);
        }

        nifti_image * BiasFieldCorrected=NULL;

        if(segment_param->flag_bc_out)
        {
            if(segment_param->verbose_level>0)
            {
                cout << "Saving Bias Field Corrected Image to file: " << segment_param->filename_bias <<endl;
            }
            BiasFieldCorrected = SEG.GetBiasCorrected(segment_param->filename_bias);
            nifti_image_write(BiasFieldCorrected);
            nifti_image_free(BiasFieldCorrected);
        }

        nifti_image_free(InputImage);
        nifti_image_free(Mask);

        free(segment_param->filename_priors);

        delete [] segment_param;
    }

    catch(std::exception & e)
    {
        std::cerr << "Standard exception: " << e.what() << std::endl;
    }

    catch(...)
    {
        std::cerr << "Unhandled Exception: Something went wrong! Please report the error to mjorgecardoso"<<(char) 64<<"gmail.com" << std::endl;
    }
    return 0;
}

