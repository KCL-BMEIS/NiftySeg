/**
 * @file seg_LoAd.cpp
 * @author M. Jorge Cardoso
 * @date 01/01/2014
 *
 * Copyright (c) 2014, University College London. All rights reserved.
 * Centre for Medical Image Computing (CMIC)
 * See the LICENSE.txt file in the nifty_seg root folder
 *
 */

#include "_seg_LoAd.h"
#include "_seg_EM.h"
#include <iostream>
#include <time.h>
using namespace std;
#define SegPrecisionTYPE float


void Usage(char *exec)
{
    printf("\nLOAD Statistical Segmentation:\nUsage:\t%s -in <filename> [OPTIONS].\n\n",exec);
    printf("\t* * Mandatory * *\n");
    printf("\t-in <filename>\t\tFilename of the input image image\n");
    printf("\t-mask <filename>\tFilename of the brainmask of the input image\n\n");
    printf("\t-out <filename>\t\tFilename of for the segmented images\n");
    printf("\t* * Options * *\n");
    printf("\t-priors <filenames> \tThe 5 priors in this order: WM,GM,CSF,dGM,iCSF\n");
    printf("\t-max_iter <int>\t\tMaximum number of iterations (default = 100)\n");
    printf("\t-rf <int>\t\tRelaxation factor 0<RF<1 (default = 1)\n");
    printf("\t-v <int>\t\tVerbose level [0 = off, 1 = verbose, 2 = debug] (default = 0)\n");
    printf("\t-mrf_beta <float>\tMRF prior strength [off = 0] (default = 0.1) \n");
    printf("\t-bc_order <int>\t\tPolinomial order for the bias field [off = 0, max = 6] (default = 5) \n");
    printf("\t-pv_off \t\tDo not preform the PV modeling \n");
    printf("\t-sg_off <int>\t\tDo not improve sulci and gyri deliniation\n");
#ifdef _GIT_HASH
    printf("\t--version\t\t|Print current source code git hash key and exit\n\t\t\t\t(%s)\n",_GIT_HASH);
#endif
    printf("\n\t* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    return;
}

void Merge_Priors(nifti_image * Priors, nifti_image ** Priors_temp)
{
    long img_size= Priors->nx * Priors->ny * Priors->nz;
    SegPrecisionTYPE * Prior_ptr_start = static_cast<SegPrecisionTYPE *>(Priors->data);
    for(long cl=0; cl<(nonPVNumClass); cl++)
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

int main(int argc, char **argv)
{
    try
    {
        set_new_handler(no_memory);



        if (argc < 2)
        {
            Usage(argv[0]);
            return 0;
        }

        seg_EM_Params * segment_param = new seg_EM_Params [1]();
        //SEG_PARAM segment_param;

        // Set defaults (SEG_PARAM constructor sets everything to zero)
        segment_param->MRF_strength=0.1f;
        segment_param->bias_order=5;
        segment_param->flag_Bias=1;
        segment_param->flag_MRF=1;
        segment_param->flag_PV_model=1;
        segment_param->flag_SG_deli=1;
        segment_param->maxIteration=100;
        segment_param->relax_factor=1.0f;

        /* read the input parameter */
        for(int i=1; i<argc; i++)
        {
            if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
                    strcmp(argv[i], "-HELP")==0 || strcmp(argv[i], "-h")==0 ||
                    strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0)
            {
                Usage(argv[0]);
                return 0;
            }

            else if(strcmp(argv[i], "-in") == 0)
            {
                segment_param->filename_T1 = argv[++i];
                segment_param->flag_T1=1;
            }
            else if(strcmp(argv[i], "-out") == 0)
            {
                segment_param->filename_out = argv[++i];
                segment_param->flag_out=1;
            }
            else if(strcmp(argv[i], "-mask") == 0)
            {
                segment_param->filename_mask=argv[++i];
                segment_param->flag_mask=1;
            }
            else if(strcmp(argv[i], "-priors") == 0)
            {
                segment_param->numb_classes=5;
                segment_param->filename_priors= (char **) calloc(segment_param->numb_classes,sizeof(char **));
                for(int classnum=0; classnum<segment_param->numb_classes; classnum++)
                {
                    segment_param->filename_priors[classnum]=argv[++i];
                }
                segment_param->flag_manual_priors=1;
                segment_param->flag_manual_priors=1;
            }
            else if(strcmp(argv[i], "-bc_order") == 0)
            {
                segment_param->bias_order=(int)(atof(argv[++i]));
                if(segment_param->bias_order==0)
                {
                    segment_param->flag_Bias=0;
                }
            }
            else if(strcmp(argv[i], "-pv_off") == 0)
            {
                segment_param->flag_PV_model=0;
            }
            else if(strcmp(argv[i], "-v") == 0)
            {
                segment_param->verbose_level=(int)atoi(argv[++i]);
            }
            /* else if(strcmp(argv[i], "-bc_out") == 0){
                segment_param->filename_bc_T1=argv[++i];
                segment_param->flag_bc_out=1;
            }*/
            else if(strcmp(argv[i], "-mrf_beta") == 0)
            {
                segment_param->MRF_strength=(SegPrecisionTYPE)atof(argv[++i]);
                if(segment_param->MRF_strength==0)
                {
                    segment_param->flag_MRF=0;
                }
            }
            else if(strcmp(argv[i], "-sg_off") == 0)
            {
                segment_param->flag_SG_deli=0;
            }
            else if(strcmp(argv[i], "-max_iter") == 0)
            {
                segment_param->maxIteration=atoi(argv[++i]);
            }
            else if(strcmp(argv[i], "-rf") == 0)
            {
                segment_param->relax_factor=atof(argv[++i]);
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
                fprintf(stderr,"Err:\tParameter %s unknown->\n",argv[i]);
                Usage(argv[0]);
                return 1;
            }
        }

        if(segment_param->flag_PV_model==0)
        {
            segment_param->flag_SG_deli=0;
        }

        if(segment_param->flag_out==0)
        {
            fprintf(stderr,"Err:\tThe output image name has to be defined.\n");
            return 1;
        }

        if(!segment_param->flag_T1)
        {
            fprintf(stderr,"Err:\tThe T1 image name has to be defined.\n");
            Usage(argv[0]);
            return 1;
        }

        // READING T1
        nifti_image * T1=nifti_image_read(segment_param->filename_T1,true);
        if(T1 == NULL)
        {
            fprintf(stderr,"* Error when reading the T1 image: %s\n",segment_param->filename_T1);
            return 1;
        }
        seg_changeDatatype<SegPrecisionTYPE>(T1);


        nifti_image * Mask = nifti_image_read(segment_param->filename_mask,true);
        if(Mask == NULL)
        {
            fprintf(stderr,"* Error when reading the mask image: %s\n",segment_param->filename_mask);
            return 1;
        }
        seg_changeDatatype<SegPrecisionTYPE>(Mask);


        nifti_image ** Priors_temp=new nifti_image * [nonPVNumClass];

        if(segment_param->flag_manual_priors)
        {

            for(int i=0; i<segment_param->numb_classes; i++)
            {
                Priors_temp[i] = nifti_image_read(segment_param->filename_priors[i],true);

                if(Priors_temp[i] == NULL)
                {
                    fprintf(stderr,"* Error when reading the WM Prior image: %s\n",segment_param->filename_priors[i]);
                    return 1;
                }
                seg_changeDatatype<SegPrecisionTYPE>(Priors_temp[i]);
            }
        }
        else
        {
            fprintf(stderr,"* Manual Priors option has to be ON for now...\n");
            return 1;
        }

        if ( segment_param->verbose_level>1 )
            segment_param->Print( cout );

        nifti_image * Priors=nifti_copy_nim_info(T1);
        Priors->dim[0]=4;
        Priors->dim[4]=nonPVNumClass;
        Priors->datatype=DT_FLOAT32;
        Priors->cal_max=1;

        nifti_update_dims_from_array(Priors);
        nifti_datatype_sizes(Priors->datatype,&Priors->nbyper,&Priors->swapsize);
        Priors->data = (void *) calloc(Priors->nvox, sizeof(SegPrecisionTYPE));

        Merge_Priors(Priors,Priors_temp);
        for(int i=0; i<nonPVNumClass; i++)
        {
            nifti_image_free(Priors_temp[i]);
            Priors_temp[i]=NULL;
        }

        nifti_image * Result = LoAd_Segment(T1,Mask,Priors,segment_param);
        Result->scl_inter=0;
        Result->scl_slope=1;
        nifti_image_write(Result);
        nifti_image_free(Result);
        nifti_image_free(Priors);
        nifti_image_free(T1);
        nifti_image_free(Mask);

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
