/**
 * @file seg_CalcTopNCC.cpp
 * @author M. Jorge Cardoso
 * @date 01/01/2014
 *
 * Copyright (c) 2014, University College London. All rights reserved.
 * Centre for Medical Image Computing (CMIC)
 * See the LICENSE.txt file in the nifty_seg root folder
 *
 */

#include "_seg_tools.h"
#include "_seg_common.h"
#include <iostream>
#include <time.h>
#include <stdlib.h>

using namespace std;
#define SegPrecisionTYPE float


void Usage(char *exec)
{
    printf("\nUsage:\t%s -target <filename> -templates <Number of templates> <Template Names> -n <Number of Top Templates> <OPTIONS>\n\n",exec);
    printf("\t* * Options * *\n");
    printf("\t-mask <filename>\tFilename of the ROI mask\n");
#ifdef _GIT_HASH
    printf("\t--version\t\tPrint current source code git hash key and exit\n\t\t\t\t(%s)\n",_GIT_HASH);
#endif
    return;
}

int main(int argc, char **argv)
{


    char * filename_target=NULL;
    char * filename_mask=NULL;

    int numb_templates=1;
    int numb_out=1;
    char ** filename_templates=NULL;

    if(argc==1)
    {
        Usage(argv[0]);
        return 1;
    }
    for(int i=1; i<argc; i++)
    {
        if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
                strcmp(argv[i], "-HELP")==0 || strcmp(argv[i], "-h")==0 ||
                strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0)
        {
            Usage(argv[0]);
            return 0;
        }

        else if(strcmp(argv[i], "-target") == 0)
        {
            filename_target = argv[++i];
        }
        else if(strcmp(argv[i], "-templates") == 0)
        {
            numb_templates = atoi(argv[++i]);
            filename_templates = new char * [numb_templates] ;
            for (int j=0; j<numb_templates; j++)
            {
                filename_templates[j]=argv[++i];
            }
        }
        else if(strcmp(argv[i], "-n") == 0)
        {
            numb_out=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-mask") == 0)
        {
            filename_mask = argv[++i];
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


    nifti_image * Image_Target=nifti_image_read(filename_target,true);
    if(filename_target==NULL)
    {
        fprintf(stderr, "This image can not be read: %s\n", filename_target);
        return 0;
    }
    seg_changeDatatype<float>(Image_Target);

    nifti_image * Mask=NULL;

    if(filename_mask!=NULL)
    {
        Mask = nifti_image_read(filename_mask,true);
        if(Mask->datatype!=DT_BINARY)
        {
            seg_convert2binary(Mask,0.5f);
        }

        if(Mask == NULL)
        {
            fprintf(stderr,"* Error when reading the mask image: %s\n",filename_mask);
            return 1;
        }
    }



    nifti_image * Image_template=NULL;
    float * nccvalues = new float [numb_templates];
    float * nccvaluesold = new float [numb_templates];
    for(int i=0; i<numb_templates; i++)
    {
        Image_template=nifti_image_read(filename_templates[i],true);
        if(Image_template==NULL)
        {
            fprintf(stderr, "This image can not be read: %s\n", filename_templates[i]);
            return 0;
        }
        seg_changeDatatype<float>(Image_template);
        if(Image_Target->nx != Image_template->nx ||
                Image_Target->ny != Image_template->ny ||
                Image_Target->nz != Image_template->nz)
        {
            fprintf(stderr, "This image does not have the correct size: %s\n", filename_templates[i]);
            return 0;
        }

        nccvalues[i]=estimateNCC3D(Image_Target,Image_template,Mask,0);
        nccvaluesold[i]=nccvalues[i];
        nifti_image_free(Image_template);
        Image_template=NULL;
    }

    int * Final_indeces=quickSort_order(nccvalues, numb_templates);

    for(int i=(numb_templates-1); i>(numb_templates-numb_out-1); i--)
    {
        cout<< filename_templates[Final_indeces[i]]<< " ";

    }
    cout << endl;
    flush(cout);
    delete [] Final_indeces;
    nifti_image_free(Image_Target);
    return 0;
}
