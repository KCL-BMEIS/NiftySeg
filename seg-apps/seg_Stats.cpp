/**
 * @file seg_Stats.cpp
 * @author M. Jorge Cardoso
 * @date 01/01/2014
 *
 * Copyright (c) 2014, University College London. All rights reserved.
 * Centre for Medical Image Computing (CMIC)
 * See the LICENSE.txt file in the nifty_seg root folder
 *
 */

#include "_seg_tools.h"
#include <limits>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <new>
#include <exception>
#include <cmath>

using namespace std;
#define SegPrecisionTYPE float

void Usage(char *exec)
{
    printf("\nStat tools:\nUsage:\t%s <in> [constrains] [statistics]\n\n",exec);
    printf("\t* * Constrains and configuration options (optional) * *\n");
    printf("\t-m <mask> \t| Only estimate statistics within the masked area.\n");
    printf("\t-t <float> \t| Only estimate statistics if voxel is larger than <float>.\n");
    printf("\t-p <int> \t| Set output precision (number of digits), by default is 6.\n");
    printf("\n\t  Note: All NaN or Inf are ignored for all stats. \n\t The -m and -t options can be used in conjusction.\n\n");
    printf("\n\t* * Statistics (at least one option is mandatory) * *\n");
    printf("\tRange operations (datatype: all)\n");
    printf("\t-r \t\t| The range <min max> of all voxels.\n");
    printf("\t-R \t\t| The robust range (assuming 2%% outliers on both sides) of all voxels\n");
    printf("\t-p <float> \t| The <float>th percentile of all voxels intensity (float=[0,100])\n");
    printf("\n\tClassical statistics (datatype: all)\n");
    printf("\t-a  \t\t| Average of all voxels \n");
    printf("\t-s  \t\t| Standard deviation of all voxels \n");
    printf("\t-v  \t\t| Volume of all voxels above 0 (<# voxels> * <volume per voxel>)\n");
    printf("\t-vl \t\t| Volume of each integer label (<# voxels per label> * <volume per voxel>)\n");
    printf("\t-vp \t\t| Volume of all probabilistic voxels (sum(<in>) * <volume per voxel>)\n");
    printf("\t-n  \t\t| Count of all voxels above 0 (<# voxels>)\n");
    printf("\t-np \t\t| Sum of all fuzzy voxels (sum(<in>))\n");
    printf("\t-e \t\t| Entropy of all voxels\n");
    printf("\t-ne \t\t| Normalized entropy of all voxels\n");

    printf("\n\tClassical statistics per slice along axis <ax> (ax=1,2,3)\n");
    printf("\t-sa  <ax> \t| Average of all voxels \n");
    printf("\t-ss  <ax> \t| Standard deviation of all voxels \n");
//    printf("\t-sv  <ax> \t\t| Volume of all voxels above 0 (<# voxels> * <volume per voxel>)\n");
//    printf("\t-svl <ax>\t\t| Volume of each integer label (<# voxels per label> * <volume per voxel>)\n");
    printf("\t-svp <ax>\t| Volume of all probabilistic voxels (sum(<in>) * <volume per voxel>)\n");
//    printf("\t-sn  <ax> \t\t| Count of all voxels above 0 (<# voxels>)\n");
//    printf("\t-snp <ax>\t\t| Sum of all fuzzy voxels (sum(<in>))\n");

    printf("\n\tImage similarities (datatype: all)\n");
    printf("\t-ncc <in2>\t| Normalized cross correlation between <in> and <in2>\n");
    printf("\t-nmi <in2>\t| Normalized Mutual Information between <in> and <in2>\n");
    printf("\n\tCoordinates operations (datatype: all)\n");
    printf("\t-x \t\t| Location (i j k x y z) of the smallest value in the image\n");
    printf("\t-X \t\t| Location (i j k x y z) of the largest value in the image\n");
    printf("\t-c \t\t| Location (i j k x y z) of the centre of mass of the object\n");
    printf("\t-B \t\t| Bounding box of all nonzero voxels [ xmin xsize ymin ysize zmin zsize ]\n");
    printf("\n\tHeader info (datatype: all)\n");
    printf("\t-xvox \t\t| Output the number of voxels in the x direction. Replace x with y/z for other directions.\n");
    printf("\t-xdim \t\t|  Output the voxel dimention in the x direction. Replace x with y/z for other directions. \n");
    printf("\n\tLabel attribute operations (datatype: char or uchar)\n");
    printf("\t-Vl <csv> \t| Volume of each integer label <in>. Save to <csv> file.\n");
    printf("\t-Nl <csv> \t| Count of each label <in>. Save to <csv> file.\n");
    printf("\t-al <in2> \t| Average value in <in> for each label in <in2> \n");
    printf("\t-Al <in2> <csv>\t| Average value in <in> for each label in <in2>. Save to <csv> file\n");
    printf("\t-d <in2>\t| Calculate the Dice score between all classes in <in> and <in2>\n");
    printf("\t-D <in2> <csv>\t| Calculate the Dice score between all classes in <in> and <in2>. Save to <csv> file\n");
#ifdef _GIT_HASH
    printf("\t--version\t| Print current source code git hash key and exit\n\t\t\t\t(%s)\n",_GIT_HASH);
#endif
    printf("\n\t* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
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

        char * filenames[100];
        nifti_image * Images[100];
        int numbimg=1;
        float * imgsort=NULL;
        if(argc<3)
        {
            Usage(argv[0]);
            return 1;
        }


        filenames[0] = argv[1];
        Images[0]=nifti_image_read(filenames[0],true);
        if(Images[0]==NULL)
        {
            fprintf(stderr, "This image %s can not be read\n", filenames[0]);
            return 0;
        }


        if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
        {
            seg_changeDatatype<float>(Images[0]);
        }
        float * Img1prt = static_cast<float *>(Images[0]->data);

        cout.precision(6);
        nifti_image * Mask=nifti_copy_nim_info(Images[0]);
        Mask->dim[0]=3;
        Mask->nt=Mask->dim[4]=0;
        Mask->datatype=NIFTI_TYPE_UINT8;
        Mask->cal_max=1;
        Mask->cal_min=0;
        nifti_update_dims_from_array(Mask);
        nifti_datatype_sizes(Mask->datatype,&Mask->nbyper,&Mask->swapsize);
        Mask->data = (void *) calloc(Mask->nvox, sizeof(unsigned char));
        unsigned char * mask = static_cast<unsigned char *>(Mask->data);

        unsigned int maskcount=0;
        for(unsigned int index=0; index<Images[0]->nvox; index++)
        {
            if(isnan(Img1prt[index])==0){
            mask[index]=1;
            maskcount++;
            }
            else{
                mask[index]=0;
            }
        }


        for(int i=2; i<argc; i++)
        {
            if(argc>50){
                fprintf(stderr, "ERROR: too many arguments\n");
                return 0;
            }
            if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
                    strcmp(argv[i], "-HELP")==0 || strcmp(argv[i], "-h")==0 ||
                    strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0)
            {
                Usage(argv[0]);
                return 0;
            }
	    else if(strcmp(argv[i], "-p") == 0 && (i+1)<argc)
            {
	        int pres = atoi(argv[++i]);
		cout.precision(pres);
	    }
            // **************************            ---------          *****************************
            // **************************            Mask Stats         *****************************
            // **************************            ---------          *****************************
            else if(strcmp(argv[i], "-m") == 0 && (i+1)<argc)
            {
                int oldnumbimg=numbimg;
                numbimg=numbimg+1;
                filenames[oldnumbimg] = argv[++i];

                Images[oldnumbimg]=nifti_image_read(filenames[oldnumbimg],true);
                if(Images[oldnumbimg]==NULL)
                {
                    fprintf(stderr, "This image can not be read: %s\n", filenames[oldnumbimg]);
                    return 0;
                }
                if(Images[oldnumbimg]->nvox!=Images[0]->nvox)
                {
                    cout<<"ERROR: The mask is not the same size as the image <in>."<<endl;
                    return 0;
                }
                if(Images[oldnumbimg]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[oldnumbimg]);
                }
                float * maskptr = static_cast<float *>(Images[oldnumbimg]->data);
                maskcount=0;
                for(unsigned int index=0; index<Images[0]->nvox; index++)
                {
                    mask[index]=(maskptr[index]>0)?mask[index]:0;
                    maskcount+=mask[index];
                }
                if(maskcount==0)
                {
                    cout<<"ERROR: Because of the choice of mask, no samples are available for further calculations."<<endl;
                    return 0;
                }
            }
            // **************************            ---------          *****************************
            // **************************         Threshold Stats       *****************************
            // **************************            ---------          *****************************
            else if(strcmp(argv[i], "-t") == 0 && (i+1)<argc)
            {
                float threshold = atof(argv[++i]);
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);

                maskcount=0;
                for(unsigned int index=0; index<Images[0]->nvox; index++)
                {
                    if(Img1prt[index]>threshold)
                    {
                        maskcount++;
                    }
                    else
                    {
                        mask[index]=0;
                    }
                }
                if(maskcount==0)
                {
                    cout<<"ERROR: Because of threshold choice, no samples are available for further calculations."<<endl;
                    return 0;
                }
            }
            // **************************            ---------          *****************************
            // **************************            CALC NCC           *****************************
            // **************************            ---------          *****************************
            else if(strcmp(argv[i], "-ncc") == 0 && (i+1)<argc)
            {
                int oldnumbimg=numbimg;
                numbimg=numbimg+1;
                filenames[oldnumbimg] = argv[++i];
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                for(int j=oldnumbimg; j<numbimg; j++)
                {
                    Images[j]=nifti_image_read(filenames[j],true);
                    if(Images[j]==NULL)
                    {
                        fprintf(stderr, "This image can not be read: %s\n", filenames[j]);
                        return 0;
                    }
                    if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                    {
                        seg_changeDatatype<float>(Images[j]);
                    }
                    if(Images[j]->nx==Images[0]->nx && Images[j]->ny==Images[0]->ny && Images[j]->nz==Images[0]->nz && Images[j]->nt==Images[0]->nt){
                    cout<< estimateNCC3D(Images[0],Images[j],Mask,0)<<endl;
                    }
                    else{
                        fprintf(stderr, "The images %s and %s have different sizes\n", filenames[0], filenames[j]);
                    }
                }
            }

            // **************************            ---------          *****************************
            // **************************            CALC NMI          *****************************
            // **************************            ---------          *****************************
            else if(strcmp(argv[i], "-nmi") == 0 && (i+1)<argc)
            {
                int oldnumbimg=numbimg;
                numbimg=numbimg+1;
                filenames[oldnumbimg] = argv[++i];
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                for(int j=oldnumbimg; j<numbimg; j++)
                {
                    Images[j]=nifti_image_read(filenames[j],true);
                    if(Images[j]==NULL)
                    {
                        fprintf(stderr, "This image can not be read: %s\n", filenames[j]);
                        return 0;
                    }
                    if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                    {
                        seg_changeDatatype<float>(Images[j]);
                    }
                    if(Images[j]->nx==Images[0]->nx && Images[j]->ny==Images[0]->ny && Images[j]->nz==Images[0]->nz && Images[j]->nt==Images[0]->nt){
                        cout<< seg_getNMIValue(Images[0],Images[j],mask)<<endl;
                    }
                    else{
                        fprintf(stderr, "The images %s and %s have different sizes\n", filenames[0], filenames[j]);
                    }
                }
            }
            // **************************            ---------          *****************************
            // **************************            CALC DICE          *****************************
            // **************************            ---------          *****************************
            else if(strcmp(argv[i], "-d") == 0 && (i+1)<argc)
            {
                int oldnumbimg=numbimg;
                numbimg=numbimg+1;
                filenames[oldnumbimg] = argv[++i];
                if(Images[0]->datatype!=NIFTI_TYPE_UINT8)
                {
                    seg_changeDatatype<unsigned char>(Images[0]);
                }
                for(int j=oldnumbimg; j<numbimg; j++)
                {
                    Images[j]=nifti_image_read(filenames[j],true);
                    if(Images[j]==NULL)
                    {
                        fprintf(stderr, "This image can not be read: %s\n", filenames[j]);
                        return 0;
                    }
                    if(Images[j]->datatype!=NIFTI_TYPE_UINT8)
                    {
                        seg_changeDatatype<unsigned char>(Images[j]);
                    }
                }
                int  CountIMG1[1000]= {0};
                unsigned char * Img1prt = static_cast<unsigned char *>(Images[0]->data);
                int  CountIMG2[1000]= {0};
                unsigned char * Img2prt = static_cast<unsigned char *>(Images[oldnumbimg]->data);
                int  CountINTERSECT[1000]= {0};
                int maxclass=0;
                for(unsigned int index=0; index<Images[0]->nvox; index++)
                {
                    CountIMG1[(int)(Img1prt[index])]++;
                    maxclass=(int)(Img1prt[index])>maxclass?(int)(Img1prt[index]):maxclass;
                    CountIMG2[(int)(Img2prt[index])]++;
                    maxclass=(int)(Img2prt[index])>maxclass?(int)(Img2prt[index]):maxclass;
                    if((int)(Img1prt[index])==(int)(Img2prt[index]))
                    {
                        CountINTERSECT[(int)(Img1prt[index])]++;
                    }
                }
                float meanDice=0;
                int meanDiceCount=0;
                for(int curtclass=0; curtclass<=maxclass; curtclass++)
                {
                    float curval=(float)2.0*(float)CountINTERSECT[curtclass]/((float)CountIMG1[curtclass]+(float)CountIMG2[curtclass]);
                    cout<< "Label["<<curtclass<<"] = "<< curval<<endl;
                    if(curval==curval)
                    {
                        meanDice+=(float)2.0*(float)CountINTERSECT[curtclass]/((float)CountIMG1[curtclass]+(float)CountIMG2[curtclass]);
                        meanDiceCount++;
                    }

                }
                if(maxclass>1)
                {
                    cout<< "Mean Dice = "<< meanDice/meanDiceCount<<"\n"<<endl;
                    flush(cout);
                }
            }
            // **************************            ---------          *****************************
            // **************************            CALC DICE          *****************************
            // **************************            ---------          *****************************
            else if(strcmp(argv[i], "-di") == 0 && (i+1)<argc)
            {
                int oldnumbimg=numbimg;
                numbimg=numbimg+1;
                float tmpthresh = atof(argv[++i]);
                filenames[oldnumbimg] = argv[++i];
                if(Images[0]->datatype!=NIFTI_TYPE_UINT8)
                {
                    seg_changeDatatype<unsigned char>(Images[0]);
                }
                for(int j=oldnumbimg; j<numbimg; j++)
                {
                    Images[j]=nifti_image_read(filenames[j],true);
                    if(Images[j]==NULL)
                    {
                        fprintf(stderr, "This image can not be read: %s\n", filenames[j]);
                        return 0;
                    }
                    if(Images[j]->datatype!=NIFTI_TYPE_UINT8)
                    {
                        seg_changeDatatype<unsigned char>(Images[j]);
                    }
                }
                int  CountIMG1[1000]= {0};
                unsigned char * Img1prt = static_cast<unsigned char *>(Images[0]->data);
                int  CountIMG2[1000]= {0};
                unsigned char * Img2prt = static_cast<unsigned char *>(Images[oldnumbimg]->data);
                int  CountINTERSECT[1000]= {0};
                int maxclass=0;
                for(unsigned int index=0; index<Images[0]->nvox; index++)
                {
                    CountIMG1[(int)(Img1prt[index])]++;
                    maxclass=(int)(Img1prt[index])>maxclass?(int)(Img1prt[index]):maxclass;
                    CountIMG2[(int)(Img2prt[index])]++;
                    maxclass=(int)(Img2prt[index])>maxclass?(int)(Img2prt[index]):maxclass;
                    if((int)(Img1prt[index])==(int)(Img2prt[index]))
                    {
                        CountINTERSECT[(int)(Img1prt[index])]++;
                    }
                }
                float meanDice=0;
                int meanDiceCount=0;
                for(int curtclass=0; curtclass<=maxclass; curtclass++)
                {
                    float curval=(float)2.0*(float)CountINTERSECT[curtclass]/((float)CountIMG1[curtclass]+(float)CountIMG2[curtclass]);
                    cout<< "Label["<<curtclass<<"] = "<< curval<<endl;
                    if(curval==curval && curval>tmpthresh)
                    {
                        meanDice+=(float)2.0*(float)CountINTERSECT[curtclass]/((float)CountIMG1[curtclass]+(float)CountIMG2[curtclass]);
                        meanDiceCount++;
                    }

                }
                if(maxclass>1)
                {
                    cout<< "Mean Dice = "<< meanDice/meanDiceCount<<"\n"<<endl;
                    flush(cout);
                }
            }
            // **************************            ---------          *****************************
            // **************************            CSV  DICE          *****************************
            // **************************            ---------          *****************************
            else if(strcmp(argv[i], "-D") == 0 && (i+2)<argc)
            {

                int oldnumbimg=numbimg;
                numbimg=numbimg+1;
                filenames[oldnumbimg] = argv[++i];
                filenames[oldnumbimg+1] = argv[++i];

                if(Images[0]->datatype!=NIFTI_TYPE_UINT8)
                {
                    seg_changeDatatype<unsigned char>(Images[0]);
                }

                cout<<filenames[oldnumbimg]<<' '<<filenames[oldnumbimg+1]<<endl;
                Images[oldnumbimg]=nifti_image_read(filenames[oldnumbimg+1],true);
                if(Images[oldnumbimg]==NULL)
                {
                    fprintf(stderr, "This image can not be read: %s\n", filenames[oldnumbimg]);
                    return 0;
                }
                if(Images[oldnumbimg]->datatype!=NIFTI_TYPE_UINT8)
                {
                    seg_changeDatatype<unsigned char>(Images[oldnumbimg]);
                }

                int  CountIMG1[1000]= {0};
                unsigned char * Img1prt = static_cast<unsigned char *>(Images[0]->data);
                int  CountIMG2[1000]= {0};
                unsigned char * Img2prt = static_cast<unsigned char *>(Images[oldnumbimg]->data);
                int  CountINTERSECT[1000]= {0};
                int maxclass=0;
                for(unsigned int index=0; index<Images[0]->nvox; index++)
                {
                    if(mask[index])
                    {
                        CountIMG1[(int)(Img1prt[index])]++;
                        maxclass=(int)(Img1prt[index])>maxclass?(int)(Img1prt[index]):maxclass;
                        CountIMG2[(int)(Img2prt[index])]++;
                        maxclass=(int)(Img2prt[index])>maxclass?(int)(Img2prt[index]):maxclass;
                        if((int)(Img1prt[index])==(int)(Img2prt[index]))
                        {
                            CountINTERSECT[(int)(Img1prt[index])]++;
                        }
                    }
                }
                filenames[oldnumbimg+1] = argv[++i];
                ofstream myfile;
                myfile.open(filenames[oldnumbimg]);

                flush(cout);
                for(int curtclass=0; curtclass<=maxclass; curtclass++)
                {
                    myfile<< (float)2.0*(float)CountINTERSECT[curtclass]/((float)CountIMG1[curtclass]+(float)CountIMG2[curtclass]);

                    if(curtclass!=maxclass)
                    {
                        myfile<<",";
                    }
                }
                myfile.close();

            }
            // **************************            ---------          *****************************
            // **************************            Lab Mean          *****************************
            // **************************            ---------          *****************************
            else if(strcmp(argv[i], "-al") == 0 && (i+1)<argc)
            {
                int oldnumbimg=numbimg;
                numbimg=numbimg+1;
                filenames[oldnumbimg] = argv[++i];
                for(int j=oldnumbimg; j<numbimg; j++)
                {
                    Images[j]=nifti_image_read(filenames[j],true);
                    if(Images[j]==NULL)
                    {
                        fprintf(stderr, "This image can not be read: %s\n", filenames[j]);
                        return 0;
                    }
                    if(Images[j]->datatype!=NIFTI_TYPE_UINT8)
                    {
                        seg_changeDatatype<unsigned char>(Images[j]);
                    }
                }

                float * Img1prt = static_cast<float *>(Images[0]->data);
                double  Count1[1000]= {0};
                double  Count2[1000]= {0};
                for(unsigned int index=0; index<1000; index++){
                   Count1[i]=0;
                   Count2[i]=0;
                }
                unsigned char * Img2prt = static_cast<unsigned char *>(Images[oldnumbimg]->data);
                int maxclass=0;
                for(unsigned int index=0; index<Images[0]->nvox; index++)
                {
                    Count1[(int)(Img2prt[index])]+=float(Img1prt[index]);
                    Count2[(int)(Img2prt[index])]+=1;
                    maxclass=(int)(Img2prt[index])>maxclass?(int)(Img2prt[index]):maxclass;
                }
                for(int curtclass=0; curtclass<=maxclass; curtclass++)
                {
                    cout<< "Label["<<curtclass<<"] = "<<Count1[curtclass]/Count2[curtclass]<<endl;
                }
                if(maxclass>1)
                {
                    flush(cout);
                }
            }
            // **************************            ---------          *****************************
            // **************************           CSV  Lab Mean       *****************************
            // **************************            ---------          *****************************
            else if(strcmp(argv[i], "-Al") == 0 && (i+2)<argc)
            {

                int oldnumbimg=numbimg;
                numbimg=numbimg+1;
                filenames[oldnumbimg] = argv[++i];
                filenames[oldnumbimg+1] = argv[++i];

                for(int j=oldnumbimg; j<numbimg; j++)
                {
                    Images[j]=nifti_image_read(filenames[j],true);
                    if(Images[j]==NULL)
                    {
                        fprintf(stderr, "This image can not be read: %s\n", filenames[j]);
                        return 0;
                    }
                    if(Images[j]->datatype!=NIFTI_TYPE_UINT8)
                    {
                        seg_changeDatatype<unsigned char>(Images[j]);
                    }
                }

                float * Img1prt = static_cast<float *>(Images[0]->data);
                float Count1[1000]= {0};
                float  Count2[1000]= {0};
                for(unsigned int index=0; index<1000; index++){
                   Count1[i]=0;
                   Count2[i]=0;
                }
                unsigned char * Img2prt = static_cast<unsigned char *>(Images[oldnumbimg]->data);
                int maxclass=0;
                for(unsigned int index=0; index<Images[0]->nvox; index++)
                {
                    Count1[(int)(Img2prt[index])]+=(Img1prt[index]);
                    Count2[(int)(Img2prt[index])]+=1;
                    maxclass=(int)(Img2prt[index])>maxclass?(int)(Img2prt[index]):maxclass;
                }

                ofstream myfile;
                myfile.open(filenames[oldnumbimg+1]);

                flush(cout);
                for(int curtclass=0; curtclass<=maxclass; curtclass++)
                {
                    myfile<<(float)(Count1[curtclass]/Count2[curtclass]);

                    if(curtclass!=maxclass)
                    {
                        myfile<<",";
                    }
                }
                myfile.close();

            }
            // **************************            ---------          *****************************
            // **************************            Fuzzy Vol          *****************************
            // **************************            ---------          *****************************

            else if(strcmp(argv[i], "-vp") == 0 && (i)<argc)
            {
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);
                double calcvol=0;
                for(unsigned int index=0; index<Images[0]->nvox; index++)
                {
                    if(mask[index] && isnan(Img1prt[index])==0)
                    {
                        calcvol += Img1prt[index];
                    }
                }

                cout << (double)(calcvol)*(double)(Images[0]->dx)*(double)(Images[0]->dy)*(double)(Images[0]->dz)<<endl;
                flush(cout);
            }
            // **************************            ---------          *****************************
            // **************************            Bin   Vol          *****************************
            // **************************            ---------          *****************************
            else if(strcmp(argv[i], "-v") == 0 && (i)<argc)
            {
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);
                double calcvol=0;
                for(unsigned int index=0; index<Images[0]->nvox; index++)
                {
                    if(mask[index])
                    {
                        calcvol += Img1prt[index]>0;
                    }
                }
                cout <<(float)(calcvol)*(double)(Images[0]->dx)*(double)(Images[0]->dy)*(double)(Images[0]->dz)<<endl;
                flush(cout);
            }

            // **************************            ---------          *****************************
            // **************************            Bin l Vol          *****************************
            // **************************            ---------          *****************************
            else if(strcmp(argv[i], "-vl") == 0 && (i)<argc)
            {
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);
                int  CountIMG1[1024]= {0};
                int maxclass=0;
                int curlab=0;
                for(unsigned int index=0; index<Images[0]->nvox; index++)
                {
                    if((round(Img1prt[index]))>1024  || (round(Img1prt[index]))<0)
                    {
                        cout<<"Too many labels... the code only handles up to 1024 labels"<<endl;
                        exit(0);
                    }
                    else
                    {
                        curlab=(int)(round(Img1prt[index]));
                    }

                    CountIMG1[curlab]++;
                    maxclass=curlab>maxclass?curlab:maxclass;
                }
                cout<< "Volumes in mm3:"<<endl;
                for(int curtclass=0; curtclass<=maxclass; curtclass++)
                {
                    cout<< "Label["<<curtclass<<"] = "<< (double)CountIMG1[curtclass]*(double)(Images[0]->dx)*(double)(Images[0]->dy)*(double)(Images[0]->dz) <<endl;
                }

                flush(cout);
            }
            // **************************            ---------          *****************************
            // **************************            Bin l Vol          *****************************
            // **************************            ---------          *****************************
            else if(strcmp(argv[i], "-Vl") == 0 && (i+1)<argc)
            {
                char * filenameCSVoutput = argv[++i];
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);
                int  CountIMG1[1024]= {0};
                int maxclass=0;
                int curlab=0;
                for(unsigned int index=0; index<Images[0]->nvox; index++)
                {
                    if((round(Img1prt[index]))>1024  || (round(Img1prt[index]))<0)
                    {
                        cout<<"Too many labels... the code only handles up to 1024 labels"<<endl;
                        exit(0);
                    }
                    else
                    {
                        curlab=(int)(round(Img1prt[index]));
                    }

                    CountIMG1[curlab]++;
                    maxclass=curlab>maxclass?curlab:maxclass;
                }
                ofstream myfile;
                myfile.open(filenameCSVoutput);

                flush(cout);
                for(int curtclass=0; curtclass<=maxclass; curtclass++)
                {
                    myfile<< (double)CountIMG1[curtclass]*(double)(Images[0]->dx)*(double)(Images[0]->dy)*(double)(Images[0]->dz);

                    if(curtclass!=maxclass)
                    {
                        myfile<<",";
                    }
                }
                myfile.close();

                flush(cout);
            }

            // **************************            ---------          *****************************
            // **************************            Bin l count          *****************************
            // **************************            ---------          *****************************
            else if(strcmp(argv[i], "-Nl") == 0 && (i+1)<argc)
            {
                char * filenameCSVoutput = argv[++i];
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);
                int  CountIMG1[1024]= {0};
                int maxclass=0;
                int curlab=0;
                for(unsigned int index=0; index<Images[0]->nvox; index++)
                {
                    if((round(Img1prt[index]))>1024  || (round(Img1prt[index]))<0)
                    {
                        cout<<"Too many labels... the code only handles up to 1024 labels"<<endl;
                        exit(0);
                    }
                    else
                    {
                        curlab=(int)(round(Img1prt[index]));
                    }

                    CountIMG1[curlab]++;
                    maxclass=curlab>maxclass?curlab:maxclass;
                }
                ofstream myfile;
                myfile.open(filenameCSVoutput);

                flush(cout);
                for(int curtclass=0; curtclass<=maxclass; curtclass++)
                {
                    myfile<< (double)CountIMG1[curtclass];

                    if(curtclass!=maxclass)
                    {
                        myfile<<",";
                    }
                }
                myfile.close();

                flush(cout);
            }

            // **************************            ---------          *****************************
            // **************************            Bounding Box       *****************************
            // **************************            ---------          *****************************
            else if(strcmp(argv[i], "-B")==0 && (i)<argc)
            {
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);
                int nx=Images[0]->nx;
                int ny=Images[0]->ny;
                int nz=Images[0]->nz;
                int locX[2]= {nx,0};
                int locY[2]= {ny,0};
                int locZ[2]= {nz,0};
                int index=0;
                for(int Zindex=0; Zindex<Images[0]->nz; Zindex++)
                {
                    for(int Yindex=0; Yindex<Images[0]->ny; Yindex++)
                    {
                        for(int Xindex=0; Xindex<Images[0]->nx; Xindex++)
                        {
                            if(mask[index] && Img1prt[index]>0)
                            {
                                if(locX[0]>Xindex)
                                    locX[0]=Xindex;
                                if(locX[1]<Xindex)
                                    locX[1]=Xindex;
                                if(locY[0]>Yindex)
                                    locY[0]=Yindex;
                                if(locY[1]<Yindex)
                                    locY[1]=Yindex;
                                if(locZ[0]>Zindex)
                                    locZ[0]=Zindex;
                                if(locZ[1]<Zindex)
                                    locZ[1]=Zindex;
                            }
                            index++;
                        }
                    }
                }
                //for(int j=0; j<2; j++){
                //   cout <<"("<<locX[j]<<","<<locY[j]<<","<<locZ[j]<<") ";
                //}
                //cout<<endl;
                cout <<locX[0]<<" "<<locX[1]-locX[0]<<" "<<locY[0]<<" "<<locY[1]-locY[0]<<" "<<locZ[0]<<" "<<locZ[1]-locZ[0]<<endl;
                flush(cout);
            }
            // **************************            ---------          *****************************
            // **************************           Centre gravity      *****************************
            // **************************            ---------          *****************************

            else if(strcmp(argv[i], "-c")==0 && (i)<argc)
            {
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float *Img1prt = static_cast<float *>(Images[0]->data);
                float intensitySum=0;
                float locVox[3]={0,0,0};
                int index=0;
                for(int Zindex=0; Zindex<Images[0]->nz; Zindex++)
                {
                    for(int Yindex=0; Yindex<Images[0]->ny; Yindex++)
                    {
                        for(int Xindex=0; Xindex<Images[0]->nx; Xindex++)
                        {
                           float intensity=Img1prt[index];
                            if(mask[index] && intensity>0)
                            {
                                locVox[0]+=Xindex*intensity;
                                locVox[1]+=Yindex*intensity;
                                locVox[2]+=Zindex*intensity;
                                intensitySum+=intensity;
                            }
                            index++;
                        }
                    }
                }
                locVox[0]/=intensitySum;
                locVox[1]/=intensitySum;
                locVox[2]/=intensitySum;
                // Compute the coordinate in mm
                float locMil[3]={0,0,0};
                if(Images[0]->sform_code!=0)
                   seg_mat44_mul(&Images[0]->sto_xyz,locVox,locMil);
                else
                   seg_mat44_mul(&Images[0]->qto_xyz,locVox,locMil);
                // Display the coordinates
                cout <<locVox[0]<<" "<<locVox[1]<<" "<<locVox[2]<<" "<<locMil[0]<<" "<<locMil[1]<<" "<<locMil[2]<<endl;
                flush(cout);
            }

            // **************************            ---------          *****************************
            // **************************            Max location       *****************************
            // **************************            ---------          *****************************
            else if(strcmp(argv[i], "-X")==0 && (i)<argc)
            {
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);
                float maxval=-1e32;
                float locVox[3]={0,0,0};
                int index=0;
                for(int Zindex=0; Zindex<Images[0]->nz; Zindex++)
                {
                    for(int Yindex=0; Yindex<Images[0]->ny; Yindex++)
                    {
                        for(int Xindex=0; Xindex<Images[0]->nx; Xindex++)
                        {
                            if(mask[index] && Img1prt[index]>maxval)
                            {
                                maxval=Img1prt[index];
                                locVox[0]=Xindex;
                                locVox[1]=Yindex;
                                locVox[2]=Zindex;
                            }
                            index++;
                        }
                    }
                }
                // Compute the coordinate in mm
                float locMil[3]={0,0,0};
                if(Images[0]->sform_code!=0)
                   seg_mat44_mul(&Images[0]->sto_xyz,locVox,locMil);
                else
                   seg_mat44_mul(&Images[0]->qto_xyz,locVox,locMil);
                // Display the coordinates
                cout <<locVox[0]<<" "<<locVox[1]<<" "<<locVox[2]<<" "<<locMil[0]<<" "<<locMil[1]<<" "<<locMil[2]<<endl;
                flush(cout);
            }

            // **************************            ---------          *****************************
            // **************************            Min location       *****************************
            // **************************            ---------          *****************************
            else if(strcmp(argv[i], "-x")==0 && (i)<argc)
            {
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);
                float minval=1e32;
                float locVox[3]={0,0,0};
                int index=0;
                for(int Zindex=0; Zindex<Images[0]->nz; Zindex++)
                {
                    for(int Yindex=0; Yindex<Images[0]->ny; Yindex++)
                    {
                        for(int Xindex=0; Xindex<Images[0]->nx; Xindex++)
                        {
                            if(mask[index] && Img1prt[index]<minval)
                            {
                                minval=Img1prt[index];
                                locVox[0]=Xindex;
                                locVox[1]=Yindex;
                                locVox[2]=Zindex;
                            }
                            index++;
                        }
                    }
                }
                // Compute the coordinate in mm
                float locMil[3]={0,0,0};
                if(Images[0]->sform_code!=0)
                   seg_mat44_mul(&Images[0]->sto_xyz,locVox,locMil);
                else
                   seg_mat44_mul(&Images[0]->qto_xyz,locVox,locMil);
                // Display the coordinates
                cout <<locVox[0]<<" "<<locVox[1]<<" "<<locVox[2]<<" "<<locMil[0]<<" "<<locMil[1]<<" "<<locMil[2]<<endl;
                flush(cout);
            }

            // **************************            ---------          *****************************
            // **************************              Average          *****************************
            // **************************            ---------          *****************************

            else if(strcmp(argv[i], "-a") == 0 && (i)<argc)
            {
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);
                double calcvol=0;
                double calcvolcount=0;
                for(unsigned int index=0; index<Images[0]->nvox; index++)
                {
                    if(mask[index])
                    {
                        calcvol += (double)Img1prt[index];
                        calcvolcount+=1;
                    }
                }

                cout << (double)(calcvol)/(double)(calcvolcount)<<endl;
                flush(cout);
            }

            // **************************            ---------          *****************************
            // **************************           Average Lab         *****************************
            // **************************            ---------          *****************************

            else if(strcmp(argv[i], "-a") == 0 && (i)<argc)
            {
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);
                float calcvol=0;
                float calcvolcount=0;
                for(unsigned int index=0; index<Images[0]->nvox; index++)
                {
                    if(mask[index])
                    {
                        calcvol += Img1prt[index];
                        calcvolcount+=1;
                    }
                }

                cout << (double)(calcvol)/(double)(calcvolcount)<<endl;
                flush(cout);
            }

            // **************************            ---------          *****************************
            // **************************             Max/min           *****************************
            // **************************            ---------          *****************************

            else if(strcmp(argv[i], "-r") == 0 && (i)<argc)
            {
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);
                float maxval=-1e-32;
                float minval=1e32;
                for(unsigned int index=0; index<Images[0]->nvox; index++)
                {
                    if(mask[index])
                    {
                        maxval=Img1prt[index]>maxval?Img1prt[index]:maxval;
                        minval=Img1prt[index]<minval?Img1prt[index]:minval;
                    }
                }

                cout << minval<<" "<<maxval<<endl;
                flush(cout);
            }

            // **************************            ---------          *****************************
            // **************************          Robust min/max       *****************************
            // **************************            ---------          *****************************

            else if(strcmp(argv[i], "-R") == 0 && (i)<argc)
            {
                float outlier=0.02f;
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);
                if(imgsort==NULL)
                {
                    imgsort=new float [maskcount];
                    int curindex=0;
                    for(unsigned int index=0; index<Images[0]->nvox; index++)
                    {
                        if(mask[index])
                        {
                            imgsort[curindex]=Img1prt[index];
                            curindex++;
                        }
                    }
                    HeapSort(imgsort,maskcount-1);
                }
                cout << imgsort[(int)(round(outlier*(maskcount-1)))]<<" "<<imgsort[(int)(round((1-outlier)*(maskcount-1)))]<<endl;
                flush(cout);
            }

            // **************************            ---------          *****************************
            // **************************          Percentile XX%       *****************************
            // **************************            ---------          *****************************

            else if(strcmp(argv[i], "-p") == 0 && (i)<argc)
            {
                string parser=argv[i+1];
                if(strtod(parser.c_str(),NULL)==0 )
                {
                    cout<<"ERROR: The <float> range in option -P is not a number or is not within the range."<<endl;
                    return 0;
                }
                float percentile = atof(argv[++i])/100.0f;
                percentile=percentile>1?1:percentile;
                percentile=percentile<0?0:percentile;
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);
                if(imgsort==NULL)
                {
                    imgsort=new float [maskcount];
                    int curindex=0;
                    for(unsigned int index=0; index<Images[0]->nvox; index++)
                    {
                        if(mask[index])
                        {
                            imgsort[curindex]=Img1prt[index];
                            curindex++;
                        }
                    }
                    HeapSort(imgsort,maskcount-1);
                }
                cout << imgsort[(int)(floor((double)(percentile)*(double)(maskcount-1)))]<<endl;

                flush(cout);
            }

            // **************************            ---------          *****************************
            // **************************               std             *****************************
            // **************************            ---------          *****************************

            else if(strcmp(argv[i], "-s") == 0 && (i)<argc)
            {
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);
                float calc=0;
                float calccount=0;
                for(unsigned int index=0; index<Images[0]->nvox; index++)
                {
                    if(mask[index])
                    {
                        calc += Img1prt[index];
                        calccount+=1;
                    }
                }
                double mean=(double)(calc)/(double)(calccount);
                calc=0;
                calccount=0;
                for(unsigned int index=0; index<Images[0]->nvox; index++)
                {
                    if(mask[index])
                    {
                        calc += powf(mean-Img1prt[index],2);
                        calccount+=1;
                    }
                }

                cout <<sqrt((double)(calc)/(double)(calccount))<<endl;
                flush(cout);
            }

            // **************************            ---------          *****************************
            // **************************            Fuzzy Numb          *****************************
            // **************************            ---------          *****************************

            else if(strcmp(argv[i], "-np") == 0 && (i)<argc)
            {
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);
                float calcvol=0;
                for(unsigned int index=0; index<Images[0]->nvox; index++)
                {
                    if(mask[index])
                    {
                        calcvol += Img1prt[index];
                    }
                }

                cout << (double)(calcvol)<<endl;
                flush(cout);
            }

            // **************************            ---------          *****************************
            // **************************            Fuzzy Numb          *****************************
            // **************************            ---------          *****************************

            else if(strcmp(argv[i], "-Np") == 0 && (i)<argc)
            {
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);
                float calcvol=0;
                for(unsigned int index=0; index<Images[0]->nvox; index++)
                {
                    if(mask[index])
                    {
                        calcvol += Img1prt[index];
                    }
                }

                cout << (double)(calcvol)<<endl;
                flush(cout);
            }
            // **************************            ---------          *****************************
            // **************************            Bin   Numb          *****************************
            // **************************            ---------          *****************************
            else if(strcmp(argv[i], "-n") == 0 && (i)<argc)
            {
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);
                float calcvol=0;
                for(unsigned int index=0; index<Images[0]->nvox; index++)
                {
                    if(mask[index])
                    {
                        calcvol += Img1prt[index]>0;
                    }
                }

                cout <<(double)(calcvol)<<endl;
                flush(cout);
            }

            // **************************            ---------          *****************************
            // **************************            Bin l count          *****************************
            // **************************            ---------          *****************************
            else if(strcmp(argv[i], "-nl") == 0 && (i)<argc)
            {
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);
                int  CountIMG1[1024]= {0};
                int maxclass=0;
                int curlab=0;
                for(unsigned int index=0; index<Images[0]->nvox; index++)
                {
                    if((round(Img1prt[index]))>1024  || (round(Img1prt[index]))<0)
                    {
                        cout<<"Too many labels... the code only handles up to 1024 labels"<<endl;
                        exit(0);
                    }
                    else
                    {
                        curlab=(int)(round(Img1prt[index]));
                    }

                    CountIMG1[curlab]++;
                    maxclass=curlab>maxclass?curlab:maxclass;
                }
                cout<< "Label counts:"<<endl;
                for(int curtclass=0; curtclass<=maxclass; curtclass++)
                {
                    cout<< "Label["<<curtclass<<"] = "<< (double)CountIMG1[curtclass]<<endl;
                }

                flush(cout);
            }
	    // **************************            ---------          *****************************
            // ********************        Entropy and normalized entropy   *****************************
            // **************************            ---------          *****************************
	    else if((strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "-ne") == 0) && (i)<argc)
            {
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);
                float max=std::numeric_limits<float>::min();
                for(unsigned int index=0; index<Images[0]->nvox; index++)
                {
                    if(mask[index])
                    {
                        if(max<Img1prt[index]) {
                            max=Img1prt[index];
                        }
                    }
                }
                float ent=0;
                float probability=0;
                int count=0;
                for(unsigned int index=0; index<Images[0]->nvox; index++)
                {
                    if(mask[index])
                    {
                        probability=Img1prt[index]/max;
                        if(probability>0) {
                            ent += probability*log(probability);
                        }
                        count++;
                    }
                }
                ent=-ent;
                if(strcmp(argv[i], "-ne") == 0) ent=ent/(float)count;
                cout <<(double)(ent)<<endl;
                flush(cout);
            }
            // **************************            ---------          *****************************
            // **************************             PerSlice Average          *****************************
            // **************************            ---------          *****************************

            else if(strcmp(argv[i], "-sa") == 0 && (i+1)<argc)
            {
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);

                unsigned int direction = atoi(argv[++i]);
                unsigned int dim1=0;
                unsigned int dim2=0;
                unsigned int dirmain=0;
                if(direction>0 && direction<=3){
                    if(direction==1){
                        dirmain=Images[0]->nx;
                        dim1=Images[0]->ny;
                        dim2=Images[0]->nz;
                    }
                    else if(direction==2){
                        dim1=Images[0]->nx;
                        dirmain=Images[0]->ny;
                        dim2=Images[0]->nz;
                      }
                    else if(direction==3){
                        dim1=Images[0]->nx;
                        dim2=Images[0]->ny;
                        dirmain=Images[0]->nz;
                      }
                }
                else{
                   cout<<"Error: Direction " <<direction<<" unknown"<<endl;
                }

                unsigned int index=0;
                double calcvol=0;
                double calcvolcount=0;
                for(unsigned int dirmain_ind=0; dirmain_ind<dirmain; dirmain_ind++)
                {
                    calcvol=calcvolcount=0;
                    for(unsigned int dim1_ind=0; dim1_ind<dim1; dim1_ind++){
                        for(unsigned int dim2_ind=0; dim2_ind<dim2; dim2_ind++){
                            if(direction==1){
                                index=dirmain_ind+dim1_ind*Images[0]->nx+dim2_ind*Images[0]->nx*Images[0]->ny;
                            }
                            else if(direction==2){
                                index=dim1_ind+dirmain_ind*Images[0]->nx+dim2_ind*Images[0]->nx*Images[0]->ny;
                              }
                            else if(direction==3){
                                index=dim1_ind+dim2_ind*Images[0]->nx+dirmain_ind*Images[0]->nx*Images[0]->ny;
                              }
                            else{
                                cout<<"Error: Direction is problematic"<<endl;
                            }
                            if(mask[index]){
                                calcvol += (double)Img1prt[index];
                                calcvolcount+=1;
                            }
                        }
                    }
                    cout << (double)(calcvol)/(double)(calcvolcount)<<endl;
                    flush(cout);
                }
            }
            // **************************            ---------          *****************************
            // **************************             PerSlice std          *****************************
            // **************************            ---------          *****************************

            else if(strcmp(argv[i], "-ss") == 0 && (i+1)<argc)
            {
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);

                unsigned int direction = atoi(argv[++i]);
                unsigned int dim1=0;
                unsigned int dim2=0;
                unsigned int dirmain=0;
                if(direction>0 && direction<=3){
                    if(direction==1){
                        dirmain=Images[0]->nx;
                        dim1=Images[0]->ny;
                        dim2=Images[0]->nz;
                    }
                    else if(direction==2){
                        dim1=Images[0]->nx;
                        dirmain=Images[0]->ny;
                        dim2=Images[0]->nz;
                      }
                    else if(direction==3){
                        dim1=Images[0]->nx;
                        dim2=Images[0]->ny;
                        dirmain=Images[0]->nz;
                      }
                }
                else{
                   cout<<"Error: Direction " <<direction<<" unknown"<<endl;
                }

                unsigned int index=0;
                double calcval=0;
                double calcvalcount=0;
                for(unsigned int dirmain_ind=0; dirmain_ind<dirmain; dirmain_ind++)
                {
                    calcval=calcvalcount=0;
                    for(unsigned int dim1_ind=0; dim1_ind<dim1; dim1_ind++){
                        for(unsigned int dim2_ind=0; dim2_ind<dim2; dim2_ind++){
                            if(direction==1){
                                index=dirmain_ind+dim1_ind*Images[0]->nx+dim2_ind*Images[0]->nx*Images[0]->ny;
                            }
                            else if(direction==2){
                                index=dim1_ind+dirmain_ind*Images[0]->nx+dim2_ind*Images[0]->nx*Images[0]->ny;
                              }
                            else if(direction==3){
                                index=dim1_ind+dim2_ind*Images[0]->nx+dirmain_ind*Images[0]->nx*Images[0]->ny;
                              }
                            else{
                                cout<<"Error: Direction is problematic"<<endl;
                            }
                            if(mask[index]){
                                calcval += (double)Img1prt[index];
                                calcvalcount+=1;
                            }
                        }
                    }
                    float mean=(double)(calcval)/(double)(calcvalcount);
                    calcval=calcvalcount=0;
                    for(unsigned int dim1_ind=0; dim1_ind<dim1; dim1_ind++){
                        for(unsigned int dim2_ind=0; dim2_ind<dim2; dim2_ind++){
                            if(direction==1){
                                index=dirmain_ind+dim1_ind*Images[0]->nx+dim2_ind*Images[0]->nx*Images[0]->ny;
                            }
                            else if(direction==2){
                                index=dim1_ind+dirmain_ind*Images[0]->nx+dim2_ind*Images[0]->nx*Images[0]->ny;
                              }
                            else if(direction==3){
                                index=dim1_ind+dim2_ind*Images[0]->nx+dirmain_ind*Images[0]->nx*Images[0]->ny;
                              }
                            else{
                                cout<<"Error: Direction is problematic"<<endl;
                            }
                            if(mask[index]){
                                calcval += (double)powf(mean-Img1prt[index],2);
                                calcvalcount+=1;
                            }
                        }
                    }
                    cout << sqrt((double)(calcval)/(double)(calcvalcount))<<endl;
                    flush(cout);
                }
            }
            // **************************            ---------          *****************************
            // ************************** PerSlice Average, output image *****************************
            // **************************            ---------          *****************************

            else if(strcmp(argv[i], "-sai") == 0 && (i+2)<argc)
            {
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);

                unsigned int direction = atoi(argv[++i]);
                unsigned int dim1=0;
                unsigned int dim2=0;
                unsigned int dirmain=0;

                nifti_image * OutImage=nifti_copy_nim_info(Images[0]);
                OutImage->dim[0]=3;
                OutImage->nt=Mask->dim[4]=0;
                OutImage->datatype=NIFTI_TYPE_FLOAT32;
                OutImage->cal_max=1;
                OutImage->cal_min=0;
                nifti_update_dims_from_array(OutImage);
                nifti_datatype_sizes(OutImage->datatype,&OutImage->nbyper,&OutImage->swapsize);
                OutImage->data = (void *) calloc(OutImage->nvox, sizeof(float));
                float * OutImagePtr = static_cast<float *>(OutImage->data);
                string OutFilename=argv[++i];

                if(direction>0 && direction<=3){
                    if(direction==1){
                        dirmain=Images[0]->nx;
                        dim1=Images[0]->ny;
                        dim2=Images[0]->nz;
                    }
                    else if(direction==2){
                        dim1=Images[0]->nx;
                        dirmain=Images[0]->ny;
                        dim2=Images[0]->nz;
                      }
                    else if(direction==3){
                        dim1=Images[0]->nx;
                        dim2=Images[0]->ny;
                        dirmain=Images[0]->nz;
                      }
                }
                else{
                   cout<<"Error: Direction " <<direction<<" unknown"<<endl;
                }

                unsigned int index=0;
                double calcvol=0;
                double calcvolcount=0;

                for(unsigned int dirmain_ind=0; dirmain_ind<dirmain; dirmain_ind++)
                {
                    calcvol=calcvolcount=0;
                    for(unsigned int dim1_ind=0; dim1_ind<dim1; dim1_ind++){
                        for(unsigned int dim2_ind=0; dim2_ind<dim2; dim2_ind++){
                            if(direction==1){
                                index=dirmain_ind+dim1_ind*Images[0]->nx+dim2_ind*Images[0]->nx*Images[0]->ny;
                            }
                            else if(direction==2){
                                index=dim1_ind+dirmain_ind*Images[0]->nx+dim2_ind*Images[0]->nx*Images[0]->ny;
                              }
                            else if(direction==3){
                                index=dim1_ind+dim2_ind*Images[0]->nx+dirmain_ind*Images[0]->nx*Images[0]->ny;
                              }
                            else{
                                cout<<"Error: Direction is problematic"<<endl;
                            }
                            if(mask[index]){
                                calcvol += (double)Img1prt[index];
                                calcvolcount+=1;
                            }
                        }
                    }
                    calcvol= (double)(calcvol)/(double)(calcvolcount);
                    for(unsigned int dim1_ind=0; dim1_ind<dim1; dim1_ind++){
                        for(unsigned int dim2_ind=0; dim2_ind<dim2; dim2_ind++){
                            if(direction==1){
                                index=dirmain_ind+dim1_ind*Images[0]->nx+dim2_ind*Images[0]->nx*Images[0]->ny;
                            }
                            else if(direction==2){
                                index=dim1_ind+dirmain_ind*Images[0]->nx+dim2_ind*Images[0]->nx*Images[0]->ny;
                              }
                            else if(direction==3){
                                index=dim1_ind+dim2_ind*Images[0]->nx+dirmain_ind*Images[0]->nx*Images[0]->ny;
                              }
                            else{
                                cout<<"Error: Direction is problematic"<<endl;
                            }
                            OutImagePtr[index]=calcvol;

                        }
                    }
                }
                nifti_set_filenames(OutImage, OutFilename.c_str(),0,0);
                nifti_image_write(OutImage);
                nifti_image_free(OutImage);
            }
            // **************************            ---------          *****************************
            // **************************  PerSlice std, output image   *****************************
            // **************************            ---------          *****************************

            else if(strcmp(argv[i], "-ssi") == 0 && (i+2)<argc)
            {
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);

                unsigned int direction = atoi(argv[++i]);
                unsigned int dim1=0;
                unsigned int dim2=0;
                unsigned int dirmain=0;

                nifti_image * OutImage=nifti_copy_nim_info(Images[0]);
                OutImage->dim[0]=3;
                OutImage->nt=Mask->dim[4]=0;
                OutImage->datatype=NIFTI_TYPE_FLOAT32;
                OutImage->cal_max=1;
                OutImage->cal_min=0;
                nifti_update_dims_from_array(OutImage);
                nifti_datatype_sizes(OutImage->datatype,&OutImage->nbyper,&OutImage->swapsize);
                OutImage->data = (void *) calloc(OutImage->nvox, sizeof(float));
                float * OutImagePtr = static_cast<float *>(OutImage->data);
                string OutFilename=argv[++i];

                if(direction>0 && direction<=3){
                    if(direction==1){
                        dirmain=Images[0]->nx;
                        dim1=Images[0]->ny;
                        dim2=Images[0]->nz;
                    }
                    else if(direction==2){
                        dim1=Images[0]->nx;
                        dirmain=Images[0]->ny;
                        dim2=Images[0]->nz;
                      }
                    else if(direction==3){
                        dim1=Images[0]->nx;
                        dim2=Images[0]->ny;
                        dirmain=Images[0]->nz;
                      }
                }
                else{
                   cout<<"Error: Direction " <<direction<<" unknown"<<endl;
                }

                unsigned int index=0;
                double calcval=0;
                double calcvalcount=0;
                for(unsigned int dirmain_ind=0; dirmain_ind<dirmain; dirmain_ind++)
                {
                    calcval=calcvalcount=0;
                    for(unsigned int dim1_ind=0; dim1_ind<dim1; dim1_ind++){
                        for(unsigned int dim2_ind=0; dim2_ind<dim2; dim2_ind++){
                            if(direction==1){
                                index=dirmain_ind+dim1_ind*Images[0]->nx+dim2_ind*Images[0]->nx*Images[0]->ny;
                            }
                            else if(direction==2){
                                index=dim1_ind+dirmain_ind*Images[0]->nx+dim2_ind*Images[0]->nx*Images[0]->ny;
                              }
                            else if(direction==3){
                                index=dim1_ind+dim2_ind*Images[0]->nx+dirmain_ind*Images[0]->nx*Images[0]->ny;
                              }
                            else{
                                cout<<"Error: Direction is problematic"<<endl;
                            }
                            if(mask[index]){
                                calcval += (double)Img1prt[index];
                                calcvalcount+=1;
                            }
                        }
                    }
                    float mean=(double)(calcval)/(double)(calcvalcount);
                    calcval=calcvalcount=0;
                    for(unsigned int dim1_ind=0; dim1_ind<dim1; dim1_ind++){
                        for(unsigned int dim2_ind=0; dim2_ind<dim2; dim2_ind++){
                            if(direction==1){
                                index=dirmain_ind+dim1_ind*Images[0]->nx+dim2_ind*Images[0]->nx*Images[0]->ny;
                            }
                            else if(direction==2){
                                index=dim1_ind+dirmain_ind*Images[0]->nx+dim2_ind*Images[0]->nx*Images[0]->ny;
                              }
                            else if(direction==3){
                                index=dim1_ind+dim2_ind*Images[0]->nx+dirmain_ind*Images[0]->nx*Images[0]->ny;
                              }
                            else{
                                cout<<"Error: Direction is problematic"<<endl;
                            }
                            if(mask[index]){
                                calcval += (double)powf(mean-Img1prt[index],2);
                                calcvalcount+=1;
                            }
                        }
                    }
                    calcval= sqrt((double)(calcval)/(double)(calcvalcount));
                    for(unsigned int dim1_ind=0; dim1_ind<dim1; dim1_ind++){
                        for(unsigned int dim2_ind=0; dim2_ind<dim2; dim2_ind++){
                            if(direction==1){
                                index=dirmain_ind+dim1_ind*Images[0]->nx+dim2_ind*Images[0]->nx*Images[0]->ny;
                            }
                            else if(direction==2){
                                index=dim1_ind+dirmain_ind*Images[0]->nx+dim2_ind*Images[0]->nx*Images[0]->ny;
                              }
                            else if(direction==3){
                                index=dim1_ind+dim2_ind*Images[0]->nx+dirmain_ind*Images[0]->nx*Images[0]->ny;
                              }
                            else{
                                cout<<"Error: Direction is problematic"<<endl;
                            }
                                OutImagePtr[index]=calcval;
                        }
                    }
                }
                nifti_set_filenames(OutImage, OutFilename.c_str(),0,0);
                nifti_image_write(OutImage);
                nifti_image_free(OutImage);
            }
            // **************************            ---------          *****************************
            // **************************             PerSlice volume          *****************************
            // **************************            ---------          *****************************

            else if(strcmp(argv[i], "-svp") == 0 && (i+1)<argc)
            {
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);

                unsigned int direction = atoi(argv[++i]);
                unsigned int dim1=0;
                unsigned int dim2=0;
                unsigned int dirmain=0;
                if(direction>0 && direction<=3){
                    if(direction==1){
                        dirmain=Images[0]->nx;
                        dim1=Images[0]->ny;
                        dim2=Images[0]->nz;
                    }
                    else if(direction==2){
                        dim1=Images[0]->nx;
                        dirmain=Images[0]->ny;
                        dim2=Images[0]->nz;
                      }
                    else if(direction==3){
                        dim1=Images[0]->nx;
                        dim2=Images[0]->ny;
                        dirmain=Images[0]->nz;
                      }
                }
                else{
                   cout<<"Error: Direction " <<direction<<" unknown"<<endl;
                }

                unsigned int index=0;
                double calcvol=0;
                for(unsigned int dirmain_ind=0; dirmain_ind<dirmain; dirmain_ind++)
                {
                    calcvol=0;
                    for(unsigned int dim1_ind=0; dim1_ind<dim1; dim1_ind++){
                        for(unsigned int dim2_ind=0; dim2_ind<dim2; dim2_ind++){
                            if(direction==1){
                                index=dirmain_ind+dim1_ind*Images[0]->nx+dim2_ind*Images[0]->nx*Images[0]->ny;
                            }
                            else if(direction==2){
                                index=dim1_ind+dirmain_ind*Images[0]->nx+dim2_ind*Images[0]->nx*Images[0]->ny;
                              }
                            else if(direction==3){
                                index=dim1_ind+dim2_ind*Images[0]->nx+dirmain_ind*Images[0]->nx*Images[0]->ny;
                              }
                            else{
                                cout<<"Error: Direction is problematic"<<endl;
                            }
                            if(mask[index]){
                                calcvol += (double)Img1prt[index];
                            }
                        }
                    }
                    cout << (double)(calcvol)*(double)(Images[0]->dx)*(double)(Images[0]->dy)*(double)(Images[0]->dz)<<endl;
                    flush(cout);
                }
            }

            // **************************            ---------          *****************************
            // **************************             PerSlice volume, output image          *****************************
            // **************************            ---------          *****************************

            else if(strcmp(argv[i], "-svpi") == 0 && (i+2)<argc)
            {
                nifti_image * OutImage=nifti_copy_nim_info(Images[0]);
                OutImage->dim[0]=3;
                OutImage->nt=Mask->dim[4]=0;
                OutImage->datatype=NIFTI_TYPE_FLOAT32;
                OutImage->cal_max=1;
                OutImage->cal_min=0;
                nifti_update_dims_from_array(OutImage);
                nifti_datatype_sizes(OutImage->datatype,&OutImage->nbyper,&OutImage->swapsize);
                OutImage->data = (void *) calloc(OutImage->nvox, sizeof(float));
                float * OutImagePtr = static_cast<float *>(OutImage->data);

                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);

                unsigned int direction = atoi(argv[++i]);
                unsigned int dim1=0;
                unsigned int dim2=0;
                unsigned int dirmain=0;
                string OutFilename=argv[++i];
                if(direction>0 && direction<=3){
                    if(direction==1){
                        dirmain=Images[0]->nx;
                        dim1=Images[0]->ny;
                        dim2=Images[0]->nz;
                    }
                    else if(direction==2){
                        dim1=Images[0]->nx;
                        dirmain=Images[0]->ny;
                        dim2=Images[0]->nz;
                      }
                    else if(direction==3){
                        dim1=Images[0]->nx;
                        dim2=Images[0]->ny;
                        dirmain=Images[0]->nz;
                      }
                }
                else{
                   cout<<"Error: Direction " <<direction<<" unknown"<<endl;
                }

                unsigned int index=0;
                double calcvol=0;
                for(unsigned int dirmain_ind=0; dirmain_ind<dirmain; dirmain_ind++)
                {
                    calcvol=0;
                    for(unsigned int dim1_ind=0; dim1_ind<dim1; dim1_ind++){
                        for(unsigned int dim2_ind=0; dim2_ind<dim2; dim2_ind++){
                            if(direction==1){
                                index=dirmain_ind+dim1_ind*Images[0]->nx+dim2_ind*Images[0]->nx*Images[0]->ny;
                            }
                            else if(direction==2){
                                index=dim1_ind+dirmain_ind*Images[0]->nx+dim2_ind*Images[0]->nx*Images[0]->ny;
                              }
                            else if(direction==3){
                                index=dim1_ind+dim2_ind*Images[0]->nx+dirmain_ind*Images[0]->nx*Images[0]->ny;
                              }
                            else{
                                cout<<"Error: Direction is problematic"<<endl;
                            }
                            if(mask[index]){
                                calcvol += (double)Img1prt[index];
                            }
                        }
                    }
                    calcvol=(double)(calcvol)*(double)(Images[0]->dx)*(double)(Images[0]->dy)*(double)(Images[0]->dz);
                    for(unsigned int dim1_ind=0; dim1_ind<dim1; dim1_ind++){
                        for(unsigned int dim2_ind=0; dim2_ind<dim2; dim2_ind++){
                            if(direction==1){
                                index=dirmain_ind+dim1_ind*Images[0]->nx+dim2_ind*Images[0]->nx*Images[0]->ny;
                            }
                            else if(direction==2){
                                index=dim1_ind+dirmain_ind*Images[0]->nx+dim2_ind*Images[0]->nx*Images[0]->ny;
                            }
                            else if(direction==3){
                                index=dim1_ind+dim2_ind*Images[0]->nx+dirmain_ind*Images[0]->nx*Images[0]->ny;
                            }
                            else{
                                cout<<"Error: Direction is problematic"<<endl;
                            }
                            OutImagePtr[index]=calcvol;
                        }
                    }
                }
                nifti_set_filenames(OutImage, OutFilename.c_str(),0,0);
                nifti_image_write(OutImage);
                nifti_image_free(OutImage);
            }
            // **************************            ---------          *****************************
            // **************************            Vox dim X/Y/Z          *****************************
            // **************************            ---------          *****************************
            else if(strcmp(argv[i], "-xvox") == 0 && (i)<argc)
            {
                cout <<(double)(Images[0]->dx)<<endl;
                flush(cout);
            }
            else if(strcmp(argv[i], "-yvox") == 0 && (i)<argc)
            {
                cout <<(double)(Images[0]->dy)<<endl;
                flush(cout);
            }
            else if(strcmp(argv[i], "-zvox") == 0 && (i)<argc)
            {
                cout <<(double)(Images[0]->dz)<<endl;
                flush(cout);
            }

            // **************************            ---------          *****************************
            // **************************            Im dim X/Y/Z           *****************************
            // **************************            ---------          *****************************
            else if(strcmp(argv[i], "-xdim") == 0 && (i)<argc)
            {
                cout <<(double)(Images[0]->nx)<<endl;
                flush(cout);
            }
            else if(strcmp(argv[i], "-ydim") == 0 && (i)<argc)
            {
                cout <<(double)(Images[0]->ny)<<endl;
                flush(cout);
            }
            else if(strcmp(argv[i], "-zdim") == 0 && (i)<argc)
            {
                cout <<(double)(Images[0]->nz)<<endl;
                flush(cout);
            }



#ifdef _GIT_HASH
            else if( strcmp(argv[i], "--version")==0)
            {
                printf("%s\n",_GIT_HASH);
                return 0;
            }
#endif
            // **************************            ---------          *****************************
            // **************************               HELP            *****************************
            // **************************            ---------          *****************************
            else
            {
                printf("Err:\tParameter %s unknown or incomplete\n\n",argv[i]);
                flush(cout);
                Usage(argv[0]);
                return 1;
            }
        }

        if(imgsort!=NULL)
            delete [] imgsort;

        for(int i=0; i<numbimg; i++)
        {
            nifti_image_free(Images[i]);
        }
        nifti_image_free(Mask);
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

