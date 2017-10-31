/**
 * @file seg_maths.cpp
 * @author M. Jorge Cardoso
 * @date 01/01/2014
 *
 * Copyright (c) 2014, University College London. All rights reserved.
 * Centre for Medical Image Computing (CMIC)
 * See the LICENSE.txt file in the nifty_seg root folder
 *
 */

#include <iostream>
#include <time.h>
#include "_seg_common.h"
#include "_seg_tools.h"
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <cfloat>

using namespace std;
#define SegPrecisionTYPE float

void Usage(char *exec)
{
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    printf("\nMath tools:\nUsage:\t%s <input> <operation> <output>.\n\n",exec);
    printf("\t* * Operations on 3-D and 4-D images* *\n");
    printf("\t-mul\t<float/file>\tMultiply image <float> value or by other image.\n");
    printf("\t-div\t<float/file>\tDivide image by <float> or by other image.\n");
    printf("\t-add\t<float/file>\tAdd image by <float> or by other image.\n");
    printf("\t-sub\t<float/file>\tSubtract image by <float> or by other image.\n");
    printf("\t-pow\t<float>\t\tImage to the power of <float>.\n");
    printf("\t-thr\t<float>\t\tThreshold the image below <float>.\n");
    printf("\t-uthr\t<float>\t\tThreshold image above <float>.\n");
    printf("\t-smo\t<float>\t\tGaussian smoothing by std <float> (in voxels and up to 4-D).\n");
    printf("\t-equal\t<int>\t\tGet voxels equal to <int>\n");
    printf("\t-replace <int1> <int2>\tReplaces voxels equal to <int1> with <int2>\n");
    printf("\t-sqrt \t\t\tSquare root of the image.\n");
    printf("\t-exp \t\t\tExponential root of the image.\n");
    printf("\t-log \t\t\tLog of the image.\n");
    printf("\t-recip \t\t\tReciprocal (1/I) of the image.\n");
    printf("\t-abs \t\t\tAbsolute value of the image.\n");
    printf("\t-bin \t\t\tBinarise the image.\n");
    printf("\t-otsu \t\t\tOtsu thresholding of the current image.\n");
    printf("\t-edge\t<float>\t\tCalculate the edges of the image using a threshold <float>.\n");
    printf("\t-sobel3\t<float>\t\tCalculate the edges of all timepoints using a Sobel filter with a 3x3x3 kernel and applying <float> gaussian smoothing.\n");
    printf("\t-sobel5\t<float>\t\tCalculate the edges of all timepoints using a Sobel filter with a 5x5x5 kernel and applying <float> gaussian smoothing.\n");
    printf("\t-min\t<file>\t\tGet the min per voxel between <current> and <file>.\n");
    printf("\n\t* * Operations on 3-D images * *\n");
    printf("\t-smol\t<float>\t\tGaussian smoothing of a 3D label image.\n");
    printf("\t-dil\t<int>\t\tDilate the image <int> times (in voxels).\n");
    printf("\t-ero\t<int>\t\tErode the image <int> times (in voxels).\n");
    printf("\t-pad\t<int>\t\tPad <int> voxels with NaN value around each 3D volume.\n");
    printf("\t-crop\t<int>\t\tCrop <int> voxels around each 3D volume.\n");    
    printf("\n\t* * Operations binary 3-D images * *\n");
    printf("\t-lconcomp\t\tTake the largest connected component\n");
    printf("\t-concomp6\t\tLabel the different connected components with a 6NN kernel\n");
    printf("\t-concomp26\t\tLabel the different connected components with a 26NN kernel\n");
    printf("\t-fill\t\t\tFill holes in binary object (e.g. fill ventricle in brain mask).\n");
    printf("\t-euc\t\t\tEuclidean distance trasnform\n");
    printf("\t-geo <float/file>\tGeodesic distance according to the speed function <float/file>\n");
    printf("\n\t* * Dimensionality reduction operations: from 4-D to 3-D * *\n");
    printf("\t-tp <int>\t\tExtract time point <int>\n");
    printf("\t-tpmax\t\t\tGet the time point with the highest value (binarise 4D probabilities)\n");
    printf("\t-tmean\t\t\tMean value of all time points.\n");
    printf("\t-tmax\t\t\tMax value of all time points.\n");
    printf("\t-tmin\t\t\tMean value of all time points.\n");
    printf("\n\t* * Dimensionality increase operations: from 3-D to 4-D * *\n");
    printf("\t-merge\t<i> <d> <files>\tMerge <i> images and the working image in the <d> dimension \n");
    printf("\t-splitlab\t\tSplit the integer labels into multiple timepoints\n");
    printf("\t-splitinter <x/y/z>\t\tSplit interleaved slices in direction <x/y/z> into separate time points\n");
    printf("\n\t* * Image similarity: Local metrics * *\n");
    printf("\t-lncc\t<file> <std>\tLocal CC between current img and <file> on a kernel with <std>\n");
    printf("\t-lssd\t<file> <std>\tLocal SSD between current img and <file> on a kernel with <std>\n");
    printf("\n\t* * Normalisation * *\n");
    printf("\t-llsnorm\t<file_norm>\t\t Linear LS normalisation between current and <file_norm>\n");
    printf("\t-lltsnorm\t<file_norm> <float>\t Linear LTS normalisation assuming <float> percent outliers\n");
    printf("\t-qlsnorm\t<order> <file_norm>\t LS normalisation of <order> between current and <file_norm>\n");
    printf("\n\t* * NaN handling * *\n");
    printf("\t-removenan\t\tRemove all NaNs and replace then with 0\n");
    printf("\t-isnan\t\t\tBinary image equal to 1 if the value is NaN and 0 otherwise\n");
    printf("\t-masknan <file_norm>\tAssign everything outside the mask (mask==0) with NaNs \n");
    printf("\n\t* * Sampling * *\n");
    printf("\t-subsamp2\t\tSubsample the image by 2 using NN sampling (qform and sform scaled) \n");
    printf("\n\t* * Image header operations * *\n");
    printf("\t-hdr_copy <file> \tCopy header from working image to <file> and save in <output>.\n");
    printf("\t-scl\t\t\tReset scale and slope info.\n");
    printf("\t-4to5\t\t\tFlip the 4th and 5th dimension.\n");
    printf("\n\t* * Output * *\n");
    printf("\t-odt <datatype> \tSet output <datatype> (char, short, int, uchar, ushort, uint, float, double).\n");
    printf("\t-range\t\t\tReset the image range to the min max\n");
    printf("\t-v\t\t\tVerbose.\n");
#if defined (_OPENMP)
    printf("\t-omp <int>\t\tNumber of openmp threads [%d]\n",omp_get_max_threads());
#endif
#ifdef _GIT_HASH
    printf("\t--version\t\tPrint current source code git hash key and exit\n\t\t\t\t(%s)\n",_GIT_HASH);
#endif
    printf("\n\t* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    return;
}

bool isEdge(float a,float b,double treshold) {
    float max=a>b?a:b;
    return (fabs(a-b)/max>treshold);
}

int isNumeric (const char *s)
{
    if(s==NULL || *s=='\0' || isspace(*s))
        return 0;
    char * p;
    strtod (s, &p);
    return *p == '\0';
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
        if (argc <= 2)
        {
            Usage(argv[0]);
            return 0;
        }
        if(strcmp(argv[1], "-help")==0 || strcmp(argv[1], "-Help")==0 ||
                strcmp(argv[1], "-HELP")==0 || strcmp(argv[1], "-h")==0 ||
                strcmp(argv[1], "--h")==0 || strcmp(argv[1], "--help")==0)
        {
            Usage(argv[0]);
            return 0;
        }


        char * filename_in=argv[1];
        nifti_image * InputImage=nifti_image_read(filename_in,true);
        if(InputImage == NULL)
        {
            fprintf(stderr,"* Error when reading the input image\n");
            return 1;
        }
        if(InputImage->datatype!=NIFTI_TYPE_FLOAT32)
        {
            seg_changeDatatype<SegPrecisionTYPE>(InputImage);
        }
        SegPrecisionTYPE * InputImagePtr = static_cast<SegPrecisionTYPE *>(InputImage->data);
        ImageSize * CurrSize = new ImageSize [1]();
        CurrSize->numel=(long)(InputImage->nx*InputImage->ny*InputImage->nz);
        CurrSize->xsize=InputImage->nx;
        CurrSize->ysize=InputImage->ny;
        CurrSize->zsize=InputImage->nz;
        CurrSize->usize=(InputImage->nu>1)?InputImage->nu:1;
        CurrSize->tsize=(InputImage->nt>1)?InputImage->nt:1;
        float Scalling[4]= { 1.0f, 1.0f, 1.0f, 1.0f };
        bool verbose=0;
        int datatypeoutput=NIFTI_TYPE_FLOAT32;

        SegPrecisionTYPE ** bufferImages = new SegPrecisionTYPE * [2];
        bufferImages[0] = new SegPrecisionTYPE [CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize];
        bufferImages[1] = new SegPrecisionTYPE [CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize];
        for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
        {
            bufferImages[0][i]=InputImagePtr[i];
        }
        int current_buffer=0;

        for(long i=2; i<(argc-1); i++)
        {
            if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
                    strcmp(argv[i], "-HELP")==0 || strcmp(argv[i], "-h")==0 ||
                    strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0)
            {
                Usage(argv[0]);
                return 0;
            }
#if defined (_OPENMP)
            else if(strcmp(argv[i], "-omp")==0 || strcmp(argv[i], "--omp")==0)
            {
                omp_set_num_threads(atoi(argv[++i]));
            }
#endif
            // *********************  MUTIPLY  *************************
            else if(strcmp(argv[i], "-mul") == 0)
            {
                string parser=argv[++i];
                if(parser.find_first_not_of("1234567890.-+")== string::npos)
                {
                    double multfactor=strtod(parser.c_str(),NULL);
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                    {
                        bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i]*multfactor;
                    }
                    current_buffer=current_buffer?0:1;
                }
                else
                {
                    nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
                    NewImage->nu=(NewImage->nu>1)?NewImage->nu:1;
                    NewImage->nt=(NewImage->nt>1)?NewImage->nt:1;
                    if(NewImage->datatype!=DT_FLOAT32)
                    {
                        seg_changeDatatype<SegPrecisionTYPE>(NewImage);
                    }
                    SegPrecisionTYPE * NewImagePtr = static_cast<SegPrecisionTYPE *>(NewImage->data);
                    if(NewImage->nx==CurrSize->xsize&&NewImage->ny==CurrSize->ysize&&NewImage->nz==CurrSize->zsize&&NewImage->nt==CurrSize->tsize&&NewImage->nu==CurrSize->usize)
                    {
                        for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                            bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i]*NewImagePtr[i];
                        current_buffer=current_buffer?0:1;
                    }
                    else
                    {
                        nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
                        NewImage->nu=(NewImage->nu>1)?NewImage->nu:1;
                        NewImage->nt=(NewImage->nt>1)?NewImage->nt:1;
                        if(NewImage->datatype!=DT_FLOAT32)
                        {
                            seg_changeDatatype<SegPrecisionTYPE>(NewImage);
                        }
                        SegPrecisionTYPE * NewImagePtr = static_cast<SegPrecisionTYPE *>(NewImage->data);
                        if(NewImage->nx==CurrSize->xsize&&NewImage->ny==CurrSize->ysize&&NewImage->nz==CurrSize->zsize&&NewImage->nt==CurrSize->tsize&&NewImage->nu==CurrSize->usize)
                        {
                            for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                                bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i]*NewImagePtr[i];
                            current_buffer=current_buffer?0:1;
                        }
                        else
                        {
                            cout << "ERROR: Image "<< parser << " is the wrong size  -  original = ( "<<CurrSize->xsize<<","
                                 <<CurrSize->ysize<<","<<CurrSize->zsize<<","<<CurrSize->tsize<<","<<CurrSize->usize<<" ) New image = ( "
                                <<NewImage->nx<<","<<NewImage->ny<<","<<NewImage->nz<<","<<NewImage->nt<<","<<NewImage->nu<<" )"<<endl;
                            i=argc;
                        }
                        nifti_image_free(NewImage);
                    }
                }
            }
            // *********************  ADD  *************************
            else if( strcmp(argv[i], "-add") == 0)
            {
                string parser=argv[++i];
                if(parser.find_first_not_of("1234567890.-+")== string::npos)
                {
                    double addfactor=strtod(parser.c_str(),NULL);
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                        bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i]+addfactor;
                    current_buffer=current_buffer?0:1;
                }
                else
                {
                    nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
                    NewImage->nu=(NewImage->nu>1)?NewImage->nu:1;
                    NewImage->nt=(NewImage->nt>1)?NewImage->nt:1;
                    if(NewImage->datatype!=DT_FLOAT32)
                    {
                        seg_changeDatatype<SegPrecisionTYPE>(NewImage);
                    }
                    SegPrecisionTYPE * NewImagePtr = static_cast<SegPrecisionTYPE *>(NewImage->data);
                    if(NewImage->nx==CurrSize->xsize&&NewImage->ny==CurrSize->ysize&&NewImage->nz==CurrSize->zsize&&NewImage->nt==CurrSize->tsize&&NewImage->nu==CurrSize->usize)
                    {
                        for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                            bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i]+NewImagePtr[i];
                        current_buffer=current_buffer?0:1;
                    }
                    else
                    {
                        cout << "ERROR: Image "<< parser << " is the wrong size  -  original = ( "<<CurrSize->xsize<<","
                             <<CurrSize->ysize<<","<<CurrSize->zsize<<","<<CurrSize->tsize<<","<<CurrSize->usize<<" )  New image = ( "<<NewImage->nx<<","
                            <<NewImage->ny<<","<<NewImage->nz<<","<<NewImage->nt<<","<<NewImage->nu<<" )"<<endl;
                        i=argc;
                    }
                    nifti_image_free(NewImage);
                }
            }
            // *********************  SUBTRACT  *************************
            else if(strcmp(argv[i], "-sub") == 0)
            {
                string parser=argv[++i];
                if(parser.find_first_not_of("1234567890.-+")== string::npos)
                {
                    double factor=strtod(parser.c_str(),NULL);
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                        bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i]-factor;
                    current_buffer=current_buffer?0:1;
                }
                else
                {
                    nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
                    NewImage->nu=(NewImage->nu>1)?NewImage->nu:1;
                    NewImage->nt=(NewImage->nt>1)?NewImage->nt:1;
                    if(NewImage->datatype!=DT_FLOAT32)
                    {
                        seg_changeDatatype<SegPrecisionTYPE>(NewImage);
                    }
                    SegPrecisionTYPE * NewImagePtr = static_cast<SegPrecisionTYPE *>(NewImage->data);
                    if(NewImage->nx==CurrSize->xsize&&NewImage->ny==CurrSize->ysize&&NewImage->nz==CurrSize->zsize&&NewImage->nt==CurrSize->tsize&&NewImage->nu==CurrSize->usize)
                    {
                        for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                            bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i]-NewImagePtr[i];
                        current_buffer=current_buffer?0:1;
                    }
                    else
                    {
                        cout << "ERROR: Image "<< parser << " is the wrong size  -  original = ( "<<CurrSize->xsize<<","
                             <<CurrSize->ysize<<","<<CurrSize->zsize<<","<<CurrSize->tsize<<","<<CurrSize->usize<<" )  New image = ( "<<NewImage->nx<<","
                            <<NewImage->ny<<","<<NewImage->nz<<","<<NewImage->nt<<","<<NewImage->nu<<" )"<<endl;
                        i=argc;
                    }
                    nifti_image_free(NewImage);
                }
            }
            // *********************  mask  *************************
            else if(strcmp(argv[i], "-masknan") == 0)
            {
                string parser=argv[++i];

                nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
                NewImage->nu=(NewImage->nu>1)?NewImage->nu:1;
                NewImage->nt=(NewImage->nt>1)?NewImage->nt:1;
                if(NewImage->datatype!=DT_FLOAT32)
                {
                    seg_changeDatatype<SegPrecisionTYPE>(NewImage);
                }
                SegPrecisionTYPE * NewImagePtr = static_cast<SegPrecisionTYPE *>(NewImage->data);
                if(NewImage->nx==CurrSize->xsize&&NewImage->ny==CurrSize->ysize&&NewImage->nz==CurrSize->zsize&&NewImage->nt==CurrSize->tsize&&NewImage->nu==CurrSize->usize)
                {
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                        bufferImages[current_buffer?0:1][i]=(NewImagePtr[i]>0)?bufferImages[current_buffer][i]:std::numeric_limits<float>::quiet_NaN();
                    current_buffer=current_buffer?0:1;
                }
                else
                {
                    cout << "ERROR: Image "<< parser << " is the wrong size  -  original = ( "<<CurrSize->xsize<<","
                         <<CurrSize->ysize<<","<<CurrSize->zsize<<","<<CurrSize->tsize<<","<<CurrSize->usize<<" )  New image = ( "<<NewImage->nx<<","
                        <<NewImage->ny<<","<<NewImage->nz<<","<<NewImage->nt<<","<<NewImage->nu<<" )"<<endl;
                    i=argc;
                }
                nifti_image_free(NewImage);

            }
            // *********************  mask  *************************
            else if(strcmp(argv[i], "-removenan") == 0)
            {

                for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                    bufferImages[current_buffer?0:1][i]=isnan(bufferImages[current_buffer][i])==1?0:bufferImages[current_buffer][i];
                current_buffer=current_buffer?0:1;

            }
	    // *********************  pad voxels  *************************
            else if(strcmp(argv[i], "-pad") == 0)
            {
                string parser=argv[++i];

                if(parser.find_first_not_of("1234567890-+")== string::npos)
                {
                    int padding=(int)strtod(parser.c_str(),NULL)*2;
		    long new_size=((CurrSize->xsize+padding)*(CurrSize->ysize+padding)*(CurrSize->zsize+padding)*CurrSize->tsize*CurrSize->usize);
		    bufferImages[current_buffer?0:1]=new SegPrecisionTYPE [new_size];
		           
		    for(long ii=0; ii<new_size; ii++)
                        bufferImages[current_buffer?0:1][ii]=std::numeric_limits<double>::quiet_NaN();

		    long old_volume=CurrSize->xsize*CurrSize->ysize*CurrSize->zsize;
		    long new_volume=(CurrSize->xsize+padding)*(CurrSize->ysize+padding)*(CurrSize->zsize+padding);
		    for (long t=0;t<CurrSize->tsize*CurrSize->usize;t++) {
			    for(long z=0; z<CurrSize->zsize; z++) {
				for(long y=0; y<CurrSize->ysize; y++) {
					for(long x=0; x<CurrSize->xsize; x++) {
						long big=t*new_volume+x+(padding/2)+(y+(padding/2))*(CurrSize->xsize+padding)+(z+(padding/2))*(CurrSize->xsize+padding)*(CurrSize->ysize+padding);
						long small=t*old_volume+x+y*CurrSize->xsize+z*(CurrSize->xsize*CurrSize->ysize);
        	                		bufferImages[current_buffer?0:1][big]=bufferImages[current_buffer][small];
					}
				}
			    }
		    }
                    current_buffer=current_buffer?0:1;
		    bufferImages[current_buffer?0:1]=new SegPrecisionTYPE [new_size];
		    for(long ii=0; ii<new_size; ii++)
                        bufferImages[current_buffer?0:1][ii]=0;
		    CurrSize->xsize+=padding;
		    CurrSize->ysize+=padding;
		    CurrSize->zsize+=padding;
    		    CurrSize->numel=CurrSize->xsize*CurrSize->ysize*CurrSize->zsize;
                }
		else
                {
                    cout << "ERROR: "<< parser << " is not a valid number"<<endl;
                    i=argc;
                }
            }
	    // *********************  crop voxels  *************************
            else if(strcmp(argv[i], "-crop") == 0)
            {
                string parser=argv[++i];

                if(parser.find_first_not_of("1234567890-+")== string::npos)
                {
                    int cropping=(int)strtod(parser.c_str(),NULL)*2;
		    long new_size=((CurrSize->xsize-cropping)*(CurrSize->ysize-cropping)*(CurrSize->zsize-cropping)*CurrSize->tsize*CurrSize->usize);
		    bufferImages[current_buffer?0:1]=new SegPrecisionTYPE [new_size];
		           
		    for(long ii=0; ii<new_size; ii++)
                        bufferImages[current_buffer?0:1][ii]=0;
		
		    long old_volume=CurrSize->xsize*CurrSize->ysize*CurrSize->zsize;
		    long new_volume=(CurrSize->xsize-cropping)*(CurrSize->ysize-cropping)*(CurrSize->zsize-cropping);
		    for (long t=0;t<CurrSize->tsize*CurrSize->usize;t++) {
                    	for(long x=cropping/2; x<CurrSize->xsize-cropping/2; x++) {
				for(long y=cropping/2; y<CurrSize->ysize-cropping/2; y++) {
					for(long z=cropping/2; z<CurrSize->zsize-cropping/2; z++) {
						long small=t*new_volume+x-(cropping/2)+(y-(cropping/2))*(CurrSize->xsize-cropping)+(z-(cropping/2))*((CurrSize->xsize-cropping)*(CurrSize->ysize-cropping));
						long big=t*old_volume+x+y*CurrSize->xsize+z*(CurrSize->xsize*CurrSize->ysize);
                    		    		bufferImages[current_buffer?0:1][small]=bufferImages[current_buffer][big];
					}
				}
		    	}
		    }
                    current_buffer=current_buffer?0:1;
		    bufferImages[current_buffer?0:1]=new SegPrecisionTYPE [new_size];
		    for(long ii=0; ii<new_size; ii++)
                        bufferImages[current_buffer?0:1][ii]=0;
		    CurrSize->xsize-=cropping;
		    CurrSize->ysize-=cropping;
		    CurrSize->zsize-=cropping;
		    CurrSize->numel=CurrSize->xsize*CurrSize->ysize*CurrSize->zsize;
                }
		else
                {
                    cout << "ERROR: "<< parser << " is not a valid number"<<endl;
                    i=argc;
                }
            }
	    //  *********************  mask edge  *************************
            else if(strcmp(argv[i], "-edge") == 0)
            {
                string parser=argv[++i];
                if(((strtod(parser.c_str(),NULL)!=0) || (parser.length()==1 && parser.find("0")!=string::npos)))
                {
                    double treshold=strtod(parser.c_str(),NULL);;
                    float * Img1prt = bufferImages[current_buffer];
                    for(int index=0; index<CurrSize->xsize*CurrSize->ysize*CurrSize->zsize; index++)
                    {
                        bool edge=false;
                        if((index+1)<CurrSize->xsize*CurrSize->ysize*CurrSize->zsize) {
                            if(isEdge(Img1prt[index],Img1prt[index+1],treshold)) edge=true;
                        }
                        if((index-1)>0) {
                            if(isEdge(Img1prt[index],Img1prt[index-1],treshold)) edge=true;
                        }
                        if((index+CurrSize->xsize)<CurrSize->xsize*CurrSize->ysize*CurrSize->zsize) {
                            if(isEdge(Img1prt[index],Img1prt[index+CurrSize->xsize],treshold)) edge=true;
                        }
                        if((index-CurrSize->xsize)>0) {
                            if(isEdge(Img1prt[index],Img1prt[index-CurrSize->xsize],treshold)) edge=true;
                        }
                        if((index+CurrSize->xsize*CurrSize->ysize)<CurrSize->xsize*CurrSize->ysize*CurrSize->zsize) {
                            if(isEdge(Img1prt[index],Img1prt[index+CurrSize->xsize*CurrSize->ysize],treshold)) edge=true;
                        }
                        if((index-CurrSize->xsize*CurrSize->ysize)>0) {
                            if(isEdge(Img1prt[index],Img1prt[index-CurrSize->xsize*CurrSize->ysize],treshold)) edge=true;
                        }
                        if(edge)
                        {
                            bufferImages[current_buffer?0:1][index]=bufferImages[current_buffer][index];
                        }
                        else {
                            bufferImages[current_buffer?0:1][index]=0;
                        }
                    }
                    current_buffer=current_buffer?0:1;
                }
                else
                {
                    cout << "ERROR: "<< parser << " is not a valid number"<<endl;
                    i=argc;
                }
            }
	    else if(strcmp(argv[i], "-sobel3") == 0)
            {
                string parser=argv[++i];
                if(((strtod(parser.c_str(),NULL)!=0) || (parser.length()==1 && parser.find("0")!=string::npos)))
                {
                    double factor=strtod(parser.c_str(),NULL);
                    float * Img1prt = bufferImages[current_buffer];
                    long tp=0;
                    #ifdef _OPENMP
                    #pragma omp parallel for \
                        private(tp)\
                        shared(CurrSize,bufferImages,Img1prt,factor,InputImage)
                    #endif
                    for(tp=0; tp<(long)(CurrSize->tsize*CurrSize->usize); tp++){
                        //create dummy nii
                        nifti_image * TMPnii = nifti_copy_nim_info(InputImage);
                        TMPnii->dim[1]=CurrSize->xsize;
                        TMPnii->dim[2]=CurrSize->ysize;
                        TMPnii->dim[3]=CurrSize->zsize;
                        TMPnii->dim[4]=TMPnii->nt=1;
                        TMPnii->dim[5]=TMPnii->nu=1;
                        nifti_update_dims_from_array(TMPnii);
                        //copy pointer, run gaussian, and set to null
                        TMPnii->data=static_cast<void*>(&Img1prt[CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*tp]);
                        if(factor>0) GaussianSmoothing5D_nifti(TMPnii,NULL,factor);
                        TMPnii->data=NULL;
                        //As TMPnii->data=NULL, the free will not cause any harm
                        nifti_image_free(TMPnii);

                        float *imgsort=new float [CurrSize->xsize*CurrSize->ysize*CurrSize->zsize];
                        for(long i=0; i<CurrSize->xsize*CurrSize->ysize*CurrSize->zsize; i++) {
                            imgsort[i]=Img1prt[i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize];
                        }
                        HeapSort(imgsort,CurrSize->xsize*CurrSize->ysize*CurrSize->zsize-1);
                        float max=imgsort[(int)(round((1-0.02)*(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize-1)))];
                        float min=imgsort[(int)(round(0.02*(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize-1)))];
                        float newMax=1,newMin=0;
                        for(long i=0; i<CurrSize->xsize*CurrSize->ysize*CurrSize->zsize; i++) {
                            if(min>Img1prt[i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize]) Img1prt[i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize]=min;
                            if(max<Img1prt[i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize]) Img1prt[i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize]=max;
                            Img1prt[i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize]=newMin+(Img1prt[i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize]-min)*(newMax-newMin)/(max-min);
                        }
                        int inz=0;
                        float xkernel[3][3][3]={
                                            {{-1,-2,-1},{-2,-4,-2},{-1,-2,-1}},
                                            {{ 0, 0, 0},{ 0, 0, 0},{ 0, 0, 0}},
                                            {{ 1, 2, 1},{ 2, 4, 2},{ 1, 2, 1}}
                                         };
                        float ykernel[3][3][3]={
                                            {{ 1, 2, 1},{ 0, 0, 0},{-1,-2,-1}},
                                            {{ 2, 4, 2},{ 0, 0, 0},{-2,-4,-2}},
                                            {{ 1, 2, 1},{ 0, 0, 0},{-1,-2,-1}}
                                         };
                        float zkernel[3][3][3]={
                                            {{-1, 0, 1},{-2, 0, 2},{-1, 0, 1}},
                                            {{-2, 0, 2},{-4, 0, 4},{-2, 0, 2}},
                                            {{-1, 0, 1},{-2, 0, 2},{-1, 0, 1}}
                                         };
                        #ifdef _OPENMP
                        #pragma omp parallel for \
                            private(inz)\
                            shared(CurrSize,bufferImages,Img1prt,xkernel,ykernel,zkernel)
                        #endif
                        for(inz=0; inz<CurrSize->zsize; inz++) {
                            for(int iny=0; iny<CurrSize->ysize; iny++) {
                                for(int inx=0; inx<CurrSize->xsize; inx++) {
                                    float sumx=0,sumy=0,sumz=0;
                                    for(int i=-1;i<=1;i++) {
                                        for(int j=-1;j<=1;j++) {
                                            for(int k=-1;k<=1;k++) {
                                                if(inx+k>=0 && iny+j>=0 && inz+i>=0 &&
                                                        inx+k<CurrSize->xsize && iny+j<CurrSize->ysize && inz+i<CurrSize->zsize) {
                                                    int index=(inx+k)+(iny+j)*CurrSize->xsize+(inz+i)*(CurrSize->xsize*CurrSize->ysize);
                                                    sumx+=xkernel[k+1][j+1][i+1]*Img1prt[index+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize];
                                                    sumy+=ykernel[j+1][k+1][i+1]*Img1prt[index+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize];
                                                    sumz+=zkernel[i+1][k+1][j+1]*Img1prt[index+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize];
                                                }
                                            }
                                        }
                                    }

                                    int index=inx+iny*CurrSize->xsize+inz*(CurrSize->xsize*CurrSize->ysize);
                                    float val=sqrt(sumx*sumx+sumy*sumy+sumz*sumz);
                                    bufferImages[current_buffer?0:1][index+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize]=val;
                                }
                            }
                        }
                        for(long i=0; i<CurrSize->xsize*CurrSize->ysize*CurrSize->zsize; i++) {
                            imgsort[i]=bufferImages[current_buffer?0:1][i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize];
                        }
                        HeapSort(imgsort,CurrSize->xsize*CurrSize->ysize*CurrSize->zsize-1);
                        max=imgsort[(int)(round((1-0.02)*(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize-1)))];
                        min=imgsort[(int)(round(0.02*(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize-1)))];
                        for(long i=0; i<CurrSize->xsize*CurrSize->ysize*CurrSize->zsize; i++) {
                            if(min>bufferImages[current_buffer?0:1][i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize]) bufferImages[current_buffer?0:1][i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize]=min;
                            if(max<bufferImages[current_buffer?0:1][i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize]) bufferImages[current_buffer?0:1][i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize]=max;
                            bufferImages[current_buffer?0:1][i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize]=newMin+(bufferImages[current_buffer?0:1][i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize]-min)*(newMax-newMin)/(max-min);
                        }
                    }
                }
                else {
                    cout << "ERROR: "<< parser << " is not a valid number"<<endl;
                    i=argc;
                }
                current_buffer=current_buffer?0:1;
            }
            else if(strcmp(argv[i], "-sobel5") == 0)
            {
                string parser=argv[++i];
                if(((strtod(parser.c_str(),NULL)!=0) || (parser.length()==1 && parser.find("0")!=string::npos)))
                {
                    double factor=strtod(parser.c_str(),NULL);
                    float * Img1prt = bufferImages[current_buffer];
                    long tp=0;
                    #ifdef _OPENMP
                    #pragma omp parallel for \
                        private(tp)\
                        shared(CurrSize,bufferImages,Img1prt,factor,InputImage)
                    #endif
                    for(tp=0; tp<(long)(CurrSize->tsize*CurrSize->usize); tp++){
                        //create dummy nii
                        nifti_image * TMPnii = nifti_copy_nim_info(InputImage);
                        TMPnii->dim[1]=CurrSize->xsize;
                        TMPnii->dim[2]=CurrSize->ysize;
                        TMPnii->dim[3]=CurrSize->zsize;
                        TMPnii->dim[4]=TMPnii->nt=1;
                        TMPnii->dim[5]=TMPnii->nu=1;
                        nifti_update_dims_from_array(TMPnii);
                        //copy pointer, run gaussian, and set to null
                        TMPnii->data=static_cast<void*>(&Img1prt[CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*tp]);
                        if(factor>0) GaussianSmoothing5D_nifti(TMPnii,NULL,factor);
                        TMPnii->data=NULL;
                        //As TMPnii->data=NULL, the free will not cause any harm
                        nifti_image_free(TMPnii);

                        float *imgsort=new float [CurrSize->xsize*CurrSize->ysize*CurrSize->zsize];
                        for(long i=0; i<CurrSize->xsize*CurrSize->ysize*CurrSize->zsize; i++) {
                            imgsort[i]=Img1prt[i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize];
                        }
                        HeapSort(imgsort,CurrSize->xsize*CurrSize->ysize*CurrSize->zsize-1);
                        float max=imgsort[(int)(round((1-0.02)*(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize-1)))];
                        float min=imgsort[(int)(round(0.02*(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize-1)))];
                        float newMax=1,newMin=0;
                        for(long i=0; i<CurrSize->xsize*CurrSize->ysize*CurrSize->zsize; i++) {
                            if(min>Img1prt[i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize]) Img1prt[i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize]=min;
                            if(max<Img1prt[i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize]) Img1prt[i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize]=max;
                            Img1prt[i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize]=newMin+(Img1prt[i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize]-min)*(newMax-newMin)/(max-min);
                        }
                        float xkernel[5][5][5]={
                                            {{-1,-4, -6,-4,-1},{-2, -8,-12, -8,-2},{-4,-16,-24,-16,-4},{-2, -8,-12, -8,-2},{-1,-4, -6,-4,-1}},
                                            {{-2,-8,-12,-8,-2},{-4,-16,-24,-16,-4},{-8,-32,-48,-32,-8},{-4,-16,-24,-16,-4},{-2,-8,-12,-8,-2}},
                                            {{ 0, 0,  0, 0, 0},{ 0,  0,  0,  0, 0},{ 0,  0,  0,  0, 0},{ 0,  0,  0,  0, 0},{ 0, 0,  0, 0, 0}},
                                            {{ 2, 8, 12, 8, 2},{ 4, 16, 24, 16, 4},{ 8, 32, 48, 32, 8},{ 4, 16, 24, 16, 4},{ 2, 8, 12, 8, 2}},
                                            {{ 1, 4,  6, 4, 1},{ 2,  8, 12,  8, 2},{ 4, 16, 24, 16, 4},{ 2,  8, 12,  8, 2},{ 1, 4,  6, 4, 1}}
                                         };
                        float ykernel[5][5][5]={
                                            {{ 1,  4,  6,  4, 1},{ 2,  8, 12,  8, 2},{ 0, 0, 0, 0, 0},{-2, -8,-12, -8,-2},{-1, -4, -6, -4,-1}},
                                            {{ 2,  8, 12,  8, 2},{ 4, 16, 24, 16, 4},{ 0, 0, 0, 0, 0},{-4,-16,-24,-16,-4},{-2, -8,-12, -8,-2}},
                                            {{ 4, 16, 24, 16, 4},{ 8, 32, 48, 32, 8},{ 0, 0, 0, 0, 0},{-8,-32,-48,-32,-8},{-4,-16,-24,-16,-4}},
                                            {{ 2,  8, 12,  8, 2},{ 4, 16, 24, 16, 4},{ 0, 0, 0, 0, 0},{-4,-16,-24,-16,-4},{-2, -8,-12, -8,-2}},
                                            {{ 1,  4,  6,  4, 1},{ 2,  8, 12,  8, 2},{ 0, 0, 0, 0, 0},{-2, -8,-12, -8,-2},{-1, -4, -6, -4,-1}}
                                         };
                        float zkernel[5][5][5]={
                                            {{-1, -2,  0,  2, 1},{ -2, -4, 0,  4,  2},{ -4, -8, 0,  8,  4},{ -2, -4, 0,  4,  2},{-1, -2,  0,  2, 1}},
                                            {{-4, -8,  0,  8, 4},{ -8,-16, 0, 16,  8},{-16,-32, 0, 32, 16},{ -8,-16, 0, 16,  8},{-4, -8,  0,  8, 4}},
                                            {{-6,-12,  0, 12, 6},{-12,-24, 0, 24, 12},{-24,-48, 0, 48, 24},{-12,-24, 0, 24, 12},{-6,-12,  0, 12, 6}},
                                            {{-4, -8,  0,  8, 4},{ -8,-16, 0, 16,  8},{-16,-32, 0, 32, 16},{ -8,-16, 0, 16,  8},{-4, -8,  0,  8, 4}},
                                            {{-1, -2,  0,  2, 1},{- 2, -4, 0,  4,  2},{ -4,  8, 0,  8,  4},{ -2, -4, 0,  4,  2},{-1, -2,  0,  2, 1}}
                                         };
                        int inz=0;
                        #ifdef _OPENMP
                        #pragma omp parallel for \
                            private(inz)\
                            shared(CurrSize,bufferImages,Img1prt,xkernel,ykernel,zkernel)
                        #endif
                        for(inz=0; inz<CurrSize->zsize; inz++) {
                            for(int iny=0; iny<CurrSize->ysize; iny++) {
                                for(int inx=0; inx<CurrSize->xsize; inx++) {
                                    float sumx=0,sumy=0,sumz=0;
                                    for(int i=-2;i<=2;i++) {
                                        for(int j=-2;j<=2;j++) {
                                            for(int k=-2;k<=2;k++) {
                                                if(inx+k>=0 && iny+j>=0 && inz+i>=0 &&
                                                        inx+k<CurrSize->xsize && iny+j<CurrSize->ysize && inz+i<CurrSize->zsize) {
                                                    int index=(inx+k)+(iny+j)*CurrSize->xsize+(inz+i)*(CurrSize->xsize*CurrSize->ysize);
                                                    sumx+=xkernel[k+2][j+2][i+2]*Img1prt[index+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize];
                                                    sumy+=ykernel[j+2][k+2][i+2]*Img1prt[index+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize];
                                                    sumz+=zkernel[i+2][k+2][j+2]*Img1prt[index+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize];
                                                }
                                            }
                                        }
                                    }
                                    int index=inx+iny*CurrSize->xsize+inz*(CurrSize->xsize*CurrSize->ysize);
                                    float val=sqrt(sumx*sumx+sumy*sumy+sumz*sumz);
                                    bufferImages[current_buffer?0:1][index+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize]=val;
                                }
                            }
                        }
                        for(long i=0; i<CurrSize->xsize*CurrSize->ysize*CurrSize->zsize; i++) {
                            imgsort[i]=bufferImages[current_buffer?0:1][i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize];
                        }
                        HeapSort(imgsort,CurrSize->xsize*CurrSize->ysize*CurrSize->zsize-1);
                        max=imgsort[(int)(round((1-0.02)*(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize-1)))];
                        min=imgsort[(int)(round(0.02*(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize-1)))];
                        for(long i=0; i<CurrSize->xsize*CurrSize->ysize*CurrSize->zsize; i++) {
                            if(min>bufferImages[current_buffer?0:1][i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize]) bufferImages[current_buffer?0:1][i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize]=min;
                            if(max<bufferImages[current_buffer?0:1][i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize]) bufferImages[current_buffer?0:1][i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize]=max;
                            bufferImages[current_buffer?0:1][i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize]=newMin+(bufferImages[current_buffer?0:1][i+tp*CurrSize->xsize*CurrSize->ysize*CurrSize->zsize]-min)*(newMax-newMin)/(max-min);
                        }
                    }
                }
                else {
                    cout << "ERROR: "<< parser << " is not a valid number"<<endl;
                    i=argc;
                }
                current_buffer=current_buffer?0:1;
            }
	    
            // *********************  ADD  *************************
            else if( strcmp(argv[i], "-div") == 0)
            {
                string parser=argv[++i];
                if(parser.find_first_not_of("1234567890.-+")== string::npos)
                {
                    double divfactor=strtod(parser.c_str(),NULL);
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                        bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i]/divfactor;
                    current_buffer=current_buffer?0:1;
                }
                else
                {
                    nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
                    NewImage->nu=(NewImage->nu>1)?NewImage->nu:1;
                    NewImage->nt=(NewImage->nt>1)?NewImage->nt:1;
                    if(NewImage->datatype!=DT_FLOAT32)
                    {
                        seg_changeDatatype<SegPrecisionTYPE>(NewImage);
                    }
                    SegPrecisionTYPE * NewImagePtr = static_cast<SegPrecisionTYPE *>(NewImage->data);
                    if(NewImage->nx==CurrSize->xsize&&NewImage->ny==CurrSize->ysize&&NewImage->nz==CurrSize->zsize&&NewImage->nt==CurrSize->tsize&&NewImage->nu==CurrSize->usize)
                    {
                        for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                            bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i]/NewImagePtr[i];
                        current_buffer=current_buffer?0:1;
                    }
                    else
                    {
                        cout << "ERROR: Image "<< parser << " is the wrong size  -  original = ( "<<CurrSize->xsize<<","
                             <<CurrSize->ysize<<","<<CurrSize->zsize<<","<<CurrSize->tsize<<","<<CurrSize->usize<<" )  New image = ( "<<NewImage->nx<<","
                            <<NewImage->ny<<","<<NewImage->nz<<","<<NewImage->nt<<","<<NewImage->nu<<" )"<<endl;
                        i=argc;
                    }
                    nifti_image_free(NewImage);
                }
            }
            // *********************  POWER  *************************
            else if(strcmp(argv[i], "-pow") == 0)
            {
                string parser=argv[++i];
                if(((strtod(parser.c_str(),NULL)!=0) || (parser.length()==1 && parser.find("0")!=string::npos)))
                {
                    float factor=strtof(parser.c_str(),NULL);
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                        bufferImages[current_buffer?0:1][i]=powf(bufferImages[current_buffer][i],factor);
                    current_buffer=current_buffer?0:1;
                }
                else
                {
                    cout << "ERROR: "<< parser << " is not a valid number"<<endl;
                    i=argc;
                }
            }
            // *********************  Is NAN  *************************
            else if(strcmp(argv[i], "-isnan") == 0)
            {

                for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                    bufferImages[current_buffer?0:1][i]=isnan(bufferImages[current_buffer][i]);
                current_buffer=current_buffer?0:1;
            }
            // *********************  square_root  *************************
            else if(strcmp(argv[i], "-sqrt") == 0)
            {
                for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                    bufferImages[current_buffer?0:1][i]=sqrtf(bufferImages[current_buffer][i]);
                current_buffer=current_buffer?0:1;
            }
            // *********************  Exponential  *************************
            else if(strcmp(argv[i], "-exp") == 0)
            {
                for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                    bufferImages[current_buffer?0:1][i]=expf(bufferImages[current_buffer][i]);
                current_buffer=current_buffer?0:1;
            }
            // *********************  Exponential  *************************
            else if(strcmp(argv[i], "-log") == 0)
            {
                for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                    bufferImages[current_buffer?0:1][i]=logf(bufferImages[current_buffer][i]);
                current_buffer=current_buffer?0:1;
            }
            // *********************  reciprocal  *************************
            else if(strcmp(argv[i], "-recip") == 0)
            {
                for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                    bufferImages[current_buffer?0:1][i]=1/(bufferImages[current_buffer][i]);
                current_buffer=current_buffer?0:1;
            }
            // *********************  absolute value  *************************
            else if(strcmp(argv[i], "-abs") == 0)
            {
                for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                    bufferImages[current_buffer?0:1][i]=fabs(bufferImages[current_buffer][i]);
                current_buffer=current_buffer?0:1;

            }
            // *********************  bin value  *************************
            else if(strcmp(argv[i], "-bin") == 0)
            {
                for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                    bufferImages[current_buffer?0:1][i]=(bufferImages[current_buffer][i]>0?1.0f:0.0f);
                current_buffer=current_buffer?0:1;
            }
            // *********************  THRESHOLD below  *************************
            else if(strcmp(argv[i], "-thr") == 0)
            {
                string parser=argv[++i];
                if(((strtod(parser.c_str(),NULL)!=0 ) || (parser.length()==1 && parser.find("0")!=string::npos)))
                {
                    double factor=strtod(parser.c_str(),NULL);
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                        bufferImages[current_buffer?0:1][i]=(bufferImages[current_buffer][i]>factor)?bufferImages[current_buffer][i]:0;
                    current_buffer=current_buffer?0:1;
                }
                else
                {
                    cout << "ERROR: "<< parser << " is not a valid number"<<endl;
                    i=argc;
                }
            }
            // *********************  THRESHOLD below  *************************
            else if(strcmp(argv[i], "-equal") == 0)
            {
                string parser=argv[++i];
                if(((strtod(parser.c_str(),NULL)!=0 ) || (parser.length()==1 && parser.find("0")!=string::npos)))
                {
                    double factor=strtod(parser.c_str(),NULL);
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                        bufferImages[current_buffer?0:1][i]=(bufferImages[current_buffer][i]==factor)?1:0;
                    current_buffer=current_buffer?0:1;
                }
                else
                {
                    cout << "ERROR: "<< parser << " is not a valid number"<<endl;
                    i=argc;
                }
            }
            // *********************  Replace below  *************************
            else if(strcmp(argv[i], "-replace") == 0)
            {
                string parser=argv[++i];
                string parser2=argv[++i];
                if(((strtod(parser.c_str(),NULL)!=0 ) || (parser.length()==1 && parser.find("0")!=string::npos)))
                {
                    double factor=strtod(parser.c_str(),NULL);
                    double factor2=strtod(parser2.c_str(),NULL);
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                        bufferImages[current_buffer?0:1][i]=(bufferImages[current_buffer][i]==factor)?factor2:bufferImages[current_buffer][i];
                    current_buffer=current_buffer?0:1;
                }
                else
                {
                    cout << "ERROR: "<< parser << " is not a valid number"<<endl;
                    i=argc;
                }
            }
            // *********************  THRESHOLD ABOVE  *************************
            else if(strcmp(argv[i], "-uthr") == 0)
            {
                string parser=argv[++i];
                if(((strtod(parser.c_str(),NULL)!=0) || (parser.length()==1 && parser.find("0")!=string::npos)))
                {
                    double factor=strtod(parser.c_str(),NULL);
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                        bufferImages[current_buffer?0:1][i]=(bufferImages[current_buffer][i]<factor)?bufferImages[current_buffer][i]:0;
                    current_buffer=current_buffer?0:1;
                }
                else
                {
                    cout << "ERROR: "<< parser << " is not a valid number"<<endl;
                    i=argc;
                }
            }
            // *********************  Dilate   *************************
            else if(strcmp(argv[i], "-dil") == 0)
            {
                string parser=argv[++i];
                if(parser.find_first_not_of("1234567890.-+")== string::npos)
                {
                    double factor=strtod(parser.c_str(),NULL);
                    Dillate(bufferImages[current_buffer],(int)round(factor),CurrSize);
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                        bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i];
                    current_buffer=current_buffer?0:1;
                }
                else
                {
                    cout << "ERROR: "<< parser << " has to be an integer > 0"<<endl;
                    i=argc;
                }
            }
            // *********************  Erosion   *************************
            else if(strcmp(argv[i], "-ero") == 0)
            {
                string parser=argv[++i];
                if(parser.find_first_not_of("1234567890.-+")== string::npos)
                {
                    double factor=strtod(parser.c_str(),NULL);
                    Erosion(bufferImages[current_buffer],(int)round(factor),CurrSize);
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                        bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i];
                    current_buffer=current_buffer?0:1;
                }
                else
                {
                    cout << "ERROR: "<< parser << " has to be an integer > 0"<<endl;
                    i=argc;
                }
            }
//            // *********************  Erosion   *************************
//            else if(strcmp(argv[i], "-eroT") == 0)
//            {
//                string parser=argv[++i];
//                if(parser.find_first_not_of("1234567890.-+")== string::npos)
//                {
//                    double factor=strtod(parser.c_str(),NULL);
//                    //TopologicalErosion(bufferImages[current_buffer],(int)round(factor),CurrSize);
//                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
//                        bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i];
//                    current_buffer=current_buffer?0:1;
//                }
//                else
//                {
//                    cout << "ERROR: "<< parser << " has to be an integer > 0"<<endl;
//                    i=argc;
//                }
//            }
            // *********************  Erosion   *************************
            else if(strcmp(argv[i], "-erot") == 0)
            {
                string parser=argv[++i];
                if(parser.find_first_not_of("1234567890.-+")== string::npos)
                {
                    double factor=strtod(parser.c_str(),NULL);
                    Erosion(bufferImages[current_buffer],(int)round(factor),CurrSize);
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                        bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i];
                    current_buffer=current_buffer?0:1;
                }
                else
                {
                    cout << "ERROR: "<< parser << " has to be an integer > 0"<<endl;
                    i=argc;
                }
            }

            // *********************  Smooth Label   *************************
            else if(strcmp(argv[i], "-smol") == 0)
            {
                string parser=argv[++i];
                if(parser.find_first_not_of("1234567890.-+")== string::npos)
                {
                    double factor=strtod(parser.c_str(),NULL);
                    SmoothLab(bufferImages[current_buffer],factor,CurrSize);
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                        bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i];
                    current_buffer=current_buffer?0:1;
                }
                else
                {
                    cout << "ERROR: "<< parser << " has to be an integer > 0"<<endl;
                    i=argc;
                }
            }
            // *********************  Euclidean Distance Transform   *************************
            else if(strcmp(argv[i], "-euc") == 0)
            {

                bool * Lable= new bool [CurrSize->numel];
                float * Speed= new float [CurrSize->numel];
                for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                {
                    Lable[i]=bufferImages[current_buffer][i];
                    Speed[i]=1.0f;
                }
                float * Distance = DoubleEuclideanDistance_3D(Lable,Speed,CurrSize);

                for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                    bufferImages[current_buffer?0:1][i]=Distance[i];
                current_buffer=current_buffer?0:1;
                delete [] Distance;
                delete [] Lable;
                delete [] Speed;

            }
            // *********************  Geodesic Distance Transform   *************************
            else if(strcmp(argv[i], "-geo") == 0)
            {


                string parser=argv[++i];
                if(parser.find_first_not_of("1234567890.-+")== string::npos)
                {
                    if(strtod(parser.c_str(),NULL)<=0)
                    {
                        cout<< "ERROR: -geo speed should be larger than zero"<<endl;
                        return 1;
                    }
                    bool * Lable= new bool [CurrSize->numel];
                    float * Speed= new float [CurrSize->numel];
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                    {
                        Lable[i]=bufferImages[current_buffer][i];
                        Speed[i]=strtod(parser.c_str(),NULL);
                    }
                    float * Distance = DoubleEuclideanDistance_3D(Lable,Speed,CurrSize);
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                        bufferImages[current_buffer?0:1][i]=Distance[i];
                    current_buffer=current_buffer?0:1;
                    delete [] Distance;
                    delete [] Lable;
                    delete [] Speed;
                }
                else
                {

                    if(   (strtod(parser.c_str(),NULL)!=0 && (parser.find(".nii")==string::npos ||parser.find(".img")==string::npos ||parser.find(".hdr")==string::npos )) ||(parser.length()==1 && parser.find("0")!=string::npos))
                    {
                        cerr<<"ERROR: "<<argv[i]<<"  has to be an image"<<endl;
                        exit(1);
                    }

                    bool * Lable= new bool [CurrSize->numel];
                    float * Speed= new float [CurrSize->numel];
                    nifti_image * SpeedImage=nifti_image_read(parser.c_str(),true);
                    SpeedImage->nu=(SpeedImage->nu>1)?SpeedImage->nu:1;
                    SpeedImage->nt=(SpeedImage->nt>1)?SpeedImage->nt:1;
                    if(SpeedImage->datatype!=DT_FLOAT32)
                    {
                        seg_changeDatatype<float>(SpeedImage);
                    }
                    float * SpeedImagePtr = static_cast<float *>(SpeedImage->data);

                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                    {
                        Lable[i]=bufferImages[current_buffer][i];
                        Speed[i]=SpeedImagePtr[i]>0.0001?SpeedImagePtr[i]:0.0001;
                    }
                    float * Distance = DoubleEuclideanDistance_3D(Lable,Speed,CurrSize);
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                        bufferImages[current_buffer?0:1][i]=Distance[i];
                    current_buffer=current_buffer?0:1;
                    delete [] Distance;
                    delete [] Lable;
                    delete [] Speed;
                    nifti_image_free(SpeedImage);
                }
            }

            // *********************  linear LS Normlise  *************************
            else if(strcmp(argv[i], "-llsnorm") == 0)
            {
                string parser=argv[++i];
                if(   (strtod(parser.c_str(),NULL)!=0 && (parser.find(".nii")==string::npos ||parser.find(".img")==string::npos ||parser.find(".hdr")==string::npos ))
                      ||(parser.length()==1 && parser.find("0")!=string::npos))
                {
                    cerr<<"ERROR: "<<argv[i]<<"  has to be an image"<<endl;
                    exit(1);
                }
                else
                {
                    nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
                    NewImage->nu=(NewImage->nu>1)?NewImage->nu:1;
                    NewImage->nt=(NewImage->nt>1)?NewImage->nt:1;
                    if(NewImage->datatype!=DT_FLOAT32)
                    {
                        seg_changeDatatype<float>(NewImage);
                    }
                    float * NewImagePtr = static_cast<float *>(NewImage->data);

                    // Y=a*X+b
                    float a=0;
                    float b=0;


                    if(NewImage->nx==CurrSize->xsize&&NewImage->ny==CurrSize->ysize&&NewImage->nz==CurrSize->zsize&&NewImage->nt==CurrSize->tsize&&NewImage->nu==CurrSize->usize)
                    {

                        LS_Vecs(bufferImages[current_buffer],NewImagePtr,NULL, (CurrSize->xsize*CurrSize->ysize*CurrSize->zsize),&a, &b);
                        for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                            bufferImages[current_buffer?0:1][i]=a*NewImagePtr[i]+b;
                        current_buffer=current_buffer?0:1;
                    }
                    else
                    {
                        cout << "ERROR: Image "<< parser << " is the wrong size  -  original = ( "<<CurrSize->xsize<<","
                             <<CurrSize->ysize<<","<<CurrSize->zsize<<","<<CurrSize->tsize<<","<<CurrSize->usize<<" ) New image = ( "
                            <<NewImage->nx<<","<<NewImage->ny<<","<<NewImage->nz<<","<<NewImage->nt<<","<<NewImage->nu<<" )"<<endl;
                        i=argc;
                    }
                    nifti_image_free(NewImage);
                }
            }

            // ********************* linear LTS Normlise  *************************
            else if(strcmp(argv[i], "-lltsnorm") == 0)
            {
                string parser=argv[++i];

                string parserout=argv[++i];
                float percent_outlier=strtod(parserout.c_str(),NULL);
                percent_outlier=percent_outlier>0.5?0.5:(percent_outlier<0?0:percent_outlier);
                if(   (strtod(parser.c_str(),NULL)!=0 && (parser.find(".nii")==string::npos ||parser.find(".img")==string::npos ||parser.find(".hdr")==string::npos ))
                      ||(parser.length()==1 && parser.find("0")!=string::npos))
                {
                    cerr<<"ERROR: "<<argv[i]<<"  has to be an image"<<endl;
                    exit(1);
                }
                else
                {
                    nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
                    NewImage->nu=(NewImage->nu>1)?NewImage->nu:1;
                    NewImage->nt=(NewImage->nt>1)?NewImage->nt:1;
                    if(NewImage->datatype!=DT_FLOAT32)
                    {
                        seg_changeDatatype<float>(NewImage);
                    }
                    float * NewImagePtr = static_cast<float *>(NewImage->data);

                    // Y=a*X+b
                    float a=0;
                    float b=0;


                    if(NewImage->nx==CurrSize->xsize&&NewImage->ny==CurrSize->ysize&&NewImage->nz==CurrSize->zsize&&NewImage->nt==CurrSize->tsize&&NewImage->nu==CurrSize->usize)
                    {

                        LTS_Vecs(bufferImages[current_buffer],NewImagePtr,NULL,percent_outlier,20, 0.001, (CurrSize->xsize*CurrSize->ysize*CurrSize->zsize),&a, &b);

                        for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                            bufferImages[current_buffer?0:1][i]=a*NewImagePtr[i]+b;
                        current_buffer=current_buffer?0:1;
                    }
                    else
                    {
                        cout << "ERROR: Image "<< parser << " is the wrong size  -  original = ( "<<CurrSize->xsize<<","
                             <<CurrSize->ysize<<","<<CurrSize->zsize<<","<<CurrSize->tsize<<","<<CurrSize->usize<<" ) New image = ( "
                            <<NewImage->nx<<","<<NewImage->ny<<","<<NewImage->nz<<","<<NewImage->nt<<","<<NewImage->nu<<" )"<<endl;
                        i=argc;
                    }
                    nifti_image_free(NewImage);
                }
            }

            // *********************  QuadraticLS Normlise  *************************
            else if(strcmp(argv[i], "-qlsnorm") == 0)
            {

                string order_str=argv[++i];
                int order=(int)round(strtod(order_str.c_str(),NULL));

                if(order>4){
                    cout << "ERROR: Order is too high... using order 5"<<endl;
                    order=4;
                }
                if(order<1){
                    cout << "ERROR: Order is too low... using order 1"<<endl;
                    order=1;
                }

                string parser=argv[++i];
                nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
                NewImage->nu=(NewImage->nu>1)?NewImage->nu:1;
                NewImage->nt=(NewImage->nt>1)?NewImage->nt:1;
                if(NewImage->datatype!=DT_FLOAT32)
                {
                    seg_changeDatatype<float>(NewImage);
                }
                float * NewImagePtr = static_cast<float *>(NewImage->data);


                const long nvox=CurrSize->xsize*CurrSize->ysize*CurrSize->zsize;

                Eigen::MatrixXf Img1(nvox,order+1);
                Eigen::VectorXf Img2(nvox,1);

                for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize); i++)
                    Img2(i)=bufferImages[current_buffer][i];

                for(int j=0; j<(order+1); j++)
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize); i++)
                        Img1(i,j)=pow(NewImagePtr[i],j);

                Eigen::MatrixXf Img1TransImg1=Img1.transpose()*Img1;
                Eigen::VectorXf Img1TransImg2=Img1.transpose()*Img2;

                Eigen::VectorXf x;
                x=Img1TransImg1.lu().solve(Img1TransImg2); // using a LU factorization

                cout<<x;

                for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize); i++){
                    bufferImages[current_buffer?0:1][i]=x(0);
                }
                for(int j=1; j<(order+1); j++){
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize); i++){
                        bufferImages[current_buffer?0:1][i]+=x(j)*pow(NewImagePtr[i],j);
                    }
                }

                current_buffer=current_buffer?0:1;
                nifti_image_free(NewImage);

            }

            else if(strcmp(argv[i], "-qlsnorm_mask") == 0)
            {

                string order_str=argv[++i];
                int order=(int)round(strtod(order_str.c_str(),NULL));

                if(order>4){
                    cout << "ERROR: Order is too high... using order 5"<<endl;
                    order=4;
                }
                if(order<1){
                    cout << "ERROR: Order is too low... using order 1"<<endl;
                    order=1;
                }

                string parser=argv[++i];
                nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
                NewImage->nu=(NewImage->nu>1)?NewImage->nu:1;
                NewImage->nt=(NewImage->nt>1)?NewImage->nt:1;
                if(NewImage->datatype!=DT_FLOAT32)
                {
                    seg_changeDatatype<float>(NewImage);
                }
                float * NewImagePtr = static_cast<float *>(NewImage->data);

                parser=argv[++i];
                nifti_image * MaskImage=nifti_image_read(parser.c_str(),true);
                MaskImage->nu=(MaskImage->nu>1)?MaskImage->nu:1;
                MaskImage->nt=(MaskImage->nt>1)?MaskImage->nt:1;
                if(MaskImage->datatype!=DT_FLOAT32)
                {
                    seg_changeDatatype<float>(MaskImage);
                }
                float * MaskImagePtr = static_cast<float *>(MaskImage->data);

                size_t nvoxmax=0;
                for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize); i++){
                    if(MaskImagePtr[i]>0 && isnan(bufferImages[current_buffer][i])==0 && isnan(NewImagePtr[i])==0)
                    {
                        nvoxmax++;
                    }
                }

                Eigen::MatrixXf Img1(nvoxmax+1,order+1);
                Eigen::VectorXf Img2(nvoxmax+1);

                size_t nvox=0;
                for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize); i++)
                {
                    if(MaskImagePtr[i]>0 && isnan(bufferImages[current_buffer][i])==0 && isnan(NewImagePtr[i])==0)
                    {
                        Img2(nvox)=bufferImages[current_buffer][i];
                        nvox++;
                    }
                }

                nvox=0;
                for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize); i++){
                    if(MaskImagePtr[i]>0 && isnan(bufferImages[current_buffer][i])==0 && isnan(NewImagePtr[i])==0)
                    {
                        for(int j=0; j<(order+1); j++){
                            Img1(nvox,j)= (j==0)? 1 : pow(NewImagePtr[i],j) ;
                        }
                        nvox++;
                    }
                }
                cout<<nvox<<endl;


                Eigen::MatrixXf Img1TransImg1=Img1.transpose()*Img1;
                Eigen::VectorXf Img1TransImg2=Img1.transpose()*Img2;

                Eigen::VectorXf x;
                x=Img1TransImg1.lu().solve(Img1TransImg2); // using a LU factorization

                cout<<x<<endl;

                cout <<"ui\n"<<x(0)<<endl;
                for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize); i++){
                    bufferImages[current_buffer?0:1][i]=x(0);
                }
                for(int j=1; j<(order+1); j++){
                    cout <<x(j)<<endl;
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize); i++){
                        bufferImages[current_buffer?0:1][i]+=x(j)*pow(NewImagePtr[i],j);
                    }
                }

                current_buffer=current_buffer?0:1;
                nifti_image_free(NewImage);
                nifti_image_free(MaskImage);

            }
            else if(strcmp(argv[i], "-qlsnorm2_mask") == 0)
            {

                string order_str=argv[++i];
                int order=(int)round(strtod(order_str.c_str(),NULL));

                if(order>4){
                    cout << "ERROR: Order is too high... using order 5"<<endl;
                    order=4;
                }
                if(order<1){
                    cout << "ERROR: Order is too low... using order 1"<<endl;
                    order=1;
                }

                string parser=argv[++i];
                nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
                NewImage->nu=(NewImage->nu>1)?NewImage->nu:1;
                NewImage->nt=(NewImage->nt>1)?NewImage->nt:1;
                if(NewImage->datatype!=DT_FLOAT32)
                {
                    seg_changeDatatype<float>(NewImage);
                }
                float * NewImagePtr = static_cast<float *>(NewImage->data);

                parser=argv[++i];
                nifti_image * MaskImage=nifti_image_read(parser.c_str(),true);
                MaskImage->nu=(MaskImage->nu>1)?MaskImage->nu:1;
                MaskImage->nt=(MaskImage->nt>1)?MaskImage->nt:1;
                if(MaskImage->datatype!=DT_FLOAT32)
                {
                    seg_changeDatatype<float>(MaskImage);
                }
                float * MaskImagePtr = static_cast<float *>(MaskImage->data);

                size_t nvoxmax=0;
                for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize); i++){
                    if(MaskImagePtr[i]>0 && isnan(bufferImages[current_buffer][i])==0 && isnan(NewImagePtr[i])==0)
                    {
                        nvoxmax++;
                    }
                }

                Eigen::MatrixXf Img1(nvoxmax+1,order);
                Eigen::VectorXf Img2(nvoxmax+1);

                size_t nvox=0;
                for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize); i++)
                {
                    if(MaskImagePtr[i]>0 && isnan(bufferImages[current_buffer][i])==0 && isnan(NewImagePtr[i])==0)
                    {
                        Img2(nvox)=bufferImages[current_buffer][i];
                        nvox++;
                    }
                }

                nvox=0;
                for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize); i++){
                    if(MaskImagePtr[i]>0 && isnan(bufferImages[current_buffer][i])==0 && isnan(NewImagePtr[i])==0)
                    {
                        for(int j=1; j<(order+1); j++){
                            Img1(nvox,j-1)= pow(NewImagePtr[i],j) ;
                        }
                        nvox++;
                    }
                }


                Eigen::MatrixXf Img1TransImg1=Img1.transpose()*Img1;
                Eigen::VectorXf Img1TransImg2=Img1.transpose()*Img2;

                Eigen::VectorXf x;
                x=Img1TransImg1.lu().solve(Img1TransImg2); // using a LU factorization

                cout<<x<<endl;

                for(int j=1; j<(order+1); j++){
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize); i++){
                        bufferImages[current_buffer?0:1][i]+=x(j-1)*pow(NewImagePtr[i],j);
                    }
                }

                current_buffer=current_buffer?0:1;
                nifti_image_free(NewImage);
                nifti_image_free(MaskImage);

            }

            // *********************  QuadraticLS Normlise  *************************
            else if(strcmp(argv[i], "-qlshnorm") == 0)
            {

                string order_str=argv[++i];
                int order=(int)round(strtod(order_str.c_str(),NULL));

                if(order>5){
                    cout << "ERROR: Order is too high... using order 5"<<endl;
                    order=4;
                }
                if(order<1){
                    cout << "ERROR: Order is too low... using order 1"<<endl;
                    order=1;
                }

                string parser=argv[++i];
                nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
                NewImage->nu=(NewImage->nu>1)?NewImage->nu:1;
                NewImage->nt=(NewImage->nt>1)?NewImage->nt:1;
                if(NewImage->datatype!=DT_FLOAT32)
                {
                    seg_changeDatatype<float>(NewImage);
                }
                float * NewImagePtr = static_cast<float *>(NewImage->data);



                // copy image, sort and fill vector
                size_t img3Dsize=(NewImage->nx*NewImage->ny*NewImage->nz);
                size_t countnan=0;
                for(size_t index=0; index<img3Dsize; index++)
                    countnan+=isnan(NewImagePtr[index])?0:1;
                float * imgsort=new float [countnan];
                size_t countindex=0;
                for(size_t index=0; index<countnan; index++)
                    if(isnan(NewImagePtr[index])==0){
                        imgsort[countindex]=NewImagePtr[index];
                        countindex++;
                    }
                HeapSort(imgsort,countnan-1);
                Eigen::VectorXf Img2(1000,1);
                for(int percentile=0; percentile<1000; percentile++)
                    Img2(percentile)=imgsort[(long)(floor(( (float)(percentile) / 1000.0f ) * (float)( countnan-1 )))];
                delete [] imgsort;


                // copy image, sort and fill vector
                img3Dsize=(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize);
                countnan=0;
                for(size_t index=0; index<img3Dsize; index++)
                    countnan+=isnan(bufferImages[current_buffer][index])?0:1;
                imgsort=new float [countnan];
                countindex=0;
                for(size_t index=0; index<countnan; index++)
                    if(isnan(bufferImages[current_buffer][index])==0){
                        imgsort[countindex]=bufferImages[current_buffer][index];
                        countindex++;
                    }
                HeapSort(imgsort,countnan-1);
                Eigen::MatrixXf Img1(1000,order+1);
                for(int j=0; j<(order+1); j++)
                    for(int percentile=0; percentile<1000; percentile++){
                        Img1(percentile,j)=pow(imgsort[(long)(floor(( (float)(percentile) / 1000.0f ) * (float)( countnan-1 )))] , j );
                    }
                delete [] imgsort;

                Eigen::MatrixXf Img1TransImg1=Img1.transpose()*Img1;
                Eigen::VectorXf Img1TransImg2=Img1.transpose()*Img2;

                Eigen::VectorXf x;
                x=Img1TransImg1.lu().solve(Img1TransImg2); // using a LU factorization

                cout<<x<<endl;

                for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize); i++){
                    bufferImages[current_buffer?0:1][i]=x(0);
                }
                for(int j=1; j<(order+1); j++)
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize); i++){
                        bufferImages[current_buffer?0:1][i]+=x(j)*pow(bufferImages[current_buffer][i],j);
                    }

                current_buffer=current_buffer?0:1;
                nifti_image_free(NewImage);

            }

            else if(strcmp(argv[i], "-qlshnorm_mask") == 0)
            {

                string order_str=argv[++i];
                int order=(int)round(strtod(order_str.c_str(),NULL));

                if(order>4){
                    cout << "ERROR: Order is too high... using order 5"<<endl;
                    order=4;
                }
                if(order<1){
                    cout << "ERROR: Order is too low... using order 1"<<endl;
                    order=1;
                }

                string parser=argv[++i];
                nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
                NewImage->nu=(NewImage->nu>1)?NewImage->nu:1;
                NewImage->nt=(NewImage->nt>1)?NewImage->nt:1;
                if(NewImage->datatype!=DT_FLOAT32)
                {
                    seg_changeDatatype<float>(NewImage);
                }
                float * NewImagePtr = static_cast<float *>(NewImage->data);

                parser=argv[++i];
                nifti_image * MaskImage=nifti_image_read(parser.c_str(),true);
                MaskImage->nu=(MaskImage->nu>1)?MaskImage->nu:1;
                MaskImage->nt=(MaskImage->nt>1)?MaskImage->nt:1;
                if(MaskImage->datatype!=DT_FLOAT32)
                {
                    seg_changeDatatype<float>(MaskImage);
                }
                float * MaskImagePtr = static_cast<float *>(MaskImage->data);




                // copy image, sort and fill vector
                size_t img3Dsize=(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize);
                size_t countnan=0;
                size_t numbsamples=1000;

                for(size_t index=0; index<img3Dsize; index++){
                    if(isnan(NewImagePtr[index])==0&&isnan(bufferImages[current_buffer][index])==0&&MaskImagePtr[index]>0){
                        countnan++;
                    }
                }
                float * imgsort=new float [countnan];
                size_t countindex=0;
                for(size_t index=0; index<img3Dsize; index++){
                    if(isnan(NewImagePtr[index])==0&&isnan(bufferImages[current_buffer][index])==0&&MaskImagePtr[index]>0){
                        imgsort[countindex]=NewImagePtr[index];
                        countindex++;
                    }
                }
                //cout<<countnan<<endl;
                //cout<<countindex<<endl;
                HeapSort(imgsort,countnan-1);
                Eigen::VectorXf Img2(numbsamples,1);
                for(size_t percentile=0; percentile<numbsamples; percentile++){
                    Img2(percentile)=imgsort[(long)(floor(( (float)(percentile) / (float)(numbsamples) ) * (float)( countnan-1 )))];
                    // cout<<percentile<<" - "<<Img2(percentile)<<endl;
                }


                // copy image, sort and fill vector

                countindex=0;
                for(size_t index=0; index<img3Dsize; index++){
                    if(isnan(NewImagePtr[index])==0&&isnan(bufferImages[current_buffer][index])==0&&MaskImagePtr[index]>0){
                        imgsort[countindex]=bufferImages[current_buffer][index];
                        countindex++;
                    }
                }
                //cout<<countnan<<endl;
                //cout<<countindex<<endl;
                HeapSort(imgsort,countnan-1);
                Eigen::MatrixXf Img1(numbsamples,order+1);
                for(size_t percentile=0; percentile<numbsamples; percentile++){
                    for(int j=0; j<(order+1); j++){
                        Img1(percentile,j)=pow(imgsort[(long)(floor(( (float)(percentile) / (float)(numbsamples) ) * (float)( countnan-1 )))] , j );

                    }
                    // cout<<percentile<<" - "<<Img1(percentile,1)<<endl;
                }
                delete [] imgsort;




                Eigen::MatrixXf Img1TransImg1=Img1.transpose()*Img1;
                Eigen::VectorXf Img1TransImg2=Img1.transpose()*Img2;

                Eigen::VectorXf x;
                x=Img1TransImg1.lu().solve(Img1TransImg2); // using a LU factorization

                cout<<x<<endl;
                for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize); i++){
                    bufferImages[current_buffer?0:1][i]=x(0);
                }
                for(int j=1; j<(order+1); j++){
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize); i++){
                        bufferImages[current_buffer?0:1][i]+=x(j)*pow(bufferImages[current_buffer][i],j);
                    }
                }

                current_buffer=current_buffer?0:1;
                nifti_image_free(NewImage);
                nifti_image_free(MaskImage);

            }
            // *********************  GAUSSIAN SMOTHING *************************
            else if(strcmp(argv[i], "-smo") == 0)
            {
                string parser=argv[++i];
                if((strtod(parser.c_str(),NULL)!=0 ))
                {
                    float factor=strtof(parser.c_str(),NULL);
                    for(long tp=0; tp<(long)(CurrSize->tsize*CurrSize->usize); tp++){
                        //create dummy nii
                        nifti_image * TMPnii = nifti_copy_nim_info(InputImage);
                        TMPnii->dim[1]=CurrSize->xsize;
                        TMPnii->dim[2]=CurrSize->ysize;
                        TMPnii->dim[3]=CurrSize->zsize;
                        TMPnii->dim[4]=TMPnii->nt=1;
                        TMPnii->dim[5]=TMPnii->nu=1;
                        nifti_update_dims_from_array(TMPnii);
                        //copy pointer, run gaussian, and set to null
                        TMPnii->data=static_cast<void*>(&bufferImages[current_buffer][CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*tp]);
                        GaussianSmoothing5D_nifti(TMPnii,NULL,factor);
                        TMPnii->data=NULL;
                        //As TMPnii->data=NULL, the free will not cause any harm
                        nifti_image_free(TMPnii);
                    }


                    //current_buffer=current_buffer?0:1;
                }
                else
                {
                    cout << "ERROR: "<< parser << " has to be a number > 0"<<endl;
                    i=argc;
                }
            }
            // *********************  GAUSSIAN SMOTHING *************************
            else if(strcmp(argv[i], "-smoNaN") == 0)
            {
                string filename=argv[++i];
                nifti_image * MaskImage=nifti_image_read(filename.c_str(),true);
                MaskImage->nu=(MaskImage->nu>1)?MaskImage->nu:1;
                MaskImage->nt=(MaskImage->nt>1)?MaskImage->nt:1;
                if(MaskImage->datatype!=DT_FLOAT32)
                {
                    seg_changeDatatype<SegPrecisionTYPE>(MaskImage);
                }


                for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                    bufferImages[current_buffer][i]=(bufferImages[current_buffer][i])?
                                bufferImages[current_buffer][i]:
                                std::numeric_limits<float>::quiet_NaN();

                for(long tp=0; tp<(long)(CurrSize->tsize*CurrSize->usize); tp++){


                    //create dummy nii
                    nifti_image * TMPnii = nifti_copy_nim_info(InputImage);
                    TMPnii->dim[1]=CurrSize->xsize;
                    TMPnii->dim[2]=CurrSize->ysize;
                    TMPnii->dim[3]=CurrSize->zsize;
                    TMPnii->dim[4]=TMPnii->nt=1;
                    TMPnii->dim[5]=TMPnii->nu=1;
                    nifti_update_dims_from_array(TMPnii);
                    //copy pointer, run gaussian, and set to null
                    TMPnii->data=static_cast<void*>(&bufferImages[current_buffer][CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*tp]);
                    GaussianSmoothing4D_Nan_nifti(TMPnii,MaskImage);
                    TMPnii->data=NULL;
                    //As TMPnii->data=NULL, the free will not cause any harm
                    nifti_image_free(TMPnii);

                    nifti_image_free(MaskImage);
                }

                    //current_buffer=current_buffer?0:1;

            }
            // *********************  GAUSSIAN sharpening  (NOT WORKING) *************************
            else if(strcmp(argv[i], "-sharp") == 0)
            {
                string parser=argv[++i];
                if((strtod(parser.c_str(),NULL)!=0 ))
                {
                    float factor=strtof(parser.c_str(),NULL);
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                        bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i];

                    GaussianFilter4D_cArray(&bufferImages[current_buffer][0], factor, CurrSize);
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                        bufferImages[current_buffer?0:1][i]=(bufferImages[current_buffer?0:1][i]-bufferImages[current_buffer][i]);

                    current_buffer=current_buffer?0:1;
                }
                else
                {
                    cout << "ERROR: "<< parser << " has to be a number > 0"<<endl;
                    i=argc;
                }
            }
            // *********************  Min  *************************
            else if(strcmp(argv[i], "-min") == 0)
            {
                string parser=argv[++i];
                if(!(parser.find_first_not_of("1234567890.-+")== string::npos))
                {

                    nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
                    NewImage->nu=(NewImage->nu>1)?NewImage->nu:1;
                    NewImage->nt=(NewImage->nt>1)?NewImage->nt:1;
                    if(NewImage->datatype!=DT_FLOAT32)
                    {
                        seg_changeDatatype<SegPrecisionTYPE>(NewImage);
                    }
                    SegPrecisionTYPE * NewImagePtr = static_cast<SegPrecisionTYPE *>(NewImage->data);
                    if(NewImage->nx==CurrSize->xsize&&NewImage->ny==CurrSize->ysize&&NewImage->nz==CurrSize->zsize&&NewImage->nt==CurrSize->tsize&&NewImage->nu==CurrSize->usize)
                    {
                        for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                            bufferImages[current_buffer?0:1][i]=min(bufferImages[current_buffer][i],NewImagePtr[i]);
                        current_buffer=current_buffer?0:1;
                    }
                }
            }

            // *********************  Otsu thresholding *************************
            else if(strcmp(argv[i], "-otsu") == 0)
            {

                otsu(bufferImages[current_buffer],NULL,CurrSize);
                for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                    bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i];

                current_buffer=current_buffer?0:1;
            }

            // *********************  Fill  *************************
            else if(strcmp(argv[i], "-fill") == 0)
            {
                if(CurrSize->tsize==1)
                {
                    Close_Forground_ConnectComp<float,float>(static_cast<void*>(bufferImages[current_buffer]),static_cast<void*>(bufferImages[current_buffer?0:1]),CurrSize);
                    current_buffer=current_buffer?0:1;
                }
                else
                {
                    cout << "ERROR: Image to -fill is not 3D"<<endl;
                    i=argc;
                }
            }
            // *********************  Largest Connected Component  *************************
            else if(strcmp(argv[i], "-lconcomp") == 0)
            {
                if(CurrSize->tsize==1)
                {
                    Largest_ConnectComp<float,float>(static_cast<void*>(bufferImages[current_buffer]),static_cast<void*>(bufferImages[current_buffer?0:1]),CurrSize);
                    current_buffer=current_buffer?0:1;
                }
                else
                {
                    cout << "ERROR: Image to -lconcomp is not 3D"<<endl;
                    i=argc;
                }
            }
            // *********************  Connected Components 6NN  *************************
            else if(strcmp(argv[i], "-concomp6") == 0)
            {
                if(CurrSize->tsize==1)
                {
                    ConnectComp6NN<float,float>(static_cast<void*>(bufferImages[current_buffer]),static_cast<void*>(bufferImages[current_buffer?0:1]),CurrSize);
                    current_buffer=current_buffer?0:1;
                }
                else
                {
                    cout << "ERROR: Image to -concomp6 is not 3D"<<endl;
                    i=argc;
                }
            }

            // *********************  Connected Components 6NN  *************************
            else if(strcmp(argv[i], "-concomp26") == 0)
            {
                if(CurrSize->tsize==1)
                {
                    ConnectComp26NN<float,float>(static_cast<void*>(bufferImages[current_buffer]),static_cast<void*>(bufferImages[current_buffer?0:1]),CurrSize);
                    current_buffer=current_buffer?0:1;
                }
                else
                {
                    cout << "ERROR: Image to -concomp26 is not 3D"<<endl;
                    i=argc;
                }
            }

            // *********************  Range  *************************
            else if(strcmp(argv[i], "-range") == 0)
            {
                float min=FLT_MAX;
                float max=-FLT_MAX;
                for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                {
                    max=bufferImages[current_buffer][i]>max?bufferImages[current_buffer][i]:max;
                    min=bufferImages[current_buffer][i]<min?bufferImages[current_buffer][i]:min;
                }
                InputImage->cal_max=max;
                InputImage->cal_min=min;
            }
            // *********************  Extract time point  *************************
            else if(strcmp(argv[i], "-tp") == 0)
            {
                string parser=argv[++i];
                if(((strtod(parser.c_str(),NULL)!=0) || (parser.length()==1 && parser.find("0")!=string::npos && parser.find("0")!=string::npos) )&& strtod(parser.c_str(),NULL)<=CurrSize->tsize )
                {
                    float factor=strtof(parser.c_str(),NULL);
                    InputImage->dim[4]=InputImage->nt=CurrSize->tsize=1;
                    InputImage->dim[0]=3;
                    InputImage->dim[5]=InputImage->nu=CurrSize->usize=1;
                    for(long i=0; i<CurrSize->numel; i++)
                        bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i+(int)round(factor)*CurrSize->numel];

                    current_buffer=current_buffer?0:1;
                }
                else
                {
                    cout << "ERROR: "<< parser << " is not an integer"<<endl;
                    i=argc;
                }
            }

            // *********************  Split Lables  *************************
            else if(strcmp(argv[i], "-splitlab") == 0)
            {
                int maxlab=0;
                for(long index=0; index<(CurrSize->numel*(CurrSize->tsize*CurrSize->usize)); index++)
                    maxlab=(round(bufferImages[current_buffer][index])>maxlab)?(int)round(bufferImages[current_buffer][index]):maxlab;
                maxlab=maxlab+1;
                if(maxlab>0 && CurrSize->tsize<=1&& CurrSize->usize<=1)
                {
                    CurrSize->tsize=maxlab;
                    CurrSize->usize=1;

                    delete [] bufferImages[current_buffer?0:1];
                    bufferImages[current_buffer?0:1]= new SegPrecisionTYPE [CurrSize->numel*maxlab];
                    for(long index=0; index<(CurrSize->numel*maxlab); index++)
                        bufferImages[current_buffer?0:1][index]=0.0f;
                    for(long index=0; index<(CurrSize->numel); index++)
                        bufferImages[current_buffer?0:1][index+(int)round(bufferImages[current_buffer][index])*CurrSize->numel]=1.0f;
                    delete [] bufferImages[current_buffer];
                    bufferImages[current_buffer]= new SegPrecisionTYPE [CurrSize->numel*maxlab];
                    for(long index=0; index<(CurrSize->numel*maxlab); index++)
                        bufferImages[current_buffer][index]=0;
                    current_buffer=current_buffer?0:1;
                }
                else
                {
                    if(CurrSize->tsize<=1&& CurrSize->usize<=1)
                    {
                        cout << "ERROR: Working image is not 3D"<<endl;
                    }
                    else
                    {
                        cout << "ERROR: Found only "<< maxlab << " labels"<<endl;
                    }
                    i=argc;
                }
            }
            // *********************  Split Lables  *************************
            else if(strcmp(argv[i], "-splitinter") == 0)
            {
                string direction=argv[++i];
                if(CurrSize->tsize<=1&& CurrSize->usize<=1){
                    CurrSize->tsize=2;
                    CurrSize->usize=1;
                    int oldxsize=CurrSize->xsize;
                    int oldysize=CurrSize->ysize;
                    int xincrement=1;
                    int yincrement=1;
                    int zincrement=1;

                    if(direction==string("x")){
                        CurrSize->xsize=round(CurrSize->xsize/2);
                        xincrement=2;
                        Scalling[0]= 0.5f;
                    }
                    else if(direction==string("y")){
                        CurrSize->ysize=round(CurrSize->ysize/2);
                        yincrement=2;
                        Scalling[1]= 0.5f;
                    }
                    else if(direction==string("z")){
                        CurrSize->zsize=round(CurrSize->zsize/2);
                        zincrement=2;
                        Scalling[2]= 0.5f;
                    }
                    else{
                        cout << "ERROR: Direction "<< direction << " is not x, y or z"<<endl;
                        exit(1);
                    }

                    CurrSize->numel=(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize);
                    for(long indexZ=0, indexZold=0; indexZ<CurrSize->zsize; indexZ++, indexZold+=zincrement){
                        for(long indexY=0, indexYold=0; indexY<CurrSize->ysize; indexY++, indexYold+=yincrement){
                            for(long indexX=0, indexXold=0; indexX<CurrSize->xsize; indexX++, indexXold+=xincrement){
                                bufferImages[current_buffer?0:1][indexX+indexY*CurrSize->xsize+indexZ*CurrSize->ysize*CurrSize->xsize]=
                                        bufferImages[current_buffer][indexXold+indexYold*oldxsize+indexZold*oldysize*oldxsize];
                                bufferImages[current_buffer?0:1][indexX+indexY*CurrSize->xsize+indexZ*CurrSize->ysize*CurrSize->xsize+CurrSize->numel]=
                                        bufferImages[current_buffer][(indexXold+xincrement-1)+(indexYold+yincrement-1)*oldxsize+(indexZold+zincrement-1)*oldysize*oldxsize];
                            }
                        }
                    }
                    current_buffer=current_buffer?0:1;

                }

            }
            // *********************  Split Lables  *************************
            else if(strcmp(argv[i], "-splitnorm") == 0)
            {
                string direction=argv[++i];
                if(CurrSize->tsize<=1&& CurrSize->usize<=1){
                    CurrSize->tsize=1;
                    CurrSize->usize=1;
                    int xincrement=1;
                    int yincrement=1;
                    int zincrement=1;
                    //bool isdirectionsizeodd=0;
                    if(direction==string("x") || direction==string("1")){
                        xincrement=2;
                        //isdirectionsizeodd=(CurrSize->xsize%2)==0;
                    }
                    else if(direction==string("y") || direction==string("2")){
                        yincrement=2;
                        //isdirectionsizeodd=(CurrSize->ysize%2)==0;
                    }
                    else if(direction==string("z") || direction==string("3")){
                        zincrement=2;
                        //isdirectionsizeodd=(CurrSize->zsize%2)==0;
                    }
                    else{
                        cout << "ERROR: Direction "<< direction << " is not x, y or z"<<endl;
                        exit(1);
                    }

                    //double regul=5.0f;
                    std::vector<float> sortedimg;
                    for(long indexZ=0; indexZ<(CurrSize->zsize); indexZ+=zincrement){
                        for(long indexY=0; indexY<(CurrSize->ysize); indexY+=yincrement){
                            for(long indexX=0; indexX<(CurrSize->xsize); indexX+=xincrement){
                                sortedimg.push_back(bufferImages[current_buffer][indexX+indexY*CurrSize->xsize+indexZ*CurrSize->ysize*CurrSize->xsize]);
                            }
                        }
                    }
                    std::sort(sortedimg.begin(), sortedimg.end());
                    float thresh=0.5f*sortedimg.at(round(sortedimg.size()*0.5f)); // Find a rough background threshold to ignore non-brain tissues
                    sortedimg.clear();

                    std::vector<float> sortedvec;
                    CurrSize->numel=(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize);
                    for(long indexZ=(zincrement-1); indexZ<(CurrSize->zsize-(zincrement-1)); indexZ++){
                        for(long indexY=(yincrement-1); indexY<(CurrSize->ysize-(yincrement-1)); indexY++){
                            for(long indexX=(xincrement-1); indexX<(CurrSize->xsize-(xincrement-1)); indexX++){

                                double previous_next_mean_val=(bufferImages[current_buffer][(indexX+xincrement-1)+(indexY+yincrement-1)*CurrSize->xsize+(indexZ+zincrement-1)*CurrSize->ysize*CurrSize->xsize]+
                                        bufferImages[current_buffer][(indexX-xincrement+1)+(indexY-yincrement+1)*CurrSize->xsize+(indexZ-zincrement+1)*CurrSize->ysize*CurrSize->xsize])/(2.0f);
                                double current_val=bufferImages[current_buffer][indexX+indexY*CurrSize->xsize+indexZ*CurrSize->ysize*CurrSize->xsize];

                                bool oddeven=(xincrement>1?indexX%2==0:(yincrement>1?indexY%2==0:(zincrement>1?indexZ%2==0:0)));

                                float curRat=(current_val*previous_next_mean_val)/(oddeven?(current_val*current_val):(previous_next_mean_val*previous_next_mean_val));
                                if(!(curRat!=curRat) && current_val>thresh){
                                    sortedvec.push_back(curRat);
                                }
                            }
                        }
                    }
                    std::sort(sortedvec.begin(), sortedvec.end());
                    double compensation_ratio=sortedvec.at(round(sortedvec.size()/2.0f)); // Get the median ratio
                    cout<<compensation_ratio<<endl;

                    sortedvec.clear();
                    for(long indexZ=0; indexZ<(CurrSize->zsize); indexZ+=zincrement){
                        for(long indexY=0; indexY<(CurrSize->ysize); indexY+=yincrement){
                            for(long indexX=0; indexX<(CurrSize->xsize); indexX+=xincrement){
                                bufferImages[current_buffer?0:1][indexX+indexY*CurrSize->xsize+indexZ*CurrSize->ysize*CurrSize->xsize]=
                                        compensation_ratio*bufferImages[current_buffer][indexX+indexY*CurrSize->xsize+indexZ*CurrSize->ysize*CurrSize->xsize];
                                if((indexX+xincrement-1)<CurrSize->xsize && (indexY+yincrement-1)<CurrSize->ysize && (indexZ+zincrement-1)<CurrSize->zsize)
                                {
                                    bufferImages[current_buffer?0:1][(indexX+xincrement-1)+(indexY+yincrement-1)*CurrSize->xsize+(indexZ+zincrement-1)*CurrSize->ysize*CurrSize->xsize]=
                                            bufferImages[current_buffer][(indexX+xincrement-1)+(indexY+yincrement-1)*CurrSize->xsize+(indexZ+zincrement-1)*CurrSize->ysize*CurrSize->xsize];
                                }
                            }
                        }
                    }
                    current_buffer=current_buffer?0:1;

                }

            }
//            // *********************  Split Lables  *************************
//            else if(strcmp(argv[i], "-splitnorm2") == 0)
//            {
//                string direction=argv[++i];
//                if(CurrSize->tsize<=1&& CurrSize->usize<=1){
//                    int cur_dims[8]={3,CurrSize->xsize,CurrSize->ysize,CurrSize->zsize,1,1,1,1};
//                    nifti_image * NewImage1=nifti_copy_nim_info(InputImage);
//                    NewImage1->data= (void *) calloc(InputImage->nvox, sizeof(float));

//                    nifti_image * NewImage2=nifti_copy_nim_info(InputImage);
//                    NewImage2->data=(void *) calloc(InputImage->nvox, sizeof(float));

//                    float* NewImage1_ptr=static_cast<SegPrecisionTYPE *>(NewImage1->data);
//                    float* NewImage2_ptr=static_cast<SegPrecisionTYPE *>(NewImage2->data);

//                    int xincrement=1;
//                    int yincrement=1;
//                    int zincrement=1;

//                    if(direction==string("x")){
//                        xincrement=2;
//                    }
//                    else if(direction==string("y")){
//                        yincrement=2;
//                    }
//                    else if(direction==string("z")){
//                        zincrement=2;
//                    }
//                    else{
//                        cout << "ERROR: Direction "<< direction << " is not x, y or z"<<endl;
//                        exit(1);
//                    }
////                    double regul=1.0e-15f;

//                    CurrSize->numel=(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize);
//                    for(long indexZ=(zincrement-1); indexZ<(CurrSize->zsize-(zincrement-1)); indexZ++){
//                        for(long indexY=(yincrement-1); indexY<(CurrSize->ysize-(yincrement-1)); indexY++){
//                            for(long indexX=(xincrement-1); indexX<(CurrSize->xsize-(xincrement-1)); indexX++){

//                                if(xincrement>1?indexX%2==0:(yincrement>1?indexY%2==0:(zincrement>1?indexZ%2==0:0))){
//                                    NewImage1_ptr[indexX+indexY*CurrSize->xsize+indexZ*CurrSize->ysize*CurrSize->xsize]=
//                                            (bufferImages[current_buffer][(indexX+xincrement-1)+(indexY+yincrement-1)*CurrSize->xsize+(indexZ+zincrement-1)*CurrSize->ysize*CurrSize->xsize]+
//                                            bufferImages[current_buffer][(indexX-xincrement+1)+(indexY-yincrement+1)*CurrSize->xsize+(indexZ-zincrement+1)*CurrSize->ysize*CurrSize->xsize])/(2.0f);
//                                    NewImage2_ptr[indexX+indexY*CurrSize->xsize+indexZ*CurrSize->ysize*CurrSize->xsize]=
//                                            bufferImages[current_buffer][indexX+indexY*CurrSize->xsize+indexZ*CurrSize->ysize*CurrSize->xsize];
//                                }
//                                else{
//                                    NewImage2_ptr[indexX+indexY*CurrSize->xsize+indexZ*CurrSize->ysize*CurrSize->xsize]=
//                                            (bufferImages[current_buffer][(indexX+xincrement-1)+(indexY+yincrement-1)*CurrSize->xsize+(indexZ+zincrement-1)*CurrSize->ysize*CurrSize->xsize]+
//                                            bufferImages[current_buffer][(indexX-xincrement+1)+(indexY-yincrement+1)*CurrSize->xsize+(indexZ-zincrement+1)*CurrSize->ysize*CurrSize->xsize])/(2.0f);
//                                    NewImage1_ptr[indexX+indexY*CurrSize->xsize+indexZ*CurrSize->ysize*CurrSize->xsize]=
//                                            bufferImages[current_buffer][indexX+indexY*CurrSize->xsize+indexZ*CurrSize->ysize*CurrSize->xsize];
//                                }

//                            }
//                        }
//                    }
//                    nifti_set_filenames(NewImage1,"img1.nii.gz",0,0);
//                    nifti_image_write(NewImage1);
//                    nifti_set_filenames(NewImage2,"img2.nii.gz",0,0);
//                    nifti_image_write(NewImage2);
//                    // Y=a*X+b
//                    float a=0;
//                    float b=0;
//                    LS_Vecs(NewImage1_ptr,NewImage2_ptr,NULL, (CurrSize->xsize*CurrSize->ysize*CurrSize->zsize),&a, &b);
//                    cout<<a<<"  "<<b<<endl;

//                    for(long indexZ=0; indexZ<(CurrSize->zsize); indexZ+=zincrement){
//                        for(long indexY=0; indexY<(CurrSize->ysize); indexY+=yincrement){
//                            for(long indexX=0; indexX<(CurrSize->xsize); indexX+=xincrement){
//                                bufferImages[current_buffer?0:1][indexX+indexY*CurrSize->xsize+indexZ*CurrSize->ysize*CurrSize->xsize]=
//                                        a*bufferImages[current_buffer][indexX+indexY*CurrSize->xsize+indexZ*CurrSize->ysize*CurrSize->xsize];
//                                if((indexX+xincrement-1)<CurrSize->xsize && (indexY+yincrement-1)<CurrSize->ysize && (indexZ+zincrement-1)<CurrSize->zsize)
//                                {
//                                    bufferImages[current_buffer?0:1][(indexX+xincrement-1)+(indexY+yincrement-1)*CurrSize->xsize+(indexZ+zincrement-1)*CurrSize->ysize*CurrSize->xsize]=
//                                            bufferImages[current_buffer][(indexX+xincrement-1)+(indexY+yincrement-1)*CurrSize->xsize+(indexZ+zincrement-1)*CurrSize->ysize*CurrSize->xsize];
//                                }
//                            }
//                        }
//                    }
//                    current_buffer=current_buffer?0:1;
//                }
//            }
            // *********************  Split Lables  *************************
            else if(strcmp(argv[i], "-joininter") == 0)
            {
                string direction=argv[++i];
                if(CurrSize->tsize==2&& CurrSize->usize<=1){
                    CurrSize->tsize=1;
                    CurrSize->usize=1;
                    int oldxsize=CurrSize->xsize;
                    int oldysize=CurrSize->ysize;
                    int oldzsize=CurrSize->zsize;
                    int xincrement=1;
                    int yincrement=1;
                    int zincrement=1;

                    if(direction==string("x")){
                        CurrSize->xsize=round(CurrSize->xsize*2);
                        xincrement=2;
                        Scalling[0]= 2.0f;
                    }
                    else if(direction==string("y")){
                        CurrSize->ysize=round(CurrSize->ysize*2);
                        yincrement=2;
                        Scalling[1]= 2.0f;
                    }
                    else if(direction==string("z")){
                        CurrSize->zsize=round(CurrSize->zsize*2);
                        zincrement=2;
                        Scalling[2]= 2.0f;
                    }
                    else{
                        cout << "ERROR: Direction "<< direction << " is not x, y or z"<<endl;
                        exit(1);
                    }

                    CurrSize->numel=(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize);
                    long oldnumel=(long)(oldxsize*oldysize*oldzsize);
                    for(long indexZ=0, indexZold=0; indexZ<CurrSize->zsize; indexZ+=zincrement, indexZold++){
                        for(long indexY=0, indexYold=0; indexY<CurrSize->ysize; indexY+=yincrement, indexYold++){
                            for(long indexX=0, indexXold=0; indexX<CurrSize->xsize; indexX+=xincrement, indexXold++){
                                bufferImages[current_buffer?0:1][(indexX)+(indexY)*CurrSize->xsize+(indexZ)*CurrSize->ysize*CurrSize->xsize]=
                                        bufferImages[current_buffer][indexXold+indexYold*oldxsize+indexZold*oldysize*oldxsize];
                                bufferImages[current_buffer?0:1][(indexX+xincrement-1)+(indexY+yincrement-1)*CurrSize->xsize+(indexZ+zincrement-1)*CurrSize->ysize*CurrSize->xsize]=
                                        bufferImages[current_buffer][indexXold+indexYold*oldxsize+indexZold*oldysize*oldxsize+oldnumel];
                            }
                        }
                    }
                    current_buffer=current_buffer?0:1;

                }
                else{
                    cout << "ERROR: Number of time points is not 2"<<endl;
                    exit(1);
                }

            }
            // *********************  merge time points  *************************
            else if(strcmp(argv[i], "-merge") == 0)
            {
                string parser=argv[++i];
                string parsertp=argv[++i];
                if(strtod(parser.c_str(),NULL) && (strtod(parser.c_str(),NULL)!=0 ))
                {
                    long numberof_new_images=(int)strtof(parser.c_str(),NULL);
                    long dim=(int)strtof(parsertp.c_str(),NULL);

                    long old_tsize=CurrSize->tsize;
                    long old_usize=CurrSize->usize;

                    long new_tsize=CurrSize->tsize;
                    long new_usize=CurrSize->usize;
                    if(dim==4)
                    {
                        new_tsize=CurrSize->tsize+(int)numberof_new_images;
                    }
                    else if(dim==5)
                    {
                        new_usize=CurrSize->usize+(int)numberof_new_images;
                    }
                    else{
                        cout<< "ERROR: dim has to be 4 or 5"<<endl;
                        return 1;
                    }


                    delete [] bufferImages[current_buffer?0:1];
                    bufferImages[current_buffer?0:1]= new SegPrecisionTYPE [CurrSize->numel*(new_tsize*new_usize)];

                    for(long index=0; index<(CurrSize->numel*(old_tsize*old_usize)); index++)
                        bufferImages[current_buffer?0:1][index]=bufferImages[current_buffer][index];

                    delete [] bufferImages[current_buffer];
                    bufferImages[current_buffer]= new SegPrecisionTYPE [CurrSize->numel*(new_tsize*new_usize)];

                    for(long index=0; index<(CurrSize->numel*(old_tsize*old_usize)); index++)
                        bufferImages[current_buffer][index]=bufferImages[current_buffer?0:1][index];

                    current_buffer=current_buffer?0:1;

                    CurrSize->usize=new_usize;
                    CurrSize->tsize=new_tsize;

                    for(long tp=0; tp<(long)numberof_new_images; tp++)
                    {
                        string parser_image_name=argv[++i];
                        if(parser_image_name.find(string(".nii"))>0 || parser_image_name.find(string(".img")) ||parser_image_name.find(string(".hdr"))>0)
                        {
                            nifti_image * NewImage=nifti_image_read(parser_image_name.c_str(),true);
                            if(NewImage == NULL)
                            {
                                cout<< "ERROR: When reading the image"<<parser_image_name<<endl;
                                return 1;
                            }
                            if(dim==4){
                                if(NewImage->nx==InputImage->nx&&NewImage->ny==InputImage->ny&&NewImage->nz==InputImage->nz && NewImage->nt<=1)
                                {
                                    if(NewImage->datatype!=DT_FLOAT32)
                                    {
                                        seg_changeDatatype<SegPrecisionTYPE>(NewImage);
                                    }
                                    SegPrecisionTYPE * NewImagePtr = static_cast<SegPrecisionTYPE *>(NewImage->data);
                                    for(long index=0; index<(long)CurrSize->numel; index++)
                                        bufferImages[current_buffer?0:1][index+(old_tsize+tp)*CurrSize->numel]=NewImagePtr[index];
                                }
                                else
                                {
                                    cout<< "ERROR: Image "<<parser_image_name<<" [nx,ny,nz] do not match or nt>1"<<endl;
                                    return 1;
                                }
                            }
                            else if(dim==5){
                                if(NewImage->nx==InputImage->nx&&NewImage->ny==InputImage->ny&&NewImage->nz==InputImage->nz&&NewImage->nt==InputImage->nt && NewImage->nu<=1)
                                {
                                    if(NewImage->datatype!=DT_FLOAT32)
                                    {
                                        seg_changeDatatype<SegPrecisionTYPE>(NewImage);
                                    }
                                    SegPrecisionTYPE * NewImagePtr = static_cast<SegPrecisionTYPE *>(NewImage->data);
                                    for(long index=0; index<(long)CurrSize->numel*old_tsize; index++)
                                        bufferImages[current_buffer?0:1][index+(old_tsize+old_tsize*tp)*CurrSize->numel]=NewImagePtr[index];
                                }
                                else
                                {
                                    cout<< "ERROR: Image "<<parser_image_name<<" [nx,ny,nz,nt] do not match or nu>1"<<endl;
                                    return 1;
                                }
                            }

                        }
                    }
                    current_buffer=current_buffer?0:1;
                }
                else
                {
                    cout << "ERROR: "<< parser << " has to be an integer > 0"<<endl;
                    i=argc;
                }
            }
            // *********************  merge time points  *************************
            else if(strcmp(argv[i], "-subsamp2") == 0)
            {

                int newx=(int)floor(CurrSize->xsize/2.0f);
                int newy=(int)floor(CurrSize->ysize/2.0f);
                int newz=(int)floor(CurrSize->zsize/2.0f);
                int newnumel=newx*newy*newz;


                for(long tp=0; tp<(long)(CurrSize->tsize*CurrSize->usize); tp++){
                    //create dummy nii
                    nifti_image * TMPnii = nifti_copy_nim_info(InputImage);
                    TMPnii->dim[1]=CurrSize->xsize;
                    TMPnii->dim[2]=CurrSize->ysize;
                    TMPnii->dim[3]=CurrSize->zsize;
                    TMPnii->dim[4]=TMPnii->nt=1;
                    TMPnii->dim[5]=TMPnii->nu=1;
                    nifti_update_dims_from_array(TMPnii);
                    //copy pointer, run gaussian, and set to null
                    TMPnii->data=static_cast<void*>(&bufferImages[current_buffer][CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*tp]);
                    GaussianSmoothing5D_nifti(TMPnii,NULL,0.5);
                    TMPnii->data=NULL;
                    //As TMPnii->data=NULL, the free will not cause any harm
                    nifti_image_free(TMPnii);
                }

                delete [] bufferImages[current_buffer?0:1];
                bufferImages[current_buffer?0:1]= new SegPrecisionTYPE [newnumel*CurrSize->tsize];
                Scalling[0]=0.5;
                Scalling[1]=0.5;
                Scalling[2]=0.5;

                for(long indexT=0; indexT<CurrSize->tsize; indexT++)
                    for(long indexZ=0; indexZ<newz; indexZ++)
                        for(long indexY=0; indexY<newy; indexY++)
                            for(long indexX=0; indexX<newx; indexX++)
                                bufferImages[current_buffer?0:1][indexX+indexY*newx+indexZ*newy*newx+indexT*newnumel]=bufferImages[current_buffer][indexX*2+indexY*2*CurrSize->xsize+indexZ*2*CurrSize->xsize*CurrSize->ysize+indexT*CurrSize->numel];




                delete [] bufferImages[current_buffer];
                bufferImages[current_buffer]= new SegPrecisionTYPE [newnumel*CurrSize->tsize];
                current_buffer=current_buffer?0:1;
                CurrSize->xsize=newx;
                CurrSize->ysize=newy;
                CurrSize->zsize=newz;
                CurrSize->numel=newnumel;
            }
            // *********************  merge time points  *************************
            else if(strcmp(argv[i], "-subsamp2xy") == 0)
            {

                int newx=(int)floor(static_cast<float>(CurrSize->xsize)/2.0f);
                int newy=(int)floor(static_cast<float>(CurrSize->ysize)/2.0f);
                int newz=(int)floor(static_cast<float>(CurrSize->zsize));
                int newnumel=newx*newy*newz;


                delete [] bufferImages[current_buffer?0:1];
                bufferImages[current_buffer?0:1]= new SegPrecisionTYPE [newnumel*CurrSize->tsize];
                Scalling[0]=0.5;
                Scalling[1]=0.5;
                Scalling[2]=1;

                for(long indexT=0; indexT<CurrSize->tsize; indexT++)
                    for(long indexZ=0; indexZ<newz; indexZ++)
                        for(long indexY=0; indexY<newy; indexY++)
                            for(long indexX=0; indexX<newx; indexX++)
                                bufferImages[current_buffer?0:1][indexX+indexY*newx+indexZ*newy*newx+indexT*newnumel]=bufferImages[current_buffer][indexX*2+indexY*2*CurrSize->xsize+indexZ*CurrSize->xsize*CurrSize->ysize+indexT*CurrSize->numel];




                delete [] bufferImages[current_buffer];
                bufferImages[current_buffer]= new SegPrecisionTYPE [newnumel*CurrSize->tsize];
                current_buffer=current_buffer?0:1;
                CurrSize->xsize=newx;
                CurrSize->ysize=newy;
                CurrSize->zsize=newz;
                CurrSize->numel=newnumel;
            }
            // *********************  Get max TP  *************************
            else if(strcmp(argv[i], "-tmax") == 0)
            {
                for(long i=0; i<CurrSize->numel; i++)
                {
                    float tmax=(float)-FLT_MAX;
                    for(long tp=0; tp<(long)CurrSize->tsize; tp++)
                    {
                        if(tmax<bufferImages[current_buffer][i+(long)(tp)*(long)CurrSize->numel])
                            tmax=bufferImages[current_buffer][i+(long)(tp)*(long)CurrSize->numel];
                    }
                    bufferImages[current_buffer?0:1][i]=tmax;
                }
                CurrSize->tsize=1;
                current_buffer=current_buffer?0:1;
            }
            // *********************  Get TP with maxval  *************************
            else if(strcmp(argv[i], "-tpmax") == 0)
            {
                for(long i=0; i<CurrSize->numel; i++)
                {
                    float tmax=(float)-FLT_MAX;
                    float tmaxindex=-1;
                    for(long tp=0; tp<(long)CurrSize->tsize; tp++)
                    {
                        if(bufferImages[current_buffer][i+(long)(tp)*(long)CurrSize->numel]>tmax)
                        {
                            tmax=bufferImages[current_buffer][i+(long)(tp)*(long)CurrSize->numel];
                            tmaxindex=(float)tp;
                        }
                    }
                    bufferImages[current_buffer?0:1][i]=(float)tmaxindex;
                }
                InputImage->cal_max=CurrSize->tsize;
                CurrSize->tsize=1;

                current_buffer=current_buffer?0:1;
            }
            // *********************  Get mean TP  *************************
            else if(strcmp(argv[i], "-tmean") == 0)
            {
                for(long i=0; i<CurrSize->numel; i++)
                {
                    float tmean=0;
                    for(long tp=0; tp<(long)CurrSize->tsize; tp++)
                    {
                        tmean+=bufferImages[current_buffer][i+(long)(tp)*CurrSize->numel];
                    }
                    bufferImages[current_buffer?0:1][i]=tmean/CurrSize->tsize;
                }
                CurrSize->tsize=1;
                current_buffer=current_buffer?0:1;
            }
            // *********************  Get min TP  *************************
            else if(strcmp(argv[i], "-tmin") == 0)
            {
                for(long i=0; i<CurrSize->numel; i++)
                {
                    float tmin=(float)FLT_MAX;
                    for(long tp=0; tp<CurrSize->tsize; tp++)
                    {
                        if(tmin>bufferImages[current_buffer][i+(int)(tp)*CurrSize->numel])
                            tmin=bufferImages[current_buffer][i+(int)(tp)*CurrSize->numel];
                    }
                    bufferImages[current_buffer?0:1][i]=tmin;
                }
                CurrSize->tsize=1;
                current_buffer=current_buffer?0:1;
            }
            // *********************  Reset SCL  *************************
            else if(strcmp(argv[i], "-scl") == 0)
            {
                InputImage->scl_inter=0;
                InputImage->scl_slope=1;

            }
            // *********************  Copy Header  *************************
            else if(strcmp(argv[i], "-hdr_copy") == 0)
            {
                string parser=argv[++i];

                nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
                if(NewImage->datatype!=DT_FLOAT32)
                {
                    seg_changeDatatype<SegPrecisionTYPE>(NewImage);
                }
                NewImage->nu=(NewImage->nu>1)?NewImage->nu:1;
                NewImage->nt=(NewImage->nt>1)?NewImage->nt:1;
                if(NewImage->nu<1)
                    NewImage->dim[5]=1;
                if(NewImage->nt<1)
                    NewImage->dim[4]=1;
                nifti_update_dims_from_array(NewImage);
                SegPrecisionTYPE * NewImagePtr = static_cast<SegPrecisionTYPE *>(NewImage->data);
                if(NewImage->nx==CurrSize->xsize&&NewImage->ny==CurrSize->ysize&&NewImage->nz==CurrSize->zsize&&NewImage->nt==CurrSize->tsize&&NewImage->nu==CurrSize->usize)
                {
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                        bufferImages[current_buffer?0:1][i]=NewImagePtr[i];
                    current_buffer=current_buffer?0:1;
                }
                else
                {
                    cout << "ERROR: Image "<< parser << " is the wrong size  -  original = ( "<<CurrSize->xsize<<","
                         <<CurrSize->ysize<<","<<CurrSize->zsize<<","<<CurrSize->tsize<<","<<CurrSize->usize<<" )  New image = ( "<<NewImage->nx<<","
                        <<NewImage->ny<<","<<NewImage->nz<<","<<NewImage->nt<<","<<NewImage->nu<<" )"<<endl;
                    exit(1);
                    i=argc;
                }
                nifti_image_free(NewImage);

            }
            // *********************  Copy Header  *************************
            else if(strcmp(argv[i], "-4to5") == 0)
            {

                int tempT=CurrSize->tsize;
                int tempU=CurrSize->usize;

                 InputImage->dim[4]=InputImage->nt=CurrSize->tsize=tempU;
                 InputImage->dim[5]=InputImage->nu=CurrSize->usize=tempT;


            }
            // *********************  Get LSSD  *************************
            else if(strcmp(argv[i], "-lssd") == 0)
            {
                string parser=argv[++i];
                nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
                NewImage->nu=(NewImage->nu>1)?NewImage->nu:1;
                NewImage->nt=(NewImage->nt>1)?NewImage->nt:1;
                if(NewImage->datatype!=DT_FLOAT32)
                {
                    seg_changeDatatype<SegPrecisionTYPE>(NewImage);
                }
                SegPrecisionTYPE * NewImagePtr = static_cast<SegPrecisionTYPE *>(NewImage->data);

                string parserstd=argv[++i];
                if(strtod(parserstd.c_str(),NULL)>0)
                {
                    if(NewImage->nt<2&&NewImage->nx==InputImage->nx&&NewImage->ny==InputImage->ny&&NewImage->nz==InputImage->nz)
                    {
                        float * NewImageMean=new float [NewImage->nx*NewImage->ny*NewImage->nz];
                        float * NewImageStd=new float [NewImage->nx*NewImage->ny*NewImage->nz];
                        float * InputImageMean=new float [InputImage->nx*InputImage->ny*InputImage->nz];
                        float * InputImageStd=new float [InputImage->nx*InputImage->ny*InputImage->nz];
                        float allmeanNew=0;
                        float allmeanInput=0;
                        float allstdNew=0;
                        float allstdInput=0;
                        for(long index=0; index<InputImage->nx*InputImage->ny*InputImage->nz; index++)
                        {
                            allmeanNew+=NewImagePtr[index];
                            allmeanInput+=bufferImages[current_buffer][index];
                            NewImageMean[index]=NewImagePtr[index];
                            NewImageStd[index]=NewImagePtr[index]*NewImagePtr[index];
                            InputImageMean[index]=bufferImages[current_buffer][index];
                            InputImageStd[index]=bufferImages[current_buffer][index]*bufferImages[current_buffer][index];
                        }
                        allmeanNew=allmeanNew/(InputImage->nx*InputImage->ny*InputImage->nz);
                        allmeanInput=allmeanInput/(InputImage->nx*InputImage->ny*InputImage->nz);

                        GaussianFilter4D_cArray(NewImageMean,strtod(parserstd.c_str(),NULL),CurrSize);
                        GaussianFilter4D_cArray(NewImageStd,strtod(parserstd.c_str(),NULL),CurrSize);
                        GaussianFilter4D_cArray(InputImageMean,strtod(parserstd.c_str(),NULL),CurrSize);
                        GaussianFilter4D_cArray(InputImageStd,strtod(parserstd.c_str(),NULL),CurrSize);
                        for(long index=0; index<InputImage->nx*InputImage->ny*InputImage->nz; index++)
                        {
                            allstdNew+=(NewImagePtr[index]-allmeanNew)*(NewImagePtr[index]-allmeanNew);
                            allstdInput+=(bufferImages[current_buffer][index]-allmeanInput)*(bufferImages[current_buffer][index]-allmeanInput);
                        }
                        allstdNew=allstdNew/(InputImage->nx*InputImage->ny*InputImage->nz);
                        allstdInput=allstdInput/(InputImage->nx*InputImage->ny*InputImage->nz);
                        for(long index=0; index<InputImage->nx*InputImage->ny*InputImage->nz; index++)
                        {
                            NewImageStd[index]=NewImageStd[index]-NewImageMean[index]*NewImageMean[index];
                            InputImageStd[index]=InputImageStd[index]-InputImageMean[index]*InputImageMean[index];
                            bufferImages[current_buffer?0:1][index]=(bufferImages[current_buffer][index]-InputImageMean[index])/(sqrt(InputImageStd[index]+0.01*allstdInput))-(NewImagePtr[index]-NewImageMean[index])/(sqrt(NewImageStd[index]+0.01*allstdNew));
                        }
                        GaussianFilter4D_cArray(bufferImages[current_buffer?0:1],strtod(parserstd.c_str(),NULL),CurrSize);
                        for(long index=0; index<InputImage->nx*InputImage->ny*InputImage->nz; index++)
                        {
                            bufferImages[current_buffer?0:1][index]=bufferImages[current_buffer?0:1][index]*bufferImages[current_buffer?0:1][index];
                        }

                        current_buffer=current_buffer?0:1;
                        delete [] NewImageMean;
                        delete [] NewImageStd;
                        delete [] InputImageMean;
                        delete [] InputImageStd;
                    }
                    else
                    {
                        cout << "ERROR: Image "<< parser << " is the wrong size  -  original = ( "<<CurrSize->xsize<<","
                             <<CurrSize->ysize<<","<<CurrSize->zsize<<","<<CurrSize->tsize<<","<<CurrSize->usize<<" )  New image = ( "<<NewImage->nx<<","
                            <<NewImage->ny<<","<<NewImage->nz<<","<<NewImage->nt<<","<<NewImage->nu<<" )"<<endl;
                        i=argc;
                    }
                }
            }
            // *********************  Get LNCC  *************************
            else if(strcmp(argv[i], "-lncc") == 0)
            {
                string parser=argv[++i];
                nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
                NewImage->nu=(NewImage->nu>1)?NewImage->nu:1;
                NewImage->nt=(NewImage->nt>1)?NewImage->nt:1;
                if(NewImage->datatype!=DT_FLOAT32)
                {
                    seg_changeDatatype<float>(NewImage);
                }
                SegPrecisionTYPE * NewImagePtr = static_cast<SegPrecisionTYPE *>(NewImage->data);

                string parserstd=argv[++i];
                if(strtod(parserstd.c_str(),NULL))
                {
                    if(NewImage->nt<2&&NewImage->nx==InputImage->nx&&NewImage->ny==InputImage->ny&&NewImage->nz==InputImage->nz)
                    {
                        float * NewImageMean=new float [NewImage->nx*NewImage->ny*NewImage->nz];
                        float * NewImageStd=new float [NewImage->nx*NewImage->ny*NewImage->nz];
                        float * InputImageMean=new float [InputImage->nx*InputImage->ny*InputImage->nz];
                        float * InputImageStd=new float [InputImage->nx*InputImage->ny*InputImage->nz];
                        float allmeanNew=0;
                        float allmeanInput=0;
                        float allstdNew=0;
                        float allstdInput=0;
                        for(long index=0; index<InputImage->nx*InputImage->ny*InputImage->nz; index++)
                        {
                            if(!isnan(NewImagePtr[index]) && !isnan(bufferImages[current_buffer][index])){
                            allmeanNew+=NewImagePtr[index];
                            NewImageMean[index]=NewImagePtr[index];
                            NewImageStd[index]=NewImagePtr[index]*NewImagePtr[index];
                            allmeanInput+=bufferImages[current_buffer][index];
                            InputImageMean[index]=bufferImages[current_buffer][index];
                            InputImageStd[index]=bufferImages[current_buffer][index]*bufferImages[current_buffer][index];
                            }
                        }
                        allmeanNew=allmeanNew/(InputImage->nx*InputImage->ny*InputImage->nz);
                        allmeanInput=allmeanInput/(InputImage->nx*InputImage->ny*InputImage->nz);
                        for(long index=0; index<InputImage->nx*InputImage->ny*InputImage->nz; index++)
                        {
                            if(!isnan(NewImagePtr[index]) && !isnan(bufferImages[current_buffer][index])){
                            allstdNew+=(NewImagePtr[index]-allmeanNew)*(NewImagePtr[index]-allmeanNew);
                            allstdInput+=(bufferImages[current_buffer][index]-allmeanInput)*(bufferImages[current_buffer][index]-allmeanInput);
                            bufferImages[current_buffer][index]=NewImagePtr[index]*bufferImages[current_buffer][index];
                            }
                        }
                        allstdNew=allstdNew/(InputImage->nx*InputImage->ny*InputImage->nz);
                        allstdInput=allstdInput/(InputImage->nx*InputImage->ny*InputImage->nz);
                        cout << allstdInput <<"  "<< allstdNew<<endl;
                        GaussianFilter4D_cArray(bufferImages[current_buffer],strtod(parserstd.c_str(),NULL),CurrSize);
                        GaussianFilter4D_cArray(NewImageMean,strtod(parserstd.c_str(),NULL),CurrSize);
                        GaussianFilter4D_cArray(NewImageStd,strtod(parserstd.c_str(),NULL),CurrSize);
                        GaussianFilter4D_cArray(InputImageMean,strtod(parserstd.c_str(),NULL),CurrSize);
                        GaussianFilter4D_cArray(InputImageStd,strtod(parserstd.c_str(),NULL),CurrSize);
                        for(long index=0; index<InputImage->nx*InputImage->ny*InputImage->nz; index++)
                        {
                            if(!isnan(NewImagePtr[index]) && !isnan(bufferImages[current_buffer][index])){
                            NewImageStd[index]=NewImageStd[index]-NewImageMean[index]*NewImageMean[index];
                            InputImageStd[index]=InputImageStd[index]-InputImageMean[index]*InputImageMean[index];
                            bufferImages[current_buffer?0:1][index]=(bufferImages[current_buffer][index]-InputImageMean[index]*NewImageMean[index])/(sqrt(NewImageStd[index]*InputImageStd[index])+sqrt(0.01*(allstdNew+allstdInput)));
                            }
                            else{
                                bufferImages[current_buffer?0:1][index]=std::numeric_limits<float>::quiet_NaN();
                            }
                        }
                        current_buffer=current_buffer?0:1;
                        delete [] NewImageMean;
                        delete [] NewImageStd;
                        delete [] InputImageMean;
                        delete [] InputImageStd;
                    }
                    else
                    {
                        cout << "ERROR: Image "<< parser << " is the wrong size  -  original = ( "<<CurrSize->xsize<<","
                             <<CurrSize->ysize<<","<<CurrSize->zsize<<","<<CurrSize->tsize<<","<<CurrSize->usize<<" )  New image = ( "<<NewImage->nx<<","
                            <<NewImage->ny<<","<<NewImage->nz<<","<<NewImage->nt<<","<<NewImage->nu<<" )"<<endl;
                        i=argc;
                    }
                }
                else
                {
                    cout << "ERROR: "<< string() << " is not a float"<<endl;
                    i=argc;
                }
                nifti_image_free(NewImage);
            }
            // ********************* z score ****************************

            else if(strcmp(argv[i], "-z") == 0)
            {
                for (long tup=0; tup<(CurrSize->tsize*CurrSize->usize); tup++)
                {
                    float mean=0;
                    int img3Dsize=(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize);
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize); i++)
                    {
                        mean+=bufferImages[current_buffer][i+img3Dsize*tup];
                    }
                    mean/=(float)(img3Dsize);
                    float std=0;
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize); i++)
                    {
                        std+=powf((bufferImages[current_buffer][i+img3Dsize*tup]-mean),2);
                    }
                    std/=(float)img3Dsize;
                    std=sqrt(std);
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize); i++)
                    {
                        bufferImages[current_buffer?0:1][i]=(bufferImages[current_buffer][i+img3Dsize*tup]-mean)/std;
                    }
                    current_buffer=current_buffer?0:1;
                }
            }

            // ********************* z score ****************************

            else if(strcmp(argv[i], "-zr") == 0)
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


                for (long tup=0; tup<(CurrSize->tsize*CurrSize->usize); tup++)
                {

                    long img3Dsize=(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize);

                    float * imgsort=new float [img3Dsize];
                    long curindex=0;
                    for(long index=0; index<img3Dsize; index++)
                    {
                        imgsort[curindex]=bufferImages[current_buffer][index+img3Dsize*tup];
                        curindex++;
                    }
                    HeapSort(imgsort,img3Dsize-1);
                    float lowThresh=imgsort[(long)(round(percentile*(img3Dsize-1)))];
                    float highThresh=imgsort[(long)(round((1-percentile)*(img3Dsize-1)))];
                    delete [] imgsort;


                    float mean=0;
                    long count=0;
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize); i++)
                    {
                        if(bufferImages[current_buffer][i+img3Dsize*tup]<highThresh && bufferImages[current_buffer][i+img3Dsize*tup]>lowThresh)
                        {
                            mean+=bufferImages[current_buffer][i+img3Dsize*tup];
                            count++;
                        }
                    }
                    mean/=(float)(count);
                    float std=0;
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize); i++)
                    {
                        if(bufferImages[current_buffer][i+img3Dsize*tup]<highThresh && bufferImages[current_buffer][i+img3Dsize*tup]>lowThresh)
                        {
                            std+=powf((bufferImages[current_buffer][i+img3Dsize*tup]-mean),2);
                        }
                    }
                    std/=(float)(count);
                    std=sqrt(std);
                    for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize); i++)
                    {
                        bufferImages[current_buffer?0:1][i]=(bufferImages[current_buffer][i+img3Dsize*tup]-mean)/std;
                    }
                    current_buffer=current_buffer?0:1;
                }
            }

            else if(strcmp(argv[i], "-flipNM") == 0) // Neuromorphometric Lab Flip for version 3
            {
                for(long i=0; i<(long)(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize); i++)
                {
                    switch((int)floor(bufferImages[current_buffer][i])){
                    case 0: bufferImages[current_buffer?0:1][i]=0; break; break; // Background and skull
                    case 1: bufferImages[current_buffer?0:1][i]=1; break; // NonBrain Low
                    case 2: bufferImages[current_buffer?0:1][i]=2; break; // NonBrain Mid
                    case 3: bufferImages[current_buffer?0:1][i]=3; break; // NonBrain High
                    case 4: bufferImages[current_buffer?0:1][i]=4; break; // Non-ventricular CSF
                    case 5: bufferImages[current_buffer?0:1][i]=5; break; // 3rd Ventricle
                    case 12: bufferImages[current_buffer?0:1][i]=12; break; // 4th Ventricle
                    case 16: bufferImages[current_buffer?0:1][i]=16; break; // 5th Ventricle
                    case 24: bufferImages[current_buffer?0:1][i]=31; break; // Right to Left Accumbens Area
                    case 31: bufferImages[current_buffer?0:1][i]=24; break; // Left to Right Accumbens Area
                    case 32: bufferImages[current_buffer?0:1][i]=33; break; // Right to Left Amygdala
                    case 33: bufferImages[current_buffer?0:1][i]=32; break; // Left to Right Amygdala
                    case 35: bufferImages[current_buffer?0:1][i]=35; break; // Pons
                    case 36: bufferImages[current_buffer?0:1][i]=36; break; // Brain Stem
                    case 37: bufferImages[current_buffer?0:1][i]=38; break; // Right to Left Caudate
                    case 38: bufferImages[current_buffer?0:1][i]=37; break; // Left to Right Caudate
                    case 39: bufferImages[current_buffer?0:1][i]=40; break; // Right to Left Cerebellum Exterior
                    case 40: bufferImages[current_buffer?0:1][i]=39; break; // Left to Right Cerebellum Exterior
                    case 41: bufferImages[current_buffer?0:1][i]=42; break; // Right to Left Cerebellum White Matter
                    case 42: bufferImages[current_buffer?0:1][i]=41; break; // Left to Right Cerebellum White Matter
                    case 43: bufferImages[current_buffer?0:1][i]=44; break; // Right to Left Cerebral Exterior
                    case 44: bufferImages[current_buffer?0:1][i]=43; break; // Left to Right Cerebral Exterior

                    case 47: bufferImages[current_buffer?0:1][i]=47; break; // 3rd Ventricle (Posterior part)
                    case 48: bufferImages[current_buffer?0:1][i]=49; break; // Right to Left Hippocampus
                    case 49: bufferImages[current_buffer?0:1][i]=48; break; // Left to Right Hippocampus
                    case 50: bufferImages[current_buffer?0:1][i]=51; break; // Right to Left Inf Lat Vent
                    case 51: bufferImages[current_buffer?0:1][i]=50; break; // Left to Right Inf Lat Vent
                    case 52: bufferImages[current_buffer?0:1][i]=53; break; // Right to Left Lateral Ventricle
                    case 53: bufferImages[current_buffer?0:1][i]=52; break; // Left to Right Lateral Ventricle

                    case 54: bufferImages[current_buffer?0:1][i]=55; break; // Right to Left Lesion
                    case 55: bufferImages[current_buffer?0:1][i]=54; break; // Left to Right Lesion

                    case 56: bufferImages[current_buffer?0:1][i]=57; break; // Right to Left Pallidum
                    case 57: bufferImages[current_buffer?0:1][i]=56; break; // Left to Right Pallidum
                    case 58: bufferImages[current_buffer?0:1][i]=59; break; // Right to Left Putamen
                    case 59: bufferImages[current_buffer?0:1][i]=58; break; // Left to Right Putamen
                    case 60: bufferImages[current_buffer?0:1][i]=61; break; // Right to Left Thalamus Proper
                    case 61: bufferImages[current_buffer?0:1][i]=60; break; // Left to Right Thalamus Proper
                    case 62: bufferImages[current_buffer?0:1][i]=63; break; // Right to Left Ventral DC
                    case 63: bufferImages[current_buffer?0:1][i]=62; break; // Left to Right Ventral DC
                    case 64: bufferImages[current_buffer?0:1][i]=65; break; // Right to Left vessel
                    case 65: bufferImages[current_buffer?0:1][i]=64; break; // Left to Right vessel

                    case 66: bufferImages[current_buffer?0:1][i]=67; break; // Right to Left Ventricular Lining
                    case 67: bufferImages[current_buffer?0:1][i]=66; break; // Left to Right Ventricular Lining

                    case 70: bufferImages[current_buffer?0:1][i]=70; break; // Optic Chiasm
                    case 72: bufferImages[current_buffer?0:1][i]=72; break; // Cerebellar Vermal Lobules I-V
                    case 73: bufferImages[current_buffer?0:1][i]=73; break; // Cerebellar Vermal Lobules VI-VII
                    case 74: bufferImages[current_buffer?0:1][i]=74; break; // Cerebellar Vermal Lobules VIII-X
                    case 77: bufferImages[current_buffer?0:1][i]=76; break; // Right to Left Basal Forebrain
                    case 76: bufferImages[current_buffer?0:1][i]=77; break; // Left to Right Basal Forebrain

                    case 81: bufferImages[current_buffer?0:1][i]=89; break; // WM region flips
                    case 82: bufferImages[current_buffer?0:1][i]=90; break;
                    case 83: bufferImages[current_buffer?0:1][i]=91; break;
                    case 84: bufferImages[current_buffer?0:1][i]=92; break;
                    case 85: bufferImages[current_buffer?0:1][i]=93; break;
                    case 86: bufferImages[current_buffer?0:1][i]=94; break;

                    case 87: bufferImages[current_buffer?0:1][i]=87; break;

                    case 89: bufferImages[current_buffer?0:1][i]=81; break;
                    case 90: bufferImages[current_buffer?0:1][i]=82; break;
                    case 91: bufferImages[current_buffer?0:1][i]=83; break;
                    case 92: bufferImages[current_buffer?0:1][i]=84; break;
                    case 93: bufferImages[current_buffer?0:1][i]=85; break;
                    case 94: bufferImages[current_buffer?0:1][i]=86; break;

                    case 96: bufferImages[current_buffer?0:1][i]=97; break; // Right to Left claustrum
                    case 97: bufferImages[current_buffer?0:1][i]=96; break; // Left to Right claustrum

                    case 101: bufferImages[current_buffer?0:1][i]=102; break; // Right to Left ACgG anterior cingulate gyrus
                    case 102: bufferImages[current_buffer?0:1][i]=101; break; // Left to Right ACgG anterior cingulate gyrus
                    case 103: bufferImages[current_buffer?0:1][i]=104; break; // Right to Left AIns anterior insula
                    case 104: bufferImages[current_buffer?0:1][i]=103; break; // Left to Right AIns anterior insula
                    case 105: bufferImages[current_buffer?0:1][i]=106; break; // Right to Left AOrG anterior orbital gyrus
                    case 106: bufferImages[current_buffer?0:1][i]=105; break; // Left to Right AOrG anterior orbital gyrus
                    case 107: bufferImages[current_buffer?0:1][i]=108; break; // Right to Left AnG angular gyrus
                    case 108: bufferImages[current_buffer?0:1][i]=107; break; // Left to Right AnG angular gyrus
                    case 109: bufferImages[current_buffer?0:1][i]=110; break; // Right to Left Calc calcarine cortex
                    case 110: bufferImages[current_buffer?0:1][i]=109; break; // Left to Right Calc calcarine cortex
                    case 113: bufferImages[current_buffer?0:1][i]=114; break; // Right to Left CO central operculum
                    case 114: bufferImages[current_buffer?0:1][i]=113; break; // Left to Right CO central operculum
                    case 115: bufferImages[current_buffer?0:1][i]=116; break; // Right to Left Cun cuneus
                    case 116: bufferImages[current_buffer?0:1][i]=115; break; // Left to Right Cun cuneus
                    case 117: bufferImages[current_buffer?0:1][i]=118; break; // Right to Left Ent entorhinal area
                    case 118: bufferImages[current_buffer?0:1][i]=117; break; // Left to Right Ent entorhinal area
                    case 119: bufferImages[current_buffer?0:1][i]=120; break; // Right to Left FO frontal operculum
                    case 120: bufferImages[current_buffer?0:1][i]=119; break; // Left to Right FO frontal operculum
                    case 121: bufferImages[current_buffer?0:1][i]=122; break; // Right to Left FRP frontal pole
                    case 122: bufferImages[current_buffer?0:1][i]=121; break; // Left to Right FRP frontal pole
                    case 123: bufferImages[current_buffer?0:1][i]=124; break; // Right to Left FuG fusiform gyrus
                    case 124: bufferImages[current_buffer?0:1][i]=123; break; // Left to Right FuG fusiform gyrus
                    case 125: bufferImages[current_buffer?0:1][i]=126; break; // Right to Left GRe gyrus rectus
                    case 126: bufferImages[current_buffer?0:1][i]=125; break; // Left to Right GRe gyrus rectus
                    case 129: bufferImages[current_buffer?0:1][i]=130; break; // Right to Left IOG inferior occipital gyrus
                    case 130: bufferImages[current_buffer?0:1][i]=129; break; // Left to Right IOG inferior occipital gyrus
                    case 133: bufferImages[current_buffer?0:1][i]=134; break; // Right to Left ITG inferior temporal gyrus
                    case 134: bufferImages[current_buffer?0:1][i]=133; break; // Left to Right ITG inferior temporal gyrus
                    case 135: bufferImages[current_buffer?0:1][i]=136; break; // Right to Left LiG lingual gyrus
                    case 136: bufferImages[current_buffer?0:1][i]=135; break; // Left to Right LiG lingual gyrus
                    case 137: bufferImages[current_buffer?0:1][i]=138; break; // Right to Left LOrG lateral orbital gyrus
                    case 138: bufferImages[current_buffer?0:1][i]=137; break; // Left to Right LOrG lateral orbital gyrus
                    case 139: bufferImages[current_buffer?0:1][i]=140; break; // Right to Left MCgG middle cingulate gyrus
                    case 140: bufferImages[current_buffer?0:1][i]=139; break; // Left to Right MCgG middle cingulate gyrus
                    case 141: bufferImages[current_buffer?0:1][i]=142; break; // Right to Left MFC medial frontal cortex
                    case 142: bufferImages[current_buffer?0:1][i]=141; break; // Left to Right MFC medial frontal cortex
                    case 143: bufferImages[current_buffer?0:1][i]=144; break; // Right to Left MFG middle frontal gyrus
                    case 144: bufferImages[current_buffer?0:1][i]=143; break; // Left to Right MFG middle frontal gyrus
                    case 145: bufferImages[current_buffer?0:1][i]=146; break; // Right to Left MOG middle occipital gyrus
                    case 146: bufferImages[current_buffer?0:1][i]=145; break; // Left to Right MOG middle occipital gyrus
                    case 147: bufferImages[current_buffer?0:1][i]=148; break; // Right to Left MOrG medial orbital gyrus
                    case 148: bufferImages[current_buffer?0:1][i]=147; break; // Left to Right MOrG medial orbital gyrus
                    case 149: bufferImages[current_buffer?0:1][i]=150; break; // Right to Left MPoG postcentral gyrus medial segment
                    case 150: bufferImages[current_buffer?0:1][i]=149; break; // Left to Right MPoG postcentral gyrus medial segment
                    case 151: bufferImages[current_buffer?0:1][i]=152; break; // Right to Left MPrG precentral gyrus medial segment
                    case 152: bufferImages[current_buffer?0:1][i]=151; break; // Left to Right MPrG precentral gyrus medial segment
                    case 153: bufferImages[current_buffer?0:1][i]=154; break; // Right to Left MSFG superior frontal gyrus medial segment
                    case 154: bufferImages[current_buffer?0:1][i]=153; break; // Left to Right MSFG superior frontal gyrus medial segment
                    case 155: bufferImages[current_buffer?0:1][i]=156; break; // Right to Left MTG middle temporal gyrus
                    case 156: bufferImages[current_buffer?0:1][i]=155; break; // Left to Right MTG middle temporal gyrus
                    case 157: bufferImages[current_buffer?0:1][i]=158; break; // Right to Left OCP occipital pole
                    case 158: bufferImages[current_buffer?0:1][i]=157; break; // Left to Right OCP occipital pole
                    case 161: bufferImages[current_buffer?0:1][i]=162; break; // Right to Left OFuG occipital fusiform gyrus
                    case 162: bufferImages[current_buffer?0:1][i]=161; break; // Left to Right OFuG occipital fusiform gyrus
                    case 163: bufferImages[current_buffer?0:1][i]=164; break; // Right to Left OpIFG opercular part of the inferior frontal gyrus
                    case 164: bufferImages[current_buffer?0:1][i]=163; break; // Left to Right OpIFG opercular part of the inferior frontal gyrus
                    case 165: bufferImages[current_buffer?0:1][i]=166; break; // Right to Left OrIFG orbital part of the inferior frontal gyrus
                    case 166: bufferImages[current_buffer?0:1][i]=165; break; // Left to Right OrIFG orbital part of the inferior frontal gyrus
                    case 167: bufferImages[current_buffer?0:1][i]=168; break; // Right to Left PCgG posterior cingulate gyrus
                    case 168: bufferImages[current_buffer?0:1][i]=167; break; // Left to Right PCgG posterior cingulate gyrus
                    case 169: bufferImages[current_buffer?0:1][i]=170; break; // Right to Left PCu precuneus
                    case 170: bufferImages[current_buffer?0:1][i]=169; break; // Left to Right PCu precuneus
                    case 171: bufferImages[current_buffer?0:1][i]=172; break; // Right to Left PHG parahippocampal gyrus
                    case 172: bufferImages[current_buffer?0:1][i]=171; break; // Left to Right PHG parahippocampal gyrus
                    case 173: bufferImages[current_buffer?0:1][i]=174; break; // Right to Left PIns posterior insula
                    case 174: bufferImages[current_buffer?0:1][i]=173; break; // Left to Right PIns posterior insula
                    case 175: bufferImages[current_buffer?0:1][i]=176; break; // Right to Left PO parietal operculum
                    case 176: bufferImages[current_buffer?0:1][i]=175; break; // Left to Right PO parietal operculum
                    case 177: bufferImages[current_buffer?0:1][i]=178; break; // Right to Left PoG postcentral gyrus
                    case 178: bufferImages[current_buffer?0:1][i]=177; break; // Left to Right PoG postcentral gyrus
                    case 179: bufferImages[current_buffer?0:1][i]=180; break; // Right to Left POrG posterior orbital gyrus
                    case 180: bufferImages[current_buffer?0:1][i]=179; break; // Left to Right POrG posterior orbital gyrus
                    case 181: bufferImages[current_buffer?0:1][i]=182; break; // Right to Left PP planum polare
                    case 182: bufferImages[current_buffer?0:1][i]=181; break; // Left to Right PP planum polare
                    case 183: bufferImages[current_buffer?0:1][i]=184; break; // Right to Left PrG precentral gyrus
                    case 184: bufferImages[current_buffer?0:1][i]=183; break; // Left to Right PrG precentral gyrus
                    case 185: bufferImages[current_buffer?0:1][i]=186; break; // Right to Left PT planum temporale
                    case 186: bufferImages[current_buffer?0:1][i]=185; break; // Left to Right PT planum temporale
                    case 187: bufferImages[current_buffer?0:1][i]=188; break; // Right to Left SCA subcallosal area
                    case 188: bufferImages[current_buffer?0:1][i]=187; break; // Left to Right SCA subcallosal area
                    case 191: bufferImages[current_buffer?0:1][i]=192; break; // Right to Left SFG superior frontal gyrus
                    case 192: bufferImages[current_buffer?0:1][i]=191; break; // Left to Right SFG superior frontal gyrus
                    case 193: bufferImages[current_buffer?0:1][i]=194; break; // Right to Left SMC supplementary motor cortex
                    case 194: bufferImages[current_buffer?0:1][i]=193; break; // Left to Right SMC supplementary motor cortex
                    case 195: bufferImages[current_buffer?0:1][i]=196; break; // Right to Left SMG supramarginal gyrus
                    case 196: bufferImages[current_buffer?0:1][i]=195; break; // Left to Right SMG supramarginal gyrus
                    case 197: bufferImages[current_buffer?0:1][i]=198; break; // Right to Left SOG superior occipital gyrus
                    case 198: bufferImages[current_buffer?0:1][i]=197; break; // Left to Right SOG superior occipital gyrus
                    case 199: bufferImages[current_buffer?0:1][i]=200; break; // Right to Left SPL superior parietal lobule
                    case 200: bufferImages[current_buffer?0:1][i]=199; break; // Left to Right SPL superior parietal lobule
                    case 201: bufferImages[current_buffer?0:1][i]=202; break; // Right to Left STG superior temporal gyrus
                    case 202: bufferImages[current_buffer?0:1][i]=201; break; // Left to Right STG superior temporal gyrus
                    case 203: bufferImages[current_buffer?0:1][i]=204; break; // Right to Left TMP temporal pole
                    case 204: bufferImages[current_buffer?0:1][i]=203; break; // Left to Right TMP temporal pole
                    case 205: bufferImages[current_buffer?0:1][i]=206; break; // Right to Left TrIFG triangular part of the inferior frontal gyrus
                    case 206: bufferImages[current_buffer?0:1][i]=205; break; // Left to Right TrIFG triangular part of the inferior frontal gyrus
                    case 207: bufferImages[current_buffer?0:1][i]=208; break; // Right to Left TTG transverse temporal gyrus
                    case 208: bufferImages[current_buffer?0:1][i]=207; break; // Left to Right TTG transverse temporal gyrus

                    }
                }

                current_buffer=current_buffer?0:1;
                for(long indexZ=0; indexZ<CurrSize->zsize; indexZ++)
                    for(long indexY=0; indexY<CurrSize->ysize; indexY++)
                        for(long indexX=0; indexX<CurrSize->xsize; indexX++)
                            bufferImages[current_buffer?0:1][((CurrSize->xsize-1-indexX)+indexY*CurrSize->xsize+indexZ*CurrSize->ysize*CurrSize->xsize)]=bufferImages[current_buffer][indexX+indexY*CurrSize->xsize+indexZ*CurrSize->ysize*CurrSize->xsize];

                current_buffer=current_buffer?0:1;


            }

            else if(strcmp(argv[i], "-fliplab") == 0)
            {
                for(long indexZ=1; indexZ<(CurrSize->zsize-1); indexZ++){
                    for(long indexY=1; indexY<(CurrSize->ysize-1); indexY++){
                        for(long indexX=1; indexX<(CurrSize->xsize-1); indexX++){
                            int indexcur=indexX+indexY*CurrSize->xsize+indexZ*CurrSize->ysize*CurrSize->xsize;
                            float curval=bufferImages[current_buffer][indexcur];
                            if(     curval!= 52 &&
                                    curval!= 53 &&
                                    curval!= 47 &&
                                    curval!= 50 &&
                                    curval!= 51 &&
                                    curval!= 46 &&
                                    curval!= 45){
                                int shiftrealsize=1;
                                int shiftspacing=1;

                                int stop=0;
                                for(int shiftz=-shiftrealsize; shiftz<=shiftrealsize; shiftz+=shiftspacing){
                                    for(int shifty=-shiftrealsize; shifty<=shiftrealsize; shifty+=shiftspacing){
                                        for(int shiftx=-shiftrealsize; shiftx<=shiftrealsize; shiftx+=shiftspacing){
                                            int index1=(indexX+shiftx)+CurrSize->xsize*(indexY+shifty)+CurrSize->xsize*CurrSize->ysize*(indexZ+shiftz);
                                            int index2=(indexX-shiftx)+CurrSize->xsize*(indexY-shifty)+CurrSize->xsize*CurrSize->ysize*(indexZ-shiftz);
                                            float curval1=bufferImages[current_buffer][index1];
                                            float curval2=bufferImages[current_buffer][index2];
                                            if(stop==0 && (fabs(shiftx)+fabs(shifty)+fabs(shiftz))<2 ){
                                                if(curval1==46){
                                                    if(curval2==47|| curval2==51|| curval2==53){
                                                        bufferImages[current_buffer?0:1][indexcur]=46;
                                                        stop=1;
                                                    }
                                                    else{
                                                        bufferImages[current_buffer?0:1][indexcur]=bufferImages[current_buffer][indexcur];

                                                    }

                                                }
                                                else if(curval1==45 ){
                                                    if(curval2==52 || curval2==50 || curval2==47 ){
                                                        bufferImages[current_buffer?0:1][indexcur]=45;
                                                        stop=1;
                                                    }
                                                    else{
                                                        bufferImages[current_buffer?0:1][indexcur]=bufferImages[current_buffer][indexcur];

                                                    }

                                                }
                                                else  if(curval1==47|| curval1==51|| curval1==53){
                                                    if(curval2==46 ){
                                                        bufferImages[current_buffer?0:1][indexcur]=46;
                                                        stop=1;
                                                    }
                                                    else{
                                                        bufferImages[current_buffer?0:1][indexcur]=bufferImages[current_buffer][indexcur];

                                                    }

                                                }
                                                else if(curval1==52 || curval1==50 || curval1==47 ){
                                                    if(curval2==45 ){
                                                        bufferImages[current_buffer?0:1][indexcur]=45;
                                                        stop=1;
                                                    }
                                                    else{
                                                        bufferImages[current_buffer?0:1][indexcur]=bufferImages[current_buffer][indexcur];
                                                    }

                                                }
                                                else{
                                                    bufferImages[current_buffer?0:1][indexcur]=bufferImages[current_buffer][indexcur];

                                                }
                                            }

                                        }
                                    }
                                }

                            }
                            else{
                                bufferImages[current_buffer?0:1][indexX+indexY*CurrSize->xsize+indexZ*CurrSize->ysize*CurrSize->xsize]=bufferImages[current_buffer][indexX+indexY*CurrSize->xsize+indexZ*CurrSize->ysize*CurrSize->xsize];

                            }
                        }
                    }
                }

                current_buffer=current_buffer?0:1;

            }
            else if(strcmp(argv[i], "-fliplab2") == 0)
              {
                for(long indexZ=1; indexZ<(CurrSize->zsize-1); indexZ++){
                    for(long indexY=1; indexY<(CurrSize->ysize-1); indexY++){
                        for(long indexX=1; indexX<(CurrSize->xsize-1); indexX++){
                            int indexcur=indexX+indexY*CurrSize->xsize+indexZ*CurrSize->ysize*CurrSize->xsize;
                            float curval=bufferImages[current_buffer][indexcur];
                            if( curval==46 || curval== 45){
                                int shiftrealsize=1;
                                int shiftspacing=1;

                                int stop=0;
                                for(int shiftz=-shiftrealsize; shiftz<=shiftrealsize; shiftz+=shiftspacing){
                                    for(int shifty=-shiftrealsize; shifty<=shiftrealsize; shifty+=shiftspacing){
                                        for(int shiftx=-shiftrealsize; shiftx<=shiftrealsize; shiftx+=shiftspacing){
                                            if( stop==0 && (fabs(shiftz)+fabs(shifty)+fabs(shiftx))<2 ){

                                                int index2=(indexX-shiftx)+CurrSize->xsize*(indexY-shifty)+CurrSize->xsize*CurrSize->ysize*(indexZ-shiftz);
                                                float curval2=bufferImages[current_buffer][index2];

                                                    if(curval2==47|| curval2==51|| curval2==53){
                                                        bufferImages[current_buffer?0:1][indexcur]=67;
                                                        stop=1;
                                                        //cout<<"hit"<<endl;
                                                    }
                                                    else if(curval2==52 || curval2==50 || curval2==47 ){
                                                        bufferImages[current_buffer?0:1][indexcur]=66;
                                                        stop=1;
                                                        //cout<<"hat"<<endl;
                                                    }
                                                    else{
                                                        bufferImages[current_buffer?0:1][indexcur]=bufferImages[current_buffer][indexcur];
                                                    }
                                            }
                                        }
                                    }
                                }
                            }

                            else{
                                bufferImages[current_buffer?0:1][indexcur]=bufferImages[current_buffer][indexcur];

                            }
                        }
                    }
                }

                  current_buffer=current_buffer?0:1;

              }
            else if(strcmp(argv[i], "-flipimgx") == 0) // X flip image
            {

                for(long indexZ=0; indexZ<CurrSize->zsize; indexZ++)
                    for(long indexY=0; indexY<CurrSize->ysize; indexY++)
                        for(long indexX=1; indexX<(CurrSize->xsize-1); indexX++)
                            bufferImages[current_buffer?0:1][((CurrSize->xsize-1-indexX)+indexY*CurrSize->xsize+indexZ*CurrSize->ysize*CurrSize->xsize)]=bufferImages[current_buffer][indexX+indexY*CurrSize->xsize+indexZ*CurrSize->ysize*CurrSize->xsize];

                current_buffer=current_buffer?0:1;

            }
            else if(strcmp(argv[i], "-flipimgy") == 0) // Y flip image
            {

                for(long indexZ=0; indexZ<CurrSize->zsize; indexZ++)
                    for(long indexY=0; indexY<CurrSize->ysize; indexY++)
                        for(long indexX=1; indexX<(CurrSize->xsize-1); indexX++)
                            bufferImages[current_buffer?0:1][((indexX)+(CurrSize->ysize-1-indexY)*CurrSize->xsize+indexZ*CurrSize->ysize*CurrSize->xsize)]=bufferImages[current_buffer][indexX+indexY*CurrSize->xsize+indexZ*CurrSize->ysize*CurrSize->xsize];

                current_buffer=current_buffer?0:1;

            }
            else if(strcmp(argv[i], "-flipimgz") == 0) // Z flip image
            {

                for(long indexZ=0; indexZ<CurrSize->zsize; indexZ++)
                    for(long indexY=0; indexY<CurrSize->ysize; indexY++)
                        for(long indexX=1; indexX<(CurrSize->xsize-1); indexX++)
                            bufferImages[current_buffer?0:1][((indexX)+indexY*CurrSize->xsize+(CurrSize->zsize-1-indexZ)*CurrSize->ysize*CurrSize->xsize)]=bufferImages[current_buffer][indexX+indexY*CurrSize->xsize+indexZ*CurrSize->ysize*CurrSize->xsize];

                current_buffer=current_buffer?0:1;

            }
            // *********************  output data type  *************************
            else if(strcmp(argv[i], "-v") == 0)
            {
                verbose=1;
            }
            else if(strcmp(argv[i], "-odt") == 0)
            {
                string parser=argv[++i];
                if(parser.find("uchar")!=string::npos)
                {
                    datatypeoutput=NIFTI_TYPE_UINT8;
                }
                else if(parser.find("ushort")!=string::npos)
                {
                    datatypeoutput=NIFTI_TYPE_UINT16;
                }
                else if(parser.find("uint")!=string::npos)
                {
                    datatypeoutput=NIFTI_TYPE_UINT32;
                }
                else if(parser.find("char")!=string::npos)
                {
                    datatypeoutput=NIFTI_TYPE_INT8;
                }
                else if(parser.find("short")!=string::npos)
                {
                    datatypeoutput=NIFTI_TYPE_INT16;
                }
                else if(parser.find("int")!=string::npos)
                {
                    datatypeoutput=NIFTI_TYPE_INT32;
                }
                else if(parser.find("float")!=string::npos)
                {
                    datatypeoutput=NIFTI_TYPE_FLOAT32;
                }
                else if(parser.find("double")!=string::npos)
                {
                    datatypeoutput=NIFTI_TYPE_FLOAT64;
                }
                else
                {
                    cout << "ERROR: Datatype "<< parser << " is unknown"<<endl;
                    i=argc;
                }
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
                cout << "Option "<< string(argv[i]) << " unkown"<<endl;
                i=argc;
                return 0;
            }

        }
        string parser=argv[argc-1];
        if(parser.find(string(".nii"))>0 || parser.find(string(".img")) ||parser.find(string(".hdr"))>0)
        {
            // saving output
            char * filename_out=argv[argc-1];
            nifti_image * OutputImage = nifti_copy_nim_info(InputImage);
            OutputImage->datatype=datatypeoutput;
            nifti_set_filenames(OutputImage,filename_out,0,0);
            OutputImage->dim[1]=OutputImage->nx=CurrSize->xsize;
            OutputImage->dim[2]=OutputImage->ny=CurrSize->ysize;
            OutputImage->dim[3]=OutputImage->nz=CurrSize->zsize;
            OutputImage->dim[4]=OutputImage->nt=CurrSize->tsize;
            OutputImage->dim[5]=OutputImage->nu=CurrSize->usize;
            OutputImage->dim[6]=OutputImage->nv=1;
            OutputImage->dim[7]=OutputImage->nw=1;
            OutputImage->dim[0]=2;
            OutputImage->dim[0]=(OutputImage->dim[3]>1?3:OutputImage->dim[0]);
            OutputImage->dim[0]=(OutputImage->dim[4]>1?4:OutputImage->dim[0]);
            OutputImage->dim[0]=(OutputImage->dim[5]>1?5:OutputImage->dim[0]);
            OutputImage->dim[0]=(OutputImage->dim[6]>1?6:OutputImage->dim[0]);
            OutputImage->dim[0]=(OutputImage->dim[7]>1?7:OutputImage->dim[0]);

            //mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
            bool scalingdiff=false;
            for(long i=0; i<4; i++)
            {
                OutputImage->sto_xyz.m[i][i]/=Scalling[i];
                OutputImage->pixdim[i+1]/=Scalling[i];
                if(Scalling[i]!=1)
                {
                    scalingdiff=true;
                }
            }
            if(scalingdiff)
            {

                cout << "A scaling factor is present. Removing Sform"<<endl;
                OutputImage->sform_code=0;
            }
            //        OutputImage->qoffset_x=translation[0];
            //        OutputImage->qoffset_y=translation[1];
            //        OutputImage->qoffset_z=translation[2];


            if(verbose)
            {
                cout << "Output Dim = [ ";
                for(long i=0; i<8; i++)
                {
                    cout<<(float)OutputImage->dim[i];
                    if(i<7)
                    {
                        cout<<" , ";
                    }
                }
                cout<<" ] "<<endl;
                flush(cout);
            }
            nifti_update_dims_from_array(OutputImage);
            nifti_datatype_sizes(OutputImage->datatype,&OutputImage->nbyper,&OutputImage->swapsize);
            if(datatypeoutput==NIFTI_TYPE_UINT8)
            {
                OutputImage->data = (void *) calloc(CurrSize->numel*CurrSize->tsize*CurrSize->usize, sizeof(unsigned char));
                unsigned char * OutputImagePtr = static_cast<unsigned char *>(OutputImage->data);
                for(long i=0; i<(long)(CurrSize->numel*CurrSize->tsize*CurrSize->usize); i++)
                {
                    OutputImagePtr[i]=(unsigned char)round(bufferImages[current_buffer][i]);
                }
            }
            else if(datatypeoutput==NIFTI_TYPE_UINT16)
            {
                OutputImage->data = (void *) calloc(OutputImage->nvox, sizeof(unsigned short));
                unsigned short * OutputImagePtr = static_cast<unsigned short *>(OutputImage->data);
                for(long i=0; i<(long)(CurrSize->numel*CurrSize->tsize*CurrSize->usize); i++)
                {
                    OutputImagePtr[i]=(unsigned short)round(bufferImages[current_buffer][i]);
                }
            }
            else if(datatypeoutput==NIFTI_TYPE_UINT32)
            {
                OutputImage->data = (void *) calloc(CurrSize->numel*CurrSize->tsize*CurrSize->usize, sizeof(unsigned int));
                unsigned int * OutputImagePtr = static_cast<unsigned int *>(OutputImage->data);
                for(long i=0; i<(long)(CurrSize->numel*CurrSize->tsize*CurrSize->usize); i++)
                {
                    OutputImagePtr[i]=(unsigned int)round(bufferImages[current_buffer][i]);
                }
            }
            else if(datatypeoutput==NIFTI_TYPE_INT8)
            {
                OutputImage->data = (void *) calloc(CurrSize->numel*CurrSize->tsize*CurrSize->usize, sizeof(char));
                char * OutputImagePtr = static_cast<char *>(OutputImage->data);
                for(long i=0; i<(long)(CurrSize->numel*CurrSize->tsize*CurrSize->usize); i++)
                {
                    OutputImagePtr[i]=(char)round(bufferImages[current_buffer][i]);
                }
            }
            else if(datatypeoutput==NIFTI_TYPE_INT16)
            {
                OutputImage->data = (void *) calloc(CurrSize->numel*CurrSize->tsize*CurrSize->usize, sizeof(short));
                short * OutputImagePtr = static_cast<short *>(OutputImage->data);
                for(long i=0; i<(long)(CurrSize->numel*CurrSize->tsize*CurrSize->usize); i++)
                {
                    OutputImagePtr[i]=(short)round(bufferImages[current_buffer][i]);
                }
            }
            else if(datatypeoutput==NIFTI_TYPE_INT32)
            {
                OutputImage->data = (void *) calloc(CurrSize->numel*CurrSize->tsize*CurrSize->usize, sizeof(int));
                int * OutputImagePtr = static_cast<int *>(OutputImage->data);
                for(long i=0; i<(long)(CurrSize->numel*CurrSize->tsize*CurrSize->usize); i++)
                {
                    OutputImagePtr[i]=(int)round(bufferImages[current_buffer][i]);
                }
            }
            else if(datatypeoutput==NIFTI_TYPE_FLOAT32)
            {
                OutputImage->data = (void *) calloc(CurrSize->numel*CurrSize->tsize*CurrSize->usize, sizeof(float));
                float * OutputImagePtr = static_cast<float *>(OutputImage->data);
                for(long i=0; i<(long)(CurrSize->numel*CurrSize->tsize*CurrSize->usize); i++)
                {
                    OutputImagePtr[i]=(float)bufferImages[current_buffer][i];
                }
            }
            else if(datatypeoutput==NIFTI_TYPE_FLOAT64)
            {
                OutputImage->data = (void *) calloc(CurrSize->numel*CurrSize->tsize*CurrSize->usize, sizeof(double));
                double * OutputImagePtr = static_cast<double *>(OutputImage->data);
                for(long i=0; i<(long)(CurrSize->numel*CurrSize->tsize*CurrSize->usize); i++)
                {
                    OutputImagePtr[i]=(double)round(bufferImages[current_buffer][i]);
                }
            }
            nifti_image_write(OutputImage);
            nifti_image_free(OutputImage);
        }

        delete [] bufferImages[0];
        delete [] bufferImages[1];
        delete [] bufferImages;
        delete [] CurrSize;

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

