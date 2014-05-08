#include "_seg_tools.h"
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <fstream>
#include <new>
#include <exception>
using namespace std;
#define SegPrecisionTYPE float


void Usage(char *exec)
{
    printf("\nUsage:\t%s <in> [constrains] [statistics]\n\n",exec);

    printf("\t* * Constrains (optional) * *\n");
    printf("\t  -m <mask> \t| Only estimate statistics within the masked area.\n");
    printf("\t  -t <float> \t| Only estimate statistics if voxel is larger than <float>.\n");
    printf("\n\t  Note: All NaN or Inf are ignored for all stats. \n\t        The -m and -t options can be used in conjusction.\n\n");

    printf("\n\t* * Statistics (at least one option is mandatory) * *\n");
    printf("\tRange operations (datatype: all)\n");
    printf("\t  -r \t\t| The range <min max> of all voxels.\n");
    printf("\t  -R \t\t| The robust range (assuming 2%% outliers on both sides) of all voxels\n");
    printf("\t  -p <float> \t| The <float>th percentile of all voxels intensity (float=[0,100])\n");
    printf("\n\tClassical operations (datatype: all)\n");
    printf("\t  -a  \t\t| Average of all voxels \n");
    printf("\t  -s  \t\t| Standard deviation of all voxels \n");
    printf("\t  -v  \t\t| Volume of all voxels above 0 (<# voxels> * <volume per voxel>)\n");
    printf("\t  -vl \t\t| Volume of each integer label (<# voxels per label> * <volume per voxel>)\n");
    printf("\t  -vp \t\t| Volume of all probabilsitic voxels (sum(<in>) * <volume per voxel>)\n");
    printf("\t  -n  \t\t| Count of all voxels above 0 (<# voxels>)\n");
    printf("\t  -np \t\t| Sum of all fuzzy voxels (sum(<in>))\n");
    printf("\n\tCoordinates operations (datatype: all)\n");
    printf("\t  -x \t\t| Location (in vox) of the smallest value in the image\n");
    printf("\t  -X \t\t| Location (in vox) of the largest value in the image\n");
    printf("\t  -c \t\t| Location (in vox) of the centre of mass of the object\n");
    printf("\t  -B \t\t| Bounding box of all nonzero voxels [ xmin xsize ymin ysize zmin zsize ]\n");
    printf("\n\tLabel attribute operations (datatype: char or uchar)\n");
    printf("\t  -Vl <csv> \t\t| Volume of each integer label <in>. Save to <csv> file.\n");
    printf("\t  -Nl <csv> \t\t| Count of each label <in>. Save to <csv> file.\n");
    printf("\t  -al <in2> \t\t| Average value in <in> for each label in <in2> \n");
    printf("\t  -Al <in2> <csv>\t| Average value in <in> for each label in <in2>. Save to <csv> file\n");
    printf("\t  -d <in2>\t\t| Calculate the Dice score between all classes in <in> and <in2>\n");
    printf("\t  -D <in2> <csv>\t| Calculate the Dice score between all classes in <in> and <in2>. Save to <csv> file\n");
    printf("\t\n");
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

        char * filenames[10];
        nifti_image * Images[10];
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

        bool * mask=new bool [Images[0]->nvox];
        unsigned int maskcount=0;
        if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
        {
            seg_changeDatatype<float>(Images[0]);
        }
        float * Img1prt = static_cast<float *>(Images[0]->data);

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
            if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
                    strcmp(argv[i], "-HELP")==0 || strcmp(argv[i], "-h")==0 ||
                    strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0)
            {
                Usage(argv[0]);
                return 0;
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
                float  Count1[1000]= {0};
                float  Count2[1000]= {0};
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
                float calcvol=0;
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
                float calcvol=0;
                for(unsigned int index=0; index<Images[0]->nvox; index++)
                {
                    if(mask[index])
                    {
                        calcvol += Img1prt[index]>0;
                    }
                }
                cout <<(double)(calcvol)*(double)(Images[0]->dx)*(double)(Images[0]->dy)*(double)(Images[0]->dz)<<endl;
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
                int oldnumbimg=numbimg;
                numbimg=numbimg+1;
                filenames[oldnumbimg] = argv[++i];
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
                myfile.open(filenames[oldnumbimg]);

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
                int oldnumbimg=numbimg;
                numbimg=numbimg+1;
                filenames[oldnumbimg] = argv[++i];
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
                myfile.open(filenames[oldnumbimg]);

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
                float * Img1prt = static_cast<float *>(Images[0]->data);
                float count=0;
                float locX=0;
                float locY=0;
                float locZ=0;
                int index=0;
                for(int Zindex=0; Zindex<Images[0]->nz; Zindex++)
                {
                    for(int Yindex=0; Yindex<Images[0]->ny; Yindex++)
                    {
                        for(int Xindex=0; Xindex<Images[0]->nx; Xindex++)
                        {
                            if(mask[index] && Img1prt[index]>0)
                            {
                                count++;
                                locX+=Xindex;
                                locY+=Yindex;
                                locZ+=Zindex;
                            }
                            index++;
                        }
                    }
                }

                cout <<locX/count<<" "<<locY/count<<" "<<locZ/count<<endl;
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
                int locX=0;
                int locY=0;
                int locZ=0;
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
                                locX=Xindex;
                                locY=Yindex;
                                locZ=Zindex;
                            }
                            index++;
                        }
                    }
                }

                cout <<locX<<" "<<locY<<" "<<locZ<<endl;
                flush(cout);
            }

            // **************************            ---------          *****************************
            // **************************            Min location       *****************************
            // **************************            ---------          *****************************
            else if(strcmp(argv[i], "-X")==0 && (i)<argc)
            {
                if(Images[0]->datatype!=NIFTI_TYPE_FLOAT32)
                {
                    seg_changeDatatype<float>(Images[0]);
                }
                float * Img1prt = static_cast<float *>(Images[0]->data);
                float minval=1e32;
                int locX=0;
                int locY=0;
                int locZ=0;
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
                                locX=Xindex;
                                locY=Yindex;
                                locZ=Zindex;
                            }
                            index++;
                        }
                    }
                }

                cout <<locX<<" "<<locY<<" "<<locZ<<endl;
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
                cout << imgsort[(int)(round(percentile*(maskcount-1)))]<<endl;

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
                float mean=(double)(calc)/(double)(calccount);
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

