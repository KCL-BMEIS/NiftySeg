/**
 * @file seg_LabFusion.cpp
 * @author M. Jorge Cardoso
 * @date 01/01/2014
 *
 * Copyright (c) 2014, University College London. All rights reserved.
 * Centre for Medical Image Computing (CMIC)
 * See the LICENSE.txt file in the nifty_seg root folder
 *
 */

#include "_seg_LabFusion.h"
#include "seg_LabFusion_CLIxml.h"

#include <iostream>
#include <time.h>
#include <stdlib.h>

using namespace std;
#define SegPrecisionTYPE float


void SmallUsage(char *exec)
{
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    printf("Label Fusion:\nUsage ->\t%s -in <filename> -<Type of Label Fusion> [OPTIONS]\n\n",exec);
    printf("\tSee the help for more details (-h).\n");
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    return;
}

void Usage(char *exec)
{

    printf("\nLabel Fusion:\nUsage ->\t%s -in <filename> -<Type of Label Fusion> [OPTIONS]\n\n",exec);
    printf("\t* * * * * * * * * * * * * * * * * * Mandatory * * * * * * * * * * * * * * * * *\n\n");
    printf("\t-in <filename>\t\t\t| Filename of the 4D integer label image\n\n");
    printf("\t\t* * * Type of Classifier Fusion (mutually exclusive) * * *\n\n");
    printf("\t-STEPS <k> <n> <i> <t> \t\t| STEPS algorithm\n");
    printf("\t\t\t\t\t| Size of the kernel (k), number of local labels to use (n),\n");
    printf("\t\t\t\t\t| Original image to segment (3D Image), registered templates (4D Image).\n");
    printf("\t-MLSTEPS <k> <l> <n> <i> <t>\t| Multi-level STEPS algorithm (Beta testing. Do not use!)\n");
    printf("\t\t\t\t\t| Size of the kernel (k), number of levels (l) and local labels to use (n),\n");
    printf("\t\t\t\t\t| Original image to segment (3D Image), registered templates (4D Image).\n");
    printf("\t-STAPLE \t\t\t| STAPLE algorithm\n");
    printf("\t-MV \t\t\t\t| Majority Vote algorithm\n");
    printf("\t-SBA \t\t\t\t| Shape Based Averaging algorithm (Beta)\n\n");

    printf("\t* * * * * * * * * * * * * * * * General Options * * * * * * * * * * * * * * * * *\n\n");
    printf("\t-v <int>\t\t\t| Verbose level [0 = off, 1 = on, 2 = debug] (default = 0)\n");
    printf("\t-unc \t\t\t\t| Only consider non-consensus voxels to calculate statistics\n");
    printf("\t-out <filename>\t\t\t| Filename of the integer segmented image (default=LabFusion.nii.gz)\n");
    printf("\t-mask <filename>\t\t| Filename of the ROI for label fusion (greatly reduces memory requirements)\n");
    printf("\t-outProb \t\t\t| Probabilistic/Fuzzy segmented image (only for 1 label)\n");
#ifdef _GIT_HASH
    printf("\t--version\t\t\t|Print current source code git hash key and exit\n\t\t\t\t\t(%s)\n",_GIT_HASH);
#endif

    printf("\n\t* * * * * * * * * * * * * * STAPLE and STEPS options * * * * * * * * * * * * * *\n\n");
    printf("\t-prop <proportion> \t\t| Proportion of the label (only for single labels)\n");
    printf("\t-prop_update \t\t\t| Update label proportions at each iteration.\n");
    //printf("\t-dil_unc <int> \t\t\t| Dilate uncertainty region by <int>.\n");
    printf("\t-setPQ <P> <Q> \t\t\t| Value of P and Q [ 0 < (P,Q) < 1 ] (default = 0.99 0.99) \n");
    printf("\t-MRF_beta <float>\t\t| MRF prior strength [ 0 < beta < 5 ] \n");
    printf("\t-max_iter <int>\t\t\t| Maximum number of iterations (default = 15)\n");
    printf("\t-uncthres <float>\t\t| If <float> percent of labels agree, then area is not uncertain \n");
    printf("\t-conv <float>\t\t\t| Ratio for convergence (default epsilon = 10^-5)\n\n");

    printf("\t* * * * * * * * Ranking for STAPLE and MV (mutually exclusive) * * * * * * * * *\n\n");
    printf("\t-ALL  (default)\t\t\t| Use all labels with no Ranking\n");
    printf("\t-GNCC <n> <img> <tmpl>\t\t| Global Normalized Cross Correlation Ranking (Calculated in the full image):\n");
    printf("\t\t\t\t\t| Number of sorted classifiers to use (n),\n");
    printf("\t\t\t\t\t| Original image to segment (3D image), registered templates (4D Image).\n");
    printf("\t-ROINCC <d> <n> <img> <tmpl>\t| ROI Normalized Cross Correlation Ranking (On the registered label ROI):\n");
    printf("\t\t\t\t\t| Dilation of the ROI ( <int> d>=1 ), Number of sorted classifiers to use (n),\n");
    printf("\t\t\t\t\t| Original image to segment (3D image), registered templates (4D Image).\n");
    printf("\t-LNCC <k> <n> <img> <tmpl>\t| Locally Normalized Cross Correlation Ranking (On a local gaussian kernel):\n");
    printf("\t\t\t\t\t| Size of the kernel (k), number of local classifiers to use (n),\n");
    printf("\t\t\t\t\t| Original image to segment (3D Image), registered templates (4D Image).\n");
    printf("\t\t\t\t\t| LNCC is only available for STAPLE and MV.\n\n");
    //printf("\t-LM <n> <metric> \t| Any voxelwise local metric (higher metric value is more similar):\n");
    //printf("\t\t\t\t\t| number of local classifiers to use (n), similarity <metric> as a 4D Image.\n");

    printf("\n\t* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n\n");
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


        char * filename_LABELS=NULL;
        char * filename_OUT=NULL;
        float prop=0.0f;
        bool propflag=0;
        float MRF_strength=0;
        int maxIteration=15;
        int verbose_level=0;
        bool PropUpdate=false;
        float tmpP=0;
        float tmpQ=0;
        float conv=0.01;
        int dilunc=0;
        float uncthres=-1;

        float LNCC_kernel=3;

        int Numb_Neigh=1;


        // STAPLE - 0
        // MV - 1
        // SBA - 2
        int LabFusType=10;
        bool LNCCflag=0;
        bool ML_LNCCflag=0;
        //bool LMflag=0;
        bool GNCCflag=0;
        bool ROINCCflag=0;
        bool UNCERTAINflag=0;
        bool UseMask=0;
        char * filename_Mask=NULL;
        int ProbOutput=0;
        char * filename_LNCC=NULL;
        char * filename_GNCC=NULL;
        char * filename_ROINCC=NULL;
        char * filename_BaseImage=NULL;
        int DilSize=0;
        int MLSTEPS_levels=0;

        float Thresh_IMG_value;
        bool Do_Thresh_IMG=false;
        /* read the input parameter */

        for(int i=1; i<argc; i++)
        {
            if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
                    strcmp(argv[i], "-HELP")==0 || strcmp(argv[i], "-h")==0 ||
                    strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0 ||
                    strcmp(argv[i], "--HELP")==0)
            {
                Usage(argv[0]);
                return 0;
            }
            else if( (strcmp(argv[i], "-xml")==0) || (strcmp(argv[i], "--xml")==0) )
            {
                cout << xml_segLabFusion;
                return 0;
            }
            else if(strcmp(argv[i], "-in") == 0 && (i+1)<argc)
            {
                filename_LABELS = argv[++i];
            }
            else if(strcmp(argv[i], "-STEPS") == 0&& (i+4)<argc)
            {
                if(LabFusType==10)
                {
                    //PropUpdate = true;
                    LabFusType = 0;
                    UNCERTAINflag = true;
                    LNCC_kernel= atof(argv[++i]);
                    Numb_Neigh = atoi(argv[++i]);
                    filename_BaseImage= argv[++i];
                    filename_LNCC = argv[++i];
                    LNCCflag=1;
                }
                else
                {
                    fprintf(stderr,"* Error: Multiple fusion algorithms selected");
                    flush(cout);
                    return 1;
                }
            }
            else if(strcmp(argv[i], "-MLSTEPS") == 0&& (i+4)<argc)
            {
                if(LabFusType==10)
                {
                    LabFusType = 0;
                    UNCERTAINflag = true;
                    LNCC_kernel= atof(argv[++i]);
                    MLSTEPS_levels= atoi(argv[++i]);
                    Numb_Neigh = atoi(argv[++i]);
                    filename_BaseImage= argv[++i];
                    filename_LNCC = argv[++i];
                    ML_LNCCflag=1;
                    LNCCflag=1;
                }
                else
                {
                    fprintf(stderr,"* Error: Multiple fusion algorithms selected");
                    flush(cout);
                    return 1;
                }
            }
            else if(strcmp(argv[i], "-unc") == 0)
            {
                UNCERTAINflag = true;
            }
            else if(strcmp(argv[i], "-uncthres") == 0)
            {
                uncthres =atof(argv[++i]);
            }
            else if(strcmp(argv[i], "-mask") == 0)
            {
                UseMask = true;
                filename_Mask=argv[++i];
            }
            else if(strcmp(argv[i], "-STAPLE") == 0)
            {
                if(LabFusType==10)
                {
                    LabFusType = 1;
                    Thresh_IMG_value=0.999;
                }
                else
                {
                    fprintf(stderr,"* Error: Multiple fusion algorithms selected");
                    flush(cout);
                    return 1;
                }
            }
            else if(strcmp(argv[i], "-MV") == 0)
            {
                if(LabFusType==10)
                {
                    LabFusType = 2;
                    Thresh_IMG_value=0.5;
                }
                else
                {
                    fprintf(stderr,"* Error: Multiple fusion algorithms selected\n");
                    flush(cout);
                    return 1;
                }

            }
            else if(strcmp(argv[i], "-SBA") == 0)
            {
                if(LabFusType==10)
                {
                    LabFusType = 3;
                    Thresh_IMG_value=0.0;
                }
                else
                {
                    fprintf(stderr,"* Error: Multiple fusion algorithms selected\n");
                    flush(cout);
                    return 1;
                }

            }
            else if(strcmp(argv[i], "-ALL") == 0)
            {
                if(!UNCERTAINflag)
                {
                    LNCCflag = false;
                    GNCCflag=false;
                    ROINCCflag=false;
                }
                else
                {
                    fprintf(stderr,"* Error: Type of ranking not allowed in STEPS\n");
                    flush(cout);
                    return 1;
                }
            }

            else if(strcmp(argv[i], "-LNCC") == 0 && (i+4)<argc)
            {
                if(!UNCERTAINflag)
                {
                    LNCC_kernel= atof(argv[++i]);
                    Numb_Neigh = atoi(argv[++i]);
                    filename_BaseImage= argv[++i];
                    filename_LNCC = argv[++i];
                    LNCCflag=1;
                }
                else
                {
                    fprintf(stderr,"* Error: Type of ranking not allowed in STEPS\n");
                    flush(cout);
                    return 1;
                }
            }
            //      else if(strcmp(argv[i], "-LM") == 0 && (i+4)<argc){
            //          if(!UNCERTAINflag){
            //              Numb_Neigh = atoi(argv[++i]);
            //              filename_LNCC = argv[++i];
            //              LMflag=1;
            //            }
            //          else{
            //              fprintf(stderr,"* Error: Type of ranking not allowed in STEPS\n");
            //              flush(cout); return 1;
            //            }
            //        }
            else if(strcmp(argv[i], "-GNCC") == 0 && (i+3)<argc)
            {
                if(!UNCERTAINflag)
                {
                    Numb_Neigh = atoi(argv[++i]);
                    filename_BaseImage= argv[++i];
                    filename_GNCC = argv[++i];
                    GNCCflag=1;
                }
                else
                {
                    fprintf(stderr,"* Error: Type of ranking not allowed in STEPS\n");
                    flush(cout);
                    return 1;
                }
            }
            else if(strcmp(argv[i], "-ROINCC") == 0 && (i+4)<argc)
            {
                if(!UNCERTAINflag)
                {
                    DilSize = atoi(argv[++i]);
                    Numb_Neigh = atoi(argv[++i]);
                    filename_BaseImage= argv[++i];
                    filename_ROINCC = argv[++i];
                    ROINCCflag=1;
                }
                else
                {
                    fprintf(stderr,"* Error: Type of ranking not allowed in STEPS\n");
                    flush(cout);
                    return 1;
                }

            }
            else if(strcmp(argv[i], "-out") == 0)
            {
                filename_OUT = argv[++i];
            }
            else if(strcmp(argv[i], "-outProb") == 0)
            {
                ProbOutput=1;
            }
            /*else if(strcmp(argv[i], "-dil_unc") == 0){
                dilunc = atoi(argv[++i]);
              }*/
            else if(strcmp(argv[i], "-conv") == 0)
            {
                conv = atof(argv[++i]);
            }
            else if(strcmp(argv[i], "-prop_update") == 0)
            {
                PropUpdate = true;
            }
            else if(strcmp(argv[i], "-setPQ") == 0)
            {
                tmpP = atof(argv[++i]);
                tmpQ = atof(argv[++i]);
            }

            else if(strcmp(argv[i], "-prop") == 0)
            {
                prop=atof(argv[++i]);
                propflag=true;
            }
            else if(strcmp(argv[i], "-v") == 0 && (i+1)<argc)
            {
                verbose_level=(int)atoi(argv[++i]);
            }
            else if(strcmp(argv[i], "-MRF_beta") == 0 && (i+1)<argc)
            {
                MRF_strength=(SegPrecisionTYPE)atof(argv[++i]);
                if(MRF_strength>4)
                {
                    cout << "WARNING: MRF Beta strenght should be less than 4. Setting MRF_beta=4."<< endl;
                    MRF_strength=4;

                }
                else if(MRF_strength<0)
                {
                    cout << "WARNING: MRF Beta strenght should be more than 0. MRF_beta will be off."<< endl;
                    MRF_strength=0;

                }
            }
            else if(strcmp(argv[i], "-max_iter") == 0 && (i+1)<argc)
            {
                maxIteration=atoi(argv[++i]);
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
                printf("Err:\tParameter %s unknown or incomplete\n\n",argv[i]);
                flush(cout);
                SmallUsage(argv[0]);
                flush(cout);
                return 1;
            }
        }

        // READING Labels
        nifti_image * CLASSIFIER=nifti_image_read(filename_LABELS,1);
        if(CLASSIFIER == NULL || CLASSIFIER->data==NULL)
        {
            fprintf(stderr,"* Error when reading the Label image: %s\n",filename_LABELS);
            flush(cout);
            return 1;
        }
        categoricalLabelType MaxLab=0;
        seg_changeDatatype<categoricalLabelType>(CLASSIFIER);


        categoricalLabelType * CLASSIFIERptr = static_cast<categoricalLabelType *>(CLASSIFIER->data);

        int * NumberOfDifferentClassesHistogram=new int [10000];
        for(int i=0; i<10000; i++)
        {
            NumberOfDifferentClassesHistogram[i]=0;
        }
        for(int i=0; i<(int)CLASSIFIER->nvox; i++)
        {
            NumberOfDifferentClassesHistogram[CLASSIFIERptr[i]]++;
        }
        for(int i=0; i<10000; i++)
        {
            MaxLab+=(NumberOfDifferentClassesHistogram[i]>0)?1:0;
        }
        delete [] NumberOfDifferentClassesHistogram;



        nifti_image * Mask=NULL;
        if(UseMask)
        {
            Mask=nifti_image_read(filename_Mask,1);
            if(Mask == NULL)
            {
                fprintf(stderr,"* Error when reading the Label image: %s\n",filename_Mask);
                flush(cout);
                return 1;
            }
            seg_convert2binary(Mask,0);
        }


        if(LabFusType==10)
        {
            LabFusType = 1;
            Thresh_IMG_value=0.999;
        }

        // *****************************
        //   CALCULATING REQUIRED SIZE
        // *****************************
        if(verbose_level>0)
        {
            cout << "Merging "<<(int)MaxLab<<" different labels from "<<CLASSIFIER->nt<<" classifiers";
            if(LabFusType==0)
            {
                if(ML_LNCCflag==0)
                {
                    cout<<" using STEPS";
                }
                else
                {
                    cout<<" using Multi Level STEPS";
                }
            }
            if(LabFusType==1)cout<<" using STAPLE";
            if(LabFusType==2)cout<<" using MV";
            if(LabFusType==3)cout<<" using SBA";
            if(LabFusType!=0)
            {
                if(LNCCflag)cout<<" ranked by LNCC";
                if(GNCCflag)cout<<" ranked by GNCC";
                if(ROINCCflag)cout<<" ranked by ROINCC";
            }
            flush(cout);
            cout << endl;
            flush(cout);
            float Number_Atlases=(float)CLASSIFIER->nt;
            float Number_Labels=(float)MaxLab;
            float tempImgSize=CLASSIFIER->nx*CLASSIFIER->ny*CLASSIFIER->nz;


            int maskcount=(int)tempImgSize;
            if(UseMask)
            {
                maskcount=0;
                bool * inputMaskPtr=static_cast<bool *>(Mask->data);
                for(int i=0; i<tempImgSize; i++)
                {
                    maskcount+=(inputMaskPtr[i]>0)?1:0;
                }
            }


            float sizeLabelImage=Number_Atlases*tempImgSize;
            float sizeTargetImage=(LNCCflag==false && GNCCflag==false &&ROINCCflag==false)?0:(tempImgSize*4);
            float sizeAtlasImage=(LNCCflag==false && GNCCflag==false &&ROINCCflag==false)?0:(tempImgSize*Number_Atlases*4);
            float sizeLNCC=(LNCCflag==true)?(tempImgSize*((float)Numb_Neigh+3.0f)*4):0;
            float sizeProbabilities=(LabFusType==2)?(tempImgSize):((tempImgSize+Number_Labels*maskcount)*4);
            float sizeMRF=(MRF_strength==0)?0:(Number_Labels*tempImgSize*4);

            float sizeConfMatrix=Number_Labels*Number_Labels*Number_Atlases*4;
            float sizeMrfMatrix=Number_Labels*Number_Atlases*4;
            float sizeresult_Image=(ProbOutput==0)?(tempImgSize*4):(tempImgSize*4);
            float sizeSum=sizeLabelImage + sizeLNCC+ sizeTargetImage + sizeAtlasImage + sizeProbabilities + sizeMRF + sizeConfMatrix + sizeMrfMatrix +sizeresult_Image;

            cout <<"Will require less than "<<setprecision (3) << (float)(1.1*(sizeSum))/powf(1024.0f,3) <<"Gb of memory"<<endl;
            if(LabFusType==0)
            {
                cout<<"Possibly much less, depending on how uncertain your labels are"<<endl<<endl;
            }
            flush(cout);
        }
        cout<<(float)MaxLab<<" "<<LabFusType<<endl;
        if(ProbOutput==1 && (MaxLab>9))
        {
            fprintf(stderr,"* Due to memory limitations, Probabilistic output only available for less than 10 labels\n");
            flush(cout);
            return 1;
        }

        nifti_image * LNCC=NULL;
        nifti_image * BaseImage=NULL;
        if(LNCCflag || ML_LNCCflag)
        {
            // READING LNCC images

            if(verbose_level>1)cout << "Read LNCC";
            LNCC=nifti_image_read(filename_LNCC,1);
            if(verbose_level>1)
            {
                cout << " - done"<<endl;
                flush(cout);
            }
            if(LNCC == NULL)
            {
                fprintf(stderr,"* Error when reading the LNCC image: %s\n",filename_LNCC);
                flush(cout);
                return 1;
            }



            if(verbose_level>1)cout << "Read BaseImage";
            BaseImage=nifti_image_read(filename_BaseImage,1);
            if(verbose_level>1)
            {
                cout << " - done"<<endl;
                flush(cout);
            }
            if(BaseImage == NULL)
            {
                fprintf(stderr,"* Error when reading the BaseImage image: %s\n",filename_BaseImage);
                flush(cout);
                return 1;
            }


            if(verbose_level>1)cout << "seg_changeDatatype BaseImage";
            seg_changeDatatype<float>(BaseImage);
            if(verbose_level>1)
            {
                cout << " - done"<<endl;
                flush(cout);
            }

            if(verbose_level>1)cout << "seg_changeDatatype LNCC";
            seg_changeDatatype<float>(LNCC);
            if(verbose_level>1)
            {
                cout << " - done"<<endl;
                flush(cout);
            }


            if(!UNCERTAINflag)
            {
                if(LNCC->nt!=CLASSIFIER->nt)
                {
                    cout << "Number of lables in the -in image is different from the number of registered templates in -STEPS.";
                    flush(cout);
                    return 1;
                }
            }
            else
            {
                if(LNCC->nt!=CLASSIFIER->nt)
                {
                    cout << "Number of lables in the -in image is different from the number of registered templates in -LNCC.";
                    flush(cout);
                    return 1;
                }
            }
        }

        //  nifti_image * METRIC=NULL;
        //  if(LMflag){
        //      // READING Metric image
        //      if(verbose_level>1)cout << "Read LNCC";
        //      METRIC=nifti_image_read(filename_LNCC,1);
        //      if(verbose_level>1){cout << " - done"<<endl;flush(cout);}
        //      if(METRIC == NULL){
        //          fprintf(stderr,"* Error when reading the metric image: %s\n",filename_LNCC);
        //          flush(cout); return 1;
        //        }
        //      if(verbose_level>1)cout << "seg_changeDatatype METRIC";
        //      if(METRIC->datatype!=NIFTI_TYPE_FLOAT32)
        //          seg_changeDatatype<float>(METRIC);
        //      if(verbose_level>1)
        //        cout << " - done"<<endl;flush(cout);
        //      if(METRIC->nt!=CLASSIFIER->nt){
        //          cout << "Number of lables in the -in image is different from the number of images in the 4d metric image.";
        //          flush(cout); return 1;
        //        }

        //    }
        //

        nifti_image * GNCC=NULL;
        if(GNCCflag)
        {
            // READING LNCC images
            GNCC=nifti_image_read(filename_GNCC,true);
            if(GNCC == NULL)
            {
                fprintf(stderr,"* Error when reading the GNCC image: %s\n",filename_GNCC);
                flush(cout);
                return 1;
            }
            if(verbose_level>1)
            {
                cout << "Read GNCC done"<<endl;
                flush(cout);
            }
            BaseImage=nifti_image_read(filename_BaseImage,true);
            if(BaseImage == NULL)
            {
                fprintf(stderr,"* Error when reading the BaseImage image: %s\n",filename_BaseImage);
                flush(cout);
                return 1;
            }
            if(verbose_level>1)
            {
                cout << "Read BaseImage done"<<endl;
                flush(cout);
            }

            seg_changeDatatype<float>(BaseImage);
            seg_changeDatatype<float>(GNCC);
            if(verbose_level>1)
            {
                cout << "seg_changeDatatype GNCC done"<<endl;
                flush(cout);
            }
            seg_changeDatatype<float>(BaseImage);

            if(verbose_level>1)
            {
                cout << "seg_changeDatatype BaseImage done"<<endl;
                flush(cout);
            }
            if(GNCC->nt!=CLASSIFIER->nt)
            {
                cout << "Number of lables in the -in image is different from the number of registered templates in -GNCC.";
                flush(cout);
                return 1;
            }
        }
        nifti_image * ROINCC=NULL;
        if(ROINCCflag)
        {
            // READING LNCC images
            ROINCC=nifti_image_read(filename_ROINCC,true);
            if(ROINCC == NULL)
            {
                fprintf(stderr,"* Error when reading the ROINCC image: %s\n",filename_ROINCC);
                flush(cout);
                return 1;
            }
            BaseImage=nifti_image_read(filename_BaseImage,true);
            if(BaseImage == NULL)
            {
                fprintf(stderr,"* Error when reading the BaseImage image: %s\n",filename_BaseImage);
                flush(cout);
                return 1;
            }
            seg_changeDatatype<float>(BaseImage);
            seg_changeDatatype<float>(ROINCC);
            if(ROINCC->nt!=CLASSIFIER->nt)
            {
                cout << "Number of lables in the -in image is different from the number of registered templates in -ROINCC.";
                flush(cout);
                return 1;
            }
        }

        time_t start,end;
        time(&start);
        if(verbose_level>1)
        {
            cout << endl<<"Creating Object";
        }
        seg_LabFusion LabFusion(CLASSIFIER->nt,MaxLab,Numb_Neigh,LNCC!=NULL?std::max(LNCC->nu,1):1);

        if(verbose_level>1)
        {
            cout << " - Done"<<endl;
            flush(cout);
        }
        if(verbose_level>1)
        {
            cout << "Setting Verbose";
        }
        LabFusion.SetVerbose(verbose_level);
        if(verbose_level>1)
        {
            cout << " - Done"<<endl;
            flush(cout);
        }
        if(dilunc>0)
        {
            LabFusion.SetDilUnc(dilunc);
        }

        if(verbose_level>1)
        {
            cout << "Setting Input";
        }
        LabFusion.SetinputCLASSIFIER(CLASSIFIER,UNCERTAINflag);
        if(verbose_level>1)
        {
            cout << " - Done"<<endl;
            flush(cout);
        }

        if(UseMask)
        {
            LabFusion.SetMask(Mask);
        }
        nifti_image_free(Mask);

        if(LNCCflag)
        {
            if(ML_LNCCflag)
            {
                if(verbose_level>1)
                {
                    cout << "Calling SetMLLNCC"<<endl;
                    flush(cout);
                }
                LabFusion.SetMLLNCC(LNCC,BaseImage,LNCC_kernel,MLSTEPS_levels,Numb_Neigh);
                nifti_image_free(LNCC);
            }
            else
            {
                if(verbose_level>1)
                {
                    cout << "Calling SetLNCC"<<endl;
                    flush(cout);
                }
                LabFusion.SetLNCC(LNCC,BaseImage,LNCC_kernel,Numb_Neigh);
                nifti_image_free(LNCC);

            }
        }

        if(GNCCflag)
        {
            if(verbose_level>1)
            {
                cout << "Calling SetGNCC"<<endl;
                flush(cout);
            }
            LabFusion.SetGNCC(GNCC,BaseImage,Numb_Neigh);
            nifti_image_free(GNCC);
        }
        if(ROINCCflag)
        {
            if(verbose_level>1)
            {
                cout << "Calling SetROINCC"<<endl;
                flush(cout);
            }
            LabFusion.SetROINCC(ROINCC,BaseImage,Numb_Neigh,DilSize);
            nifti_image_free(ROINCC);
        }
        //  if(LMflag){
        //      LabFusion.SetLMETRIC(METRIC,Numb_Neigh);
        //      nifti_image_free(METRIC);
        //    }


        if(filename_OUT!=NULL)
        {
            LabFusion.SetFilenameOut(filename_OUT);
        }

        if(Do_Thresh_IMG)
        {
            LabFusion.SetImgThresh(Thresh_IMG_value);
        }

        if (LabFusType == 0 || LabFusType == 1 )
        {
            if(conv>0)
            {
                LabFusion.SetConv(conv);
            }
            if(propflag)
            {
                LabFusion.SetProp(prop);
            }
            if(MRF_strength>0.0f)
            {
                LabFusion.Turn_MRF_ON(MRF_strength);
            }
            if(PropUpdate>0.0f)
            {
                LabFusion.Turn_Prop_Update_ON();
            }
            if(uncthres>0.5)
            {
                LabFusion.SetUncThresh(uncthres);
            }
            if(tmpP>0 && tmpQ>0 && tmpP<1 && tmpQ<1)
            {
                LabFusion.SetPQ(tmpP,tmpQ);
            }
            LabFusion.SetMaximalIterationNumber(maxIteration);
            LabFusion.Run_STAPLE_or_STEPS();
        }
        else if(LabFusType == 2)
        {
            LabFusion.Run_MV();
        }
        else
        {
            LabFusion.Run_SBA();

        }

        if(ProbOutput==0)
        {
            nifti_image * Result = LabFusion.GetResult_label();
            nifti_image_write(Result);
            nifti_image_free(Result);
        }
        else
        {
            nifti_image * Result = LabFusion.GetResult_probability();
            nifti_image_write(Result);
            nifti_image_free(Result);
        }






        time(&end);
        if(verbose_level>0)
        {
            cout << "Finished in "<<difftime(end,start)<<"sec"<< endl;
        }


        nifti_image_free(CLASSIFIER);
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
