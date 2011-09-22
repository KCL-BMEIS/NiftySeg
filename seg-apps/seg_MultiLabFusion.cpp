#include "_seg_MultiLabFusion.h"
#include <iostream>
#include <time.h>
#include <stdlib.h>

using namespace std;
#define PrecisionTYPE float


void Usage(char *exec)
{

    printf("\n  Multi Label Fusion (integer labels):\n  Usage:\t%s -in <filename> -<Type of Label Fusion> [OPTIONS]\n\n",exec);
    printf("  * * * * * * * * * * * * * * * * * * Mandatory * * * * * * * * * * * * * * * * *\n\n");
    printf("  -in <filename>\t\t\tFilename of the 4D integer label image\n\n");
    printf("  \t\t- Type of Classifier Fusion (mutualy exclusive) -\n\n");
//    printf("  -STEPS <k> <n> <img> <tmpl> \tSTEPS algorithm\n");
    printf("  \t\t\t\tSize of the kernel (k), number of local labels to use (n),\n");
    printf("  \t\t\t\tOriginal image to segment (3D Image), registered templates (4D Image).\n");
    printf("  -STAPLE \t\t\tSTAPLE algorithm\n");
    printf("  -MV \t\t\t\tMajority Vote algorithm\n\n");
    printf("  -SBA \t\t\t\tShape Based Averaging algorithm\n\n");

    printf("  * * * * * * * * * * * * * * * * * * Options * * * * * * * * * * * * * * * * * *\n\n");
    printf("  -v <int>\t\t\tVerbose level [0 = off, 1 = on, 2 = debug] (default = 0)\n");
//    printf("  -unc <int>\t\t\tOnly consider non-consensus voxels\n");
    printf("  -out <filename>\t\tFilename of the segmented image (default=LabFusion.nii.gz)\n");
    //printf("  -thr <float> \t\t\tThreshold the final result at a different level (only with binary lables) \n\n");

    printf("  * * * * * * * * * * * * * * STAPLE and STEPS options * * * * * * * * * * * * * *\n\n");
    printf("  -prop <proportion> \t\tProportion of the classifier (automaticaly estimated by default)\n");
    printf("  -prop_update \t\t\tUpdate proportions at each iteration.\n");
    printf("  -setPQ <P> <Q> \t\tValue of P and Q [ 0 < (P,Q) < 1 ] (default = 0.99 0.99) \n");
    printf("  -MRF_beta <float>\t\tMRF prior strength [ 0 < beta < 5 ] \n");
    printf("  -max_iter <int>\t\tMaximum number of iterations (default = 50)\n");
    printf("  -conv <float>\t\t\tRatio for convergence (default epsilon = 10^-6)\n\n");

    printf("  * * * * * * * * Ranking for STAPLE and MV (mutualy exclusive) * * * * * * * * *\n\n");
    printf("  -ALL  (default)\t\tUse all labels with no Ranking\n");
    printf("  -GNCC <n> <img> <tmpl>\tGlobal Normalized Cross Correlation Ranking (Calculated in the full image):\n");
    printf("  \t\t\t\tNumber of sorted classifiers to use (n),\n");
    printf("  \t\t\t\tOriginal image to segment (3D image), registered templates (4D Image).\n");
    printf("  -ROINCC <d> <n> <img> <tmpl>\tROI Normalized Cross Correlation Ranking (On the registered label ROI):\n");
    printf("  \t\t\t\tDilation of the ROI ( <int> d>=1 ), Number of sorted classifiers to use (n),\n");
    printf("  \t\t\t\tOriginal image to segment (3D image), registered templates (4D Image).\n");
    printf("  -LNCC <k> <n> <img> <tmpl>\tLocally Normalized Cross Correlation Ranking (On a local gaussian kernel):\n");
    printf("  \t\t\t\tSize of the kernel (k), number of local classifiers to use (n),\n");
    printf("  \t\t\t\tOriginal image to segment (3D Image), registered templates (4D Image).\n\n");
    printf("  \t\t\t\tLNCC is only availabel for STAPLE and MV.\n\n");

    printf("  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n\n");
    return;
}

int main(int argc, char **argv)
{



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
    int maxIteration=50;
    int verbose_level=0;
    bool PropUpdate=false;
    float tmpP=0;
    float tmpQ=0;
    float conv=0.000001;

    float LNCC_kernel=3;

    int Numb_Neigh=1;


    // STAPLE - 0
    // MV - 1
    // SBA - 2
    int LabFusType=10;
    bool LNCCflag=0;
    bool GNCCflag=0;
    bool ROINCCflag=0;
    bool UNCERTAINflag=0;
    char * filename_LNCC=NULL;
    char * filename_GNCC=NULL;
    char * filename_ROINCC=NULL;
    char * filename_BaseImage=NULL;
    int DilSize=0;

    float Thresh_IMG_value;
    bool Do_Thresh_IMG=false;
    /* read the input parameter */

    for(int i=1;i<argc;i++){
        if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
           strcmp(argv[i], "-HELP")==0 || strcmp(argv[i], "-h")==0 ||
           strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0 ||
           strcmp(argv[i], "--HELP")==0){
            Usage(argv[0]);
            return 0;
        }

        else if(strcmp(argv[i], "-in") == 0){
            filename_LABELS = argv[++i];
        }
        else if(strcmp(argv[i], "-STEPS") == 0){
            if(LabFusType==10){
                LabFusType = 0;
                Thresh_IMG_value=0.999;
                UNCERTAINflag = true;
                LNCC_kernel= atof(argv[++i]);
                Numb_Neigh = atoi(argv[++i]);
                filename_BaseImage= argv[++i];
                filename_LNCC = argv[++i];
                LNCCflag=1;
            }
            else{
                fprintf(stderr,"* Error: Multiple fusion algorithms selected");
                return 1;
            }
        }
        else if(strcmp(argv[i], "-unc") == 0){
            UNCERTAINflag = true;
        }
        else if(strcmp(argv[i], "-STAPLE") == 0){
            if(LabFusType==10){
                LabFusType = 1;
                Thresh_IMG_value=0.999;
            }
            else{
                fprintf(stderr,"* Error: Multiple fusion algorithms selected");
                return 1;
            }
        }
        else if(strcmp(argv[i], "-MV") == 0){
            if(LabFusType==10){
                LabFusType = 2;
                Thresh_IMG_value=0.5;
            }
            else{
                fprintf(stderr,"* Error: Multiple fusion algorithms selected\n");
                return 1;
            }

        }
        else if(strcmp(argv[i], "-SBA") == 0){
            if(LabFusType==10){
                LabFusType = 3;
                Thresh_IMG_value=0.0;
            }
            else{
                fprintf(stderr,"* Error: Multiple fusion algorithms selected\n");
                return 1;
            }

        }
        else if(strcmp(argv[i], "-ALL") == 0){
            if(!UNCERTAINflag){
                LNCCflag = false;
                GNCCflag=false;
                ROINCCflag=false;
            }
            else{
                fprintf(stderr,"* Error: Type of ranking not allowed in STEPS\n");
                return 1;
            }
        }

        else if(strcmp(argv[i], "-LNCC") == 0){
            if(!UNCERTAINflag){
                LNCC_kernel= atof(argv[++i]);
                Numb_Neigh = atoi(argv[++i]);
                filename_BaseImage= argv[++i];
                filename_LNCC = argv[++i];
                LNCCflag=1;
            }
            else{
                fprintf(stderr,"* Error: Type of ranking not allowed in STEPS\n");
                return 1;
            }
        }
        else if(strcmp(argv[i], "-GNCC") == 0){
            if(!UNCERTAINflag){
                Numb_Neigh = atoi(argv[++i]);
                filename_BaseImage= argv[++i];
                filename_GNCC = argv[++i];
                GNCCflag=1;
            }
            else{
                fprintf(stderr,"* Error: Type of ranking not allowed in STEPS\n");
                return 1;
            }
        }
        else if(strcmp(argv[i], "-ROINCC") == 0){
            if(!UNCERTAINflag){
                DilSize = atoi(argv[++i]);
                Numb_Neigh = atoi(argv[++i]);
                filename_BaseImage= argv[++i];
                filename_ROINCC = argv[++i];
                ROINCCflag=1;
            }
            else{
                fprintf(stderr,"* Error: Type of ranking not allowed in STEPS\n");
                return 1;
            }

        }
        else if(strcmp(argv[i], "-out") == 0){
            filename_OUT = argv[++i];
        }
        else if(strcmp(argv[i], "-conv") == 0){
            conv = atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-prop_update") == 0){
            PropUpdate = true;
        }
        else if(strcmp(argv[i], "-setPQ") == 0){
            tmpP = atof(argv[++i]);
            tmpQ = atof(argv[++i]);
        }

        else if(strcmp(argv[i], "-prop") == 0){
            prop=atof(argv[++i]);
            propflag=true;
        }
        else if(strcmp(argv[i], "-v") == 0){
            verbose_level=(int)atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-MRF_beta") == 0){
            MRF_strength=(PrecisionTYPE)atof(argv[++i]);
            if(MRF_strength>4){
                cout << "WARNING: MRF Beta strenght should be less than 4. Setting MRF_beta=4."<< endl;
                MRF_strength=4;

            }
            else if(MRF_strength<0){
                cout << "WARNING: MRF Beta strenght should be more than 0. MRF_beta will be off."<< endl;
                MRF_strength=0;

            }
        }
        else if(strcmp(argv[i], "-max_iter") == 0){
            maxIteration=atoi(argv[++i]);
        }
        else{
            fprintf(stderr,"Err:\tParameter %s unknown->\n",argv[i]);
            Usage(argv[0]);
            return 1;
        }
    }

    // READING Labels
    nifti_image * CLASSIFIER=nifti_image_read(filename_LABELS,1);
    if(CLASSIFIER == NULL){
        fprintf(stderr,"* Error when reading the Label image: %s\n",filename_LABELS);
        return 1;
    }

    classifier_datatype MaxLab=0;
    seg_changeDatatype<classifier_datatype>(CLASSIFIER);

    classifier_datatype * CLASSIFIERptr = static_cast<classifier_datatype *>(CLASSIFIER->data);

    int * NumberOfDifferentClassesHistogram=new int [10000];
    for(int i=0;i<10000;i++){
        NumberOfDifferentClassesHistogram[i]=0;
    }
    for(int i=0;i<(int)CLASSIFIER->nvox;i++){
        NumberOfDifferentClassesHistogram[CLASSIFIERptr[i]]++;
    }
    for(int i=0;i<10000;i++){
        MaxLab+=(NumberOfDifferentClassesHistogram[i]>0)?1:0;
    }
    delete [] NumberOfDifferentClassesHistogram;

    if(verbose_level>0){
        cout << "Merging "<<(int)MaxLab<<" different labels from "<<CLASSIFIER->nt+1<<" classifiers";
        if(LabFusType==0)cout<<" using STEPS";
        if(LabFusType==1)cout<<" using STAPLE";
        if(LabFusType==2)cout<<" using MV";
        if(LabFusType==3)cout<<" using SBA";
        if(LNCCflag)cout<<" ranked by LNCC" <<endl;
        if(GNCCflag)cout<<" ranked by GNCC" <<endl;
        if(ROINCCflag)cout<<" ranked by ROINCC" <<endl;
        flush(cout);
    }
    nifti_image * LNCC=NULL;
    nifti_image * BaseImage=NULL;
    if(LNCCflag){
        // READING LNCC images
        LNCC=nifti_image_read(filename_LNCC,1);
        BaseImage=nifti_image_read(filename_BaseImage,1);

        seg_changeDatatype<float>(BaseImage);
        seg_changeDatatype<float>(LNCC);
        if(!UNCERTAINflag){
            if(LNCC->nt!=CLASSIFIER->nt){
                cout << "Number of lables in the -in image is different from the number of registered templates in -STEPS.";
                return 1;
            }
        }
        else{
            if(LNCC->nt!=CLASSIFIER->nt){
                cout << "Number of lables in the -in image is different from the number of registered templates in -LNCC.";
                return 1;
            }
        }
    }

    //

    nifti_image * GNCC=NULL;
    if(GNCCflag){
        // READING LNCC images
        GNCC=nifti_image_read(filename_GNCC,true);
        BaseImage=nifti_image_read(filename_BaseImage,true);

        seg_changeDatatype<float>(BaseImage);
        seg_changeDatatype<float>(GNCC);
        if(GNCC->nt!=CLASSIFIER->nt){
            cout << "Number of lables in the -in image is different from the number of registered templates in -GNCC.";
            return 1;
        }
    }
    nifti_image * ROINCC=NULL;
    if(ROINCCflag){
        // READING LNCC images
        ROINCC=nifti_image_read(filename_ROINCC,true);
        BaseImage=nifti_image_read(filename_BaseImage,true);

        seg_changeDatatype<float>(BaseImage);
        seg_changeDatatype<float>(ROINCC);
        if(ROINCC->nt!=CLASSIFIER->nt){
            cout << "Number of lables in the -in image is different from the number of registered templates in -ROINCC.";
            return 1;
        }
    }

    time_t start,end;
    time(&start);

    seg_LabFusion LabFusion(CLASSIFIER->nt,MaxLab,Numb_Neigh);
    LabFusion.SetVerbose(verbose_level);
    LabFusion.SetinputCLASSIFIER(CLASSIFIER,UNCERTAINflag);

    if(LNCCflag){
        LabFusion.SetLNCC(LNCC,BaseImage,LNCC_kernel,Numb_Neigh);
        nifti_image_free(LNCC);
    }
    if(GNCCflag){
        LabFusion.SetGNCC(GNCC,BaseImage,Numb_Neigh);
        nifti_image_free(GNCC);
    }
    if(ROINCCflag){
        LabFusion.SetROINCC(ROINCC,BaseImage,Numb_Neigh,DilSize);
        nifti_image_free(ROINCC);
    }


    if(filename_OUT!=NULL){LabFusion.SetFilenameOut(filename_OUT);}

    if(Do_Thresh_IMG){
        LabFusion.SetImgThresh(Thresh_IMG_value);
    }

    if (LabFusType == 0 || LabFusType == 1 ){
        if(conv>0){LabFusion.SetConv(conv);}
        if(propflag){LabFusion.SetProp(prop);}
        if(MRF_strength>0.0f){LabFusion.Turn_MRF_ON(MRF_strength);}
        if(PropUpdate>0.0f){LabFusion.Turn_Prop_Update_ON();}
        if(tmpP>0 && tmpQ>0 && tmpP<1 && tmpQ<1){LabFusion.SetPQ(tmpP,tmpQ);}
        LabFusion.SetMaximalIterationNumber(maxIteration);
        LabFusion.Run_STAPLE();
    }
    else if(LabFusType == 2){
        LabFusion.Run_MV();
    }
    else{
        LabFusion.Run_SBA();

    }

    if(verbose_level>0){
        cout << "Saving Fused Label"<<endl;
    }


    nifti_image * Result = LabFusion.GetResult();
    nifti_image_write(Result);




    time(&end);
    if(verbose_level>0){
        cout << "Finished in "<<difftime(end,start)<<"sec"<< endl;
    }

    nifti_image_free(Result);
    nifti_image_free(CLASSIFIER);

    return 0;
}
