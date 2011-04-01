#include "_seg_STAPLE.h"
#include <iostream>
#include <time.h>
#include <stdlib.h>

using namespace std;
#define PrecisionTYPE float


void Usage(char *exec)
{
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    printf("Usage:\t%s -in <filename> [OPTIONS].\n\n",exec);
    printf("\t* * Mandatory * *\n");
    printf("\t-in <filename>\t\t\tFilename of the lable image\n\n");

    printf("\t* * Model Options * *\n");
    printf("\t-out <filename>\t\t\tFilename of the segmented image (default=STAPLE.nii.gz)\n");
    printf("\t-prop <proportion> \tProportion of the lable\n");
    printf("\t-prop_update \tDynamically update proportions\n");
    printf("\t-setPQ <P> <Q> \tValue of P and Q [ 0 < (P,Q) < 1 ] (default = 0.99 0.99) \n");
    printf("\t-mrf_beta <float>\t\tMRF prior strength [off = 0] (default = 0) \n");

    printf("\t* * General Options * *\n");
    printf("\t-max_iter <int>\t\t\tMaximum number of iterations (default = 100)\n");
    printf("\t-conv <float>\t\tRatio for convergence (default = 0)\n");
    printf("\t-v <int>\t\t\tVerbose level [0 = off, 1 = on, 2 = debug] (default = 0)\n");

    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    return;
}

int main(int argc, char **argv)
{



    if (argc < 2)
    {
        Usage(argv[0]);
        return 0;
    }


    char * filename_LABLES=NULL;
    char * filename_OUT=NULL;
    float prop=0.0f;
    bool propflag=0;
    float MRF_strength=0;
    int maxIteration=100;
    int verbose_level=0;
    bool PropUpdate=false;
    float tmpP=0;
    float tmpQ=0;
    float conv=0;
    int LNCCdistance=3;
    bool saveLNCC=0;
    float lnccthresh=0.05;

    char * filename_PRIOR=NULL;
    bool LoadPrior=false;


    bool LNCCflag=false;
    char * filename_LNCC=NULL;
    char * filename_BaseImage=NULL;
    /* read the input parameter */

    for(int i=1;i<argc;i++){
        if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
           strcmp(argv[i], "-HELP")==0 || strcmp(argv[i], "-h")==0 ||
           strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0){
            Usage(argv[0]);
            return 0;
        }

        else if(strcmp(argv[i], "-in") == 0){
            filename_LABLES = argv[++i];
        }
        else if(strcmp(argv[i], "-lncc") == 0){
            LNCCdistance= atoi(argv[++i]);
            filename_BaseImage= argv[++i];
            filename_LNCC = argv[++i];
            LNCCflag=true;
        }
        else if(strcmp(argv[i], "-savelncc") == 0){
            saveLNCC=1;
        }
        else if(strcmp(argv[i], "-lnccthresh") == 0){
            lnccthresh = atof(argv[++i]);
        }

        else if(strcmp(argv[i], "-prior") == 0){
            filename_PRIOR= argv[++i];
            LoadPrior=true;
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
        else if(strcmp(argv[i], "-mrf_beta") == 0){
            MRF_strength=(PrecisionTYPE)atof(argv[++i]);
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

    // READING Lables
    nifti_image * LABLE=nifti_image_read(filename_LABLES,true);
    if(LABLE == NULL){
        fprintf(stderr,"* Error when reading the Lable image: %s\n",filename_LABLES);
        return 1;
    }
    seg_changeDatatype<unsigned char>(LABLE);

    nifti_image * LNCC=NULL;
    nifti_image * BaseImage=NULL;
    if(LNCCflag){
        // READING LNCC images
        LNCC=nifti_image_read(filename_LNCC,true);
        BaseImage=nifti_image_read(filename_BaseImage,true);

        seg_changeDatatype<float>(BaseImage);
        seg_changeDatatype<float>(LNCC);
    }

    nifti_image * Prior=NULL;
    if(LoadPrior){
        // READING LNCC images
        Prior=nifti_image_read(filename_PRIOR,true);
    }

    if(LABLE->nt>MaxSTAPLElable){
        cout << "There is a predefined maximum of"<<MaxSTAPLElable<<" lables. Please change the ''#define MaxSTAPLElable 40'' in _seg_common.h to a higher number and recompile. "<< endl;
        return 0;
    }

    seg_STAPLE STAPLE(LABLE->nt);
    STAPLE.SetVerbose(verbose_level);
    STAPLE.SetInputLables(LABLE);
    if(LoadPrior){
        STAPLE.SetPrior(Prior,lnccthresh);
    }
    else{
        if(LNCCflag){STAPLE.SetLNCC(LNCC,BaseImage,LNCCdistance,lnccthresh,saveLNCC);}
    }
    if(conv>0){STAPLE.SetConv(conv);}
    if(propflag){STAPLE.SetProp(prop);}
    if(filename_OUT!=NULL){STAPLE.SetFilenameOut(filename_OUT);}
    if(MRF_strength>0.0f){STAPLE.Turn_MRF_ON(MRF_strength);}
    if(PropUpdate>0.0f){STAPLE.Turn_Prop_Update_ON();}
    if(tmpP>0 && tmpQ>0 && tmpP<1 && tmpQ<1){STAPLE.SetPQ(tmpP,tmpQ);}

    STAPLE.SetMaximalIterationNumber(maxIteration);

    STAPLE.Run_EM();

    if(verbose_level>0){
        cout << "Saving Segmentation"<<endl;
    }
    nifti_image * Result = STAPLE.GetResult();
    nifti_image_write(Result);
    nifti_image_free(Result);
    nifti_image_free(LABLE);

    return 0;
}
