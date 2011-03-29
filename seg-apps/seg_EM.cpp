#include "_seg_EM.h"
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
    printf("\t-in <filename>\t\t| Filename of the input image image\n");
    printf("\t\t\t\t| The input image should be 2D, 3D or 4D images. 2D images should be on the XY plane.\n");
    printf("\t\t\t\t| 4D images are segmented as if they were multimodal.\n");
    printf("\t-priors <n> <filenames>\t| The number of priors (n>0) and their filenames\n");
    printf("\t\t\t\t| Priors are mandatory for now. Won't be in the next release.\n\n");
    printf("\t* * Options * *\n");
    printf("\t-out <filename>\t\t| Filename of the segmented image ( default = Segmentation.nii.gz)\n");
    printf("\t-mask <filename>\t| Filename of the brainmask of the input image\n");
    printf("\t-max_iter <int>\t\t| Maximum number of iterations (default = 100)\n");
    printf("\t-v <int>\t\t| Verbose level [0 = off, 1 = on, 2 = debug] (default = 0)\n");
    printf("\t-mrf_beta <float>\t| MRF prior strength [off = 0, max = 1] (default = 0.25) \n");
    printf("\t-bc_order <int>\t\t| Polinomial order for the bias field [off = 0, max = 5] (default = 5) \n");
    printf("\t-bc_thresh <float>\t| Bias field correction will run only if the ratio of improvment is below bc_thresh (default=0 [OFF]) \n");
    printf("\t-bc_out <filename>\t| Output the bias corrected image\n");
    //printf("\t-aprox_off \t\t| All aproximations off\n");
    printf("\t-rf <rel> <gstd>\t| Relax Priors [relaxation factor: 0<rel<1 (default = 0.5), gaussian regularization: gstd>0 (default=2.0)] /only 3D/\n");
    printf("\t-MAP <M V M V ...> \t| MAP formulation: M and V are the a priori mean and variance for each class (default=off) /only 3D/\n");
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    return;
}

void Merge_Priors(nifti_image * Priors, nifti_image ** Priors_temp, SEG_PARAM * segment_param)
{
    long img_size= Priors->nx * Priors->ny * Priors->nz;
    PrecisionTYPE * Prior_ptr_start = static_cast<PrecisionTYPE *>(Priors->data);
    for(long cl=0;cl<(segment_param->numb_classes);cl++){
        PrecisionTYPE * Prior_tmp_ptr = static_cast<PrecisionTYPE *>(Priors_temp[cl]->data);
        PrecisionTYPE * Prior_ptr = &Prior_ptr_start[cl*img_size];
        for(long i=0; i<img_size;i++){
            (*Prior_ptr)=(*Prior_tmp_ptr);
            Prior_ptr++;
            Prior_tmp_ptr++;
        }
    }
    return;
}

int main(int argc, char **argv)
{




    if (argc < 2)
    {
        Usage(argv[0]);
        return 0;
    }


    SEG_PARAM * segment_param = new SEG_PARAM [1]();
    segment_param->maxIteration=100;
    segment_param->flag_T1=0;
    segment_param->flag_out=0;
    segment_param->flag_mask=0;
    segment_param->flag_MRF=1;
    segment_param->flag_Bias=1;
    segment_param->flag_SG_deli=1;
    segment_param->flag_bc_out=0;
    segment_param->relax_factor=0;
    segment_param->relax_gauss_kernel=0;
    segment_param->flag_PV_model=1;
    segment_param->verbose_level=0;
    segment_param->flag_manual_priors=0;
    segment_param->bias_order=5;
    segment_param->MRF_strength=0.25f;
    segment_param->Bias_threshold=0;
    segment_param->numb_classes=0;
    segment_param->aprox=false;
    /* read the input parameter */
    int To_do_MAP_index_argv=0;
    for(int i=1;i<argc;i++){
        if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
           strcmp(argv[i], "-HELP")==0 || strcmp(argv[i], "-h")==0 ||
           strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0){
            Usage(argv[0]);
            return 0;
        }

        else if(strcmp(argv[i], "-in") == 0){
            segment_param->filename_T1 = argv[++i];
            segment_param->flag_T1=1;
        }
        else if(strcmp(argv[i], "-out") == 0){
            segment_param->filename_out = argv[++i];
            segment_param->flag_out=1;
        }
        else if(strcmp(argv[i], "-mask") == 0){
            segment_param->filename_mask=argv[++i];
            segment_param->flag_mask=1;
        }
        else if(strcmp(argv[i], "-priors") == 0){
            segment_param->numb_classes=atoi(argv[++i]);
            if(segment_param->numb_classes<1){
                cout<<"Number of classes has to be bigger than 0";
                return 0;
            }
            segment_param->filename_priors= (char **) calloc(segment_param->numb_classes,sizeof(char *));
            for(int classnum=0; classnum<segment_param->numb_classes; classnum++){
                segment_param->filename_priors[classnum]=argv[++i];
            }
            segment_param->flag_manual_priors=1;
        }
        else if(strcmp(argv[i], "-MAP") == 0){
            if(segment_param->numb_classes>0){
                segment_param->MAP_M= (float*) calloc(segment_param->numb_classes,sizeof(float));
                segment_param->MAP_V= (float*) calloc(segment_param->numb_classes,sizeof(float));
                for(int classnum=0; classnum<segment_param->numb_classes; classnum++){
                    segment_param->MAP_M[classnum]=atof(argv[++i]);
                    segment_param->MAP_V[classnum]=atof(argv[++i]);
                }
                segment_param->flag_MAP=true;
            }
            else{
                To_do_MAP_index_argv=i;
            }
        }
        else if(strcmp(argv[i], "-bc_order") == 0){
            segment_param->bias_order=(int)(atof(argv[++i]));
            if(segment_param->bias_order==0){
                segment_param->flag_Bias=0;
            }
        }
        else if(strcmp(argv[i], "-bc_thresh") == 0){
            segment_param->Bias_threshold=atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-rf") == 0){
            segment_param->relax_factor=atof(argv[++i]);
            segment_param->relax_gauss_kernel=atof(argv[++i]);
        }
        //else if(strcmp(argv[i], "-aprox_off") == 0){
        //    segment_param->aprox=false;
        //}

        else if(strcmp(argv[i], "-v") == 0){
            segment_param->verbose_level=(int)atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-bc_out") == 0){
            segment_param->filename_bias = argv[++i];
            segment_param->flag_bc_out=1;
        }
        else if(strcmp(argv[i], "-mrf_beta") == 0){
            segment_param->MRF_strength=(PrecisionTYPE)atof(argv[++i]);
            if(segment_param->MRF_strength==0){
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

    string envVarName="NIFTYSEG";
    char * currpath_char=NULL;
    string currpath;
    currpath_char=getenv(envVarName.c_str());

    if(To_do_MAP_index_argv>0){
        int i=To_do_MAP_index_argv;
        if(segment_param->numb_classes>0){
            segment_param->MAP_M= (float*) calloc(segment_param->numb_classes,sizeof(float));
            segment_param->MAP_V= (float*) calloc(segment_param->numb_classes,sizeof(float));
            for(int classnum=0; classnum<segment_param->numb_classes; classnum++){
                segment_param->MAP_M[classnum]=atof(argv[++i]);
                segment_param->MAP_V[classnum]=atof(argv[++i]);
            }
            segment_param->flag_MAP=true;
        }
    }
    if(segment_param->flag_out==0){
        fprintf(stderr,"Err:\tThe output image name has to be defined.\n");
        return 1;
    }

    if(!segment_param->flag_T1){
        fprintf(stderr,"Err:\tThe T1 image name has to be defined.\n");
        Usage(argv[0]);
        return 1;
    }

    // READING T1
    nifti_image * T1=nifti_image_read(segment_param->filename_T1,true);
    if(T1 == NULL){
        fprintf(stderr,"* Error when reading the T1 image: %s\n",segment_param->filename_T1);
        return 1;
    }
    reg_changeDatatype<PrecisionTYPE>(T1);

    // READING MASK - mandatory for now
    nifti_image * Mask=NULL;
    if(segment_param->flag_mask){
        Mask = nifti_image_read(segment_param->filename_mask,true);
        reg_changeDatatype<PrecisionTYPE>(Mask);

        if(Mask == NULL){
            fprintf(stderr,"* Error when reading the mask image: %s\n",segment_param->filename_mask);
            return 1;
        }
    }

    // READING PRIORS - mandatory for now
    nifti_image ** Priors_temp=new nifti_image * [segment_param->numb_classes];

    if(segment_param->flag_manual_priors){

        for(int i=0; i<segment_param->numb_classes; i++){
            Priors_temp[i] = nifti_image_read(segment_param->filename_priors[i],true);

            if(Priors_temp[i] == NULL){
                fprintf(stderr,"* Error when reading the WM Prior image: %s\n",segment_param->filename_priors[i]);
                return 1;
            }
            reg_changeDatatype<PrecisionTYPE>(Priors_temp[i]);
        }
    }
    else{
        fprintf(stderr,"* Priors option have to be ON for now...\n");
        return 1;
    }

    nifti_image * Priors=nifti_copy_nim_info(T1);
    Priors->dim[0]=4;
    Priors->dim[4]=segment_param->numb_classes;
    Priors->datatype=DT_FLOAT32;
    Priors->cal_max=1;

    nifti_update_dims_from_array(Priors);
    nifti_datatype_sizes(Priors->datatype,&Priors->nbyper,&Priors->swapsize);
    Priors->data = (void *) calloc(Priors->nvox, sizeof(PrecisionTYPE));

    Merge_Priors(Priors,Priors_temp,segment_param);
    for(int i=0;i<segment_param->numb_classes;i++){
        nifti_image_free(Priors_temp[i]);
        Priors_temp[i]=NULL;
    }
    delete [] Priors_temp;


    seg_EM<float> SEG(segment_param->numb_classes,T1->dim[4],1);
    SEG.SetInputImage(T1);
    if(segment_param->flag_mask) {
        SEG.SetMaskImage(Mask);
    }
    SEG.SetPriorImage(Priors);
    SEG.SetVerbose(segment_param->verbose_level);
    SEG.SetFilenameOut(segment_param->filename_out);
    if(segment_param->flag_Bias)
        SEG.Turn_BiasField_ON(segment_param->bias_order,segment_param->Bias_threshold);
    if(segment_param->flag_MRF)
        SEG.Turn_MRF_ON(segment_param->MRF_strength);
    if(segment_param->relax_factor>0)
        SEG.Turn_Relaxation_ON(segment_param->relax_factor,segment_param->relax_gauss_kernel);
    if(segment_param->flag_MAP)
        SEG.SetMAP(segment_param->MAP_M,segment_param->MAP_V);
    SEG.SetAprox(segment_param->aprox);
    SEG.SetMaximalIterationNumber(segment_param->maxIteration);


    SEG.Run_EM();

    if(segment_param->verbose_level>0){
        cout << "Saving Segmentation"<<endl;
    }
    nifti_image * Result = SEG.GetResult();
    nifti_image_write(Result);
    nifti_image_free(Result);
    nifti_image_free(Priors);

    nifti_image * BiasFieldCorrected=NULL;

    if(segment_param->flag_bc_out){
        if(segment_param->verbose_level>0){
            cout << "Saving Bias Field Corrected Image"<<endl;
        }
        BiasFieldCorrected = SEG.GetBiasCorrected(segment_param->filename_bias);
        nifti_image_write(BiasFieldCorrected);
        nifti_image_free(BiasFieldCorrected);
    }



    nifti_image_free(T1);
    nifti_image_free(Mask);

    free(segment_param->filename_priors);

    delete [] segment_param;
    return 0;
}

