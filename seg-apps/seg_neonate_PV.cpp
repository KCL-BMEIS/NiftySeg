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
    printf("\t-in <filename>\t\t\tFilename of the input image image\n\n");

    printf("\t* * Options * *\n");
    printf("\t-out <filename>\t\t\tFilename of the segmented image (default=Segmentation.nii.gz)\n");
    printf("\t-mask <filename>\t\tFilename of the brainmask of the input image\n");
    printf("\t-EMseg <filename> \tSegmentation output from the EM algorithm\n");
    printf("\t-max_iter <int>\t\t\tMaximum number of iterations (default = 100)\n");
    printf("\t-v <int>\t\t\tVerbose level [0 = off, 1 = on, 2 = debug] (default = 0)\n");
    printf("\t-mrf_beta <float>\t\tMRF prior strength [off = 0, max = 1] (default = 0.25) \n");
    printf("\t-bc_order <int>\t\t\tPolinomial order for the bias field [off = 0, max = 5] (default = 5) \n");
    printf("\t-bc_thresh <float>\t\t\tBias field correction will run only if the ratio of improvment is below bc_thresh (default=0 [OFF]) \n");
    printf("\t-rf <relax_factor> <gaussian>\tRelax Priors [0<relax_factor<1 (default = 0.5), gaussian>0 (default=2.0)] \n");
    printf("\t-MAP <M V M V ...> \t\tMAP formulation: M and V are the a priori mean and variance for each class (default=off)\n");
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    return;
}

nifti_image * Transform_Priors_EM(nifti_image * PriorsEM,nifti_image * Mask, float thresh)
{
    nifti_image * Priors=nifti_copy_nim_info(PriorsEM);
    Priors->dim[0]=4;
    Priors->dim[4]=7;
    Priors->datatype=DT_FLOAT32;
    Priors->cal_max=1;

    nifti_update_dims_from_array(Priors);
    nifti_datatype_sizes(Priors->datatype,&Priors->nbyper,&Priors->swapsize);
    Priors->data = (void *) calloc(Priors->nvox, sizeof(PrecisionTYPE));
    char filename[] = "input_priors.nii";
    nifti_set_filenames(Priors,filename,0,0);
    int img_size= Priors->nx * Priors->ny * Priors->nz;
    PrecisionTYPE * PriorsEM_ptr= static_cast<PrecisionTYPE *>(PriorsEM->data);
    bool * Mask_ptr= static_cast<bool *>(Mask->data);

    bool * GM_to_smooth= new bool [img_size]();
    bool * CSF_to_smooth= new bool [img_size]();

    for(int i=0; i<img_size;i++){
        GM_to_smooth[i]=(PriorsEM_ptr[i+img_size]+PriorsEM_ptr[i+img_size*3]+PriorsEM_ptr[i+img_size*4])>0.3;
        CSF_to_smooth[i]=((PriorsEM_ptr[i+img_size*2]>0.1)+(!Mask_ptr[i]));
    }

    ImageSize * CurrSizes = new ImageSize [1]();
    CurrSizes->numel=(int)(PriorsEM->nx*PriorsEM->ny*PriorsEM->nz);
    CurrSizes->xsize=PriorsEM->nx;
    CurrSizes->ysize=PriorsEM->ny;
    CurrSizes->zsize=PriorsEM->nz;
    CurrSizes->usize=1;
    CurrSizes->tsize=1;
    CurrSizes->numclass=1;
    CurrSizes->numelmasked=0;
    CurrSizes->numelbias=0;



    int dimensions[3];
    dimensions[0]=CurrSizes->xsize;
    dimensions[1]=CurrSizes->ysize;
    dimensions[2]=CurrSizes->zsize;


    Dillate(GM_to_smooth,1,dimensions);
    Dillate(CSF_to_smooth,1,dimensions);



    float * PVprior= new float [img_size]();
    for(int i=0; i<img_size;i++){
        if(Mask_ptr[i]>0){
            PVprior[i]=((float)((CSF_to_smooth[i]&&GM_to_smooth[i])))*PriorsEM_ptr[i];
        }
    }




    PrecisionTYPE * Priors_ptr = static_cast<PrecisionTYPE *>(Priors->data);
    for(long cl=0;cl<6;cl++){
        PrecisionTYPE * Priors_ptr_class = &Priors_ptr[cl*img_size];
        PrecisionTYPE * PriorsEM_ptr_class = &PriorsEM_ptr[cl*img_size];
        for(long i=0; i<img_size;i++,Priors_ptr_class++,PriorsEM_ptr_class++){
            if(Mask_ptr[i]>0){
                if(cl==0){
                    (*Priors_ptr_class)=(*PriorsEM_ptr_class)-PVprior[i];
                }
                else if(cl==2){
                    (*Priors_ptr_class)=(*PriorsEM_ptr_class);
                }
                else{
                    (*Priors_ptr_class)=(*PriorsEM_ptr_class);
                }
            }
        }
    }

    PrecisionTYPE * Priors_ptr_class = &Priors_ptr[6*img_size];
    for(long i=0; i<img_size;i++){
        if(Mask_ptr[i]){
            (*Priors_ptr_class)=PVprior[i];
        }
        Priors_ptr_class++;
    }

    float * Tmpimage=NULL;
    for(long cl=0;cl<6;cl++){
        PrecisionTYPE * Priors_ptr_class = &Priors_ptr[cl*img_size];
        Tmpimage=Gaussian_Filter_4D_inside_mask(Priors_ptr_class,Mask_ptr,1.0,CurrSizes);
        for(int i=0; i<img_size; i++){
            if(Mask_ptr[i]){
                Priors_ptr_class[i]=Tmpimage[i];
            }
        }
        delete [] Tmpimage;
    }

    int ups=0;
    for(int i=0; i<img_size; i++){
        float sumexp=0;
        if(Mask_ptr[i]>0){
            for(int currclass=0; currclass<Priors->dim[4];currclass++){
                if((Priors_ptr[i+currclass*img_size])>0.01){
                    sumexp+=Priors_ptr[i+currclass*img_size];
                }
            }
            if(sumexp>0){
                for(int currclass=0; currclass<Priors->dim[4];currclass++){
                    if((Priors_ptr[i+currclass*img_size])>thresh){
                        Priors_ptr[i+currclass*img_size]=Priors_ptr[i+currclass*img_size]/sumexp;
                    }
                    else{
                        Priors_ptr[i+currclass*img_size]=0;
                    }
                }
            }
            else{
                ups++;
                for(int currclass=0;currclass<Priors->dim[4];currclass++){
                    Priors_ptr[i+currclass*img_size]=1/7;
                }
            }
        }
    }

    delete [] GM_to_smooth;
    delete [] CSF_to_smooth;
    delete [] PVprior;
    //nifti_image_write(Priors);


    return Priors;
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
    segment_param->numb_classes=7;
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
            segment_param->filename_priors= (char **) calloc(1,sizeof(char *));
            segment_param->filename_priors[0]=argv[++i];
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
        else if(strcmp(argv[i], "-v") == 0){
            segment_param->verbose_level=(int)atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-bc_out") == 0){
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
    seg_changeDatatype<PrecisionTYPE>(T1);

    // READING MASK - mandatory for now
    nifti_image * Mask=NULL;
    if(segment_param->flag_mask){
        Mask = nifti_image_read(segment_param->filename_mask,true);
        seg_changeDatatype<PrecisionTYPE>(Mask);
        if(Mask->datatype!=DT_BINARY){
            seg_convert2binary(Mask,0.0f);
        }
        if(Mask == NULL){
            fprintf(stderr,"* Error when reading the mask image: %s\n",segment_param->filename_mask);
            return 1;
        }
    }

    cout << "\n\n***********************\n***********************\n      Stage 1 \n***********************\n***********************\n\n";
    flush(cout);
    // READING PRIORS - mandatory for now
    nifti_image * Priors_EM=NULL;
    nifti_image * Priors=NULL;
    if(segment_param->flag_manual_priors){
        Priors_EM = nifti_image_read(segment_param->filename_priors[0],true);
        Priors=Transform_Priors_EM(Priors_EM,Mask,0.005);
    }
    else{
        fprintf(stderr,"* Priors option have to be ON for now...\n");
        return 1;
    }
    nifti_image_free(Priors_EM);

    seg_EM<float>* SEG = new seg_EM<float>(segment_param->numb_classes,1,1);
    //seg_EM<float> SEG(segment_param->numb_classes,T1->dim[4],1);
    SEG->SetInputImage(T1);
    if(segment_param->flag_mask) {
        SEG->SetMaskImage(Mask);
    }
    SEG->SetPriorImage(Priors);
    SEG->SetVerbose(segment_param->verbose_level);
    SEG->SetFilenameOut(segment_param->filename_out);
    if(segment_param->flag_Bias)
        SEG->Turn_BiasField_ON(segment_param->bias_order,segment_param->Bias_threshold);
    if(segment_param->flag_MRF)
        SEG->Turn_MRF_ON(segment_param->MRF_strength);
    if(segment_param->relax_factor>0)
        SEG->Turn_Relaxation_ON(segment_param->relax_factor,segment_param->relax_gauss_kernel);
    if(segment_param->flag_MAP)
        SEG->SetMAP(segment_param->MAP_M,segment_param->MAP_V);
    SEG->SetAprox(segment_param->aprox);
    SEG->SetMaximalIterationNumber(segment_param->maxIteration);


    SEG->Run_EM();
    nifti_image * Result1 = SEG->GetResultNeonate();
    delete SEG;
    nifti_image_free(Priors);
    //nifti_image_write(Result1);

    cout << "\n\n***********************\n***********************\n      Stage 2 \n***********************\n***********************\n\n";
flush(cout);
    nifti_image * Priors2=NULL;
    Priors2=Transform_Priors_EM(Result1,Mask,0.005);
    nifti_image_free(Result1);

    seg_EM<float>* SEG2 = new seg_EM<float>(segment_param->numb_classes,T1->dim[4],1);
    SEG2->SetInputImage(T1);
    if(segment_param->flag_mask) {
        SEG2->SetMaskImage(Mask);
    }
    SEG2->SetPriorImage(Priors2);
    SEG2->SetVerbose(segment_param->verbose_level);
    SEG2->SetFilenameOut(segment_param->filename_out);
    if(segment_param->flag_Bias)
        SEG2->Turn_BiasField_ON(segment_param->bias_order,segment_param->Bias_threshold);
    if(segment_param->flag_MRF)
        SEG2->Turn_MRF_ON(segment_param->MRF_strength);
    if(segment_param->relax_factor>0)
        SEG2->Turn_Relaxation_ON(segment_param->relax_factor,segment_param->relax_gauss_kernel);
    if(segment_param->flag_MAP)
        SEG2->SetMAP(segment_param->MAP_M,segment_param->MAP_V);
    SEG2->SetAprox(segment_param->aprox);
    SEG2->SetMaximalIterationNumber(segment_param->maxIteration);
    SEG2->Run_EM();

    nifti_image * Result2 = SEG2->GetResultNeonate();
    delete SEG2;
    nifti_image_free(Priors2);
    //nifti_image_write(Result2);
    cout << "\n\n***********************\n***********************\n      Stage 3 \n***********************\n***********************\n\n";
  flush(cout);
    nifti_image * Priors3=NULL;
    if(segment_param->flag_manual_priors){
        Priors3=Transform_Priors_EM(Result2,Mask,0.01);
    }
    else{
        fprintf(stderr,"* Priors option have to be ON for now...\n");
        return 1;
    }
    nifti_image_free(Result2);

    seg_EM<float>* SEG3 = new seg_EM<float>(segment_param->numb_classes,T1->dim[4],1);
    SEG3->SetInputImage(T1);
    if(segment_param->flag_mask) {
        SEG3->SetMaskImage(Mask);
    }
    SEG3->SetPriorImage(Priors3);
    SEG3->SetVerbose(segment_param->verbose_level);
    SEG3->SetFilenameOut(segment_param->filename_out);
    if(segment_param->flag_Bias)
        SEG3->Turn_BiasField_ON(segment_param->bias_order,segment_param->Bias_threshold);
    if(segment_param->flag_MRF)
        SEG3->Turn_MRF_ON(segment_param->MRF_strength);
    if(segment_param->relax_factor>0)
        SEG3->Turn_Relaxation_ON(segment_param->relax_factor,segment_param->relax_gauss_kernel);
    if(segment_param->flag_MAP)
        SEG3->SetMAP(segment_param->MAP_M,segment_param->MAP_V);
    SEG3->SetAprox(segment_param->aprox);
    SEG3->SetMaximalIterationNumber(segment_param->maxIteration);
    SEG3->Run_EM();

    nifti_image * Result3 = SEG3->GetResultNeonate();
    delete SEG3;
    nifti_image_free(Priors3);
    nifti_image_write(Result3);
    nifti_image_free(Result3);
/*
    cout << "\n\n***********************\n***********************\n      Stage 4 \n***********************\n***********************\n\n";
    flush(cout);
    nifti_image * Priors4=NULL;
    if(segment_param->flag_manual_priors){
        Priors4=Transform_Priors_EM(Result3,Mask,0.02);
    }
    else{
        fprintf(stderr,"* Priors option have to be ON for now...\n");
        return 1;
    }
    nifti_image_free(Result3);

    seg_EM<float>* SEG4 = new seg_EM<float>(segment_param->numb_classes,T1->dim[4],1);
    SEG4->SetInputImage(T1);
    if(segment_param->flag_mask) {
        SEG4->SetMaskImage(Mask);
    }
    SEG4->SetPriorImage(Priors4);
    SEG4->SetVerbose(segment_param->verbose_level);
    SEG4->SetFilenameOut(segment_param->filename_out);
    if(segment_param->flag_Bias)
        SEG4->Turn_BiasField_ON(segment_param->bias_order,segment_param->Bias_threshold);
    if(segment_param->flag_MRF)
        SEG4->Turn_MRF_ON(segment_param->MRF_strength);
    if(segment_param->relax_factor>0)
        SEG4->Turn_Relaxation_ON(segment_param->relax_factor,segment_param->relax_gauss_kernel);
    if(segment_param->flag_MAP)
        SEG4->SetMAP(segment_param->MAP_M,segment_param->MAP_V);
    SEG4->SetAprox(segment_param->aprox);
    SEG4->SetMaximalIterationNumber(segment_param->maxIteration);
    SEG4->Run_EM();

    nifti_image * Result4 = SEG4->GetResultNeonate();
    delete SEG4;
    nifti_image_free(Priors4);
    nifti_image_write(Result4);
    nifti_image_free(Result4);
    */
    nifti_image_free(T1);
    nifti_image_free(Mask);

    free(segment_param->filename_priors);

    delete [] segment_param;
    return 0;
}

