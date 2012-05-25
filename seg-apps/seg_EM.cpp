#include "_seg_EM.h"
#include <iostream>
#include <time.h>
#include <stdlib.h>

using namespace std;
#define SegPrecisionTYPE float


void Usage(char *exec)
{
  printf("\n  EM Statistical Segmentation:\n  Usage ->\t%s -in <filename> [OPTIONS]\n\n",exec);
  printf("  * * * * * * * * * * * * * * * * * Mandatory * * * * * * * * * * * * * * * * * *\n\n");
  printf("\t-in <filename>\t\t| Filename of the input image\n");
  printf("\t-out <filename>\t\t| Filename of the segmented image\n");
  printf("\t\t\t\t| The input image should be 2D, 3D or 4D images. 2D images should be on the XY plane.\n");
  printf("\t\t\t\t| 4D images are segmented as if they were multimodal.\n\n");
  printf("  \t\t- Select one of the following (mutually exclusive) -\n\n");
  printf("\t-priors <n> <fnames>\t| The number of priors (n>0) and their filenames. Priors should be registerd to the input image\n");
  printf("\t-priors4D <fname>\t| 4D image with the piors stacked in the 4th dimension. Priors should be registerd to the input image\n");
  printf("\t-nopriors <n>\t\t| The number of classes (n>0) \n\n");
  printf("  * * * * * * * * * * * * * * * * General Options * * * * * * * * * * * * * * * * *\n\n");
  printf("\t-mask <filename>\t| Filename of the brain-mask of the input image\n");
  printf("\t-max_iter <int>\t\t| Maximum number of iterations (default = 100)\n");
  printf("\t-v <int>\t\t| Verbose level [0 = off, 1 = on, 2 = debug] (default = 0)\n");
  printf("\t-mrf_beta <float>\t| MRF prior strength [off = 0, max = 1] (default = 0.4) \n");
  printf("\t-bc_order <int>\t\t| Polynomial order for the bias field [off = 0, max = 5] (default = 3) \n");
  printf("\t-bc_thresh <float>\t| Bias field correction will run only if the ratio of improvement is below bc_thresh (default=0 [OFF]) \n");
  printf("\t-bc_out <filename>\t| Output the bias corrected image\n");
  printf("\t-reg <filename>\t\t| Amount of regularization over the diagonal of the covariance matrix [above 1]\n");
  printf("\t-outlier <fl1> <fl2>\t| Outlier detection as in (Van Leemput TMI 2003). <fl1> is the Mahalanobis threshold [recommended between 3 and 7] \n");
  printf("\t\t\t\t| <fl2> is a convergence ratio below wich the outlier detection is going to be done [recommended 0.001].\n");
  printf("\t-out_outlier <filename>\t| Output outlierness image \n");
  printf("\t-rf <rel> <gstd>\t| Relax Priors [relaxation factor: 0<rf<1 (recommended=0.5), gaussian regularization: gstd>0 (recommended=2.0)] /only 3D/\n");
  printf("\t-MAP <M V M V ...> \t| MAP formulation: M and V are the parameters (mean & variance) of the semiconjugate prior over the class mean\n\n");
  printf(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
  return;
}

void Merge_Priors(nifti_image * Priors, nifti_image ** Priors_temp, SEG_PARAM * segment_param)
{
  long img_size= Priors->nx * Priors->ny * Priors->nz;
  SegPrecisionTYPE * Prior_ptr_start = static_cast<SegPrecisionTYPE *>(Priors->data);
  for(long cl=0;cl<(segment_param->numb_classes);cl++){
      SegPrecisionTYPE * Prior_tmp_ptr = static_cast<SegPrecisionTYPE *>(Priors_temp[cl]->data);
      SegPrecisionTYPE * Prior_ptr = &Prior_ptr_start[cl*img_size];
      for(long i=0; i<img_size;i++){
          (*Prior_ptr)=(*Prior_tmp_ptr);
          Prior_ptr++;
          Prior_tmp_ptr++;
        }
    }
  return;
}

void no_memory () {
  cout << "Failed to allocate memory!\n";
  exit (1);
}

int main(int argc, char **argv)
{
  try
  {

    set_new_handler(no_memory);

    if (argc <=1)
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
    segment_param->bias_order=3;
    segment_param->MRF_strength=0.4f;
    segment_param->Bias_threshold=0;
    segment_param->numb_classes=1;
    segment_param->aprox=false;
    segment_param->flag_Outlierness=0;
    segment_param->OutliernessThreshold=0;
    segment_param->flag_out_outlier=0;
    float OutliernessRatio=0.01;
    int To_do_MAP_index_argv=0;
    float regularization_amount=1;
    bool priors4D=false;



    /* read the input parameter */
    for(int i=1;i<argc;i++){
        if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
           strcmp(argv[i], "-HELP")==0 || strcmp(argv[i], "-h")==0 ||
           strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0){
            Usage(argv[0]);
            return 0;
          }

        else if(strcmp(argv[i], "-in") == 0 && (i+1)<argc){
            segment_param->filename_T1 = argv[++i];
            segment_param->flag_T1=1;
          }
        else if(strcmp(argv[i], "-out") == 0 && (i+1)<argc){
            segment_param->filename_out = argv[++i];
            segment_param->flag_out=1;
          }
        else if(strcmp(argv[i], "-mask") == 0 && (i+1)<argc){
            segment_param->filename_mask=argv[++i];
            segment_param->flag_mask=1;
          }
        else if(strcmp(argv[i], "-priors") == 0 && (i+1)<argc){
            segment_param->numb_classes=atoi(argv[++i]);
            if(segment_param->numb_classes<2){
                cout<<"Number of classes has to be bigger than 1";
                return 0;
              }
            if((i+segment_param->numb_classes)<argc){
                segment_param->filename_priors= (char **) calloc(segment_param->numb_classes,sizeof(char *));
                for(int classnum=0; classnum<segment_param->numb_classes; classnum++){
                    segment_param->filename_priors[classnum]=argv[++i];
                  }
              }
            else{
                fprintf(stderr,"Err:\tParameter -priors are incomplete\n");
                Usage(argv[0]);
                return 1;
              }
            segment_param->flag_manual_priors=1;
          }
        else if(strcmp(argv[i], "-priors4D") == 0 && (i+1)<argc){
            segment_param->filename_priors= (char **) calloc(1,sizeof(char *));
            segment_param->filename_priors[0]=argv[++i];
            segment_param->flag_manual_priors=1;
            nifti_image * tmpread = nifti_image_read(segment_param->filename_priors[0],false);
            segment_param->numb_classes=tmpread->nt;
            if(segment_param->numb_classes<2){
                cout<<"Number of classes has to be bigger than 1";
                return 0;
              }
            nifti_image_free(tmpread);
            priors4D=true;
          }
        else if(strcmp(argv[i], "-nopriors") == 0 && (i+1)<argc){
            segment_param->flag_manual_priors=0;
            segment_param->numb_classes=atoi(argv[++i]);
          }
        else if(strcmp(argv[i], "-outlier") == 0 && (i+2)<argc){
            segment_param->flag_Outlierness=1;
            segment_param->OutliernessThreshold=atof(argv[++i]);
            OutliernessRatio=atof(argv[++i]);
          }
        else if(strcmp(argv[i], "-out_outlier") == 0 && (i+1)<argc){
            segment_param->filename_out_outlier = argv[++i];
            segment_param->flag_out_outlier=1;
          }
        else if(strcmp(argv[i], "-MAP") == 0){
            if(segment_param->numb_classes>0 && (i+segment_param->numb_classes)<argc){
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
        else if(strcmp(argv[i], "-bc_order") == 0 && (i+1)<argc){
            segment_param->bias_order=(int)(atof(argv[++i]));
            if(segment_param->bias_order==0){
                segment_param->flag_Bias=0;
              }
          }
        else if(strcmp(argv[i], "-bc_thresh") == 0 && (i+1)<argc){
            segment_param->Bias_threshold=atof(argv[++i]);
          }
        else if(strcmp(argv[i], "-reg") == 0 && (i+1)<argc){
            regularization_amount=atof(argv[++i]);
          }
        else if(strcmp(argv[i], "-rf") == 0 && (i+1)<argc){
            segment_param->relax_factor=atof(argv[++i]);
            segment_param->relax_gauss_kernel=atof(argv[++i]);
          }
        //else if(strcmp(argv[i], "-aprox_off") == 0){
        //    segment_param->aprox=false;
        //}

        else if(strcmp(argv[i], "-v") == 0 &&(i+1)<argc){
            segment_param->verbose_level=(int)atoi(argv[++i]);
          }
        else if(strcmp(argv[i], "-bc_out") == 0 && (i+1)<argc){
            segment_param->filename_bias = argv[++i];
            segment_param->flag_bc_out=1;
          }
        else if(strcmp(argv[i], "-mrf_beta") == 0 && (i+1)<argc){
            segment_param->MRF_strength=(SegPrecisionTYPE)atof(argv[++i]);

          }
        else if(strcmp(argv[i], "-max_iter") == 0 && (i+1)<argc){
            segment_param->maxIteration=atoi(argv[++i]);
          }
        else{
            fprintf(stderr,"Err:\tParameter %s unknown or incomplete \n",argv[i]);
            Usage(argv[0]);
            return 1;
          }
      }
    if(segment_param->MRF_strength<=0){
        segment_param->flag_MRF=0;
      }

    //string envVarName="NIFTYSEG";
    //char * currpath_char=NULL;
    //currpath_char=getenv(envVarName.c_str());


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

    if(!segment_param->flag_T1){
        fprintf(stderr,"Err:\tThe T1 image name has to be defined.\n");
        Usage(argv[0]);
        return 1;
      }

    if(segment_param->flag_out==0){
        fprintf(stderr,"Err:\tThe output image name has to be defined.\n");
        Usage(argv[0]);
        return 1;
      }



    // READING T1
    nifti_image * InputImage=nifti_image_read(segment_param->filename_T1,true);
    if(InputImage == NULL){
        fprintf(stderr,"* Error when reading the T1 image: %s\n",segment_param->filename_T1);
        return 1;
      }
    if(InputImage->datatype!=NIFTI_TYPE_FLOAT32)
      seg_changeDatatype<SegPrecisionTYPE>(InputImage);

    InputImage->dim[4]=InputImage->nt=(InputImage->nt<1)?1:InputImage->nt;
    InputImage->dim[5]=InputImage->nu=(InputImage->nu<1)?1:InputImage->nu;
    if(InputImage->nu>1){
        InputImage->dim[5]=InputImage->nu;
        InputImage->dim[4]=1;
      }
    else if(InputImage->nt>1){
        InputImage->dim[5]=InputImage->nt;
        InputImage->dim[4]=1;
      }
    else if(InputImage->nt>1){
        InputImage->dim[5]=1;
        InputImage->dim[4]=1;
      }
    InputImage->dim[0]=5;
    nifti_update_dims_from_array(InputImage);


    nifti_image * Mask=NULL;
    if(segment_param->flag_mask){
        Mask = nifti_image_read(segment_param->filename_mask,true);
        if(Mask->datatype!=NIFTI_TYPE_FLOAT32)
          seg_changeDatatype<SegPrecisionTYPE>(Mask);

        if(Mask == NULL){
            fprintf(stderr,"* Error when reading the mask image: %s\n",segment_param->filename_mask);
            return 1;
          }
        if( (Mask->nx != InputImage->nx) || (Mask->ny != InputImage->ny) || (Mask->nz != InputImage->nz) ){
            fprintf(stderr,"* Error: Mask image not the same size as input\n");
            return 1;
          }
        Mask->dim[4]=Mask->nt=(Mask->nt<1)?1:Mask->nt;
        Mask->dim[5]=Mask->nu=(Mask->nu<1)?1:Mask->nu;
        nifti_update_dims_from_array(Mask);
      }

    nifti_image ** Priors_temp=new nifti_image * [segment_param->numb_classes];
    nifti_image * Priors=NULL;

    if(segment_param->flag_manual_priors){

        if(priors4D){
            int i=0;
            Priors = nifti_image_read(segment_param->filename_priors[i],true);
            if(Priors == NULL){
                fprintf(stderr,"* Error when reading the Prior image: %s\n",segment_param->filename_priors[i]);
                return 1;
              }
            if( (Priors->nx != InputImage->nx) || (Priors->ny != InputImage->ny) || (Priors->nz != InputImage->nz) ){
                fprintf(stderr,"* Error: Prior image ( %s ) not the same size as input\n",segment_param->filename_priors[i]);
                return 1;
              }
            seg_changeDatatype<SegPrecisionTYPE>(Priors);

          }
        else{

            for(int i=0; i<segment_param->numb_classes; i++){
                Priors_temp[i] = nifti_image_read(segment_param->filename_priors[i],true);

                if(Priors_temp[i] == NULL){
                    fprintf(stderr,"* Error when reading the Prior image: %s\n",segment_param->filename_priors[i]);
                    return 1;
                  }
                if( (Priors_temp[i]->nx != InputImage->nx) || (Priors_temp[i]->ny != InputImage->ny) || (Priors_temp[i]->nz != InputImage->nz) ){
                    fprintf(stderr,"* Error: Prior image ( %s ) not the same size as input\n",segment_param->filename_priors[i]);
                    return 1;
                  }
                seg_changeDatatype<SegPrecisionTYPE>(Priors_temp[i]);
              }

            Priors=nifti_copy_nim_info(InputImage);
            Priors->dim[0]=4;
            Priors->nt=Priors->dim[4]=segment_param->numb_classes;
            Priors->datatype=NIFTI_TYPE_FLOAT32;
            Priors->cal_max=1;

            nifti_update_dims_from_array(Priors);
            nifti_datatype_sizes(Priors->datatype,&Priors->nbyper,&Priors->swapsize);
            Priors->data = (void *) calloc(Priors->nvox, sizeof(SegPrecisionTYPE));

            Merge_Priors(Priors,Priors_temp,segment_param);
            for(int i=0;i<segment_param->numb_classes;i++){
                nifti_image_free(Priors_temp[i]);
                Priors_temp[i]=NULL;
              }
          }
      }
    else{
        if(segment_param->numb_classes<2){
            fprintf(stderr,"Please use either the -priors or -nopriors options with at least 2 classes\n");
            return 1;
          }
        else{
            if(segment_param->verbose_level>0)
              printf("WARNING: No priors selected. Results might be inconsistent\n");
          }

      }

    delete [] Priors_temp;


    seg_EM SEG(segment_param->numb_classes,InputImage->dim[5],InputImage->dim[4]);
    SEG.SetInputImage(InputImage);
    if(segment_param->flag_mask)
      SEG.SetMaskImage(Mask);
    if(segment_param->flag_manual_priors)
      SEG.SetPriorImage(Priors);

    SEG.SetVerbose(segment_param->verbose_level);
    SEG.SetFilenameOut(segment_param->filename_out);
    SEG.SetAprox(segment_param->aprox);
    SEG.SetMaximalIterationNumber(segment_param->maxIteration);

    if(segment_param->flag_Outlierness)
      SEG.OutliernessON(segment_param->OutliernessThreshold,OutliernessRatio);
    if(segment_param->flag_Bias)
      SEG.Turn_BiasField_ON(segment_param->bias_order,segment_param->Bias_threshold);
    if(segment_param->flag_MRF)
      SEG.Turn_MRF_ON(segment_param->MRF_strength);
    if(segment_param->relax_factor>0)
      SEG.Turn_Relaxation_ON(segment_param->relax_factor,segment_param->relax_gauss_kernel);
    if(segment_param->flag_MAP)
      SEG.SetMAP(segment_param->MAP_M,segment_param->MAP_V);
    if(regularization_amount>0)
      SEG.SetRegValue(regularization_amount);

    SEG.Run_EM();

    if(segment_param->verbose_level>0){
        cout << "Saving Segmentation"<<endl;
      }
    nifti_image * Result = SEG.GetResult();
    nifti_image_write(Result);
    nifti_image_free(Result);
    nifti_image_free(Priors);


    nifti_image * OutliernessImage=NULL;
    if(segment_param->flag_out_outlier && segment_param->flag_Outlierness){
        if(segment_param->verbose_level>0){
            cout << "Saving Outlierness Image"<<endl;
          }
        OutliernessImage=SEG.GetOutlierness(segment_param->filename_out_outlier);
        nifti_image_write(OutliernessImage);
        nifti_image_free(OutliernessImage);
      }

    nifti_image * BiasFieldCorrected=NULL;

    if(segment_param->flag_bc_out){
        if(segment_param->verbose_level>0){
            cout << "Saving Bias Field Corrected Image"<<endl;
          }
        BiasFieldCorrected = SEG.GetBiasCorrected(segment_param->filename_bias);
        nifti_image_write(BiasFieldCorrected);
        nifti_image_free(BiasFieldCorrected);
      }



    nifti_image_free(InputImage);
    nifti_image_free(Mask);

    free(segment_param->filename_priors);

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

