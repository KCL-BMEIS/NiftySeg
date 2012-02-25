#ifndef _SEG_EM_CPP
#define _SEG_EM_CPP
#include "_seg_EM.h"

seg_EM::seg_EM(int _numb_classes, int _nu,int _nt)
{

  this->inputImage=NULL; // pointer to external

  this->FilenameOut="Segmentation.nii.gz";

  this->dimentions=1;
  this->nx=1;
  this->ny=0;
  this->nz=0;
  this->nu=_nu;
  this->nt=_nt;
  this->dx=0;
  this->dy=0;
  this->dz=0;
  this->numel=0;
  this->aprox=true;
  this->iter=0;
  this->checkpoint_iter=0;
  this->ratio=0;

  this->numb_classes=_numb_classes;
  this->M= new float [MaxMultispectalSize*max_numbclass];
  for(int i=0; i<(MaxMultispectalSize*max_numbclass); i++){
      this->M[i]=0.0f;
    }
  this->V= new float [MaxMultispectalSize*MaxMultispectalSize*max_numbclass];
  for(int i=0; i<(MaxMultispectalSize*MaxMultispectalSize*max_numbclass); i++){
      this->V[i]=0.0f;
    }

  this->Expec=NULL;
  this->ShortPrior=NULL;
  this->Short_2_Long_Indices=NULL;
  this->Long_2_Short_Indices=NULL;
  this->CurrSizes=NULL;
  this->reg_factor=1.1f;

  this->maxIteration=100;
  this->verbose_level=0;
  this->loglik=2.0;
  this->oldloglik=1.0;


  this->maskImage_status=false;
  this->Mask=NULL; // pointer to external
  this->numelmasked=0;

  this->Priors_status=false;
  this->Priors=NULL;

  this->MRF_status=false;
  this->MRF_strength=0.0f;
  this->MRF=NULL;
  this->MRF_beta=NULL;
  this->MRF_transitionMatrix=NULL;

  this->Outlierness=NULL;
  this->OutliernessUSE=NULL;
  this->OutliernessFlag=false;
  this->OutliernessThreshold=0.0f;
  this->Outlierness_ratio=0.01f;

  this->BiasField_status=false;
  this->BiasField_order=0;
  this->BiasField=NULL;
  this->BiasField_coeficients=NULL;
  this->BiasField_ratio=0;

  this->LoAd_phase=-1;
  this->PV_model_status=false;
  this->SG_deli_status=false;
  this->Relax_status=false;
  this->Relax_factor=0;
  this->RelaxGaussKernelSize=0;
  this->MAP_status=0;
  this->MAP_M=NULL;
  this->MAP_V=NULL;
}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

seg_EM::~seg_EM()
{
  if(this->Expec!=NULL){
      delete[] this->Expec;
    }
  this->Expec=NULL;

  if(this->ShortPrior!=NULL){
      delete [] this->ShortPrior;
    }
  this->ShortPrior=NULL;

  if(this->BiasField!=NULL){
      delete [] this->BiasField;
    }
  this->BiasField=NULL;

  if(this->BiasField_coeficients!=NULL){
      delete [] this->BiasField_coeficients;
    }
  this->BiasField_coeficients=NULL;

  if(this->Outlierness!=NULL){
      delete [] this->Outlierness;
    }
  if(this->MRF_status){
      delete [] this->MRF;
      delete [] this->MRF_transitionMatrix;
    }
  this->MRF=NULL;
  this->MRF_transitionMatrix=NULL;

  if(this->maskImage_status){
      delete [] this->Long_2_Short_Indices;
      delete [] this->Short_2_Long_Indices;

    }

  if(this->CurrSizes!=NULL)
    delete [] this->CurrSizes;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM::SetInputImage(nifti_image *r)
{
  this->inputImage = r;
  this->inputImage_status = true;
  // Size
  this->dimentions=(int)((r->nx)>1)+(int)((r->ny)>1)+(int)((r->nz)>1)+(int)((r->nt)>1)+(int)((r->nu)>1);
  this->nx=r->nx;
  this->ny=r->ny;
  this->nz=r->nz;
  this->nt=r->nt;
  this->nu=r->nu;
  this->dx=r->dx;
  this->dy=r->dy;
  this->dz=r->dz;
  this->numel=r->nz*r->ny*r->nx;
  if(this->nx==1 ||this->ny==1){
      cout<<"Error: The segmentation algorithm only takes 2D, 3D and 5D images. 2D images should be on the XY plane"<<endl;
      return 1;
    }

  return 0;
}


int seg_EM::SetMAP( float *M, float* V)
{
  this->MAP_status=true;
  this->MAP_M=new float[this->numb_classes];
  this->MAP_V=new float[this->numb_classes];


  for(int i=0;i<this->numb_classes;i++){
      this->MAP_M[i]=M[i];
      this->MAP_V[i]=V[i];
    }

  return 0;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM::SetPriorImage(nifti_image *r)
{
  this->Priors = r;
  this->Priors_status = true;
  // Size
  this->dimentions=(int)((r->nx)>1)+(int)((r->ny)>1)+(int)((r->nz)>1);
  if(this->nx==r->nx && this->ny==r->ny && this->nz==r->nz){
      return 0;
    }
  else{
      cout << "ERROR: Priors have wrong dimentions" << endl;
      return 1;
    }


}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM::SetFilenameOut(char *f)
{
  this->FilenameOut = f;
  return 0;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM::SetMaskImage(nifti_image *f)
{
  this->Mask = f;
  this->maskImage_status = true;

  if(Mask->datatype!=DT_BINARY){
      seg_convert2binary(Mask,0.0f);
    }

  bool * MaskDataPtr = static_cast<bool *>(Mask->data);
  this->numelmasked=0;
  for(int i=0; i<this->numel; i++,MaskDataPtr++){
      if((*MaskDataPtr)>0){
          this->numelmasked++;
        }
    }
  if(this->nx==f->nx && this->ny==f->ny && this->nz==f->nz){
      return 0;
    }
  else{
      cout << "ERROR: Mask has wrong dimentions" << endl;
      return 1;
    }

  return 0;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM::SetVerbose(unsigned int verblevel)
{
  if(verblevel<0){
      this->verbose_level=0;
    }
  else if(verblevel>2){
      this->verbose_level=2;
    }
  else{
      this->verbose_level=verblevel;
    }
  return 0;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM::SetRegValue(float reg)
{
  if(reg>0){
      this->reg_factor=reg;
    }
  else{
      this->reg_factor=1;
      cout << "Non valid regularization value. Will assume -reg 1."<<endl;
    }
  return 0;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM::SetMaximalIterationNumber(unsigned int numberiter)
{
  if(numberiter<1){
      this->maxIteration=1;
      cout << "Warning: It will only stop at iteration 1"<< endl;
    }
  else{
      this->maxIteration=numberiter;
    }
  return 0;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM::SetAprox(bool aproxval)
{
  this->aprox=aproxval;
  return 0;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM::Turn_MRF_ON(float strength)
{
  this->MRF_status=true;
  this->MRF_strength=strength;
  this->MRF_transitionMatrix=new float[this->numb_classes*this->numb_classes]();
  Create_diagonal_MRF_transitionMatrix();


  if(this->maskImage_status){
      this->MRF=new float[this->numelmasked*this->numb_classes*this->nu*this->nt]();
      for(int i=0; i<(this->numelmasked*this->numb_classes); i++){
          MRF[i]=(float)(1.0);
        }
    }
  else{
      this->MRF=new float[this->numel*this->numb_classes*this->nu*this->nt]();
      for(int i=0; i<(this->numel*this->numb_classes); i++){
          MRF[i]=(float)(1.0);
        }
    }

  Create_diagonal_MRF_transitionMatrix();
  return 0;
}


int seg_EM::Turn_Relaxation_ON(float relax_factor,float relax_gauss_kernel){
  this->Relax_status=true;
  this->Relax_factor=relax_factor;
  this->RelaxGaussKernelSize=relax_gauss_kernel;
  return 1;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM::Turn_BiasField_ON(int powerOrder,float ratiothresh)
{
  this->BiasField_status=true;
  this->BiasField_order=powerOrder;
  this->BiasField_coeficients = new SegPrecisionTYPE[((powerOrder+1)*(powerOrder+2)/2*(powerOrder+3)/3)*this->nu*this->nt]();
  this->BiasField_ratio=ratiothresh;
  if(this->maskImage_status){
      this->BiasField = new SegPrecisionTYPE[this->numelmasked*this->nu*this->nt]();
    }
  else{
      this->BiasField = new SegPrecisionTYPE[this->numel*this->nu*this->nt]();
    }
  return 0;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM::Create_diagonal_MRF_transitionMatrix()
{
  for(int i=0;i<this->numb_classes;i++){
      for(int j=0;j<this->numb_classes;j++){
          if(j==i){
              this->MRF_transitionMatrix[i+j*this->numb_classes]=0;
            }
          else{
              this->MRF_transitionMatrix[i+j*this->numb_classes]=this->MRF_strength;
            }
        }
    }
  return 0;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM::Create_CurrSizes()
{
  this->CurrSizes = new ImageSize [1]();
  CurrSizes->numel=(int)(this->nx*this->ny*this->nz);
  CurrSizes->xsize=this->nx;
  CurrSizes->ysize=this->ny;
  CurrSizes->zsize=this->nz;
  CurrSizes->usize=(this->nu>1)?this->nu:1;
  CurrSizes->tsize=(this->nt>1)?this->nt:1;
  CurrSizes->numclass=this->numb_classes;
  CurrSizes->numelmasked=0;
  CurrSizes->numelbias=0;
  return 0;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
int seg_EM::OutliernessON(float in_OutliernessThreshold, float ratio){

  this->OutliernessFlag=true;
  this->OutliernessThreshold=in_OutliernessThreshold;
  this->Outlierness_ratio=ratio;
  return 0;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
int seg_EM::Maximization()
{
  if(this->MAP_status){
      if(this->maskImage_status){
          calcM_mask(this->inputImage,this->Expec,this->BiasField,this->Outlierness,this->Short_2_Long_Indices,M,V,this->MAP_M,this->MAP_V,this->reg_factor,CurrSizes,this->verbose_level);
        }
      else{
          calcM(this->inputImage,this->Expec,this->BiasField,this->Outlierness,M,V,this->MAP_M,this->MAP_V,this->reg_factor,CurrSizes,this->verbose_level);
        }
    }
  else{
      if(this->maskImage_status){
          calcM_mask(this->inputImage,this->Expec,this->BiasField,this->Outlierness,this->Short_2_Long_Indices,M,V,NULL,NULL,this->reg_factor,CurrSizes,this->verbose_level);
        }
      else{
          calcM(this->inputImage,this->Expec,this->BiasField,this->Outlierness,M,V,NULL,NULL,this->reg_factor,CurrSizes,this->verbose_level);
        }
    }
  return 1;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM::Expectation()
{

  if(this->ratio<(SegPrecisionTYPE)(this->Outlierness_ratio) && this->iter>3 && this->OutliernessUSE==NULL && this->OutliernessFlag){
      this->OutliernessUSE=this->Outlierness;
      if(this->verbose_level>0){
          cout << "Updating Outlierness - LogRatio = "<<ratio<<endl;
        }
      this->loglik=2;
      this->oldloglik=1;
      this->checkpoint_iter=this->iter+2;
    }

  //update
  if(this->MRF_status){
      if(this->maskImage_status){
          calcE_mask(this->inputImage,this->MRF,this->Expec,&this->loglik,this->BiasField,this->OutliernessUSE,this->OutliernessThreshold,this->Short_2_Long_Indices,this->M,this->V,CurrSizes,this->verbose_level);
        }
      else{
          calcE(this->inputImage,this->MRF,this->Expec,&this->loglik,this->BiasField,this->OutliernessUSE,this->OutliernessThreshold,this->M,this->V,CurrSizes,this->verbose_level);
        }
    }
  else{

      if(this->maskImage_status){
          calcE_mask(this->inputImage,this->ShortPrior,this->Expec,&this->loglik,this->BiasField,this->OutliernessUSE,this->OutliernessThreshold,this->Short_2_Long_Indices,this->M,this->V,CurrSizes,this->verbose_level);
        }
      else{
          calcE(this->inputImage,this->ShortPrior,this->Expec,&this->loglik,this->BiasField,this->OutliernessUSE,this->OutliernessThreshold,this->M,this->V,CurrSizes,this->verbose_level);
        }
    }
  //}
  return 0;

}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM::UpdateMRF()
{
  if(this->nz>1){
      if(this->MRF_status){
          if(this->maskImage_status){
              MRFregularization_mask(this->Expec,this->MRF_transitionMatrix,this->MRF_transitionMatrix,this->MRF_beta,this->MRF,this->ShortPrior,this->Long_2_Short_Indices,this->Short_2_Long_Indices,this->CurrSizes,this->MRF_status, this->verbose_level);
            }
          else{
              MRFregularization(this->Expec,this->MRF_transitionMatrix,this->MRF_transitionMatrix,this->MRF_beta,this->MRF,this->ShortPrior,this->CurrSizes,this->MRF_status, this->verbose_level);
            }
        }
    }
  else if(this->nz==1){
      if(this->MRF_status){
          if(this->maskImage_status){
              MRFregularization_mask2D(this->Expec,this->MRF_transitionMatrix,this->MRF_transitionMatrix,this->MRF_beta,this->MRF,this->ShortPrior,this->Long_2_Short_Indices,this->Short_2_Long_Indices,this->CurrSizes,this->MRF_status, this->verbose_level);
            }
          else{
              MRFregularization2D(this->Expec,this->MRF_transitionMatrix,this->MRF_transitionMatrix,this->MRF_beta,this->MRF,this->ShortPrior,this->CurrSizes,this->MRF_status, this->verbose_level);
            }
        }
    }
  return 1;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM::UpdateBiasField()
{
  if(this->BiasField_status){
      if((((this->loglik-this->oldloglik)/fabs(this->oldloglik))<(SegPrecisionTYPE)(this->BiasField_ratio) && this->iter>3)||((SegPrecisionTYPE)(this->BiasField_ratio)==0.0f)){
          if(this->maskImage_status){
              if(this->nz>1){

                  BiasCorrection_mask(this->BiasField,this->BiasField_coeficients,this->inputImage,this->Long_2_Short_Indices,Expec,this->OutliernessUSE,this->M,this->V,this->BiasField_order,CurrSizes,this->BiasField_status,this->verbose_level);
                }
              else{
                  BiasCorrection_mask2D(this->BiasField,this->BiasField_coeficients,this->inputImage,this->Long_2_Short_Indices,Expec,this->OutliernessUSE,this->M,this->V,this->BiasField_order,CurrSizes,this->BiasField_status,this->verbose_level);
                }
            }
          else{
              BiasCorrection(this->BiasField,this->BiasField_coeficients,this->inputImage,Expec,this->OutliernessUSE,this->M,this->V,this->BiasField_order,CurrSizes,this->BiasField_status,this->verbose_level);
            }
        }
    }
  return 1;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM::UpdatePriorWeight()
{
  if(this->Relax_status){

      if(this->maskImage_status){
          if((int)(this->verbose_level)>(int)(0)){
              cout << "Relaxing Priors"<< endl;
            }
          PriorWeight_mask(this->ShortPrior,this->Priors,this->Expec,this->RelaxGaussKernelSize,this->Relax_factor,this->Short_2_Long_Indices,this->Long_2_Short_Indices,CurrSizes,this->verbose_level);
        }
      else{
          if((int)(this->verbose_level)>(int)(0)){
              cout << "Relaxing Priors only available on masked images"<< endl;
            }
          // PriorWeight(this->ShortPrior,this->Priors,this->Expec,CurrSizes,this->verbose_level);
        }
    }
  return 1;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM::Normalize_Image_and_Priors()
{
  if(this->maskImage_status){
      if(this->Priors_status)Normalize_NaN_Priors_mask(this->Priors,this->Mask,this->verbose_level);
      Normalize_Image_mask(this->inputImage,this->Mask,CurrSizes,this->verbose_level);

      this->Short_2_Long_Indices = Create_Short_2_Long_Matrix_from_NII(this->Mask,&(CurrSizes->numelmasked));
      this->Long_2_Short_Indices = Create_Long_2_Short_Matrix_from_NII(this->Mask);
      this->numelmasked=(CurrSizes->numelmasked);
    }
  else{
      if(this->Priors_status)Normalize_NaN_Priors(this->Priors,this->verbose_level);

      Normalize_Image(this->inputImage,CurrSizes,this->verbose_level);
    }

  if(this->MAP_status){
      for(int i=0;i<this->numb_classes;i++){
          this->MAP_M[i]=logf(((this->MAP_M[i]-CurrSizes->rescale_min[0])/(CurrSizes->rescale_max[0]-CurrSizes->rescale_min[0]))+1)/0.693147181;;
          if(this->verbose_level>0){
              cout << "MAP_M["<<i<<"] = "<<this->MAP_M[i]<< endl;
            }
          this->M[i]=this->MAP_M[i];
          this->V[i]=1.0/this->numb_classes;
        }
    }
  return 0;

}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM::Allocate_and_Initialize()
{


  if(this->Priors_status){
      if(this->maskImage_status){
          this->Expec = Create_cArray_from_Prior_mask(this->Mask,this->Priors,CurrSizes->numclass,this->PV_model_status);
          this->ShortPrior = Create_cArray_from_Prior_mask(this->Mask,this->Priors,CurrSizes->numclass,this->PV_model_status);
        }
      else{
          this->Expec = Create_cArray_from_Prior(this->Priors,CurrSizes->numclass,this->PV_model_status);
          this->ShortPrior = Create_cArray_from_Prior(this->Priors,CurrSizes->numclass,this->PV_model_status);
        }
    }
  else{
      Intensity_Based_Inisitalization_of_Means();
      int tmpnumb_elem=0;
      if(this->maskImage_status){
          tmpnumb_elem=(this->numelmasked*(this->numb_classes+(int)(this->PV_model_status)*2));
        }
      else{
          tmpnumb_elem=(numel*(this->numb_classes+(int)(this->PV_model_status)*2));
        }

      float tmpnumb_clas=((this->numb_classes+(int)(this->PV_model_status)*2));
      this->Expec=new SegPrecisionTYPE [tmpnumb_elem] ();
      this->ShortPrior=new SegPrecisionTYPE [tmpnumb_elem] ();
      for(int i=0; i<tmpnumb_elem; i++){
          this->Expec[i]=1.0/tmpnumb_clas;
          this->ShortPrior[i]=1.0/tmpnumb_clas;
        }
      if(this->maskImage_status){
          calcE_mask(this->inputImage,this->ShortPrior,this->Expec,&this->loglik,this->BiasField,NULL,0,this->Short_2_Long_Indices,this->M,this->V,CurrSizes,this->verbose_level);
        }
      else{
          calcE(this->inputImage,this->ShortPrior,this->Expec,&this->loglik,this->BiasField,NULL,0,this->M,this->V,CurrSizes,this->verbose_level);
        }

    }

  if(this->OutliernessFlag){
      int tmpnumb_elem=0;
      if(this->maskImage_status){
          tmpnumb_elem=(this->numelmasked*(this->numb_classes+(int)(this->PV_model_status)*2));
        }
      else{
          tmpnumb_elem=(numel*(this->numb_classes+(int)(this->PV_model_status)*2));
        }
      this->OutliernessUSE=NULL;
      this->Outlierness=new SegPrecisionTYPE [tmpnumb_elem] ();
      for(int i=0; i<tmpnumb_elem; i++){
          this->Outlierness[i]=1.0;
        }
    }

  return 0;

}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
int seg_EM::Intensity_Based_Inisitalization_of_Means()
{

  for(int ms=0; ms<this->nu; ms++){
      SegPrecisionTYPE * Intensity_PTR = static_cast<SegPrecisionTYPE *>(this->inputImage->data);
      bool * MaskDataPtr=NULL;
      if(this->maskImage_status){
          MaskDataPtr = static_cast<bool *>(this->Mask->data);
        }
      int mycounter=0;
      float meanval=0.0;
      float variance=0.0;

      for(int i=0; i<this->numel; i++){
          if(!this->maskImage_status || MaskDataPtr[i]>0){
              mycounter++;
              meanval+=(Intensity_PTR[i+ms*this->numel]);
            }
        }
      meanval=meanval/mycounter;

      for(int i=0; i<this->numel; i++){
          if(!this->maskImage_status || MaskDataPtr[i]>0 ){
              variance+=pow((meanval-Intensity_PTR[i+ms*this->numel]),2);
            }
        }
      variance=variance/mycounter;
      int histogram[1001];
      for(int i=0;i<1000;i++){
          histogram[i]=0;
        }
      float tmpmax=-1e32f;
      float tmpmin=1e32f;



      for(int i=0; i<this->numel; i++){
          if(!this->maskImage_status || MaskDataPtr[i]>0){
              if(tmpmax<(int)(Intensity_PTR[i+ms*this->numel])){
                  tmpmax=(int)(Intensity_PTR[i+ms*this->numel]);
                }
              if(tmpmin>(int)(Intensity_PTR[i+ms*this->numel])){
                  tmpmin=(int)(Intensity_PTR[i+ms*this->numel]);
                }
            }
        }

      for(int i=0; i<this->numel; i++){
          if(!this->maskImage_status || MaskDataPtr[i]>0){

              int index4hist=(int)(1000.0f*(float)(Intensity_PTR[i+ms*this->numel]-tmpmin)/(float)(tmpmax-tmpmin));
              if((index4hist>1000) & (index4hist<0)){
                  cout<< "error"<<endl;
                }
              histogram[(int)(1000.0*(float)(Intensity_PTR[i+ms*this->numel]-tmpmin)/(float)(tmpmax-tmpmin))]++;
            }
        }


      for(int clas=0; clas<this->numb_classes; clas++){
          float tmpsum=0;
          int tmpindex=0;
          float percentile=((float)clas+1)/(this->numb_classes+1);
          for(int i=999;i>0;i--){
              tmpsum+=histogram[i];
              tmpindex=i;
              if((float)(tmpsum)>((1.0f-percentile)*(float)(mycounter))){
                  i=0;
                }
            }
          M[clas*CurrSizes->usize+ms]=float(tmpindex)*(tmpmax-tmpmin)/1000.0f+(tmpmin);
          V[clas*CurrSizes->usize*CurrSizes->usize+ms*CurrSizes->usize+ms]=variance/this->numb_classes/2;
        }
    }


  for (int cl=0; cl<this->numb_classes; cl++) {
      if(this->verbose_level>0){
          if(CurrSizes->usize==1){
              cout.fill('0');
              cout<< "M["<<(int)(cl)<<"]= "<<setw(10)<<setprecision(7)<<left<<(SegPrecisionTYPE)(M[cl])<<"\tV["<<(int)(cl)<<"]="<<setw(10)<<setprecision(7)<<left<<(SegPrecisionTYPE)(V[cl])<< endl;
              flush(cout);
            }
          else{

              cout<< "M["<<(int)(cl)<<"]= ";
              for(int Multispec=0; Multispec<CurrSizes->usize; Multispec++) {
                  cout<< setw(10)<<setprecision(7)<<left<<(SegPrecisionTYPE)(M[cl*CurrSizes->usize+Multispec])<<"\t";
                }
              cout<< endl;
              flush(cout);
              cout<< "V["<<(int)(cl)<<"]= ";
              for(int Multispec=0; Multispec<CurrSizes->usize; Multispec++) {
                  if(Multispec>0){
                      cout<< "      ";
                    }
                  for(int Multispec2=0; Multispec2<CurrSizes->usize; Multispec2++) {
                      cout<< setw(10)<<setprecision(7)<<left<<(SegPrecisionTYPE)(V[cl*CurrSizes->usize*CurrSizes->usize+Multispec*CurrSizes->usize+Multispec2])<<"\t";
                    }
                  cout<< endl;
                }
              cout<< endl;
              flush(cout);
            }
        }
    }

  return 0;
}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

nifti_image * seg_EM::GetResult()
{
  nifti_image * Result=NULL;

  if(this->maskImage_status){
      Result = Copy_Expec_to_Result_mask(this->Expec,this->Short_2_Long_Indices,this->inputImage,(char*)this->FilenameOut.c_str(),this->CurrSizes);
    }
  else{
      Result = Copy_Expec_to_Result(this->Expec,this->inputImage,(char*)this->FilenameOut.c_str(),this->CurrSizes);
    }
  return Result;

}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

nifti_image * seg_EM::GetResultNeonate()
{
  nifti_image * Result=NULL;

  Result = Copy_Expec_to_Result_Neonate_mask(this->Expec,this->Short_2_Long_Indices,this->Long_2_Short_Indices,this->inputImage,this->BiasField,this->M,(char*)this->FilenameOut.c_str(),this->CurrSizes);

  return Result;

}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

nifti_image * seg_EM::GetBiasCorrected(char * filename)
{
  nifti_image * Result=NULL;

  if(this->maskImage_status){
      Result = Get_Bias_Corrected_mask(this->BiasField_coeficients,this->inputImage,filename,this->CurrSizes,this->BiasField_order);
    }
  else{
      Result = Get_Bias_Corrected(this->BiasField,this->inputImage,filename,this->CurrSizes);
    }
  return Result;

}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

nifti_image * seg_EM::GetOutlierness(char * filename)
{
  nifti_image * Result = nifti_copy_nim_info(this->inputImage);
  Result->dim[0]=3;
  Result->dim[4]=1;
  Result->dim[5]=1;
  Result->datatype=DT_FLOAT32;
  Result->cal_max=1;
  nifti_set_filenames(Result,filename,0,0);
  nifti_update_dims_from_array(Result);
  nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
  Result->data = (void *) calloc(Result->nvox, sizeof(SegPrecisionTYPE));
  SegPrecisionTYPE * Resultdata = static_cast<SegPrecisionTYPE *>(Result->data);
  for(unsigned int i=0; i<Result->nvox; i++){Resultdata[i]=0;}



  if(this->maskImage_status){
      for(int i=0; i<CurrSizes->numelmasked; i++){
          float currsum=0;
          for(int currclass=0; currclass<CurrSizes->numclass;currclass++){
              currsum+=this->Outlierness[i+(currclass)*CurrSizes->numelmasked]*Expec[i+(currclass)*CurrSizes->numelmasked];
            }
          Resultdata[Short_2_Long_Indices[i]]=1-currsum;
        }

    }
  else{
      int class_nvox=Result->nx*Result->ny*Result->nz;
      for(int i=0; i<CurrSizes->numel; i++){
          float currsum=0;
          for(int currclass=0; currclass<CurrSizes->numclass;currclass++){
              currsum+=this->Outlierness[i+(currclass)*class_nvox];
            }
          Resultdata[i]=1-currsum;
        }
    }


  return Result;

}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */


int *  seg_EM::Run_EM()
{
  time_t start,end;
  time(&start);
  if((int)(this->verbose_level)>(int)(0)){
      cout << "EM: Verbose level " << this->verbose_level << endl;
    }
  if(!Priors_status){

    }

  if(this->CurrSizes==NULL) {
      Create_CurrSizes();
    }

  Normalize_Image_and_Priors();
  Allocate_and_Initialize();
  if((int)(this->verbose_level)>(int)(0)){
      cout << "Number of voxels inside the mask = " << this->numelmasked << endl;
    }
  //**************
  // EM Algorithm
  //**************
  //bool MRFreset=0;
  this->iter=0;
  bool out= true;

  while (out) {

      if(this->verbose_level>0){
          cout << endl << "*******************************" << endl;
          cout << "Iteration " << iter << endl;
        }

      // Iterative Components - EM, MRF, Bias Correction

      //Maximization
      Maximization();
      //Expectation
      Expectation();
      //MRF
      UpdateMRF();
      //Bias Correction
      UpdateBiasField();
      //Update Weight
      UpdatePriorWeight();

      // Print LogLik depending on the verbose level
      if(this->verbose_level>0 && this->iter>0){
          printloglik(iter,this->loglik,this->oldloglik);
        }
      // Preform MRF reset or Exit
      if((((this->loglik-this->oldloglik)/fabs(this->oldloglik))<(SegPrecisionTYPE)(0.0005) && this->iter>3 && this->iter>this->checkpoint_iter) || iter>=this->maxIteration || (isinf(this->loglik) && this->iter>3)){
          out=false;
        }
      this->ratio=((this->loglik-this->oldloglik)/fabs(this->oldloglik));
      // Update LogLik
      this->oldloglik=this->loglik;
      iter++;
    }

  time(&end);

  if(this->verbose_level>0){
      int minutes = (int)floorf(float(end-start)/60.0f);
      int seconds = (int)(end-start - 60*minutes);
      cout << "Finished in "<<minutes<<"min "<<seconds<<"sec"<< endl;
    }
  return 0;
}

#endif
