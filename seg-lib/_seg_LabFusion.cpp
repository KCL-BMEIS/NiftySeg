#include "_seg_LabFusion.h"




seg_LabFusion::seg_LabFusion(int _numb_classif, int numbclasses, int _Numb_Neigh)
{
  TYPE_OF_FUSION=0;
  NumberOfLabels=numbclasses;
  inputCLASSIFIER = NULL; // pointer to external
  inputImage_status=0;
  verbose_level=0;

  // Size
  this->dimentions=4;
  this->nx=0;
  this->ny=0;
  this->nz=0;
  this->nt=0;
  this->nu=0;
  this->dx=0;
  this->dy=0;
  this->dz=0;
  this->Conv=0.000001f;
  this->numel=0;
  this->iter=0;
  this->CurrSizes=NULL;
  this->Numb_Neigh=_Numb_Neigh;

  this->Thresh_IMG_value=0;
  this->Thresh_IMG_DO=false;


  // SegParameters
  for (int i=0;i<(5000);i++){
      LableCorrespondences_big_to_small[i]=0;
      LableCorrespondences_small_to_big[i]=0;
    }
  ConfusionMatrix=new LabFusion_datatype [_numb_classif*numbclasses*numbclasses];
  if(ConfusionMatrix == NULL){
      fprintf(stderr,"* Error when alocating ConfusionMatrix: Not enough memory\n");
      exit(1);
    }

  for (int i=0;i<(_numb_classif);i++){
      for (int j=0;j<(numbclasses);j++){
          for (int k=0;k<(numbclasses);k++){
              if(k==j){
                  ConfusionMatrix[k+j*numbclasses+i*numbclasses*numbclasses]=0.80;
                }
              else{
                  ConfusionMatrix[k+j*numbclasses+i*numbclasses*numbclasses]=0.2/(numbclasses-1);
                }
            }
        }
    }

  this->W=NULL;
  this->FinalSeg=NULL;
  this->maskAndUncertainIndeces=NULL;
  this->sizeAfterMaskingAndUncertainty=0;
  this->uncertainflag=false;
  this->dilunc=0;
  this->Prop=new LabFusion_datatype [numbclasses];
  this->Fixed_Prop_status=false;
  this->PropUpdate=false;
  this->numb_classif=_numb_classif;
  this->loglik=1;
  this->oldloglik=0;
  this->maxIteration=100;
  this->LNCC=NULL;
  this->NCC=NULL;
  this->Numb_Neigh=3;
  this->LNCC_status=false;
  this->NCC_status=false;


  // MRF Specific
  this->MRF_status=0;
  this->MRF_strength=0.0f;
  this->MRF_matrix=new LabFusion_datatype [numbclasses*numbclasses];
  if(MRF_matrix == NULL){
      fprintf(stderr,"* Error when alocating MRF_matrix: Not enough memory\n");
      exit(1);
    }

  for (int j=0;j<(numbclasses);j++){
      for (int k=0;k<(numbclasses);k++){
          if(k==j){
              this->MRF_matrix[k+j*numbclasses]=0;
            }
          else{
              this->MRF_matrix[k+j*numbclasses]=1;
            }
        }
    }

  MRF=NULL;

}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

seg_LabFusion::~seg_LabFusion()
{

  if(this->W!=NULL)
    delete[] this->W;
  if(this->MRF!=NULL)
    delete[] this->MRF;
  if(this->CurrSizes!=NULL)
    delete [] this->CurrSizes;
  if(this->NCC!=NULL)
    delete [] this->NCC;
  if(this->LNCC!=NULL)
    delete [] this->LNCC;
  if(this->Prop!=NULL)
    delete [] this->Prop;
  if(this->ConfusionMatrix!=NULL)
    delete [] this->ConfusionMatrix;
  if(this->MRF_matrix!=NULL)
    delete [] this->MRF_matrix;
  if(this->maskAndUncertainIndeces!=NULL)
    delete [] this->maskAndUncertainIndeces;
  if(this->FinalSeg!=NULL)
    delete [] this->FinalSeg;

  this->W=NULL;
  this->MRF=NULL;
  this->NCC=NULL;
  this->LNCC=NULL;
  this->Prop=NULL;
  this->CurrSizes=NULL;
  this->ConfusionMatrix=NULL;
  this->MRF_matrix=NULL;
  this->FinalSeg=NULL;
  this->maskAndUncertainIndeces=NULL;

}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::SetMask(nifti_image *Mask)
{

  bool * inputMaskPtr=static_cast<bool *>(Mask->data);
  this->sizeAfterMaskingAndUncertainty=0;
  for(int i=0;i<(this->numel);i++){
      if(inputMaskPtr[i]>0){
          this->maskAndUncertainIndeces[i]=this->sizeAfterMaskingAndUncertainty;
          this->sizeAfterMaskingAndUncertainty++;
        }
      else{
          this->maskAndUncertainIndeces[i]=-1;
        }
    }
  return 0;
}

int seg_LabFusion::SetinputCLASSIFIER(nifti_image *r,bool UNCERTAINflag)
{

  //if(r->datatype!=DT_BINARY){
  //    seg_convert2binary(r,0.5f);
  //}

  this->inputCLASSIFIER = r;
  this->inputImage_status = true;
  this->uncertainflag=UNCERTAINflag;
  // Size
  this->dimentions=(int)((r->nx)>1)+(int)((r->ny)>1)+(int)((r->nz)>1)+(int)((r->nt)>1)+(int)((r->nu)>1);
  this->nx=r->nx;
  this->ny=r->ny;
  this->nz=r->nz;
  this->dx=r->dx;
  this->dy=r->dy;
  this->dz=r->dz;
  this->Numb_Neigh=r->nt;
  this->numel=r->nz*r->ny*r->nx;
  if(this->CurrSizes==NULL) Create_CurrSizes();

  classifier_datatype * CLASSIFIERptr = static_cast<classifier_datatype *>(this->inputCLASSIFIER->data);
  int * NumberOfDifferentClassesHistogram=new int [5000];
  if(NumberOfDifferentClassesHistogram == NULL){
      fprintf(stderr,"* Error when alocating NumberOfDifferentClassesHistogram: Not enough memory\n");
      exit(1);
    }

  for(int i=0;i<5000;i++){
      NumberOfDifferentClassesHistogram[i]=0;
    }
  for(int i=0;i<(int)this->inputCLASSIFIER->nvox;i++){
      NumberOfDifferentClassesHistogram[CLASSIFIERptr[i]]++;
    }

  int currlabindex=0;
  for(int i=0;i<5000;i++){
      if(NumberOfDifferentClassesHistogram[i]>0){
          LableCorrespondences_small_to_big[currlabindex]=i;
          LableCorrespondences_big_to_small[i]=currlabindex;
          currlabindex++;
        }
      else{
          LableCorrespondences_big_to_small[i]=-1;
        }
    }
  for(int i=0;i<(int)this->inputCLASSIFIER->nvox;i++){
      CLASSIFIERptr[i]=LableCorrespondences_big_to_small[CLASSIFIERptr[i]];
    }


  this->maskAndUncertainIndeces=new int [this->numel];
  if(maskAndUncertainIndeces == NULL){
      fprintf(stderr,"* Error when alocating mask_and_uncertain_indexes: Not enough memory\n");
      exit(1);
    }

  this->sizeAfterMaskingAndUncertainty=this->numel;
  for(int i=0;i<(this->numel);i++){
      this->maskAndUncertainIndeces[i]=i;
    }
  delete [] NumberOfDifferentClassesHistogram;

  return 0;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::SetLNCC(nifti_image * _LNCC,nifti_image * BaseImage,LabFusion_datatype distance,int Numb_Neigh)
{
  if((Numb_Neigh<_LNCC->nt) & (Numb_Neigh>0)){
      this->Numb_Neigh=(int)Numb_Neigh;
    }
  else{
      this->Numb_Neigh=_LNCC->nt;
    }

  if(this->nx!=_LNCC->nx || this->ny!=_LNCC->ny || this->nz!=_LNCC->nz){
      fprintf(stderr,"* The image size of the images do not match");
      exit(1);
    }
  if(this->nx!=BaseImage->nx || this->ny!=BaseImage->ny || this->nz!=BaseImage->nz){
      fprintf(stderr,"* The image size of the images do not match");
      exit(1);
    }

  if(_LNCC->nt==this->numb_classif){
      if(_LNCC->datatype==DT_FLOAT){
          this->LNCC=seg_norm4LNCC(BaseImage,_LNCC,distance,Numb_Neigh,CurrSizes,this->verbose_level);
          this->LNCC_status = true;

        }
      else{
          cout << "LNCC is not LabFusion_datatype"<< endl;
        }
    }
  return 0;
}



/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

/*int seg_LabFusion::SetLMETRIC(nifti_image * _InMetric,int Numb_Neigh)
{
  if((Numb_Neigh<_InMetric->nt) & (Numb_Neigh>0)){
      this->Numb_Neigh=(int)Numb_Neigh;
    }
  else{
      this->Numb_Neigh=_InMetric->nt;
    }

  if(this->nx!=_InMetric->nx || this->ny!=_InMetric->ny || this->nz!=_InMetric->nz){
      fprintf(stderr,"* The image size of the images do not match");
      exit(1);
    }


  if(_InMetric->nt==this->numb_classif){
      if(_InMetric->datatype==DT_FLOAT){
          this->LNCC=new char [Numb_Neigh*this->nx*this->ny*this->nz];
          if(this->LNCC == NULL){
              fprintf(stderr,"* Error when alocating this->LNCC in function SetLMETRIC");
              exit(-1);
            }

          LabFusion_datatype * InputMetricPtr = static_cast<LabFusion_datatype *>(_InMetric->data);

#ifdef _OPENMP
#pragma omp parallel for default(none) \
  shared(LNCC,BaseImage,LNCCptr,CurrSizes,numberordered,LNCC_ordered)
#endif
          for(int i=0; i<this->nx*this->ny*this->nz;i++){
              LabFusion_datatype * LNCCvalue_tmp = new LabFusion_datatype [_InMetric->nt];
              for(int currlable=0;currlable<_InMetric->nt; currlable++){
                  LNCCvalue_tmp[currlable]=InputMetricPtr[i+currlable*this->nx*this->ny*this->nz];
                }
              int * ordertmp=quickSort_order(&LNCCvalue_tmp[0],_InMetric->nt);
              for(int lable_order=0;lable_order<Numb_Neigh;lable_order++){
                  this->LNCC[i+lable_order*this->nx*this->ny*this->nz]=(char)ordertmp[_InMetric->nt-lable_order-1];
                }
              delete [] ordertmp;
            }
          if (this->verbose_level>0){
              cout << "Finished sorting Metric"<< endl;
              flush(cout);
            }
        }
      else{
          cout << "Metric is not float"<< endl;
        }
    }
  return 0;
}
*/

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
int seg_LabFusion::SetMLLNCC(nifti_image * _LNCC,nifti_image * BaseImage,LabFusion_datatype distance,int levels,int Numb_Neigh)
{
  if((Numb_Neigh<_LNCC->nt) & (Numb_Neigh>0)){
      this->Numb_Neigh=(int)Numb_Neigh;
    }
  else{
      this->Numb_Neigh=_LNCC->nt;
    }

  if(this->nx!=_LNCC->nx || this->ny!=_LNCC->ny || this->nz!=_LNCC->nz){
      fprintf(stderr,"* The image size of the images do not match");
      exit(1);
    }
  if(this->nx!=BaseImage->nx || this->ny!=BaseImage->ny || this->nz!=BaseImage->nz){
      fprintf(stderr,"* The image size of the images do not match");
      exit(1);
    }

  if(_LNCC->nt==this->numb_classif){
      if(_LNCC->datatype==DT_FLOAT){
          this->LNCC=seg_norm4MLLNCC(BaseImage,_LNCC,distance,levels,Numb_Neigh,CurrSizes,this->verbose_level);

          this->LNCC_status = true;

        }
      else{
          cout << "LNCC is not LabFusion_datatype"<< endl;
        }
    }
  return 0;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
int seg_LabFusion::SetGNCC(nifti_image * _GNCC,nifti_image * BaseImage,int Numb_Neigh)
{
  if((Numb_Neigh<_GNCC->nt) & (Numb_Neigh>0)){
      this->Numb_Neigh=(int)Numb_Neigh;
    }
  else{
      this->Numb_Neigh=_GNCC->nt;
    }

  if(this->nx!=_GNCC->nx || this->ny!=_GNCC->ny || this->nz!=_GNCC->nz){
      fprintf(stderr,"* The image size of the images do not match");
      exit(1);
    }
  if(this->nx!=BaseImage->nx || this->ny!=BaseImage->ny || this->nz!=BaseImage->nz){
      fprintf(stderr,"* The image size of the images do not match");
      exit(1);
    }

  if(_GNCC->nt==this->numb_classif){
      if(_GNCC->datatype==DT_FLOAT){
          this->NCC=seg_norm_4D_GNCC(BaseImage,_GNCC,Numb_Neigh,CurrSizes,this->verbose_level);
          this->NCC_status = true;
        }
      else{
          cout << "GNCC is not LabFusion_datatype"<< endl;
        }
    }

  return 0;
}

int seg_LabFusion::SetROINCC(nifti_image * _ROINCC,nifti_image * BaseImage,int Numb_Neigh,int DilSize)
{
  if((Numb_Neigh<_ROINCC->nt) & (Numb_Neigh>0)){
      this->Numb_Neigh=(int)Numb_Neigh;
    }
  else{
      this->Numb_Neigh=_ROINCC->nt;
    }


  if(this->nx!=_ROINCC->nx || this->ny!=_ROINCC->ny || this->nz!=_ROINCC->nz){
      fprintf(stderr,"* The image size of the images do not match");
      exit(1);
    }
  if(this->nx!=BaseImage->nx || this->ny!=BaseImage->ny || this->nz!=BaseImage->nz){
      fprintf(stderr,"* The image size of the images do not match");
      exit(1);
    }

  if(_ROINCC->nt==this->numb_classif){
      if(_ROINCC->datatype==DT_FLOAT){
          this->NCC=seg_norm4ROINCC(this->inputCLASSIFIER,BaseImage,_ROINCC,Numb_Neigh,CurrSizes,DilSize,this->verbose_level);
          this->NCC_status = true;
        }
      else{
          cout << "ROINCC is not LabFusion_datatype"<< endl;
        }
    }
  return 0;
}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::SetProp(LabFusion_datatype r)
{
  this->Prop[0] = (LabFusion_datatype)r;
  this->Fixed_Prop_status = true;
  return 0;
}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::SetConv(LabFusion_datatype r)
{
  if(this->verbose_level){
      cout<< "Convergence Ratio = " << r <<endl;
      flush(cout);
    }
  this->Conv = (LabFusion_datatype)r;
  return 0;
}

int seg_LabFusion::SetImgThresh(LabFusion_datatype _Thresh_IMG_value)
{
  this->Thresh_IMG_value=_Thresh_IMG_value;
  this->Thresh_IMG_DO=true;
  return 0;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::SetFilenameOut(char *f)
{
  this->FilenameOut = f;
  return 0;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::SetVerbose(unsigned int verblevel)
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

int seg_LabFusion::SetDilUnc(int _dilunc)
{
  this->dilunc=_dilunc;
  return 0;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::SetMaximalIterationNumber(unsigned int numberiter)
{
  if(numberiter<3){
      this->maxIteration=1;
      cout << "Warning: It will only stop at iteration 1. For less than 1 iteration, use majority voting."<< endl;
    }
  else{
      this->maxIteration=numberiter;
    }
  return 0;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::Turn_MRF_ON(LabFusion_datatype strength)
{
  this->MRF_status=true;
  this->MRF_strength=(LabFusion_datatype)strength;
  return 0;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::Create_CurrSizes()
{
  this->CurrSizes = new ImageSize [1]();
  if(CurrSizes == NULL){
      fprintf(stderr,"* Error when alocating CurrSizes: Not enough memory\n");
      exit(1);
    }
  CurrSizes->numel=(int)(this->nx*this->ny*this->nz);
  CurrSizes->xsize=this->nx;
  CurrSizes->ysize=this->ny;
  CurrSizes->zsize=this->nz;
  CurrSizes->usize=1;
  CurrSizes->tsize=1;
  CurrSizes->numclass=this->numb_classif;
  CurrSizes->numelmasked=0;
  CurrSizes->numelbias=0;
  return 0;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */



int seg_LabFusion::STAPLE_STEPS_Multiclass_Expectation_Maximization()
{

  LabFusion_datatype *tmpW=new LabFusion_datatype [this->NumberOfLabels];
  if(tmpW == NULL){
      fprintf(stderr,"* Error when alocating tmpW: Not enough memory\n");
      exit(1);
    }

  classifier_datatype * inputCLASSIFIERptr = static_cast<classifier_datatype *>(this->inputCLASSIFIER->data);

  bool * nccexists_once= new bool [this->CurrSizes->numclass];
  if(nccexists_once == NULL){
      fprintf(stderr,"* Error when alocating nccexists_once: Not enough memory\n");
      exit(1);
    }

  for(int classifier=0; classifier<this->CurrSizes->numclass;classifier++)nccexists_once[classifier]=false;

  LabFusion_datatype * ConfusionMatrix2=new LabFusion_datatype [this->numb_classif*this->NumberOfLabels*this->NumberOfLabels];
  if(ConfusionMatrix2 == NULL){
      fprintf(stderr,"* Error when alocating ConfusionMatrix2: Not enough memory\n");
      exit(1);
    }

  for(int i=0; i<(this->CurrSizes->numclass*this->NumberOfLabels*this->NumberOfLabels);i++){
      ConfusionMatrix2[i]=0;
    }

  if(this->verbose_level>0){
      cout << "Updating the Posteriors and the Performance Parameters"<< endl;
      flush(cout);
    }
  if(this->verbose_level>1){
      cout << "[lab]  =  P \t\tQ"<< endl;
      flush(cout);
    }

  this->loglik=0;
  for(int classifier=0; classifier<this->CurrSizes->numclass;classifier++)nccexists_once[classifier]=false;

  for(int i=0;i<(this->numel);i++){
      if(this->maskAndUncertainIndeces[i]>=0){
          if(this->MRF_status){
              // **************************
              //   Expectation with MRF
              // **************************
              for(int currclass=0; currclass<this->NumberOfLabels; currclass++){
                  tmpW[currclass]=1;
                }

              if(this->LNCC_status){
                  for(int classifier=0; classifier<this->Numb_Neigh;classifier++){
                      int LNCCvalue=(int)(this->LNCC[i+classifier*this->CurrSizes->numel]);
                      if(LNCCvalue>=0 && LNCCvalue<this->CurrSizes->numclass){
                          for(int currclass=0; currclass<this->NumberOfLabels;currclass++){
                              tmpW[currclass]*=((LabFusion_datatype)(this->ConfusionMatrix[(int)inputCLASSIFIERptr[i+LNCCvalue*this->CurrSizes->numel]+currclass*this->NumberOfLabels+LNCCvalue*this->NumberOfLabels*this->NumberOfLabels]));
                            }
                        }
                    }
                }
              else if(this->NCC_status){
                  for(int classifier=0; classifier<this->Numb_Neigh;classifier++){
                      int NCCvalue=(int)(this->NCC[classifier]);
                      if(NCCvalue>=0 && NCCvalue<this->CurrSizes->numclass){
                          for(int currclass=0; currclass<this->NumberOfLabels;currclass++){
                              tmpW[currclass]*=((LabFusion_datatype)(this->ConfusionMatrix[(int)inputCLASSIFIERptr[i+NCCvalue*this->CurrSizes->numel]+currclass*this->NumberOfLabels+NCCvalue*this->NumberOfLabels*this->NumberOfLabels]));
                            }
                        }
                    }
                }
              else{
                  for(int classifier=0; classifier<this->CurrSizes->numclass;classifier++){
                      for(int currclass=0; currclass<this->NumberOfLabels;currclass++){
                          tmpW[currclass]*=((LabFusion_datatype)(this->ConfusionMatrix[(int)inputCLASSIFIERptr[i+classifier*this->CurrSizes->numel]+currclass*this->NumberOfLabels+classifier*this->NumberOfLabels*this->NumberOfLabels]));
                        }
                    }
                }

              LabFusion_datatype sumW=0;
              for(int currclass=0; currclass<this->NumberOfLabels; currclass++){
                  tmpW[currclass]=((LabFusion_datatype)this->Prop[currclass])*tmpW[currclass]*this->MRF[this->maskAndUncertainIndeces[i]+currclass*this->sizeAfterMaskingAndUncertainty];
                  sumW+=tmpW[currclass];
                }

              if(sumW<=0){
                  if(verbose_level>1)cout << "Normalize by zero at voxel - "<<i<<endl;
//                  for(int currclass=0; currclass<this->NumberOfLabels; currclass++){
//                      W[this->maskAndUncertainIndeces[i]+currclass*this->sizeAfterMaskingAndUncertainty]=1/this->NumberOfLabels;
//                    }
                }
              else{
                  for(int currclass=0; currclass<this->NumberOfLabels; currclass++){
                      W[this->maskAndUncertainIndeces[i]+currclass*this->sizeAfterMaskingAndUncertainty]=tmpW[currclass]/sumW;
                    }
                }
              for(int currclass=0; currclass<this->NumberOfLabels; currclass++){
                  if(this->W[this->maskAndUncertainIndeces[i]+currclass*this->sizeAfterMaskingAndUncertainty]<0){
                      this->W[this->maskAndUncertainIndeces[i]+currclass*this->sizeAfterMaskingAndUncertainty]=0;
                    }
                  if(this->W[this->maskAndUncertainIndeces[i]+currclass*this->sizeAfterMaskingAndUncertainty]>1){
                      this->W[this->maskAndUncertainIndeces[i]+currclass*this->sizeAfterMaskingAndUncertainty]=1;
                    }
                  //this->loglik+=W[this->mask_and_uncertain_indexes[i]+currclass*this->size_after_masking_and_uncertainty];
                }
              this->loglik+=sumW;
            }

          else{
              // **************************
              //  Expectation without MRF
              // **************************
              for(int currclass=0; currclass<this->NumberOfLabels; currclass++){
                  tmpW[currclass]=1;
                }

              if(this->LNCC_status){
                  for(int classifier=0; classifier<this->Numb_Neigh;classifier++){
                      int LNCCvalue=(int)(this->LNCC[i+classifier*this->CurrSizes->numel]);
                      if(LNCCvalue>=0 && LNCCvalue<this->CurrSizes->numclass){
                          for(int currclass=0; currclass<this->NumberOfLabels;currclass++){
                              tmpW[currclass]*=((LabFusion_datatype)(this->ConfusionMatrix[(int)inputCLASSIFIERptr[i+LNCCvalue*this->CurrSizes->numel]+currclass*this->NumberOfLabels+LNCCvalue*this->NumberOfLabels*this->NumberOfLabels]));
                            }
                        }
                    }
                }
              else if(this->NCC_status){
                  for(int classifier=0; classifier<this->Numb_Neigh;classifier++){
                      int NCCvalue=(int)(this->NCC[classifier]);
                      if(NCCvalue>=0 && NCCvalue<this->CurrSizes->numclass){
                          for(int currclass=0; currclass<this->NumberOfLabels;currclass++){
                              tmpW[currclass]*=((LabFusion_datatype)(this->ConfusionMatrix[(int)inputCLASSIFIERptr[i+NCCvalue*this->CurrSizes->numel]+currclass*this->NumberOfLabels+NCCvalue*this->NumberOfLabels*this->NumberOfLabels]));
                            }
                        }
                    }
                }
              else{
                  for(int classifier=0; classifier<this->CurrSizes->numclass;classifier++){
                      for(int currclass=0; currclass<this->NumberOfLabels;currclass++){
                          tmpW[currclass]*=((LabFusion_datatype)(this->ConfusionMatrix[(int)inputCLASSIFIERptr[i+classifier*this->CurrSizes->numel]+currclass*this->NumberOfLabels+classifier*this->NumberOfLabels*this->NumberOfLabels]));
                        }
                    }
                }
              LabFusion_datatype sumW=0;
              for(int currclass=0; currclass<this->NumberOfLabels; currclass++){
                  tmpW[currclass]=((LabFusion_datatype)this->Prop[currclass])*tmpW[currclass];
                  sumW+=tmpW[currclass];
                }

              if(sumW<=0){
                  if(verbose_level>1)cout << "Normalize by zero at voxel - "<<i<<endl;
                }
              else{
                  for(int currclass=0; currclass<this->NumberOfLabels; currclass++){
                      W[this->maskAndUncertainIndeces[i]+currclass*this->sizeAfterMaskingAndUncertainty]=tmpW[currclass]/sumW;
                    }
                }



              for(int currclass=0; currclass<this->NumberOfLabels; currclass++){
                  if(this->W[this->maskAndUncertainIndeces[i]+currclass*this->sizeAfterMaskingAndUncertainty]<0){
                      this->W[this->maskAndUncertainIndeces[i]+currclass*this->sizeAfterMaskingAndUncertainty]=0;
                    }
                  if(this->W[this->maskAndUncertainIndeces[i]+currclass*this->sizeAfterMaskingAndUncertainty]>1){
                      this->W[this->maskAndUncertainIndeces[i]+currclass*this->sizeAfterMaskingAndUncertainty]=1;
                    }
                  //this->loglik+=W[uncertainindex+currclass*this->size_after_masking_and_uncertainty];
                }
              this->loglik+=sumW;
            }

          // **************************
          //       MAXIMIZATION
          // **************************

          if(this->LNCC_status){

              for(int classifier=0; classifier<this->CurrSizes->numclass;classifier++){
                  bool lnccexists=false;
                  for(classifier_datatype lnccindex=0;lnccindex<this->Numb_Neigh;lnccindex++){
                      if(classifier==(this->LNCC[i+lnccindex*this->CurrSizes->numel])){
                          lnccexists=true;
                          nccexists_once[classifier]=true;
                          lnccindex=this->Numb_Neigh;
                        }
                    }

                  if(lnccexists){
                      for(int currclass=0; currclass<this->NumberOfLabels; currclass++){

                          ConfusionMatrix2[(int)inputCLASSIFIERptr[i+this->CurrSizes->numel*classifier]
                              +currclass*this->NumberOfLabels+classifier*
                              this->NumberOfLabels*this->NumberOfLabels]+=
                              W[this->maskAndUncertainIndeces[i]+currclass*this->sizeAfterMaskingAndUncertainty];
                        }
                    }
                }
            }
          else if(this->NCC_status){
              bool nccexists=false;
              for(int classifier=0; classifier<this->CurrSizes->numclass;classifier++){
                  for(int nccindex=0;nccindex<this->Numb_Neigh;nccindex++){
                      if(classifier==(this->NCC[nccindex])){
                          nccexists=true;
                          nccindex=this->Numb_Neigh;
                          nccexists_once[classifier]=true;
                        }
                    }

                  if(nccexists){
                      for(int currclass=0; currclass<this->NumberOfLabels; currclass++){

                          ConfusionMatrix2[(int)inputCLASSIFIERptr[i+this->CurrSizes->numel*classifier]
                              +currclass*this->NumberOfLabels+classifier*
                              this->NumberOfLabels*this->NumberOfLabels]+=W[this->maskAndUncertainIndeces[i]+currclass*this->sizeAfterMaskingAndUncertainty];
                        }
                    }
                }
            }
          else{
              for(int classifier=0; classifier<this->CurrSizes->numclass;classifier++){
                  for(int currclass=0; currclass<this->NumberOfLabels; currclass++){

                      ConfusionMatrix2[(int)inputCLASSIFIERptr[i+this->CurrSizes->numel*classifier]
                          +currclass*this->NumberOfLabels+classifier*
                          this->NumberOfLabels*this->NumberOfLabels]+=W[this->maskAndUncertainIndeces[i]+currclass*this->sizeAfterMaskingAndUncertainty];
                    }
                }
            }
        }
    }

  // **************************
  // Normalize Confusion Matrix
  // **************************
  for(int classifier=0; classifier<this->CurrSizes->numclass;classifier++){
      for(int currclass=0; currclass<this->NumberOfLabels; currclass++){
          double sumConf=0;
          for(int currclass2=0; currclass2<this->NumberOfLabels; currclass2++){
              sumConf+=ConfusionMatrix2[currclass2+currclass*this->NumberOfLabels+classifier*
                  this->NumberOfLabels*this->NumberOfLabels];
            }

          if(sumConf<=0){
              if(verbose_level>1){
                  cout << "NOMALISE BY ZERO - "<<classifier<<" , "<< currclass<<endl;
                }
              for(int currclass2=0; currclass2<this->NumberOfLabels; currclass2++){
                  ConfusionMatrix2[currclass2+currclass*this->NumberOfLabels+classifier*
                      this->NumberOfLabels*this->NumberOfLabels]=ConfusionMatrix[currclass2+currclass*this->NumberOfLabels+classifier*
                      this->NumberOfLabels*this->NumberOfLabels];
                }
            }
          else{
              for(int currclass2=0; currclass2<this->NumberOfLabels; currclass2++){
                  ConfusionMatrix2[currclass2+currclass*this->NumberOfLabels+classifier*
                      this->NumberOfLabels*this->NumberOfLabels]=
                      ConfusionMatrix2[currclass2+currclass*this->NumberOfLabels+classifier*
                      this->NumberOfLabels*this->NumberOfLabels]/sumConf;
                }
            }
        }

      if(this->verbose_level>1){
          if(this->NCC_status){
              if(nccexists_once[classifier]){
                  cout <<endl<< "["<<classifier+1<<"]=";
                  for(int currclass=0; currclass<this->NumberOfLabels; currclass++){
                      for(int currclass2=0; currclass2<this->NumberOfLabels; currclass2++){
                          cout <<"\t"<< setprecision(4)<<((ConfusionMatrix2[currclass2+currclass*this->NumberOfLabels+classifier*this->NumberOfLabels*this->NumberOfLabels]<0.0001)?0:ConfusionMatrix2[currclass2+currclass*this->NumberOfLabels+classifier*this->NumberOfLabels*this->NumberOfLabels]);
                        }
                      cout << endl;
                    }
                }
            }
          else if(this->LNCC_status){
              if(nccexists_once[classifier]){
                  cout <<endl<< "["<<classifier+1<<"]=";
                  for(int currclass=0; currclass<this->NumberOfLabels; currclass++){
                      for(int currclass2=0; currclass2<this->NumberOfLabels; currclass2++){
                          cout <<"\t"<< setprecision(4)<<((ConfusionMatrix2[currclass2+currclass*this->NumberOfLabels+classifier*this->NumberOfLabels*this->NumberOfLabels]<0.0001)?0:ConfusionMatrix2[currclass2+currclass*this->NumberOfLabels+classifier*this->NumberOfLabels*this->NumberOfLabels]);
                        }
                      cout << endl;
                    }
                }
            }
          else{
              cout <<endl<< "[ "<<classifier+1<<" ]  =  ";
              cout << endl;
              for(int currclass=0; currclass<this->NumberOfLabels; currclass++){
                  cout << "["<<classifier+1<<"]=";
                  for(int currclass=0; currclass<this->NumberOfLabels; currclass++){
                      for(int currclass2=0; currclass2<this->NumberOfLabels; currclass2++){
                          cout <<"\t"<< setprecision(4)<<((ConfusionMatrix2[currclass2+currclass*this->NumberOfLabels+classifier*this->NumberOfLabels*this->NumberOfLabels]<0.0001)?0:ConfusionMatrix2[currclass2+currclass*this->NumberOfLabels+classifier*this->NumberOfLabels*this->NumberOfLabels]);
                        }
                      cout << endl;
                    }
                }
            }
        }
    }

  for(int i=0; i<(this->numb_classif*this->NumberOfLabels*this->NumberOfLabels);i++){
      ConfusionMatrix[i]=ConfusionMatrix2[i];
    }



  delete [] tmpW;
  delete [] ConfusionMatrix2;
  delete [] nccexists_once;
  return 1;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
int seg_LabFusion::SBA_Estimate()
{

  classifier_datatype * inputCLASSIFIERptr = static_cast<classifier_datatype *>(this->inputCLASSIFIER->data);
  //LabFusion_datatype eps=0.1;
  bool * CurrLableImage=new bool [this->numel];
  if(CurrLableImage == NULL){
      fprintf(stderr,"* Error when alocating CurrLableImage: Not enough memory\n");
      exit(1);
    }
  if(this->verbose_level>0){
      cout<<"Calculating Euclidean Distances"<<endl;
      flush(cout);
    }

  //#ifdef _OPENMP
  //#pragma omp parallel for
  //#endif
  for(int classifier=0; classifier<this->numb_classif;classifier++){
      LabFusion_datatype * geotime=NULL;
      LabFusion_datatype * speedfunc=NULL;
      if(this->LNCC_status){
          speedfunc= (LabFusion_datatype *) calloc(this->numel, sizeof(LabFusion_datatype));
          for(int i=0;i<(this->CurrSizes->numel);i++){
              bool lnccexists=false;
              for(classifier_datatype lnccindex=0;lnccindex<this->Numb_Neigh;lnccindex++){
                  if(classifier==(this->LNCC[i+lnccindex*this->CurrSizes->numel])){
                      lnccexists=true;
                    }
                }
              speedfunc[i]=((LabFusion_datatype)(lnccexists)+0.1);
            }
          for(int currclass=0; currclass<this->NumberOfLabels;currclass++){
              if(this->verbose_level>0){
                  cout<<"Classifier "<<classifier+1<<" - Lable "<<currclass<<endl;
                  flush(cout);
                }

              for(int i=0;i<(this->CurrSizes->numel);i++){
                  CurrLableImage[i]=(inputCLASSIFIERptr[i+this->CurrSizes->numel*classifier]==currclass);
                }
              geotime=DoubleEuclideanDistance_3D(CurrLableImage,speedfunc,CurrSizes);
              for(int i=0;i<(this->CurrSizes->numel);i++){
                  this->W[i+currclass*this->numel]+=geotime[i];
                }
              free(geotime);
            }
          free(speedfunc);
        }
      else if(this->NCC_status){
          int NCCvalue=(int)(this->NCC[classifier]);
          if(NCCvalue>=0 && NCCvalue<this->CurrSizes->numclass){
              for(int currclass=0; currclass<this->NumberOfLabels;currclass++){
                  if(this->verbose_level>0){
                      cout<<"Classifier "<<classifier+1<<" - Lable "<<currclass<<endl;
                      flush(cout);
                    }
                  for(int i=0;i<(this->CurrSizes->numel);i++){
                      CurrLableImage[i]=(inputCLASSIFIERptr[i+this->CurrSizes->numel*classifier]==currclass);
                    }
                  geotime=DoubleEuclideanDistance_3D(CurrLableImage,NULL,CurrSizes);

                  for(int i=0;i<(this->CurrSizes->numel);i++){
                      this->W[i+currclass*this->numel]+=geotime[i];
                    }
                  free(geotime);
                }
            }

        }
      else{
          for( int currclass=0; currclass<this->NumberOfLabels; currclass++ ){
              if(this->verbose_level>0){
                  cout<<"Classifier "<<classifier+1<<" - Lable "<<currclass<<endl;
                  flush(cout);
                }
              for(int i=0;i<(this->CurrSizes->numel);i++){
                  CurrLableImage[i]=(inputCLASSIFIERptr[i+this->CurrSizes->numel*classifier]==currclass);
                }
              geotime=DoubleEuclideanDistance_3D(CurrLableImage,NULL,CurrSizes);

              for(int i=0;i<(this->CurrSizes->numel);i++){
                  this->W[i+currclass*this->numel]+=geotime[i];
                }
              free(geotime);
            }
        }
    }

  delete [] CurrLableImage;


  return 1;

}


int seg_LabFusion::MV_Estimate()
{

  LabFusion_datatype * tmpW=new LabFusion_datatype [this->NumberOfLabels];
  if(tmpW == NULL){
      fprintf(stderr,"* Error when alocating tmpW: Not enough memory\n");
      exit(1);
    }

  classifier_datatype * inputCLASSIFIERptr = static_cast<classifier_datatype *>(this->inputCLASSIFIER->data);

  this->loglik=0;
  for(int i=0;i<(this->CurrSizes->numel);i++){
      if(this->LNCC_status){
          for(int currclass=0; currclass<(this->NumberOfLabels);currclass++){
              tmpW[currclass]=0.0f;
            }
          for(int classifier=0; classifier<this->Numb_Neigh;classifier++){
              int LNCCvalue=(int)(this->LNCC[i+classifier*this->CurrSizes->numel]);
              if(LNCCvalue>=0 && LNCCvalue<this->CurrSizes->numclass){
                  tmpW[inputCLASSIFIERptr[i+LNCCvalue*this->CurrSizes->numel]]++;
                }
            }

          LabFusion_datatype tmpmaxval=-1;
          LabFusion_datatype tmpmaxindex=-1;

          if(this->NumberOfLabels>2){
              for(int currclass=0; currclass<(this->NumberOfLabels);currclass++){
                  if(tmpmaxval<tmpW[currclass]){
                      tmpmaxval=tmpW[currclass];
                      tmpmaxindex=currclass;
                    }
                }
              this->W[i]=tmpmaxindex;
            }
          else{
              this->W[i]=tmpW[1]/(float)this->Numb_Neigh;
            }

        }
      else if(this->NCC_status){
          for(int currclass=0; currclass<(this->NumberOfLabels);currclass++){
              tmpW[currclass]=0.0f;
            }
          for(int classifier=0; classifier<this->Numb_Neigh;classifier++){
              int NCCvalue=(int)(this->NCC[classifier]);
              if(NCCvalue>=0 && NCCvalue<this->CurrSizes->numclass){
                  tmpW[inputCLASSIFIERptr[i+NCCvalue*this->CurrSizes->numel]]++;
                }
            }
          LabFusion_datatype tmpmaxval=0;
          LabFusion_datatype tmpmaxindex=0;
          if(this->NumberOfLabels>2){
              for(int currclass=0; currclass<(this->NumberOfLabels);currclass++){
                  if(tmpmaxval<tmpW[currclass]){
                      tmpmaxval=tmpW[currclass];
                      tmpmaxindex=currclass;
                    }
                }
              this->W[i]=tmpmaxindex;
            }
          else{
              this->W[i]=tmpW[1]/(float)this->Numb_Neigh;
            }
        }
      else{
          for(int currclass=0; currclass<(this->NumberOfLabels);currclass++){
              tmpW[currclass]=0.0f;
            }
          for(int classifier=0; classifier<this->CurrSizes->numclass;classifier++){
              tmpW[inputCLASSIFIERptr[i+classifier*this->CurrSizes->numel]]++;

            }
          LabFusion_datatype tmpmaxval=-1;
          LabFusion_datatype tmpmaxindex=-1;
          if(this->NumberOfLabels>2){
              for(int currclass=0; currclass<(this->NumberOfLabels);currclass++){
                  if(tmpmaxval<tmpW[currclass]){
                      tmpmaxval=tmpW[currclass];
                      tmpmaxindex=currclass;
                    }
                }
              this->W[i]=tmpmaxindex;
            }
          else{
              this->W[i]=tmpW[1]/(float)this->CurrSizes->numclass;
            }
        }
    }
  delete [] tmpW;
  return 1;

}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::UpdateMRF()
{

  if(this->MRF_status){
      for(int i=0;i<(this->sizeAfterMaskingAndUncertainty*this->NumberOfLabels);i++){
          MRF[i]=0;
        }

      if(this->verbose_level>0){
          cout<<"Updating MRF"<<endl;
          flush(cout);
        }
      int col_size, plane_size,indexCentre, indexWest, indexEast, indexSouth, indexNorth, indexTop, indexBottom;
      int indexWestSmall, indexEastSmall, indexSouthSmall, indexNorthSmall, indexTopSmall, indexBottomSmall;
      int ix, iy, iz,maxiy, maxix, maxiz, neighbourclass;
      SegPrecisionTYPE Sum_Temp_MRF_Class_Expect;
      col_size = (int)(CurrSizes->xsize);
      plane_size = (int)(CurrSizes->xsize)*(CurrSizes->ysize);
      maxix = (int)(CurrSizes->xsize);
      maxiy = (int)(CurrSizes->ysize);
      maxiz = (int)(CurrSizes->zsize);
      SegPrecisionTYPE Clique[MaxMultiLableClass]={0};
      if(Clique == NULL){
          fprintf(stderr,"* The variable Clique was not allocated: OUT OF MEMORY!");
          exit(1);
        }

      SegPrecisionTYPE Temp_MRF_Class_Expect[MaxMultiLableClass]={0};
      if(Temp_MRF_Class_Expect == NULL){
          fprintf(stderr,"* The variable Temp_MRF_Class_Expect was not allocated: OUT OF MEMORY!");
          exit(1);
        }
      register int currclass;


      indexCentre=0;
      int index_within_mask=0;
      for (iz=0; iz<maxiz; iz++) {
          for (iy=0; iy<maxiy; iy++) {
              for (ix=0; ix<maxix; ix++) {
                  if(this->maskAndUncertainIndeces[indexCentre]>=0){

                      indexWest=(indexCentre-col_size)>=0?(indexCentre-col_size):-1;
                      indexEast=(indexCentre+col_size)<this->numel?(indexCentre+col_size):-1;
                      indexNorth=(indexCentre-1)>=0?(indexCentre-1):-1;
                      indexSouth=(indexCentre+1)<this->numel?(indexCentre+1):-1;
                      indexBottom=(indexCentre-plane_size)>=0?(indexCentre-plane_size):-1;
                      indexTop=(indexCentre+plane_size)<this->numel?(indexCentre+plane_size):-1;


                      indexWestSmall=this->maskAndUncertainIndeces[indexWest];
                      indexEastSmall=this->maskAndUncertainIndeces[indexEast];
                      indexNorthSmall=this->maskAndUncertainIndeces[indexNorth];
                      indexSouthSmall=this->maskAndUncertainIndeces[indexSouth];
                      indexBottomSmall=this->maskAndUncertainIndeces[indexBottom];
                      indexTopSmall=this->maskAndUncertainIndeces[indexTop];


                      for (currclass=0; currclass<this->NumberOfLabels; currclass++){

                          Clique[currclass] = W[index_within_mask];

                          Clique[currclass]+=((indexWestSmall>=0)?W[indexWestSmall]:(this->FinalSeg[indexWest]==currclass?1.0f:0.0f));
                          Clique[currclass]+=((indexEastSmall>=0)?W[indexEastSmall]:(this->FinalSeg[indexEast]==currclass?1.0f:0.0f));
                          Clique[currclass]+=((indexNorthSmall>=0)?W[indexNorthSmall]:(this->FinalSeg[indexNorth]==currclass?1.0f:0.0f));
                          Clique[currclass]+=((indexSouthSmall>=0)?W[indexSouthSmall]:(this->FinalSeg[indexSouth]==currclass?1.0f:0.0f));
                          Clique[currclass]+=((indexTopSmall>=0)?W[indexTopSmall]:(this->FinalSeg[indexTop]==currclass?1.0f:0.0f));
                          Clique[currclass]+=((indexBottomSmall>=0)?W[indexBottomSmall]:(this->FinalSeg[indexBottom]==currclass?1.0f:0.0f));

                          if(currclass<this->NumberOfLabels){
                              indexWestSmall+=(indexWestSmall>=0)?this->sizeAfterMaskingAndUncertainty:0;
                              indexEastSmall+=(indexEastSmall>=0)?this->sizeAfterMaskingAndUncertainty:0;
                              indexNorthSmall+=(indexNorthSmall>=0)?this->sizeAfterMaskingAndUncertainty:0;
                              indexSouthSmall+=(indexSouthSmall>=0)?this->sizeAfterMaskingAndUncertainty:0;
                              indexTopSmall+=(indexTopSmall>=0)?this->sizeAfterMaskingAndUncertainty:0;
                              indexBottomSmall+=(indexBottomSmall>=0)?this->sizeAfterMaskingAndUncertainty:0;
                            }
                        }
                      Sum_Temp_MRF_Class_Expect = 0;
                      for (currclass=0; currclass<this->NumberOfLabels; currclass++){
                          Temp_MRF_Class_Expect[currclass]=0;
                          for (neighbourclass=0; neighbourclass<this->NumberOfLabels; neighbourclass++){
                              Temp_MRF_Class_Expect[currclass]-=this->MRF_matrix[currclass+(this->NumberOfLabels)*neighbourclass]*Clique[neighbourclass];
                            }

                          Temp_MRF_Class_Expect[currclass] = exp(this->MRF_strength*Temp_MRF_Class_Expect[currclass]);
                          Sum_Temp_MRF_Class_Expect += Temp_MRF_Class_Expect[currclass];
                        }
                      if(Sum_Temp_MRF_Class_Expect<=0){
                          Sum_Temp_MRF_Class_Expect=0.0001;
                      }
                      for (currclass=0; currclass<this->NumberOfLabels; currclass++) {

                          MRF[index_within_mask+currclass*this->sizeAfterMaskingAndUncertainty]=(Temp_MRF_Class_Expect[currclass]/Sum_Temp_MRF_Class_Expect);

                        }
                      index_within_mask++;
                    }
                  indexCentre++;
                }
            }
        }
      if(this->verbose_level>0){
          cout<<"MRF updated"<<endl;
          flush(cout);
        }
    }


  return 1;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::EstimateInitialDensity()
{
  classifier_datatype * inputCLASSIFIERptr = static_cast<classifier_datatype *>(this->inputCLASSIFIER->data);


  int tempsum=0;
  for(int currclass=0;currclass<this->NumberOfLabels;currclass++){
      this->Prop[currclass]=0;
    }
  for( int i=0; i<this->numel; i++){
      if(this->maskAndUncertainIndeces[i]>=0){
          tempsum+=1;
          for(int classifier=0; classifier<this->numb_classif; classifier++){
              this->Prop[inputCLASSIFIERptr[i+this->numel*classifier]]+=1.0;
            }
        }
    }
  for(int currclass=0;currclass<this->NumberOfLabels;currclass++){
      this->Prop[currclass]=this->Prop[currclass]/(LabFusion_datatype)(tempsum*this->numb_classif);
    }


  if(this->verbose_level>0){
      cout << "Estimating initial proportion"<<endl;
      if(this->verbose_level>1){
          for(int currclass=0;currclass<this->NumberOfLabels;currclass++){
              cout<<"\tlabel["<<currclass<<"] - "<<this->Prop[currclass]<<endl;
            }
        }
      flush(cout);

    }
  return 1;

}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::SetPQ(LabFusion_datatype tmpP,LabFusion_datatype tmpQ)
{
  for(int i=0; i<this->numb_classif;i++){
      this->ConfusionMatrix[i]=(LabFusion_datatype)tmpP;
    }
  return 1;

}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::UpdateDensity()
{

  if(this->PropUpdate){

      LabFusion_datatype tempsum=0;
      for(int currclass=0;currclass<this->NumberOfLabels;currclass++){
          LabFusion_datatype tempprop=0;
          for( int i=0; i<this->numel; i++){
              if(this->maskAndUncertainIndeces[i]){
                  tempprop+=this->W[this->maskAndUncertainIndeces[i]+currclass*this->sizeAfterMaskingAndUncertainty];
                  tempsum+=this->W[this->maskAndUncertainIndeces[i]+currclass*this->sizeAfterMaskingAndUncertainty];
                }
            }
          this->Prop[currclass]=tempprop+0.0000001;
        }
      for(int currclass=0;currclass<this->NumberOfLabels;currclass++){
          this->Prop[currclass]=(this->Prop[currclass])/(LabFusion_datatype)(tempsum);
        }


      if(this->verbose_level>0){
          cout << "Estimating proportion"<<endl;
        }
      if(this->verbose_level>1){
          for(int currclass=0;currclass<this->NumberOfLabels;currclass++){
              if(this->verbose_level>1){
                  cout<<"\t"<<this->Prop[currclass]<<endl;
                }
            }
          flush(cout);

        }
    }
  return 1;

}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::Turn_Prop_Update_ON()
{
  this->PropUpdate=true;
  return 1;

}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::Allocate_Stuff_MV()
{


  this->W=new LabFusion_datatype [this->numel];
  if(W == NULL){
      fprintf(stderr,"* Error when alocating W: Not enough memory\n");
      exit(1);
    }

  for(int i=0;i<(this->numel);i++){
      this->W[i]=0;
    }

  return 0;

}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::Allocate_Stuff_SBA()
{


  this->W=new LabFusion_datatype [this->numel*this->NumberOfLabels];
  if(W == NULL){
      fprintf(stderr,"* Error when alocating W: Not enough memory\n");
      exit(1);
    }

  for(int i=0;i<(this->numel*this->NumberOfLabels);i++){
      this->W[i]=0;
    }

  return 0;

}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::Allocate_Stuff_STAPLE()
{
  if(this->verbose_level>1){
      cout<< "Allocating this->num_true";
      flush(cout);
    }



  this->FinalSeg = new LabFusion_datatype [this->numel];
  if(this->FinalSeg == NULL){
      fprintf(stderr,"\n* Error when alocating FinalSeg: Not enough memory\n");
      exit(1);
    }


  classifier_datatype * inputCLASSIFIERptr = static_cast<classifier_datatype *>(this->inputCLASSIFIER->data);

  if(this->uncertainflag){
      int * num_true=new int [this->NumberOfLabels+1];
      if(num_true == NULL){
          fprintf(stderr,"\n* Error when alocating num_true: Not enough memory\n");
          exit(1);
        }
      int current_size_after_masking_and_uncertainty=0;
      for(int i=0;i<(this->numel);i++){
          if(this->maskAndUncertainIndeces[i]>=0){
              for(int currClass=0; currClass<this->NumberOfLabels;currClass++){
                  num_true[currClass]=0;
                }
              for(int currClassifier=0; currClassifier<this->Numb_Neigh;currClassifier++){
                  if(this->LNCC_status){
                      num_true[(int)(inputCLASSIFIERptr[i+(int)LNCC[i+currClassifier*(this->numel)]*(this->numel)])]++;
                    }
                  else if(this->NCC_status){
                      num_true[(int)(inputCLASSIFIERptr[i+(int)NCC[currClassifier]*(this->numel)])]++;
                    }
                  else{
                      num_true[(int)(inputCLASSIFIERptr[i+(int)currClassifier*(this->numel)])]++;
                    }
                }
              bool is_uncertain=true;
              for(int currClass=0; currClass<this->NumberOfLabels;currClass++){
                  if(num_true[currClass]>=(this->Numb_Neigh)){
                      is_uncertain=false;
                    }
                }

              if(is_uncertain){
                  this->maskAndUncertainIndeces[i]=current_size_after_masking_and_uncertainty;
                  current_size_after_masking_and_uncertainty++;
                }
              else{
                  this->maskAndUncertainIndeces[i]=-1;
                }
            }
        }
      this->sizeAfterMaskingAndUncertainty=current_size_after_masking_and_uncertainty;
      delete [] num_true;
    }


  if(this->verbose_level>0 && this->uncertainflag){
      cout<<"Percentage of 'Active Area' = "<<(float)this->sizeAfterMaskingAndUncertainty/(float)this->numel*100.0f<<"%" <<endl;
    }

  if(this->verbose_level>1){
      cout<< "Allocating this->W ( "<<(float)(((float)this->sizeAfterMaskingAndUncertainty*(float)this->NumberOfLabels)/1024.0f/1024.0f/1024.0f) << " Gb )";
      flush(cout);
    }

  this->W=new LabFusion_datatype [this->sizeAfterMaskingAndUncertainty*this->NumberOfLabels];
  if(this->W == NULL){
      fprintf(stderr,"\n* Error when alocating this->W: Not enough memory\n");
      exit(1);
    }
  if(this->verbose_level>1){
      cout<< " - Done"<<endl;
      flush(cout);
    }

  if(this->verbose_level>1){
      cout<< "Initializing this->W";
      flush(cout);
    }
  for(int i=0;i<(this->sizeAfterMaskingAndUncertainty*this->NumberOfLabels);i++){
      this->W[i]=1.0f/(float)(this->NumberOfLabels);
    }
  if(this->verbose_level>1){
      cout<< " - Done"<<endl;
      flush(cout);
    }


  if(this->verbose_level>1){
      cout<< " - Done"<<endl;
      flush(cout);
    }


  if(this->verbose_level>1){
      cout<< "Calc Initial this->W";
      flush(cout);
    }


  int * num_true=new int [this->NumberOfLabels+1];
  if(num_true == NULL){
      fprintf(stderr,"\n* Error when alocating num_true: Not enough memory\n");
      exit(1);
    }
  for(int i=0;i<(this->numel);i++){

      // If it is uncertain area
      if(this->maskAndUncertainIndeces[i]>=0){
          for(int currLabel=0; currLabel<this->NumberOfLabels;currLabel++){
              num_true[currLabel]=0;
            }
          for(int currClassifier=0; currClassifier<this->Numb_Neigh;currClassifier++){
              if(this->LNCC_status){
                  num_true[(int)(inputCLASSIFIERptr[i+(int)LNCC[i+currClassifier*(this->numel)]*(this->numel)])]++;
                }
              else if(this->NCC_status){
                  num_true[(int)(inputCLASSIFIERptr[i+(int)NCC[currClassifier]*(this->numel)])]++;
                }
              else{
                  num_true[(int)(inputCLASSIFIERptr[i+(int)currClassifier*(this->numel)])]++;
                }
            }


          for(int currLabel=0; currLabel<this->NumberOfLabels;currLabel++){
              this->W[this->maskAndUncertainIndeces[i]+currLabel*(this->sizeAfterMaskingAndUncertainty)]=(float)(num_true[currLabel]/this->Numb_Neigh);
            }
        }
      // If it is not uncertain area
      else{
          for(int currLabel=0; currLabel<this->NumberOfLabels;currLabel++){
              num_true[currLabel]=0;
            }
          for(int currClassifier=0; currClassifier<this->Numb_Neigh;currClassifier++){
              if(this->LNCC_status){
                  num_true[(int)(inputCLASSIFIERptr[i+(int)LNCC[i+currClassifier*(this->numel)]*(this->numel)])]++;
                }
              else if(this->NCC_status){
                  num_true[(int)(inputCLASSIFIERptr[i+(int)NCC[currClassifier]*(this->numel)])]++;
                }
              else{
                  num_true[(int)(inputCLASSIFIERptr[i+(int)currClassifier*(this->numel)])]++;
                }
            }
          int maxClass=0;
          int maxCount=0;
          for(int currLabel=0; currLabel<this->NumberOfLabels;currLabel++){
              if(num_true[currLabel]>maxCount){
                  maxClass=currLabel;
                  maxCount=num_true[currLabel];
                }
            }
          this->FinalSeg[i]=(float)maxClass;
        }
    }
  delete [] num_true;

  if(this->verbose_level>1){
      cout<< " - Done"<<endl;
      flush(cout);
    }

  int *dim_array=new int[10];
  dim_array[0]=(int)inputCLASSIFIER->nx;
  dim_array[1]=(int)inputCLASSIFIER->ny;
  dim_array[2]=(int)inputCLASSIFIER->nz;

  /*
  if(this->verbose_level>1){
      cout<< "Dilating uncertainarea";
      flush(cout);
    }

  if(this->uncertainflag){
      if(dilunc>0){
          Dillate(this->uncertainarea,dilunc,dim_array,this->verbose_level);
        }
    }

  if(this->verbose_level>1){
      cout<< " - Done"<<endl;
      flush(cout);
    }
    */


  if(this->verbose_level>1){
      cout<< "Allocating MRF";
      flush(cout);
    }
  if(this->MRF_status){
      this->MRF=new LabFusion_datatype [this->sizeAfterMaskingAndUncertainty*this->NumberOfLabels];
      if(MRF == NULL){
          fprintf(stderr,"* The variable MRF was not allocated: OUT OF MEMORY!");
          exit(1);
        }
      for(int i=0; i<(this->sizeAfterMaskingAndUncertainty*this->NumberOfLabels); i++)
        this->MRF[i]=1.0f/this->NumberOfLabels;
    }
  if(this->verbose_level>1){
      cout<< " - Done"<<endl;
      flush(cout);
    }

  return 0;

}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

nifti_image * seg_LabFusion::GetResult_probability()
{
  nifti_image * Result = nifti_copy_nim_info(this->inputCLASSIFIER);

  Result->cal_min=10.0e30;
  Result->cal_max=-10.0e30;
  Result->datatype=DT_FLOAT;

  Result->cal_min=0;
  Result->cal_max=this->LableCorrespondences_small_to_big[this->NumberOfLabels];

  if(this->FilenameOut.empty()){
      this->FilenameOut.assign("LabFusion.nii.gz");
    }

  nifti_set_filenames(Result,(char*)this->FilenameOut.c_str(),0,0);
  nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
  Result->cal_max=(this->NumberOfLabels-1);

  if(TYPE_OF_FUSION==1 || TYPE_OF_FUSION==3){



      Result->dim[0]=4;
      Result->dim[4]=this->NumberOfLabels;
      nifti_update_dims_from_array(Result);
      Result->data = (void *)W;
      if(this->verbose_level>0){
          cout << "Saving Probabilistic Fused Label. One label per time point and all labels will be projected down in value."<<endl;
        }
    }
  else{
      Result->dim[0]=4;
      Result->dim[4]=1;
      nifti_update_dims_from_array(Result);

      if(this->NumberOfLabels>2){
          for(int i=0;i<(this->numel);i++){
              W[i]=this->LableCorrespondences_small_to_big[(int)W[i]];
            }
        }
      else{
          for(int i=0;i<(this->numel);i++){
              W[i]=this->LableCorrespondences_small_to_big[(int)(W[i]>0.5)];
            }
        }
      Result->data = (void *)W;
      if(this->verbose_level>0){
          cout << "Saving Integer Fused Label instead of Probabilistic. Probabilistic output is not available for MV"<<endl;
        }
    }

  this->W=NULL;

  return Result;
}

nifti_image * seg_LabFusion::GetResult_label()
{
  nifti_image * Result = nifti_copy_nim_info(this->inputCLASSIFIER);
  Result->cal_min=0;
  Result->cal_max=this->LableCorrespondences_small_to_big[this->NumberOfLabels];
  Result->datatype=DT_FLOAT;
  if(this->FilenameOut.empty()){
      this->FilenameOut.assign("LabFusion.nii.gz");
    }

  nifti_set_filenames(Result,(char*)this->FilenameOut.c_str(),0,0);
  nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
  Result->cal_max=(this->NumberOfLabels-1);


  if(this->verbose_level>0){
      cout << "Saving Integer Fused Label"<<endl;
    }

  if(TYPE_OF_FUSION==1){
      Result->data = this->FinalSeg;
      this->FinalSeg=NULL;
    }
  else{
      Result->data = (void *) calloc(Result->nvox, sizeof(float));
    }
  Result->dim[0]=4;
  Result->dim[4]=1;
  nifti_update_dims_from_array(Result);
  float * Resultdata = static_cast<float *>(Result->data);
  //cout << "TYPE_OF_FUSION = "<<TYPE_OF_FUSION<<endl;

  if(TYPE_OF_FUSION==3){
      for(int i=0;i<(this->numel);i++){
          LabFusion_datatype wmax=-1;
          int wmaxindex=0;
          for(int currclass=0; currclass<this->NumberOfLabels;currclass++){
              if(wmax<=W[i+currclass*this->numel]){
                  wmax=W[i+currclass*this->numel];
                  wmaxindex=currclass;
                }
            }
          Resultdata[i]=this->LableCorrespondences_small_to_big[wmaxindex];


        }
    }
  else if(TYPE_OF_FUSION==2){
      for(int i=0;i<(this->numel);i++){
          Resultdata[i]= this->LableCorrespondences_small_to_big[(int)W[i]];
        }

    }
  else if( TYPE_OF_FUSION==1){
      int uncertainindex=0;
      for(int i=0;i<(this->numel);i++){
          if(this->maskAndUncertainIndeces[i]>=0){
              LabFusion_datatype wmax=-1;
              int wmaxindex=0;
              for(int currclass=0; currclass<this->NumberOfLabels;currclass++){
                  if(wmax<=W[uncertainindex+currclass*this->sizeAfterMaskingAndUncertainty]){
                      wmax=W[uncertainindex+currclass*this->sizeAfterMaskingAndUncertainty];
                      wmaxindex=currclass;
                    }
                }
              Resultdata[i]=this->LableCorrespondences_small_to_big[wmaxindex];
              uncertainindex++;
            }
          else{
              Resultdata[i]=this->LableCorrespondences_small_to_big[(int)Resultdata[i]];
            }
        }
    }
  else{
      cout << "ERROR: Variable TYPE_OF_FUSION is not set. Cannot save."<<endl;
    }

  return Result;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
int  seg_LabFusion::Run_MV()
{

  TYPE_OF_FUSION=2;
  if((int)(this->verbose_level)>(int)(0)){
      cout << "Majority Vote: Verbose level " << this->verbose_level << endl;
      if(this->NCC_status ){
          cout << "Number of labels used = "<< this->Numb_Neigh<<endl;
        }
      else if(this->LNCC_status){
          cout << "Number of local labels used = "<< this->Numb_Neigh<<endl;
        }
      else{
          cout << "Using all classifiers"<<endl;
        }
    }


  Allocate_Stuff_MV();
  MV_Estimate();


  return 0;
}

int  seg_LabFusion::Run_SBA()
{
  TYPE_OF_FUSION=3;
  if((int)(this->verbose_level)>(int)(0)){
      cout << "SBA: Verbose level " << this->verbose_level << endl;
      if(this->NCC_status ){
          cout << "Number of labels used = "<< this->Numb_Neigh<<endl;
        }
      else if(this->LNCC_status ){
          cout << "Number of labels used = "<< this->Numb_Neigh<<endl;
        }
      else{
          cout << "Using all classifiers"<<endl;
        }
    }


  Allocate_Stuff_SBA();
  SBA_Estimate();
  //Find_WMax();
  return 0;
}

int  seg_LabFusion::Run_STAPLE_or_STEPS()
{

  this->TYPE_OF_FUSION=1;
  if((int)(this->verbose_level)>(int)(0)){
      cout << "STAPLE: Verbose level " << this->verbose_level << endl;
      if(this->NCC_status ){
          cout << "Number of labels used = "<< this->Numb_Neigh<<endl;
          cout << "MRF_beta = " << this->MRF_strength<<endl;
        }
      else if(this->LNCC_status){
          cout << "Number of local labels used = "<< this->Numb_Neigh<<endl;
          cout << "MRF_beta = " << this->MRF_strength<<endl;
        }
      else{
          cout << "Using all classifiers"<<endl;
          cout << "MRF_beta = " << this->MRF_strength<<endl;
        }
    }



  Allocate_Stuff_STAPLE();

  if(!(this->Fixed_Prop_status)){
      //UpdateDensity();
      EstimateInitialDensity();

    }
  //**************
  // EM Algorithm
  //**************



  if(this->verbose_level>0){
      cout << endl << "*******************************" << endl;
      cout << "Initialising " << endl;
    }
  this->iter=0;
  STAPLE_STEPS_Multiclass_Expectation_Maximization();
  if(this->verbose_level>0)printloglik(this->iter,this->loglik,this->oldloglik);


  this->iter=1;
  bool out= true;
  while (out) {

      if(this->verbose_level>0){
          cout << endl << "*******************************" << endl;
          cout << "Iteration " << iter << endl;
        }

      // Iterative Components - EM, MRF



      //MRF
      //out=false;
      UpdateMRF();
      //Update Density
      UpdateDensity();
      //Expectation
      //STAPLE_STEPS_Multiclass_Expectation();
      //Maximization
      //STAPLE_STEPS_Multiclass_Maximization();
      STAPLE_STEPS_Multiclass_Expectation_Maximization();
      // Print LogLik depending on the verbose level
      if(this->verbose_level>0)printloglik(this->iter,this->loglik,this->oldloglik);

      // EXIT CHECKS
      // Check convergence
      if( (((this->loglik-this->oldloglik)/this->oldloglik)<=this->Conv && this->iter>3) || iter>=this->maxIteration || isnan(this->loglik) )out=false;
      // Exit if this->Numb_Neigh==1
      if(this->Numb_Neigh==1)out=false;


      // Update LogLik
      this->oldloglik=this->loglik;
      iter++;
    }

  return 0;
}
