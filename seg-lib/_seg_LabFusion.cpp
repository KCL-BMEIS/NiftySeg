#include "_seg_LabFusion.h"




seg_LabFusion::seg_LabFusion(int _numb_classif, int numbclasses, int _Numb_Neigh, int _numb_modalities)
{
    TYPE_OF_FUSION=0;
    NumberOfLabels=numbclasses;
    inputCLASSIFIER = NULL; // pointer to external
    inputImage_status=0;
    verbose_level=0;
    loglik=0.01;
    oldloglik=0.001;

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
    this->Conv=0.01;
    this->numel=0;
    this->iter=0;
    this->CurrSizes=NULL;
    this->Numb_Neigh=_Numb_Neigh;

    this->Thresh_IMG_value=0;
    this->Thresh_IMG_DO=false;


    // SegParameters
    for (int i=0; i<(5000); i++)
    {
        LableCorrespondences_big_to_small[i]=0;
        LableCorrespondences_small_to_big[i]=0;
    }
    ConfusionMatrix=new segPrecisionTYPE [_numb_classif*numbclasses*numbclasses];
    if(ConfusionMatrix == NULL)
    {
        fprintf(stderr,"* Error when alocating ConfusionMatrix: Not enough memory\n");
        exit(1);
    }

    for (int i=0; i<(_numb_classif); i++)
    {
        for (int j=0; j<(numbclasses); j++)
        {
            for (int k=0; k<(numbclasses); k++)
            {
                if(k==j)
                {
                    ConfusionMatrix[k+j*numbclasses+i*numbclasses*numbclasses]=0.95;
                }
                else
                {
                    ConfusionMatrix[k+j*numbclasses+i*numbclasses*numbclasses]=0.05/(numbclasses-1);
                }
            }
        }
    }

    this->W=NULL;
    this->FinalSeg=NULL;
    this->maskAndUncertainIndeces=NULL;
    this->sizeAfterMaskingAndUncertainty=0;
    this->uncertainflag=false;
    this->uncertainthresh=0.999;
    this->dilunc=0;
    this->Prop=new segPrecisionTYPE [numbclasses];
    this->Fixed_Prop_status=false;
    this->PropUpdate=false;
    this->numb_classif=_numb_classif;
    this->numb_nummod=_numb_modalities;
    this->tracePQ=1;
    this->oldTracePQ=0;
    this->maxIteration=100;
    this->LNCC=NULL;
    this->NCC=NULL;
    this->Numb_Neigh=3;
    this->LNCC_status=false;
    this->NCC_status=false;


    // MRF Specific
    this->MRF_status=0;
    this->MRF_strength=0.0f;
    this->MRF_matrix=new segPrecisionTYPE [numbclasses*numbclasses];
    if(MRF_matrix == NULL)
    {
        fprintf(stderr,"* Error when alocating MRF_matrix: Not enough memory\n");
        exit(1);
    }

    for (int j=0; j<(numbclasses); j++)
    {
        for (int k=0; k<(numbclasses); k++)
        {
            if(k==j)
            {
                this->MRF_matrix[k+j*numbclasses]=0;
            }
            else
            {
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
    for(long i=0; i<(this->numel); i++)
    {
        if(inputMaskPtr[i]>0)
        {
            this->maskAndUncertainIndeces[i]=this->sizeAfterMaskingAndUncertainty;
            this->sizeAfterMaskingAndUncertainty++;
        }
        else
        {
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

    categoricalLabelType * CLASSIFIERptr = static_cast<categoricalLabelType *>(this->inputCLASSIFIER->data);
    int * NumberOfDifferentClassesHistogram=new int [5000];
    if(NumberOfDifferentClassesHistogram == NULL)
    {
        fprintf(stderr,"* Error when alocating NumberOfDifferentClassesHistogram: Not enough memory\n");
        exit(1);
    }

    for(long i=0; i<5000; i++)
    {
        NumberOfDifferentClassesHistogram[i]=0;
    }
    for(long i=0; i<(int)this->inputCLASSIFIER->nvox; i++)
    {
        NumberOfDifferentClassesHistogram[CLASSIFIERptr[i]]++;
    }

    int currlabindex=0;
    for(long i=0; i<5000; i++)
    {
        if(NumberOfDifferentClassesHistogram[i]>0)
        {
            LableCorrespondences_small_to_big[currlabindex]=i;
            LableCorrespondences_big_to_small[i]=currlabindex;
            currlabindex++;
        }
        else
        {
            LableCorrespondences_big_to_small[i]=-1;
        }
    }

    //this->NumberOfLabels=currlabindex;
    for(long i=0; i<(int)this->inputCLASSIFIER->nvox; i++)
    {
        CLASSIFIERptr[i]=LableCorrespondences_big_to_small[CLASSIFIERptr[i]];
    }


    this->maskAndUncertainIndeces=new long [this->numel];
    if(maskAndUncertainIndeces == NULL)
    {
        fprintf(stderr,"* Error when alocating mask_and_uncertain_indexes: Not enough memory\n");
        exit(1);
    }

    this->sizeAfterMaskingAndUncertainty=this->numel;
    for(long i=0; i<(this->numel); i++)
    {
        this->maskAndUncertainIndeces[i]=i;
    }
    delete [] NumberOfDifferentClassesHistogram;

    return 0;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::SetLNCC(nifti_image * _LNCC,nifti_image * BaseImage,segPrecisionTYPE distance,int Numb_Neigh)
{
    if((Numb_Neigh<_LNCC->nt) & (Numb_Neigh>0))
    {
        this->Numb_Neigh=(int)Numb_Neigh;
    }
    else
    {
        this->Numb_Neigh=_LNCC->nt;
    }

    if(this->nx!=_LNCC->nx || this->ny!=_LNCC->ny || this->nz!=_LNCC->nz)
    {
        fprintf(stderr,"* The image size of the images do not match");
        exit(1);
    }
    if(this->nx!=BaseImage->nx || this->ny!=BaseImage->ny || this->nz!=BaseImage->nz)
    {
        fprintf(stderr,"* The image size of the images do not match");
        exit(1);
    }

    if(_LNCC->nt==this->numb_classif)
    {
        if(_LNCC->datatype==DT_FLOAT)
        {
            this->LNCC=estimateLNCC5D(BaseImage,_LNCC,distance,Numb_Neigh,this->CurrSizes,this->verbose_level);
            this->LNCC_status = true;

        }
        else
        {
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
          for(long i=0; i<this->nx*this->ny*this->nz;i++){
              LabFusion_datatype * LNCCvalue_tmp = new LabFusion_datatype [_InMetric->nt];
              for(long currlable=0;currlable<_InMetric->nt; currlable++){
                  LNCCvalue_tmp[currlable]=InputMetricPtr[i+currlable*this->nx*this->ny*this->nz];
                }
              int * ordertmp=quickSort_order(&LNCCvalue_tmp[0],_InMetric->nt);
              for(long lable_order=0;lable_order<Numb_Neigh;lable_order++){
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
int seg_LabFusion::SetMLLNCC(nifti_image * _LNCC,nifti_image * BaseImage,segPrecisionTYPE distance,int levels,int Numb_Neigh)
{
    if((Numb_Neigh<_LNCC->nt) & (Numb_Neigh>0))
    {
        this->Numb_Neigh=(int)Numb_Neigh;
    }
    else
    {
        this->Numb_Neigh=_LNCC->nt;
    }

    if(this->nx!=_LNCC->nx || this->ny!=_LNCC->ny || this->nz!=_LNCC->nz)
    {
        fprintf(stderr,"* The image size of the images do not match");
        exit(1);
    }
    if(this->nx!=BaseImage->nx || this->ny!=BaseImage->ny || this->nz!=BaseImage->nz)
    {
        fprintf(stderr,"* The image size of the images do not match");
        exit(1);
    }

    if(_LNCC->nt==this->numb_classif)
    {
        if(_LNCC->datatype==DT_FLOAT)
        {
            this->LNCC=estimateMLNCC4D(BaseImage,_LNCC,distance,levels,Numb_Neigh,CurrSizes,this->verbose_level);

            this->LNCC_status = true;

        }
        else
        {
            cout << "LNCC is not LabFusion_datatype"<< endl;
        }
    }
    return 0;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
int seg_LabFusion::SetGNCC(nifti_image * _GNCC,nifti_image * BaseImage,int Numb_Neigh)
{
    if((Numb_Neigh<_GNCC->nt) & (Numb_Neigh>0))
    {
        this->Numb_Neigh=(int)Numb_Neigh;
    }
    else
    {
        this->Numb_Neigh=_GNCC->nt;
    }

    if(this->nx!=_GNCC->nx || this->ny!=_GNCC->ny || this->nz!=_GNCC->nz)
    {
        fprintf(stderr,"* The image size of the images do not match");
        exit(1);
    }
    if(this->nx!=BaseImage->nx || this->ny!=BaseImage->ny || this->nz!=BaseImage->nz)
    {
        fprintf(stderr,"* The image size of the images do not match");
        exit(1);
    }

    if(_GNCC->nt==this->numb_classif)
    {
        if(_GNCC->datatype==DT_FLOAT)
        {
            this->NCC=estimateNCC4D(BaseImage,_GNCC,Numb_Neigh,CurrSizes,this->verbose_level);
            this->NCC_status = true;
        }
        else
        {
            cout << "GNCC is not LabFusion_datatype"<< endl;
        }
    }

    return 0;
}

int seg_LabFusion::SetROINCC(nifti_image * _ROINCC,nifti_image * BaseImage,int Numb_Neigh,int DilSize)
{
    if((Numb_Neigh<_ROINCC->nt) & (Numb_Neigh>0))
    {
        this->Numb_Neigh=(int)Numb_Neigh;
    }
    else
    {
        this->Numb_Neigh=_ROINCC->nt;
    }


    if(this->nx!=_ROINCC->nx || this->ny!=_ROINCC->ny || this->nz!=_ROINCC->nz)
    {
        fprintf(stderr,"* The image size of the images do not match");
        exit(1);
    }
    if(this->nx!=BaseImage->nx || this->ny!=BaseImage->ny || this->nz!=BaseImage->nz)
    {
        fprintf(stderr,"* The image size of the images do not match");
        exit(1);
    }

    if(_ROINCC->nt==this->numb_classif)
    {
        if(_ROINCC->datatype==DT_FLOAT)
        {
            this->NCC=estimateROINCC4D(this->inputCLASSIFIER,BaseImage,_ROINCC,Numb_Neigh,CurrSizes,DilSize,this->verbose_level);
            this->NCC_status = true;
        }
        else
        {
            cout << "ROINCC is not LabFusion_datatype"<< endl;
        }
    }
    return 0;
}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::SetProp(segPrecisionTYPE r)
{
    this->Prop[0] = (segPrecisionTYPE)r;
    this->Fixed_Prop_status = true;
    return 0;
}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::SetConv(segPrecisionTYPE r)
{
    if(this->verbose_level)
    {
        cout<< "Convergence Ratio = " << r <<endl;
        flush(cout);
    }
    this->Conv = (segPrecisionTYPE)r;
    return 0;
}

int seg_LabFusion::SetImgThresh(segPrecisionTYPE _Thresh_IMG_value)
{
    this->Thresh_IMG_value=_Thresh_IMG_value;
    this->Thresh_IMG_DO=true;
    return 0;
}

int seg_LabFusion::SetUncThresh(float _uncthresh)
{
    this->uncertainthresh=(_uncthresh>=0.5)?((_uncthresh<=1)?_uncthresh:1):0.5;

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

    if(verblevel>2)
    {
        this->verbose_level=2;
    }
    else
    {
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
    if(numberiter<1)
    {
        this->maxIteration=1;
        cout << "Warning: It will only stop at iteration 1. For less than 1 iteration, use majority voting."<< endl;
    }
    else
    {
        this->maxIteration=numberiter;
    }
    return 0;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::Turn_MRF_ON(segPrecisionTYPE strength)
{
    this->MRF_status=true;
    this->MRF_strength=(segPrecisionTYPE)strength;
    return 0;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::Create_CurrSizes()
{
    this->CurrSizes = new ImageSize [1]();
    if(CurrSizes == NULL)
    {
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
    CurrSizes->nummod=this->numb_nummod;
    CurrSizes->numelmasked=0;
    CurrSizes->numelbias=0;
    return 0;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */



int seg_LabFusion::STAPLE_STEPS_Multiclass_Expectation_Maximization()
{

    segPrecisionTYPE *tmpW=new segPrecisionTYPE [this->NumberOfLabels];
    if(tmpW == NULL)
    {
        fprintf(stderr,"* Error when alocating tmpW: Not enough memory\n");
        exit(1);
    }

    categoricalLabelType * inputHumanRater = static_cast<categoricalLabelType *>(this->inputCLASSIFIER->data);

    bool * nccexists_once= new bool [this->CurrSizes->numclass];
    if(nccexists_once == NULL)
    {
        fprintf(stderr,"* Error when alocating nccexists_once: Not enough memory\n");
        exit(1);
    }

    for(long humanRater=0; humanRater<this->CurrSizes->numclass; humanRater++)nccexists_once[humanRater]=false;

    segPrecisionTYPE * ConfusionMatrix2=new segPrecisionTYPE [this->numb_classif*this->NumberOfLabels*this->NumberOfLabels];
    if(ConfusionMatrix2 == NULL)
    {
        fprintf(stderr,"* Error when alocating ConfusionMatrix2: Not enough memory\n");
        exit(1);
    }

    for(long i=0; i<(this->CurrSizes->numclass*this->NumberOfLabels*this->NumberOfLabels); i++)
    {
        ConfusionMatrix2[i]=0;
    }

    if(this->verbose_level>0)
    {
        cout << "Updating the Posteriors and the Performance Parameters"<< endl;
        flush(cout);
    }
    if(this->verbose_level>1)
    {
        cout << "[lab]  =  P \t\tQ"<< endl;
        flush(cout);
    }

    this->tracePQ=0;
    for(long humanRater=0; humanRater<this->CurrSizes->numclass; humanRater++)nccexists_once[humanRater]=false;

    for(long i=0; i<(this->numel); i++)
    {
        if(this->maskAndUncertainIndeces[i]>=0)
        {
            if(this->MRF_status)
            {
                // **************************
                //   Expectation with MRF
                // **************************
                // Set all tmpW to 0. If the label exits at least once, then set to 1.
                for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                {
                    tmpW[currlabelnumb]=0;
                }
                for(long humanRater=0; humanRater<this->Numb_Neigh; humanRater++)
                {
                    int LNCCvalue=(int)(this->LNCC[i+humanRater*this->CurrSizes->numel]);
                    tmpW[(int)inputHumanRater[i+LNCCvalue*this->CurrSizes->numel]]=1;
                }

                if(this->LNCC_status)
                {
                    for(long humanRater=0; humanRater<this->Numb_Neigh; humanRater++)
                    {
                        int LNCCvalue=(int)(this->LNCC[i+humanRater*this->CurrSizes->numel]);
                        if(LNCCvalue>=0 && LNCCvalue<this->CurrSizes->numclass)
                        {
                            for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                            {
                                tmpW[currlabelnumb]*=((segPrecisionTYPE)(this->ConfusionMatrix[(int)inputHumanRater[i+LNCCvalue*this->CurrSizes->numel]+currlabelnumb*this->NumberOfLabels+LNCCvalue*this->NumberOfLabels*this->NumberOfLabels]));
                            }
                        }
                    }
                }
                else if(this->NCC_status)
                {
                    for(long humanater=0; humanater<this->Numb_Neigh; humanater++)
                    {
                        int NCCvalue=(int)(this->NCC[humanater]);
                        if(NCCvalue>=0 && NCCvalue<this->CurrSizes->numclass)
                        {
                            for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                            {
                                tmpW[currlabelnumb]*=((segPrecisionTYPE)(this->ConfusionMatrix[(int)inputHumanRater[i+NCCvalue*this->CurrSizes->numel]+currlabelnumb*this->NumberOfLabels+NCCvalue*this->NumberOfLabels*this->NumberOfLabels]));
                            }
                        }
                    }
                }
                else
                {
                    for(long humanater=0; humanater<this->CurrSizes->numclass; humanater++)
                    {
                        for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                        {
                            tmpW[currlabelnumb]*=((segPrecisionTYPE)(this->ConfusionMatrix[(int)inputHumanRater[i+humanater*this->CurrSizes->numel]+currlabelnumb*this->NumberOfLabels+humanater*this->NumberOfLabels*this->NumberOfLabels]));
                        }
                    }
                }
                //cout << "1- "<<this->0] << "  "  <<this->Prop[1] <<endl;
                segPrecisionTYPE sumW=0;
                for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                {
                    tmpW[currlabelnumb]=((segPrecisionTYPE)this->Prop[currlabelnumb])*tmpW[currlabelnumb]*this->MRF[this->maskAndUncertainIndeces[i]+currlabelnumb*this->sizeAfterMaskingAndUncertainty];
                    sumW+=tmpW[currlabelnumb];
                }
                if(sumW<=0)
                {
                    if(verbose_level>1)
                    {
                        cout << "Normalize by zero at voxel - "<<i<<endl;
                    }
                }
                else
                {
                    for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                    {
                        W[this->maskAndUncertainIndeces[i]+currlabelnumb*this->sizeAfterMaskingAndUncertainty]=tmpW[currlabelnumb]/sumW;
                    }
                }
                for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                {
                    if(this->W[this->maskAndUncertainIndeces[i]+currlabelnumb*this->sizeAfterMaskingAndUncertainty]<0)
                    {
                        this->W[this->maskAndUncertainIndeces[i]+currlabelnumb*this->sizeAfterMaskingAndUncertainty]=0;
                    }
                    if(this->W[this->maskAndUncertainIndeces[i]+currlabelnumb*this->sizeAfterMaskingAndUncertainty]>1)
                    {
                        this->W[this->maskAndUncertainIndeces[i]+currlabelnumb*this->sizeAfterMaskingAndUncertainty]=1;
                    }
                }
            }

            else
            {
                // **************************
                //  Expectation without MRF
                // **************************
                for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                {
                    tmpW[currlabelnumb]=0;
                }
                for(long humanRater=0; humanRater<this->Numb_Neigh; humanRater++)
                {
                    int SimilarityIndex=this->LNCC==NULL?humanRater:(int)(this->LNCC[i+humanRater*this->CurrSizes->numel]);
                    tmpW[(int)inputHumanRater[i+SimilarityIndex*this->CurrSizes->numel]]=1;
                }

                if(this->LNCC_status)
                {

                    for(long humanater=0; humanater<this->Numb_Neigh; humanater++)
                    {

                        int LNCCvalue=(int)(this->LNCC[i+humanater*this->CurrSizes->numel]);
                        if(LNCCvalue>=0 && LNCCvalue<this->CurrSizes->numclass)
                        {
                            for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                            {

                                tmpW[currlabelnumb]*=((segPrecisionTYPE)(this->ConfusionMatrix[(int)inputHumanRater[i+LNCCvalue*this->CurrSizes->numel]+currlabelnumb*this->NumberOfLabels+LNCCvalue*this->NumberOfLabels*this->NumberOfLabels]));
                            }
                        }


                    }
                }
                else if(this->NCC_status)
                {
                    for(long humanater=0; humanater<this->Numb_Neigh; humanater++)
                    {
                        int NCCvalue=(int)(this->NCC[humanater]);
                        if(NCCvalue>=0 && NCCvalue<this->CurrSizes->numclass)
                        {
                            for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                            {
                                tmpW[currlabelnumb]*=((segPrecisionTYPE)(this->ConfusionMatrix[(int)inputHumanRater[i+NCCvalue*this->CurrSizes->numel]+
                                        currlabelnumb*this->NumberOfLabels+
                                        NCCvalue*this->NumberOfLabels*this->NumberOfLabels]));
                            }
                        }
                    }
                }
                else
                {
                    for(long humanater=0; humanater<this->CurrSizes->numclass; humanater++)
                    {
                        for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                        {
                            tmpW[currlabelnumb]*=((segPrecisionTYPE)(this->ConfusionMatrix[(int)inputHumanRater[i+humanater*this->CurrSizes->numel]+
                                    currlabelnumb*this->NumberOfLabels+
                                    humanater*this->NumberOfLabels*this->NumberOfLabels]));
                        }
                    }
                }
                segPrecisionTYPE sumW=0;
                for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                {
                    tmpW[currlabelnumb]=((segPrecisionTYPE)this->Prop[currlabelnumb])*tmpW[currlabelnumb];
                    sumW+=tmpW[currlabelnumb];
                }

                if(sumW<=0)
                {
                    if(verbose_level>1)cout << "Normalize by zero at voxel - "<<i<<endl;
                }
                else
                {

                    for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                    {
                        W[this->maskAndUncertainIndeces[i]+currlabelnumb*this->sizeAfterMaskingAndUncertainty]=tmpW[currlabelnumb]/sumW;

                    }

                }



                for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                {
                    if(this->W[this->maskAndUncertainIndeces[i]+currlabelnumb*this->sizeAfterMaskingAndUncertainty]<0)
                    {
                        this->W[this->maskAndUncertainIndeces[i]+currlabelnumb*this->sizeAfterMaskingAndUncertainty]=0;
                    }
                    if(this->W[this->maskAndUncertainIndeces[i]+currlabelnumb*this->sizeAfterMaskingAndUncertainty]>1)
                    {
                        this->W[this->maskAndUncertainIndeces[i]+currlabelnumb*this->sizeAfterMaskingAndUncertainty]=1;
                    }

                }

            }

            // **************************
            //       MAXIMIZATION
            // **************************

            if(this->LNCC_status)
            {

                for(long humanater=0; humanater<this->CurrSizes->numclass; humanater++)
                {
                    bool lnccexists=false;
                    for(categoricalLabelType lnccindex=0; lnccindex<this->Numb_Neigh; lnccindex++)
                    {
                        if(humanater==(this->LNCC[i+lnccindex*this->CurrSizes->numel]))
                        {
                            lnccexists=true;
                            nccexists_once[humanater]=true;
                            lnccindex=this->Numb_Neigh;
                        }
                    }

                    if(lnccexists)
                    {
                        for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                        {

                            ConfusionMatrix2[(int)inputHumanRater[i+this->CurrSizes->numel*humanater]
                                    +currlabelnumb*this->NumberOfLabels+
                                    humanater*this->NumberOfLabels*this->NumberOfLabels]+=
                                    W[this->maskAndUncertainIndeces[i]+currlabelnumb*this->sizeAfterMaskingAndUncertainty];
                        }
                    }
                }
            }
            else if(this->NCC_status)
            {
                bool nccexists=false;
                for(long humanater=0; humanater<this->CurrSizes->numclass; humanater++)
                {
                    for(long nccindex=0; nccindex<this->Numb_Neigh; nccindex++)
                    {
                        if(humanater==(this->NCC[nccindex]))
                        {
                            nccexists=true;
                            nccindex=this->Numb_Neigh;
                            nccexists_once[humanater]=true;
                        }
                    }

                    if(nccexists)
                    {
                        for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                        {

                            ConfusionMatrix2[(int)inputHumanRater[i+this->CurrSizes->numel*humanater]
                                    +currlabelnumb*this->NumberOfLabels+humanater*
                                    this->NumberOfLabels*this->NumberOfLabels]+=W[this->maskAndUncertainIndeces[i]+currlabelnumb*this->sizeAfterMaskingAndUncertainty];
                        }
                    }
                }
            }
            else
            {
                for(long humanater=0; humanater<this->CurrSizes->numclass; humanater++)
                {
                    for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                    {

                        ConfusionMatrix2[(long)inputHumanRater[i+this->CurrSizes->numel*humanater]
                                +currlabelnumb*this->NumberOfLabels+humanater*
                                this->NumberOfLabels*this->NumberOfLabels]+=W[this->maskAndUncertainIndeces[i]+currlabelnumb*this->sizeAfterMaskingAndUncertainty];
                    }
                }
            }
        }
    }

    // **************************
    // Normalize Confusion Matrix
    // **************************
    this->tracePQ=0;
    for(long humanater=0; humanater<this->CurrSizes->numclass; humanater++)
    {

        for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
        {
            double sumConf=0;
            for(long currlabelnumb2=0; currlabelnumb2<this->NumberOfLabels; currlabelnumb2++)
            {
                sumConf+=ConfusionMatrix2[currlabelnumb2+currlabelnumb*this->NumberOfLabels+humanater*
                        this->NumberOfLabels*this->NumberOfLabels]+0.0001f; // the +0.0001 is an EPS to avoid a binary case
            }

            if(sumConf<=0)
            {
                if(verbose_level>1)
                {
                    cout << "NOMALISE BY ZERO - "<<humanater<<" , "<< currlabelnumb<<endl;
                }
                for(long currlabelnumb2=0; currlabelnumb2<this->NumberOfLabels; currlabelnumb2++)
                {
                    ConfusionMatrix2[currlabelnumb2+currlabelnumb*this->NumberOfLabels+humanater*
                            this->NumberOfLabels*this->NumberOfLabels]=ConfusionMatrix[currlabelnumb2+currlabelnumb*this->NumberOfLabels+humanater*
                            this->NumberOfLabels*this->NumberOfLabels]+0.0001f; // same as above
                }
            }
            else
            {
                for(long currlabelnumb2=0; currlabelnumb2<this->NumberOfLabels; currlabelnumb2++)
                {
                    ConfusionMatrix2[currlabelnumb2+currlabelnumb*this->NumberOfLabels+humanater*
                            this->NumberOfLabels*this->NumberOfLabels]=(
                                ConfusionMatrix2[currlabelnumb2+currlabelnumb*this->NumberOfLabels+humanater*
                            this->NumberOfLabels*this->NumberOfLabels]+0.0001f)/sumConf; // again, same as above
                    this->tracePQ+=(currlabelnumb2==currlabelnumb)?ConfusionMatrix2[currlabelnumb2+currlabelnumb*this->NumberOfLabels+humanater*
                            this->NumberOfLabels*this->NumberOfLabels]:0;
                }
            }
        }


        if(this->verbose_level>1)
        {
            if(this->NCC_status)
            {
                if(nccexists_once[humanater])
                {
                    cout <<endl<< "["<<humanater+1<<"]=";
                    for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                    {
                        for(long currlabelnumb2=0; currlabelnumb2<this->NumberOfLabels; currlabelnumb2++)
                        {
                            cout <<"\t"<< setprecision(4)<<((ConfusionMatrix2[currlabelnumb2+currlabelnumb*this->NumberOfLabels+humanater*this->NumberOfLabels*this->NumberOfLabels]<0.0001)?0:ConfusionMatrix2[currlabelnumb2+currlabelnumb*this->NumberOfLabels+humanater*this->NumberOfLabels*this->NumberOfLabels]);
                        }
                        cout << endl;
                    }
                }
            }
            else if(this->LNCC_status)
            {
                if(nccexists_once[humanater])
                {
                    cout <<endl<< "["<<humanater+1<<"]=";
                    for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                    {
                        for(long currlabelnumb2=0; currlabelnumb2<this->NumberOfLabels; currlabelnumb2++)
                        {
                            cout <<"\t"<< setprecision(4)<<((ConfusionMatrix2[currlabelnumb2+currlabelnumb*this->NumberOfLabels+humanater*this->NumberOfLabels*this->NumberOfLabels]<0.0001)?0:ConfusionMatrix2[currlabelnumb2+currlabelnumb*this->NumberOfLabels+humanater*this->NumberOfLabels*this->NumberOfLabels]);
                        }
                        cout << endl;
                    }
                }
            }
            else
            {
                cout <<endl<< "[ "<<humanater+1<<" ]  =  ";
                cout << endl;
                for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                {
                    cout << "["<<humanater+1<<"]=";
                    for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                    {
                        for(long currlabelnumb2=0; currlabelnumb2<this->NumberOfLabels; currlabelnumb2++)
                        {
                            cout <<"\t"<< setprecision(4)<<((ConfusionMatrix2[currlabelnumb2+currlabelnumb*this->NumberOfLabels+humanater*this->NumberOfLabels*this->NumberOfLabels]<0.0001)?0:ConfusionMatrix2[currlabelnumb2+currlabelnumb*this->NumberOfLabels+humanater*this->NumberOfLabels*this->NumberOfLabels]);
                        }
                        cout << endl;
                    }
                }
            }
        }
    }

    for(long i=0; i<(this->numb_classif*this->NumberOfLabels*this->NumberOfLabels); i++)
    {
        ConfusionMatrix[i]=ConfusionMatrix2[i];
    }

    this->tracePQ=this->tracePQ/(this->numb_classif*this->NumberOfLabels);
    delete [] tmpW;
    delete [] ConfusionMatrix2;
    delete [] nccexists_once;
    return 1;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
int seg_LabFusion::SBA_Estimate()
{

    categoricalLabelType * inputCLASSIFIERptr = static_cast<categoricalLabelType *>(this->inputCLASSIFIER->data);
    //LabFusion_datatype eps=0.1;
    bool * CurrLableImage=new bool [this->numel];
    if(CurrLableImage == NULL)
    {
        fprintf(stderr,"* Error when alocating CurrLableImage: Not enough memory\n");
        exit(1);
    }
    if(this->verbose_level>0)
    {
        cout<<"Calculating Euclidean Distances"<<endl;
        flush(cout);
    }


    for(long classifier=0; classifier<this->numb_classif; classifier++)
    {
        segPrecisionTYPE * geotime=NULL;
        segPrecisionTYPE * speedfunc=NULL;
        if(this->LNCC_status)
        {
            speedfunc= (segPrecisionTYPE *) calloc(this->numel, sizeof(segPrecisionTYPE));
            for(long i=0; i<(this->CurrSizes->numel); i++)
            {
                bool lnccexists=false;
                for(long lnccindex=0; lnccindex<this->Numb_Neigh; lnccindex++)
                {
                    if(classifier==(this->LNCC[i+lnccindex*this->CurrSizes->numel]))
                    {
                        lnccexists=true;
                    }
                }
                speedfunc[i]=((segPrecisionTYPE)(lnccexists)+0.1);
            }
            for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
            {
                if(this->verbose_level>0)
                {
                    cout<<"Classifier "<<classifier+1<<" - Lable "<<currlabelnumb<<endl;
                    flush(cout);
                }

                for(long i=0; i<(this->CurrSizes->numel); i++)
                {
                    CurrLableImage[i]=(inputCLASSIFIERptr[i+this->CurrSizes->numel*classifier]==currlabelnumb);
                }
                geotime=DoubleEuclideanDistance_3D(CurrLableImage,speedfunc,CurrSizes);
                for(long i=0; i<(this->CurrSizes->numel); i++)
                {
                    this->W[i+currlabelnumb*this->numel]+=geotime[i]*fabs(geotime[i]);
                }
                free(geotime);
            }
            free(speedfunc);
        }
        else if(this->NCC_status)
        {
            int NCCvalue=(int)(this->NCC[classifier]);
            if(NCCvalue>=0 && NCCvalue<this->CurrSizes->numclass)
            {
                for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                {
                    if(this->verbose_level>0)
                    {
                        cout<<"Classifier "<<classifier+1<<" - Lable "<<currlabelnumb<<endl;
                        flush(cout);
                    }
                    for(long i=0; i<(this->CurrSizes->numel); i++)
                    {
                        CurrLableImage[i]=(inputCLASSIFIERptr[i+this->CurrSizes->numel*classifier]==currlabelnumb);
                    }
                    geotime=DoubleEuclideanDistance_3D(CurrLableImage,NULL,CurrSizes);

                    for(long i=0; i<(this->CurrSizes->numel); i++)
                    {
                        this->W[i+currlabelnumb*this->numel]+=geotime[i]*fabs(geotime[i]);
                    }
                    free(geotime);
                }
            }

        }
        else
        {
            for( long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++ )
            {
                if(this->verbose_level>0)
                {
                    cout<<"Classifier "<<classifier+1<<" - Lable "<<currlabelnumb<<endl;
                    flush(cout);
                }
                for(long i=0; i<(this->CurrSizes->numel); i++)
                {
                    CurrLableImage[i]=(inputCLASSIFIERptr[i+this->CurrSizes->numel*classifier]==currlabelnumb);
                }
                geotime=DoubleEuclideanDistance_3D(CurrLableImage,NULL,CurrSizes);

                for(long i=0; i<(this->CurrSizes->numel); i++)
                {
                    this->W[i+currlabelnumb*this->numel]+=geotime[i]*fabs(geotime[i]);
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

    segPrecisionTYPE * tmpW=new segPrecisionTYPE [this->NumberOfLabels];
    if(tmpW == NULL)
    {
        fprintf(stderr,"* Error when alocating tmpW: Not enough memory\n");
        exit(1);
    }

    categoricalLabelType * inputCLASSIFIERptr = static_cast<categoricalLabelType *>(this->inputCLASSIFIER->data);

    this->tracePQ=0;
    for(long i=0; i<(this->CurrSizes->numel); i++)
    {
        if(this->LNCC_status)
        {
            for(long currlabelnumb=0; currlabelnumb<(this->NumberOfLabels); currlabelnumb++)
            {
                tmpW[currlabelnumb]=0.0f;
            }
            for(long classifier=0; classifier<this->Numb_Neigh; classifier++)
            {
                int LNCCvalue=(int)(this->LNCC[i+classifier*this->CurrSizes->numel]);
                if(LNCCvalue>=0 && LNCCvalue<this->CurrSizes->numclass)
                {
                    tmpW[inputCLASSIFIERptr[i+LNCCvalue*this->CurrSizes->numel]]++;
                }
            }

            segPrecisionTYPE tmpmaxval=-1;
            segPrecisionTYPE tmpmaxindex=-1;

            if(this->NumberOfLabels>2)
            {
                for(long currlabelnumb=0; currlabelnumb<(this->NumberOfLabels); currlabelnumb++)
                {
                    if(tmpmaxval<tmpW[currlabelnumb])
                    {
                        tmpmaxval=tmpW[currlabelnumb];
                        tmpmaxindex=currlabelnumb;
                    }
                }
                this->W[i]=tmpmaxindex;
            }
            else
            {
                this->W[i]=tmpW[1]/(float)this->Numb_Neigh;
            }

        }
        else if(this->NCC_status)
        {
            for(long currlabelnumb=0; currlabelnumb<(this->NumberOfLabels); currlabelnumb++)
            {
                tmpW[currlabelnumb]=0.0f;
            }
            for(long classifier=0; classifier<this->Numb_Neigh; classifier++)
            {
                int NCCvalue=(int)(this->NCC[classifier]);
                if(NCCvalue>=0 && NCCvalue<this->CurrSizes->numclass)
                {
                    tmpW[inputCLASSIFIERptr[i+NCCvalue*this->CurrSizes->numel]]++;
                }
            }
            segPrecisionTYPE tmpmaxval=0;
            segPrecisionTYPE tmpmaxindex=0;
            if(this->NumberOfLabels>2)
            {
                for(long currlabelnumb=0; currlabelnumb<(this->NumberOfLabels); currlabelnumb++)
                {
                    if(tmpmaxval<tmpW[currlabelnumb])
                    {
                        tmpmaxval=tmpW[currlabelnumb];
                        tmpmaxindex=currlabelnumb;
                    }
                }
                this->W[i]=tmpmaxindex;
            }
            else
            {
                this->W[i]=tmpW[1]/(float)this->Numb_Neigh;
            }
        }
        else
        {
            for(long currlabelnumb=0; currlabelnumb<(this->NumberOfLabels); currlabelnumb++)
            {
                tmpW[currlabelnumb]=0.0f;
            }
            for(long classifier=0; classifier<this->CurrSizes->numclass; classifier++)
            {
                tmpW[inputCLASSIFIERptr[i+classifier*this->CurrSizes->numel]]++;

            }
            segPrecisionTYPE tmpmaxval=-1;
            segPrecisionTYPE tmpmaxindex=-1;
            if(this->NumberOfLabels>2)
            {
                for(long currlabelnumb=0; currlabelnumb<(this->NumberOfLabels); currlabelnumb++)
                {
                    if(tmpmaxval<tmpW[currlabelnumb])
                    {
                        tmpmaxval=tmpW[currlabelnumb];
                        tmpmaxindex=currlabelnumb;
                    }
                }
                this->W[i]=tmpmaxindex;
            }
            else
            {
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

    if(this->MRF_status)
    {
        for(long i=0; i<(this->sizeAfterMaskingAndUncertainty*this->NumberOfLabels); i++)
        {
            MRF[i]=0;
        }

        if(this->verbose_level>0)
        {
            cout<<"Updating MRF"<<endl;
            flush(cout);
        }
        int col_size, plane_size,indexCentre, indexWest, indexEast, indexSouth, indexNorth, indexTop, indexBottom;
        int indexWestSmall, indexEastSmall, indexSouthSmall, indexNorthSmall, indexTopSmall, indexBottomSmall;
        int ix, iy, iz, maxiy, maxix, maxiz, neighbourclass;
        segPrecisionTYPE Sum_Temp_MRF_Class_Expect;
        col_size = (int)(CurrSizes->xsize);
        plane_size = (int)(CurrSizes->xsize)*(CurrSizes->ysize);
        maxix = (int)(CurrSizes->xsize);
        maxiy = (int)(CurrSizes->ysize);
        maxiz = (int)(CurrSizes->zsize);
        segPrecisionTYPE Clique[maxMultiLableClass]= {0};
        segPrecisionTYPE Temp_MRF_Class_Expect[maxMultiLableClass]= {0};

        register int currlabelnumb;
        indexCentre=0;
        int index_within_mask=0;
        for (iz=0; iz<maxiz; iz++)
        {
            for (iy=0; iy<maxiy; iy++)
            {
                for (ix=0; ix<maxix; ix++)
                {
                    if(this->maskAndUncertainIndeces[indexCentre]>=0)
                    {

                        indexWest=(indexCentre-col_size)>=0?(indexCentre-col_size):-1;
                        indexEast=(indexCentre+col_size)<this->numel?(indexCentre+col_size):-1;
                        indexNorth=(indexCentre-1)>=0?(indexCentre-1):-1;
                        indexSouth=(indexCentre+1)<this->numel?(indexCentre+1):-1;
                        indexBottom=(indexCentre-plane_size)>=0?(indexCentre-plane_size):-1;
                        indexTop=(indexCentre+plane_size)<this->numel?(indexCentre+plane_size):-1;


                        indexWestSmall=indexWest>=0?this->maskAndUncertainIndeces[indexWest]:-1;
                        indexEastSmall=indexEast>=0?this->maskAndUncertainIndeces[indexEast]:-1;
                        indexNorthSmall=indexNorth>=0?this->maskAndUncertainIndeces[indexNorth]:-1;
                        indexSouthSmall=indexSouth>=0?this->maskAndUncertainIndeces[indexSouth]:-1;
                        indexBottomSmall=indexBottom>=0?this->maskAndUncertainIndeces[indexBottom]:-1;
                        indexTopSmall=indexTop>=0?this->maskAndUncertainIndeces[indexTop]:-1;


                        for (currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                        {

                            Clique[currlabelnumb] = W[index_within_mask];

                            Clique[currlabelnumb]+=((indexWestSmall>=0)?W[indexWestSmall]:0.0f);
                            Clique[currlabelnumb]+=((indexEastSmall>=0)?W[indexEastSmall]:0.0f);
                            Clique[currlabelnumb]+=((indexNorthSmall>=0)?W[indexNorthSmall]:0.0f);
                            Clique[currlabelnumb]+=((indexSouthSmall>=0)?W[indexSouthSmall]:0.0f);
                            Clique[currlabelnumb]+=((indexTopSmall>=0)?W[indexTopSmall]:0.0f);
                            Clique[currlabelnumb]+=((indexBottomSmall>=0)?W[indexBottomSmall]:0.0f);

                            if(currlabelnumb<this->NumberOfLabels)
                            {
                                indexWestSmall+=(indexWestSmall>=0)?this->sizeAfterMaskingAndUncertainty:0;
                                indexEastSmall+=(indexEastSmall>=0)?this->sizeAfterMaskingAndUncertainty:0;
                                indexNorthSmall+=(indexNorthSmall>=0)?this->sizeAfterMaskingAndUncertainty:0;
                                indexSouthSmall+=(indexSouthSmall>=0)?this->sizeAfterMaskingAndUncertainty:0;
                                indexTopSmall+=(indexTopSmall>=0)?this->sizeAfterMaskingAndUncertainty:0;
                                indexBottomSmall+=(indexBottomSmall>=0)?this->sizeAfterMaskingAndUncertainty:0;
                            }
                        }
                        Sum_Temp_MRF_Class_Expect = 0;
                        for (currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                        {
                            //cout <<  Clique[currlabelnumb] <<" ";
                            Temp_MRF_Class_Expect[currlabelnumb]=0;
                            for (neighbourclass=0; neighbourclass<this->NumberOfLabels; neighbourclass++)
                            {
                                Temp_MRF_Class_Expect[currlabelnumb]-=this->MRF_matrix[currlabelnumb+(this->NumberOfLabels)*neighbourclass]*Clique[neighbourclass];
                            }

                            Temp_MRF_Class_Expect[currlabelnumb] = exp(this->MRF_strength*Temp_MRF_Class_Expect[currlabelnumb]);
                            Sum_Temp_MRF_Class_Expect += Temp_MRF_Class_Expect[currlabelnumb];
                        }
                        //cout <<endl;

                        if(Sum_Temp_MRF_Class_Expect<=0)
                        {
                            Sum_Temp_MRF_Class_Expect=0.0001;
                        }
                        //cout<<Temp_MRF_Class_Expect[0]/Sum_Temp_MRF_Class_Expect<<" "<<Temp_MRF_Class_Expect[1]/Sum_Temp_MRF_Class_Expect<<endl;
                        for (currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                        {

                            MRF[index_within_mask+currlabelnumb*this->sizeAfterMaskingAndUncertainty]=(Temp_MRF_Class_Expect[currlabelnumb]/Sum_Temp_MRF_Class_Expect);

                        }
                        index_within_mask++;
                    }
                    indexCentre++;
                }
            }
        }
        if(this->verbose_level>0)
        {
            cout<<"MRF updated"<<endl;
            flush(cout);
        }
    }


    return 1;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::EstimateInitialDensity()
{
//    classifier_datatype * inputCLASSIFIERptr = static_cast<classifier_datatype *>(this->inputCLASSIFIER->data);
//    int tempsum=0;
//    int * propcount = new int [this->NumberOfLabels];
//    for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
//    {
//        propcount[currlabelnumb]=0;
//        this->Prop[currlabelnumb]=0;
//    }
//    for( int i=0; i<this->numel; i++)
//    {
//        if(this->maskAndUncertainIndeces[i]>=0)
//        {
//            tempsum+=1;
//            for(long classifier=0; classifier<this->numb_classif; classifier++)
//            {
//                propcount[inputCLASSIFIERptr[i+this->numel*classifier]]++;
//
//                //cout << (float)inputCLASSIFIERptr[i+this->numel*classifier];
//            }
//            //cout << "\n";
//        }
//   }
    for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
    {
        //cout<<"\tlabel["<<currlabelnumb<<"] - "<<propcount[currlabelnumb]<<" - "<<tempsum<<endl;
        //this->Prop[currlabelnumb]=propcount[currlabelnumb]/(LabFusion_datatype)((float)tempsum*(float)this->numb_classif);
        this->Prop[currlabelnumb]=1.0f/this->NumberOfLabels;
    }

    // Adding the 0.0001 trace amount and renormalise works like an EPS for numerical precision.
//    float tempsum2=0.0f;
//    for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
//    {
//        this->Prop[currlabelnumb]=(this->Prop[currlabelnumb]+0.0001f);
//        tempsum2+=this->Prop[currlabelnumb];
//    }
//    for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
//    {
//        //cout<<"\tlabel["<<currlabelnumb<<"] - "<<this->Prop[currlabelnumb]<<endl;
//        this->Prop[currlabelnumb]=(this->Prop[currlabelnumb])/(LabFusion_datatype)(tempsum2);;
//    }
//
//    delete [] propcount;
//    if(this->verbose_level>0)
//    {
//        cout << "Estimating initial proportion"<<endl;
//        if(this->verbose_level>0)
//        {
//            for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
//            {
//                cout<<"\tlabel["<<currlabelnumb<<"] - "<<this->Prop[currlabelnumb]<<endl;
//            }
//        }
//        flush(cout);
//    }

    return 1;

}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::SetPQ(segPrecisionTYPE tmpP,segPrecisionTYPE tmpQ)
{
    for(long i=0; i<this->numb_classif; i++)
    {
        this->ConfusionMatrix[i]=(segPrecisionTYPE)tmpP;
    }
    return 1;

}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::UpdateDensity()
{

    if(this->PropUpdate)
    {

        segPrecisionTYPE tempsum=0;
        for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
        {
            segPrecisionTYPE tempprop=0;
            for( int i=0; i<this->numel; i++)
            {
                if(this->maskAndUncertainIndeces[i]>0)
                {
                    tempprop+=this->W[this->maskAndUncertainIndeces[i]+currlabelnumb*this->sizeAfterMaskingAndUncertainty];
                    tempsum+=this->W[this->maskAndUncertainIndeces[i]+currlabelnumb*this->sizeAfterMaskingAndUncertainty];
                }
            }
            this->Prop[currlabelnumb]=tempprop;
        }
        for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
        {
            this->Prop[currlabelnumb]=(this->Prop[currlabelnumb])/(segPrecisionTYPE)(tempsum);
        }

        // Adding the 0.0001 trace amount and renormalise works like an EPS for numerical precision.
        tempsum=0;
        for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
        {
            this->Prop[currlabelnumb]=(this->Prop[currlabelnumb]+0.0001f);
            tempsum+=this->Prop[currlabelnumb];
        }
        for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
        {
            this->Prop[currlabelnumb]=(this->Prop[currlabelnumb])/(segPrecisionTYPE)(tempsum);;
        }


        if(this->verbose_level>0)
        {
            cout << "Estimating proportion"<<endl;
        }
        if(this->verbose_level>1)
        {
            for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
            {
                if(this->verbose_level>1)
                {
                    cout<<"\t"<<this->Prop[currlabelnumb]<<endl;
                }
            }
            flush(cout);

        }
    }
    return 1;

}

int seg_LabFusion::UpdateDensity_noTest()
{



    segPrecisionTYPE tempsum=0;
    for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
    {
        segPrecisionTYPE tempprop=0;
        for( int i=0; i<this->numel; i++)
        {
            if(this->maskAndUncertainIndeces[i])
            {
                tempprop+=this->W[this->maskAndUncertainIndeces[i]+currlabelnumb*this->sizeAfterMaskingAndUncertainty];
                tempsum+=this->W[this->maskAndUncertainIndeces[i]+currlabelnumb*this->sizeAfterMaskingAndUncertainty];
            }
        }
        this->Prop[currlabelnumb]=tempprop+0.0001f;
    }
    for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
    {
        this->Prop[currlabelnumb]=(this->Prop[currlabelnumb])/(segPrecisionTYPE)(tempsum);
    }


    if(this->verbose_level>0)
    {
        cout << "Estimating proportion"<<endl;
    }
    if(this->verbose_level>1)
    {
        for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
        {
            if(this->verbose_level>1)
            {
                cout<<"\t"<<this->Prop[currlabelnumb]<<endl;
            }
        }
        flush(cout);

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


    this->W=new segPrecisionTYPE [this->numel];
    if(W == NULL)
    {
        fprintf(stderr,"* Error when alocating W: Not enough memory\n");
        exit(1);
    }

    for(long i=0; i<(this->numel); i++)
    {
        this->W[i]=0;
    }

    return 0;

}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::Allocate_Stuff_SBA()
{


    this->W=new segPrecisionTYPE [this->numel*this->NumberOfLabels];
    if(W == NULL)
    {
        fprintf(stderr,"* Error when alocating W: Not enough memory\n");
        exit(1);
    }

    for(long i=0; i<(this->numel*this->NumberOfLabels); i++)
    {
        this->W[i]=0;
    }

    return 0;

}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::Allocate_Stuff_STAPLE()
{
    if(this->verbose_level>1)
    {
        cout<< "Allocating this->FinalSeg";
        flush(cout);
    }
    this->FinalSeg = new segPrecisionTYPE [this->numel];
    if(this->verbose_level>1)
    {
        cout<< " -> Done"<<endl;
        flush(cout);
    }
    if(this->FinalSeg == NULL)
    {
        fprintf(stderr,"\n* Error when alocating FinalSeg: Not enough memory\n");
        exit(1);
    }

    categoricalLabelType * inputCLASSIFIERptr = static_cast<categoricalLabelType *>(this->inputCLASSIFIER->data);

    if(this->uncertainflag)
    {
        int * num_true=new int [this->NumberOfLabels+1];
        if(num_true == NULL)
        {
            fprintf(stderr,"\n* Error when alocating num_true: Not enough memory\n");
            exit(1);
        }
        long current_size_after_masking_and_uncertainty=0;
        for(long i=0; i<(this->numel); i++)
        {
            if(this->maskAndUncertainIndeces[i]>=0)
            {
                for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                {
                    num_true[currlabelnumb]=0;
                }
                for(long currlabelnumbifier=0; currlabelnumbifier<this->Numb_Neigh; currlabelnumbifier++)
                {
                    if(this->LNCC_status)
                    {
                        num_true[(int)(inputCLASSIFIERptr[i+(int)LNCC[i+currlabelnumbifier*(this->numel)]*(this->numel)])]++;
                    }
                    else if(this->NCC_status)
                    {
                        num_true[(int)(inputCLASSIFIERptr[i+(int)NCC[currlabelnumbifier]*(this->numel)])]++;
                    }
                    else
                    {
                        num_true[(int)(inputCLASSIFIERptr[i+(int)currlabelnumbifier*(this->numel)])]++;
                    }
                }
                bool is_uncertain=true;
                for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                {
                    float thisthresh=ceil(this->uncertainthresh*(double)this->Numb_Neigh);
                    //float thisthresh=7;
                    if(num_true[currlabelnumb]>=(thisthresh))
                    {
                        is_uncertain=false;
                    }
                }

                if(is_uncertain)
                {
                    this->maskAndUncertainIndeces[i]=current_size_after_masking_and_uncertainty;
                    current_size_after_masking_and_uncertainty++;
                }
                else
                {
                    this->maskAndUncertainIndeces[i]=-1;
                }
            }
        }
        this->sizeAfterMaskingAndUncertainty=current_size_after_masking_and_uncertainty;
        delete [] num_true;
    }


    if(this->verbose_level>0 && this->uncertainflag)
    {
        cout<<"Percentage of 'Active Area' = "<<(float)this->sizeAfterMaskingAndUncertainty/(float)this->numel*100.0f<<"%" <<endl;
    }

    if(this->verbose_level>1)
    {
        cout<< "Allocating this->W ( "<<(float)(((float)this->sizeAfterMaskingAndUncertainty*(float)this->NumberOfLabels)/1024.0f/1024.0f/1024.0f) << " Gb )";
        flush(cout);
    }
    this->W=new segPrecisionTYPE [this->sizeAfterMaskingAndUncertainty*this->NumberOfLabels];
    if(this->W == NULL)
    {
        fprintf(stderr,"\n* Error when alocating this->W: Not enough memory\n");
        exit(1);
    }
    if(this->verbose_level>1)
    {
        cout<< " - Done"<<endl;
        flush(cout);
    }

    if(this->verbose_level>1)
    {
        cout<< "Initializing this->W";
        flush(cout);
    }
    for(long i=0; i<(this->sizeAfterMaskingAndUncertainty*this->NumberOfLabels); i++)
    {
        this->W[i]=1.0f/(float)(this->NumberOfLabels);
    }
    if(this->verbose_level>1)
    {
        cout<< " - Done"<<endl;
        flush(cout);
    }

    if(this->verbose_level>1)
    {
        cout<< "Calc Initial this->W";
        flush(cout);
    }

    int * num_true=new int [this->NumberOfLabels+1];
    if(num_true == NULL)
    {
        fprintf(stderr,"\n* Error when alocating num_true: Not enough memory\n");
        exit(1);
    }
    for(long i=0; i<(this->numel); i++)
    {

        // If it is uncertain area
        if(this->maskAndUncertainIndeces[i]>=0)
        {
            for(long currLabel=0; currLabel<this->NumberOfLabels; currLabel++)
            {
                num_true[currLabel]=0;
            }
            for(long currlabelnumbifier=0; currlabelnumbifier<this->Numb_Neigh; currlabelnumbifier++)
            {
                if(this->LNCC_status)
                {
                    num_true[(int)(inputCLASSIFIERptr[i+(int)LNCC[i+currlabelnumbifier*(this->numel)]*(this->numel)])]++;
                }
                else if(this->NCC_status)
                {
                    num_true[(int)(inputCLASSIFIERptr[i+(int)NCC[currlabelnumbifier]*(this->numel)])]++;
                }
                else
                {
                    num_true[(int)(inputCLASSIFIERptr[i+(int)currlabelnumbifier*(this->numel)])]++;
                }
            }


            for(long currLabel=0; currLabel<this->NumberOfLabels; currLabel++)
            {
                this->W[this->maskAndUncertainIndeces[i]+currLabel*(this->sizeAfterMaskingAndUncertainty)]=(float)(num_true[currLabel]/this->Numb_Neigh);
            }
            this->FinalSeg[i]=0;
        }
        // If it is not uncertain area
        else
        {
            for(long currLabel=0; currLabel<this->NumberOfLabels; currLabel++)
            {
                num_true[currLabel]=0;
            }
            for(long currlabelnumbifier=0; currlabelnumbifier<this->Numb_Neigh; currlabelnumbifier++)
            {
                if(this->LNCC_status)
                {
                    num_true[(int)(inputCLASSIFIERptr[i+(int)LNCC[i+currlabelnumbifier*(this->numel)]*(this->numel)])]++;
                }
                else if(this->NCC_status)
                {
                    num_true[(int)(inputCLASSIFIERptr[i+(int)NCC[currlabelnumbifier]*(this->numel)])]++;
                }
                else
                {
                    num_true[(int)(inputCLASSIFIERptr[i+(int)currlabelnumbifier*(this->numel)])]++;
                }
            }
            long maxClass=0;
            long maxCount=0;
            for(long currLabel=0; currLabel<this->NumberOfLabels; currLabel++)
            {
                if(num_true[currLabel]>maxCount)
                {
                    maxClass=currLabel;
                    maxCount=num_true[currLabel];
                }
            }
            this->FinalSeg[i]=(float)maxClass;
        }
    }
    delete [] num_true;

    if(this->verbose_level>1)
    {
        cout<< " - Done"<<endl;
        flush(cout);
    }

    if(this->verbose_level>1)
    {
        cout<< "Allocating MRF";
        flush(cout);
    }
    if(this->MRF_status)
    {
        this->MRF=new segPrecisionTYPE [this->sizeAfterMaskingAndUncertainty*this->NumberOfLabels];
        if(MRF == NULL)
        {
            fprintf(stderr,"* The variable MRF was not allocated: OUT OF MEMORY!");
            exit(1);
        }
        for(long i=0; i<(this->sizeAfterMaskingAndUncertainty*this->NumberOfLabels); i++)
            this->MRF[i]=1.0f/this->NumberOfLabels;
    }
    if(this->verbose_level>1)
    {
        cout<< " - Done"<<endl;
        flush(cout);
    }

    return 0;

}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

nifti_image * seg_LabFusion::GetResult_probability()
{
    nifti_image * Result = nifti_copy_nim_info(this->inputCLASSIFIER);
    Result->dim[0]=4;
    Result->dim[4]=this->NumberOfLabels;
    Result->dim[5]=1;
    nifti_update_dims_from_array(Result);

    Result->cal_min=10.0e30;
    Result->cal_max=-10.0e30;
    Result->datatype=DT_FLOAT;

    Result->cal_min=0;
    Result->cal_max=this->LableCorrespondences_small_to_big[this->NumberOfLabels];
    //classifier_datatype * inputCLASSIFIERptr = static_cast<classifier_datatype *>(this->inputCLASSIFIER->data);

    if(this->FilenameOut.empty())
    {
        this->FilenameOut.assign("LabFusion.nii.gz");
    }

    nifti_set_filenames(Result,(char*)this->FilenameOut.c_str(),0,0);
    nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
    Result->cal_max=(1);
    Result->cal_min=(1);
    Result->data = (void *) calloc(Result->nx*Result->ny*Result->nz*Result->nt, sizeof(float));
    float * Resultdata = static_cast<float *>(Result->data);
    if(TYPE_OF_FUSION==1 || TYPE_OF_FUSION==3 )
    {
        int uncertainindex=0;
        for(long i=0; i<(this->numel); i++)
        {
            if(this->maskAndUncertainIndeces[i]>=0)
            {

                for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                {
                Resultdata[i+currlabelnumb*Result->nx*Result->ny*Result->nz]=
                        W[uncertainindex+currlabelnumb*this->sizeAfterMaskingAndUncertainty];
                }
                uncertainindex++;
            }
            else
            {
                Resultdata[i]=(int)this->FinalSeg[i]>1?1:(int)this->FinalSeg[i];
            }
        }

    }
    else if(TYPE_OF_FUSION==2)
    {
        if(this->NumberOfLabels>2)
            fprintf(stderr,"* Warning: Probabilistic output for MV is only implemented for binary fusion (2 lab max)");
        for(long i=0; i<(this->numel); i++)
        {
            Resultdata[i]=W[i];
        }

    }
    return Result;
}

nifti_image * seg_LabFusion::GetResult_label()
{
    nifti_image * Result = nifti_copy_nim_info(this->inputCLASSIFIER);
    Result->cal_min=0;
    Result->cal_max=this->LableCorrespondences_small_to_big[this->NumberOfLabels];
    Result->datatype=DT_FLOAT;
    if(this->FilenameOut.empty())
    {
        this->FilenameOut.assign("LabFusion.nii.gz");
    }

    nifti_set_filenames(Result,(char*)this->FilenameOut.c_str(),0,0);
    nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
    Result->cal_max=(this->NumberOfLabels-1);


    if(this->verbose_level>0)
    {
        cout << "Saving Integer Fused Label"<<endl;
    }

    if(TYPE_OF_FUSION==1)
    {
        Result->data = this->FinalSeg;
        this->FinalSeg=NULL;
    }
    else
    {
        Result->data = (void *) calloc(Result->nvox, sizeof(float));
    }
    Result->dim[0]=4;
    Result->dim[4]=1;
    nifti_update_dims_from_array(Result);
    float * Resultdata = static_cast<float *>(Result->data);
    //cout << "TYPE_OF_FUSION = "<<TYPE_OF_FUSION<<endl;

    if(TYPE_OF_FUSION==3)
    {
        for(long i=0; i<(this->numel); i++)
        {
            segPrecisionTYPE wmax=-1;
            int wmaxindex=0;
            for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
            {
                if(wmax<=W[i+currlabelnumb*this->numel])
                {
                    wmax=W[i+currlabelnumb*this->numel];
                    wmaxindex=currlabelnumb;
                }
            }
            Resultdata[i]=this->LableCorrespondences_small_to_big[wmaxindex];


        }
    }
    else if(TYPE_OF_FUSION==2)
    {
        for(long i=0; i<(this->numel); i++)
        {
            Resultdata[i]=this->LableCorrespondences_small_to_big[(int)round(W[i])];
        }

    }
    else if( TYPE_OF_FUSION==1)
    {
        int uncertainindex=0;
        for(long i=0; i<(this->numel); i++)
        {
            if(this->maskAndUncertainIndeces[i]>=0)
            {
                segPrecisionTYPE wmax=-1;
                int wmaxindex=0;
                for(long currlabelnumb=0; currlabelnumb<this->NumberOfLabels; currlabelnumb++)
                {
                    if(wmax<=W[uncertainindex+currlabelnumb*this->sizeAfterMaskingAndUncertainty])
                    {
                        wmax=W[uncertainindex+currlabelnumb*this->sizeAfterMaskingAndUncertainty];
                        wmaxindex=currlabelnumb;
                    }
                }
                Resultdata[i]=this->LableCorrespondences_small_to_big[wmaxindex];
                uncertainindex++;
            }
            else
            {
                Resultdata[i]=this->LableCorrespondences_small_to_big[(int)Resultdata[i]];
            }
        }
    }
    else
    {
        cout << "ERROR: Variable TYPE_OF_FUSION is not set. Cannot save."<<endl;
    }

    return Result;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
int  seg_LabFusion::Run_MV()
{

    TYPE_OF_FUSION=2;
    if((int)(this->verbose_level)>(int)(0))
    {
        cout << "Majority Vote: Verbose level " << this->verbose_level << endl;
        if(this->NCC_status )
        {
            cout << "Number of raters used = "<< this->Numb_Neigh<<endl;
        }
        else if(this->LNCC_status)
        {
            cout << "Number of local raters used = "<< this->Numb_Neigh<<endl;
        }
        else
        {
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
    if((int)(this->verbose_level)>(int)(0))
    {
        cout << "SBA: Verbose level " << this->verbose_level << endl;
        if(this->NCC_status )
        {
            cout << "Number of labels used = "<< this->Numb_Neigh<<endl;
        }
        else if(this->LNCC_status )
        {
            cout << "Number of labels used = "<< this->Numb_Neigh<<endl;
        }
        else
        {
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
    if((int)(this->verbose_level)>(int)(0))
    {
        cout << "STAPLE: Verbose level " << this->verbose_level << endl;
        if(this->NCC_status )
        {
            cout << "Number of raters used = "<< this->Numb_Neigh<<endl;
            cout << "MRF_beta = " << this->MRF_strength<<endl;
        }
        else if(this->LNCC_status)
        {
            cout << "Number of local raters used = "<< this->Numb_Neigh<<endl;
            cout << "MRF_beta = " << this->MRF_strength<<endl;
        }
        else
        {
            cout << "Using all classifiers"<<endl;
            cout << "MRF_beta = " << this->MRF_strength<<endl;
        }
    }
    if(this->verbose_level>0)
    {
        cout <<"Uncertainty area threshold is "<<this->uncertainthresh<<"percent"<<endl;
    }

    Allocate_Stuff_STAPLE();

    if(!(this->Fixed_Prop_status))
    {
        //UpdateDensity();
        EstimateInitialDensity();
        // UpdateDensity_noTest();

    }
    //**************
    // EM Algorithm
    //**************



    if(this->verbose_level>0)
    {
        cout << endl << "*******************************" << endl;
        cout << "Initialising " << endl;
    }
    this->iter=0;
    STAPLE_STEPS_Multiclass_Expectation_Maximization();
    //if(this->verbose_level>0)printloglik(this->iter,this->tracePQ,this->oldTracePQ);
    this->oldTracePQ=this->tracePQ;

    this->iter=1;
    bool out= true;
    while (out)
    {
        //out= false;
        if(this->verbose_level>0)
        {
            cout << endl << "*******************************" << endl;
            cout << "Iteration " << iter << endl;
        }

        // Iterative Components - EM, MRF



        //MRF
        //out=false;
        UpdateMRF();
        //Update Density
        UpdateDensity();
        //Expectation & Maximization
        STAPLE_STEPS_Multiclass_Expectation_Maximization();
        // Print Trace depending on the verbose level
        if(this->verbose_level>0)
        {
            if(this->iter>0)
            {
                cout<< "Normalized Trace Change = " << fabs((this->tracePQ-this->oldTracePQ)/fabsf(this->oldTracePQ)) << endl;
            }
            else
            {
                cout<< "Initial Normalized Trace Value = " << setprecision(7) << this->tracePQ << endl;
            }
        }

        // EXIT CHECKS
        // Check convergence
        if( (fabs(((this->tracePQ-this->oldTracePQ)/this->oldTracePQ))<=this->Conv && this->iter>1) || iter>=this->maxIteration || isnan(this->tracePQ) )out=false;
        // Exit if this->Numb_Neigh==1
        if(this->Numb_Neigh==1)out=false;


        // Update LogLik
        this->oldTracePQ=this->tracePQ;
        iter++;
    }

    return 0;
}
