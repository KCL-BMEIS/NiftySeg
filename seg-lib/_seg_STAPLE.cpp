#ifndef _SEG_STAPLE_CPP
#define _SEG_STAPLE_CPP
#include "_seg_STAPLE.h"


template <class T>
        seg_STAPLE<T>::seg_STAPLE(int _numb_lables)
{

    inputLABLES = NULL; // pointer to external
    inputImage_status=0;
    FilenameOut="STAPLE.nii.gz";
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
    this->Conv=0;
    this->numel=0;
    this->iter=0;
    this->CurrSizes=NULL;

    // SegParameters
    for (int i=0;i<_numb_lables;i++){
        this->Q[i]=0.99;
        this->P[i]=0.99;
    }

    this->W=NULL;
    this->Prop=0;
    this->Fixed_Prop_status=false;
    this->PropUpdate=false;
    this->numb_lables=_numb_lables;
    this->loglik=1;
    this->oldloglik=0;
    this->maxIteration=100;
    this->LNCC=NULL;
    this->lnccthresh=0.05;
    this->LNCC_status=false;


    // MRF Specific
    this->MRF_status=0;
    this->MRF_strength=0.0f;
    MRF=NULL;

}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
template <class T>
        seg_STAPLE<T>::~seg_STAPLE()
{

    if(this->W!=NULL){
        delete[] this->W;
    }
    this->W=NULL;

    if(this->MRF!=NULL){
        delete[] this->MRF;
    }
    this->MRF=NULL;

    if(this->CurrSizes!=NULL)
        delete [] this->CurrSizes;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
template<class T>
int seg_STAPLE<T>::SetInputLables(nifti_image *r)
{

    if(r->datatype!=DT_BINARY){
        seg_convert2binary(r,0.0f);
    }

    this->inputLABLES = r;
    this->inputImage_status = true;
    // Size
    this->dimentions=(int)((r->nx)>1)+(int)((r->ny)>1)+(int)((r->nz)>1)+(int)((r->nt)>1)+(int)((r->nu)>1);
    this->nx=r->nx;
    this->ny=r->ny;
    this->nz=r->nz;
    this->dx=r->dx;
    this->dy=r->dy;
    this->dz=r->dz;
    this->numel=r->nz*r->ny*r->nx;
    if(this->CurrSizes==NULL) Create_CurrSizes();
    return 0;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
template<class T>
int seg_STAPLE<T>::SetPrior(nifti_image * _Prior,float _lnccthresh)
{
    PrecisionTYPE * _PriorPTR=NULL;
    if(_lnccthresh<1& _lnccthresh>0){
        this->lnccthresh=_lnccthresh;
    }
    else{
        this->lnccthresh=0.05;
    }
    if(_Prior->datatype==DT_FLOAT){
        _PriorPTR = static_cast<float *>(_Prior->data);
        this->LNCC=_PriorPTR;
        this->LNCC_status = true;
        // Size

    }
    else{
        cout << "Prior is not float"<< endl;
    }



    return 0;
}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
template<class T>
int seg_STAPLE<T>::SetLNCC(nifti_image * _LNCC,nifti_image * BaseImage,int distance,float _lnccthresh,bool saveLNCC)
{
    if(_lnccthresh<1& _lnccthresh>0){
        this->lnccthresh=_lnccthresh;
    }
    else{
        this->lnccthresh=0.05;
    }
    string filenameLNCC="LNCC.nii";
    if(_LNCC->nt==this->numb_lables){
        if(_LNCC->datatype==DT_FLOAT){
            this->LNCC=seg_norm4NCC(BaseImage,_LNCC,this->inputLABLES,distance,CurrSizes);
            this->LNCC_status = true;
            if(saveLNCC){
                nifti_set_filenames(_LNCC,filenameLNCC.c_str(),0,0);
                nifti_image_write(_LNCC);
            }
        }
        else{
            cout << "LNCC is not float"<< endl;
        }
    }

    return 0;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
template<class T>
int seg_STAPLE<T>::SetProp(float r)
{
    this->Prop = r;
    this->Fixed_Prop_status = true;
}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
template<class T>
int seg_STAPLE<T>::SetConv(float r)
{
    if(this->verbose_level){
        cout<< "Convergence Ratio = " << r <<endl;
        flush(cout);
    }
    this->Conv = r;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
template<class T>
int seg_STAPLE<T>::SetFilenameOut(char *f)
{
    this->FilenameOut = f;
    return 0;
}



/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
template<class T>
int seg_STAPLE<T>::SetVerbose(unsigned int verblevel)
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
template<class T>
int seg_STAPLE<T>::SetMaximalIterationNumber(unsigned int numberiter)
{
    if(numberiter<3){
        this->maxIteration=3;
        cout << "Warning: It will only stop at iteration 3"<< endl;
    }
    else{
        this->maxIteration=numberiter;
    }
    return 0;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
template<class T>
int seg_STAPLE<T>::Turn_MRF_ON(float strength)
{
    this->MRF_status=true;
    this->MRF_strength=strength;
    return 0;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
template<class T>
int seg_STAPLE<T>::Create_CurrSizes()
{
    this->CurrSizes = new ImageSize [1]();
    CurrSizes->numel=(int)(this->nx*this->ny*this->nz);
    CurrSizes->xsize=this->nx;
    CurrSizes->ysize=this->ny;
    CurrSizes->zsize=this->nz;
    CurrSizes->usize=1;
    CurrSizes->tsize=1;
    CurrSizes->numclass=this->numb_lables;
    CurrSizes->numelmasked=0;
    CurrSizes->numelbias=0;
    return 0;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
template<class T>
int seg_STAPLE<T>::Maximization()
{


    float sumW=0;
    float sumWinv=0;

    // Get pointers to lables
    bool ** Lablepointer=new bool * [this->CurrSizes->numclass];
    bool * inputLABLESptr = static_cast<bool *>(this->inputLABLES->data);

    for(int lable=0; lable<this->CurrSizes->numclass;lable++){
        Lablepointer[lable]=&inputLABLESptr[this->CurrSizes->numel*lable];
    }


    float sumWclass=0;
    float sumWinvclass=0;
    if(this->verbose_level>0){
        cout << "[lab]  =  P \t\tQ"<< endl;
        flush(cout);
    }
    for(int lable=0; lable<this->CurrSizes->numclass;lable++){
        sumWclass=0;
        sumWinvclass=0;
        sumW=0;
        sumWinv=0;
        if(this->LNCC_status){
            for(int i=0;i<(this->CurrSizes->numel);i++,Lablepointer[lable]++){
                sumWclass+=(*Lablepointer[lable]&&this->LNCC[i+lable*this->CurrSizes->numel]>this->lnccthresh)?W[i]:0;
                sumWinvclass+=(!(*Lablepointer[lable])&&this->LNCC[i+lable*this->CurrSizes->numel]>this->lnccthresh)?(1-W[i]):0;
                sumW+=this->LNCC[i+lable*this->CurrSizes->numel]>this->lnccthresh?W[i]:0;
                sumWinv+=this->LNCC[i+lable*this->CurrSizes->numel]>this->lnccthresh?(1-W[i]):0;
            }
        }
        else{

            for(int i=0;i<(this->CurrSizes->numel);i++,Lablepointer[lable]++){
                sumWclass+=(*Lablepointer[lable])?W[i]:0;
                sumWinvclass+=(!(*Lablepointer[lable]))?(1-W[i]):0;
                sumW+=W[i];
                sumWinv+=(1-W[i]);
            }
        }

        P[lable]=sumWclass/sumW;
        Q[lable]=sumWinvclass/sumWinv;

        if(this->verbose_level>0){
            cout << "[ "<<lable<<" ]  =  "<< P[lable] << "\t" << Q[lable]<< endl;
        }
    }

    return 1;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
template<class T>
int seg_STAPLE<T>::Expectation()
{
    float DIJ1=1;
    float DIJ0=1;
    float DIQ1=1;
    float DIQ0=1;
    bool ** Lablepointer=new bool * [this->CurrSizes->numclass];
    bool * inputLABLESptr = static_cast<bool *>(this->inputLABLES->data);

    for(int lable=0; lable<this->CurrSizes->numclass;lable++){
        Lablepointer[lable]=&inputLABLESptr[this->CurrSizes->numel*lable];
    }
    this->loglik=0;
    float maxW=0;
    float minW=0;
    if(this->MRF_status){
        for(int i=0;i<(this->CurrSizes->numel);i++){
            DIJ1=1;
            DIJ0=1;
            DIQ1=1;
            DIQ0=1;

            if(!this->LNCC_status){
                for(int lable=0; lable<this->CurrSizes->numclass;lable++){
                    if(*Lablepointer[lable]){
                        DIJ1*=P[lable];
                        DIQ1*=(1-Q[lable]);
                    }
                    else{
                        DIJ0*=(1-P[lable]);
                        DIQ0*=(Q[lable]);
                    }
                    Lablepointer[lable]++;
                }
                W[i]=((this->Prop)*MRF[i]*DIJ0*DIJ1)/(this->Prop*MRF[i]*DIJ0*DIJ1+(1-this->Prop)*MRF[i+this->numel]*DIQ0*DIQ1);
            }
            else{
                for(int lable=0; lable<this->CurrSizes->numclass;lable++){
                    float LNCCvalue=(this->LNCC[i+lable*this->CurrSizes->numel]);
                    if(LNCCvalue>this->lnccthresh){
                        if(*Lablepointer[lable]){
                            DIJ1*=(P[lable]);
                            DIQ1*=(1-Q[lable]);
                        }
                        else{
                            DIJ0*=(1-P[lable]);
                            DIQ0*=(Q[lable]);
                        }
                    }
                    Lablepointer[lable]++;
                }
                W[i]=((this->Prop)*MRF[i]*DIJ0*DIJ1)/((this->Prop)*MRF[i]*DIJ0*DIJ1+(1-this->Prop)*MRF[i+this->numel]*DIQ0*DIQ1);
            }
            if(this->W[i]<0){
                this->W[i]=0;
            }
            if(this->W[i]>1){
                this->W[i]=1;
            }
            this->loglik+=W[i];

        }

    }
    else{
        for(int i=0;i<(this->CurrSizes->numel);i++){
            DIJ1=1;
            DIJ0=1;
            DIQ1=1;
            DIQ0=1;
            if(!this->LNCC_status){
                for(int lable=0; lable<this->CurrSizes->numclass;lable++){
                    if(*Lablepointer[lable]){
                        DIJ1*=P[lable];
                        DIQ1*=(1-Q[lable]);
                    }
                    else{
                        DIJ0*=(1-P[lable]);
                        DIQ0*=(Q[lable]);
                    }
                    Lablepointer[lable]++;
                }
                W[i]=((this->Prop)*DIJ0*DIJ1)/(this->Prop*DIJ0*DIJ1+(1-this->Prop)*DIQ0*DIQ1);
            }
            else{
                for(int lable=0; lable<this->CurrSizes->numclass;lable++){
                    float LNCCvalue=(this->LNCC[i+lable*this->CurrSizes->numel]);
                    if(LNCCvalue>this->lnccthresh){
                        if(*Lablepointer[lable]){
                            DIJ1*=(P[lable]);
                            DIQ1*=(1-Q[lable]);
                        }
                        else{
                            DIJ0*=(1-P[lable]);
                            DIQ0*=(Q[lable]);
                        }
                    }
                    Lablepointer[lable]++;
                }
                W[i]=((this->Prop)*DIJ0*DIJ1)/((this->Prop)*DIJ0*DIJ1+(1-this->Prop)*DIQ0*DIQ1);
            }

            if(this->W[i]<0){
                this->W[i]=0;
            }
            if(this->W[i]>1){
                this->W[i]=1;
            }
            this->loglik+=W[i];
        }
    }


    return 1;

}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
template<class T>
int seg_STAPLE<T>::UpdateMRF()
{
    if(this->MRF_status){
        if(this->verbose_level>0){
            cout<<"Updating MRF"<<endl;
            flush(cout);
        }
        int maxix = (int)(CurrSizes->xsize);
        int maxiy = (int)(CurrSizes->ysize);
        int maxiz = (int)(CurrSizes->zsize);
        unsigned int numel_currclass_shift[2];
        float U[2];

        for(int i=0; i<this->numb_lables; i++){
            numel_currclass_shift[i]=i*this->numel;
        }
        int indexCentre=0;
        int planesize=(int)CurrSizes->xsize*(int)CurrSizes->ysize;
        float * Wptr=NULL;

        for(int i=0;i<(this->numel*2);i++){
            MRF[i]=0;
        }

        int neighShift[6];
        neighShift[0]=1;
        neighShift[1]=-1;
        neighShift[2]=CurrSizes->xsize;
        neighShift[3]=-CurrSizes->xsize;
        neighShift[4]=planesize;
        neighShift[5]=-planesize;
        float * MRFptr2=NULL;
        float * MRFptr1=NULL;
        float tmpW=0;
        for(int neighbour=0;neighbour<6;neighbour++){
            for (int iz=1; iz<maxiz-1; iz++) {
                for (int iy=1; iy<maxiy-1; iy++) {
                    indexCentre=1+iy*(int)CurrSizes->xsize+iz*planesize;
                    Wptr=&W[indexCentre-neighShift[neighbour]];
                    MRFptr2=&MRF[indexCentre+this->numel];
                    MRFptr1=&MRF[indexCentre];
                    for (int ix=1; ix<(maxix-1); ix++,Wptr++,MRFptr2++,MRFptr1++){
                        *MRFptr2+=(-this->MRF_strength)*(*Wptr);
                        *MRFptr1+=(-this->MRF_strength)*(1-(*Wptr));
                    }
                }
            }
        }

        MRFptr1=&MRF[0];
        MRFptr2=&MRF[this->numel];
        for(int i=0;i<this->numel;i++,MRFptr1++,MRFptr2++){
            (*MRFptr2) = exp((*MRFptr2));
            (*MRFptr1) = exp((*MRFptr1));
            (*MRFptr1)=(*MRFptr1)/((*MRFptr1)+(*MRFptr2));
            (*MRFptr2)=(*MRFptr2)/((*MRFptr1)+(*MRFptr2));
        }
        if(this->verbose_level>0){
            cout<<"MRF Updated"<<endl;
            flush(cout);
        }
    }
    return 1;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
template<class T>
int seg_STAPLE<T>::EstimateInitialDensity()
{
    bool * inputLABLESptr = static_cast<bool *>(this->inputLABLES->data);

    this->Prop=0;
    for(int i=0; i<this->inputLABLES->nvox; i++){
        this->Prop+=(float)inputLABLESptr[i];
    }
    this->Prop=this->Prop/(float)(this->inputLABLES->nvox);

    if(this->verbose_level>0){
        cout << "Estimated initial proportion = "<<this->Prop<<endl;
        flush(cout);

    }
    return 1;

}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
template<class T>
int seg_STAPLE<T>::SetPQ(float tmpP,float tmpQ)
{
    for(int i=0; i<this->numb_lables;i++){
        this->P[i]=tmpP;
        this->Q[i]=tmpQ;
    }
    return 1;

}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
template<class T>
int seg_STAPLE<T>::UpdateDensity()
{
    if(this->PropUpdate && !this->LNCC_status){
        this->Prop=0;
        for(int i=0; i<this->numel; i++){
            this->Prop+=W[i];
        }
        this->Prop=this->Prop/(float)(this->numel);

        if(this->verbose_level>0){
            cout << "Proportion = " << this->Prop << endl;
            flush(cout);

        }
    }
    return 1;

}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
template<class T>
int seg_STAPLE<T>::Turn_Prop_Update_ON()
{
    this->PropUpdate=true;
    return 1;

}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
template<class T>
int seg_STAPLE<T>::Allocate_Expec_and_MRF()
{

    this->W=new float [this->numel];
    if(this->MRF_status){
        this->MRF=new float [this->numel*2];
        for(int i=0; i<this->numel*2; i++)
            this->MRF[i]=0.5;
    }

    return 0;

}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
template<class T>
nifti_image * seg_STAPLE<T>::GetResult()
{
    nifti_image * Result = Copy_single_image_to_Result(this->W,this->inputLABLES,(char*)this->FilenameOut.c_str());
    return Result;

}



/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

template<class T>
int *  seg_STAPLE<T>::Run_EM()
{
    time_t start,end;
    time(&start);
    if((int)(this->verbose_level)>(int)(0)){
        cout << "STAPLE: Verbose level " << this->verbose_level << endl;
        if(this->LNCC_status){
            cout << "LNCC thresh = "<< this->lnccthresh<<endl;
        }
    }

    if(!(this->Fixed_Prop_status)){
        EstimateInitialDensity();
    }

    Allocate_Expec_and_MRF();

    //**************
    // EM Algorithm
    //**************
    this->iter=0;
    bool out= true;

    while (out) {

        if(this->verbose_level>0){
            cout << endl << "*******************************" << endl;
            cout << "Iteration " << iter << endl;
        }

        // Iterative Components - EM, MRF
        //Expectation
        Expectation();
        //Maximization
        Maximization();
        //MRF
        UpdateMRF();
        //Update Density
        UpdateDensity();

        // Print LogLik depending on the verbose level
        if(this->verbose_level>0){
            printloglik(this->iter,this->loglik,this->oldloglik);
        }
        // Preform Segmentation Refinement Steps or Exit
        if( (((this->loglik-this->oldloglik)/this->oldloglik)<=this->Conv && this->iter>3) || iter>this->maxIteration || isnan(this->loglik) ){
            out=false;
        }

        // Update LogLik
        this->oldloglik=this->loglik;
        iter++;
    }

    time(&end);

    if(this->verbose_level>0){
        cout << "Finished in "<<difftime(end,start)<<"sec"<< endl;
    }
    return 0;
}


#endif
