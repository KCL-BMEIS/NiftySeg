#include "_seg_ProbLabFusion.h"




seg_LabFusion::seg_LabFusion(int _numb_lables)
{

    inputLABELS = NULL; // pointer to external
    inputImage_status=0;
    FilenameOut="LabFusion.nii.gz";
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


    this->Thresh_IMG_value=0;
    this->Thresh_IMG_DO=false;


    // SegParameters
    Q=new LabFusion_datatype [_numb_lables];
    P=new LabFusion_datatype [_numb_lables];
    for (int i=0;i<_numb_lables;i++){
        this->Q[i]=0.95;
        this->P[i]=0.95;
    }

    this->W=NULL;
    this->uncertainarea=NULL;
    this->Prop=0;
    this->Fixed_Prop_status=false;
    this->PropUpdate=false;
    this->numb_lables=_numb_lables;
    this->loglik=1;
    this->oldloglik=0;
    this->maxIteration=100;
    this->LNCC=NULL;
    this->NCC=NULL;
    this->Numb_Neigh=5;
    this->LNCC_status=false;
    this->NCC_status=false;


    // MRF Specific
    this->MRF_status=0;
    this->MRF_strength=0.0f;
    MRF=NULL;

}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

seg_LabFusion::~seg_LabFusion()
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

int seg_LabFusion::SetInputLabels(nifti_image *r,bool UNCERTAINflag)
{

    if(r->datatype!=DT_BINARY){
        seg_convert2binary(r,0.5f);
    }

    this->inputLABELS = r;
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

    int dim_array[3];
    dim_array[0]=(int)r->nx;
    dim_array[1]=(int)r->ny;
    dim_array[2]=(int)r->nz;

    this->uncertainarea=new bool [this->numel];
    this->W=new LabFusion_datatype [this->numel];

    for(int i=0;i<(this->numel);i++){
        this->uncertainarea[i]=true;
    }
    if(UNCERTAINflag){
        bool * inputLABELSptr = static_cast<bool *>(this->inputLABELS->data);

        for(int i=0;i<(this->numel);i++){
            int num_true=0;
            int num_false=0;
            for(int lable=0; lable<this->numb_lables;lable++){
                if(inputLABELSptr[i+lable*(this->numel)]){
                    num_true++;
                }
                else{
                    num_false++;
                }
            }
            if (num_true==this->numb_lables ){
                this->uncertainarea[i]=false;
                this->W[i]=1;
            }
            else if(num_false==this->numb_lables){
                this->uncertainarea[i]=false;
                this->W[i]=0;
            }
            else{
                this->uncertainarea[i]=true;
            }
        }
    }
    //Dillate(this->uncertainarea,1,dim_array);

    return 0;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::SetLNCC(nifti_image * _LNCC,nifti_image * BaseImage,float distance,int Numb_Neigh)
{
    if((Numb_Neigh<_LNCC->nt) & (Numb_Neigh>0)){
        this->Numb_Neigh=(int)Numb_Neigh;
    }
    else{
        this->Numb_Neigh=_LNCC->nt;
    }

    if(_LNCC->nt==this->numb_lables){
        if(_LNCC->datatype==DT_FLOAT){
            this->LNCC=seg_norm4LNCC(BaseImage,_LNCC,distance,Numb_Neigh,CurrSizes,this->verbose_level);
            this->LNCC_status = true;

        }
        else{
            cout << "LNCC is not float"<< endl;
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


    if(_GNCC->nt==this->numb_lables){
        if(_GNCC->datatype==DT_FLOAT){
            this->NCC=seg_norm_4D_GNCC(BaseImage,_GNCC,Numb_Neigh,CurrSizes,this->verbose_level);
            this->NCC_status = true;
        }
        else{
            cout << "GNCC is not float"<< endl;
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


    if(_ROINCC->nt==this->numb_lables){
        if(_ROINCC->datatype==DT_FLOAT){
            this->NCC=seg_norm4ROINCC(this->inputLABELS,BaseImage,_ROINCC,Numb_Neigh,CurrSizes,DilSize,this->verbose_level);
            this->NCC_status = true;
        }
        else{
            cout << "ROINCC is not float"<< endl;
        }
    }


    return 0;
}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::SetProp(float r)
{
    this->Prop = (LabFusion_datatype)r;
    this->Fixed_Prop_status = true;
    return 0;
}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::SetConv(float r)
{
    if(this->verbose_level){
        cout<< "Convergence Ratio = " << r <<endl;
        flush(cout);
    }
    this->Conv = (LabFusion_datatype)r;
    return 0;
}

int seg_LabFusion::SetImgThresh(float _Thresh_IMG_value)
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

int seg_LabFusion::SetMaximalIterationNumber(unsigned int numberiter)
{
    if(numberiter<3){
        this->maxIteration=3;
        cout << "Warning: It will only stop at iteration 3. For less than 3 iterations, use majority voting."<< endl;
    }
    else{
        this->maxIteration=numberiter;
    }
    return 0;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::Turn_MRF_ON(float strength)
{
    this->MRF_status=true;
    this->MRF_strength=(LabFusion_datatype)strength;
    return 0;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::Create_CurrSizes()
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

int seg_LabFusion::STAPLE_Maximization()
{


    LabFusion_datatype sumW=0;
    LabFusion_datatype sumWinv=0;

    // Get pointers to lables
    bool * inputLABELSptr = static_cast<bool *>(this->inputLABELS->data);

    LabFusion_datatype sumWclass=0;
    LabFusion_datatype sumWinvclass=0;
    if(this->verbose_level>0){
        cout << "[lab]  =  P \t\tQ"<< endl;
        flush(cout);
    }
    for(int lable=0; lable<this->CurrSizes->numclass;lable++){
        sumWclass=0;
        sumWinvclass=0;
        sumW=0;
        sumWinv=0;
        bool nccexists=false;
        if(this->LNCC_status){
            for(int i=0;i<(this->CurrSizes->numel);i++){
                if(this->uncertainarea[i]){
                    bool lnccexists=false;
                    for(char lnccindex=0;lnccindex<this->Numb_Neigh;lnccindex++){
                        if(lable==(this->LNCC[i+lnccindex*this->CurrSizes->numel])){
                            lnccexists=true;
                            lnccindex=this->Numb_Neigh;
                        }
                    }

                    if(lnccexists){
                        sumWclass+=(inputLABELSptr[i+lable*this->numel])?W[i]:0;
                        sumWinvclass+=(!(inputLABELSptr[i+lable*this->numel]))?(1-W[i]):0;
                        sumW+=W[i];
                        sumWinv+=(1-W[i]);
                    }
                }
            }
        }
        else if(this->NCC_status){

            for(int nccindex=0;nccindex<this->Numb_Neigh;nccindex++){
                if(lable==(this->NCC[nccindex])){
                    nccexists=true;
                    nccindex=this->Numb_Neigh;
                }
            }

            for(int i=0;i<(this->CurrSizes->numel);i++,inputLABELSptr[i+lable*this->numel]++){
                if(this->uncertainarea[i]){

                    if(nccexists){
                        sumWclass+=(inputLABELSptr[i+lable*this->numel])?W[i]:0;
                        sumWinvclass+=(!(inputLABELSptr[i+lable*this->numel]))?(1-W[i]):0;
                        sumW+=W[i];
                        sumWinv+=(1-W[i]);
                    }
                }
            }
        }
        else{
            for(int i=0;i<(this->CurrSizes->numel);i++,inputLABELSptr[i+lable*this->numel]++){
                if(this->uncertainarea[i]){
                    sumWclass+=(inputLABELSptr[i+lable*this->numel])?W[i]:0;
                    sumWinvclass+=(!(inputLABELSptr[i+lable*this->numel]))?(1-W[i]):0;
                    sumW+=W[i];
                    sumWinv+=(1-W[i]);
                }
            }
        }
        Q[lable]=sumWinvclass/sumWinv;
        P[lable]=sumWclass/sumW;
        if(this->verbose_level>0){
            if(this->NCC_status){
                if(nccexists){
                    cout << "[ "<<lable+1<<" ]  =  "<<  setprecision(4)<<P[lable] << "\t" << Q[lable]<< endl;
                }
            }
            else{
                cout << "[ "<<lable+1<<" ]  =  "<<  setprecision(4)<<P[lable] << "\t" << Q[lable]<< endl;
            }
        }
    }

    return 1;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::STAPLE_Expectation()
{
    LabFusion_datatype DIJ1=1;
    LabFusion_datatype DIJ0=1;
    LabFusion_datatype DIQ1=1;
    LabFusion_datatype DIQ0=1;
    LabFusion_datatype tmpW=0;

    bool * inputLABELSptr = static_cast<bool *>(this->inputLABELS->data);


    this->loglik=0;
    if(this->MRF_status){
        for(int i=0;i<(this->CurrSizes->numel);i++){
            if(this->uncertainarea[i]){

                DIJ1=1;
                DIJ0=1;
                DIQ1=1;
                DIQ0=1;
                if(this->LNCC_status){
                    for(int lable=0; lable<this->Numb_Neigh;lable++){
                        int LNCCvalue=(int)(this->LNCC[i+lable*this->CurrSizes->numel]);
                        if(LNCCvalue>=0 && LNCCvalue<this->CurrSizes->numclass){
                            if(inputLABELSptr[i+LNCCvalue*this->CurrSizes->numel] ){
                                DIJ1*=(LabFusion_datatype)(P[LNCCvalue]);
                                DIQ1*=(LabFusion_datatype)(1-Q[LNCCvalue]);
                            }
                            else{
                                DIJ0*=(LabFusion_datatype)(1-P[LNCCvalue]);
                                DIQ0*=(LabFusion_datatype)(Q[LNCCvalue]);
                            }
                        }
                    }
                }
                else if(this->NCC_status){
                    for(int lable=0; lable<this->Numb_Neigh;lable++){
                        int NCCvalue=(int)(this->NCC[lable]);
                        if(NCCvalue>=0 && NCCvalue<this->CurrSizes->numclass){
                            if(inputLABELSptr[i+NCCvalue*this->CurrSizes->numel] ){
                                DIJ1*=(LabFusion_datatype)(P[NCCvalue]);
                                DIQ1*=(LabFusion_datatype)(1-Q[NCCvalue]);
                            }
                            else{
                                DIJ0*=(LabFusion_datatype)(1-P[NCCvalue]);
                                DIQ0*=(LabFusion_datatype)(Q[NCCvalue]);
                            }
                        }
                    }
                }
                else{
                    for(int lable=0; lable<this->CurrSizes->numclass;lable++){
                        if(inputLABELSptr[i+lable*this->numel]){
                            DIJ1*=(LabFusion_datatype)P[lable];
                            DIQ1*=(LabFusion_datatype)(1.0-Q[lable]);
                        }
                        else{
                            DIJ0*=(LabFusion_datatype)(1-P[lable]);
                            DIQ0*=(LabFusion_datatype)(Q[lable]);
                        }
                    }
                }
                tmpW=(((LabFusion_datatype)this->Prop)*(LabFusion_datatype)MRF[i]*DIJ0*DIJ1)/(((LabFusion_datatype)this->Prop)*(LabFusion_datatype)MRF[i]*DIJ0*DIJ1+(1-(LabFusion_datatype)this->Prop)*(LabFusion_datatype)MRF[i+this->numel]*DIQ0*DIQ1);
                if(tmpW<0)tmpW=0;
                if(tmpW>1)tmpW=1;
                this->loglik+=W[i];
                this->W[i]=tmpW;
            }
        }

    }
    else{
        for(int i=0;i<(this->CurrSizes->numel);i++){
            if(this->uncertainarea[i]){
                DIJ1=1;
                DIJ0=1;
                DIQ1=1;
                DIQ0=1;
                if(this->LNCC_status){
                    for(int lable=0; lable<this->Numb_Neigh;lable++){
                        int LNCCvalue=(int)(this->LNCC[i+lable*this->CurrSizes->numel]);
                        if(LNCCvalue>=0 && LNCCvalue<this->CurrSizes->numclass){
                            if(inputLABELSptr[i+LNCCvalue*this->CurrSizes->numel] ){
                                DIJ1*=(LabFusion_datatype)(P[LNCCvalue]);
                                DIQ1*=(LabFusion_datatype)(1-Q[LNCCvalue]);
                            }
                            else{
                                DIJ0*=(LabFusion_datatype)(1-P[LNCCvalue]);
                                DIQ0*=(LabFusion_datatype)(Q[LNCCvalue]);
                            }
                        }
                    }
                }
                else if(this->NCC_status){
                    for(int lable=0; lable<this->Numb_Neigh;lable++){
                        int NCCvalue=(int)(this->NCC[lable]);
                        if(NCCvalue>=0 && NCCvalue<this->CurrSizes->numclass){
                            if(inputLABELSptr[i+NCCvalue*this->CurrSizes->numel] ){
                                DIJ1*=(LabFusion_datatype)(P[NCCvalue]);
                                DIQ1*=(LabFusion_datatype)(1-Q[NCCvalue]);
                            }
                            else{
                                DIJ0*=(LabFusion_datatype)(1-P[NCCvalue]);
                                DIQ0*=(LabFusion_datatype)(Q[NCCvalue]);
                            }
                        }
                    }
                }
                else{
                    for(int lable=0; lable<this->CurrSizes->numclass;lable++){
                        if(inputLABELSptr[i+lable*this->numel]){
                            DIJ1*=(LabFusion_datatype)P[lable];
                            DIQ1*=(LabFusion_datatype)(1.0-Q[lable]);
                        }
                        else{
                            DIJ0*=(LabFusion_datatype)(1-P[lable]);
                            DIQ0*=(LabFusion_datatype)(Q[lable]);
                        }
                    }
                }
                W[i]=(((LabFusion_datatype)this->Prop)*DIJ0*DIJ1)/((LabFusion_datatype)this->Prop*DIJ0*DIJ1+(1-(LabFusion_datatype)this->Prop)*DIQ0*DIQ1);

                if(this->W[i]<0){
                    this->W[i]=0;
                }
                if(this->W[i]>1){
                    this->W[i]=1;
                }
                this->loglik+=W[i];
            }
        }
    }
    return 1;
}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::MV_Estimate()
{

    LabFusion_datatype tmpW=0;

    bool * inputLABELSptr = static_cast<bool *>(this->inputLABELS->data);

    this->loglik=0;
    for(int i=0;i<(this->CurrSizes->numel);i++){
        if(this->LNCC_status){
            tmpW=0.0f;
            for(int lable=0; lable<this->Numb_Neigh;lable++){
                int LNCCvalue=(int)(this->LNCC[i+lable*this->CurrSizes->numel]);
                if(LNCCvalue>=0 && LNCCvalue<this->CurrSizes->numclass){
                    if(inputLABELSptr[i+LNCCvalue*this->CurrSizes->numel] ){
                        tmpW++;
                    }
                }
            }
            this->W[i]=tmpW/this->Numb_Neigh;

        }
        else if(this->NCC_status){
            tmpW=0.0f;
            for(int lable=0; lable<this->Numb_Neigh;lable++){
                int NCCvalue=(int)(this->NCC[lable]);
                if(NCCvalue>=0 && NCCvalue<this->CurrSizes->numclass){
                    if(inputLABELSptr[i+NCCvalue*this->CurrSizes->numel] ){
                        tmpW++;
                    }
                }
            }
            this->W[i]=tmpW/this->Numb_Neigh;
        }
        else{
            tmpW=0.0f;
            for(int lable=0; lable<this->CurrSizes->numclass;lable++){
                if(inputLABELSptr[i+lable*this->numel]){
                    tmpW++;
                }
            }
            this->W[i]=tmpW/this->Numb_Neigh;
        }
    }
    return 1;

}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::UpdateMRF()
{
    if(this->MRF_status){
        if(this->verbose_level>0){
            cout<<"Updating MRF"<<endl;
            flush(cout);
        }
        int maxix = (int)(CurrSizes->xsize);
        int maxiy = (int)(CurrSizes->ysize);
        int maxiz = (int)(CurrSizes->zsize);     
        int indexCentre=0;
        int planesize=(int)CurrSizes->xsize*(int)CurrSizes->ysize;
        LabFusion_datatype * Wptr=NULL;
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
        LabFusion_datatype * MRFptr2=NULL;
        LabFusion_datatype * MRFptr1=NULL;
        for(int neighbour=0;neighbour<6;neighbour++){
            for (int iz=1; iz<maxiz-1; iz++) {
                for (int iy=1; iy<maxiy-1; iy++) {
                    indexCentre=1+iy*(int)CurrSizes->xsize+iz*planesize;
                    Wptr=&W[indexCentre+neighShift[neighbour]];
                    MRFptr2=&MRF[indexCentre+this->numel];
                    MRFptr1=&MRF[indexCentre];
                    for (int ix=1; ix<(maxix-1); ix++,Wptr++,MRFptr2++,MRFptr1++){
                        *MRFptr2+=(*Wptr);
                        *MRFptr1+=(1-(*Wptr));
                    }
                }
            }
        }

        Wptr=&W[0];
        MRFptr1=&MRF[0];
        MRFptr2=&MRF[this->numel];
        Wptr=&W[0];
        for(int i=0;i<this->numel;i++,MRFptr1++,MRFptr2++,Wptr++){
            (*MRFptr2) = exp((-this->MRF_strength)*(*MRFptr2)*(*Wptr));
            (*MRFptr1) = exp((-this->MRF_strength)*(*MRFptr1)*(*Wptr));
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

int seg_LabFusion::EstimateInitialDensity()
{
    bool * inputLABELSptr = static_cast<bool *>(this->inputLABELS->data);

    this->Prop=0;
    int tempprop=0;
    int tempsum=0;
    for(int clas=0; clas<this->numb_lables; clas++){
        for( int i=0; i<this->numel; i++){

            if(this->uncertainarea[i]){
                if(inputLABELSptr[i+this->numel*clas]){
                tempprop++;
                }
                tempsum++;
            }
        }
    }
    this->Prop=(float)(tempprop)/(float)(tempsum);

    if(this->verbose_level>0){
        cout << "Estimated initial proportion = "<<this->Prop<<endl;
        flush(cout);

    }
    return 1;

}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::SetPQ(float tmpP,float tmpQ)
{
    for(int i=0; i<this->numb_lables;i++){
        this->P[i]=(LabFusion_datatype)tmpP;
        this->Q[i]=(LabFusion_datatype)tmpQ;
    }
    return 1;

}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_LabFusion::UpdateDensity()
{
    if(this->PropUpdate){
        this->Prop=0;
        int tempsum=0;
        for(int i=0; i<this->numel; i++){
            if(this->uncertainarea[i]){
                tempsum++;
                this->Prop+=this->W[i];
            }
        }
        this->Prop=this->Prop/(tempsum);

        if(this->verbose_level>0){
            cout << "Proportion = " << this->Prop << endl;
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

int seg_LabFusion::Allocate_Stuff()
{

    if(this->MRF_status){
        this->MRF=new LabFusion_datatype [this->numel*2];
        for(int i=0; i<this->numel*2; i++)
            this->MRF[i]=0.5;
    }

    return 0;

}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

nifti_image * seg_LabFusion::GetResult()
{
    nifti_image * Result = Copy_single_image_to_Result(this->W,this->inputLABELS,(char*)this->FilenameOut.c_str());
    if(this->Thresh_IMG_DO){
        float * Resultdata = static_cast<float *>(Result->data);
        for(unsigned int i=0; i<Result->nvox; i++){
            Resultdata[i]=(float)(Resultdata[i]>=this->Thresh_IMG_value);
        }

    }
    return Result;
    return 0;
}



/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
int  seg_LabFusion::Run_MV()
{

    if((int)(this->verbose_level)>(int)(0)){
        cout << "Majority Vote: Verbose level " << this->verbose_level << endl;
        if(this->NCC_status ){
            cout << "Number of labels used = "<< this->Numb_Neigh<<endl;
        }
        else if(this->LNCC_status){
            cout << "Number of local labels used = "<< this->Numb_Neigh<<endl;
        }
        else{
            cout << "Using all lables"<<endl;
        }
    }


    Allocate_Stuff();
    MV_Estimate();


    return 0;
}

int  seg_LabFusion::Run_SBA()
{
    cout << "Not Implemented" << this->verbose_level << endl;
    return 0;
}

int  seg_LabFusion::Run_STAPLE()
{

    if((int)(this->verbose_level)>(int)(0)){
        cout << "STAPLE: Verbose level " << this->verbose_level << endl;
        if(this->NCC_status ){
            cout << "Number of labels used = "<< this->Numb_Neigh<<endl;
        }
        else if(this->LNCC_status){
            cout << "Number of local labels used = "<< this->Numb_Neigh<<endl;
        }
        else{
            cout << "Using all lables"<<endl;
        }
    }

    if(!(this->Fixed_Prop_status)){
        EstimateInitialDensity();
    }

    Allocate_Stuff();
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
        STAPLE_Expectation();
        //Maximization
        STAPLE_Maximization();
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

    return 0;
}
