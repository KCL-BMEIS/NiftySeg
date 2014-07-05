#ifndef _SEG_EM_CPP
#define _SEG_EM_CPP
#include "_seg_EM.h"


/// @brief Class constructor that takes in the number of classes and the number of multimodal data (nu*nt).
/// As multiple time points (nt) are treated in the same way as multimodal data (nu), the value of this->nt is set to 1,
/// and the value of this->nu is set to nu*nt
/// @param _numb_classes The number of classes
/// @param _nu The number of modalities per time point
/// @param _nt The number of time points
seg_EM::seg_EM(int _numb_classes, int _nu,int _nt)
{

    this->InputImage=NULL;  // pointer to external

    this->filenameOut="Segmentation.nii.gz";

    this->dimentions=1;
    this->nx=1;
    this->ny=0;
    this->nz=0;
    this->nu=_nu*_nt;
    this->nt=1;

    this->dx=0;
    this->dy=0;
    this->dz=0;

    this->numel=0;

    this->iter=0;
    this->ratio=1000;
    this->aprox=true;

    this->numb_classes=_numb_classes;
    this->M= new float [maxMultispectalSize*maxNumbClass];
    for(int i=0; i<(maxMultispectalSize*maxNumbClass); i++)
    {
        this->M[i]=0.0f;
    }
    this->V= new float [maxMultispectalSize*maxMultispectalSize*maxNumbClass];
    for(int i=0; i<(maxMultispectalSize*maxMultispectalSize*maxNumbClass); i++)
    {
        this->V[i]=0.0f;
    }

    this->Expec=NULL;
    this->ShortPrior=NULL;
    this->S2L=NULL;
    this->L2S=NULL;
    this->CurrSizes=NULL;
    this->reg_factor=1.1f;

    this->maxIteration=30;
    this->minIteration=4;
    this->verbose_level=0;
    this->loglik=2.0;
    this->oldloglik=1.0;


    this->maskImageStatus=0;
    this->Mask=NULL;        // pointer to external
    this->numelmasked=0;

    this->priorsStatus=false;
    this->Priors=NULL;      // pointer to external

    this->mrfStatus=false;
    this->mrfStrength=0.0f;
    this->MRF=NULL;
    this->MRFBeta=NULL;
    this->MRFTransitionMatrix=NULL;

    this->Outlierness=NULL;
    this->OutliernessUSE=NULL;
    this->outliernessStatus=false;
    this->OutliernessThreshold=0.0f;
    this->Outlierness_ratio=0.01f;

    this->biasFieldStatus=false;
    this->biasFieldOrder=0;
    this->BiasField=NULL;
    this->BiasField_coeficients=NULL;
    this->biasFieldRatio=0;

    this->stageLoAd=-1;
    this->pvModelStatus=false;
    this->sgDelineationStatus=false;
    this->relaxStatus=false;
    this->relaxFactor=0;
    this->relaxGaussKernelSize=0;
    this->mapStatus=0;
    this->MAP_M=NULL;
    this->MAP_V=NULL;
}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

/// @brief Class destructor which deletes all pointers that are different from null (deffiend) and owned by the seg_EM object.
/// The destructor does not delete the *inputImage, *Priors and *Mask pointers as they are not owned by the object.
seg_EM::~seg_EM()
{
    if(this->Expec!=NULL)
    {
        delete [] this->Expec;
    }
    this->Expec=NULL;

    if(this->ShortPrior!=NULL)
    {
        delete [] this->ShortPrior;
    }
    this->ShortPrior=NULL;

    if(this->BiasField!=NULL)
    {
        delete [] this->BiasField;
    }
    this->BiasField=NULL;

    if(this->BiasField_coeficients!=NULL)
    {
        delete [] this->BiasField_coeficients;
    }
    this->BiasField_coeficients=NULL;

    if(this->Outlierness!=NULL)
    {
        delete [] this->Outlierness;
    }

    if(this->mrfStatus)
    {
        delete [] this->MRF;
        delete [] this->MRFTransitionMatrix;
    }
    this->MRF=NULL;
    this->MRFTransitionMatrix=NULL;

    if(this->maskImageStatus>0)
    {
        delete [] this->L2S;
        delete [] this->S2L;
    }

    if(this->CurrSizes!=NULL)
        delete [] this->CurrSizes;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Defines the input image and all the image sizes.
/// @param _r The image pointer. This pointer is owned by the caller.
///
/// Note that the nifti_image pointer is owned by the caller, but the called should not free the pointer while seg_EM is still using it. More specifically, seg_EM will use _r internally, but will not free it when clearing the object.
void seg_EM::SetInputImage(nifti_image *_r)
{
    this->InputImage = _r;
    this->inputImage_status = true;
    // Count the number of dimentions with size above 1.
    this->dimentions=(int)((_r->nx)>1)+(int)((_r->ny)>1)+(int)((_r->nz)>1)+(int)((_r->nt)>1)+(int)((_r->nu)>1);
    this->nx=_r->nx;
    this->ny=_r->ny;
    this->nz=_r->nz;
    this->nt=(_r->nt>1)?_r->nt:1;
    this->nu=(_r->nu>1)?_r->nu:1;
    this->dx=_r->dx;
    this->dy=_r->dy;
    this->dz=_r->dz;
    this->numel=_r->nz*_r->ny*_r->nx;
    if(this->nx==1 ||this->ny==1)
    {
        cout<<"Error: The segmentation algorithm only takes 2D, 3D and 5D images. 2D images should be on the XY plane"<<endl;
        return;
    }

    return;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Defines the the mean and variances of a semiconjugate prior distribution over the actual model mean this->M
/// @param _M Pointer to the prior Gaussian Means.
/// @param _V Pointer to the prior Gaussian Variances.
///
/// Note that the data is copied internally. Thus, the user can delete the input pointers.
void seg_EM::SetMAP( float *_M, float* _V)
{
    this->mapStatus=true;
    this->MAP_M=new float[this->numb_classes];
    this->MAP_V=new float[this->numb_classes];


    for(int i=0; i<this->numb_classes; i++)
    {
        this->MAP_M[i]=_M[i];
        this->MAP_V[i]=_V[i];
    }

    return;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Defines the input anatomical priors.
/// @param _r The prior pointer. This pointer is owned by the caller.
///
/// Note that the nifti_image pointer is owned by the caller, but the called should not free the pointer while seg_EM is still using it. More specifically, seg_EM will use _r internally, but will not free it when clearing the object.
void seg_EM::SetPriorImage(nifti_image *_r)
{
    this->Priors = _r;
    this->priorsStatus = true;
    // Size
    this->dimentions=(int)((_r->nx)>1)+(int)((_r->ny)>1)+(int)((_r->nz)>1);
    if(this->nx==_r->nx && this->ny==_r->ny && this->nz==_r->nz)
    {
        return;
    }
    else
    {
        cout << "ERROR: Priors have wrong dimentions" << endl;
        return;
    }


}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Sets the filename of the output
/// @param f A char pointer to the output file name.
///

void seg_EM::SetFilenameOut(char *f)
{
    this->filenameOut = f;
    return;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Defines the ROI mask.
/// @param _r The mask pointer. This pointer is owned by the caller.
///
/// Note that the nifti_image pointer is owned by the caller, but the called should not free the pointer while seg_EM is still using it. More specifically, seg_EM will use _r internally, but will not free it when clearing the object. Also, note that if(Mask->datatype!=DT_BINARY), then the data in the input pointer will be modified, as it has to be converted to bool. Thus, the caller has to be carefull to delete the input, and has to be aware that the data was changed. This will be fixed in the near future.
void seg_EM::SetMaskImage(nifti_image *_r)
{
    this->Mask = _r;
    this->maskImageStatus = 1;

    if(Mask->datatype!=DT_BINARY)
    {
        seg_convert2binary(Mask,0.5f);
    }

    if(Mask->nt>1 ||Mask->nt>2)
    {
        cout << "ERROR: Mask has wrong dimentions" << endl;
        this->Mask->nt=1;
        this->Mask->nu=1;
        this->Mask->ndim=3;
    }

    bool * MaskDataPtr = static_cast<bool *>(Mask->data);
    this->numelmasked=0;
    for(int i=0; i<this->numel; i++,MaskDataPtr++)
    {
        if((*MaskDataPtr)>0)
        {
            this->numelmasked++;
        }
    }
    if(this->nx==_r->nx && this->ny==_r->ny && this->nz==_r->nz)
    {
        return;
    }
    else
    {
        cout << "ERROR: Mask has wrong dimentions" << endl;
        return;
    }
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Sets the filename of the output
/// @param verblevel Am unsigned int defining the verbose level.
///
/// In this definition of verbose, verblevel==0 means no verbose, verblevel==1 means normal verbose, verblevel==2 is 'almost' debug level verbose.
///
void seg_EM::SetVerbose(unsigned int verblevel)
{
    if(verblevel>2)
    {
        this->verbose_level=2;
    }
    else
    {
        this->verbose_level=verblevel;
    }
    return;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

void seg_EM::SetRegValue(float reg)
{
    if(reg>=1)
    {
        this->reg_factor=reg;
    }
    else
    {
        this->reg_factor=1;
        cout << "Non valid regularization value. Will assume -reg 1."<<endl;
    }
    return;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

void seg_EM::SetMaximalIterationNumber(unsigned int numberiter)
{
    this->maxIteration=numberiter;
    return;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

void seg_EM::SetMinIterationNumber(unsigned int numberiter)
{
    this->minIteration=numberiter;
    return;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

void seg_EM::SetAprox(bool aproxval)
{
    this->aprox=aproxval;
    return;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

void seg_EM::SetMRF(float strength)
{
    this->mrfStatus=true;
    this->mrfStrength=strength;
    this->MRFTransitionMatrix=new float[this->numb_classes*this->numb_classes]();
    CreateDiagonalMRFTransitionMatrix();


    if(this->maskImageStatus>0)
    {
        this->MRF=new float[this->numelmasked*this->numb_classes*this->nu*this->nt]();
        for(int i=0; i<(this->numelmasked*this->numb_classes); i++)
        {
            MRF[i]=(float)(1.0);
        }
    }
    else
    {
        this->MRF=new float[this->numel*this->numb_classes*this->nu*this->nt]();
        for(int i=0; i<(this->numel*this->numb_classes); i++)
        {
            MRF[i]=(float)(1.0);
        }
    }

    CreateDiagonalMRFTransitionMatrix();
    return;
}


void seg_EM::SetRelaxation(float relax_factor,float relax_gauss_kernel)
{
    this->relaxStatus=true;
    this->relaxFactor=relax_factor;
    this->relaxGaussKernelSize=relax_gauss_kernel;
    return;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

void seg_EM::SetBiasField(int _BiasField_order, float _BiasField_ratio)
{
    this->biasFieldStatus=true;
    this->biasFieldOrder=_BiasField_order;
    this->BiasField_coeficients = new segPrecisionTYPE[((_BiasField_order+1)*(_BiasField_order+2)/2*(_BiasField_order+3)/3)*this->nu*this->nt]();
    this->biasFieldRatio=_BiasField_ratio;
    if(this->maskImageStatus>0)
    {
        this->BiasField = new segPrecisionTYPE[this->numelmasked*this->nu*this->nt]();
    }
    else
    {
        this->BiasField = new segPrecisionTYPE[this->numel*this->nu*this->nt]();
    }
    return;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
void seg_EM::SetOutlierness(float _OutliernessThreshold, float _Outlierness_ratio)
{

    this->outliernessStatus=true;
    this->OutliernessThreshold=_OutliernessThreshold;
    this->Outlierness_ratio=_Outlierness_ratio;
    return;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

void seg_EM::InitializeAndAllocate()
{


    if(this->priorsStatus)
    {
        this->CreateExpectationAndShortPriors();
    }
    else
    {
        this->InitializeMeansUsingIntensity();
        int tmpnumb_elem=0;
        if(this->maskImageStatus>0)
        {
            tmpnumb_elem=(this->numelmasked*(this->numb_classes+(int)(this->pvModelStatus)*2));
        }
        else
        {
            tmpnumb_elem=(numel*(this->numb_classes+(int)(this->pvModelStatus)*2));
        }

        float tmpnumb_clas=((this->numb_classes+(int)(this->pvModelStatus)*2));
        this->Expec=new segPrecisionTYPE [tmpnumb_elem] ();
        this->ShortPrior=new segPrecisionTYPE [tmpnumb_elem] ();
        for(int i=0; i<tmpnumb_elem; i++)
        {
            this->Expec[i]=1.0/tmpnumb_clas;
            this->ShortPrior[i]=1.0/tmpnumb_clas;
        }
        this->RunExpectation();

    }

    for (int cl=0; cl<this->numb_classes; cl++)
    {
        if(this->outliernessStatus)
        {
            int tmpnumb_elem=0;
            if(this->maskImageStatus>0)
            {
                tmpnumb_elem=(this->numelmasked*(this->numb_classes+(int)(this->pvModelStatus)*2));
            }
            else
            {
                tmpnumb_elem=(numel*(this->numb_classes+(int)(this->pvModelStatus)*2));
            }
            this->OutliernessUSE=NULL;
            this->Outlierness=new segPrecisionTYPE [tmpnumb_elem] ();
            for(int i=0; i<tmpnumb_elem; i++)
            {
                this->Outlierness[i]=1.0;
            }
        }
    }

    return;

}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

void seg_EM::InitializeAndNormalizeImageAndPriors()
{
    if(this->maskImageStatus==0)
    {
        this->maskImageStatus=2;
        this->Mask=nifti_copy_nim_info(this->InputImage);
        Mask->dim[0]=3;
        Mask->dim[4]=1;
        Mask->dim[5]=1;
        Mask->datatype=DT_BINARY;
        Mask->cal_max=1;
        nifti_set_filenames(Mask,(const char *) "tmpInternalMask.nii.gz",0,0);
        nifti_update_dims_from_array(Mask);
        nifti_datatype_sizes(Mask->datatype,&Mask->nbyper,&Mask->swapsize);
        Mask->data = (void *) calloc(Mask->nvox, sizeof(bool));
        bool * Maskdata = static_cast<bool *>(Mask->data);
        for(unsigned int i=0; i<Mask->nvox; i++)
        {
            Maskdata[i]=1;
        }
        this->numelmasked=this->numel;
    }

    this->InitializeAndNormalizeNaNPriors();
    this->InitializeAndNormalizeImage();
    this->CreateShort2LongMatrix();
    this->CreateLong2ShortMatrix();


    if(this->mapStatus)
    {
        for(int i=0; i<this->numb_classes; i++)
        {
            this->MAP_M[i]=logf(((this->MAP_M[i]-this->rescale_min[0])/(this->rescale_max[0]-this->rescale_min[0]))+1)/0.693147181;;
            if(this->verbose_level>0)
            {
                cout << "MAP_M["<<i<<"] = "<<this->MAP_M[i]<< endl;
            }
            this->M[i]=this->MAP_M[i];
            this->V[i]=1.0/this->numb_classes;
        }
    }

    return;

}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

void seg_EM::InitializeAndNormalizeImage()
{

    if(this->verbose_level>0)
    {
        cout<< "Normalizing Input Image" << endl;
    }
    int numel=(int)(this->InputImage->nx*this->InputImage->ny*this->InputImage->nz);

    if(Mask->datatype!=DT_BINARY)
    {
        seg_convert2binary(Mask,0.0f);
    }
    if(this->InputImage->datatype!=NIFTI_TYPE_FLOAT32)
    {
        seg_changeDatatype<segPrecisionTYPE>(this->InputImage);
    }
    for(long udir=0; udir<this->nu; udir++) // Per Multispectral Image
    {
        bool * brainmaskptr = static_cast<bool *> (this->Mask->data);
        segPrecisionTYPE * Inputptrtmp = static_cast<segPrecisionTYPE *>(this->InputImage->data);
        segPrecisionTYPE * Inputptr=&Inputptrtmp[numel*udir];

        float tempmax=-(1e32);
        float tempmin=1e32;

        for (int i=0; i<numel; i++)
        {
            if(brainmaskptr[i])
            {
                if (Inputptr[i]<tempmin)
                {
                    tempmin=Inputptr[i];
                }
                if (Inputptr[i]>tempmax)
                {
                    tempmax=Inputptr[i];
                }
            }
        }
        this->rescale_max[udir]=tempmax;
        this->rescale_min[udir]=tempmin;
        if(this->verbose_level>0)
        {
            cout << "Normalization["<<udir<<"] = ["<<tempmin<<","<<tempmax<<"]"<<endl;
        }
        Inputptr=&Inputptrtmp[numel*udir];
        brainmaskptr = static_cast<bool *> (Mask->data);
        bool nanflag=false;
        for (int i=0; i<numel; i++)
        {
            Inputptr[i]=logf((((Inputptr[i])-tempmin)/(tempmax-tempmin))+1)/0.693147181;
            if(Inputptr[i]!=Inputptr[i])
            {
                if(nanflag==0){
                    cout<< "Warning: The image at timepoint="<<udir<<" has NaNs. This can cause problems." << endl;
                    nanflag=true;
                }
            }
        }
    }
    return;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

void seg_EM::InitializeAndNormalizeNaNPriors()
{
    if(this->priorsStatus){
        register int numel = Mask->nx*Mask->ny*Mask->nz;
        register int ups=0;
        register int good=0;
        if(this->verbose_level>0)
        {
            cout<< "Normalizing Priors" << endl;
        }
        if(Mask->datatype==DT_BINARY)
        {
            if(Priors->datatype==NIFTI_TYPE_FLOAT32)
            {
                segPrecisionTYPE * priorsptr = static_cast<segPrecisionTYPE *>(this->Priors->data);
                bool * brainmaskptr = static_cast<bool *> (this->Mask->data);

                for (int i=0; i<numel; i++)
                {
                    if(brainmaskptr[i])
                    {
                        float tempsum=0;
                        for (int j=0; j<this->Priors->nt; j++)
                        {
                            int tempind=i+numel*j;
                            if( priorsptr[tempind]<0.0 || priorsptr[tempind]!=priorsptr[tempind] || priorsptr[tempind]>1000 )
                            {
                                priorsptr[tempind]=0.0;
                            }
                            tempsum+=priorsptr[tempind];
                        }
                        if (tempsum>0 && tempsum<1000)
                        {
                            for (int j=0; j<Priors->nt; j++)
                            {
                                int tempind=i+numel*j;
                                priorsptr[tempind]=priorsptr[tempind]/tempsum;
                            }
                            good++;
                        }
                        else
                        {
                            for (int j=0; j<Priors->nt; j++)
                            {
                                int tempind=i+numel*j;
                                priorsptr[tempind]=1.0f/(Priors->nt);
                            }
                            ups++;

                        }
                    }
                    else
                    {
                        for (int j=0; j<Priors->nt; j++)
                        {
                            int tempind=i+numel*j;
                            priorsptr[tempind]=0;
                        }
                    }
                }
            }
            else
            {

                printf("Error:\tNormalize_NaN_Priors\tWrong Image datatype\n");

            }
        }
        else
        {

            printf("Error:\tNormalize_NaN_Priors\tWrong mask datatype\n");

        }

        if(this->verbose_level>0)
        {
            cout<<"Priors: "<< good<<" good voxels and "<<ups<<" bad voxels" << endl;
            flush(cout);
        }
    }
    return;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
void seg_EM::InitializeMeansUsingIntensity()
{

    for(int ms=0; ms<this->nu; ms++)
    {
        segPrecisionTYPE * Intensity_PTR = static_cast<segPrecisionTYPE *>(this->InputImage->data);
        bool * MaskDataPtr=NULL;
        if(this->maskImageStatus>0)
        {
            MaskDataPtr = static_cast<bool *>(this->Mask->data);
        }
        int mycounter=0;
        float meanval=0.0;
        float variance=0.0;

        for(int i=0; i<this->numel; i++)
        {
            if(this->maskImageStatus==0 || MaskDataPtr[i]>0)
            {
                mycounter++;
                meanval+=(Intensity_PTR[i+ms*this->numel]);
            }
        }
        meanval=meanval/mycounter;

        for(int i=0; i<this->numel; i++)
        {
            if(this->maskImageStatus==0 || MaskDataPtr[i]>0 )
            {
                variance+=pow((meanval-Intensity_PTR[i+ms*this->numel]),2);
            }
        }
        variance=variance/mycounter;
        int histogram[1001];
        for(int i=0; i<1000; i++)
        {
            histogram[i]=0;
        }
        float tmpmax=-1e32f;
        float tmpmin=1e32f;



        for(int i=0; i<this->numel; i++)
        {
            if(this->maskImageStatus==0 || MaskDataPtr[i]>0)
            {
                if(tmpmax<(int)(Intensity_PTR[i+ms*this->numel]))
                {
                    tmpmax=(int)(Intensity_PTR[i+ms*this->numel]);
                }
                if(tmpmin>(int)(Intensity_PTR[i+ms*this->numel]))
                {
                    tmpmin=(int)(Intensity_PTR[i+ms*this->numel]);
                }
            }
        }

        for(int i=0; i<this->numel; i++)
        {
            if(this->maskImageStatus==0 || MaskDataPtr[i]>0)
            {

                int index4hist=(int)(1000.0f*(float)(Intensity_PTR[i+ms*this->numel]-tmpmin)/(float)(tmpmax-tmpmin));
                if((index4hist>1000) & (index4hist<0))
                {
                    cout<< "error"<<endl;
                }
                histogram[(int)(1000.0*(float)(Intensity_PTR[i+ms*this->numel]-tmpmin)/(float)(tmpmax-tmpmin))]++;
            }
        }


        for(int clas=0; clas<this->numb_classes; clas++)
        {
            float tmpsum=0;
            int tmpindex=0;
            float percentile=((float)clas+1)/(this->numb_classes+1);
            for(int i=999; i>0; i--)
            {
                tmpsum+=histogram[i];
                tmpindex=i;
                if((float)(tmpsum)>((1.0f-percentile)*(float)(mycounter)))
                {
                    i=0;
                }
            }
            M[clas*this->nu+ms]=float(tmpindex)*(tmpmax-tmpmin)/1000.0f+(tmpmin);
            V[clas*this->nu*this->nu+ms*this->nu+ms]=variance/this->numb_classes/2;
        }
    }


    for (int cl=0; cl<this->numb_classes; cl++)
    {
        if(this->verbose_level>0)
        {
            if(this->nu==1)
            {
                cout.fill('0');
                cout<< "M["<<(int)(cl)<<"]= "<<setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(M[cl])<<"\tV["<<(int)(cl)<<"]="<<setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(V[cl])<< endl;
                flush(cout);
            }
            else
            {

                cout<< "M["<<(int)(cl)<<"]= ";
                for(int Multispec=0; Multispec<this->nu; Multispec++)
                {
                    cout<< setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(M[cl*this->nu+Multispec])<<"\t";
                }
                cout<< endl;
                flush(cout);
                cout<< "V["<<(int)(cl)<<"]= ";
                for(int Multispec=0; Multispec<this->nu; Multispec++)
                {
                    if(Multispec>0)
                    {
                        cout<< "      ";
                    }
                    for(int Multispec2=0; Multispec2<this->nu; Multispec2++)
                    {
                        cout<< setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(V[cl*this->nu*this->nu+Multispec*this->nu+Multispec2])<<"\t";
                    }
                    cout<< endl;
                }
                cout<< endl;
                flush(cout);
            }
        }
    }

    return;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

void seg_EM::CreateDiagonalMRFTransitionMatrix()
{
    for(int i=0; i<this->numb_classes; i++)
    {
        for(int j=0; j<this->numb_classes; j++)
        {
            if(j==i)
            {
                this->MRFTransitionMatrix[i+j*this->numb_classes]=0;
            }
            else
            {
                this->MRFTransitionMatrix[i+j*this->numb_classes]=this->mrfStrength;
            }
        }
    }
    return;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

void seg_EM::CreateCurrSizes()
{

    this->CurrSizes = new ImageSize [1]();
    CurrSizes->numel=(int)(this->nx*this->ny*this->nz);
    CurrSizes->xsize=this->nx;
    CurrSizes->ysize=this->ny;
    CurrSizes->zsize=this->nz;
    CurrSizes->usize=(this->nu>1)?this->nu:1;
    CurrSizes->tsize=(this->nt>1)?this->nt:1;
    CurrSizes->numclass=this->numb_classes;
    CurrSizes->numelmasked=this->numelmasked;
    CurrSizes->numelbias=0;

    return;

}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */


void seg_EM::CreateShort2LongMatrix()
{
    int numel_masked=0;
    int numel = this->Mask->nvox;
    if(this->Mask->datatype==DT_BINARY)
    {
        bool * Maskptr = static_cast<bool *> (this->Mask->data);
        bool * Maskptrtmp = Maskptr;
        for (int i=0; i<numel; i++, Maskptrtmp++)
        {
            (*Maskptrtmp)>0?numel_masked++:0;
        }
        this->numelmasked=numel_masked;

        this->S2L= new int [numel_masked]();
        int * Short_2_Long_Indices_PTR = (int *)(this->S2L);

        Maskptrtmp = Maskptr;
        int tempindex=0;
        for (int i=0; i<numel; i++)
        {
            if ((*Maskptrtmp)>0)
            {
                Short_2_Long_Indices_PTR[tempindex]=i;
                tempindex++;
            }
            Maskptrtmp++;
        }
        return;
    }
    else
    {
        printf("Error:\tCreate_Short_2_Long_Matrix\tWrong Mask datatype\n");
        this->S2L=NULL;
        return;

    }

}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

void seg_EM::CreateLong2ShortMatrix()
{
    int numel = this->Mask->nvox;
    this->L2S= new int [numel]();
    if(this->Mask->datatype==DT_BINARY)
    {
        bool * Maskptr = static_cast<bool *> (this->Mask->data);
        bool * Maskptrtmp = Maskptr;
        int * Long_2_Short_Indices_PTR = (int *) L2S;

        Maskptrtmp = Maskptr;
        int tempindex=0;
        for (int i=0; i<numel; i++,Maskptrtmp++,Long_2_Short_Indices_PTR++)
        {
            if ((*Maskptrtmp)>0)
            {
                (*Long_2_Short_Indices_PTR)=tempindex;
                tempindex++;
            }
            else
            {
                (*Long_2_Short_Indices_PTR)=-1;
            }
        }
        return;
    }
    else
    {
        cout<< "Error:\tCreate_Correspondace_Matrices\tWrong Mask datatype\n" << endl;
    }
    return;
}




/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
void seg_EM::CreateExpectationAndShortPriors()
{
    register long numel=(int)(this->Mask->nx*this->Mask->ny*this->Mask->nz);
    register long numel_masked=0;

    bool * Maskptrtmp = static_cast<bool *> (this->Mask->data);;
    for (long i=0; i<numel; i++, Maskptrtmp++)
    {
        *Maskptrtmp?numel_masked++:0;
    }
    int pluspv=(int)(this->pvModelStatus)*2;

    this->Expec = new segPrecisionTYPE [numel_masked*(this->numb_classes+pluspv)] ();
    segPrecisionTYPE * tempExpec= (segPrecisionTYPE *) Expec;
    this->ShortPrior = new segPrecisionTYPE [numel_masked*(this->numb_classes+pluspv)] ();
    segPrecisionTYPE * tempShortPrior= (segPrecisionTYPE *) ShortPrior;
    segPrecisionTYPE * PriorPTR = static_cast<segPrecisionTYPE *>(this->Priors->data);
    for(long cl=0; cl<this->numb_classes; cl++)
    {
        Maskptrtmp = static_cast<bool *> (this->Mask->data);;
        for (int i=numel; i--; Maskptrtmp++,PriorPTR++)
        {
            if(*Maskptrtmp)
            {
                *tempExpec = *PriorPTR;
                *tempShortPrior= *PriorPTR;
                tempExpec++;
                tempShortPrior++;
            }
        }
    }

    return;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

float * seg_EM::GetMeans()
{
    float * OutM= new float [this->nu*this->numb_classes];
    int index=0;
    for(int j=0; j<(this->numb_classes); j++)
    {
        for(int i=0; i<(this->nu); i++)
        {
            //cout << "M["<<j<<"]="<<this->M[j+i*this->numb_classes]<<endl;
            float resize=exp((this->M[index++])*0.693147181)-1;
            OutM[j+i*this->numb_classes]=(resize*(this->rescale_max[i]-this->rescale_min[i])+this->rescale_min[i]);
        }
    }

    return OutM;
}

float * seg_EM::GetSTD()
{
    float * OutV= new float [this->nu *this->numb_classes*this->numb_classes];
    for(int i=0; i<(this->nu); i++)
    {
        for(int j=0; j<(this->numb_classes*this->numb_classes); j++)
        {
            float resize=exp((sqrt(this->V[j+i*this->numb_classes*this->numb_classes]))*0.693147181)-1;
            OutV[j+i*this->numb_classes*this->numb_classes]=(resize*(this->rescale_max[i]-this->rescale_min[i]));
        }
    }
    return OutV;
}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

nifti_image * seg_EM::GetResult()
{

    //Result = Copy_Expec_to_Result_mask(this->Expec,this->S2L,this->InputImage,(char*)this->silenameOut.c_str(),this->CurrSizes);

    nifti_image * Result = nifti_copy_nim_info(this->InputImage);
    Result->dim[0]=4;
    Result->dim[4]=this->numb_classes;
    Result->dim[5]=1;
    Result->scl_inter=0;
    Result->scl_slope=1;
    Result->datatype=DT_FLOAT32;
    Result->cal_max=1;
    nifti_set_filenames(Result,(char*)this->filenameOut.c_str(),0,0);
    nifti_update_dims_from_array(Result);
    nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
    Result->data = (void *) calloc(Result->nvox, sizeof(segPrecisionTYPE));
    segPrecisionTYPE * Resultdata = static_cast<segPrecisionTYPE *>(Result->data);
    for(unsigned int i=0; i<Result->nvox; i++)
    {
        Resultdata[i]=0;
    }

    int class_nvox=Result->nx*Result->ny*Result->nz;
    for(long currclass=0; currclass<this->numb_classes; currclass++)
    {

        segPrecisionTYPE * Resultdata_class = &Resultdata[(currclass)*class_nvox];
        segPrecisionTYPE * Expec_PTR = &Expec[(currclass)*this->numelmasked];

        for(long i=0; i<(long)this->numelmasked; i++,Expec_PTR++)
        {
            Resultdata_class[this->S2L[i]]=*Expec_PTR;
        }
    }

    return Result;

}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

nifti_image * seg_EM::GetBiasCorrected(char * filename)
{

    int UsedBasisFunctions=(int)((this->biasFieldOrder+1) * (this->biasFieldOrder+2)/2 *(this->biasFieldOrder+3)/3);
    segPrecisionTYPE * InputImageData = static_cast<segPrecisionTYPE *>(this->InputImage->data);

    nifti_image * Result = nifti_copy_nim_info(this->InputImage);
    Result->dim[0]=4;
    Result->dim[4]=this->nu;
    Result->datatype=DT_FLOAT32;
    Result->cal_max=(this->rescale_max[0]);
    Result->scl_inter=0;
    Result->scl_slope=1;

    float * brainmask= new float [this->numel];
    for(long i=0; i<(long)this->numel; i++)
    {

        brainmask[i]=(InputImageData[i]!=InputImageData[i])?0:InputImageData[i];
    }
    otsu(brainmask,NULL,CurrSizes);
    Dillate(brainmask,5,CurrSizes);
    Erosion(brainmask,4,CurrSizes);
    bool* Maskptr = static_cast<bool * >(Mask->data);
    for(long i=0; i<(long)this->numel; i++)
    {
        brainmask[i]*=Maskptr[i];
    }
    GaussianFilter4D_cArray(brainmask, 3.0f, CurrSizes);


    nifti_set_filenames(Result,filename,0,0);
    nifti_update_dims_from_array(Result);
    nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
    Result->data = (void *) calloc(Result->nvox, sizeof(segPrecisionTYPE));
    segPrecisionTYPE * BiasCorrected_PTR = static_cast<segPrecisionTYPE *>(Result->data);

    float BiasField=0;
    segPrecisionTYPE currxpower[maxAllowedBCPowerOrder];
    segPrecisionTYPE currypower[maxAllowedBCPowerOrder];
    segPrecisionTYPE currzpower[maxAllowedBCPowerOrder];
    float xpos=0.0f;
    float ypos=0.0f;
    float zpos=0.0f;
    segPrecisionTYPE not_point_five_times_dims_x=(0.5f*(segPrecisionTYPE)this->nx);
    segPrecisionTYPE not_point_five_times_dims_y=(0.5f*(segPrecisionTYPE)this->ny);
    segPrecisionTYPE not_point_five_times_dims_z=(0.5f*(segPrecisionTYPE)this->nz);
    segPrecisionTYPE inv_not_point_five_times_dims_x=1.0f/(0.5f*(segPrecisionTYPE)this->nx);
    segPrecisionTYPE inv_not_point_five_times_dims_y=1.0f/(0.5f*(segPrecisionTYPE)this->ny);
    segPrecisionTYPE inv_not_point_five_times_dims_z=1.0f/(0.5f*(segPrecisionTYPE)this->nz);
    int ind=0;

    for(long multispec=0; multispec<this->nu; multispec++)
    {

        BiasCorrected_PTR = static_cast<segPrecisionTYPE *>(Result->data);
        BiasCorrected_PTR = &BiasCorrected_PTR[multispec*this->numel];
        InputImageData = static_cast<segPrecisionTYPE *>(this->InputImage->data);
        InputImageData = &InputImageData[multispec*this->numel];

        float * BiasFieldCoefs_multispec = &this->BiasField_coeficients[multispec*UsedBasisFunctions];


        for(long i=0; i<(long)this->numel; i++)
        {
            BiasCorrected_PTR[i]=0;
        }

        float to_resize=0;
        int index_full=0;
        for (int iz=0; iz<this->nz; iz++)
        {
            for (int iy=0; iy<this->ny; iy++)
            {
                for (int ix=0; ix<this->nz; ix++)
                {
                    BiasField=0.0f;
                    xpos=(((segPrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                    ypos=(((segPrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                    zpos=(((segPrecisionTYPE)iz-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z);

                    // Get the polynomial basis order
                    int order=1;
                    currxpower[0]=1;
                    currypower[0]=1;
                    currzpower[0]=1;
                    int orderminusone=0;
                    int maxorderplusone=this->biasFieldOrder+1;
                    while (order<maxorderplusone)
                    {
                        currxpower[order]=currxpower[orderminusone]*xpos;
                        currypower[order]=currypower[orderminusone]*ypos;
                        currzpower[order]=currzpower[orderminusone]*zpos;
                        order++;
                        orderminusone++;
                    }

                    // Estimate the basis
                    ind=0;
                    for(long order=0; order<=this->biasFieldOrder; order++)
                    {
                        for(long xorder=0; xorder<=order; xorder++)
                        {
                            for(long yorder=0; yorder<=(order-xorder); yorder++)
                            {
                                int zorder=order-yorder-xorder;
                                BiasField-=BiasFieldCoefs_multispec[ind]*currxpower[xorder]*currypower[yorder]*currzpower[zorder];
                                ind++;
                            }
                        }
                    }
                    BiasField*=brainmask[index_full];

                    to_resize=exp((BiasField+InputImageData[index_full])*0.693147181)-1;
                    BiasCorrected_PTR[index_full]=(to_resize*(this->rescale_max[multispec]-this->rescale_min[multispec])+this->rescale_min[multispec]);
                    index_full++;
                }
            }
        }
    }
    return Result;
}



/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

nifti_image * seg_EM::GetOutlierness(char * filename)
{
    nifti_image * Result = nifti_copy_nim_info(this->InputImage);
    Result->dim[0]=3;
    Result->dim[4]=1;
    Result->dim[5]=1;
    Result->datatype=DT_FLOAT32;
    Result->cal_max=1;
    nifti_set_filenames(Result,filename,0,0);
    nifti_update_dims_from_array(Result);
    nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
    Result->data = (void *) calloc(Result->nvox, sizeof(segPrecisionTYPE));
    segPrecisionTYPE * Resultdata = static_cast<segPrecisionTYPE *>(Result->data);
    for(unsigned int i=0; i<Result->nvox; i++)
    {
        Resultdata[i]=0;
    }



    if(this->maskImageStatus>0)
    {
        for(int i=0; i<this->numelmasked; i++)
        {
            float currsum=0;
            for(int currclass=0; currclass<this->numb_classes; currclass++)
            {
                currsum+=this->Outlierness[i+(currclass)*this->numelmasked]*Expec[i+(currclass)*this->numelmasked];
            }
            Resultdata[S2L[i]]=1-currsum;
        }

    }
    else
    {
        int class_nvox=Result->nx*Result->ny*Result->nz;
        for(int i=0; i<this->numel; i++)
        {
            float currsum=0;
            for(int currclass=0; currclass<this->numb_classes; currclass++)
            {
                currsum+=this->Outlierness[i+(currclass)*class_nvox];
            }
            Resultdata[i]=1-currsum;
        }
    }


    return Result;

}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */


void  seg_EM::Run_EM()
{
    time_t start,end;
    time(&start);
    if((int)(this->verbose_level)>(int)(0))
    {
        cout << "EM: Verbose level " << this->verbose_level << endl;
    }
    if(!priorsStatus)
    {

    }

    if(this->CurrSizes==NULL)
    {
        CreateCurrSizes();
    }

    InitializeAndNormalizeImageAndPriors();
    InitializeAndAllocate();
    if((int)(this->verbose_level)>(int)(0))
    {
        cout << "Number of voxels inside the mask = " << this->numelmasked << endl;
    }

    //**************
    // EM Algorithm
    //**************
    this->iter=0;
    bool out= true;

    while (out)
    {
        if(this->verbose_level>0)
        {
            cout << endl << "*******************************" << endl;
            cout << "Iteration " << iter << endl;
        }

        // Iterative Components - EM, MRF, Bias Correction

        //Maximization
        this->RunMaximization();
        //Expectation
        this->RunExpectation();
        //MRF
        this->RunMRF();
        //Bias Correction
        RunBiasField();
        //Update Weight
        this->RunPriorRelaxation();

        // Print LogLik depending on the verbose level
        if(this->verbose_level>0 && this->iter>0)
        {
            if(iter>0)
            {
                if ((this->loglik-this->oldloglik)/fabsf(this->oldloglik)>0 && (this->loglik-this->oldloglik)/fabsf(this->oldloglik)<100)
                {
                    cout<< "Loglik = " << setprecision(7)<<this->loglik << " : Ratio = " << (this->loglik-this->oldloglik)/fabsf(this->oldloglik) << endl;
                }
                else
                {
                    cout<< "Loglik = " << setprecision(7)<<this->loglik << endl;
                }
            }
            else
            {
                cout<< "Initial Loglik = " << setprecision(7) <<this->loglik << endl;
            }
        }
        // Preform Exit
        if((((this->loglik-this->oldloglik)/(fabs(this->oldloglik+this->loglik)/2.0f))<(segPrecisionTYPE)(0.0005) && this->iter>this->minIteration) || iter>=this->maxIteration || (isinf(this->loglik) && this->iter>3))
        {
            out=false;
        }
        this->ratio=((this->loglik-this->oldloglik)/fabs(this->oldloglik));
        // Update LogLik
        this->oldloglik=this->loglik;
        iter++;
    }

    time(&end);

    if(this->verbose_level>0)
    {
        int minutes = (int)floorf(float(end-start)/60.0f);
        int seconds = (int)(end-start - 60*minutes);
        cout << "Finished in "<<minutes<<"min "<<seconds<<"sec"<< endl;
    }
    return;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
void seg_EM::RunMaximization()
{

    int verbose=this->verbose_level;
    if(this->verbose_level>0)
    {
        cout<< "Optimising Gaussian Parameters" << endl;
        flush(cout);
    }
    bool OutliernessFlag=(Outlierness==NULL)?0:1;

    int numel_masked=this->numelmasked;
    int num_class=this->numb_classes;
    int Expec_offset[maxNumbClass];
    for (int cl=0; cl<num_class; cl++)
    {
        Expec_offset[cl]=cl*numel_masked;
    }

#ifdef _OPENMP
#pragma omp parallel for shared(T1,BiasField,Outlierness)
#endif
    // ***********
    for( int cl=0; cl<num_class; cl++)
    {
        int * S2L_PTR = (int *) this->S2L;
        segPrecisionTYPE * Expec_PTR = (segPrecisionTYPE *) this->Expec;
        segPrecisionTYPE * OutliernessPTR = (segPrecisionTYPE *) Outlierness;
        segPrecisionTYPE * T1_PTR = static_cast<segPrecisionTYPE *>(this->InputImage->data);
        segPrecisionTYPE * BiasField_PTR= (segPrecisionTYPE *) this->BiasField;
        segPrecisionTYPE tempsum= (segPrecisionTYPE) 0.0;
        segPrecisionTYPE SumPriors= (segPrecisionTYPE) 0.0;
        segPrecisionTYPE * T1_PTR2= static_cast<segPrecisionTYPE *>(this->InputImage->data);
        segPrecisionTYPE * BiasField_PTR2= (segPrecisionTYPE *) this->BiasField;

        // MEAN
        for(long Multispec=0; Multispec<this->nu; Multispec++)
        {
            Expec_PTR=(segPrecisionTYPE *) &Expec[Expec_offset[cl]];
            OutliernessPTR=(segPrecisionTYPE *) &Outlierness[Expec_offset[cl]];
            S2L_PTR = (int *) this->S2L;

            T1_PTR = static_cast<segPrecisionTYPE *>(this->InputImage->data);
            T1_PTR =&T1_PTR[Multispec*this->numel];
            tempsum=(segPrecisionTYPE)0.0;
            SumPriors=(segPrecisionTYPE)0.0;

            if(OutliernessFlag)
            {
                if(BiasField!=NULL)
                {
                    BiasField_PTR= &BiasField[Multispec*numel_masked];
                    for (int i=0; i<numel_masked; i++, Expec_PTR++,OutliernessPTR++,BiasField_PTR++,S2L_PTR++)
                    {
                        float current_value=(*Expec_PTR)*(*OutliernessPTR)*(T1_PTR[(*S2L_PTR)]+(*BiasField_PTR));
                        if(current_value==current_value)
                        {
                            tempsum+=current_value;
                            SumPriors+=(*Expec_PTR)*(*OutliernessPTR);
                        }
                    }
                }
                else
                {
                    for (int i=0; i<numel_masked; i++, Expec_PTR++,S2L_PTR++)
                    {
                        float current_value=(*Expec_PTR)*(*OutliernessPTR)*(T1_PTR[(*S2L_PTR)]);
                        if(current_value==current_value)
                        {
                            tempsum+=current_value;
                            SumPriors+=(*Expec_PTR)*(*OutliernessPTR);
                        }
                    }
                }
            }
            else
            {
                if(BiasField!=NULL)
                {
                    BiasField_PTR= &BiasField[Multispec*numel_masked];
                    for (int i=0; i<numel_masked; i++, Expec_PTR++,BiasField_PTR++,S2L_PTR++)
                    {
                        float current_value=(*Expec_PTR)*(T1_PTR[(*S2L_PTR)]+(*BiasField_PTR));
                        if(current_value==current_value)
                        {
                            tempsum+=current_value;
                            SumPriors+=(*Expec_PTR);
                        }
                    }
                }
                else
                {
                    for (int i=0; i<numel_masked; i++, Expec_PTR++,S2L_PTR++)
                    {
                        float current_value=(*Expec_PTR)*(T1_PTR[(*S2L_PTR)]);
                        if(current_value==current_value)
                        {
                            tempsum+=current_value;
                            SumPriors+=(*Expec_PTR);
                        }
                    }
                }

            }
            if(SumPriors==SumPriors && SumPriors>0)
            {
                if(this->mapStatus && this->MAP_M!=NULL)
                {
                    M[cl*(this->nu)+Multispec]=(tempsum/SumPriors/powf(V[cl*(this->nu)+Multispec],2)+this->MAP_M[cl*(this->nu)+Multispec]/powf(this->MAP_V[cl*(this->nu)+Multispec],2))/(1/powf(V[cl*(this->nu)+Multispec],2)+1/powf(this->MAP_V[cl*(this->nu)+Multispec],2));
                }
                else
                {
                    M[cl*(this->nu)+Multispec]=tempsum/SumPriors;

                }

                for(long Multispec2=Multispec; Multispec2<this->nu; Multispec2++)
                {
                    S2L_PTR = (int *) this->S2L;;

                    T1_PTR = static_cast<segPrecisionTYPE *>(this->InputImage->data);
                    T1_PTR =&T1_PTR[Multispec*this->numel];

                    T1_PTR2 = static_cast<segPrecisionTYPE *>(this->InputImage->data);
                    T1_PTR2 =&T1_PTR2[Multispec2*this->numel];
                    float tmpM=this->M[cl*this->nu+Multispec];
                    float tmpM2=this->M[cl*this->nu+Multispec2];
                    //STD
                    tempsum=0;
                    Expec_PTR=&Expec[Expec_offset[cl]];
                    OutliernessPTR=(segPrecisionTYPE *) &Outlierness[Expec_offset[cl]];
                    if(BiasField!=NULL)
                    {
                        BiasField_PTR=&BiasField[Multispec*numel_masked];
                        BiasField_PTR2=&BiasField[Multispec2*numel_masked];
                        if(OutliernessFlag)
                        {
                            for (int i=0; i<numel_masked; i++,Expec_PTR++,BiasField_PTR++,OutliernessPTR++,BiasField_PTR2++,S2L_PTR++)
                            {
                                float current_vaue=(*Expec_PTR) * (*OutliernessPTR)*(T1_PTR[(*S2L_PTR)]+(*BiasField_PTR)-tmpM) * (T1_PTR2[(*S2L_PTR)]+(*BiasField_PTR2)-tmpM2);
                                if(current_vaue==current_vaue)
                                {
                                    tempsum+=current_vaue;
                                }
                            }
                        }
                        else
                        {
                            for (int i=0; i<numel_masked; i++,Expec_PTR++,BiasField_PTR++,BiasField_PTR2++,S2L_PTR++)
                            {
                                float current_vaue=(*Expec_PTR) * (T1_PTR[(*S2L_PTR)]+(*BiasField_PTR)-tmpM) * (T1_PTR2[(*S2L_PTR)]+(*BiasField_PTR2)-tmpM2);
                                if(current_vaue==current_vaue)
                                {
                                    tempsum+=current_vaue;
                                }
                            }
                        }
                    }
                    else
                    {
                        for (int i=0; i<numel_masked; i++,Expec_PTR++,S2L_PTR++)
                        {
                            float current_vaue=(*Expec_PTR) * (T1_PTR[(*S2L_PTR)]-tmpM) * (T1_PTR2[(*S2L_PTR)]-tmpM2);
                            if(current_vaue==current_vaue)
                            {
                                tempsum+=current_vaue;
                            }
                        }

                    }
                    if( (tempsum/SumPriors>0) && SumPriors>0  && (!isnan(tempsum/SumPriors)))
                    {
                        V[cl*this->nu*this->nu+Multispec+Multispec2*this->nu]=tempsum/SumPriors;
                        if(Multispec2!=Multispec)
                        {
                            V[cl*this->nu*this->nu+Multispec2+Multispec*this->nu]=V[cl*this->nu*this->nu+Multispec+Multispec2*this->nu];
                            V[cl*this->nu*this->nu+Multispec+Multispec2*this->nu]/=reg_factor;
                            V[cl*this->nu*this->nu+Multispec2+Multispec*this->nu]/=reg_factor;
                        }
                    }
                }
            }
        }
    }

    for (int cl=0; cl<num_class; cl++)
    {
        if(verbose>0)
        {
            if(this->nu==1)
            {
                cout.fill('0');
                cout<< "M["<<(int)(cl)<<"]= "<<setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(M[cl])<<"\tV["<<(int)(cl)<<"]="<<setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(V[cl])<< endl;
                flush(cout);
            }
            else
            {

                cout<< "M["<<(int)(cl)<<"]= ";
                for(long Multispec=0; Multispec<this->nu; Multispec++)
                {
                    cout<< setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(M[cl*this->nu+Multispec])<<"\t";
                }
                cout<< endl;
                flush(cout);
                cout<< "V["<<(int)(cl)<<"]= ";
                for(long Multispec=0; Multispec<this->nu; Multispec++)
                {
                    if(Multispec>0)
                    {
                        cout<< "      ";
                    }
                    for(long Multispec2=0; Multispec2<this->nu; Multispec2++)
                    {
                        cout<< setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(V[cl*this->nu*this->nu+Multispec*this->nu+Multispec2])<<"\t";
                    }
                    cout<< endl;
                }
                cout<< endl;
                flush(cout);
            }
        }
    }
    return;
}



/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

void seg_EM::RunExpectation()
{


    if(this->ratio<(segPrecisionTYPE)(this->Outlierness_ratio) && this->outliernessStatus)
    {
        this->OutliernessUSE=this->Outlierness;
        if(this->verbose_level>0)
        {
            cout << "Updating Outlierness - LogRatio = "<<ratio<<endl;
        }
    }

    segPrecisionTYPE * IterPrior=NULL;
    if(this->mrfStatus)
    {
        IterPrior=this->MRF;
    }
    else
    {
        IterPrior=this->ShortPrior;
    }

    int numel_masked=this->numelmasked;
    int num_class=this->numb_classes;
    bool OutliernessFlag=(this->OutliernessUSE==NULL)?0:1;
    segPrecisionTYPE inv_v [maxNumbClass*maxMultispectalSize*maxMultispectalSize]= {0.0f};
    segPrecisionTYPE inv_sqrt_V_2pi [maxNumbClass]= {0.0f};

    int Expec_offset [maxNumbClass]= {0};

    for (int cl=0; cl<num_class; cl++)
    {
        Expec_offset[cl]=(int) cl*numel_masked;
        if(this->nu>1)
        {
            seg_Matrix <double> Vmat(this->nu,this->nu);

            for(long j2=0; j2<this->nu; j2++)
            {
                for(long i2=j2; i2<this->nu; i2++)
                {
                    Vmat.setvalue(i2,j2,(double)(this->V[i2+j2*this->nu+cl*this->nu*this->nu]));
                    Vmat.setvalue(j2,i2,(double)(this->V[i2+j2*this->nu+cl*this->nu*this->nu]));
                }
            }
            inv_sqrt_V_2pi[cl]=1/(sqrtf(2*M_PI*Vmat.determinant()));
            if(this->verbose_level>1)
            {
                cout<<endl<<"inv_sqrt_V_2pi["<< cl <<"]= "<< inv_sqrt_V_2pi[cl] << endl;
                flush(cout);
            }
            Vmat.invert();
            double cvalue=0.0f;
            bool success;
            if(this->verbose_level>1)
            {
                cout<<"inv_V["<< cl <<"]= ";
                flush(cout);
            }
            for(long j2=0; j2<this->nu; j2++)
            {
                if(this->verbose_level>1)
                {
                    if(j2!=0)
                    {
                        cout<< endl << "          ";
                    }
                }
                for(long i2=0; i2<this->nu; i2++)
                {
                    Vmat.getvalue(i2,j2,cvalue,success);
                    inv_v[i2+j2*this->nu+cl*this->nu*this->nu]=(segPrecisionTYPE)(cvalue);
                    if(this->verbose_level>1)
                    {
                        cout<<inv_v[i2+j2*this->nu+cl*this->nu*this->nu]<< "\t";
                        flush(cout);
                    }
                }

            }
            if(this->verbose_level>1)
            {
                cout<< endl;
            }
        }
        else
        {
            inv_sqrt_V_2pi[cl]=1/(sqrtf(2*M_PI*V[cl]));
            inv_v[cl]=1/V[cl];
        }
    }
    this->loglik=0;

    float logliktmp=0.0f;


#ifdef _OPENMP
    float * loglikthread = new float [omp_get_max_threads()]();
    for(long i=0; i<(long)omp_get_max_threads(); i++)
        loglikthread[i]=0;

#pragma omp parallel for shared(Expec,loglikthread,T1,BiasField,Outlierness,IterPrior)
#endif
    for (int i=0; i<numel_masked; i++)
    {
        segPrecisionTYPE * T1_PTR = static_cast<segPrecisionTYPE *>(this->InputImage->data);
        segPrecisionTYPE T1_Bias_corr[maxMultispectalSize];
        segPrecisionTYPE SumExpec=0.0f;

        for(long Multispec=0; Multispec<this->nu; Multispec++)
            T1_Bias_corr[Multispec]=(BiasField!=NULL)?(T1_PTR[this->S2L[i]+Multispec*this->numel] + BiasField[i+Multispec*numel_masked]):(T1_PTR[this->S2L[i]+Multispec*this->numel]);

        //Expec_offset_PTR=Expec_offset;

        for (int cl=0; cl<num_class; cl++)
        {
            segPrecisionTYPE mahal=0.0f;
            for(long Multispec=0; Multispec<this->nu; Multispec++)
            {
                segPrecisionTYPE tmpT1_BC_minusM=(T1_Bias_corr[Multispec] - this->M[cl*(this->nu)+Multispec]);
                for(long Multispec2=0; Multispec2<this->nu; Multispec2++)
                {
                    mahal-=(0.5f)*(T1_Bias_corr[Multispec2] - M[cl*(this->nu)+Multispec2])*inv_v[cl*this->nu*this->nu+Multispec+Multispec2*this->nu]*tmpT1_BC_minusM;
                }
            }

            if(OutliernessFlag)
            {
                float outvalue=(expf(mahal)+0.01)/(expf(mahal)+expf(-0.5*(this->OutliernessThreshold*this->OutliernessThreshold))+0.01);
                Outlierness[i+Expec_offset[cl]]=outvalue;
            }
            Expec[i+Expec_offset[cl]]=IterPrior[i+Expec_offset[cl]] * expf(mahal) * inv_sqrt_V_2pi[cl];
            SumExpec+=Expec[i+Expec_offset[cl]];
        }

        if (SumExpec<=0.0 || SumExpec!=SumExpec)
        {
            for (int cl=0; cl<num_class; cl++)
            {
                Expec[i+Expec_offset[cl]]=(float)(1)/(float)(num_class);
            }

        }
        else
        {

            for (int cl=0; cl<num_class; cl++)
            {
                Expec[i+Expec_offset[cl]]=Expec[i+Expec_offset[cl]]/SumExpec;
            }
#ifdef _OPENMP
            loglikthread[omp_get_thread_num()]+=logf(SumExpec);
#else
            logliktmp+=logf(SumExpec);
#endif
        }
    }

#ifdef _OPENMP
    for(long i =0; i<(long)omp_get_max_threads(); i++)
        logliktmp+=loglikthread[i];
#endif

    loglik=logliktmp;
    return;

}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

void seg_EM::RunMRF()
{
    if(this->nz>1)
    {

        RunMRF3D();

    }
    else if(this->nz==1)
    {
        RunMRF2D();
    }
    return;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

void seg_EM::RunMRF2D()
{

    int numelmasked=this->numelmasked;
    int numclass=this->numb_classes;
    segPrecisionTYPE * G =this->MRFTransitionMatrix;

    if(this->mrfStatus)
    {
        segPrecisionTYPE * MRFpriorPtr = (segPrecisionTYPE *)this->MRF;
        int * Long_2_Short_IndicesPtr = (int *)this->L2S;
        int col_size, indexCentre, indexWest, indexEast, indexSouth, indexNorth;
        int ix, iy,maxiy, maxix, neighbourclass;
        segPrecisionTYPE Sum_Temp_MRF_Class_Expect;
        col_size = (int)(this->nx);
        maxix = (int)(this->nx);
        maxiy = (int)(this->ny);
        segPrecisionTYPE Gplane[maxNumbClass];
        segPrecisionTYPE Temp_MRF_Class_Expect[maxNumbClass];
        if(verbose_level>0)
        {
            cout << "Optimising MRF"<<endl;
            flush(cout);
        }
        register int currclass;

        unsigned int numelmasked_currclass_shift[maxNumbClass];
        for(int i=0; i<numclass; i++)
        {
            numelmasked_currclass_shift[i]=i*numelmasked;
        }
        int curr_short_centreindex;
        for (iy=1; iy<maxiy-1; iy++)
        {
            indexCentre=(col_size*iy);
            for (ix=1; ix<maxix-1; ix++)
            {
                indexCentre++;
                Sum_Temp_MRF_Class_Expect = 0;
                curr_short_centreindex=Long_2_Short_IndicesPtr[indexCentre];
                if (curr_short_centreindex>=0)
                {
                    indexWest=Long_2_Short_IndicesPtr[indexCentre-col_size]>-1?Long_2_Short_IndicesPtr[indexCentre-col_size]:0;
                    indexEast=Long_2_Short_IndicesPtr[indexCentre+col_size]>-1?Long_2_Short_IndicesPtr[indexCentre+col_size]:0;
                    indexNorth=Long_2_Short_IndicesPtr[indexCentre-1]>-1?Long_2_Short_IndicesPtr[indexCentre-1]:0;
                    indexSouth=Long_2_Short_IndicesPtr[indexCentre+1]>-1?Long_2_Short_IndicesPtr[indexCentre+1]:0;
                    for (currclass=0; currclass<numclass; currclass++)
                    {
                        Gplane[currclass] = 0.0;
                        Temp_MRF_Class_Expect[currclass] = 0.0;
                        Gplane[currclass]+=this->Expec[indexWest];
                        Gplane[currclass]+=this->Expec[indexEast];
                        Gplane[currclass]+=this->Expec[indexNorth];
                        Gplane[currclass]+=this->Expec[indexSouth];
                        if(currclass<numclass)
                        {
                            indexWest+=numelmasked;
                            indexEast+=numelmasked;
                            indexNorth+=numelmasked;
                            indexSouth+=numelmasked;
                        }
                    }
                    for (currclass=0; currclass<numclass; currclass++)
                    {
                        for (neighbourclass=0; neighbourclass<numclass; neighbourclass++)
                        {
                            Temp_MRF_Class_Expect[currclass]-=G[currclass+(numclass)*neighbourclass]*Gplane[neighbourclass];
                        }
                        if(this->MRFBeta==NULL)
                        {
                            Temp_MRF_Class_Expect[currclass] = exp(Temp_MRF_Class_Expect[currclass])*this->ShortPrior[curr_short_centreindex+numelmasked_currclass_shift[currclass]];
                        }
                        else
                        {
                            Temp_MRF_Class_Expect[currclass] = exp(this->MRFBeta[curr_short_centreindex]*Temp_MRF_Class_Expect[currclass])*this->ShortPrior[curr_short_centreindex+numelmasked_currclass_shift[currclass]];
                        }
                        Sum_Temp_MRF_Class_Expect += Temp_MRF_Class_Expect[currclass];
                    }
                    for (currclass=0; currclass<numclass; currclass++)
                    {
                        MRFpriorPtr[curr_short_centreindex+numelmasked_currclass_shift[currclass]]=(Temp_MRF_Class_Expect[currclass]/Sum_Temp_MRF_Class_Expect);
                    }
                }
            }
        }

    }
    else
    {
        for(int currclass=0; currclass<numclass; currclass++)
        {
            for(int i=0; i<(numelmasked); i++)
            {
                this->MRF[i+currclass*numelmasked]=this->ShortPrior[i+currclass*numelmasked];
            }
        }
    }
    return;

}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

void seg_EM::RunMRF3D()
{
    int numelmasked=this->numelmasked;
    int numclass=this->numb_classes;

    segPrecisionTYPE * G =this->MRFTransitionMatrix;
    segPrecisionTYPE * H =this->MRFTransitionMatrix;

    if(this->mrfStatus)
    {
        segPrecisionTYPE * MRFpriorPtr = (segPrecisionTYPE *)this->MRF;
        int * Long_2_Short_IndicesPtr = (int *)this->L2S;
        int col_size, plane_size;
        int maxiy, maxix, maxiz;
        col_size = (int)(this->nx);
        plane_size = (int)(this->nx)*(this->ny);

        maxix = (int)(this->nx);
        maxiy = (int)(this->ny);
        maxiz = (int)(this->nz);

        if(verbose_level>0)
        {
            cout << "Optimising MRF"<<endl;
            flush(cout);
        }

        unsigned int numelmasked_currclass_shift[maxNumbClass];
        //unsigned int image_size_currclass_shift[max_numbclass];
        for(int i=0; i<numclass; i++)
        {
            numelmasked_currclass_shift[i]=i*numelmasked;

        }

#ifdef _OPENMP
#pragma omp parallel for shared(Expec,MRFprior,Long_2_Short_IndicesPtr)
#endif
        for (int iz=1; iz<maxiz-1; iz++)
        {
            segPrecisionTYPE Sum_Temp_MRF_Class_Expect;
            register int currclass;
            segPrecisionTYPE Temp_MRF_Class_Expect[maxNumbClass];
            segPrecisionTYPE Gplane[maxNumbClass];
            segPrecisionTYPE Hplane[maxNumbClass];
            int indexCentre, indexWest, indexEast, indexSouth, indexNorth, indexTop, indexBottom;

            for (int iy=1; iy<maxiy-1; iy++)
            {
                indexCentre=(col_size*iy)+(plane_size*iz);
                for (int ix=1; ix<maxix-1; ix++)
                {
                    indexCentre++;
                    Sum_Temp_MRF_Class_Expect = 0;
                    int curr_short_centreindex=Long_2_Short_IndicesPtr[indexCentre];
                    if (curr_short_centreindex>=0)
                    {
                        indexWest=Long_2_Short_IndicesPtr[indexCentre-col_size]>-1?Long_2_Short_IndicesPtr[indexCentre-col_size]:0;
                        indexEast=Long_2_Short_IndicesPtr[indexCentre+col_size]>-1?Long_2_Short_IndicesPtr[indexCentre+col_size]:0;
                        indexNorth=Long_2_Short_IndicesPtr[indexCentre-1]>-1?Long_2_Short_IndicesPtr[indexCentre-1]:0;
                        indexSouth=Long_2_Short_IndicesPtr[indexCentre+1]>-1?Long_2_Short_IndicesPtr[indexCentre+1]:0;
                        indexBottom=Long_2_Short_IndicesPtr[indexCentre+plane_size]>-1?Long_2_Short_IndicesPtr[indexCentre+plane_size]:0;
                        indexTop=Long_2_Short_IndicesPtr[indexCentre-plane_size]>-1?Long_2_Short_IndicesPtr[indexCentre-plane_size]:0;
                        for (currclass=0; currclass<numclass; currclass++)
                        {
                            Gplane[currclass] = 0.0;
                            Hplane[currclass] = 0.0;
                            Temp_MRF_Class_Expect[currclass] = 0.0;
                            Gplane[currclass]+=this->Expec[indexWest];
                            Gplane[currclass]+=this->Expec[indexEast];
                            Gplane[currclass]+=this->Expec[indexNorth];
                            Gplane[currclass]+=this->Expec[indexSouth];
                            Hplane[currclass]+=this->Expec[indexTop];
                            Hplane[currclass]+=this->Expec[indexBottom];
                            if(currclass<numclass)
                            {
                                indexWest+=numelmasked;
                                indexEast+=numelmasked;
                                indexNorth+=numelmasked;
                                indexSouth+=numelmasked;
                                indexTop+=numelmasked;
                                indexBottom+=numelmasked;
                            }
                        }
                        for (currclass=0; currclass<numclass; currclass++)
                        {
                            for (int neighbourclass=0; neighbourclass<numclass; neighbourclass++)
                            {
                                Temp_MRF_Class_Expect[currclass]-=G[currclass+(numclass)*neighbourclass]*Gplane[neighbourclass]+H[currclass+(numclass)*neighbourclass]*Hplane[neighbourclass];
                            }
                            if(this->MRFBeta==NULL)
                            {
                                Temp_MRF_Class_Expect[currclass] = exp(Temp_MRF_Class_Expect[currclass])*this->ShortPrior[curr_short_centreindex+numelmasked_currclass_shift[currclass]];
                            }
                            else
                            {
                                Temp_MRF_Class_Expect[currclass] = exp(this->MRFBeta[curr_short_centreindex]*Temp_MRF_Class_Expect[currclass])*this->ShortPrior[curr_short_centreindex+numelmasked_currclass_shift[currclass]];
                            }
                            Sum_Temp_MRF_Class_Expect+=Temp_MRF_Class_Expect[currclass];
                        }
                        for (currclass=0; currclass<numclass; currclass++)
                        {
                            MRFpriorPtr[curr_short_centreindex+numelmasked_currclass_shift[currclass]]=(Temp_MRF_Class_Expect[currclass]/Sum_Temp_MRF_Class_Expect);
                        }
                    }
                }
            }
        }
    }
    else
    {
        for(int currclass=0; currclass<numclass; currclass++)
        {
            for(int i=0; i<(numelmasked); i++)
            {
                this->MRF[i+currclass*numelmasked]=this->ShortPrior[i+currclass*numelmasked];
            }
        }
    }
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

void seg_EM::RunBiasField()
{

    if(this->biasFieldStatus && ((((this->loglik-this->oldloglik)/fabs(this->oldloglik))<(this->biasFieldRatio)
                                  && this->iter>3)||((this->biasFieldRatio)==0.0f)))
    {
        if(this->nz>1)
        {
            this->RunBiasField3D();
        }
        else
        {
            this->RunBiasField2D();
        }
    }
    return;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

void seg_EM::RunBiasField3D()
{

    if(!this->biasFieldStatus)
    {
        return;
    }
    if(verbose_level>0)
    {
        cout << "Optimising the Bias Field with order " << this->biasFieldOrder<< endl;
        if(this->nu>1)
        {
            cout<< "Assuming fully decoupled bias-fields" << endl;
        }
        flush(cout);
    }
    int reduxfactor=reduxFactorForBias;
    int nrOfClasses = this->numb_classes;
    int TotalLength = this->numelmasked;
    int UsedBasisFunctions=(int)((this->biasFieldOrder+1) * (this->biasFieldOrder+2)/2 *(this->biasFieldOrder+3)/3);
    segPrecisionTYPE * sampledData = static_cast<segPrecisionTYPE *>(this->InputImage->data);


    // Precompute Powers depending on the current BiasOrder
    int PowerOrder [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3))]= {0};
    int ind=0;
    for(int order=0; order<=this->biasFieldOrder; order++)
    {
        for(int xorder=0; xorder<=order; xorder++)
        {
            for(int yorder=0; yorder<=(order-xorder); yorder++)
            {
                int zorder=order-yorder-xorder;
                PowerOrder[ind] =xorder;
                PowerOrder[ind+1] =yorder;
                PowerOrder[ind+2] =zorder;
                ind += 3;
            }
        }
    }
    segPrecisionTYPE invV[maxNumbClass];
    segPrecisionTYPE currM[maxNumbClass];

    for(long multispec=0; multispec<this->nu; multispec++)
    {
        sampledData = static_cast<segPrecisionTYPE *>(this->InputImage->data);
        sampledData = &sampledData[multispec*this->numel];
        // Precompute the M and V  inverses
        for(int i=0; i<nrOfClasses; i++)
        {
            invV[i]=1.0f/(V[i*this->nu*this->nu+multispec+multispec*this->nu]);
            currM[i]=M[i*this->nu+multispec];
        }
        segPrecisionTYPE A [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)*((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)]= {0.0f};
        segPrecisionTYPE B [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)]= {0.0f};
        segPrecisionTYPE C [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)]= {0.0f};


        // Precompute sizes
        int col_size = (int)(this->nx);
        int plane_size = (int)(this->nx)*(this->ny);
        int maxix = (int)(this->nx);
        int maxiy = (int)(this->ny);
        int maxiz = (int)(this->nz);
        int Dims3d[3]= {0};
        Dims3d[0]=maxix;
        Dims3d[1]=maxiy;
        Dims3d[2]=maxiz;

        // Precompute number of samples as it was never computed
        int samplecount=0;
        int linearindexes=0;
        int currshortindex=0;
        if(this->numelbias==0)
        {
            for (int iz=0; iz<maxiz; iz+=reduxfactor)
            {
                for (int iy=0; iy<maxiy; iy+=reduxfactor)
                {
                    for (int ix=0; ix<maxix; ix+=reduxfactor)
                    {
                        linearindexes=iz*plane_size+iy*col_size+ix;
                        if(this->L2S[linearindexes]>=0)
                        {
                            samplecount++;
                        }
                    }
                }
            }
            if(verbose_level>0)
            {
                cout << "Number of samples for BiasField = " << samplecount<<"\n";
                flush(cout);
            }
            this->numelbias=samplecount;
        }
        else
        {
            samplecount=this->numelbias;
        }



        // CALC MATRIX A

        // Calc W (Van Leemput 1999 eq 7)
        //cout << "Calculating C = inv(A'WA) WR";
        //flush(cout);

        segPrecisionTYPE * Tempvar= new segPrecisionTYPE [samplecount] ();
        segPrecisionTYPE Tempvar_tmp=0;
        currshortindex=0;
        int tempvarindex=0;
        for (int iz=0; iz<maxiz; iz+=reduxfactor)
        {
            for (int iy=0; iy<maxiy; iy+=reduxfactor)
            {
                for (int ix=0; ix<maxix; ix+=reduxfactor)
                {
                    currshortindex=this->L2S[iz*plane_size+iy*col_size+ix];
                    if(currshortindex>=0)
                    {
                        Tempvar_tmp=0;
                        if(Outlierness==NULL)
                        {
                            for(int j=0; j<nrOfClasses; j++)
                            {
                                Tempvar_tmp+=Expec[currshortindex+TotalLength*j]*invV[j];
                            }
                        }
                        else
                        {
                            for(int j=0; j<nrOfClasses; j++)
                            {
                                Tempvar_tmp+=Expec[currshortindex+TotalLength*j]*Outlierness[currshortindex+TotalLength*j]*invV[j];
                            }

                        }
                        Tempvar[tempvarindex]=Tempvar_tmp;
                        tempvarindex++;
                    }
                }
            }
        }

        // Precompute shifts
        segPrecisionTYPE not_point_five_times_dims_x=(0.5f*(segPrecisionTYPE)Dims3d[0]);
        segPrecisionTYPE not_point_five_times_dims_y=(0.5f*(segPrecisionTYPE)Dims3d[1]);
        segPrecisionTYPE not_point_five_times_dims_z=(0.5f*(segPrecisionTYPE)Dims3d[2]);

        segPrecisionTYPE inv_not_point_five_times_dims_x=1.0f/(0.5f*(segPrecisionTYPE)Dims3d[0]);
        segPrecisionTYPE inv_not_point_five_times_dims_y=1.0f/(0.5f*(segPrecisionTYPE)Dims3d[1]);
        segPrecisionTYPE inv_not_point_five_times_dims_z=1.0f/(0.5f*(segPrecisionTYPE)Dims3d[2]);

        segPrecisionTYPE * Basis= new segPrecisionTYPE[UsedBasisFunctions]();
        segPrecisionTYPE xpos=0.0f;
        segPrecisionTYPE ypos=0.0f;
        segPrecisionTYPE zpos=0.0f;
        int x_bias_index_shift=0;
        int y_bias_index_shift=1;
        int z_bias_index_shift=2;
        segPrecisionTYPE current_Tempvar=0.0f;

        // Calc A'WA (Van Leemput 1999 eq 7)
        tempvarindex=0;
        segPrecisionTYPE * Basisptr1= (segPrecisionTYPE *) Basis;
        segPrecisionTYPE * Basisptr2= (segPrecisionTYPE *) Basis;
        segPrecisionTYPE * Aptr= (segPrecisionTYPE *) A;
        for (int iz=0; iz<maxiz; iz+=reduxfactor)
        {
            for (int iy=0; iy<maxiy; iy+=reduxfactor)
            {
                for (int ix=0; ix<maxix; ix+=reduxfactor)
                {
                    linearindexes=(iz)*(this->nx)*(this->ny)+(iy)*(this->nx)+ix;
                    currshortindex=this->L2S[linearindexes];
                    if(currshortindex>=0)
                    {
                        Basisptr1= (segPrecisionTYPE *) Basis;
                        current_Tempvar=Tempvar[tempvarindex];
                        xpos=(((segPrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                        ypos=(((segPrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                        zpos=(((segPrecisionTYPE)iz-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z);
                        x_bias_index_shift=0;
                        y_bias_index_shift=1;
                        z_bias_index_shift=2;
                        for(int j2=0; j2<UsedBasisFunctions; j2++,x_bias_index_shift+=3,y_bias_index_shift+=3,z_bias_index_shift+=3, Basisptr1++)
                        {
                            // Because Powerorder is always int, use a special power function (faster)
                            *Basisptr1=(pow_int(xpos,PowerOrder[x_bias_index_shift])*pow_int(ypos,PowerOrder[y_bias_index_shift])*pow_int(zpos,PowerOrder[z_bias_index_shift]));
                        }


                        Basisptr1= (segPrecisionTYPE *) Basis;
                        Aptr= (segPrecisionTYPE *) A;
                        for(int j2=0; j2<UsedBasisFunctions; j2++, Basisptr1++)
                        {
                            Basisptr2= &Basis[j2];
                            Aptr= &A[j2+j2*UsedBasisFunctions];
                            for(int i2=j2; i2<UsedBasisFunctions; i2++, Aptr++, Basisptr2++)
                            {
                                (*Aptr)+=(*Basisptr2)*(current_Tempvar)*(*Basisptr1);
                            }
                        }
                        tempvarindex++;
                    }
                }
            }
        }



        seg_Matrix <double> RealA(UsedBasisFunctions,UsedBasisFunctions);

        for(int j2=0; j2<UsedBasisFunctions; j2++)
        {
            for(int i2=j2; i2<UsedBasisFunctions; i2++)
            {
                RealA.setvalue(i2,j2,(double)(A[i2+j2*UsedBasisFunctions]));
                RealA.setvalue(j2,i2,(double)(A[i2+j2*UsedBasisFunctions]));
            }
        }

        seg_Matrix <double> RealA_inv(UsedBasisFunctions);
        RealA_inv.copymatrix(RealA);
        RealA_inv.invert();

        if(verbose_level>1)
        {
            seg_Matrix <double> RealA_test(UsedBasisFunctions);
            RealA_test.settoproduct(RealA,RealA_inv);
            RealA_test.comparetoidentity();
            //RealA.dumpmatrix();
        }



        // CALC MATRIX B

        //Precompute WR (Van Leemput 1999 eq 7)
        segPrecisionTYPE Wi;
        segPrecisionTYPE Wij;
        segPrecisionTYPE Yest;
        segPrecisionTYPE Ysum;
        tempvarindex=0;
        for (int iz=0; iz<maxiz; iz+=reduxfactor)
        {
            for (int iy=0; iy<maxiy; iy+=reduxfactor)
            {
                for (int ix=0; ix<maxix; ix+=reduxfactor)
                {
                    linearindexes=(iz)*(this->nx)*(this->ny)+(iy)*(this->nx)+ix;
                    currshortindex=this->L2S[linearindexes];
                    if(currshortindex>=0)
                    {
                        Wi=0;
                        Wij=0;
                        Yest=0;
                        Ysum=0;
                        for(int j=0; j<nrOfClasses; j++)
                        {
                            segPrecisionTYPE tmpexpec = (segPrecisionTYPE)Expec[currshortindex+TotalLength*j];
                            Wij=tmpexpec*(invV[j]);
                            Wi+=Wij;
                            Yest+=Wij*(currM[j]);
                            Ysum+=Wij;
                        }
                        Tempvar[tempvarindex]=Wi*(sampledData[linearindexes]-(Yest/Ysum));
                        tempvarindex++;
                    }
                }
            }
        }

        //#ifdef _OPENMP
        //#pragma omp parallel shared(Tempvar,PowerOrder,CurrSizes,B,Long_2_Short_Indices) private(UsedBasisFunctions,maxiz,maxiy,maxix, not_point_five_times_dims_x, not_point_five_times_dims_y, not_point_five_times_dims_z,reduxfactor,inv_not_point_five_times_dims_x, inv_not_point_five_times_dims_y, inv_not_point_five_times_dims_z)
        //#endif
        for(int bfindex=0; bfindex<UsedBasisFunctions; bfindex++)
        {
            int tempvarindex2=0;
            int linearindexes2=0;
            segPrecisionTYPE * BPTR=&B[bfindex];
            *BPTR=0;
            for (int iz2=0; iz2<maxiz; iz2+=reduxfactor)
            {
                for (int iy2=0; iy2<maxiy; iy2+=reduxfactor)
                {
                    for (int ix2=0; ix2<maxix; ix2+=reduxfactor)
                    {
                        linearindexes2=(iz2)*(this->nx)*(this->ny)+(iy2)*(this->nx)+ix2;
                        if(this->L2S[linearindexes2]>=0)
                        {
                            *BPTR+=pow_int((((segPrecisionTYPE)ix2-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x),PowerOrder[0+bfindex*3])*
                                    pow_int((((segPrecisionTYPE)iy2-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y),PowerOrder[1+bfindex*3])*
                                    pow_int((((segPrecisionTYPE)iz2-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z),PowerOrder[2+bfindex*3])*
                                    Tempvar[tempvarindex2];


                            tempvarindex2++;
                        }
                    }
                }
            }
            if(*BPTR!=*BPTR)
            {
                *BPTR=1;
            }
        }

        //#ifdef _OPENMP
        //#pragma omp barrier
        //#endif
        seg_Matrix <double> RealB(UsedBasisFunctions,1);

        for(int i2=0; i2<UsedBasisFunctions; i2++)
        {
            RealB.setvalue(i2,0,(double)(B[i2]));
        }

        seg_Matrix <double> RealC(UsedBasisFunctions,1);

        RealC.settoproduct(RealA_inv,RealB);
        if(verbose_level>1)
        {
            cout << "C= " << endl;
            RealC.dumpmatrix();
        }

        double cvalue=0.0f;
        bool success;
        for(int i2=0; i2<UsedBasisFunctions; i2++)
        {
            RealC.getvalue(i2,0,cvalue,success);
            C[i2]=(segPrecisionTYPE)(cvalue);
        }

        //#ifdef _OPENMP
        //#pragma omp parallel
        //#endif
        for (int iz=0; iz<maxiz; iz++)
        {
            segPrecisionTYPE currxpower[maxAllowedBCPowerOrder];
            segPrecisionTYPE currypower[maxAllowedBCPowerOrder];
            segPrecisionTYPE currzpower[maxAllowedBCPowerOrder];
            segPrecisionTYPE tmpbiasfield=0.0f;
            for (int iy=0; iy<maxiy; iy++)
            {
                for (int ix=0; ix<maxix; ix++)
                {
                    linearindexes=(iz)*(this->nx)*(this->ny)+(iy)*(this->nx)+ix;
                    currshortindex=this->L2S[linearindexes];
                    if(currshortindex>=0)
                    {
                        tmpbiasfield=0.0f;
                        xpos=(((segPrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                        ypos=(((segPrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                        zpos=(((segPrecisionTYPE)iz-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z);

                        // Get the polynomial basis order
                        int order=1;
                        currxpower[0]=1;
                        currypower[0]=1;
                        currzpower[0]=1;
                        int orderminusone=0;
                        int maxorderplusone=this->biasFieldOrder+1;
                        while (order<maxorderplusone)
                        {
                            currxpower[order]=currxpower[orderminusone]*xpos;
                            currypower[order]=currypower[orderminusone]*ypos;
                            currzpower[order]=currzpower[orderminusone]*zpos;
                            order++;
                            orderminusone++;
                        }
                        // Estimate the basis
                        int ind2=0;
                        for(int order=0; order<=this->biasFieldOrder; order++)
                        {
                            for(int xorder=0; xorder<=order; xorder++)
                            {
                                for(int yorder=0; yorder<=(order-xorder); yorder++)
                                {
                                    int zorder=order-yorder-xorder;
                                    tmpbiasfield-=C[ind2]*currxpower[xorder]*currypower[yorder]*currzpower[zorder];
                                    ind2++;
                                }
                            }
                        }
                        this->BiasField[currshortindex+multispec*this->numelmasked]=tmpbiasfield;
                    }
                }
            }
        }

        for( int i=0; i<UsedBasisFunctions; i++)
        {
            this->BiasField_coeficients[i+multispec*UsedBasisFunctions]=C[i];
        }
        delete [] Basis;
        delete [] Tempvar;
    }
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

void seg_EM::RunBiasField2D()
{
    if(!this->biasFieldStatus)
    {
        return;
    }

    if(this->verbose_level>0)
    {
        cout << "Optimising the Bias Field with order " << this->biasFieldOrder<< endl;
        if(this->nu>1)
        {
            cout<< "Assuming fully decoupled bias-fields" << endl;
        }
        flush(cout);
    }
    int reduxfactor=reduxFactorForBias;
    long nrOfClasses = this->numb_classes;
    //nrOfClasses = 1;
    //int nrOfClasses = non_PV_numclass;
    long TotalLength = this->numelmasked;
    int UsedBasisFunctions=(int)((this->biasFieldOrder+1) * (this->biasFieldOrder+2)/2);
    segPrecisionTYPE * sampledData = static_cast<segPrecisionTYPE *>(this->InputImage->data);


    // Precompute Powers depending on the current BiasOrder
    int PowerOrder [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2)]= {0};
    int ind=0;
    for(int order=0; order<=this->biasFieldOrder; order++)
    {
        for(int xorder=0; xorder<=order; xorder++)
        {
            int yorder=order-xorder;
            PowerOrder[ind] =xorder;
            PowerOrder[ind+1] =yorder;
            ind += 2;
        }
    }
    segPrecisionTYPE invV[maxNumbClass];
    segPrecisionTYPE currM[maxNumbClass];

    for(long multispec=0; multispec<this->nu; multispec++)
    {
        sampledData = static_cast<segPrecisionTYPE *>(this->InputImage->data);
        sampledData = &sampledData[multispec*this->numel];
        // Precompute the M and V  inverses
        for(long i=0; i<nrOfClasses; i++)
        {
            invV[i]=1.0f/(V[i*this->nu*this->nu+multispec+multispec*this->nu]);
            currM[i]=M[i*this->nu+multispec];
        }
        segPrecisionTYPE A [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2))/2*((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2)]= {0.0f};
        segPrecisionTYPE B [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2))/2]= {0.0f};
        segPrecisionTYPE C [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2))/2]= {0.0f};


        // Precompute sizes
        int col_size = (int)(this->nx);
        int maxix = (int)(this->nx);
        int maxiy = (int)(this->ny);
        int Dims3d[2]= {0};
        Dims3d[0]=maxix;
        Dims3d[1]=maxiy;

        // Precompute number of samples as it was never computed
        int samplecount=0;
        int linearindexes=0;
        int currshortindex=0;
        if(this->numelbias==0)
        {
            for (int iy=0; iy<maxiy; iy+=reduxfactor)
            {
                for (int ix=0; ix<maxix; ix+=reduxfactor)
                {
                    linearindexes=iy*col_size+ix;
                    if(this->L2S[linearindexes]>=0)
                    {
                        samplecount++;
                    }
                }
            }
            if(verbose_level>0)
            {
                cout << "Number of samples for BiasField = " << samplecount<<"\n";
                flush(cout);
            }
            this->numelbias=samplecount;
        }
        else
        {
            samplecount=this->numelbias;
        }



        // CALC MATRIX A

        // Calc W (Van Leemput 1999 eq 7)
        //cout << "Calculating C = inv(A'WA) WR";
        //flush(cout);
        segPrecisionTYPE * Tempvar= new segPrecisionTYPE [samplecount] ();
        segPrecisionTYPE Tempvar_tmp=0;
        currshortindex=0;
        int tempvarindex=0;
        for (int iy=0; iy<maxiy; iy+=reduxfactor)
        {
            for (int ix=0; ix<maxix; ix+=reduxfactor)
            {
                currshortindex=this->L2S[iy*col_size+ix];
                if(currshortindex>=0)
                {
                    Tempvar_tmp=0;
                    for(long j=0; j<nrOfClasses; j++)
                    {
                        Tempvar_tmp+=Expec[currshortindex+TotalLength*j]*invV[j];
                    }
                    Tempvar[tempvarindex]=Tempvar_tmp;
                    tempvarindex++;
                }
            }
        }
        // Precompute shifts
        segPrecisionTYPE not_point_five_times_dims_x=(0.5f*(segPrecisionTYPE)Dims3d[0]);
        segPrecisionTYPE not_point_five_times_dims_y=(0.5f*(segPrecisionTYPE)Dims3d[1]);

        segPrecisionTYPE inv_not_point_five_times_dims_x=1.0f/(0.5f*(segPrecisionTYPE)Dims3d[0]);
        segPrecisionTYPE inv_not_point_five_times_dims_y=1.0f/(0.5f*(segPrecisionTYPE)Dims3d[1]);


        segPrecisionTYPE * Basis= new segPrecisionTYPE[UsedBasisFunctions]();
        segPrecisionTYPE xpos=0.0f;
        segPrecisionTYPE ypos=0.0f;
        int x_bias_index_shift=0;
        int y_bias_index_shift=1;
        segPrecisionTYPE current_Tempvar=0.0f;
        // Calc A'WA (Van Leemput 1999 eq 7)
        tempvarindex=0;
        segPrecisionTYPE * Basisptr1= (segPrecisionTYPE *) Basis;
        segPrecisionTYPE * Basisptr2= (segPrecisionTYPE *) Basis;
        segPrecisionTYPE * Aptr= (segPrecisionTYPE *) A;
        for (int iy=0; iy<maxiy; iy+=reduxfactor)
        {
            for (int ix=0; ix<maxix; ix+=reduxfactor)
            {
                linearindexes=(iy)*(this->nx)+ix;
                currshortindex=this->L2S[linearindexes];
                if(currshortindex>=0)
                {
                    Basisptr1= (segPrecisionTYPE *) Basis;
                    current_Tempvar=Tempvar[tempvarindex];
                    xpos=(((segPrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                    ypos=(((segPrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                    x_bias_index_shift=0;
                    y_bias_index_shift=1;
                    for(int j2=0; j2<UsedBasisFunctions; j2++,x_bias_index_shift+=2,y_bias_index_shift+=2, Basisptr1++)
                    {
                        // Because Powerorder is always int, use a special power function (faster)
                        *Basisptr1=(pow_int(xpos,PowerOrder[x_bias_index_shift])*pow_int(ypos,PowerOrder[y_bias_index_shift]));
                    }
                    Basisptr1= (segPrecisionTYPE *) Basis;
                    Aptr= (segPrecisionTYPE *) A;
                    for(int j2=0; j2<UsedBasisFunctions; j2++, Basisptr1++)
                    {
                        Basisptr2= &Basis[j2];
                        Aptr= &A[j2+j2*UsedBasisFunctions];
                        for(int i2=j2; i2<UsedBasisFunctions; i2++, Aptr++, Basisptr2++)
                        {
                            (*Aptr)+=(*Basisptr2)*(current_Tempvar)*(*Basisptr1);
                        }
                    }
                    tempvarindex++;
                }
            }
        }

        seg_Matrix <double> RealA(UsedBasisFunctions,UsedBasisFunctions);

        for(int j2=0; j2<UsedBasisFunctions; j2++)
        {
            for(int i2=j2; i2<UsedBasisFunctions; i2++)
            {
                RealA.setvalue(i2,j2,(double)(A[i2+j2*UsedBasisFunctions]));
                RealA.setvalue(j2,i2,(double)(A[i2+j2*UsedBasisFunctions]));
            }
        }

        seg_Matrix <double> RealA_inv(UsedBasisFunctions);
        RealA_inv.copymatrix(RealA);
        RealA_inv.invert();

        if(verbose_level>1)
        {
            seg_Matrix <double> RealA_test(UsedBasisFunctions);
            RealA_test.settoproduct(RealA,RealA_inv);
            RealA_test.comparetoidentity();
            //RealA.dumpmatrix();
        }



        // CALC MATRIX B

        //Precompute WR (Van Leemput 1999 eq 7)
        segPrecisionTYPE Wi;
        segPrecisionTYPE Wij;
        segPrecisionTYPE Yest;
        segPrecisionTYPE Ysum;
        tempvarindex=0;
        for (int iy=0; iy<maxiy; iy+=reduxfactor)
        {
            for (int ix=0; ix<maxix; ix+=reduxfactor)
            {
                linearindexes=(iy)*(this->nx)+ix;
                currshortindex=this->L2S[linearindexes];
                if(currshortindex>=0)
                {
                    Wi=0;
                    Wij=0;
                    Yest=0;
                    Ysum=0;
                    for(long j=0; j<nrOfClasses; j++)
                    {
                        segPrecisionTYPE tmpexpec = (segPrecisionTYPE)this->Expec[currshortindex+TotalLength*j];
                        Wij=tmpexpec*(invV[j]);
                        Wi+=Wij;
                        Yest+=Wij*(currM[j]);
                        Ysum+=Wij;
                    }
                    Tempvar[tempvarindex]=Wi*(sampledData[linearindexes]-(Yest/Ysum));
                    tempvarindex++;
                }
            }
        }

        for(int i2=0; i2<UsedBasisFunctions; i2++)
        {
            tempvarindex=0;
            B[i2]=0;
            for (int iy=0; iy<maxiy; iy+=reduxfactor)
            {

                for (int ix=0; ix<maxix; ix+=reduxfactor)
                {
                    linearindexes=(iy)*(this->nx)+ix;
                    currshortindex=this->L2S[linearindexes];
                    if(currshortindex>=0)
                    {
                        B[i2]+=pow_int((((segPrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x),PowerOrder[0+i2*2])*
                                pow_int((((segPrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y),PowerOrder[1+i2*2])*
                                Tempvar[tempvarindex];

                        if(B[i2]!=B[i2])
                        {
                            B[i2]=1;
                        }
                        tempvarindex++;
                    }
                }
            }
        }

        seg_Matrix <double> RealB(UsedBasisFunctions,1);

        for(int i2=0; i2<UsedBasisFunctions; i2++)
        {
            RealB.setvalue(i2,0,(double)(B[i2]));
        }

        seg_Matrix <double> RealC(UsedBasisFunctions,1);

        RealC.settoproduct(RealA_inv,RealB);
        if(verbose_level>1)
        {
            cout << "C= " << endl;
            RealC.dumpmatrix();
        }

        double cvalue=0.0f;
        bool success;
        for(int i2=0; i2<UsedBasisFunctions; i2++)
        {
            RealC.getvalue(i2,0,cvalue,success);
            C[i2]=(segPrecisionTYPE)(cvalue);
        }

        segPrecisionTYPE currxpower[maxAllowedBCPowerOrder];
        segPrecisionTYPE currypower[maxAllowedBCPowerOrder];
        segPrecisionTYPE tmpbiasfield=0.0f;
        for (int iy=0; iy<maxiy; iy++)
        {
            for (int ix=0; ix<maxix; ix++)
            {
                linearindexes=(iy)*(this->nx)+ix;
                currshortindex=this->L2S[linearindexes];
                if(currshortindex>=0)
                {
                    tmpbiasfield=0.0f;
                    xpos=(((segPrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                    ypos=(((segPrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);

                    // Get the polynomial order power
                    int order=1;
                    currxpower[0]=1;
                    currypower[0]=1;
                    int orderminusone=0;
                    int maxorderplusone=this->biasFieldOrder+1;
                    while (order<maxorderplusone)
                    {
                        currxpower[order]=currxpower[orderminusone]*xpos;
                        currypower[order]=currypower[orderminusone]*ypos;
                        order++;
                        orderminusone++;
                    }

                    ind=0;
                    for(int order=0; order<=this->biasFieldOrder; order++)
                    {
                        for(int xorder=0; xorder<=order; xorder++)
                        {
                            int yorder=order-xorder;
                            tmpbiasfield-=C[ind]*currxpower[xorder]*currypower[yorder];
                            ind++;
                        }
                    }
                    this->BiasField[currshortindex+multispec*this->numelmasked]=tmpbiasfield;
                }
            }

        }

        for( int i=0; i<UsedBasisFunctions; i++)
        {
            this->BiasField_coeficients[i+multispec*UsedBasisFunctions]=C[i];
        }
        delete [] Basis;
        delete [] Tempvar;
    }
}



/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

void seg_EM::RunPriorRelaxation()
{
    if(this->relaxStatus)
    {
        if((int)(this->verbose_level)>(int)(0))
        {
            cout << "Relaxing Priors"<< endl;
        }
        for(long i=0; i<(this->numb_classes*this->numelmasked); i++)this->ShortPrior[i]=Expec[i];
        GaussianFilter4D_cArray(this->ShortPrior,this->S2L,this->L2S,this->relaxGaussKernelSize,this->CurrSizes,CSFclass);
        float * PriorsPtr = static_cast<float *>(this->Priors->data);
        long currindex=0;
        for(long k=0; k<this->numb_classes; k++)
        {
            currindex=k*this->numelmasked;
            for(long i=0; i<(this->numelmasked); i++)
            {
                this->ShortPrior[currindex]*=(1-this->relaxFactor);
                this->ShortPrior[currindex]+=(this->relaxFactor)*PriorsPtr[S2L[i]+k*this->numel];
                currindex++;
            }
        }
    }
    return;
}
#endif
