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

    this->dimensions=1;
    this->nx=1;
    this->ny=0;
    this->nz=0;
    this->nu=_nu*_nt;
    this->nt=1;

    this->dx=0;
    this->dy=0;
    this->dz=0;


    this->numElements=0;


    this->iter=0;
    this->ratio=1000;

    this->numberOfClasses=_numb_classes;
    this->M= new segPrecisionTYPE [maxMultispectalSize*maxNumbClass];
    for(int i=0; i<(maxMultispectalSize*maxNumbClass); i++)
    {
        this->M[i]=0.0f;
    }
    this->V= new segPrecisionTYPE [maxMultispectalSize*maxMultispectalSize*maxNumbClass];
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
        this->convCrit=0.0005;

    this->verbose_level=0;
    this->loglik=2.0;
    this->oldloglik=1.0;


    this->maskImageStatus=0;
    this->Mask=NULL;        // pointer to external
    this->numElementsMasked=0;

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
    this->outliernessThreshold=0.0f;
    this->outliernessRatio=0.01f;

    this->biasFieldStatus=false;
    this->biasFieldOrder=0;
    this->BiasField=NULL;
    this->biasFieldCoeficients=NULL;
    this->biasFieldRatio=0;
    this->numelBias=0;

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

    if(this->biasFieldCoeficients!=NULL)
    {
        delete [] this->biasFieldCoeficients;
    }
    this->biasFieldCoeficients=NULL;

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

    if(this->M!=NULL)
	delete [] this->M;

    if(this->V!=NULL)
	delete [] this->V;
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
    // Count the number of dimensions with size above 1.
    this->dimensions=(int)((_r->nx)>1)+(int)((_r->ny)>1)+(int)((_r->nz)>1)+(int)((_r->nt)>1)+(int)((_r->nu)>1);
    this->nx=_r->nx;
    this->ny=_r->ny;
    this->nz=_r->nz;
    this->nt=(_r->nt>1)?_r->nt:1;
    this->nu=(_r->nu>1)?_r->nu:1;
    this->dx=_r->dx;
    this->dy=_r->dy;
    this->dz=_r->dz;
    this->numElements=_r->nz*_r->ny*_r->nx;
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
void seg_EM::SetMAP( segPrecisionTYPE *_M, segPrecisionTYPE* _V)
{
    this->mapStatus=true;
    this->MAP_M=new segPrecisionTYPE[this->numberOfClasses];
    this->MAP_V=new segPrecisionTYPE[this->numberOfClasses];


    for(int i=0; i<this->numberOfClasses; i++)
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
    this->dimensions=(int)((_r->nx)>1)+(int)((_r->ny)>1)+(int)((_r->nz)>1);
    if(this->nx==_r->nx && this->ny==_r->ny && this->nz==_r->nz)
    {
        return;
    }
    else
    {
        cout << "ERROR: Priors have wrong dimensions" << endl;
        return;
    }


}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Sets the filename of the output
/// @param _f A char pointer to the output file name.
///

void seg_EM::SetFilenameOut(char * _f)
{
    this->filenameOut = _f;
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
        cout << "ERROR: Mask has wrong dimensions" << endl;
        this->Mask->nt=1;
        this->Mask->nu=1;
        this->Mask->ndim=3;
    }

    bool * MaskDataPtr = static_cast<bool *>(Mask->data);
    this->numElementsMasked=0;
    for(int i=0; i<this->numElements; i++,MaskDataPtr++)
    {
        if((*MaskDataPtr)>0)
        {
            this->numElementsMasked++;
        }
    }
    if(this->nx==_r->nx && this->ny==_r->ny && this->nz==_r->nz)
    {
        return;
    }
    else
    {
        cout << "ERROR: Mask has wrong dimensions" << endl;
        return;
    }
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Sets the filename of the output
/// @param verblevel An unsigned int defining the verbose level.
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
/// @brief Sets the regularisation of the -ffdiagonal componentes of the covariance matrix
/// @param reg A segPrecisionTYPE defining the regularisation value.
///
/// The off-diagonal componentes of the covariance matrix are divided by this->reg_factor in order to improve the rank of the matrix (needs to be positive-semi-definite covariance).
/// Thus, reg>1, making the off-diagonal smaller than it should be. This will reduce the model fit quality, but it will help the stability of the algorithm. It also makes the covariances more isotropic, reducing the amount of covariance between chanels.
///
void seg_EM::SetRegValue(segPrecisionTYPE reg)
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
/// @brief Sets the maximum number of iterations.
/// @param verblevel An unsigned int defining the maximum number of iterations.
///
/// this->maxIteration overrides the convergence criteria.
///
void seg_EM::SetMaximalIterationNumber(unsigned int numberiter)
{
    if(numberiter<1)
    {
        this->maxIteration=1;
    }
    else
    {
        this->maxIteration=numberiter;
    }
    return;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Sets the minimum number of iterations
/// @param verblevel An unsigned int defining the minimum number of iterations.
///
/// this->minIteration overrides the convergence criteria.
///
void seg_EM::SetMinIterationNumber(unsigned int numberiter)
{
    this->minIteration=numberiter;
    return;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Sets the convergence criteria threshold
/// @param inConvCrit An float defining the minimum log lik ration increase
///
///
void seg_EM::SetConvergenceCriteria(float inConvCrit)
{
    this->convCrit=inConvCrit;
    return;
}




/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Allocates and sets the beta parameter of the MRF energy function.
/// @param strength The MRF beta parameter
///
/// The beta parameter is defined as \f$\exp(-\beta U_{MRF})\f$
///
void seg_EM::SetMRF(segPrecisionTYPE strength)
{
    this->mrfStatus=true;
    this->mrfStrength=strength;
    this->MRFTransitionMatrix=new segPrecisionTYPE[this->numberOfClasses*this->numberOfClasses]();
    CreateDiagonalMRFTransitionMatrix();


    if(this->maskImageStatus>0)
    {
        this->MRF=new segPrecisionTYPE[this->numElementsMasked*this->numberOfClasses*this->nu*this->nt]();
        for(int i=0; i<(this->numElementsMasked*this->numberOfClasses); i++)
        {
            MRF[i]=(segPrecisionTYPE)(1.0);
        }
    }
    else
    {
        this->MRF=new segPrecisionTYPE[this->numElements*this->numberOfClasses*this->nu*this->nt]();
        for(int i=0; i<(this->numElements*this->numberOfClasses); i++)
        {
            MRF[i]=(segPrecisionTYPE)(1.0);
        }
    }

    CreateDiagonalMRFTransitionMatrix();
    return;
}


/// @brief Sets the flag for using the Derrichlet prior, the associated relaxation factor and the Gaussian standard deviation of the regularisation, as deffined in the AdaPT paper (Cardoso et al. Neuroimage 2012).
/// @param relax_factor The Derrichlet prior relaxation factor
/// @param relax_gauss_kernel The Gaussian standard deviation of the regularisation
///
void seg_EM::SetRelaxation(segPrecisionTYPE relax_factor,segPrecisionTYPE relax_gauss_kernel)
{
    this->relaxStatus=true;
    this->relaxFactor=relax_factor;
    this->relaxGaussKernelSize=relax_gauss_kernel;
    return;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Sets the flag for using the Bias field correction, the associated relaxation factor and the Gaussian standard deviation of the regularisation, as defined in the AdaPT paper (Cardoso et al. Neuroimage 2012).
/// @param _BiasFieldOrder The order of the polynomial bias field
/// @param _BiasFieldRatio The convergence ratio below which the bias field correction is used.
///
/// If the _BiasFieldRatio value is large (>0.1), the bias field can "go" in the wrong direction.
///
void seg_EM::SetBiasField(int _BiasFieldOrder, segPrecisionTYPE _BiasFieldRatio)
{
    if(_BiasFieldOrder<=0){
        this->biasFieldStatus=true;
        this->biasFieldOrder=0;
        this->biasFieldRatio=0.0f;
    }
    else{
        this->biasFieldStatus=true;
        this->biasFieldOrder=_BiasFieldOrder;
        this->biasFieldCoeficients = new segPrecisionTYPE[((int)(((float)(this->biasFieldOrder)+1.0f) * ((float)(this->biasFieldOrder)+2.0f)/2.0f *((float)(this->biasFieldOrder)+3.0f)/3.0f))*this->nu*this->nt]();
        this->biasFieldRatio=_BiasFieldRatio;
        if(this->maskImageStatus>0)
        {
            this->BiasField = new segPrecisionTYPE[this->numElementsMasked*this->nu*this->nt]();
        }
        else
        {
            this->BiasField = new segPrecisionTYPE[this->numElements*this->nu*this->nt]();
        }
    }
    return;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Sets the flag for using the Outlier detection, the associated relaxation factor and the Gaussian standard deviation of the regularisation, as defined in the AdaPT paper (Cardoso et al. Neuroimage 2012).
/// @param _OutliernessThreshold The Mahalanobis distance threshold for classifying an observation as an outlier.
/// @param _OutliernessRatio The convergence ration below which the outlier detection is used.
///
/// If the _OutliernessRatio value is large (>0.1), the outlierness can "go" in the wrong direction.
///
void seg_EM::SetOutlierness(segPrecisionTYPE _OutliernessThreshold, segPrecisionTYPE _OutliernessRatio)
{

    this->outliernessStatus=true;
    this->outliernessThreshold=_OutliernessThreshold;
    this->outliernessRatio=_OutliernessRatio;
    return;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Initialises an allocates all the vectors.
///
/// If priors are not used, it initialises the means and variances using the intensity distributions. It also allocates the this->ShortPrior and this->Outlierness vectors when not previously deffined, even though they are all equal to 1/K and 1 respectively (thus not affecting the convergence);
///
void seg_EM::InitializeAndAllocate()
{


    if(this->priorsStatus)
    {
        register long numel=(int)(this->Mask->nx*this->Mask->ny*this->Mask->nz);
        register long numel_masked=0;

        bool * Maskptrtmp = static_cast<bool *> (this->Mask->data);
        for (long i=0; i<numel; i++, Maskptrtmp++)
        {
            *Maskptrtmp?numel_masked++:0;
        }
        int pluspv=(int)(this->pvModelStatus)*2;

        this->Expec = new segPrecisionTYPE [numel_masked*(this->numberOfClasses+pluspv)] ();
        segPrecisionTYPE * tempExpec= (segPrecisionTYPE *) Expec;
        this->ShortPrior = new segPrecisionTYPE [numel_masked*(this->numberOfClasses+pluspv)] ();
        segPrecisionTYPE * tempShortPrior= (segPrecisionTYPE *) ShortPrior;
        segPrecisionTYPE * PriorPTR = static_cast<segPrecisionTYPE *>(this->Priors->data);
        for(long cl=0; cl<this->numberOfClasses; cl++)
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
    }
    else
    {
        this->InitializeMeansUsingIntensity();
        int tmpnumb_elem=0;
        if(this->maskImageStatus>0)
        {
            tmpnumb_elem=(this->numElementsMasked*(this->numberOfClasses+(int)(this->pvModelStatus)*2));
        }
        else
        {
            tmpnumb_elem=(numElements*(this->numberOfClasses+(int)(this->pvModelStatus)*2));
        }

        segPrecisionTYPE tmpnumb_clas=((this->numberOfClasses+(int)(this->pvModelStatus)*2));
        this->Expec=new segPrecisionTYPE [tmpnumb_elem] ();
        this->ShortPrior=new segPrecisionTYPE [tmpnumb_elem] ();
        for(int i=0; i<tmpnumb_elem; i++)
        {
            this->Expec[i]=1.0/tmpnumb_clas;
            this->ShortPrior[i]=1.0/tmpnumb_clas;
        }
        this->RunExpectation();

    }

    for (int cl=0; cl<this->numberOfClasses; cl++)
    {
        if(this->outliernessStatus)
        {
            int tmpnumb_elem=0;
            if(this->maskImageStatus>0)
            {
                tmpnumb_elem=(this->numElementsMasked*(this->numberOfClasses+(int)(this->pvModelStatus)*2));
            }
            else
            {
                tmpnumb_elem=(numElements*(this->numberOfClasses+(int)(this->pvModelStatus)*2));
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
/// @brief Initialises and normalises the image and the Population Priors. The MAP priors are also normalised using the same transformation.
///
/// If the mask is not used, it also creates a binary mask image which is 1 everywhere. This mask has to be created before using this->CreateShort2LongMatrix() and this->CreateLong2ShortMatrix(), as they will use the this->Mask to estimate the L2S and S2L mappings.

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
        this->numElementsMasked=this->numElements;
    }

    this->InitializeAndNormalizeNaNPriors();
    this->InitializeAndNormalizeImage();

    this->CreateShort2LongMatrix();
    this->CreateLong2ShortMatrix();


    if(this->mapStatus)
    {
        for(int i=0; i<this->numberOfClasses; i++)
        {
            this->MAP_M[i]=logf(((this->MAP_M[i]-this->rescale_min[0])/(this->rescale_max[0]-this->rescale_min[0]))+1)/0.693147181;;
            if(this->verbose_level>0)
            {
                cout << "MAP_M["<<i<<"] = "<<this->MAP_M[i]<< endl;
            }
            this->M[i]=this->MAP_M[i];
            this->V[i]=1.0/this->numberOfClasses;
        }
    }

    return;

}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Initialises, normalises and log transforms (for the bias field computation) the image itself.
///
/// The log transformation is used to make the bias field aditive.

void seg_EM::InitializeAndNormalizeImage()
{

    if(this->verbose_level>0)
    {
        cout<< "Normalizing Input Image" << endl;
    }
    int numel=(int)(this->InputImage->nx*this->InputImage->ny*this->InputImage->nz);


    seg_changeDatatype<segPrecisionTYPE>(this->InputImage);


    for(long udir=0; udir<this->nu; udir++) // Per Multispectral Image
    {
        bool * brainmaskptr = static_cast<bool *> (this->Mask->data);
        segPrecisionTYPE * Inputptrtmp = static_cast<segPrecisionTYPE *>(this->InputImage->data);
        segPrecisionTYPE * Inputptr=&Inputptrtmp[numel*udir];

        segPrecisionTYPE tempmax=-(1.0e32);
        segPrecisionTYPE tempmin=1.0e32;

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
/// @brief Initialises, normalises and checks for NaN's in the priors.
///
/// It also does an empirical check of the priors.
/// If the prior sum (over the classes) per voxel is <0 or >1000, then something is probably wrong or the user is using it wrong.
/// If this happens, output an Error.
///
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
            if(this->InputImage->datatype!=NIFTI_TYPE_FLOAT32 || this->InputImage->datatype!=NIFTI_TYPE_FLOAT64)
            {
                segPrecisionTYPE * priorsptr = static_cast<segPrecisionTYPE *>(this->Priors->data);
                bool * brainmaskptr = static_cast<bool *> (this->Mask->data);

                for (int i=0; i<numel; i++)
                {
                    if(brainmaskptr[i])
                    {
                        segPrecisionTYPE tempsum=0;
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
/// @brief Initialises the mean and variance parameters using the image intensity.
///
/// Sets the initial value of the class means by first estimating the image histogram, and the partitionin the histogram using percentiles.
/// It uses only the first image chanel if the input image is multimodal.
/// It also sets the variance to the overall variance devided by K.
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
        segPrecisionTYPE meanval=0.0;
        segPrecisionTYPE variance=0.0;

        for(int i=0; i<this->numElements; i++)
        {
            if(this->maskImageStatus==0 || MaskDataPtr[i]>0)
            {
                mycounter++;
                meanval+=(Intensity_PTR[i+ms*this->numElements]);
            }
        }
        meanval=meanval/mycounter;

        for(int i=0; i<this->numElements; i++)
        {
            if(this->maskImageStatus==0 || MaskDataPtr[i]>0 )
            {
                variance+=pow((meanval-Intensity_PTR[i+ms*this->numElements]),2);
            }
        }
        variance=variance/mycounter;
        int histogram[1001];
        for(int i=0; i<1000; i++)
        {
            histogram[i]=0;
        }
        segPrecisionTYPE tmpmax=-1e32f;
        segPrecisionTYPE tmpmin=1e32f;



        for(int i=0; i<this->numElements; i++)
        {
            if(this->maskImageStatus==0 || MaskDataPtr[i]>0)
            {
                if(tmpmax<(int)(Intensity_PTR[i+ms*this->numElements]))
                {
                    tmpmax=(int)(Intensity_PTR[i+ms*this->numElements]);
                }
                if(tmpmin>(int)(Intensity_PTR[i+ms*this->numElements]))
                {
                    tmpmin=(int)(Intensity_PTR[i+ms*this->numElements]);
                }
            }
        }

        for(int i=0; i<this->numElements; i++)
        {
            if(this->maskImageStatus==0 || MaskDataPtr[i]>0)
            {

                int index4hist=(int)(1000.0f*(segPrecisionTYPE)(Intensity_PTR[i+ms*this->numElements]-tmpmin)/(segPrecisionTYPE)(tmpmax-tmpmin));
                if((index4hist>1000) & (index4hist<0))
                {
                    cout<< "error"<<endl;
                }
                histogram[(int)(1000.0*(segPrecisionTYPE)(Intensity_PTR[i+ms*this->numElements]-tmpmin)/(segPrecisionTYPE)(tmpmax-tmpmin))]++;
            }
        }


        for(int clas=0; clas<this->numberOfClasses; clas++)
        {
            segPrecisionTYPE tmpsum=0;
            int tmpindex=0;
            segPrecisionTYPE percentile=((segPrecisionTYPE)clas+1)/(this->numberOfClasses+1);
            for(int i=999; i>0; i--)
            {
                tmpsum+=histogram[i];
                tmpindex=i;
                if((segPrecisionTYPE)(tmpsum)>((1.0f-percentile)*(segPrecisionTYPE)(mycounter)))
                {
                    i=0;
                }
            }
            M[clas*this->nu+ms]=segPrecisionTYPE(tmpindex)*(tmpmax-tmpmin)/1000.0f+(tmpmin);
            V[clas*this->nu*this->nu+ms*this->nu+ms]=variance/this->numberOfClasses/2;
        }
    }


    for (int cl=0; cl<this->numberOfClasses; cl++)
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
/// @brief Initialises the MRF transtition matri.
///
/// Sets MRF transition matrix to 0 in the diagonal and to this->mrfStrength in the offdiagonal.
///
void seg_EM::CreateDiagonalMRFTransitionMatrix()
{
    for(int i=0; i<this->numberOfClasses; i++)
    {
        for(int j=0; j<this->numberOfClasses; j++)
        {
            if(j==i)
            {
                this->MRFTransitionMatrix[i+j*this->numberOfClasses]=0;
            }
            else
            {
                this->MRFTransitionMatrix[i+j*this->numberOfClasses]=this->mrfStrength;
            }
        }
    }
    return;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Creates the CurrSizes variable used when relaxin the priors or when outputing the bias field.
///
/// The CurrSizes is used as a bridge to stransport the image dimensions between the seg_EM object and basic c functions
///
void seg_EM::CreateCurrSizes()
{
    this->CurrSizes = new ImageSize [1]();
    CurrSizes->numel=(int)(this->nx*this->ny*this->nz);
    CurrSizes->xsize=this->nx;
    CurrSizes->ysize=this->ny;
    CurrSizes->zsize=this->nz;
    CurrSizes->usize=(this->nu>1)?this->nu:1;
    CurrSizes->tsize=(this->nt>1)?this->nt:1;
    CurrSizes->numclass=this->numberOfClasses;
    CurrSizes->numelmasked=this->numElementsMasked;
    CurrSizes->numelbias=0;

    return;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Creates a mapping to go from the masked vector to the full image.
///
/// This saves a lot of memory as the brain occupies 1/3 of the image size. This mapping is the inverse of the this->CreateLong2ShortMatrix() mapping.
/// \sa this->CreateLong2ShortMatrix()
///
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
        this->numElementsMasked=numel_masked;

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
/// @brief Creates a mapping to go from the full image to the masked vector.
///
/// This saves a lot of memory as the brain occupies 1/3 of the image size. This mapping is the inverse of the this->CreateShort2LongMatrix() mapping.
/// \sa this->CreateShort2LongMatrix()
///
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
/// @brief Gets the Mean vector. This is normally used by the caller after the algorithm converges to obtain the model parameters.
///
/// This is normally used by the caller after the algorithm converges to obtain the model parameters. Note that the means have to be pub back into a non-log non-normalised form.
///
segPrecisionTYPE * seg_EM::GetMeans()
{
    segPrecisionTYPE * OutM= new segPrecisionTYPE [this->nu*this->numberOfClasses];
    int index=0;
    for(int j=0; j<(this->numberOfClasses); j++)
    {
        for(int i=0; i<(this->nu); i++)
        {
            //cout << "M["<<j<<"]="<<this->M[j+i*this->numb_classes]<<endl;
            segPrecisionTYPE resize=exp((this->M[index++])*0.693147181)-1;
            OutM[j+i*this->numberOfClasses]=(resize*(this->rescale_max[i]-this->rescale_min[i])+this->rescale_min[i]);
        }
    }

    return OutM;
}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Gets the Covariance matrix.
///
/// This is normally used by the caller after the algorithm converges to obtain the model parameters. Note that the variances are in the log-transformed space and are not rescalled.
///
segPrecisionTYPE * seg_EM::GetSTD()
{
    segPrecisionTYPE * OutV= new segPrecisionTYPE [this->nu*this->numberOfClasses*this->numberOfClasses];
    for(int i=0; i<(this->nu); i++)
    {
        for(int j=0; j<(this->numberOfClasses*this->numberOfClasses); j++)
        {
            OutV[j+i*this->numberOfClasses*this->numberOfClasses] = this->V[j+i*this->numberOfClasses*this->numberOfClasses];
        }
    }

    return OutV;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Gets the final segmentation results.
///
/// Allocates the output result nifti_image structure and fills it with the already normalised results from the expectation vector.
/// As Expec_PTR is only deffined within the mask, the this->S2L[i] is used to go from a masked index to the image index.
///
nifti_image * seg_EM::GetResult()
{

    nifti_image * Result = nifti_copy_nim_info(this->InputImage);
    Result->dim[0]=4;
    Result->dim[4]=this->numberOfClasses;
    Result->dim[5]=1;
    Result->scl_inter=0;
    Result->scl_slope=1;
    Result->datatype=this->InputImage->datatype;
    Result->cal_max=1;
    Result->cal_min=0;
    nifti_set_filenames(Result,(char*)this->filenameOut.c_str(),0,0);
    nifti_update_dims_from_array(Result);
    nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
    Result->data = (void *) calloc(Result->nvox, sizeof(segPrecisionTYPE));
    segPrecisionTYPE * Resultdata = static_cast<segPrecisionTYPE *>(Result->data);
    // First, zero all the elements of Resultdata
    for(unsigned int i=0; i<Result->nvox; i++)
    {
        Resultdata[i]=0;
    }

    int class_nvox=Result->nx*Result->ny*Result->nz;
    // Then, for each class
    for(long currclass=0; currclass<this->numberOfClasses; currclass++)
    {

        // Get the vector pointer for the class
        segPrecisionTYPE * Resultdata_class = &Resultdata[(currclass)*class_nvox];
        segPrecisionTYPE * Expec_PTR = &Expec[(currclass)*this->numElementsMasked];
        // Copy the probability from the masked vector Expec_PTR to the full image vector Resultdata_class.
        // The copying requires using the mapping S2L to go from the masked vector index to the full image index
        for(long i=0; i<(long)this->numElementsMasked; i++,Expec_PTR++)
        {
            Resultdata_class[this->S2L[i]]=*Expec_PTR;
        }
    }

    return Result;

}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Gets the final biasfield corrected result.
///
/// Allocates the output bias field corrected image.
/// To do so, it needs to estimate the bias field, subtract it from the observations and then un-log transform and un-nomalise the image.
/// The bias field correction is only done within the image.
/// The boundaries of the bias field are then smoothed out using a gaussian filter in order to avoid bias field extrapolation problems.
///
nifti_image * seg_EM::GetBiasCorrected(char * filename)
{

    if(this->biasFieldOrder==0)
    {
        return NULL;
    }
    int UsedBasisFunctions=(int)(((float)(this->biasFieldOrder)+1.0f) * ((float)(this->biasFieldOrder)+2.0f)/2.0f *((float)(this->biasFieldOrder)+3.0f)/3.0f);
    segPrecisionTYPE * InputImageData = static_cast<segPrecisionTYPE *>(this->InputImage->data);

    nifti_image * Result = nifti_copy_nim_info(this->InputImage);
    Result->dim[0]=4;
    Result->dim[4]=this->nu;
    Result->datatype=this->InputImage->datatype;
    Result->cal_max=(this->rescale_max[0]);
    Result->scl_inter=0;
    Result->scl_slope=1;

    segPrecisionTYPE * brainmask= new segPrecisionTYPE [this->numElements];
    bool * Maskptrtmp = static_cast<bool *> (this->Mask->data);;

    for(long i=0; i<(long)this->numElements; i++)
    {
        if(InputImageData[i]!=InputImageData[i] )
        {
            brainmask[i]=0.0f;
        }
        else{
            brainmask[i]=Maskptrtmp[i];
        }
    }
    Dillate(brainmask,7,CurrSizes);
    Erosion(brainmask,3,CurrSizes);
    GaussianFilter4D_cArray(brainmask, 3.0f, CurrSizes);


    nifti_set_filenames(Result,filename,0,0);
    nifti_update_dims_from_array(Result);
    nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
    Result->data = (void *) calloc(Result->nvox, sizeof(segPrecisionTYPE));
    segPrecisionTYPE * BiasCorrected_PTR = static_cast<segPrecisionTYPE *>(Result->data);

    segPrecisionTYPE BiasField=0;
    segPrecisionTYPE currxpower[maxAllowedBCPowerOrder]={0};
    segPrecisionTYPE currypower[maxAllowedBCPowerOrder]={0};
    segPrecisionTYPE currzpower[maxAllowedBCPowerOrder]={0};
    segPrecisionTYPE xpos=0.0f;
    segPrecisionTYPE ypos=0.0f;
    segPrecisionTYPE zpos=0.0f;
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
        BiasCorrected_PTR = &BiasCorrected_PTR[multispec*this->numElements];
        InputImageData = static_cast<segPrecisionTYPE *>(this->InputImage->data);
        InputImageData = &InputImageData[multispec*this->numElements];

        segPrecisionTYPE * BiasFieldCoefs_multispec = &this->biasFieldCoeficients[multispec*UsedBasisFunctions];


        for(long i=0; i<(long)this->numElements; i++)
        {
            BiasCorrected_PTR[i]=0;
        }


        segPrecisionTYPE to_resize=0;
        int index_full=0;
        for (int iz=0; iz<this->nz; iz++)
        {
            for (int iy=0; iy<this->ny; iy++)
            {
                for (int ix=0; ix<this->nx; ix++)
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
    delete [] brainmask;

    return Result;
}



/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Gets the final outlier segmentation.
///
/// Estimates 1 - (the sum of the Outlierness vector over all classes).
///
nifti_image * seg_EM::GetOutlierness(char * filename)
{
    nifti_image * Result = nifti_copy_nim_info(this->InputImage);
    Result->dim[0]=3;
    Result->dim[4]=1;
    Result->dim[5]=1;
    Result->datatype=this->InputImage->datatype;
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


    for(int i=0; i<this->numElementsMasked; i++)
    {
        segPrecisionTYPE currsum=0;
        for(int currclass=0; currclass<this->numberOfClasses; currclass++)
        {
            currsum+=this->Outlierness[i+(currclass)*this->numElementsMasked]*Expec[i+(currclass)*this->numElementsMasked];
        }
        Resultdata[S2L[i]]=1-currsum;
    }

    return Result;

}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief This is the main function of the seg_EM object. It starts the segmentation process.
///
/// This should only be called after all the parameters have been set.
///
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
        cout << "Number of voxels inside the mask = " << this->numElementsMasked << endl;
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
        this->RunBiasField();
        //Update Weight
        this->RunPriorRelaxation();

        // Print LogLik depending on the verbose level
        if(this->verbose_level>0 && this->iter>0)
        {
            if(iter>0)
            {
                if ((this->loglik-this->oldloglik)/fabs(this->oldloglik)>0 && (this->loglik-this->oldloglik)/fabs(this->oldloglik)<100)
                {
                    cout<< "Loglik = " << setprecision(7)<<this->loglik <<
                           " : Ratio = " << (this->loglik-this->oldloglik)/fabs(this->oldloglik) << endl;
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
        if((((this->loglik-this->oldloglik)/(fabs(this->oldloglik+this->loglik)/2.0f))<(segPrecisionTYPE)(this->convCrit)
            && this->iter>this->minIteration)
                || iter>=this->maxIteration
	        || (std::isinf(this->loglik) && this->iter>3))
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
        int minutes = (int)floorf(segPrecisionTYPE(end-start)/60.0f);
        int seconds = (int)(end-start - 60*minutes);
        cout << "Finished in "<<minutes<<"min "<<seconds<<"sec"<< endl;
    }
    return;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Optimisation function taking care of the MAximisation step
///
/// Estimates the mean and covariance matrix parameters given the current responsabilities.
///
void seg_EM::RunMaximization()
{

    int verbose=this->verbose_level;
    if(this->verbose_level>0)
    {
        cout<< "Optimising Mixture Model Parameters" << endl;
        flush(cout);
    }
    bool OutliernessFlag=(Outlierness==NULL)?0:1;

    int numel_masked=this->numElementsMasked;
    int num_class=this->numberOfClasses;
    int Expec_offset[maxNumbClass];
    for (int cl=0; cl<num_class; cl++)
    {
        Expec_offset[cl]=cl*numel_masked;
    }

    segPrecisionTYPE * ExpectationTmpPTR = (segPrecisionTYPE *) this->Expec;
    segPrecisionTYPE * OutliernessTmpPTR = (segPrecisionTYPE *) this->Outlierness;
    segPrecisionTYPE * InputImageTmpPtr = static_cast<segPrecisionTYPE *>(this->InputImage->data);
    segPrecisionTYPE * BiasFieldTmpPTR= (segPrecisionTYPE *) this->BiasField;

#ifdef _OPENMP
#pragma omp parallel for shared(InputImageTmpPtr,BiasFieldTmpPTR,OutliernessTmpPTR)
#endif
    // ***********
    // For each class, get all the temporary pointers
    for( int cl=0; cl<num_class; cl++)
    {
        int * S2L_PTR = (int *) this->S2L;
        segPrecisionTYPE * ExpectationPTR = (segPrecisionTYPE *) ExpectationTmpPTR;
        segPrecisionTYPE * OutliernessPTR = (segPrecisionTYPE *) OutliernessTmpPTR;
        segPrecisionTYPE * InputImagePtr = (segPrecisionTYPE *) InputImageTmpPtr;
        segPrecisionTYPE * BiasFieldPTR= (segPrecisionTYPE *) BiasFieldTmpPTR;
        segPrecisionTYPE * T1_PTR2= (segPrecisionTYPE *) InputImageTmpPtr;
        segPrecisionTYPE * BiasField_PTR2= (segPrecisionTYPE *) BiasFieldTmpPTR;

        // MEAN
        // For each multispectral data (or each time point), estimate the mean vector. This involves doing a weighted sum (tempsum/SumPriors) of the observed intensities, weighted by the responsabilities Expec_PTR, the OutliernessPTR. The Estimated mean uses the bias field corrected intensities (T1_PTR-BiasField_PTR)
        for(long Multispec=0; Multispec<this->nu; Multispec++)
        {
            ExpectationPTR=(segPrecisionTYPE *) &Expec[Expec_offset[cl]];
            OutliernessPTR=(segPrecisionTYPE *) &Outlierness[Expec_offset[cl]];
            S2L_PTR = (int *) this->S2L;

            InputImagePtr = static_cast<segPrecisionTYPE *>(this->InputImage->data);
            InputImagePtr = &InputImagePtr[Multispec*this->numElements];
            segPrecisionTYPE tempsum=(segPrecisionTYPE)0.0;
            segPrecisionTYPE SumPriors=(segPrecisionTYPE)0.0;
            // First, it estimates the weights tempsum and SumPriors.
            if(OutliernessFlag)
            {
                if(BiasField!=NULL)
                {
                    BiasFieldPTR= &BiasField[Multispec*numel_masked];
                    for (int i=0; i<numel_masked; i++, ExpectationPTR++,OutliernessPTR++,BiasFieldPTR++,S2L_PTR++)
                    {
                        segPrecisionTYPE current_value=(*ExpectationPTR)*(*OutliernessPTR)*(InputImagePtr[(*S2L_PTR)]+(*BiasFieldPTR));
                        if(current_value==current_value)
                        {
                            tempsum+=current_value;
                            SumPriors+=(*ExpectationPTR)*(*OutliernessPTR);
                        }
                    }
                }
                else
                {
                    for (int i=0; i<numel_masked; i++, ExpectationPTR++,S2L_PTR++)
                    {
                        segPrecisionTYPE current_value=(*ExpectationPTR)*(*OutliernessPTR)*(InputImagePtr[(*S2L_PTR)]);
                        if(current_value==current_value)
                        {
                            tempsum+=current_value;
                            SumPriors+=(*ExpectationPTR)*(*OutliernessPTR);
                        }
                    }
                }
            }
            else
            {
                if(BiasField!=NULL)
                {
                    BiasFieldPTR= &BiasField[Multispec*numel_masked];
                    for (int i=0; i<numel_masked; i++, ExpectationPTR++,BiasFieldPTR++,S2L_PTR++)
                    {
                        segPrecisionTYPE current_value=(*ExpectationPTR)*(InputImagePtr[(*S2L_PTR)]+(*BiasFieldPTR));
                        if(current_value==current_value)
                        {
                            tempsum+=current_value;
                            SumPriors+=(*ExpectationPTR);
                        }
                    }
                }
                else
                {
                    for (int i=0; i<numel_masked; i++, ExpectationPTR++,S2L_PTR++)
                    {
                        segPrecisionTYPE current_value=(*ExpectationPTR)*(InputImagePtr[(*S2L_PTR)]);
                        if(current_value==current_value)
                        {
                            tempsum+=current_value;
                            SumPriors+=(*ExpectationPTR);
                        }
                    }
                }

            }

            // Second, estimates the actual mean M, the variance V
            if(SumPriors==SumPriors && SumPriors>0)
            {
                if(this->mapStatus && this->MAP_M!=NULL)
                {
                    // The mean M
                    M[cl*(this->nu)+Multispec]=(tempsum/SumPriors/powf(V[cl*(this->nu)+Multispec],2)+this->MAP_M[cl*(this->nu)+Multispec]/powf(this->MAP_V[cl*(this->nu)+Multispec],2))/(1/powf(V[cl*(this->nu)+Multispec],2)+1/powf(this->MAP_V[cl*(this->nu)+Multispec],2));
                }
                else
                {
                    M[cl*(this->nu)+Multispec]=tempsum/SumPriors;

                }

                // The Covariances V, which require looping through all the multispectral chanels again.
                for(long Multispec2=Multispec; Multispec2<this->nu; Multispec2++)
                {
                    S2L_PTR = (int *) this->S2L;

                    InputImagePtr = static_cast<segPrecisionTYPE *>(this->InputImage->data);
                    InputImagePtr =&InputImagePtr[Multispec*this->numElements];

                    T1_PTR2 = static_cast<segPrecisionTYPE *>(this->InputImage->data);
                    T1_PTR2 =&T1_PTR2[Multispec2*this->numElements];
                    segPrecisionTYPE tmpM=this->M[cl*this->nu+Multispec];
                    segPrecisionTYPE tmpM2=this->M[cl*this->nu+Multispec2];
                    tempsum=0;
                    ExpectationPTR=&Expec[Expec_offset[cl]];
                    OutliernessPTR=(segPrecisionTYPE *) &Outlierness[Expec_offset[cl]];

                    if(BiasField!=NULL)
                    {
                        BiasFieldPTR=&BiasField[Multispec*numel_masked];
                        BiasField_PTR2=&BiasField[Multispec2*numel_masked];
                        if(OutliernessFlag)
                        {
                            for (int i=0; i<numel_masked; i++,ExpectationPTR++,BiasFieldPTR++,OutliernessPTR++,BiasField_PTR2++,S2L_PTR++)
                            {

                                segPrecisionTYPE currentValue=(*ExpectationPTR) * (*OutliernessPTR)*(InputImagePtr[(*S2L_PTR)]+(*BiasFieldPTR)-tmpM) * (T1_PTR2[(*S2L_PTR)]+(*BiasField_PTR2)-tmpM2);
                                if(currentValue==currentValue)
                                {
                                    tempsum+=currentValue;
                                }
                            }
                        }
                        else
                        {
                            for (int i=0; i<numel_masked; i++,ExpectationPTR++,BiasFieldPTR++,BiasField_PTR2++,S2L_PTR++)
                            {
                                segPrecisionTYPE currentValue=(*ExpectationPTR) * (InputImagePtr[(*S2L_PTR)]+(*BiasFieldPTR)-tmpM) * (T1_PTR2[(*S2L_PTR)]+(*BiasField_PTR2)-tmpM2);
                                if(currentValue==currentValue)
                                {
                                    tempsum+=currentValue;
                                }
                            }
                        }
                    }
                    else
                    {
                        for (int i=0; i<numel_masked; i++,ExpectationPTR++,S2L_PTR++)
                        {
                            segPrecisionTYPE current_vaue=(*ExpectationPTR) * (InputImagePtr[(*S2L_PTR)]-tmpM) * (T1_PTR2[(*S2L_PTR)]-tmpM2);
                            if(current_vaue==current_vaue)
                            {
                                tempsum+=current_vaue;
                            }
                        }

                    }
                    if( (tempsum/SumPriors>0) && SumPriors>0  && (!isnan(tempsum/SumPriors)))
                    {
                        // assign tempsum/SumPriors to the uper triangular part
                        V[cl*this->nu*this->nu+Multispec+Multispec2*this->nu]=tempsum/SumPriors;
                        //if the value is not in a diagonal
                        if(Multispec2!=Multispec)
                        {
                            // then copy to the botom triangular part
                            V[cl*this->nu*this->nu+Multispec2+Multispec*this->nu]=V[cl*this->nu*this->nu+Multispec+Multispec2*this->nu];
                            // and regularise the value by dividing by the reg_factor.
                            V[cl*this->nu*this->nu+Multispec+Multispec2*this->nu]/=reg_factor;
                            V[cl*this->nu*this->nu+Multispec2+Multispec*this->nu]/=reg_factor;
                        }
                    }
                }
            }
        }
    }

    // This section is only for printing.
    if(verbose>0)
    {
        for (int cl=0; cl<num_class; cl++)
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
/// @brief Optimisation function taking care of the Expectation step and also the outlierness estimation.
///
/// Estimates the responsabilities given the current estimates of the mean and covariance matrix parameters.
/// Also estimates the outlierness, given the Mahalanobis threshold.
///
void seg_EM::RunExpectation()
{


    if(this->ratio<(segPrecisionTYPE)(this->outliernessRatio) && this->outliernessStatus)
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

    int numel_masked=this->numElementsMasked;
    int num_class=this->numberOfClasses;
    bool OutliernessFlag=(this->OutliernessUSE==NULL)?0:1;
    segPrecisionTYPE inv_v [maxNumbClass*maxMultispectalSize*maxMultispectalSize]= {0.0f};
    segPrecisionTYPE inv_sqrt_V_2pi [maxNumbClass]= {0.0f};

    int Expec_offset [maxNumbClass]= {0};

    // for each class, do
    for (int cl=0; cl<num_class; cl++)
    {
        Expec_offset[cl]=(int) cl*numel_masked;
        // if multimodal, do
        if(this->nu>1)
        {
            seg_Matrix <double> Vmat(this->nu,this->nu);
            for (long j2=0; j2<this->nu; j2++)
            {
                for (long i2=j2; i2<this->nu; i2++)
                {
                    Vmat.setvalue(i2,j2,(double)(this->V[i2+j2*this->nu+cl*this->nu*this->nu]));
                    Vmat.setvalue(j2,i2,(double)(this->V[i2+j2*this->nu+cl*this->nu*this->nu]));
                }
            }
            // Get the Gaussian normaliser
            inv_sqrt_V_2pi[cl]=1/(sqrtf(2*M_PI*Vmat.determinant()));
            if (this->verbose_level>1)
            {
                cout<<endl<<"inv_sqrt_V_2pi["<< cl <<"]= "<< inv_sqrt_V_2pi[cl] << endl;
                flush(cout);
            }
            // Get the inverted covariance matrix
            Vmat.invert();
            double covarianceValue=0.0f;
            bool success;
            // Print if in debug, i.e. verbose_level==2
            if (this->verbose_level>1)
            {
                cout<<"inv_V["<< cl <<"]= ";
                flush(cout);
            }
            for (long j2=0; j2<this->nu; j2++)
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
                    // copy data from Vmat to inv_v
                    Vmat.getvalue(i2,j2,covarianceValue,success);
                    inv_v[i2+j2*this->nu+cl*this->nu*this->nu]=(segPrecisionTYPE)(covarianceValue);
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
        // else, just get the gaussian normaliser and the inverse covariance determinant
        else
        {
            inv_sqrt_V_2pi[cl]=1/(sqrtf(2*M_PI*V[cl]));
            inv_v[cl]=1/V[cl];
        }
    }
    this->loglik=0;
    segPrecisionTYPE logliktmp=0.0f;


    segPrecisionTYPE * ExpectationTmpPTR = (segPrecisionTYPE *) this->Expec;
    segPrecisionTYPE * OutliernessTmpPTR = (segPrecisionTYPE *) this->Outlierness;
    segPrecisionTYPE * InputImageTmpPtr = static_cast<segPrecisionTYPE *>(this->InputImage->data);
    segPrecisionTYPE * BiasFieldTmpPTR= (segPrecisionTYPE *) this->BiasField;

#ifdef _OPENMP
    segPrecisionTYPE * loglikthread = new segPrecisionTYPE [omp_get_max_threads()]();
    for(long i=0; i<(long)omp_get_max_threads(); i++)
        loglikthread[i]=0;

#pragma omp parallel for shared(ExpectationTmpPTR,loglikthread,InputImageTmpPtr,BiasFieldTmpPTR,OutliernessTmpPTR,IterPrior)
#endif
    // Now that we have the Gaussian normaliser and the inverted determinant of the covariance, estimate the Gaussian PDF. We do that independently per voxel.
    for (int i=0; i<numel_masked; i++)
    {
        segPrecisionTYPE * T1_PTR = (segPrecisionTYPE *)(InputImageTmpPtr);
        segPrecisionTYPE T1_Bias_corr[maxMultispectalSize];
        segPrecisionTYPE SumExpec=0.0f;

        // For each modality, estimate the bias field corrected intensity first
        for(long Multispec=0; Multispec<this->nu; Multispec++)
            T1_Bias_corr[Multispec]=(BiasFieldTmpPTR!=NULL)?(T1_PTR[this->S2L[i]+Multispec*this->numElements] + BiasFieldTmpPTR[i+Multispec*numel_masked]):(T1_PTR[this->S2L[i]+Multispec*this->numElements]);

        // Then, for each class and for each modality, iterate over all modalities and subtract their (Y-\Mu)/inv(|\Sigma|) in order to get the Mahalanobis distance.
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

            // If the outlierness threshold is used, estimate here the new outlierness values. Note that this value is different per class. The 0.01 are there for stability reasons.
            if(OutliernessFlag)
            {
                segPrecisionTYPE outvalue=(expf(mahal)+0.01)/(expf(mahal)+expf(-0.5*(this->outliernessThreshold*this->outliernessThreshold))+0.01);
                OutliernessTmpPTR[i+Expec_offset[cl]]=outvalue;
            }
            // Update the Expectation by calculating, Prior*GaussianPDF
            ExpectationTmpPTR[i+Expec_offset[cl]]=IterPrior[i+Expec_offset[cl]] * ( inv_sqrt_V_2pi[cl] * expf(mahal) ) ;
            // Update the normaliser
            SumExpec+=ExpectationTmpPTR[i+Expec_offset[cl]];
        }

        // If something went wrong, just set the expectation to 1/K
        if (SumExpec<=0.0 || SumExpec!=SumExpec)
        {
            for (int cl=0; cl<num_class; cl++)
            {
                ExpectationTmpPTR[i+Expec_offset[cl]]=(segPrecisionTYPE)(1)/(segPrecisionTYPE)(num_class);
            }

        }
        // If it worked, then normalise the expectations, thus obtaining a per voxel responsability
        else
        {

            for (int cl=0; cl<num_class; cl++)
            {
                ExpectationTmpPTR[i+Expec_offset[cl]]=ExpectationTmpPTR[i+Expec_offset[cl]]/SumExpec;
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

    // update the log likelihood
    loglik=logliktmp;

    return;

}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Optimisation function taking care of the MRF optimisation step
///
/// This function will run the MRF in 2D or 3D depending on the clique structure. A ND version is coming soon, providing a temporal MRF.
///
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
/// @brief Optimisation function taking care of the 2D MRF
///
/// Estimates the 2D MRF energy given the MRFTransitionMatrix. This is similar to this->RunMRF3D() but in 2D. \sa this->RunMRF3D()
///
void seg_EM::RunMRF2D()
{

    int numelmasked=this->numElementsMasked;
    int numclass=this->numberOfClasses;
    segPrecisionTYPE * G =this->MRFTransitionMatrix;

    // if the MRF optimisation is ON
    if(this->mrfStatus)
    {
        // Define all the pointers and varibles
        segPrecisionTYPE * MRFpriorPtr = (segPrecisionTYPE *) this->MRF;
        int * Long_2_Short_IndicesPtr = (int *) this->L2S;
        int col_size, indexCentre, indexWest, indexEast, indexSouth, indexNorth;
        int ix, iy,maxiy, maxix, neighbourclass;
        segPrecisionTYPE Sum_Temp_MRF_Class_Expect;
        segPrecisionTYPE Gplane[maxNumbClass];
        segPrecisionTYPE Temp_MRF_Class_Expect[maxNumbClass];
        register int currclass;
        unsigned int numelmasked_currclass_shift[maxNumbClass];
        col_size = (int)(this->nx);
        maxix = (int)(this->nx);
        maxiy = (int)(this->ny);


        if(verbose_level>0)
        {
            cout << "Optimising MRF"<<endl;
            flush(cout);
        }

        // precompute the index shifts between modalities (performance reasons)
        for(int i=0; i<numclass; i++)
        {
            numelmasked_currclass_shift[i]=i*numelmasked;
        }

        // As it is 2D, iterate over all Y and X's
        int curr_short_centreindex;
        for (iy=1; iy<maxiy-1; iy++)
        {
            // This updates the indexCentre to the correct column
            indexCentre=(col_size*iy);
            for (ix=1; ix<maxix-1; ix++)
            {
                // indexcentre is incremented before because ix starts at 1
                indexCentre++;
                Sum_Temp_MRF_Class_Expect = 0;
                curr_short_centreindex=Long_2_Short_IndicesPtr[indexCentre];

                if (curr_short_centreindex>=0)
                {
                    // Get the index of all the 4 neighbours
                    indexWest=Long_2_Short_IndicesPtr[indexCentre-col_size]>-1?Long_2_Short_IndicesPtr[indexCentre-col_size]:0;
                    indexEast=Long_2_Short_IndicesPtr[indexCentre+col_size]>-1?Long_2_Short_IndicesPtr[indexCentre+col_size]:0;
                    indexNorth=Long_2_Short_IndicesPtr[indexCentre-1]>-1?Long_2_Short_IndicesPtr[indexCentre-1]:0;
                    indexSouth=Long_2_Short_IndicesPtr[indexCentre+1]>-1?Long_2_Short_IndicesPtr[indexCentre+1]:0;
                    for (currclass=0; currclass<numclass; currclass++)
                    {
                        // Get the sum of the expectations for each class over all the neighbours
                        Gplane[currclass] = 0.0;
                        Temp_MRF_Class_Expect[currclass] = 0.0;
                        Gplane[currclass]+=this->Expec[indexWest];
                        Gplane[currclass]+=this->Expec[indexEast];
                        Gplane[currclass]+=this->Expec[indexNorth];
                        Gplane[currclass]+=this->Expec[indexSouth];
                        // increment the indexes to shift them to the next class
                        if(currclass<numclass)
                        {
                            indexWest+=numelmasked;
                            indexEast+=numelmasked;
                            indexNorth+=numelmasked;
                            indexSouth+=numelmasked;
                        }
                    }
                    // After you have the probabilities, estimate exp(-Beta*U_MRF), with U_MRF beeing the sum over all classes of G*Gplane
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
                        // Also estimate the normaliser
                        Sum_Temp_MRF_Class_Expect += Temp_MRF_Class_Expect[currclass];
                    }
                    // Normalise the MRF prob using the MRF_Exp/Sum_MRF_Exp
                    for (currclass=0; currclass<numclass; currclass++)
                    {
                        MRFpriorPtr[curr_short_centreindex+numelmasked_currclass_shift[currclass]]=(Temp_MRF_Class_Expect[currclass]/Sum_Temp_MRF_Class_Expect);
                    }
                }
            }
        }

    }
    return;

}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Optimisation function taking care of the 3D MRF
///
/// Estimates the 3D MRF energy given the MRFTransitionMatrix. This is similar to this->RunMRF2D() but in 3D. \sa this->RunMRF2D()
///
void seg_EM::RunMRF3D()
{
    int numelmasked=this->numElementsMasked;
    int numclass=this->numberOfClasses;

    segPrecisionTYPE * G =this->MRFTransitionMatrix;
    segPrecisionTYPE * H =this->MRFTransitionMatrix;

    // if the MRF optimisation is ON
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

        // precompute the index shifts between modalities (performance reasons)
        for(int i=0; i<numclass; i++)
        {
            numelmasked_currclass_shift[i]=i*numelmasked;

        }

        segPrecisionTYPE * ExpectationTmpPTR = (segPrecisionTYPE *) this->Expec;

#ifdef _OPENMP
#pragma omp parallel for shared(ExpectationTmpPTR,MRFpriorPtr,Long_2_Short_IndicesPtr)
#endif
        // As it is 3D, iterate over all Z, Y and X's
        for (int iz=1; iz<maxiz-1; iz++)
        {
            // Define all the pointers and varibles
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
                    // indexcentre is incremented before because ix starts at 1
                    indexCentre++;
                    Sum_Temp_MRF_Class_Expect = 0;
                    int curr_short_centreindex=Long_2_Short_IndicesPtr[indexCentre];
                    if (curr_short_centreindex>=0)
                    {
                        // Get the index of all the 6 neighbours
                        indexWest=Long_2_Short_IndicesPtr[indexCentre-col_size]>-1?Long_2_Short_IndicesPtr[indexCentre-col_size]:0;
                        indexEast=Long_2_Short_IndicesPtr[indexCentre+col_size]>-1?Long_2_Short_IndicesPtr[indexCentre+col_size]:0;
                        indexNorth=Long_2_Short_IndicesPtr[indexCentre-1]>-1?Long_2_Short_IndicesPtr[indexCentre-1]:0;
                        indexSouth=Long_2_Short_IndicesPtr[indexCentre+1]>-1?Long_2_Short_IndicesPtr[indexCentre+1]:0;
                        indexBottom=Long_2_Short_IndicesPtr[indexCentre+plane_size]>-1?Long_2_Short_IndicesPtr[indexCentre+plane_size]:0;
                        indexTop=Long_2_Short_IndicesPtr[indexCentre-plane_size]>-1?Long_2_Short_IndicesPtr[indexCentre-plane_size]:0;
                        for (currclass=0; currclass<numclass; currclass++)
                        {
                            // Get the sum of the expectations for each class over all the neighbours
                            Gplane[currclass] = 0.0;
                            Hplane[currclass] = 0.0;
                            Temp_MRF_Class_Expect[currclass] = 0.0;
                            Gplane[currclass]+=ExpectationTmpPTR[indexWest];
                            Gplane[currclass]+=ExpectationTmpPTR[indexEast];
                            Gplane[currclass]+=ExpectationTmpPTR[indexNorth];
                            Gplane[currclass]+=ExpectationTmpPTR[indexSouth];
                            Hplane[currclass]+=ExpectationTmpPTR[indexTop];
                            Hplane[currclass]+=ExpectationTmpPTR[indexBottom];
                            // increment the indexes to shift them to the next class

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
                        // After you have the probabilities, estimate exp(-Beta*U_MRF), with U_MRF beeing the sum over all classes of ( G*Gplane + H*Hplane )

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
                            // Estimate the normaliser
                            Sum_Temp_MRF_Class_Expect+=Temp_MRF_Class_Expect[currclass];
                        }
                        // Normalise the MRF prob using the MRF_Exp/Sum_MRF_Exp
                        for (currclass=0; currclass<numclass; currclass++)
                        {
                            MRFpriorPtr[curr_short_centreindex+numelmasked_currclass_shift[currclass]]=(Temp_MRF_Class_Expect[currclass]/Sum_Temp_MRF_Class_Expect);
                        }
                    }
                }
            }
        }
    }
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Function taking care of the bias field optimisation
///
/// This function calls the this->RunBiasField3D() and this->RunBiasField2D(), depending on the dimention of the data.
///

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
/// @brief Optimisation function taking care of the 3D Bias Field
///
/// Estimates the 3D MRF energy given the numbe of basis functions. This is similar to this->RunBiasField2D() but in 3D. \sa this->RunBiasField2D()
///

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
            cout<< "Assuming independent bias-fields between time points" << endl;
        }
        flush(cout);
    }

    int reduxfactor=reduxFactorForBias;
    segPrecisionTYPE * sampledData = static_cast<segPrecisionTYPE *>(this->InputImage->data);

    // Get number of basis functions
    int UsedBasisFunctions=(int)(((float)(this->biasFieldOrder)+1.0f) * ((float)(this->biasFieldOrder)+2.0f)/2.0f *((float)(this->biasFieldOrder)+3.0f)/3.0f);
    // Precompute Powers depending on the current BiasOrder. The power order to sum to this->biasFieldOrder max.
    int PowerOrder [(int)((((float)(maxAllowedBCPowerOrder)+1.0f) * ((float)(maxAllowedBCPowerOrder)+2.0f)/2.0f *((float)(maxAllowedBCPowerOrder)+3.0f)/3.0f))*3]= {0};
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
    for(long i=0; i<maxNumbClass; i++)
    {
        invV[i]=0;
        currM[i]=0;
    }


    // Solve the problem separetly per modality, as the Bias Field is assumed to be indpendent between modalities.
    for(long multispec=0; multispec<this->nu; multispec++)
    {

        // Get pointers to the current modality
        sampledData = static_cast<segPrecisionTYPE *>(this->InputImage->data);
        sampledData = &sampledData[multispec*this->numElements];
        // Precompute the M and V  inverses
        for(int i=0; i<this->numberOfClasses; i++)
        {
            invV[i]=1.0f/(V[i*this->nu*this->nu+multispec+multispec*this->nu]);
            currM[i]=M[i*this->nu+multispec];
        }

        //
        segPrecisionTYPE AWA [(int)((((float)(maxAllowedBCPowerOrder)+1.0f) * ((float)(maxAllowedBCPowerOrder)+2.0f)/2.0f *((float)(maxAllowedBCPowerOrder)+3.0f)/3.0f))*(int)((((float)(maxAllowedBCPowerOrder)+1.0f) * ((float)(maxAllowedBCPowerOrder)+2.0f)/2.0f *((float)(maxAllowedBCPowerOrder)+3.0f)/3.0f))]= {0.0f};
        segPrecisionTYPE WR [(int)((((float)(maxAllowedBCPowerOrder)+1.0f) * ((float)(maxAllowedBCPowerOrder)+2.0f)/2.0f *((float)(maxAllowedBCPowerOrder)+3.0f)/3.0f))]= {0.0f};
        segPrecisionTYPE FinalCoefs [(int)((((float)(maxAllowedBCPowerOrder)+1.0f) * ((float)(maxAllowedBCPowerOrder)+2.0f)/2.0f *((float)(maxAllowedBCPowerOrder)+3.0f)/3.0f))]= {0.0f};


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
        if(this->numelBias==0)
        {
            for (int iz=0; iz<this->nz; iz+=reduxfactor)
            {
                for (int iy=0; iy<this->ny; iy+=reduxfactor)
                {
                    for (int ix=0; ix<this->nx; ix+=reduxfactor)
                    {
                        linearindexes = iz*this->nx*this->ny + iy*this->nx + ix;

                        if(this->L2S[linearindexes]>=0)
                        {
                            samplecount++;
                        }
                    }
                }
            }
            if(verbose_level>0)
            {
                cout << "Number of samples for BiasField = " << samplecount<<", with "<<UsedBasisFunctions<<" basis functions\n";
                flush(cout);
            }
            this->numelBias=samplecount;
        }
        else
        {
            samplecount=this->numelBias;
        }



        // CALC A'WA

        segPrecisionTYPE * W;
        try
        {
            W = new segPrecisionTYPE [samplecount] ();
        }
        catch (std::bad_alloc& ba)
        {
            std::cerr << "ERROR:Low memory, could not allocate Q (" <<samplecount<<" float elements): "<< ba.what() <<std::endl;
            exit(1);
        }

        segPrecisionTYPE W_tmp=0;
        currshortindex=0;
        int Windex=0;
        // Calc W from the A'WA part of the equations (Van Leemput 1999 eq 7)
        for (int iz=0; iz<maxiz; iz+=reduxfactor)
        {
            for (int iy=0; iy<maxiy; iy+=reduxfactor)
            {
                for (int ix=0; ix<maxix; ix+=reduxfactor)
                {
                    currshortindex=this->L2S[iz*plane_size+iy*col_size+ix];
                    if(currshortindex>=0)
                    {
                        W_tmp=0;
                        if(Outlierness==NULL)
                        {
                            for(int j=0; j<this->numberOfClasses; j++)
                            {
                                W_tmp+=Expec[currshortindex+this->numElementsMasked*j]*invV[j];
                            }
                        }
                        else
                        {
                            for(int j=0; j<this->numberOfClasses; j++)
                            {
                                W_tmp+=Expec[currshortindex+this->numElementsMasked*j]*Outlierness[currshortindex+this->numElementsMasked*j]*invV[j];
                            }

                        }
                        // This is the actual W. Note that W is stored as a vector, as W is a diagonal matrix
                        W[Windex]=W_tmp;
                        Windex++;
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

        segPrecisionTYPE * Basis;
        try
        {
            Basis= new segPrecisionTYPE[UsedBasisFunctions]();
        }
        catch (std::bad_alloc& ba)
        {
            std::cerr << "ERROR:Low memory, could not allocate Basis (" <<UsedBasisFunctions<<" float elements): "<< ba.what() <<std::endl;
            exit(1);
        }

        segPrecisionTYPE xpos=0.0f;
        segPrecisionTYPE ypos=0.0f;
        segPrecisionTYPE zpos=0.0f;
        int x_bias_index_shift=0;
        int y_bias_index_shift=1;
        int z_bias_index_shift=2;
        segPrecisionTYPE currentW=0.0f;
        // Calc A'WA (Van Leemput 1999 eq 7)
        Windex=0;
        segPrecisionTYPE * Basisptr1= (segPrecisionTYPE *) Basis;
        segPrecisionTYPE * Basisptr2= (segPrecisionTYPE *) Basis;
        segPrecisionTYPE * AWAptr= (segPrecisionTYPE *) AWA;
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
                        currentW=W[Windex];
                        xpos=(((segPrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                        ypos=(((segPrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                        zpos=(((segPrecisionTYPE)iz-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z);
                        x_bias_index_shift=0;
                        y_bias_index_shift=1;
                        z_bias_index_shift=2;
                        for(int j2=0; j2<UsedBasisFunctions; j2++,x_bias_index_shift+=3,y_bias_index_shift+=3,z_bias_index_shift+=3, Basisptr1++)
                        {
                            // Because Powerorder is always int, use a special power function (faster) to estimate the basis.
                            *Basisptr1=(pow_int(xpos,PowerOrder[x_bias_index_shift])*pow_int(ypos,PowerOrder[y_bias_index_shift])*pow_int(zpos,PowerOrder[z_bias_index_shift]));
                        }


                        Basisptr1= (segPrecisionTYPE *) Basis;
                        AWAptr= (segPrecisionTYPE *) AWA;
                        for(int j2=0; j2<UsedBasisFunctions; j2++, Basisptr1++)
                        {
                            Basisptr2= &Basis[j2];
                            AWAptr= &AWA[j2+j2*UsedBasisFunctions];
                            for(int i2=j2; i2<UsedBasisFunctions; i2++, AWAptr++, Basisptr2++)
                            {
                                // Estimate A'WA, i.e. Basis*W*Basis
                                (*AWAptr)+=(*Basisptr2)*(currentW)*(*Basisptr1);
                            }
                        }
                        Windex++;
                    }
                }
            }
        }

        // Transfer from the AWA matrix to an actual matrix structure.
        seg_Matrix <double> RealAWA(UsedBasisFunctions,UsedBasisFunctions);
        for(int j2=0; j2<UsedBasisFunctions; j2++)
        {
            for(int i2=j2; i2<UsedBasisFunctions; i2++)
            {
                RealAWA.setvalue(i2,j2,(double)(AWA[i2+j2*UsedBasisFunctions]));
                RealAWA.setvalue(j2,i2,(double)(AWA[i2+j2*UsedBasisFunctions]));
            }
        }

        // Invert AWA
        seg_Matrix <double> RealAWA_inv(UsedBasisFunctions);
        RealAWA_inv.copymatrix(RealAWA);
        RealAWA_inv.invert();

        // This is only from printing reasons, it does nothing to the vectors
        if(verbose_level>1)
        {
            seg_Matrix <double> RealAWA_test(UsedBasisFunctions);
            RealAWA_test.settoproduct(RealAWA,RealAWA_inv);
            RealAWA_test.comparetoidentity();
        }



        // CALC MATRIX WR

        //Precompute W (Van Leemput 1999 eq 7)
        segPrecisionTYPE Wi;
        segPrecisionTYPE Wij;
        segPrecisionTYPE Yest;
        segPrecisionTYPE Ysum;
        Windex=0;
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
                        for(int j=0; j<this->numberOfClasses; j++)
                        {
                            segPrecisionTYPE tmpexpec = (segPrecisionTYPE)Expec[currshortindex+this->numElementsMasked*j];
                            Wij=tmpexpec*(invV[j]);
                            Wi+=Wij;
                            Yest+=Wij*(currM[j]);
                            Ysum+=Wij;
                        }
                        W[Windex]=Wi*(sampledData[linearindexes]-(Yest/Ysum));
                        Windex++;
                    }
                }
            }
        }

        //Compute WR (Van Leemput 1999 eq 7)
        for(int bfindex=0; bfindex<UsedBasisFunctions; bfindex++)
        {
            int tempvarindex2=0;
            int linearindexes2=0;
            segPrecisionTYPE * WRptr=&WR[bfindex];
            *WRptr=0;
            for (int iz2=0; iz2<maxiz; iz2+=reduxfactor)
            {
                for (int iy2=0; iy2<maxiy; iy2+=reduxfactor)
                {
                    for (int ix2=0; ix2<maxix; ix2+=reduxfactor)
                    {
                        linearindexes2=(iz2)*(this->nx)*(this->ny)+(iy2)*(this->nx)+ix2;
                        if(this->L2S[linearindexes2]>=0)
                        {
                            *WRptr+=pow_int((((segPrecisionTYPE)ix2-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x),PowerOrder[0+bfindex*3])*
                                    pow_int((((segPrecisionTYPE)iy2-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y),PowerOrder[1+bfindex*3])*
                                    pow_int((((segPrecisionTYPE)iz2-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z),PowerOrder[2+bfindex*3])*
                                    W[tempvarindex2];


                            tempvarindex2++;
                        }
                    }
                }
            }
            // if nan, set to 1
            if(*WRptr!=*WRptr)
            {
                *WRptr=1;
            }
        }
        // Copy data from WR to the matrix structure
        seg_Matrix <double> RealWR(UsedBasisFunctions,1);

        for(int i2=0; i2<UsedBasisFunctions; i2++)
        {
            RealWR.setvalue(i2,0,(double)(WR[i2]));
        }
        // Get coeficients by multiplying RealAWA_inv with RealWR
        seg_Matrix <double> RealCoefs(UsedBasisFunctions,1);
        RealCoefs.settoproduct(RealAWA_inv,RealWR);
        if(verbose_level>1)
        {
            cout << "C= " << endl;
            RealCoefs.dumpmatrix();
        }
        double coefValue=0.0f;
        bool success;
        for(int i2=0; i2<UsedBasisFunctions; i2++)
        {
            RealCoefs.getvalue(i2,0,coefValue,success);
            FinalCoefs[i2]=(segPrecisionTYPE)(coefValue);
        }


        // Now that we have the Coefs, we can get the actual bias field by multiplying with the Basis functions.
        for (int iz=0; iz<maxiz; iz++)
        {
            segPrecisionTYPE currxpower[maxAllowedBCPowerOrder]={0};
            segPrecisionTYPE currypower[maxAllowedBCPowerOrder]={0};
            segPrecisionTYPE currzpower[maxAllowedBCPowerOrder]={0};
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
                                    tmpbiasfield-=FinalCoefs[ind2]*currxpower[xorder]*currypower[yorder]*currzpower[zorder];
                                    ind2++;
                                }
                            }
                        }
                        this->BiasField[currshortindex+multispec*this->numElementsMasked]=tmpbiasfield;
                    }
                }
            }
        }

        // Save the coeficients back
        for( int i=0; i<UsedBasisFunctions; i++)
        {
            this->biasFieldCoeficients[i+multispec*UsedBasisFunctions]=FinalCoefs[i];
        }

        if(Basis!=NULL)
        {
            delete [] Basis;
        }
        if(W!=NULL)
        {
            delete [] W;
        }
    }
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Optimisation function taking care of the 2D Bias Field
///
/// Estimates the 2D MRF energy given the numbe of basis functions. This is similar to this->RunBiasField3D() but in 2D. \sa this->RunBiasField3D()
///
void seg_EM::RunBiasField2D()
{
    // For more comments on this function, see the RunBiasField3D() function, as it is the same but with only a 2D basis.
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
    long nrOfClasses = this->numberOfClasses;
    //nrOfClasses = 1;
    //int nrOfClasses = non_PV_numclass;
    long TotalLength = this->numElementsMasked;
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
    for(long i=0; i<maxNumbClass; i++)
    {
        invV[i]=0;
        currM[i]=0;
    }


    for(long multispec=0; multispec<this->nu; multispec++)
    {
        sampledData = static_cast<segPrecisionTYPE *>(this->InputImage->data);
        sampledData = &sampledData[multispec*this->numElements];
        // Precompute the M and V  inverses
        for(long i=0; i<nrOfClasses; i++)
        {
            invV[i]=1.0f/(V[i*this->nu*this->nu+multispec+multispec*this->nu]);
            currM[i]=M[i*this->nu+multispec];
        }
        segPrecisionTYPE AWA [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2))/2*((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2)]= {0.0f};
        segPrecisionTYPE WR [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2))/2]= {0.0f};
        segPrecisionTYPE FinalCoefs [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2))/2]= {0.0f};


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
        if(this->numelBias==0)
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
            this->numelBias=samplecount;
        }
        else
        {
            samplecount=this->numelBias;
        }



        // CALC MATRIX A

        // Calc W (Van Leemput 1999 eq 7)
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
        segPrecisionTYPE * AWAptr= (segPrecisionTYPE *) AWA;
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
                    AWAptr= (segPrecisionTYPE *) AWA;
                    for(int j2=0; j2<UsedBasisFunctions; j2++, Basisptr1++)
                    {
                        Basisptr2= &Basis[j2];
                        AWAptr= &AWA[j2+j2*UsedBasisFunctions];
                        for(int i2=j2; i2<UsedBasisFunctions; i2++, AWAptr++, Basisptr2++)
                        {
                            (*AWAptr)+=(*Basisptr2)*(current_Tempvar)*(*Basisptr1);
                        }
                    }
                    tempvarindex++;
                }
            }
        }

        seg_Matrix <double> RealAWA(UsedBasisFunctions,UsedBasisFunctions);

        for(int j2=0; j2<UsedBasisFunctions; j2++)
        {
            for(int i2=j2; i2<UsedBasisFunctions; i2++)
            {
                RealAWA.setvalue(i2,j2,(double)(AWA[i2+j2*UsedBasisFunctions]));
                RealAWA.setvalue(j2,i2,(double)(AWA[i2+j2*UsedBasisFunctions]));
            }
        }

        seg_Matrix <double> RealAWA_inv(UsedBasisFunctions);
        RealAWA_inv.copymatrix(RealAWA);
        RealAWA_inv.invert();

        if(verbose_level>1)
        {
            seg_Matrix <double> RealAWA_test(UsedBasisFunctions);
            RealAWA_test.settoproduct(RealAWA,RealAWA_inv);
            RealAWA_test.comparetoidentity();
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
            WR[i2]=0;
            for (int iy=0; iy<maxiy; iy+=reduxfactor)
            {

                for (int ix=0; ix<maxix; ix+=reduxfactor)
                {
                    linearindexes=(iy)*(this->nx)+ix;
                    currshortindex=this->L2S[linearindexes];
                    if(currshortindex>=0)
                    {
                        WR[i2]+=pow_int((((segPrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x),PowerOrder[0+i2*2])*
                                pow_int((((segPrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y),PowerOrder[1+i2*2])*
                                Tempvar[tempvarindex];

                        if(WR[i2]!=WR[i2])
                        {
                            WR[i2]=1;
                        }
                        tempvarindex++;
                    }
                }
            }
        }

        seg_Matrix <double> RealWR(UsedBasisFunctions,1);

        for(int i2=0; i2<UsedBasisFunctions; i2++)
        {
            RealWR.setvalue(i2,0,(double)(WR[i2]));
        }

        seg_Matrix <double> RealFinalCoefs(UsedBasisFunctions,1);

        RealFinalCoefs.settoproduct(RealAWA_inv,RealWR);
        if(verbose_level>1)
        {
            cout << "C= " << endl;
            RealFinalCoefs.dumpmatrix();
        }

        double cvalue=0.0f;
        bool success;
        for(int i2=0; i2<UsedBasisFunctions; i2++)
        {
            RealFinalCoefs.getvalue(i2,0,cvalue,success);
            FinalCoefs[i2]=(segPrecisionTYPE)(cvalue);
        }

        segPrecisionTYPE currxpower[maxAllowedBCPowerOrder]={0};
        segPrecisionTYPE currypower[maxAllowedBCPowerOrder]={0};
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
                            tmpbiasfield-=FinalCoefs[ind]*currxpower[xorder]*currypower[yorder];
                            ind++;
                        }
                    }
                    this->BiasField[currshortindex+multispec*this->numElementsMasked]=tmpbiasfield;
                }
            }

        }

        for( int i=0; i<UsedBasisFunctions; i++)
        {
            this->biasFieldCoeficients[i+multispec*UsedBasisFunctions]=FinalCoefs[i];
        }
        delete [] Basis;
        delete [] Tempvar;
    }
}



/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Optimisation function taking care of the Prior relaxation
///
/// Estimates the Gaussian-regularised Derrichlet prior relaxation as described in the AdaPT (Cardoso et al. Neuroimage) paper
///
void seg_EM::RunPriorRelaxation()
{
    if(this->relaxStatus)
    {
        if((int)(this->verbose_level)>(int)(0))
        {
            cout << "Relaxing Priors"<< endl;
        }

        // Copy expectation to the this->ShortPrior and convolve with a gaussian filter
        for(long i=0; i<(this->numberOfClasses*this->numElementsMasked); i++)
        {
            this->ShortPrior[i]=Expec[i];
        }

        int oldTSize=this->CurrSizes->tsize;
        this->CurrSizes->tsize=this->CurrSizes->numclass;
        GaussianFilter4D_cArray(this->ShortPrior,this->S2L,this->L2S,this->relaxGaussKernelSize,this->CurrSizes);
        this->CurrSizes->tsize=oldTSize;

        // Estimate the new priors, i.e. (1-alpha)(Gaussian(Expec)) + (alpha)(Prior)
        segPrecisionTYPE * PriorsPtr = static_cast<segPrecisionTYPE *>(this->Priors->data);
        long currindex=0;
        for(long k=0; k<this->numberOfClasses; k++)
        {
            currindex=k*this->numElementsMasked;
            for(long i=0; i<(this->numElementsMasked); i++, currindex++)
            {
                // gets (1-alpha)(Gaussian(Expec))
                this->ShortPrior[currindex]*=(1-this->relaxFactor);
                // gets  (alpha)(Prior)
                this->ShortPrior[currindex]+=(this->relaxFactor)*PriorsPtr[S2L[i]+k*this->numElements];

            }
        }
    }
    return;
}
#endif
