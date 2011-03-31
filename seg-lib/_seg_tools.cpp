#include "_seg_tools.h"

static union
{
    double d;
    struct
    {
#ifdef LITTLE_ENDIAN
        int j, i;
#else
        int i, j;
#endif
    } n;
} _eco;
#define EXP_A (1048576 /M_LN2)
#define EXP_C 60801
/* use 1512775 for integer version */
/* see text for choice of c values */
#define EXP(y) (_eco.n.i = EXP_A*(y) + (1072693248 - EXP_C), _eco.d)


void Create_diagonal_GH_Nclass(PrecisionTYPE * G,
                               PrecisionTYPE * H,
                               PrecisionTYPE ratio,
                               SEG_PARAM * segment_param)
{
    //WM - neighbours GM and dGM
    int numclass=segment_param->numb_classes;
    for(int i=0; i<segment_param->numb_classes; i++){
        for(int j=0; j<segment_param->numb_classes; j++){
            if(i!=j){
                G[i+j*numclass]=segment_param->MRF_strength;
                H[i+j*numclass]=ratio*G[i+j*numclass];
            }
            else{
                G[i+j*numclass]=0.0;
                H[i+j*numclass]=ratio*G[i+j*numclass];
            }

        }
    }

    //Print G
    if(segment_param->verbose_level>1){
        cout<<"G=" << endl;
        for (int i=0; i<numclass; i++) {
            for (int  j=0; j<numclass; j++) {
                cout<< (PrecisionTYPE)G[i+j*numclass] << '\t';
            }
            cout<< endl;
        }
        cout<< endl << "H=" << endl;
        for (int i=0; i<numclass; i++) {
            for (int j=0; j<numclass; j++) {
                cout<< G[i+j*numclass] << '\t';
            }
            cout<< endl;
        }
        cout<< endl;
    }

    return;
}


void Create_GH_5class(PrecisionTYPE * G,
                      PrecisionTYPE * H,
                      PrecisionTYPE ba,
                      PrecisionTYPE be,
                      PrecisionTYPE ratio,
                      SEG_PARAM * segment_param)
{
    //WM - neighbours GM and dGM
    int numclass=5;
    G[WMclass+WMclass*numclass]=0.0;
    G[GMclass+WMclass*numclass]=be;
    G[CSFclass+WMclass*numclass]=ba;
    G[dGMclass+WMclass*numclass]=be;
    G[iCSFclass+WMclass*numclass]=ba;
    //GM  - neighbours WM and CSF
    G[GMclass+GMclass*numclass]=0.0;
    G[CSFclass+GMclass*numclass]=be;
    G[dGMclass+GMclass*numclass]=ba;
    G[iCSFclass+GMclass*numclass]=ba;
    //CSF - neighbours GM
    G[CSFclass+CSFclass*numclass]=0.0;
    G[dGMclass+CSFclass*numclass]=ba;
    G[iCSFclass+CSFclass*numclass]=ba;
    //dGM - neighbours WM and iCSF
    G[dGMclass+dGMclass*numclass]=0.0;
    G[iCSFclass+dGMclass*numclass]=be;
    //iCSF - neighbours dGM
    G[iCSFclass+iCSFclass*numclass]=0.0;
    int i,j;

    //Copy to lower triangle
    for (j=0; j<numclass; j++) { for (i=1+j; i<numclass; i++) { G[j+i*numclass]=G[i+j*numclass]; } }
    //Copy G to H
    for (int i=0; i<(numclass*numclass); i++) { H[i]=ratio*G[i]; }
    //Print G
    if(segment_param->verbose_level>1){
        cout<<"G=" << endl;
        for (i=0; i<numclass; i++) {
            for (j=0; j<numclass; j++) {
                cout<< (PrecisionTYPE)G[i+j*numclass] << '\t';
            }
            cout<< endl;
        }
        cout<< endl << "H=" << endl;
        for (i=0; i<numclass; i++) {
            for (j=0; j<numclass; j++) {
                cout<< G[i+j*numclass] << '\t';
            }
            cout<< endl;
        }
        cout<< endl;
    }

    return;
}


void Create_GH_7class(PrecisionTYPE * G,
                      PrecisionTYPE * H,
                      PrecisionTYPE ba,
                      PrecisionTYPE be,
                      PrecisionTYPE ratio,
                      SEG_PARAM * segment_param)
{
    //WM - neighbours WMGMpv and dGM
    int numclass=7;
    G[WMclass+WMclass*numclass]=0.0;
    G[GMclass+WMclass*numclass]=ba;
    G[CSFclass+WMclass*numclass]=ba;
    G[dGMclass+WMclass*numclass]=be;
    G[iCSFclass+WMclass*numclass]=ba;
    G[WMGMpvclass+WMclass*numclass]=be;
    G[GMCSFpvclass+WMclass*numclass]=ba;
    //GM  - neighbours WMGMpv and GMCSFpv
    G[GMclass+GMclass*numclass]=0.0;
    G[CSFclass+GMclass*numclass]=ba;
    G[dGMclass+GMclass*numclass]=ba;
    G[iCSFclass+GMclass*numclass]=ba;
    G[WMGMpvclass+GMclass*numclass]=be;
    G[GMCSFpvclass+GMclass*numclass]=be;
    //CSF - neighbours GMCSFpv
    G[CSFclass+CSFclass*numclass]=0.0;
    G[dGMclass+CSFclass*numclass]=ba;
    G[iCSFclass+CSFclass*numclass]=ba;
    G[WMGMpvclass+CSFclass*numclass]=ba;
    G[GMCSFpvclass+CSFclass*numclass]=be;
    //dGM - neighbours WM and iCSF
    G[dGMclass+dGMclass*numclass]=0.0;
    G[iCSFclass+dGMclass*numclass]=be;
    G[WMGMpvclass+dGMclass*numclass]=ba;
    G[GMCSFpvclass+dGMclass*numclass]=ba;
    //iCSF - neighbours dGM
    G[iCSFclass+iCSFclass*numclass]=0.0;
    G[WMGMpvclass+iCSFclass*numclass]=ba;
    G[GMCSFpvclass+iCSFclass*numclass]=ba;
    //WMGMpv - neighbours WM and GM
    G[WMGMpvclass+WMGMpvclass*numclass]=0.0;
    G[GMCSFpvclass+WMGMpvclass*numclass]=ba;
    //GMCSFpv - neighbours GM and CSF
    G[GMCSFpvclass+GMCSFpvclass*numclass]=0.0;
    int i,j;

    //Copy to lower triangle
    for (j=0; j<numclass; j++) { for (i=1+j; i<numclass; i++) { G[j+i*numclass]=G[i+j*numclass]; } }
    //Copy G to H
    for (int i=0; i<(numclass*numclass); i++) { H[i]=ratio*G[i]; }
    //Print G
    if(segment_param->verbose_level>1){
        cout<<"G=" << endl;
        for (i=0; i<numclass; i++) {
            for (j=0; j<numclass; j++) {
                cout<< G[i+j*numclass] << '\t';
            }
            cout<< endl;
        }
        cout<< endl << "H=" << endl;
        for (i=0; i<numclass; i++) {
            for (j=0; j<numclass; j++) {
                cout<< G[i+j*numclass] << '\t';
            }
            cout<< endl;
        }
        cout<< endl;
    }


    return;
}

void seg_convert2binary(nifti_image *image,
                        float thresh)
{
    switch(image->datatype){
    case NIFTI_TYPE_UINT8:
        seg_convert2binary_data<unsigned char>(image,thresh);
        break;
    case NIFTI_TYPE_INT8:
        seg_convert2binary_data<char>(image,thresh);
        break;
    case NIFTI_TYPE_UINT16:
        seg_convert2binary_data<unsigned short>(image,thresh);
        break;
    case NIFTI_TYPE_INT16:
        seg_convert2binary_data<short>(image,thresh);
        break;
    case NIFTI_TYPE_UINT32:
        seg_convert2binary_data<unsigned int>(image,thresh);
        break;
    case NIFTI_TYPE_INT32:
        seg_convert2binary_data<int>(image,thresh);
        break;
    case NIFTI_TYPE_FLOAT32:
        seg_convert2binary_data<float>(image,thresh);
        break;
    case NIFTI_TYPE_FLOAT64:
        seg_convert2binary_data<PrecisionTYPE>(image,thresh);
        break;
    default:
        printf("err\tseg_changeDatatype\tThe initial image data type is not supported\n");
        return;
    }
}


template <class DTYPE>
        void seg_convert2binary_data(nifti_image *image,
                                     float thresh)
{
    // the initial array is saved and freeed
    DTYPE *initialValue = (DTYPE *)malloc(image->nvox*sizeof(DTYPE));
    memcpy(initialValue, image->data, image->nvox*sizeof(DTYPE));

    // the new array is allocated and then filled
    image->datatype = DT_BINARY;
    free(image->data);
    image->nbyper = sizeof(bool);
    image->data=(void *)calloc(image->nvox,sizeof(bool));
    bool *dataPtr = static_cast<bool *>(image->data);

    for(unsigned int i=0; i<image->nvox; i++){
        dataPtr[i] = (bool)((initialValue[i])>thresh);
    }
    free(initialValue);
    return;
}


void Normalize_NaN_Priors(nifti_image * Priors,
                          bool verbose)
{
    register int numel = Priors->nx*Priors->ny*Priors->nz;
    register int ups=0;
    register int good=0;
    if(verbose>0){
        cout<< "Normalizing Priors" << endl;
    }
    if(Priors->datatype==NIFTI_TYPE_FLOAT32){
        PrecisionTYPE * priorsptr = static_cast<PrecisionTYPE *>(Priors->data);
        for (int i=0; i<numel; i++) {
            float tempsum=0;
            for (int j=0; j<Priors->nt; j++) {
                int tempind=i+numel*j;
                if( priorsptr[tempind]<0.0 || priorsptr[tempind]!=priorsptr[tempind] || priorsptr[tempind]>1000 ){
                    priorsptr[tempind]=0.0;
                }
                tempsum+=priorsptr[tempind];
            }
            if (tempsum>0 && tempsum<1000) {
                for (int j=0; j<Priors->nt; j++) {
                    int tempind=i+numel*j;
                    priorsptr[tempind]=priorsptr[tempind]/tempsum;
                }
                good++;
            }
            else{
                for (int j=0; j<Priors->nt; j++) {
                    int tempind=i+numel*j;
                    priorsptr[tempind]=0.2f;
                }
                ups++;
            }
        }
    }
    else{

        printf("err\tNormalize_NaN_Priors\tWrong Image datatype\n");

    }

    if(verbose>0){
        cout<<"good="<<good<<"\t ups="<<ups<<"\n";
        flush(cout);
    }

    return;
}

void PriorWeight_mask(float * ShortPrior,nifti_image * Priors, float * Expec,float GaussKernelSize,float RelaxFactor, int * S2L, int * L2S,ImageSize * CurrSizes,int verbose_level){

    for(int i=0; i<(CurrSizes->numclass*CurrSizes->numelmasked);i++)ShortPrior[i]=Expec[i];
    Gaussian_Filter_Short_4D(ShortPrior,S2L,L2S,GaussKernelSize,CurrSizes,CSFclass);
    float * PriorsPtr = static_cast<float *>(Priors->data);
    int currindex=0;
    for(int k=0; k<CurrSizes->numclass; k++){
        currindex=k*CurrSizes->numelmasked;
        for(int i=0; i<(CurrSizes->numelmasked);i++){
            ShortPrior[currindex]*=(1-RelaxFactor);
            ShortPrior[currindex]+=(RelaxFactor)*PriorsPtr[S2L[i]+k*CurrSizes->numel];
            currindex++;
        }
    }

}




void Normalize_NaN_Priors_mask(nifti_image * Priors,
                               nifti_image * Mask,
                               bool verbose)
{
    register int numel = Mask->nvox;
    register int ups=0;
    register int good=0;
    if(verbose>0){
        cout<< "Normalizing Priors" << endl;
    }
    if(Mask->datatype==DT_BINARY){
        if(Priors->datatype==NIFTI_TYPE_FLOAT32){
            PrecisionTYPE * priorsptr = static_cast<PrecisionTYPE *>(Priors->data);
            bool * brainmaskptr = static_cast<bool *> (Mask->data);

            for (int i=0; i<numel; i++) {
                float tempsum=0;
                if(brainmaskptr[i]){
                    for (int j=0; j<Priors->nt; j++) {
                        int tempind=i+numel*j;
                        if( priorsptr[tempind]<0.0 || priorsptr[tempind]!=priorsptr[tempind] || priorsptr[tempind]>1000 ){
                            priorsptr[tempind]=0.0;
                        }
                        tempsum+=priorsptr[tempind];
                    }
                    if (tempsum>0 && tempsum<1000) {
                        for (int j=0; j<Priors->nt; j++) {
                            int tempind=i+numel*j;
                            priorsptr[tempind]=priorsptr[tempind]/tempsum;
                        }
                        good++;
                    }
                    else{
                        for (int j=0; j<Priors->nt; j++) {
                            int tempind=i+numel*j;
                            priorsptr[tempind]=0.2f;
                        }
                        //brainmaskptr[i]=0.0;
                        ups++;

                    }
                }
                else{
                    for (int j=0; j<Priors->nt; j++) {
                        int tempind=i+numel*j;
                        priorsptr[tempind]=0;
                    }
                }
            }
        }
        else{

            printf("err\tNormalize_NaN_Priors\tWrong Image datatype\n");

        }
    }
    else{

        printf("err\tNormalize_NaN_Priors\tWrong mask datatype\n");

    }

    if(verbose>0){
        cout<<"Priors: "<< ups<<" bad voxels" << endl;
        flush(cout);
    }

    return;
}


void Normalize_Image_mask(nifti_image * input,
                          nifti_image * Mask,
                          ImageSize * CurrSizes,
                          bool verbose)
{
    if(input->datatype!=NIFTI_TYPE_FLOAT32){
        seg_changeDatatype<float>(input);
    }
    if(verbose>0){
        cout<< "Normalizing Input Image" << endl;
    }
    int numel=(int)(rowsize(input)*colsize(input)*depth(input));
    if(Mask->datatype==DT_BINARY){
        if(input->datatype==NIFTI_TYPE_FLOAT32){
            for(int udir=0; udir<CurrSizes->usize;udir++){ // Per Multispectral Image
                for(int tdir=0; tdir<CurrSizes->tsize;tdir++){ // Per Time point Image
                    bool * brainmaskptr = static_cast<bool *> (Mask->data);
                    PrecisionTYPE * Inputptrtmp = static_cast<PrecisionTYPE *>(input->data);
                    PrecisionTYPE * Inputptr=&Inputptrtmp[numel*tdir+(CurrSizes->tsize)*numel*udir];

                    float tempmax=-1000000000;
                    float tempmin=1000000000;
                    for (int i=0; i<numel; i++) {
                        if(*brainmaskptr){
                            if (*Inputptr<tempmin) {
                                tempmin=*Inputptr;
                            }
                            if (*Inputptr>tempmax) {
                                tempmax=*Inputptr;
                            }
                        }
                        brainmaskptr++;
                        Inputptr++;
                    }
                    CurrSizes->rescale_max[udir]=tempmax;
                    CurrSizes->rescale_min[udir]=tempmin;
                    Inputptr=&Inputptrtmp[numel*tdir+(CurrSizes->tsize)*numel*udir];
                    brainmaskptr = static_cast<bool *> (Mask->data);
                    for (int i=0; i<numel; i++) {
                        //log(number_between_0_and_1 + 1)/log(2)

                        *Inputptr=logf((((*Inputptr)-tempmin)/(tempmax-tempmin))+1)/0.693147181;

                        brainmaskptr++;
                        Inputptr++;
                    }
                }
            }
        }
        else{
            {
                printf("err\tNormalize_T1\tWrong Mask datatype\n");
            }
        }
    }
    return;
}

void Normalize_Image(nifti_image * input,
                     ImageSize * CurrSizes,
                     bool verbose)
{
    if(input->datatype!=NIFTI_TYPE_FLOAT32){
        seg_changeDatatype<float>(input);
    }
    if(verbose>0){
        cout<< "Normalizing Input Image" << endl;
    }
    if(input->datatype==NIFTI_TYPE_FLOAT32){
        // if mask is not set up
        for(int udir=0; udir<CurrSizes->usize;udir++){ // Per Multispectral Image
            for(int tdir=0; tdir<CurrSizes->tsize;tdir++){ // Per Time point Image
                int numel=(int)(rowsize(input)*colsize(input)*depth(input));
                PrecisionTYPE * Inputptrtmp = static_cast<PrecisionTYPE *>(input->data);
                PrecisionTYPE * Inputptr=&Inputptrtmp[numel*tdir+(CurrSizes->tsize)*numel*udir];

                float tempmax=0;
                float tempmin=10000000000;
                for (int i=0; i<numel; i++) {
                    if (*Inputptr<tempmin) {
                        tempmin=*Inputptr;
                    }
                    if (*Inputptr>tempmax) {
                        tempmax=*Inputptr;
                    }
                    Inputptr++;
                }
                CurrSizes->rescale_max[udir]=tempmax;
                CurrSizes->rescale_min[udir]=tempmin;
                Inputptr=&Inputptrtmp[numel*tdir+(CurrSizes->tsize)*numel*udir];
                for (int i=0; i<numel; i++) {
                    //log(number_between_0_and_1 + 1)/log(2)
                    *Inputptr=logf((((*Inputptr)-tempmin)/(tempmax-tempmin))+1)/0.693147181;
                    Inputptr++;
                }
            }
        }
    }
    return;
}

void Normalize_T1_and_MV(nifti_image * T1,
                         nifti_image * Mask,
                         PrecisionTYPE * M,
                         PrecisionTYPE * V,
                         ImageSize * CurrSizes)
{

    PrecisionTYPE * T1ptrtmp = static_cast<PrecisionTYPE *>(T1->data);
    bool * brainmaskptr = static_cast<bool *>(Mask->data);
    PrecisionTYPE * T1ptr=T1ptrtmp;
    int numel=(int)(rowsize(T1)*colsize(T1)*depth(T1));
    float tempmax=0;
    float tempmin=100000;
    for (int i=0; i<numel; i++) {
        if(*brainmaskptr){
            if (*T1ptr<tempmin) {
                tempmin=*T1ptr;
            }
            if (*T1ptr>tempmax) {
                tempmax=*T1ptr;
            }
        }
        brainmaskptr++;
        T1ptr++;
    }
    CurrSizes->rescale_max[0]=tempmax;
    CurrSizes->rescale_min[0]=tempmin;
    T1ptr=T1ptrtmp;
    brainmaskptr = static_cast<bool *>(Mask->data);
    for (int i=0; i<numel; i++) {
        if(*brainmaskptr){
            *T1ptr=logf((((*T1ptr)-tempmin)/(tempmax-tempmin))+1)/0.693147181;
        }else{
            *T1ptr=0;
        }
        brainmaskptr++;
        T1ptr++;
    }
    for(int cl=0; cl<CurrSizes->numclass; cl++){
        M[cl]=log(1+(M[cl]-tempmin)/(tempmax-tempmin));
        V[cl]=log(1+(V[cl])/(tempmax-tempmin));
    }

    return;
}

int * Create_Short_2_Long_Matrix_from_NII(nifti_image * Mask,
                                          int * shortsize){
    int numel_masked=0;
    int numel = Mask->nvox;
    if(Mask->datatype==DT_BINARY){
        bool * Maskptr = static_cast<bool *> (Mask->data);
        bool * Maskptrtmp = Maskptr;
        for (int i=0; i<numel; i++, Maskptrtmp++) {
            (*Maskptrtmp)>0?numel_masked++:0;
        }
        shortsize[0]=numel_masked;

        int * Short_2_Long_Indices= new int [numel_masked]();
        int * Short_2_Long_Indices_PTR = (int *)(Short_2_Long_Indices);

        Maskptrtmp = Maskptr;
        int tempindex=0;
        for (int i=0; i<numel; i++) {
            if ((*Maskptrtmp)>0) {
                Short_2_Long_Indices_PTR[tempindex]=i;
                tempindex++;
            }
            Maskptrtmp++;
        }
        return Short_2_Long_Indices;
    }
    else {
        printf("err\tCreate_Short_2_Long_Matrix\tWrong Mask datatype\n");
        return NULL;

    }
}

int *  Create_Long_2_Short_Matrix_from_NII(nifti_image * Mask){
    int numel = Mask->nvox;
    if(Mask->datatype==DT_BINARY){


        bool * Maskptr = static_cast<bool *> (Mask->data);
        bool * Maskptrtmp = Maskptr;
        int * Long_2_Short_Indices= new int [rowsize(Mask)*colsize(Mask)*depth(Mask)]();
        int * Long_2_Short_Indices_PTR = (int *) Long_2_Short_Indices;

        Maskptrtmp = Maskptr;
        int tempindex=0;
        for (int i=0; i<numel; i++,Maskptrtmp++,Long_2_Short_Indices_PTR++) {
            if ((*Maskptrtmp)>0) {
                (*Long_2_Short_Indices_PTR)=tempindex;
                tempindex++;
            }
            else{
                (*Long_2_Short_Indices_PTR)=-1;
            }
        }
        return Long_2_Short_Indices;
    }
    else {
        cout<< "err\tCreate_Correspondace_Matrices\tWrong Mask datatype\n" << endl;
    }
}


int * Create_Short_2_Long_Matrix_from_Carray(bool * Mask,
                                             int * shortsize,
                                             int nvox){
    int numel_masked=0;
    int numel = nvox;
    bool * Maskptrtmp =(bool *) Mask;
    for (int i=0; i<numel; i++) {
        Mask[i]?numel_masked++:0;
    }
    shortsize[0]=numel_masked;

    int * Short_2_Long_Indices= new int [numel_masked]();
    int * Short_2_Long_Indices_PTR = (int *)(Short_2_Long_Indices);

    int tempindex=0;
    for (int i=0; i<numel; i++) {
        if ((Mask[i])>0) {
            Short_2_Long_Indices_PTR[tempindex]=i;
            tempindex++;
        }

    }
    return Short_2_Long_Indices;

}

int *  Create_Long_2_Short_Matrix_from_Carray(bool * Mask,
                                              int nvox){
    int * Long_2_Short_Indices= new int [nvox]();
    int * Long_2_Short_Indices_PTR = (int *) Long_2_Short_Indices;
    int tempindex=0;
    for (int i=0; i<nvox; i++,Long_2_Short_Indices_PTR++) {
        if ((Mask[i])>0) {
            (*Long_2_Short_Indices_PTR)=tempindex;
            tempindex++;
        }
        else{
            (*Long_2_Short_Indices_PTR)=-1;
        }
    }
    return Long_2_Short_Indices;

}

PrecisionTYPE * Create_cArray_from_Prior_mask(nifti_image * Mask,
                                              nifti_image * Priors,
                                              int numclass,
                                              bool PV_ON)
{
    register int numel=(int)(rowsize(Mask)*colsize(Mask)*depth(Mask));
    register int numel_masked=0;

    bool * Maskptrtmp = static_cast<bool *> (Mask->data);;
    for (int i=0; i<numel; i++, Maskptrtmp++) {
        *Maskptrtmp?numel_masked++:0;
    }
    int pluspv=(int)(PV_ON)*2;

    PrecisionTYPE * Expec = new PrecisionTYPE [numel_masked*(numclass+pluspv)] ();
    PrecisionTYPE * tempExpec= (PrecisionTYPE *) Expec;
    PrecisionTYPE * PriorPTR = static_cast<PrecisionTYPE *>(Priors->data);
    for(int cl=0; cl<numclass;cl++){
        Maskptrtmp = static_cast<bool *> (Mask->data);;
        for (int i=numel; i--; Maskptrtmp++,PriorPTR++) {
            if(*Maskptrtmp){
                *tempExpec = *PriorPTR;
                tempExpec++;
            }
        }
    }

    return Expec;
}

PrecisionTYPE * Create_cArray_from_Prior(nifti_image * Priors,
                                         int numclass,
                                         bool PV_ON)
{
    register int numel=(int)(rowsize(Priors)*colsize(Priors)*depth(Priors));
    int pluspv=(int)(PV_ON)*2;
    PrecisionTYPE * Expec = new PrecisionTYPE [numel*(numclass+pluspv)] ();
    PrecisionTYPE * Expec_PTR= Expec;
    PrecisionTYPE * PriorPTR = static_cast<PrecisionTYPE *>(Priors->data);
    for(int cl=0; cl<numclass;cl++){
        for (int i=numel; i--; PriorPTR++,Expec_PTR++) {
            *Expec_PTR = *PriorPTR;
        }
    }
    return Expec;
}

PrecisionTYPE * Create_cArray_from_3D_image(nifti_image * Mask,
                                            nifti_image * SourceImage)
{
    register int numel=(int)(rowsize(Mask)*colsize(Mask)*depth(Mask));
    register int numel_masked=0;

    bool * Maskptrtmp = static_cast<bool *> (Mask->data);;
    for (int i=0; i<numel; i++, Maskptrtmp++) {
        *Maskptrtmp?numel_masked++:0;
    }

    PrecisionTYPE * outimage = new PrecisionTYPE [numel_masked] ();
    PrecisionTYPE * outimage_ptr= outimage;
    PrecisionTYPE * SourceImagePTR = static_cast<PrecisionTYPE *>(SourceImage->data);
    Maskptrtmp = static_cast<bool *> (Mask->data);

    for (int i=numel; i--; Maskptrtmp++,SourceImagePTR++) {
        if(*Maskptrtmp){
            *outimage_ptr = *SourceImagePTR;
            outimage_ptr++;
        }
    }


    return outimage;
}


void calcE_mask(nifti_image * T1,
                PrecisionTYPE * IterPrior,
                PrecisionTYPE * Expec,
                PrecisionTYPE * loglik,
                PrecisionTYPE * BiasField,
                int * S2L,
                PrecisionTYPE * M,
                PrecisionTYPE * V,
                ImageSize * CurrSizes,
                int verbose)
{
    int numel_masked=CurrSizes->numelmasked;
    int num_class=CurrSizes->numclass;

    PrecisionTYPE * IterPrior_PTR= (PrecisionTYPE *) IterPrior;
    PrecisionTYPE * Expec_PTR= (PrecisionTYPE *) Expec;
    PrecisionTYPE * T1_PTR = static_cast<PrecisionTYPE *>(T1->data);
    PrecisionTYPE SumExpec=0.0f;
    PrecisionTYPE expectmp=0.0f;



    PrecisionTYPE inv_v [max_numbclass*MaxMultispectalSize*MaxMultispectalSize]={0.0f};
    PrecisionTYPE inv_sqrt_V_2pi [max_numbclass]={0.0f};
    PrecisionTYPE newM [max_numbclass*MaxMultispectalSize]={0.0f};
    PrecisionTYPE T1_Bias_corr[MaxMultispectalSize]={0.0f};
    PrecisionTYPE tmpT1_BC_minusM=0;

    int Expec_offset [max_numbclass]={0};

    for (int cl=0; cl<num_class; cl++) {
        Expec_offset[cl]=(int) cl*numel_masked;
        if(CurrSizes->usize>1){
            matrix <double> Vmat(CurrSizes->usize,CurrSizes->usize);

            for(int j2=0; j2<CurrSizes->usize; j2++){
                for(int i2=j2; i2<CurrSizes->usize; i2++){
                    Vmat.setvalue(i2,j2,(double)(V[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]));
                    Vmat.setvalue(j2,i2,(double)(V[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]));
                }
            }
            inv_sqrt_V_2pi[cl]=1/(sqrtf(2*M_PI*Vmat.determinant()));
            if(verbose>1){
                cout<<endl<<"inv_sqrt_V_2pi["<< cl <<"]= "<< inv_sqrt_V_2pi[cl] << endl;
                flush(cout);
            }
            Vmat.invert();
            double cvalue;
            bool success;
            if(verbose>1){
                cout<<"inv_V["<< cl <<"]= ";
                flush(cout);
            }
            for(int j2=0; j2<CurrSizes->usize; j2++){
                if(verbose>1){
                    if(j2!=0){
                        cout<< endl << "          ";
                    }
                }
                for(int i2=0; i2<CurrSizes->usize; i2++){
                    Vmat.getvalue(i2,j2,cvalue,success);
                    inv_v[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]=(PrecisionTYPE)(cvalue);
                    if(verbose>1){
                        cout<<inv_v[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]<< "\t";
                        flush(cout);
                    }
                }

            }
            if(verbose>1){
                cout<< endl;
            }
        }
        else{
            inv_sqrt_V_2pi[cl]=1/(sqrtf(2*M_PI* V[cl]));
            inv_v[cl]=1/V[cl];
            newM[cl]=M[cl];

        }
    }
    loglik[0]=0;
    SumExpec=0.0f;
    //int * Expec_offset_PTR= (int *) Expec_offset;
    register PrecisionTYPE tempvar=0.0f;
    float mahal=0.0f;
    float logliktmp=0.0f;
    for (int i=0; i<numel_masked;i++, Expec_PTR++, IterPrior_PTR++) {
        for(int Multispec=0; Multispec<CurrSizes->usize; Multispec++) {
            T1_Bias_corr[Multispec]=(BiasField!=NULL)?(T1_PTR[S2L[i]+Multispec*CurrSizes->numel] + BiasField[i+Multispec*numel_masked]):(T1_PTR[S2L[i]+Multispec*CurrSizes->numel]);
        }
        SumExpec=0.0f;
        expectmp=0.0f;
        mahal=0.0f;
        //Expec_offset_PTR=Expec_offset;
        for (int cl=0; cl<num_class; cl++) {
            mahal=0.0f;
            for(int Multispec=0; Multispec<CurrSizes->usize; Multispec++) {
                tmpT1_BC_minusM=(T1_Bias_corr[Multispec] - M[cl*(CurrSizes->usize)+Multispec]);
                for(int Multispec2=0; Multispec2<CurrSizes->usize; Multispec2++) {
                    mahal-=(0.5f)*(T1_Bias_corr[Multispec2] - M[cl*(CurrSizes->usize)+Multispec2])*inv_v[cl*CurrSizes->usize*CurrSizes->usize+Multispec+Multispec2*CurrSizes->usize]*tmpT1_BC_minusM;
                }
            }
            Expec_PTR[Expec_offset[cl]]=IterPrior_PTR[Expec_offset[cl]] * expf(mahal) * inv_sqrt_V_2pi[cl];

            SumExpec+=Expec_PTR[Expec_offset[cl]];
        }

        logliktmp+=logf(SumExpec);



        if (SumExpec<=0.0 || SumExpec!=SumExpec){
            for (int cl=0; cl<num_class; cl++) {
                Expec_PTR[Expec_offset[cl]]=(float)(1)/(float)(num_class);
            }
        }
        else{
            logliktmp+=logf(SumExpec);
            for (int cl=0; cl<num_class; cl++) {
                Expec_PTR[Expec_offset[cl]]=Expec_PTR[Expec_offset[cl]]/SumExpec;
            }
        }
    }
    loglik[0]=logliktmp;
    return;
}


void calcE(nifti_image * T1,
           PrecisionTYPE * IterPrior,
           PrecisionTYPE * Expec,
           PrecisionTYPE * loglik,
           PrecisionTYPE * BiasField,
           PrecisionTYPE * M,
           PrecisionTYPE * V,
           ImageSize * CurrSizes,
           int verbose)
{
    int numel=CurrSizes->numel;
    int num_class=CurrSizes->numclass;

    PrecisionTYPE * IterPrior_PTR= (PrecisionTYPE *) IterPrior;
    PrecisionTYPE * Expec_PTR= (PrecisionTYPE *) Expec;
    PrecisionTYPE * T1_PTR = static_cast<PrecisionTYPE *>(T1->data);
    PrecisionTYPE SumExpec=0.0f;
    PrecisionTYPE expectmp=0.0f;



    PrecisionTYPE inv_v [max_numbclass*MaxMultispectalSize*MaxMultispectalSize]={0.0f};
    PrecisionTYPE inv_sqrt_V_2pi [max_numbclass]={0.0f};
    PrecisionTYPE newM [max_numbclass*MaxMultispectalSize]={0.0f};
    PrecisionTYPE T1_Bias_corr[MaxMultispectalSize]={0.0f};
    PrecisionTYPE tmpT1_BC_minusM=0;

    int Expec_offset [max_numbclass]={0};

    for (int cl=0; cl<num_class; cl++) {
        Expec_offset[cl]=(int) cl*numel;
        if(CurrSizes->usize>1){
            matrix <double> Vmat(CurrSizes->usize,CurrSizes->usize);

            for(int j2=0; j2<CurrSizes->usize; j2++){
                for(int i2=j2; i2<CurrSizes->usize; i2++){
                    Vmat.setvalue(i2,j2,(double)(V[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]));
                    Vmat.setvalue(j2,i2,(double)(V[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]));
                }
            }
            inv_sqrt_V_2pi[cl]=1/(sqrtf(2*M_PI)*Vmat.determinant());
            if(verbose>1){
                cout<<endl<<"inv_sqrt_V_2pi["<< cl <<"]= "<< inv_sqrt_V_2pi[cl] << endl;
                flush(cout);
            }
            Vmat.invert();
            double cvalue;
            bool success;
            if(verbose>1){
                cout<<"inv_V["<< cl <<"]= ";
                flush(cout);
            }
            for(int j2=0; j2<CurrSizes->usize; j2++){
                if(verbose>1){
                    if(j2!=0){
                        cout<< endl << "          ";
                    }
                }
                for(int i2=0; i2<CurrSizes->usize; i2++){
                    Vmat.getvalue(i2,j2,cvalue,success);
                    inv_v[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]=(PrecisionTYPE)(cvalue);
                    if(verbose>1){
                        cout<<inv_v[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]<< "\t";
                        flush(cout);
                    }
                }

            }
            if(verbose>1){
                cout<< endl;
            }
        }
        else{
            inv_sqrt_V_2pi[cl]=1/(sqrtf(2*M_PI)* V[cl]);
            inv_v[cl]=1/V[cl];
            newM[cl]=M[cl];

        }
    }
    loglik[0]=0;
    SumExpec=0.0f;
    //int * Expec_offset_PTR= (int *) Expec_offset;
    register PrecisionTYPE tempvar=0.0f;
    float mahal=0.0f;
    float logliktmp=0.0f;
    for (int i=0; i<numel;i++, Expec_PTR++, IterPrior_PTR++) {
        for(int Multispec=0; Multispec<CurrSizes->usize; Multispec++) {
            T1_Bias_corr[Multispec]=(BiasField!=NULL)?(T1_PTR[i+Multispec*numel] + BiasField[i+Multispec*numel]):(T1_PTR[i+Multispec*numel]);
        }
        SumExpec=0.0f;
        expectmp=0.0f;
        mahal=0.0f;
        //Expec_offset_PTR=Expec_offset;
        for (int cl=0; cl<num_class; cl++) {
            mahal=0.0f;
            for(int Multispec=0; Multispec<CurrSizes->usize; Multispec++) {
                tmpT1_BC_minusM=(T1_Bias_corr[Multispec] - M[cl*(CurrSizes->usize)+Multispec]);
                for(int Multispec2=0; Multispec2<CurrSizes->usize; Multispec2++) {
                    mahal-=(0.5f)*(T1_Bias_corr[Multispec2] - M[cl*(CurrSizes->usize)+Multispec2])*inv_v[cl*CurrSizes->usize*CurrSizes->usize+Multispec+Multispec2*CurrSizes->usize]*tmpT1_BC_minusM;
                }
            }
            Expec_PTR[Expec_offset[cl]]=IterPrior_PTR[Expec_offset[cl]] * expf(mahal) * inv_sqrt_V_2pi[cl];
            SumExpec+=Expec_PTR[Expec_offset[cl]];
        }
        if (SumExpec<=0.0 || SumExpec!=SumExpec){
            SumExpec=0.01;
        }
        logliktmp+=logf(SumExpec);

        if (SumExpec<=0.0 || SumExpec!=SumExpec){
            for (int cl=0; cl<num_class; cl++) {
                Expec_PTR[Expec_offset[cl]]=(float)(1)/(float)(num_class);
            }
        }
        else{
            logliktmp+=logf(SumExpec);
            for (int cl=0; cl<num_class; cl++) {
                Expec_PTR[Expec_offset[cl]]=Expec_PTR[Expec_offset[cl]]/SumExpec;
            }
        }
    }
    loglik[0]=logliktmp;
    return;
}




/*
void calcE_mask_aprox(nifti_image * T1,
                      PrecisionTYPE * IterPrior,
                      PrecisionTYPE * Expec,
                      PrecisionTYPE * loglik,
                      PrecisionTYPE * BiasField,
                      int * S2L,
                      PrecisionTYPE * M,
                      PrecisionTYPE * V,
                      ImageSize * CurrSizes,
                      int verbose)
{
    int numel_masked=CurrSizes->numelmasked;
    int num_class=CurrSizes->numclass;

    PrecisionTYPE * IterPrior_PTR= (PrecisionTYPE *) IterPrior;
    PrecisionTYPE * Expec_PTR= (PrecisionTYPE *) Expec;
    PrecisionTYPE * T1_PTR = static_cast<PrecisionTYPE *>(T1->data);
    PrecisionTYPE SumExpec=0.0f;
    PrecisionTYPE expectmp=0.0f;



    PrecisionTYPE inv_v [max_numbclass*MaxMultispectalSize*MaxMultispectalSize]={0.0f};
    PrecisionTYPE inv_sqrt_V_2pi [max_numbclass]={0.0f};
    PrecisionTYPE newM [max_numbclass*MaxMultispectalSize]={0.0f};
    PrecisionTYPE T1_Bias_corr[MaxMultispectalSize]={0.0f};
    PrecisionTYPE tmpT1_BC_minusM=0;

    int Expec_offset [max_numbclass]={0};

    for (int cl=0; cl<num_class; cl++) {
        Expec_offset[cl]=(int) cl*numel_masked;
        if(CurrSizes->usize>1){
            matrix <double> Vmat(CurrSizes->usize,CurrSizes->usize);

            for(int j2=0; j2<CurrSizes->usize; j2++){
                for(int i2=j2; i2<CurrSizes->usize; i2++){
                    Vmat.setvalue(i2,j2,(double)(V[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]));
                    Vmat.setvalue(j2,i2,(double)(V[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]));
                }
            }
            inv_sqrt_V_2pi[cl]=1/(sqrtf(2*M_PI)*Vmat.determinant());
            if(verbose>1){
                cout<<endl<<"inv_sqrt_V_2pi["<< cl <<"]= "<< inv_sqrt_V_2pi[cl] << endl;
                flush(cout);
            }
            Vmat.invert();
            double cvalue;
            bool success;
            if(verbose>1){
                cout<<"inv_V["<< cl <<"]= ";
                flush(cout);
            }
            for(int j2=0; j2<CurrSizes->usize; j2++){
                if(verbose>1){
                    if(j2!=0){
                        cout<< endl << "          ";
                    }
                }
                for(int i2=0; i2<CurrSizes->usize; i2++){
                    Vmat.getvalue(i2,j2,cvalue,success);


                    inv_v[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]=(PrecisionTYPE)(cvalue);
                    if(verbose>1){
                        cout<<inv_v[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]<< "\t";
                        flush(cout);
                    }
                }

            }
            if(verbose>1){
                cout<< endl;
            }
        }
        else{
            inv_sqrt_V_2pi[cl]=1/(sqrtf(2*M_PI)* V[cl]);
            inv_v[cl]=1/V[cl];
            newM[cl]=M[cl];

        }
    }
    loglik[0]=0;
    SumExpec=0.0f;
    //int * Expec_offset_PTR= (int *) Expec_offset;
    register PrecisionTYPE tempvar=0.0f;
    float mahal=0.0f;
    float logliktmp=0.0f;
    for (int i=0; i<numel_masked;i++, Expec_PTR++, IterPrior_PTR++) {
        for(int Multispec=0; Multispec<CurrSizes->usize; Multispec++) {
            T1_Bias_corr[Multispec]=(BiasField!=NULL)?(T1_PTR[S2L[i]+Multispec*CurrSizes->numel] + BiasField[i+Multispec*numel_masked]):(T1_PTR[S2L[i]+Multispec*CurrSizes->numel]);
        }
        SumExpec=0.0f;
        expectmp=0.0f;
        mahal=0.0f;
        //Expec_offset_PTR=Expec_offset;
        for (int cl=0; cl<num_class; cl++) {
            mahal=0.0f;
            for(int Multispec=0; Multispec<CurrSizes->usize; Multispec++) {
                tmpT1_BC_minusM=(T1_Bias_corr[Multispec] - M[cl*(CurrSizes->usize)+Multispec]);
                for(int Multispec2=0; Multispec2<CurrSizes->usize; Multispec2++) {
                    mahal-=(0.5f)*(T1_Bias_corr[Multispec2] - M[cl*(CurrSizes->usize)+Multispec2])*inv_v[cl*CurrSizes->usize*CurrSizes->usize+Multispec2+Multispec*CurrSizes->usize]*tmpT1_BC_minusM;
                }
            }
            Expec_PTR[Expec_offset[cl]]=IterPrior_PTR[Expec_offset[cl]] * (float)(EXP(mahal)) * inv_sqrt_V_2pi[cl];
            SumExpec+=Expec_PTR[Expec_offset[cl]];
        }
        if (SumExpec<=0.0 || SumExpec!=SumExpec){
            SumExpec=0.01;
        }
        logliktmp+=logf(SumExpec);



        if (SumExpec<=0.0 || SumExpec!=SumExpec){
            for (int cl=0; cl<num_class; cl++) {
                Expec_PTR[Expec_offset[cl]]=(float)(1)/(float)(num_class);
            }
        }
        else{
            logliktmp+=logf(SumExpec);
            for (int cl=0; cl<num_class; cl++) {

                Expec_PTR[Expec_offset[cl]]=Expec_PTR[Expec_offset[cl]]/SumExpec;
            }
        }
    }
    loglik[0]=logliktmp;
    return;
}


*/
/*
void calcE_aprox(nifti_image * T1,
                 PrecisionTYPE * IterPrior,
                 PrecisionTYPE * Expec,
                 PrecisionTYPE * loglik,
                 PrecisionTYPE * BiasField,
                 PrecisionTYPE * M,
                 PrecisionTYPE * V,
                 ImageSize * CurrSizes,
                 int verbose)
{
    int numel=CurrSizes->numel;
    int num_class=CurrSizes->numclass;

    PrecisionTYPE * IterPrior_PTR= (PrecisionTYPE *) IterPrior;
    PrecisionTYPE * Expec_PTR= (PrecisionTYPE *) Expec;
    PrecisionTYPE * T1_PTR = static_cast<PrecisionTYPE *>(T1->data);
    PrecisionTYPE SumExpec=0.0f;
    PrecisionTYPE expectmp=0.0f;



    PrecisionTYPE inv_v [max_numbclass*MaxMultispectalSize*MaxMultispectalSize]={0.0f};
    PrecisionTYPE inv_sqrt_V_2pi [max_numbclass]={0.0f};
    PrecisionTYPE newM [max_numbclass*MaxMultispectalSize]={0.0f};
    PrecisionTYPE T1_Bias_corr[MaxMultispectalSize]={0.0f};
    PrecisionTYPE tmpT1_BC_minusM=0;

    int Expec_offset [max_numbclass]={0};

    for (int cl=0; cl<num_class; cl++) {
        Expec_offset[cl]=(int) cl*numel;
        if(CurrSizes->usize>1){
            matrix <double> Vmat(CurrSizes->usize,CurrSizes->usize);

            for(int j2=0; j2<CurrSizes->usize; j2++){
                for(int i2=j2; i2<CurrSizes->usize; i2++){
                    Vmat.setvalue(i2,j2,(double)(V[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]));
                    Vmat.setvalue(j2,i2,(double)(V[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]));
                }
            }
            inv_sqrt_V_2pi[cl]=1/(sqrtf(2*M_PI)*Vmat.determinant());
            if(verbose>1){
                cout<<endl<<"inv_sqrt_V_2pi["<< cl <<"]= "<< inv_sqrt_V_2pi[cl] << endl;
                flush(cout);
            }
            Vmat.invert();
            double cvalue;
            bool success;
            if(verbose>1){
                cout<<"inv_V["<< cl <<"]= ";
                flush(cout);
            }
            for(int j2=0; j2<CurrSizes->usize; j2++){
                if(verbose>1){
                    if(j2!=0){
                        cout<< endl << "          ";
                    }
                }
                for(int i2=0; i2<CurrSizes->usize; i2++){
                    Vmat.getvalue(i2,j2,cvalue,success);
                    inv_v[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]=(PrecisionTYPE)(cvalue);
                    if(verbose>1){
                        cout<<inv_v[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]<< "\t";
                        flush(cout);
                    }
                }

            }
            if(verbose>1){
                cout<< endl;
            }
        }
        else{
            inv_sqrt_V_2pi[cl]=1/(sqrtf(2*M_PI)* V[cl]);
            inv_v[cl]=1/V[cl];
            newM[cl]=M[cl];

        }
    }
    loglik[0]=0;
    SumExpec=0.0f;
    //int * Expec_offset_PTR= (int *) Expec_offset;
    register PrecisionTYPE tempvar=0.0f;
    float mahal=0.0f;
    float logliktmp=0.0f;
    IterPrior_PTR= (PrecisionTYPE *) IterPrior;
    Expec_PTR= (PrecisionTYPE *) Expec;
    for (int i=0; i<numel;i++, Expec_PTR++, IterPrior_PTR++) {
        for(int Multispec=0; Multispec<CurrSizes->usize; Multispec++) {
            T1_Bias_corr[Multispec]=(BiasField!=NULL)?(T1_PTR[i+Multispec*numel] + BiasField[i+Multispec*numel]):(T1_PTR[i+Multispec*numel]);
        }
        SumExpec=0.0f;
        expectmp=0.0f;
        mahal=0.0f;
        //Expec_offset_PTR=Expec_offset;
        for (int cl=0; cl<num_class; cl++) {
            mahal=0.0f;
            for(int Multispec=0; Multispec<CurrSizes->usize; Multispec++) {
                tmpT1_BC_minusM=(T1_Bias_corr[Multispec] - M[cl*(CurrSizes->usize)+Multispec]);
                for(int Multispec2=0; Multispec2<CurrSizes->usize; Multispec2++) {
                    mahal-=(0.5f)*(T1_Bias_corr[Multispec2] - M[cl*(CurrSizes->usize)+Multispec2])*inv_v[cl*CurrSizes->usize*CurrSizes->usize+Multispec+Multispec2*CurrSizes->usize]*tmpT1_BC_minusM;
                }
            }

            Expec_PTR[Expec_offset[cl]]=IterPrior_PTR[Expec_offset[cl]] * (float)(EXP(mahal)) * inv_sqrt_V_2pi[cl];
            SumExpec+=Expec_PTR[Expec_offset[cl]];
        }
        if (SumExpec<=0.0 || SumExpec!=SumExpec){
            SumExpec=0.01;
        }
        logliktmp+=logf(SumExpec);

        if (SumExpec<=0.0 || SumExpec!=SumExpec){
            for (int cl=0; cl<num_class; cl++) {
                Expec_PTR[Expec_offset[cl]]=(float)(1)/(float)(num_class);
            }
        }
        else{
            logliktmp+=logf(SumExpec);
            for (int cl=0; cl<num_class; cl++) {
                Expec_PTR[Expec_offset[cl]]=Expec_PTR[Expec_offset[cl]]/SumExpec;
            }
        }
    }
    loglik[0]=logliktmp;
    return;
}



*/
void calcM(nifti_image * T1,
           PrecisionTYPE * Expec,
           PrecisionTYPE * BiasField,
           PrecisionTYPE * M,
           PrecisionTYPE * V,
           PrecisionTYPE * M_MAP,
           PrecisionTYPE * V_MAP,
           ImageSize * CurrSizes,
           int verbose)
{

    if(verbose>0){
        cout<< "Optimising Gaussian Parameters" << endl;
        flush(cout);
    }
    int numel=CurrSizes->numel;
    int currentnum_class=CurrSizes->numclass;
    PrecisionTYPE * Expec_PTR = (PrecisionTYPE *) Expec;
    PrecisionTYPE * T1_PTR = static_cast<PrecisionTYPE *>(T1->data);
    PrecisionTYPE * BiasField_PTR= (PrecisionTYPE *) BiasField;
    PrecisionTYPE * T1_PTR2= static_cast<PrecisionTYPE *>(T1->data);
    PrecisionTYPE * BiasField_PTR2= (PrecisionTYPE *) BiasField;

    int Expec_offset[PV_numbclass];
    for (int cl=0; cl<currentnum_class; cl++) {
        Expec_offset[cl]=cl*numel;
    }
    PrecisionTYPE tempsum= (PrecisionTYPE) 0.0;
    PrecisionTYPE SumPriors= (PrecisionTYPE) 0.0;
    PrecisionTYPE tempexpec= (PrecisionTYPE) 0.0;

    // ***********

    for (int cl=0; cl<currentnum_class; cl++) {
        // MEAN
        tempexpec=0;
        for(int Multispec=0; Multispec<CurrSizes->usize; Multispec++) {
            Expec_PTR=(PrecisionTYPE *) &Expec[Expec_offset[cl]];

            T1_PTR = static_cast<PrecisionTYPE *>(T1->data);
            T1_PTR =&T1_PTR[Multispec*CurrSizes->numel];
            tempsum=(PrecisionTYPE)0.0;
            SumPriors=(PrecisionTYPE)0.0;

            if(BiasField!=NULL){
                BiasField_PTR= &BiasField[Multispec*numel];
                for (int i=0; i<numel; i++, Expec_PTR++,BiasField_PTR++) {
                    tempsum+=(*Expec_PTR) * (T1_PTR[i]+(*BiasField_PTR));
                    SumPriors+=(*Expec_PTR);
                }
            }else{
                for (int i=0; i<numel; i++, Expec_PTR++) {
                    tempsum+=(*Expec_PTR)*(T1_PTR[i]);
                    SumPriors+=(*Expec_PTR);
                }
            }
            if(M_MAP==NULL){
                M[cl*(CurrSizes->usize)+Multispec]=tempsum/SumPriors;
            }
            else{
                M[cl*(CurrSizes->usize)+Multispec]=(tempsum/SumPriors/powf(V[cl*(CurrSizes->usize)+Multispec],2)+M_MAP[cl*(CurrSizes->usize)+Multispec]/powf(V_MAP[cl*(CurrSizes->usize)+Multispec],2))/(1/powf(V[cl*(CurrSizes->usize)+Multispec],2)+1/powf(V_MAP[cl*(CurrSizes->usize)+Multispec],2));
            }

            for(int Multispec2=Multispec; Multispec2<CurrSizes->usize; Multispec2++) {

                T1_PTR = static_cast<PrecisionTYPE *>(T1->data);
                T1_PTR =&T1_PTR[Multispec*CurrSizes->numel];

                T1_PTR2 = static_cast<PrecisionTYPE *>(T1->data);
                T1_PTR2 =&T1_PTR2[Multispec2*CurrSizes->numel];
                float tmpM=M[cl*CurrSizes->usize+Multispec];
                float tmpM2=M[cl*CurrSizes->usize+Multispec2];
                //STD
                tempsum=0;
                Expec_PTR=&Expec[Expec_offset[cl]];
                if(BiasField!=NULL){
                    BiasField_PTR=&BiasField[Multispec*numel];
                    BiasField_PTR2=&BiasField[Multispec2*numel];
                    for (int i=0; i<numel;i++,Expec_PTR++,BiasField_PTR++,BiasField_PTR2++) {
                        tempsum+=(*Expec_PTR) * (T1_PTR[i]+(*BiasField_PTR)-tmpM) * (T1_PTR2[i]+(*BiasField_PTR2)-tmpM2);
                    }
                }else{
                    for (int i=0; i<numel;i++,Expec_PTR++) {
                        tempsum+=(*Expec_PTR) * (T1_PTR[i]-tmpM) * (T1_PTR2[i]-tmpM2);
                    }

                }
                V[cl*CurrSizes->usize*CurrSizes->usize+Multispec+Multispec2*CurrSizes->usize]=tempsum/SumPriors+0.00001f;
                V[cl*CurrSizes->usize*CurrSizes->usize+Multispec2+Multispec*CurrSizes->usize]=V[cl*CurrSizes->usize*CurrSizes->usize+Multispec+Multispec2*CurrSizes->usize];
            }
        }
        if(verbose>0){
            if(CurrSizes->usize==1){
                cout.fill('0');
                cout<< "M["<<(int)(cl)<<"]= "<<setw(10)<<setprecision(7)<<left<<(PrecisionTYPE)(M[cl])<<"\tV["<<(int)(cl)<<"]="<<setw(10)<<setprecision(7)<<left<<(PrecisionTYPE)(V[cl])<< endl;
                flush(cout);
            }
            else{

                cout<< "M["<<(int)(cl)<<"]= ";
                for(int Multispec=0; Multispec<CurrSizes->usize; Multispec++) {
                    cout<< setw(10)<<setprecision(7)<<left<<(PrecisionTYPE)(M[cl*CurrSizes->usize+Multispec])<<"\t";
                }
                cout<< endl;
                flush(cout);
                cout<< "V["<<(int)(cl)<<"]= ";
                for(int Multispec=0; Multispec<CurrSizes->usize; Multispec++) {
                    if(Multispec>0){
                        cout<< "      ";
                    }
                    for(int Multispec2=0; Multispec2<CurrSizes->usize; Multispec2++) {
                        cout<<setw(10)<<setprecision(7)<<left<<(PrecisionTYPE)(V[cl*CurrSizes->usize*CurrSizes->usize+Multispec*CurrSizes->usize+Multispec2])<<"\t";
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




void calcM_mask(nifti_image * T1,
                PrecisionTYPE * Expec,
                PrecisionTYPE * BiasField,
                int * S2L,
                PrecisionTYPE * M,
                PrecisionTYPE * V,
                PrecisionTYPE * M_MAP,
                PrecisionTYPE * V_MAP,
                ImageSize * CurrSizes,
                int verbose)
{

    if(verbose>0){
        cout<< "Optimising Gaussian Parameters" << endl;
        flush(cout);
    }
    int numel_masked=CurrSizes->numelmasked;
    int num_class=CurrSizes->numclass;
    int * S2L_PTR = (int *) S2L;
    PrecisionTYPE * Expec_PTR = (PrecisionTYPE *) Expec;
    PrecisionTYPE * T1_PTR = static_cast<PrecisionTYPE *>(T1->data);
    PrecisionTYPE * BiasField_PTR= (PrecisionTYPE *) BiasField;
    PrecisionTYPE * T1_PTR2= static_cast<PrecisionTYPE *>(T1->data);
    PrecisionTYPE * BiasField_PTR2= (PrecisionTYPE *) BiasField;

    int Expec_offset[PV_numbclass];
    for (int cl=0; cl<num_class; cl++) {
        Expec_offset[cl]=cl*numel_masked;
    }
    PrecisionTYPE tempsum= (PrecisionTYPE) 0.0;
    PrecisionTYPE SumPriors= (PrecisionTYPE) 0.0;
    PrecisionTYPE tempexpec= (PrecisionTYPE) 0.0;

    // ***********

    for (int cl=0; cl<num_class; cl++) {
        // MEAN
        tempexpec=0;
        for(int Multispec=0; Multispec<CurrSizes->usize; Multispec++) {
            Expec_PTR=(PrecisionTYPE *) &Expec[Expec_offset[cl]];
            S2L_PTR = (int *) S2L;

            T1_PTR = static_cast<PrecisionTYPE *>(T1->data);
            T1_PTR =&T1_PTR[Multispec*CurrSizes->numel];
            tempsum=(PrecisionTYPE)0.0;
            SumPriors=(PrecisionTYPE)0.0;

            if(BiasField!=NULL){
                BiasField_PTR= &BiasField[Multispec*numel_masked];
                for (int i=0; i<numel_masked; i++, Expec_PTR++,BiasField_PTR++,S2L_PTR++) {
                    tempsum+=(*Expec_PTR) * (T1_PTR[(*S2L_PTR)]+(*BiasField_PTR));
                    SumPriors+=(*Expec_PTR);
                }
            }else{
                for (int i=0; i<numel_masked; i++, Expec_PTR++,S2L_PTR++) {
                    tempsum+=(*Expec_PTR)*(T1_PTR[(*S2L_PTR)]);
                    SumPriors+=(*Expec_PTR);
                }
            }
            if(M_MAP==NULL){

                M[cl*(CurrSizes->usize)+Multispec]=tempsum/SumPriors;
            }
            else{
                M[cl*(CurrSizes->usize)+Multispec]=(tempsum/SumPriors/powf(V[cl*(CurrSizes->usize)+Multispec],2)+M_MAP[cl*(CurrSizes->usize)+Multispec]/powf(V_MAP[cl*(CurrSizes->usize)+Multispec],2))/(1/powf(V[cl*(CurrSizes->usize)+Multispec],2)+1/powf(V_MAP[cl*(CurrSizes->usize)+Multispec],2));
            }

            for(int Multispec2=Multispec; Multispec2<CurrSizes->usize; Multispec2++) {
                S2L_PTR = (int *) S2L;

                T1_PTR = static_cast<PrecisionTYPE *>(T1->data);
                T1_PTR =&T1_PTR[Multispec*CurrSizes->numel];

                T1_PTR2 = static_cast<PrecisionTYPE *>(T1->data);
                T1_PTR2 =&T1_PTR2[Multispec2*CurrSizes->numel];
                float tmpM=M[cl*CurrSizes->usize+Multispec];
                float tmpM2=M[cl*CurrSizes->usize+Multispec2];
                //STD
                tempsum=0;
                Expec_PTR=&Expec[Expec_offset[cl]];
                if(BiasField!=NULL){
                    BiasField_PTR=&BiasField[Multispec*numel_masked];
                    BiasField_PTR2=&BiasField[Multispec2*numel_masked];
                    for (int i=0; i<numel_masked;i++,Expec_PTR++,BiasField_PTR++,BiasField_PTR2++,S2L_PTR++) {
                        tempsum+=(*Expec_PTR) * (T1_PTR[(*S2L_PTR)]+(*BiasField_PTR)-tmpM) * (T1_PTR2[(*S2L_PTR)]+(*BiasField_PTR2)-tmpM2);
                    }
                }else{
                    for (int i=0; i<numel_masked;i++,Expec_PTR++,S2L_PTR++) {
                        tempsum+=(*Expec_PTR) * (T1_PTR[(*S2L_PTR)]-tmpM) * (T1_PTR2[(*S2L_PTR)]-tmpM2);
                    }

                }

                V[cl*CurrSizes->usize*CurrSizes->usize+Multispec+Multispec2*CurrSizes->usize]=tempsum/SumPriors+0.0001f;
                if(Multispec2!=Multispec){
                    V[cl*CurrSizes->usize*CurrSizes->usize+Multispec2+Multispec*CurrSizes->usize]=V[cl*CurrSizes->usize*CurrSizes->usize+Multispec+Multispec2*CurrSizes->usize];
                }
            }
        }
        if(verbose>0){
            if(CurrSizes->usize==1){
                cout.fill('0');
                cout<< "M["<<(int)(cl)<<"]= "<<setw(10)<<setprecision(7)<<left<<(PrecisionTYPE)(M[cl])<<"\tV["<<(int)(cl)<<"]="<<setw(10)<<setprecision(7)<<left<<(PrecisionTYPE)(V[cl])<< endl;
                flush(cout);
            }
            else{

                cout<< "M["<<(int)(cl)<<"]= ";
                for(int Multispec=0; Multispec<CurrSizes->usize; Multispec++) {
                    cout<< setw(10)<<setprecision(7)<<left<<(PrecisionTYPE)(M[cl*CurrSizes->usize+Multispec])<<"\t";
                }
                cout<< endl;
                flush(cout);
                cout<< "V["<<(int)(cl)<<"]= ";
                for(int Multispec=0; Multispec<CurrSizes->usize; Multispec++) {
                    if(Multispec>0){
                        cout<< "      ";
                    }
                    for(int Multispec2=0; Multispec2<CurrSizes->usize; Multispec2++) {
                        cout<< setw(10)<<setprecision(7)<<left<<(PrecisionTYPE)(V[cl*CurrSizes->usize*CurrSizes->usize+Multispec*CurrSizes->usize+Multispec2])<<"\t";
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





void calcM_mask_LoAd(nifti_image * T1,
                     PrecisionTYPE * Expec,
                     PrecisionTYPE * BiasField,
                     int * S2L,
                     PrecisionTYPE * M,
                     PrecisionTYPE * V,
                     ImageSize * CurrSizes,
                     int verbose,
                     bool PVon)
{

    if(verbose>0){
        cout<< "Optimising Gaussian Parameters" << endl;
        flush(cout);
    }
    int numel_masked=CurrSizes->numelmasked;
    int currentnum_class=CurrSizes->numclass;
    int * S2L_PTR = (int *) S2L;
    PrecisionTYPE * Expec_PTR = (PrecisionTYPE *) Expec;
    PrecisionTYPE * T1_PTR = static_cast<PrecisionTYPE *>(T1->data);
    PrecisionTYPE * BiasField_PTR= (PrecisionTYPE *) BiasField;
    int Expec_offset[PV_numbclass];
    for (int cl=0; cl<currentnum_class; cl++) {
        Expec_offset[cl]=cl*numel_masked;
    }
    PrecisionTYPE tempsum= (PrecisionTYPE) 0.0;
    PrecisionTYPE SumPriors= (PrecisionTYPE) 0.0;
    PrecisionTYPE tempexpec= (PrecisionTYPE) 0.0;


    for (int cl=0; cl<5; cl++) {
        Expec_PTR=(PrecisionTYPE *) &Expec[Expec_offset[cl]];
        BiasField_PTR= (PrecisionTYPE *) BiasField;
        S2L_PTR = (int *) S2L;
        tempsum=(PrecisionTYPE)0.0;
        SumPriors=(PrecisionTYPE)0.0;

        // MEAN
        tempexpec=0;
        if(BiasField!=NULL){
            for (int i=0; i<numel_masked; i++, Expec_PTR++,S2L_PTR++,BiasField_PTR++) {
                tempsum+=(*Expec_PTR) * (T1_PTR[(*S2L_PTR)]+(*BiasField_PTR));
                SumPriors+=(*Expec_PTR);
            }
        }else{
            for (int i=0; i<numel_masked; i++, Expec_PTR++,S2L_PTR++) {
                tempsum+=(*Expec_PTR) * (T1_PTR[(*S2L_PTR)]);
                SumPriors+=(*Expec_PTR);
            }
        }

        M[cl]=tempsum/SumPriors;

        //STD
        tempsum=0;
        Expec_PTR=&Expec[Expec_offset[cl]];
        BiasField_PTR=BiasField;
        S2L_PTR = (int *) S2L;
        if(BiasField!=NULL){
            for (int i=0; i<numel_masked;i++,Expec_PTR++,S2L_PTR++,BiasField_PTR++) {
                tempsum+=(*Expec_PTR) * (T1_PTR[(*S2L_PTR)]+(*BiasField_PTR)-M[cl]) * (T1_PTR[(*S2L_PTR)]+(*BiasField_PTR)-M[cl]);
            }
        }else{
            for (int i=0; i<numel_masked;i++,Expec_PTR++,S2L_PTR++) {
                tempsum+=(*Expec_PTR) * (T1_PTR[(*S2L_PTR)]-M[cl]) * (T1_PTR[(*S2L_PTR)]-M[cl]);
            }

        }
        V[cl]=tempsum/SumPriors;

        if(V[cl]<0.01)V[cl]=0.01;

        if(verbose>0){
            cout.fill('0');
            cout<< "M["<<(int)(cl)<<"]= "<<setw(10)<<setprecision(7)<<left<<(PrecisionTYPE)(M[cl])<<"\tV["<<(int)(cl)<<"]="<<setw(10)<<setprecision(7)<<left<<(PrecisionTYPE)(V[cl])<< endl;
            flush(cout);
        }
    }

    if(PVon){
        M[WMGMpvclass]=(M[WMclass]+M[GMclass])/(float)(2.0);
        M[GMCSFpvclass]=(M[GMclass]+M[CSFclass])/(float)(2.0);
        V[WMGMpvclass]=sqrt(0.5*0.5*pow(V[WMclass],2)+(0.5*0.5)*powf(V[GMclass],2));
        if(V[WMGMpvclass]<0.003)V[WMGMpvclass]=0.003;
        V[GMCSFpvclass]=sqrt(0.5*0.5*pow(V[GMclass],2)+(0.5*0.5)*powf(V[CSFclass],2));
        if(V[GMCSFpvclass]<0.003)V[GMCSFpvclass]=0.003;
        cout<< "M["<<(int)(WMGMpvclass)<<"]= "<<setw(10)<<setprecision(7)<<left<<(PrecisionTYPE)(M[WMGMpvclass])<<"\tV["<<(int)(WMGMpvclass)<<"]="<<setw(10)<<setprecision(7)<<left<<(PrecisionTYPE)(V[WMGMpvclass])<< endl;
        cout<< "M["<<(int)(GMCSFpvclass)<<"]= "<<setw(10)<<setprecision(7)<<left<<(PrecisionTYPE)(M[GMCSFpvclass])<<"\tV["<<(int)(GMCSFpvclass)<<"]="<<setw(10)<<setprecision(7)<<left<<(PrecisionTYPE)(V[GMCSFpvclass])<< endl;
        flush(cout);
    }

    return;
}

void printloglik(int iter,
                 PrecisionTYPE loglik,
                 PrecisionTYPE oldloglik){
    if(iter>1){
        if ((loglik-oldloglik)/fabsf(oldloglik)>0 && (loglik-oldloglik)/fabsf(oldloglik)<100){
            cout<< "Loglik = " << loglik << " : Ratio = " << (loglik-oldloglik)/fabsf(oldloglik) << endl;
        }
        else{
            cout<< "Loglik = " << loglik << endl;
        }
    }
    else {
        cout<< "Loglik = " << loglik << endl;
    }
    return;
}

void Relax_Priors(PrecisionTYPE * Priors,
                  PrecisionTYPE * Expec,
                  PrecisionTYPE * MRF,
                  int * S2L,
                  int * L2S,
                  float RelaxFactor,
                  PrecisionTYPE * G,
                  PrecisionTYPE ba,
                  PrecisionTYPE be,
                  ImageSize *  CurrSizes,
                  SEG_PARAM * segment_param){

    if(RelaxFactor<0){
        RelaxFactor=0;
        if(segment_param->verbose_level>0){
            cout<< "RelaxFactor was bellow 0 (min). Assuming RelaxFactor=0." << endl;
        }
    }

    if(RelaxFactor>1){
        RelaxFactor=1;
        if(segment_param->verbose_level>0){
            cout<< "RelaxFactor was above 1 (max). Assuming RelaxFactor=1." << endl;
        }
    }


    if(RelaxFactor>=0){
        if(segment_param->verbose_level>0){
            cout<< "Relaxing Priors with Rf = "<< RelaxFactor << endl;
            flush(cout);
        }
        Gaussian_Filter_Short_4D(Expec,S2L,L2S,2.0, CurrSizes,CSFclass);
        Relax_Priors_Share(Priors,Expec,RelaxFactor,G,be,CurrSizes);
        for(int i=0; i<(CurrSizes->numclass*CurrSizes->numelmasked);i++)MRF[i]=Priors[i];
        memcpy(Priors,Expec,CurrSizes->numelmasked*CurrSizes->numclass*sizeof(PrecisionTYPE));
    }

    return;

}

void Relax_Priors_Share(PrecisionTYPE * Priors,
                        PrecisionTYPE * Expec,
                        float RelaxFactor,
                        PrecisionTYPE * G,
                        PrecisionTYPE be,
                        ImageSize * CurrSizes){
    int numclass=CurrSizes->numclass;

    int N[PV_numbclass*PV_numbclass];
    PrecisionTYPE Currexpec[PV_numbclass];
    PrecisionTYPE Currexpec_updated[PV_numbclass];
    PrecisionTYPE Currexpec_sumall;
    PrecisionTYPE Currexpec_class;
    int Class_shift[PV_numbclass];
    // Precompute Class Shift
    for (int cl=0; cl<numclass; cl++) {
        Class_shift[cl]=cl*CurrSizes->numelmasked;
    }

    //Find neighbouring rules
    for(int i=0; i<numclass; i++){
        for(int j=0; j<numclass; j++){
            N[i+j*numclass]=(G[i+j*numclass]<=be)?1:0;
        }
    }

    for (int i=0; i<CurrSizes->numelmasked;i++) {
        //Copy current expec for all classes
        for (int cl=0; cl<numclass; cl++) {
            Currexpec[cl]=(1-RelaxFactor)*Expec[i+Class_shift[cl]]+(RelaxFactor)*Priors[i+Class_shift[cl]];
        }
        //Compute relaxed expec
        Currexpec_sumall=0;
        for (int cl=0; cl<numclass; cl++) {
            Currexpec_updated[cl]=0;
            Currexpec_class=Currexpec[cl];
            for (int cl2=0; cl2<numclass; cl2++) {
                Currexpec_updated[cl]+=Currexpec_class*Currexpec[cl2]*double(N[cl+cl2*numclass]);
            }
            Currexpec_sumall+=Currexpec_updated[cl];
        }
        Currexpec_sumall=Currexpec_sumall==0?1:Currexpec_sumall;
        //Copy current expec back all classes
        for (int cl=0; cl<numclass; cl++) {
            Expec[i+Class_shift[cl]]=Currexpec_updated[cl]/Currexpec_sumall;
        }

    }
}

void Gaussian_Filter_Short_4D(PrecisionTYPE * ShortData,
                              int * S2L,
                              int * L2S,
                              PrecisionTYPE gauss_std,
                              ImageSize * CurrSizes,
                              int class_with_CSF){

    int kernelsize=0;
    int kernelsizemin=(int)floorf(gauss_std*3.0);
    int kernelsizemax=(int)ceilf(gauss_std*3.0);

    if((kernelsizemin/2.0)==(double)(kernelsizemin/2) && kernelsizemin!=kernelsizemax){ // Which one is odd? kernelsizemin or kernelsizemax?
        kernelsize=kernelsizemax;}
    else if((kernelsizemax/2.0)==(double)(kernelsizemax/2) && kernelsizemin!=kernelsizemax){
        kernelsize=kernelsizemin;}
    else{
        kernelsize=kernelsizemin+1;}

    int kernelshift=(int)floorf(kernelsize/2);
    PrecisionTYPE GaussKernel [100]= {0};

    for(int i=0; i<kernelsize; i++){
        float kernelvalue=expf((float)(-0.5*powf((i-kernelshift)/gauss_std, 2)))/(sqrtf(2*3.14159265*powf(gauss_std, 2)));
        GaussKernel[i]=kernelvalue;
    }

    PrecisionTYPE * Buffer= new PrecisionTYPE [CurrSizes->numel]();
    PrecisionTYPE * LongData= new PrecisionTYPE [CurrSizes->numel]();


    int shiftdirection[3];
    shiftdirection[0]=1;
    shiftdirection[1]=(int)CurrSizes->xsize;
    shiftdirection[2]=(int)CurrSizes->xsize*(int)CurrSizes->ysize;
    int dim_array[3];
    dim_array[0]=(int)CurrSizes->xsize;
    dim_array[1]=(int)CurrSizes->ysize;
    dim_array[2]=(int)CurrSizes->zsize;


    //Outside the mask is considered Pure CSF
    int outsiderangevalue[10];
    for(int i=0; i<10; i++){
        outsiderangevalue[i]=0;
    }
    outsiderangevalue[class_with_CSF]=1;

    for(int curr4d=0; curr4d<CurrSizes->numclass; curr4d++){ //For Each Class
        int current_4dShift_short=curr4d*CurrSizes->numelmasked;
        for(int index=0;index<CurrSizes->numelmasked;index++){ //Copy Class to Buffer in LongFormat
            Buffer[S2L[index]]=ShortData[index+current_4dShift_short];
        }

        PrecisionTYPE tempsum=0;
        int xyzpos[3];
        for(int currentdirection=0;currentdirection<3;currentdirection++){ //Blur Buffer along each direction
            int index=0;
            for(xyzpos[2]=0;xyzpos[2]<(int)CurrSizes->zsize;xyzpos[2]++){
                for(xyzpos[1]=0;xyzpos[1]<(int)CurrSizes->ysize;xyzpos[1]++){
                    for(xyzpos[0]=0;xyzpos[0]<(int)CurrSizes->xsize;xyzpos[0]++){
                        PrecisionTYPE tmpvalue=0.0f;
                        PrecisionTYPE tmpkernelsum=0.0f;
                        LongData[index]=0.0f;
                        if(L2S[index]>=0){
                            for(int shift=((xyzpos[currentdirection]<kernelshift)?-xyzpos[currentdirection]:-kernelshift);shift<=((xyzpos[currentdirection]>=(dim_array[currentdirection]-kernelshift))?(int)dim_array[currentdirection]-xyzpos[currentdirection]-1:kernelshift) ; shift++){
                                tmpvalue+=(L2S[index+shift*shiftdirection[currentdirection]]==-1)?GaussKernel[shift+kernelshift]*outsiderangevalue[curr4d]:GaussKernel[shift+kernelshift]*Buffer[index+shift*shiftdirection[currentdirection]];
                                tmpkernelsum+=GaussKernel[shift+kernelshift];
                            }
                            LongData[index]=tmpvalue/tmpkernelsum;
                        }
                        index++;
                    }
                }
            }
            if(currentdirection<2){
                for(int index2=0;index2<CurrSizes->numel;index2++){
                    Buffer[index2]=LongData[index2];
                }
            }
        }

        for(int index=0;index<CurrSizes->numelmasked;index++){ //Copy Class to Buffer in LongFormat
            ShortData[index+current_4dShift_short]=LongData[S2L[index]];
        }


    }
    delete [] LongData;
    delete [] Buffer;
    return;
}

PrecisionTYPE * Gaussian_Filter_4D(PrecisionTYPE * LongData,
                                   PrecisionTYPE gauss_std,
                                   ImageSize * CurrSizes){

    int kernelsize=0;
    int kernelsizemin=(int)floorf(gauss_std*3.0);
    int kernelsizemax=(int)ceilf(gauss_std*3.0);

    if((kernelsizemin/2.0)==(double)(kernelsizemin/2) && kernelsizemin!=kernelsizemax){ // Which one is odd? kernelsizemin or kernelsizemax?
        kernelsize=kernelsizemax;}
    else if((kernelsizemax/2.0)==(double)(kernelsizemax/2) && kernelsizemin!=kernelsizemax){
        kernelsize=kernelsizemin;}
    else{
        kernelsize=kernelsizemin+1;}

    int kernelshift=(int)floorf(kernelsize/2);
    PrecisionTYPE GaussKernel [100]= {0};

    for(int i=0; i<kernelsize; i++){
        float kernelvalue=expf((float)(-0.5*powf((i-kernelshift)/gauss_std, 2)))/(sqrtf(2*3.14159265*powf(gauss_std, 2)));
        GaussKernel[i]=kernelvalue;
    }

    PrecisionTYPE * Buffer= new PrecisionTYPE [CurrSizes->numel]();



    int shiftdirection[3];
    shiftdirection[0]=1;
    shiftdirection[1]=(int)CurrSizes->xsize;
    shiftdirection[2]=(int)CurrSizes->xsize*(int)CurrSizes->ysize;
    int dim_array[3];
    dim_array[0]=(int)CurrSizes->xsize;
    dim_array[1]=(int)CurrSizes->ysize;
    dim_array[2]=(int)CurrSizes->zsize;


    //Outside the mask is considered Pure CSF
    int outsiderangevalue[10];
    for(int i=0; i<10; i++){
        outsiderangevalue[i]=0;
    }
    for(int curr4d=0; curr4d<CurrSizes->numclass; curr4d++){ //For Each Class

        int current_4dShift_short=curr4d*CurrSizes->numel;
        PrecisionTYPE * longdataptr=&LongData[current_4dShift_short];
        for(int index=0;index<CurrSizes->numel;index++){ //Copy Class to Buffer in LongFormat
            Buffer[index]=longdataptr[index];
        }

        PrecisionTYPE tempsum=0;
        int xyzpos[3];
        for(int currentdirection=0;currentdirection<3;currentdirection++){ //Blur Buffer along each direction
            int index=0;
            for(xyzpos[2]=0;xyzpos[2]<(int)CurrSizes->zsize;xyzpos[2]++){
                for(xyzpos[1]=0;xyzpos[1]<(int)CurrSizes->ysize;xyzpos[1]++){
                    for(xyzpos[0]=0;xyzpos[0]<(int)CurrSizes->xsize;xyzpos[0]++){
                        PrecisionTYPE tmpvalue=0.0f;
                        PrecisionTYPE tmpkernelsum=0.0f;
                        longdataptr[index]=0.0f;
                        for(int shift=((xyzpos[currentdirection]<kernelshift)?-xyzpos[currentdirection]:-kernelshift);shift<=((xyzpos[currentdirection]>=(dim_array[currentdirection]-kernelshift))?(int)dim_array[currentdirection]-xyzpos[currentdirection]-1:kernelshift) ; shift++){
                            tmpvalue+=GaussKernel[shift+kernelshift]*Buffer[index+shift*shiftdirection[currentdirection]];
                            tmpkernelsum+=GaussKernel[shift+kernelshift];
                        }
                        longdataptr[index]=tmpvalue/tmpkernelsum;
                        index++;
                    }
                }
            }
            if(currentdirection<2){
                for(int index2=0;index2<CurrSizes->numel;index2++){
                    Buffer[index2]=longdataptr[index2];
                }
            }
        }
    }
    return Buffer;
}

PrecisionTYPE * Gaussian_Filter_4D_inside_mask(PrecisionTYPE * LongData,
                                               bool * mask,
                                               PrecisionTYPE gauss_std,
                                               ImageSize * CurrSizes){

    int kernelsize=0;
    int kernelsizemin=(int)floorf(gauss_std*3.0);
    int kernelsizemax=(int)ceilf(gauss_std*3.0);

    if((kernelsizemin/2.0)==(double)(kernelsizemin/2) && kernelsizemin!=kernelsizemax){ // Which one is odd? kernelsizemin or kernelsizemax?
        kernelsize=kernelsizemax;}
    else if((kernelsizemax/2.0)==(double)(kernelsizemax/2) && kernelsizemin!=kernelsizemax){
        kernelsize=kernelsizemin;}
    else{
        kernelsize=kernelsizemin+1;}

    int kernelshift=(int)floorf(kernelsize/2);
    PrecisionTYPE GaussKernel [100]= {0};

    for(int i=0; i<kernelsize; i++){
        float kernelvalue=expf((float)(-0.5*powf((i-kernelshift)/gauss_std, 2)))/(sqrtf(2*3.14159265*powf(gauss_std, 2)));
        GaussKernel[i]=kernelvalue;
    }

    PrecisionTYPE * Buffer= new PrecisionTYPE [CurrSizes->numel]();



    int shiftdirection[3];
    shiftdirection[0]=1;
    shiftdirection[1]=(int)CurrSizes->xsize;
    shiftdirection[2]=(int)CurrSizes->xsize*(int)CurrSizes->ysize;
    int dim_array[3];
    dim_array[0]=(int)CurrSizes->xsize;
    dim_array[1]=(int)CurrSizes->ysize;
    dim_array[2]=(int)CurrSizes->zsize;


    //Outside the mask is considered Pure CSF
    int outsiderangevalue[10];
    for(int i=0; i<10; i++){
        outsiderangevalue[i]=0;
    }
    for(int curr4d=0; curr4d<CurrSizes->numclass; curr4d++){ //For Each Class
        int current_4dShift_short=curr4d*CurrSizes->numel;
        for(int index=0;index<CurrSizes->numel;index++){ //Copy Class to Buffer in LongFormat
            Buffer[index]=LongData[index+current_4dShift_short];
        }

        PrecisionTYPE tempsum=0;
        int xyzpos[3];
        for(int currentdirection=0;currentdirection<3;currentdirection++){ //Blur Buffer along each direction
            int index=0;
            for(xyzpos[2]=0;xyzpos[2]<(int)CurrSizes->zsize;xyzpos[2]++){
                for(xyzpos[1]=0;xyzpos[1]<(int)CurrSizes->ysize;xyzpos[1]++){
                    for(xyzpos[0]=0;xyzpos[0]<(int)CurrSizes->xsize;xyzpos[0]++){
                        PrecisionTYPE tmpvalue=0.0f;
                        PrecisionTYPE tmpkernelsum=0.0f;
                        LongData[index]=0.0f;
                        if(mask[index]>0){
                            for(int shift=((xyzpos[currentdirection]<kernelshift)?-xyzpos[currentdirection]:-kernelshift);shift<=((xyzpos[currentdirection]>=(dim_array[currentdirection]-kernelshift))?(int)dim_array[currentdirection]-xyzpos[currentdirection]-1:kernelshift) ; shift++){
                                tmpvalue+=GaussKernel[shift+kernelshift]*Buffer[index+shift*shiftdirection[currentdirection]];
                                tmpkernelsum+=GaussKernel[shift+kernelshift];
                            }
                            LongData[index]=tmpvalue/tmpkernelsum;
                        }
                        else{
                            if(curr4d==2){
                                LongData[index]=1;
                            }
                            else{
                                LongData[index]=0;
                            }
                        }
                        index++;
                    }
                }
            }
            if(currentdirection<2){
                for(int index2=0;index2<CurrSizes->numel;index2++){
                    Buffer[index2]=LongData[index2];
                }
            }
        }
    }
    return Buffer;
}

void Convert_to_PV(nifti_image * T1,
                   PrecisionTYPE * BiasField,
                   PrecisionTYPE * ShortPrior,
                   PrecisionTYPE * Expec,
                   PrecisionTYPE * MRF,
                   PrecisionTYPE * M,
                   PrecisionTYPE * V,
                   int * S2L,
                   int * L2S,
                   ImageSize * CurrSizes,
                   SEG_PARAM * segment_param){
    if(segment_param->verbose_level>0){
        cout<< "Covert to Explicit PV modeling" << endl;
    }
    Gaussian_Filter_Short_4D(Expec,S2L,L2S,1.0,CurrSizes,CSFclass);
    Gaussian_Filter_Short_4D(MRF,S2L,L2S,1.0,CurrSizes,CSFclass);
    Convert_WM_and_GM_to_PV(T1,BiasField,ShortPrior,Expec,S2L,M,V,CurrSizes);
    CurrSizes->numclass=7;
    for(int i=0; i<(CurrSizes->numclass*CurrSizes->numelmasked);i++)MRF[i]==1.0f;

}

void Sulci_and_gyri_correction(PrecisionTYPE * MRF_Beta,
                               PrecisionTYPE * ShortPrior,
                               PrecisionTYPE * Expec,
                               PrecisionTYPE *MRF,
                               int * S2L,
                               int * L2S,
                               ImageSize *CurrSizes,
                               SEG_PARAM *segment_param){

    // Deep Sulci ->  Seed=WM+WMGMpv+dGM+iCSF
    cout<< "Sucli and Gyri correction" << endl;
    bool * Seed_Mask= new bool [CurrSizes->numelmasked]();
    PrecisionTYPE * SpeedFunc= new PrecisionTYPE [CurrSizes->numelmasked]();
    PrecisionTYPE * wSulci= new PrecisionTYPE [CurrSizes->numelmasked]();
    PrecisionTYPE * wGyri= new PrecisionTYPE [CurrSizes->numelmasked]();

    for(int i=0; i<CurrSizes->numelmasked;i++){
        Seed_Mask[i]=(Expec[i+WMclass*CurrSizes->numelmasked]+
                      Expec[i+WMGMpvclass*CurrSizes->numelmasked]+
                      Expec[i+dGMclass*CurrSizes->numelmasked]+
                      Expec[i+iCSFclass*CurrSizes->numelmasked])>0.5f;
        SpeedFunc[i]=(Expec[i+CSFclass*CurrSizes->numelmasked]+
                      Expec[i+GMCSFpvclass*CurrSizes->numelmasked])*2;
    }

    FMM(Seed_Mask, SpeedFunc, wSulci,30, L2S, S2L, CurrSizes);
    TransformGeoTime(wSulci,30, L2S, S2L, CurrSizes);

    for(int i=0; i<CurrSizes->numelmasked;i++){
        Seed_Mask[i]=(Expec[i+CSFclass*CurrSizes->numelmasked]+
                      Expec[i+GMCSFpvclass*CurrSizes->numelmasked])>0.5f;
        SpeedFunc[i]=(Expec[i+WMclass*CurrSizes->numelmasked]+
                      Expec[i+WMGMpvclass*CurrSizes->numelmasked]+
                      Expec[i+dGMclass*CurrSizes->numelmasked]+
                      Expec[i+iCSFclass*CurrSizes->numelmasked])*2;
    }

    FMM(Seed_Mask, SpeedFunc, wGyri,30, L2S, S2L, CurrSizes);
    TransformGeoTime(wGyri,30, L2S, S2L, CurrSizes);


    float sumed_all=0;
    float toremove_Gyri=0;
    float toremove_Sulci=0;
    float wSulci_tmp=0;
    float wGyri_tmp=0;
    float MRFbeta_tmp=0;

    for(int i=0; i<CurrSizes->numelmasked;i++){
        wSulci_tmp=wSulci[i];
        wGyri_tmp=wGyri[i];
        MRFbeta_tmp=(1-wSulci_tmp)*(1-wGyri_tmp);
        MRF_Beta[i]=MRFbeta_tmp;
        sumed_all=0;
        if((wSulci_tmp+wGyri_tmp)>0){
            toremove_Gyri=Expec[i+GMclass*CurrSizes->numelmasked]*(1-MRFbeta_tmp)*(wGyri_tmp/(wSulci_tmp+wGyri_tmp));
            toremove_Sulci=Expec[i+GMclass*CurrSizes->numelmasked]*(1-MRFbeta_tmp)*(wSulci_tmp/(wSulci_tmp+wGyri_tmp));

            Expec[i+WMGMpvclass*CurrSizes->numelmasked]=Expec[i+WMGMpvclass*CurrSizes->numelmasked]+wGyri[i]*Expec[i+GMclass*CurrSizes->numelmasked];
            Expec[i+GMCSFpvclass*CurrSizes->numelmasked]=Expec[i+GMCSFpvclass*CurrSizes->numelmasked]+wSulci[i]*Expec[i+GMclass*CurrSizes->numelmasked];
            Expec[i+GMclass*CurrSizes->numelmasked]=Expec[i+GMclass*CurrSizes->numelmasked]*(MRFbeta_tmp);

            for(int c=0;c<7;c++){
                sumed_all+=Expec[i+c*CurrSizes->numelmasked];
            }
            if(sumed_all>0){
                for(int c=0;c<7;c++){
                    Expec[i+c*CurrSizes->numelmasked]=Expec[i+c*CurrSizes->numelmasked]/sumed_all;
                }
            }else{
                for(int c=0;c<7;c++){
                    Expec[i+c*CurrSizes->numelmasked]=1.0/7.0;
                }

            }
        }

    }
    Gaussian_Filter_Short_4D(Expec,S2L,L2S,1.0, CurrSizes,CSFclass);

    for(int i=0; i<CurrSizes->numelmasked;i++){
        for(int c=0;c<7;c++){
            ShortPrior[i+c*CurrSizes->numelmasked]=Expec[i+c*CurrSizes->numelmasked];
            MRF[i+c*CurrSizes->numelmasked]=Expec[i+c*CurrSizes->numelmasked];
        }
    }

    delete [] Seed_Mask;
    delete [] SpeedFunc;
    delete [] wSulci;
    delete [] wGyri;
    return;
}

void Convert_WM_and_GM_to_PV(nifti_image * T1,
                             PrecisionTYPE * BiasField,
                             PrecisionTYPE * ShortPrior,
                             PrecisionTYPE * Expec,
                             int * S2L,
                             PrecisionTYPE * M,
                             PrecisionTYPE * V,
                             ImageSize * CurrSizes){

    int WMindex=WMclass*CurrSizes->numelmasked;
    int GMindex=GMclass*CurrSizes->numelmasked;
    int CSFindex=CSFclass*CurrSizes->numelmasked;
    int dGMindex=dGMclass*CurrSizes->numelmasked;
    int iCSFindex=iCSFclass*CurrSizes->numelmasked;
    int WMGMpvindex=WMGMpvclass*CurrSizes->numelmasked;
    int GMCSFpvindex=GMCSFpvclass*CurrSizes->numelmasked;

    float averageFC_WMGM=0;
    float averageFCtmp_WMGM=0;
    float averageFC_GMCSF=0;
    float averageFCtmp_GMCSF=0;
    float currentT1=0;
    int count1=0;
    int count2=0;

    PrecisionTYPE * T1ptrtmp = static_cast<PrecisionTYPE *>(T1->data);
    for (int i=0; i<CurrSizes->numelmasked;i++){
        currentT1=T1ptrtmp[S2L[i]]+BiasField[i];

        averageFCtmp_WMGM=(currentT1-M[GMclass])/(M[WMclass]-M[GMclass]);

        if(averageFCtmp_WMGM<1 && averageFCtmp_WMGM>0){
            averageFC_WMGM+=averageFCtmp_WMGM;
            count1++;
        }

        averageFCtmp_GMCSF=(currentT1-M[CSFclass])/(M[GMclass]-M[CSFclass]);

        if(averageFCtmp_GMCSF<1 && averageFCtmp_GMCSF>0){
            averageFC_GMCSF+=averageFCtmp_GMCSF;
            count2++;
        }
    }

    averageFC_WMGM=(count1>0)?averageFC_WMGM/(float)(count1):0.5;
    averageFC_GMCSF=(count2>0)?averageFC_GMCSF/(float)(count2):0.5;


    cout<< "avr_WMGMfc=" << averageFC_WMGM << endl;
    cout<< "avr_GMCSFfc=" << averageFC_GMCSF << endl;

    for (int i=0; i<CurrSizes->numelmasked;i++,WMindex++,GMindex++,CSFindex++,dGMindex++,iCSFindex++,WMGMpvindex++,GMCSFpvindex++) {
        //Copy current expec for all classes
        ShortPrior[WMindex]=(ShortPrior[WMindex]*Expec[WMindex]);
        Expec[WMindex]=ShortPrior[WMindex];
        ShortPrior[GMindex]=(ShortPrior[GMindex]*Expec[GMindex]);
        Expec[GMindex]=ShortPrior[GMindex];
        ShortPrior[CSFindex]=(ShortPrior[CSFindex]*Expec[CSFindex]);
        Expec[CSFindex]=ShortPrior[CSFindex];
        ShortPrior[dGMindex]=(ShortPrior[dGMindex]*Expec[dGMindex]);
        Expec[dGMindex]=ShortPrior[dGMindex];
        ShortPrior[iCSFindex]=(ShortPrior[iCSFindex]*Expec[iCSFindex]);
        Expec[iCSFindex]=ShortPrior[iCSFindex];
        ShortPrior[WMGMpvindex]=sqrtf(ShortPrior[WMindex]*ShortPrior[GMindex]);
        Expec[WMGMpvindex]=(ShortPrior[WMGMpvindex]==ShortPrior[WMGMpvindex])?ShortPrior[WMGMpvindex]:0;
        ShortPrior[WMGMpvindex]=Expec[WMGMpvindex];
        ShortPrior[GMCSFpvindex]=sqrtf(ShortPrior[CSFindex]*ShortPrior[GMindex]);
        Expec[GMCSFpvindex]=(ShortPrior[GMCSFpvindex]==ShortPrior[GMCSFpvindex])?ShortPrior[GMCSFpvindex]:0;
        ShortPrior[GMCSFpvindex]=Expec[GMCSFpvindex];
    }



    M[WMGMpvclass]=(averageFC_WMGM*M[WMclass]+(1-averageFC_WMGM)*M[GMclass]);
    M[GMCSFpvclass]=(averageFC_GMCSF*M[GMclass]+(1-averageFC_GMCSF)*M[CSFclass]);
    V[WMGMpvclass]=sqrtf(averageFC_WMGM*averageFC_WMGM*pow(V[WMclass],2)+(1-averageFC_WMGM)*(1-averageFC_WMGM)*V[GMclass]);
    V[GMCSFpvclass]=sqrtf(averageFC_GMCSF*averageFC_GMCSF*pow(V[GMclass],2)+(1-averageFC_GMCSF)*(1-averageFC_GMCSF)*V[CSFclass]);

}



nifti_image * Copy_ShortExpec_to_Result(nifti_image * T1,
                                        PrecisionTYPE * Expec,
                                        PrecisionTYPE * BiasField,
                                        PrecisionTYPE * BiasFieldCoefs,
                                        int * S2L,
                                        nifti_image * Priors,
                                        SEG_PARAM * segment_param,
                                        PrecisionTYPE * M,
                                        ImageSize * CurrSizes){

    nifti_image * Result = nifti_copy_nim_info(Priors);
    if(segment_param->flag_PV_model){

        PrecisionTYPE * T1ptrtmp = static_cast<PrecisionTYPE *>(T1->data);

        Result->dim[0]=4;
        Result->dim[4]=6;
        Result->datatype=DT_FLOAT32;
        Result->cal_max=1;
        nifti_set_filenames(Result,segment_param->filename_out,0,0);
        nifti_update_dims_from_array(Result);
        nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
        Result->data = (void *) calloc(Result->nvox, sizeof(PrecisionTYPE));
        PrecisionTYPE * Result_PTR = static_cast<PrecisionTYPE *>(Result->data);
        for(int i=0; i<Result->nvox; i++){Result_PTR[i]=0;}

        int Short_2_Long_Indices_tmp = 0;
        int class_nvox=Result->nx*Result->ny*Result->nz;
        PrecisionTYPE * Resultdata= static_cast<PrecisionTYPE *>(Result->data);
        PrecisionTYPE * Expec_tmp = new PrecisionTYPE [CurrSizes->numclass]();
        PrecisionTYPE Expec_tmp_sum=0;
        PrecisionTYPE Fractional_content=0;


        PrecisionTYPE * T1data = static_cast<PrecisionTYPE *>(T1->data);
        float to_resize=0;
        if(BiasField!=NULL){
            for(int i=0; i<CurrSizes->numelmasked; i++){
                to_resize=exp((BiasField[i]+T1data[S2L[i]])*0.693147181)-1;
                Resultdata[S2L[i]]=(to_resize*(CurrSizes->rescale_max[0]-CurrSizes->rescale_min[0])+CurrSizes->rescale_min[0]);
            }
        }else{
            for(int i=0; i<CurrSizes->numelmasked; i++){
                to_resize=exp((T1data[S2L[i]])*0.693147181)-1;
                Resultdata[S2L[i]]=(to_resize*(CurrSizes->rescale_max[0]-CurrSizes->rescale_min[0])+CurrSizes->rescale_min[0]);
            }
        }




        if(BiasField!=NULL){
            for(int i=0; i<CurrSizes->numelmasked; i++){
                Short_2_Long_Indices_tmp=S2L[i];
                for(int currclass=0; currclass<CurrSizes->numclass;currclass++){
                    Expec_tmp[currclass]=Expec[i+currclass*CurrSizes->numelmasked];
                }
                float classthreshold=0.1;
                Resultdata[Short_2_Long_Indices_tmp+(WMclass+1)*class_nvox]=Expec_tmp[WMclass]>classthreshold?Expec_tmp[WMclass]:0;
                Resultdata[Short_2_Long_Indices_tmp+(GMclass+1)*class_nvox]=Expec_tmp[GMclass]>classthreshold?Expec_tmp[GMclass]:0;
                Resultdata[Short_2_Long_Indices_tmp+(CSFclass+1)*class_nvox]=Expec_tmp[CSFclass]>classthreshold?Expec_tmp[CSFclass]:0;
                Resultdata[Short_2_Long_Indices_tmp+(dGMclass+1)*class_nvox]=Expec_tmp[dGMclass]>classthreshold?Expec_tmp[dGMclass]:0;
                Resultdata[Short_2_Long_Indices_tmp+(iCSFclass+1)*class_nvox]=Expec_tmp[iCSFclass]>classthreshold?Expec_tmp[iCSFclass]:0;


                // Estimating WM/GM fractional content from the T1 bias corrected data
                if(Expec_tmp[WMGMpvclass]>classthreshold){
                    if((Expec_tmp[WMclass])<(Expec_tmp[WMGMpvclass]) & Expec_tmp[GMclass]<(Expec_tmp[WMGMpvclass])){
                        Expec_tmp_sum=Expec_tmp[WMclass]+Expec_tmp[WMGMpvclass]+Expec_tmp[GMclass];
                        Fractional_content=(M[WMclass]-(T1ptrtmp[Short_2_Long_Indices_tmp]+BiasField[i]))/(M[WMclass]-M[GMclass]);
                        Fractional_content=(Fractional_content<0)?0:Fractional_content;
                        Fractional_content=(Fractional_content>1)?1:Fractional_content;
                        Resultdata[Short_2_Long_Indices_tmp+(WMclass+1)*class_nvox]=(1-Fractional_content)/Expec_tmp_sum;
                        Resultdata[Short_2_Long_Indices_tmp+(GMclass+1)*class_nvox]=(Fractional_content)/Expec_tmp_sum;
                    }
                    else if(Expec_tmp[WMclass]>Expec_tmp[GMclass]){
                        Resultdata[Short_2_Long_Indices_tmp+(WMclass+1)*class_nvox]=1;
                        Resultdata[Short_2_Long_Indices_tmp+(GMclass+1)*class_nvox]=0;
                    }
                    else{
                        Resultdata[Short_2_Long_Indices_tmp+(WMclass+1)*class_nvox]=0;
                        Resultdata[Short_2_Long_Indices_tmp+(GMclass+1)*class_nvox]=1;
                    }
                }

                // Estimating CSF/GM fractional content from the T1 bias corrected data
                if(Expec_tmp[GMCSFpvclass]>classthreshold){
                    if(Expec_tmp[CSFclass]<(Expec_tmp[GMCSFpvclass]) & Expec_tmp[GMclass]<(Expec_tmp[GMCSFpvclass])){
                        Expec_tmp_sum=Expec_tmp[CSFclass]+Expec_tmp[GMCSFpvclass]+Expec_tmp[GMclass];
                        Fractional_content=(M[CSFclass]-(T1ptrtmp[Short_2_Long_Indices_tmp]+BiasField[i]))/(M[CSFclass]-M[GMclass]);
                        Fractional_content=(Fractional_content<0)?0:Fractional_content;
                        Fractional_content=(Fractional_content>1)?1:Fractional_content;
                        Resultdata[Short_2_Long_Indices_tmp+(CSFclass+1)*class_nvox]=(1-Fractional_content)/Expec_tmp_sum;
                        Resultdata[Short_2_Long_Indices_tmp+(GMclass+1)*class_nvox]=(Fractional_content)/Expec_tmp_sum;
                    }
                    else if(Expec_tmp[CSFclass]>Expec_tmp[GMclass]){
                        Resultdata[Short_2_Long_Indices_tmp+(CSFclass+1)*class_nvox]=1;
                        Resultdata[Short_2_Long_Indices_tmp+(GMclass+1)*class_nvox]=0;
                    }
                    else{
                        Resultdata[Short_2_Long_Indices_tmp+(CSFclass+1)*class_nvox]=0;
                        Resultdata[Short_2_Long_Indices_tmp+(GMclass+1)*class_nvox]=1;
                    }
                }

                // Normalizing the fractional contents
                Expec_tmp_sum=0;
                for(int currclass=0; currclass<non_PV_numclass;currclass++){
                    Expec_tmp_sum+=Resultdata[Short_2_Long_Indices_tmp+(currclass+1)*class_nvox];
                }
                for(int currclass=0; currclass<non_PV_numclass;currclass++){
                    Resultdata[Short_2_Long_Indices_tmp+(currclass+1)*class_nvox]=Resultdata[Short_2_Long_Indices_tmp+(currclass+1)*class_nvox]/Expec_tmp_sum;
                }

            }
        }
        else{

            for(int i=0; i<CurrSizes->numelmasked; i++){
                Short_2_Long_Indices_tmp=S2L[i];
                for(int currclass=0; currclass<CurrSizes->numclass;currclass++){
                    Expec_tmp[currclass]=Expec[i+currclass*CurrSizes->numelmasked];
                }
                float classthreshold=0.1;
                Resultdata[Short_2_Long_Indices_tmp+(WMclass+1)*class_nvox]=Expec_tmp[WMclass]>classthreshold?Expec_tmp[WMclass]:0;
                Resultdata[Short_2_Long_Indices_tmp+(GMclass+1)*class_nvox]=Expec_tmp[GMclass]>classthreshold?Expec_tmp[GMclass]:0;
                Resultdata[Short_2_Long_Indices_tmp+(CSFclass+1)*class_nvox]=Expec_tmp[CSFclass]>classthreshold?Expec_tmp[CSFclass]:0;
                Resultdata[Short_2_Long_Indices_tmp+(dGMclass+1)*class_nvox]=Expec_tmp[dGMclass]>classthreshold?Expec_tmp[dGMclass]:0;
                Resultdata[Short_2_Long_Indices_tmp+(iCSFclass+1)*class_nvox]=Expec_tmp[iCSFclass]>classthreshold?Expec_tmp[iCSFclass]:0;


                // Estimating WM/GM fractional content from the T1 data
                if(Expec_tmp[WMGMpvclass]>classthreshold){
                    if((Expec_tmp[WMclass])<(Expec_tmp[WMGMpvclass]) & Expec_tmp[GMclass]<(Expec_tmp[WMGMpvclass])){
                        Expec_tmp_sum=Expec_tmp[WMclass]+Expec_tmp[WMGMpvclass]+Expec_tmp[GMclass];
                        Fractional_content=(M[WMclass]-(T1ptrtmp[Short_2_Long_Indices_tmp]))/(M[WMclass]-M[GMclass]);
                        Fractional_content=(Fractional_content<0)?0:Fractional_content;
                        Fractional_content=(Fractional_content>1)?1:Fractional_content;
                        Resultdata[Short_2_Long_Indices_tmp+(WMclass+1)*class_nvox]=(1-Fractional_content)/Expec_tmp_sum;
                        Resultdata[Short_2_Long_Indices_tmp+(GMclass+1)*class_nvox]=(Fractional_content)/Expec_tmp_sum;
                    }
                    else if(Expec_tmp[WMclass]>Expec_tmp[GMclass]){
                        Resultdata[Short_2_Long_Indices_tmp+(WMclass+1)*class_nvox]=1;
                        Resultdata[Short_2_Long_Indices_tmp+(GMclass+1)*class_nvox]=0;
                    }
                    else{
                        Resultdata[Short_2_Long_Indices_tmp+(WMclass+1)*class_nvox]=0;
                        Resultdata[Short_2_Long_Indices_tmp+(GMclass+1)*class_nvox]=1;
                    }
                }

                // Estimating CSF/GM fractional content from the T1 data
                if(Expec_tmp[GMCSFpvclass]>classthreshold){
                    if(Expec_tmp[CSFclass]<Expec_tmp[GMCSFpvclass] & Expec_tmp[GMclass]<Expec_tmp[GMCSFpvclass]){
                        Expec_tmp_sum=Expec_tmp[CSFclass]+Expec_tmp[GMCSFpvclass]+Expec_tmp[GMclass];
                        Fractional_content=(M[CSFclass]-(T1ptrtmp[Short_2_Long_Indices_tmp]))/(M[CSFclass]-M[GMclass]);
                        Fractional_content=(Fractional_content<0)?0:Fractional_content;
                        Fractional_content=(Fractional_content>1)?1:Fractional_content;
                        Resultdata[Short_2_Long_Indices_tmp+(CSFclass+1)*class_nvox]=(1-Fractional_content)/Expec_tmp_sum;
                        Resultdata[Short_2_Long_Indices_tmp+(GMclass+1)*class_nvox]=(Fractional_content)/Expec_tmp_sum;
                    }
                    else if(Expec_tmp[CSFclass]>Expec_tmp[GMclass]){
                        Resultdata[Short_2_Long_Indices_tmp+(CSFclass+1)*class_nvox]=1;
                        Resultdata[Short_2_Long_Indices_tmp+(GMclass+1)*class_nvox]=0;
                    }
                    else{
                        Resultdata[Short_2_Long_Indices_tmp+(CSFclass+1)*class_nvox]=0;
                        Resultdata[Short_2_Long_Indices_tmp+(GMclass+1)*class_nvox]=1;
                    }
                }

                // Normalizing the fractional contents
                Expec_tmp_sum=0;
                for(int currclass=0; currclass<non_PV_numclass;currclass++){
                    Expec_tmp_sum+=Resultdata[Short_2_Long_Indices_tmp+(currclass+1)*class_nvox];
                }
                for(int currclass=0; currclass<non_PV_numclass;currclass++){
                    Resultdata[Short_2_Long_Indices_tmp+(currclass+1)*class_nvox]=Resultdata[Short_2_Long_Indices_tmp+(currclass+1)*class_nvox]/Expec_tmp_sum;
                }



            }
            delete [] Expec_tmp;


        }
    }
    else{
        Result->dim[0]=4;
        Result->dim[4]=(CurrSizes->numclass+1);
        Result->datatype=DT_FLOAT32;
        Result->cal_max=1;
        nifti_set_filenames(Result,segment_param->filename_out,0,0);
        nifti_update_dims_from_array(Result);
        nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
        Result->data = (void *) calloc(Result->nvox, sizeof(PrecisionTYPE));
        PrecisionTYPE * Result_PTR = static_cast<PrecisionTYPE *>(Result->data);
        for(int i=0; i<Result->nvox; i++){Result_PTR[i]=0;}

        int * S2L_PRT = (int *) S2L;
        PrecisionTYPE * Resultdata= static_cast<PrecisionTYPE *>(Result->data);


        int class_nvox=CurrSizes->numel;
        for(int currclass=0; currclass<CurrSizes->numclass;currclass++){

            PrecisionTYPE * Resultdata_class = &Resultdata[(currclass+1)*class_nvox];
            PrecisionTYPE * Expec_PTR = &Expec[(currclass)*CurrSizes->numelmasked];
            S2L_PRT= (int *) S2L;

            for(int i=0; i<CurrSizes->numelmasked; i++,S2L_PRT++,Expec_PTR++){
                Resultdata_class[*S2L_PRT]=*Expec_PTR;
            }
        }
    }
    return Result;
}

bool * binarise_image(PrecisionTYPE * SingleImage,
                      PrecisionTYPE Threshold,
                      ImageSize * CurrSizes){

    bool * Result = new bool [CurrSizes->numelmasked];

    for(int i=0; i<CurrSizes->numelmasked; i++){
        Result[i]=(SingleImage[i]>Threshold);
    }
    return Result;
}

nifti_image * Copy_Single_ShortImage_to_Result(PrecisionTYPE * SingleImage,
                                               int * Short_2_Long_Indices,
                                               nifti_image * Sourceimage,
                                               char * filename,
                                               ImageSize * CurrSizes){

    nifti_image * Result = nifti_copy_nim_info(Sourceimage);
    nifti_set_filenames(Result,filename,0,0);
    Result->data = (void *) calloc(Result->nvox, sizeof(PrecisionTYPE));
    PrecisionTYPE * Result_PTR = static_cast<PrecisionTYPE *>(Result->data);
    for(int i=0; i<Result->nvox; i++){Result_PTR[i]=0;}

    int * Short_2_Long_Indices_PRT = (int *) Short_2_Long_Indices;
    PrecisionTYPE * Resultdata= static_cast<PrecisionTYPE *>(Result->data);

    PrecisionTYPE * SingleImage_PTR =SingleImage;
    Short_2_Long_Indices_PRT= (int *) Short_2_Long_Indices;

    for(int i=0; i<CurrSizes->numelmasked; i++,Short_2_Long_Indices_PRT++,SingleImage_PTR++){
        Resultdata[*Short_2_Long_Indices_PRT]=(PrecisionTYPE)(*SingleImage_PTR);
    }
    return Result;
}



nifti_image * Copy_BiasCorrected_to_Result_mask(PrecisionTYPE * BiasField,
                                                int * Short_2_Long_Indices,
                                                nifti_image * T1,
                                                char * filename,
                                                ImageSize * CurrSizes){

    nifti_image * Result = nifti_copy_nim_info(T1);
    Result->dim[0]=4;
    Result->dim[4]=1;
    Result->datatype=DT_FLOAT32;
    Result->cal_max=1;
    nifti_set_filenames(Result,filename,0,0);
    nifti_update_dims_from_array(Result);
    nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
    Result->data = (void *) calloc(Result->nvox, sizeof(PrecisionTYPE));
    PrecisionTYPE * Resultdata = static_cast<PrecisionTYPE *>(Result->data);
    PrecisionTYPE * T1data = static_cast<PrecisionTYPE *>(T1->data);
    for(int i=0; i<Result->nvox; i++){Resultdata[i]=0;}

    int * Short_2_Long_Indices_PRT = (int *) Short_2_Long_Indices;

    int class_nvox=Result->nx*Result->ny*Result->nz;

    Short_2_Long_Indices_PRT= (int *) Short_2_Long_Indices;
    PrecisionTYPE to_resize=0;
    if(BiasField!=NULL){
        for(int i=0; i<CurrSizes->numelmasked; i++,Short_2_Long_Indices_PRT++){
            to_resize=exp((BiasField[i]+T1data[*Short_2_Long_Indices_PRT])*0.693147181)-1;
            Resultdata[*Short_2_Long_Indices_PRT]=(to_resize*(CurrSizes->rescale_max[0]-CurrSizes->rescale_min[0])+CurrSizes->rescale_min[0]);
        }
    }else{
        for(int i=0; i<CurrSizes->numelmasked; i++,Short_2_Long_Indices_PRT++){
            to_resize=exp((T1data[*Short_2_Long_Indices_PRT])*0.693147181)-1;
            Resultdata[*Short_2_Long_Indices_PRT]=(to_resize*(CurrSizes->rescale_max[0]-CurrSizes->rescale_min[0])+CurrSizes->rescale_min[0]);
        }
    }
    return Result;
}

nifti_image * Copy_BiasCorrected_to_Result(PrecisionTYPE * BiasField,
                                           nifti_image * T1,
                                           char * filename,
                                           ImageSize * CurrSizes){

    nifti_image * Result = nifti_copy_nim_info(T1);
    Result->dim[0]=4;
    Result->dim[4]=1;
    Result->datatype=DT_FLOAT32;
    Result->cal_max=1;
    nifti_set_filenames(Result,filename,0,0);
    nifti_update_dims_from_array(Result);
    nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
    Result->data = (void *) calloc(Result->nvox, sizeof(PrecisionTYPE));
    PrecisionTYPE * Resultdata = static_cast<PrecisionTYPE *>(Result->data);
    PrecisionTYPE * T1data = static_cast<PrecisionTYPE *>(T1->data);
    for(int i=0; i<Result->nvox; i++){Resultdata[i]=0;}
    int class_nvox=Result->nx*Result->ny*Result->nz;
    PrecisionTYPE to_resize=0;
    if(BiasField!=NULL){
        for(int i=0; i<CurrSizes->numel; i++){
            to_resize=exp((BiasField[i]+T1data[i])*0.693147181)-1;
            Resultdata[i]=(to_resize*(CurrSizes->rescale_max[0]-CurrSizes->rescale_min[0])+CurrSizes->rescale_min[0]);
        }
    }else{
        for(int i=0; i<CurrSizes->numel; i++){
            to_resize=exp((T1data[i])*0.693147181)-1;
            Resultdata[i]=(to_resize*(CurrSizes->rescale_max[0]-CurrSizes->rescale_min[0])+CurrSizes->rescale_min[0]);
        }
    }
    return Result;
}



nifti_image * Copy_Expec_to_Result_Neonate_mask(PrecisionTYPE * Expec,
                                                int * Short_2_Long_Indices,
                                                int * Long_2_Short_Indices,
                                                nifti_image * T1,
                                                float * BiasField,
                                                float * M,
                                                char * filename,
                                                ImageSize * CurrSizes){

    nifti_image * Result = nifti_copy_nim_info(T1);
    Result->dim[0]=4;
    Result->dim[4]=6;
    Result->dim[5]=1;
    Result->datatype=DT_FLOAT32;
    Result->cal_max=1;
    nifti_set_filenames(Result,filename,0,0);
    nifti_update_dims_from_array(Result);
    nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
    Result->data = (void *) calloc(Result->nvox, sizeof(PrecisionTYPE));
    PrecisionTYPE * Resultdata = static_cast<PrecisionTYPE *>(Result->data);
    PrecisionTYPE * T1ptrtmp = static_cast<PrecisionTYPE *>(T1->data);
    for(int i=0; i<Result->nvox; i++){Resultdata[i]=0;}

    int * Short_2_Long_Indices_PRT = (int *) Short_2_Long_Indices;

    int class_nvox=Result->nx*Result->ny*Result->nz;

    Short_2_Long_Indices_PRT= (int *) Short_2_Long_Indices;
    for(int currclass=0; currclass<6;currclass++){
        PrecisionTYPE * Resultdata_class = &Resultdata[(currclass)*class_nvox];
        PrecisionTYPE * Expec_PTR = &Expec[(currclass)*CurrSizes->numelmasked];
        Short_2_Long_Indices_PRT= (int *) Short_2_Long_Indices;
        for(int i=0; i<CurrSizes->numelmasked; i++,Short_2_Long_Indices_PRT++,Expec_PTR++){
            Resultdata_class[(*Short_2_Long_Indices_PRT)]=(*Expec_PTR);
        }
    }


    int xyzpos[3]={0};
    int distance=4;
    PrecisionTYPE * Resultdata_class = &Resultdata[(6)*class_nvox];
    PrecisionTYPE * Expec_PTR = &Expec[(6)*CurrSizes->numelmasked];
    Short_2_Long_Indices_PRT= (int *) Short_2_Long_Indices;
    float Fractional_content=0;
    int biggestclass=-1;
    float biggestclass_prob=-1;
    int biggestclass2=-1;
    float biggestclass2_prob=-1;
    int neigh_index=0;

    for(int i=0; i<CurrSizes->numelmasked; i++,Short_2_Long_Indices_PRT++){
        if(Expec[(6)*CurrSizes->numelmasked+i]>0.1){
            biggestclass=-1;
            biggestclass_prob=0.1;
            biggestclass2=-2;
            biggestclass2_prob=0.1;
            for(xyzpos[2]=-distance;xyzpos[2]<distance;xyzpos[2]++){
                for(xyzpos[1]=-distance;xyzpos[1]<distance;xyzpos[1]++){
                    for(xyzpos[0]=-distance;xyzpos[0]<distance;xyzpos[0]++){
                        neigh_index=Long_2_Short_Indices[(*Short_2_Long_Indices_PRT)+xyzpos[0]+xyzpos[1]*CurrSizes->xsize+xyzpos[2]*CurrSizes->xsize*CurrSizes->ysize];
                        if(neigh_index>0 && neigh_index<CurrSizes->numel){
                            for(int neigh_class=1; neigh_class<6;neigh_class++){
                                if(Expec[neigh_index+(neigh_class)*CurrSizes->numelmasked]>biggestclass_prob){
                                    if(biggestclass!=biggestclass2){
                                        biggestclass2_prob=biggestclass_prob;
                                        biggestclass2=biggestclass;
                                    }
                                    biggestclass_prob=Expec[neigh_index+(neigh_class)*CurrSizes->numelmasked];
                                    biggestclass=neigh_class;
                                }
                                if(Expec[neigh_index+(neigh_class)*CurrSizes->numelmasked]>biggestclass2_prob && Expec[neigh_index+(neigh_class)*CurrSizes->numelmasked]<biggestclass2_prob){
                                    if(biggestclass!=neigh_class){
                                        biggestclass2_prob=Expec[neigh_index+(neigh_class)*CurrSizes->numelmasked];
                                        biggestclass2=neigh_class;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if(biggestclass>0 && biggestclass2>0){
                Fractional_content=(M[biggestclass]-(T1ptrtmp[(*Short_2_Long_Indices_PRT)]+BiasField[i]))/(M[biggestclass]-M[biggestclass2]);
                Fractional_content=(Fractional_content<0)?0:Fractional_content;
                Fractional_content=(Fractional_content>1)?1:Fractional_content;
                Resultdata[(*Short_2_Long_Indices_PRT)+(biggestclass)*class_nvox]=(1-Fractional_content);
                Resultdata[(*Short_2_Long_Indices_PRT)+(biggestclass2)*class_nvox]=(Fractional_content);
                for(int otherclasses=0;otherclasses<6;otherclasses++){
                    if(otherclasses!=biggestclass && otherclasses!=biggestclass2){
                        Resultdata[(*Short_2_Long_Indices_PRT)+(otherclasses)*class_nvox]=0;
                    }
                }
            }
            else{
                for(int otherclasses=0;otherclasses<6;otherclasses++){
                    Resultdata[(*Short_2_Long_Indices_PRT)+(otherclasses)*class_nvox]=0;
                }
                Resultdata[(*Short_2_Long_Indices_PRT)+(2)*class_nvox]=1;
            }
        }
    }


    Short_2_Long_Indices_PRT= (int *) Short_2_Long_Indices;
    float sumexp=0;

    for(int i=0; i<CurrSizes->numel; i++){
        if(Long_2_Short_Indices[i]>0){
            sumexp=0;
            for(int currclass=0; currclass<6;currclass++){
                if((Resultdata[i+currclass*CurrSizes->numel])>0.01){
                    sumexp+=Resultdata[i+currclass*CurrSizes->numel];
                }
            }
            if(sumexp>0){
                for(int currclass=0; currclass<6;currclass++){
                    if((Resultdata[i+currclass*CurrSizes->numel])>0.01){
                        Resultdata[i+currclass*CurrSizes->numel]=Resultdata[i+currclass*CurrSizes->numel]/sumexp;
                    }
                    else{
                        Resultdata[i+currclass*CurrSizes->numel]=0;
                    }
                }
            }
        }
    }

    return Result;
}

nifti_image * Copy_Expec_to_Result_mask(PrecisionTYPE * Expec,
                                        int * Short_2_Long_Indices,
                                        nifti_image * T1,
                                        char * filename,
                                        ImageSize * CurrSizes){

    nifti_image * Result = nifti_copy_nim_info(T1);
    Result->dim[0]=4;
    Result->dim[4]=CurrSizes->numclass;
    Result->dim[5]=1;
    Result->datatype=DT_FLOAT32;
    Result->cal_max=1;
    nifti_set_filenames(Result,filename,0,0);
    nifti_update_dims_from_array(Result);
    nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
    Result->data = (void *) calloc(Result->nvox, sizeof(PrecisionTYPE));
    PrecisionTYPE * Resultdata = static_cast<PrecisionTYPE *>(Result->data);
    for(int i=0; i<Result->nvox; i++){Resultdata[i]=0;}

    int * Short_2_Long_Indices_PRT = (int *) Short_2_Long_Indices;

    int class_nvox=Result->nx*Result->ny*Result->nz;

    Short_2_Long_Indices_PRT= (int *) Short_2_Long_Indices;
    for(int currclass=0; currclass<CurrSizes->numclass;currclass++){

        PrecisionTYPE * Resultdata_class = &Resultdata[(currclass)*class_nvox];
        PrecisionTYPE * Expec_PTR = &Expec[(currclass)*CurrSizes->numelmasked];
        Short_2_Long_Indices_PRT= (int *) Short_2_Long_Indices;

        for(int i=0; i<CurrSizes->numelmasked; i++,Short_2_Long_Indices_PRT++,Expec_PTR++){
            Resultdata_class[Short_2_Long_Indices[i]]=*Expec_PTR;
        }
    }
    return Result;
}


nifti_image * Copy_Expec_to_Result(PrecisionTYPE * Expec,
                                   nifti_image * T1,
                                   char * filename,
                                   ImageSize * CurrSizes){

    nifti_image * Result = nifti_copy_nim_info(T1);
    Result->dim[0]=4;
    Result->dim[4]=CurrSizes->numclass;
    Result->datatype=DT_FLOAT32;
    Result->cal_max=1;
    nifti_set_filenames(Result,filename,0,0);
    nifti_update_dims_from_array(Result);
    nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
    Result->data = (void *) calloc(Result->nvox, sizeof(PrecisionTYPE));
    PrecisionTYPE * Resultdata = static_cast<PrecisionTYPE *>(Result->data);
    PrecisionTYPE * T1data = static_cast<PrecisionTYPE *>(T1->data);
    for(int i=0; i<Result->nvox; i++){Resultdata[i]=0;}
    int class_nvox=Result->nx*Result->ny*Result->nz;
    PrecisionTYPE to_resize=0;

    for(int currclass=0; currclass<CurrSizes->numclass;currclass++){

        PrecisionTYPE * Resultdata_class = &Resultdata[(currclass)*class_nvox];
        PrecisionTYPE * Expec_PTR = &Expec[(currclass)*CurrSizes->numel];

        for(int i=0; i<CurrSizes->numel; i++,Expec_PTR++){
            Resultdata_class[i]=*Expec_PTR;
        }
    }
    return Result;
}


nifti_image * Copy_single_image_to_Result(bool * Mask,
                                          nifti_image * Original,
                                          char * filename){

    nifti_image * Result = nifti_copy_nim_info(Original);
    Result->dim[0]=4;
    Result->dim[4]=1;

    Result->datatype=DT_INT32;
    Result->cal_max=1;
    nifti_set_filenames(Result,filename,0,0);
    nifti_update_dims_from_array(Result);

    nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);

    Result->data = (void *) calloc(Result->nvox, sizeof(int));

    int * Resultdata = static_cast<int *>(Result->data);

    bool * Expec_PTR = static_cast<bool *>(Mask);

    for(int i=0; i<Result->nvox; i++,Expec_PTR++,Resultdata++){
        *Resultdata=(int)(*Expec_PTR);
    }

    return Result;
}


nifti_image * Copy_single_image_to_Result(float * W,
                                          nifti_image * Original,
                                          char * filename){

    nifti_image * Result = nifti_copy_nim_info(Original);
    Result->dim[0]=4;
    Result->dim[4]=1;

    Result->datatype=DT_FLOAT32;
    Result->cal_max=1;
    nifti_set_filenames(Result,filename,0,0);
    nifti_update_dims_from_array(Result);

    nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);

    Result->data = (void *) calloc(Result->nvox, sizeof(float));

    float * Resultdata = static_cast<float *>(Result->data);

    float * Expec_PTR = static_cast<float *>(W);

    for(int i=0; i<Result->nvox; i++,Expec_PTR++,Resultdata++){
        *Resultdata=(float)(*Expec_PTR);
    }

    return Result;
}

nifti_image * Copy_single_image_to_Result(int * Image,
                                          nifti_image * Original,
                                          char * filename){

    nifti_image * Result = nifti_copy_nim_info(Original);
    Result->dim[0]=4;
    Result->dim[4]=1;

    Result->datatype=DT_INT32;

    Result->cal_max=1;
    nifti_set_filenames(Result,filename,0,0);
    nifti_update_dims_from_array(Result);

    nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);

    Result->data = (void *) calloc(Result->nvox, sizeof(int));

    int * Resultdata = static_cast<int *>(Result->data);
    int * Image_PTR = Image;

    for(int i=0; i<Result->nvox; i++,Image_PTR++,Resultdata++){
        *Resultdata=(int)(*Image_PTR);
    }

    return Result;
}


void quickSort(int *arr, int elements) {

    int  piv, beg[300], end[300], i=0, L, R, swap ;

    beg[0]=0; end[0]=elements;
    while (i>=0) {
        L=beg[i]; R=end[i]-1;
        if (L<R) {
            piv=arr[L];
            while (L<R) {
                while (arr[R]>=piv && L<R) R--; if (L<R) arr[L++]=arr[R];
                while (arr[L]<=piv && L<R) L++; if (L<R) arr[R--]=arr[L];
            }
            arr[L]=piv; beg[i+1]=L+1; end[i+1]=end[i]; end[i++]=L;
            if (end[i]-beg[i]>end[i-1]-beg[i-1]) {
                swap=beg[i];
                beg[i]=beg[i-1];
                beg[i-1]=swap;
                swap=end[i];
                end[i]=end[i-1];
                end[i-1]=swap;
            }
        }
        else {
            i--;
        }
    }
}



nifti_image * Get_Bias_Corrected(float * BiasField,
                                 nifti_image * T1,
                                 char * filename,
                                 ImageSize * CurrSizes){

    nifti_image * BiasCorrected = nifti_copy_nim_info(T1);
    BiasCorrected->dim[0]=4;
    BiasCorrected->dim[4]=CurrSizes->usize;
    BiasCorrected->datatype=DT_FLOAT32;
    BiasCorrected->cal_max=(CurrSizes->rescale_max[0]);

    nifti_set_filenames(BiasCorrected,filename,0,0);
    nifti_update_dims_from_array(BiasCorrected);
    nifti_datatype_sizes(BiasCorrected->datatype,&BiasCorrected->nbyper,&BiasCorrected->swapsize);
    BiasCorrected->data = (void *) calloc(BiasCorrected->nvox, sizeof(PrecisionTYPE));
    PrecisionTYPE * BiasCorrected_PTR = static_cast<PrecisionTYPE *>(BiasCorrected->data);
    PrecisionTYPE * T1data = static_cast<PrecisionTYPE *>(T1->data);
    for(int multispec=0;multispec<CurrSizes->usize;multispec++){

        BiasCorrected_PTR = static_cast<PrecisionTYPE *>(BiasCorrected->data);
        BiasCorrected_PTR = &BiasCorrected_PTR[multispec*BiasCorrected->nvox];

        T1data = static_cast<PrecisionTYPE *>(T1->data);
        T1data = &T1data[multispec*BiasCorrected->nvox];

        for(int i=0; i<BiasCorrected->nvox; i++){
            BiasCorrected_PTR[i]=0;
        }

        float to_resize=0;

        for(int i=0; i<CurrSizes->numel; i++){
            to_resize=exp((BiasField[i]+T1data[i])*0.693147181)-1;
            BiasCorrected_PTR[i]=(to_resize*(CurrSizes->rescale_max[multispec]-CurrSizes->rescale_min[multispec])+CurrSizes->rescale_min[multispec]);
        }
    }
    return BiasCorrected;
}


nifti_image * Get_Bias_Corrected_mask(float * BiasFieldCoefs,
                                      nifti_image * T1,
                                      char * filename,
                                      ImageSize * CurrSizes,
                                      int biasOrder){

    int UsedBasisFunctions=(int)((biasOrder+1) * (biasOrder+2)/2 *(biasOrder+3)/3);

    nifti_image * BiasCorrected = nifti_copy_nim_info(T1);
    BiasCorrected->dim[0]=4;
    BiasCorrected->dim[4]=CurrSizes->usize;
    BiasCorrected->datatype=DT_FLOAT32;
    BiasCorrected->cal_max=(CurrSizes->rescale_max[0]);

    nifti_set_filenames(BiasCorrected,filename,0,0);
    nifti_update_dims_from_array(BiasCorrected);
    nifti_datatype_sizes(BiasCorrected->datatype,&BiasCorrected->nbyper,&BiasCorrected->swapsize);
    BiasCorrected->data = (void *) calloc(BiasCorrected->nvox, sizeof(PrecisionTYPE));
    PrecisionTYPE * BiasCorrected_PTR = static_cast<PrecisionTYPE *>(BiasCorrected->data);
    PrecisionTYPE * T1data = static_cast<PrecisionTYPE *>(T1->data);
    float BiasField=0;
    PrecisionTYPE currxpower[maxallowedpowerorder];
    PrecisionTYPE currypower[maxallowedpowerorder];
    PrecisionTYPE currzpower[maxallowedpowerorder];
    float xpos=0.0f;
    float ypos=0.0f;
    float zpos=0.0f;
    PrecisionTYPE not_point_five_times_dims_x=(0.5f*(PrecisionTYPE)CurrSizes->xsize);
    PrecisionTYPE not_point_five_times_dims_y=(0.5f*(PrecisionTYPE)CurrSizes->ysize);
    PrecisionTYPE not_point_five_times_dims_z=(0.5f*(PrecisionTYPE)CurrSizes->zsize);
    PrecisionTYPE inv_not_point_five_times_dims_x=1.0f/(0.5f*(PrecisionTYPE)CurrSizes->xsize);
    PrecisionTYPE inv_not_point_five_times_dims_y=1.0f/(0.5f*(PrecisionTYPE)CurrSizes->ysize);
    PrecisionTYPE inv_not_point_five_times_dims_z=1.0f/(0.5f*(PrecisionTYPE)CurrSizes->zsize);
    int ind=0;

    for(int multispec=0;multispec<CurrSizes->usize;multispec++){

        BiasCorrected_PTR = static_cast<PrecisionTYPE *>(BiasCorrected->data);
        BiasCorrected_PTR = &BiasCorrected_PTR[multispec*CurrSizes->numel];
        T1data = static_cast<PrecisionTYPE *>(T1->data);
        T1data = &T1data[multispec*CurrSizes->numel];

        float * BiasFieldCoefs_multispec = &BiasFieldCoefs[multispec*UsedBasisFunctions];


        for(int i=0; i<BiasCorrected->nvox; i++){
            BiasCorrected_PTR[i]=0;
        }

        float to_resize=0;
        int index_full=0;
        for (int iz=0; iz<CurrSizes->zsize; iz++) {
            for (int iy=0; iy<CurrSizes->ysize; iy++) {
                for (int ix=0; ix<CurrSizes->xsize; ix++) {
                    BiasField=0.0f;
                    xpos=(((PrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                    ypos=(((PrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                    zpos=(((PrecisionTYPE)iz-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z);
                    get_xyz_pow_int(xpos, ypos, zpos, currxpower, currypower, currzpower, biasOrder);
                    ind=0;
                    for(int order=0; order<=biasOrder; order++){
                        for(int xorder=0; xorder<=order; xorder++){
                            for(int yorder=0; yorder<=(order-xorder); yorder++){
                                int zorder=order-yorder-xorder;
                                BiasField-=BiasFieldCoefs_multispec[ind]*currxpower[xorder]*currypower[yorder]*currzpower[zorder];
                                ind++;
                            }
                        }
                    }

                    to_resize=exp((BiasField+T1data[index_full])*0.693147181)-1;
                    BiasCorrected_PTR[index_full]=(to_resize*(CurrSizes->rescale_max[multispec]-CurrSizes->rescale_min[multispec])+CurrSizes->rescale_min[multispec]);
                    index_full++;
                }
            }
        }
    }
    return BiasCorrected;
}

float * seg_norm4NCC(nifti_image * BaseImage, nifti_image * LNCC, nifti_image * Lables,int distance,ImageSize * CurrSizes){

    PrecisionTYPE * LNCCptr = static_cast<float *>(LNCC->data);
    PrecisionTYPE * BaseImageptr = static_cast<float *>(BaseImage->data);
    bool * Lablesptr = static_cast<bool *>(Lables->data);

    int distance2=ceil((float)(distance));
    int distance3=ceil((float)(distance));
    int directionShift[3];
    int dim_array[3];
    dim_array[0]=(int)BaseImage->nx;
    dim_array[1]=(int)BaseImage->ny;
    dim_array[2]=(int)BaseImage->nz;
    directionShift[0]=1;
    directionShift[1]=BaseImage->nx;
    directionShift[2]=BaseImage->nx*BaseImage->ny;
    float * BaseMean=new float [BaseImage->nx*BaseImage->ny*BaseImage->nz];
    float * BaseSTD=new float [BaseImage->nx*BaseImage->ny*BaseImage->nz];


    float * bufferMean=new float [BaseImage->nx*BaseImage->ny*BaseImage->nz];
    float * bufferSTD=new float [BaseImage->nx*BaseImage->ny*BaseImage->nz];
    float * bufferDATA=new float [BaseImage->nx*BaseImage->ny*BaseImage->nz];

    float tmpsum=0;
    float tmpsumshift=0;

    // CALC MEAN AND STD OF THE BASE
    cout << "Calculating LNCC"<<endl;
    cout << "Local Mean and STD of the base image"<<endl;
    flush(cout);
    for(int i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz;i++){
        BaseMean[i]=0;
        BaseSTD[i]=0;
    }
    // Calc Mean

    int xyzpos[3];
    int index=0;
    for(xyzpos[2]=0;xyzpos[2]<(int)BaseImage->nz;xyzpos[2]++){
        for(xyzpos[1]=0;xyzpos[1]<(int)BaseImage->ny;xyzpos[1]++){
            for(xyzpos[0]=0;xyzpos[0]<(int)BaseImage->nx;xyzpos[0]++){
                tmpsum=0.0f;
                tmpsumshift=0;
                int xyzShift[3];
                for(xyzShift[2]=((xyzpos[2]<distance)?-xyzpos[2]:-distance);xyzShift[2]<=((xyzpos[2]>=(dim_array[2]-distance))?(int)dim_array[2]-xyzpos[2]-1:distance);xyzShift[2]++){
                    for(xyzShift[1]=((xyzpos[1]<distance)?-xyzpos[1]:-distance);xyzShift[1]<=((xyzpos[1]>=(dim_array[1]-distance))?(int)dim_array[1]-xyzpos[1]-1:distance);xyzShift[1]++){
                        for(xyzShift[0]=((xyzpos[0]<distance)?-xyzpos[0]:-distance);xyzShift[0]<=((xyzpos[0]>=(dim_array[0]-distance))?(int)dim_array[0]-xyzpos[0]-1:distance);xyzShift[0]++){
                            tmpsum+=(BaseImageptr[index+xyzShift[0]*directionShift[0]+xyzShift[1]*directionShift[1]+xyzShift[2]*directionShift[2]]);
                            tmpsumshift+=1;
                        }
                    }
                }
                BaseMean[index]=tmpsum/tmpsumshift;
                index++;
            }
        }
    }



    // Calc STD
    float tmpmean;
    float tmpsumstep;
    index=0;
    for(xyzpos[2]=0;xyzpos[2]<(int)BaseImage->nz;xyzpos[2]++){
        for(xyzpos[1]=0;xyzpos[1]<(int)BaseImage->ny;xyzpos[1]++){
            for(xyzpos[0]=0;xyzpos[0]<(int)BaseImage->nx;xyzpos[0]++){
                tmpsum=0.0f;
                tmpsumstep=0;
                tmpmean=BaseMean[index];
                tmpsumshift=0;
                int xyzShift[3];
                for(xyzShift[2]=((xyzpos[2]<distance2)?-xyzpos[2]:-distance2);xyzShift[2]<=((xyzpos[2]>=(dim_array[2]-distance2))?(int)dim_array[2]-xyzpos[2]-1:distance2);xyzShift[2]++){
                    for(xyzShift[1]=((xyzpos[1]<distance2)?-xyzpos[1]:-distance2);xyzShift[1]<=((xyzpos[1]>=(dim_array[1]-distance2))?(int)dim_array[1]-xyzpos[1]-1:distance2);xyzShift[1]++){
                        for(xyzShift[0]=((xyzpos[0]<distance2)?-xyzpos[0]:-distance2);xyzShift[0]<=((xyzpos[0]>=(dim_array[0]-distance2))?(int)dim_array[0]-xyzpos[0]-1:distance2);xyzShift[0]++){
                            tmpsumstep=(BaseImageptr[index+xyzShift[0]*directionShift[0]+xyzShift[1]*directionShift[1]+xyzShift[2]*directionShift[2]]-tmpmean);
                            tmpsum+=(tmpsumstep*tmpsumstep);
                            tmpsumshift+=1;
                        }
                    }
                }
                BaseSTD[index]=sqrt(tmpsum/tmpsumshift);
                index++;
            }
        }
    }



    // CALC LNCC FOR EACH LABLE
    cout << "Local Mean and STD of the Template images"<<endl;
    flush(cout);
    for(int currlable=0;currlable<LNCC->nt; currlable++){
        //for(int currlable=0;currlable<3; currlable++){

        cout << currlable+1 << "/" << LNCC->nt<< endl;
        flush(cout);
        float * currLNCCptr=&LNCCptr[currlable*BaseImage->nx*BaseImage->ny*BaseImage->nz];
        for(int i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz;i++){
            bufferDATA[i]=currLNCCptr[i];
        }
        // Calc Mean

        int xyzpos[3];
        int index=0;
        for(xyzpos[2]=0;xyzpos[2]<(int)BaseImage->nz;xyzpos[2]++){
            for(xyzpos[1]=0;xyzpos[1]<(int)BaseImage->ny;xyzpos[1]++){
                for(xyzpos[0]=0;xyzpos[0]<(int)BaseImage->nx;xyzpos[0]++){
                    tmpsum=0.0f;
                    tmpsumshift=0;
                    int xyzShift[3];
                    for(xyzShift[2]=((xyzpos[2]<distance)?-xyzpos[2]:-distance);xyzShift[2]<=((xyzpos[2]>=(dim_array[2]-distance))?(int)dim_array[2]-xyzpos[2]-1:distance);xyzShift[2]++){
                        for(xyzShift[1]=((xyzpos[1]<distance)?-xyzpos[1]:-distance);xyzShift[1]<=((xyzpos[1]>=(dim_array[1]-distance))?(int)dim_array[1]-xyzpos[1]-1:distance);xyzShift[1]++){
                            for(xyzShift[0]=((xyzpos[0]<distance)?-xyzpos[0]:-distance);xyzShift[0]<=((xyzpos[0]>=(dim_array[0]-distance))?(int)dim_array[0]-xyzpos[0]-1:distance);xyzShift[0]++){
                                tmpsum+=(bufferDATA[index+xyzShift[0]*directionShift[0]+xyzShift[1]*directionShift[1]+xyzShift[2]*directionShift[2]]);
                                tmpsumshift+=1;
                            }
                        }
                    }
                    bufferMean[index]=tmpsum/tmpsumshift;
                    index++;
                }
            }
        }



        // Calc STD
        float tmpmean;
        float tmpsumstep;
        index=0;
        for(xyzpos[2]=0;xyzpos[2]<(int)BaseImage->nz;xyzpos[2]++){
            for(xyzpos[1]=0;xyzpos[1]<(int)BaseImage->ny;xyzpos[1]++){
                for(xyzpos[0]=0;xyzpos[0]<(int)BaseImage->nx;xyzpos[0]++){
                    tmpsum=0.0f;
                    tmpsumstep=0;
                    tmpmean=bufferMean[index];
                    tmpsumshift=0;
                    int xyzShift[3];
                    for(xyzShift[2]=((xyzpos[2]<distance)?-xyzpos[2]:-distance);xyzShift[2]<=((xyzpos[2]>=(dim_array[2]-distance))?(int)dim_array[2]-xyzpos[2]-1:distance);xyzShift[2]++){
                        for(xyzShift[1]=((xyzpos[1]<distance)?-xyzpos[1]:-distance);xyzShift[1]<=((xyzpos[1]>=(dim_array[1]-distance))?(int)dim_array[1]-xyzpos[1]-1:distance);xyzShift[1]++){
                            for(xyzShift[0]=((xyzpos[0]<distance)?-xyzpos[0]:-distance);xyzShift[0]<=((xyzpos[0]>=(dim_array[0]-distance))?(int)dim_array[0]-xyzpos[0]-1:distance);xyzShift[0]++){
                                tmpsumstep=(bufferDATA[index+xyzShift[0]*directionShift[0]+xyzShift[1]*directionShift[1]+xyzShift[2]*directionShift[2]]-tmpmean);
                                tmpsum+=(tmpsumstep*tmpsumstep);
                                tmpsumshift+=1;
                            }
                        }
                    }
                    bufferSTD[index]=sqrt(tmpsum/tmpsumshift);
                    index++;
                }
            }
        }
        index=0;
        float ncctmp=0;
        for(xyzpos[2]=0;xyzpos[2]<(int)BaseImage->nz;xyzpos[2]++){
            for(xyzpos[1]=0;xyzpos[1]<(int)BaseImage->ny;xyzpos[1]++){
                for(xyzpos[0]=0;xyzpos[0]<(int)BaseImage->nx;xyzpos[0]++){
                    ncctmp=0;
                    tmpsumshift=0;
                    int xyzShift[3];
                    for(xyzShift[2]=((xyzpos[2]<distance2)?-xyzpos[2]:-distance2);xyzShift[2]<=((xyzpos[2]>=(dim_array[2]-distance2))?(int)dim_array[2]-xyzpos[2]-1:distance2);xyzShift[2]++){
                        for(xyzShift[1]=((xyzpos[1]<distance2)?-xyzpos[1]:-distance2);xyzShift[1]<=((xyzpos[1]>=(dim_array[1]-distance2))?(int)dim_array[1]-xyzpos[1]-1:distance2);xyzShift[1]++){
                            for(xyzShift[0]=((xyzpos[0]<distance2)?-xyzpos[0]:-distance2);xyzShift[0]<=((xyzpos[0]>=(dim_array[0]-distance2))?(int)dim_array[0]-xyzpos[0]-1:distance2);xyzShift[0]++){
                                int currindex=index+xyzShift[0]*directionShift[0]+xyzShift[1]*directionShift[1]+xyzShift[2]*directionShift[2];
                                ncctmp+=((bufferDATA[currindex]-bufferMean[index])*(BaseImageptr[currindex]-BaseMean[index]))/(bufferSTD[index]*BaseSTD[index]);
                                tmpsumshift+=1;
                            }
                        }
                    }
                    currLNCCptr[index]=(ncctmp>0?ncctmp:0)/(tmpsumshift+1);
                    currLNCCptr[index]=(currLNCCptr[index]>1?1:currLNCCptr[index]);
                    index++;
                }
            }
        }
    }

    cout << "Filtering LNCC"<< endl;
    Gaussian_Filter_4D(LNCCptr,(float)(distance3),CurrSizes);


    delete [] bufferSTD;
    delete [] bufferMean;
    delete [] BaseSTD;
    delete [] BaseMean;
    delete [] bufferDATA;

    float currSUM;


    cout << "Normalizing LNCC"<< endl;
    for(int i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz;i++){
        currSUM=0;
        for(int currlable=0;currlable<LNCC->nt; currlable++){
            currSUM+=LNCCptr[i+currlable*BaseImage->nx*BaseImage->ny*BaseImage->nz];
        }
        for(int currlable=0;currlable<LNCC->nt; currlable++){
            if(currSUM>0.01){
                LNCCptr[i+currlable*BaseImage->nx*BaseImage->ny*BaseImage->nz]=(float)(LNCCptr[i+currlable*BaseImage->nx*BaseImage->ny*BaseImage->nz])/currSUM;
            }
            else{
                LNCCptr[i+currlable*BaseImage->nx*BaseImage->ny*BaseImage->nz]=(float)1.0f/(float)LNCC->nt;
            }
        }
    }


    return LNCCptr;
}

/* *************************************************************** */
template <class NewTYPE, class DTYPE>
void seg_changeDatatype1(nifti_image *image)
{
        // the initial array is saved and freeed
        DTYPE *initialValue = (DTYPE *)malloc(image->nvox*sizeof(DTYPE));
        memcpy(initialValue, image->data, image->nvox*sizeof(DTYPE));

        // the new array is allocated and then filled
        if(sizeof(NewTYPE)==sizeof(unsigned char)) image->datatype = NIFTI_TYPE_UINT8;
        else if(sizeof(NewTYPE)==sizeof(float)) image->datatype = NIFTI_TYPE_FLOAT32;
        else if(sizeof(NewTYPE)==sizeof(double)) image->datatype = NIFTI_TYPE_FLOAT64;
        else{
        fprintf(stderr,"[NiftyReg ERROR] seg_changeDatatype\tOnly change to unsigned char, float or double are supported\n");
        exit(1);
        }
        free(image->data);
        image->nbyper = sizeof(NewTYPE);
        image->data = (void *)calloc(image->nvox,sizeof(NewTYPE));
        NewTYPE *dataPtr = static_cast<NewTYPE *>(image->data);
        for(unsigned int i=0; i<image->nvox; i++) dataPtr[i] = (NewTYPE)(initialValue[i]);

        free(initialValue);
        return;
}
/* *************************************************************** */
template <class NewTYPE>
void seg_changeDatatype(nifti_image *image)
{
        switch(image->datatype){
                case NIFTI_TYPE_UINT8:
                        seg_changeDatatype1<NewTYPE,unsigned char>(image);
                        break;
                case NIFTI_TYPE_INT8:
                        seg_changeDatatype1<NewTYPE,char>(image);
                        break;
                case NIFTI_TYPE_UINT16:
                        seg_changeDatatype1<NewTYPE,unsigned short>(image);
                        break;
                case NIFTI_TYPE_INT16:
                        seg_changeDatatype1<NewTYPE,short>(image);
                        break;
                case NIFTI_TYPE_UINT32:
                        seg_changeDatatype1<NewTYPE,unsigned int>(image);
                        break;
                case NIFTI_TYPE_INT32:
                        seg_changeDatatype1<NewTYPE,int>(image);
                        break;
                case NIFTI_TYPE_FLOAT32:
                        seg_changeDatatype1<NewTYPE,float>(image);
                        break;
                case NIFTI_TYPE_FLOAT64:
                        seg_changeDatatype1<NewTYPE,double>(image);
                        break;
                default:
            fprintf(stderr,"[NiftyReg ERROR] seg_changeDatatype\tThe initial image data type is not supported\n");
            exit(1);
        }
}
/* *************************************************************** */
template void seg_changeDatatype<unsigned char>(nifti_image *);
template void seg_changeDatatype<float>(nifti_image *);
template void seg_changeDatatype<double>(nifti_image *);
/* *************************************************************** */

