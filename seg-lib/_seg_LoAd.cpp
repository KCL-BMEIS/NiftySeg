#include "_seg_LoAd.h"

/// @brief The LoAd segment tool
/// \deprecated
nifti_image * LoAd_Segment(nifti_image * T1, nifti_image * Mask, nifti_image * Priors, seg_EM_Params * segment_param)
{

    time_t start,end;
    time(&start);
    if((int)(segment_param->verbose_level)>(int)(0))
    {
        int verboselevel=segment_param->verbose_level;
        cout << "LoAd: Verbose level " << verboselevel << endl;
    }
    // TEMPORARY int NumClass=Priors->nt;
    ImageSize * CurrSizes = new ImageSize [1]();
    CurrSizes->numel=(int)(T1->nx*T1->ny*T1->nz);
    CurrSizes->xsize=T1->nx;
    CurrSizes->ysize=T1->ny;
    CurrSizes->zsize=T1->nz;
    CurrSizes->usize=1;
    CurrSizes->tsize=1;
    CurrSizes->numclass=5;
    CurrSizes->numelmasked=0;
    CurrSizes->numelbias=0;

    seg_EM_Flags * flags = new seg_EM_Flags [1]();
    flags->improv_phase=0;
    flags->out=1;
    flags->oldloglik=(segPrecisionTYPE)1.0;
    flags->loglik=(segPrecisionTYPE)2.0;
    flags->prior_relax=!((segment_param->relax_factor)>0);
    flags->do_pv_modeling=segment_param->flag_PV_model;
    flags->pv_modeling_on=0;
    flags->sg_delineation=!(segment_param->flag_SG_deli);

    // Creating G and H matrixes
    segPrecisionTYPE be=segment_param->MRF_strength*1.0f;
    segPrecisionTYPE ba=segment_param->MRF_strength*5.0f;
    segPrecisionTYPE ratioGH=1.0f;

    // Alocate more than enough space for G and H even if the PV option is active
    segPrecisionTYPE G[maxNumbClass*maxNumbClass]= {0.0};
    segPrecisionTYPE H[maxNumbClass*maxNumbClass]= {0.0};
    Create_GH_5class(G,H,ba,be,ratioGH,segment_param);

    // Sanity check... Normalizing Priors and removing NaN
    if(segment_param->verbose_level>0)
    {
        cout << "Normalizing Data" << endl;
    }
    if(Mask->datatype!=DT_BINARY)
    {
        seg_convert2binary(Mask,0.0f);
    }
    Normalize_NaN_Priors_mask(Priors,Mask,segment_param->verbose_level);
    Normalize_Image_mask(T1,Mask,CurrSizes,segment_param->verbose_level);
    int * Short_2_Long_Indices = Create_Short_2_Long_Matrix_from_NII(Mask,&(CurrSizes->numelmasked));
    int * Long_2_Short_Indices = Create_Long_2_Short_Matrix_from_NII(Mask);

    // if PV modeling is on, prealocate extra memory for the 2 PV classes

    segPrecisionTYPE * Expec = Create_cArray_from_Prior_mask(Mask,Priors,CurrSizes->numclass,segment_param->flag_PV_model);
    segPrecisionTYPE * ShortPrior = Create_cArray_from_Prior_mask(Mask,Priors,CurrSizes->numclass,segment_param->flag_PV_model);
    segPrecisionTYPE * MRF = new segPrecisionTYPE[CurrSizes->numelmasked*(CurrSizes->numclass+(int)(segment_param->flag_PV_model*2))]();
    for(int i=0; i<(CurrSizes->numelmasked*(CurrSizes->numclass+(int)(segment_param->flag_PV_model*2))); i++)
    {
        MRF[i]=(segPrecisionTYPE)(1.0);
    }
    segPrecisionTYPE * BiasField = new segPrecisionTYPE[CurrSizes->numelmasked]();
    segPrecisionTYPE * BiasFieldCoefs = new segPrecisionTYPE[((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)]();
    for(int i=0; i<CurrSizes->numelmasked; i++)
    {
        BiasField[i]=0;
    }
    segPrecisionTYPE M [maxNumbClass]= {0.0f};
    segPrecisionTYPE V [maxNumbClass]= {0.0f};
    segPrecisionTYPE * MRF_Beta =NULL;


    // Initialize
    calcM_mask(T1,Expec,BiasField,NULL,Short_2_Long_Indices,M, V,NULL,NULL,1.0,CurrSizes,segment_param->verbose_level);
    bool BFupdate=!(segment_param->flag_Bias);
    bool MRFupdate=!(segment_param->flag_MRF);
    int CurrentBF=segment_param->bias_order<=1?segment_param->bias_order:1;
    //**************
    // EM Algorithm
    //**************
    CurrSizes->numclass=5;
    int iter=(int)0;
    int iterchange=0;
    float last_SGcorrect_loglik=1;
    while (flags->out && iter < segment_param->maxIteration)
    {
        iter++;
        iterchange++;
        if(segment_param->verbose_level>0)
        {
            cout << "\n*******************************" << endl;
            cout << "Iteration " << iter << endl;
        }

        // Iterative Components - EM, MRF, Bias Correction
        //Expectation and LogLik
        calcE_mask(T1,MRF,Expec,&flags->loglik,BiasField,NULL,0,Short_2_Long_Indices,M,V,CurrSizes,segment_param->verbose_level);
        //Bias Correction
        if(BFupdate)BiasCorrection_mask(BiasField,BiasFieldCoefs,T1,Long_2_Short_Indices,Expec,NULL,M,V,CurrentBF,CurrSizes,segment_param->flag_Bias,segment_param->verbose_level);
        //MRF
        if(MRFupdate)MRFregularization_mask(Expec,G,H,MRF_Beta,MRF,ShortPrior,Long_2_Short_Indices,Short_2_Long_Indices,CurrSizes,segment_param->flag_MRF,segment_param->verbose_level);
        // Print LogLik depending on the verbose level
        if(segment_param->verbose_level>0)printloglik(iter,flags->loglik,flags->oldloglik);
        //Maximization
        calcM_mask_LoAd(T1,Expec,BiasField,Short_2_Long_Indices,M,V,CurrSizes,segment_param->verbose_level,flags->pv_modeling_on);

        // Preform Segmentation Refinement Steps or Exit
        if((((flags->loglik-flags->oldloglik)/fabs(flags->oldloglik))<(segPrecisionTYPE)(0.005) && iterchange>3) || iter>segment_param->maxIteration)
        {
            switch(flags->improv_phase)
            {
            case 0: // Activate BC
                if(BFupdate==0 || CurrentBF<segment_param->bias_order)
                {
                    flags->loglik=1;
                    CurrentBF++;
                    BFupdate=true;
                    if(CurrentBF>=segment_param->bias_order)flags->improv_phase++;
                    break;
                }
                else
                {
                    flags->improv_phase++;
                }
            case 1: // Activate MRF
                if(MRFupdate==0)
                {
                    flags->improv_phase++;
                    MRFupdate=true;
                    flags->loglik=1;
                    iterchange=0;
                    break;
                }
                else
                {
                    flags->improv_phase++;
                }
            case 2: // Run Prior Relaxation
                if((segment_param->relax_factor)>0 && flags->prior_relax==false)
                {
                    flags->improv_phase++;
                    flags->prior_relax=true;
                    Relax_Priors(ShortPrior,Expec,MRF,Short_2_Long_Indices,Long_2_Short_Indices,segment_param->relax_factor,G,ba,be,CurrSizes,segment_param);
                    flags->loglik=1;
                    iterchange=0;
                    break;
                }
                else
                {
                    flags->improv_phase++;
                }
            case 3: // Run Implicit PV modeling
                if(flags->do_pv_modeling==true)
                {
                    flags->improv_phase++;
                    flags->pv_modeling_on=true;
                    Convert_to_PV(T1,BiasField,ShortPrior,Expec,MRF,M,V,Short_2_Long_Indices,Long_2_Short_Indices,CurrSizes,segment_param);
                    Create_GH_7class(G,H,ba,be,ratioGH,segment_param);
                    CurrSizes->numclass=7;
                    //segment_param->flag_Bias=false;
                    flags->loglik=1;
                    iterchange=0;
                    break;
                }
                else
                {
                    flags->improv_phase++;
                }
            case 4: // Run Sulci and Gyri deliniation
                if(flags->sg_delineation==false && flags->do_pv_modeling==true)
                {
                    if(MRF_Beta==NULL)
                    {
                        MRF_Beta=new segPrecisionTYPE [CurrSizes->numelmasked]();
                    }

                    if((flags->loglik-last_SGcorrect_loglik)/fabs(last_SGcorrect_loglik)<(segPrecisionTYPE)(0.05) )
                    {
                        flags->sg_delineation=true;
                        flags->out=false;
                    }
                    else
                    {
                        Sulci_and_gyri_correction(MRF_Beta,ShortPrior,Expec,MRF,Short_2_Long_Indices,Long_2_Short_Indices,CurrSizes);
                        last_SGcorrect_loglik=flags->loglik;
                        flags->loglik=1;
                        iterchange=0;
                    }

                    break;
                }
                else
                {
                    flags->improv_phase++;
                }
            case 5: // Exit EM
                flags->out=false;
                break;
            }

        }
        // Update LogLik
        flags->oldloglik=flags->loglik;
    }


    // ******************* OUT IMAGE
    nifti_image * Result=NULL;
    //temporary save single image
    if(true)
    {
        if(flags->do_pv_modeling)
        {
            Result = Copy_ShortExpec_to_Result(T1,Expec,BiasField,BiasFieldCoefs,Short_2_Long_Indices,Priors,segment_param,M,CurrSizes);
        }
        else
        {
            CurrSizes->numclass=5;
            Result = Copy_Expec_to_Result_mask(Expec,Short_2_Long_Indices,T1,segment_param->filename_out,CurrSizes);
        }
    }
    else
    {
        Result = Copy_Single_ShortImage_to_Result(MRF_Beta,Short_2_Long_Indices,Priors,segment_param->filename_out,CurrSizes);
    }

    delete [] Short_2_Long_Indices;
    delete [] Expec;
    delete [] MRF;
    delete [] BiasField;
    delete [] Long_2_Short_Indices;
    delete [] ShortPrior;
    delete [] BiasFieldCoefs;
    delete [] CurrSizes;
    delete [] flags;

    time(&end);

    if(segment_param->verbose_level>0)
    {
        cout << "Finished in "<<difftime(end,start)<<"sec"<< endl;
    }
    return Result;
}




int Create_diagonal_GH_Nclass(segPrecisionTYPE * G,
                              segPrecisionTYPE * H,
                              segPrecisionTYPE ratio,
                              seg_EM_Params * segment_param)
{
    //WM - neighbours GM and dGM
    int numclass=segment_param->numb_classes;
    for(long i=0; i<segment_param->numb_classes; i++)
    {
        for(long j=0; j<segment_param->numb_classes; j++)
        {
            if(i!=j)
            {
                G[i+j*numclass]=segment_param->MRF_strength;
                H[i+j*numclass]=ratio*G[i+j*numclass];
            }
            else
            {
                G[i+j*numclass]=0.0;
                H[i+j*numclass]=ratio*G[i+j*numclass];
            }

        }
    }

    //Print G
    if(segment_param->verbose_level>1)
    {
        cout<<"G=" << endl;
        for (int i=0; i<numclass; i++)
        {
            for (int  j=0; j<numclass; j++)
            {
                cout<< (segPrecisionTYPE)G[i+j*numclass] << '\t';
            }
            cout<< endl;
        }
        cout<< endl << "H=" << endl;
        for (int i=0; i<numclass; i++)
        {
            for (int j=0; j<numclass; j++)
            {
                cout<< G[i+j*numclass] << '\t';
            }
            cout<< endl;
        }
        cout<< endl;
    }

    return 1;
}


int Create_GH_5class(segPrecisionTYPE * G,
                     segPrecisionTYPE * H,
                     segPrecisionTYPE ba,
                     segPrecisionTYPE be,
                     segPrecisionTYPE ratio,
                     seg_EM_Params * segment_param)
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
    for (j=0; j<numclass; j++)
    {
        for (i=1+j; i<numclass; i++)
        {
            G[j+i*numclass]=G[i+j*numclass];
        }
    }
    //Copy G to H
    for (int i=0; i<(numclass*numclass); i++)
    {
        H[i]=ratio*G[i];
    }
    //Print G
    if(segment_param->verbose_level>1)
    {
        cout<<"G=" << endl;
        for (i=0; i<numclass; i++)
        {
            for (j=0; j<numclass; j++)
            {
                cout<< (segPrecisionTYPE)G[i+j*numclass] << '\t';
            }
            cout<< endl;
        }
        cout<< endl << "H=" << endl;
        for (i=0; i<numclass; i++)
        {
            for (j=0; j<numclass; j++)
            {
                cout<< G[i+j*numclass] << '\t';
            }
            cout<< endl;
        }
        cout<< endl;
    }

    return 1;
}


int Create_GH_7class(segPrecisionTYPE * G,
                     segPrecisionTYPE * H,
                     segPrecisionTYPE ba,
                     segPrecisionTYPE be,
                     segPrecisionTYPE ratio,
                     seg_EM_Params * segment_param)
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
    for (j=0; j<numclass; j++)
    {
        for (i=1+j; i<numclass; i++)
        {
            G[j+i*numclass]=G[i+j*numclass];
        }
    }
    //Copy G to H
    for (int i=0; i<(numclass*numclass); i++)
    {
        H[i]=ratio*G[i];
    }
    //Print G
    if(segment_param->verbose_level>1)
    {
        cout<<"G=" << endl;
        for (i=0; i<numclass; i++)
        {
            for (j=0; j<numclass; j++)
            {
                cout<< G[i+j*numclass] << '\t';
            }
            cout<< endl;
        }
        cout<< endl << "H=" << endl;
        for (i=0; i<numclass; i++)
        {
            for (j=0; j<numclass; j++)
            {
                cout<< G[i+j*numclass] << '\t';
            }
            cout<< endl;
        }
        cout<< endl;
    }


    return 1;
}

int Normalize_NaN_Priors(nifti_image * Priors,
                         bool verbose)
{
    register int numel = Priors->nx*Priors->ny*Priors->nz;
    register int ups=0;
    register int good=0;
    if(verbose>0)
    {
        cout<< "Normalizing Priors" << endl;
    }
    if(Priors->datatype==NIFTI_TYPE_FLOAT32)
    {
        segPrecisionTYPE * priorsptr = static_cast<segPrecisionTYPE *>(Priors->data);
        for (int i=0; i<numel; i++)
        {
            float tempsum=0;
            for (int j=0; j<Priors->nt; j++)
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
                    priorsptr[tempind]=0.2f;
                }
                ups++;
            }
        }
    }
    else
    {

        printf("err\tNormalize_NaN_Priors\tWrong Image datatype\n");

    }

    if(verbose>0)
    {
        cout<<"Priors: "<< good<<" good voxels and "<<ups<<" bad voxels" << endl;
        flush(cout);
    }

    return 1;
}

int PriorWeight_mask(float * ShortPrior,nifti_image * Priors, float * Expec,float GaussKernelSize,float RelaxFactor, int * S2L, int * L2S,ImageSize * CurrSizes,int verbose_level)
{

    for(long i=0; i<(CurrSizes->numclass*CurrSizes->numelmasked); i++)ShortPrior[i]=Expec[i];
    GaussianFilter4D_cArray(ShortPrior,S2L,L2S,GaussKernelSize,CurrSizes);
    float * PriorsPtr = static_cast<float *>(Priors->data);
    long currindex=0;
    for(long k=0; k<CurrSizes->numclass; k++)
    {
        currindex=k*CurrSizes->numelmasked;
        for(long i=0; i<(CurrSizes->numelmasked); i++)
        {
            ShortPrior[currindex]*=(1-RelaxFactor);
            ShortPrior[currindex]+=(RelaxFactor)*PriorsPtr[S2L[i]+k*CurrSizes->numel];
            currindex++;
        }
    }

    return 1;
}




int Normalize_NaN_Priors_mask(nifti_image * Priors,
                              nifti_image * Mask,
                              bool verbose)
{
    register int numel = Mask->nx*Mask->ny*Mask->nz;
    register int ups=0;
    register int good=0;
    if(verbose>0)
    {
        cout<< "Normalizing Priors" << endl;
    }
    if(Mask->datatype==DT_BINARY)
    {
        if(Priors->datatype==NIFTI_TYPE_FLOAT32)
        {
            segPrecisionTYPE * priorsptr = static_cast<segPrecisionTYPE *>(Priors->data);
            bool * brainmaskptr = static_cast<bool *> (Mask->data);

            for (int i=0; i<numel; i++)
            {
                if(brainmaskptr[i])
                {
                    float tempsum=0;
                    for (int j=0; j<Priors->nt; j++)
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

            printf("err\tNormalize_NaN_Priors\tWrong Image datatype\n");

        }
    }
    else
    {

        printf("err\tNormalize_NaN_Priors\tWrong mask datatype\n");

    }

    if(verbose>0)
    {
        cout<<"Priors: "<< good<<" good voxels and "<<ups<<" bad voxels" << endl;
        flush(cout);
    }

    return 1;
}


int Normalize_Image_mask(nifti_image * input,
                         nifti_image * Mask,
                         ImageSize * CurrSizes,
                         bool verbose)
{
    if(input->datatype!=NIFTI_TYPE_FLOAT32)
    {
        seg_changeDatatype<float>(input);
    }
    if(verbose>0)
    {
        cout<< "Normalizing Input Image" << endl;
    }
    int numel=(int)(input->nx*input->ny*input->nz);
    if(Mask->datatype!=DT_BINARY)
    {
        seg_convert2binary(Mask,0.0f);
    }
    if(input->datatype!=NIFTI_TYPE_FLOAT32)
    {
        seg_changeDatatype<segPrecisionTYPE>(input);
    }
    for(long udir=0; udir<CurrSizes->usize; udir++) // Per Multispectral Image
    {
        bool * brainmaskptr = static_cast<bool *> (Mask->data);
        segPrecisionTYPE * Inputptrtmp = static_cast<segPrecisionTYPE *>(input->data);
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
        CurrSizes->rescale_max[udir]=tempmax;
        CurrSizes->rescale_min[udir]=tempmin;
        if(verbose>0)
        {
            cout << "Normalization["<<udir<<"] = ["<<tempmin<<","<<tempmax<<"]"<<endl;
        }
        Inputptr=&Inputptrtmp[numel*udir];
        brainmaskptr = static_cast<bool *> (Mask->data);
        bool nanflag=false;
        for (int i=0; i<numel; i++)
        {
            //if(brainmaskptr[i]>0){
            //log(number_between_0_and_1 + 1)/log(2)
            Inputptr[i]=logf((((Inputptr[i])-tempmin)/(tempmax-tempmin))+1)/0.693147181;
            if(Inputptr[i]!=Inputptr[i])
            {
                if(nanflag==0){
                    cout<< "Warning: The image at timepoint="<<udir<<" has NaNs. This can cause problems." << endl;
                    nanflag=true;
                }
            }
            /*}
            else{
                Inputptr[i]=0;
            }*/
        }
    }
    return 1;
}

int Normalize_Image(nifti_image * input,
                    ImageSize * CurrSizes,
                    bool verbose)
{
    if(input->datatype!=NIFTI_TYPE_FLOAT32)
    {
        seg_changeDatatype<float>(input);
    }
    if(verbose>0)
    {
        cout<< "Normalizing Input Image" << endl;
    }
    if(input->datatype==NIFTI_TYPE_FLOAT32)
    {
        // if mask is not set up
        for(long udir=0; udir<CurrSizes->usize; udir++) // Per Multispectral Image
        {
            int numel=(int)(input->nx*input->ny*input->nz);
            segPrecisionTYPE * Inputptrtmp = static_cast<segPrecisionTYPE *>(input->data);
            segPrecisionTYPE * Inputptr=&Inputptrtmp[numel*udir];

            float tempmax=0;
            float tempmin=1000000.f;
            for (int i=0; i<numel; i++)
            {
                if (*Inputptr<tempmin)
                {
                    tempmin=*Inputptr;
                }
                if (*Inputptr>tempmax)
                {
                    tempmax=*Inputptr;
                }
                Inputptr++;
            }
            CurrSizes->rescale_max[udir]=tempmax;
            CurrSizes->rescale_min[udir]=tempmin;
            Inputptr=&Inputptrtmp[numel*udir];
            bool nanflag=false;
            for (int i=0; i<numel; i++)
            {
                //log(number_between_0_and_1 + 1)/log(2)
                *Inputptr=logf((((*Inputptr)-tempmin)/(tempmax-tempmin))+1)/0.693147181;
                if(*Inputptr!=*Inputptr)
                {
                    if(nanflag==0){
                        cout<< "Warning: The image at timepoint="<<udir<<" has NaNs. This can cause problems." << endl;
                        nanflag=true;
                    }
                }
                Inputptr++;
            }
        }
    }

    return 1;
}

int Normalize_T1_and_MV(nifti_image * T1,
                        nifti_image * Mask,
                        segPrecisionTYPE * M,
                        segPrecisionTYPE * V,
                        ImageSize * CurrSizes)
{

    segPrecisionTYPE * T1ptrtmp = static_cast<segPrecisionTYPE *>(T1->data);
    bool * brainmaskptr = static_cast<bool *>(Mask->data);
    segPrecisionTYPE * T1ptr=T1ptrtmp;
    int numel=(int)(T1->nx*T1->ny*T1->nz);
    float tempmax=0;
    float tempmin=100000;
    for (int i=0; i<numel; i++)
    {
        if(*brainmaskptr)
        {
            if (*T1ptr<tempmin)
            {
                tempmin=*T1ptr;
            }
            if (*T1ptr>tempmax)
            {
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
    for (int i=0; i<numel; i++)
    {
        if(*brainmaskptr)
        {
            *T1ptr=logf((((*T1ptr)-tempmin)/(tempmax-tempmin))+1)/0.693147181;
        }
        else
        {
            *T1ptr=0;
        }
        brainmaskptr++;
        T1ptr++;
    }
    for(long cl=0; cl<CurrSizes->numclass; cl++)
    {
        M[cl]=log(1+(M[cl]-tempmin)/(tempmax-tempmin));
        V[cl]=log(1+(V[cl])/(tempmax-tempmin));
    }

    return 1;
}

int * Create_Short_2_Long_Matrix_from_NII(nifti_image * Mask,
                                          long * shortsize)
{
    int numel_masked=0;
    int numel = Mask->nvox;
    if(Mask->datatype==DT_BINARY)
    {
        bool * Maskptr = static_cast<bool *> (Mask->data);
        bool * Maskptrtmp = Maskptr;
        for (int i=0; i<numel; i++, Maskptrtmp++)
        {
            (*Maskptrtmp)>0?numel_masked++:0;
        }
        shortsize[0]=numel_masked;

        int * Short_2_Long_Indices= new int [numel_masked]();
        int * Short_2_Long_Indices_PTR = (int *)(Short_2_Long_Indices);

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
        return Short_2_Long_Indices;
    }
    else
    {
        printf("err\tCreate_Short_2_Long_Matrix\tWrong Mask datatype\n");
        return NULL;

    }
}

int *  Create_Long_2_Short_Matrix_from_NII(nifti_image * Mask)
{
    int numel = Mask->nvox;
    int * Long_2_Short_Indices= new int [numel]();
    if(Mask->datatype==DT_BINARY)
    {
        bool * Maskptr = static_cast<bool *> (Mask->data);
        bool * Maskptrtmp = Maskptr;
        int * Long_2_Short_Indices_PTR = (int *) Long_2_Short_Indices;

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
        return Long_2_Short_Indices;
    }
    else
    {
        cout<< "err\tCreate_Correspondace_Matrices\tWrong Mask datatype\n" << endl;
    }
    return Long_2_Short_Indices;
}


int * Create_Short_2_Long_Matrix_from_Carray(bool * Mask,
                                             int * shortsize,
                                             int nvox)
{
    int numel_masked=0;
    int numel = nvox;
    for (int i=0; i<numel; i++)
    {
        Mask[i]?numel_masked++:0;
    }
    shortsize[0]=numel_masked;

    int * Short_2_Long_Indices= new int [numel_masked]();
    int * Short_2_Long_Indices_PTR = (int *)(Short_2_Long_Indices);

    int tempindex=0;
    for (int i=0; i<numel; i++)
    {
        if ((Mask[i])>0)
        {
            Short_2_Long_Indices_PTR[tempindex]=i;
            tempindex++;
        }

    }
    return Short_2_Long_Indices;

}

int *  Create_Long_2_Short_Matrix_from_Carray(bool * Mask,
                                              int nvox)
{
    int * Long_2_Short_Indices= new int [nvox]();
    int * Long_2_Short_Indices_PTR = (int *) Long_2_Short_Indices;
    int tempindex=0;
    for (int i=0; i<nvox; i++,Long_2_Short_Indices_PTR++)
    {
        if ((Mask[i])>0)
        {
            (*Long_2_Short_Indices_PTR)=tempindex;
            tempindex++;
        }
        else
        {
            (*Long_2_Short_Indices_PTR)=-1;
        }
    }
    return Long_2_Short_Indices;

}

segPrecisionTYPE * Create_cArray_from_Prior_mask(nifti_image * Mask,
                                                 nifti_image * Priors,
                                                 long numclass,
                                                 bool PV_ON)
{
    register long numel=(int)(Mask->nx*Mask->ny*Mask->nz);
    register long numel_masked=0;

    bool * Maskptrtmp = static_cast<bool *> (Mask->data);;
    for (long i=0; i<numel; i++, Maskptrtmp++)
    {
        *Maskptrtmp?numel_masked++:0;
    }
    int pluspv=(int)(PV_ON)*2;

    segPrecisionTYPE * Expec = new segPrecisionTYPE [numel_masked*(numclass+pluspv)] ();
    segPrecisionTYPE * tempExpec= (segPrecisionTYPE *) Expec;
    segPrecisionTYPE * PriorPTR = static_cast<segPrecisionTYPE *>(Priors->data);
    for(long cl=0; cl<numclass; cl++)
    {
        Maskptrtmp = static_cast<bool *> (Mask->data);;
        for (int i=numel; i--; Maskptrtmp++,PriorPTR++)
        {
            if(*Maskptrtmp)
            {
                *tempExpec = *PriorPTR;
                tempExpec++;
            }
        }
    }

    return Expec;
}

segPrecisionTYPE * Create_cArray_from_Prior(nifti_image * Priors,
                                            long numclass,
                                            bool PV_ON)
{
    register long numel=(int)(Priors->nx*Priors->ny*Priors->nz);
    long pluspv=(int)(PV_ON)*2;
    segPrecisionTYPE * Expec = new segPrecisionTYPE [numel*(numclass+pluspv)] ();
    segPrecisionTYPE * Expec_PTR= Expec;
    segPrecisionTYPE * PriorPTR = static_cast<segPrecisionTYPE *>(Priors->data);
    for(long cl=0; cl<numclass; cl++)
    {
        for (int i=numel; i--; PriorPTR++,Expec_PTR++)
        {
            *Expec_PTR = *PriorPTR;
        }
    }
    return Expec;
}

segPrecisionTYPE * Create_cArray_from_3D_image(nifti_image * Mask,
                                               nifti_image * SourceImage)
{
    register int numel=(int)(Mask->nx*Mask->ny*Mask->nz);
    register int numel_masked=0;

    bool * Maskptrtmp = static_cast<bool *> (Mask->data);;
    for (int i=0; i<numel; i++, Maskptrtmp++)
    {
        *Maskptrtmp?numel_masked++:0;
    }

    segPrecisionTYPE * outimage = new segPrecisionTYPE [numel_masked] ();
    segPrecisionTYPE * outimage_ptr= outimage;
    segPrecisionTYPE * SourceImagePTR = static_cast<segPrecisionTYPE *>(SourceImage->data);
    Maskptrtmp = static_cast<bool *> (Mask->data);

    for (int i=numel; i--; Maskptrtmp++,SourceImagePTR++)
    {
        if(*Maskptrtmp)
        {
            *outimage_ptr = *SourceImagePTR;
            outimage_ptr++;
        }
    }


    return outimage;
}


int calcE_mask(nifti_image * T1,
               segPrecisionTYPE * IterPrior,
               segPrecisionTYPE * Expec,
               double * loglik,
               segPrecisionTYPE * BiasField,
               segPrecisionTYPE * Outlierness,
               segPrecisionTYPE OutliernessThreshold,
               int * S2L,
               segPrecisionTYPE * M,
               segPrecisionTYPE * V,
               ImageSize * CurrSizes,
               int verbose)
{
    int numel_masked=CurrSizes->numelmasked;
    int num_class=CurrSizes->numclass;
    bool OutliernessFlag=(Outlierness==NULL)?0:1;
    segPrecisionTYPE inv_v [maxNumbClass*maxMultispectalSize*maxMultispectalSize]= {0.0f};
    segPrecisionTYPE inv_sqrt_V_2pi [maxNumbClass]= {0.0f};

    int Expec_offset [maxNumbClass]= {0};

    for (int cl=0; cl<num_class; cl++)
    {
        Expec_offset[cl]=(int) cl*numel_masked;
        if(CurrSizes->usize>1)
        {
            seg_Matrix <double> Vmat(CurrSizes->usize,CurrSizes->usize);

            for(long j2=0; j2<CurrSizes->usize; j2++)
            {
                for(long i2=j2; i2<CurrSizes->usize; i2++)
                {
                    Vmat.setvalue(i2,j2,(double)(V[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]));
                    Vmat.setvalue(j2,i2,(double)(V[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]));
                }
            }
            inv_sqrt_V_2pi[cl]=1/(sqrtf(2*M_PI*Vmat.determinant()));
            if(verbose>1)
            {
                cout<<endl<<"inv_sqrt_V_2pi["<< cl <<"]= "<< inv_sqrt_V_2pi[cl] << endl;
                flush(cout);
            }
            Vmat.invert();
            double cvalue=0.0f;
            bool success;
            if(verbose>1)
            {
                cout<<"inv_V["<< cl <<"]= ";
                flush(cout);
            }
            for(long j2=0; j2<CurrSizes->usize; j2++)
            {
                if(verbose>1)
                {
                    if(j2!=0)
                    {
                        cout<< endl << "          ";
                    }
                }
                for(long i2=0; i2<CurrSizes->usize; i2++)
                {
                    Vmat.getvalue(i2,j2,cvalue,success);
                    inv_v[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]=(segPrecisionTYPE)(cvalue);
                    if(verbose>1)
                    {
                        cout<<inv_v[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]<< "\t";
                        flush(cout);
                    }
                }

            }
            if(verbose>1)
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
    loglik[0]=0;

    //int * Expec_offset_PTR= (int *) Expec_offset;

    float logliktmp=0.0f;


#ifdef _OPENMP
    float * loglikthread = new float [omp_get_max_threads()]();
    for(long i=0; i<(long)omp_get_max_threads(); i++)
        loglikthread[i]=0;

#pragma omp parallel for shared(Expec,loglikthread,T1,BiasField,Outlierness,IterPrior)
#endif
    for (int i=0; i<numel_masked; i++)
    {
        segPrecisionTYPE * T1_PTR = static_cast<segPrecisionTYPE *>(T1->data);
        segPrecisionTYPE T1_Bias_corr[maxMultispectalSize];
        segPrecisionTYPE SumExpec=0.0f;

        for(long Multispec=0; Multispec<CurrSizes->usize; Multispec++)
            T1_Bias_corr[Multispec]=(BiasField!=NULL)?(T1_PTR[S2L[i]+Multispec*CurrSizes->numel] + BiasField[i+Multispec*numel_masked]):(T1_PTR[S2L[i]+Multispec*CurrSizes->numel]);



        //Expec_offset_PTR=Expec_offset;

        for (int cl=0; cl<num_class; cl++)
        {
            segPrecisionTYPE mahal=0.0f;
            for(long Multispec=0; Multispec<CurrSizes->usize; Multispec++)
            {
                segPrecisionTYPE tmpT1_BC_minusM=(T1_Bias_corr[Multispec] - M[cl*(CurrSizes->usize)+Multispec]);
                for(long Multispec2=0; Multispec2<CurrSizes->usize; Multispec2++)
                {
                    mahal-=(0.5f)*(T1_Bias_corr[Multispec2] - M[cl*(CurrSizes->usize)+Multispec2])*inv_v[cl*CurrSizes->usize*CurrSizes->usize+Multispec+Multispec2*CurrSizes->usize]*tmpT1_BC_minusM;
                }
            }

            if(OutliernessFlag)
            {
                float outvalue=(expf(mahal)+0.01)/(expf(mahal)+expf(-0.5*(OutliernessThreshold*OutliernessThreshold))+0.01);
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

    loglik[0]=logliktmp;
    return 1;
}


int calcE(nifti_image * T1,
          segPrecisionTYPE * IterPrior,
          segPrecisionTYPE * Expec,
          double * loglik,
          segPrecisionTYPE * BiasField,
          segPrecisionTYPE * Outlierness,
          segPrecisionTYPE OutliernessThreshold,
          segPrecisionTYPE * M,
          segPrecisionTYPE * V,
          ImageSize * CurrSizes,
          int verbose)
{
    int numel=CurrSizes->numel;
    int num_class=CurrSizes->numclass;
    bool OutliernessFlag=(Outlierness==NULL)?0:1;


    segPrecisionTYPE * IterPrior_PTR= (segPrecisionTYPE *) IterPrior;
    segPrecisionTYPE * Expec_PTR= (segPrecisionTYPE *) Expec;
    segPrecisionTYPE * Outlierness_PTR= (segPrecisionTYPE *) Outlierness;
    segPrecisionTYPE * T1_PTR = static_cast<segPrecisionTYPE *>(T1->data);
    segPrecisionTYPE inv_v [maxNumbClass*maxMultispectalSize*maxMultispectalSize]= {0.0f};
    segPrecisionTYPE inv_sqrt_V_2pi [maxNumbClass]= {0.0f};
    segPrecisionTYPE tmpT1_BC_minusM=0;

    int Expec_offset [maxNumbClass]= {0};

    for (int cl=0; cl<num_class; cl++)
    {
        Expec_offset[cl]=(int) cl*numel;
        if(CurrSizes->usize>1)
        {
            seg_Matrix <double> Vmat(CurrSizes->usize,CurrSizes->usize);

            for(long j2=0; j2<CurrSizes->usize; j2++)
            {
                for(long i2=j2; i2<CurrSizes->usize; i2++)
                {
                    Vmat.setvalue(i2,j2,(double)(V[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]));
                    Vmat.setvalue(j2,i2,(double)(V[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]));
                }
            }
            inv_sqrt_V_2pi[cl]=1/(sqrtf(2*M_PI)*Vmat.determinant());
            if(verbose>1)
            {
                cout<<endl<<"inv_sqrt_V_2pi["<< cl <<"]= "<< inv_sqrt_V_2pi[cl] << endl;
                flush(cout);
            }
            Vmat.invert();
            double cvalue=0.0f;
            bool success;
            if(verbose>1)
            {
                cout<<"inv_V["<< cl <<"]= ";
                flush(cout);
            }
            for(long j2=0; j2<CurrSizes->usize; j2++)
            {
                if(verbose>1)
                {
                    if(j2!=0)
                    {
                        cout<< endl << "          ";
                    }
                }
                for(long i2=0; i2<CurrSizes->usize; i2++)
                {
                    Vmat.getvalue(i2,j2,cvalue,success);
                    inv_v[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]=(segPrecisionTYPE)(cvalue);
                    if(verbose>1)
                    {
                        cout<<inv_v[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]<< "\t";
                        flush(cout);
                    }
                }

            }
            if(verbose>1)
            {
                cout<< endl;
            }
        }
        else
        {
            inv_sqrt_V_2pi[cl]=1/(sqrtf(2*M_PI)* V[cl]);
            inv_v[cl]=1/V[cl];
        }
    }
    loglik[0]=0;



    float logliktmp=0.0f;
#ifdef _OPENMP
    float * loglikthread = new float [omp_get_max_threads()]();
    for(long i=0; i<(long)omp_get_max_threads(); i++)
        loglikthread[i]=0;

#pragma omp parallel for shared(Expec,loglikthread,T1,BiasField,Outlierness,IterPrior)
#endif
    for (int i=0; i<numel; i++)
    {
        segPrecisionTYPE T1_Bias_corr[maxMultispectalSize]= {0.0f};
        for(long Multispec=0; Multispec<CurrSizes->usize; Multispec++)
        {
            T1_Bias_corr[Multispec]=(BiasField!=NULL)?(T1_PTR[i+Multispec*numel] + BiasField[i+Multispec*numel]):(T1_PTR[i+Multispec*numel]);
        }
        segPrecisionTYPE mahal=0.0f;
        segPrecisionTYPE SumExpec=0.0f;

        //Expec_offset_PTR=Expec_offset;

        for (int cl=0; cl<num_class; cl++)
        {
            mahal=0.0f;
            for(long Multispec=0; Multispec<CurrSizes->usize; Multispec++)
            {
                tmpT1_BC_minusM=(T1_Bias_corr[Multispec] - M[cl*(CurrSizes->usize)+Multispec]);
                for(long Multispec2=0; Multispec2<CurrSizes->usize; Multispec2++)
                {
                    mahal-=(0.5f)*(T1_Bias_corr[Multispec2] - M[cl*(CurrSizes->usize)+Multispec2])*inv_v[cl*CurrSizes->usize*CurrSizes->usize+Multispec+Multispec2*CurrSizes->usize]*tmpT1_BC_minusM;
                }
            }
            Expec_PTR[i+Expec_offset[cl]]=IterPrior_PTR[i+Expec_offset[cl]] * expf(mahal) * inv_sqrt_V_2pi[cl];
            if(OutliernessFlag)
            {
                Outlierness_PTR[i+Expec_offset[cl]]=(expf(mahal))/(expf(mahal)+expf(-0.5*(OutliernessThreshold*OutliernessThreshold)));
            }
            SumExpec+=Expec_PTR[i+Expec_offset[cl]];
        }
        if (SumExpec<=0.0 || SumExpec!=SumExpec)
        {
            SumExpec=0.01;
        }

        if (SumExpec<=0.0 || SumExpec!=SumExpec)
        {
            for (int cl=0; cl<num_class; cl++)
            {
                Expec_PTR[i+Expec_offset[cl]]=(float)(1)/(float)(num_class);
            }
        }
        else
        {
            for (int cl=0; cl<num_class; cl++)
            {
                Expec_PTR[i+Expec_offset[cl]]=Expec_PTR[i+Expec_offset[cl]]/SumExpec;
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


    loglik[0]=logliktmp;
    return 1;
}




/*
int calcE_mask_aprox(nifti_image * T1,
                      SegPrecisionTYPE * IterPrior,
                      SegPrecisionTYPE * Expec,
                      SegPrecisionTYPE * loglik,
                      SegPrecisionTYPE * BiasField,
                      int * S2L,
                      SegPrecisionTYPE * M,
                      SegPrecisionTYPE * V,
                      ImageSize * CurrSizes,
                      int verbose)
{
    int numel_masked=CurrSizes->numelmasked;
    int num_class=CurrSizes->numclass;

    SegPrecisionTYPE * IterPrior_PTR= (SegPrecisionTYPE *) IterPrior;
    SegPrecisionTYPE * Expec_PTR= (SegPrecisionTYPE *) Expec;
    SegPrecisionTYPE * T1_PTR = static_cast<SegPrecisionTYPE *>(T1->data);
    SegPrecisionTYPE SumExpec=0.0f;
    SegPrecisionTYPE expectmp=0.0f;



    SegPrecisionTYPE inv_v [max_numbclass*MaxMultispectalSize*MaxMultispectalSize]={0.0f};
    SegPrecisionTYPE inv_sqrt_V_2pi [max_numbclass]={0.0f};
    SegPrecisionTYPE newM [max_numbclass*MaxMultispectalSize]={0.0f};
    SegPrecisionTYPE T1_Bias_corr[MaxMultispectalSize]={0.0f};
    SegPrecisionTYPE tmpT1_BC_minusM=0;

    int Expec_offset [max_numbclass]={0};

    for (int cl=0; cl<num_class; cl++) {
        Expec_offset[cl]=(int) cl*numel_masked;
        if(CurrSizes->usize>1){
            matrix <double> Vmat(CurrSizes->usize,CurrSizes->usize);

            for(long j2=0; j2<CurrSizes->usize; j2++){
                for(long i2=j2; i2<CurrSizes->usize; i2++){
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
            for(long j2=0; j2<CurrSizes->usize; j2++){
                if(verbose>1){
                    if(j2!=0){
                        cout<< endl << "          ";
                    }
                }
                for(long i2=0; i2<CurrSizes->usize; i2++){
                    Vmat.getvalue(i2,j2,cvalue,success);


                    inv_v[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]=(SegPrecisionTYPE)(cvalue);
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
    register SegPrecisionTYPE tempvar=0.0f;
    float mahal=0.0f;
    float logliktmp=0.0f;
    for (int i=0; i<numel_masked;i++, Expec_PTR++, IterPrior_PTR++) {
        for(long Multispec=0; Multispec<CurrSizes->usize; Multispec++) {
            T1_Bias_corr[Multispec]=(BiasField!=NULL)?(T1_PTR[S2L[i]+Multispec*CurrSizes->numel] + BiasField[i+Multispec*numel_masked]):(T1_PTR[S2L[i]+Multispec*CurrSizes->numel]);
        }
        SumExpec=0.0f;
        expectmp=0.0f;
        mahal=0.0f;
        //Expec_offset_PTR=Expec_offset;
        for (int cl=0; cl<num_class; cl++) {
            mahal=0.0f;
            for(long Multispec=0; Multispec<CurrSizes->usize; Multispec++) {
                tmpT1_BC_minusM=(T1_Bias_corr[Multispec] - M[cl*(CurrSizes->usize)+Multispec]);
                for(long Multispec2=0; Multispec2<CurrSizes->usize; Multispec2++) {
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
    return 1;
}


*/
/*
int calcE_aprox(nifti_image * T1,
                 SegPrecisionTYPE * IterPrior,
                 SegPrecisionTYPE * Expec,
                 SegPrecisionTYPE * loglik,
                 SegPrecisionTYPE * BiasField,
                 SegPrecisionTYPE * M,
                 SegPrecisionTYPE * V,
                 ImageSize * CurrSizes,
                 int verbose)
{
    int numel=CurrSizes->numel;
    int num_class=CurrSizes->numclass;

    SegPrecisionTYPE * IterPrior_PTR= (SegPrecisionTYPE *) IterPrior;
    SegPrecisionTYPE * Expec_PTR= (SegPrecisionTYPE *) Expec;
    SegPrecisionTYPE * T1_PTR = static_cast<SegPrecisionTYPE *>(T1->data);
    SegPrecisionTYPE SumExpec=0.0f;
    SegPrecisionTYPE expectmp=0.0f;



    SegPrecisionTYPE inv_v [max_numbclass*MaxMultispectalSize*MaxMultispectalSize]={0.0f};
    SegPrecisionTYPE inv_sqrt_V_2pi [max_numbclass]={0.0f};
    SegPrecisionTYPE newM [max_numbclass*MaxMultispectalSize]={0.0f};
    SegPrecisionTYPE T1_Bias_corr[MaxMultispectalSize]={0.0f};
    SegPrecisionTYPE tmpT1_BC_minusM=0;

    int Expec_offset [max_numbclass]={0};

    for (int cl=0; cl<num_class; cl++) {
        Expec_offset[cl]=(int) cl*numel;
        if(CurrSizes->usize>1){
            matrix <double> Vmat(CurrSizes->usize,CurrSizes->usize);

            for(long j2=0; j2<CurrSizes->usize; j2++){
                for(long i2=j2; i2<CurrSizes->usize; i2++){
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
            for(long j2=0; j2<CurrSizes->usize; j2++){
                if(verbose>1){
                    if(j2!=0){
                        cout<< endl << "          ";
                    }
                }
                for(long i2=0; i2<CurrSizes->usize; i2++){
                    Vmat.getvalue(i2,j2,cvalue,success);
                    inv_v[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]=(SegPrecisionTYPE)(cvalue);
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
    register SegPrecisionTYPE tempvar=0.0f;
    float mahal=0.0f;
    float logliktmp=0.0f;
    IterPrior_PTR= (SegPrecisionTYPE *) IterPrior;
    Expec_PTR= (SegPrecisionTYPE *) Expec;
    for (int i=0; i<numel;i++, Expec_PTR++, IterPrior_PTR++) {
        for(long Multispec=0; Multispec<CurrSizes->usize; Multispec++) {
            T1_Bias_corr[Multispec]=(BiasField!=NULL)?(T1_PTR[i+Multispec*numel] + BiasField[i+Multispec*numel]):(T1_PTR[i+Multispec*numel]);
        }
        SumExpec=0.0f;
        expectmp=0.0f;
        mahal=0.0f;
        //Expec_offset_PTR=Expec_offset;
        for (int cl=0; cl<num_class; cl++) {
            mahal=0.0f;
            for(long Multispec=0; Multispec<CurrSizes->usize; Multispec++) {
                tmpT1_BC_minusM=(T1_Bias_corr[Multispec] - M[cl*(CurrSizes->usize)+Multispec]);
                for(long Multispec2=0; Multispec2<CurrSizes->usize; Multispec2++) {
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
    return 1;
}



*/
int calcM(nifti_image * T1,
          segPrecisionTYPE * Expec,
          segPrecisionTYPE * BiasField,
          segPrecisionTYPE * Outlierness,
          segPrecisionTYPE * M,
          segPrecisionTYPE * V,
          segPrecisionTYPE * M_MAP,
          segPrecisionTYPE * V_MAP,
          segPrecisionTYPE reg_factor,
          ImageSize * CurrSizes,
          int verbose)
{

    bool OutliernessFlag=(Outlierness==NULL)?0:1;

    if(verbose>0)
    {
        cout<< "Optimising Gaussian Parameters" << endl;
        flush(cout);
    }
    int numel=CurrSizes->numel;
    int currentnum_class=CurrSizes->numclass;
    segPrecisionTYPE * Expec_PTR = (segPrecisionTYPE *) Expec;
    segPrecisionTYPE * OutliernessPTR = (segPrecisionTYPE *) Outlierness;
    segPrecisionTYPE * T1_PTR = static_cast<segPrecisionTYPE *>(T1->data);
    segPrecisionTYPE * BiasField_PTR= (segPrecisionTYPE *) BiasField;
    segPrecisionTYPE * T1_PTR2= static_cast<segPrecisionTYPE *>(T1->data);
    segPrecisionTYPE * BiasField_PTR2= (segPrecisionTYPE *) BiasField;

    int Expec_offset[maxNumbClass];
    for (int cl=0; cl<currentnum_class; cl++)
    {
        Expec_offset[cl]=cl*numel;
    }
    segPrecisionTYPE tempsum= (segPrecisionTYPE) 0.0;
    segPrecisionTYPE SumPriors= (segPrecisionTYPE) 0.0;

    // ***********

    for (int cl=0; cl<currentnum_class; cl++)
    {
        // MEAN
        for(long Multispec=0; Multispec<CurrSizes->usize; Multispec++)
        {
            Expec_PTR=(segPrecisionTYPE *) &Expec[Expec_offset[cl]];
            OutliernessPTR=(segPrecisionTYPE *) &Outlierness[Expec_offset[cl]];
            T1_PTR = static_cast<segPrecisionTYPE *>(T1->data);
            T1_PTR =&T1_PTR[Multispec*CurrSizes->numel];
            tempsum=(segPrecisionTYPE)0.0;
            SumPriors=(segPrecisionTYPE)0.0;


            if(OutliernessFlag)
            {

                if(BiasField!=NULL)
                {
                    BiasField_PTR= &BiasField[Multispec*numel];
                    for (int i=0; i<numel; i++, Expec_PTR++,OutliernessPTR++,BiasField_PTR++)
                    {
                        float current_value=(*Expec_PTR) * (*OutliernessPTR)*(T1_PTR[i]+(*BiasField_PTR));
                        if(current_value==current_value)
                        {
                            tempsum+=current_value;
                            SumPriors+=(*Expec_PTR)*(*OutliernessPTR);
                        }
                    }
                }
                else
                {
                    for (int i=0; i<numel; i++, Expec_PTR++,OutliernessPTR++)
                    {
                        float current_value=(*Expec_PTR) * (*OutliernessPTR)*(T1_PTR[i]);
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
                    BiasField_PTR= &BiasField[Multispec*numel];
                    for (int i=0; i<numel; i++, Expec_PTR++,BiasField_PTR++)
                    {
                        float current_value=(*Expec_PTR) * (T1_PTR[i]+(*BiasField_PTR));
                        if(current_value==current_value)
                        {
                            tempsum+=current_value;
                            SumPriors+=(*Expec_PTR);
                        }
                    }
                }
                else
                {
                    for (int i=0; i<numel; i++, Expec_PTR++)
                    {
                        float current_value=(*Expec_PTR) * (T1_PTR[i]);
                        if(current_value==current_value)
                        {
                            tempsum+=current_value;
                            SumPriors+=(*Expec_PTR);
                        }
                    }
                }

            }


            if(M_MAP==NULL)
            {
                M[cl*(CurrSizes->usize)+Multispec]=tempsum/SumPriors;
            }
            else
            {
                M[cl*(CurrSizes->usize)+Multispec]=(tempsum/SumPriors/powf(V[cl*(CurrSizes->usize)+Multispec],2)+M_MAP[cl*(CurrSizes->usize)+Multispec]/powf(V_MAP[cl*(CurrSizes->usize)+Multispec],2))/(1/powf(V[cl*(CurrSizes->usize)+Multispec],2)+1/powf(V_MAP[cl*(CurrSizes->usize)+Multispec],2));
            }

            for(long Multispec2=Multispec; Multispec2<CurrSizes->usize; Multispec2++)
            {

                T1_PTR = static_cast<segPrecisionTYPE *>(T1->data);
                T1_PTR =&T1_PTR[Multispec*CurrSizes->numel];

                T1_PTR2 = static_cast<segPrecisionTYPE *>(T1->data);
                T1_PTR2 =&T1_PTR2[Multispec2*CurrSizes->numel];
                float tmpM=M[cl*CurrSizes->usize+Multispec];
                float tmpM2=M[cl*CurrSizes->usize+Multispec2];
                //STD
                tempsum=0;
                Expec_PTR=&Expec[Expec_offset[cl]];
                OutliernessPTR=(segPrecisionTYPE *) &Outlierness[Expec_offset[cl]];
                if(OutliernessFlag)
                {
                    if(BiasField!=NULL)
                    {
                        BiasField_PTR=&BiasField[Multispec*numel];
                        BiasField_PTR2=&BiasField[Multispec2*numel];
                        for (int i=0; i<numel; i++,Expec_PTR++,BiasField_PTR++,OutliernessPTR++,BiasField_PTR2++)
                        {
                            tempsum+=(*Expec_PTR) *(*OutliernessPTR)* (T1_PTR[i]+(*BiasField_PTR)-tmpM) * (T1_PTR2[i]+(*BiasField_PTR2)-tmpM2);
                        }
                    }
                    else
                    {
                        for (int i=0; i<numel; i++,Expec_PTR++,OutliernessPTR++)
                        {
                            tempsum+=(*Expec_PTR) *(*OutliernessPTR)* (T1_PTR[i]-tmpM) * (T1_PTR2[i]-tmpM2);
                        }

                    }
                }
                else
                {
                    if(BiasField!=NULL)
                    {
                        BiasField_PTR=&BiasField[Multispec*numel];
                        BiasField_PTR2=&BiasField[Multispec2*numel];
                        for (int i=0; i<numel; i++,Expec_PTR++,BiasField_PTR++,BiasField_PTR2++)
                        {
                            tempsum+=(*Expec_PTR) * (T1_PTR[i]+(*BiasField_PTR)-tmpM) * (T1_PTR2[i]+(*BiasField_PTR2)-tmpM2);
                        }
                    }
                    else
                    {
                        for (int i=0; i<numel; i++,Expec_PTR++)
                        {
                            tempsum+=(*Expec_PTR) * (T1_PTR[i]-tmpM) * (T1_PTR2[i]-tmpM2);
                        }

                    }
                }

                V[cl*CurrSizes->usize*CurrSizes->usize+Multispec+Multispec2*CurrSizes->usize]=tempsum/SumPriors+0.00001f;
                V[cl*CurrSizes->usize*CurrSizes->usize+Multispec2+Multispec*CurrSizes->usize]=V[cl*CurrSizes->usize*CurrSizes->usize+Multispec+Multispec2*CurrSizes->usize];
                if(Multispec==Multispec2)
                    V[cl*CurrSizes->usize*CurrSizes->usize+Multispec+Multispec2*CurrSizes->usize]*=reg_factor;

            }
        }
    }

    for (int cl=0; cl<currentnum_class; cl++)
    {
        for(long Multispec=0; Multispec<CurrSizes->usize; Multispec++)
        {
            for(long Multispec2=0; Multispec2<CurrSizes->usize; Multispec2++)
            {
                V[cl*CurrSizes->usize*CurrSizes->usize+Multispec+Multispec2*CurrSizes->usize]/=reg_factor;
            }
        }
    }

    for (int cl=0; cl<currentnum_class; cl++)
    {
        if(verbose>0)
        {
            if(CurrSizes->usize==1)
            {
                cout.fill('0');
                cout<< "M["<<(int)(cl)<<"]= "<<setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(M[cl])<<"\tV["<<(int)(cl)<<"]="<<setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(V[cl])<< endl;
                flush(cout);
            }
            else
            {

                cout<< "M["<<(int)(cl)<<"]= ";
                for(long Multispec=0; Multispec<CurrSizes->usize; Multispec++)
                {
                    cout<< setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(M[cl*CurrSizes->usize+Multispec])<<"\t";
                }
                cout<< endl;
                flush(cout);
                cout<< "V["<<(int)(cl)<<"]= ";
                for(long Multispec=0; Multispec<CurrSizes->usize; Multispec++)
                {
                    if(Multispec>0)
                    {
                        cout<< "      ";
                    }
                    for(long Multispec2=0; Multispec2<CurrSizes->usize; Multispec2++)
                    {
                        cout<< setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(V[cl*CurrSizes->usize*CurrSizes->usize+Multispec*CurrSizes->usize+Multispec2])<<"\t";
                    }
                    cout<< endl;
                }
                cout<< endl;
                flush(cout);
            }
        }
    }

    return 1;
}




int calcM_mask(nifti_image * T1,
               segPrecisionTYPE * Expec,
               segPrecisionTYPE * BiasField,
               segPrecisionTYPE * Outlierness,
               int * S2L,
               segPrecisionTYPE * M,
               segPrecisionTYPE * V,
               segPrecisionTYPE * M_MAP,
               segPrecisionTYPE * V_MAP,
               segPrecisionTYPE reg_factor,
               ImageSize * CurrSizes,
               int verbose)
{

    if(verbose>0)
    {
        cout<< "Optimising Gaussian Parameters" << endl;
        flush(cout);
    }
    bool OutliernessFlag=(Outlierness==NULL)?0:1;

    int numel_masked=CurrSizes->numelmasked;
    int num_class=CurrSizes->numclass;
    int Expec_offset[maxNumbClass];
    for (int cl=0; cl<num_class; cl++)
    {
        Expec_offset[cl]=cl*numel_masked;
    }

#ifdef _OPENMP
#pragma omp parallel for shared(T1,BiasField,Outlierness)
#endif
    // ***********
    for (int cl=0; cl<num_class; cl++)
    {
        int * S2L_PTR = (int *) S2L;
        segPrecisionTYPE * Expec_PTR = (segPrecisionTYPE *) Expec;
        segPrecisionTYPE * OutliernessPTR = (segPrecisionTYPE *) Outlierness;
        segPrecisionTYPE * T1_PTR = static_cast<segPrecisionTYPE *>(T1->data);
        segPrecisionTYPE * BiasField_PTR= (segPrecisionTYPE *) BiasField;
        segPrecisionTYPE tempsum= (segPrecisionTYPE) 0.0;
        segPrecisionTYPE SumPriors= (segPrecisionTYPE) 0.0;
        segPrecisionTYPE * T1_PTR2= static_cast<segPrecisionTYPE *>(T1->data);
        segPrecisionTYPE * BiasField_PTR2= (segPrecisionTYPE *) BiasField;

        // MEAN
        for(long Multispec=0; Multispec<CurrSizes->usize; Multispec++)
        {
            Expec_PTR=(segPrecisionTYPE *) &Expec[Expec_offset[cl]];
            OutliernessPTR=(segPrecisionTYPE *) &Outlierness[Expec_offset[cl]];
            S2L_PTR = (int *) S2L;

            T1_PTR = static_cast<segPrecisionTYPE *>(T1->data);
            T1_PTR =&T1_PTR[Multispec*CurrSizes->numel];
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
                if(M_MAP==NULL)
                {
                    M[cl*(CurrSizes->usize)+Multispec]=tempsum/SumPriors;
                }
                else
                {
                    M[cl*(CurrSizes->usize)+Multispec]=(tempsum/SumPriors/powf(V[cl*(CurrSizes->usize)+Multispec],2)+M_MAP[cl*(CurrSizes->usize)+Multispec]/powf(V_MAP[cl*(CurrSizes->usize)+Multispec],2))/(1/powf(V[cl*(CurrSizes->usize)+Multispec],2)+1/powf(V_MAP[cl*(CurrSizes->usize)+Multispec],2));
                }

                for(long Multispec2=Multispec; Multispec2<CurrSizes->usize; Multispec2++)
                {
                    S2L_PTR = (int *) S2L;

                    T1_PTR = static_cast<segPrecisionTYPE *>(T1->data);
                    T1_PTR =&T1_PTR[Multispec*CurrSizes->numel];

                    T1_PTR2 = static_cast<segPrecisionTYPE *>(T1->data);
                    T1_PTR2 =&T1_PTR2[Multispec2*CurrSizes->numel];
                    float tmpM=M[cl*CurrSizes->usize+Multispec];
                    float tmpM2=M[cl*CurrSizes->usize+Multispec2];
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
                        V[cl*CurrSizes->usize*CurrSizes->usize+Multispec+Multispec2*CurrSizes->usize]=tempsum/SumPriors;
                        if(Multispec2!=Multispec)
                        {
                            V[cl*CurrSizes->usize*CurrSizes->usize+Multispec2+Multispec*CurrSizes->usize]=V[cl*CurrSizes->usize*CurrSizes->usize+Multispec+Multispec2*CurrSizes->usize];
                            V[cl*CurrSizes->usize*CurrSizes->usize+Multispec+Multispec2*CurrSizes->usize]/=reg_factor;
                            V[cl*CurrSizes->usize*CurrSizes->usize+Multispec2+Multispec*CurrSizes->usize]/=reg_factor;
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
            if(CurrSizes->usize==1)
            {
                cout.fill('0');
                cout<< "M["<<(int)(cl)<<"]= "<<setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(M[cl])<<"\tV["<<(int)(cl)<<"]="<<setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(V[cl])<< endl;
                flush(cout);
            }
            else
            {

                cout<< "M["<<(int)(cl)<<"]= ";
                for(long Multispec=0; Multispec<CurrSizes->usize; Multispec++)
                {
                    cout<< setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(M[cl*CurrSizes->usize+Multispec])<<"\t";
                }
                cout<< endl;
                flush(cout);
                cout<< "V["<<(int)(cl)<<"]= ";
                for(long Multispec=0; Multispec<CurrSizes->usize; Multispec++)
                {
                    if(Multispec>0)
                    {
                        cout<< "      ";
                    }
                    for(long Multispec2=0; Multispec2<CurrSizes->usize; Multispec2++)
                    {
                        cout<< setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(V[cl*CurrSizes->usize*CurrSizes->usize+Multispec*CurrSizes->usize+Multispec2])<<"\t";
                    }
                    cout<< endl;
                }
                cout<< endl;
                flush(cout);
            }
        }
    }
    return 1;
}





int calcM_mask_LoAd(nifti_image * T1,
                    segPrecisionTYPE * Expec,
                    segPrecisionTYPE * BiasField,
                    int * S2L,
                    segPrecisionTYPE * M,
                    segPrecisionTYPE * V,
                    ImageSize * CurrSizes,
                    int verbose,
                    bool PVon)
{

    if(verbose>0)
    {
        cout<< "Optimising Gaussian Parameters" << endl;
        flush(cout);
    }
    int numel_masked=CurrSizes->numelmasked;
    int currentnum_class=CurrSizes->numclass;
    int * S2L_PTR = (int *) S2L;
    segPrecisionTYPE * Expec_PTR = (segPrecisionTYPE *) Expec;
    segPrecisionTYPE * T1_PTR = static_cast<segPrecisionTYPE *>(T1->data);
    segPrecisionTYPE * BiasField_PTR= (segPrecisionTYPE *) BiasField;
    int Expec_offset[maxNumbClass];
    for (int cl=0; cl<currentnum_class; cl++)
    {
        Expec_offset[cl]=cl*numel_masked;
    }
    segPrecisionTYPE tempsum= (segPrecisionTYPE) 0.0;
    segPrecisionTYPE SumPriors= (segPrecisionTYPE) 0.0;


    for (int cl=0; cl<5; cl++)
    {
        Expec_PTR=(segPrecisionTYPE *) &Expec[Expec_offset[cl]];
        BiasField_PTR= (segPrecisionTYPE *) BiasField;
        S2L_PTR = (int *) S2L;
        tempsum=(segPrecisionTYPE)0.0;
        SumPriors=(segPrecisionTYPE)0.0;

        // MEAN
        if(BiasField!=NULL)
        {
            for (int i=0; i<numel_masked; i++, Expec_PTR++,S2L_PTR++,BiasField_PTR++)
            {
                tempsum+=(*Expec_PTR) * (T1_PTR[(*S2L_PTR)]+(*BiasField_PTR));
                SumPriors+=(*Expec_PTR);
            }
        }
        else
        {
            for (int i=0; i<numel_masked; i++, Expec_PTR++,S2L_PTR++)
            {
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
        if(BiasField!=NULL)
        {
            for (int i=0; i<numel_masked; i++,Expec_PTR++,S2L_PTR++,BiasField_PTR++)
            {
                tempsum+=(*Expec_PTR) * (T1_PTR[(*S2L_PTR)]+(*BiasField_PTR)-M[cl]) * (T1_PTR[(*S2L_PTR)]+(*BiasField_PTR)-M[cl]);
            }
        }
        else
        {
            for (int i=0; i<numel_masked; i++,Expec_PTR++,S2L_PTR++)
            {
                tempsum+=(*Expec_PTR) * (T1_PTR[(*S2L_PTR)]-M[cl]) * (T1_PTR[(*S2L_PTR)]-M[cl]);
            }

        }
        V[cl]=tempsum/SumPriors;

        if(V[cl]<0.0001 || V[cl]!=V[cl])
        {
            V[cl]=0.0001;
        }

        if(verbose>0)
        {
            cout.fill('0');
            cout << "M[" << (int)(cl) << "]= " << setw(10) << setprecision(7) << left << (segPrecisionTYPE)(M[cl]) << "\tV[" << (int)(cl) << "]=" << setw(10) << setprecision(7) <<left << (segPrecisionTYPE)(V[cl]) << endl;
            flush(cout);
        }
    }

    if(PVon)
    {
        M[WMGMpvclass]=(M[WMclass]*V[GMclass]+M[GMclass]*V[WMclass])/(float)(V[WMclass]+V[GMclass]);
        M[GMCSFpvclass]=(M[CSFclass]*V[GMclass]+M[GMclass]*V[CSFclass])/(float)(V[CSFclass]+V[GMclass]);
        V[WMGMpvclass]=0.5*V[WMclass]+0.5*V[GMclass];
        if(V[WMGMpvclass]<0.0003 || V[WMGMpvclass]!=V[WMGMpvclass])V[WMGMpvclass]=0.0003;
        V[GMCSFpvclass]=0.5*V[GMclass]+(0.5)*V[CSFclass];
        if(V[GMCSFpvclass]<0.0003 || V[GMCSFpvclass]!=V[GMCSFpvclass])V[GMCSFpvclass]=0.0003;
        cout<< "M["<<(int)(WMGMpvclass)<<"]= "<<setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(M[WMGMpvclass])<<"\tV["<<(int)(WMGMpvclass)<<"]="<<setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(V[WMGMpvclass])<< endl;
        cout<< "M["<<(int)(GMCSFpvclass)<<"]= "<<setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(M[GMCSFpvclass])<<"\tV["<<(int)(GMCSFpvclass)<<"]="<<setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(V[GMCSFpvclass])<< endl;
        flush(cout);
    }

    return 1;
}

int Relax_Priors(segPrecisionTYPE * Priors,
                 segPrecisionTYPE * Expec,
                 segPrecisionTYPE * MRF,
                 int * S2L,
                 int * L2S,
                 float RelaxFactor,
                 segPrecisionTYPE * G,
                 segPrecisionTYPE ba,
                 segPrecisionTYPE be,
                 ImageSize *  CurrSizes,
                 seg_EM_Params * segment_param)
{

    if(RelaxFactor<0)
    {
        RelaxFactor=0;
        if(segment_param->verbose_level>0)
        {
            cout<< "RelaxFactor was bellow 0 (min). Assuming RelaxFactor=0." << endl;
        }
    }

    if(RelaxFactor>1)
    {
        RelaxFactor=1;
        if(segment_param->verbose_level>0)
        {
            cout<< "RelaxFactor was above 1 (max). Assuming RelaxFactor=1." << endl;
        }
    }


    if(RelaxFactor>=0)
    {
        if(segment_param->verbose_level>0)
        {
            cout<< "Relaxing Priors with Rf = "<< RelaxFactor << endl;
            flush(cout);
        }
        GaussianFilter4D_cArray(Expec,S2L,L2S,2.0, CurrSizes);
        Relax_Priors_Share(Priors,Expec,RelaxFactor,G,be,CurrSizes);
        for(long i=0; i<(CurrSizes->numclass*CurrSizes->numelmasked); i++)MRF[i]=Priors[i];
        memcpy(Priors,Expec,CurrSizes->numelmasked*CurrSizes->numclass*sizeof(segPrecisionTYPE));
    }

    return 1;

}

int Relax_Priors_Share(segPrecisionTYPE * Priors,
                       segPrecisionTYPE * Expec,
                       float RelaxFactor,
                       segPrecisionTYPE * G,
                       segPrecisionTYPE be,
                       ImageSize * CurrSizes)
{
    long numclass=CurrSizes->numclass;

    int N[maxNumbClass*maxNumbClass];
    segPrecisionTYPE Currexpec[maxNumbClass];
    segPrecisionTYPE Currexpec_updated[maxNumbClass];
    segPrecisionTYPE Currexpec_sumall;
    segPrecisionTYPE Currexpec_class;
    int Class_shift[maxNumbClass];
    // Precompute Class Shift
    for (long cl=0; cl<numclass; cl++)
    {
        Class_shift[cl]=cl*CurrSizes->numelmasked;
    }

    //Find neighbouring rules
    for(long i=0; i<numclass; i++)
    {
        for(long j=0; j<numclass; j++)
        {
            N[i+j*numclass]=(G[i+j*numclass]<=be)?1:0;
        }
    }

    for (long i=0; i<(long)CurrSizes->numelmasked; i++)
    {
        //Copy current expec for all classes
        for (long cl=0; cl<numclass; cl++)
        {
            Currexpec[cl]=(1-RelaxFactor)*Expec[i+Class_shift[cl]]+(RelaxFactor)*Priors[i+Class_shift[cl]];
        }
        //Compute relaxed expec
        Currexpec_sumall=0;
        for (long cl=0; cl<numclass; cl++)
        {
            Currexpec_updated[cl]=0;
            Currexpec_class=Currexpec[cl];
            for (long cl2=0; cl2<numclass; cl2++)
            {
                Currexpec_updated[cl]+=Currexpec_class*Currexpec[cl2]*double(N[cl+cl2*numclass]);
            }
            Currexpec_sumall+=Currexpec_updated[cl];
        }
        Currexpec_sumall=Currexpec_sumall==0?1:Currexpec_sumall;
        //Copy current expec back all classes
        for (long cl=0; cl<numclass; cl++)
        {
            Expec[i+Class_shift[cl]]=Currexpec_updated[cl]/Currexpec_sumall;
        }

    }
    return 1;
}


int Convert_to_PV(nifti_image * T1,
                  segPrecisionTYPE * BiasField,
                  segPrecisionTYPE * ShortPrior,
                  segPrecisionTYPE * Expec,
                  segPrecisionTYPE * MRF,
                  segPrecisionTYPE * M,
                  segPrecisionTYPE * V,
                  int * S2L,
                  int * L2S,
                  ImageSize * CurrSizes,
                  seg_EM_Params * segment_param)
{
    if(segment_param->verbose_level>0)
    {
        cout<< "Covert to Explicit PV modeling" << endl;
    }
    GaussianFilter4D_cArray(Expec,S2L,L2S,1.0,CurrSizes);
    GaussianFilter4D_cArray(MRF,S2L,L2S,1.0,CurrSizes);
    Convert_WM_and_GM_to_PV(T1,BiasField,ShortPrior,Expec,S2L,M,V,CurrSizes);
    CurrSizes->numclass=7;
    for(long i=0; i<(CurrSizes->numclass*CurrSizes->numelmasked); i++)
    {
        MRF[i]=1.0f;
    }
    return 1;
}

int Sulci_and_gyri_correction(segPrecisionTYPE * MRF_Beta,
                              segPrecisionTYPE * ShortPrior,
                              segPrecisionTYPE * Expec,
                              segPrecisionTYPE *MRF,
                              int * S2L,
                              int * L2S,
                              ImageSize *CurrSizes)
{

    // Deep Sulci ->  Seed=WM+WMGMpv+dGM+iCSF
    cout<< "Sucli and Gyri correction" << endl;
    bool * Seed_Mask= new bool [CurrSizes->numelmasked]();
    segPrecisionTYPE * SpeedFunc= new segPrecisionTYPE [CurrSizes->numelmasked]();
    segPrecisionTYPE * wSulci= new segPrecisionTYPE [CurrSizes->numelmasked]();
    segPrecisionTYPE * wGyri= new segPrecisionTYPE [CurrSizes->numelmasked]();

    for(long i=0; i<(long)CurrSizes->numelmasked; i++)
    {
        Seed_Mask[i]=(Expec[i+WMclass*CurrSizes->numelmasked]+
                Expec[i+WMGMpvclass*CurrSizes->numelmasked]+
                Expec[i+dGMclass*CurrSizes->numelmasked]+
                Expec[i+iCSFclass*CurrSizes->numelmasked])>0.5f;
        SpeedFunc[i]=(Expec[i+CSFclass*CurrSizes->numelmasked]+
                Expec[i+GMCSFpvclass*CurrSizes->numelmasked])*2;
    }

    FMM(Seed_Mask, SpeedFunc, wSulci,30, L2S, S2L, CurrSizes);
    TransformGeoTime(wSulci,30, L2S, S2L, CurrSizes);

    for(long i=0; i<(long)CurrSizes->numelmasked; i++)
    {
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
    float wSulci_tmp=0;
    float wGyri_tmp=0;
    float MRFbeta_tmp=0;

    for(long i=0; i<(long)CurrSizes->numelmasked; i++)
    {
        wSulci_tmp=wSulci[i];
        wGyri_tmp=wGyri[i];
        MRFbeta_tmp=(1-wSulci_tmp)*(1-wGyri_tmp);
        MRF_Beta[i]=MRFbeta_tmp;
        sumed_all=0;
        if((wSulci_tmp+wGyri_tmp)>0)
        {
            Expec[i+WMGMpvclass*CurrSizes->numelmasked]=Expec[i+WMGMpvclass*CurrSizes->numelmasked]+wGyri[i]*Expec[i+GMclass*CurrSizes->numelmasked];
            Expec[i+GMCSFpvclass*CurrSizes->numelmasked]=Expec[i+GMCSFpvclass*CurrSizes->numelmasked]+wSulci[i]*Expec[i+GMclass*CurrSizes->numelmasked];
            Expec[i+GMclass*CurrSizes->numelmasked]=Expec[i+GMclass*CurrSizes->numelmasked]*(MRFbeta_tmp);

            for(long c=0; c<7; c++)
            {
                sumed_all+=Expec[i+c*CurrSizes->numelmasked];
            }
            if(sumed_all>0)
            {
                for(long c=0; c<7; c++)
                {
                    Expec[i+c*CurrSizes->numelmasked]=Expec[i+c*CurrSizes->numelmasked]/sumed_all;
                }
            }
            else
            {
                for(long c=0; c<7; c++)
                {
                    Expec[i+c*CurrSizes->numelmasked]=1.0/7.0;
                }

            }
        }

    }
    GaussianFilter4D_cArray(Expec,S2L,L2S,1.0, CurrSizes);

    for(long i=0; i<(long)CurrSizes->numelmasked; i++)
    {
        for(long c=0; c<7; c++)
        {
            ShortPrior[i+c*CurrSizes->numelmasked]=Expec[i+c*CurrSizes->numelmasked];
            MRF[i+c*CurrSizes->numelmasked]=Expec[i+c*CurrSizes->numelmasked];
        }
    }

    delete [] Seed_Mask;
    delete [] SpeedFunc;
    delete [] wSulci;
    delete [] wGyri;
    return 1;

}

int Convert_WM_and_GM_to_PV(nifti_image * T1,
                            segPrecisionTYPE * BiasField,
                            segPrecisionTYPE * ShortPrior,
                            segPrecisionTYPE * Expec,
                            int * S2L,
                            segPrecisionTYPE * M,
                            segPrecisionTYPE * V,
                            ImageSize * CurrSizes)
{

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

    segPrecisionTYPE * T1ptrtmp = static_cast<segPrecisionTYPE *>(T1->data);
    for (long i=0; i<(long)CurrSizes->numelmasked; i++)
    {
        currentT1=T1ptrtmp[S2L[i]]+BiasField[i];

        averageFCtmp_WMGM=(currentT1-M[GMclass])/(M[WMclass]-M[GMclass]);

        if(averageFCtmp_WMGM<1 && averageFCtmp_WMGM>0)
        {
            averageFC_WMGM+=averageFCtmp_WMGM;
            count1++;
        }

        averageFCtmp_GMCSF=(currentT1-M[CSFclass])/(M[GMclass]-M[CSFclass]);

        if(averageFCtmp_GMCSF<1 && averageFCtmp_GMCSF>0)
        {
            averageFC_GMCSF+=averageFCtmp_GMCSF;
            count2++;
        }
    }

    averageFC_WMGM=(count1>0)?averageFC_WMGM/(float)(count1):0.5;
    averageFC_GMCSF=(count2>0)?averageFC_GMCSF/(float)(count2):0.5;


    cout<< "avr_WMGMfc=" << averageFC_WMGM << endl;
    cout<< "avr_GMCSFfc=" << averageFC_GMCSF << endl;

    for (long i=0; i<(long)CurrSizes->numelmasked; i++,WMindex++,GMindex++,CSFindex++,dGMindex++,iCSFindex++,WMGMpvindex++,GMCSFpvindex++)
    {
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

    return 1;
}



nifti_image * Copy_ShortExpec_to_Result(nifti_image * T1,
                                        segPrecisionTYPE * Expec,
                                        segPrecisionTYPE * BiasField,
                                        segPrecisionTYPE * BiasFieldCoefs,
                                        int * S2L,
                                        nifti_image * Priors,
                                        seg_EM_Params * segment_param,
                                        segPrecisionTYPE * M,
                                        ImageSize * CurrSizes)
{

    nifti_image * Result = nifti_copy_nim_info(Priors);
    if(segment_param->flag_PV_model)
    {

        segPrecisionTYPE * T1ptrtmp = static_cast<segPrecisionTYPE *>(T1->data);

        Result->dim[0]=4;
        Result->dim[4]=6;
        Result->datatype=DT_FLOAT32;
        Result->cal_max=1;
        nifti_set_filenames(Result,segment_param->filename_out,0,0);
        nifti_update_dims_from_array(Result);
        nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
        Result->data = (void *) calloc(Result->nvox, sizeof(segPrecisionTYPE));
        segPrecisionTYPE * Result_PTR = static_cast<segPrecisionTYPE *>(Result->data);
        for(unsigned int i=0; i<Result->nvox; i++)
        {
            Result_PTR[i]=0;
        }

        int Short_2_Long_Indices_tmp = 0;
        int class_nvox=Result->nx*Result->ny*Result->nz;
        segPrecisionTYPE * Resultdata= static_cast<segPrecisionTYPE *>(Result->data);
        segPrecisionTYPE * Expec_tmp = new segPrecisionTYPE [CurrSizes->numclass]();
        segPrecisionTYPE Expec_tmp_sum=0;
        segPrecisionTYPE Fractional_content=0;

        if(BiasField!=NULL)
        {
            for(long i=0; i<(long)CurrSizes->numelmasked; i++)
            {
                Short_2_Long_Indices_tmp=S2L[i];
                for(long currclass=0; currclass<CurrSizes->numclass; currclass++)
                {
                    Expec_tmp[currclass]=Expec[i+currclass*CurrSizes->numelmasked];
                }
                float classthreshold=0.1;
                Resultdata[Short_2_Long_Indices_tmp+(WMclass)*class_nvox]=Expec_tmp[WMclass]>classthreshold?Expec_tmp[WMclass]:0;
                Resultdata[Short_2_Long_Indices_tmp+(GMclass)*class_nvox]=Expec_tmp[GMclass]>classthreshold?Expec_tmp[GMclass]:0;
                Resultdata[Short_2_Long_Indices_tmp+(CSFclass)*class_nvox]=Expec_tmp[CSFclass]>classthreshold?Expec_tmp[CSFclass]:0;
                Resultdata[Short_2_Long_Indices_tmp+(dGMclass)*class_nvox]=Expec_tmp[dGMclass]>classthreshold?Expec_tmp[dGMclass]:0;
                Resultdata[Short_2_Long_Indices_tmp+(iCSFclass)*class_nvox]=Expec_tmp[iCSFclass]>classthreshold?Expec_tmp[iCSFclass]:0;


                // Estimating WM/GM fractional content from the T1 bias corrected data
                if(Expec_tmp[WMGMpvclass]>classthreshold)
                {
                    if(((Expec_tmp[WMclass])<(Expec_tmp[WMGMpvclass])) & (Expec_tmp[GMclass]<(Expec_tmp[WMGMpvclass])))
                    {
                        Expec_tmp_sum=Expec_tmp[WMclass]+Expec_tmp[WMGMpvclass]+Expec_tmp[GMclass];
                        Fractional_content=(M[WMclass]-(T1ptrtmp[Short_2_Long_Indices_tmp]+BiasField[i]))/(M[WMclass]-M[GMclass]);
                        Fractional_content=(Fractional_content<0)?0:Fractional_content;
                        Fractional_content=(Fractional_content>1)?1:Fractional_content;
                        Resultdata[Short_2_Long_Indices_tmp+(WMclass)*class_nvox]=(1-Fractional_content)/Expec_tmp_sum;
                        Resultdata[Short_2_Long_Indices_tmp+(GMclass)*class_nvox]=(Fractional_content)/Expec_tmp_sum;
                    }
                    else if(Expec_tmp[WMclass]>Expec_tmp[GMclass])
                    {
                        Resultdata[Short_2_Long_Indices_tmp+(WMclass)*class_nvox]=1;
                        Resultdata[Short_2_Long_Indices_tmp+(GMclass)*class_nvox]=0;
                    }
                    else
                    {
                        Resultdata[Short_2_Long_Indices_tmp+(WMclass)*class_nvox]=0;
                        Resultdata[Short_2_Long_Indices_tmp+(GMclass)*class_nvox]=1;
                    }
                }

                // Estimating CSF/GM fractional content from the T1 bias corrected data
                if(Expec_tmp[GMCSFpvclass]>classthreshold)
                {
                    if((Expec_tmp[CSFclass]<(Expec_tmp[GMCSFpvclass])) & (Expec_tmp[GMclass]<(Expec_tmp[GMCSFpvclass])))
                    {
                        Expec_tmp_sum=Expec_tmp[CSFclass]+Expec_tmp[GMCSFpvclass]+Expec_tmp[GMclass];
                        Fractional_content=(M[CSFclass]-(T1ptrtmp[Short_2_Long_Indices_tmp]+BiasField[i]))/(M[CSFclass]-M[GMclass]);
                        Fractional_content=(Fractional_content<0)?0:Fractional_content;
                        Fractional_content=(Fractional_content>1)?1:Fractional_content;
                        Resultdata[Short_2_Long_Indices_tmp+(CSFclass)*class_nvox]=(1-Fractional_content)/Expec_tmp_sum;
                        Resultdata[Short_2_Long_Indices_tmp+(GMclass)*class_nvox]=(Fractional_content)/Expec_tmp_sum;
                    }
                    else if(Expec_tmp[CSFclass]>Expec_tmp[GMclass])
                    {
                        Resultdata[Short_2_Long_Indices_tmp+(CSFclass)*class_nvox]=1;
                        Resultdata[Short_2_Long_Indices_tmp+(GMclass)*class_nvox]=0;
                    }
                    else
                    {
                        Resultdata[Short_2_Long_Indices_tmp+(CSFclass)*class_nvox]=0;
                        Resultdata[Short_2_Long_Indices_tmp+(GMclass)*class_nvox]=1;
                    }
                }

                // Normalizing the fractional contents
                Expec_tmp_sum=0;
                for(long currclass=0; currclass<nonPVNumClass; currclass++)
                {
                    Expec_tmp_sum+=Resultdata[Short_2_Long_Indices_tmp+(currclass)*class_nvox];
                }
                for(long currclass=0; currclass<nonPVNumClass; currclass++)
                {
                    Resultdata[Short_2_Long_Indices_tmp+(currclass)*class_nvox]=Resultdata[Short_2_Long_Indices_tmp+(currclass)*class_nvox]/Expec_tmp_sum;
                }

            }
        }
        else
        {

            for(long i=0; i<(long)CurrSizes->numelmasked; i++)
            {
                Short_2_Long_Indices_tmp=S2L[i];
                for(long currclass=0; currclass<CurrSizes->numclass; currclass++)
                {
                    Expec_tmp[currclass]=Expec[i+currclass*CurrSizes->numelmasked];
                }
                float classthreshold=0.1;
                Resultdata[Short_2_Long_Indices_tmp+(WMclass)*class_nvox]=Expec_tmp[WMclass]>classthreshold?Expec_tmp[WMclass]:0;
                Resultdata[Short_2_Long_Indices_tmp+(GMclass)*class_nvox]=Expec_tmp[GMclass]>classthreshold?Expec_tmp[GMclass]:0;
                Resultdata[Short_2_Long_Indices_tmp+(CSFclass)*class_nvox]=Expec_tmp[CSFclass]>classthreshold?Expec_tmp[CSFclass]:0;
                Resultdata[Short_2_Long_Indices_tmp+(dGMclass)*class_nvox]=Expec_tmp[dGMclass]>classthreshold?Expec_tmp[dGMclass]:0;
                Resultdata[Short_2_Long_Indices_tmp+(iCSFclass)*class_nvox]=Expec_tmp[iCSFclass]>classthreshold?Expec_tmp[iCSFclass]:0;


                // Estimating WM/GM fractional content from the T1 data
                if(Expec_tmp[WMGMpvclass]>classthreshold)
                {
                    if(((Expec_tmp[WMclass])<(Expec_tmp[WMGMpvclass])) & (Expec_tmp[GMclass]<(Expec_tmp[WMGMpvclass])))
                    {
                        Expec_tmp_sum=Expec_tmp[WMclass]+Expec_tmp[WMGMpvclass]+Expec_tmp[GMclass];
                        Fractional_content=(M[WMclass]-(T1ptrtmp[Short_2_Long_Indices_tmp]))/(M[WMclass]-M[GMclass]);
                        Fractional_content=(Fractional_content<0)?0:Fractional_content;
                        Fractional_content=(Fractional_content>1)?1:Fractional_content;
                        Resultdata[Short_2_Long_Indices_tmp+(WMclass)*class_nvox]=(1-Fractional_content)/Expec_tmp_sum;
                        Resultdata[Short_2_Long_Indices_tmp+(GMclass)*class_nvox]=(Fractional_content)/Expec_tmp_sum;
                    }
                    else if(Expec_tmp[WMclass]>Expec_tmp[GMclass])
                    {
                        Resultdata[Short_2_Long_Indices_tmp+(WMclass)*class_nvox]=1;
                        Resultdata[Short_2_Long_Indices_tmp+(GMclass)*class_nvox]=0;
                    }
                    else
                    {
                        Resultdata[Short_2_Long_Indices_tmp+(WMclass)*class_nvox]=0;
                        Resultdata[Short_2_Long_Indices_tmp+(GMclass)*class_nvox]=1;
                    }
                }

                // Estimating CSF/GM fractional content from the T1 data
                if(Expec_tmp[GMCSFpvclass]>classthreshold)
                {
                    if((Expec_tmp[CSFclass]<Expec_tmp[GMCSFpvclass]) & (Expec_tmp[GMclass]<Expec_tmp[GMCSFpvclass]))
                    {
                        Expec_tmp_sum=Expec_tmp[CSFclass]+Expec_tmp[GMCSFpvclass]+Expec_tmp[GMclass];
                        Fractional_content=(M[CSFclass]-(T1ptrtmp[Short_2_Long_Indices_tmp]))/(M[CSFclass]-M[GMclass]);
                        Fractional_content=(Fractional_content<0)?0:Fractional_content;
                        Fractional_content=(Fractional_content>1)?1:Fractional_content;
                        Resultdata[Short_2_Long_Indices_tmp+(CSFclass)*class_nvox]=(1-Fractional_content)/Expec_tmp_sum;
                        Resultdata[Short_2_Long_Indices_tmp+(GMclass)*class_nvox]=(Fractional_content)/Expec_tmp_sum;
                    }
                    else if(Expec_tmp[CSFclass]>Expec_tmp[GMclass])
                    {
                        Resultdata[Short_2_Long_Indices_tmp+(CSFclass)*class_nvox]=1;
                        Resultdata[Short_2_Long_Indices_tmp+(GMclass)*class_nvox]=0;
                    }
                    else
                    {
                        Resultdata[Short_2_Long_Indices_tmp+(CSFclass)*class_nvox]=0;
                        Resultdata[Short_2_Long_Indices_tmp+(GMclass)*class_nvox]=1;
                    }
                }

                // Normalizing the fractional contents
                Expec_tmp_sum=0;
                for(long currclass=0; currclass<nonPVNumClass; currclass++)
                {
                    Expec_tmp_sum+=Resultdata[Short_2_Long_Indices_tmp+(currclass)*class_nvox];
                }
                for(long currclass=0; currclass<nonPVNumClass; currclass++)
                {
                    Resultdata[Short_2_Long_Indices_tmp+(currclass)*class_nvox]=Resultdata[Short_2_Long_Indices_tmp+(currclass)*class_nvox]/Expec_tmp_sum;
                }



            }
            delete [] Expec_tmp;


        }
    }
    else
    {
        Result->dim[0]=4;
        Result->dim[4]=(CurrSizes->numclass);
        Result->datatype=DT_FLOAT32;
        Result->cal_max=1;
        nifti_set_filenames(Result,segment_param->filename_out,0,0);
        nifti_update_dims_from_array(Result);
        nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
        Result->data = (void *) calloc(Result->nvox, sizeof(segPrecisionTYPE));
        segPrecisionTYPE * Result_PTR = static_cast<segPrecisionTYPE *>(Result->data);
        for(unsigned int i=0; i<Result->nvox; i++)
        {
            Result_PTR[i]=0;
        }

        int * S2L_PRT = (int *) S2L;
        segPrecisionTYPE * Resultdata= static_cast<segPrecisionTYPE *>(Result->data);


        int class_nvox=CurrSizes->numel;
        for(long currclass=0; currclass<CurrSizes->numclass; currclass++)
        {

            segPrecisionTYPE * Resultdata_class = &Resultdata[(currclass)*class_nvox];
            segPrecisionTYPE * Expec_PTR = &Expec[(currclass)*CurrSizes->numelmasked];
            S2L_PRT= (int *) S2L;

            for(long i=0; i<(long)CurrSizes->numelmasked; i++,S2L_PRT++,Expec_PTR++)
            {
                Resultdata_class[*S2L_PRT]=*Expec_PTR;
            }
        }
    }
    return Result;
}

bool * binarise_image(segPrecisionTYPE * SingleImage,
                      segPrecisionTYPE Threshold,
                      ImageSize * CurrSizes)
{

    bool * Result = new bool [CurrSizes->numelmasked];

    for(long i=0; i<(long)CurrSizes->numelmasked; i++)
    {
        Result[i]=(SingleImage[i]>Threshold);
    }
    return Result;
}

nifti_image * Copy_Single_ShortImage_to_Result(segPrecisionTYPE * SingleImage,
                                               int * Short_2_Long_Indices,
                                               nifti_image * Sourceimage,
                                               char * filename,
                                               ImageSize * CurrSizes)
{

    nifti_image * Result = nifti_copy_nim_info(Sourceimage);
    nifti_set_filenames(Result,filename,0,0);
    Result->data = (void *) calloc(Result->nvox, sizeof(segPrecisionTYPE));
    segPrecisionTYPE * Result_PTR = static_cast<segPrecisionTYPE *>(Result->data);
    for(unsigned int i=0; i<Result->nvox; i++)
    {
        Result_PTR[i]=0;
    }

    int * Short_2_Long_Indices_PRT = (int *) Short_2_Long_Indices;
    segPrecisionTYPE * Resultdata= static_cast<segPrecisionTYPE *>(Result->data);

    segPrecisionTYPE * SingleImage_PTR =SingleImage;
    Short_2_Long_Indices_PRT= (int *) Short_2_Long_Indices;

    for(long i=0; i<(long)CurrSizes->numelmasked; i++,Short_2_Long_Indices_PRT++,SingleImage_PTR++)
    {
        Resultdata[*Short_2_Long_Indices_PRT]=(segPrecisionTYPE)(*SingleImage_PTR);
    }
    return Result;
}



nifti_image * Copy_BiasCorrected_to_Result_mask(segPrecisionTYPE * BiasField,
                                                int * Short_2_Long_Indices,
                                                nifti_image * T1,
                                                char * filename,
                                                ImageSize * CurrSizes)
{

    nifti_image * Result = nifti_copy_nim_info(T1);
    Result->dim[0]=4;
    Result->dim[4]=1;
    Result->datatype=DT_FLOAT32;
    Result->cal_max=1;
    nifti_set_filenames(Result,filename,0,0);
    nifti_update_dims_from_array(Result);
    nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
    Result->data = (void *) calloc(Result->nvox, sizeof(segPrecisionTYPE));
    segPrecisionTYPE * Resultdata = static_cast<segPrecisionTYPE *>(Result->data);
    segPrecisionTYPE * T1data = static_cast<segPrecisionTYPE *>(T1->data);
    for(unsigned int i=0; i<Result->nvox; i++)
    {
        Resultdata[i]=0;
    }

    int * Short_2_Long_Indices_PRT = (int *) Short_2_Long_Indices;

    Short_2_Long_Indices_PRT= (int *) Short_2_Long_Indices;
    segPrecisionTYPE to_resize=0;
    if(BiasField!=NULL)
    {
        for(long i=0; i<(long)CurrSizes->numelmasked; i++,Short_2_Long_Indices_PRT++)
        {
            to_resize=exp((BiasField[i]+T1data[*Short_2_Long_Indices_PRT])*0.693147181)-1;
            Resultdata[*Short_2_Long_Indices_PRT]=(to_resize*(CurrSizes->rescale_max[0]-CurrSizes->rescale_min[0])+CurrSizes->rescale_min[0]);
        }
    }
    else
    {
        for(long i=0; i<(long)CurrSizes->numelmasked; i++,Short_2_Long_Indices_PRT++)
        {
            to_resize=exp((T1data[*Short_2_Long_Indices_PRT])*0.693147181)-1;
            Resultdata[*Short_2_Long_Indices_PRT]=(to_resize*(CurrSizes->rescale_max[0]-CurrSizes->rescale_min[0])+CurrSizes->rescale_min[0]);
        }
    }
    return Result;
}

nifti_image * Copy_BiasCorrected_to_Result(segPrecisionTYPE * BiasField,
                                           nifti_image * T1,
                                           char * filename,
                                           ImageSize * CurrSizes)
{

    nifti_image * Result = nifti_copy_nim_info(T1);
    Result->dim[0]=4;
    Result->dim[4]=1;
    Result->datatype=DT_FLOAT32;
    Result->cal_max=1;
    nifti_set_filenames(Result,filename,0,0);
    nifti_update_dims_from_array(Result);
    nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
    Result->data = (void *) calloc(Result->nvox, sizeof(segPrecisionTYPE));
    segPrecisionTYPE * Resultdata = static_cast<segPrecisionTYPE *>(Result->data);
    segPrecisionTYPE * T1data = static_cast<segPrecisionTYPE *>(T1->data);
    for(unsigned int i=0; i<Result->nvox; i++)
    {
        Resultdata[i]=0;
    }
    segPrecisionTYPE to_resize=0;
    if(BiasField!=NULL)
    {
        for(long i=0; i<(long)CurrSizes->numel; i++)
        {
            to_resize=exp((BiasField[i]+T1data[i])*0.693147181)-1;
            Resultdata[i]=(to_resize*(CurrSizes->rescale_max[0]-CurrSizes->rescale_min[0])+CurrSizes->rescale_min[0]);
        }
    }
    else
    {
        for(long i=0; i<(long)CurrSizes->numel; i++)
        {
            to_resize=exp((T1data[i])*0.693147181)-1;
            Resultdata[i]=(to_resize*(CurrSizes->rescale_max[0]-CurrSizes->rescale_min[0])+CurrSizes->rescale_min[0]);
        }
    }
    return Result;
}



nifti_image * Copy_Expec_to_Result_Neonate_mask(segPrecisionTYPE * Expec,
                                                int * Short_2_Long_Indices,
                                                int * Long_2_Short_Indices,
                                                nifti_image * T1,
                                                float * BiasField,
                                                float * M,
                                                char * filename,
                                                ImageSize * CurrSizes)
{

    nifti_image * Result = nifti_copy_nim_info(T1);
    Result->dim[0]=4;
    Result->dim[4]=6;
    Result->dim[5]=1;
    Result->datatype=DT_FLOAT32;
    Result->cal_max=1;
    nifti_set_filenames(Result,filename,0,0);
    nifti_update_dims_from_array(Result);
    nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
    Result->data = (void *) calloc(Result->nvox, sizeof(segPrecisionTYPE));
    segPrecisionTYPE * Resultdata = static_cast<segPrecisionTYPE *>(Result->data);

    segPrecisionTYPE * T1ptrtmp = static_cast<segPrecisionTYPE *>(T1->data);
    for(unsigned int i=0; i<Result->nvox; i++)
    {
        Resultdata[i]=0;
    }

    int * Short_2_Long_Indices_PRT = (int *) Short_2_Long_Indices;

    int class_nvox=Result->nx*Result->ny*Result->nz;

    Short_2_Long_Indices_PRT= (int *) Short_2_Long_Indices;
    for(long currclass=0; currclass<6; currclass++)
    {
        segPrecisionTYPE * Resultdata_class = &Resultdata[(currclass)*class_nvox];
        segPrecisionTYPE * Expec_PTR = &Expec[(currclass)*CurrSizes->numelmasked];
        Short_2_Long_Indices_PRT= (int *) Short_2_Long_Indices;
        for(long i=0; i<(long)CurrSizes->numelmasked; i++,Short_2_Long_Indices_PRT++,Expec_PTR++)
        {
            Resultdata_class[(*Short_2_Long_Indices_PRT)]=(*Expec_PTR);
        }
    }


    int xyzpos[3]= {0};
    int distance=4;
    Short_2_Long_Indices_PRT= (int *) Short_2_Long_Indices;
    float Fractional_content=0;
    int biggestclass=-1;
    float biggestclass_prob=-1;
    int biggestclass2=-1;
    float biggestclass2_prob=-1;
    int neigh_index=0;

    for(long i=0; i<(long)CurrSizes->numelmasked; i++,Short_2_Long_Indices_PRT++)
    {
        if(Expec[(6)*CurrSizes->numelmasked+i]>0.1)
        {
            biggestclass=-1;
            biggestclass_prob=0.1;
            biggestclass2=-2;
            biggestclass2_prob=0.1;
            for(xyzpos[2]=-distance; xyzpos[2]<distance; xyzpos[2]++)
            {
                for(xyzpos[1]=-distance; xyzpos[1]<distance; xyzpos[1]++)
                {
                    for(xyzpos[0]=-distance; xyzpos[0]<distance; xyzpos[0]++)
                    {
                        neigh_index=Long_2_Short_Indices[(*Short_2_Long_Indices_PRT)+xyzpos[0]+xyzpos[1]*CurrSizes->xsize+xyzpos[2]*CurrSizes->xsize*CurrSizes->ysize];
                        if(neigh_index>0 && neigh_index<(int)CurrSizes->numel)
                        {
                            for(long neigh_class=1; neigh_class<6; neigh_class++)
                            {
                                if(Expec[neigh_index+(neigh_class)*CurrSizes->numelmasked]>biggestclass_prob)
                                {
                                    if(biggestclass!=biggestclass2)
                                    {
                                        biggestclass2_prob=biggestclass_prob;
                                        biggestclass2=biggestclass;
                                    }
                                    biggestclass_prob=Expec[neigh_index+(neigh_class)*CurrSizes->numelmasked];
                                    biggestclass=neigh_class;
                                }
                                if(Expec[neigh_index+(neigh_class)*CurrSizes->numelmasked]>biggestclass2_prob && Expec[neigh_index+(neigh_class)*CurrSizes->numelmasked]<biggestclass2_prob)
                                {
                                    if(biggestclass!=neigh_class)
                                    {
                                        biggestclass2_prob=Expec[neigh_index+(neigh_class)*CurrSizes->numelmasked];
                                        biggestclass2=neigh_class;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if(biggestclass>0 && biggestclass2>0)
            {
                Fractional_content=(M[biggestclass]-(T1ptrtmp[(*Short_2_Long_Indices_PRT)]+BiasField[i]))/(M[biggestclass]-M[biggestclass2]);
                Fractional_content=(Fractional_content<0)?0:Fractional_content;
                Fractional_content=(Fractional_content>1)?1:Fractional_content;
                Resultdata[(*Short_2_Long_Indices_PRT)+(biggestclass)*class_nvox]=(1-Fractional_content);
                Resultdata[(*Short_2_Long_Indices_PRT)+(biggestclass2)*class_nvox]=(Fractional_content);
                for(long otherclasses=0; otherclasses<6; otherclasses++)
                {
                    if(otherclasses!=biggestclass && otherclasses!=biggestclass2)
                    {
                        Resultdata[(*Short_2_Long_Indices_PRT)+(otherclasses)*class_nvox]=0;
                    }
                }
            }
            else
            {
                for(long otherclasses=0; otherclasses<6; otherclasses++)
                {
                    Resultdata[(*Short_2_Long_Indices_PRT)+(otherclasses)*class_nvox]=0;
                }
                Resultdata[(*Short_2_Long_Indices_PRT)+(2)*class_nvox]=1;
            }
        }
    }


    Short_2_Long_Indices_PRT= (int *) Short_2_Long_Indices;
    float sumexp=0;

    for(long i=0; i<(long)CurrSizes->numel; i++)
    {
        if(Long_2_Short_Indices[i]>0)
        {
            sumexp=0;
            for(long currclass=0; currclass<6; currclass++)
            {
                if((Resultdata[i+currclass*CurrSizes->numel])>0.01)
                {
                    sumexp+=Resultdata[i+currclass*CurrSizes->numel];
                }
            }
            if(sumexp>0)
            {
                for(long currclass=0; currclass<6; currclass++)
                {
                    if((Resultdata[i+currclass*CurrSizes->numel])>0.01)
                    {
                        Resultdata[i+currclass*CurrSizes->numel]=Resultdata[i+currclass*CurrSizes->numel]/sumexp;
                    }
                    else
                    {
                        Resultdata[i+currclass*CurrSizes->numel]=0;
                    }
                }
            }
        }
    }

    return Result;
}

nifti_image * Copy_Expec_to_Result_mask(segPrecisionTYPE * Expec,
                                        int * Short_2_Long_Indices,
                                        nifti_image * T1,
                                        char * filename,
                                        ImageSize * CurrSizes)
{

    nifti_image * Result = nifti_copy_nim_info(T1);
    Result->dim[0]=4;
    Result->dim[4]=CurrSizes->numclass;
    Result->dim[5]=1;
    Result->scl_inter=0;
    Result->scl_slope=1;
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

    int * Short_2_Long_Indices_PRT = (int *) Short_2_Long_Indices;

    int class_nvox=Result->nx*Result->ny*Result->nz;

    Short_2_Long_Indices_PRT= (int *) Short_2_Long_Indices;
    for(long currclass=0; currclass<CurrSizes->numclass; currclass++)
    {

        segPrecisionTYPE * Resultdata_class = &Resultdata[(currclass)*class_nvox];
        segPrecisionTYPE * Expec_PTR = &Expec[(currclass)*CurrSizes->numelmasked];
        Short_2_Long_Indices_PRT= (int *) Short_2_Long_Indices;

        for(long i=0; i<(long)CurrSizes->numelmasked; i++,Short_2_Long_Indices_PRT++,Expec_PTR++)
        {
            Resultdata_class[Short_2_Long_Indices[i]]=*Expec_PTR;
        }
    }
    return Result;
}


nifti_image * Copy_Expec_to_Result(segPrecisionTYPE * Expec,
                                   nifti_image * T1,
                                   char * filename,
                                   ImageSize * CurrSizes)
{

    nifti_image * Result = nifti_copy_nim_info(T1);
    Result->dim[0]=4;
    Result->dim[4]=CurrSizes->numclass;
    Result->datatype=DT_FLOAT32;
    Result->cal_max=1;
    Result->scl_inter=0;
    Result->scl_slope=1;
    nifti_set_filenames(Result,filename,0,0);
    nifti_update_dims_from_array(Result);
    nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
    Result->data = (void *) calloc(Result->nvox, sizeof(segPrecisionTYPE));
    segPrecisionTYPE * Resultdata = static_cast<segPrecisionTYPE *>(Result->data);
    for(unsigned int i=0; i<Result->nvox; i++)
    {
        Resultdata[i]=0;
    }
    int class_nvox=Result->nx*Result->ny*Result->nz;

    for(long currclass=0; currclass<CurrSizes->numclass; currclass++)
    {

        segPrecisionTYPE * Resultdata_class = &Resultdata[(currclass)*class_nvox];
        segPrecisionTYPE * Expec_PTR = &Expec[(currclass)*CurrSizes->numel];

        for(long i=0; i<(long)CurrSizes->numel; i++,Expec_PTR++)
        {
            Resultdata_class[i]=*Expec_PTR;
        }
    }
    return Result;
}

nifti_image * Get_Bias_Corrected(float * BiasField,
                                 nifti_image * T1,
                                 char * filename,
                                 ImageSize * CurrSizes)
{

    nifti_image * BiasCorrected = nifti_copy_nim_info(T1);
    BiasCorrected->dim[0]=4;
    BiasCorrected->dim[4]=CurrSizes->usize;
    BiasCorrected->datatype=DT_FLOAT32;
    BiasCorrected->cal_max=(CurrSizes->rescale_max[0]);
    BiasCorrected->scl_inter=0;
    BiasCorrected->scl_slope=1;
    segPrecisionTYPE * T1data = static_cast<segPrecisionTYPE *>(T1->data);



    nifti_set_filenames(BiasCorrected,filename,0,0);
    nifti_update_dims_from_array(BiasCorrected);
    nifti_datatype_sizes(BiasCorrected->datatype,&BiasCorrected->nbyper,&BiasCorrected->swapsize);
    BiasCorrected->data = (void *) calloc(BiasCorrected->nvox, sizeof(segPrecisionTYPE));
    segPrecisionTYPE * BiasCorrected_PTR = static_cast<segPrecisionTYPE *>(BiasCorrected->data);

    for(long multispec=0; multispec<CurrSizes->usize; multispec++)
    {

        BiasCorrected_PTR = static_cast<segPrecisionTYPE *>(BiasCorrected->data);
        BiasCorrected_PTR = &BiasCorrected_PTR[multispec*BiasCorrected->nvox];

        T1data = static_cast<segPrecisionTYPE *>(T1->data);
        T1data = &T1data[multispec*BiasCorrected->nvox];

        for( unsigned int i=0; i<BiasCorrected->nvox; i++)
        {
            BiasCorrected_PTR[i]=0;
        }

        float to_resize=0;

        for(long i=0; i<(long)CurrSizes->numel; i++)
        {
            to_resize=exp((BiasField[i]+T1data[i])*0.693147181)-1;
            BiasCorrected_PTR[i]=(to_resize*(CurrSizes->rescale_max[multispec]-CurrSizes->rescale_min[multispec])+CurrSizes->rescale_min[multispec]);
        }
    }
    return BiasCorrected;
}


nifti_image * Get_Bias_Corrected_mask(float * BiasFieldCoefs,
                                      nifti_image * T1,
                                      nifti_image * Mask,
                                      char * filename,
                                      ImageSize * CurrSizes,
                                      int biasOrder)
{

    int UsedBasisFunctions=(int)((biasOrder+1) * (biasOrder+2)/2 *(biasOrder+3)/3);
    segPrecisionTYPE * T1data = static_cast<segPrecisionTYPE *>(T1->data);

    nifti_image * BiasCorrected = nifti_copy_nim_info(T1);
    BiasCorrected->dim[0]=4;
    BiasCorrected->dim[4]=CurrSizes->usize;
    BiasCorrected->datatype=DT_FLOAT32;
    BiasCorrected->cal_max=(CurrSizes->rescale_max[0]);
    BiasCorrected->scl_inter=0;
    BiasCorrected->scl_slope=1;

    float * brainmask= new float [CurrSizes->numel];
    for(long i=0; i<(long)CurrSizes->numel; i++)
    {

        brainmask[i]=(T1data[i]!=T1data[i])?0:T1data[i];
    }
    otsu(brainmask,NULL,CurrSizes);
    Dillate(brainmask,5,CurrSizes);
    Erosion(brainmask,4,CurrSizes);
    bool* Maskptr = static_cast<bool * >(Mask->data);
    for(long i=0; i<(long)CurrSizes->numel; i++)
    {
        brainmask[i]*=Maskptr[i];
    }
    GaussianFilter4D_cArray(brainmask, 3.0f, CurrSizes);


    nifti_set_filenames(BiasCorrected,filename,0,0);
    nifti_update_dims_from_array(BiasCorrected);
    nifti_datatype_sizes(BiasCorrected->datatype,&BiasCorrected->nbyper,&BiasCorrected->swapsize);
    BiasCorrected->data = (void *) calloc(BiasCorrected->nvox, sizeof(segPrecisionTYPE));
    segPrecisionTYPE * BiasCorrected_PTR = static_cast<segPrecisionTYPE *>(BiasCorrected->data);

    float BiasField=0;
    segPrecisionTYPE currxpower[maxAllowedBCPowerOrder];
    segPrecisionTYPE currypower[maxAllowedBCPowerOrder];
    segPrecisionTYPE currzpower[maxAllowedBCPowerOrder];
    float xpos=0.0f;
    float ypos=0.0f;
    float zpos=0.0f;
    segPrecisionTYPE not_point_five_times_dims_x=(0.5f*(segPrecisionTYPE)CurrSizes->xsize);
    segPrecisionTYPE not_point_five_times_dims_y=(0.5f*(segPrecisionTYPE)CurrSizes->ysize);
    segPrecisionTYPE not_point_five_times_dims_z=(0.5f*(segPrecisionTYPE)CurrSizes->zsize);
    segPrecisionTYPE inv_not_point_five_times_dims_x=1.0f/(0.5f*(segPrecisionTYPE)CurrSizes->xsize);
    segPrecisionTYPE inv_not_point_five_times_dims_y=1.0f/(0.5f*(segPrecisionTYPE)CurrSizes->ysize);
    segPrecisionTYPE inv_not_point_five_times_dims_z=1.0f/(0.5f*(segPrecisionTYPE)CurrSizes->zsize);
    int ind=0;

    for(long multispec=0; multispec<CurrSizes->usize; multispec++)
    {

        BiasCorrected_PTR = static_cast<segPrecisionTYPE *>(BiasCorrected->data);
        BiasCorrected_PTR = &BiasCorrected_PTR[multispec*CurrSizes->numel];
        T1data = static_cast<segPrecisionTYPE *>(T1->data);
        T1data = &T1data[multispec*CurrSizes->numel];

        float * BiasFieldCoefs_multispec = &BiasFieldCoefs[multispec*UsedBasisFunctions];


        for(long i=0; i<(long)CurrSizes->numel; i++)
        {
            BiasCorrected_PTR[i]=0;
        }

        float to_resize=0;
        int index_full=0;
        for (int iz=0; iz<CurrSizes->zsize; iz++)
        {
            for (int iy=0; iy<CurrSizes->ysize; iy++)
            {
                for (int ix=0; ix<CurrSizes->xsize; ix++)
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

                    for (order=1; order<maxAllowedBCPowerOrder; order++, orderminusone++)
                    {
                        currxpower[order]=currxpower[orderminusone]*xpos;
                        currypower[order]=currypower[orderminusone]*ypos;
                        currzpower[order]=currzpower[orderminusone]*zpos;
                    }

                    // Estimate the basis
                    ind=0;
                    for(long order=0; order<=biasOrder; order++)
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

                    to_resize=exp((BiasField+T1data[index_full])*0.693147181)-1;
                    BiasCorrected_PTR[index_full]=(to_resize*(CurrSizes->rescale_max[multispec]-CurrSizes->rescale_min[multispec])+CurrSizes->rescale_min[multispec]);
                    index_full++;
                }
            }
        }
    }
    return BiasCorrected;
}

int printloglik(int iter,
                segPrecisionTYPE loglik,
                segPrecisionTYPE oldloglik)
{
    if(iter>0)
    {
        if ((loglik-oldloglik)/fabsf(oldloglik)>0 && (loglik-oldloglik)/fabsf(oldloglik)<100)
        {
            cout<< "Loglik = " << setprecision(7)<<loglik << " : Ratio = " << (loglik-oldloglik)/fabsf(oldloglik) << endl;
        }
        else
        {
            cout<< "Loglik = " << setprecision(7)<<loglik << endl;
        }
    }
    else
    {
        cout<< "Initial Loglik = " << setprecision(7) << loglik << endl;
    }
    return 1;
}



void MRFregularization_mask(const segPrecisionTYPE * Expec,
                            const segPrecisionTYPE * G,
                            const segPrecisionTYPE * H,
                            segPrecisionTYPE * MRFbeta,
                            segPrecisionTYPE * MRFprior,
                            segPrecisionTYPE * AtlasPrior,
                            int * Long_2_Short_Indices,
                            int * Short_2_Long_Indices,
                            ImageSize * CurrSizes,
                            bool MRFflag,
                            int verbose_level)
{
    int numelmasked=CurrSizes->numelmasked;
    int numclass=CurrSizes->numclass;

    if(MRFflag)
    {
        segPrecisionTYPE * MRFpriorPtr = (segPrecisionTYPE *)MRFprior;
        int * Long_2_Short_IndicesPtr = (int *)Long_2_Short_Indices;
        int col_size, plane_size;
        int maxiy, maxix, maxiz;
        col_size = (int)(CurrSizes->xsize);
        plane_size = (int)(CurrSizes->xsize)*(CurrSizes->ysize);

        maxix = (int)(CurrSizes->xsize);
        maxiy = (int)(CurrSizes->ysize);
        maxiz = (int)(CurrSizes->zsize);

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
                            Gplane[currclass]+=Expec[indexWest];
                            Gplane[currclass]+=Expec[indexEast];
                            Gplane[currclass]+=Expec[indexNorth];
                            Gplane[currclass]+=Expec[indexSouth];
                            Hplane[currclass]+=Expec[indexTop];
                            Hplane[currclass]+=Expec[indexBottom];
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
                            if(MRFbeta==NULL)
                            {
                                Temp_MRF_Class_Expect[currclass] = exp(Temp_MRF_Class_Expect[currclass])*AtlasPrior[curr_short_centreindex+numelmasked_currclass_shift[currclass]];
                            }
                            else
                            {
                                Temp_MRF_Class_Expect[currclass] = exp(MRFbeta[curr_short_centreindex]*Temp_MRF_Class_Expect[currclass])*AtlasPrior[curr_short_centreindex+numelmasked_currclass_shift[currclass]];
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
                MRFprior[i+currclass*numelmasked]=AtlasPrior[i+currclass*numelmasked];
            }
        }
    }
}


void MRFregularization(const segPrecisionTYPE * Expec,
                       const segPrecisionTYPE * G,
                       const segPrecisionTYPE * H,
                       segPrecisionTYPE * MRFbeta,
                       segPrecisionTYPE * MRFprior,
                       segPrecisionTYPE * AtlasPrior,
                       ImageSize * CurrSizes,
                       bool MRFflag,
                       int verbose_level)
{

    int numel=CurrSizes->numel;
    int numclass=CurrSizes->numclass;

    if(MRFflag)
    {
        segPrecisionTYPE * MRFpriorPtr = (segPrecisionTYPE *)MRFprior;
        int col_size, plane_size, indexCentre, indexWest, indexEast, indexSouth, indexNorth, indexTop, indexBottom;
        int ix, iy, iz,maxiy, maxix, maxiz, neighbourclass;
        segPrecisionTYPE Sum_Temp_MRF_Class_Expect;
        col_size = (int)(CurrSizes->xsize);
        plane_size = (int)(CurrSizes->xsize)*(CurrSizes->ysize);

        maxix = (int)(CurrSizes->xsize);
        maxiy = (int)(CurrSizes->ysize);
        maxiz = (int)(CurrSizes->zsize);
        segPrecisionTYPE Gplane[maxNumbClass];
        segPrecisionTYPE Hplane[maxNumbClass];
        segPrecisionTYPE Temp_MRF_Class_Expect[maxNumbClass];
        if(verbose_level>0)
        {
            cout << "Optimising MRF"<<endl;
            flush(cout);
        }
        register int currclass;

        unsigned int numel_currclass_shift[maxNumbClass];
        //unsigned int image_size_currclass_shift[max_numbclass];
        for(int i=0; i<numclass; i++)
        {
            numel_currclass_shift[i]=i*numel;
        }
        indexCentre=0;
        for (iz=0; iz<maxiz-0; iz++)
        {
            for (iy=0; iy<maxiy-0; iy++)
            {
                for (ix=0; ix<maxix-0; ix++)
                {
                    Sum_Temp_MRF_Class_Expect = 0;
                    indexWest=(indexCentre-col_size)>=0?(indexCentre-col_size):-1;
                    indexEast=(indexCentre+col_size)<numel?(indexCentre+col_size):-1;
                    indexNorth=(indexCentre-1)>=0?(indexCentre-1):-1;
                    indexSouth=(indexCentre+1)<numel?(indexCentre+1):-1;
                    indexBottom=(indexCentre-plane_size)>=0?(indexCentre-plane_size):-1;
                    indexTop=(indexCentre+plane_size)<numel?(indexCentre+plane_size):-1;
                    for (currclass=0; currclass<numclass; currclass++)
                    {
                        Gplane[currclass] = 0.0;
                        Hplane[currclass] = 0.0;
                        Temp_MRF_Class_Expect[currclass] = 0.0;
                        Gplane[currclass]+=((indexWest>=0)?Expec[indexWest]:0);
                        Gplane[currclass]+=((indexEast>=0)?Expec[indexEast]:0);
                        Gplane[currclass]+=((indexNorth>=0)?Expec[indexNorth]:0);
                        Gplane[currclass]+=((indexSouth>=0)?Expec[indexSouth]:0);
                        Hplane[currclass]+=((indexTop>=0)?Expec[indexTop]:0);
                        Hplane[currclass]+=((indexBottom>=0)?Expec[indexBottom]:0);
                        if(currclass<numclass)
                        {
                            indexWest+=(indexWest>=0)?numel:0;
                            indexEast+=(indexEast>=0)?numel:0;
                            indexNorth+=(indexNorth>=0)?numel:0;
                            indexSouth+=(indexSouth>=0)?numel:0;
                            indexTop+=(indexTop>=0)?numel:0;
                            indexBottom+=(indexBottom>=0)?numel:0;
                        }
                    }
                    for (currclass=0; currclass<numclass; currclass++)
                    {
                        for (neighbourclass=0; neighbourclass<numclass; neighbourclass++)
                        {
                            Temp_MRF_Class_Expect[currclass]-=G[currclass+(numclass)*neighbourclass]*Gplane[neighbourclass]+H[currclass+(numclass)*neighbourclass]*Hplane[neighbourclass];
                        }
                        Temp_MRF_Class_Expect[currclass] = exp(Temp_MRF_Class_Expect[currclass]) * AtlasPrior[indexCentre+numel_currclass_shift[currclass]];
                        Sum_Temp_MRF_Class_Expect += Temp_MRF_Class_Expect[currclass];
                    }
                    for (currclass=0; currclass<numclass; currclass++)
                    {
                        MRFpriorPtr[indexCentre+numel_currclass_shift[currclass]]=(Temp_MRF_Class_Expect[currclass]/Sum_Temp_MRF_Class_Expect);
                    }
                    indexCentre++;
                }
            }
        }
    }
    else
    {
        for(int currclass=0; currclass<numclass; currclass++)
        {
            for(int i=0; i<(numel); i++)
            {
                MRFprior[i+currclass*numel]=AtlasPrior[i+currclass*numel];
            }
        }

    }

}




void MRFregularization_mask2D(const segPrecisionTYPE * Expec,
                              const segPrecisionTYPE * G,
                              const segPrecisionTYPE * H,
                              segPrecisionTYPE * MRFbeta,
                              segPrecisionTYPE * MRFprior,
                              segPrecisionTYPE * AtlasPrior,
                              int * Long_2_Short_Indices,
                              int * Short_2_Long_Indices,
                              ImageSize * CurrSizes,
                              bool MRFflag,
                              int verbose_level)
{

    int numelmasked=CurrSizes->numelmasked;
    int numclass=CurrSizes->numclass;
    if(MRFflag)
    {
        segPrecisionTYPE * MRFpriorPtr = (segPrecisionTYPE *)MRFprior;
        int * Long_2_Short_IndicesPtr = (int *)Long_2_Short_Indices;
        int col_size, indexCentre, indexWest, indexEast, indexSouth, indexNorth;
        int ix, iy,maxiy, maxix, neighbourclass;
        segPrecisionTYPE Sum_Temp_MRF_Class_Expect;
        col_size = (int)(CurrSizes->xsize);
        maxix = (int)(CurrSizes->xsize);
        maxiy = (int)(CurrSizes->ysize);
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
                        Gplane[currclass]+=Expec[indexWest];
                        Gplane[currclass]+=Expec[indexEast];
                        Gplane[currclass]+=Expec[indexNorth];
                        Gplane[currclass]+=Expec[indexSouth];
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
                        if(MRFbeta==NULL)
                        {
                            Temp_MRF_Class_Expect[currclass] = exp(Temp_MRF_Class_Expect[currclass])*AtlasPrior[curr_short_centreindex+numelmasked_currclass_shift[currclass]];
                        }
                        else
                        {
                            Temp_MRF_Class_Expect[currclass] = exp(MRFbeta[curr_short_centreindex]*Temp_MRF_Class_Expect[currclass])*AtlasPrior[curr_short_centreindex+numelmasked_currclass_shift[currclass]];
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
                MRFprior[i+currclass*numelmasked]=AtlasPrior[i+currclass*numelmasked];
            }
        }

    }

}


void MRFregularization2D(const segPrecisionTYPE * Expec,
                         const segPrecisionTYPE * G,
                         const segPrecisionTYPE * H,
                         segPrecisionTYPE * MRFbeta,
                         segPrecisionTYPE * MRFprior,
                         segPrecisionTYPE * AtlasPrior,
                         ImageSize * CurrSizes,
                         bool MRFflag,
                         int verbose_level)
{

    int numel=CurrSizes->numel;
    int numclass=CurrSizes->numclass;
    if(MRFflag)
    {
        segPrecisionTYPE * MRFpriorPtr = (segPrecisionTYPE *)MRFprior;
        int col_size, indexCentre, indexWest, indexEast, indexSouth, indexNorth;
        int ix, iy,maxiy, maxix, neighbourclass;
        segPrecisionTYPE Sum_Temp_MRF_Class_Expect;
        col_size = (int)(CurrSizes->xsize);
        maxix = (int)(CurrSizes->xsize);
        maxiy = (int)(CurrSizes->ysize);
        segPrecisionTYPE Gplane[maxNumbClass];
        segPrecisionTYPE Temp_MRF_Class_Expect[maxNumbClass];
        if(verbose_level>0)
        {
            cout << "Optimising MRF"<<endl;
            flush(cout);
        }
        register int currclass;

        unsigned int numel_currclass_shift[maxNumbClass];
        for(int i=0; i<numclass; i++)
        {
            numel_currclass_shift[i]=i*numel;
        }
        indexCentre=0;
        for (iy=0; iy<maxiy-0; iy++)
        {
            for (ix=0; ix<maxix-0; ix++)
            {
                Sum_Temp_MRF_Class_Expect = 0;
                indexWest=(indexCentre-col_size)>=0?(indexCentre-col_size):-1;
                indexEast=(indexCentre+col_size)<numel?(indexCentre+col_size):-1;
                indexNorth=(indexCentre-1)>=0?(indexCentre-1):-1;
                indexSouth=(indexCentre+1)<numel?(indexCentre+1):-1;
                for (currclass=0; currclass<numclass; currclass++)
                {
                    Gplane[currclass] = 0.0;
                    Temp_MRF_Class_Expect[currclass] = 0.0;
                    Gplane[currclass]+=((indexWest>=0)?Expec[indexWest]:0);
                    Gplane[currclass]+=((indexEast>=0)?Expec[indexEast]:0);
                    Gplane[currclass]+=((indexNorth>=0)?Expec[indexNorth]:0);
                    Gplane[currclass]+=((indexSouth>=0)?Expec[indexSouth]:0);
                    if(currclass<numclass)
                    {
                        indexWest+=(indexWest>=0)?numel:0;
                        indexEast+=(indexEast>=0)?numel:0;
                        indexNorth+=(indexNorth>=0)?numel:0;
                        indexSouth+=(indexSouth>=0)?numel:0;
                    }
                }
                for (currclass=0; currclass<numclass; currclass++)
                {
                    for (neighbourclass=0; neighbourclass<numclass; neighbourclass++)
                    {
                        Temp_MRF_Class_Expect[currclass]-=G[currclass+(numclass)*neighbourclass]*Gplane[neighbourclass];
                    }
                    Temp_MRF_Class_Expect[currclass] = exp(Temp_MRF_Class_Expect[currclass]) * AtlasPrior[indexCentre+numel_currclass_shift[currclass]];
                    Sum_Temp_MRF_Class_Expect += Temp_MRF_Class_Expect[currclass];
                }
                for (currclass=0; currclass<numclass; currclass++)
                {
                    MRFpriorPtr[indexCentre+numel_currclass_shift[currclass]]=(Temp_MRF_Class_Expect[currclass]/Sum_Temp_MRF_Class_Expect);
                }
                indexCentre++;
            }
        }
    }
    else
    {
        for(int currclass=0; currclass<numclass; currclass++)
        {
            for(int i=0; i<(numel); i++)
            {
                MRFprior[i+currclass*numel]=AtlasPrior[i+currclass*numel];
            }
        }

    }

}


void BiasCorrection(segPrecisionTYPE * BiasField,
                    segPrecisionTYPE * BiasFieldCoefs,
                    nifti_image * T1,
                    segPrecisionTYPE * Expec,
                    segPrecisionTYPE * Outlierness,
                    segPrecisionTYPE * M,
                    segPrecisionTYPE * V,
                    int biasOrder,
                    ImageSize * CurrSizes,
                    bool flag_Bias,
                    int verbose_level)
{

    if(verbose_level>0)
    {
        cout << "Optimising the Bias Field with order " << biasOrder<< endl;
        if(CurrSizes->usize>1)
        {
            cout<< "Assuming fully decoupled bias-fields" << endl;
        }
        flush(cout);
    }
    int reduxfactor=reduxFactorForBias;
    int nrOfClasses = CurrSizes->numclass;
    //nrOfClasses = 1;
    //int nrOfClasses = non_PV_numclass;
    int TotalLength = CurrSizes->numel;
    int UsedBasisFunctions=(int)((biasOrder+1) * (biasOrder+2)/2 *(biasOrder+3)/3);
    segPrecisionTYPE * sampledData = static_cast<segPrecisionTYPE *>(T1->data);


    // Precompute Powers depending on the current BiasOrder
    int PowerOrder [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3))]= {0};
    int ind=0;
    for(int order=0; order<=biasOrder; order++)
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
    for(int i=0; i<nrOfClasses; i++)
    {
        invV[i]=0;
        currM[i]=0;
    }

    for(long multispec=0; multispec<CurrSizes->usize; multispec++)
    {
        sampledData = &sampledData[multispec*CurrSizes->numel];
        // Precompute the M and V  inverses

        for(int i=0; i<nrOfClasses; i++)
        {
            invV[i]=1.0f/(V[i*CurrSizes->usize*CurrSizes->usize+multispec+multispec*CurrSizes->usize]);
            currM[i]=M[i*CurrSizes->usize+multispec];
        }



        segPrecisionTYPE A [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)*((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)]= {0.0f};
        segPrecisionTYPE B [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)]= {0.0f};
        segPrecisionTYPE C [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)]= {0.0f};

//#ifdef _OPENMP
//        SegPrecisionTYPE * Athread =new SegPrecisionTYPE [omp_get_max_threads()*((maxallowedpowerorder+1)*(maxallowedpowerorder+2)/2*(maxallowedpowerorder+3)/3)*((maxallowedpowerorder+1)*(maxallowedpowerorder+2)/2*(maxallowedpowerorder+3)/3)]();
//#endif
// Precompute sizes
        int col_size = (int)(CurrSizes->xsize);
        int plane_size = (int)(CurrSizes->xsize)*(CurrSizes->ysize);
        int maxix = (int)(CurrSizes->xsize);
        int maxiy = (int)(CurrSizes->ysize);
        int maxiz = (int)(CurrSizes->zsize);
        int Dims3d[3]= {0};
        Dims3d[0]=maxix;
        Dims3d[1]=maxiy;
        Dims3d[2]=maxiz;

// Precompute number of samples as it was never computed
        int samplecount=0;
        int linearindexes=0;
        int currindex=0;
        if(CurrSizes->numelbias==0)
        {
            for (int iz=0; iz<maxiz; iz+=reduxfactor)
            {
                for (int iy=0; iy<maxiy; iy+=reduxfactor)
                {
                    for (int ix=0; ix<maxix; ix+=reduxfactor)
                    {
                        samplecount++;
                    }
                }
            }
            if(verbose_level>0)
            {
                cout << "Samplecount = " << samplecount<<"\n";
                flush(cout);
            }
            CurrSizes->numelbias=samplecount;
        }
        else
        {
            samplecount=CurrSizes->numelbias;
        }

// CALC MATRIX A

// Calc W (Van Leemput 1999 eq 7)

        segPrecisionTYPE * Tempvar= new segPrecisionTYPE [samplecount] ();
        segPrecisionTYPE Tempvar_tmp=0;
        currindex=0;
        int tempvarindex=0;

        for (int iz=0; iz<maxiz; iz+=reduxfactor)
        {
            for (int iy=0; iy<maxiy; iy+=reduxfactor)
            {
                for (int ix=0; ix<maxix; ix+=reduxfactor)
                {
                    currindex=iz*plane_size+iy*col_size+ix;
                    Tempvar_tmp=0;
                    for(int j=0; j<nrOfClasses; j++)
                    {
                        Tempvar_tmp+=Expec[currindex+TotalLength*j]*invV[j];
                    }
                    Tempvar[tempvarindex]=Tempvar_tmp;
                    tempvarindex++;
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
//#ifdef _OPENMP
//#pragma omp parallel
//#endif
        for (int iz=0; iz<maxiz; iz+=reduxfactor)
        {
            for (int iy=0; iy<maxiy; iy+=reduxfactor)
            {
                for (int ix=0; ix<maxix; ix+=reduxfactor)
                {
                    linearindexes=(iz)*(CurrSizes->xsize)*(CurrSizes->ysize)+(iy)*(CurrSizes->xsize)+ix;

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
                        // Because Powerorder is alwais int, use a special power function (faster)
                        *Basisptr1=(pow_int(xpos,PowerOrder[x_bias_index_shift])*pow_int(ypos,PowerOrder[y_bias_index_shift])*pow_int(zpos,PowerOrder[z_bias_index_shift]));
                        // Instead, although slower, one can use
                        // TmpA=(pow(xpos,PowerOrder[0+j2*3])*pow(ypos,PowerOrder[1+j2*3])*pow(zpos,PowerOrder[2+j2*3]));
                    }
                    Basisptr1= (segPrecisionTYPE *) Basis;
                    //#ifdef _OPENMP
                    //                    Aptr=&Athread[omp_get_thread_num()*((maxallowedpowerorder+1)*(maxallowedpowerorder+2)/2*(maxallowedpowerorder+3)/3)*((maxallowedpowerorder+1)*(maxallowedpowerorder+2)/2*(maxallowedpowerorder+3)/3)];
                    //                    for(int j2=0; j2<UsedBasisFunctions; j2++, Basisptr1++){
                    //                        Basisptr2= &Basis[j2];
                    //                        SegPrecisionTYPE * Aptr2= &Aptr[j2+j2*UsedBasisFunctions];
                    //                        for(int i2=j2; i2<UsedBasisFunctions; i2++, Aptr2++, Basisptr2++){
                    //                            (*Aptr2)+=(*Basisptr2)*(current_Tempvar)*(*Basisptr1);
                    //                        }
                    //                    }
                    //#else
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
                    linearindexes=(iz)*(CurrSizes->xsize)*(CurrSizes->ysize)+(iy)*(CurrSizes->xsize)+ix;

                    Wi=0;
                    Wij=0;
                    Yest=0;
                    Ysum=0;
                    for(int j=0; j<nrOfClasses; j++)
                    {
                        segPrecisionTYPE tmpexpec = (segPrecisionTYPE)Expec[linearindexes+TotalLength*j];
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
            for (int iz=0; iz<maxiz; iz+=reduxfactor)
            {
                for (int iy=0; iy<maxiy; iy+=reduxfactor)
                {
                    for (int ix=0; ix<maxix; ix+=reduxfactor)
                    {
                        linearindexes=(iz)*(CurrSizes->xsize)*(CurrSizes->ysize)+(iy)*(CurrSizes->xsize)+ix;
                        B[i2]+=pow_int((((segPrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x),PowerOrder[0+i2*3])*
                               pow_int((((segPrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y),PowerOrder[1+i2*3])*
                               pow_int((((segPrecisionTYPE)iz-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z),PowerOrder[2+i2*3])*
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



        for (int iz=0; iz<maxiz; iz++)
        {
            for (int iy=0; iy<maxiy; iy++)
            {
                for (int ix=0; ix<maxix; ix++)
                {
                    segPrecisionTYPE tmpbiasfield=0.0f;
                    segPrecisionTYPE currxpower[maxAllowedBCPowerOrder]={0};
                    segPrecisionTYPE currypower[maxAllowedBCPowerOrder]={0};
                    segPrecisionTYPE currzpower[maxAllowedBCPowerOrder]={0};
                    linearindexes=(iz)*(CurrSizes->xsize)*(CurrSizes->ysize)+(iy)*(CurrSizes->xsize)+ix;
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
                    int maxorderplusone=biasOrder+1;
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
                    for(int order=0; order<=biasOrder; order++)
                    {
                        for(int xorder=0; xorder<=order; xorder++)
                        {
                            for(int yorder=0; yorder<=(order-xorder); yorder++)
                            {
                                int zorder=order-yorder-xorder;
                                tmpbiasfield-=C[ind]*currxpower[xorder]*currypower[yorder]*currzpower[zorder];
                                ind++;
                            }
                        }
                    }
                    BiasField[linearindexes+multispec*CurrSizes->numel]=tmpbiasfield;
                }
            }
        }

        for( int i=0; i<UsedBasisFunctions; i++)
        {
            BiasFieldCoefs[i+multispec*UsedBasisFunctions]=C[i];
        }

        delete [] Basis;
        delete [] Tempvar;
    }


}


void BiasCorrection_SPARCS(float * BiasField,
                           float * T1,
                           float * Expec,
                           float * Mask,
                           float * M,
                           float * V,
                           int biasOrder,
                           int nrOfClasses,
                           int aceletation_factor,
                           int xyzsize[3])
{


    if(aceletation_factor<1)
    {
        cout<<"ERROR: The acceleration factor has to be above or equal to 1"<<endl;
        flush(cout);
        return;
    }
    //int aceletation_factor=redux_factor_for_bias;
    int TotalLength = xyzsize[0]*xyzsize[1]*xyzsize[2];
    int UsedBasisFunctions=(int)((biasOrder+1) * (biasOrder+2)/2 *(biasOrder+3)/3);
    // Precompute Powers depending on the current BiasOrder
    int PowerOrder [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3))]= {0};
    int ind=0;
    for(int order=0; order<=biasOrder; order++)
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

    float invV[maxNumbClass];
    float currM[maxNumbClass];

    float * sampledData = &T1[0];
    // Precompute the M and V  inverses

    for(int i=0; i<nrOfClasses; i++)
    {
        invV[i]=1.0f/V[i];

        currM[i]=M[i];
    }

    float A [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)*((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)]= {0.0f};
    float B [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)]= {0.0f};
    float C [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)]= {0.0f};

    // Precompute sizes
    int maxix = (int)(xyzsize[0]);
    int maxiy = (int)(xyzsize[1]);
    int maxiz = (int)(xyzsize[2]);
    int Dims3d[3]= {0};
    Dims3d[0]=maxix;
    Dims3d[1]=maxiy;
    Dims3d[2]=maxiz;

    // Precompute number of samples as it was never computed
    int samplecount=0;
    int linearindexes=0;

    for (int iz=0; iz<maxiz; iz+=aceletation_factor)
    {
        for (int iy=0; iy<maxiy; iy+=aceletation_factor)
        {
            for (int ix=0; ix<maxix; ix+=aceletation_factor)
            {
                linearindexes=(iz)*(xyzsize[0])*(xyzsize[1])+(iy)*(xyzsize[0])+ix;
                if((Mask[linearindexes])>0)
                {
                    samplecount++;
                }
            }
        }
    }


    // CALC MATRIX A

    // Calc W (Van Leemput 1999 eq 7)

    float * Tempvar= new float [samplecount] ();
    float Tempvar_tmp=0;
    int tempvarindex=0;
    for (int iz=0; iz<maxiz; iz+=aceletation_factor)
    {
        for (int iy=0; iy<maxiy; iy+=aceletation_factor)
        {
            for (int ix=0; ix<maxix; ix+=aceletation_factor)
            {
                linearindexes=(iz)*(xyzsize[0])*(xyzsize[1])+(iy)*(xyzsize[0])+ix;
                if((Mask[linearindexes])>0)
                {
                    Tempvar_tmp=0;
                    for(int j=0; j<nrOfClasses; j++)
                    {
                        Tempvar_tmp+=Expec[linearindexes+TotalLength*j]*invV[j];
                    }
                    Tempvar[tempvarindex]=Tempvar_tmp;
                    tempvarindex++;
                }
            }
        }
    }

    // Precompute shifts
    float not_point_five_times_dims_x=(0.5f*(float)Dims3d[0]);
    float not_point_five_times_dims_y=(0.5f*(float)Dims3d[1]);
    float not_point_five_times_dims_z=(0.5f*(float)Dims3d[2]);
    float inv_not_point_five_times_dims_x=1.0f/(0.5f*(float)Dims3d[0]);
    float inv_not_point_five_times_dims_y=1.0f/(0.5f*(float)Dims3d[1]);
    float inv_not_point_five_times_dims_z=1.0f/(0.5f*(float)Dims3d[2]);
    float * Basis= new float[UsedBasisFunctions]();
    float xpos=0.0f;
    float ypos=0.0f;
    float zpos=0.0f;
    int x_bias_index_shift=0;
    int y_bias_index_shift=1;
    int z_bias_index_shift=2;
    float current_Tempvar=0.0f;

    // Calc A'WA (Van Leemput 1999 eq 7)
    tempvarindex=0;
    float * Basisptr1= (float *) Basis;
    float * Basisptr2= (float *) Basis;
    float * Aptr= (float *) A;
    for (int iz=0; iz<maxiz; iz+=aceletation_factor)
    {
        for (int iy=0; iy<maxiy; iy+=aceletation_factor)
        {
            for (int ix=0; ix<maxix; ix+=aceletation_factor)
            {
                linearindexes=(iz)*(xyzsize[0])*(xyzsize[1])+(iy)*(xyzsize[0])+ix;
                if((Mask[linearindexes])>0)
                {
                    current_Tempvar=Tempvar[tempvarindex];
                    xpos=(((float)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                    ypos=(((float)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                    zpos=(((float)iz-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z);
                    x_bias_index_shift=0;
                    y_bias_index_shift=1;
                    z_bias_index_shift=2;
                    Basisptr1= (float *) Basis;
                    for(int j2=0; j2<UsedBasisFunctions; j2++,x_bias_index_shift+=3,y_bias_index_shift+=3,z_bias_index_shift+=3, Basisptr1++)
                    {
                        *Basisptr1=(pow_int(xpos,PowerOrder[x_bias_index_shift])*pow_int(ypos,PowerOrder[y_bias_index_shift])*pow_int(zpos,PowerOrder[z_bias_index_shift]));
                    }
                    Basisptr1= (float *) Basis;
                    Aptr= (float *) A;
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


    // CALC MATRIX B
    //Precompute WR (Van Leemput 1999 eq 7)
    float Wi;
    float Wij;
    float Yest;
    float Ysum;
    tempvarindex=0;

    for (int iz=0; iz<maxiz; iz+=aceletation_factor)
    {
        for (int iy=0; iy<maxiy; iy+=aceletation_factor)
        {
            for (int ix=0; ix<maxix; ix+=aceletation_factor)
            {
                linearindexes=(iz)*(xyzsize[0])*(xyzsize[1])+(iy)*(xyzsize[0])+ix;
                if((Mask[linearindexes])>0)
                {
                    Wi=0;
                    Wij=0;
                    Yest=0;
                    Ysum=0;
                    for(int j=0; j<nrOfClasses; j++)
                    {
                        float tmpexpec = (float)Expec[linearindexes+TotalLength*j];
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


    for(int i2=0; i2<UsedBasisFunctions; i2++)
    {
        tempvarindex=0;
        B[i2]=0;
        for (int iz=0; iz<maxiz; iz+=aceletation_factor)
        {
            for (int iy=0; iy<maxiy; iy+=aceletation_factor)
            {
                for (int ix=0; ix<maxix; ix+=aceletation_factor)
                {
                    linearindexes=(iz)*(xyzsize[0])*(xyzsize[1])+(iy)*(xyzsize[0])+ix;
                    if((Mask[linearindexes])>0)
                    {
                        linearindexes=(iz)*(xyzsize[0])*(xyzsize[1])+(iy)*(xyzsize[0])+ix;
                        B[i2]+=pow_int((((float)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x),PowerOrder[0+i2*3])*
                               pow_int((((float)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y),PowerOrder[1+i2*3])*
                               pow_int((((float)iz-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z),PowerOrder[2+i2*3])*
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
    }

    seg_Matrix <double> RealB(UsedBasisFunctions,1);

    for(int i2=0; i2<UsedBasisFunctions; i2++)
    {
        RealB.setvalue(i2,0,(double)(B[i2]));
    }

    seg_Matrix <double> RealC(UsedBasisFunctions,1);

    RealC.settoproduct(RealA_inv,RealB);

    double cvalue=0.0f;
    bool success;
    for(int i2=0; i2<UsedBasisFunctions; i2++)
    {
        RealC.getvalue(i2,0,cvalue,success);
        C[i2]=(float)(cvalue);
    }

    for(int j2=0; j2<UsedBasisFunctions; j2++)
    {
        for(int i2=0; i2<UsedBasisFunctions; i2++)
        {
            double cvalue=0.0f;
            bool success;
            RealB.getvalue(i2,j2,cvalue,success);

        }
    }

    for (int iz=0; iz<maxiz; iz++)
    {
        for (int iy=0; iy<maxiy; iy++)
        {
            for (int ix=0; ix<maxix; ix++)
            {
                linearindexes=(iz)*(xyzsize[0])*(xyzsize[1])+(iy)*(xyzsize[0])+ix;
                if((Mask[linearindexes])>0)
                {
                    float tmpbiasfield=0.0f;
                    float currxpower[maxAllowedBCPowerOrder]={0};
                    float currypower[maxAllowedBCPowerOrder]={0};
                    float currzpower[maxAllowedBCPowerOrder]={0};
                    tmpbiasfield=0.0f;
                    xpos=(((float)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                    ypos=(((float)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                    zpos=(((float)iz-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z);

                    // Get the polynomial basis order
                    int order=1;
                    currxpower[0]=1;
                    currypower[0]=1;
                    currzpower[0]=1;
                    int orderminusone=0;
                    int maxorderplusone=biasOrder+1;
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
                    for(int order=0; order<=biasOrder; order++)
                    {
                        for(int xorder=0; xorder<=order; xorder++)
                        {
                            for(int yorder=0; yorder<=(order-xorder); yorder++)
                            {
                                int zorder=order-yorder-xorder;
                                tmpbiasfield-=C[ind]*currxpower[xorder]*currypower[yorder]*currzpower[zorder];
                                ind++;
                            }
                        }
                    }
                    BiasField[linearindexes]=tmpbiasfield;

                }
                else
                {
                    BiasField[linearindexes]=1;
                }
            }
        }
    }

    delete [] Basis;
    delete [] Tempvar;
}


void BiasCorrection_mask(segPrecisionTYPE * BiasField,
                         segPrecisionTYPE * BiasFieldCoefs,
                         nifti_image * T1,
                         int * Long_2_Short_Indices,
                         segPrecisionTYPE * Expec,
                         segPrecisionTYPE * Outlierness,
                         segPrecisionTYPE * M,
                         segPrecisionTYPE * V,
                         int biasOrder,
                         ImageSize * CurrSizes,
                         bool flag_Bias,
                         int verbose_level)
{

    if(verbose_level>0)
    {
        cout << "Optimising the Bias Field with order " << biasOrder<< endl;
        if(CurrSizes->usize>1)
        {
            cout<< "Assuming fully decoupled bias-fields" << endl;
        }
        flush(cout);
    }
    int reduxfactor=reduxFactorForBias;
    int nrOfClasses = CurrSizes->numclass;
    //nrOfClasses = 1;
    //int nrOfClasses = non_PV_numclass;
    int TotalLength = CurrSizes->numelmasked;
    int UsedBasisFunctions=(int)((biasOrder+1) * (biasOrder+2)/2 *(biasOrder+3)/3);
    segPrecisionTYPE * sampledData = static_cast<segPrecisionTYPE *>(T1->data);


    // Precompute Powers depending on the current BiasOrder
    int PowerOrder [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3))]= {0};
    int ind=0;
    for(int order=0; order<=biasOrder; order++)
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

    for(long multispec=0; multispec<CurrSizes->usize; multispec++)
    {
        sampledData = static_cast<segPrecisionTYPE *>(T1->data);
        sampledData = &sampledData[multispec*CurrSizes->numel];
        // Precompute the M and V  inverses
        for(int i=0; i<nrOfClasses; i++)
        {
            invV[i]=1.0f/(V[i*CurrSizes->usize*CurrSizes->usize+multispec+multispec*CurrSizes->usize]);
            currM[i]=M[i*CurrSizes->usize+multispec];
        }
        segPrecisionTYPE A [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)*((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)]= {0.0f};
        segPrecisionTYPE B [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)]= {0.0f};
        segPrecisionTYPE C [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)]= {0.0f};


// Precompute sizes
        int col_size = (int)(CurrSizes->xsize);
        int plane_size = (int)(CurrSizes->xsize)*(CurrSizes->ysize);
        int maxix = (int)(CurrSizes->xsize);
        int maxiy = (int)(CurrSizes->ysize);
        int maxiz = (int)(CurrSizes->zsize);
        int Dims3d[3]= {0};
        Dims3d[0]=maxix;
        Dims3d[1]=maxiy;
        Dims3d[2]=maxiz;

// Precompute number of samples as it was never computed
        int samplecount=0;
        int linearindexes=0;
        int currshortindex=0;
        if(CurrSizes->numelbias==0)
        {
            for (int iz=0; iz<maxiz; iz+=reduxfactor)
            {
                for (int iy=0; iy<maxiy; iy+=reduxfactor)
                {
                    for (int ix=0; ix<maxix; ix+=reduxfactor)
                    {
                        linearindexes=iz*plane_size+iy*col_size+ix;
                        if(Long_2_Short_Indices[linearindexes]>=0)
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
            CurrSizes->numelbias=samplecount;
        }
        else
        {
            samplecount=CurrSizes->numelbias;
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
                    currshortindex=Long_2_Short_Indices[iz*plane_size+iy*col_size+ix];
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
                    linearindexes=(iz)*(CurrSizes->xsize)*(CurrSizes->ysize)+(iy)*(CurrSizes->xsize)+ix;
                    currshortindex=Long_2_Short_Indices[linearindexes];
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
                    linearindexes=(iz)*(CurrSizes->xsize)*(CurrSizes->ysize)+(iy)*(CurrSizes->xsize)+ix;
                    currshortindex=Long_2_Short_Indices[linearindexes];
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
                        linearindexes2=(iz2)*(CurrSizes->xsize)*(CurrSizes->ysize)+(iy2)*(CurrSizes->xsize)+ix2;
                        if(Long_2_Short_Indices[linearindexes2]>=0)
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
            segPrecisionTYPE currxpower[maxAllowedBCPowerOrder]={0};
            segPrecisionTYPE currypower[maxAllowedBCPowerOrder]={0};
            segPrecisionTYPE currzpower[maxAllowedBCPowerOrder]={0};
            segPrecisionTYPE tmpbiasfield=0.0f;
            for (int iy=0; iy<maxiy; iy++)
            {
                for (int ix=0; ix<maxix; ix++)
                {
                    linearindexes=(iz)*(CurrSizes->xsize)*(CurrSizes->ysize)+(iy)*(CurrSizes->xsize)+ix;
                    currshortindex=Long_2_Short_Indices[linearindexes];
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
                        int maxorderplusone=biasOrder+1;
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
                        for(int order=0; order<=biasOrder; order++)
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
                        BiasField[currshortindex+multispec*CurrSizes->numelmasked]=tmpbiasfield;
                    }
                }
            }
        }

        for( int i=0; i<UsedBasisFunctions; i++)
        {
            BiasFieldCoefs[i+multispec*UsedBasisFunctions]=C[i];
        }
        delete [] Basis;
        delete [] Tempvar;
    }
}


void BiasCorrection2D(segPrecisionTYPE * BiasField,
                      segPrecisionTYPE * BiasFieldCoefs,
                      nifti_image * T1,
                      segPrecisionTYPE * Expec,
                      segPrecisionTYPE * Outlierness,
                      segPrecisionTYPE * M,
                      segPrecisionTYPE * V,
                      int biasOrder,
                      ImageSize * CurrSizes,
                      bool flag_Bias,
                      int verbose_level)
{

    if(verbose_level>0)
    {
        cout << "Optimising the Bias Field with order " << biasOrder<< endl;
        if(CurrSizes->usize>1)
        {
            cout<< "Assuming fully decoupled bias-fields" << endl;
        }
        flush(cout);
    }
    int reduxfactor=reduxFactorForBias;
    int nrOfClasses = CurrSizes->numclass;

    int TotalLength = CurrSizes->numel;
    int UsedBasisFunctions=(int)((biasOrder+1) * (biasOrder+2)/2);
    segPrecisionTYPE * sampledData = static_cast<segPrecisionTYPE *>(T1->data);


    // Precompute Powers depending on the current BiasOrder
    int PowerOrder [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2)]= {0};
    int ind=0;
    for(int order=0; order<=min(biasOrder,maxAllowedBCPowerOrder); order++)
    {
        for(int xorder=0; xorder<=order; xorder++)
        {
            int yorder=order-xorder;
            PowerOrder[ind] =xorder;
            PowerOrder[ind+1] =yorder;
            ind += 3;
        }
    }
    segPrecisionTYPE invV[maxNumbClass];
    segPrecisionTYPE currM[maxNumbClass];
    for(long i=0; i<maxNumbClass; i++)
    {
        invV[i]=0;
        currM[i]=0;
    }

    for(long multispec=0; multispec<CurrSizes->usize; multispec++)
    {
        sampledData = &sampledData[multispec*CurrSizes->numel];
        // Precompute the M and V  inverses

        for(int i=0; i<nrOfClasses; i++)
        {
            invV[i]=1.0f/(V[i*CurrSizes->usize*CurrSizes->usize+multispec+multispec*CurrSizes->usize]);
            currM[i]=M[i*CurrSizes->usize+multispec];
        }



        segPrecisionTYPE A [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2)*((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2)]= {0.0f};
        segPrecisionTYPE B [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2)]= {0.0f};
        segPrecisionTYPE C [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2)]= {0.0f};


// Precompute sizes
        int col_size = (int)(CurrSizes->xsize);
        int maxix = (int)(CurrSizes->xsize);
        int maxiy = (int)(CurrSizes->ysize);
        int Dims3d[2]= {0};
        Dims3d[0]=maxix;
        Dims3d[1]=maxiy;

// Precompute number of samples as it was never computed
        int samplecount=0;
        int linearindexes=0;
        int currindex=0;
        if(CurrSizes->numelbias==0)
        {
            for (int iy=0; iy<maxiy; iy+=reduxfactor)
            {
                for (int ix=0; ix<maxix; ix+=reduxfactor)
                {
                    samplecount++;
                }
            }
            if(verbose_level>0)
            {
                cout << "Samplecount = " << samplecount<<"\n";
                flush(cout);
            }
            CurrSizes->numelbias=samplecount;
        }
        else
        {
            samplecount=CurrSizes->numelbias;
        }



// CALC MATRIX A

// Calc W (Van Leemput 1999 eq 7)

        segPrecisionTYPE * Tempvar= new segPrecisionTYPE [samplecount] ();
        segPrecisionTYPE Tempvar_tmp=0;
        currindex=0;
        int tempvarindex=0;
        for (int iy=0; iy<maxiy; iy+=reduxfactor)
        {
            for (int ix=0; ix<maxix; ix+=reduxfactor)
            {
                currindex=iy*col_size+ix;
                Tempvar_tmp=0;
                for(int j=0; j<nrOfClasses; j++)
                {
                    Tempvar_tmp+=Expec[currindex+TotalLength*j]*invV[j];
                }
                Tempvar[tempvarindex]=Tempvar_tmp;
                tempvarindex++;
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
                linearindexes=(iy)*(CurrSizes->xsize)+ix;
                Basisptr1= (segPrecisionTYPE *) Basis;
                current_Tempvar=Tempvar[tempvarindex];
                xpos=(((segPrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                ypos=(((segPrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                x_bias_index_shift=0;
                y_bias_index_shift=1;
                for(int j2=0; j2<UsedBasisFunctions; j2++,x_bias_index_shift+=2,y_bias_index_shift+=2, Basisptr1++)
                {
                    // Because Powerorder is alwais int, use a special power function (faster)
                    *Basisptr1=(pow_int(xpos,PowerOrder[x_bias_index_shift])*pow_int(ypos,PowerOrder[y_bias_index_shift]));
                    // Instead, although slower, one can use
                    // TmpA=(pow(xpos,PowerOrder[0+j2*3])*pow(ypos,PowerOrder[1+j2*3])*pow(zpos,PowerOrder[2+j2*3]));
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
                linearindexes=(iy)*(CurrSizes->xsize)+ix;

                Wi=0;
                Wij=0;
                Yest=0;
                Ysum=0;
                for(int j=0; j<nrOfClasses; j++)
                {
                    segPrecisionTYPE tmpexpec = (segPrecisionTYPE)Expec[linearindexes+TotalLength*j];
                    Wij=tmpexpec*(invV[j]);
                    Wi+=Wij;
                    Yest+=Wij*(currM[j]);
                    Ysum+=Wij;
                }
                Tempvar[tempvarindex]=Wi*(sampledData[linearindexes]-(Yest/Ysum));
                tempvarindex++;

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
                    linearindexes=(iy)*(CurrSizes->xsize)+ix;
                    B[i2]+=pow_int((((segPrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x),PowerOrder[0+i2*3])*
                           pow_int((((segPrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y),PowerOrder[1+i2*3])*
                           Tempvar[tempvarindex];
                    if(B[i2]!=B[i2])
                    {
                        B[i2]=1;
                    }
                    tempvarindex++;
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

        segPrecisionTYPE currxpower[maxAllowedBCPowerOrder]={0};
        segPrecisionTYPE currypower[maxAllowedBCPowerOrder]={0};
        segPrecisionTYPE tmpbiasfield=0.0f;
        for (int iy=0; iy<maxiy; iy++)
        {
            for (int ix=0; ix<maxix; ix++)
            {
                linearindexes=(iy)*(CurrSizes->xsize)+ix;
                tmpbiasfield=0.0f;
                xpos=(((segPrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                ypos=(((segPrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);

                // Get the polynomial order power
                int order=1;
                currxpower[0]=1;
                currypower[0]=1;
                int orderminusone=0;
                int maxorderplusone=biasOrder+1;
                while (order<maxorderplusone)
                {
                    currxpower[order]=currxpower[orderminusone]*xpos;
                    currypower[order]=currypower[orderminusone]*ypos;
                    order++;
                    orderminusone++;
                }


                ind=0;
                for(int order=0; order<=biasOrder; order++)
                {
                    for(int xorder=0; xorder<=order; xorder++)
                    {
                        int yorder=order-xorder;
                        tmpbiasfield-=C[ind]*currxpower[xorder]*currypower[yorder];
                        ind++;
                    }
                }
                BiasField[linearindexes+multispec*CurrSizes->numel]=tmpbiasfield;
            }
        }

        for( int i=0; i<UsedBasisFunctions; i++)
        {
            BiasFieldCoefs[i+multispec*UsedBasisFunctions]=C[i];
        }

        delete [] Basis;
        delete [] Tempvar;
    }


}





void BiasCorrection_mask2D(segPrecisionTYPE * BiasField,
                           segPrecisionTYPE * BiasFieldCoefs,
                           nifti_image * T1,
                           int * Long_2_Short_Indices,
                           segPrecisionTYPE * Expec,
                           segPrecisionTYPE * Outlierness,
                           segPrecisionTYPE * M,
                           segPrecisionTYPE * V,
                           int biasOrder,
                           ImageSize * CurrSizes,
                           bool flag_Bias,
                           int verbose_level)
{
    if(verbose_level>0)
    {
        cout << "Optimising the Bias Field with order " << biasOrder<< endl;
        if(CurrSizes->usize>1)
        {
            cout<< "Assuming fully decoupled bias-fields" << endl;
        }
        flush(cout);
    }
    int reduxfactor=reduxFactorForBias;
    long nrOfClasses = CurrSizes->numclass;
    //nrOfClasses = 1;
    //int nrOfClasses = non_PV_numclass;
    long TotalLength = CurrSizes->numelmasked;
    int UsedBasisFunctions=(int)((biasOrder+1) * (biasOrder+2)/2);
    segPrecisionTYPE * sampledData = static_cast<segPrecisionTYPE *>(T1->data);


    // Precompute Powers depending on the current BiasOrder
    int PowerOrder [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2)]= {0};
    int ind=0;
    for(int order=0; order<=biasOrder; order++)
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

    for(long multispec=0; multispec<CurrSizes->usize; multispec++)
    {
        sampledData = static_cast<segPrecisionTYPE *>(T1->data);
        sampledData = &sampledData[multispec*CurrSizes->numel];
        // Precompute the M and V  inverses
        for(long i=0; i<nrOfClasses; i++)
        {
            invV[i]=1.0f/(V[i*CurrSizes->usize*CurrSizes->usize+multispec+multispec*CurrSizes->usize]);
            currM[i]=M[i*CurrSizes->usize+multispec];
        }
        segPrecisionTYPE A [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2))/2*((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2)]= {0.0f};
        segPrecisionTYPE B [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2))/2]= {0.0f};
        segPrecisionTYPE C [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2))/2]= {0.0f};


// Precompute sizes
        int col_size = (int)(CurrSizes->xsize);
        int maxix = (int)(CurrSizes->xsize);
        int maxiy = (int)(CurrSizes->ysize);
        int Dims3d[2]= {0};
        Dims3d[0]=maxix;
        Dims3d[1]=maxiy;

// Precompute number of samples as it was never computed
        int samplecount=0;
        int linearindexes=0;
        int currshortindex=0;
        if(CurrSizes->numelbias==0)
        {
            for (int iy=0; iy<maxiy; iy+=reduxfactor)
            {
                for (int ix=0; ix<maxix; ix+=reduxfactor)
                {
                    linearindexes=iy*col_size+ix;
                    if(Long_2_Short_Indices[linearindexes]>=0)
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
            CurrSizes->numelbias=samplecount;
        }
        else
        {
            samplecount=CurrSizes->numelbias;
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
                currshortindex=Long_2_Short_Indices[iy*col_size+ix];
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
                linearindexes=(iy)*(CurrSizes->xsize)+ix;
                currshortindex=Long_2_Short_Indices[linearindexes];
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
                linearindexes=(iy)*(CurrSizes->xsize)+ix;
                currshortindex=Long_2_Short_Indices[linearindexes];
                if(currshortindex>=0)
                {
                    Wi=0;
                    Wij=0;
                    Yest=0;
                    Ysum=0;
                    for(long j=0; j<nrOfClasses; j++)
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

        for(int i2=0; i2<UsedBasisFunctions; i2++)
        {
            tempvarindex=0;
            B[i2]=0;
            for (int iy=0; iy<maxiy; iy+=reduxfactor)
            {

                for (int ix=0; ix<maxix; ix+=reduxfactor)
                {
                    linearindexes=(iy)*(CurrSizes->xsize)+ix;
                    currshortindex=Long_2_Short_Indices[linearindexes];
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

        segPrecisionTYPE currxpower[maxAllowedBCPowerOrder]={0};
        segPrecisionTYPE currypower[maxAllowedBCPowerOrder]={0};
        segPrecisionTYPE tmpbiasfield=0.0f;
        for (int iy=0; iy<maxiy; iy++)
        {
            for (int ix=0; ix<maxix; ix++)
            {
                linearindexes=(iy)*(CurrSizes->xsize)+ix;
                currshortindex=Long_2_Short_Indices[linearindexes];
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
                    int maxorderplusone=biasOrder+1;
                    while (order<maxorderplusone)
                    {
                        currxpower[order]=currxpower[orderminusone]*xpos;
                        currypower[order]=currypower[orderminusone]*ypos;
                        order++;
                        orderminusone++;
                    }

                    ind=0;
                    for(int order=0; order<=biasOrder; order++)
                    {
                        for(int xorder=0; xorder<=order; xorder++)
                        {
                            int yorder=order-xorder;
                            tmpbiasfield-=C[ind]*currxpower[xorder]*currypower[yorder];
                            ind++;
                        }
                    }
                    BiasField[currshortindex+multispec*CurrSizes->numelmasked]=tmpbiasfield;
                }
            }

        }

        for( int i=0; i<UsedBasisFunctions; i++)
        {
            BiasFieldCoefs[i+multispec*UsedBasisFunctions]=C[i];
        }
        delete [] Basis;
        delete [] Tempvar;
    }
}

