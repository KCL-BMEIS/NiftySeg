#include "_seg_LoAd.h"

nifti_image * LoAd_Segment(nifti_image * T1, nifti_image * Mask, nifti_image * Priors, SEG_PARAM * segment_param){

    time_t start,end;
    time(&start);
    if((int)(segment_param->verbose_level)>(int)(0)){
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

    FLAGS * flags = new FLAGS [1]();
    flags->improv_phase=0;
    flags->out=1;
    flags->oldloglik=(PrecisionTYPE)1.0;
    flags->loglik=(PrecisionTYPE)2.0;
    flags->prior_relax=!((segment_param->relax_factor)>0);
    flags->pv_modeling=!(segment_param->flag_PV_model);;
    flags->sg_delineation=!(segment_param->flag_SG_deli);;

    // Creating G and H matrixes
    PrecisionTYPE be=segment_param->MRF_strength*1.0f;
    PrecisionTYPE ba=segment_param->MRF_strength*5.0f;
    PrecisionTYPE ratioGH=1.0f;

    // Alocate more than enough space for G and H even if the PV option is active
    PrecisionTYPE G[PV_numbclass*PV_numbclass]={0.0};
    PrecisionTYPE H[PV_numbclass*PV_numbclass]={0.0};
    Create_GH_5class(G,H,ba,be,ratioGH,segment_param);

    // Sanity check... Normalizing Priors and removing NaN
    if(segment_param->verbose_level>0){cout << "Normalizing Data" << endl;}
    if(Mask->datatype!=DT_BINARY){
        seg_convert2binary(Mask,0.0f);
    }
    Normalize_NaN_Priors_mask(Priors,Mask,segment_param->verbose_level);
    Normalize_Image_mask(T1,Mask,CurrSizes,segment_param->verbose_level);
    int * Short_2_Long_Indices = Create_Short_2_Long_Matrix_from_NII(Mask,&(CurrSizes->numelmasked));
    int * Long_2_Short_Indices = Create_Long_2_Short_Matrix_from_NII(Mask);

    // if PV modeling is on, prealocate extra memory for the 2 PV classes

    PrecisionTYPE * Expec = Create_cArray_from_Prior_mask(Mask,Priors,CurrSizes->numclass,segment_param->flag_PV_model);
    PrecisionTYPE * ShortPrior = Create_cArray_from_Prior_mask(Mask,Priors,CurrSizes->numclass,segment_param->flag_PV_model);
    PrecisionTYPE * MRF = new PrecisionTYPE[CurrSizes->numelmasked*(CurrSizes->numclass+(int)(segment_param->flag_PV_model*2))]();
    for(int i=0; i<(CurrSizes->numelmasked*(CurrSizes->numclass+(int)(segment_param->flag_PV_model*2))); i++){
        MRF[i]=(PrecisionTYPE)(1.0);
    }
    PrecisionTYPE * BiasField = new PrecisionTYPE[CurrSizes->numelmasked]();
    PrecisionTYPE * BiasFieldCoefs = new PrecisionTYPE[((maxallowedpowerorder+1)*(maxallowedpowerorder+2)/2*(maxallowedpowerorder+3)/3)]();
    for(int i=0; i<CurrSizes->numelmasked; i++){
        BiasField[i]=0;
    }
    PrecisionTYPE M [PV_numbclass]={0.0f};
    PrecisionTYPE V [PV_numbclass]={0.0f};
    PrecisionTYPE * MRF_Beta =NULL;


    // Initialize
    calcM_mask(T1,Expec,BiasField,Short_2_Long_Indices,M, V,NULL,NULL,CurrSizes,segment_param->verbose_level);
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
    while (flags->out && iter < segment_param->maxIteration) {
        iter++;
        iterchange++;
        if(segment_param->verbose_level>0){
            cout << "\n*******************************" << endl;
            cout << "Iteration " << iter << endl;
        }

        // Iterative Components - EM, MRF, Bias Correction
        //Expectation and LogLik
        calcE_mask(T1,MRF,Expec,&flags->loglik,BiasField,Short_2_Long_Indices,M,V,CurrSizes,segment_param->verbose_level);
        //Bias Correction
        if(BFupdate)BiasCorrection_mask(BiasField,BiasFieldCoefs,T1,Long_2_Short_Indices,Expec,M,V,CurrentBF,CurrSizes,segment_param->flag_Bias,segment_param->verbose_level);
        //MRF
        if(MRFupdate)MRFregularization_mask(Expec,G,H,MRF_Beta,MRF,ShortPrior,Long_2_Short_Indices,Short_2_Long_Indices,CurrSizes,segment_param->flag_MRF,segment_param->verbose_level);
        // Print LogLik depending on the verbose level
        if(segment_param->verbose_level>0)printloglik(iter,flags->loglik,flags->oldloglik);
        //Maximization
        calcM_mask_LoAd(T1,Expec,BiasField,Short_2_Long_Indices,M,V,CurrSizes,segment_param->verbose_level,flags->pv_modeling);

        // Preform Segmentation Refinement Steps or Exit
        if((((flags->loglik-flags->oldloglik)/fabs(flags->oldloglik))<(PrecisionTYPE)(0.005) && iterchange>3) || iter>segment_param->maxIteration){
            switch(flags->improv_phase){
            case 0: // Activate BC
                if(BFupdate==0 || CurrentBF<segment_param->bias_order){
                    flags->loglik=1;
                    CurrentBF++;
                    BFupdate=true;
                    if(CurrentBF>=segment_param->bias_order)flags->improv_phase++;
                    break;
                }
                else{flags->improv_phase++;}
             case 1: // Activate MRF
                 if(MRFupdate==0){
                     flags->improv_phase++;
                     MRFupdate=true;
                     flags->loglik=1;
                     iterchange=0;
                     break;
                 }
                 else{flags->improv_phase++;}
             case 2: // Run Prior Relaxation
                 if((segment_param->relax_factor)>0 && flags->prior_relax==false){
                     flags->improv_phase++;
                     flags->prior_relax=true;
                     Relax_Priors(ShortPrior,Expec,MRF,Short_2_Long_Indices,Long_2_Short_Indices,segment_param->relax_factor,G,ba,be,CurrSizes,segment_param);
                     flags->loglik=1;
                     iterchange=0;
                     break;
                 }
                 else{flags->improv_phase++;}
             case 3: // Run Implicit PV modeling
                 if(flags->pv_modeling==false){
                     flags->improv_phase++;
                     flags->pv_modeling=true;
                     Convert_to_PV(T1,BiasField,ShortPrior,Expec,MRF,M,V,Short_2_Long_Indices,Long_2_Short_Indices,CurrSizes,segment_param);
                     Create_GH_7class(G,H,ba,be,ratioGH,segment_param);
                     CurrSizes->numclass=7;
                     //segment_param->flag_Bias=false;
                     flags->loglik=1;
                     iterchange=0;
                     break;}
                 else{flags->improv_phase++;}
            case 4: // Run Sulci and Gyri deliniation
                if(flags->sg_delineation==false){
                    if(MRF_Beta==NULL){
                        MRF_Beta=new PrecisionTYPE [CurrSizes->numelmasked]();
                    }

                    if((flags->loglik-last_SGcorrect_loglik)/fabs(last_SGcorrect_loglik)<(PrecisionTYPE)(0.05) ){
                        flags->sg_delineation=true;
                        flags->out=false;
                    }
                    else{
                        Sulci_and_gyri_correction(MRF_Beta,ShortPrior,Expec,MRF,Short_2_Long_Indices,Long_2_Short_Indices,CurrSizes);
                        last_SGcorrect_loglik=flags->loglik;
                        flags->loglik=1;
                        iterchange=0;
                    }

                    break;}
                else{flags->improv_phase++;}
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
    if(true){
        if(flags->pv_modeling){
        Result = Copy_ShortExpec_to_Result(T1,Expec,BiasField,BiasFieldCoefs,Short_2_Long_Indices,Priors,segment_param,M,CurrSizes);
    }
        else{
            CurrSizes->numclass=5;
            Result = Copy_Expec_to_Result_mask(Expec,Short_2_Long_Indices,T1,segment_param->filename_out,CurrSizes);
        }
    }
    else{
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

    if(segment_param->verbose_level>0){
        cout << "Finished in "<<difftime(end,start)<<"sec"<< endl;
    }
    return Result;
}

