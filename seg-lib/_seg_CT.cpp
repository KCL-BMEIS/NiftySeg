#include "_seg_CT.h"

nifti_image * CT_Lung(nifti_image * T1, nifti_image * Mask,char* filename_out, SEG_PARAM * segment_param){


    PrecisionTYPE * LungCT_PTR = static_cast<PrecisionTYPE *>(T1->data);
    int dimensions[3];
    dimensions[0]=T1->nx;
    dimensions[1]=T1->ny;
    dimensions[2]=T1->nz;


    /* Creating Mask
    float thresh=800.0f;
    bool * Mask = new bool [dimensions[0]*dimensions[1]*dimensions[2]]();
    for(int i=0; i<T1->nvox; i++)Mask[i]=(PrecisionTYPE)(LungCT_PTR[i]<thresh & LungCT_PTR[i]>100.0f);
    cout << "Filling in the Thresholded mask" << endl;
    /*Dillate(Mask,1,dimensions);
    for(int i=0; i<T1->nvox; i++)Mask[i]=!Mask[i];
    Dillate(Mask,2,dimensions);
    for(int i=0; i<T1->nvox; i++)Mask[i]=!Mask[i];
    Dillate(Mask,1,dimensions);
    for(int i=0; i<T1->nvox; i++)Mask[i]=(int)(Mask[i]*(LungCT_PTR[i]<700));
    Dillate(Mask,1,dimensions);

    int * Buffer_Int = new int [dimensions[0]*dimensions[1]*dimensions[2]]();
    for(int i=0; i<T1->nvox; i++)Buffer_Int[i]=(int)(Mask[i]);

    int * Buffer_Int2 = new int [dimensions[0]*dimensions[1]*dimensions[2]]();
    cout << "Finding the connected components" << endl;
    ConnectComp(Buffer_Int,Buffer_Int2,dimensions);
    for(int i=0; i<T1->nvox; i++)Mask[i]=(bool)(Buffer_Int2[i]==1);
    delete [] Buffer_Int;
    delete [] Buffer_Int2;
    cout << "Re-filling in the Thresholded mask" << endl;

    Dillate(Mask,1,dimensions);
    for(int i=0; i<T1->nvox; i++)Mask[i]=(bool)(Mask[i] && LungCT_PTR[i]<thresh & LungCT_PTR[i]>100.0f);
    Dillate(Mask,1,dimensions);
    for(int i=0; i<T1->nvox; i++)Mask[i]=!Mask[i];
    Dillate(Mask,2,dimensions);
    for(int i=0; i<T1->nvox; i++)Mask[i]=!Mask[i];
*/

    // Start Segmentation
    ImageSize * CurrSizes = new ImageSize [1]();
    CurrSizes->numel=(int)(T1->nx*T1->ny*T1->nz);
    CurrSizes->xsize=T1->nx;
    CurrSizes->ysize=T1->ny;
    CurrSizes->zsize=T1->nz;
    CurrSizes->numclass=2;
    CurrSizes->numelmasked=0;
    CurrSizes->numelbias=0;

    if(Mask->datatype!=DT_BINARY){seg_convert2binary(Mask,0.0f);}
    int * S2L = Create_Short_2_Long_Matrix_from_NII(Mask,&(CurrSizes->numelmasked));
    int * L2S = Create_Long_2_Short_Matrix_from_NII(Mask);
    //int * S2L = Create_Short_2_Long_Matrix_from_Carray(Mask,&(CurrSizes->numelmasked),CurrSizes->numel);
    //int * L2S = Create_Long_2_Short_Matrix_from_Carray(Mask,CurrSizes->numel);

    PrecisionTYPE * Expec =  new PrecisionTYPE [CurrSizes->numclass*CurrSizes->numelmasked]();
    PrecisionTYPE * MRF = new PrecisionTYPE[CurrSizes->numclass*CurrSizes->numelmasked]();
    PrecisionTYPE * ShortPrior = new PrecisionTYPE[CurrSizes->numclass*CurrSizes->numelmasked]();
    for(int i=0; i<(CurrSizes->numelmasked*(CurrSizes->numclass)); i++){
        MRF[i]=1.0f;
        ShortPrior[i]=1.0f;
    }
    PrecisionTYPE * BiasField = new PrecisionTYPE[CurrSizes->numelmasked]();
    PrecisionTYPE * BiasFieldCoefs = new PrecisionTYPE[((maxallowedpowerorder+1)*(maxallowedpowerorder+2)/2*(maxallowedpowerorder+3)/3)]();
    for(int i=0; i<CurrSizes->numelmasked; i++){BiasField[i]=0.0f;}
    PrecisionTYPE * M = new PrecisionTYPE[CurrSizes->numclass]();
    PrecisionTYPE * V = new PrecisionTYPE[CurrSizes->numclass]();
    M[0]=300.0f;
    M[1]=1000.0f;
    V[0]=200.0f;
    V[1]=50.0f;
    Normalize_T1_and_MV(T1,Mask,M,V,CurrSizes);

    PrecisionTYPE * G = new PrecisionTYPE [segment_param->numb_classes*segment_param->numb_classes]();
    PrecisionTYPE * H = new PrecisionTYPE [segment_param->numb_classes*segment_param->numb_classes]();
    Create_diagonal_GH_Nclass(G,H,1,segment_param);

    float loglik=2;
    float oldloglik=1;
    // Iterative Components - EM, MRF, Bias Correction
    //Expectation and LogLik
    bool out=true;
    int iter=1;
    while (out) {
        calcE_class_mask(T1,MRF,Expec,&loglik,BiasField,S2L,M,V,CurrSizes);
        //Maximization
        calcM_class_mask(T1,Expec,BiasField,S2L,M,V,CurrSizes,segment_param->verbose_level);

        MRFregularization_mask(Expec,G,H,NULL,MRF,ShortPrior,L2S,S2L,CurrSizes,segment_param->flag_MRF,segment_param->verbose_level);
        //Bias Correction
        BiasCorrection_mask(BiasField,BiasFieldCoefs,T1,L2S,Expec,M,V,4,CurrSizes,segment_param->flag_Bias,segment_param->verbose_level);
        // Print LogLik
        if(segment_param->verbose_level>0){
            printloglik(iter,loglik,oldloglik);
        }

        if(((loglik-oldloglik)/fabs(oldloglik))<(PrecisionTYPE)(0.005) && iter>3){
            out=false;
        }
        oldloglik=loglik;
        iter++;
    }


    // OUTPUT!!
    nifti_image * Result = Copy_ShortExpec_to_Result_Lungs(Expec,BiasField,S2L,T1,filename_out,CurrSizes);
    cout << "Finished" << endl;

    delete [] Expec;
    delete [] MRF;
    delete [] BiasField;
    delete [] S2L;
    delete [] L2S;
    return Result;
}






//

nifti_image * Copy_ShortExpec_to_Result_Lungs(PrecisionTYPE * Expec,
                                              PrecisionTYPE * BiasField,
                                              int * Short_2_Long_Indices,
                                              nifti_image * T1,
                                              char * filename,
                                              ImageSize * CurrSizes){

    nifti_image * Result = nifti_copy_nim_info(T1);
    Result->dim[0]=4;
    Result->dim[4]=CurrSizes->numclass+2;
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
    PrecisionTYPE * Resultdata_class = &Resultdata[class_nvox];
    for(int i=0; i<CurrSizes->numelmasked; i++,Short_2_Long_Indices_PRT++){
        //Resultdata[*Short_2_Long_Indices_PRT]=(exp(BiasField[i]+T1data[*Short_2_Long_Indices_PRT])-1)/0.6;
        to_resize=exp((BiasField[i]+T1data[*Short_2_Long_Indices_PRT])*0.693147181)-1;
        Resultdata[*Short_2_Long_Indices_PRT]=(to_resize*(CurrSizes->rescale_max-CurrSizes->rescale_min)+CurrSizes->rescale_min);
        Resultdata_class[*Short_2_Long_Indices_PRT]=(exp(BiasField[i]*0.693147181)-1)*(CurrSizes->rescale_max-CurrSizes->rescale_min)+CurrSizes->rescale_min;
    }

    for(int currclass=0; currclass<CurrSizes->numclass;currclass++){

        PrecisionTYPE * Resultdata_class = &Resultdata[(currclass+2)*class_nvox];
        PrecisionTYPE * Expec_PTR = &Expec[(currclass)*CurrSizes->numelmasked];
        Short_2_Long_Indices_PRT= (int *) Short_2_Long_Indices;

        for(int i=0; i<CurrSizes->numelmasked; i++,Short_2_Long_Indices_PRT++,Expec_PTR++){
            Resultdata_class[*Short_2_Long_Indices_PRT]=*Expec_PTR;
        }
    }
    return Result;
}
