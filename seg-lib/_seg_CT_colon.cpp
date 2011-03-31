#include "_seg_CT_colon.h"

nifti_image * CT_colon(nifti_image * CT,char* filename_out, SEG_PARAM * segment_param){





    PrecisionTYPE * CT_PTR = static_cast<PrecisionTYPE *>(CT->data);


    PrecisionTYPE * M = new PrecisionTYPE[segment_param->numb_classes]();
    PrecisionTYPE * V = new PrecisionTYPE[segment_param->numb_classes]();
    M[0]=50.0f;
    M[1]=850.0f;
    M[2]=1350.0f;
    V[0]=100.0f;
    V[1]=300.0f;
    V[2]=35.0f;



    int dimensions[3];
    dimensions[0]=CT->nx;
    dimensions[1]=CT->ny;
    dimensions[2]=CT->nz;
    cout << dimensions[0] << ' '<<dimensions[1]  << ' '<<dimensions[2] << endl;

    if(CT->datatype!=NIFTI_TYPE_FLOAT32){
        seg_changeDatatype<float>(CT);
    }

    time_t start,end;
    time(&start);

    bool * Mask_carr = new bool [dimensions[0]*dimensions[1]*dimensions[2]]();
    for(int i=0; i<CT->nvox; i++)Mask_carr[i]=(bool)(CT_PTR[i]<(M[0]+400));
    cout << "Filling in the ROI";
    flush(cout);


    for(int i=0; i<CT->nvox; i++){
        Mask_carr[i]=!Mask_carr[i];
    }
    Dillate(Mask_carr,3,dimensions);
    for(int i=0; i<CT->nvox; i++){
        Mask_carr[i]=!Mask_carr[i];
    }
    Dillate(Mask_carr,5,dimensions);
    for(int i=0; i<CT->nvox; i++){
        Mask_carr[i]=!Mask_carr[i];
    }
    Dillate(Mask_carr,2,dimensions);
    for(int i=0; i<CT->nvox; i++){
        Mask_carr[i]=!Mask_carr[i];
    }
    cout << " - Done" << endl;





    int * Buffer_Int = new int [dimensions[0]*dimensions[1]*dimensions[2]]();
    for(int i=0; i<(dimensions[0]*dimensions[1]*dimensions[2]); i++)Buffer_Int[i]=(int)(Mask_carr[i]);

    int * Buffer_Int2 = new int [dimensions[0]*dimensions[1]*dimensions[2]]();
    cout << "Finding the connected components";
    flush(cout);
    ConnectComp(Buffer_Int,Buffer_Int2,dimensions);
    for(int i=0; i<(dimensions[0]*dimensions[1]*dimensions[2]); i++)Mask_carr[i]=(bool)(Buffer_Int2[i]==1);
    cout << " - Done" << endl;
    delete [] Buffer_Int;
    delete [] Buffer_Int2;


    cout << "Expanding the ROI";
    flush(cout);
    Dillate(Mask_carr,25,dimensions);
    cout << " - Done" << endl;

    flush(cout);

    nifti_image * Mask = Copy_single_image_to_Result(Mask_carr,CT,filename_out);



    int mycounter=0;
    for(int i=0; i<CT->nvox; i++){
        if(Mask_carr[i]>0){
            mycounter++;
        }
    }

    int * tmpvec = new int [mycounter]();
    int newindex=0;
    for(int i=0; i<CT->nvox; i++){
        if(Mask_carr[i]>0){
            tmpvec[newindex]=(int)(CT_PTR[i]);
            newindex++;
        }
    }
    cout << "Sorting - " << mycounter << " long array" << endl;
    flush(cout);
    quickSort(tmpvec,mycounter);

    cout<< "99.5 percentile - "<< tmpvec[(int)(floorf(mycounter*0.995f))] <<endl;
    if(tmpvec[(int)(floorf(mycounter*0.995f))]>M[2]){
        M[2]=tmpvec[(int)(floorf(mycounter*0.995f))]-50;
    }
    cout << "Start Segmentation";
    flush(cout);
    delete [] tmpvec;
    delete [] Mask_carr;

    ImageSize * CurrSizes = new ImageSize [1]();
    CurrSizes->numel=(int)(CT->nx*CT->ny*CT->nz);
    CurrSizes->xsize=CT->nx;
    CurrSizes->ysize=CT->ny;
    CurrSizes->zsize=CT->nz;

    CurrSizes->numclass=segment_param->numb_classes;
    CurrSizes->numelmasked=0;
    CurrSizes->numelbias=0;

    if(Mask->datatype!=DT_BINARY){seg_convert2binary(Mask,0.0f);}
    int * S2L = Create_Short_2_Long_Matrix_from_NII(Mask,&(CurrSizes->numelmasked));
    int * L2S = Create_Long_2_Short_Matrix_from_NII(Mask);

    PrecisionTYPE * Expec =  new PrecisionTYPE [segment_param->numb_classes*CurrSizes->numelmasked]();
    PrecisionTYPE * MRF = new PrecisionTYPE[segment_param->numb_classes*CurrSizes->numelmasked]();
    PrecisionTYPE * ShortPrior = new PrecisionTYPE[segment_param->numb_classes*CurrSizes->numelmasked]();
    for(int i=0; i<(CurrSizes->numelmasked*(segment_param->numb_classes)); i++){
        MRF[i]=1.0f;
        ShortPrior[i]=1.0f;
    }

    PrecisionTYPE * G = new PrecisionTYPE [segment_param->numb_classes*segment_param->numb_classes]();
    PrecisionTYPE * H = new PrecisionTYPE [segment_param->numb_classes*segment_param->numb_classes]();
    Create_diagonal_GH_Nclass(G,H,1,segment_param);

    float loglik=2;
    float oldloglik=1;
    // Iterative Components - EM, MRF
    //Expectation and LogLik
    bool out=true;
    int iter=1;
    while (out) {
        cout << "Iter = " << iter ;
        flush(cout);

        calcE_mask(CT,MRF,Expec,&loglik,NULL,S2L,M,V,CurrSizes,segment_param->verbose_level);

        //Maximization
        CurrSizes->numclass=2;
        calcM_mask(CT,Expec,NULL,S2L,M,V,NULL,NULL,CurrSizes,segment_param->verbose_level);
        CurrSizes->numclass=3;
        MRFregularization_mask(Expec,G,H,NULL,MRF,ShortPrior,L2S,S2L,CurrSizes,segment_param->flag_MRF,segment_param->verbose_level);
        //Bias Correction

        cout << "  \\ loglik = " << fabs((loglik-oldloglik)/fabs(oldloglik)) << endl;
        for(int classi=0;classi<segment_param->numb_classes;classi++){
            cout << "       M["<<classi<<"] = " <<M[classi]<<"  \\  V["<<classi<<"] = " <<V[classi] << endl;
        }
        if(fabs((loglik-oldloglik)/fabs(oldloglik))<(PrecisionTYPE)(0.005) && iter>3){
            out=false;
        }

        oldloglik=loglik;
        iter++;
    }


    for(int i=0; i<CT->nvox; i++)Mask_carr[i]=0;
    for(int i=0; i<CurrSizes->numelmasked; i++)Mask_carr[S2L[i]]=(bool)(Expec[i]>0.3);


    int * Buffer_Int_afterseg = new int [dimensions[0]*dimensions[1]*dimensions[2]]();
    for(int i=0; i<CT->nvox; i++)Buffer_Int_afterseg[i]=(int)(Mask_carr[i]);

    int * Buffer_Int2_afterseg = new int [dimensions[0]*dimensions[1]*dimensions[2]]();
    cout << "Finding biggest connected component";
    flush(cout);
    ConnectComp(Buffer_Int_afterseg,Buffer_Int2_afterseg,dimensions);
    for(int i=0; i<CT->nvox; i++)Mask_carr[i]=(Buffer_Int2_afterseg[i]==1);
    cout << " - Done" << endl;
    cout << "Adding Tagged Fluid";
    flush(cout);
    bool * Tagged = new bool [dimensions[0]*dimensions[1]*dimensions[2]]();
    for(int i=0; i<CT->nvox; i++)Tagged[i]=0;
    for(int i=0; i<CurrSizes->numelmasked; i++)Tagged[S2L[i]]=(bool)((Expec[i+2*CurrSizes->numelmasked]>Expec[i+CurrSizes->numelmasked]) & (Expec[i+2*CurrSizes->numelmasked]>0.1) );
    Dillate(Tagged,2,dimensions);
    for(int i=0; i<CT->nvox; i++)Tagged[i]=!Tagged[i];
    Dillate_const(Tagged,Mask_carr,2,dimensions);
    Dillate(Tagged,1,dimensions);
    for(int i=0; i<CT->nvox; i++)Tagged[i]=!Tagged[i];
    Dillate(Tagged,1,dimensions);
    cout << " - Done" << endl;
    flush(cout);

    for(int i=0; i<CT->nvox; i++){
        Buffer_Int_afterseg[i]=(int)(Mask_carr[i]||Tagged[i]);
        Buffer_Int2_afterseg[i]=0;
    }
    cout << "Finding biggest connected component with the tagged fluid";
    flush(cout);
    ConnectComp(Buffer_Int_afterseg,Buffer_Int2_afterseg,dimensions);
    for(int i=0; i<CT->nvox; i++)Mask_carr[i]=(Buffer_Int2_afterseg[i]==1);
    cout << " - Done" << endl;


    delete [] Expec;
    delete [] MRF;
    delete [] S2L;
    delete [] L2S;
    delete [] Buffer_Int_afterseg;
    delete [] Buffer_Int2_afterseg;

    nifti_image * Result = Copy_single_image_to_Result(Mask_carr,CT,filename_out);

    time(&end);
    cout << "Finished in "<<difftime(end,start)<<"sec"<< endl;



    return Result;
}
