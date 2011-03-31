#ifndef _SEG_BIASCORRECTION_CPP
#define _SEG_BIASCORRECTION_CPP


#include "_seg_BiasCorrection.h"

void BiasCorrection(PrecisionTYPE * BiasField,
                    PrecisionTYPE * BiasFieldCoefs,
                    nifti_image * T1,
                    PrecisionTYPE * Expec,
                    PrecisionTYPE * M,
                    PrecisionTYPE * V,
                    int biasOrder,
                    ImageSize * CurrSizes,
                    bool flag_Bias,
                    int verbose_level) {

    if(verbose_level>0){
        cout << "Optimising the Bias Field with order " << biasOrder<< endl;
        if(CurrSizes->usize>1){
            cout<< "Assuming fully decoupled bias-fields" << endl;
        }
        flush(cout);
    }
    int reduxfactor=redux_factor_for_bias;
    int nrOfClasses = CurrSizes->numclass;
    //nrOfClasses = 1;
    //int nrOfClasses = non_PV_numclass;
    int TotalLength = CurrSizes->numel;
    int UsedBasisFunctions=(int)((biasOrder+1) * (biasOrder+2)/2 *(biasOrder+3)/3);
    PrecisionTYPE * sampledData = static_cast<PrecisionTYPE *>(T1->data);


    // Precompute Powers depending on the current BiasOrder
    int PowerOrder [((maxallowedpowerorder+1)*(maxallowedpowerorder+2)/2*(maxallowedpowerorder+3))]={0};
    int ind=0;
    for(int order=0; order<=biasOrder; order++){
        for(int xorder=0; xorder<=order; xorder++){
            for(int yorder=0; yorder<=(order-xorder); yorder++){
                int zorder=order-yorder-xorder;
                PowerOrder[ind] =xorder;
                PowerOrder[ind+1] =yorder;
                PowerOrder[ind+2] =zorder;
                ind += 3;
            }
        }
    }
    PrecisionTYPE invV[nrOfClasses];
    PrecisionTYPE currM[nrOfClasses];

    for(int multispec=0; multispec<CurrSizes->usize; multispec++){
        sampledData = &sampledData[multispec*CurrSizes->numel];
        // Precompute the M and V  inverses

        for(int i=0; i<nrOfClasses; i++){
            invV[i]=1.0f/(V[i*CurrSizes->usize*CurrSizes->usize+multispec+multispec*CurrSizes->usize]);
            currM[i]=M[i*CurrSizes->usize+multispec];
        }



        PrecisionTYPE A [((maxallowedpowerorder+1)*(maxallowedpowerorder+2)/2*(maxallowedpowerorder+3)/3)*((maxallowedpowerorder+1)*(maxallowedpowerorder+2)/2*(maxallowedpowerorder+3)/3)]={0.0f};
        PrecisionTYPE B [((maxallowedpowerorder+1)*(maxallowedpowerorder+2)/2*(maxallowedpowerorder+3)/3)]={0.0f};
        PrecisionTYPE C [((maxallowedpowerorder+1)*(maxallowedpowerorder+2)/2*(maxallowedpowerorder+3)/3)]={0.0f};


        // Precompute sizes
        int col_size = (int)(CurrSizes->xsize);
        int plane_size = (int)(CurrSizes->xsize)*(CurrSizes->ysize);
        int maxix = (int)(CurrSizes->xsize);
        int maxiy = (int)(CurrSizes->ysize);
        int maxiz = (int)(CurrSizes->zsize);
        int Dims3d[3]={0};
        Dims3d[0]=maxix;
        Dims3d[1]=maxiy;
        Dims3d[2]=maxiz;

        // Precompute number of samples as it was never computed
        int samplecount=0;
        int linearindexes=0;
        int currindex=0;
        if(CurrSizes->numelbias==0){
            for (int iz=0; iz<maxiz; iz+=reduxfactor) {
                for (int iy=0; iy<maxiy; iy+=reduxfactor) {
                    for (int ix=0; ix<maxix; ix+=reduxfactor) {
                        samplecount++;
                    }
                }
            }
            if(verbose_level>0){
                cout << "Samplecount = " << samplecount<<"\n";
                flush(cout);
            }
            CurrSizes->numelbias=samplecount;
        }
        else{
            samplecount=CurrSizes->numelbias;
        }



        // CALC MATRIX A

        // Calc W (Van Leemput 1999 eq 7)

        PrecisionTYPE * Tempvar= new PrecisionTYPE [samplecount] ();
        PrecisionTYPE Tempvar_tmp=0;
        currindex=0;
        int tempvarindex=0;
        for (int iz=0; iz<maxiz; iz+=reduxfactor) {
            for (int iy=0; iy<maxiy; iy+=reduxfactor) {
                for (int ix=0; ix<maxix; ix+=reduxfactor) {
                    currindex=iz*plane_size+iy*col_size+ix;
                    Tempvar_tmp=0;
                    for(int j=0; j<nrOfClasses; j++){
                        Tempvar_tmp+=Expec[currindex+TotalLength*j]*invV[j];
                    }
                    Tempvar[tempvarindex]=Tempvar_tmp;
                    tempvarindex++;
                }
            }
        }

        // Precompute shifts
        PrecisionTYPE not_point_five_times_dims_x=(0.5f*(PrecisionTYPE)Dims3d[0]);
        PrecisionTYPE not_point_five_times_dims_y=(0.5f*(PrecisionTYPE)Dims3d[1]);
        PrecisionTYPE not_point_five_times_dims_z=(0.5f*(PrecisionTYPE)Dims3d[2]);

        PrecisionTYPE inv_not_point_five_times_dims_x=1.0f/(0.5f*(PrecisionTYPE)Dims3d[0]);
        PrecisionTYPE inv_not_point_five_times_dims_y=1.0f/(0.5f*(PrecisionTYPE)Dims3d[1]);
        PrecisionTYPE inv_not_point_five_times_dims_z=1.0f/(0.5f*(PrecisionTYPE)Dims3d[2]);



        PrecisionTYPE * Basis= new PrecisionTYPE[UsedBasisFunctions]();
        PrecisionTYPE xpos=0.0f;
        PrecisionTYPE ypos=0.0f;
        PrecisionTYPE zpos=0.0f;
        int x_bias_index_shift=0;
        int y_bias_index_shift=1;
        int z_bias_index_shift=2;
        PrecisionTYPE current_Tempvar=0.0f;

        // Calc A'WA (Van Leemput 1999 eq 7)
        tempvarindex=0;
        PrecisionTYPE * Basisptr1= (PrecisionTYPE *) Basis;
        PrecisionTYPE * Basisptr2= (PrecisionTYPE *) Basis;
        PrecisionTYPE * Aptr= (PrecisionTYPE *) A;
        for (int iz=0; iz<maxiz; iz+=reduxfactor) {
            for (int iy=0; iy<maxiy; iy+=reduxfactor) {
                for (int ix=0; ix<maxix; ix+=reduxfactor) {
                    linearindexes=(iz)*(CurrSizes->xsize)*(CurrSizes->ysize)+(iy)*(CurrSizes->xsize)+ix;

                    Basisptr1= (PrecisionTYPE *) Basis;
                    current_Tempvar=Tempvar[tempvarindex];
                    xpos=(((PrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                    ypos=(((PrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                    zpos=(((PrecisionTYPE)iz-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z);
                    x_bias_index_shift=0;
                    y_bias_index_shift=1;
                    z_bias_index_shift=2;
                    for(int j2=0; j2<UsedBasisFunctions; j2++,x_bias_index_shift+=3,y_bias_index_shift+=3,z_bias_index_shift+=3, Basisptr1++){
                        // Because Powerorder is alwais int, use a special power function (faster)
                        *Basisptr1=(pow_int(xpos,PowerOrder[x_bias_index_shift])*pow_int(ypos,PowerOrder[y_bias_index_shift])*pow_int(zpos,PowerOrder[z_bias_index_shift]));
                        // Instead, although slower, one can use
                        // TmpA=(pow(xpos,PowerOrder[0+j2*3])*pow(ypos,PowerOrder[1+j2*3])*pow(zpos,PowerOrder[2+j2*3]));
                    }
                    Basisptr1= (PrecisionTYPE *) Basis;
                    Aptr= (PrecisionTYPE *) A;
                    for(int j2=0; j2<UsedBasisFunctions; j2++, Basisptr1++){
                        Basisptr2= &Basis[j2];
                        Aptr= &A[j2+j2*UsedBasisFunctions];
                        for(int i2=j2; i2<UsedBasisFunctions; i2++, Aptr++, Basisptr2++){
                            (*Aptr)+=(*Basisptr2)*(current_Tempvar)*(*Basisptr1);
                        }
                    }
                    tempvarindex++;
                }
            }
        }



        matrix <double> RealA(UsedBasisFunctions,UsedBasisFunctions);

        for(int j2=0; j2<UsedBasisFunctions; j2++){
            for(int i2=j2; i2<UsedBasisFunctions; i2++){
                RealA.setvalue(i2,j2,(double)(A[i2+j2*UsedBasisFunctions]));
                RealA.setvalue(j2,i2,(double)(A[i2+j2*UsedBasisFunctions]));
            }
        }

        matrix <double> RealA_inv(UsedBasisFunctions);
        RealA_inv.copymatrix(RealA);
        RealA_inv.invert();

        if(verbose_level>1){
            matrix <double> RealA_test(UsedBasisFunctions);
            RealA_test.settoproduct(RealA,RealA_inv);
            RealA_test.comparetoidentity();
            //RealA.dumpmatrix();
        }



        // CALC MATRIX B

        //Precompute WR (Van Leemput 1999 eq 7)
        PrecisionTYPE Wi;
        PrecisionTYPE Wij;
        PrecisionTYPE Yest;
        PrecisionTYPE Ysum;
        tempvarindex=0;
        for (int iz=0; iz<maxiz; iz+=reduxfactor) {
            for (int iy=0; iy<maxiy; iy+=reduxfactor) {
                for (int ix=0; ix<maxix; ix+=reduxfactor) {
                    linearindexes=(iz)*(CurrSizes->xsize)*(CurrSizes->ysize)+(iy)*(CurrSizes->xsize)+ix;

                    Wi=0;
                    Wij=0;
                    Yest=0;
                    Ysum=0;
                    for(int j=0; j<nrOfClasses; j++){
                        PrecisionTYPE tmpexpec = (PrecisionTYPE)Expec[linearindexes+TotalLength*j];
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

        for(int i2=0; i2<UsedBasisFunctions; i2++){
            tempvarindex=0;
            B[i2]=0;
            for (int iz=0; iz<maxiz; iz+=reduxfactor) {
                for (int iy=0; iy<maxiy; iy+=reduxfactor) {
                    for (int ix=0; ix<maxix; ix+=reduxfactor) {
                        linearindexes=(iz)*(CurrSizes->xsize)*(CurrSizes->ysize)+(iy)*(CurrSizes->xsize)+ix;
                        B[i2]+=pow_int((((PrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x),PowerOrder[0+i2*3])*
                               pow_int((((PrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y),PowerOrder[1+i2*3])*
                               pow_int((((PrecisionTYPE)iz-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z),PowerOrder[2+i2*3])*
                               Tempvar[tempvarindex];
                        if(B[i2]!=B[i2]){
                            B[i2]=1;
                        }
                        tempvarindex++;
                    }
                }
            }
        }

        matrix <double> RealB(UsedBasisFunctions,1);

        for(int i2=0; i2<UsedBasisFunctions; i2++){
            RealB.setvalue(i2,0,(double)(B[i2]));
        }

        matrix <double> RealC(UsedBasisFunctions,1);

        RealC.settoproduct(RealA_inv,RealB);
        if(verbose_level>1){
            cout << "C= " << endl;
            RealC.dumpmatrix();
        }

        double cvalue;
        bool success;
        for(int i2=0; i2<UsedBasisFunctions; i2++){
            RealC.getvalue(i2,0,cvalue,success);
            C[i2]=(PrecisionTYPE)(cvalue);
        }

        PrecisionTYPE currxpower[maxallowedpowerorder];
        PrecisionTYPE currypower[maxallowedpowerorder];
        PrecisionTYPE currzpower[maxallowedpowerorder];
        PrecisionTYPE tmpbiasfield=0.0f;
        for (int iz=0; iz<maxiz; iz++) {
            for (int iy=0; iy<maxiy; iy++) {
                for (int ix=0; ix<maxix; ix++) {
                    linearindexes=(iz)*(CurrSizes->xsize)*(CurrSizes->ysize)+(iy)*(CurrSizes->xsize)+ix;
                    tmpbiasfield=0.0f;
                    xpos=(((PrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                    ypos=(((PrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                    zpos=(((PrecisionTYPE)iz-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z);
                    get_xyz_pow_int(xpos, ypos, zpos, currxpower, currypower, currzpower, biasOrder);
                    ind=0;
                    for(int order=0; order<=biasOrder; order++){
                        for(int xorder=0; xorder<=order; xorder++){
                            for(int yorder=0; yorder<=(order-xorder); yorder++){
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

        for( int i=0; i<UsedBasisFunctions; i++){
            BiasFieldCoefs[i+multispec*UsedBasisFunctions]=C[i];
        }

        delete [] Basis;
        delete [] Tempvar;
    }


}





void BiasCorrection_mask(PrecisionTYPE * BiasField,
                         PrecisionTYPE * BiasFieldCoefs,
                         nifti_image * T1,
                         int * Long_2_Short_Indices,
                         PrecisionTYPE * Expec,
                         PrecisionTYPE * M,
                         PrecisionTYPE * V,
                         int biasOrder,
                         ImageSize * CurrSizes,
                         bool flag_Bias,
                         int verbose_level) {

    if(verbose_level>0){
        cout << "Optimising the Bias Field with order " << biasOrder<< endl;
        if(CurrSizes->usize>1){
            cout<< "Assuming fully decoupled bias-fields" << endl;
        }
        flush(cout);
    }
    int reduxfactor=redux_factor_for_bias;
    int nrOfClasses = CurrSizes->numclass;
    //nrOfClasses = 1;
    //int nrOfClasses = non_PV_numclass;
    int TotalLength = CurrSizes->numelmasked;
    int UsedBasisFunctions=(int)((biasOrder+1) * (biasOrder+2)/2 *(biasOrder+3)/3);
    PrecisionTYPE * sampledData = static_cast<PrecisionTYPE *>(T1->data);


    // Precompute Powers depending on the current BiasOrder
    int PowerOrder [((maxallowedpowerorder+1)*(maxallowedpowerorder+2)/2*(maxallowedpowerorder+3))]={0};
    int ind=0;
    for(int order=0; order<=biasOrder; order++){
        for(int xorder=0; xorder<=order; xorder++){
            for(int yorder=0; yorder<=(order-xorder); yorder++){
                int zorder=order-yorder-xorder;
                PowerOrder[ind] =xorder;
                PowerOrder[ind+1] =yorder;
                PowerOrder[ind+2] =zorder;
                ind += 3;
            }
        }
    }
    PrecisionTYPE invV[nrOfClasses];
    PrecisionTYPE currM[nrOfClasses];

    for(int multispec=0; multispec<CurrSizes->usize; multispec++){
        sampledData = static_cast<PrecisionTYPE *>(T1->data);
        sampledData = &sampledData[multispec*CurrSizes->numel];
        // Precompute the M and V  inverses
        for(int i=0; i<nrOfClasses; i++){
            invV[i]=1.0f/(V[i*CurrSizes->usize*CurrSizes->usize+multispec+multispec*CurrSizes->usize]);
            currM[i]=M[i*CurrSizes->usize+multispec];
        }
        PrecisionTYPE A [((maxallowedpowerorder+1)*(maxallowedpowerorder+2)/2*(maxallowedpowerorder+3)/3)*((maxallowedpowerorder+1)*(maxallowedpowerorder+2)/2*(maxallowedpowerorder+3)/3)]={0.0f};
        PrecisionTYPE B [((maxallowedpowerorder+1)*(maxallowedpowerorder+2)/2*(maxallowedpowerorder+3)/3)]={0.0f};
        PrecisionTYPE C [((maxallowedpowerorder+1)*(maxallowedpowerorder+2)/2*(maxallowedpowerorder+3)/3)]={0.0f};


        // Precompute sizes
        int col_size = (int)(CurrSizes->xsize);
        int plane_size = (int)(CurrSizes->xsize)*(CurrSizes->ysize);
        int maxix = (int)(CurrSizes->xsize);
        int maxiy = (int)(CurrSizes->ysize);
        int maxiz = (int)(CurrSizes->zsize);
        int Dims3d[3]={0};
        Dims3d[0]=maxix;
        Dims3d[1]=maxiy;
        Dims3d[2]=maxiz;

        // Precompute number of samples as it was never computed
        int samplecount=0;
        int linearindexes=0;
        int currshortindex=0;
        if(CurrSizes->numelbias==0){
            for (int iz=0; iz<maxiz; iz+=reduxfactor) {
                for (int iy=0; iy<maxiy; iy+=reduxfactor) {
                    for (int ix=0; ix<maxix; ix+=reduxfactor) {
                        linearindexes=iz*plane_size+iy*col_size+ix;
                        if(Long_2_Short_Indices[linearindexes]>=0){
                            samplecount++;
                        }
                    }
                }
            }
            if(verbose_level>0){
                cout << "Number of samples for BiasField = " << samplecount<<"\n";
                flush(cout);
            }
            CurrSizes->numelbias=samplecount;
        }
        else{
            samplecount=CurrSizes->numelbias;
        }



        // CALC MATRIX A

        // Calc W (Van Leemput 1999 eq 7)
        //cout << "Calculating C = inv(A'WA) WR";
        //flush(cout);

        PrecisionTYPE * Tempvar= new PrecisionTYPE [samplecount] ();
        PrecisionTYPE Tempvar_tmp=0;
        currshortindex=0;
        int tempvarindex=0;
        for (int iz=0; iz<maxiz; iz+=reduxfactor) {
            for (int iy=0; iy<maxiy; iy+=reduxfactor) {
                for (int ix=0; ix<maxix; ix+=reduxfactor) {
                    currshortindex=Long_2_Short_Indices[iz*plane_size+iy*col_size+ix];
                    if(currshortindex>=0){
                        Tempvar_tmp=0;
                        for(int j=0; j<nrOfClasses; j++){
                            Tempvar_tmp+=Expec[currshortindex+TotalLength*j]*invV[j];
                        }
                        Tempvar[tempvarindex]=Tempvar_tmp;
                        tempvarindex++;
                    }
                }
            }
        }

        // Precompute shifts
        PrecisionTYPE not_point_five_times_dims_x=(0.5f*(PrecisionTYPE)Dims3d[0]);
        PrecisionTYPE not_point_five_times_dims_y=(0.5f*(PrecisionTYPE)Dims3d[1]);
        PrecisionTYPE not_point_five_times_dims_z=(0.5f*(PrecisionTYPE)Dims3d[2]);

        PrecisionTYPE inv_not_point_five_times_dims_x=1.0f/(0.5f*(PrecisionTYPE)Dims3d[0]);
        PrecisionTYPE inv_not_point_five_times_dims_y=1.0f/(0.5f*(PrecisionTYPE)Dims3d[1]);
        PrecisionTYPE inv_not_point_five_times_dims_z=1.0f/(0.5f*(PrecisionTYPE)Dims3d[2]);

        PrecisionTYPE * Basis= new PrecisionTYPE[UsedBasisFunctions]();
        PrecisionTYPE xpos=0.0f;
        PrecisionTYPE ypos=0.0f;
        PrecisionTYPE zpos=0.0f;
        int x_bias_index_shift=0;
        int y_bias_index_shift=1;
        int z_bias_index_shift=2;
        PrecisionTYPE current_Tempvar=0.0f;

        // Calc A'WA (Van Leemput 1999 eq 7)
        tempvarindex=0;
        PrecisionTYPE * Basisptr1= (PrecisionTYPE *) Basis;
        PrecisionTYPE * Basisptr2= (PrecisionTYPE *) Basis;
        PrecisionTYPE * Aptr= (PrecisionTYPE *) A;
        for (int iz=0; iz<maxiz; iz+=reduxfactor) {
            for (int iy=0; iy<maxiy; iy+=reduxfactor) {
                for (int ix=0; ix<maxix; ix+=reduxfactor) {
                    linearindexes=(iz)*(CurrSizes->xsize)*(CurrSizes->ysize)+(iy)*(CurrSizes->xsize)+ix;
                    currshortindex=Long_2_Short_Indices[linearindexes];
                    if(currshortindex>=0){
                        Basisptr1= (PrecisionTYPE *) Basis;
                        current_Tempvar=Tempvar[tempvarindex];
                        xpos=(((PrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                        ypos=(((PrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                        zpos=(((PrecisionTYPE)iz-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z);
                        x_bias_index_shift=0;
                        y_bias_index_shift=1;
                        z_bias_index_shift=2;
                        for(int j2=0; j2<UsedBasisFunctions; j2++,x_bias_index_shift+=3,y_bias_index_shift+=3,z_bias_index_shift+=3, Basisptr1++){
                            // Because Powerorder is always int, use a special power function (faster)
                            *Basisptr1=(pow_int(xpos,PowerOrder[x_bias_index_shift])*pow_int(ypos,PowerOrder[y_bias_index_shift])*pow_int(zpos,PowerOrder[z_bias_index_shift]));
                        }


                        Basisptr1= (PrecisionTYPE *) Basis;
                        Aptr= (PrecisionTYPE *) A;
                        for(int j2=0; j2<UsedBasisFunctions; j2++, Basisptr1++){
                            Basisptr2= &Basis[j2];
                            Aptr= &A[j2+j2*UsedBasisFunctions];
                            for(int i2=j2; i2<UsedBasisFunctions; i2++, Aptr++, Basisptr2++){
                                (*Aptr)+=(*Basisptr2)*(current_Tempvar)*(*Basisptr1);
                            }
                        }
                        tempvarindex++;
                    }
                }
            }
        }



        matrix <double> RealA(UsedBasisFunctions,UsedBasisFunctions);

        for(int j2=0; j2<UsedBasisFunctions; j2++){
            for(int i2=j2; i2<UsedBasisFunctions; i2++){
                RealA.setvalue(i2,j2,(double)(A[i2+j2*UsedBasisFunctions]));
                RealA.setvalue(j2,i2,(double)(A[i2+j2*UsedBasisFunctions]));
            }
        }

        matrix <double> RealA_inv(UsedBasisFunctions);
        RealA_inv.copymatrix(RealA);
        RealA_inv.invert();

        if(verbose_level>1){
            matrix <double> RealA_test(UsedBasisFunctions);
            RealA_test.settoproduct(RealA,RealA_inv);
            RealA_test.comparetoidentity();
            //RealA.dumpmatrix();
        }



        // CALC MATRIX B

        //Precompute WR (Van Leemput 1999 eq 7)
        PrecisionTYPE Wi;
        PrecisionTYPE Wij;
        PrecisionTYPE Yest;
        PrecisionTYPE Ysum;
        tempvarindex=0;
        for (int iz=0; iz<maxiz; iz+=reduxfactor) {
            for (int iy=0; iy<maxiy; iy+=reduxfactor) {
                for (int ix=0; ix<maxix; ix+=reduxfactor) {
                    linearindexes=(iz)*(CurrSizes->xsize)*(CurrSizes->ysize)+(iy)*(CurrSizes->xsize)+ix;
                    currshortindex=Long_2_Short_Indices[linearindexes];
                    if(currshortindex>=0){
                        Wi=0;
                        Wij=0;
                        Yest=0;
                        Ysum=0;
                        for(int j=0; j<nrOfClasses; j++){
                            PrecisionTYPE tmpexpec = (PrecisionTYPE)Expec[currshortindex+TotalLength*j];
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

        for(int i2=0; i2<UsedBasisFunctions; i2++){
            tempvarindex=0;
            B[i2]=0;
            for (int iz=0; iz<maxiz; iz+=reduxfactor) {
                for (int iy=0; iy<maxiy; iy+=reduxfactor) {
                    for (int ix=0; ix<maxix; ix+=reduxfactor) {
                        linearindexes=(iz)*(CurrSizes->xsize)*(CurrSizes->ysize)+(iy)*(CurrSizes->xsize)+ix;
                        currshortindex=Long_2_Short_Indices[linearindexes];
                        if(currshortindex>=0){

                            B[i2]+=pow_int((((PrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x),PowerOrder[0+i2*3])*
                                   pow_int((((PrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y),PowerOrder[1+i2*3])*
                                   pow_int((((PrecisionTYPE)iz-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z),PowerOrder[2+i2*3])*
                                   Tempvar[tempvarindex];

                            if(B[i2]!=B[i2]){
                                B[i2]=1;
                            }
                            tempvarindex++;
                        }
                    }
                }
            }
        }

        matrix <double> RealB(UsedBasisFunctions,1);

        for(int i2=0; i2<UsedBasisFunctions; i2++){
            RealB.setvalue(i2,0,(double)(B[i2]));
        }

        matrix <double> RealC(UsedBasisFunctions,1);

        RealC.settoproduct(RealA_inv,RealB);
        if(verbose_level>1){
            cout << "C= " << endl;
            RealC.dumpmatrix();
        }

        double cvalue;
        bool success;
        for(int i2=0; i2<UsedBasisFunctions; i2++){
            RealC.getvalue(i2,0,cvalue,success);
            C[i2]=(PrecisionTYPE)(cvalue);
        }

        PrecisionTYPE currxpower[maxallowedpowerorder];
        PrecisionTYPE currypower[maxallowedpowerorder];
        PrecisionTYPE currzpower[maxallowedpowerorder];
        PrecisionTYPE tmpbiasfield=0.0f;
        for (int iz=0; iz<maxiz; iz++) {
            for (int iy=0; iy<maxiy; iy++) {
                for (int ix=0; ix<maxix; ix++) {
                    linearindexes=(iz)*(CurrSizes->xsize)*(CurrSizes->ysize)+(iy)*(CurrSizes->xsize)+ix;
                    currshortindex=Long_2_Short_Indices[linearindexes];
                    if(currshortindex>=0){
                        tmpbiasfield=0.0f;
                        xpos=(((PrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                        ypos=(((PrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                        zpos=(((PrecisionTYPE)iz-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z);
                        get_xyz_pow_int(xpos, ypos, zpos, currxpower, currypower, currzpower, biasOrder);
                        ind=0;
                        for(int order=0; order<=biasOrder; order++){
                            for(int xorder=0; xorder<=order; xorder++){
                                for(int yorder=0; yorder<=(order-xorder); yorder++){
                                    int zorder=order-yorder-xorder;
                                    tmpbiasfield-=C[ind]*currxpower[xorder]*currypower[yorder]*currzpower[zorder];
                                    ind++;
                                }
                            }
                        }
                        BiasField[currshortindex+multispec*CurrSizes->numelmasked]=tmpbiasfield;
                    }
                }
            }
        }

        for( int i=0; i<UsedBasisFunctions; i++){
            BiasFieldCoefs[i+multispec*UsedBasisFunctions]=C[i];
        }
        delete [] Basis;
        delete [] Tempvar;
    }
}

void get_xyz_pow_int(PrecisionTYPE xpos,
                     PrecisionTYPE ypos,
                     PrecisionTYPE zpos,
                     PrecisionTYPE *currxpower,
                     PrecisionTYPE *currypower,
                     PrecisionTYPE *currzpower,
                     int maxorder){
    int order=1;
    currxpower[0]=1;
    currypower[0]=1;
    currzpower[0]=1;
    int orderminusone=0;
    int maxorderplusone=maxorder+1;
    while (order<maxorderplusone){
        currxpower[order]=currxpower[orderminusone]*xpos;
        currypower[order]=currypower[orderminusone]*ypos;
        currzpower[order]=currzpower[orderminusone]*zpos;
        order++;
        orderminusone++;
    }
    return;
}



void BiasCorrection2D(PrecisionTYPE * BiasField,
                      PrecisionTYPE * BiasFieldCoefs,
                      nifti_image * T1,
                      PrecisionTYPE * Expec,
                      PrecisionTYPE * M,
                      PrecisionTYPE * V,
                      int biasOrder,
                      ImageSize * CurrSizes,
                      bool flag_Bias,
                      int verbose_level) {

    if(verbose_level>0){
        cout << "Optimising the Bias Field with order " << biasOrder<< endl;
        if(CurrSizes->usize>1){
            cout<< "Assuming fully decoupled bias-fields" << endl;
        }
        flush(cout);
    }
    int reduxfactor=redux_factor_for_bias;
    int nrOfClasses = CurrSizes->numclass;

    int TotalLength = CurrSizes->numel;
    int UsedBasisFunctions=(int)((biasOrder+1) * (biasOrder+2)/2);
    PrecisionTYPE * sampledData = static_cast<PrecisionTYPE *>(T1->data);


    // Precompute Powers depending on the current BiasOrder
    int PowerOrder [((maxallowedpowerorder+1)*(maxallowedpowerorder+2)/2)]={0};
    int ind=0;
    for(int order=0; order<=biasOrder; order++){
        for(int xorder=0; xorder<=order; xorder++){
            for(int yorder=0; yorder<=(order-xorder); yorder++){
                int zorder=order-yorder-xorder;
                PowerOrder[ind] =xorder;
                PowerOrder[ind+1] =yorder;
                ind += 3;
            }
        }
    }
    PrecisionTYPE invV[nrOfClasses];
    PrecisionTYPE currM[nrOfClasses];

    for(int multispec=0; multispec<CurrSizes->usize; multispec++){
        sampledData = &sampledData[multispec*CurrSizes->numel];
        // Precompute the M and V  inverses

        for(int i=0; i<nrOfClasses; i++){
            invV[i]=1.0f/(V[i*CurrSizes->usize*CurrSizes->usize+multispec+multispec*CurrSizes->usize]);
            currM[i]=M[i*CurrSizes->usize+multispec];
        }



        PrecisionTYPE A [((maxallowedpowerorder+1)*(maxallowedpowerorder+2)/2)*((maxallowedpowerorder+1)*(maxallowedpowerorder+2)/2)]={0.0f};
        PrecisionTYPE B [((maxallowedpowerorder+1)*(maxallowedpowerorder+2)/2)]={0.0f};
        PrecisionTYPE C [((maxallowedpowerorder+1)*(maxallowedpowerorder+2)/2)]={0.0f};


        // Precompute sizes
        int col_size = (int)(CurrSizes->xsize);
        int maxix = (int)(CurrSizes->xsize);
        int maxiy = (int)(CurrSizes->ysize);
        int Dims3d[2]={0};
        Dims3d[0]=maxix;
        Dims3d[1]=maxiy;

        // Precompute number of samples as it was never computed
        int samplecount=0;
        int linearindexes=0;
        int currindex=0;
        if(CurrSizes->numelbias==0){
                for (int iy=0; iy<maxiy; iy+=reduxfactor) {
                    for (int ix=0; ix<maxix; ix+=reduxfactor) {
                        samplecount++;
                    }
                }
            if(verbose_level>0){
                cout << "Samplecount = " << samplecount<<"\n";
                flush(cout);
            }
            CurrSizes->numelbias=samplecount;
        }
        else{
            samplecount=CurrSizes->numelbias;
        }



        // CALC MATRIX A

        // Calc W (Van Leemput 1999 eq 7)

        PrecisionTYPE * Tempvar= new PrecisionTYPE [samplecount] ();
        PrecisionTYPE Tempvar_tmp=0;
        currindex=0;
        int tempvarindex=0;
        for (int iy=0; iy<maxiy; iy+=reduxfactor) {
            for (int ix=0; ix<maxix; ix+=reduxfactor) {
                currindex=iy*col_size+ix;
                Tempvar_tmp=0;
                for(int j=0; j<nrOfClasses; j++){
                    Tempvar_tmp+=Expec[currindex+TotalLength*j]*invV[j];
                }
                Tempvar[tempvarindex]=Tempvar_tmp;
                tempvarindex++;
            }
        }


        // Precompute shifts
        PrecisionTYPE not_point_five_times_dims_x=(0.5f*(PrecisionTYPE)Dims3d[0]);
        PrecisionTYPE not_point_five_times_dims_y=(0.5f*(PrecisionTYPE)Dims3d[1]);

        PrecisionTYPE inv_not_point_five_times_dims_x=1.0f/(0.5f*(PrecisionTYPE)Dims3d[0]);
        PrecisionTYPE inv_not_point_five_times_dims_y=1.0f/(0.5f*(PrecisionTYPE)Dims3d[1]);


        PrecisionTYPE * Basis= new PrecisionTYPE[UsedBasisFunctions]();
        PrecisionTYPE xpos=0.0f;
        PrecisionTYPE ypos=0.0f;
        int x_bias_index_shift=0;
        int y_bias_index_shift=1;
        PrecisionTYPE current_Tempvar=0.0f;

        // Calc A'WA (Van Leemput 1999 eq 7)
        tempvarindex=0;
        PrecisionTYPE * Basisptr1= (PrecisionTYPE *) Basis;
        PrecisionTYPE * Basisptr2= (PrecisionTYPE *) Basis;
        PrecisionTYPE * Aptr= (PrecisionTYPE *) A;
        for (int iy=0; iy<maxiy; iy+=reduxfactor) {
            for (int ix=0; ix<maxix; ix+=reduxfactor) {
                linearindexes=(iy)*(CurrSizes->xsize)+ix;
                Basisptr1= (PrecisionTYPE *) Basis;
                current_Tempvar=Tempvar[tempvarindex];
                xpos=(((PrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                ypos=(((PrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                x_bias_index_shift=0;
                y_bias_index_shift=1;
                for(int j2=0; j2<UsedBasisFunctions; j2++,x_bias_index_shift+=2,y_bias_index_shift+=2, Basisptr1++){
                    // Because Powerorder is alwais int, use a special power function (faster)
                    *Basisptr1=(pow_int(xpos,PowerOrder[x_bias_index_shift])*pow_int(ypos,PowerOrder[y_bias_index_shift]));
                    // Instead, although slower, one can use
                    // TmpA=(pow(xpos,PowerOrder[0+j2*3])*pow(ypos,PowerOrder[1+j2*3])*pow(zpos,PowerOrder[2+j2*3]));
                }
                Basisptr1= (PrecisionTYPE *) Basis;
                Aptr= (PrecisionTYPE *) A;
                for(int j2=0; j2<UsedBasisFunctions; j2++, Basisptr1++){
                    Basisptr2= &Basis[j2];
                    Aptr= &A[j2+j2*UsedBasisFunctions];
                    for(int i2=j2; i2<UsedBasisFunctions; i2++, Aptr++, Basisptr2++){
                        (*Aptr)+=(*Basisptr2)*(current_Tempvar)*(*Basisptr1);
                    }
                }
                tempvarindex++;
            }
        }




    matrix <double> RealA(UsedBasisFunctions,UsedBasisFunctions);

    for(int j2=0; j2<UsedBasisFunctions; j2++){
        for(int i2=j2; i2<UsedBasisFunctions; i2++){
            RealA.setvalue(i2,j2,(double)(A[i2+j2*UsedBasisFunctions]));
            RealA.setvalue(j2,i2,(double)(A[i2+j2*UsedBasisFunctions]));
        }
    }

    matrix <double> RealA_inv(UsedBasisFunctions);
    RealA_inv.copymatrix(RealA);
    RealA_inv.invert();

    if(verbose_level>1){
        matrix <double> RealA_test(UsedBasisFunctions);
        RealA_test.settoproduct(RealA,RealA_inv);
        RealA_test.comparetoidentity();
        //RealA.dumpmatrix();
    }



    // CALC MATRIX B

    //Precompute WR (Van Leemput 1999 eq 7)
    PrecisionTYPE Wi;
    PrecisionTYPE Wij;
    PrecisionTYPE Yest;
    PrecisionTYPE Ysum;
    tempvarindex=0;
    for (int iy=0; iy<maxiy; iy+=reduxfactor) {
        for (int ix=0; ix<maxix; ix+=reduxfactor) {
            linearindexes=(iy)*(CurrSizes->xsize)+ix;

            Wi=0;
            Wij=0;
            Yest=0;
            Ysum=0;
            for(int j=0; j<nrOfClasses; j++){
                PrecisionTYPE tmpexpec = (PrecisionTYPE)Expec[linearindexes+TotalLength*j];
                Wij=tmpexpec*(invV[j]);
                Wi+=Wij;
                Yest+=Wij*(currM[j]);
                Ysum+=Wij;
            }
            Tempvar[tempvarindex]=Wi*(sampledData[linearindexes]-(Yest/Ysum));
            tempvarindex++;

        }
    }

    for(int i2=0; i2<UsedBasisFunctions; i2++){
        tempvarindex=0;
        B[i2]=0;
        for (int iy=0; iy<maxiy; iy+=reduxfactor) {
            for (int ix=0; ix<maxix; ix+=reduxfactor) {
                linearindexes=(iy)*(CurrSizes->xsize)+ix;
                B[i2]+=pow_int((((PrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x),PowerOrder[0+i2*3])*
                       pow_int((((PrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y),PowerOrder[1+i2*3])*
                       Tempvar[tempvarindex];
                if(B[i2]!=B[i2]){
                    B[i2]=1;
                }
                tempvarindex++;
            }
        }
    }

    matrix <double> RealB(UsedBasisFunctions,1);

    for(int i2=0; i2<UsedBasisFunctions; i2++){
        RealB.setvalue(i2,0,(double)(B[i2]));
    }

    matrix <double> RealC(UsedBasisFunctions,1);

    RealC.settoproduct(RealA_inv,RealB);
    if(verbose_level>1){
        cout << "C= " << endl;
        RealC.dumpmatrix();
    }

    double cvalue;
    bool success;
    for(int i2=0; i2<UsedBasisFunctions; i2++){
        RealC.getvalue(i2,0,cvalue,success);
        C[i2]=(PrecisionTYPE)(cvalue);
    }

    PrecisionTYPE currxpower[maxallowedpowerorder];
    PrecisionTYPE currypower[maxallowedpowerorder];
    PrecisionTYPE tmpbiasfield=0.0f;
    for (int iy=0; iy<maxiy; iy++) {
        for (int ix=0; ix<maxix; ix++) {
            linearindexes=(iy)*(CurrSizes->xsize)+ix;
            tmpbiasfield=0.0f;
            xpos=(((PrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
            ypos=(((PrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
            get_xy_pow_int(xpos, ypos, currxpower, currypower, biasOrder);
            ind=0;
            for(int order=0; order<=biasOrder; order++){
                for(int xorder=0; xorder<=order; xorder++){
                    int yorder=order-xorder;
                    tmpbiasfield-=C[ind]*currxpower[xorder]*currypower[yorder];
                    ind++;
                }
            }
            BiasField[linearindexes+multispec*CurrSizes->numel]=tmpbiasfield;
        }
    }

    for( int i=0; i<UsedBasisFunctions; i++){
        BiasFieldCoefs[i+multispec*UsedBasisFunctions]=C[i];
    }

    delete [] Basis;
    delete [] Tempvar;
}


}





void BiasCorrection_mask2D(PrecisionTYPE * BiasField,
                           PrecisionTYPE * BiasFieldCoefs,
                           nifti_image * T1,
                           int * Long_2_Short_Indices,
                           PrecisionTYPE * Expec,
                           PrecisionTYPE * M,
                           PrecisionTYPE * V,
                           int biasOrder,
                           ImageSize * CurrSizes,
                           bool flag_Bias,
                           int verbose_level)
{
    if(verbose_level>0){
        cout << "Optimising the Bias Field with order " << biasOrder<< endl;
        if(CurrSizes->usize>1){
            cout<< "Assuming fully decoupled bias-fields" << endl;
        }
        flush(cout);
    }
    int reduxfactor=redux_factor_for_bias;
    int nrOfClasses = CurrSizes->numclass;
    //nrOfClasses = 1;
    //int nrOfClasses = non_PV_numclass;
    int TotalLength = CurrSizes->numelmasked;
    int UsedBasisFunctions=(int)((biasOrder+1) * (biasOrder+2)/2);
    PrecisionTYPE * sampledData = static_cast<PrecisionTYPE *>(T1->data);


    // Precompute Powers depending on the current BiasOrder
    int PowerOrder [((maxallowedpowerorder+1)*(maxallowedpowerorder+2)/2)]={0};
    int ind=0;
    for(int order=0; order<=biasOrder; order++){
        for(int xorder=0; xorder<=order; xorder++){
            int yorder=order-xorder;
            PowerOrder[ind] =xorder;
            PowerOrder[ind+1] =yorder;
            ind += 2;
        }
    }
    PrecisionTYPE invV[nrOfClasses];
    PrecisionTYPE currM[nrOfClasses];

    for(int multispec=0; multispec<CurrSizes->usize; multispec++){
        sampledData = static_cast<PrecisionTYPE *>(T1->data);
        sampledData = &sampledData[multispec*CurrSizes->numel];
        // Precompute the M and V  inverses
        for(int i=0; i<nrOfClasses; i++){
            invV[i]=1.0f/(V[i*CurrSizes->usize*CurrSizes->usize+multispec+multispec*CurrSizes->usize]);
            currM[i]=M[i*CurrSizes->usize+multispec];
        }
        PrecisionTYPE A [((maxallowedpowerorder+1)*(maxallowedpowerorder+2))/2*((maxallowedpowerorder+1)*(maxallowedpowerorder+2)/2)]={0.0f};
        PrecisionTYPE B [((maxallowedpowerorder+1)*(maxallowedpowerorder+2))/2]={0.0f};
        PrecisionTYPE C [((maxallowedpowerorder+1)*(maxallowedpowerorder+2))/2]={0.0f};


        // Precompute sizes
        int col_size = (int)(CurrSizes->xsize);
        int maxix = (int)(CurrSizes->xsize);
        int maxiy = (int)(CurrSizes->ysize);
        int Dims3d[2]={0};
        Dims3d[0]=maxix;
        Dims3d[1]=maxiy;

        // Precompute number of samples as it was never computed
        int samplecount=0;
        int linearindexes=0;
        int currshortindex=0;
        if(CurrSizes->numelbias==0){
            for (int iy=0; iy<maxiy; iy+=reduxfactor) {
                for (int ix=0; ix<maxix; ix+=reduxfactor) {
                    linearindexes=iy*col_size+ix;
                    if(Long_2_Short_Indices[linearindexes]>=0){
                        samplecount++;
                    }
                }
            }
            if(verbose_level>0){
                cout << "Number of samples for BiasField = " << samplecount<<"\n";
                flush(cout);
            }
            CurrSizes->numelbias=samplecount;
        }
        else{
            samplecount=CurrSizes->numelbias;
        }



        // CALC MATRIX A

        // Calc W (Van Leemput 1999 eq 7)
        //cout << "Calculating C = inv(A'WA) WR";
        //flush(cout);
        PrecisionTYPE * Tempvar= new PrecisionTYPE [samplecount] ();
        PrecisionTYPE Tempvar_tmp=0;
        currshortindex=0;
        int tempvarindex=0;
        for (int iy=0; iy<maxiy; iy+=reduxfactor) {
            for (int ix=0; ix<maxix; ix+=reduxfactor) {
                currshortindex=Long_2_Short_Indices[iy*col_size+ix];
                if(currshortindex>=0){
                    Tempvar_tmp=0;
                    for(int j=0; j<nrOfClasses; j++){
                        Tempvar_tmp+=Expec[currshortindex+TotalLength*j]*invV[j];
                    }
                    Tempvar[tempvarindex]=Tempvar_tmp;
                    tempvarindex++;
                }
            }
        }
        // Precompute shifts
        PrecisionTYPE not_point_five_times_dims_x=(0.5f*(PrecisionTYPE)Dims3d[0]);
        PrecisionTYPE not_point_five_times_dims_y=(0.5f*(PrecisionTYPE)Dims3d[1]);

        PrecisionTYPE inv_not_point_five_times_dims_x=1.0f/(0.5f*(PrecisionTYPE)Dims3d[0]);
        PrecisionTYPE inv_not_point_five_times_dims_y=1.0f/(0.5f*(PrecisionTYPE)Dims3d[1]);


        PrecisionTYPE * Basis= new PrecisionTYPE[UsedBasisFunctions]();
        PrecisionTYPE xpos=0.0f;
        PrecisionTYPE ypos=0.0f;
        int x_bias_index_shift=0;
        int y_bias_index_shift=1;
        PrecisionTYPE current_Tempvar=0.0f;
        // Calc A'WA (Van Leemput 1999 eq 7)
        tempvarindex=0;
        PrecisionTYPE * Basisptr1= (PrecisionTYPE *) Basis;
        PrecisionTYPE * Basisptr2= (PrecisionTYPE *) Basis;
        PrecisionTYPE * Aptr= (PrecisionTYPE *) A;
        for (int iy=0; iy<maxiy; iy+=reduxfactor) {
            for (int ix=0; ix<maxix; ix+=reduxfactor) {
                linearindexes=(iy)*(CurrSizes->xsize)+ix;
                currshortindex=Long_2_Short_Indices[linearindexes];
                if(currshortindex>=0){
                    Basisptr1= (PrecisionTYPE *) Basis;
                    current_Tempvar=Tempvar[tempvarindex];
                    xpos=(((PrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                    ypos=(((PrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                    x_bias_index_shift=0;
                    y_bias_index_shift=1;
                    for(int j2=0; j2<UsedBasisFunctions; j2++,x_bias_index_shift+=2,y_bias_index_shift+=2, Basisptr1++){
                        // Because Powerorder is always int, use a special power function (faster)
                        *Basisptr1=(pow_int(xpos,PowerOrder[x_bias_index_shift])*pow_int(ypos,PowerOrder[y_bias_index_shift]));
                    }
                    Basisptr1= (PrecisionTYPE *) Basis;
                    Aptr= (PrecisionTYPE *) A;
                    for(int j2=0; j2<UsedBasisFunctions; j2++, Basisptr1++){
                        Basisptr2= &Basis[j2];
                        Aptr= &A[j2+j2*UsedBasisFunctions];
                        for(int i2=j2; i2<UsedBasisFunctions; i2++, Aptr++, Basisptr2++){
                            (*Aptr)+=(*Basisptr2)*(current_Tempvar)*(*Basisptr1);
                        }
                    }
                    tempvarindex++;
                }
            }
        }

        matrix <double> RealA(UsedBasisFunctions,UsedBasisFunctions);

        for(int j2=0; j2<UsedBasisFunctions; j2++){
            for(int i2=j2; i2<UsedBasisFunctions; i2++){
                RealA.setvalue(i2,j2,(double)(A[i2+j2*UsedBasisFunctions]));
                RealA.setvalue(j2,i2,(double)(A[i2+j2*UsedBasisFunctions]));
            }
        }

        matrix <double> RealA_inv(UsedBasisFunctions);
        RealA_inv.copymatrix(RealA);
        RealA_inv.invert();

        if(verbose_level>1){
            matrix <double> RealA_test(UsedBasisFunctions);
            RealA_test.settoproduct(RealA,RealA_inv);
            RealA_test.comparetoidentity();
            //RealA.dumpmatrix();
        }



        // CALC MATRIX B

        //Precompute WR (Van Leemput 1999 eq 7)
        PrecisionTYPE Wi;
        PrecisionTYPE Wij;
        PrecisionTYPE Yest;
        PrecisionTYPE Ysum;
        tempvarindex=0;
        for (int iy=0; iy<maxiy; iy+=reduxfactor) {
            for (int ix=0; ix<maxix; ix+=reduxfactor) {
                linearindexes=(iy)*(CurrSizes->xsize)+ix;
                currshortindex=Long_2_Short_Indices[linearindexes];
                if(currshortindex>=0){
                    Wi=0;
                    Wij=0;
                    Yest=0;
                    Ysum=0;
                    for(int j=0; j<nrOfClasses; j++){
                        PrecisionTYPE tmpexpec = (PrecisionTYPE)Expec[currshortindex+TotalLength*j];
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

        for(int i2=0; i2<UsedBasisFunctions; i2++){
            tempvarindex=0;
            B[i2]=0;
            for (int iy=0; iy<maxiy; iy+=reduxfactor) {

                for (int ix=0; ix<maxix; ix+=reduxfactor) {
                    linearindexes=(iy)*(CurrSizes->xsize)+ix;
                    currshortindex=Long_2_Short_Indices[linearindexes];
                    if(currshortindex>=0){
                        B[i2]+=pow_int((((PrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x),PowerOrder[0+i2*2])*
                               pow_int((((PrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y),PowerOrder[1+i2*2])*
                               Tempvar[tempvarindex];

                        if(B[i2]!=B[i2]){
                            B[i2]=1;
                        }
                        tempvarindex++;
                    }
                }
            }
        }

        matrix <double> RealB(UsedBasisFunctions,1);

        for(int i2=0; i2<UsedBasisFunctions; i2++){
            RealB.setvalue(i2,0,(double)(B[i2]));
        }

        matrix <double> RealC(UsedBasisFunctions,1);

        RealC.settoproduct(RealA_inv,RealB);
        if(verbose_level>1){
            cout << "C= " << endl;
            RealC.dumpmatrix();
        }

        double cvalue;
        bool success;
        for(int i2=0; i2<UsedBasisFunctions; i2++){
            RealC.getvalue(i2,0,cvalue,success);
            C[i2]=(PrecisionTYPE)(cvalue);
        }

        PrecisionTYPE currxpower[maxallowedpowerorder];
        PrecisionTYPE currypower[maxallowedpowerorder];
        PrecisionTYPE tmpbiasfield=0.0f;
        for (int iy=0; iy<maxiy; iy++) {
            for (int ix=0; ix<maxix; ix++) {
                linearindexes=(iy)*(CurrSizes->xsize)+ix;
                currshortindex=Long_2_Short_Indices[linearindexes];
                if(currshortindex>=0){
                    tmpbiasfield=0.0f;
                    xpos=(((PrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                    ypos=(((PrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                    get_xy_pow_int(xpos, ypos, currxpower, currypower, biasOrder);
                    ind=0;
                    for(int order=0; order<=biasOrder; order++){
                        for(int xorder=0; xorder<=order; xorder++){
                            int yorder=order-xorder;
                            tmpbiasfield-=C[ind]*currxpower[xorder]*currypower[yorder];
                            ind++;
                        }
                    }
                    BiasField[currshortindex+multispec*CurrSizes->numelmasked]=tmpbiasfield;
                }
            }

        }

        for( int i=0; i<UsedBasisFunctions; i++){
            BiasFieldCoefs[i+multispec*UsedBasisFunctions]=C[i];
        }
        delete [] Basis;
        delete [] Tempvar;
    }
}


inline PrecisionTYPE pow_int(const PrecisionTYPE base,
                             int exp){
    if(exp==0){return 1;}
    PrecisionTYPE result = base;
    while (--exp){result *= base;}
    return result;
}


void get_xy_pow_int(PrecisionTYPE xpos,
                    PrecisionTYPE ypos,
                    PrecisionTYPE *currxpower,
                    PrecisionTYPE *currypower,
                    int maxorder){
    int order=1;
    currxpower[0]=1;
    currypower[0]=1;
    int orderminusone=0;
    int maxorderplusone=maxorder+1;
    while (order<maxorderplusone){
        currxpower[order]=currxpower[orderminusone]*xpos;
        currypower[order]=currypower[orderminusone]*ypos;
        order++;
        orderminusone++;
    }
    return;
}

#endif
