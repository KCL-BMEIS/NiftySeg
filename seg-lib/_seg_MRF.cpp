#include "_seg_MRF.h"


void MRFregularization_mask(const PrecisionTYPE * Expec,
                            const PrecisionTYPE * G,
                            const PrecisionTYPE * H,
                            PrecisionTYPE * MRFbeta,
                            PrecisionTYPE * MRFprior,
                            PrecisionTYPE * AtlasPrior,
                            int * Long_2_Short_Indices,
                            int * Short_2_Long_Indices,
                            ImageSize * CurrSizes,
                            bool MRFflag,
                            int verbose_level)
{

    int numelmasked=CurrSizes->numelmasked;
    int numclass=CurrSizes->numclass;
    int image_size = (int)(CurrSizes->xsize)*(CurrSizes->ysize)*(CurrSizes->zsize);
    if(MRFflag){
        PrecisionTYPE * MRFpriorPtr = (PrecisionTYPE *)MRFprior;
        int * Long_2_Short_IndicesPtr = (int *)Long_2_Short_Indices;
        int col_size, plane_size, indexCentre, indexWest, indexEast, indexSouth, indexNorth, indexTop, indexBottom;
        int ix, iy, iz,maxiy, maxix, maxiz, neighbourclass;
        PrecisionTYPE Sum_Temp_MRF_Class_Expect;
        col_size = (int)(CurrSizes->xsize);
        plane_size = (int)(CurrSizes->xsize)*(CurrSizes->ysize);
        image_size = (int)(CurrSizes->xsize)*(CurrSizes->ysize)*(CurrSizes->zsize);
        maxix = (int)(CurrSizes->xsize);
        maxiy = (int)(CurrSizes->ysize);
        maxiz = (int)(CurrSizes->zsize);
        PrecisionTYPE Gplane[max_numbclass];
        PrecisionTYPE Hplane[max_numbclass];
        PrecisionTYPE Temp_MRF_Class_Expect[max_numbclass];
        if(verbose_level>0){
            cout << "Optimising MRF"<<endl;
            flush(cout);
        }
        register int currclass;

        unsigned int numelmasked_currclass_shift[max_numbclass];
        unsigned int image_size_currclass_shift[max_numbclass];
        for(int i=0; i<numclass; i++){
            numelmasked_currclass_shift[i]=i*numelmasked;
            image_size_currclass_shift[i]=i*image_size;
        }
        int curr_short_centreindex;
        for (iz=1; iz<maxiz-1; iz++) {
            for (iy=1; iy<maxiy-1; iy++) {
                indexCentre=(col_size*iy)+(plane_size*iz);
                for (ix=1; ix<maxix-1; ix++) {
                    indexCentre++;
                    Sum_Temp_MRF_Class_Expect = 0;
                    curr_short_centreindex=Long_2_Short_IndicesPtr[indexCentre];
                    if (curr_short_centreindex>=0) {
                        indexWest=Long_2_Short_IndicesPtr[indexCentre-col_size]>-1?Long_2_Short_IndicesPtr[indexCentre-col_size]:0;
                        indexEast=Long_2_Short_IndicesPtr[indexCentre+col_size]>-1?Long_2_Short_IndicesPtr[indexCentre+col_size]:0;
                        indexNorth=Long_2_Short_IndicesPtr[indexCentre-1]>-1?Long_2_Short_IndicesPtr[indexCentre-1]:0;
                        indexSouth=Long_2_Short_IndicesPtr[indexCentre+1]>-1?Long_2_Short_IndicesPtr[indexCentre+1]:0;
                        indexBottom=Long_2_Short_IndicesPtr[indexCentre+plane_size]>-1?Long_2_Short_IndicesPtr[indexCentre+plane_size]:0;
                        indexTop=Long_2_Short_IndicesPtr[indexCentre-plane_size]>-1?Long_2_Short_IndicesPtr[indexCentre-plane_size]:0;
                        for (currclass=0; currclass<numclass; currclass++){
                            Gplane[currclass] = 0.0;
                            Hplane[currclass] = 0.0;
                            Temp_MRF_Class_Expect[currclass] = 0.0;
                            Gplane[currclass]+=Expec[indexWest];
                            Gplane[currclass]+=Expec[indexEast];
                            Gplane[currclass]+=Expec[indexNorth];
                            Gplane[currclass]+=Expec[indexSouth];
                            Hplane[currclass]+=Expec[indexTop];
                            Hplane[currclass]+=Expec[indexBottom];
                            if(currclass<numclass){
                                indexWest+=numelmasked;
                                indexEast+=numelmasked;
                                indexNorth+=numelmasked;
                                indexSouth+=numelmasked;
                                indexTop+=numelmasked;
                                indexBottom+=numelmasked;
                            }
                        }
                        for (currclass=0; currclass<numclass; currclass++){
                            for (neighbourclass=0; neighbourclass<numclass; neighbourclass++){
                                Temp_MRF_Class_Expect[currclass]-=G[currclass+(numclass)*neighbourclass]*Gplane[neighbourclass]+H[currclass+(numclass)*neighbourclass]*Hplane[neighbourclass];
                            }
                            if(MRFbeta==NULL){
                                Temp_MRF_Class_Expect[currclass] = exp(Temp_MRF_Class_Expect[currclass])*AtlasPrior[curr_short_centreindex+numelmasked_currclass_shift[currclass]];
                            }
                            else{
                                Temp_MRF_Class_Expect[currclass] = exp(MRFbeta[curr_short_centreindex]*Temp_MRF_Class_Expect[currclass])*AtlasPrior[curr_short_centreindex+numelmasked_currclass_shift[currclass]];
                            }
                            Sum_Temp_MRF_Class_Expect += Temp_MRF_Class_Expect[currclass];
                        }
                        for (currclass=0; currclass<numclass; currclass++) {
                            MRFpriorPtr[curr_short_centreindex+numelmasked_currclass_shift[currclass]]=(Temp_MRF_Class_Expect[currclass]/Sum_Temp_MRF_Class_Expect);
                        }
                    }
                }
            }
        }
    }
    else{
        for(int currclass=0; currclass<numclass; currclass++){
            for(int i=0; i<(numelmasked);i++){
                MRFprior[i+currclass*numelmasked]=AtlasPrior[i+currclass*numelmasked];
            }
        }

    }

}


void MRFregularization(const PrecisionTYPE * Expec,
                       const PrecisionTYPE * G,
                       const PrecisionTYPE * H,
                       PrecisionTYPE * MRFbeta,
                       PrecisionTYPE * MRFprior,
                       PrecisionTYPE * AtlasPrior,
                       ImageSize * CurrSizes,
                       bool MRFflag,
                       int verbose_level)
{

    int numel=CurrSizes->numel;
    int numclass=CurrSizes->numclass;
    int image_size = (int)(CurrSizes->xsize)*(CurrSizes->ysize)*(CurrSizes->zsize);
    if(MRFflag){
        PrecisionTYPE * MRFpriorPtr = (PrecisionTYPE *)MRFprior;
        int col_size, plane_size, indexCentre, indexWest, indexEast, indexSouth, indexNorth, indexTop, indexBottom;
        int ix, iy, iz,maxiy, maxix, maxiz, neighbourclass;
        PrecisionTYPE Sum_Temp_MRF_Class_Expect;
        col_size = (int)(CurrSizes->xsize);
        plane_size = (int)(CurrSizes->xsize)*(CurrSizes->ysize);
        image_size = (int)(CurrSizes->xsize)*(CurrSizes->ysize)*(CurrSizes->zsize);
        maxix = (int)(CurrSizes->xsize);
        maxiy = (int)(CurrSizes->ysize);
        maxiz = (int)(CurrSizes->zsize);
        PrecisionTYPE Gplane[max_numbclass];
        PrecisionTYPE Hplane[max_numbclass];
        PrecisionTYPE Temp_MRF_Class_Expect[max_numbclass];
        if(verbose_level>0){
            cout << "Optimising MRF"<<endl;
            flush(cout);
        }
        register int currclass;

        unsigned int numel_currclass_shift[max_numbclass];
        unsigned int image_size_currclass_shift[max_numbclass];
        for(int i=0; i<numclass; i++){
            numel_currclass_shift[i]=i*numel;
            image_size_currclass_shift[i]=i*image_size;
        }
        indexCentre=0;
        for (iz=0; iz<maxiz-0; iz++) {
            for (iy=0; iy<maxiy-0; iy++) {
                for (ix=0; ix<maxix-0; ix++) {
                    Sum_Temp_MRF_Class_Expect = 0;
                    indexWest=(indexCentre-col_size)>=0?(indexCentre-col_size):-1;
                    indexEast=(indexCentre+col_size)<numel?(indexCentre+col_size):-1;
                    indexNorth=(indexCentre-1)>=0?(indexCentre-1):-1;
                    indexSouth=(indexCentre+1)<numel?(indexCentre+1):-1;
                    indexBottom=(indexCentre-plane_size)>=0?(indexCentre-plane_size):-1;
                    indexTop=(indexCentre+plane_size)<numel?(indexCentre+plane_size):-1;
                    for (currclass=0; currclass<numclass; currclass++){
                        Gplane[currclass] = 0.0;
                        Hplane[currclass] = 0.0;
                        Temp_MRF_Class_Expect[currclass] = 0.0;
                        Gplane[currclass]+=((indexWest>=0)?Expec[indexWest]:0);
                        Gplane[currclass]+=((indexEast>=0)?Expec[indexEast]:0);
                        Gplane[currclass]+=((indexNorth>=0)?Expec[indexNorth]:0);
                        Gplane[currclass]+=((indexSouth>=0)?Expec[indexSouth]:0);
                        Hplane[currclass]+=((indexTop>=0)?Expec[indexTop]:0);
                        Hplane[currclass]+=((indexBottom>=0)?Expec[indexBottom]:0);
                        if(currclass<numclass){
                            indexWest+=(indexWest>=0)?numel:0;
                            indexEast+=(indexEast>=0)?numel:0;
                            indexNorth+=(indexNorth>=0)?numel:0;
                            indexSouth+=(indexSouth>=0)?numel:0;
                            indexTop+=(indexTop>=0)?numel:0;
                            indexBottom+=(indexBottom>=0)?numel:0;
                        }
                    }
                    for (currclass=0; currclass<numclass; currclass++){
                        for (neighbourclass=0; neighbourclass<numclass; neighbourclass++){
                            Temp_MRF_Class_Expect[currclass]-=G[currclass+(numclass)*neighbourclass]*Gplane[neighbourclass]+H[currclass+(numclass)*neighbourclass]*Hplane[neighbourclass];
                        }
                        Temp_MRF_Class_Expect[currclass] = exp(Temp_MRF_Class_Expect[currclass]) * AtlasPrior[indexCentre+numel_currclass_shift[currclass]];
                        Sum_Temp_MRF_Class_Expect += Temp_MRF_Class_Expect[currclass];
                    }
                    for (currclass=0; currclass<numclass; currclass++) {
                        MRFpriorPtr[indexCentre+numel_currclass_shift[currclass]]=(Temp_MRF_Class_Expect[currclass]/Sum_Temp_MRF_Class_Expect);
                    }
                    indexCentre++;
                }
            }
        }
    }
    else{
        for(int currclass=0; currclass<numclass; currclass++){
            for(int i=0; i<(numel);i++){
                MRFprior[i+currclass*numel]=AtlasPrior[i+currclass*numel];
            }
        }

    }

}




void MRFregularization_mask2D(const PrecisionTYPE * Expec,
                            const PrecisionTYPE * G,
                            const PrecisionTYPE * H,
                            PrecisionTYPE * MRFbeta,
                            PrecisionTYPE * MRFprior,
                            PrecisionTYPE * AtlasPrior,
                            int * Long_2_Short_Indices,
                            int * Short_2_Long_Indices,
                            ImageSize * CurrSizes,
                            bool MRFflag,
                            int verbose_level)
{

    int numelmasked=CurrSizes->numelmasked;
    int numclass=CurrSizes->numclass;
    int image_size = (int)(CurrSizes->xsize)*(CurrSizes->ysize)*(CurrSizes->zsize);
    if(MRFflag){
        PrecisionTYPE * MRFpriorPtr = (PrecisionTYPE *)MRFprior;
        int * Long_2_Short_IndicesPtr = (int *)Long_2_Short_Indices;
        int col_size, indexCentre, indexWest, indexEast, indexSouth, indexNorth;
        int ix, iy,maxiy, maxix, neighbourclass;
        PrecisionTYPE Sum_Temp_MRF_Class_Expect;
        col_size = (int)(CurrSizes->xsize);
        image_size = (int)(CurrSizes->xsize)*(CurrSizes->ysize);
        maxix = (int)(CurrSizes->xsize);
        maxiy = (int)(CurrSizes->ysize);
        PrecisionTYPE Gplane[max_numbclass];
        PrecisionTYPE Temp_MRF_Class_Expect[max_numbclass];
        if(verbose_level>0){
            cout << "Optimising MRF"<<endl;
            flush(cout);
        }
        register int currclass;

        unsigned int numelmasked_currclass_shift[max_numbclass];
        unsigned int image_size_currclass_shift[max_numbclass];
        for(int i=0; i<numclass; i++){
            numelmasked_currclass_shift[i]=i*numelmasked;
            image_size_currclass_shift[i]=i*image_size;
        }
        int curr_short_centreindex;
            for (iy=1; iy<maxiy-1; iy++) {
                indexCentre=(col_size*iy);
                for (ix=1; ix<maxix-1; ix++) {
                    indexCentre++;
                    Sum_Temp_MRF_Class_Expect = 0;
                    curr_short_centreindex=Long_2_Short_IndicesPtr[indexCentre];
                    if (curr_short_centreindex>=0) {
                        indexWest=Long_2_Short_IndicesPtr[indexCentre-col_size]>-1?Long_2_Short_IndicesPtr[indexCentre-col_size]:0;
                        indexEast=Long_2_Short_IndicesPtr[indexCentre+col_size]>-1?Long_2_Short_IndicesPtr[indexCentre+col_size]:0;
                        indexNorth=Long_2_Short_IndicesPtr[indexCentre-1]>-1?Long_2_Short_IndicesPtr[indexCentre-1]:0;
                        indexSouth=Long_2_Short_IndicesPtr[indexCentre+1]>-1?Long_2_Short_IndicesPtr[indexCentre+1]:0;
                        for (currclass=0; currclass<numclass; currclass++){
                            Gplane[currclass] = 0.0;
                            Temp_MRF_Class_Expect[currclass] = 0.0;
                            Gplane[currclass]+=Expec[indexWest];
                            Gplane[currclass]+=Expec[indexEast];
                            Gplane[currclass]+=Expec[indexNorth];
                            Gplane[currclass]+=Expec[indexSouth];
                            if(currclass<numclass){
                                indexWest+=numelmasked;
                                indexEast+=numelmasked;
                                indexNorth+=numelmasked;
                                indexSouth+=numelmasked;
                            }
                        }
                        for (currclass=0; currclass<numclass; currclass++){
                            for (neighbourclass=0; neighbourclass<numclass; neighbourclass++){
                                Temp_MRF_Class_Expect[currclass]-=G[currclass+(numclass)*neighbourclass]*Gplane[neighbourclass];
                            }
                            if(MRFbeta==NULL){
                                Temp_MRF_Class_Expect[currclass] = exp(Temp_MRF_Class_Expect[currclass])*AtlasPrior[curr_short_centreindex+numelmasked_currclass_shift[currclass]];
                            }
                            else{
                                Temp_MRF_Class_Expect[currclass] = exp(MRFbeta[curr_short_centreindex]*Temp_MRF_Class_Expect[currclass])*AtlasPrior[curr_short_centreindex+numelmasked_currclass_shift[currclass]];
                            }
                            Sum_Temp_MRF_Class_Expect += Temp_MRF_Class_Expect[currclass];
                        }
                        for (currclass=0; currclass<numclass; currclass++) {
                            MRFpriorPtr[curr_short_centreindex+numelmasked_currclass_shift[currclass]]=(Temp_MRF_Class_Expect[currclass]/Sum_Temp_MRF_Class_Expect);
                        }
                    }
                }
            }

    }
    else{
        for(int currclass=0; currclass<numclass; currclass++){
            for(int i=0; i<(numelmasked);i++){
                MRFprior[i+currclass*numelmasked]=AtlasPrior[i+currclass*numelmasked];
            }
        }

    }

}


void MRFregularization2D(const PrecisionTYPE * Expec,
                       const PrecisionTYPE * G,
                       const PrecisionTYPE * H,
                       PrecisionTYPE * MRFbeta,
                       PrecisionTYPE * MRFprior,
                       PrecisionTYPE * AtlasPrior,
                       ImageSize * CurrSizes,
                       bool MRFflag,
                       int verbose_level)
{

    int numel=CurrSizes->numel;
    int numclass=CurrSizes->numclass;
    int image_size = (int)(CurrSizes->xsize)*(CurrSizes->ysize)*(CurrSizes->zsize);
    if(MRFflag){
        PrecisionTYPE * MRFpriorPtr = (PrecisionTYPE *)MRFprior;
        int col_size, indexCentre, indexWest, indexEast, indexSouth, indexNorth;
        int ix, iy,maxiy, maxix, neighbourclass;
        PrecisionTYPE Sum_Temp_MRF_Class_Expect;
        col_size = (int)(CurrSizes->xsize);
        image_size = (int)(CurrSizes->xsize)*(CurrSizes->ysize);
        maxix = (int)(CurrSizes->xsize);
        maxiy = (int)(CurrSizes->ysize);
        PrecisionTYPE Gplane[max_numbclass];
        PrecisionTYPE Temp_MRF_Class_Expect[max_numbclass];
        if(verbose_level>0){
            cout << "Optimising MRF"<<endl;
            flush(cout);
        }
        register int currclass;

        unsigned int numel_currclass_shift[max_numbclass];
        unsigned int image_size_currclass_shift[max_numbclass];
        for(int i=0; i<numclass; i++){
            numel_currclass_shift[i]=i*numel;
            image_size_currclass_shift[i]=i*image_size;
        }
        indexCentre=0;
            for (iy=0; iy<maxiy-0; iy++) {
                for (ix=0; ix<maxix-0; ix++) {
                    Sum_Temp_MRF_Class_Expect = 0;
                    indexWest=(indexCentre-col_size)>=0?(indexCentre-col_size):-1;
                    indexEast=(indexCentre+col_size)<numel?(indexCentre+col_size):-1;
                    indexNorth=(indexCentre-1)>=0?(indexCentre-1):-1;
                    indexSouth=(indexCentre+1)<numel?(indexCentre+1):-1;
                    for (currclass=0; currclass<numclass; currclass++){
                        Gplane[currclass] = 0.0;
                        Temp_MRF_Class_Expect[currclass] = 0.0;
                        Gplane[currclass]+=((indexWest>=0)?Expec[indexWest]:0);
                        Gplane[currclass]+=((indexEast>=0)?Expec[indexEast]:0);
                        Gplane[currclass]+=((indexNorth>=0)?Expec[indexNorth]:0);
                        Gplane[currclass]+=((indexSouth>=0)?Expec[indexSouth]:0);
                        if(currclass<numclass){
                            indexWest+=(indexWest>=0)?numel:0;
                            indexEast+=(indexEast>=0)?numel:0;
                            indexNorth+=(indexNorth>=0)?numel:0;
                            indexSouth+=(indexSouth>=0)?numel:0;
                        }
                    }
                    for (currclass=0; currclass<numclass; currclass++){
                        for (neighbourclass=0; neighbourclass<numclass; neighbourclass++){
                            Temp_MRF_Class_Expect[currclass]-=G[currclass+(numclass)*neighbourclass]*Gplane[neighbourclass];
                        }
                        Temp_MRF_Class_Expect[currclass] = exp(Temp_MRF_Class_Expect[currclass]) * AtlasPrior[indexCentre+numel_currclass_shift[currclass]];
                        Sum_Temp_MRF_Class_Expect += Temp_MRF_Class_Expect[currclass];
                    }
                    for (currclass=0; currclass<numclass; currclass++) {
                        MRFpriorPtr[indexCentre+numel_currclass_shift[currclass]]=(Temp_MRF_Class_Expect[currclass]/Sum_Temp_MRF_Class_Expect);
                    }
                    indexCentre++;
                }
            }
    }
    else{
        for(int currclass=0; currclass<numclass; currclass++){
            for(int i=0; i<(numel);i++){
                MRFprior[i+currclass*numel]=AtlasPrior[i+currclass*numel];
            }
        }

    }

}
