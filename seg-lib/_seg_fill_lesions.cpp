#ifndef _SEG_FILL_LESIONS_CPP
#define _SEG_FILL_LESIONS_CPP

#ifdef _OPENMP
#include "omp.h"
#endif

#include "_seg_fill_lesions.h"

template <class T>
seg_fill_lesions<T>::seg_fill_lesions() {
    this->inputImage=NULL;
    this->meanImage=NULL;
    this->normImage=NULL;
    this->imgPtr=NULL;
    this->lesMaskPtr=NULL;
    this->maskPtr=NULL;
    this->normImgPtr=NULL;
    this->meanImgPtr=NULL;
    this->currSize=NULL;
    this->verbose=0;
    this->debug=0;
    this->mult=0.1;
    this->patchSearchAreaSize=4;
    this->percentage=0.5;
    this->k=2.0;
    this->patch2D=false;
}

template <class T>
seg_fill_lesions<T>::seg_fill_lesions(const seg_fill_lesions &copy) {
    this->inputImage=copy->inputImage;
    this->meanImage=copy->meanImage;
    this->normImage=copy->normImage;
    this->outputImage=copy->outputImage;
    this->imgPtr=copy->imgPtr;
    this->lesMaskPtr=copy->lesMaskPtr;
    this->maskPtr=copy->maskPtr;
    this->normImgPtr=copy->normImgPtr;
    this->meanImgPtr=copy->meanImgPtr;
    this->tmpLesMask=copy->tmpLesMask;
    this->tmpImg=copy->tmpImg;
    this->originalLesMask=copy->originalLesMask;
    this->tmpNorm=copy->tmpNorm;
    this->currSize=copy->currSize;
    this->verbose=copy->verbose;
    this->debug=copy->debug;
    this->mult=copy->mult;
    this->patchSearchAreaSize=copy->patchSearchAreaSize;
    this->percentage=copy->percentage;
    this->numTP=copy->numTP;
    this->k=copy->k;
    this->patch2D=copy->patch2D;
}

template <class T>
int seg_fill_lesions<T>::operator=(const seg_fill_lesions &copy) {
    this->inputImage=copy->inputImage;
    this->meanImage=copy->meanImage;
    this->normImage=copy->normImage;
    this->outputImage=copy->outputImage;
    this->imgPtr=copy->imgPtr;
    this->lesMaskPtr=copy->lesMaskPtr;
    this->maskPtr=copy->maskPtr;
    this->normImgPtr=copy->normImgPtr;
    this->meanImgPtr=copy->meanImgPtr;
    this->tmpLesMask=copy->tmpLesMask;
    this->tmpImg=copy->tmpImg;
    this->originalLesMask=copy->originalLesMask;
    this->tmpNorm=copy->tmpNorm;
    this->currSize=copy->currSize;
    this->verbose=copy->verbose;
    this->debug=copy->debug;
    this->mult=copy->mult;
    this->patchSearchAreaSize=copy->patchSearchAreaSize;
    this->percentage=copy->percentage;
    this->numTP=copy->numTP;
    this->k=copy->k;
    this->patch2D=copy->patch2D;
    
    return *this;
}


template<class T>
void seg_fill_lesions<T>::setInputImage(nifti_image *input) {
    this->inputImage = input;
    this->init();
}

template<class T>
void seg_fill_lesions<T>::init() {
    this->imgPtr = static_cast<T *>(this->inputImage->data);
    this->setNumTP(this->inputImage->nt);
    this->currSize = new ImageSize [1]();
    this->currSize->numel=(long)(this->inputImage->nx*this->inputImage->ny*this->inputImage->nz);
    this->currSize->xsize=this->inputImage->nx;
    this->currSize->ysize=this->inputImage->ny;
    this->currSize->zsize=this->inputImage->nz;
    this->currSize->usize=this->inputImage->nu;
    this->currSize->tsize=this->inputImage->nt;

    // We initialize mean image
    this->meanImage = nifti_copy_nim_info(this->inputImage);
    this->meanImage->data = (void *) calloc(this->meanImage->nvox, sizeof(T));
    this->meanImgPtr = static_cast<T *>(this->meanImage->data);

    // We initialize norm image
    this->normImage = nifti_copy_nim_info(this->inputImage);
    this->normImage->data = (void *) calloc(this->normImage->nvox, sizeof(T));
    this->normImgPtr = static_cast<T *>(this->normImage->data);
}

template<class T>
void seg_fill_lesions<T>::setInputMask(nifti_image *mask) {
    this->inputMask = mask;
    this->maskPtr = static_cast<T *>(this->inputMask->data);
    cout <<"Mixing input mask and lesion mask" << endl;
    long tp;
    #ifdef _OPENMP
    #pragma omp parallel for default(none)\
        shared(std::cout)\
        private(tp)
    #endif
    for(tp=0;tp<this->getNumTP();tp++) {
        if(this->getVerbose()) cout <<"TP "<<tp<<endl;
        for(long i=0; i<this->getSingleVolumSize(); i++) {
             this->maskPtr[i+tp*this->getSingleVolumSize()] = this->maskPtr[i+tp*this->getSingleVolumSize()]>0 || this->lesMaskPtr[i+tp*this->getSingleVolumSize()]>0;
        }
    }
    char filename[100];
    if(this->getDebug()) {
        sprintf(filename,"segFillLesions_mixed_mask.nii.gz");
        this->saveImagePtr(this->maskPtr,this->inputImage,filename);
    }
}

template<class T>
void seg_fill_lesions<T>::setInputLesionMask(nifti_image *mask) {
    this->inputLesionMask = mask;
    this->lesMaskPtr = static_cast<T *>(this->inputLesionMask->data);
}

template<class T>
void seg_fill_lesions<T>::setPatchSearchAreaSize(int value) {
    this->patchSearchAreaSize = value;
}

template<class T>
int seg_fill_lesions<T>::getPatchSearchAreaSize() {
    return this->patchSearchAreaSize;
}

template<class T>
void seg_fill_lesions<T>::setPatchPercentage(float value) {
    this->percentage = value;
}

template<class T>
void seg_fill_lesions<T>::setSmoothing(float value) {
    this->mult = value;
}

template<class T>
float seg_fill_lesions<T>::getSmoothing() {
    return this->mult;
}

template<class T>
void seg_fill_lesions<T>::setExpandingPercentage(float value) {
    this->expanding = value;
}

template<class T>
float seg_fill_lesions<T>::getExpandingPercentage() {
    return this->expanding;
}

template<class T>
void seg_fill_lesions<T>::setVerbose(int value) {
    this->verbose = value;
}

template<class T>
int seg_fill_lesions<T>::getVerbose() {
    return this->verbose;
}

template<class T>
void seg_fill_lesions<T>::setDebug(int value) {
    this->debug = value;
}

template<class T>
int seg_fill_lesions<T>::getDebug() {
    return this->debug;
}

template<class T>
void seg_fill_lesions<T>::setNumTP(int value) {
    this->numTP = value;
}

template<class T>
int seg_fill_lesions<T>::getNumTP() {
    return this->numTP;
}

template<class T>
void seg_fill_lesions<T>::setK(float value) {
    this->k = value;
}

template<class T>
float seg_fill_lesions<T>::getK() {
    return this->k;
}

template<class T>
long seg_fill_lesions<T>::getSingleVolumSize() {
    return this->currSize->numel;
}

template<class T>
long seg_fill_lesions<T>::getTotalVolumSize() {
    return this->currSize->numel*this->getNumTP();
}

template<class T>
long seg_fill_lesions<T>::getSliceSize() {
    return this->currSize->xsize*this->currSize->ysize;
}

template<class T>
long seg_fill_lesions<T>::getXAxisSize() {
    return this->currSize->xsize;
}

template<class T>
long seg_fill_lesions<T>::getYAxisSize() {
    return this->currSize->ysize;
}

template<class T>
long seg_fill_lesions<T>::getZAxisSize() {
    return this->currSize->zsize;
}

template<class T>
void seg_fill_lesions<T>::setDimensionality(bool is2D) {
    this->patch2D=is2D;
}

template<class T>
bool seg_fill_lesions<T>::isIt2D() {
    return this->patch2D;
}

template<class T>
float seg_fill_lesions<T>::calculateDistance(int location1, int location2,int patchSize,int countlesionvoxel){

    const int numvox=this->getSingleVolumSize();
    float totalpatchsize=0;
    int patchSizeZ=0;
    if(this->isIt2D()) {
        totalpatchsize=(patchSize*2+1)*(patchSize*2+1)*this->getNumTP();
	patchSizeZ=0;
    }
    else {
        totalpatchsize=(patchSize*2+1)*(patchSize*2+1)*(patchSize*2+1)*this->getNumTP();
	patchSizeZ=patchSize;
    }
        
    float distance=0;
    int count=0;
    int shiftx=0;
    int shifty=0;
    int shiftz=0;

    int maxshift=patchSize+this->getXAxisSize()*(patchSize+(this->getYAxisSize()*patchSizeZ));

    bool WillTouchBorder=1;
    if( (location1+maxshift)<numvox &&
            (location1-maxshift)>=0 &&
            (location2+maxshift)<numvox &&
            (location2-maxshift)>=0){
        WillTouchBorder=0;
    }

    for(shiftx=-patchSize; shiftx<=patchSize; shiftx++){
        for(shifty=-patchSize; shifty<=patchSize; shifty++){
            for(shiftz=-patchSizeZ; shiftz<=patchSizeZ; shiftz++){
                int shift=(shiftx)+this->getXAxisSize()*(shifty+(this->getYAxisSize()*shiftz));
                int index1=location1+shift;
                int index2=location2+shift;
                for(int tp=0;tp<this->getNumTP();tp++) {
                    if(WillTouchBorder)
                    {
                        if(index1<numvox   &&  // Is within the index bounds (also checks z)
                                index1>=0  &&  // Is within the index bounds (also checks z)

                                index2<numvox   &&  // Is within the index bounds (also checks z)
                                index2>=0       &&  // Is within the index bounds (also checks z)

                                this->tmpLesMask[index1+tp*numvox]==0 &&  // Isn't within the mask
                                this->tmpLesMask[index2+tp*numvox]==0 )   // Isn't within the mask
                        {
                            distance+=(this->normImgPtr[index1+tp*numvox]-this->normImgPtr[index2+tp*numvox])*(this->normImgPtr[index1+tp*numvox]-this->normImgPtr[index2+tp*numvox]);
                            count++;
                        }
                    }
                    else{
                        if(this->tmpLesMask[index1+tp*numvox]==0     &&  // Isn't within the mask
                                this->tmpLesMask[index2+tp*numvox]==0)   // Isn't within the mask
                        {
                            distance+=(this->normImgPtr[index1+tp*numvox]-this->normImgPtr[index2+tp*numvox])*(this->normImgPtr[index1+tp*numvox]-this->normImgPtr[index2+tp*numvox]);
                            count++;
                        }
                    }

                }
            }
        }
    }
    
    return (count< round(this->percentage*(float)(totalpatchsize-countlesionvoxel))) ? std::numeric_limits<float>::quiet_NaN() : (distance/(float)(pow(count,this->getK())));
}

template<class T>
void seg_fill_lesions<T>::normalizeImageIntesities(float newMin,float newMax) {
    long tp;
    #ifdef _OPENMP
    #pragma omp parallel for default(none)\
        shared(newMin,newMax,std::cout)\
        private(tp)
    #endif
    for(tp=0; tp<this->getNumTP(); tp++) {
        float max=std::numeric_limits<float>::min();
        float min=std::numeric_limits<float>::max();
        for(long i=0; i<this->getSingleVolumSize(); i++) {
            max=max>this->imgPtr[i+tp*this->getSingleVolumSize()]?max:this->imgPtr[i+tp*this->getSingleVolumSize()];
            min=min<this->imgPtr[i+tp*this->getSingleVolumSize()]?min:this->imgPtr[i+tp*this->getSingleVolumSize()];
        }
        if(this->getVerbose()) cout<<"["<<tp<<"] MIN="<<min<<" MAX="<<max<<endl;
        for(long i=0; i<this->getSingleVolumSize(); i++) {
            this->normImgPtr[i+tp*this->getSingleVolumSize()]=newMin+(this->imgPtr[i+tp*this->getSingleVolumSize()]-min)*(newMax-newMin)/(max-min);
        }
    }
}

template<class T>
void seg_fill_lesions<T>::calculateEuclideanDistance(float *Distance) {
    long tp;
    ImageSize *eucl_Size = new ImageSize [1]();
    eucl_Size->numel=(long)(this->inputImage->nx*this->inputImage->ny*this->inputImage->nz);
    eucl_Size->xsize=this->inputImage->nx;
    eucl_Size->ysize=this->inputImage->ny;
    eucl_Size->zsize=this->inputImage->nz;
    eucl_Size->usize=this->inputImage->nu;
    eucl_Size->tsize=1;
    #ifdef _OPENMP
    #pragma omp parallel for \
        shared(Distance,std::cout)\
        private(tp)
    #endif
    for(tp=0;tp<this->getNumTP();tp++) {
        if(this->getVerbose()) cout <<"TP "<<tp<<endl;
        bool *Lable= new bool[this->getSingleVolumSize()];
        float *Speed= new float[this->getSingleVolumSize()];
        float *LableFloat= new float[this->getSingleVolumSize()];
        for(long i=0; i<this->getSingleVolumSize(); i++) {
            if(this->maskPtr!=NULL) {
                Lable[i]=this->maskPtr[i+tp*this->getSingleVolumSize()]<1 || this->lesMaskPtr[i+tp*this->getSingleVolumSize()]>0;
            }
            else {
                Lable[i]=this->lesMaskPtr[i+tp*this->getSingleVolumSize()]>0;
            }
            Speed[i]=1.0f;
            LableFloat[i]=Lable[i]?1.0:0.0;
        }
        this->currSize->tsize=1;
        float * result = DoubleEuclideanDistance_3D(Lable,Speed,eucl_Size);
        for(long ii=0;ii<this->getSingleVolumSize();ii++) {
            int tpf=this->isMaskInAnyTP(ii,this->originalLesMask);
            if(tpf>=0){
                Distance[ii+tp*this->getSingleVolumSize()]=result[ii];
            }
        }
    }
}

template<class T>
void seg_fill_lesions<T>::expandPatches(float *Distance) {
    long tp;
    long change=0,total=0;
    #ifdef _OPENMP
    #pragma omp parallel for default(none)\
        shared(Distance,std::cout,change,total)\
        private(tp)
    #endif
    for(tp=0;tp<this->getNumTP();tp++) {
        if(this->getVerbose()) cout <<"TP "<<tp<<endl;
        for(long i=0; i<this->getSingleVolumSize(); i++) {
            long index=i+tp*this->getSingleVolumSize();
            if(this->tmpLesMask[index]>0) {
                int patchSize=round(Distance[index]+1);
                int oldPatchSize=patchSize;
                int patchshiftx=(patchSize*2+1);
                float totalpatchsize=patchshiftx*patchshiftx*patchshiftx*this->getNumTP();
                int countlesionvoxel=this->countLesionVoxels(index,patchSize);
                while (((totalpatchsize-countlesionvoxel)/totalpatchsize)<this->getExpandingPercentage()) {
                    patchSize++;
                    patchshiftx=(patchSize*2+1);
                    totalpatchsize=patchshiftx*patchshiftx*patchshiftx*this->getNumTP();
                    countlesionvoxel=this->countLesionVoxels(index,patchSize);
                }
                if(oldPatchSize!=patchSize) {
                    Distance[index]=patchSize-1;
                    change++;//cout<<oldPatchSize<<" "<<patchSize<<endl;
                }
                total++;
            }
        }
    }
    if(this->getVerbose()) {
        cout<<"Patch size changed in :"<<change<<"/"<<total<<endl;
    }
}

template<class T>
long seg_fill_lesions<T>::isMaskInAnyTP(long index,int *inputMask) {
    long tp=0;
    bool mask=false;
    while(tp<this->getNumTP() && !mask) {
        mask=inputMask[index+tp*this->getSingleVolumSize()]>0;
        tp++;
    }
    if(mask) tp=tp-1;
    else tp=-1;
    return tp;
}

template<class T>
bool seg_fill_lesions<T>::isValidVoxelInAnyTP(long index) {
    bool valid=false;
     if(this->maskPtr!=NULL) {
        long tp=0;
        while(tp<this->getNumTP() && !valid) {
            valid=this->maskPtr[index+tp*this->getSingleVolumSize()]>0;
            tp++;
        }
    }
    else {
        valid=true;
    }
    return valid;
}


template<class T>
void seg_fill_lesions<T>::runIt(){
    long index=0;
    char filename[100];
    long i,tp;
    long countvox=0;
    long curcountvox=0;

    if(this->getVerbose()) {
        cout << "File dimensions"<<endl;
        cout<<"Size ["<<this->getXAxisSize()<<","<<this->getYAxisSize()<<","<<this->getZAxisSize()<<"] TP="<<this->getNumTP()<<endl;
        cout << "Normalizing Image"<<endl;
    }
    // Intensity image normalization
    this->normalizeImageIntesities(0.0f,1024.0f);

    if(this->getDebug()) {
        sprintf(filename,"segFillLesions_normalized_image.nii.gz");
        this->saveImage(this->normImage,filename);
    }

    this->tmpLesMask= new int [this->getTotalVolumSize()];
    this->originalLesMask= new int [this->getTotalVolumSize()];
    this->tmpImg= new float [this->getTotalVolumSize()];
    this->tmpNorm= new float [this->getTotalVolumSize()];
    int *level= new int [this->getTotalVolumSize()];
    // For all voxels and timepoints
    for(i=0; i<this->getTotalVolumSize(); i++) {
        level[i]=0;
        this->tmpLesMask[i]=this->lesMaskPtr[i]>0;
        this->originalLesMask[i]=this->lesMaskPtr[i]>0;
        this->meanImgPtr[i]=this->lesMaskPtr[i]>0?std::numeric_limits<float>::quiet_NaN():this->normImgPtr[i];
    }
    // Total of voxels to be analyzed
    for(i=0; i<this->getSingleVolumSize(); i++) {
        int tpf=this->isMaskInAnyTP(i,this->tmpLesMask);
        if(tpf>=0) countvox++;
    }

    if(this->getVerbose()) cout<<"Calculating Euclidean distance"<<endl;
    float *Distance = new float[this->getTotalVolumSize()];
    this->calculateEuclideanDistance(Distance);
    if(this->getDebug()) {
        sprintf(filename,"segFillLesions_euclidean_distance.nii.gz");
        this->saveImagePtr(Distance,this->inputImage,filename);
    }
    float max=0;
    // For all voxels and timepoints
    for(tp=0;tp<this->getNumTP();tp++) {
        for(i=0; i<this->getSingleVolumSize(); i++) {
            int tpf=this->isMaskInAnyTP(i,this->tmpLesMask);
            if(tpf>=0) max=max>Distance[i+tp*this->getSingleVolumSize()]?max:Distance[i+tp*this->getSingleVolumSize()];
        }
    }
    max=round(max+1);
    if(this->getVerbose()) {
        cout<<"Search area size="<<max<<endl;
        cout<<"Calculating Smooth Mean Image"<<endl;
    }
    BlockSmoothing(this->meanImage,NULL,max*2+1);
    if(this->getDebug()) {
        sprintf(filename,"segFillLesions_mean_image-%d.nii.gz",0);
        this->saveImage(this->meanImage,filename);
    }

    if(this->getExpandingPercentage()>0) {
        if(this->getVerbose()) cout<<"Expanding euclidean maps to improve patches"<<endl;
        this->expandPatches(Distance);
        if(this->getDebug()) {
            sprintf(filename,"segFillLesions_expanded_distance.nii.gz");
            this->saveImagePtr(Distance,this->inputImage,filename);
        }
    }

    if(this->getVerbose()) cout<<"Getting max and min for each TP"<<endl;
    float *maxinten=new float[this->getNumTP()];
    float *mininten=new float[this->getNumTP()];
    for(tp=0;tp<this->getNumTP();tp++) {
        mininten[tp]=std::numeric_limits<float>::max();
        maxinten[tp]=std::numeric_limits<float>::min();
        for(i=0; i<this->getSingleVolumSize(); i++) {
            maxinten[tp]=(this->meanImgPtr[i+tp*this->getSingleVolumSize()]>maxinten[tp])?this->meanImgPtr[i+tp*this->getSingleVolumSize()]:maxinten[tp];
            mininten[tp]=(this->meanImgPtr[i+tp*this->getSingleVolumSize()]<mininten[tp])?this->meanImgPtr[i+tp*this->getSingleVolumSize()]:mininten[tp];
        }
        if(this->getVerbose()) cout<<"["<<tp<<"] MIN="<<mininten[tp]<<" MAX="<<maxinten[tp]<<endl;
    }
    if(this->getVerbose()) cout<<"Filling lesions"<<endl;
    long count=1;
    int iteration=0;
    int lastpercent=0;
    while(count>0 && iteration<max+1){
        iteration++;
        count=0;
        if(this->getVerbose()) {
            cout<<"Iteration: "<<iteration<<"/"<<(int)(max+1)<<endl;
            cout<<"Filled voxels: "<<(float)(curcountvox)<<"/"<<(float)(countvox)<<endl;
        }
        int inz=0;
        #ifdef _OPENMP
        #pragma omp parallel for \
            shared(inz,curcountvox,std::cout,Distance,maxinten,mininten,level,iteration,lastpercent,countvox)\
            private(index)\
            reduction(+:count)
        #endif
        for(inz=0; inz<this->getZAxisSize(); inz++) {

            if(lastpercent!=floor((float)(curcountvox)/(float)(countvox)*100.0f)){
                if(this->getVerbose()) cout<<"Filled voxels: "<<(float)(curcountvox)<<"/"<<(float)(countvox)<<endl;
                lastpercent=floor((float)(curcountvox)/(float)(countvox)*100.0f);
            }
            float mindistance=0;
            float curdist=0;
            int inx=0;
            int iny=0;
            int shiftx=0;
            int shifty=0;
            int shiftz=0;
            long index2=0;
            for(iny=0; iny<this->getYAxisSize(); iny++) {
                for(inx=0; inx<this->getXAxisSize(); inx++) {

                    index=inx+this->getXAxisSize()*iny+this->getSliceSize()*inz;
                    long tpf=this->isMaskInAnyTP(index,this->tmpLesMask); // We obtain the timepoint
                    // If it is mask in any timepoint
                    if( tpf>=0 &&
                            (long)(index+1)<this->getSingleVolumSize() &&
                            (long)(index-1)>=0 &&
                            (long)(index+this->getXAxisSize())<this->getSingleVolumSize() &&
                            (long)(index-this->getXAxisSize())>=0 &&
                            ((long)(index+this->getSliceSize())<this->getSingleVolumSize() || this->isIt2D())  &&
                            ((long)(index-this->getSliceSize())>=0 || this->isIt2D()))

                    {   // Boundary for the time point affected
                        if(level[index]<1) level[index]=1;
                        if((this->tmpLesMask[(index+tpf*this->getSingleVolumSize())+1]==0 ||
                            this->tmpLesMask[(index+tpf*this->getSingleVolumSize())-1]==0 ||
                            this->tmpLesMask[(index+tpf*this->getSingleVolumSize())+this->getXAxisSize()]==0 ||
                            this->tmpLesMask[(index+tpf*this->getSingleVolumSize())-this->getXAxisSize()]==0 ||
                            (this->tmpLesMask[(index+tpf*this->getSingleVolumSize())+this->getSliceSize()]==0 || this->isIt2D()) ||
                            (this->tmpLesMask[(index+tpf*this->getSingleVolumSize())-this->getSliceSize()]==0 || this->isIt2D()) )
                                ){
                            if(level[index]<3) level[index]=3;
                            curcountvox++;
                            count++;
                            // If it is a lesion, then search
                            mindistance=std::numeric_limits<float>::max();
                            int maxPatchsize=0;
                            // For all the timepoints we search the biggest distance to the boundary
                            for(long ii=0;ii<this->getNumTP();ii++) {
                                int ps=round(Distance[index+ii*this->getSingleVolumSize()]+1);
                                if(maxPatchsize<ps) {
                                    maxPatchsize=ps;
                                }
                            }
                            int shiftrealsize=maxPatchsize*this->getPatchSearchAreaSize(); //>15?maxPatchsize*3:15;
                            int totalLesionVoxels=this->countLesionVoxels(index,maxPatchsize);
                            for(shiftz=-shiftrealsize; shiftz<=shiftrealsize; shiftz++) {
                                for(shifty=-shiftrealsize; shifty<=shiftrealsize; shifty++) {
                                    for(shiftx=-shiftrealsize; shiftx<=shiftrealsize; shiftx++) {
                                        index2=(inx+shiftx)+this->getXAxisSize()*(iny+shifty)+this->getSliceSize()*(inz+shiftz);
                                        if(index2<this->getSingleVolumSize() && // Is within the index bounds (also checks z)
                                                index2>=0 && // Is within the index bounds (also checks z)
                                                this->isMaskInAnyTP(index2,this->originalLesMask)<0 &&//this->originalMask[index2+tpf*this->getSingleVolumSize()]==0 && // Is outside the mask
                                                (inx+shiftx) < this->getXAxisSize() && // Does it flip around the image
                                                (inx-shiftx) >= 0 && // Does it flip around the image
                                                (iny+shifty) < this->getYAxisSize() && // Does it flip around the image
                                                (iny-shifty) >= 0 &&     // Does it flip around the image
                                                index!=index2 &&
                                                this->tmpLesMask[index2+tpf*this->getSingleVolumSize()]==0 && // Is outside the mask for this timepoint
                                                this->isValidVoxelInAnyTP(index2) &&
                                                fabs(0.1f*(maxinten[tpf]-mininten[tpf]))>fabs(this->meanImgPtr[index+tpf*this->getSingleVolumSize()]-this->meanImgPtr[index2+tpf*this->getSingleVolumSize()]) // For the specific timepoint
                                                )
                                        {
                                            // Multidimensional patch (4D)					    
                                            curdist=this->calculateDistance(index,index2,maxPatchsize,totalLesionVoxels);
                                            if(mindistance > curdist && isnan(curdist)==0) {
                                                mindistance=curdist;
                                                for(int ii=0;ii<this->getNumTP();ii++) {
                                                    if(this->tmpLesMask[index+ii*this->getSingleVolumSize()]==1) { // We replace the values only for the timepoints that are masked.
                                                        this->imgPtr[index+ii*this->getSingleVolumSize()]=this->imgPtr[index2+ii*this->getSingleVolumSize()];
                                                        this->normImgPtr[index+ii*this->getSingleVolumSize()]=this->normImgPtr[index2+ii*this->getSingleVolumSize()];
                                                        this->lesMaskPtr[index+ii*this->getSingleVolumSize()]=2;
                                                        if(level[index]<6) level[index]=6;
                                                    }
                                                }
                                            }
                                            else {
                                                if(level[index]<5) level[index]=5;
                                            }
                                        }
                                        else {
                                            if(level[index]<4) level[index]=4;
                                        }
                                    }
                                }

                            }
                        }
                    }
                    else{
                        if(level[index]<2) level[index]=2;
                    }
                }
            }
        }
        if(this->getDebug()) {
            sprintf(filename,"segFillLesions_filemask-med-%d.nii.gz",iteration);
            this->saveImagePtr(level,this->inputLesionMask,filename);
        }
        for(tp=0;tp<this->getNumTP();tp++) {
            for(i=0; i<this->getSingleVolumSize(); i++) {
                level[i+tp*this->getSingleVolumSize()]=0;
            }
        }
        #ifdef _OPENMP
        #pragma omp parallel for \
            private(tp,i)
        #endif
        for(tp=0;tp<this->getNumTP();tp++) {
            for(i=0; i<this->getSingleVolumSize(); i++) {
                this->tmpLesMask[i+tp*this->getSingleVolumSize()]=this->lesMaskPtr[i+tp*this->getSingleVolumSize()];
                this->tmpImg[i+tp*this->getSingleVolumSize()]=this->imgPtr[i+tp*this->getSingleVolumSize()];
                this->tmpNorm[i+tp*this->getSingleVolumSize()]=this->normImgPtr[i+tp*this->getSingleVolumSize()];
            }
        }
        #ifdef _OPENMP
        #pragma omp parallel for \
            private(tp,i)
        #endif
        for(tp=0;tp<this->getNumTP();tp++) {
            for(i=0; i<this->getSingleVolumSize(); i++) {
                if(this->tmpLesMask[i+tp*this->getSingleVolumSize()]==2) {
                    float intensity=this->tmpImg[i+tp*this->getSingleVolumSize()];
                    float normintensity=this->tmpNorm[i+tp*this->getSingleVolumSize()];
                    float mult=this->getSmoothing();
                    float density=1;

                    if((i-1)>=0 && this->tmpLesMask[(i+tp*this->getSingleVolumSize())-1]!=1 ){
                        intensity+=mult*this->tmpImg[(i+tp*this->getSingleVolumSize())-1];
                        normintensity+=mult*this->tmpNorm[(i+tp*this->getSingleVolumSize())-1];
                        density+=mult;
                    }
                    if((i+1)<this->getSingleVolumSize() && this->tmpLesMask[(i+tp*this->getSingleVolumSize())+1]!=1 ){
                        intensity+=mult*this->tmpImg[(i+tp*this->getSingleVolumSize())+1];
                        normintensity+=mult*this->tmpNorm[(i+tp*this->getSingleVolumSize())+1];
                        density+=mult;
                    }

                    if((i-this->getXAxisSize())>=0 && this->tmpLesMask[(i+tp*this->getSingleVolumSize())-this->getXAxisSize()]!=1 ){
                        intensity+=mult*this->tmpImg[(i+tp*this->getSingleVolumSize())-this->getXAxisSize()];
                        normintensity+=mult*this->tmpNorm[(i+tp*this->getSingleVolumSize())-this->getXAxisSize()];
                        density+=mult;
                    }
                    if((i+this->getXAxisSize())<this->getSingleVolumSize() && this->tmpLesMask[(i+tp*this->getSingleVolumSize())+this->getXAxisSize()]!=1 ){
                        intensity+=mult*this->tmpImg[(i+tp*this->getSingleVolumSize())+this->getXAxisSize()];
                        normintensity+=mult*this->tmpNorm[(i+tp*this->getSingleVolumSize())+this->getXAxisSize()];
                        density+=mult;
                    }

                    if(!this->isIt2D() && (i-this->getSliceSize())>=0 && this->tmpLesMask[(i+tp*this->getSingleVolumSize())-this->getSliceSize()]!=1 ){
                        intensity+=mult*this->tmpImg[(i+tp*this->getSingleVolumSize())-this->getSliceSize()];
                        normintensity+=mult*this->tmpNorm[(i+tp*this->getSingleVolumSize())-this->getSliceSize()];
                        density+=mult;
                    }
                    if(!this->isIt2D() && (i+this->getSliceSize())<this->getSingleVolumSize() && this->tmpLesMask[(i+tp*this->getSingleVolumSize())+this->getSliceSize()]!=1 ){
                        intensity+=mult*this->tmpImg[(i+tp*this->getSingleVolumSize())+this->getSliceSize()];
                        normintensity+=mult*this->tmpNorm[(i+tp*this->getSingleVolumSize())+this->getSliceSize()];
                        density+=mult;
                    }

                    this->imgPtr[i+tp*this->getSingleVolumSize()]=intensity/density;
                    this->normImgPtr[i+tp*this->getSingleVolumSize()]=normintensity/density;
                    this->lesMaskPtr[i+tp*this->getSingleVolumSize()]=0;
                }
            }
        }
        if(this->getDebug()) {
            sprintf(filename,"segFillLesions_file-%d.nii.gz",iteration);
            this->saveImage(this->inputImage,filename);
            sprintf(filename,"segFillLesions_filemask-%d.nii.gz",iteration);
            this->saveImage(this->inputLesionMask,filename);
        }
        for(i=0;i<this->getTotalVolumSize();i++) {
            this->tmpLesMask[i]=this->lesMaskPtr[i]>0;
            this->meanImgPtr[i]=this->lesMaskPtr[i]>0?std::numeric_limits<float>::quiet_NaN():this->normImgPtr[i];
        }
        BlockSmoothing(this->meanImage,NULL,max*2+1);
        if(this->getDebug()) {
            sprintf(filename,"segFillLesions_mean_image-%d.nii.gz",iteration);
            this->saveImage(this->meanImage,filename);
        }
    }
    if(this->getVerbose()) {
        std::cout<<"Filled voxels: "<<(float)(curcountvox)<<"/"<<(float)(countvox)<<endl;
        std::cout << "Done"<<endl;
    }

    nifti_image_free(this->normImage);
    nifti_image_free(this->meanImage);
    delete [] this->tmpLesMask;
    delete [] this->tmpImg;
    delete [] this->tmpNorm;
    delete [] Distance;

    return;
}

template<class T>
void seg_fill_lesions<T>::saveImagePtr(int * ImagePtr,nifti_image *Image,char *filename_out) {
    float *tmp= new float [Image->nvox*Image->nt];
    for(size_t i=0; i<Image->nvox*(size_t)Image->nt; i++) {
        tmp[i]=(float)ImagePtr[i];
    }
    saveImagePtr(tmp,Image,filename_out);
    delete [] tmp;
}

template<class T>
void seg_fill_lesions<T>::saveImagePtr(float * ImagePtr,nifti_image *Image,char *filename_out) {
    std::cout<<"Saving intermidium file: "<<filename_out<<std::endl;
    // saving output
    nifti_image * OutputImage = nifti_copy_nim_info(Image);
    OutputImage->datatype=NIFTI_TYPE_FLOAT32;
    nifti_set_filenames(OutputImage,filename_out,0,0);
    OutputImage->dim[1]=Image->nx;
    OutputImage->dim[2]=Image->ny;
    OutputImage->dim[3]=Image->nz;
    OutputImage->dim[4]=OutputImage->nt=Image->nt;
    OutputImage->dim[5]=OutputImage->nu=Image->nu;
    OutputImage->dim[6]=OutputImage->nv=1;
    OutputImage->dim[7]=OutputImage->nw=1;
    OutputImage->dim[0]=3;
    OutputImage->dim[0]=(OutputImage->dim[4]>1?4:OutputImage->dim[0]);
    OutputImage->dim[0]=(OutputImage->dim[5]>1?5:OutputImage->dim[0]);
    OutputImage->dim[0]=(OutputImage->dim[6]>1?6:OutputImage->dim[0]);
    OutputImage->dim[0]=(OutputImage->dim[7]>1?7:OutputImage->dim[0]);
    OutputImage->data = (void *) calloc(Image->nx*Image->ny*Image->nz*Image->nt*Image->nu, sizeof(float));
    float * OutputImagePtr = static_cast<float *>(OutputImage->data);
    for(long i=0; i<(long)(Image->nx*Image->ny*Image->nz*Image->nt*Image->nu); i++)
    {
        OutputImagePtr[i]=ImagePtr[i];
    }
    nifti_image_write(OutputImage);
    nifti_image_free(OutputImage);
}

template<class T>
void seg_fill_lesions<T>::saveImage(nifti_image * Image,char *filename_out) {
    std::cout<<"Saving intermidium file: "<<filename_out<<std::endl;
    // saving output
    nifti_set_filenames(Image,filename_out,0,0);
    nifti_image_write(Image);
}

template<class T>
int seg_fill_lesions<T>::countLesionVoxels(int location1,int patchSize) {
    int countlesionvoxel=0;
    int shiftx=0;
    int shifty=0;
    int shiftz=0;
    for(shiftx=-patchSize; shiftx<=patchSize; shiftx++){
        for(shifty=-patchSize; shifty<=patchSize; shifty++){
            for(shiftz=-patchSize; shiftz<=patchSize; shiftz++){
                int shift=shiftx+this->getXAxisSize()*(shifty+(this->getYAxisSize()*shiftz));
                int index1=location1+shift;
                for(int tp=0;tp<this->getNumTP();tp++) {
                    if(index1<this->getSingleVolumSize()   &&  // Is within the index bounds (also checks z)
                            index1>=0  &&  // Is within the index bounds (also checks z)
                            this->isValidVoxelInAnyTP(index1) &&
                            this->tmpLesMask[index1+tp*this->getSingleVolumSize()]>0)
                    {
                        countlesionvoxel++;
                    }
                }
            }
        }
    }
    return countlesionvoxel;
}

#endif // _SEG_FILL_LESIONS_CPP
