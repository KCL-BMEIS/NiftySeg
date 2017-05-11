#ifndef _SEG_FILL_LESIONS_OTHER_CPP
#define _SEG_FILL_LESIONS_OTHER_CPP

#ifdef _OPENMP
#include "omp.h"
#endif

#include "_seg_fill_lesions_other.h"

template <class T>
seg_fill_lesions_other<T>::seg_fill_lesions_other() {
    this->inputImage=NULL;
    this->normImage=NULL;
    this->imgPtr=NULL;
    this->lesMaskPtr=NULL;
    this->normImgPtr=NULL;
    this->currSize=NULL;
    this->verbose=0;
    this->debug=0;
    this->mult=0.1;
    this->patchSearchAreaSize=10;
    this->patchSize=4;
}

template <class T>
seg_fill_lesions_other<T>::seg_fill_lesions_other(const seg_fill_lesions_other &copy) {
    this->inputImage=copy->inputImage;
    this->normImage=copy->normImage;
    this->outputImage=copy->outputImage;
    this->imgPtr=copy->imgPtr;
    this->lesMaskPtr=copy->lesMaskPtr;
    this->normImgPtr=copy->normImgPtr;
    this->tmpLesMask=copy->tmpLesMask;
    this->tmpImg=copy->tmpImg;
    this->originalLesMask=copy->originalLesMask;
    this->tmpNorm=copy->tmpNorm;
    this->currSize=copy->currSize;
    this->verbose=copy->verbose;
    this->debug=copy->debug;
    this->mult=copy->mult;
    this->patchSearchAreaSize=copy->patchSearchAreaSize;
    this->patchSize=copy->patchSize;
    this->numTP=copy->numTP;
}

template <class T>
int seg_fill_lesions_other<T>::operator=(const seg_fill_lesions_other &copy) {
    this->inputImage=copy->inputImage;
    this->normImage=copy->normImage;
    this->outputImage=copy->outputImage;
    this->imgPtr=copy->imgPtr;
    this->lesMaskPtr=copy->lesMaskPtr;
    this->normImgPtr=copy->normImgPtr;
    this->tmpLesMask=copy->tmpLesMask;
    this->tmpImg=copy->tmpImg;
    this->originalLesMask=copy->originalLesMask;
    this->tmpNorm=copy->tmpNorm;
    this->currSize=copy->currSize;
    this->verbose=copy->verbose;
    this->debug=copy->debug;
    this->mult=copy->mult;
    this->patchSearchAreaSize=copy->patchSearchAreaSize;
    this->patchSize=copy->patchSize;
    this->numTP=copy->numTP;

    return *this;
}


template<class T>
void seg_fill_lesions_other<T>::setInputImage(nifti_image *input) {
   this->inputImage = input;
   this->init();
}

template<class T>
void seg_fill_lesions_other<T>::init() {
    this->imgPtr = static_cast<T *>(this->inputImage->data);
    this->setNumTP(this->inputImage->nt);
    this->currSize = new ImageSize [1]();
    this->currSize->numel=(long)(this->inputImage->nx*this->inputImage->ny*this->inputImage->nz);
    this->currSize->xsize=this->inputImage->nx;
    this->currSize->ysize=this->inputImage->ny;
    this->currSize->zsize=this->inputImage->nz;
    this->currSize->usize=this->inputImage->nu;
    this->currSize->tsize=this->inputImage->nt;
 
    // We initialize norm image
    this->normImage = nifti_copy_nim_info(this->inputImage);
    this->normImage->data = (void *) calloc(this->normImage->nvox, sizeof(T));
    this->normImgPtr = static_cast<T *>(this->normImage->data);
}

template<class T>
void seg_fill_lesions_other<T>::setInputLesionMask(nifti_image *mask) {
    this->inputLesionMask = mask;
    this->lesMaskPtr = static_cast<T *>(this->inputLesionMask->data);
}

template<class T>
void seg_fill_lesions_other<T>::setSearchArea(int value) {
    this->patchSearchAreaSize = value;
}

template<class T>
int seg_fill_lesions_other<T>::getSearchArea() {
    return this->patchSearchAreaSize;
}

template<class T>
void seg_fill_lesions_other<T>::setPatchSize(int value) {
    this->patchSize = value;
}

template<class T>
int seg_fill_lesions_other<T>::getPatchSize() {
    return this->patchSize;
}

template<class T>
void seg_fill_lesions_other<T>::setVerbose(int value) {
    this->verbose = value;
}

template<class T>
int seg_fill_lesions_other<T>::getVerbose() {
    return this->verbose;
}

template<class T>
void seg_fill_lesions_other<T>::setDebug(int value) {
    this->debug = value;
}

template<class T>
int seg_fill_lesions_other<T>::getDebug() {
    return this->debug;
}

template<class T>
void seg_fill_lesions_other<T>::setNumTP(int value) {
    this->numTP = value;
}

template<class T>
int seg_fill_lesions_other<T>::getNumTP() {
    return this->numTP;
}

template<class T>
long seg_fill_lesions_other<T>::getSingleVolumSize() {
    return this->currSize->numel;
}

template<class T>
long seg_fill_lesions_other<T>::getTotalVolumSize() {
    return this->currSize->numel*this->getNumTP();
}

template<class T>
long seg_fill_lesions_other<T>::getSliceSize() {
    return this->currSize->xsize*this->currSize->ysize;
}

template<class T>
long seg_fill_lesions_other<T>::getXAxisSize() {
    return this->currSize->xsize;
}

template<class T>
long seg_fill_lesions_other<T>::getYAxisSize() {
    return this->currSize->ysize;
}

template<class T>
long seg_fill_lesions_other<T>::getZAxisSize() {
    return this->currSize->zsize;
}

template<class T>
float seg_fill_lesions_other<T>::calculateDistance(int location1, int location2){

    const int numvox=this->getSingleVolumSize();
    float distance=0;
    int count=0;
    int shiftx=0;
    int shifty=0;
    int shiftz=0;
    int maxshift=this->getPatchSize()+this->getXAxisSize()*(this->getPatchSize()+(this->getYAxisSize()*this->getPatchSize()));

    bool WillTouchBorder=1;
    if( (location1+maxshift)<numvox &&
            (location1-maxshift)>=0 &&
            (location2+maxshift)<numvox &&
            (location2-maxshift)>=0){
        WillTouchBorder=0;
    }

    for(shiftx=-this->getPatchSize(); shiftx<=this->getPatchSize(); shiftx++){
        for(shifty=-this->getPatchSize(); shifty<=this->getPatchSize(); shifty++){
            for(shiftz=-this->getPatchSize(); shiftz<=this->getPatchSize(); shiftz++){
                int shift=(shiftx)+this->getXAxisSize()*(shifty+(this->getYAxisSize()*shiftz));
                int index1=location1+shift;
                int index2=location2+shift;
                for(int tp=0;tp<this->getNumTP();tp++) {
                    if(WillTouchBorder)
                    {
                        if(index1<numvox   &&  // Is within the index bounds (also checks z)
                                index1>=0  &&  // Is within the index bounds (also checks z)

                                index2<numvox   &&  // Is within the index bounds (also checks z)
                                index2>=0) 
                        {
                            distance+=(this->tmpNorm[index1+tp*numvox]-this->tmpNorm[index2+tp*numvox])*(this->tmpNorm[index1+tp*numvox]-this->tmpNorm[index2+tp*numvox]);
                            count++;
                        }
                    }
                    else{
                        distance+=(this->tmpNorm[index1+tp*numvox]-this->tmpNorm[index2+tp*numvox])*(this->tmpNorm[index1+tp*numvox]-this->tmpNorm[index2+tp*numvox]);
                        count++;
                    }

                }
            }
        }
    }
    return distance;
}

template<class T>
void seg_fill_lesions_other<T>::calculateEuclideanDistance(float *Distance) {
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
            Lable[i]=this->lesMaskPtr[i+tp*this->getSingleVolumSize()]>0;
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
void seg_fill_lesions_other<T>::normalizeImageIntesities(float newMin,float newMax) {
    long tp;
    #ifdef _OPENMP
    #pragma omp parallel for \
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
long seg_fill_lesions_other<T>::isMaskInAnyTP(long index,int *inputMask) {
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
void seg_fill_lesions_other<T>::runIt(){
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
    this->normalizeImageIntesities(0.0f,1.0f);

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
	this->tmpImg[i]=this->imgPtr[i];
	this->tmpNorm[i]=this->normImgPtr[i];
    }
    // Total of voxels to be analyzed
    for(i=0; i<this->getSingleVolumSize(); i++) {
        int tpf=this->isMaskInAnyTP(i,this->tmpLesMask);
        if(tpf>=0) countvox++;
    }
    
    if(this->getVerbose()) cout << "Calculating Euclidean distance"<<endl;
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

    if(this->getVerbose()) cout << "Filling lesions"<<endl;
    
    float smoothing_parameter[5]={0.9,0.7,0.5,0.3,0.1};
    for(int i_smooth=0;i_smooth<5;i_smooth++) {
      if(this->getVerbose()) {
            cout<< "Smoothing factor: "<<smoothing_parameter[i_smooth]<<endl;
      }
      for(tp=0;tp<this->getNumTP();tp++) {
            for(i=0; i<this->getSingleVolumSize(); i++) {
                this->tmpLesMask[i+tp*this->getSingleVolumSize()]=this->originalLesMask[i+tp*this->getSingleVolumSize()];
                this->lesMaskPtr[i+tp*this->getSingleVolumSize()]=this->originalLesMask[i+tp*this->getSingleVolumSize()];
            }
      }
      long count=1;
      int iteration=0;
      int lastpercent=0; 
      curcountvox=0;
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
            shared(inz,lastpercent,curcountvox,std::cout,level,iteration,countvox,i_smooth,smoothing_parameter)\
            private(index)\
            reduction(+:count)
        #endif
        for(inz=0; inz<this->getZAxisSize(); inz++) {
            if(lastpercent!=floor((float)(curcountvox)/(float)(countvox)*100.0f)){
                if(this->getVerbose()) cout<<"Filled voxels: "<<(float)(curcountvox)<<"/"<<(float)(countvox)<<endl;
                lastpercent=floor((float)(curcountvox)/(float)(countvox)*100.0f);
            }
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
                            (long)(index+this->getSliceSize())<this->getSingleVolumSize() &&
                            (long)(index-this->getSliceSize())>=0)

                    {   // Boundary for the time point affected
                        if(level[index]<1) level[index]=1;
                        if((this->tmpLesMask[(index+tpf*this->getSingleVolumSize())+1]==0 ||
                            this->tmpLesMask[(index+tpf*this->getSingleVolumSize())-1]==0 ||
                            this->tmpLesMask[(index+tpf*this->getSingleVolumSize())+this->getXAxisSize()]==0 ||
                            this->tmpLesMask[(index+tpf*this->getSingleVolumSize())-this->getXAxisSize()]==0 ||
                            this->tmpLesMask[(index+tpf*this->getSingleVolumSize())+this->getSliceSize()]==0 ||
                            this->tmpLesMask[(index+tpf*this->getSingleVolumSize())-this->getSliceSize()]==0)
                                ){
                            if(level[index]<3) level[index]=3;
                            curcountvox++;
                            count++;
			    for(int ii=0;ii<this->getNumTP();ii++) {
			      if(this->tmpLesMask[index+ii*this->getSingleVolumSize()]==1) {
				// If it is a lesion, then search
				float result_raw=0;
				float result_normalised=0;
				float sum_weight=0;
				for(shiftz=-this->getSearchArea(); shiftz<=this->getSearchArea(); shiftz++) {
				    for(shifty=-this->getSearchArea(); shifty<=this->getSearchArea(); shifty++) {
					for(shiftx=-this->getSearchArea(); shiftx<=this->getSearchArea(); shiftx++) {
					    index2=(inx+shiftx)+this->getXAxisSize()*(iny+shifty)+this->getSliceSize()*(inz+shiftz);
					    if(index2<this->getSingleVolumSize() && // Is within the index bounds (also checks z)
						    index2>=0 && // Is within the index bounds (also checks z)
						    (inx+shiftx) < this->getXAxisSize() && // Does it flip around the image
						    (inx-shiftx) >= 0 && // Does it flip around the image
						    (iny+shifty) < this->getYAxisSize() && // Does it flip around the image
						    (iny-shifty) >= 0 &&     // Does it flip around the image
						    index!=index2  
					      )
					    {
						curdist=this->calculateDistance(index,index2);
						float weight=exp(-curdist/(smoothing_parameter[i_smooth])); // Because it is h^2 by itself in the paper at page 3, right column
						result_raw+=weight*this->tmpImg[index2+ii*this->getSingleVolumSize()];
						result_normalised+=weight*this->tmpNorm[index2+ii*this->getSingleVolumSize()];
						sum_weight+=weight;
					    }
					    else {
						if(level[index]<4) level[index]=4;
					    }
					}
				    }
				} 
				// We replace the values only for the timepoints that are masked.
				this->imgPtr[index+ii*this->getSingleVolumSize()]=result_raw/sum_weight;
				this->normImgPtr[index+ii*this->getSingleVolumSize()]=result_normalised/sum_weight;
				this->lesMaskPtr[index+ii*this->getSingleVolumSize()]=0;
				if(level[index]<6) level[index]=6;
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
        if(this->getDebug()) {
            sprintf(filename,"segFillLesions_file-%d.nii.gz",iteration);
            this->saveImagePtr(this->tmpImg,this->inputImage,filename);
            sprintf(filename,"segFillLesions_filemask-%d.nii.gz",iteration);
            this->saveImagePtr(this->tmpLesMask,this->inputLesionMask,filename);
        }
      }
    }
    
    if(this->getVerbose()) {
        std::cout<<"Filled voxels: "<<(float)(curcountvox)<<"/"<<(float)(countvox)<<endl;
        std::cout << "Done"<<endl;
    }

    nifti_image_free(this->normImage);
    delete [] this->tmpLesMask;
    delete [] this->tmpImg;
    delete [] this->tmpNorm;

    return;
}

template<class T>
void seg_fill_lesions_other<T>::saveImagePtr(int * ImagePtr,nifti_image *Image,char *filename_out) {
    float *tmp= new float [Image->nvox*Image->nt];
    for(size_t i=0; i<Image->nvox*(size_t)Image->nt; i++) {
        tmp[i]=(float)ImagePtr[i];
    }
    saveImagePtr(tmp,Image,filename_out);
    delete [] tmp;
}

template<class T>
void seg_fill_lesions_other<T>::saveImagePtr(float * ImagePtr,nifti_image *Image,char *filename_out) {
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
void seg_fill_lesions_other<T>::saveImage(nifti_image * Image,char *filename_out) {
    std::cout<<"Saving intermidium file: "<<filename_out<<std::endl;
    // saving output
    nifti_set_filenames(Image,filename_out,0,0);
    nifti_image_write(Image);
}


#endif // _SEG_FILL_LESIONS_CPP
