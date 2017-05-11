/**
 * @file _seg_PatchMatch.cpp
 * @author Ferran Prados
 * @date 16/06/2016
 *
 * Copyright (c) 2016, University College London. All rights reserved.
 * Centre for Medical Image Computing (CMIC)
 * See the LICENSE.txt file in the nifty_seg root folder
 *
 */
 
#ifndef _SEG_PATCHMATCH_CPP
#define _SEG_PATCHMATCH_CPP

#ifdef _OPENMP
#include "omp.h"
#endif

#include "_seg_PatchMatch.h"

template <class T>
seg_PatchMatch<T>::seg_PatchMatch() {
    this->inputImage=NULL;
    this->inputMask=NULL;
    this->imgPtr=NULL;
    this->currSize=NULL;
    this->verbose=0;
    this->debug=0;
    this->result=NULL;
    this->outputMask=NULL;
    this->patchSize=3;
    this->distance=0;
    this->outputMaskPtr=NULL;
    this->count=0;
    this->maskPtr=NULL;
    this->matchingCount=NULL;
    this->filling=false;
}

template<class T>
void seg_PatchMatch<T>::setInputImage(nifti_image *input) {
    this->inputImage = input;
    this->init();
}


template<class T>
void seg_PatchMatch<T>::setInputImageDatabase(vector<string> input) {
    this->imageDatabaseFileList = input;
}

template<class T>
void seg_PatchMatch<T>::setInputMaskDatabase(vector<string> input) {
    this->maskDatabaseFileList = input;
}

template<class T>
void seg_PatchMatch<T>::setOutputFilesDatabase(vector<string> input) {
    this->outputDatabaseFileList = input;
}

template<class T>
void seg_PatchMatch<T>::setInputMask(nifti_image *input) {
    this->inputMask = input;

    // We initialize mask
    this->maskPtr = static_cast<T *>(this->inputMask->data);
}

template<class T>
float * seg_PatchMatch<T>::getOutputResult() {
    return this->outputMaskPtr;
}

template<class T>
void seg_PatchMatch<T>::init() {
    this->imgPtr = static_cast<T *>(this->inputImage->data);
    this->setNumTP(this->inputImage->nt);
    this->currSize = new ImageSize [1]();
    this->currSize->numel=(long)(this->inputImage->nx*this->inputImage->ny*this->inputImage->nz);
    this->currSize->xsize=this->inputImage->nx;
    this->currSize->ysize=this->inputImage->ny;
    this->currSize->zsize=this->inputImage->nz;
    this->currSize->usize=this->inputImage->nu;
    this->currSize->tsize=this->inputImage->nt;
}

template<class T>
void seg_PatchMatch<T>::setPatchSize(int value) {
    this->patchSize = value;
}

template<class T>
int seg_PatchMatch<T>::getPatchSize() {
    return this->patchSize;
}
template<class T>
void seg_PatchMatch<T>::setBetterMatch(int value) {
    this->better_match= value;
}

template<class T>
int seg_PatchMatch<T>::getBetterMatch() {
    return this->better_match;
}


template<class T>
void seg_PatchMatch<T>::setPatchMatchExecutions(int value) {
    this->pmatch_executions= value;
}

template<class T>
int seg_PatchMatch<T>::getPatchMatchExecutions() {
    return this->pmatch_executions;
}


template<class T>
void seg_PatchMatch<T>::setConstrainedSearchAreaSize(int value) {
    this->csa_size= value;
}

template<class T>
int seg_PatchMatch<T>::getConstrainedSearchAreaSize() {
    return this->csa_size;
}


template<class T>
void seg_PatchMatch<T>::setPatchMatchIterations(int value) {
    this->pm_iter= value;
}

template<class T>
int seg_PatchMatch<T>::getPatchMatchIterations() {
    return this->pm_iter;
}

template<class T>
void seg_PatchMatch<T>::setDistance(int value) {
    this->distance = value;
}

template<class T>
int seg_PatchMatch<T>::getDistance() {
    return this->distance;
}

template<class T>
void seg_PatchMatch<T>::setVerbose(int value) {
    this->verbose = value;
}

template<class T>
int seg_PatchMatch<T>::getVerbose() {
    return this->verbose;
}

template<class T>
void seg_PatchMatch<T>::setDebug(int value) {
    this->debug = value;
}

template<class T>
int seg_PatchMatch<T>::getDebug() {
    return this->debug;
}

template<class T>
void seg_PatchMatch<T>::setNumTP(int value) {
    this->numTP = value;
}

template<class T>
int seg_PatchMatch<T>::getNumTP() {
    return this->numTP;
}

template<class T>
long seg_PatchMatch<T>::getSingleVolumSize() {
    return this->currSize->numel;
}

template<class T>
long seg_PatchMatch<T>::getTotalVolumSize() {
    return this->currSize->numel*this->getNumTP();
}

template<class T>
long seg_PatchMatch<T>::getSliceSize() {
    return this->currSize->xsize*this->currSize->ysize;
}

template<class T>
long seg_PatchMatch<T>::getXAxisSize() {
    return this->currSize->xsize;
}

template<class T>
long seg_PatchMatch<T>::getYAxisSize() {
    return this->currSize->ysize;
}

template<class T>
long seg_PatchMatch<T>::getZAxisSize() {
    return this->currSize->zsize;
}

template<class T>
void seg_PatchMatch<T>::setFilling(bool v) {
    this->filling=v;
}

template<class T>
bool seg_PatchMatch<T>::isFilling() {
    return this->filling;
}

template<class T>
float seg_PatchMatch<T>::calculateSSDDistance(long location1, long location2,int image,int tp,float actual_distance){

    float distance=0;
    int count=0;
    int shiftx=0;
    int shifty=0;
    int shiftz=0;
    int pSize=(int)((this->getPatchSize()-1)/2);
    int maxshift=pSize+this->getXAxisSize()*(pSize+(this->getYAxisSize()*pSize));

    bool WillTouchBorder=1;
    if( (location1+maxshift)<this->getSingleVolumSize() &&
            (location1-maxshift)>=0 &&
            (location2+maxshift)<this->getSingleVolumSize() &&
            (location2-maxshift)>=0){
        WillTouchBorder=0;
    }
    for(shiftx=-pSize; shiftx<=pSize; shiftx++){
        for(shifty=-pSize; shifty<=pSize; shifty++){
            for(shiftz=-pSize; shiftz<=pSize; shiftz++){
                int shift=shiftx+this->getXAxisSize()*shifty+this->getXAxisSize()*this->getYAxisSize()*shiftz;
                int index1=location1+shift;
                int index2=location2+shift;
                long pos=index1+tp*this->getSingleVolumSize();
                long posSTP=index1;
                long posDB=image*this->getTotalVolumSize()+index2+tp*this->getSingleVolumSize();
                long posDB_STP=image*this->getSingleVolumSize()+index2;
                if(WillTouchBorder)
                {
                    if(index1<this->getSingleVolumSize()   &&  // Is within the index bounds (also checks z)
                            index1>=0  &&  // Is within the index bounds (also checks z)

                            index2<this->getSingleVolumSize()   &&  // Is within the index bounds (also checks z)
                            index2>=0       )   // Is within the index bounds (also checks z)
                    {
                        if(this->db_list_mask[posDB_STP]>0 && this->maskPtr[posSTP]>0) {
                            distance+=(this->imgPtr[pos]-this->db_list_image[posDB])*(this->imgPtr[pos]-this->db_list_image[posDB]);
                            count++;
                        }
                    }
                }
                else{
                    if(this->db_list_mask[posDB_STP]>0 && this->maskPtr[posSTP]>0) {
                        distance+=(this->imgPtr[pos]-this->db_list_image[posDB])*(this->imgPtr[pos]-this->db_list_image[posDB]);
                        count++;
                    }
                }
                if(distance>actual_distance) {
                    return std::numeric_limits<float>::quiet_NaN();
                }
            }
        }
    }
    return count==0?std::numeric_limits<float>::quiet_NaN():distance;
}

template<class T>
float seg_PatchMatch<T>::calculateSSDistance1(long location1,long location2,int image, float dist_old){
    float distance=0;
    // We compute the SSD of the two patches
    for(int tp=0;tp<this->getNumTP();tp++) {
        float dist=this->calculateSSDDistance(location1,location2,image,tp,dist_old);
        if(isnan(dist)==1 ){
            return std::numeric_limits<float>::quiet_NaN();
        }
        distance+=dist;
        if(distance>dist_old) {
            return std::numeric_limits<float>::quiet_NaN();
        }
    }

    return distance;
}

template<class T>
float seg_PatchMatch<T>::calculateSSDDistanceSliced(long location1, long location2,int image,int tp,int *pSizeX,int *pSizeY,int *pSizeZ){

    float distance=0;
    int count=0;
    int shiftx=0;
    int shifty=0;
    int shiftz=0;
    int pSize=(int)((this->getPatchSize()-1)/2);
    int maxshift=pSize+this->getXAxisSize()*(pSize+(this->getYAxisSize()*pSize));

    bool WillTouchBorder=1;
    if( (location1+maxshift)<this->getSingleVolumSize() &&
            (location1-maxshift)>=0 &&
            (location2+maxshift)<this->getSingleVolumSize() &&
            (location2-maxshift)>=0){
        WillTouchBorder=0;
    }
    for(shiftx=pSizeX[0]; shiftx<=pSizeX[1]; shiftx++){
        for(shifty=pSizeY[0]; shifty<=pSizeY[1]; shifty++){
            for(shiftz=pSizeZ[0]; shiftz<=pSizeZ[1]; shiftz++){
                int shift=shiftx+this->getXAxisSize()*shifty+this->getXAxisSize()*this->getYAxisSize()*shiftz;
                int index1=location1+shift;
                int index2=location2+shift;
                long pos=index1+tp*this->getSingleVolumSize();
                long posSTP=index1;
                long posDB=image*this->getTotalVolumSize()+index2+tp*this->getSingleVolumSize();
                long posDB_STP=image*this->getSingleVolumSize()+index2;
                if(WillTouchBorder)
                {
                    if(index1<this->getSingleVolumSize()   &&  // Is within the index bounds (also checks z)
                            index1>=0  &&  // Is within the index bounds (also checks z)

                            index2<this->getSingleVolumSize()   &&  // Is within the index bounds (also checks z)
                            index2>=0       )   // Is within the index bounds (also checks z)
                    {
                        if(this->db_list_mask[posDB_STP]>0 && this->maskPtr[posSTP]>0) {
                            distance+=(this->imgPtr[pos]-this->db_list_image[posDB])*(this->imgPtr[pos]-this->db_list_image[posDB]);
                            count++;
                        }
                    }
                }
                else{
                    if(this->db_list_mask[posDB_STP]>0 && this->maskPtr[posSTP]>0) {
                        distance+=(this->imgPtr[pos]-this->db_list_image[posDB])*(this->imgPtr[pos]-this->db_list_image[posDB]);
                        count++;
                    }
                }
            }
        }
    }
    return count==0?std::numeric_limits<float>::quiet_NaN():distance;
}

template<class T>
float seg_PatchMatch<T>::calculateSSDistance1Sliced(long location1,long location1DB,int image1,long location2,long location2DB,int image2, float distance1,float distance2,int pSize[6][2]){
    // We compute the SSD of the two patches but just for the lesioned voxel
    for(int tp=0;tp<this->getNumTP();tp++) {
        float dist1=this->calculateSSDDistanceSliced(location1,location1DB,image1,tp,pSize[0],pSize[1],pSize[2]);
        float dist2=this->calculateSSDDistanceSliced(location2,location2DB,image2,tp,pSize[3],pSize[4],pSize[5]);
        if(isnan(dist1)==1 ||
           isnan(dist2)==1 ) {
            return std::numeric_limits<float>::quiet_NaN();
        }
        distance2-=(dist2-dist1);
    }
    if(distance2>distance1) {
        return std::numeric_limits<float>::quiet_NaN();
    }
    return distance2;
}

template<class T>
void seg_PatchMatch<T>::calculate(int location1, int location2,int image,int tp,float &sum,int &count,float &meanImg,float &stdImg,float &meanDB,float &stdDB){

    meanImg=0;
    stdImg=0;
    meanDB=0;
    stdDB=0;
    sum=0;
    count=0;
    int shiftx=0;
    int shifty=0;
    int shiftz=0;
    int pSize=(int)((this->getPatchSize()-1)/2);
    int maxshift=pSize+this->getXAxisSize()*(pSize+(this->getYAxisSize()*pSize));
    bool WillTouchBorder=1;
    if( (location1+maxshift)<this->getSingleVolumSize()  &&
            (location1-maxshift)>=0 &&
            (location2+maxshift)<this->getSingleVolumSize()  &&
            (location2-maxshift)>=0){
        WillTouchBorder=0;
    }
    for(shiftx=-pSize; shiftx<=pSize; shiftx++){
        for(shifty=-pSize; shifty<=pSize; shifty++){
            for(shiftz=-pSize; shiftz<=pSize; shiftz++){
                int shift=shiftx+this->getXAxisSize()*shifty+this->getXAxisSize()*this->getYAxisSize()*shiftz;
                int index1=location1+shift;
                int index2=location2+shift;
                long pos=index1+tp*this->getSingleVolumSize();
                long posSTP=index1;
                long posDB=image*this->getTotalVolumSize()+index2+tp*this->getSingleVolumSize();
                long posDB_STP=image*this->getSingleVolumSize()+index2;
                if(WillTouchBorder) {
                    if(index1<this->getSingleVolumSize()    &&  // Is within the index bounds (also checks z)
                            index1>=0  &&  // Is within the index bounds (also checks z)

                            index2<this->getSingleVolumSize()    &&  // Is within the index bounds (also checks z)
                            index2>=0       )   // Is within the index bounds (also checks z)
                    {
                        if(this->db_list_mask[posDB_STP]>0 &&
                           this->maskPtr[posSTP]>0) {
                            meanImg+=this->imgPtr[pos];
                            meanDB+=this->db_list_image[posDB];
                            sum+=(this->imgPtr[pos]*this->db_list_image[posDB]);
                            count++;
                        }
                    }
                }
                else {
                    if(this->db_list_mask[posDB_STP]>0 &&
                       this->maskPtr[posSTP]>0) {
                        meanImg+=this->imgPtr[pos];
                        meanDB+=this->db_list_image[posDB];
                        sum+=(this->imgPtr[pos]*this->db_list_image[posDB]);
                        count++;
                    }
                }
            }
        }
    }
    meanImg=meanImg/count;
    meanDB=meanDB/count;
    for(shiftx=-pSize; shiftx<=pSize; shiftx++){
        for(shifty=-pSize; shifty<=pSize; shifty++){
            for(shiftz=-pSize; shiftz<=pSize; shiftz++){
                int shift=shiftx+this->getXAxisSize()*shifty+this->getXAxisSize()*this->getYAxisSize()*shiftz;
                int index1=location1+shift;
                int index2=location2+shift;
                long pos=index1+tp*this->getSingleVolumSize();
                long posSTP=index1;
                long posDB=image*this->getTotalVolumSize()+index2+tp*this->getSingleVolumSize();
                long posDB_STP=image*this->getSingleVolumSize()+index2;
                if(WillTouchBorder) {
                    if(index1<this->getSingleVolumSize()   &&  // Is within the index bounds (also checks z)
                            index1>=0  &&  // Is within the index bounds (also checks z)

                            index2<this->getSingleVolumSize()   &&  // Is within the index bounds (also checks z)
                            index2>=0       )   // Is within the index bounds (also checks z)
                    {
                        if(this->db_list_mask[posDB_STP]>0 &&
                           this->maskPtr[posSTP]>0) {
                            stdImg+=(this->imgPtr[pos]-meanImg)*(this->imgPtr[pos]-meanImg);
                            stdDB+=(this->db_list_image[posDB]-meanDB)*(this->db_list_image[posDB]-meanDB);
                        }
                    }
                }
                else {
                    if(this->db_list_mask[posDB_STP]>0 &&
                       this->maskPtr[posSTP]>0) {
                        stdImg+=(this->imgPtr[pos]-meanImg)*(this->imgPtr[pos]-meanImg);
                        stdDB+=(this->db_list_image[posDB]-meanDB)*(this->db_list_image[posDB]-meanDB);
                    }
                }
            }
        }
    }
    stdImg=sqrt(stdImg/count)+0.0000001;
    stdDB=sqrt(stdDB/count)+0.0000001;
}

template<class T>
float seg_PatchMatch<T>::calculateLNCCDistance1(long location1,long location2,int image,float dist_old){

    float distance=0;
    // We compute the LNCC of the two patches but just for the lesioned voxel
    for(int tp=0;tp<this->getNumTP();tp++) {
        float sum=0;
        int numPatchVox=0;
        float meanImg=0;
        float stdImg=0;
        float meanDB=0;
        float stdDB=0;
        this->calculate(location1,location2,image,tp,sum,numPatchVox,meanImg,stdImg,meanDB,stdDB);

        float dist=((sum-meanImg*meanDB)/numPatchVox)/(stdImg*stdDB);
        distance+=dist;
        if(distance>dist_old) {
            return std::numeric_limits<float>::quiet_NaN();
        }
    }

    return distance;
}

template<class T>
void seg_PatchMatch<T>::normalizeImageIntesitiesOutliers(float newMin,float newMax,T *image) {
    long tp;
    #ifdef _OPENMP
    #pragma omp parallel for \
        shared(newMin,newMax,std::cout)\
        private(tp)
    #endif
    for(tp=0; tp<this->getNumTP(); tp++) {
        float *imgsort=new float [this->getSingleVolumSize()];
        for(long i=0; i<this->getSingleVolumSize(); i++) {
            imgsort[i]=image[i+tp*this->getSingleVolumSize()];
        }
        HeapSort(imgsort,this->getSingleVolumSize()-1);
        float max=imgsort[(int)(round((1-0.02)*(this->getSingleVolumSize()-1)))];
        float min=imgsort[(int)(round(0.02*(this->getSingleVolumSize()-1)))];
        if(this->getVerbose()) cout<<"[TP="<<tp<<"] MIN="<<min<<" MAX="<<max<<endl;
        for(long i=0; i<this->getSingleVolumSize(); i++) {
            if(min>image[i+tp*this->getSingleVolumSize()]) image[i+tp*this->getSingleVolumSize()]=min;
            if(max<image[i+tp*this->getSingleVolumSize()]) image[i+tp*this->getSingleVolumSize()]=max;
            image[i+tp*this->getSingleVolumSize()]=newMin+(image[i+tp*this->getSingleVolumSize()]-min)*(newMax-newMin)/(max-min);
        }
    }
}

template<class T>
void seg_PatchMatch<T>::normalizeImageIntesities(float newMin,float newMax,T *image) {
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
            max=max>image[i+tp*this->getSingleVolumSize()]?max:image[i+tp*this->getSingleVolumSize()];
            min=min<image[i+tp*this->getSingleVolumSize()]?min:image[i+tp*this->getSingleVolumSize()];
        }
        if(this->getVerbose()) cout<<"[TP="<<tp<<"] MIN="<<min<<" MAX="<<max<<endl;
        for(long i=0; i<this->getSingleVolumSize(); i++) {
            image[i+tp*this->getSingleVolumSize()]=newMin+(image[i+tp*this->getSingleVolumSize()]-min)*(newMax-newMin)/(max-min);
        }
    }
}

template<class T>
void seg_PatchMatch<T>::loadFile(nifti_image *&InputImage,string filename,float *&filePtr,ImageSize *&fileSize) {
    InputImage=nifti_image_read(filename.c_str(),true);
    if(InputImage == NULL)
    {
        fprintf(stderr,"* Error when reading the input database file: %s\n",filename.c_str());
        exit(-1);
    }
    if(InputImage->datatype!=NIFTI_TYPE_FLOAT32)
    {
        seg_changeDatatype<SegPrecisionTYPE>(InputImage);
    }
    filePtr = static_cast<T *>(InputImage->data);

    fileSize = new ImageSize [1]();
    fileSize->numel=(long)(InputImage->nx*InputImage->ny*InputImage->nz);
    fileSize->xsize=InputImage->nx;
    fileSize->ysize=InputImage->ny;
    fileSize->zsize=InputImage->nz;
    fileSize->usize=(InputImage->nu>1)?InputImage->nu:1;
    fileSize->tsize=(InputImage->nt>1)?InputImage->nt:1;
}

template<class T>
void seg_PatchMatch<T>::getNextRecordDatabase(int num,float *&imagePtr,float *&maskDBPtr) {
    char filename[100];

    // Load Image File
    string filename_file=this->imageDatabaseFileList[num];
    nifti_image *nifti_image_file=NULL,*nifti_mask_file=NULL;
    ImageSize *imageSize=NULL;

    if(this->getVerbose()) {
        cout <<"[DB="<<num<<"] Reading a new file from the database: "<<filename_file;
    }
    this->loadFile(nifti_image_file,filename_file,imagePtr,imageSize);

    // Load mask File
    string filename_mask=this->maskDatabaseFileList[num];
    ImageSize *maskSize=NULL;
    if(this->getVerbose()) cout << " and "<<filename_mask<<endl;
    this->loadFile(nifti_mask_file,filename_mask,maskDBPtr,maskSize);

    if(!(maskSize->xsize==imageSize->xsize && maskSize->ysize==imageSize->ysize && maskSize->zsize==imageSize->zsize))
    {
        cout << "ERROR: Database mask (2nd column from the input text file) "<< filename_mask << " is the wrong size  respect to database  image = ( "<<imageSize->xsize<<","
             <<imageSize->ysize<<","<<imageSize->zsize<<" )  mask = ( "<<maskSize->xsize<<","
            <<maskSize->ysize<<","<<maskSize->zsize<<" )"<<endl;
        exit(-1);
    }
    if(!(this->getXAxisSize()==imageSize->xsize && this->getYAxisSize()==imageSize->ysize && this->getZAxisSize()==imageSize->zsize && this->getNumTP()==imageSize->tsize))
    {
        cout << "ERROR: Database images (1st column from the input text file) have different size respect to input image = ( "<<this->getXAxisSize()<<","
             <<this->getYAxisSize()<<","<<this->getZAxisSize()<<","<<this->getNumTP()<<")  database image = ( "<<imageSize->xsize<<","
            <<imageSize->ysize<<","<<imageSize->zsize<<","<<imageSize->tsize<<" )"<<endl;
        exit(-1);
    }

    // Removing voxels outside the mask
    long nvox=imageSize->xsize*imageSize->ysize*imageSize->zsize;
    long tp=0;
    for(int i=0; i<nvox; i++) {
        if(maskDBPtr[i]<0.5) maskDBPtr[i]=0;
        else maskDBPtr[i]=1;
    }

    if(this->getVerbose()) {
        cout <<"[DB="<<num<<"] Applying mask to the Image"<<endl;
    }
    #ifdef _OPENMP
    #pragma omp parallel for \
        shared(imageSize,maskDBPtr,imagePtr,nvox)\
        private(tp)
    #endif
    for(tp=0;tp<imageSize->tsize;tp++) {
        for(int i=0; i<nvox; i++) {
            if(maskDBPtr[i]>0) imagePtr[i+tp*nvox]=imagePtr[i+tp*nvox];
            else imagePtr[i+tp*nvox]=0;
        }
    }

    // Normalizing intensities
    if(this->getVerbose()) {
        cout <<"[DB="<<num<<"] Normalizing Database Image"<<endl;
    }
    //this->normalizeImageIntesitiesOutliers(0.0f,1024.0f,imagePtr);
    if(this->getDebug()) {
        sprintf(filename,"segPatchMatch_normalized_database_image_%d.nii.gz",num);
        this->saveImagePtr(imagePtr,nifti_image_file,filename);
    }

    // Counting
    this->count++;
}

template<class T>
float seg_PatchMatch<T>::calculateDistance(long index,long patchM,int nImage,float current_distances){

    float curdist=-1;
    switch(this->getDistance()) {
        case 0: curdist=this->calculateSSDistance1(index,
                                                   patchM,
                                                   nImage,
                                                   current_distances);
                break;
        case 1: curdist=this->calculateLNCCDistance1(index,
                                                     patchM,
                                                     nImage,
                                                     current_distances);
                break;
        default: cout<<"ERROR: Distance not implemented!!!!"<<endl;
                 exit(-1);
    }
    return curdist;
}

template<class T>
void seg_PatchMatch<T>::loadInputDatabase(){
    int size=this->imageDatabaseFileList.size();
    // We will compare with the entire database
    this->db_list_image=new float[size*this->getTotalVolumSize()];
    this->db_list_mask=new float[size*this->getSingleVolumSize()];
    int num=0;
    #ifdef _OPENMP
    #pragma omp parallel for  shared(std::cout,size)\
        private(num)
    #endif
    for(num=0;num<size;num++) {
        if(this->getVerbose()) {
            cout <<"Loading image "<<(num+1)<<"/"<<size<<" from the database"<<endl;
        }
        float *databaseImagePtr=NULL;
        float *databaseMaskPtr=NULL;
        this->getNextRecordDatabase(num,databaseImagePtr,databaseMaskPtr);

        for(int tp=0;tp<this->getNumTP();tp++) {
            for(long ii=0; ii<this->getSingleVolumSize(); ii++) {
                this->db_list_image[ii+tp*this->getSingleVolumSize()+num*this->getTotalVolumSize()]=databaseImagePtr[ii+tp*this->getSingleVolumSize()];
            }
        }
        for(long ii=0;ii<this->getSingleVolumSize();ii++) {
            this->db_list_mask[ii+num*this->getSingleVolumSize()]=databaseMaskPtr[ii];
        }
    }
}

template<class T>
void seg_PatchMatch<T>::loadOuputDatabase(){
    // Load Image File
    string filename_file=this->outputDatabaseFileList[0];
    float *imagePtr=NULL;
    long size=this->outputDatabaseFileList.size();
    long num=0;
    if(this->getVerbose()) {
        cout <<"Loading output "<<(num+1)<<"/"<<size<<" from the database"<<endl;
    }
    this->loadFile(this->ouputImageSample,filename_file,imagePtr,this->outputSize);

    if(!(this->getXAxisSize()==this->outputSize->xsize && this->getYAxisSize()==this->outputSize->ysize && this->getZAxisSize()==this->outputSize->zsize))
    {
        cout << "ERROR: Output database images (3rd column from the input text file) have different size respect to input image = ( "<<this->getXAxisSize()<<","
             <<this->getYAxisSize()<<","<<this->getZAxisSize()<<" )  output database image = ( "<<this->outputSize->xsize<<","
            <<this->outputSize->ysize<<","<<this->outputSize->zsize<<" )"<<endl;
        exit(-1);
    }

    long NUMVOXTOTAL=this->getSingleVolumSize()*this->outputSize->tsize;
    this->db_list_output=new float[size*NUMVOXTOTAL];
    for(long ii=0; ii<NUMVOXTOTAL; ii++) {
        this->db_list_output[ii+num*NUMVOXTOTAL]=imagePtr[ii];
    }
    #ifdef _OPENMP
    #pragma omp parallel for \
        private(num)\
        shared(NUMVOXTOTAL,size,cout)
    #endif
    for(num=1;num<size;num++) {
        if(this->getVerbose()) {
            cout <<"Loading output "<<(num+1)<<"/"<<size<<" from the database"<<endl;
        }
        float *databaseImagePtr=NULL;
        string filename_file_DB=this->outputDatabaseFileList[num];
        nifti_image *nifti_image_file_DB=NULL;
        ImageSize *imageSize_DB=NULL;
        this->loadFile(nifti_image_file_DB,filename_file_DB,databaseImagePtr,imageSize_DB);

	/*if(this->isFilling()) {
	  this->normalizeImageIntesitiesOutliers(0.0f,1024.0f,databaseImagePtr);
	}*/
        for(long ii=0; ii<NUMVOXTOTAL; ii++) {
            this->db_list_output[ii+num*NUMVOXTOTAL]=databaseImagePtr[ii];
        }
    }
}

template<class T>
void seg_PatchMatch<T>::patchMatch(){
    int threads;
    #ifdef _OPENMP
    #pragma omp parallel for \
        private(threads)
    #endif
    for(threads=0;threads<this->getPatchMatchExecutions();threads++) {
        this->computePatchMatch(threads);
    }
}

template<class T>
void seg_PatchMatch<T>::sortKNNResults(){
    if(this->getVerbose()) {
        cout<<"Sorting the results"<<endl;
    }
    long index=0;
    #ifdef _OPENMP
    #pragma omp parallel for \
        private(index)
    #endif
    for(index=0;index<this->getSingleVolumSize();index++) {
        if(this->KNN[index]!=NULL && this->KNN[index]->getPatchMatch()>-1) {
            this->sortResults(index);
        }
    }
}

template<class T>
void seg_PatchMatch<T>::loadingInputData(){
    char filename[100];
    if(this->getVerbose()) {
        cout <<"File dimensions"<<endl;
        cout <<"Size ["<<this->getXAxisSize()<<","<<this->getYAxisSize()<<","<<this->getZAxisSize()<<"] TP="<<this->getNumTP()<<endl;
        cout <<"Normalizing Image"<<endl;
    }
    long tp;
    #ifdef _OPENMP
    #pragma omp parallel for \
        private(tp)
    #endif
    for(tp=0;tp<this->getNumTP();tp++) {
        for(int i=0; i<this->getSingleVolumSize(); i++) {
            if(this->maskPtr[i]<0.5) {
                this->maskPtr[i]=0;
                this->imgPtr[i+tp*this->getSingleVolumSize()]=0;
            }
            else {
                this->imgPtr[i+tp*this->getSingleVolumSize()]=this->imgPtr[i+tp*this->getSingleVolumSize()];
                this->maskPtr[i]=1;
            }
        }
    }
    // Intensity image normalization
    //this->normalizeImageIntesitiesOutliers(0.0f,1024.0f,this->imgPtr);
    if(this->getDebug()) {
        sprintf(filename,"segPatchMatch_normalized_image.nii.gz");
        this->saveImagePtr(this->imgPtr,this->inputImage,filename);
    }
}

template<class T>
void seg_PatchMatch<T>::initKNNVectors(){
    if(this->getVerbose()) {
        cout <<"Init KNN vector"<<endl;
    }
    this->KNN.clear();
    this->KNN.reserve(this->getTotalVolumSize()*this->getPatchMatchExecutions());
    if(this->getVerbose()) {
        cout <<"KNN vector ready"<<endl;
    }
}

template<class T>
void seg_PatchMatch<T>::labelFusion() {
    if(this->getVerbose()) {
        cout<<"Fusing the "<<this->getBetterMatch()<<" best matchings"<<endl;
    }
    long NUMVOXTOTAL=this->getSingleVolumSize()*this->outputSize->tsize;
    this->outputMaskPtr=new float[NUMVOXTOTAL];
    this->matchingCount=new int[NUMVOXTOTAL];
    long ii=0;

    #ifdef _OPENMP
    #pragma omp parallel for \
        private(ii)\
        shared(NUMVOXTOTAL)
    #endif
    for(ii=0; ii<NUMVOXTOTAL; ii++) {
        this->outputMaskPtr[ii]=0;
        this->matchingCount[ii]=0;
    }

    long matchings;
    #ifdef _OPENMP
    #pragma omp parallel for \
        private(matchings)\
        shared(NUMVOXTOTAL)
    #endif
    for(matchings=0;matchings<this->getBetterMatch();matchings++) {
        for(long i=0;i<this->getSingleVolumSize();i++) {
            long position=i+matchings*this->getTotalVolumSize();
            if(this->KNN[position]!=NULL && this->KNN[position]->getPatchMatch()>-1) {
                this->fusePatch(i,
                                this->KNN[position]->getPatchMatch(),
                                this->KNN[position]->getImage(),
                                NUMVOXTOTAL);
            }
        }
    }
}

template<class T>
void seg_PatchMatch<T>::saveResults() {
    char filename[100];
    long NUMVOXTOTAL=this->getSingleVolumSize()*this->outputSize->tsize;
    float *outMaskPtr=new float[NUMVOXTOTAL];
    long i;
    #ifdef _OPENMP
    #pragma omp parallel for \
        private(i)\
        shared(outMaskPtr,NUMVOXTOTAL)
    #endif
    for(i=0; i<NUMVOXTOTAL; i++) {
        if(this->matchingCount[i]>0) {
            outMaskPtr[i]=(float)(this->outputMaskPtr[i]/this->matchingCount[i]);
            this->outputMaskPtr[i]=outMaskPtr[i];
        }
        else {
            outMaskPtr[i]=0;
            this->outputMaskPtr[i]=0;
        }
    }
    if(this->getDebug()) {
        sprintf(filename,"segPatchMatch_intermidium_fused-patch.nii.gz");
        this->saveImagePtr(outMaskPtr,this->ouputImageSample,filename);
    }
}

template<class T>
void seg_PatchMatch<T>::fusePatch(long index,long location,int image,long NUMVOXTOTAL){
    int shiftx=0;
    int shifty=0;
    int shiftz=0;
    int pSize=(int)((this->getPatchSize()-1)/2);

    for(shiftx=-pSize; shiftx<=pSize; shiftx++){
        for(shifty=-pSize; shifty<=pSize; shifty++){
            for(shiftz=-pSize; shiftz<=pSize; shiftz++){
                for(int tp=0;tp<this->outputSize->tsize;tp++) {
                    int shift=shiftx+this->getXAxisSize()*shifty+this->getXAxisSize()*this->getYAxisSize()*shiftz;
                    long pos=index+shift+tp*this->getSingleVolumSize();
                    long posDB=location+shift+tp*this->getSingleVolumSize()+image*NUMVOXTOTAL;
                    if(index+shift>=0 &&
                       index+shift<this->getSingleVolumSize() &&
                       location+shift>=0 &&
                       location+shift<this->getSingleVolumSize() &&
                       (this->maskPtr[index+shift]>0 || this->isFilling())&&
                       this->db_list_mask[location+shift+image*this->getSingleVolumSize()]>0) {
                            this->outputMaskPtr[pos]+=this->db_list_output[posDB];
                            this->matchingCount[pos]++;
                    }
               }
            }
        }
    }
}

template<class T>
void seg_PatchMatch<T>::saveDebugResults() {
    char filename[100];
    float *ANN=new float[this->getTotalVolumSize()];
    int tp=0;
    #ifdef _OPENMP
    #pragma omp parallel for \
        private(tp)\
        shared(ANN)
    #endif
    for(tp=0;tp<this->getNumTP();tp++){
        for(long i=0;i<this->getSingleVolumSize();i++) {
            if(this->KNN[i]!=NULL) {
                if(this->KNN[i]->getANN()==std::numeric_limits<float>::max())
                    ANN[i+tp*this->getSingleVolumSize()]=0;
                else
                    ANN[i+tp*this->getSingleVolumSize()]=this->KNN[i]->getANN();
            }
            else {
                ANN[i+tp*this->getSingleVolumSize()]=0;
            }
        }
    }
    this->normalizeImageIntesities(0.0f,1024.0f,ANN);

    sprintf(filename,"segPatchMatch_intermidium_ANN-best.nii.gz");
    this->saveImagePtr(ANN,this->inputImage,filename);

    int *image=new int[this->getSingleVolumSize()];
    long i=0;
    #ifdef _OPENMP
    #pragma omp parallel for \
        private(i)\
        shared(image)
    #endif
    for(i=0;i<this->getSingleVolumSize();i++) {
        if(this->KNN[i]!=NULL) {
            if(this->KNN[i]->getANN()==std::numeric_limits<float>::max()) image[i]=0;
            else image[i]=this->KNN[i]->getImage();
         }
         else {
             image[i]=0;
         }
    }
    sprintf(filename,"segPatchMatch_intermidium_Image-best.nii.gz");
    this->saveImagePtr(image,this->inputMask,filename);
    
    sprintf(filename,"segPatchMatch_intermidium_matching_counting.nii.gz");
    this->saveImagePtr(this->matchingCount,this->inputMask,filename);
}

template<class T>
long seg_PatchMatch<T>::recomputeInputData(int iteration){
    char filename[100];
    if(this->getVerbose()) {
      cout<<"["<<iteration<<"] Filling"<<endl;
      cout<<"["<<iteration<<"] Recomputing/filling input data"<<endl;
    }
    long tp;
    long count=0;
    long painted=0;
    int average=(int)this->getPatchSize()*this->getPatchSize()*this->getPatchMatchExecutions()/2;
    for(int i=0; i<this->getSingleVolumSize(); i++) {
      if(this->maskPtr[i]<1) {
	for(tp=0;tp<this->getNumTP();tp++) {
	  if(this->matchingCount[i+tp*this->getSingleVolumSize()]>average) {
	    this->maskPtr[i]=1;
	    this->imgPtr[i+tp*this->getSingleVolumSize()]=this->outputMaskPtr[i+tp*this->getSingleVolumSize()]/this->matchingCount[i+tp*this->getSingleVolumSize()];
	  }
	  else {
	    count++;
	  }
        }
        painted++;
      }
    }
    if(this->getVerbose()) {
      cout<<"["<<iteration<<"] Number of voxels: "<<painted<<" "<<count<<endl;
    }
    if(this->getDebug()) {
        sprintf(filename,"segPatchMatch_recomputed_image_%d.nii.gz",iteration);
        this->saveImagePtr(this->imgPtr,this->inputImage,filename);
	sprintf(filename,"segPatchMatch_recomputed_mask_%d.nii.gz",iteration);
        this->saveImagePtr(this->maskPtr,this->inputMask,filename);
    }
    return count;
}

template<class T>
void seg_PatchMatch<T>::runIt(){

    // We load the input data
    this->loadingInputData();

    // Initializing results
    this->initKNNVectors();

    // We load the database
    this->loadInputDatabase();

    // We are going to compute PatchMatch
    this->patchMatch();

    // Sorting the results
    this->sortKNNResults();

    // We will load the output files
    this->loadOuputDatabase();

    // Fusing results taking into account KNN best matches and the output data
    this->labelFusion();
    
    // If we are doing filling we iterate until we fill all the voxels
    if(this->isFilling()) {
      int iteration=0;
      long voxels=this->recomputeInputData(iteration);
      while(voxels>0) {
	// Initializing results
	this->initKNNVectors();
	// We are going to compute PatchMatch
        this->patchMatch();
	// Sorting the results
	this->sortKNNResults();
	// Fusing results taking into account KNN best matches and the output data
	this->labelFusion();
	
	// We check if we need to recompute again the filling
	iteration++;
	voxels=this->recomputeInputData(iteration);
      }
    }

    // Saving results
    this->saveResults();

    // Saving debugging results
    if(this->getDebug()) {
        this->saveDebugResults();
    }

    if(this->getVerbose()) {
        cout<<"The patchmatch is done"<<endl;
    }

    return;
}

template<class T>
void seg_PatchMatch<T>::sortResults(long index){
    PatchMatchResult  *piv;
    std::vector<long> beg;
    std::vector<long> end;
    long i=0, L, R, swap ;

    beg.reserve(this->getPatchMatchExecutions());
    end.reserve(this->getPatchMatchExecutions());
    beg[0]=0;
    end[0]=this->getPatchMatchExecutions();
    while (i>=0) {
        L=beg[i];
        R=end[i]-1;
        if (L<R) {
            piv=this->KNN[index+L*this->getTotalVolumSize()];
            while (L<R) {
                while (this->KNN[index+R*this->getTotalVolumSize()]>=piv && L<R) {
                    R--;
                }
                if (L<R) {
                    this->KNN[index+L*this->getTotalVolumSize()]=this->KNN[index+R*this->getTotalVolumSize()];
                    L++;
                }
                while (this->KNN[index+L*this->getTotalVolumSize()]<=piv && L<R) {
                    L++;
                }
                if (L<R) {
                    this->KNN[index+R*this->getTotalVolumSize()]=this->KNN[index+L*this->getTotalVolumSize()];
                    R--;
               }
            }
            this->KNN[index+L*this->getTotalVolumSize()]=piv;
            beg[i+1]=L+1;
            end[i+1]=end[i];
            end[i++]=L;
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

template<class T>
void seg_PatchMatch<T>::computePatchMatch(int num){
    if(this->getVerbose()) {
        cout<<"["<<num<<"] Initializing data"<<endl;
    }
    srand(num);
    if(this->getVerbose()) {
        cout<<"["<<num<<"] Initializing Random Search"<<endl;
        cout <<"["<<num<<"] Size ["<<this->getXAxisSize()<<","<<this->getYAxisSize()<<","<<this->getZAxisSize()<<"] TP="<<this->getNumTP()<<endl;
    }
    long times=this->getConstrainedSearchAreaSize();

    int inz=0;
    #ifdef _OPENMP
    #pragma omp parallel for \
        private(inz)
    #endif
    // We init the PatchMatch Algorithm, first we need to init the random search
    for(inz=0; inz<this->getZAxisSize(); inz++) {
        for(int iny=0; iny<this->getYAxisSize(); iny++) {
            for(int inx=0; inx<this->getXAxisSize(); inx++) {
               long index=inx+this->getXAxisSize()*iny+this->getXAxisSize()*this->getYAxisSize()*inz;
               // If it is an important patch we will look for a similar patch
               if(this->maskPtr[index]>0) {
                    long position=index+num*this->getTotalVolumSize();
                    int image=-1;
                    long shift=0;
                    this->getRandomIndex(index,times,shift,image);
                    this->KNN[position]=new PatchMatchResult();
                    this->KNN[position]->setPatchMatch(index+shift);
                    this->KNN[position]->setImage(image);
                    float current_distances=this->KNN[position]->getANN();
                    float curdist=this->calculateDistance(index,
                                                          this->KNN[position]->getPatchMatch(),
                                                          this->KNN[position]->getImage(),
                                                          current_distances);
                    if(isnan(curdist)==0) {
                        this->KNN[position]->setANN(curdist);
                    }
                }
            }
        }
    }

    if(this->getVerbose()) {
        cout<<"["<<num<<"] Propagation step"<<endl;
    }
    int pSize=(int)((this->getPatchSize()-1)/2);
    const long vshift[6]={-1,1,
                          this->getXAxisSize()*-1,this->getXAxisSize()*1,
                          this->getXAxisSize()*this->getYAxisSize()*-1,this->getXAxisSize()*this->getYAxisSize()*1};
    int pSizeShift[6][6][2]={
                            {   // shift -1 -> X
                                { pSize, pSize},{-pSize, pSize},{-pSize, pSize}, //Input image
                                {-pSize,-pSize},{-pSize, pSize},{-pSize, pSize}  //Database image
                            },
                            {   // shift +1 -> X
                                {-pSize,-pSize},{-pSize, pSize},{-pSize, pSize},
                                { pSize, pSize},{-pSize, pSize},{-pSize, pSize},
                            },
                            {   // shift -1 -> Y
                                {-pSize, pSize},{ pSize, pSize},{-pSize, pSize},
                                {-pSize, pSize},{-pSize,-pSize},{-pSize, pSize}
                            },
                            {   // shift +1 -> Y
                                {-pSize, pSize},{-pSize,-pSize},{-pSize, pSize},
                                {-pSize, pSize},{ pSize, pSize},{-pSize, pSize}
                            },
                            {   // shift -1 -> Z
                                {-pSize, pSize},{-pSize, pSize},{ pSize, pSize},
                                {-pSize, pSize},{-pSize, pSize},{-pSize,-pSize}
                            },
                            {   // shift +1 -> Z
                                {-pSize, pSize},{-pSize, pSize},{-pSize,-pSize},
                                {-pSize, pSize},{-pSize, pSize},{ pSize, pSize}
                            }
    };
    for(int iterations=0;iterations<this->getPatchMatchIterations();iterations++) {
        if(this->getVerbose()) {
            cout<<"["<<num<<"] Iteration "<<(iterations+1)<<endl;
        }
        for(inz=0; inz<this->getZAxisSize(); inz++) {
            for(int iny=0; iny<this->getYAxisSize(); iny++) {
                for(int inx=0; inx<this->getXAxisSize(); inx++) {
                    long index=inx+this->getXAxisSize()*iny+this->getYAxisSize()*this->getXAxisSize()*inz;
                    long position=index+num*this->getTotalVolumSize();
                    if(this->KNN[position]!=NULL && this->KNN[position]->getPatchMatch()>-1){
                        for(int i=0;i<6;i++) {
                            this->propagationSTEP(num,index,vshift[i],pSizeShift[i]);
                        }
                        long init_times=times-1;
                        while(init_times>1) {
                            long shiftImage;
                            // Last parameter: Image is not used, we remain in
                            // the same image looking for a better patch randomly
                            // and constraining the search area
                            int image=this->KNN[position]->getImage();
                            this->getRandomIndex(this->KNN[position]->getPatchMatch(),
                                                 init_times,
                                                 shiftImage,
                                                 image);
                            this->constrainedRandomSearch(num,index,shiftImage);
                            init_times=init_times-1;
                        }
                    }
                }
            }
        }
    }
}

template<class T>
void seg_PatchMatch<T>::getRandomIndex(long index,long times,long &shift,int &image) {
    long shiftrealsize=(int)((this->getPatchSize()*times-1)/2);
    long shiftz=(rand() % (shiftrealsize*2+1))-shiftrealsize;
    long shifty=(rand() % (shiftrealsize*2+1))-shiftrealsize;
    long shiftx=(rand() % (shiftrealsize*2+1))-shiftrealsize;
    int imageNew=image;
    shift=shiftx+this->getXAxisSize()*shifty+this->getXAxisSize()*this->getYAxisSize()*shiftz;
    int total_image=this->count;
    if(image<0) {
        imageNew=(rand() % (total_image));
    }
    while(index+shift<0 ||
          index+shift>this->getSingleVolumSize() ||
          this->db_list_mask[imageNew*this->getSingleVolumSize()+index+shift]<1) {
        shiftz=(rand() % (shiftrealsize*2+1))-shiftrealsize;
        shifty=(rand() % (shiftrealsize*2+1))-shiftrealsize;
        shiftx=(rand() % (shiftrealsize*2+1))-shiftrealsize;
        shift=shiftx+this->getXAxisSize()*shifty+this->getXAxisSize()*this->getYAxisSize()*shiftz;
        if(image<0) {
            imageNew=(rand() % (total_image));
        }
    }
    image=imageNew;
}

template<class T>
void seg_PatchMatch<T>::propagationSTEP(int num,long location,long shift,int pSize[6][2]) {
    long locationNeighboor=location+shift;
    if( locationNeighboor>0 &&
        locationNeighboor<this->getSingleVolumSize() &&
        this->maskPtr[locationNeighboor]>0) {
            long posTrg=location+num*this->getTotalVolumSize();
            long posSrc=locationNeighboor+num*this->getTotalVolumSize();
            long positionNeighboor=this->KNN[posSrc]->getPatchMatch();
            long newPosition=positionNeighboor-shift;
            int newImage=this->KNN[posSrc]->getImage();
            if( newPosition>0 &&
                newPosition<this->getSingleVolumSize() &&
                this->db_list_mask[newPosition+newImage*this->getSingleVolumSize()]>0) {
                    float new_distances1=this->KNN[posTrg]->getANN();
                    float new_distances2=this->KNN[posSrc]->getANN();
                    float curdist=this->calculateSSDistance1Sliced(location,newPosition,newImage,locationNeighboor,positionNeighboor,newImage,new_distances1,new_distances2,pSize);
                    if(isnan(curdist)==0) {
                        this->KNN[posTrg]->setANN(curdist);
                        this->KNN[posTrg]->setPatchMatch(newPosition);
                        this->KNN[posTrg]->setImage(newImage);
                    }
            }
    }
}

template<class T>
void seg_PatchMatch<T>::constrainedRandomSearch(int num,long location,long shift) {
    long position=location+num*this->getTotalVolumSize();
    long newPosition=this->KNN[position]->getPatchMatch()+shift;
    int image=this->KNN[position]->getImage();
    if( newPosition>0 &&
        newPosition<this->getSingleVolumSize() &&
        this->db_list_mask[newPosition+image*this->getSingleVolumSize()]>0) {
            float new_distances=this->KNN[position]->getANN();
            float curdist=this->calculateDistance(location,newPosition,image,new_distances);
            if(isnan(curdist)==0) {
                this->KNN[position]->setANN(curdist);
                this->KNN[position]->setPatchMatch(newPosition);
                this->KNN[position]->setImage(image);
            }
    }
}
/* This function is not the most optimal */
template<class T>
void seg_PatchMatch<T>::propagationSTEP(int num,long location,long shift) {
    long locationNeighboor=location+shift;
    if( locationNeighboor>0 &&
        locationNeighboor<this->getSingleVolumSize() &&
        this->maskPtr[locationNeighboor]>0) {
            long posTrg=location+num*this->getTotalVolumSize();
            long posSrc=locationNeighboor+num*this->getTotalVolumSize();
            long newPosition=this->KNN[posSrc]->getPatchMatch()-shift;
            int newImage=this->KNN[posSrc]->getImage();
            if( newPosition>0 &&
                newPosition<this->getSingleVolumSize() &&
                this->db_list_mask[newPosition+newImage*this->getSingleVolumSize()]>0) {
                    float new_distances=this->KNN[posTrg]->getANN();
                    float curdist=this->calculateDistance(location,newPosition,newImage,new_distances);
                    if(isnan(curdist)==0) {
                        this->KNN[posTrg]->setANN(curdist);
                        this->KNN[posTrg]->setPatchMatch(newPosition);
                        this->KNN[posTrg]->setImage(newImage);
                    }
            }
    }
}

template<class T>
void seg_PatchMatch<T>::saveImagePtr(int * ImagePtr,nifti_image *Image,char *filename_out) {
    float *tmp= new float [Image->nvox];
    for(size_t i=0; i<Image->nvox; i++) {
        tmp[i]=(float)ImagePtr[i];
    }
    saveImagePtr(tmp,Image,filename_out);
    delete []tmp;
}

template<class T>
void seg_PatchMatch<T>::saveImagePtr(float * ImagePtr,nifti_image *Image,char *filename_out) {
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
    OutputImage->scl_inter=0;
    OutputImage->scl_slope=1;
    float max=std::numeric_limits<float>::min();
    float min=std::numeric_limits<float>::max();
    float * OutputImagePtr = static_cast<float *>(OutputImage->data);
    for(long i=0; i<(long)(Image->nx*Image->ny*Image->nz*Image->nt*Image->nu); i++)
    {
        max=max>ImagePtr[i]?max:ImagePtr[i];
        min=min<ImagePtr[i]?min:ImagePtr[i];
        OutputImagePtr[i]=ImagePtr[i];

    }
    OutputImage->cal_max=max;
    OutputImage->cal_min=min;
    nifti_image_write(OutputImage);
    nifti_image_free(OutputImage);
}

template<class T>
void seg_PatchMatch<T>::saveImage(nifti_image * Image,char *filename_out) {
    std::cout<<"Saving intermidium file: "<<filename_out<<std::endl;
    // saving output
    nifti_set_filenames(Image,filename_out,0,0);
    nifti_image_write(Image);
}

#endif // _seg_PatchMatch_CPP
