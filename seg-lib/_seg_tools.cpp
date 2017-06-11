#include "_seg_tools.h"
#include <cfloat>



void GaussianFilter4D_cArray(segPrecisionTYPE * ShortData,
                             int * S2L,
                             int * L2S,
                             segPrecisionTYPE gauss_std,
                             ImageSize * CurrSizes)
{

    const int dims[8]={(int)((CurrSizes->xsize>1)+(CurrSizes->ysize>1)+(CurrSizes->zsize>1)+(CurrSizes->tsize>1)+(CurrSizes->usize>1)),
                       (int)(CurrSizes->xsize),(int)(CurrSizes->ysize),(int)(CurrSizes->zsize),(int)(CurrSizes->tsize),(int)(CurrSizes->usize),0,0};
    nifti_image *TMPimg = nifti_make_new_nim(dims,NIFTI_TYPE_FLOAT32,0);
    TMPimg->dim[0]=TMPimg->ndim=dims[0];
    TMPimg->dim[4]=TMPimg->nt=CurrSizes->tsize;
    TMPimg->dim[5]=TMPimg->nu=1;
    TMPimg->pixdim[5]=TMPimg->du=1;
    nifti_update_dims_from_array(TMPimg);
    TMPimg->datatype = NIFTI_TYPE_FLOAT32;
    TMPimg->nbyper = sizeof(float);
    TMPimg->data = (void *)(malloc(CurrSizes->numel*CurrSizes->numclass*sizeof(segPrecisionTYPE)));

    segPrecisionTYPE * TMPdataPtr=static_cast<segPrecisionTYPE *>(TMPimg->data);

    for(long curr4d=0; curr4d<CurrSizes->numclass; curr4d++)  //For Each Class
    {
        for(long index=0; index<(long)CurrSizes->numel; index++) //Every voxel is initially a NaN
        {
            TMPdataPtr[index+curr4d*CurrSizes->numel]=std::numeric_limits<float>::quiet_NaN();
        }
        for(long index=0; index<(long)CurrSizes->numelmasked; index++) //Copy Class to TMPdataPtr in LongFormat if inside S2L
        {
            TMPdataPtr[S2L[index]+curr4d*CurrSizes->numel]=ShortData[index+curr4d*CurrSizes->numelmasked];
        }

    }

    GaussianSmoothing5D_nifti(TMPimg,NULL,gauss_std);


    for(long curr4d=0; curr4d<CurrSizes->numclass; curr4d++)  //For Each Class
    {
        for(long index=0; index<(long)CurrSizes->numelmasked; index++) //Copy data from TMPdataPtr back to the shortData
        {
            ShortData[index+curr4d*CurrSizes->numelmasked]=TMPdataPtr[S2L[index]+curr4d*CurrSizes->numel];
        }
    }

    nifti_image_free(TMPimg);
    return;
}

void GaussianFilter4D_cArray(segPrecisionTYPE * LongData,
                             segPrecisionTYPE gauss_std,
                             ImageSize * CurrSizes)
{

    const int dims[8]={(int)((CurrSizes->xsize>1)+(CurrSizes->ysize>1)+(CurrSizes->zsize>1)+(CurrSizes->tsize>1)+(CurrSizes->usize>1)),
                       (int)(CurrSizes->xsize),(int)(CurrSizes->ysize),(int)(CurrSizes->zsize),(int)(CurrSizes->tsize),(int)(CurrSizes->usize),0,0};
    nifti_image *TMPimg = nifti_make_new_nim(dims,NIFTI_TYPE_FLOAT32,0);
    TMPimg->dim[0]=TMPimg->ndim=3;
    TMPimg->dim[4]=TMPimg->nt=CurrSizes->tsize;
    TMPimg->dim[5]=TMPimg->nu=1;
    TMPimg->pixdim[5]=TMPimg->du=1;
    nifti_update_dims_from_array(TMPimg);
    TMPimg->datatype = NIFTI_TYPE_FLOAT32;
    TMPimg->nbyper = sizeof(float);
    TMPimg->data = (void *)LongData;
    GaussianSmoothing5D_nifti(TMPimg,NULL,gauss_std);
    TMPimg->data=NULL;
    nifti_image_free(TMPimg);
    return;
}




void GaussianSmoothing5D_nifti(nifti_image * Data,int * mask,float gauss_std_in)
{

    int nx=Data->nx;
    int ny=Data->ny;
    int nz=Data->nz;

    float * ImageBuffer= new float [nx*ny*nz]();
    float * Density= new float [nx*ny*nz]();
    float * DensityBuffer= new float [nx*ny*nz]();

    float * DataPTR=static_cast<float *>(Data->data);

    int shiftdirection[3]= {1,nx,nx*ny};
    int dim_array[3]= {nx,ny,nz};
    float dist_array[3]= {Data->dx,Data->dy,Data->dz};

    int numel=nx*ny*nz;
    long number4D5D=std::max(Data->nt,1)*std::max(Data->nu,1);
    for(long curr4D5Dindex=0; curr4D5Dindex<number4D5D; curr4D5Dindex++)  //For Each TP and modality
    {
        long current_4dShift_short=curr4D5Dindex*nx*ny*nz;

        // Masking density and area
        int i=0;
#ifdef _OPENMP
#pragma omp parallel for default(none) \
    shared(DataPTR,mask,current_4dShift_short,nx,ny,nz,Density) \
    private(i)
#endif
        for(i=0; i<nx*ny*nz; i++)
        {
            if(mask!=NULL){
                Density[i]=(mask[i]>0 && isnan(DataPTR[i+current_4dShift_short])==0);
            }
            else{
                Density[i]=(isnan(DataPTR[i+current_4dShift_short])==0);
            }
            DataPTR[i+current_4dShift_short]=(Density[i]>0)?DataPTR[i+current_4dShift_short]:0;
        }

        //Blur Buffer and density along each direction
        for(long currentdirection=0; currentdirection<3; currentdirection++)
        {

            // Setup kernel
            float gauss_std=gauss_std_in>0?gauss_std_in:fabs(gauss_std_in/(float)dist_array[currentdirection]);

            int kernelsize=(int)(gauss_std*6.0f) % 2 != 0 ? (int)(gauss_std*6.0f) : (int)(gauss_std*6.0f)+1;
            int kernelshift=(int)(kernelsize/2.0f);
            float GaussKernel[400]= {0};
            float kernelsum=0;
            for(i=0; i<kernelsize; i++)
            {
                GaussKernel[i]=expf((float)(-0.5f*powf((float)((float)i-(float)kernelshift)/gauss_std, 2.0f)))/(sqrtf(2.0f*3.14159265*powf(gauss_std, 2)));
                kernelsum+=GaussKernel[i];
            }
            for(i=0; i<kernelsize; i++)
                GaussKernel[i]/=kernelsum;

            // Updating buffers
            int index;
#ifdef _OPENMP
#pragma omp parallel for default(none) \
    shared(DataPTR,Density,ImageBuffer,DensityBuffer,current_4dShift_short,numel) \
    private(i)
#endif
            for(index=0; index<numel; index++)
            {
                ImageBuffer[index]=DataPTR[index+current_4dShift_short];
                DensityBuffer[index]=Density[index];
            }

            // openmp defines

            float TmpDataConvolution=0;
            float TmpMaskConvolution=0;
            float TmpKernDensity=0;
            int shiftstart=0;
            int shiftstop=0;
            int shift=0;
            int xyzpos[3];
            int xyzpos2;
            // For every pixel

#ifdef _OPENMP
#pragma omp parallel for default(none) \
    shared(DataPTR,GaussKernel,ImageBuffer,DensityBuffer,nx,ny,nz,kernelshift,Density,dim_array,currentdirection,shiftdirection,current_4dShift_short) \
    private(xyzpos2,TmpDataConvolution,TmpMaskConvolution,TmpKernDensity,index,shiftstart,shiftstop,shift,xyzpos)
#endif
            for(xyzpos2=0; xyzpos2<nz; xyzpos2++)
            {
                xyzpos[2]=xyzpos2;
                for(xyzpos[1]=0; xyzpos[1]<ny; xyzpos[1]++)
                {
                    for(xyzpos[0]=0; xyzpos[0]<nx; xyzpos[0]++)
                    {

                        TmpDataConvolution=0;
                        TmpMaskConvolution=0;
                        TmpKernDensity=0;
                        index=xyzpos[0]+(xyzpos[1]+xyzpos[2]*ny)*nx;
                        // Calculate allowed kernel shifts
                        shiftstart=((xyzpos[currentdirection]<kernelshift)?-xyzpos[currentdirection]:-kernelshift);
                        shiftstop=((xyzpos[currentdirection]>=(dim_array[currentdirection]-kernelshift))?(int)dim_array[currentdirection]-xyzpos[currentdirection]-1:kernelshift);

                        for(shift=shiftstart; shift<=shiftstop; shift++)
                        {
                            // Data Blur
                            TmpDataConvolution+=GaussKernel[shift+kernelshift]*ImageBuffer[index+shift*shiftdirection[currentdirection]];
                            // Mask Blur
                            TmpMaskConvolution+=GaussKernel[shift+kernelshift]*DensityBuffer[index+shift*shiftdirection[currentdirection]];
                            // Kernel density
                            TmpKernDensity+=GaussKernel[shift+kernelshift];
                        }
                        // Devide convolutions by the kernel density
                        DataPTR[index+current_4dShift_short]=TmpDataConvolution/TmpKernDensity;
                        Density[index]=TmpMaskConvolution/TmpKernDensity;

                    }
                }
            }
        }
        int index=0;
#ifdef _OPENMP
#pragma omp parallel for default(none) \
    shared(DataPTR, mask, Density,current_4dShift_short,numel) \
    private(index)
#endif
        for(index=0; index<numel; index++)
        {
            if(mask!=NULL){
                DataPTR[index+current_4dShift_short]=(mask[index]>0)?DataPTR[index+current_4dShift_short]/Density[index]:0;
            }
            else{
                DataPTR[index+current_4dShift_short]=DataPTR[index+current_4dShift_short]/Density[index];
            }
        }
    }

    delete [] ImageBuffer;
    delete [] Density;
    delete [] DensityBuffer;
    return;
}



void GaussianSmoothing4D_Nan_nifti(nifti_image * Data,nifti_image * maskNifti)
{

    float gauss_std_in=0.7f;
    if(maskNifti->datatype!=DT_FLOAT32)
    {
        seg_changeDatatype<float>(maskNifti);
    }
    float * mask = static_cast<float *>(maskNifti->data);

    int nx=Data->nx;
    int ny=Data->ny;
    int nz=Data->nz;

    float * ImageBuffer= new float [nx*ny*nz]();
    float * Density= new float [nx*ny*nz]();
    float * DensityBuffer= new float [nx*ny*nz]();


    float * DataPTR=static_cast<float *>(Data->data);

    int shiftdirection[3]= {1,nx,nx*ny};
    int dim_array[3]= {nx,ny,nz};
    float dist_array[3]= {Data->dx,Data->dy,Data->dz};

    int numel=nx*ny*nz;
    for(long curr4d=0; curr4d<(long)Data->nt; curr4d++)  //For Each TP
    {
        long current_4dShift_short=curr4d*nx*ny*nz;

        // Masking density and area
        int i=0;
#ifdef _OPENMP
#pragma omp parallel for default(none) \
    shared(DataPTR,mask,current_4dShift_short,nx,ny,nz,Density) \
    private(i)
#endif
        for(i=0; i<nx*ny*nz; i++)
        {
            if(mask!=NULL){
                Density[i]=(mask[i]>0 && isnan(DataPTR[i+current_4dShift_short])==0);
            }
            else{
                Density[i]=(isnan(DataPTR[i+current_4dShift_short])==0);
            }
            DataPTR[i+current_4dShift_short]=(Density[i]>0)?DataPTR[i+current_4dShift_short]:0;
        }

        //Blur Buffer and density along each direction
        for(long currentdirection=0; currentdirection<3; currentdirection++)
        {

            // Setup kernel
            float gauss_std=gauss_std_in>0?gauss_std_in:fabs(gauss_std_in/(float)dist_array[currentdirection]);

            int kernelsize=(int)(gauss_std*6.0f) % 2 != 0 ? (int)(gauss_std*6.0f) : (int)(gauss_std*6.0f)+1;
            int kernelshift=(int)(kernelsize/2.0f);
            float GaussKernel[400]= {0};
            float kernelsum=0;
            for(i=0; i<kernelsize; i++)
            {
                GaussKernel[i]=expf((float)(-0.5f*powf((float)((float)i-(float)kernelshift)/gauss_std, 2.0f)))/(sqrtf(2.0f*3.14159265*powf(gauss_std, 2)))+0.0000001;
                kernelsum+=GaussKernel[i];
            }
            for(i=0; i<kernelsize; i++)
                GaussKernel[i]/=kernelsum;

            // Updating buffers
            int index;
#ifdef _OPENMP
#pragma omp parallel for default(none) \
    shared(DataPTR,Density,ImageBuffer,DensityBuffer,current_4dShift_short,numel) \
    private(i)
#endif
            for(index=0; index<numel; index++)
            {
                ImageBuffer[index]=DataPTR[index+current_4dShift_short];
                DensityBuffer[index]=Density[index];
            }

            // openmp defines

            float TmpDataConvolution=0;
            float TmpMaskConvolution=0;
            float TmpKernDensity=0;
            int shiftstart=0;
            int shiftstop=0;
            int shift=0;
            int xyzpos[3];
            int xyzpos2;
            // For every pixel
#ifdef _OPENMP
#pragma omp parallel for default(none) \
    shared(DataPTR,GaussKernel,ImageBuffer,DensityBuffer,nx,ny,nz,kernelshift,Density,dim_array,currentdirection,shiftdirection,current_4dShift_short) \
    private(xyzpos2,TmpDataConvolution,TmpMaskConvolution,TmpKernDensity,index,shiftstart,shiftstop,shift,xyzpos)
#endif
            for(xyzpos2=0; xyzpos2<nz; xyzpos2++)
            {
                xyzpos[2]=xyzpos2;
                for(xyzpos[1]=0; xyzpos[1]<ny; xyzpos[1]++)
                {
                    for(xyzpos[0]=0; xyzpos[0]<nx; xyzpos[0]++)
                    {

                        TmpDataConvolution=0;
                        TmpMaskConvolution=0;
                        TmpKernDensity=0;
                        index=xyzpos[0]+(xyzpos[1]+xyzpos[2]*ny)*nx;
                        // Calculate allowed kernel shifts
                        shiftstart=((xyzpos[currentdirection]<kernelshift)?-xyzpos[currentdirection]:-kernelshift);
                        shiftstop=((xyzpos[currentdirection]>=(dim_array[currentdirection]-kernelshift))?(int)dim_array[currentdirection]-xyzpos[currentdirection]-1:kernelshift);

                        for(shift=shiftstart; shift<=shiftstop; shift++)
                        {
                            // Data Blur
                            TmpDataConvolution+=GaussKernel[shift+kernelshift]*ImageBuffer[index+shift*shiftdirection[currentdirection]];
                            // Mask Blur
                            TmpMaskConvolution+=GaussKernel[shift+kernelshift]*DensityBuffer[index+shift*shiftdirection[currentdirection]];
                            // Kernel density
                            TmpKernDensity+=GaussKernel[shift+kernelshift];
                        }
                        // Devide convolutions by the kernel density
                        DataPTR[index+current_4dShift_short]=TmpDataConvolution/TmpKernDensity;
                        Density[index]=TmpMaskConvolution/TmpKernDensity;

                    }
                }
            }
        }
        int index=0;
#ifdef _OPENMP
#pragma omp parallel for default(none) \
    shared(DataPTR, mask, Density,current_4dShift_short,numel) \
    private(index)
#endif
        for(index=0; index<numel; index++)
        {
            if(mask!=NULL){
                DataPTR[index+current_4dShift_short]=(mask[index]>0 && Density[index] >0)?DataPTR[index+current_4dShift_short]/Density[index]:std::numeric_limits<float>::quiet_NaN();
            }
            else{
                DataPTR[index+current_4dShift_short]=(Density[index] >0)?DataPTR[index+current_4dShift_short]/Density[index]:std::numeric_limits<float>::quiet_NaN();
            }
        }
    }

    delete [] ImageBuffer;
    delete [] Density;
    delete [] DensityBuffer;
    return;
}


void BlockSmoothing(nifti_image * Data,int * mask,int side_size)
{

    int nx=Data->nx;
    int ny=Data->ny;
    int nz=Data->nz;

    float * ImageBuffer= new float [nx*ny*nz]();
    float * Density= new float [nx*ny*nz]();
    float * DensityBuffer= new float [nx*ny*nz]();

    float * DataPTR=static_cast<float *>(Data->data);

    int shiftdirection[3]= {1,nx,nx*ny};
    int dim_array[3]= {nx,ny,nz};

    int numel=nx*ny*nz;
    for(long curr4d=0; curr4d<(long)Data->nt; curr4d++)  //For Each TP
    {
        long current_4dShift_short=curr4d*nx*ny*nz;

        // Masking density and area
        int i=0;
#ifdef _OPENMP
#pragma omp parallel for default(none) \
    shared(DataPTR,mask,current_4dShift_short,nx,ny,nz,Density) \
    private(i)
#endif
        for(i=0; i<nx*ny*nz; i++)
        {
            if(mask!=NULL){
                Density[i]=(mask[i]>0 && isnan(DataPTR[i+current_4dShift_short])==0);
            }
            else{
                Density[i]=(isnan(DataPTR[i+current_4dShift_short])==0);
            }
            DataPTR[i+current_4dShift_short]=(Density[i]>0)?DataPTR[i+current_4dShift_short]:0;
        }

        //Blur Buffer and density along each direction
        for(long currentdirection=0; currentdirection<3; currentdirection++)
        {

            // Setup kernel
            int kernelshift=(int)floor((float)(side_size)/2.0f);

            // Updating buffers
            int index;
#ifdef _OPENMP
#pragma omp parallel for default(none) \
    shared(DataPTR,Density,ImageBuffer,DensityBuffer,current_4dShift_short,numel) \
    private(i)
#endif
            for(index=0; index<numel; index++)
            {
                ImageBuffer[index]=DataPTR[index+current_4dShift_short];
                DensityBuffer[index]=Density[index];
            }

            // openmp defines

            float TmpDataConvolution=0;
            float TmpMaskConvolution=0;
            int shiftstart=0;
            int shiftstop=0;
            int shift=0;
            int xyzpos[3];
            int xyzpos2;
            // For every pixel
#ifdef _OPENMP
#pragma omp parallel for default(none) \
    shared(DataPTR,ImageBuffer,DensityBuffer,nx,ny,nz,kernelshift,Density,dim_array,currentdirection,shiftdirection,current_4dShift_short) \
    private(xyzpos2,TmpDataConvolution,TmpMaskConvolution,index,shiftstart,shiftstop,shift,xyzpos)
#endif
            for(xyzpos2=0; xyzpos2<nz; xyzpos2++)
            {
                xyzpos[2]=xyzpos2;
                for(xyzpos[1]=0; xyzpos[1]<ny; xyzpos[1]++)
                {
                    for(xyzpos[0]=0; xyzpos[0]<nx; xyzpos[0]++)
                    {

                        TmpDataConvolution=0;
                        TmpMaskConvolution=0;
                        index=xyzpos[0]+(xyzpos[1]+xyzpos[2]*ny)*nx;
                        // Calculate allowed kernel shifts
                        shiftstart=((xyzpos[currentdirection]<kernelshift)?-xyzpos[currentdirection]:-kernelshift);
                        shiftstop=((xyzpos[currentdirection]>=(dim_array[currentdirection]-kernelshift))?(int)dim_array[currentdirection]-xyzpos[currentdirection]-1:kernelshift);

                        for(shift=shiftstart; shift<=shiftstop; shift++)
                        {
                            // Data Blur
                            TmpDataConvolution+=ImageBuffer[index+shift*shiftdirection[currentdirection]];
                            // Mask Blur
                            TmpMaskConvolution+=DensityBuffer[index+shift*shiftdirection[currentdirection]];
                        }
                        // Devide convolutions by the kernel density
                        DataPTR[index+current_4dShift_short]=TmpDataConvolution;
                        Density[index]=TmpMaskConvolution;

                    }
                }
            }
        }
        int index=0;
#ifdef _OPENMP
#pragma omp parallel for default(none) \
    shared(DataPTR, mask, Density,current_4dShift_short,numel) \
    private(index)
#endif
        for(index=0; index<numel; index++)
        {
            if(mask!=NULL){
                DataPTR[index+current_4dShift_short]=(mask[index]>0)?DataPTR[index+current_4dShift_short]/Density[index]:0;
            }
            else{
                DataPTR[index+current_4dShift_short]=DataPTR[index+current_4dShift_short]/Density[index];
            }
        }
    }

    delete [] ImageBuffer;
    delete [] Density;
    delete [] DensityBuffer;
    return;
}



void SmoothLab(float * DataPTR,float factor, ImageSize * Currentsize){

    typedef std::map <unsigned int, float> DataPointMap;
    typedef std::pair <unsigned int, float> DataPointPair;

    int nx=Currentsize->xsize;
    int ny=Currentsize->ysize;
    int nz=Currentsize->zsize;

    float * ImageBuffer= new float [nx*ny*nz]();
    float gauss_std=factor;
    float gaussX_std=gauss_std;
    float gaussY_std=gauss_std;
    float gaussZ_std=gauss_std;

    int shiftdirection[3]= {1,nx,nx*ny};
    int dim_array[3]= {nx,ny,nz};

    int numel=nx*ny*nz;
    for(long curr4d=0; curr4d<(long)Currentsize->tsize; curr4d++)  //For Each TP
    {
        long current_4dShift_short=curr4d*nx*ny*nz;

        int index=0;
        int xyzpos[3]={0};

#ifdef _OPENMP
#pragma omp parallel for default(none) \
    shared(DataPTR,ImageBuffer,nx,ny,nz,dim_array,gauss_std,shiftdirection,current_4dShift_short,gaussX_std,gaussY_std,gaussZ_std) \
    private(index,xyzpos)
#endif
        for(int zpos=0; zpos<nz; zpos++)
        {
            xyzpos[2]=zpos;
            for(xyzpos[1]=0; xyzpos[1]<ny; xyzpos[1]++)
            {
                for(xyzpos[0]=0; xyzpos[0]<nx; xyzpos[0]++)
                {

                    index=xyzpos[0]+(xyzpos[1]+xyzpos[2]*ny)*nx;

                    // Calculate allowed kernel shifts
                    int kernelXsize=(int)(gaussX_std*6.0f) % 2 != 0 ? (int)(gaussX_std*6.0f) : (int)(gaussX_std*6.0f)+1;
                    int kernelXshift=(int)(kernelXsize/2.0f);
                    int shiftXstart=((xyzpos[0]<kernelXshift)?-xyzpos[0]:-kernelXshift);
                    int shiftXstop=((xyzpos[0]>=(dim_array[0]-kernelXshift))?(int)dim_array[0]-xyzpos[0]-1:kernelXshift);

                    int kernelYsize=(int)(gaussY_std*6.0f) % 2 != 0 ? (int)(gaussY_std*6.0f) : (int)(gaussY_std*6.0f)+1;
                    int kernelYshift=(int)(kernelYsize/2.0f);
                    int shiftYstart=((xyzpos[1]<kernelYshift)?-xyzpos[1]:-kernelYshift);
                    int shiftYstop=((xyzpos[1]>=(dim_array[1]-kernelYshift))?(int)dim_array[1]-xyzpos[1]-1:kernelYshift);

                    int kernelZsize=(int)(gaussZ_std*6.0f) % 2 != 0 ? (int)(gaussZ_std*6.0f) : (int)(gaussZ_std*6.0f)+1;
                    int kernelZshift=(int)(kernelZsize/2.0f);
                    int shiftZstart=((xyzpos[2]<kernelZshift)?-xyzpos[2]:-kernelZshift);
                    int shiftZstop=((xyzpos[2]>=(dim_array[2]-kernelZshift))?(int)dim_array[2]-xyzpos[2]-1:kernelZshift);


                    DataPointMap tmp_lab;
                    for(int shiftx=shiftXstart; shiftx<=shiftXstop; shiftx++)
                    {
                        for(int shifty=shiftYstart; shifty<=shiftYstop; shifty++)
                        {
                            for(int shiftz=shiftZstart; shiftz<=shiftZstop; shiftz++)
                            {
                                // Data Blur
                                float kernelval=expf((float)(-0.5f *(powf(shiftx,2)/powf(gaussX_std, 2.0f)
                                                                     +powf(shifty,2)/powf(gaussY_std, 2.0f)
                                                                     +powf(shiftz,2)/powf(gaussZ_std, 2.0f)
                                                                     )))/(sqrtf(2.0f*3.14159265*powf(gauss_std, 2)));
                                int index2=index+(shiftx*shiftdirection[0])+(shifty*shiftdirection[1])+(shiftz*shiftdirection[2]);

                                if( !(DataPTR[index2]!=DataPTR[index2]) ){

                                    std::map<unsigned int,float>::iterator location=tmp_lab.find(DataPTR[index2]);

                                    if(location!=tmp_lab.end())
                                    {
                                        location->second=location->second+kernelval;
                                    }
                                    else
                                    {
                                        unsigned int tmpIndex = (unsigned int)round(DataPTR[index2]);
                                        DataPointPair tmpPair = std::make_pair(tmpIndex, kernelval);
                                        tmp_lab.insert(tmpPair);
                                    }
                                }
                            }
                        }
                    }
                    std::map<unsigned int,float>::iterator currIterator = tmp_lab.begin();
                    int maxindex=std::numeric_limits<float>::quiet_NaN();
                    float maxval=-FLT_MAX;
                    if(currIterator!=tmp_lab.end()){
                        while(currIterator != tmp_lab.end())
                        {
                            if(currIterator->second>maxval)
                            {
                                maxindex=currIterator->first;
                                maxval=currIterator->second;
                            }
                            currIterator++;
                        }
                        ImageBuffer[index]=maxindex;
                    }
                    else{
                        ImageBuffer[index]=std::numeric_limits<float>::quiet_NaN();
                    }
                }
            }
        }

        for(index=0; index<numel; index++)
        {
            DataPTR[index+current_4dShift_short]=ImageBuffer[index];
        }
    }

    delete [] ImageBuffer;
    return;
}


int quickSort(int *arr, int elements)
{

    int  piv, beg[300], end[300], i=0, L, R, swap ;

    beg[0]=0;
    end[0]=elements;
    while (i>=0)
    {
        L=beg[i];
        R=end[i]-1;
        if (L<R)
        {
            piv=arr[L];
            while (L<R)
            {
                while (arr[R]>=piv && L<R) R--;
                if (L<R) arr[L++]=arr[R];
                while (arr[L]<=piv && L<R) L++;
                if (L<R) arr[R--]=arr[L];
            }
            arr[L]=piv;
            beg[i+1]=L+1;
            end[i+1]=end[i];
            end[i++]=L;
            if (end[i]-beg[i]>end[i-1]-beg[i-1])
            {
                swap=beg[i];
                beg[i]=beg[i-1];
                beg[i-1]=swap;
                swap=end[i];
                end[i]=end[i-1];
                end[i-1]=swap;
            }
        }
        else
        {
            i--;
        }
    }
    return 1;
}


int quickSort(float *arr, int elements)
{

    float  piv;
    int beg[300], end[300], i=0, L, R, swap ;

    beg[0]=0;
    end[0]=elements;

    while (i>=0)
    {
        L=beg[i];
        R=end[i]-1;
        if (L<R)
        {
            piv=arr[L];
            while (L<R)
            {
                while (arr[R]>=piv && L<R) R--;
                if (L<R) arr[L++]=arr[R];
                while (arr[L]<=piv && L<R) L++;
                if (L<R) arr[R--]=arr[L];
            }
            arr[L]=piv;
            beg[i+1]=L+1;
            end[i+1]=end[i];
            end[i++]=L;
            if (end[i]-beg[i]>end[i-1]-beg[i-1])
            {
                swap=beg[i];
                beg[i]=beg[i-1];
                beg[i-1]=swap;
                swap=end[i];
                end[i]=end[i-1];
                end[i-1]=swap;
            }
        }
        else
        {
            i--;
        }
    }
    return 1;
}



int * quickSort_order(float *arr, int elements)
{

    float  piv;
    int piv_index, beg[300], end[300], i=0, L, R, swap ;
    int * order = new int [elements];
    for(long index=0; index<elements; index++)
    {
        order[index]=index;
    }
    beg[0]=0;
    end[0]=elements;
    while (i>=0)
    {
        L=beg[i];
        R=end[i]-1;
        if (L<R)
        {
            piv=arr[L];
            piv_index=order[L];
            while (L<R)
            {
                while (arr[R]>=piv && L<R)
                {
                    R--;
                }
                if (L<R)
                {
                    arr[L]=arr[R];
                    order[L++]=order[R];
                }
                while (arr[L]<=piv && L<R)
                {
                    L++;
                }
                if (L<R)
                {
                    arr[R]=arr[L];
                    order[R--]=order[L];
                }
            }
            arr[L]=piv;
            order[L]=piv_index;
            beg[i+1]=L+1;
            end[i+1]=end[i];
            end[i++]=L;
            if (end[i]-beg[i]>end[i-1]-beg[i-1])
            {
                swap=beg[i];
                beg[i]=beg[i-1];
                beg[i-1]=swap;
                swap=end[i];
                end[i]=end[i-1];
                end[i-1]=swap;
            }
        }
        else
        {
            i--;
        }
    }
    return order;
}




void MaxHeapify(float * a,int i,int n)
{
    int l,r,lr;
    float t;
    l=2*i+1;
    r=2*i+2;
    if((l<=n)&&(a[l]>a[i]))lr=l;
    else lr=i;
    if((r<=n)&&(a[r]>a[lr]))lr=r;
    if(lr!=i)
    {
        t=a[i];
        a[i]=a[lr];
        a[lr]=t;
        MaxHeapify(a,lr,n);
    }
}

void BuildMaxHeap(float * a,int n)
{
    int i;
    for(i=(n/2); i>=0; i--)
        MaxHeapify(a,i,n);
}

void HeapSort(float * a,int n)
{
    int i;
    float t;
    BuildMaxHeap(a,n);
    for(i=n; i>0; i--)
    {
        t=a[0];
        a[0]=a[i];
        a[i]=t;
        n--;
        MaxHeapify(a,0,n);
    }
}




unsigned char * estimateNCC4D(nifti_image * BaseImage,nifti_image * NCC,int numberordered,ImageSize * CurrSizes,int verbose)
{

    segPrecisionTYPE * NCCptr = static_cast<segPrecisionTYPE *>(NCC->data);
    segPrecisionTYPE * BaseImageptr = static_cast<segPrecisionTYPE *>(BaseImage->data);

    segPrecisionTYPE BaseMean=0.0f;
    segPrecisionTYPE BaseSTD=0.0f;


    segPrecisionTYPE bufferMean=0.0f;
    segPrecisionTYPE bufferSTD=0.0f;


    // CALC MEAN AND STD OF THE BASE
    if (verbose>0)
    {
        cout << "Calculating NCC"<<endl;
        cout << "Local Mean and STD of the base image"<<endl;
        flush(cout);
    }
    for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
    {
        BaseMean+=BaseImageptr[i];

    }
    BaseMean=BaseMean/(BaseImage->nx*BaseImage->ny*BaseImage->nz);
    for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
    {
        BaseSTD+=(BaseImageptr[i]-BaseMean)*(BaseImageptr[i]-BaseMean);
    }
    BaseSTD=sqrtf(BaseSTD/(BaseImage->nx*BaseImage->ny*BaseImage->nz));
    // Calc Mean
    int realt=CurrSizes->numclass;


    // CALC NCC FOR EACH
    if (verbose>0)
    {
        cout << "Local Mean and STD of the Template images"<<endl;
        flush(cout);
    }

    segPrecisionTYPE * NCCval= new segPrecisionTYPE [NCC->nt];
    for(long currlable=0; currlable<NCC->nt; currlable++)
    {
        NCCval[currlable]=0;
        bufferMean=0;
        bufferSTD=0;
        float * currNCCptr=&NCCptr[currlable*BaseImage->nx*BaseImage->ny*BaseImage->nz];

        for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
        {
            bufferMean+=currNCCptr[i];
        }
        bufferMean=bufferMean/(BaseImage->nx*BaseImage->ny*BaseImage->nz);
        for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
        {
            bufferSTD+=(currNCCptr[i]-bufferMean)*(currNCCptr[i]-bufferMean);
        }
        bufferSTD=sqrtf(bufferSTD/(BaseImage->nx*BaseImage->ny*BaseImage->nz));
        NCCval[currlable]=0;
        for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
        {
            NCCval[currlable]+=((currNCCptr[i]-bufferMean)*(BaseImageptr[i]-BaseMean));
        }
        NCCval[currlable]=NCCval[currlable]/((bufferSTD*BaseSTD)+0.0000001)/(BaseImage->nx*BaseImage->ny*BaseImage->nz);
        if (verbose>0)
        {
            cout << currlable+1 << "/" << NCC->nt<<" - GNCC= "<< NCCval[currlable]<<endl;
            flush(cout);
        }
    }

    CurrSizes->numclass=realt;


    unsigned char * NCC_ordered=new unsigned char [numberordered];

    if (verbose>0)
    {
        cout << "Used Labels after sorting the GNCC"<<endl;
    }

    int * ordertmp=quickSort_order(&NCCval[0],NCC->nt);
    for(long lable_order=0; lable_order<numberordered; lable_order++)
    {
        NCC_ordered[lable_order]=(char)ordertmp[NCC->nt-lable_order-1];
        if (verbose>0)
        {
            cout << (int)NCC_ordered[lable_order]+1 << "/" << numberordered<<endl;
        }
    }


    return NCC_ordered;
}


float estimateNCC3D(nifti_image * BaseImage,nifti_image * Template,nifti_image * Mask,int verbose)
{

    segPrecisionTYPE * Templateptr = static_cast<segPrecisionTYPE *>(Template->data);
    segPrecisionTYPE * BaseImageptr = static_cast<segPrecisionTYPE *>(BaseImage->data);

    bool * Maskptr=NULL;
    int Maskcount=0;
    bool usemask=false;
    if (Mask!=NULL)
    {
        Maskptr = static_cast<bool * >(Mask->data);
        usemask=true;
        for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
        {
            // Account for NANs
            if(Templateptr[i]!=Templateptr[i] || BaseImageptr[i]!=BaseImageptr[i]){
                Maskptr[i]=0;
            }
            Maskcount+=(float)(Maskptr[i]>0);
        }
    }

    segPrecisionTYPE BaseMean=0.0f;
    segPrecisionTYPE BaseSTD=0.0f;


    segPrecisionTYPE bufferMean=0.0f;
    segPrecisionTYPE bufferSTD=0.0f;


    // CALC MEAN AND STD OF THE BASE
    if (verbose>0)
    {
        cout << "Calculating NCC"<<endl;
    }


    for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
    {

        if(usemask)
        {
            BaseMean+=(Maskptr[i])?BaseImageptr[i]:0;
        }
        else
        {
            BaseMean+=BaseImageptr[i];
        }
    }

    if(usemask)
    {
        BaseMean=BaseMean/Maskcount;
    }
    else
    {
        BaseMean=BaseMean/(BaseImage->nx*BaseImage->ny*BaseImage->nz);
    }

    for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
    {
        if(usemask)
        {
            BaseSTD+=(Maskptr[i])?((BaseImageptr[i]-BaseMean)*(BaseImageptr[i]-BaseMean)):0;
        }
        else
        {
            BaseSTD+=(BaseImageptr[i]-BaseMean)*(BaseImageptr[i]-BaseMean);
        }
    }
    if(usemask)
    {
        BaseSTD=sqrtf(BaseSTD/Maskcount);
    }
    else
    {
        BaseSTD=sqrtf(BaseSTD/(BaseImage->nx*BaseImage->ny*BaseImage->nz));
    }

    if (verbose>0)
    {
        cout << "Local Mean and STD of the base image = "<<BaseMean << "( " << BaseSTD<< " )"<<endl;
        flush(cout);
    }

    // Calc Mean
    // CALC NCC FOR EACH


    segPrecisionTYPE NCCval=0;

    bufferMean=0;
    bufferSTD=0;

    for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
    {
        if(usemask)
        {
            bufferMean+=(Maskptr[i])?Templateptr[i]:0;
        }
        else
        {
            bufferMean+=Templateptr[i];
        }
    }
    if(usemask)
    {
        bufferMean=bufferMean/Maskcount;
    }
    else
    {
        bufferMean=bufferMean/(BaseImage->nx*BaseImage->ny*BaseImage->nz);
    }
    for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
    {
        if(usemask)
        {
            bufferSTD+=(Maskptr[i])?((Templateptr[i]-bufferMean)*(Templateptr[i]-bufferMean)):0;
        }
        else
        {
            bufferSTD+=(Templateptr[i]-bufferMean)*(Templateptr[i]-bufferMean);
        }
    }
    if(usemask)
    {
        bufferSTD=sqrtf(bufferSTD/Maskcount);
    }
    else
    {
        bufferSTD=sqrtf(bufferSTD/(BaseImage->nx*BaseImage->ny*BaseImage->nz));
    }

    if (verbose>0)
    {
        cout << "Local Mean and STD of the Template image = "<<bufferMean << "( " << bufferSTD<< " )"<<endl;
        flush(cout);
    }

    NCCval=0;
    for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
    {
        if(usemask)
        {
            NCCval+=(Maskptr[i])?((Templateptr[i]-bufferMean)*(BaseImageptr[i]-BaseMean)):0;
        }
        else
        {
            NCCval+=((Templateptr[i]-bufferMean)*(BaseImageptr[i]-BaseMean));
        }

    }

    if(usemask)
    {
        NCCval=NCCval/((bufferSTD*BaseSTD)+0.0000001)/(Maskcount);
    }
    else
    {
        NCCval=NCCval/((bufferSTD*BaseSTD)+0.0000001)/(BaseImage->nx*BaseImage->ny*BaseImage->nz);
    }

    return NCCval;
}





unsigned char * estimateROINCC4D(nifti_image * LableImage,nifti_image * BaseImage,nifti_image * NCC,int numberordered,ImageSize * CurrSizes,int DilSize, int verbose)
{

    segPrecisionTYPE * NCCptr = static_cast<segPrecisionTYPE *>(NCC->data);
    bool * LableImageptr = static_cast<bool *>(LableImage->data);
    segPrecisionTYPE * BaseImageptr = static_cast<segPrecisionTYPE *>(BaseImage->data);
    float * ROIarea=new float [BaseImage->nx*BaseImage->ny*BaseImage->nz];

    segPrecisionTYPE BaseMean=0.0f;
    segPrecisionTYPE BaseSTD=0.0f;
    segPrecisionTYPE bufferMean=0.0f;
    segPrecisionTYPE bufferSTD=0.0f;

    // CALC MEAN AND STD OF THE BASE
    if (verbose>0)
    {
        cout << "Calculating ROI from Lables"<<endl;
        flush(cout);
    }
    for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
    {
        ROIarea[i]=0;
    }



    int ROIsize=0;
    for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
    {
        bool ROItmp=false;
        for(long currlable=0; currlable<NCC->nt; currlable++)
        {
            if(LableImageptr[i+currlable*BaseImage->nx*BaseImage->ny*BaseImage->nz])
            {
                ROItmp=true;
            }
        }
        if(ROItmp==true)
        {
            ROIsize++;
            ROIarea[i]=1;
        }
    }

    Dillate(ROIarea,DilSize,CurrSizes);


    // CALC MEAN AND STD OF THE BASE
    if (verbose>0)
    {
        cout << "Calculating NCC"<<endl;
        cout << "Local Mean and STD of the base image"<<endl;
        flush(cout);
    }
    for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
    {
        if(ROIarea[i]>0)
        {
            BaseMean+=BaseImageptr[i];
        }
    }
    BaseMean=BaseMean/(ROIsize);
    for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
    {
        if(ROIarea[i]>0)
        {
            BaseSTD+=(BaseImageptr[i]-BaseMean)*(BaseImageptr[i]-BaseMean);
        }
    }
    BaseSTD=sqrtf(BaseSTD/(ROIsize));
    // Calc Mean
    int realt=CurrSizes->numclass;


    // CALC NCC FOR EACH
    if (verbose>0)
    {
        cout << "Local Mean and STD of the Template images"<<endl;
        flush(cout);
    }

    segPrecisionTYPE * NCCval= new segPrecisionTYPE [NCC->nt];
    for(long currlable=0; currlable<NCC->nt; currlable++)
    {
        NCCval[currlable]=0;
        float * currNCCptr=&NCCptr[currlable*BaseImage->nx*BaseImage->ny*BaseImage->nz];
        bufferMean=0;
        bufferSTD=0;
        for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
        {
            if(ROIarea[i]>0)
            {
                bufferMean+=currNCCptr[i];
            }
        }
        bufferMean=bufferMean/(ROIsize);
        for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
        {
            if(ROIarea[i]>0)
            {
                bufferSTD+=(currNCCptr[i]-bufferMean)*(currNCCptr[i]-bufferMean);
            }
        }
        bufferSTD=sqrtf(bufferSTD/(ROIsize));
        NCCval[currlable]=0;
        for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
        {
            if(ROIarea[i]>0)
            {
                NCCval[currlable]+=((currNCCptr[i]-bufferMean)*(BaseImageptr[i]-BaseMean));
            }
        }
        NCCval[currlable]=NCCval[currlable]/((bufferSTD*BaseSTD)+0.0000001)/(ROIsize);
        if (verbose>0)
        {
            cout << currlable+1 << "/" << NCC->nt<<" - ROINCC= "<< NCCval[currlable]<<endl;
            flush(cout);
        }
    }

    CurrSizes->numclass=realt;


    unsigned char * NCC_ordered=new unsigned char [numberordered];

    if (verbose>0)
    {
        cout << "Used Labels after sorting the ROINCC"<<endl;
    }

    int * ordertmp=quickSort_order(&NCCval[0],NCC->nt);
    for(long lable_order=0; lable_order<numberordered; lable_order++)
    {
        NCC_ordered[lable_order]=(char)ordertmp[NCC->nt-lable_order-1];
        if (verbose>0)
        {
            cout << (int)NCC_ordered[lable_order]+1 << "/" << numberordered<<endl;
        }
    }


    return NCC_ordered;
}


unsigned char * estimateMLNCC4D(nifti_image * BaseImage, nifti_image * LNCC,float distance,int levels, int numberordered,ImageSize * CurrSizes,int verbose)
{

    segPrecisionTYPE * LNCCptr = static_cast<segPrecisionTYPE *>(LNCC->data);
    segPrecisionTYPE * BaseImageptr = static_cast<segPrecisionTYPE *>(BaseImage->data);
    unsigned char * LNCC_ordered=NULL;
    unsigned char * LNCC_ordered_save=NULL;


    int numbordered_level_old=LNCC->nt;

    for(long curlevel=levels; curlevel; curlevel--)
    {

        if(curlevel==levels)
        {
            LNCC_ordered_save=new unsigned char [numbordered_level_old*BaseImage->nx*BaseImage->ny*BaseImage->nz];
            if(LNCC_ordered_save == NULL)
            {
                fprintf(stderr,"* Error when alocating LNCC_ordered_save in function seg_norm4MLLNCC");
                exit(-1);
            }

            for(long cl=0; cl<numbordered_level_old; cl++)
            {
                for(long cl_index=0; cl_index<BaseImage->nx*BaseImage->ny*BaseImage->nz; cl_index++)
                {
                    LNCC_ordered_save[cl_index+cl*(BaseImage->nx*BaseImage->ny*BaseImage->nz)]=cl;
                }
            }
        }


        float distance_level=distance*pow(2.0,(curlevel-1));
        int numbordered_level=(curlevel*numberordered)<=LNCC->nt?(curlevel*numberordered):LNCC->nt;


        segPrecisionTYPE * BaseMean=new segPrecisionTYPE [BaseImage->nx*BaseImage->ny*BaseImage->nz];
        if(BaseMean == NULL)
        {
            fprintf(stderr,"* Error when alocating BaseMean in function seg_norm4LNCC");
            exit(-1);
        }
        segPrecisionTYPE * BaseSTD=new segPrecisionTYPE [BaseImage->nx*BaseImage->ny*BaseImage->nz];
        if(BaseSTD == NULL)
        {
            fprintf(stderr,"* Error when alocating BaseSTD in function seg_norm4LNCC");
            exit(-1);
        }

        // CALC MEAN AND STD OF THE BASE
        if (verbose>0)
        {
            cout << "Calculating LNCC at level "<<curlevel<< " ( kernel size = "<<distance_level<<" , number templates = "<<numbordered_level<<" )"<<endl;
            cout << "Local Mean and STD of the base image"<<endl;
            flush(cout);
        }
        for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
        {
            BaseMean[i]=BaseImageptr[i];
            BaseSTD[i]=BaseImageptr[i]*BaseImageptr[i];
        }
        // Calc Mean
        int realt=CurrSizes->tsize;
        CurrSizes->tsize=1;
        GaussianFilter4D_cArray(BaseMean,(float)(distance_level),CurrSizes);
        GaussianFilter4D_cArray(BaseSTD,(float)(distance_level),CurrSizes);

        for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
        {
            BaseSTD[i]=BaseSTD[i]-BaseMean[i]*BaseMean[i];
        }

        // CALC LNCC FOR EACH
        if (verbose>0)
        {
            cout << "Local Mean and STD of the Template images"<<endl;
            flush(cout);
        }
        int currlable=0;

        //#ifdef _OPENMP
        //#pragma omp parallel for shared(BaseImageptr, BaseSTD, BaseMean, LNCC,BaseImage,stderr, verbose,cout,LNCCptr,CurrSizes,distance_level)
        //#endif

        for(currlable=0; currlable<LNCC->nt; currlable++)
        {
            segPrecisionTYPE * bufferMean=new segPrecisionTYPE [BaseImage->nx*BaseImage->ny*BaseImage->nz];
            if(bufferMean == NULL)
            {
                fprintf(stderr,"* Error when alocating bufferMean in function seg_norm4LNCC");
                exit(-1);
            }
            segPrecisionTYPE * bufferSTD=new segPrecisionTYPE [ BaseImage->nx * BaseImage->ny * BaseImage->nz ];
            if(bufferSTD == NULL)
            {
                fprintf(stderr,"* Error when alocating bufferSTD in function seg_norm4LNCC");
                exit(-1);
            }
            segPrecisionTYPE * bufferDATA=new segPrecisionTYPE [BaseImage->nx*BaseImage->ny*BaseImage->nz];
            if(bufferDATA == NULL)
            {
                fprintf(stderr,"* Error when alocating bufferDATA in function seg_norm4LNCC");
                exit(-1);
            }

            if (verbose>0)
            {
                cout << currlable+1 << "/" << LNCC->nt<<"\n";
                flush(cout);
            }
            segPrecisionTYPE * currLNCCptr=&LNCCptr[currlable*BaseImage->nx*BaseImage->ny*BaseImage->nz];
            for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
            {
                bufferDATA[i]=currLNCCptr[i]*BaseImageptr[i];
                bufferMean[i]=currLNCCptr[i];
                bufferSTD[i]=currLNCCptr[i]*currLNCCptr[i];
            }

            // Calc Mean
            GaussianFilter4D_cArray(bufferMean,(float)(distance_level),CurrSizes);
            GaussianFilter4D_cArray(bufferSTD,(float)(distance_level),CurrSizes);
            GaussianFilter4D_cArray(bufferDATA,(float)(distance_level),CurrSizes);
            for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
            {
                bufferSTD[i]=bufferSTD[i]-bufferMean[i]*bufferMean[i];
                currLNCCptr[i]=(bufferDATA[i]-BaseMean[i]*bufferMean[i])/(sqrt(bufferSTD[i]*BaseSTD[i])+0.0000001);
                currLNCCptr[i]=currLNCCptr[i]>0?currLNCCptr[i]:0;
            }
            delete [] bufferSTD;
            delete [] bufferMean;
            delete [] bufferDATA;
        }

        delete [] BaseSTD;
        delete [] BaseMean;

        CurrSizes->tsize=realt;
        LNCC_ordered=new unsigned char [numbordered_level*BaseImage->nx*BaseImage->ny*BaseImage->nz];
        if(LNCC_ordered == NULL)
        {
            fprintf(stderr,"* Error when alocating LNCC_ordered in function seg_norm4LNCC");
            exit(-1);
        }


        if (verbose>0)
        {
            cout << "Sorting"<<endl;
        }

#ifdef _OPENMP
#pragma omp parallel for default(none) \
    shared(LNCC,LNCCptr,CurrSizes,numbordered_level,LNCC_ordered,numbordered_level_old,LNCC_ordered_save)
#endif
        for(long i=0; i<LNCC->nx*LNCC->ny*LNCC->nz; i++)
        {
            segPrecisionTYPE LNCCvalue_tmp[1000];
            char old_sort[1000];;
            for(long currlable=0; currlable<LNCC->nt; currlable++)
            {
                LNCCvalue_tmp[currlable]=LNCCptr[i+currlable*LNCC->nx*LNCC->ny*LNCC->nz];
            }

            for(long currlable=0; currlable<numbordered_level_old; currlable++)
            {
                old_sort[currlable]=LNCC_ordered_save[i+currlable*LNCC->nx*LNCC->ny*LNCC->nz];
            }

            int * ordertmp=quickSort_order(&LNCCvalue_tmp[0],LNCC->nt);

            int label_order_index=0;
            for(long lable_order=0; lable_order<LNCC->nt; lable_order++)
            {
                char cur_LNCC_ordered_label=(char)ordertmp[LNCC->nt-lable_order-1];
                for(long currlable=0; currlable<numbordered_level_old; currlable++)
                {
                    if(cur_LNCC_ordered_label==old_sort[currlable])
                    {
                        LNCC_ordered[i+label_order_index*LNCC->nx*LNCC->ny*LNCC->nz]=cur_LNCC_ordered_label;
                        label_order_index++;
                        currlable=numbordered_level_old;
                    }
                }
                if((label_order_index)>=numbordered_level)
                {
                    lable_order=LNCC->nt;
                }

            }
            delete [] ordertmp;

        }

        if(curlevel>1)
        {
            numbordered_level_old=numbordered_level;
            delete [] LNCC_ordered_save;
            LNCC_ordered_save=new unsigned char [numbordered_level_old*LNCC->nx*LNCC->ny*LNCC->nz];
            if(LNCC_ordered_save == NULL)
            {
                fprintf(stderr,"* Error when alocating LNCC_ordered_save in function seg_norm4MLLNCC");
                exit(-1);
            }
            for(long cl=0; cl<numbordered_level_old*LNCC->nx*LNCC->ny*LNCC->nz; cl++)
            {
                LNCC_ordered_save[cl]=LNCC_ordered[cl];
            }
            delete [] LNCC_ordered;
        }

        if (verbose>0)
        {
            cout << "Finished sorting"<< endl;
            flush(cout);
        }
    }

    return LNCC_ordered;
}

/* *************************************************************** */
/* *************************************************************** */
template<class PrecisionTYPE>
PrecisionTYPE GetBasisSplineValue(PrecisionTYPE x)
{
    x=fabs(x);
    PrecisionTYPE value=0.0;
    if(x<2.0)
    {
        if(x<1.0)
            value = (PrecisionTYPE)(2.0f/3.0f + (0.5f*x-1.0)*x*x);
        else
        {
            x-=2.0f;
            value = -x*x*x/6.0f;
        }
    }
    return value;
}
/* *************************************************************** */
/* *************************************************************** */
float seg_getNMIValue(nifti_image *img1,
                      nifti_image *img2,
                      unsigned char *referenceMask
                      )
{
    // Create pointers to the image data arrays
    float *img1Ptr = static_cast<float *>(img1->data);
    float *img2Ptr = static_cast<float *>(img2->data);
    // Useful variable
    size_t voxelNumber = (size_t)img1->nx *
            img1->ny *
            img1->nz;

    unsigned short referenceBinNumber=68;
    unsigned short floatingBinNumber=68;
    unsigned short totalBinNumber=floatingBinNumber*referenceBinNumber+referenceBinNumber+floatingBinNumber;
    float entropyValues[4]={0};

    // Empty the joint histogram
    double * jointHistogramPro=new double [totalBinNumber];
    double * jointHistogramProPtr=(double *)jointHistogramPro;
    memset((double *)jointHistogramProPtr,0,totalBinNumber*sizeof(double));
    //cout<<jointHistogramPro<<endl;

    double * jointHistogramLog=new double [totalBinNumber];
    double * jointHistogramLogPtr=(double *)jointHistogramLog;
    memset((double *)jointHistogramLogPtr,0,totalBinNumber*sizeof(double));

    // Fill the joint histograms using an approximation
    float *refPtr = &img1Ptr[0];
    float *warPtr = &img2Ptr[0];
    for(size_t voxel=0; voxel<voxelNumber; ++voxel)
    {
        if(referenceMask[voxel]>0)
        {
            float refValue=refPtr[voxel];
            float warValue=warPtr[voxel];
            if(refValue==refValue && warValue==warValue &&
                    refValue>=0 && warValue>=0 &&
                    refValue<referenceBinNumber &&
                    warValue<floatingBinNumber)
            {
                ++jointHistogramProPtr[static_cast<int>(refValue) +
                        static_cast<int>(warValue) * referenceBinNumber];
            }
        }
    }
    // Convolve the histogram with a cubic B-spline kernel
    double kernel[3];
    kernel[0]=kernel[2]=GetBasisSplineValue(-1.);
    kernel[1]=GetBasisSplineValue(0.);
    // Histogram is first smooth along the reference axis
    memset((double*)jointHistogramLogPtr,0,totalBinNumber*sizeof(double));
    for(int f=0; f<floatingBinNumber; ++f)
    {
        for(int r=0; r<referenceBinNumber; ++r)
        {
            double value=0.0;
            int index = r-1;
            double *ptrHisto = &jointHistogramProPtr[index+referenceBinNumber*f];

            for(int it=0; it<3; it++)
            {
                if(-1<index && index<referenceBinNumber)
                {
                    value += *ptrHisto * kernel[it];
                }
                ++ptrHisto;
                ++index;
            }
            jointHistogramLogPtr[r+referenceBinNumber*f] = value;
        }
    }
    // Histogram is then smooth along the warped floating axis
    for(int r=0; r<referenceBinNumber; ++r)
    {
        for(int f=0; f<floatingBinNumber; ++f)
        {
            double value=0.;
            int index = f-1;
            double *ptrHisto = &jointHistogramLogPtr[r+referenceBinNumber*index];

            for(int it=0; it<3; it++)
            {
                if(-1<index && index<floatingBinNumber)
                {
                    value += *ptrHisto * kernel[it];
                }
                ptrHisto+=referenceBinNumber;
                ++index;
            }
            jointHistogramProPtr[r+referenceBinNumber*f] = value;
        }
    }
    // Normalise the histogram
    double activeVoxel=0.f;
    for(int i=0; i<totalBinNumber; ++i)
        activeVoxel+=jointHistogramProPtr[i];
    entropyValues[3]=activeVoxel;
    for(int i=0; i<totalBinNumber; ++i)
        jointHistogramProPtr[i]/=activeVoxel;
    // Marginalise over the reference axis
    for(int r=0; r<referenceBinNumber; ++r)
    {
        double sum=0.;
        int index=r;
        for(int f=0; f<floatingBinNumber; ++f)
        {
            sum+=jointHistogramProPtr[index];
            index+=referenceBinNumber;
        }
        jointHistogramProPtr[referenceBinNumber*
                floatingBinNumber+r]=sum;
    }
    // Marginalise over the warped floating axis
    for(int f=0; f<floatingBinNumber; ++f)
    {
        double sum=0.;
        int index=referenceBinNumber*f;
        for(int r=0; r<referenceBinNumber; ++r)
        {
            sum+=jointHistogramProPtr[index];
            ++index;
        }
        jointHistogramProPtr[referenceBinNumber*
                floatingBinNumber+referenceBinNumber+f]=sum;
    }
    // Set the log values to zero
    memset(jointHistogramLogPtr,0,totalBinNumber*sizeof(double));
    // Compute the entropy of the reference image
    double referenceEntropy=0.;
    for(int r=0; r<referenceBinNumber; ++r)
    {
        double valPro=jointHistogramProPtr[referenceBinNumber*floatingBinNumber+r];
        if(valPro>0)
        {
            double valLog=log(valPro);
            referenceEntropy -= valPro * valLog;
            jointHistogramLogPtr[referenceBinNumber*floatingBinNumber+r]=valLog;
        }
    }
    entropyValues[0]=referenceEntropy;
    // Compute the entropy of the warped floating image
    double warpedEntropy=0.;
    for(int f=0; f<floatingBinNumber; ++f)
    {
        double valPro=jointHistogramProPtr[referenceBinNumber*floatingBinNumber+
                referenceBinNumber+f];
        if(valPro>0)
        {
            double valLog=log(valPro);
            warpedEntropy -= valPro * valLog;
            jointHistogramLogPtr[referenceBinNumber*floatingBinNumber+
                    referenceBinNumber+f]=valLog;
        }
    }
    entropyValues[1]=warpedEntropy;
    // Compute the joint entropy
    double jointEntropy=0.;
    for(int i=0; i<referenceBinNumber*floatingBinNumber; ++i)
    {
        double valPro=jointHistogramProPtr[i];
        if(valPro>0)
        {
            double valLog=log(valPro);
            jointEntropy -= valPro * valLog;
            jointHistogramLogPtr[i]=valLog;
        }
    }
    entropyValues[2]=jointEntropy;

    delete [] jointHistogramLog;
    //cout<<jointHistogramPro<<endl;
    delete [] jointHistogramPro;


    return (entropyValues[0]+entropyValues[1])/entropyValues[2];
}

/* *************************************************************** */

unsigned char * estimateLNCC5D(nifti_image * BaseImage, nifti_image * LNCC,float distance,int numberordered,ImageSize * CurrSizes,int verbose)
{

    segPrecisionTYPE * LNCCptr = static_cast<segPrecisionTYPE *>(LNCC->data);
    segPrecisionTYPE * BaseImageptr = static_cast<segPrecisionTYPE *>(BaseImage->data);
    unsigned char * LNCC_ordered=NULL;
    long BaseImageXYZsize=BaseImage->nx*BaseImage->ny*BaseImage->nz;
    segPrecisionTYPE * BaseMean=new segPrecisionTYPE [BaseImageXYZsize*BaseImage->nt];
    if(BaseMean == NULL)
    {
        fprintf(stderr,"* Error when alocating BaseMean in function seg_norm4LNCC");
        exit(-1);
    }
    segPrecisionTYPE * BaseSTD=new segPrecisionTYPE [BaseImageXYZsize*BaseImage->nt];
    if(BaseSTD == NULL)
    {
        fprintf(stderr,"* Error when alocating BaseSTD in function seg_norm4LNCC");
        exit(-1);
    }

    // CALC MEAN AND STD OF THE BASE
    if (verbose>0)
    {
        cout << "Calculating LNCC"<<endl;
        cout << "Local Mean and STD of the base image"<<endl;
        flush(cout);
    }
    for(long i=0; i<BaseImageXYZsize*BaseImage->nt; i++)
    {
        BaseMean[i]=BaseImageptr[i];
        BaseSTD[i]=BaseImageptr[i]*BaseImageptr[i];
    }
    // Calc Mean
    int realt=CurrSizes->tsize;
    CurrSizes->tsize=1;
    for(long currmodality=0; currmodality<BaseImage->nt; currmodality++)
    {
        if (verbose>0)
        {
            cout << "Modality "<<currmodality+1<<"/"<<BaseImage->nt<<endl;
            flush(cout);
        }
        GaussianFilter4D_cArray(&BaseMean[currmodality*BaseImageXYZsize],(float)(distance),CurrSizes);
        GaussianFilter4D_cArray(&BaseSTD[currmodality*BaseImageXYZsize],(float)(distance),CurrSizes);
    }

    for(long i=0; i<BaseImageXYZsize*BaseImage->nt; i++)
    {
        BaseSTD[i]=BaseSTD[i]-BaseMean[i]*BaseMean[i];
    }

    // CALC LNCC FOR EACH
    if (verbose>0)
    {
        cout << "Local Mean and STD of the Template images"<<endl;
        flush(cout);
    }

#ifdef _OPENMP
#pragma omp parallel for shared(BaseImageptr, BaseSTD, BaseMean, LNCC, BaseImage, verbose, cout, LNCCptr, CurrSizes, distance)
#endif

    for(long currlable=0; currlable<CurrSizes->numclass; currlable++)
    {
        for(long currmodality=0; currmodality<CurrSizes->nummod; currmodality++)
        {
            segPrecisionTYPE * bufferMean=new segPrecisionTYPE [BaseImageXYZsize];
            if(bufferMean == NULL)
            {
                fprintf(stderr,"* Error when alocating bufferMean in function seg_norm4LNCC");
                exit(-1);
            }
            segPrecisionTYPE * bufferSTD=new segPrecisionTYPE [BaseImageXYZsize];
            if(bufferSTD == NULL)
            {
                fprintf(stderr,"* Error when alocating bufferSTD in function seg_norm4LNCC");
                exit(-1);
            }
            segPrecisionTYPE * bufferDATA=new segPrecisionTYPE [BaseImageXYZsize];
            if(bufferDATA == NULL)
            {
                fprintf(stderr,"* Error when alocating bufferDATA in function seg_norm4LNCC");
                exit(-1);
            }
            //for(long currlable=0;currlable<3; currlable++){
            if (verbose>0)
            {
                std::stringstream stream;
                stream << "Template "<<currlable+1 << "/" << LNCC->nt<<"  - Modality "<<currmodality+1<<"/"<<CurrSizes->nummod<<"\n";
                std::cout << stream.str();
            }
            segPrecisionTYPE * currLNCCptr = &LNCCptr[currlable*BaseImageXYZsize + currmodality*BaseImageXYZsize*CurrSizes->numclass];
            for(long i=0; i<BaseImageXYZsize; i++)
            {
                bufferDATA[i]=currLNCCptr[i]*BaseImageptr[i+currmodality*BaseImageXYZsize];
                bufferMean[i]=currLNCCptr[i];
                bufferSTD[i]=currLNCCptr[i]*currLNCCptr[i];
            }

            // Calc Mean
            GaussianFilter4D_cArray(bufferMean,(float)(distance),CurrSizes);
            GaussianFilter4D_cArray(bufferSTD,(float)(distance),CurrSizes);
            GaussianFilter4D_cArray(bufferDATA,(float)(distance),CurrSizes);

            currLNCCptr=&LNCCptr[currlable*BaseImageXYZsize];
            for(long i=0; i<BaseImageXYZsize; i++)
            {
                bufferSTD[i]=bufferSTD[i]-bufferMean[i]*bufferMean[i];
                float tmpLNCC=(bufferDATA[i]-BaseMean[i+currmodality*BaseImageXYZsize]*bufferMean[i])/(sqrt(bufferSTD[i]*BaseSTD[i+currmodality*BaseImageXYZsize])+0.000001);
                if(currmodality==0){
                    currLNCCptr[i]=(tmpLNCC>0?tmpLNCC:0);
                }
                else{
                    LNCCptr[i+currlable*BaseImageXYZsize]+=(tmpLNCC>0?tmpLNCC:0);
                }

            }
            delete [] bufferSTD;
            delete [] bufferMean;
            delete [] bufferDATA;
        }
    }
    delete [] BaseSTD;
    delete [] BaseMean;

    CurrSizes->tsize=realt;
    // cout << "Filtering LNCC"<< endl;
    //  Gaussian_Filter_4D(LNCCptr,(float)(distance),CurrSizes);

    LNCC_ordered=new unsigned char [numberordered*BaseImage->nx*BaseImage->ny*BaseImage->nz];
    if(LNCC_ordered == NULL)
    {
        fprintf(stderr,"* Error when alocating LNCC_ordered in function seg_norm4LNCC");
        exit(-1);
    }


    if (verbose>0)
    {
        cout << "Sorting LNCC"<<endl;
    }

#ifdef _OPENMP
#pragma omp parallel for default(none) \
    shared(LNCC,BaseImage,LNCCptr,CurrSizes,numberordered,LNCC_ordered)
#endif

    for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
    {
        segPrecisionTYPE * LNCCvalue_tmp = new segPrecisionTYPE [LNCC->nt];
        for(long currlable=0; currlable<LNCC->nt; currlable++)
        {
            LNCCvalue_tmp[currlable]=LNCCptr[i+currlable*BaseImage->nx*BaseImage->ny*BaseImage->nz];
        }

        int * ordertmp=quickSort_order(&LNCCvalue_tmp[0],LNCC->nt);

        for(long lable_order=0; lable_order<numberordered; lable_order++)
        {
            LNCC_ordered[i+lable_order*BaseImage->nx*BaseImage->ny*BaseImage->nz]=(unsigned char)ordertmp[LNCC->nt-lable_order-1];
        }

        delete [] ordertmp;
        delete [] LNCCvalue_tmp;
    }
    if (verbose>0)
    {
        cout << "Finished sorting LNCC"<< endl;
        flush(cout);
    }

    return LNCC_ordered;
}
/* *************************************************************** */

unsigned char * estimateLNCC4D(nifti_image * BaseImage, nifti_image * LNCC,float distance,int numberordered,ImageSize * CurrSizes,int verbose)
{

    segPrecisionTYPE * LNCCptr = static_cast<segPrecisionTYPE *>(LNCC->data);
    segPrecisionTYPE * BaseImageptr = static_cast<segPrecisionTYPE *>(BaseImage->data);
    unsigned char * LNCC_ordered=NULL;
    segPrecisionTYPE * BaseMean=new segPrecisionTYPE [BaseImage->nx*BaseImage->ny*BaseImage->nz];
    if(BaseMean == NULL)
    {
        fprintf(stderr,"* Error when alocating BaseMean in function seg_norm4LNCC");
        exit(-1);
    }
    segPrecisionTYPE * BaseSTD=new segPrecisionTYPE [BaseImage->nx*BaseImage->ny*BaseImage->nz];
    if(BaseSTD == NULL)
    {
        fprintf(stderr,"* Error when alocating BaseSTD in function seg_norm4LNCC");
        exit(-1);
    }

    // CALC MEAN AND STD OF THE BASE
    if (verbose>0)
    {
        cout << "Calculating LNCC"<<endl;
        cout << "Local Mean and STD of the base image"<<endl;
        flush(cout);
    }
    for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
    {
        BaseMean[i]=BaseImageptr[i];
        BaseSTD[i]=BaseImageptr[i]*BaseImageptr[i];
    }
    // Calc Mean
    int realt=CurrSizes->tsize;
    CurrSizes->tsize=1;
    GaussianFilter4D_cArray(BaseMean,(float)(distance),CurrSizes);
    GaussianFilter4D_cArray(BaseSTD,(float)(distance),CurrSizes);

    for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
    {
        BaseSTD[i]=BaseSTD[i]-BaseMean[i]*BaseMean[i];
    }

    // CALC LNCC FOR EACH
    if (verbose>0)
    {
        cout << "Local Mean and STD of the Template images"<<endl;
        flush(cout);
    }

#ifdef _OPENMP
#pragma omp parallel for shared(BaseImageptr, BaseSTD, BaseMean, LNCC, BaseImage, verbose, cout, LNCCptr, CurrSizes, distance)
#endif

    for(long currlable=0; currlable<LNCC->nt; currlable++)
    {
        segPrecisionTYPE * bufferMean=new segPrecisionTYPE [BaseImage->nx*BaseImage->ny*BaseImage->nz];
        if(bufferMean == NULL)
        {
            fprintf(stderr,"* Error when alocating bufferMean in function seg_norm4LNCC");
            exit(-1);
        }
        segPrecisionTYPE * bufferSTD=new segPrecisionTYPE [ BaseImage->nx * BaseImage->ny * BaseImage->nz ];
        if(bufferSTD == NULL)
        {
            fprintf(stderr,"* Error when alocating bufferSTD in function seg_norm4LNCC");
            exit(-1);
        }
        segPrecisionTYPE * bufferDATA=new segPrecisionTYPE [BaseImage->nx*BaseImage->ny*BaseImage->nz];
        if(bufferDATA == NULL)
        {
            fprintf(stderr,"* Error when alocating bufferDATA in function seg_norm4LNCC");
            exit(-1);
        }
        //for(long currlable=0;currlable<3; currlable++){
        if (verbose>0)
        {
            std::stringstream stream;
            stream << "Template "<<currlable+1 << "/" << LNCC->nt<<"\n";
            std::cout << stream.str();
        }
        segPrecisionTYPE * currLNCCptr=&LNCCptr[currlable*BaseImage->nx*BaseImage->ny*BaseImage->nz];
        for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
        {
            bufferDATA[i]=currLNCCptr[i]*BaseImageptr[i];
            bufferMean[i]=currLNCCptr[i];
            bufferSTD[i]=currLNCCptr[i]*currLNCCptr[i];
        }

        // Calc Mean
        GaussianFilter4D_cArray(bufferMean,(float)(distance),CurrSizes);
        GaussianFilter4D_cArray(bufferSTD,(float)(distance),CurrSizes);
        GaussianFilter4D_cArray(bufferDATA,(float)(distance),CurrSizes);


        for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
        {
            bufferSTD[i]=bufferSTD[i]-bufferMean[i]*bufferMean[i];
            currLNCCptr[i]=(bufferDATA[i]-BaseMean[i]*bufferMean[i])/(sqrt(bufferSTD[i]*BaseSTD[i])+0.000001);
            currLNCCptr[i]=currLNCCptr[i]>0?currLNCCptr[i]:0;

        }


        delete [] bufferSTD;
        delete [] bufferMean;
        delete [] bufferDATA;
    }
    delete [] BaseSTD;
    delete [] BaseMean;


    // cout << "Filtering LNCC"<< endl;
    CurrSizes->tsize=realt;
    //  Gaussian_Filter_4D(LNCCptr,(float)(distance),CurrSizes);

    LNCC_ordered=new unsigned char [numberordered*BaseImage->nx*BaseImage->ny*BaseImage->nz];
    if(LNCC_ordered == NULL)
    {
        fprintf(stderr,"* Error when alocating LNCC_ordered in function seg_norm4LNCC");
        exit(-1);
    }


    if (verbose>0)
    {
        cout << "Sorting LNCC"<<endl;
    }

#ifdef _OPENMP
#pragma omp parallel for default(none) \
    shared(LNCC,BaseImage,LNCCptr,CurrSizes,numberordered,LNCC_ordered)
#endif

    for(long i=0; i<BaseImage->nx*BaseImage->ny*BaseImage->nz; i++)
    {
        segPrecisionTYPE * LNCCvalue_tmp = new segPrecisionTYPE [LNCC->nt];
        for(long currlable=0; currlable<LNCC->nt; currlable++)
        {
            LNCCvalue_tmp[currlable]=LNCCptr[i+currlable*BaseImage->nx*BaseImage->ny*BaseImage->nz];
        }

        int * ordertmp=quickSort_order(&LNCCvalue_tmp[0],LNCC->nt);

        for(long lable_order=0; lable_order<numberordered; lable_order++)
        {
            LNCC_ordered[i+lable_order*BaseImage->nx*BaseImage->ny*BaseImage->nz]=(unsigned char)ordertmp[LNCC->nt-lable_order-1];
        }

        delete [] ordertmp;
        delete [] LNCCvalue_tmp;
    }
    if (verbose>0)
    {
        cout << "Finished sorting LNCC"<< endl;
        flush(cout);
    }

    return LNCC_ordered;
}

/* *************************************************************** */

int seg_convert2binary(nifti_image *image,
                       float thresh)
{
    switch(image->datatype)
    {
    case DT_BINARY:
        break;
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
        seg_convert2binary_data<segPrecisionTYPE>(image,thresh);
        break;
    default:
        printf("err\tseg_convert2binary\tThe initial image data type (%d) is not supported\n",image->datatype);
        return 1;
    }
    return 1;
}


template <class DTYPE>
int seg_convert2binary_data(nifti_image *image,
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

    for(unsigned int i=0; i<image->nvox; i++)
    {
        dataPtr[i] = (bool)((initialValue[i])>thresh);
    }
    free(initialValue);
    return 1;
}


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
    else
    {
        fprintf(stderr,"[NiftyReg ERROR] reg_tools_changeDatatype\tOnly change to unsigned char, float or double are supported\n");
        exit(1);
    }
    free(image->data);
    image->nbyper = sizeof(NewTYPE);
    image->data = (void *)calloc(image->nvox,sizeof(NewTYPE));
    NewTYPE *dataPtr = static_cast<NewTYPE *>(image->data);
    for(size_t i=0; i<image->nvox; i++)
        dataPtr[i] = (NewTYPE)(initialValue[i]);

    free(initialValue);
    return;
}
/* *************************************************************** */
template <class NewTYPE>
int seg_changeDatatype(nifti_image *image)
{
    switch(image->datatype)
    {
    case DT_BINARY:
        seg_changeDatatype1<NewTYPE,bool>(image);
        break;
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
        fprintf(stderr,"[NiftyReg ERROR] seg_changeDatatype\tThe initial image data type (%d) is not supported\n",image->datatype);
        exit(1);
    }
    return 1;
}
/* *************************************************************** */
template NIFTYSEG_WINEXPORT int seg_changeDatatype<unsigned char>(nifti_image *);
template NIFTYSEG_WINEXPORT int seg_changeDatatype<float>(nifti_image *);
template NIFTYSEG_WINEXPORT int seg_changeDatatype<double>(nifti_image *);

/* *************************************************************** */
template<class SourceTYPE, class FieldTYPE>
void TrilinearResampleSourceImage_for_GIF(  nifti_image *sourceImage,
                                            nifti_image *deformationField,
                                            nifti_image *resultImage,
                                            int *mask,
                                            FieldTYPE bgValue)
{
    // The resampling scheme is applied along each time
    SourceTYPE *sourceIntensityPtr = static_cast<SourceTYPE *>(sourceImage->data);
    SourceTYPE *resultIntensityPtr = static_cast<SourceTYPE *>(resultImage->data);
    FieldTYPE *deformationFieldPtrX = static_cast<FieldTYPE *>(deformationField->data);
    int targetVoxelNumber = resultImage->nx*resultImage->ny*resultImage->nz;
    int sourceVoxelNumber = sourceImage->nx*sourceImage->ny*sourceImage->nz;
    FieldTYPE *deformationFieldPtrY = &deformationFieldPtrX[targetVoxelNumber];
    FieldTYPE *deformationFieldPtrZ = &deformationFieldPtrY[targetVoxelNumber];

    int *maskPtr = &mask[0];
    mat44 *sourceIJKMatrix;
    if(sourceImage->sform_code>0)
        sourceIJKMatrix=&(sourceImage->sto_ijk);
    else sourceIJKMatrix=&(sourceImage->qto_ijk);

    for(long t=0; t<resultImage->nt; t++)
    {
#ifndef NDEBUG
        printf("[NiftyReg DEBUG] 3D linear resampling of volume number %i\n",t);
#endif

        SourceTYPE *resultIntensity = &resultIntensityPtr[t*targetVoxelNumber];
        SourceTYPE *sourceIntensity = &sourceIntensityPtr[t*sourceVoxelNumber];

        FieldTYPE xBasis[2], yBasis[2], zBasis[2], relative;
        int a, b, c, Y, Z, previous[3], index;
        SourceTYPE *zPointer, *xyzPointer;
        FieldTYPE xTempNewValue, yTempNewValue, xTempNewBasis, yTempNewBasis, intensity, basis, world[3], position[3];
#ifdef _OPENMP
#pragma omp parallel for default(none) \
    private(index, intensity, world, position, previous, xBasis, yBasis, zBasis, relative, \
    a, b, c, Y, Z, zPointer, xyzPointer, xTempNewValue, yTempNewValue) \
    shared(sourceIntensity, resultIntensity, targetVoxelNumber, sourceVoxelNumber, \
    deformationFieldPtrX, deformationFieldPtrY, deformationFieldPtrZ, maskPtr, \
    sourceIJKMatrix, sourceImage, bgValue)
#endif // _OPENMP
        for(index=0; index<targetVoxelNumber; index++)
        {

            intensity=0.0;
            basis=0.0;

            if(maskPtr[index]>-1)
            {

                world[0]=(FieldTYPE) deformationFieldPtrX[index];
                world[1]=(FieldTYPE) deformationFieldPtrY[index];
                world[2]=(FieldTYPE) deformationFieldPtrZ[index];

                /* real -> voxel; source space */
                reg_mat44_mul(sourceIJKMatrix, world, position);

                if( position[0]>=0.f && position[0]<(FieldTYPE)(sourceImage->nx-1) &&
                        position[1]>=0.f && position[1]<(FieldTYPE)(sourceImage->ny-1) &&
                        position[2]>=0.f && position[2]<(FieldTYPE)(sourceImage->nz-1) )
                {

                    previous[0] = (int)position[0];
                    previous[1] = (int)position[1];
                    previous[2] = (int)position[2];
                    // basis values along the x axis
                    relative=position[0]-(FieldTYPE)previous[0];
                    if(relative<0) relative=0.0; // rounding error correction
                    xBasis[0]= (FieldTYPE)(1.0-relative);
                    xBasis[1]= relative;
                    // basis values along the y axis
                    relative=position[1]-(FieldTYPE)previous[1];
                    if(relative<0) relative=0.0; // rounding error correction
                    yBasis[0]= (FieldTYPE)(1.0-relative);
                    yBasis[1]= relative;
                    // basis values along the z axis
                    relative=position[2]-(FieldTYPE)previous[2];
                    if(relative<0) relative=0.0; // rounding error correction
                    zBasis[0]= (FieldTYPE)(1.0-relative);
                    zBasis[1]= relative;

                    for(c=0; c<2; c++)
                    {
                        Z= previous[2]+c;
                        zPointer = &sourceIntensity[Z*sourceImage->nx*sourceImage->ny];
                        yTempNewValue=0.0;
                        yTempNewBasis=0.0;
                        for(b=0; b<2; b++)
                        {
                            Y= previous[1]+b;
                            xyzPointer = &zPointer[Y*sourceImage->nx+previous[0]];
                            xTempNewValue=0.0;
                            xTempNewBasis=0.0;
                            for(a=0; a<2; a++)
                            {
                                if((FieldTYPE)*xyzPointer>=0)
                                {
                                    xTempNewValue +=  (FieldTYPE)*xyzPointer * xBasis[a];
                                    xTempNewBasis += xBasis[a];
                                }
                                xyzPointer++;
                            }
                            yTempNewValue += (xTempNewValue * yBasis[b]);
                            yTempNewBasis += xTempNewBasis;
                        }
                        intensity += yTempNewValue * zBasis[c];
                        basis += yTempNewBasis;
                    }

                    if(basis<0.99 && basis>0)
                        intensity/=basis;
                }
                else intensity = -1.0f;
            }

            switch(sourceImage->datatype)
            {
            case NIFTI_TYPE_FLOAT32:
                resultIntensity[index]=(SourceTYPE)intensity;
                break;
            case NIFTI_TYPE_FLOAT64:
                resultIntensity[index]=(SourceTYPE)intensity;
                break;
            case NIFTI_TYPE_UINT8:
                resultIntensity[index]=(SourceTYPE)(intensity>0?round(intensity):0);
                break;
            case NIFTI_TYPE_UINT16:
                resultIntensity[index]=(SourceTYPE)(intensity>0?round(intensity):0);
                break;
            case NIFTI_TYPE_UINT32:
                resultIntensity[index]=(SourceTYPE)(intensity>0?round(intensity):0);
                break;
            default:
                resultIntensity[index]=(SourceTYPE)round(intensity);
                break;
            }
        }
    }
}


int get_all_files_and_folders_in_dir (string dir, vector<string> &files , vector<string> &folders)
{
    DIR *dp;
    struct dirent64 *dirp;
    if((dp  = opendir(dir.c_str())) == NULL)
    {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return errno;
    }

    while ((dirp = readdir64(dp)) != NULL)
    {
        if(dirp->d_name[0]!='.')
        {
            if(dirp->d_type==8 ||dirp->d_type==10 ||dirp->d_type==0 )
            {
                if((&files)!=NULL)
                {
                    string curstring=dir;
                    curstring.append(SEP);
                    curstring.append(dirp->d_name);
                    files.push_back(curstring);
                }
            }
            else
            {
                if((&folders)!=NULL)
                {
                    string curstring=dir;
                    curstring.append(SEP);
                    curstring.append(dirp->d_name);
                    folders.push_back(curstring);
                }
            }
        }
    }
    closedir(dp);
    return 0;
}
int get_all_files_that_match_string (string dir, vector<string> &files , string string_to_match)
{
    DIR *dp;
    struct dirent64 *dirp;
    if((dp  = opendir(dir.c_str())) == NULL)
    {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return errno;
    }

    while ((dirp = readdir64(dp)) != NULL)
    {
        if(dirp->d_name[0]!='.')
        {
            if(dirp->d_type==8 ||dirp->d_type==10 ||dirp->d_type==0 )
            {
                string curstring=dirp->d_name;
                if(curstring.find(string_to_match)!=string::npos)
                {
                    files.push_back(dir+string(SEP)+string(dirp->d_name));
                }
            }
        }
    }
    closedir(dp);
    return 0;
}
int get_all_files_that_match_2_strings(string dir, vector<string> &files , string string_to_match, string string_to_match2)
{
    DIR *dp;
    struct dirent64 *dirp;
    if((dp  = opendir(dir.c_str())) == NULL)
    {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return errno;
    }
    while ((dirp = readdir64(dp)) != NULL)
    {
        string curstring=dirp->d_name;
        if(dirp->d_name[0]!='.' ||  curstring.size()>2)
        {
            if(dirp->d_type==8 ||dirp->d_type==10 ||dirp->d_type==0 )
            {
                if((bool)(curstring.find(string_to_match)!=string::npos) && (bool)(curstring.find(string_to_match2)!=string::npos))
                {
                    files.push_back(dir+string(SEP)+string(dirp->d_name));
                }
            }
            //            else{
            //                folders.push_back(dir+string("/")+string(dirp->d_name));
            //            }
        }
    }
    closedir(dp);
    return 0;
}
int get_all_files_in_dir_without_extension(string dir, vector<string> &files)
{
    DIR *dp;
    struct dirent64 *dirp;
    if((dp  = opendir(dir.c_str())) == NULL)
    {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return errno;
    }
    while ((dirp = readdir64(dp)) != NULL)
    {
        if(dirp->d_name[0]!='.')
        {
            if(dirp->d_type==8 ||dirp->d_type==10 ||dirp->d_type==0 )
            {
                string tmpstring=dir+string(SEP)+string(dirp->d_name);
                tmpstring.erase(tmpstring.find_first_of(".",tmpstring.find_last_of(SEP)),tmpstring.length());
                files.push_back(tmpstring);
            }
        }
    }
    closedir(dp);
    return 0;
}


void LTS_Vecs(float * Y, float * X,int * mask, float percentOutliers,int maxNumbIter, float convergenceRatio, unsigned int size, float *a, float *b)
{

    percentOutliers=percentOutliers<0?0:(percentOutliers>0.499)?0.499:percentOutliers;
    float distance_threshold=FLT_MAX;
    int iteration=maxNumbIter;
    float Aval=1;
    float Bval=0;
    float olddistance_threshold=distance_threshold;


    while(iteration)
    {
        float sumX=0;
        float sumY=0;
        float sumXY=0;
        float sumXsquared=0;
        float sizefloat=0;



        for(unsigned int i=0; i<size; i++)
        {
            if(mask==NULL || (mask!=NULL && mask[i]==true))
            {
                if( distance_threshold > fabs(Y[i]-Aval*X[i]+Bval) )
                {
                    sumX+=X[i];
                    sumY+=Y[i];
                    sumXY+=X[i]*Y[i];
                    sumXsquared+=X[i]*X[i];
                    sizefloat++;
                }
            }
        }

        Aval=(sizefloat*sumXY-sumX*sumY)/(sizefloat*sumXsquared-sumX*sumX);
        Bval=(sumY*sumXsquared-sumX*sumXY)/(sizefloat*sumXsquared-sumX*sumX);



        float * distance=new float [(int)(sizefloat)];
        float maxdistance=0;
        int indexmask=0;
        for(unsigned int i=0; i<size; i++)
        {
            if(mask==NULL || (mask!=NULL && mask[i]==true))
            {
                distance[indexmask]=fabs(Y[i]-Aval*X[i]+Bval);
                maxdistance=(maxdistance<distance[indexmask])?distance[indexmask]:maxdistance;
                indexmask++;
            }
        }


        // Fill histogram
        float histsize=200.0f;
        float histo[201]= {0};
        for(long i=0; i<(int)histsize; i++) histo[i]=0;

        for(long i=0; i<indexmask; i++)
        {
            float location=histsize*distance[i]/maxdistance;
            float weight=location-floorf(location);
            histo[(int)floor(location)]+=(1-weight);
            histo[(int)ceil(location)]+=(weight);
        }

        // Normalise histogram
        float sumhisto=0;
        for(long i=0; i<(int)histsize; i++)
            sumhisto+=histo[i];
        for(long i=0; i<(int)histsize; i++)
            histo[i]=histo[i]/sumhisto;


        // Find location of percentile
        float targetcount=(1-percentOutliers*(1.0f/(float)iteration));
        float after_currentcount=0;
        float before_currentcount=0;
        float index_count=histsize;
        for(long i=0; i<(int)histsize; i++)
        {
            before_currentcount=after_currentcount;
            after_currentcount+=histo[i];
            if(after_currentcount>targetcount)
            {
                index_count=i;
                i=(int)histsize;
            }
        }

        // Interpolate Val
        float valbefore=((index_count-1)>=0?(index_count-1):0)*maxdistance/histsize;
        float valafter=(index_count)*maxdistance/histsize;
        float InterpWeight=(targetcount-before_currentcount)/(after_currentcount-before_currentcount);

        // New threshold
        distance_threshold=valafter*InterpWeight+(1-InterpWeight)*valbefore;
        delete [] distance;
        iteration--;
        if( (((olddistance_threshold-distance_threshold)/distance_threshold)<convergenceRatio) && ( fabs(maxNumbIter-iteration)>3 || maxNumbIter<=3))
        {
            iteration=0;
            a[0]=Aval;
            b[0]=Bval;
        }
        else
        {
            a[0]=Aval;
            b[0]=Bval;
        }
        olddistance_threshold=distance_threshold;
    }
}

void LS_Vecs(float * Y, float * X,int * mask, unsigned int size, float *a, float *b)
{

    float sumX=0;
    float sumY=0;
    float sumXY=0;
    float sumXsquared=0;
    float sizefloat=0;
    for(unsigned int i=0; i<size; i++)
    {
        if(mask==NULL || (mask!=NULL && mask[i]==true))
        {
            sumX+=X[i];
            sumY+=Y[i];
            sumXY+=X[i]*Y[i];
            sumXsquared+=X[i]*X[i];
            sizefloat++;
        }
    }

    *a=(sizefloat*sumXY-sumX*sumY)/(sizefloat*sumXsquared-sumX*sumX);
    *b=(sumY*sumXsquared-sumX*sumXY)/(sizefloat*sumXsquared-sumX*sumX);
    //cout << "LS: VAL_A = "<<*a<<" , VAL_B = "<<*b<<endl;

}



void otsu(float * Image,
          int * mask,
          ImageSize *Currentsize )
{

    //find image max and min
    float tempmax=-(1e32);
    float tempmin=1e32;

    for (int i=0; i<Currentsize->numel; i++)
    {
        if(mask==NULL || mask[i]>0)
        {
            if (Image[i]<tempmin)
            {
                tempmin=Image[i];
            }
            if (Image[i]>tempmax)
            {
                tempmax=Image[i];
            }
        }
    }
    //cout<< "min/max = "<<tempmin<<"/"<<tempmax<<endl;


    // Fill histogram
    float histsize=1002.0f;
    float histo[1003]= {0};
    for(long i=0; i<(int)histsize; i++) histo[i]=0;


    for(long i=0; i<Currentsize->numel; i++)
    {
        float location=(histsize-2)*(Image[i]-tempmin)/(tempmax-tempmin)+1;
        float weight=location-floorf(location);
        histo[(int)floor(location)]+=(1-weight);
        histo[(int)ceil(location)]+=(weight);
    }

    // Normalise histogram
    float sumhisto=0;
    for(long i=0; i<(int)histsize; i++)
        sumhisto+=histo[i];
    for(long i=0; i<(int)histsize; i++)
        histo[i]=histo[i]/sumhisto;


    float  w = 0;                // first order cumulative
    float  u = 0;                // second order cumulative
    float  uT = 0;               // total mean level
    float  work1, work2;		// working variables
    double work3 = 0.0;
    float threshold=0;



    // Calculate total mean level
    for (int i=1; i<=histsize; i++)
        uT+=(i*histo[i-1]);


    // Find optimal threshold value
    for (int i=1; i<histsize; i++)
    {
        w+=histo[i-1];
        u+=(i*histo[i-1]);
        work1 = (uT * w - u);
        work2 = (work1 * work1) / ( w * (1.0f-w) );
        if (work2>work3 && work2==work2 && isinf(work2)==0 ){
            work3=work2;
            threshold=(i-1)/(histsize-2)*(tempmax-tempmin)+tempmin;
        }

    }


    // Convert the final value to an integer
    for(long i=0; i<Currentsize->numel; i++)
    {
        Image[i]=Image[i]>threshold;
    }
    return;
}

void ConnectComp(int * Old,
                 int * New,
                 int dimensions[3],int varin)
{

    // Create Counter image
    int index;
    int CCcounter=1;
    int NumElements=((int)dimensions[0]*(int)dimensions[1]*(int)dimensions[2]);

    //  **********    Foreground    ***********
    index=0;
    for(int z=0; z<((int)dimensions[2]); z++)
    {
        for(int y=0; y<((int)dimensions[1]); y++)
        {
            for(int x=0; x<((int)dimensions[0]); x++)
            {
                New[index]=0;
                if(Old[index]>0)
                {
                    New[index]=CCcounter;
                    CCcounter++;
                }
                index++;
            }
        }
    }
    //cout << "Test [" << indexprob<<"] = "<< New[indexprob]<< " "<<Old[indexprob]<<endl;

    int tempmin;
    int numbchanges=1;


    int *CClist = (int *) calloc(CCcounter, sizeof(int));
    for(int i=1; i<CCcounter; i++)
    {
        CClist[i]=i;
    }

    int iter=0;
    while(numbchanges!=0)
    {
        //while(iter<3){

        flush(cout);
        iter++;
        numbchanges=0;
        int currindex=0;
        for(int z=1; z<((int)dimensions[2]-1); z++)
        {
            for(int y=1; y<((int)dimensions[1]-1); y++)
            {
                for(int x=1; x<((int)dimensions[0]-1); x++)
                {
                    index = z*dimensions[1]*dimensions[0]+y*dimensions[0]+x;
                    if(Old[index]>0 && CClist[New[index]]>0)
                    {
                        tempmin=CClist[New[index]];

                        for(int deltaZ=-1; deltaZ<=1; deltaZ+=2)
                        {
                            currindex=index+deltaZ*dimensions[0]*dimensions[1];
                            if(Old[currindex]>0 && tempmin>CClist[New[currindex]])
                            {
                                tempmin=CClist[New[currindex]];
                            }
                        }
                        for(int deltaY=-1; deltaY<=1; deltaY+=2)
                        {
                            currindex=index+deltaY*dimensions[0];
                            if(Old[currindex]>0 && tempmin>CClist[New[currindex]])
                            {
                                tempmin=CClist[New[currindex]];
                            }
                        }
                        for(int deltaX=-1; deltaX<=1; deltaX+=2)
                        {
                            currindex=index+deltaX;
                            if(Old[currindex]>0 && tempmin>CClist[New[currindex]])
                            {
                                tempmin=CClist[New[currindex]];
                            }
                        }
                        if(tempmin>0 && tempmin<CClist[New[index]])
                        {
                            CClist[New[index]]=tempmin;
                            numbchanges++;
                        }
                    }
                }
            }
        }

        for(int index=0; index<((int)dimensions[0]*(int)dimensions[1]*(int)dimensions[2]); index++)
        {
            if(Old[index]>0)
            {
                New[index]=CClist[New[index]];
            }
        }

        //cout << numbchanges<<" "<<New[index]<<" "<<CClist[New[index]]<< endl;
        //probarea=New[index];
    }


    // Update C. Components


    //Find lable counts
    int *Pixelcounter = (int *) calloc((int)CCcounter, sizeof(int));
    int maxForground=0, maxForgroundIndex=0;

    for(index=0; index<NumElements; index++)
    {
        if(New[index]>0 && Old[index]>0)
        {
            Pixelcounter[(int)New[index]]++;
        }
    }

    //cout<<Pixelcounter[probarea]<<endl;
    // If lable touches the edge, errase that lable
    index=0;
    for(int iz=0; iz<dimensions[2]; iz++)
    {
        for(int iy=0; iy<dimensions[1]; iy++)
        {
            for(int ix=0; ix<dimensions[0]; ix++)
            {
                if(((ix==0) || (iy==0) || (iz==0) || (ix==(dimensions[0]-1)) || (iy==(dimensions[1]-1)) || (iz==(dimensions[2]-1))) && (Old[index]>0))
                {
                    Pixelcounter[(int)New[index]]=0;
                }
                index++;
            }
        }
    }
    // Find Biggest Component
    for(index=0; index<CCcounter; index++)
    {
        if(maxForground<Pixelcounter[index])
        {
            maxForgroundIndex=(int)index;
            maxForground=(int)Pixelcounter[index];
        }
    }
    cout << maxForground <<" "<< maxForgroundIndex<<endl;

    //Rassign to oposit class
    for(index=0; index<NumElements; index++)
    {
        if(Old[index]>0)
        {
            if(varin==0)
            {
                if(New[index]==maxForgroundIndex)
                {
                    New[index]=1;
                }
                else
                {
                    New[index]=0;
                }
            }
            else
            {
                if(Pixelcounter[New[index]]>varin)
                {
                    New[index]=1;
                }

                else
                {
                    New[index]=0;
                }
            }
        }
        else
        {
            New[index]=0;
        }
    }
    delete [] Pixelcounter;


    return;
}

template <class OldType, class NewType>
void Close_Forground_ConnectComp(void * Old_void, void * New_void, ImageSize * Currentsize)
{

    OldType * Old = static_cast<OldType *>(Old_void);
    NewType * New = static_cast<NewType *>(New_void);
    int dimensions[3];
    dimensions[0]=Currentsize->xsize;
    dimensions[1]=Currentsize->ysize;
    dimensions[2]=Currentsize->zsize;

    // Create Counter image
    int index;
    int CCcounter=1;
    int NumElements=((int)dimensions[0]*(int)dimensions[1]*(int)dimensions[2]);
    int* tempimg = new int[NumElements];

    //  **********    Foreground    ***********
    index=0;
    for(int z=0; z<((int)dimensions[2]); z++)
    {
        for(int y=0; y<((int)dimensions[1]); y++)
        {
            for(int x=0; x<((int)dimensions[0]); x++)
            {
                tempimg[index]=0;
                if(Old[index]==0)
                {
                    Old[index]=1;
                }
                else
                {
                    Old[index]=0;
                }
                index++;
            }
        }
    }

    index=0;
    for(int z=0; z<((int)dimensions[2]); z++)
    {
        for(int y=0; y<((int)dimensions[1]); y++)
        {
            for(int x=0; x<((int)dimensions[0]); x++)
            {
                tempimg[index]=0;
                if(Old[index]>0)
                {
                    tempimg[index]=CCcounter;
                    CCcounter++;
                }
                index++;
            }
        }
    }
    //cout << "Test [" << indexprob<<"] = "<< New[indexprob]<< " "<<Old[indexprob]<<endl;

    int tempmin;
    int numbchanges=1;


    int *CClist = (int *) calloc(CCcounter, sizeof(int));
    for(int i=1; i<CCcounter; i++)
    {
        CClist[i]=i;
    }

    int iter=0;
    while(numbchanges!=0)
    {
        //while(iter<3){

        flush(cout);
        iter++;
        numbchanges=0;
        int currindex=0;
        for(int z=1; z<((int)dimensions[2]-1); z++)
        {
            for(int y=1; y<((int)dimensions[1]-1); y++)
            {
                for(int x=1; x<((int)dimensions[0]-1); x++)
                {
                    index = z*dimensions[1]*dimensions[0]+y*dimensions[0]+x;
                    if(Old[index]>0 && CClist[tempimg[index]]>0)
                    {
                        tempmin=CClist[tempimg[index]];

                        for(int deltaZ=-1; deltaZ<=1; deltaZ+=2)
                        {
                            currindex=index+deltaZ*dimensions[0]*dimensions[1];
                            if(Old[currindex]>0 && tempmin>CClist[tempimg[currindex]])
                            {
                                tempmin=CClist[tempimg[currindex]];
                            }
                        }
                        for(int deltaY=-1; deltaY<=1; deltaY+=2)
                        {
                            currindex=index+deltaY*dimensions[0];
                            if(Old[currindex]>0 && tempmin>CClist[tempimg[currindex]])
                            {
                                tempmin=CClist[tempimg[currindex]];
                            }
                        }
                        for(int deltaX=-1; deltaX<=1; deltaX+=2)
                        {
                            currindex=index+deltaX;
                            if(Old[currindex]>0 && tempmin>CClist[tempimg[currindex]])
                            {
                                tempmin=CClist[tempimg[currindex]];
                            }
                        }
                        if(tempmin>0 && tempmin<CClist[tempimg[index]])
                        {
                            CClist[tempimg[index]]=tempmin;
                            numbchanges++;
                        }
                    }
                }
            }
        }

        for(int index=0; index<((int)dimensions[0]*(int)dimensions[1]*(int)dimensions[2]); index++)
        {
            if(Old[index]>0)
            {
                tempimg[index]=CClist[tempimg[index]];
            }
        }

        //cout << numbchanges<<" "<<New[index]<<" "<<CClist[New[index]]<< endl;
        //probarea=New[index];
    }

    free(CClist);

    // Update C. Components


    //Find lable counts
    int *Pixelcounter = new int[CCcounter];

    for(index=0; index<NumElements; index++)
    {
        if(tempimg[index]>0 && Old[index]>0)
        {
            Pixelcounter[(int)tempimg[index]]++;
        }
    }

    //cout<<Pixelcounter[probarea]<<endl;
    // If lable touches the edge, errase that lable
    index=0;
    for(int iz=0; iz<dimensions[2]; iz++)
    {
        for(int iy=0; iy<dimensions[1]; iy++)
        {
            for(int ix=0; ix<dimensions[0]; ix++)
            {
                if(((ix==0) || (iy==0) || (iz==0) || (ix==(dimensions[0]-1)) || (iy==(dimensions[1]-1)) || (iz==(dimensions[2]-1))) && (Old[index]>0))
                {
                    Pixelcounter[(int)tempimg[index]]=0;
                }
                index++;
            }
        }
    }
    // Find Biggest Component


    //Rassign to oposit class
    for(index=0; index<NumElements; index++)
    {
        if(Old[index]>0)
        {
            if(Pixelcounter[tempimg[index]]>0)
            {
                New[index]=1;
            }
            else
            {
                New[index]=0;
            }
        }
        else
        {
            New[index]=1;
        }
    }
    delete [] tempimg;
    delete [] Pixelcounter;


    return;
}
template NIFTYSEG_WINEXPORT
void Close_Forground_ConnectComp<unsigned char, unsigned char>(void * Old, void * New, ImageSize * Currentsize);

template NIFTYSEG_WINEXPORT
void Close_Forground_ConnectComp<float, float>(void * Old, void * New, ImageSize * Currentsize);


template <class OldType, class NewType>
void Largest_ConnectComp(void * Old_void, void * New_void, ImageSize * Currentsize)
{

    OldType * Old = static_cast<OldType *>(Old_void);
    NewType * New = static_cast<NewType *>(New_void);
    int dimensions[3];
    dimensions[0]=Currentsize->xsize;
    dimensions[1]=Currentsize->ysize;
    dimensions[2]=Currentsize->zsize;

    // Create Counter image
    int index;
    int CCcounter=1;
    int NumElements=((int)dimensions[0]*(int)dimensions[1]*(int)dimensions[2]);
    int * tempimg= (int *) calloc(NumElements, sizeof(int));

    //  **********    Foreground    ***********
    index=0;
    for(int z=0; z<((int)dimensions[2]); z++)
    {
        for(int y=0; y<((int)dimensions[1]); y++)
        {
            for(int x=0; x<((int)dimensions[0]); x++)
            {
                tempimg[index]=0;
                if(Old[index]==0)
                {
                    Old[index]=0;
                }
                else
                {
                    Old[index]=1;
                }
                index++;
            }
        }
    }

    index=0;
    for(int z=0; z<((int)dimensions[2]); z++)
    {
        for(int y=0; y<((int)dimensions[1]); y++)
        {
            for(int x=0; x<((int)dimensions[0]); x++)
            {
                tempimg[index]=0;
                if(Old[index]>0)
                {
                    tempimg[index]=CCcounter;
                    CCcounter++;
                }
                index++;
            }
        }
    }
    //cout << "Test [" << indexprob<<"] = "<< New[indexprob]<< " "<<Old[indexprob]<<endl;

    int tempmin;
    int numbchanges=1;


    int *CClist = (int *) calloc(CCcounter, sizeof(int));
    for(int i=1; i<CCcounter; i++)
    {
        CClist[i]=i;
    }

    int iter=0;
    while(numbchanges!=0)
    {
        //while(iter<3){

        flush(cout);
        iter++;
        numbchanges=0;
        int currindex=0;
        if((int)dimensions[2]>1){
            for(int z=1; z<((int)dimensions[2]-1); z++)
            {
                for(int y=1; y<((int)dimensions[1]-1); y++)
                {
                    for(int x=1; x<((int)dimensions[0]-1); x++)
                    {
                        index = z*dimensions[1]*dimensions[0]+y*dimensions[0]+x;
                        if(Old[index]>0 && CClist[tempimg[index]]>0)
                        {
                            tempmin=CClist[tempimg[index]];

                            for(int deltaZ=-1; deltaZ<=1; deltaZ+=2)
                            {
                                currindex=index+deltaZ*dimensions[0]*dimensions[1];
                                if(Old[currindex]>0 && tempmin>CClist[tempimg[currindex]])
                                {
                                    tempmin=CClist[tempimg[currindex]];
                                }
                            }
                            for(int deltaY=-1; deltaY<=1; deltaY+=2)
                            {
                                currindex=index+deltaY*dimensions[0];
                                if(Old[currindex]>0 && tempmin>CClist[tempimg[currindex]])
                                {
                                    tempmin=CClist[tempimg[currindex]];
                                }
                            }
                            for(int deltaX=-1; deltaX<=1; deltaX+=2)
                            {
                                currindex=index+deltaX;
                                if(Old[currindex]>0 && tempmin>CClist[tempimg[currindex]])
                                {
                                    tempmin=CClist[tempimg[currindex]];
                                }
                            }
                            if(tempmin>0 && tempmin<CClist[tempimg[index]])
                            {
                                CClist[tempimg[index]]=tempmin;
                                numbchanges++;
                            }
                        }
                    }
                }
            }
        }
        else{
            int z=0;
            for(int y=1; y<((int)dimensions[1]-1); y++)
            {
                for(int x=1; x<((int)dimensions[0]-1); x++)
                {
                    index = z*dimensions[1]*dimensions[0]+y*dimensions[0]+x;
                    if(Old[index]>0 && CClist[tempimg[index]]>0)
                    {
                        tempmin=CClist[tempimg[index]];
                        for(int deltaY=-1; deltaY<=1; deltaY+=2)
                        {
                            currindex=index+deltaY*dimensions[0];
                            if(Old[currindex]>0 && tempmin>CClist[tempimg[currindex]])
                            {
                                tempmin=CClist[tempimg[currindex]];
                            }
                        }
                        for(int deltaX=-1; deltaX<=1; deltaX+=2)
                        {
                            currindex=index+deltaX;
                            if(Old[currindex]>0 && tempmin>CClist[tempimg[currindex]])
                            {
                                tempmin=CClist[tempimg[currindex]];
                            }
                        }
                        if(tempmin>0 && tempmin<CClist[tempimg[index]])
                        {
                            CClist[tempimg[index]]=tempmin;
                            numbchanges++;
                        }
                    }
                }
            }
        }

        for(int index=0; index<((int)dimensions[0]*(int)dimensions[1]*(int)dimensions[2]); index++)
        {
            if(Old[index]>0)
            {
                tempimg[index]=CClist[tempimg[index]];
            }
        }

        //cout << numbchanges<<" "<<New[index]<<" "<<CClist[New[index]]<< endl;
        //probarea=New[index];
    }

    free(CClist);

    // Update C. Components


    //Find lable counts
    int *Pixelcounter = (int *) calloc((int)CCcounter, sizeof(int));
    int maxForground=0;

    for(index=0; index<NumElements; index++)
    {
        if(tempimg[index]>0 && Old[index]>0)
        {
            Pixelcounter[(int)tempimg[index]]++;
        }
    }

    // Find Biggest Component
    for(index=0; index<CCcounter; index++)
    {
        if(maxForground<Pixelcounter[index])
        {
            maxForground=(int)Pixelcounter[index];
        }
    }


    //Rassign to oposit class
    for(index=0; index<NumElements; index++)
    {
        if(Old[index]>0)
        {

            if(Pixelcounter[tempimg[index]]==maxForground)
            {
                New[index]=1;
            }
            else
            {
                New[index]=0;
            }
        }
        else
        {
            New[index]=0;
        }
    }
    delete [] tempimg;
    delete [] Pixelcounter;


    return;
}
template NIFTYSEG_WINEXPORT void Largest_ConnectComp<unsigned char, unsigned char>(void * Old, void * New, ImageSize * Currentsize);
template NIFTYSEG_WINEXPORT void Largest_ConnectComp<float, float>(void * Old, void * New, ImageSize * Currentsize);


template <class OldType, class NewType>
void ConnectComp26NN(void * Old_void, void * New_void, ImageSize * Currentsize)
{

    OldType * Old = static_cast<OldType *>(Old_void);
    NewType * New = static_cast<NewType *>(New_void);
    int dimensions[3];
    dimensions[0]=Currentsize->xsize;
    dimensions[1]=Currentsize->ysize;
    dimensions[2]=Currentsize->zsize;

    // Create Counter image
    int index;
    int CCcounter=1;
    int NumElements=((int)dimensions[0]*(int)dimensions[1]*(int)dimensions[2]);
    int * tempimg= (int *) calloc(NumElements, sizeof(int));

    //  **********    Foreground    ***********
    index=0;
    for(int z=0; z<((int)dimensions[2]); z++)
    {
        for(int y=0; y<((int)dimensions[1]); y++)
        {
            for(int x=0; x<((int)dimensions[0]); x++)
            {
                tempimg[index]=0;
                Old[index]=Old[index]>0;
                index++;
            }
        }
    }

    index=0;
    for(int z=0; z<((int)dimensions[2]); z++)
    {
        for(int y=0; y<((int)dimensions[1]); y++)
        {
            for(int x=0; x<((int)dimensions[0]); x++)
            {
                tempimg[index]=0;
                if(Old[index]>0)
                {
                    tempimg[index]=CCcounter;
                    CCcounter++;
                }
                index++;
            }
        }
    }
    //cout << "Test [" << indexprob<<"] = "<< New[indexprob]<< " "<<Old[indexprob]<<endl;

    int tempmin;
    int numbchanges=1;


    int *CClist = (int *) calloc(CCcounter, sizeof(int));
    for(int i=1; i<CCcounter; i++)
    {
        CClist[i]=i;
    }

    int iter=0;
    while(numbchanges!=0)
    {
        //while(iter<3){
        flush(cout);
        iter++;
        numbchanges=0;
        int currindex=0;
        if((int)dimensions[2]>1){
            for(int z=1; z<((int)dimensions[2]-1); z++)
            {
                for(int y=1; y<((int)dimensions[1]-1); y++)
                {
                    for(int x=1; x<((int)dimensions[0]-1); x++)
                    {
                        index =z*dimensions[1]*dimensions[0]+y*dimensions[0]+x;
                        if(Old[index]>0 && CClist[tempimg[index]]>0)
                        {
                            tempmin=CClist[tempimg[index]];
                            for(int deltaZ=-1; deltaZ<2; deltaZ++)
                            {
                                for(int deltaY=-1; deltaY<2; deltaY++)
                                {
                                    for(int deltaX=-1; deltaX<2; deltaX++)
                                    {
                                        currindex=index+ (deltaZ*dimensions[0]*dimensions[1]) + (deltaY*dimensions[0]) +(deltaX);
                                        if(currindex!=index)
                                        {
                                            if(Old[currindex]>0 && tempmin>CClist[tempimg[currindex]])
                                            {
                                                tempmin=CClist[tempimg[currindex]];
                                            }
                                        }
                                    }
                                }
                            }

                            if(tempmin>0 && tempmin<CClist[tempimg[index]])
                            {
                                CClist[tempimg[index]]=tempmin;
                                numbchanges++;
                            }
                        }
                    }
                }
            }
        }
        else{
            int z=0;
            for(int y=1; y<((int)dimensions[1]-1); y++)
            {
                for(int x=1; x<((int)dimensions[0]-1); x++)
                {
                    index = z*dimensions[1]*dimensions[0]+y*dimensions[0]+x;
                    if(Old[index]>0 && CClist[tempimg[index]]>0)
                    {
                        tempmin=CClist[tempimg[index]];
                        for(int deltaY=-1; deltaY<=1; deltaY+=2)
                        {
                            currindex=index+deltaY*dimensions[0];
                            if(Old[currindex]>0 && tempmin>CClist[tempimg[currindex]])
                            {
                                tempmin=CClist[tempimg[currindex]];
                            }
                        }
                        for(int deltaX=-1; deltaX<=1; deltaX+=2)
                        {
                            currindex=index+deltaX;
                            if(Old[currindex]>0 && tempmin>CClist[tempimg[currindex]])
                            {
                                tempmin=CClist[tempimg[currindex]];
                            }
                        }
                        if(tempmin>0 && tempmin<CClist[tempimg[index]])
                        {
                            CClist[tempimg[index]]=tempmin;
                            numbchanges++;
                        }
                    }
                }
            }
        }

        for(int index=0; index<((int)dimensions[0]*(int)dimensions[1]*(int)dimensions[2]); index++)
        {
            if(Old[index]>0)
            {
                tempimg[index]=CClist[tempimg[index]];
            }
        }
        //probarea=New[index];
    }



    //Find lable counts
    float *Pixelcounter = (float *) calloc((int)CCcounter, sizeof(float));
    int *PixelcounterOrder;

    for(index=0; index<NumElements; index++)
    {
        if(tempimg[index]>0 && Old[index]>0)
        {
            Pixelcounter[(int)tempimg[index]]++;
        }
    }
    PixelcounterOrder=quickSort_order(Pixelcounter,CCcounter);

    //Rassign to oposit class
    for(index=0; index<NumElements; index++)
    {
        if(tempimg[index]>0 && Old[index]>0)
        {
            New[index]=CCcounter-PixelcounterOrder[(int)tempimg[index]];
        }
        else{
            New[index]=0;
        }
    }
    delete [] tempimg;

    return;
}
template NIFTYSEG_WINEXPORT void ConnectComp26NN<unsigned char, unsigned char>(void * Old, void * New, ImageSize * Currentsize);
template NIFTYSEG_WINEXPORT void ConnectComp26NN<float, float>(void * Old, void * New, ImageSize * Currentsize);

template <class OldType, class NewType>
void ConnectComp6NN(void * Old_void, void * New_void, ImageSize * Currentsize)
{

    OldType * Old = static_cast<OldType *>(Old_void);
    NewType * New = static_cast<NewType *>(New_void);
    int dimensions[3];
    dimensions[0]=Currentsize->xsize;
    dimensions[1]=Currentsize->ysize;
    dimensions[2]=Currentsize->zsize;

    // Create Counter image
    int index;
    int CCcounter=1;
    int NumElements=((int)dimensions[0]*(int)dimensions[1]*(int)dimensions[2]);
    int * tempimg= (int *) calloc(NumElements, sizeof(int));

    //  **********    Foreground    ***********
    index=0;
    for(int z=0; z<((int)dimensions[2]); z++)
    {
        for(int y=0; y<((int)dimensions[1]); y++)
        {
            for(int x=0; x<((int)dimensions[0]); x++)
            {
                tempimg[index]=0;
                if(Old[index]==0)
                {
                    Old[index]=0;
                }
                else
                {
                    Old[index]=1;
                }
                index++;
            }
        }
    }

    index=0;
    for(int z=0; z<((int)dimensions[2]); z++)
    {
        for(int y=0; y<((int)dimensions[1]); y++)
        {
            for(int x=0; x<((int)dimensions[0]); x++)
            {
                tempimg[index]=0;
                if(Old[index]>0)
                {
                    tempimg[index]=CCcounter;
                    CCcounter++;
                }
                index++;
            }
        }
    }


    int tempmin;
    int numbchanges=1;


    int *CClist = (int *) calloc(CCcounter, sizeof(int));
    for(int i=1; i<CCcounter; i++)
    {
        CClist[i]=i;
    }

    int iter=0;
    while(numbchanges!=0)
    {
        //while(iter<3){

        flush(cout);
        iter++;
        numbchanges=0;
        int currindex=0;
        if((int)dimensions[2]>1){
            for(int z=1; z<((int)dimensions[2]-1); z++)
            {
                for(int y=1; y<((int)dimensions[1]-1); y++)
                {
                    for(int x=1; x<((int)dimensions[0]-1); x++)
                    {
                        index = z*dimensions[1]*dimensions[0]+y*dimensions[0]+x;
                        if(Old[index]>0 && CClist[tempimg[index]]>0)
                        {
                            tempmin=CClist[tempimg[index]];

                            for(int deltaZ=-1; deltaZ<=1; deltaZ+=2)
                            {
                                currindex=index+deltaZ*dimensions[0]*dimensions[1];
                                if(Old[currindex]>0 && tempmin>CClist[tempimg[currindex]])
                                {
                                    tempmin=CClist[tempimg[currindex]];
                                }
                            }
                            for(int deltaY=-1; deltaY<=1; deltaY+=2)
                            {
                                currindex=index+deltaY*dimensions[0];
                                if(Old[currindex]>0 && tempmin>CClist[tempimg[currindex]])
                                {
                                    tempmin=CClist[tempimg[currindex]];
                                }
                            }
                            for(int deltaX=-1; deltaX<=1; deltaX+=2)
                            {
                                currindex=index+deltaX;
                                if(Old[currindex]>0 && tempmin>CClist[tempimg[currindex]])
                                {
                                    tempmin=CClist[tempimg[currindex]];
                                }
                            }
                            if(tempmin>0 && tempmin<CClist[tempimg[index]])
                            {
                                CClist[tempimg[index]]=tempmin;
                                numbchanges++;
                            }
                        }
                    }
                }
            }
        }
        else{
            int z=0;
            for(int y=1; y<((int)dimensions[1]-1); y++)
            {
                for(int x=1; x<((int)dimensions[0]-1); x++)
                {
                    index = z*dimensions[1]*dimensions[0]+y*dimensions[0]+x;
                    if(Old[index]>0 && CClist[tempimg[index]]>0)
                    {
                        tempmin=CClist[tempimg[index]];
                        for(int deltaY=-1; deltaY<=1; deltaY+=2)
                        {
                            currindex=index+deltaY*dimensions[0];
                            if(Old[currindex]>0 && tempmin>CClist[tempimg[currindex]])
                            {
                                tempmin=CClist[tempimg[currindex]];
                            }
                        }
                        for(int deltaX=-1; deltaX<=1; deltaX+=2)
                        {
                            currindex=index+deltaX;
                            if(Old[currindex]>0 && tempmin>CClist[tempimg[currindex]])
                            {
                                tempmin=CClist[tempimg[currindex]];
                            }
                        }
                        if(tempmin>0 && tempmin<CClist[tempimg[index]])
                        {
                            CClist[tempimg[index]]=tempmin;
                            numbchanges++;
                        }
                    }
                }
            }
        }

        for(int index=0; index<((int)dimensions[0]*(int)dimensions[1]*(int)dimensions[2]); index++)
        {
            if(Old[index]>0)
            {
                tempimg[index]=CClist[tempimg[index]];
            }
        }

        //cout << numbchanges<<" "<<New[index]<<" "<<CClist[New[index]]<< endl;
        //probarea=New[index];
    }



    //Find lable counts
    float *Pixelcounter = (float *) calloc((int)CCcounter, sizeof(float));
    int *PixelcounterOrder;

    for(index=0; index<NumElements; index++)
    {
        if(tempimg[index]>0 && Old[index]>0)
        {
            Pixelcounter[(int)tempimg[index]]++;
        }
    }

    PixelcounterOrder=quickSort_order(Pixelcounter,CCcounter);


    //Rassign to oposit class
    for(index=0; index<NumElements; index++)
    {
        if(tempimg[index]>0 && Old[index]>0)
        {
            New[index]=CCcounter-PixelcounterOrder[(int)tempimg[index]];
        }
        else{
            New[index]=0;
        }
    }
    delete [] tempimg;

    return;
}
template NIFTYSEG_WINEXPORT void ConnectComp6NN<unsigned char, unsigned char>(void * Old, void * New, ImageSize * Currentsize);
template NIFTYSEG_WINEXPORT void ConnectComp6NN<float, float>(void * Old, void * New, ImageSize * Currentsize);


void Dillate(float * Image,
             int kernel,
             ImageSize * Currentsize)
{

    int dimensions[3];
    dimensions[0]=Currentsize->xsize;
    dimensions[1]=Currentsize->ysize;
    dimensions[2]=Currentsize->zsize;
    int xyzpos[3];
    int shiftdirection[3];
    shiftdirection[0]=1;
    shiftdirection[1]=dimensions[0];
    shiftdirection[2]=dimensions[1]*dimensions[0];
    float * Buffer = new float [(dimensions[1])*(dimensions[0])*(dimensions[2])];

    if(Buffer == NULL)
    {
        fprintf(stderr,"* Error when alocating Buffer in void Dillate(): Not enough memory\n");
        exit(1);
    }

    float tmpvalue=0;
    for(int tp=0; tp<Currentsize->tsize; tp++){
        for(int currentdirection=0; currentdirection<3; currentdirection++) //Buffer aint each direction
        {
            int index=0;
            for(xyzpos[2]=0; xyzpos[2]<dimensions[2]; xyzpos[2]++)
            {
                for(xyzpos[1]=0; xyzpos[1]<dimensions[1]; xyzpos[1]++)
                {
                    for(xyzpos[0]=0; xyzpos[0]<dimensions[0]; xyzpos[0]++)
                    {
                        tmpvalue=-1e32;
                        for(int shift=((xyzpos[currentdirection]<kernel)?-xyzpos[currentdirection]:-kernel); shift<=((xyzpos[currentdirection]>=(dimensions[currentdirection]-kernel))?(int)dimensions[currentdirection]-xyzpos[currentdirection]-1:kernel); shift++)
                        {
                            if(Image[index+shift*shiftdirection[currentdirection]+tp*dimensions[0]*dimensions[1]*dimensions[2]]>tmpvalue)
                            {
                                tmpvalue=Image[index+shift*shiftdirection[currentdirection]+tp*dimensions[0]*dimensions[1]*dimensions[2]];
                            }

                        }
                        Buffer[index]=tmpvalue;
                        index++;
                    }
                }
            }
            for(int i=0; i<(dimensions[1]*dimensions[0]*dimensions[2]); i++)
            {
                Image[i+tp*dimensions[0]*dimensions[1]*dimensions[2]]=Buffer[i];
            }
        }
    }

    delete [] Buffer;

    return;
}


void Erosion(float * Image,
             int kernel,
             ImageSize *Currentsize )
{

    int dimensions[3];
    dimensions[0]=Currentsize->xsize;
    dimensions[1]=Currentsize->ysize;
    dimensions[2]=Currentsize->zsize;
    int xyzpos[3];
    int shiftdirection[3];
    shiftdirection[0]=1;
    shiftdirection[1]=dimensions[0];
    shiftdirection[2]=dimensions[1]*dimensions[0];
    float * Buffer = new float [(dimensions[1])*(dimensions[0])*(dimensions[2])];
    if(Buffer == NULL)
    {
        fprintf(stderr,"* Error when alocating Buffer in void Dillate(): Not enough memory\n");
        exit(1);
    }

    float tmpvalue=0;
    for(int tp=0; tp<Currentsize->tsize; tp++){
        for(int currentdirection=0; currentdirection<3; currentdirection++) //Buffer aint each direction
        {
            int index=0;
            for(xyzpos[2]=0; xyzpos[2]<dimensions[2]; xyzpos[2]++)
            {
                for(xyzpos[1]=0; xyzpos[1]<dimensions[1]; xyzpos[1]++)
                {
                    for(xyzpos[0]=0; xyzpos[0]<dimensions[0]; xyzpos[0]++)
                    {
                        tmpvalue=1e32;
                        for(int shift=((xyzpos[currentdirection]<kernel)?-xyzpos[currentdirection]:-kernel); shift<=((xyzpos[currentdirection]>=(dimensions[currentdirection]-kernel))?(int)dimensions[currentdirection]-xyzpos[currentdirection]-1:kernel); shift++)
                        {
                            if(Image[index+shift*shiftdirection[currentdirection]+tp*dimensions[0]*dimensions[1]*dimensions[2]]<tmpvalue)
                            {
                                tmpvalue=Image[index+shift*shiftdirection[currentdirection]+tp*dimensions[0]*dimensions[1]*dimensions[2]];
                            }

                        }
                        Buffer[index]=tmpvalue;
                        index++;
                    }
                }
            }
            for(int i=0; i<(dimensions[1]*dimensions[0]*dimensions[2]); i++)
            {
                Image[i+tp*dimensions[0]*dimensions[1]*dimensions[2]]=Buffer[i];
            }
        }
    }
    delete [] Buffer;

    return;
}

bool isSimplePoint(bool * SimplePointTestBlock){

  bool CounterBlock[27]={0};
  char CCcounter=0;
  int index=0;
  for(int z=0; z<3; z++)
  {
    for(int y=0; y<3; y++)
    {
      for(int x=0; x<3; x++)
      {
        if(SimplePointTestBlock[index]>0){
          CounterBlock[index]=CCcounter;
          CCcounter++;
          index++;
        }
        else{
          CounterBlock[index]=0;
        }
      }
    }
  }
  int CClist[27]={0};
  for(int i=0; i<27; i++)
  {
    CClist[i]=CounterBlock[i]?i:-1;
  }

  int iter=0;
  int numbchanges=1;
  int totalnumbchanges=0;
  while(numbchanges!=0)
  {
    //while(iter<3){

    flush(cout);
    iter++;
    numbchanges=0;
    int currindex=0;

    index=0;
    int tempmin;
    for(int z=1; z<3; z++)
    {
      for(int y=1; y<3; y++)
      {
        for(int x=1; x<3; x++)
        {
          if(SimplePointTestBlock[index]>0 && CClist[CounterBlock[index]]>0)
          {
            tempmin=CClist[CounterBlock[index]];

            for(int deltaZ=-1; deltaZ<=1; deltaZ+=2)
            {
              currindex=index+deltaZ*9;
              if(SimplePointTestBlock[currindex]>0 && tempmin>CClist[CounterBlock[currindex]])
              {
                tempmin=CClist[CounterBlock[currindex]];
              }
            }
            for(int deltaY=-1; deltaY<=1; deltaY+=2)
            {
              currindex=index+deltaY*3;
              if(SimplePointTestBlock[currindex]>0 && tempmin>CClist[CounterBlock[currindex]])
              {
                tempmin=CClist[CounterBlock[currindex]];
              }
            }
            for(int deltaX=-1; deltaX<=1; deltaX+=2)
            {
              currindex=index+deltaX;
              if(SimplePointTestBlock[currindex]>0 && tempmin>CClist[CounterBlock[currindex]])
              {
                tempmin=CClist[CounterBlock[currindex]];
              }
            }
            if(tempmin>0 && tempmin<CClist[CounterBlock[index]])
            {
              CClist[CounterBlock[index]]=tempmin;
              numbchanges++;
            }
          }
          index++;
        }
      }
    }

    for(int i=0; i<27; i++)
    {
      if(CClist[i]==i && CClist[i]>=0)
      {
      totalnumbchanges += 1;
      }
    }

  }


  bool SimplePointTestBlockMod[27]={0};
  CCcounter=0;
  index=0;
  for(int z=0; z<3; z++)
  {
    for(int y=0; y<3; y++)
    {
      for(int x=0; x<3; x++)
      {
        SimplePointTestBlockMod[index]=(z==1 && y==1 && x==1)?0:SimplePointTestBlock[index];

        if(SimplePointTestBlockMod[index]>0){
          CounterBlock[index]=CCcounter;
          CCcounter++;
          index++;
        }
        else{
          CounterBlock[index]=0;
        }
      }
    }
  }
  for(int i=0; i<27; i++)
  {
    CClist[i]=CounterBlock[i]?i:-1;
  }

  iter=0;
  numbchanges=1;
  int totalnumbchanges2=0;
  while(numbchanges!=0)
  {

    flush(cout);
    iter++;
    numbchanges=0;
    int currindex=0;

    index=0;
    int tempmin;
    for(int z=1; z<3; z++)
    {
      for(int y=1; y<3; y++)
      {
        for(int x=1; x<3; x++)
        {
          if(SimplePointTestBlockMod[index]>0 && CClist[CounterBlock[index]]>0)
          {
            tempmin=CClist[CounterBlock[index]];

            for(int deltaZ=-1; deltaZ<=1; deltaZ+=2)
            {
              currindex=index+deltaZ*9;
              if(SimplePointTestBlockMod[currindex]>0 && tempmin>CClist[CounterBlock[currindex]])
              {
                tempmin=CClist[CounterBlock[currindex]];
              }
            }
            for(int deltaY=-1; deltaY<=1; deltaY+=2)
            {
              currindex=index+deltaY*3;
              if(SimplePointTestBlockMod[currindex]>0 && tempmin>CClist[CounterBlock[currindex]])
              {
                tempmin=CClist[CounterBlock[currindex]];
              }
            }
            for(int deltaX=-1; deltaX<=1; deltaX+=2)
            {
              currindex=index+deltaX;
              if(SimplePointTestBlockMod[currindex]>0 && tempmin>CClist[CounterBlock[currindex]])
              {
                tempmin=CClist[CounterBlock[currindex]];
              }
            }
            if(tempmin>0 && tempmin<CClist[CounterBlock[index]])
            {
              CClist[CounterBlock[index]]=tempmin;
              numbchanges++;
            }
          }
          index++;
        }
      }
    }

    for(int i=0; i<27; i++)
    {
      if(CClist[i]==i && CClist[i]>=0)
      {
      totalnumbchanges2 += 1;
      }
    }

  }

  if(totalnumbchanges2>0 && totalnumbchanges>0){
    cout << totalnumbchanges2<< "  "<< totalnumbchanges<<endl;
  }
  return (totalnumbchanges2==totalnumbchanges);
}

//void TopologicalErosion(float * Image,
//             int kernel,
//             ImageSize *Currentsize )
//{

//    int dimensions[3];
//    dimensions[0]=Currentsize->xsize;
//    dimensions[1]=Currentsize->ysize;
//    dimensions[2]=Currentsize->zsize;
//    int xyzpos[3];
//    int shiftdirection[3];
//    shiftdirection[0]=1;
//    shiftdirection[1]=dimensions[0];
//    shiftdirection[2]=dimensions[1]*dimensions[0];
//    float * Buffer = new float [(dimensions[1])*(dimensions[0])*(dimensions[2])];

//    if(Buffer == NULL)
//    {
//        fprintf(stderr,"* Error when alocating Buffer in void Dillate(): Not enough memory\n");
//        exit(1);
//    }

//    for(int tp=0; tp<Currentsize->tsize; tp++){
//      int index=0;
//      for(xyzpos[2]=0; xyzpos[2]<dimensions[2]; xyzpos[2]++)
//      {
//        for(xyzpos[1]=0; xyzpos[1]<dimensions[1]; xyzpos[1]++)
//        {
//          for(xyzpos[0]=0; xyzpos[0]<dimensions[0]; xyzpos[0]++)
//          {
//            if(Image[index]){
//              bool SimplePointTestBlock[27]={0};
//              for(char dx=-1; dx<2; dx++){
//                for(char dy=-1; dy<2; dy++){
//                  for(char dz=-1; dz<2; dz++){
//                    if(Image[(xyzpos[0]+dx)+((xyzpos[1]+dy)+(xyzpos[2]+dz)*dimensions[1])*dimensions[0]]>0 &&
//                       (xyzpos[0]+dx)>0 && (xyzpos[0]+dx)<(dimensions[0]-1)
//                       && (xyzpos[1]+dy)>0 && (xyzpos[1]+dy)<(dimensions[1]-1)
//                       && (xyzpos[2]+dz)>0 && (xyzpos[2]+dz)<(dimensions[2]-1)){
//                      SimplePointTestBlock[(dx+1)+(dy+1)*3+(dz+1)*9]=Image[index+dx+(dy+(dz*dimensions[1]))*dimensions[0]];
//                    }
//                  }
//                }
//              }
//              Buffer[index]=isSimplePoint(SimplePointTestBlock);
//            }
//            else{
//              Buffer[index]=0;
//            }
//            index++;
//          }
//        }
//      }
//    }
//    for(int tp=0; tp<Currentsize->tsize; tp++){
//      int index=0;
//      for(xyzpos[2]=0; xyzpos[2]<dimensions[2]; xyzpos[2]++)
//      {
//        for(xyzpos[1]=0; xyzpos[1]<dimensions[1]; xyzpos[1]++)
//        {
//          for(xyzpos[0]=0; xyzpos[0]<dimensions[0]; xyzpos[0]++)
//          {
//            Image[index]=Buffer[index];
//            index++;
//          }
//        }
//      }
//    }
//    delete [] Buffer;

//    return;
//}


void Dillate_const(bool * Image,
                   bool * Const,
                   int kernel,
                   int dimensions[3],
int direction)
{
    int xyzpos[3];
    int shiftdir=0;
    if(fabs(direction)==1)
    {
        shiftdir=1;
    }
    if(fabs(direction)==2)
    {
        shiftdir=dimensions[0];
    }
    if(fabs(direction)==3)
    {
        shiftdir=dimensions[1]*dimensions[0];
    }

    bool * Buffer = new bool [dimensions[1]*dimensions[0]*dimensions[2]]();
    bool tmpvalue=0;

    int index=0;
    for(xyzpos[2]=0; xyzpos[2]<dimensions[2]; xyzpos[2]++)
    {
        for(xyzpos[1]=0; xyzpos[1]<dimensions[1]; xyzpos[1]++)
        {
            for(xyzpos[0]=0; xyzpos[0]<dimensions[0]; xyzpos[0]++)
            {
                tmpvalue=false;
                if(direction>0)
                {
                    for(int shift=0; shift<=((xyzpos[(int)fabs(direction)]>=(dimensions[(int)fabs(direction)]-kernel))?(int)dimensions[(int)fabs(direction)]-xyzpos[(int)fabs(direction)]-1:kernel); shift++)
                    {
                        if(Image[index+shift*shiftdir])
                        {
                            tmpvalue=true;
                            shift=10000;
                        }
                    }
                }
                if(direction>0)
                {
                    for(int shift=((xyzpos[(int)fabs(direction)]<kernel)?-xyzpos[(int)fabs(direction)]:-kernel); shift<=0; shift++)
                    {
                        if(Image[index+shift*shiftdir])
                        {
                            tmpvalue=true;
                            shift=10000;
                        }
                    }
                }
                Buffer[index]=tmpvalue;
                index++;
            }
        }
    }
    for(int i=0; i<(dimensions[1]*dimensions[0]*dimensions[2]); i++)
    {
        Image[i]=Buffer[i];
    }
    index=0;
    for(xyzpos[2]=0; xyzpos[2]<dimensions[2]; xyzpos[2]++)
    {
        for(xyzpos[1]=0; xyzpos[1]<dimensions[1]; xyzpos[1]++)
        {
            for(xyzpos[0]=0; xyzpos[0]<dimensions[0]; xyzpos[0]++)
            {
                if(Image[index]>0)
                {
                    tmpvalue=true;
                    if(direction>0)
                    {
                        for(int shift=0; shift<=((xyzpos[(int)fabs(direction)]>=(dimensions[(int)fabs(direction)]-kernel))?(int)dimensions[(int)fabs(direction)]-xyzpos[(int)fabs(direction)]-1:kernel); shift++)
                        {
                            if(Image[index+shift*shiftdir]==0 && Const[index+shift*shiftdir]==0)
                            {
                                tmpvalue=false;
                                shift=10000;
                            }
                        }
                    }
                    if(direction>0)
                    {
                        for(int shift=((xyzpos[(int)fabs(direction)]<kernel)?-xyzpos[(int)fabs(direction)]:-kernel); shift<=0; shift++)
                        {
                            if(Image[index+shift*shiftdir]==0 && Const[index+shift*shiftdir]==0)
                            {
                                tmpvalue=false;
                                shift=10000;
                            }
                        }
                    }
                    Buffer[index]=tmpvalue;
                }
                else
                {
                    Buffer[index]=0;
                }
                index++;
            }
        }
    }
    for(int i=0; i<(dimensions[1]*dimensions[0]*dimensions[2]); i++)
    {
        Image[i]=Buffer[i];
    }


    delete [] Buffer;

    return;
}


/* *************************************************************** */
/* *************************************************************** */
template <class DTYPE>
void seg_mat44_mul(mat44 const* mat,
                   DTYPE const* in,
                   DTYPE *out)
{
    out[0]= (DTYPE)mat->m[0][0]*in[0] +
            (DTYPE)mat->m[0][1]*in[1] +
            (DTYPE)mat->m[0][2]*in[2] +
            (DTYPE)mat->m[0][3];
    out[1]= (DTYPE)mat->m[1][0]*in[0] +
            (DTYPE)mat->m[1][1]*in[1] +
            (DTYPE)mat->m[1][2]*in[2] +
            (DTYPE)mat->m[1][3];
    out[2]= (DTYPE)mat->m[2][0]*in[0] +
            (DTYPE)mat->m[2][1]*in[1] +
            (DTYPE)mat->m[2][2]*in[2] +
            (DTYPE)mat->m[2][3];
    return;
}
template NIFTYSEG_WINEXPORT void seg_mat44_mul<float>(mat44 const*, float const*, float*);
template NIFTYSEG_WINEXPORT void seg_mat44_mul<double>(mat44 const*, double const*, double*);
