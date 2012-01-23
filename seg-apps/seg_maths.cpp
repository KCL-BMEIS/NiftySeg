#include <iostream>
#include <time.h>
#include "_seg_common.h"
#include "_seg_tools.h"
#include "_seg_Topo.h"

using namespace std;
#define PrecisionTYPE float

void Usage(char *exec)
{
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    printf("Usage:\t%s <input> <operation> <output>.\n\n",exec);
    printf("\t* * Operations on 3-D and 4-D images* *\n");
    printf("\t-mul\t<float/file>\tMultiply image <float> value or by other image.\n");
    printf("\t-div\t<float/file>\tDivide image by <float> or by other image.\n");
    printf("\t-add\t<float/file>\tAdd image by <float> or by other image.\n");
    printf("\t-sub\t<float/file>\tSubtract image by <float> or by other image.\n");
    printf("\t-pow\t<float>\t\tImage to the power of <float>.\n");
    printf("\t-thr\t<float>\t\tThreshold the image below <float>.\n");
    printf("\t-uthr\t<float>\t\tThreshold image above <float>.\n");
    printf("\t-smo\t<float>\t\tGaussian smoothing by std <float> (in voxels and up to 4-D).\n");
    printf("\t-sqrt \t\t\tSquare root of the image.\n");
    printf("\t-exp \t\t\tExponential root of the image.\n");
    printf("\t-recip \t\t\tReciprocal (1/I) of the image.\n");
    printf("\t-abs \t\t\tAbsolute value of the image.\n");
    printf("\t-bin \t\t\tBinarise the image.\n");
    printf("\n\t* * Operations on 3-D images * *\n");
    printf("\t-dil\t<int>\t\tDilate the image <int> times (in voxels).\n");
    printf("\t-ero\t<int>\t\tErode the image <int> times (in voxels).\n");
    printf("\t-lconcomp\t\tTake the largest connected component\n");
    printf("\t-fill\t\t\tFill holes in binary object (e.g. fill ventricle in brain mask).\n");
    printf("\n\t* * Dimensionality reduction operations: from 4-D to 3-D * *\n");
    printf("\t-tp\t<int>\t\tExtract time point <int>\n");
    printf("\t-tmean\t<int>\t\tMean value of all time points.\n");
    printf("\t-tmax\t<int>\t\tMax value of all time points.\n");
    printf("\t-tmin\t\t\tMean value of all time points.\n");
    printf("\n\t* * Dimensionality increase operations: from 3-D to 4-D * *\n");
    printf("\t-merge\t<int> <files>\tMerge <int> images and the working image in time (4-D)\n");
    printf("\n\t* * Image similarity: Local metrics * *\n");
    printf("\t-lncc\t<file> <std>\tLocal CC between current img and <int> on a kernel with <std>\n");
    printf("\t-lssd\t<file> <std>\tLocal SSD between current img and <int> on a kernel with <std>\n");
    //printf("\t-lmi\t<file> <std>\tLocal MI between current img and <int> on a kernel with <std>\n");
    printf("\n\t* * Image header operations * *\n");
    printf("\t-hdr_copy <file> \tCopy header from working image to <file> and save in <output>.\n");
    printf("\n\t* * Datatype output * *\n");
    printf("\t-odt <datatype> \tSet output <datatype> (char, short, int, uchar, ushort, uint, float, double).\n");
    printf("\n\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    return;
}
int main(int argc, char **argv)
{

    if (argc <= 2)
    {
        Usage(argv[0]);
        return 0;
    }
    if(strcmp(argv[1], "-help")==0 || strcmp(argv[1], "-Help")==0 ||
            strcmp(argv[1], "-HELP")==0 || strcmp(argv[1], "-h")==0 ||
            strcmp(argv[1], "--h")==0 || strcmp(argv[1], "--help")==0){
        Usage(argv[0]);
        return 0;
    }


    char * filename_in=argv[1];
    nifti_image * InputImage=nifti_image_read(filename_in,true);
    if(InputImage == NULL){
        fprintf(stderr,"* Error when reading the input Segmentation image\n");
        return 1;
    }
    if(InputImage->datatype!=DT_FLOAT32){
        seg_changeDatatype<PrecisionTYPE>(InputImage);
    }
    PrecisionTYPE * InputImagePtr = static_cast<PrecisionTYPE *>(InputImage->data);
    ImageSize * CurrSize = new ImageSize [1]();
    CurrSize->numel=(int)(InputImage->nx*InputImage->ny*InputImage->nz);
    CurrSize->xsize=InputImage->nx;
    CurrSize->ysize=InputImage->ny;
    CurrSize->zsize=InputImage->nz;
    CurrSize->usize=(InputImage->nu>1)?InputImage->nu:1;
    CurrSize->tsize=(InputImage->nt>1)?InputImage->nt:1;
    int datatypeoutput=NIFTI_TYPE_FLOAT32;

    PrecisionTYPE ** bufferImages = new PrecisionTYPE * [2];
    bufferImages[0] = new PrecisionTYPE [InputImage->nvox];
    bufferImages[1] = new PrecisionTYPE [InputImage->nvox];
    for(unsigned int i=0; i<InputImage->nvox; i++){
        bufferImages[0][i]=InputImagePtr[i];
    }
    int current_buffer=0;



    for(int i=2;i<(argc-1);i++){
        if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
                strcmp(argv[i], "-HELP")==0 || strcmp(argv[i], "-h")==0 ||
                strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0){
            Usage(argv[0]);
            return 0;
        }
        // *********************  MUTIPLY  *************************
        else if(strcmp(argv[i], "-mul") == 0){
            string parser=argv[++i];
            if(strtod(parser.c_str(),NULL)){
                double multfactor=strtod(parser.c_str(),NULL);
                for(unsigned int i=0; i<InputImage->nvox; i++)
                    bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i]*multfactor;
                current_buffer=current_buffer?0:1;
            }
            else{
                nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
                if(InputImage->datatype!=DT_FLOAT32){
                    seg_changeDatatype<PrecisionTYPE>(NewImage);
                }
                PrecisionTYPE * NewImagePtr = static_cast<PrecisionTYPE *>(NewImage->data);
                if(NewImage->nx==InputImage->nx&&NewImage->ny==InputImage->ny&&NewImage->nz==InputImage->nz&&NewImage->nt==InputImage->nt&&NewImage->nu==InputImage->nu&&NewImage->nv==InputImage->nv&&NewImage->nw==InputImage->nw){
                    for(unsigned int i=0; i<InputImage->nvox; i++)
                        bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i]*NewImagePtr[i];

                    current_buffer=current_buffer?0:1;
                }
                else{
                    cout << "ERROR: Image "<< parser << " is the wrong size"<<endl;
                    i=argc;
                }
                nifti_image_free(NewImage);
            }
        }
        // *********************  ADD  *************************
        else if(strcmp(argv[i], "-add") == 0){
            string parser=argv[++i];
            if(strtod(parser.c_str(),NULL)){
                double addfactor=strtod(parser.c_str(),NULL);
                for(unsigned int i=0; i<InputImage->nvox; i++)
                    bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i]+addfactor;
                current_buffer=current_buffer?0:1;
            }
            else{
                nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
                if(NewImage->datatype!=DT_FLOAT32){
                    seg_changeDatatype<PrecisionTYPE>(NewImage);
                }
                PrecisionTYPE * NewImagePtr = static_cast<PrecisionTYPE *>(NewImage->data);
                if(NewImage->nx==InputImage->nx&&NewImage->ny==InputImage->ny&&NewImage->nz==InputImage->nz&&NewImage->nt==InputImage->nt&&NewImage->nu==InputImage->nu&&NewImage->nv==InputImage->nv&&NewImage->nw==InputImage->nw){
                    for(unsigned int i=0; i<InputImage->nvox; i++)
                        bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i]+NewImagePtr[i];
                    current_buffer=current_buffer?0:1;
                }
                else{
                    cout << "ERROR: Image "<< parser << " is the wrong size"<<endl;
                    i=argc;
                }
                nifti_image_free(NewImage);
            }
        }
        // *********************  SUBTRACT  *************************
        else if(strcmp(argv[i], "-sub") == 0){
            string parser=argv[++i];
            if(strtod(parser.c_str(),NULL)){
                double factor=strtod(parser.c_str(),NULL);
                for(unsigned int i=0; i<InputImage->nvox; i++)
                    bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i]-factor;
                current_buffer=current_buffer?0:1;
            }
            else{
                nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
                if(NewImage->datatype!=DT_FLOAT32){
                    seg_changeDatatype<PrecisionTYPE>(NewImage);
                }
                PrecisionTYPE * NewImagePtr = static_cast<PrecisionTYPE *>(NewImage->data);
                if(NewImage->nx==InputImage->nx&&NewImage->ny==InputImage->ny&&NewImage->nz==InputImage->nz&&NewImage->nt==InputImage->nt&&NewImage->nu==InputImage->nu&&NewImage->nv==InputImage->nv&&NewImage->nw==InputImage->nw){
                    for(unsigned int i=0; i<InputImage->nvox; i++)
                        bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i]-NewImagePtr[i];
                    current_buffer=current_buffer?0:1;
                }
                else{
                    cout << "ERROR: Image "<< parser << " is the wrong size"<<endl;
                    i=argc;
                }
                nifti_image_free(NewImage);
            }
        }
        // *********************  DIV  *************************
        else if(strcmp(argv[i], "-sub") == 0){
            string parser=argv[++i];
            if(strtod(parser.c_str(),NULL)){
                double factor=strtod(parser.c_str(),NULL);
                for(unsigned int i=0; i<InputImage->nvox; i++)
                    bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i]/factor;
                current_buffer=current_buffer?0:1;
            }
            else{
                nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
                if(NewImage->datatype!=DT_FLOAT32){
                    seg_changeDatatype<PrecisionTYPE>(NewImage);
                }
                PrecisionTYPE * NewImagePtr = static_cast<PrecisionTYPE *>(NewImage->data);
                if(NewImage->nx==InputImage->nx&&NewImage->ny==InputImage->ny&&NewImage->nz==InputImage->nz&&NewImage->nt==InputImage->nt&&NewImage->nu==InputImage->nu&&NewImage->nv==InputImage->nv&&NewImage->nw==InputImage->nw){
                    for(unsigned int i=0; i<InputImage->nvox; i++)
                        bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i]/NewImagePtr[i];
                    current_buffer=current_buffer?0:1;
                }
                else{
                    cout << "ERROR: Image "<< parser << " is the wrong size"<<endl;
                    i=argc;
                }
                nifti_image_free(NewImage);
            }
        }
        // *********************  POWER  *************************
        else if(strcmp(argv[i], "-pow") == 0){
            string parser=argv[++i];
            if(strtod(parser.c_str(),NULL)){
                float factor=strtof(parser.c_str(),NULL);
                for(unsigned int i=0; i<InputImage->nvox; i++)
                    bufferImages[current_buffer?0:1][i]=powf(bufferImages[current_buffer][i],factor);
                current_buffer=current_buffer?0:1;
            }
            else{
                cout << "ERROR: "<< parser << " is not a number"<<endl;
                i=argc;
            }
        }
        // *********************  square_root  *************************
        else if(strcmp(argv[i], "-sqrt") == 0){
            string parser=argv[++i];
            if(strtod(parser.c_str(),NULL)){
                for(unsigned int i=0; i<InputImage->nvox; i++)
                    bufferImages[current_buffer?0:1][i]=sqrtf(bufferImages[current_buffer][i]);
                current_buffer=current_buffer?0:1;
            }
            else{
                cout << "ERROR: "<< parser << " is not a number"<<endl;
                i=argc;
            }
        }
        // *********************  square_root  *************************
        else if(strcmp(argv[i], "-exp") == 0){
            string parser=argv[++i];
            if(strtod(parser.c_str(),NULL)){
                for(unsigned int i=0; i<InputImage->nvox; i++)
                    bufferImages[current_buffer?0:1][i]=expf(bufferImages[current_buffer][i]);
                current_buffer=current_buffer?0:1;
            }
            else{
                cout << "ERROR: "<< parser << " is not a number"<<endl;
                i=argc;
            }
        }
        // *********************  reciprocal  *************************
        else if(strcmp(argv[i], "-recip") == 0){
            string parser=argv[++i];
            if(strtod(parser.c_str(),NULL)){
                for(unsigned int i=0; i<InputImage->nvox; i++)
                    bufferImages[current_buffer?0:1][i]=1/(bufferImages[current_buffer][i]);
                current_buffer=current_buffer?0:1;
            }
            else{
                cout << "ERROR: "<< parser << " is not a number"<<endl;
                i=argc;
            }
        }
        // *********************  absolute value  *************************
        else if(strcmp(argv[i], "-abs") == 0){
            for(unsigned int i=0; i<InputImage->nvox; i++)
                bufferImages[current_buffer?0:1][i]=fabs(bufferImages[current_buffer][i]);
            current_buffer=current_buffer?0:1;

        }
        // *********************  bin value  *************************
        else if(strcmp(argv[i], "-bin") == 0){
            for(unsigned int i=0; i<InputImage->nvox; i++)
                bufferImages[current_buffer?0:1][i]=(bufferImages[current_buffer][i]>0?1.0f:0.0f);
            current_buffer=current_buffer?0:1;
        }
        // *********************  THRESHOLD below  *************************
        else if(strcmp(argv[i], "-thr") == 0){
            string parser=argv[++i];
            if(strtod(parser.c_str(),NULL)){
                double factor=strtod(parser.c_str(),NULL);
                for(unsigned int i=0; i<InputImage->nvox; i++)
                    bufferImages[current_buffer?0:1][i]=(bufferImages[current_buffer][i]>factor)?bufferImages[current_buffer][i]:0;
                current_buffer=current_buffer?0:1;
            }
            else{
                cout << "ERROR: "<< parser << " is not a number"<<endl;
                i=argc;
            }
        }
        // *********************  THRESHOLD ABOVE  *************************
        else if(strcmp(argv[i], "-uthr") == 0){
            string parser=argv[++i];
            if(strtod(parser.c_str(),NULL)){
                double factor=strtod(parser.c_str(),NULL);
                for(unsigned int i=0; i<InputImage->nvox; i++)
                    bufferImages[current_buffer?0:1][i]=(bufferImages[current_buffer][i]<factor)?bufferImages[current_buffer][i]:0;
                current_buffer=current_buffer?0:1;
            }
            else{
                cout << "ERROR: "<< parser << " is not a number"<<endl;
                i=argc;
            }
        }
        // *********************  Dilate   *************************
        else if(strcmp(argv[i], "-dil") == 0){
            string parser=argv[++i];
            if(strtod(parser.c_str(),NULL)){
                double factor=strtod(parser.c_str(),NULL);
                Dillate(bufferImages[current_buffer],(int)round(factor),CurrSize);
                for(unsigned int i=0; i<InputImage->nvox; i++)
                    bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i];
                current_buffer=current_buffer?0:1;
            }
            else{
                cout << "ERROR: "<< parser << " is not a number"<<endl;
                i=argc;
            }
        }
        // *********************  Erosion   *************************
        else if(strcmp(argv[i], "-ero") == 0){
            string parser=argv[++i];
            if(strtod(parser.c_str(),NULL)){
                double factor=strtod(parser.c_str(),NULL);
                Erosion(bufferImages[current_buffer],(int)round(factor),CurrSize);
                for(unsigned int i=0; i<InputImage->nvox; i++)
                    bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i];
                current_buffer=current_buffer?0:1;
            }
            else{
                cout << "ERROR: "<< parser << " is not a number"<<endl;
                i=argc;
            }
        }

        // *********************  GAUSSIAN SMOTHING *************************
        else if(strcmp(argv[i], "-smo") == 0){
            string parser=argv[++i];
            if(strtod(parser.c_str(),NULL)){
                float factor=strtof(parser.c_str(),NULL);

                Gaussian_Filter_4D(&bufferImages[current_buffer][0], factor, CurrSize);
                for(unsigned int i=0; i<InputImage->nvox; i++)
                    bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i];

                current_buffer=current_buffer?0:1;
            }
            else{
                cout << "ERROR: "<< parser << " is not a number"<<endl;
                i=argc;
            }
        }
        // *********************  Fill  *************************
        else if(strcmp(argv[i], "-fill") == 0){
            if(CurrSize->tsize==1){
                Close_Forground_ConnectComp<float,float>(static_cast<void*>(bufferImages[current_buffer]),static_cast<void*>(bufferImages[current_buffer?0:1]),CurrSize);
                current_buffer=current_buffer?0:1;
            }
            else{
                cout << "ERROR: Image to -fill is not 3D"<<endl;
                i=argc;
            }
        }
        // *********************  Largest Connected Component  *************************
        else if(strcmp(argv[i], "-lconcomp") == 0){
            if(CurrSize->tsize==1){
                Largest_ConnectComp<float,float>(static_cast<void*>(bufferImages[current_buffer]),static_cast<void*>(bufferImages[current_buffer?0:1]),CurrSize);
                current_buffer=current_buffer?0:1;
            }
            else{
                cout << "ERROR: Image to -lconcomp is not 3D"<<endl;
                i=argc;
            }
        }

        // *********************  Extract time point  *************************
        else if(strcmp(argv[i], "-tp") == 0){
            string parser=argv[++i];
            if(strtod(parser.c_str(),NULL) && strtod(parser.c_str(),NULL)<CurrSize->tsize){
                float factor=strtof(parser.c_str(),NULL);
                CurrSize->tsize=1;
                for(int i=0; i<CurrSize->numel; i++)
                    bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i+(int)round(factor)*CurrSize->numel];

                current_buffer=current_buffer?0:1;
            }
            else{
                cout << "ERROR: "<< parser << " is not a number"<<endl;
                i=argc;
            }
        }
        // *********************  merge time points  *************************
        else if(strcmp(argv[i], "-merge") == 0){
            string parser=argv[++i];
            if(strtod(parser.c_str(),NULL) && strtod(parser.c_str(),NULL)>0){
                int numberofTP=(int)strtof(parser.c_str(),NULL);
                int oldnumbTP=CurrSize->tsize;
                delete [] bufferImages[current_buffer?0:1];
                bufferImages[current_buffer?0:1]= new PrecisionTYPE [CurrSize->numel*(CurrSize->tsize+(int)numberofTP)];
                for(int index=0; index<(CurrSize->numel*CurrSize->tsize); index++)
                    bufferImages[current_buffer?0:1][index]=bufferImages[current_buffer][index];
                delete [] bufferImages[current_buffer];
                bufferImages[current_buffer]= new PrecisionTYPE [CurrSize->numel*(CurrSize->tsize+(int)numberofTP)];
                for(int index=0; index<(CurrSize->numel*CurrSize->tsize); index++)
                    bufferImages[current_buffer][index]=bufferImages[current_buffer?0:1][index];
                current_buffer=current_buffer?0:1;
                CurrSize->tsize=CurrSize->tsize+numberofTP;
                for(int tp=oldnumbTP; tp<CurrSize->tsize;tp++){
                    string parser_image_name=argv[++i];
                    if(parser_image_name.find(string(".nii"))>0 || parser_image_name.find(string(".img")) ||parser_image_name.find(string(".hdr"))>0){
                        nifti_image * NewImage=nifti_image_read(parser_image_name.c_str(),true);
                        if(NewImage == NULL){
                            cout<< "ERROR: When reading the image"<<parser_image_name<<endl;
                            return 1;
                        }
                        if(NewImage->nx==InputImage->nx&&NewImage->ny==InputImage->ny&&NewImage->nz==InputImage->nz&&NewImage->nt==1){
                            if(NewImage->datatype!=DT_FLOAT32){
                                seg_changeDatatype<PrecisionTYPE>(NewImage);
                            }
                            PrecisionTYPE * NewImagePtr = static_cast<PrecisionTYPE *>(NewImage->data);
                            for(int index=0; index<CurrSize->numel; index++)
                                bufferImages[current_buffer?0:1][index+tp*CurrSize->numel]=NewImagePtr[index];
                        }
                        else{
                            cout<< "ERROR: Image "<<parser_image_name<<" is not single time point or [nx,ny,nz] do not match"<<endl;
                            return 1;
                        }
                    }
                }
                current_buffer=current_buffer?0:1;
            }
            else{
                cout << "ERROR: "<< parser << " is not a number"<<endl;
                i=argc;
            }
        }
        // *********************  Get max TP  *************************
        else if(strcmp(argv[i], "-tmax") == 0){
            for(int i=0; i<CurrSize->numel; i++){
                float tmax=(float)-1.0e32;
                for(int tp=0; tp<CurrSize->tsize; tp++){
                    if(tmax<bufferImages[current_buffer][i+(int)(tp)*CurrSize->numel])
                        tmax=bufferImages[current_buffer][i+(int)(tp)*CurrSize->numel];
                }
                bufferImages[current_buffer?0:1][i]=tmax;
            }
            CurrSize->tsize=1;
            current_buffer=current_buffer?0:1;
        }
        // *********************  Get mean TP  *************************
        else if(strcmp(argv[i], "-tmean") == 0){
            for(int i=0; i<CurrSize->numel; i++){
                float tmean=0;
                for(int tp=0; tp<CurrSize->tsize; tp++){
                    tmean+=bufferImages[current_buffer][i+(int)(tp)*CurrSize->numel];
                }
                bufferImages[current_buffer?0:1][i]=tmean/CurrSize->tsize;
            }
            CurrSize->tsize=1;
            current_buffer=current_buffer?0:1;
        }
        // *********************  Get min TP  *************************
        else if(strcmp(argv[i], "-tmin") == 0){
            for(int i=0; i<CurrSize->numel; i++){
                float tmin=(float)1.0e32;
                for(int tp=0; tp<CurrSize->tsize; tp++){
                    if(tmin>bufferImages[current_buffer][i+(int)(tp)*CurrSize->numel])
                        tmin=bufferImages[current_buffer][i+(int)(tp)*CurrSize->numel];
                }
                bufferImages[current_buffer?0:1][i]=tmin;
            }
            CurrSize->tsize=1;
            current_buffer=current_buffer?0:1;
        }
        // *********************  Copy Header  *************************
        else if(strcmp(argv[i], "-hdr_copy") == 0){
            string parser=argv[++i];

            nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
            if(NewImage->datatype!=DT_FLOAT32){
                seg_changeDatatype<PrecisionTYPE>(NewImage);
            }
            PrecisionTYPE * NewImagePtr = static_cast<PrecisionTYPE *>(NewImage->data);
            if(NewImage->nx==InputImage->nx&&NewImage->ny==InputImage->ny&&NewImage->nz==InputImage->nz&&NewImage->nt==InputImage->nt&&NewImage->nu==InputImage->nu&&NewImage->nv==InputImage->nv&&NewImage->nw==InputImage->nw){
                for(unsigned int i=0; i<InputImage->nvox; i++)
                    bufferImages[current_buffer?0:1][i]=NewImagePtr[i];
                current_buffer=current_buffer?0:1;
            }
            else{
                cout << "ERROR: Image "<< parser << " is the wrong size"<<endl;
                i=argc;
            }
            nifti_image_free(NewImage);

        }
        // *********************  Get LSSD  *************************
        else if(strcmp(argv[i], "-lssd") == 0){
            string parser=argv[++i];
            nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
            if(InputImage->datatype!=DT_FLOAT32){
                seg_changeDatatype<PrecisionTYPE>(NewImage);
            }
            PrecisionTYPE * NewImagePtr = static_cast<PrecisionTYPE *>(NewImage->data);

            string parserstd=argv[++i];
            if(strtod(parserstd.c_str(),NULL)){
                if(NewImage->nt<2&&NewImage->nx==InputImage->nx&&NewImage->ny==InputImage->ny&&NewImage->nz==InputImage->nz){
                    float * NewImageMean=new float [NewImage->nx*NewImage->ny*NewImage->nz];
                    //float * NewImageStd=new float [NewImage->nx*NewImage->ny*NewImage->nz];
                    float * InputImageMean=new float [InputImage->nx*InputImage->ny*InputImage->nz];
                    //float * InputImageStd=new float [InputImage->nx*InputImage->ny*InputImage->nz];
                    float allmeanNew=0;
                    float allmeanInput=0;
                    float allstdNew=0;
                    float allstdInput=0;
                    for(int i=0; i<InputImage->nx*InputImage->ny*InputImage->nz;i++){
                        allmeanNew+=NewImagePtr[i];
                        NewImageMean[i]=NewImagePtr[i];
                        //NewImageStd[i]=NewImagePtr[i]*NewImagePtr[i];
                        allmeanInput+=bufferImages[current_buffer][i];
                        InputImageMean[i]=bufferImages[current_buffer][i];
                        //InputImageStd[i]=bufferImages[current_buffer][i]*bufferImages[current_buffer][i];
                    }
                    Gaussian_Filter_4D(NewImageMean,strtod(parserstd.c_str(),NULL)*3,CurrSize);
                    //Gaussian_Filter_4D(NewImageStd,strtod(parserstd.c_str(),NULL)*3,CurrSize);
                    Gaussian_Filter_4D(InputImageMean,strtod(parserstd.c_str(),NULL)*3,CurrSize);
                    //Gaussian_Filter_4D(InputImageStd,strtod(parserstd.c_str(),NULL)*3,CurrSize);


                    allmeanNew=allmeanNew/(InputImage->nx*InputImage->ny*InputImage->nz);
                    allmeanInput=allmeanInput/(InputImage->nx*InputImage->ny*InputImage->nz);
                    for(int i=0; i<InputImage->nx*InputImage->ny*InputImage->nz;i++){
                        allstdNew+=(NewImagePtr[i]-allmeanNew)*(NewImagePtr[i]-allmeanNew);
                        allstdInput+=(bufferImages[current_buffer][i]-allmeanInput)*(bufferImages[current_buffer][i]-allmeanInput);

                    }
                    allstdNew=allstdNew/(InputImage->nx*InputImage->ny*InputImage->nz);
                    allstdInput=allstdInput/(InputImage->nx*InputImage->ny*InputImage->nz);
                    //cout << allstdInput <<"  "<< allstdNew<<endl;
                    for(int i=0; i<InputImage->nx*InputImage->ny*InputImage->nz;i++){
                        bufferImages[current_buffer][i]=((bufferImages[current_buffer][i]-InputImageMean[i])/allstdInput)-((NewImagePtr[i]-NewImageMean[i])/allstdNew);
                        bufferImages[current_buffer?0:1][i]=(bufferImages[current_buffer][i]*bufferImages[current_buffer][i]);
                    }
                    current_buffer=current_buffer?0:1;
                    Gaussian_Filter_4D(bufferImages[current_buffer],strtod(parserstd.c_str(),NULL),CurrSize);
                    delete [] NewImageMean;
                    //delete [] NewImageStd;
                    delete [] InputImageMean;
                    //delete [] InputImageStd;
                }
                else{
                    cout << "ERROR: Image "<< parser << " is the wrong size"<<endl;
                    i=argc;
                }
            }
        }
        // *********************  Get LNCC  *************************
        else if(strcmp(argv[i], "-lncc") == 0){
            string parser=argv[++i];
            nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
            if(InputImage->datatype!=DT_FLOAT32){
                seg_changeDatatype<PrecisionTYPE>(NewImage);
            }
            PrecisionTYPE * NewImagePtr = static_cast<PrecisionTYPE *>(NewImage->data);

            string parserstd=argv[++i];
            if(strtod(parserstd.c_str(),NULL)){
                if(NewImage->nt<2&&NewImage->nx==InputImage->nx&&NewImage->ny==InputImage->ny&&NewImage->nz==InputImage->nz){
                    float * NewImageMean=new float [NewImage->nx*NewImage->ny*NewImage->nz];
                    float * NewImageStd=new float [NewImage->nx*NewImage->ny*NewImage->nz];
                    float * InputImageMean=new float [InputImage->nx*InputImage->ny*InputImage->nz];
                    float * InputImageStd=new float [InputImage->nx*InputImage->ny*InputImage->nz];
                    float allmeanNew=0;
                    float allmeanInput=0;
                    float allstdNew=0;
                    float allstdInput=0;
                    for(int i=0; i<InputImage->nx*InputImage->ny*InputImage->nz;i++){
                        allmeanNew+=NewImagePtr[i];
                        NewImageMean[i]=NewImagePtr[i];
                        NewImageStd[i]=NewImagePtr[i]*NewImagePtr[i];
                        allmeanInput+=bufferImages[current_buffer][i];
                        InputImageMean[i]=bufferImages[current_buffer][i];
                        InputImageStd[i]=bufferImages[current_buffer][i]*bufferImages[current_buffer][i];
                    }
                    allmeanNew=allmeanNew/(InputImage->nx*InputImage->ny*InputImage->nz);
                    allmeanInput=allmeanInput/(InputImage->nx*InputImage->ny*InputImage->nz);
                    for(int i=0; i<InputImage->nx*InputImage->ny*InputImage->nz;i++){
                        allstdNew+=(NewImagePtr[i]-allmeanNew)*(NewImagePtr[i]-allmeanNew);
                        allstdInput+=(bufferImages[current_buffer][i]-allmeanInput)*(bufferImages[current_buffer][i]-allmeanInput);
                        bufferImages[current_buffer][i]=NewImagePtr[i]*bufferImages[current_buffer][i];
                    }
                    allstdNew=allstdNew/(InputImage->nx*InputImage->ny*InputImage->nz);
                    allstdInput=allstdInput/(InputImage->nx*InputImage->ny*InputImage->nz);
                    //cout << allstdInput <<"  "<< allstdNew<<endl;
                    Gaussian_Filter_4D(bufferImages[current_buffer],strtod(parserstd.c_str(),NULL),CurrSize);
                    Gaussian_Filter_4D(NewImageMean,strtod(parserstd.c_str(),NULL),CurrSize);
                    Gaussian_Filter_4D(NewImageStd,strtod(parserstd.c_str(),NULL),CurrSize);
                    Gaussian_Filter_4D(InputImageMean,strtod(parserstd.c_str(),NULL),CurrSize);
                    Gaussian_Filter_4D(InputImageStd,strtod(parserstd.c_str(),NULL),CurrSize);
                    for(int i=0; i<InputImage->nx*InputImage->ny*InputImage->nz;i++){
                        NewImageStd[i]=NewImageStd[i]-NewImageMean[i]*NewImageMean[i];
                        InputImageStd[i]=InputImageStd[i]-InputImageMean[i]*InputImageMean[i];
                        bufferImages[current_buffer?0:1][i]=(bufferImages[current_buffer][i]-InputImageMean[i]*NewImageMean[i])/(sqrt(NewImageStd[i]*InputImageStd[i])+(0.01*(allstdNew+allstdInput)));
                    }
                    current_buffer=current_buffer?0:1;
                    delete [] NewImageMean;
                    delete [] NewImageStd;
                    delete [] InputImageMean;
                    delete [] InputImageStd;
                }
                else{
                    cout << "ERROR: Image "<< parser << " is the wrong size"<<endl;
                    i=argc;
                }
            }
            else{
                cout << "ERROR: "<< string() << " is not a float"<<endl;
                i=argc;
            }
            nifti_image_free(NewImage);
        }
        // *********************  output data type  *************************
        else if(strcmp(argv[i], "-odt") == 0){
            string parser=argv[++i];
            cout << parser<<" "<<parser.find("char")<<endl;
            if(parser.find("uchar")>=0){
                datatypeoutput=NIFTI_TYPE_UINT8;
            }
            else if(parser.find("ushort")>=0){
                datatypeoutput=NIFTI_TYPE_UINT16;
            }
            else if(parser.find("uint")>=0){
                datatypeoutput=NIFTI_TYPE_UINT32;
            }
            else if(parser.find("char")>=0){
                datatypeoutput=NIFTI_TYPE_INT8;
            }
            else if(parser.find("short")>=0){
                datatypeoutput=NIFTI_TYPE_INT16;
            }
            else if(parser.find("int")>=0){
                datatypeoutput=NIFTI_TYPE_INT32;
            }
            else if(parser.find("float")>=0){
                datatypeoutput=NIFTI_TYPE_FLOAT32;
            }
            else if(parser.find("double")>=0){
                datatypeoutput=NIFTI_TYPE_FLOAT64;
            }
            else{
                cout << "ERROR: Datatype "<< parser << " is unknown"<<endl;
                i=argc;
            }
        }
        else{
            cout << "Option "<< string(argv[i]) << " unkown"<<endl;
            i=argc;
            return 0;
        }

    }
    string parser=argv[argc-1];
    if(parser.find(string(".nii"))>0 || parser.find(string(".img")) ||parser.find(string(".hdr"))>0){
        // saving output
        char * filename_out=argv[argc-1];
        nifti_image * OutputImage = nifti_copy_nim_info(InputImage);
        OutputImage->datatype=datatypeoutput;
        nifti_set_filenames(OutputImage,filename_out,0,0);
        OutputImage->dim[4]=OutputImage->nt=CurrSize->tsize;
        OutputImage->dim[5]=OutputImage->nu=CurrSize->usize;
        OutputImage->dim[6]=OutputImage->nv=1;
        OutputImage->dim[7]=OutputImage->nw=1;
        OutputImage->dim[0]=(int)(OutputImage->dim[1]>0)+(int)(OutputImage->dim[2]>0)+(int)(OutputImage->dim[3]>0)+(int)(OutputImage->dim[4]>0)+(int)(OutputImage->dim[5]>0)+(int)(OutputImage->dim[6]>0)+(int)(OutputImage->dim[7]>0);
        nifti_update_dims_from_array(OutputImage);
        nifti_datatype_sizes(OutputImage->datatype,&OutputImage->nbyper,&OutputImage->swapsize);
        if(datatypeoutput==NIFTI_TYPE_UINT8){
            OutputImage->data = (void *) calloc(OutputImage->nvox, sizeof(unsigned char));
            unsigned char * OutputImagePtr = static_cast<unsigned char *>(OutputImage->data);
            for(int i=0; i<(CurrSize->numel*CurrSize->tsize); i++){OutputImagePtr[i]=(unsigned char)round(bufferImages[current_buffer][i]);}
        }
        else if(datatypeoutput==NIFTI_TYPE_UINT16){
            OutputImage->data = (void *) calloc(OutputImage->nvox, sizeof(unsigned short));
            unsigned short * OutputImagePtr = static_cast<unsigned short *>(OutputImage->data);
            for(int i=0; i<(CurrSize->numel*CurrSize->tsize); i++){OutputImagePtr[i]=(unsigned short)round(bufferImages[current_buffer][i]);}
        }
        else if(datatypeoutput==NIFTI_TYPE_UINT32){
            OutputImage->data = (void *) calloc(OutputImage->nvox, sizeof(unsigned int));
            unsigned int * OutputImagePtr = static_cast<unsigned int *>(OutputImage->data);
            for(int i=0; i<(CurrSize->numel*CurrSize->tsize); i++){OutputImagePtr[i]=(unsigned int)round(bufferImages[current_buffer][i]);}
        }
        else if(datatypeoutput==NIFTI_TYPE_INT8){
            OutputImage->data = (void *) calloc(OutputImage->nvox, sizeof(char));
            char * OutputImagePtr = static_cast<char *>(OutputImage->data);
            for(int i=0; i<(CurrSize->numel*CurrSize->tsize); i++){OutputImagePtr[i]=(char)round(bufferImages[current_buffer][i]);}
        }
        else if(datatypeoutput==NIFTI_TYPE_INT16){
            OutputImage->data = (void *) calloc(OutputImage->nvox, sizeof(short));
            short * OutputImagePtr = static_cast<short *>(OutputImage->data);
            for(int i=0; i<(CurrSize->numel*CurrSize->tsize); i++){OutputImagePtr[i]=(short)round(bufferImages[current_buffer][i]);}
        }
        else if(datatypeoutput==NIFTI_TYPE_INT32){
            OutputImage->data = (void *) calloc(OutputImage->nvox, sizeof(int));
            int * OutputImagePtr = static_cast<int *>(OutputImage->data);
            for(int i=0; i<(CurrSize->numel*CurrSize->tsize); i++){OutputImagePtr[i]=(int)round(bufferImages[current_buffer][i]);}
        }
        else if(datatypeoutput==NIFTI_TYPE_FLOAT32){
            OutputImage->data = (void *) calloc(OutputImage->nvox, sizeof(float));
            float * OutputImagePtr = static_cast<float *>(OutputImage->data);
            for(int i=0; i<(CurrSize->numel*CurrSize->tsize); i++){OutputImagePtr[i]=(float)bufferImages[current_buffer][i];}
        }
        else if(datatypeoutput==NIFTI_TYPE_FLOAT64){
            OutputImage->data = (void *) calloc(OutputImage->nvox, sizeof(double));
            double * OutputImagePtr = static_cast<double *>(OutputImage->data);
            for(int i=0; i<(CurrSize->numel*CurrSize->tsize); i++){OutputImagePtr[i]=(double)round(bufferImages[current_buffer][i]);}
        }
        nifti_image_write(OutputImage);
        nifti_image_free(OutputImage);
    }

    delete [] bufferImages[0];
    delete [] bufferImages[1];
    delete [] bufferImages;
    delete [] CurrSize;
    return 0;
}

