#include <iostream>
#include <time.h>
#include "_seg_common.h"
#include "_seg_tools.h"
#include "_seg_Topo.h"

using namespace std;
#define SegPrecisionTYPE float

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
  printf("\n\t* * Operations binary 3-D images * *\n");
  printf("\t-lconcomp\t\tTake the largest connected component\n");
  printf("\t-fill\t\t\tFill holes in binary object (e.g. fill ventricle in brain mask).\n");
  printf("\t-euc\t\t\tEuclidean distance trasnform\n");
  printf("\t-geo <float/file>\tGeodesic distance according to the speed function <float/file>\n");
  printf("\n\t* * Dimensionality reduction operations: from 4-D to 3-D * *\n");
  printf("\t-tp <int>\t\tExtract time point <int>\n");
  printf("\t-tpmax\t\t\tGet the time point with the highest value (binarise 4D probabilities)\n");
  printf("\t-tmean\t\t\tMean value of all time points.\n");
  printf("\t-tmax\t\t\tMax value of all time points.\n");
  printf("\t-tmin\t\t\tMean value of all time points.\n");
  printf("\n\t* * Dimensionality increase operations: from 3-D to 4-D * *\n");
  printf("\t-merge\t<i> <d> <files>\tMerge <i> images and the working image in the <d> dimension \n");
  printf("\t-splitlab\t\tSplit the integer labels into multiple timepoints\n");
  printf("\n\t* * Image similarity: Local metrics * *\n");
  printf("\t-lncc\t<file> <std>\tLocal CC between current img and <int> on a kernel with <std>\n");
  printf("\t-lssd\t<file> <std>\tLocal SSD between current img and <int> on a kernel with <std>\n");
  //printf("\t-lmi\t<file> <std>\tLocal MI between current img and <int> on a kernel with <std>\n");
  printf("\n\t* * Sampling * *\n");
  printf("\t-subsamp2\t\tSubsample the image by 2 using NN sampling (qform and sform scaled) \n");
  printf("\n\t* * Image header operations * *\n");
  printf("\t-hdr_copy <file> \tCopy header from working image to <file> and save in <output>.\n");
  printf("\t-scl\t\t\tReset scale and slope info.\n");
  printf("\n\t* * Output * *\n");
  printf("\t-odt <datatype> \tSet output <datatype> (char, short, int, uchar, ushort, uint, float, double).\n");
  printf("\t-range\t\t\tReset the image range to the min max\n");
  printf("\t-v\t\t\tVerbose.\n");
  printf("\n\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
  return;
}

void no_memory () {
  cout << "Failed to allocate memory!\n";
  exit (1);
}

int main(int argc, char **argv)
{
  try
  {
    set_new_handler(no_memory);
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
        fprintf(stderr,"* Error when reading the input image\n");
        return 1;
      }
    if(InputImage->datatype!=NIFTI_TYPE_FLOAT32){
        seg_changeDatatype<SegPrecisionTYPE>(InputImage);
      }
    SegPrecisionTYPE * InputImagePtr = static_cast<SegPrecisionTYPE *>(InputImage->data);
    ImageSize * CurrSize = new ImageSize [1]();
    CurrSize->numel=(int)(InputImage->nx*InputImage->ny*InputImage->nz);
    CurrSize->xsize=InputImage->nx;
    CurrSize->ysize=InputImage->ny;
    CurrSize->zsize=InputImage->nz;
    CurrSize->usize=(InputImage->nu>1)?InputImage->nu:1;
    CurrSize->tsize=(InputImage->nt>1)?InputImage->nt:1;
    float Scalling=1;
    bool verbose=0;
    int datatypeoutput=NIFTI_TYPE_FLOAT32;

    SegPrecisionTYPE ** bufferImages = new SegPrecisionTYPE * [2];
    bufferImages[0] = new SegPrecisionTYPE [CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize];
    bufferImages[1] = new SegPrecisionTYPE [CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize];
    for(long i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++){
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
            if(   (strtod(parser.c_str(),NULL)!=0 && (parser.find(".nii")==string::npos ||parser.find(".img")==string::npos ||parser.find(".hdr")==string::npos ))    ||     (parser.length()==1 && parser.find("0")!=string::npos)   ){
                double multfactor=strtod(parser.c_str(),NULL);
                for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                  bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i]*multfactor;
                current_buffer=current_buffer?0:1;
              }
            else{
                nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
                NewImage->nu=(NewImage->nu>1)?NewImage->nu:1;
                NewImage->nt=(NewImage->nt>1)?NewImage->nt:1;
                if(NewImage->datatype!=DT_FLOAT32){
                    seg_changeDatatype<SegPrecisionTYPE>(NewImage);
                  }
                SegPrecisionTYPE * NewImagePtr = static_cast<SegPrecisionTYPE *>(NewImage->data);
                if(NewImage->nx==CurrSize->xsize&&NewImage->ny==CurrSize->ysize&&NewImage->nz==CurrSize->zsize&&NewImage->nt==CurrSize->tsize&&NewImage->nu==CurrSize->usize){
                    for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                      bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i]*NewImagePtr[i];
                    current_buffer=current_buffer?0:1;
                  }
                else{
                    cout << "ERROR: Image "<< parser << " is the wrong size  -  original = ( "<<CurrSize->xsize<<","
                         <<CurrSize->ysize<<","<<CurrSize->ysize<<","<<CurrSize->tsize<<","<<CurrSize->usize<<" ) New image = ( "
                        <<NewImage->nx<<","<<NewImage->ny<<","<<NewImage->nz<<","<<NewImage->nt<<","<<NewImage->nu<<" )"<<endl;
                    i=argc;
                  }
                nifti_image_free(NewImage);
              }
          }
        // *********************  ADD  *************************
        else if(strcmp(argv[i], "-add") == 0){
            string parser=argv[++i];
            if(((strtod(parser.c_str(),NULL)!=0 && (parser.find(".nii")==string::npos ||parser.find(".img")==string::npos ||parser.find(".hdr")==string::npos ))|| (parser.length()==1 && parser.find("0")!=string::npos))){
                double addfactor=strtod(parser.c_str(),NULL);
                for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                  bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i]+addfactor;
                current_buffer=current_buffer?0:1;
              }
            else{
                nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
                NewImage->nu=(NewImage->nu>1)?NewImage->nu:1;
                NewImage->nt=(NewImage->nt>1)?NewImage->nt:1;
                if(NewImage->datatype!=DT_FLOAT32){
                    seg_changeDatatype<SegPrecisionTYPE>(NewImage);
                  }
                SegPrecisionTYPE * NewImagePtr = static_cast<SegPrecisionTYPE *>(NewImage->data);
                if(NewImage->nx==CurrSize->xsize&&NewImage->ny==CurrSize->ysize&&NewImage->nz==CurrSize->zsize&&NewImage->nt==CurrSize->tsize&&NewImage->nu==CurrSize->usize){
                    for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                      bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i]+NewImagePtr[i];
                    current_buffer=current_buffer?0:1;
                  }
                else{
                    cout << "ERROR: Image "<< parser << " is the wrong size  -  original = ( "<<CurrSize->xsize<<","
                         <<CurrSize->ysize<<","<<CurrSize->ysize<<","<<CurrSize->tsize<<","<<CurrSize->usize<<" )  New image = ( "<<NewImage->nx<<","
                        <<NewImage->ny<<","<<NewImage->nz<<","<<NewImage->nt<<","<<NewImage->nu<<" )"<<endl;
                    i=argc;
                  }
                nifti_image_free(NewImage);
              }
          }
        // *********************  SUBTRACT  *************************
        else if(strcmp(argv[i], "-sub") == 0){
            string parser=argv[++i];
            if(((strtod(parser.c_str(),NULL)!=0 && (parser.find(".nii")==string::npos ||parser.find(".img")==string::npos ||parser.find(".hdr")==string::npos ))|| (parser.length()==1 && parser.find("0")!=string::npos))){
                double factor=strtod(parser.c_str(),NULL);
                for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                  bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i]-factor;
                current_buffer=current_buffer?0:1;
              }
            else{
                nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
                NewImage->nu=(NewImage->nu>1)?NewImage->nu:1;
                NewImage->nt=(NewImage->nt>1)?NewImage->nt:1;
                if(NewImage->datatype!=DT_FLOAT32){
                    seg_changeDatatype<SegPrecisionTYPE>(NewImage);
                  }
                SegPrecisionTYPE * NewImagePtr = static_cast<SegPrecisionTYPE *>(NewImage->data);
                if(NewImage->nx==CurrSize->xsize&&NewImage->ny==CurrSize->ysize&&NewImage->nz==CurrSize->zsize&&NewImage->nt==CurrSize->tsize&&NewImage->nu==CurrSize->usize){
                    for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                      bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i]-NewImagePtr[i];
                    current_buffer=current_buffer?0:1;
                  }
                else{
                    cout << "ERROR: Image "<< parser << " is the wrong size  -  original = ( "<<CurrSize->xsize<<","
                         <<CurrSize->ysize<<","<<CurrSize->ysize<<","<<CurrSize->tsize<<","<<CurrSize->usize<<" )  New image = ( "<<NewImage->nx<<","
                        <<NewImage->ny<<","<<NewImage->nz<<","<<NewImage->nt<<","<<NewImage->nu<<" )"<<endl;
                    i=argc;
                  }
                nifti_image_free(NewImage);
              }
          }
        // *********************  DIV  *************************
        else if(strcmp(argv[i], "-div") == 0){
            string parser=argv[++i];
            if(((strtod(parser.c_str(),NULL)!=0 && (parser.find(".nii")==string::npos ||parser.find(".img")==string::npos ||parser.find(".hdr")==string::npos ))|| (parser.length()==1 && parser.find("0")!=string::npos))){
                double factor=strtod(parser.c_str(),NULL);
                for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                  bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i]/factor;
                current_buffer=current_buffer?0:1;
              }
            else{
                nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
                NewImage->nu=(NewImage->nu>1)?NewImage->nu:1;
                NewImage->nt=(NewImage->nt>1)?NewImage->nt:1;
                if(NewImage->datatype!=DT_FLOAT32){
                    seg_changeDatatype<SegPrecisionTYPE>(NewImage);
                  }
                SegPrecisionTYPE * NewImagePtr = static_cast<SegPrecisionTYPE *>(NewImage->data);
                if(NewImage->nx==CurrSize->xsize&&NewImage->ny==CurrSize->ysize&&NewImage->nz==CurrSize->zsize&&NewImage->nt==CurrSize->tsize&&NewImage->nu==CurrSize->usize){
                    for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                      bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i]/NewImagePtr[i];
                    current_buffer=current_buffer?0:1;
                  }
                else{
                    cout << "ERROR: Image "<< parser << " is the wrong size  -  original = ( "<<CurrSize->xsize<<","
                         <<CurrSize->ysize<<","<<CurrSize->ysize<<","<<CurrSize->tsize<<","<<CurrSize->usize<<" )  New image = ( "<<NewImage->nx<<","
                        <<NewImage->ny<<","<<NewImage->nz<<","<<NewImage->nt<<","<<NewImage->nu<<" )"<<endl;
                    i=argc;
                  }
                nifti_image_free(NewImage);
              }
          }
        // *********************  POWER  *************************
        else if(strcmp(argv[i], "-pow") == 0){
            string parser=argv[++i];
            if(((strtod(parser.c_str(),NULL)!=0) || (parser.length()==1 && parser.find("0")!=string::npos))){
                float factor=strtof(parser.c_str(),NULL);
                for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                  bufferImages[current_buffer?0:1][i]=powf(bufferImages[current_buffer][i],factor);
                current_buffer=current_buffer?0:1;
              }
            else{
                cout << "ERROR: "<< parser << " is not a valid number"<<endl;
                i=argc;
              }
          }
        // *********************  square_root  *************************
        else if(strcmp(argv[i], "-sqrt") == 0){
            for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
              bufferImages[current_buffer?0:1][i]=sqrtf(bufferImages[current_buffer][i]);
            current_buffer=current_buffer?0:1;
          }
        // *********************  Exponential  *************************
        else if(strcmp(argv[i], "-exp") == 0){
            for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
              bufferImages[current_buffer?0:1][i]=expf(bufferImages[current_buffer][i]);
            current_buffer=current_buffer?0:1;
          }
        // *********************  reciprocal  *************************
        else if(strcmp(argv[i], "-recip") == 0){
            for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
              bufferImages[current_buffer?0:1][i]=1/(bufferImages[current_buffer][i]);
            current_buffer=current_buffer?0:1;
          }
        // *********************  absolute value  *************************
        else if(strcmp(argv[i], "-abs") == 0){
            for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
              bufferImages[current_buffer?0:1][i]=fabs(bufferImages[current_buffer][i]);
            current_buffer=current_buffer?0:1;

          }
        // *********************  bin value  *************************
        else if(strcmp(argv[i], "-bin") == 0){
            for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
              bufferImages[current_buffer?0:1][i]=(bufferImages[current_buffer][i]>0?1.0f:0.0f);
            current_buffer=current_buffer?0:1;
          }
        // *********************  THRESHOLD below  *************************
        else if(strcmp(argv[i], "-thr") == 0){
            string parser=argv[++i];
            if(((strtod(parser.c_str(),NULL)!=0 ) || (parser.length()==1 && parser.find("0")!=string::npos))){
                double factor=strtod(parser.c_str(),NULL);
                for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                  bufferImages[current_buffer?0:1][i]=(bufferImages[current_buffer][i]>factor)?bufferImages[current_buffer][i]:0;
                current_buffer=current_buffer?0:1;
              }
            else{
                cout << "ERROR: "<< parser << " is not a valid number"<<endl;
                i=argc;
              }
          }
        // *********************  THRESHOLD ABOVE  *************************
        else if(strcmp(argv[i], "-uthr") == 0){
            string parser=argv[++i];
            if(((strtod(parser.c_str(),NULL)!=0) || (parser.length()==1 && parser.find("0")!=string::npos))){
                double factor=strtod(parser.c_str(),NULL);
                for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                  bufferImages[current_buffer?0:1][i]=(bufferImages[current_buffer][i]<factor)?bufferImages[current_buffer][i]:0;
                current_buffer=current_buffer?0:1;
              }
            else{
                cout << "ERROR: "<< parser << " is not a valid number"<<endl;
                i=argc;
              }
          }
        // *********************  Dilate   *************************
        else if(strcmp(argv[i], "-dil") == 0){
            string parser=argv[++i];
            if((strtod(parser.c_str(),NULL)>0)){
                double factor=strtod(parser.c_str(),NULL);
                Dillate(bufferImages[current_buffer],(int)round(factor),CurrSize);
                for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                  bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i];
                current_buffer=current_buffer?0:1;
              }
            else{
                cout << "ERROR: "<< parser << " has to be an integer > 0"<<endl;
                i=argc;
              }
          }
        // *********************  Erosion   *************************
        else if(strcmp(argv[i], "-ero") == 0){
            string parser=argv[++i];
            if((strtod(parser.c_str(),NULL)>0 )){
                double factor=strtod(parser.c_str(),NULL);
                Erosion(bufferImages[current_buffer],(int)round(factor),CurrSize);
                for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                  bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i];
                current_buffer=current_buffer?0:1;
              }
            else{
                cout << "ERROR: "<< parser << " has to be an integer > 0"<<endl;
                i=argc;
              }
          }
        // *********************  Euclidean Distance Transform   *************************
        else if(strcmp(argv[i], "-euc") == 0){

            bool * Lable= new bool [CurrSize->numel];
            float * Speed= new float [CurrSize->numel];
            for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++){
                Lable[i]=bufferImages[current_buffer][i];
                Speed[i]=1.0f;
              }
            float * Distance = DoubleEuclideanDistance_3D(Lable,Speed,CurrSize);

            for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
              bufferImages[current_buffer?0:1][i]=Distance[i];
            current_buffer=current_buffer?0:1;
            delete [] Distance;
            delete [] Lable;
            delete [] Speed;

          }
        // *********************  Geodesic Distance Transform   *************************
        else if(strcmp(argv[i], "-geo") == 0){


            string parser=argv[++i];
            if(((strtod(parser.c_str(),NULL)!=0 && (parser.find(".nii")==string::npos ||parser.find(".img")==string::npos ||parser.find(".hdr")==string::npos ))|| (parser.length()==1 && parser.find("0")!=string::npos))){
                if(strtod(parser.c_str(),NULL)<=0){
                    cout<< "ERROR: -geo speed should be larger than zero"<<endl;
                    return 1;
                  }
                bool * Lable= new bool [CurrSize->numel];
                float * Speed= new float [CurrSize->numel];
                for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++){
                    Lable[i]=bufferImages[current_buffer][i];
                    Speed[i]=strtod(parser.c_str(),NULL);
                  }
                float * Distance = DoubleEuclideanDistance_3D(Lable,Speed,CurrSize);

                for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                  bufferImages[current_buffer?0:1][i]=Distance[i];
                current_buffer=current_buffer?0:1;
                delete [] Distance;
                delete [] Lable;
                delete [] Speed;
              }
            else{
                nifti_image * Speed=nifti_image_read(parser.c_str(),true);
                if(Speed->datatype!=DT_FLOAT32){
                    seg_changeDatatype<SegPrecisionTYPE>(Speed);
                  }
                SegPrecisionTYPE * SpeedPtr = static_cast<SegPrecisionTYPE *>(Speed->data);

                bool * Lable= new bool [CurrSize->numel];
                for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                  Lable[i]=bufferImages[current_buffer][i];

                float * Distance = DoubleEuclideanDistance_3D(Lable,SpeedPtr,CurrSize);
                for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                  bufferImages[current_buffer?0:1][i]=Distance[i];
                current_buffer=current_buffer?0:1;

                delete [] Distance;
                delete [] Lable;
                nifti_image_free(Speed);
              }
          }

        // *********************  GAUSSIAN SMOTHING *************************
        else if(strcmp(argv[i], "-smo") == 0){
            string parser=argv[++i];
            if((strtod(parser.c_str(),NULL)!=0 )){
                float factor=strtof(parser.c_str(),NULL);
                Gaussian_Filter_4D(&bufferImages[current_buffer][0], factor, CurrSize);
                for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                  bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i];

                current_buffer=current_buffer?0:1;
              }
            else{
                cout << "ERROR: "<< parser << " has to be a number > 0"<<endl;
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
        // *********************  Range  *************************
        else if(strcmp(argv[i], "-range") == 0){
            float min=1e32;
            float max=-1e32;
            for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++){
                max=bufferImages[current_buffer][i]>max?bufferImages[current_buffer][i]:max;
                min=bufferImages[current_buffer][i]<min?bufferImages[current_buffer][i]:min;
              }
            InputImage->cal_max=max;
            InputImage->cal_min=min;
          }
        // *********************  Extract time point  *************************
        else if(strcmp(argv[i], "-tp") == 0){
            string parser=argv[++i];
            if(((strtod(parser.c_str(),NULL)!=0) || (parser.length()==1 && parser.find("0")!=string::npos && parser.find("0")!=string::npos) )&& strtod(parser.c_str(),NULL)<=CurrSize->tsize ){
                float factor=strtof(parser.c_str(),NULL);
                InputImage->dim[4]=InputImage->nt=CurrSize->tsize=1;
                InputImage->dim[0]=3;
                InputImage->dim[5]=InputImage->nu=CurrSize->usize=1;
                for(int i=0; i<CurrSize->numel; i++)
                  bufferImages[current_buffer?0:1][i]=bufferImages[current_buffer][i+(int)round(factor)*CurrSize->numel];

                current_buffer=current_buffer?0:1;
              }
            else{
                cout << "ERROR: "<< parser << " is not an integer"<<endl;
                i=argc;
              }
          }

        // *********************  Split Lables  *************************
        else if(strcmp(argv[i], "-splitlab") == 0){
            int maxlab=0;
            for(int index=0; index<(CurrSize->numel*(CurrSize->tsize*CurrSize->usize)); index++)
              maxlab=(round(bufferImages[current_buffer][index])>maxlab)?(int)round(bufferImages[current_buffer][index]):maxlab;
            maxlab=maxlab+1;
            if(maxlab>0 && CurrSize->tsize<=1&& CurrSize->usize<=1){
                CurrSize->tsize=maxlab;
                CurrSize->usize=1;

                delete [] bufferImages[current_buffer?0:1];
                bufferImages[current_buffer?0:1]= new SegPrecisionTYPE [CurrSize->numel*maxlab];
                for(int index=0; index<(CurrSize->numel*maxlab); index++)
                  bufferImages[current_buffer?0:1][index]=0.0f;
                for(int index=0; index<(CurrSize->numel); index++)
                  bufferImages[current_buffer?0:1][index+(int)round(bufferImages[current_buffer][index])*CurrSize->numel]=1.0f;
                delete [] bufferImages[current_buffer];
                bufferImages[current_buffer]= new SegPrecisionTYPE [CurrSize->numel*maxlab];
                for(int index=0; index<(CurrSize->numel*maxlab); index++)
                  bufferImages[current_buffer][index]=0;
                current_buffer=current_buffer?0:1;
              }
            else{
                if(CurrSize->tsize<=1&& CurrSize->usize<=1){
                    cout << "ERROR: Working image is not 3D"<<endl;
                  }
                else{
                    cout << "ERROR: Found only "<< maxlab << " labels"<<endl;
                  }
                i=argc;
              }
          }
        // *********************  merge time points  *************************
        else if(strcmp(argv[i], "-merge") == 0){
            string parser=argv[++i];
            string parsertp=argv[++i];
            if(strtod(parser.c_str(),NULL) && (strtod(parser.c_str(),NULL)!=0 )){
                int numberofTP=(int)strtof(parser.c_str(),NULL);
                int dim=(int)strtof(parsertp.c_str(),NULL);
                int oldnumbTP=0;
                if(dim==4){
                    oldnumbTP=CurrSize->tsize;
                  }
                else if(dim==5){
                    oldnumbTP=CurrSize->usize;
                  }
                delete [] bufferImages[current_buffer?0:1];
                bufferImages[current_buffer?0:1]= new SegPrecisionTYPE [CurrSize->numel*(oldnumbTP+(int)numberofTP)];
                for(int index=0; index<(CurrSize->numel*oldnumbTP); index++)
                  bufferImages[current_buffer?0:1][index]=bufferImages[current_buffer][index];
                delete [] bufferImages[current_buffer];
                bufferImages[current_buffer]= new SegPrecisionTYPE [CurrSize->numel*(oldnumbTP+(int)numberofTP)];
                for(int index=0; index<(CurrSize->numel*oldnumbTP); index++)
                  bufferImages[current_buffer][index]=bufferImages[current_buffer?0:1][index];
                current_buffer=current_buffer?0:1;
                if(dim==4){
                    CurrSize->usize=1;
                    CurrSize->tsize=oldnumbTP+numberofTP;
                  }
                else if(dim==5){
                    CurrSize->tsize=1;
                    CurrSize->usize=oldnumbTP+numberofTP;
                  }
                for(int tp=0; tp<numberofTP;tp++){
                    string parser_image_name=argv[++i];
                    if(parser_image_name.find(string(".nii"))>0 || parser_image_name.find(string(".img")) ||parser_image_name.find(string(".hdr"))>0){
                        nifti_image * NewImage=nifti_image_read(parser_image_name.c_str(),true);
                        if(NewImage == NULL){
                            cout<< "ERROR: When reading the image"<<parser_image_name<<endl;
                            return 1;
                          }
                        if(NewImage->nx==InputImage->nx&&NewImage->ny==InputImage->ny&&NewImage->nz==InputImage->nz){
                            if(NewImage->datatype!=DT_FLOAT32){
                                seg_changeDatatype<SegPrecisionTYPE>(NewImage);
                              }
                            SegPrecisionTYPE * NewImagePtr = static_cast<SegPrecisionTYPE *>(NewImage->data);
                            for(int index=0; index<CurrSize->numel; index++)
                              bufferImages[current_buffer?0:1][index+(oldnumbTP+tp)*CurrSize->numel]=NewImagePtr[index];
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
                cout << "ERROR: "<< parser << " has to be an integer > 0"<<endl;
                i=argc;
              }
          }
        // *********************  merge time points  *************************
        else if(strcmp(argv[i], "-subsamp2") == 0){

            int newx=floor(CurrSize->xsize/2.0f);
            int newy=floor(CurrSize->ysize/2.0f);
            int newz=floor(CurrSize->zsize/2.0f);
            int newnumel=newx*newy*newz;


            delete [] bufferImages[current_buffer?0:1];
            bufferImages[current_buffer?0:1]= new SegPrecisionTYPE [newnumel*CurrSize->tsize];
            Scalling=0.5;

            for(int indexT=0; indexT<CurrSize->tsize; indexT++)
              for(int indexZ=0; indexZ<newz; indexZ++)
                for(int indexY=0; indexY<newy; indexY++)
                  for(int indexX=0; indexX<newx; indexX++)
                    bufferImages[current_buffer?0:1][indexX+indexY*newx+indexZ*newy*newx+indexT*newnumel]=bufferImages[current_buffer][indexX*2+indexY*2*CurrSize->xsize+indexZ*2*CurrSize->xsize*CurrSize->ysize+indexT*CurrSize->numel];




            delete [] bufferImages[current_buffer];
            bufferImages[current_buffer]= new SegPrecisionTYPE [newnumel*CurrSize->tsize];
            current_buffer=current_buffer?0:1;
            CurrSize->xsize=newx;
            CurrSize->ysize=newy;
            CurrSize->zsize=newz;
            CurrSize->numel=newnumel;
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
        // *********************  Get TP with maxval  *************************
        else if(strcmp(argv[i], "-tpmax") == 0){
            for(int i=0; i<CurrSize->numel; i++){
                float tmax=(float)-1.0e32;
                float tmaxindex=-1;
                for(int tp=0; tp<CurrSize->tsize; tp++){
                    if(bufferImages[current_buffer][i+(int)(tp)*CurrSize->numel]>tmax){
                        tmax=bufferImages[current_buffer][i+(int)(tp)*CurrSize->numel];
                        tmaxindex=(float)tp;
                      }
                  }
                bufferImages[current_buffer?0:1][i]=(float)tmaxindex;
              }
            CurrSize->tsize=1;
            InputImage->cal_max=CurrSize->tsize;

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
        // *********************  Reset SCL  *************************
        else if(strcmp(argv[i], "-scl") == 0){
            InputImage->scl_inter=0;
            InputImage->scl_slope=1;

         }
        // *********************  Copy Header  *************************
        else if(strcmp(argv[i], "-hdr_copy") == 0){
            string parser=argv[++i];

            nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
            if(NewImage->datatype!=DT_FLOAT32){
                seg_changeDatatype<SegPrecisionTYPE>(NewImage);
              }
            NewImage->nu=(NewImage->nu>1)?NewImage->nu:1;
            NewImage->nt=(NewImage->nt>1)?NewImage->nt:1;
            if(NewImage->nu<1)
              NewImage->dim[5]=1;
            if(NewImage->nt<1)
              NewImage->dim[4]=1;
            nifti_update_dims_from_array(NewImage);
            SegPrecisionTYPE * NewImagePtr = static_cast<SegPrecisionTYPE *>(NewImage->data);
            if(NewImage->nx==CurrSize->xsize&&NewImage->ny==CurrSize->ysize&&NewImage->nz==CurrSize->zsize&&NewImage->nt==CurrSize->tsize&&NewImage->nu==CurrSize->usize){
                for(int i=0; i<(CurrSize->xsize*CurrSize->ysize*CurrSize->zsize*CurrSize->tsize*CurrSize->usize); i++)
                  bufferImages[current_buffer?0:1][i]=NewImagePtr[i];
                current_buffer=current_buffer?0:1;
              }
            else{
                cout << "ERROR: Image "<< parser << " is the wrong size  -  original = ( "<<CurrSize->xsize<<","
                     <<CurrSize->ysize<<","<<CurrSize->zsize<<","<<CurrSize->tsize<<","<<CurrSize->usize<<" )  New image = ( "<<NewImage->nx<<","
                    <<NewImage->ny<<","<<NewImage->nz<<","<<NewImage->nt<<","<<NewImage->nu<<" )"<<endl;
                exit(1);
                i=argc;
              }
            nifti_image_free(NewImage);

          }
        // *********************  Get LSSD  *************************
        else if(strcmp(argv[i], "-lssd") == 0){
            string parser=argv[++i];
            nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
            NewImage->nu=(NewImage->nu>1)?NewImage->nu:1;
            NewImage->nt=(NewImage->nt>1)?NewImage->nt:1;
            if(NewImage->datatype!=DT_FLOAT32){
                seg_changeDatatype<SegPrecisionTYPE>(NewImage);
              }
            SegPrecisionTYPE * NewImagePtr = static_cast<SegPrecisionTYPE *>(NewImage->data);

            string parserstd=argv[++i];
            if(strtod(parserstd.c_str(),NULL)>0){
                if(NewImage->nt<2&&NewImage->nx==InputImage->nx&&NewImage->ny==InputImage->ny&&NewImage->nz==InputImage->nz){
                    float * NewImageMean=new float [NewImage->nx*NewImage->ny*NewImage->nz];
                    float * NewImageStd=new float [NewImage->nx*NewImage->ny*NewImage->nz];
                    float * InputImageMean=new float [InputImage->nx*InputImage->ny*InputImage->nz];
                    float * InputImageStd=new float [InputImage->nx*InputImage->ny*InputImage->nz];
                    float allmeanNew=0;
                    float allmeanInput=0;
                    float allstdNew=0;
                    float allstdInput=0;
                    for(int index=0; index<InputImage->nx*InputImage->ny*InputImage->nz;index++){
                        allmeanNew+=NewImagePtr[index];
                        allmeanInput+=bufferImages[current_buffer][index];
                        NewImageMean[index]=NewImagePtr[index];
                        NewImageStd[index]=NewImagePtr[index]*NewImagePtr[index];
                        InputImageMean[index]=bufferImages[current_buffer][index];
                        InputImageStd[index]=bufferImages[current_buffer][index]*bufferImages[current_buffer][index];
                      }
                    allmeanNew=allmeanNew/(InputImage->nx*InputImage->ny*InputImage->nz);
                    allmeanInput=allmeanInput/(InputImage->nx*InputImage->ny*InputImage->nz);

                    Gaussian_Filter_4D(NewImageMean,strtod(parserstd.c_str(),NULL),CurrSize);
                    Gaussian_Filter_4D(NewImageStd,strtod(parserstd.c_str(),NULL),CurrSize);
                    Gaussian_Filter_4D(InputImageMean,strtod(parserstd.c_str(),NULL),CurrSize);
                    Gaussian_Filter_4D(InputImageStd,strtod(parserstd.c_str(),NULL),CurrSize);
                    for(int index=0; index<InputImage->nx*InputImage->ny*InputImage->nz;index++){
                        allstdNew+=(NewImagePtr[index]-allmeanNew)*(NewImagePtr[index]-allmeanNew);
                        allstdInput+=(bufferImages[current_buffer][index]-allmeanInput)*(bufferImages[current_buffer][index]-allmeanInput);
                      }
                    allstdNew=allstdNew/(InputImage->nx*InputImage->ny*InputImage->nz);
                    allstdInput=allstdInput/(InputImage->nx*InputImage->ny*InputImage->nz);
                    for(int index=0; index<InputImage->nx*InputImage->ny*InputImage->nz;index++){
                        NewImageStd[index]=NewImageStd[index]-NewImageMean[index]*NewImageMean[index];
                        InputImageStd[index]=InputImageStd[index]-InputImageMean[index]*InputImageMean[index];
                        bufferImages[current_buffer?0:1][index]=(bufferImages[current_buffer][index]-InputImageMean[index])/(sqrt(InputImageStd[index]+0.01*allstdInput))-(NewImagePtr[index]-NewImageMean[index])/(sqrt(NewImageStd[index]+0.01*allstdNew));
                      }
                    Gaussian_Filter_4D(bufferImages[current_buffer?0:1],strtod(parserstd.c_str(),NULL),CurrSize);
                    for(int index=0; index<InputImage->nx*InputImage->ny*InputImage->nz;index++){
                        bufferImages[current_buffer?0:1][index]=bufferImages[current_buffer?0:1][index]*bufferImages[current_buffer?0:1][index];
                      }

                    current_buffer=current_buffer?0:1;
                    delete [] NewImageMean;
                    delete [] NewImageStd;
                    delete [] InputImageMean;
                    delete [] InputImageStd;
                  }
                else{
                    cout << "ERROR: Image "<< parser << " is the wrong size  -  original = ( "<<CurrSize->xsize<<","
                         <<CurrSize->ysize<<","<<CurrSize->ysize<<","<<CurrSize->tsize<<","<<CurrSize->usize<<" )  New image = ( "<<NewImage->nx<<","
                        <<NewImage->ny<<","<<NewImage->nz<<","<<NewImage->nt<<","<<NewImage->nu<<" )"<<endl;
                    i=argc;
                  }
              }
          }
        // *********************  Get LNCC  *************************
        else if(strcmp(argv[i], "-lncc") == 0){
            string parser=argv[++i];
            nifti_image * NewImage=nifti_image_read(parser.c_str(),true);
            NewImage->nu=(NewImage->nu>1)?NewImage->nu:1;
            NewImage->nt=(NewImage->nt>1)?NewImage->nt:1;
            if(NewImage->datatype!=DT_FLOAT32){
                seg_changeDatatype<float>(NewImage);
              }
            SegPrecisionTYPE * NewImagePtr = static_cast<SegPrecisionTYPE *>(NewImage->data);

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
                    for(int index=0; index<InputImage->nx*InputImage->ny*InputImage->nz;index++){
                        allmeanNew+=NewImagePtr[index];
                        NewImageMean[index]=NewImagePtr[index];
                        NewImageStd[index]=NewImagePtr[index]*NewImagePtr[index];
                        allmeanInput+=bufferImages[current_buffer][index];
                        InputImageMean[index]=bufferImages[current_buffer][index];
                        InputImageStd[index]=bufferImages[current_buffer][index]*bufferImages[current_buffer][index];
                      }
                    allmeanNew=allmeanNew/(InputImage->nx*InputImage->ny*InputImage->nz);
                    allmeanInput=allmeanInput/(InputImage->nx*InputImage->ny*InputImage->nz);
                    for(int index=0; index<InputImage->nx*InputImage->ny*InputImage->nz;index++){
                        allstdNew+=(NewImagePtr[index]-allmeanNew)*(NewImagePtr[index]-allmeanNew);
                        allstdInput+=(bufferImages[current_buffer][index]-allmeanInput)*(bufferImages[current_buffer][index]-allmeanInput);
                        bufferImages[current_buffer][index]=NewImagePtr[index]*bufferImages[current_buffer][index];
                      }
                    allstdNew=allstdNew/(InputImage->nx*InputImage->ny*InputImage->nz);
                    allstdInput=allstdInput/(InputImage->nx*InputImage->ny*InputImage->nz);
                    //cout << allstdInput <<"  "<< allstdNew<<endl;
                    Gaussian_Filter_4D(bufferImages[current_buffer],strtod(parserstd.c_str(),NULL),CurrSize);
                    Gaussian_Filter_4D(NewImageMean,strtod(parserstd.c_str(),NULL),CurrSize);
                    Gaussian_Filter_4D(NewImageStd,strtod(parserstd.c_str(),NULL),CurrSize);
                    Gaussian_Filter_4D(InputImageMean,strtod(parserstd.c_str(),NULL),CurrSize);
                    Gaussian_Filter_4D(InputImageStd,strtod(parserstd.c_str(),NULL),CurrSize);
                    for(int index=0; index<InputImage->nx*InputImage->ny*InputImage->nz;index++){
                        NewImageStd[index]=NewImageStd[index]-NewImageMean[index]*NewImageMean[index];
                        InputImageStd[index]=InputImageStd[index]-InputImageMean[index]*InputImageMean[index];
                        bufferImages[current_buffer?0:1][index]=(bufferImages[current_buffer][index]-InputImageMean[index]*NewImageMean[index])/(sqrt(NewImageStd[index]*InputImageStd[index])+sqrt(0.01*(allstdNew+allstdInput)));
                      }
                    current_buffer=current_buffer?0:1;
                    delete [] NewImageMean;
                    delete [] NewImageStd;
                    delete [] InputImageMean;
                    delete [] InputImageStd;
                  }
                else{
                    cout << "ERROR: Image "<< parser << " is the wrong size  -  original = ( "<<CurrSize->xsize<<","
                         <<CurrSize->ysize<<","<<CurrSize->ysize<<","<<CurrSize->tsize<<","<<CurrSize->usize<<" )  New image = ( "<<NewImage->nx<<","
                        <<NewImage->ny<<","<<NewImage->nz<<","<<NewImage->nt<<","<<NewImage->nu<<" )"<<endl;
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
        else if(strcmp(argv[i], "-v") == 0){
            verbose=1;
          }
        else if(strcmp(argv[i], "-odt") == 0){
            string parser=argv[++i];
            if(parser.find("uchar")!=string::npos){
                datatypeoutput=NIFTI_TYPE_UINT8;
              }
            else if(parser.find("ushort")!=string::npos){
                datatypeoutput=NIFTI_TYPE_UINT16;
              }
            else if(parser.find("uint")!=string::npos){
                datatypeoutput=NIFTI_TYPE_UINT32;
              }
            else if(parser.find("char")!=string::npos){
                datatypeoutput=NIFTI_TYPE_INT8;
              }
            else if(parser.find("short")!=string::npos){
                datatypeoutput=NIFTI_TYPE_INT16;
              }
            else if(parser.find("int")!=string::npos){
                datatypeoutput=NIFTI_TYPE_INT32;
              }
            else if(parser.find("float")!=string::npos){
                datatypeoutput=NIFTI_TYPE_FLOAT32;
              }
            else if(parser.find("double")!=string::npos){
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
        OutputImage->dim[1]=CurrSize->xsize;
        OutputImage->dim[2]=CurrSize->ysize;
        OutputImage->dim[3]=CurrSize->zsize;
        OutputImage->dim[4]=OutputImage->nt=CurrSize->tsize;
        OutputImage->dim[5]=OutputImage->nu=CurrSize->usize;
        OutputImage->dim[6]=OutputImage->nv=1;
        OutputImage->dim[7]=OutputImage->nw=1;
        OutputImage->dim[0]=3;
        OutputImage->dim[0]=(OutputImage->dim[4]>1?4:OutputImage->dim[0]);
        OutputImage->dim[0]=(OutputImage->dim[5]>1?5:OutputImage->dim[0]);
        OutputImage->dim[0]=(OutputImage->dim[6]>1?6:OutputImage->dim[0]);
        OutputImage->dim[0]=(OutputImage->dim[7]>1?7:OutputImage->dim[0]);

        if(Scalling!=1.0f){
            OutputImage->dx=OutputImage->dx/Scalling;
            OutputImage->dy=OutputImage->dy/Scalling;
            OutputImage->dz=OutputImage->dz/Scalling;
            OutputImage->pixdim[1]=OutputImage->dx;
            OutputImage->pixdim[2]=OutputImage->dy;
            OutputImage->pixdim[3]=OutputImage->dz;
            if(OutputImage->sform_code>0){
                OutputImage->sto_xyz.m[0][0]/=Scalling;
                OutputImage->sto_xyz.m[0][1]/=Scalling;
                OutputImage->sto_xyz.m[0][2]/=Scalling;
                OutputImage->sto_xyz.m[1][0]/=Scalling;
                OutputImage->sto_xyz.m[1][1]/=Scalling;
                OutputImage->sto_xyz.m[1][2]/=Scalling;
                OutputImage->sto_xyz.m[2][0]/=Scalling;
                OutputImage->sto_xyz.m[2][1]/=Scalling;
                OutputImage->sto_xyz.m[2][2]/=Scalling;
              }
            //OutputImage->s
          }

        if(verbose){
            cout << "Output Dim = [ ";
            for(int i=0; i<8; i++){
                cout<<(float)OutputImage->dim[i];
                if(i<7){
                    cout<<" , ";
                  }
              }
            cout<<" ] "<<endl;
            flush(cout);
          }
        nifti_update_dims_from_array(OutputImage);
        nifti_datatype_sizes(OutputImage->datatype,&OutputImage->nbyper,&OutputImage->swapsize);
        if(datatypeoutput==NIFTI_TYPE_UINT8){
            OutputImage->data = (void *) calloc(CurrSize->numel*CurrSize->tsize*CurrSize->usize, sizeof(unsigned char));
            unsigned char * OutputImagePtr = static_cast<unsigned char *>(OutputImage->data);
            for(int i=0; i<(CurrSize->numel*CurrSize->tsize*CurrSize->usize); i++){OutputImagePtr[i]=(unsigned char)round(bufferImages[current_buffer][i]);}
          }
        else if(datatypeoutput==NIFTI_TYPE_UINT16){
            OutputImage->data = (void *) calloc(OutputImage->nvox, sizeof(unsigned short));
            unsigned short * OutputImagePtr = static_cast<unsigned short *>(OutputImage->data);
            for(int i=0; i<(CurrSize->numel*CurrSize->tsize*CurrSize->usize); i++){OutputImagePtr[i]=(unsigned short)round(bufferImages[current_buffer][i]);}
          }
        else if(datatypeoutput==NIFTI_TYPE_UINT32){
            OutputImage->data = (void *) calloc(CurrSize->numel*CurrSize->tsize*CurrSize->usize, sizeof(unsigned int));
            unsigned int * OutputImagePtr = static_cast<unsigned int *>(OutputImage->data);
            for(int i=0; i<(CurrSize->numel*CurrSize->tsize*CurrSize->usize); i++){OutputImagePtr[i]=(unsigned int)round(bufferImages[current_buffer][i]);}
          }
        else if(datatypeoutput==NIFTI_TYPE_INT8){
            OutputImage->data = (void *) calloc(CurrSize->numel*CurrSize->tsize*CurrSize->usize, sizeof(char));
            char * OutputImagePtr = static_cast<char *>(OutputImage->data);
            for(int i=0; i<(CurrSize->numel*CurrSize->tsize*CurrSize->usize); i++){OutputImagePtr[i]=(char)round(bufferImages[current_buffer][i]);}
          }
        else if(datatypeoutput==NIFTI_TYPE_INT16){
            OutputImage->data = (void *) calloc(CurrSize->numel*CurrSize->tsize*CurrSize->usize, sizeof(short));
            short * OutputImagePtr = static_cast<short *>(OutputImage->data);
            for(int i=0; i<(CurrSize->numel*CurrSize->tsize*CurrSize->usize); i++){OutputImagePtr[i]=(short)round(bufferImages[current_buffer][i]);}
          }
        else if(datatypeoutput==NIFTI_TYPE_INT32){
            OutputImage->data = (void *) calloc(CurrSize->numel*CurrSize->tsize*CurrSize->usize, sizeof(int));
            int * OutputImagePtr = static_cast<int *>(OutputImage->data);
            for(int i=0; i<(CurrSize->numel*CurrSize->tsize*CurrSize->usize); i++){OutputImagePtr[i]=(int)round(bufferImages[current_buffer][i]);}
          }
        else if(datatypeoutput==NIFTI_TYPE_FLOAT32){
            OutputImage->data = (void *) calloc(CurrSize->numel*CurrSize->tsize*CurrSize->usize, sizeof(float));
            float * OutputImagePtr = static_cast<float *>(OutputImage->data);
            for(int i=0; i<(CurrSize->numel*CurrSize->tsize*CurrSize->usize); i++){OutputImagePtr[i]=(float)bufferImages[current_buffer][i];}
          }
        else if(datatypeoutput==NIFTI_TYPE_FLOAT64){
            OutputImage->data = (void *) calloc(CurrSize->numel*CurrSize->tsize*CurrSize->usize, sizeof(double));
            double * OutputImagePtr = static_cast<double *>(OutputImage->data);
            for(int i=0; i<(CurrSize->numel*CurrSize->tsize*CurrSize->usize); i++){OutputImagePtr[i]=(double)round(bufferImages[current_buffer][i]);}
          }
        nifti_image_write(OutputImage);
        nifti_image_free(OutputImage);
      }

    delete [] bufferImages[0];
    delete [] bufferImages[1];
    delete [] bufferImages;
    delete [] CurrSize;

  }

  catch(std::exception & e)
  {
    std::cerr << "Standard exception: " << e.what() << std::endl;
  }

  catch(...)
  {
    std::cerr << "Unhandled Exception: Something went wrong! Please report the error to mjorgecardoso"<<(char) 64<<"gmail.com" << std::endl;
  }
  return 0;
}

