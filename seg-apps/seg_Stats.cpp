#include "_seg_tools.h"
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <fstream>
using namespace std;
#define SegPrecisionTYPE float


void Usage(char *exec)
{
  printf("\nUsage:\t%s <img1> <options>\n\n",exec);
  printf("\t -Dice <img2>\t\t| Calculate the Dice score between all classes in <img1> and <img2>\n");
  printf("\t -BinVol <bin> \t\t| Calculate volume of the image binarized at <bin> [sum(img>bin)*voxel_size]\n");
  printf("\t -FuzVol \t\t| Calculate the fractional content volume of a structure [sum(img)*voxel_size]\n");
  printf("\t -DiceCSV <csv_file> <img2>\t| Calculate the Dice score between all classes in <img1> and <img2> and save to csf file\n");
  printf("\t\n");
  return;
}

int main(int argc, char **argv)
{


  char * filenames[10];
  nifti_image * Images[10];
  int numbimg=0;
  if(argc<3){
      Usage(argv[0]);
      return 1;
    }


  filenames[0] = argv[1];
  Images[0]=nifti_image_read(filenames[0],true);
  if(Images[0]==NULL){
      fprintf(stderr, "This image can not be read: %s\n", filenames[0]);
      return 0;
    }

  for(int i=2;i<argc;i++){
      if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
         strcmp(argv[i], "-HELP")==0 || strcmp(argv[i], "-h")==0 ||
         strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0){
          Usage(argv[0]);
          return 0;
        }
      // **************************            ---------          *****************************
      // **************************            CALC DICE          *****************************
      // **************************            ---------          *****************************
      if(strcmp(argv[i], "-Dice") == 0 && (i+1)<argc){
          filenames[1] = argv[++i];
          nifti_image * Images[2];
          for(int i=0; i<2; i++){
              Images[i]=nifti_image_read(filenames[i],true);
              if(Images[i]==NULL){
                  fprintf(stderr, "This image can not be read: %s\n", filenames[i]);
                  return 0;
                }
            }
          for(int i=0; i<numbimg; i++){
              seg_changeDatatype<unsigned char>(Images[i]);
            }
          int  CountIMG1[1000]={0};
          unsigned char * Img1prt = static_cast<unsigned char *>(Images[0]->data);
          int  CountIMG2[1000]={0};
          unsigned char * Img2prt = static_cast<unsigned char *>(Images[1]->data);
          int  CountINTERSECT[1000]={0};
          int maxclass=0;
          for(unsigned int index=0; index<Images[0]->nvox; index++){
              CountIMG1[(int)(Img1prt[index])]++;
              maxclass=(int)(Img1prt[index])>maxclass?(int)(Img1prt[index]):maxclass;
              CountIMG2[(int)(Img2prt[index])]++;
              maxclass=(int)(Img2prt[index])>maxclass?(int)(Img2prt[index]):maxclass;
              if((int)(Img1prt[index])==(int)(Img2prt[index])){
                  CountINTERSECT[(int)(Img1prt[index])]++;
                }
            }
          float meanDice=0;
          float meanDiceSelected=0;
          int numbmeanDiceSelected=0;
          cout << "Calculating Dice from "<<maxclass << " Labels"<< endl;
          flush(cout);
          for(int curtclass=1; curtclass<=maxclass; curtclass++){
              cout<< "class["<<curtclass<<"] = "<< (float)2.0*(float)CountINTERSECT[curtclass]/((float)CountIMG1[curtclass]+(float)CountIMG2[curtclass])<<endl;
              meanDice+=(float)2.0*(float)CountINTERSECT[curtclass]/((float)CountIMG1[curtclass]+(float)CountIMG2[curtclass]);

              if(curtclass==39 ||curtclass==38 ||curtclass==43 ||curtclass==42 ||curtclass==34 ||curtclass==35 ||curtclass==40 ||curtclass==41 ){
                  numbmeanDiceSelected++;
                  meanDiceSelected+=(float)2.0*(float)CountINTERSECT[curtclass]/((float)CountIMG1[curtclass]+(float)CountIMG2[curtclass]);
                }

            }
          if(maxclass>1){
              cout<< "MEAN DICE = "<< meanDice/(maxclass-1)<<endl;
              flush(cout);
            }
          for(int i=0; i<numbimg; i++){
              nifti_image_free(Images[i]);
            }
        }
      // **************************            ---------          *****************************
      // **************************            Test Vol          *****************************
      // **************************            ---------          *****************************

      else if(strcmp(argv[i], "-Test") == 0 && (i)<argc){

          time_t start,end;
          time(&start);
          filenames[1] = argv[++i];
          nifti_image * Images[2];
          for(int i=0; i<(90); i++){
              cout << i << endl;
              Images[0]=nifti_image_read(filenames[1],true);
              if(Images[0]==NULL){
                  fprintf(stderr, "This image can not be read: %s\n", filenames[0]);
                  return 0;
                }
              short * Img1prt = static_cast<short *>(Images[0]->data);
              cout << "in  "<<Images[0]->nvox<< endl;
              for (unsigned int k=0; k<100; k++){
                  Images[1]=nifti_image_read(filenames[1],true);
                  for (unsigned int j=0; j<Images[1]->nvox; j++){
                      Img1prt[j]=i;
                    }
                  nifti_image_free(Images[1]);
                }

              cout << "out"<< endl;
              nifti_image_write(Images[0]);
              nifti_image_free(Images[0]);

            }

          time(&end);


          int minutes = (int)floorf(float(end-start)/60.0f);
          int seconds = (int)(end-start - 60*minutes);
          cout << "Finished in "<<minutes<<"min "<<seconds<<"sec"<< endl;


        }
      // **************************            ---------          *****************************
      // **************************            Fuzzy Vol          *****************************
      // **************************            ---------          *****************************

      else if(strcmp(argv[i], "-FuzVol") == 0 && (i)<argc){
          seg_changeDatatype<float>(Images[0]);
          float * Img1prt = static_cast<float *>(Images[0]->data);
          float calcvol=0;
          for(unsigned int index=0; index<Images[0]->nvox; index++){
              calcvol += Img1prt[index];
            }

          cout << calcvol<<(double)(calcvol)*(double)(Images[0]->dx)*(double)(Images[0]->dy)*(double)(Images[0]->dz)<<endl;
          flush(cout);
        }
      // **************************            ---------          *****************************
      // **************************            Bin   Vol          *****************************
      // **************************            ---------          *****************************
      else if(strcmp(argv[i], "-BinVol") == 0 && (i+1)<argc){
          float binval=atof(argv[++i]);
          seg_changeDatatype<float>(Images[0]);
          float * Img1prt = static_cast<float *>(Images[0]->data);
          float calcvol=0;
          for(unsigned int index=0; index<Images[0]->nvox; index++){
              calcvol += Img1prt[index]>binval;
            }

          cout <<(double)(calcvol)*(double)(Images[0]->dx)*(double)(Images[0]->dy)*(double)(Images[0]->dz)<<endl;
          flush(cout);
        }
      // **************************            ---------          *****************************
      // **************************            CSV  DICE          *****************************
      // **************************            ---------          *****************************
      else if(strcmp(argv[i], "-DiceCSV") == 0 && (i+2)<argc){
          filenames[2] = argv[++i];
          filenames[1] = argv[++i];
          nifti_image * Images[2];
          for(int i=0; i<2; i++){
              Images[i]=nifti_image_read(filenames[i],true);
              if(Images[i]==NULL){
                  fprintf(stderr, "This image can not be read: %s\n", filenames[i]);
                  return 0;
                }
            }
          for(int i=0; i<numbimg; i++){
              seg_changeDatatype<unsigned char>(Images[i]);
            }
          int  CountIMG1[1000]={0};
          unsigned char * Img1prt = static_cast<unsigned char *>(Images[0]->data);
          int  CountIMG2[1000]={0};
          unsigned char * Img2prt = static_cast<unsigned char *>(Images[1]->data);
          int  CountINTERSECT[1000]={0};
          int maxclass=0;
          for(unsigned int index=0; index<Images[0]->nvox; index++){
              CountIMG1[(int)(Img1prt[index])]++;
              maxclass=(int)(Img1prt[index])>maxclass?(int)(Img1prt[index]):maxclass;
              CountIMG2[(int)(Img2prt[index])]++;
              maxclass=(int)(Img2prt[index])>maxclass?(int)(Img2prt[index]):maxclass;
              if((int)(Img1prt[index])==(int)(Img2prt[index])){
                  CountINTERSECT[(int)(Img1prt[index])]++;
                }
            }
          ofstream myfile;
          myfile.open(filenames[2]);

          flush(cout);
          for(int curtclass=1; curtclass<=maxclass; curtclass++){
              myfile<< (float)2.0*(float)CountINTERSECT[curtclass]/((float)CountIMG1[curtclass]+(float)CountIMG2[curtclass]);

              if(curtclass!=maxclass){
                  myfile<<",";
                }
            }
          myfile.close();
          for(int i=0; i<numbimg; i++){
              nifti_image_free(Images[i]);
            }
        }

      else{
          printf("Err:\tParameter %s unknown or incomplete\n\n",argv[i]);
          flush(cout);
          Usage(argv[0]);
          return 1;
        }
    }
  nifti_image_free(Images[0]);
  return 0;
}
