#include "_seg_tools.h"
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <fstream>
using namespace std;
#define PrecisionTYPE float


void Usage(char *exec)
{
    printf("\nUsage:\t%s <options>\n\n",exec);

    printf("\t * * * * Options * * * *\n\n");
    printf("\t -Dice <img1> <img2>\n");
    printf("\t -DiceCSV <csv_file_name> <img1> <img2>\n");
    printf("\t\n");
    return;
}

int main(int argc, char **argv)
{


    char * filenames[10];
    int numbimg=0;

    int option=0;

    if(argc==1){
        Usage(argv[0]);
        return 1;
    }
    for(int i=1;i<argc;i++){
        if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
                strcmp(argv[i], "-HELP")==0 || strcmp(argv[i], "-h")==0 ||
                strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0){
            Usage(argv[0]);
            return 0;
        }
        else if(strcmp(argv[i], "-Dice") == 0){
            filenames[0] = argv[++i];
            filenames[1] = argv[++i];
            numbimg=2;
            option=1;
        }
        else if(strcmp(argv[i], "-DiceCSV") == 0){
            filenames[2] = argv[++i];
            filenames[0] = argv[++i];
            filenames[1] = argv[++i];
            numbimg=2;
            option=2;
        }
        else{
            fprintf(stderr,"Err:\tParameter %s unknown->\n",argv[i]);
            Usage(argv[0]);
            return 1;
        }
    }


    nifti_image * Images[10];

    for(int i=0; i<numbimg; i++){
        Images[i]=nifti_image_read(filenames[i],true);
        if(Images[i]==NULL){
            fprintf(stderr, "This image can not be read: %s\n", filenames[i]);
            return 0;
        }
    }


    //CALC DICE
    if(option==1){

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
        cout<< "MEAN DICE Internal Areas in Hammer's= "<< meanDiceSelected/(numbmeanDiceSelected)<<endl;
        flush(cout);
        }
    }

    //CALC DICE and save to CSV file
    if(option==2){

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

    }




    for(int i=0; i<numbimg; i++){
        nifti_image_free(Images[i]);
    }
    return 0;
}
