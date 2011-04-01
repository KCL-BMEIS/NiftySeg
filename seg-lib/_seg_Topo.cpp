#include "_seg_Topo.h"

void ConnectComp(int * Old,
                 int * New,
                 int dimensions[3]){

    // Create Counter image
    int index;
    int CCcounter=1;
    int flag=1;
    int NumElements=((int)dimensions[0]*(int)dimensions[1]*(int)dimensions[2]);


    //  **********    Foreground    ***********
    for(int z=0; z<((int)dimensions[2]); z++){
        for(int y=0; y<((int)dimensions[1]); y++){
            flag=0;
            for(int x=0; x<((int)dimensions[0]); x++){
                index = (z*dimensions[1]+y)*dimensions[0]+x;
                New[index]=0;
                if(Old[index]>0){
                    New[index]=CCcounter;
                    flag=1;
                }
                if(Old[index]==0 && flag==1){
                    CCcounter++;
                    flag=0;
                }
            }
        }
    }

    int tempmin;
    int numbchanges=1;


    int *CClist = (int *) calloc(CCcounter, sizeof(int));
    for(int i=1; i<CCcounter; i++){
        CClist[i]=i;
    }

    int iter=0;
    while(numbchanges!=0){
    //while(iter<3){

        flush(cout);
        iter++;
        numbchanges=0;
        int currindex=0;
        for(int z=1; z<((int)dimensions[2]-1); z++){
            for(int y=1; y<((int)dimensions[1]-1); y++){
                for(int x=1; x<((int)dimensions[0]-1); x++){
                    index = (z*dimensions[1]+y)*dimensions[0]+x;
                    if(Old[index]==1 && New[index]>0){
                        tempmin=CClist[New[index]];

                        for(int deltaZ=-1;deltaZ<=1;deltaZ+=2){
                            currindex=index+deltaZ*dimensions[0]*dimensions[1];
                            if(Old[currindex]>0 && tempmin>CClist[New[currindex]]){tempmin=CClist[New[currindex]];}
                        }
                        for(int deltaY=-1;deltaY<=1;deltaY+=2){
                            currindex=index+deltaY*dimensions[0];
                            if(Old[currindex]>0 && tempmin>CClist[New[currindex]]){tempmin=CClist[New[currindex]];}
                        }
                        for(int deltaX=-1;deltaX<=1;deltaX+=2){
                            currindex=index+deltaX;
                            if(Old[currindex]>0 && tempmin>CClist[New[currindex]]){tempmin=CClist[New[currindex]];}
                        }
                        if(tempmin<CClist[New[index]]){
                        CClist[New[index]]=tempmin;
                        numbchanges++;
                        }
                    }
                }
            }
        }
        // Update C. Components
        for(int index=0;index<((int)dimensions[0]*(int)dimensions[1]*(int)dimensions[2]);index++){
            New[index]=CClist[New[index]];
        }
//        cout << numbchanges<< endl;
    }





    int numblab=0;
    for(index=0; index<CCcounter; index++){
        if(index==CClist[index]){
            CClist[index]=numblab;
            numblab++;
        }
    }

    for(index=0;index<NumElements;index++){
        New[index]=CClist[New[index]];
    }

    //Find lable counts
    int *Pixelcounter = (int *) calloc((int)numblab, sizeof(int));
    int maxForground=0, maxForgroundIndex=0;

    for(index=0;index<NumElements;index++){
        if((int)New[index]!=0){
            Pixelcounter[(int)New[index]]++;
        }
    }

    // If lable touches the edge, errase that lable
    index=0;
    for(int iz=0;iz<dimensions[2]; iz++){
        for(int iy=0;iy<dimensions[1]; iy++){
            for(int ix=0;ix<dimensions[0]; ix++){
                if(ix==0 || iy==0 || iz==0 || ix==dimensions[0] || iy==dimensions[1] || iz==dimensions[2]){
                    Pixelcounter[(int)New[index]]=0;
                }
                index++;
            }
        }
    }
    // Find Biggest Component
    for(index=0;index<numblab;index++){
        if(maxForground<Pixelcounter[index]){
            maxForgroundIndex=(int)index;
            maxForground=(int)Pixelcounter[index];
        }
    }

    //Rassign to oposit class
    for(index=0;index<((int)dimensions[0]*(int)dimensions[1]*(int)dimensions[2]);index++){
        if(Old[index]==1){
            if(New[index]==maxForgroundIndex){
                New[index]=1;
            }
            else{New[index]=0;}
        }
        else{
            New[index]=0;
        }
    }



    return;
}


void Dillate(bool * Image,
             int kernel,
             int dimensions[3] ){

    int xyzpos[3];
    int shiftdirection[3];
    shiftdirection[0]=1;
    shiftdirection[1]=dimensions[0];
    shiftdirection[2]=dimensions[1]*dimensions[0];
    bool * Buffer = new bool [dimensions[1]*dimensions[0]*dimensions[2]]();
    bool tmpvalue=0;
    for(int currentdirection=0;currentdirection<3;currentdirection++){ //Buffer aint each direction
        int index=0;
        for(xyzpos[2]=0;xyzpos[2]<dimensions[2];xyzpos[2]++){
            for(xyzpos[1]=0;xyzpos[1]<dimensions[1];xyzpos[1]++){
                for(xyzpos[0]=0;xyzpos[0]<dimensions[0];xyzpos[0]++){
                    tmpvalue=false;
                    for(int shift=((xyzpos[currentdirection]<kernel)?-xyzpos[currentdirection]:-kernel); shift<=((xyzpos[currentdirection]>=(dimensions[currentdirection]-kernel))?(int)dimensions[currentdirection]-xyzpos[currentdirection]-1:kernel);shift++){
                        if(Image[index+shift*shiftdirection[currentdirection]]){
                            tmpvalue=true;
                            shift=10000;
                        }

                    }
                    Buffer[index]=tmpvalue;
                    index++;
                }
            }
        }
        for(int i=0; i<(dimensions[1]*dimensions[0]*dimensions[2]); i++){
            Image[i]=Buffer[i];
        }
    }

    delete [] Buffer;

    return;
}



void Dillate_const(bool * Image,
                   bool * Const,
                   int kernel,
                   int dimensions[3] ){

    int xyzpos[3];
    int shiftdirection[3];
    shiftdirection[0]=1;
    shiftdirection[1]=dimensions[0];
    shiftdirection[2]=dimensions[1]*dimensions[0];
    bool * Buffer = new bool [dimensions[1]*dimensions[0]*dimensions[2]]();
    bool tmpvalue=0;
    for(int currentdirection=0;currentdirection<3;currentdirection++){ //Buffer aint each direction
        int index=0;
        for(xyzpos[2]=0;xyzpos[2]<dimensions[2];xyzpos[2]++){
            for(xyzpos[1]=0;xyzpos[1]<dimensions[1];xyzpos[1]++){
                for(xyzpos[0]=0;xyzpos[0]<dimensions[0];xyzpos[0]++){
                    tmpvalue=false;
                    for(int shift=((xyzpos[currentdirection]<kernel)?-xyzpos[currentdirection]:-kernel); shift<=((xyzpos[currentdirection]>=(dimensions[currentdirection]-kernel))?(int)dimensions[currentdirection]-xyzpos[currentdirection]-1:kernel);shift++){
                        if(Image[index+shift*shiftdirection[currentdirection]] && !Const[index+shift*shiftdirection[currentdirection]]){
                            tmpvalue=true;
                            shift=10000;
                        }
                    }
                    Buffer[index]=tmpvalue;
                    index++;
                }
            }
        }
        for(int i=0; i<(dimensions[1]*dimensions[0]*dimensions[2]); i++){
            Image[i]=Buffer[i];
        }
    }

    delete [] Buffer;

    return;
}
