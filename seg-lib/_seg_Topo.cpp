#include "_seg_Topo.h"

void ConnectComp(int * Old,
                 int * New,
                 int dimensions[3],int varin){

    // Create Counter image
    int index;
    int CCcounter=1;
    int NumElements=((int)dimensions[0]*(int)dimensions[1]*(int)dimensions[2]);

    //  **********    Foreground    ***********
    index=0;
    for(int z=0; z<((int)dimensions[2]); z++){
        for(int y=0; y<((int)dimensions[1]); y++){
            for(int x=0; x<((int)dimensions[0]); x++){
                New[index]=0;
                if(Old[index]>0){
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
                    index = z*dimensions[1]*dimensions[0]+y*dimensions[0]+x;
                    if(Old[index]>0 && CClist[New[index]]>0){
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
                        if(tempmin>0 && tempmin<CClist[New[index]]){
                            CClist[New[index]]=tempmin;
                            numbchanges++;
                        }
                    }
                }
            }
        }

        for(int index=0;index<((int)dimensions[0]*(int)dimensions[1]*(int)dimensions[2]);index++){
            if(Old[index]>0){
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

    for(index=0;index<NumElements;index++){
        if(New[index]>0 && Old[index]>0){
            Pixelcounter[(int)New[index]]++;
        }
    }

    //cout<<Pixelcounter[probarea]<<endl;
    // If lable touches the edge, errase that lable
    index=0;
    for(int iz=0;iz<dimensions[2]; iz++){
        for(int iy=0;iy<dimensions[1]; iy++){
            for(int ix=0;ix<dimensions[0]; ix++){
                if(((ix==0) || (iy==0) || (iz==0) || (ix==(dimensions[0]-1)) || (iy==(dimensions[1]-1)) || (iz==(dimensions[2]-1))) && (Old[index]>0)){
                    Pixelcounter[(int)New[index]]=0;
                }
                index++;
            }
        }
    }
    // Find Biggest Component
    for(index=0;index<CCcounter;index++){
        if(maxForground<Pixelcounter[index]){
            maxForgroundIndex=(int)index;
            maxForground=(int)Pixelcounter[index];
        }
    }
    cout << maxForground <<" "<< maxForgroundIndex<<endl;

    //Rassign to oposit class
    for(index=0;index<NumElements;index++){
        if(Old[index]>0){
            if(varin==0){
                if(New[index]==maxForgroundIndex){
                    New[index]=1;
                }
                else{New[index]=0;
                }
            }
            else{
                if(Pixelcounter[New[index]]>varin){
                    New[index]=1;
                }

                else{New[index]=0;}
            }
        }
        else{
            New[index]=0;
        }
    }
    delete [] Pixelcounter;


    return;
}


void Close_Forground_ConnectComp(int * Old,
                                 int * New,
                                 int dimensions[3]){

    // Create Counter image
    int index;
    int CCcounter=1;
    int NumElements=((int)dimensions[0]*(int)dimensions[1]*(int)dimensions[2]);

    //  **********    Foreground    ***********
    index=0;
    for(int z=0; z<((int)dimensions[2]); z++){
        for(int y=0; y<((int)dimensions[1]); y++){
            for(int x=0; x<((int)dimensions[0]); x++){
                New[index]=0;
                if(Old[index]==0){
                    Old[index]=1;
                }
                else{
                    Old[index]=0;
                }
                index++;
            }
        }
    }

    index=0;
    for(int z=0; z<((int)dimensions[2]); z++){
        for(int y=0; y<((int)dimensions[1]); y++){
            for(int x=0; x<((int)dimensions[0]); x++){
                New[index]=0;
                if(Old[index]>0){
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
                    index = z*dimensions[1]*dimensions[0]+y*dimensions[0]+x;
                    if(Old[index]>0 && CClist[New[index]]>0){
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
                        if(tempmin>0 && tempmin<CClist[New[index]]){
                            CClist[New[index]]=tempmin;
                            numbchanges++;
                        }
                    }
                }
            }
        }

        for(int index=0;index<((int)dimensions[0]*(int)dimensions[1]*(int)dimensions[2]);index++){
            if(Old[index]>0){
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

    for(index=0;index<NumElements;index++){
        if(New[index]>0 && Old[index]>0){
            Pixelcounter[(int)New[index]]++;
        }
    }

    //cout<<Pixelcounter[probarea]<<endl;
    // If lable touches the edge, errase that lable
    index=0;
    for(int iz=0;iz<dimensions[2]; iz++){
        for(int iy=0;iy<dimensions[1]; iy++){
            for(int ix=0;ix<dimensions[0]; ix++){
                if(((ix==0) || (iy==0) || (iz==0) || (ix==(dimensions[0]-1)) || (iy==(dimensions[1]-1)) || (iz==(dimensions[2]-1))) && (Old[index]>0)){
                    Pixelcounter[(int)New[index]]=0;
                }
                index++;
            }
        }
    }
    // Find Biggest Component
    for(index=0;index<CCcounter;index++){
        if(maxForground<Pixelcounter[index]){
            maxForgroundIndex=(int)index;
            maxForground=(int)Pixelcounter[index];
        }
    }


    //Rassign to oposit class
    for(index=0;index<NumElements;index++){
        if(Old[index]>0){

            if(Pixelcounter[New[index]]>0){
                New[index]=1;
            }
            else{New[index]=0;}
        }
        else{
            New[index]=1;
        }
    }
    delete [] Pixelcounter;


    return;
}

void Dillate(bool * Image,
             int kernel,
             int dimensions[3], int verbose ){

    int xyzpos[3];
    int shiftdirection[3];
    shiftdirection[0]=1;
    shiftdirection[1]=dimensions[0];
    shiftdirection[2]=dimensions[1]*dimensions[0];
    if(verbose>1){
        cout<<"Buffer allocating size:" << dimensions[1]*dimensions[0]*dimensions[2] << " - " << dimensions[0]<<" - " << dimensions[1]<<" - " << dimensions[2]<<endl;
        flush(cout);
    }
    bool * Buffer = new bool [(dimensions[1])*(dimensions[0])*(dimensions[2])];

    if(verbose>1){
        cout<<"Buffer allocated sucessfully"<<endl;
        flush(cout);
    }
    if(Buffer == NULL){
        fprintf(stderr,"* Error when alocating Buffer in void Dillate(): Not enough memory\n");
        exit(1);
    }

    if(verbose>1){
        cout<<"Dilating image now";
        flush(cout);
    }
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

    if(verbose>1){
        cout<<" - Done"<<endl;
        flush(cout);
    }
    delete [] Buffer;

    return;
}



void Dillate_const(bool * Image,
                   bool * Const,
                   int kernel,
                   int dimensions[3],
                   int direction){
    int xyzpos[3];
    int shiftdir=0;
    if(fabs(direction)==1){
        shiftdir=1;
    }
    if(fabs(direction)==2){
        shiftdir=dimensions[0];
    }
    if(fabs(direction)==3){
        shiftdir=dimensions[1]*dimensions[0];
    }

    bool * Buffer = new bool [dimensions[1]*dimensions[0]*dimensions[2]]();
    bool tmpvalue=0;

    int index=0;
    for(xyzpos[2]=0;xyzpos[2]<dimensions[2];xyzpos[2]++){
        for(xyzpos[1]=0;xyzpos[1]<dimensions[1];xyzpos[1]++){
            for(xyzpos[0]=0;xyzpos[0]<dimensions[0];xyzpos[0]++){
                tmpvalue=false;
                if(direction>0){
                    for(int shift=0; shift<=((xyzpos[(int)fabs(direction)]>=(dimensions[(int)fabs(direction)]-kernel))?(int)dimensions[(int)fabs(direction)]-xyzpos[(int)fabs(direction)]-1:kernel);shift++){
                        if(Image[index+shift*shiftdir]){
                            tmpvalue=true;
                            shift=10000;
                        }
                    }
                }
                if(direction>0){
                    for(int shift=((xyzpos[(int)fabs(direction)]<kernel)?-xyzpos[(int)fabs(direction)]:-kernel); shift<=0;shift++){
                        if(Image[index+shift*shiftdir]){
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
    for(int i=0; i<(dimensions[1]*dimensions[0]*dimensions[2]); i++){
        Image[i]=Buffer[i];
    }
    index=0;
    for(xyzpos[2]=0;xyzpos[2]<dimensions[2];xyzpos[2]++){
        for(xyzpos[1]=0;xyzpos[1]<dimensions[1];xyzpos[1]++){
            for(xyzpos[0]=0;xyzpos[0]<dimensions[0];xyzpos[0]++){
                if(Image[index]>0){
                    tmpvalue=true;
                    if(direction>0){
                        for(int shift=0; shift<=((xyzpos[(int)fabs(direction)]>=(dimensions[(int)fabs(direction)]-kernel))?(int)dimensions[(int)fabs(direction)]-xyzpos[(int)fabs(direction)]-1:kernel);shift++){
                            if(Image[index+shift*shiftdir]==0 && Const[index+shift*shiftdir]==0){
                                tmpvalue=false;
                                shift=10000;
                            }
                        }
                    }
                    if(direction>0){
                        for(int shift=((xyzpos[(int)fabs(direction)]<kernel)?-xyzpos[(int)fabs(direction)]:-kernel); shift<=0;shift++){
                            if(Image[index+shift*shiftdir]==0 && Const[index+shift*shiftdir]==0){
                                tmpvalue=false;
                                shift=10000;
                            }
                        }
                    }
                    Buffer[index]=tmpvalue;
                }
                else{
                    Buffer[index]=0;
                }
                index++;
            }
        }
    }
    for(int i=0; i<(dimensions[1]*dimensions[0]*dimensions[2]); i++){
        Image[i]=Buffer[i];
    }


    delete [] Buffer;

    return;
}
