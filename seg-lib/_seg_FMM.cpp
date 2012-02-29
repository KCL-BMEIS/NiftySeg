# include "_seg_FMM.h"

float * DoubleEuclideanDistance_3D(bool *LablePtr, float * speedptr,
                                   ImageSize * CurrSizes){

  int i, NumElements, centre_index;
  NumElements=CurrSizes->xsize*CurrSizes->ysize*CurrSizes->zsize;
  int MaxGeoTime=100;
  i=0;
  int neighbour6[6]={0};
  neighbour6[0]=-1;
  neighbour6[1]=1;
  neighbour6[2]=-1*CurrSizes->xsize;
  neighbour6[3]=CurrSizes->xsize;
  neighbour6[4]=-1*CurrSizes->xsize*CurrSizes->ysize;
  neighbour6[5]=CurrSizes->xsize*CurrSizes->ysize;
  float speed=1;

  typedef std::multimap <SegPrecisionTYPE, int> MapLengthType;
  typedef std::pair <SegPrecisionTYPE, int> PairType;


  MapLengthType MapLength;
  MapLengthType::iterator first_element=MapLength.begin();

  char * Labels = (char *) calloc(NumElements, sizeof(char));
  int mycounter=0;
  bool *Border= (bool *) calloc(NumElements, sizeof(bool));
  float *GeoTime= (float *) calloc(NumElements, sizeof(float));

  for(centre_index=0;centre_index<NumElements;centre_index++){
      if(speedptr!=NULL){
          if(speedptr[centre_index]<=0.000001)
            speedptr[centre_index]=speedptr[centre_index]+0.000001;
        }
      if(LablePtr[centre_index]>0.5){
          Labels[centre_index]=0;
          GeoTime[centre_index]=0;
        }
      else{
          mycounter++;
          Labels[centre_index]=2;
          GeoTime[centre_index]=MaxGeoTime;
        }
    }


  centre_index=0;
  for(int zpos=0; zpos<CurrSizes->zsize; zpos++){
      for(int ypos=0; ypos<CurrSizes->ysize; ypos++){
          for(int xpos=0; xpos<CurrSizes->xsize; xpos++){
              Border[centre_index]=0;
              if(xpos==0 || xpos==(CurrSizes->xsize-1) || ypos==0 || ypos==(CurrSizes->ysize-1) || zpos==0 || zpos==(CurrSizes->zsize-1)){
                  Border[centre_index]=1;
                }
              centre_index++;
            }
        }
    }


  bool flagadd=false;
  int index_neighbour;
  for(centre_index=0;centre_index<NumElements;centre_index++){
      if(!Border[centre_index] && Labels[centre_index]==2){
          flagadd=false;
          for(i=0; i<6; i++){
              index_neighbour=centre_index+neighbour6[i];
              if(!flagadd && index_neighbour>=0 && Labels[index_neighbour]==0){
                  flagadd=true;
                }
            }
          if(flagadd){
              Labels[centre_index]=1;
              if(speedptr!=NULL)speed=speedptr[centre_index];
              MapLength.insert(PairType(speed*0.5, centre_index));
              GeoTime[centre_index]=speed*0.5;
            }
        }
    }

  SegPrecisionTYPE oldGeoTime, newGeoTime;
  while(MapLength.empty()==false){
      first_element=MapLength.begin();
      centre_index=first_element->second;
      if(first_element->first==GeoTime[centre_index]){
          Labels[centre_index]=0;
          if(!Border[centre_index]){
              for(i=0; i<6; i++){
                  index_neighbour=centre_index+neighbour6[i];
                  if(index_neighbour>=0){
                      if (Labels[index_neighbour]>0 && !Border[index_neighbour]){
                          if(speedptr!=NULL)speed=speedptr[index_neighbour];
                          oldGeoTime=GeoTime[index_neighbour];
                          newGeoTime=CalcGeoTime_long(index_neighbour, GeoTime, speed, neighbour6, MaxGeoTime);
                          if(newGeoTime<oldGeoTime && newGeoTime<MaxGeoTime){
                              GeoTime[index_neighbour]=newGeoTime;
                              MapLength.insert(PairType(newGeoTime, index_neighbour));
                              Labels[index_neighbour]=1;
                            }
                        }
                    }
                }
            }
        }
      MapLength.erase(first_element);
    }


  centre_index=0;
  for(int zpos=0; zpos<CurrSizes->zsize; zpos++){
      for(int ypos=0; ypos<CurrSizes->ysize; ypos++){
          for(int xpos=0; xpos<CurrSizes->xsize; xpos++){
              if(xpos==0 || xpos==(CurrSizes->xsize-1) || ypos==0 || ypos==(CurrSizes->ysize-1) || zpos==0 || zpos==(CurrSizes->zsize-1)){
                  GeoTime[centre_index]=MaxGeoTime;
                }
              centre_index++;
            }
        }
    }


  float *GeoTime2= (float *) calloc(NumElements, sizeof(float));
  for(centre_index=0;centre_index<NumElements;centre_index++){
      if(LablePtr[centre_index]<0.5){
          Labels[centre_index]=0;
          GeoTime2[centre_index]=0;
        }
      else{
          mycounter++;
          Labels[centre_index]=2;
          GeoTime2[centre_index]=MaxGeoTime;
        }
    }

  for(centre_index=0;centre_index<NumElements;centre_index++){
      if(!Border[centre_index] && Labels[centre_index]==2){
          flagadd=false;
          for(i=0; i<6; i++){
              index_neighbour=centre_index+neighbour6[i];
              if(!flagadd && index_neighbour>=0 && Labels[index_neighbour]==0){
                  flagadd=true;
                }
            }
          if(flagadd){

              Labels[centre_index]=1;
              if(speedptr!=NULL)speed=speedptr[centre_index];

              MapLength.insert(PairType(speed*0.5, centre_index));
              GeoTime2[centre_index]=speed*0.5;
            }
        }
    }

  while(MapLength.empty()==false){
      first_element=MapLength.begin();
      centre_index=first_element->second;
      if(first_element->first==GeoTime2[centre_index]){
          Labels[centre_index]=0;
          if(!Border[centre_index]){
              for(i=0; i<6; i++){
                  index_neighbour=centre_index+neighbour6[i];
                  if(index_neighbour>=0){
                      if (Labels[index_neighbour]>0 && !Border[index_neighbour]){
                          oldGeoTime=GeoTime2[index_neighbour];
                          if(speedptr!=NULL)speed=speedptr[index_neighbour];
                          newGeoTime=CalcGeoTime_long(index_neighbour, GeoTime2, speed, neighbour6, MaxGeoTime);
                          if(newGeoTime<oldGeoTime && newGeoTime<MaxGeoTime){
                              GeoTime2[index_neighbour]=newGeoTime;
                              MapLength.insert(PairType(newGeoTime, index_neighbour));
                              Labels[index_neighbour]=1;
                            }
                        }
                    }
                }
            }
        }
      MapLength.erase(first_element);
    }



  for(centre_index=0;centre_index<NumElements;centre_index++){
      GeoTime[centre_index]=LablePtr[centre_index]<0.5?(-1*GeoTime[centre_index]):(GeoTime2[centre_index]);
    }

  centre_index=0;
  for(int zpos=0; zpos<CurrSizes->zsize; zpos++){
      for(int ypos=0; ypos<CurrSizes->ysize; ypos++){
          for(int xpos=0; xpos<CurrSizes->xsize; xpos++){
              if(Border[centre_index]){
                  int xshift=(xpos==0)?1:(xpos==(CurrSizes->xsize-1)?-1:0);
                  int yshift=(ypos==0)?1:(ypos==(CurrSizes->ysize-1)?-1:0);
                  int zshift=(zpos==0)?1:(zpos==(CurrSizes->zsize-1)?-1:0);
                  GeoTime[centre_index]=GeoTime[(centre_index+neighbour6[1]*xshift+neighbour6[3]*yshift+neighbour6[5]*zshift)];
                }
              centre_index++;
            }
        }
    }

  free(Labels);
  free(Border);
  free(GeoTime2);
  return GeoTime;
}

void FMM(bool *Seeds,
         SegPrecisionTYPE *SpeedI,
         SegPrecisionTYPE *GeoTime,
         SegPrecisionTYPE MaxGeoTime,
         int * L2S,
         int * S2L,
         ImageSize * CurrSizes){

  int i, NumElements, centre_short_index;
  int index_long;
  int index_neighbour_short;
  NumElements=CurrSizes->numelmasked;
  i=0;
  int neighbour6[6]={0};
  neighbour6[0]=-1;
  neighbour6[1]=1;
  neighbour6[2]=-1*CurrSizes->xsize;
  neighbour6[3]=CurrSizes->xsize;
  neighbour6[4]=-1*CurrSizes->xsize*CurrSizes->ysize;
  neighbour6[5]=CurrSizes->xsize*CurrSizes->ysize;



  typedef std::multimap <SegPrecisionTYPE, int> MapLengthType;
  typedef std::pair <SegPrecisionTYPE, int> PairType;


  MapLengthType MapLength;
  MapLengthType::iterator first_element=MapLength.begin();

  char * Labels = (char *) calloc(NumElements, sizeof(char));
  int mycounter=0;
  bool *Border= (bool *) calloc(NumElements, sizeof(bool));

  for(centre_short_index=0;centre_short_index<NumElements;centre_short_index++){
      if(Seeds[centre_short_index]){
          Labels[centre_short_index]=0;
          GeoTime[centre_short_index]=0;
        }
      else{
          mycounter++;
          Labels[centre_short_index]=2;
          GeoTime[centre_short_index]=MaxGeoTime;
        }
    }

  bool curr_res=false;
  for(centre_short_index=0; centre_short_index<NumElements; centre_short_index++){
      index_long=S2L[centre_short_index];
      curr_res=false;
      for(i=0; i<6; i++){
          if(!curr_res){
              curr_res=(L2S[index_long+neighbour6[i]]<=0)?true:false;
            }
        }
      Border[centre_short_index]=curr_res;
    }


  bool flagadd=false;
  for(centre_short_index=0;centre_short_index<NumElements;centre_short_index++){
      if(!Border[centre_short_index] && Labels[centre_short_index]==2){
          index_long=S2L[centre_short_index];
          flagadd=false;
          for(i=0; i<6; i++){
              index_neighbour_short=L2S[index_long+neighbour6[i]];
              if(!flagadd && index_neighbour_short>=0 && Labels[index_neighbour_short]==0){
                  flagadd=true;
                }
            }
          if(flagadd){
              Labels[centre_short_index]=1;
              MapLength.insert(PairType((1+SpeedI[centre_short_index])*0.5, centre_short_index));
              GeoTime[centre_short_index]=0.5*(1+SpeedI[centre_short_index]);
            }
        }
    }

  SegPrecisionTYPE oldGeoTime, newGeoTime;
  while(MapLength.empty()==false){
      first_element=MapLength.begin();
      centre_short_index=first_element->second;
      if(first_element->first==GeoTime[centre_short_index]){
          Labels[centre_short_index]=0;
          if(!Border[centre_short_index]){
              index_long=S2L[centre_short_index];
              for(i=0; i<6; i++){
                  index_neighbour_short=L2S[index_long+neighbour6[i]];
                  if(index_neighbour_short>=0){
                      if (Labels[index_neighbour_short]>0 && !Border[index_neighbour_short]){
                          oldGeoTime=GeoTime[index_neighbour_short];
                          newGeoTime=CalcGeoTime(index_long+neighbour6[i], GeoTime, SpeedI, neighbour6, L2S, MaxGeoTime);
                          if(newGeoTime<oldGeoTime && newGeoTime<MaxGeoTime){
                              GeoTime[index_neighbour_short]=newGeoTime;
                              MapLength.insert(PairType(newGeoTime, index_neighbour_short));
                              Labels[index_neighbour_short]=1;
                            }
                        }
                    }
                }
            }
        }
      MapLength.erase(first_element);
    }

  for(centre_short_index=0;centre_short_index<NumElements;centre_short_index++){
      if(Border[centre_short_index]){
          index_long=S2L[centre_short_index];
          float curr_max=0;
          for(i=0; i<6; i++){
              if(!Border[L2S[index_long+neighbour6[i]]] && GeoTime[L2S[index_long+neighbour6[i]]]>curr_max){
                  curr_max=GeoTime[L2S[index_long+neighbour6[i]]];
                }
            }
          GeoTime[centre_short_index]=curr_max;
        }
    }
  free(Labels);
  free(Border);
}

// **************************************************************************************
// ********************************** OTHER FUNCTIONS  **********************************
// **************************************************************************************

SegPrecisionTYPE CalcGeoTime(int index,
                             SegPrecisionTYPE *GeoTime,
                             SegPrecisionTYPE * SpeedI,
                             int * neighbour,
                             int * L2S,
                             SegPrecisionTYPE Max){

  SegPrecisionTYPE average=0.0f;
  int index_neighbour_short;
  SegPrecisionTYPE counter=0.0f;
  SegPrecisionTYPE a=0.0f, c=0.0f;

  int neighbour_index[2]={0};
  for(int i=0; i<6; i+=2){
      neighbour_index[0]=L2S[index+neighbour[i]];
      neighbour_index[1]=L2S[index+neighbour[i+1]];
      if((GeoTime[neighbour_index[0]]<Max || GeoTime[neighbour_index[1]]<Max)){
          index_neighbour_short=(GeoTime[neighbour_index[0]]<GeoTime[neighbour_index[1]])?neighbour_index[0]:neighbour_index[1];
          average+=GeoTime[index_neighbour_short];
          counter+=1.0f;

        }
    }
  average=average/counter;

  for(int i=0; i<6; i+=2){
      neighbour_index[0]=L2S[index+neighbour[i]];
      neighbour_index[1]=L2S[index+neighbour[i+1]];
      if((GeoTime[neighbour_index[0]]<Max || GeoTime[neighbour_index[1]]<Max)){
          a+=1.0f;
          index_neighbour_short=(GeoTime[neighbour_index[0]]<GeoTime[neighbour_index[1]])?neighbour_index[0]:neighbour_index[1];
          c += powf(GeoTime[index_neighbour_short]-average,2);
        }
    }

  c -= 1+SpeedI[L2S[index]];
  if((-4.0f*a*c)>0){
      return sqrtf(-4.0f*a*c)/(2.0f*a)+average;
    }
  else{
      return average;
    }
}



SegPrecisionTYPE CalcGeoTime_long(int index,
                                  SegPrecisionTYPE *GeoTime,
                                  SegPrecisionTYPE SpeedI,
                                  int * neighbour,
                                  SegPrecisionTYPE Max){

  SegPrecisionTYPE average=0.0f;
  int neighbour_index_true;
  SegPrecisionTYPE counter=0.0f;
  SegPrecisionTYPE a=0.0f, c=0.0f;

  int neighbour_index[2]={0};
  for(int i=0; i<6; i+=2){
      neighbour_index[0]=index+neighbour[i];
      neighbour_index[1]=index+neighbour[i+1];
      if(((GeoTime[neighbour_index[0]]<Max &&GeoTime[neighbour_index[0]]>0) || (GeoTime[neighbour_index[1]]<Max &&GeoTime[neighbour_index[1]]>0))){
          neighbour_index_true=(GeoTime[neighbour_index[0]]<GeoTime[neighbour_index[1]])?neighbour_index[0]:neighbour_index[1];
          average+=GeoTime[neighbour_index_true];
          counter+=1.0f;

        }
    }
  average=average/counter;

  for(int i=0; i<6; i+=2){
      neighbour_index[0]=index+neighbour[i];
      neighbour_index[1]=index+neighbour[i+1];
      if(((GeoTime[neighbour_index[0]]<Max &&GeoTime[neighbour_index[0]]>0) || (GeoTime[neighbour_index[1]]<Max &&GeoTime[neighbour_index[1]]>0))){
          a+=1.0f;
          neighbour_index_true=(GeoTime[neighbour_index[0]]<GeoTime[neighbour_index[1]])?neighbour_index[0]:neighbour_index[1];
          c += powf(GeoTime[neighbour_index_true]-average,2);
        }
    }

  c -= SpeedI;
  if((-4.0f*a*c)>0){
      return sqrtf(-4.0f*a*c)/(2.0f*a)+average;
    }
  else{
      return average;
    }
}

void TransformGeoTime(SegPrecisionTYPE *GeoTime,
                      SegPrecisionTYPE MaxGeoTime,
                      int * L2S,
                      int * S2L,
                      ImageSize * CurrSizes){



  float trace_xgrad,trace_ygrad,trace_zgrad,normgrad,trace,sqrtgrad;
  int neighbour6[6];
  neighbour6[0]=-1;
  neighbour6[1]=1;
  neighbour6[2]=-1*CurrSizes->xsize;
  neighbour6[3]=CurrSizes->xsize;
  neighbour6[4]=-1*CurrSizes->xsize*CurrSizes->ysize;
  neighbour6[5]=CurrSizes->xsize*CurrSizes->ysize;

  float *Buffer= (float *) calloc(CurrSizes->numelmasked, sizeof(float));
  float *Gradxyz= (float *) calloc(CurrSizes->numelmasked*3, sizeof(float));
  bool isBroder=false;
  for(int centre_short_index=0; centre_short_index<CurrSizes->numelmasked; centre_short_index++){
      int index_long=S2L[centre_short_index];
      isBroder=false;
      for(int i=0; i<6; i++){
          if(!isBroder){
              isBroder=(L2S[index_long+neighbour6[i]]<=0)?true:false;
            }
        }
      if(!isBroder){
          /*Gradxyz[centre_short_index*3]=GeoTime[L2S[index_long+neighbour6[0]]]-GeoTime[L2S[index_long+neighbour6[1]]];

            Gradxyz[centre_short_index*3+1]=GeoTime[L2S[index_long+neighbour6[2]]]-GeoTime[L2S[index_long+neighbour6[3]]];

            Gradxyz[centre_short_index*3+2]=GeoTime[L2S[index_long+neighbour6[4]]]-GeoTime[L2S[index_long+neighbour6[5]]];*/


          Gradxyz[centre_short_index*3]=GeoTime[L2S[index_long+neighbour6[0]]]-GeoTime[L2S[index_long+neighbour6[1]]]+
              0.3*GeoTime[L2S[index_long+neighbour6[0]+neighbour6[2]]]-0.3*GeoTime[L2S[index_long+neighbour6[1]+neighbour6[3]]]+
              0.3*GeoTime[L2S[index_long+neighbour6[0]-neighbour6[2]]]-0.3*GeoTime[L2S[index_long+neighbour6[1]-neighbour6[3]]]+
              0.3*GeoTime[L2S[index_long+neighbour6[0]+neighbour6[4]]]-0.3*GeoTime[L2S[index_long+neighbour6[1]+neighbour6[5]]]+
              0.3*GeoTime[L2S[index_long+neighbour6[0]-neighbour6[4]]]-0.3*GeoTime[L2S[index_long+neighbour6[1]-neighbour6[5]]];

          Gradxyz[centre_short_index*3+1]=GeoTime[L2S[index_long+neighbour6[2]]]-GeoTime[L2S[index_long+neighbour6[3]]]+
              0.3*GeoTime[L2S[index_long+neighbour6[2]+neighbour6[0]]]-0.3*GeoTime[L2S[index_long+neighbour6[3]+neighbour6[1]]]+
              0.3*GeoTime[L2S[index_long+neighbour6[2]-neighbour6[0]]]-0.3*GeoTime[L2S[index_long+neighbour6[3]-neighbour6[1]]]+
              0.3*GeoTime[L2S[index_long+neighbour6[2]+neighbour6[4]]]-0.3*GeoTime[L2S[index_long+neighbour6[3]+neighbour6[5]]]+
              0.3*GeoTime[L2S[index_long+neighbour6[2]-neighbour6[4]]]-0.3*GeoTime[L2S[index_long+neighbour6[3]-neighbour6[5]]];

          Gradxyz[centre_short_index*3+2]=GeoTime[L2S[index_long+neighbour6[4]]]-GeoTime[L2S[index_long+neighbour6[5]]]+
              0.3*GeoTime[L2S[index_long+neighbour6[4]+neighbour6[2]]]-0.3*GeoTime[L2S[index_long+neighbour6[5]+neighbour6[3]]]+
              0.3*GeoTime[L2S[index_long+neighbour6[4]-neighbour6[2]]]-0.3*GeoTime[L2S[index_long+neighbour6[5]-neighbour6[3]]]+
              0.3*GeoTime[L2S[index_long+neighbour6[4]+neighbour6[0]]]-0.3*GeoTime[L2S[index_long+neighbour6[5]+neighbour6[1]]]+
              0.3*GeoTime[L2S[index_long+neighbour6[4]-neighbour6[0]]]-0.3*GeoTime[L2S[index_long+neighbour6[5]-neighbour6[1]]];
        }
    }

  for(int centre_short_index=0; centre_short_index<CurrSizes->numelmasked; centre_short_index++){
      int index_long=S2L[centre_short_index];
      isBroder=false;
      for(int i=0; i<6; i++){
          if(!isBroder){
              isBroder=(L2S[index_long+neighbour6[i]]<=0)?true:false;
            }
        }
      if(!isBroder){
          sqrtgrad=sqrt(Gradxyz[centre_short_index*3]*Gradxyz[centre_short_index*3]+Gradxyz[centre_short_index*3+1]*Gradxyz[centre_short_index*3+1]+Gradxyz[centre_short_index*3+2]*Gradxyz[centre_short_index*3+2]);
          normgrad=(sqrtgrad/7)<1?(1-(sqrtgrad/7)):0;

          trace_xgrad = Gradxyz[L2S[index_long+neighbour6[0]]*3]   -  Gradxyz[L2S[index_long+neighbour6[1]]*3];
          trace_ygrad = Gradxyz[L2S[index_long+neighbour6[2]]*3+1] -  Gradxyz[L2S[index_long+neighbour6[3]]*3+1];
          trace_zgrad = Gradxyz[L2S[index_long+neighbour6[4]]*3+2] -  Gradxyz[L2S[index_long+neighbour6[5]]*3+2];
          trace = (trace_xgrad+trace_ygrad+trace_zgrad);
          trace=trace<0?fabs(trace/7):0;

          //Buffer[centre_short_index]=sqrtgrad;
          Buffer[centre_short_index]=((trace*normgrad)>1)?1:((trace*normgrad)>0.2?(trace*normgrad):0);
          Buffer[centre_short_index]=pow(Buffer[centre_short_index],4);
        }
    }


  for(int centre_short_index=0; centre_short_index<CurrSizes->numelmasked; centre_short_index++){
      GeoTime[centre_short_index]=(GeoTime[centre_short_index]>2 && GeoTime[centre_short_index]<10)?Buffer[centre_short_index]:0;
      //GeoTime[centre_short_index]=Buffer[centre_short_index];
    }

  free(Gradxyz);
  free(Buffer);
}


