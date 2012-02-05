#pragma once

#include "_seg_common.h"

void ConnectComp(int * Old, int * New, int dimensions[3],int varin);

template <class OldType, class NewType>
void Largest_ConnectComp(void * Old, void * New, ImageSize * Currentsize);

template <class OldType, class NewType>
void Close_Forground_ConnectComp(void * Old, void * New, ImageSize * Currentsize);

void Dillate(bool * Image, int kernel, int dimensions[3],int verbose);

void Dillate(float * Image,int kernel,ImageSize * Currentsize );

void Erosion(float * Image,int kernel,ImageSize * Currentsize );

void Dillate_const(bool * Image,
                   bool * Const,
                   int kernel,
                   int dimensions[3],
                   int direction);
