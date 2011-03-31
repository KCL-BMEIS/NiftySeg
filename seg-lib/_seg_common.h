# pragma once

//Global includes
#include <stdio.h>
#include <new>
#include <fstream>
#include <limits>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <iostream>
#include <iomanip>
#include "nifti1_io.h"

// LoAd Defines
#define PrecisionTYPE float
#define non_PV_numclass 5
#define PV_numbclass 7
#define max_numbclass 10
#define maxallowedpowerorder 6
#define redux_factor_for_bias 3
#define MaxMultispectalSize 6
#define MaxSTAPLElable 40

// Class Define
#define WMclass 0
#define GMclass 1
#define CSFclass 2
#define dGMclass 3
#define iCSFclass 4
#define WMGMpvclass 5
#define GMCSFpvclass 6

// Macros
#define colsize(I)          ((I)->nx)
#define rowsize(I)          ((I)->ny)
#define depth(I)            ((I)->nz)
#define numclass(I)         ((I)->nt)
#define numelem(I)          ((I)->nx*(I)->ny*(I)->nz)
#define numelemtotal(I)     ((I)->nvox)

using namespace std;



typedef struct{
    int xsize;
    int ysize;
    int zsize;
    int usize;
    int tsize;
    int numclass;
    int numel;
    int numelmasked;
    int numelbias;
    float rescale_max[MaxMultispectalSize];
    float rescale_min[MaxMultispectalSize];
}ImageSize;


typedef struct{
    PrecisionTYPE  relax_factor;
    PrecisionTYPE relax_gauss_kernel;
    int  bias_order;
    float Bias_threshold;
    PrecisionTYPE  MRF_strength;
    int  maxIteration;
    int  verbose_level;
    int numb_classes;
    bool flag_T1;
    bool flag_mask;
    bool flag_MRF;
    bool flag_out;
    bool flag_Bias;
    bool flag_PV_model;
    bool flag_SG_deli;
    bool flag_bc_out;
    bool flag_manual_priors;
    float * MAP_M;
    float * MAP_V;
    bool flag_MAP;
    bool aprox;
    char * filename_T1;
    char * filename_out;
    char * filename_bias;
    char * filename_mask;
    char ** filename_priors;
}SEG_PARAM;

typedef struct{
    PrecisionTYPE loglik;
    PrecisionTYPE oldloglik;
    int improv_phase;
    bool prior_relax;
    bool pv_modeling;
    bool sg_delineation;
    bool out;
}FLAGS;


