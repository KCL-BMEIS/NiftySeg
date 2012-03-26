# pragma once

//Global includes
#include <sys/types.h>


#ifndef _WINDOWS
#include <dirent.h>
#endif

#if defined(_WIN32) && !defined(__CYGWIN__)
#include <dirent_win.h>
#include <float.h>
#include <time.h>
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

template<typename T> inline bool isinf(T value) { return std::numeric_limits<T>::has_infinity && value == std::numeric_limits<T>::infinity(); }
#ifndef isnan(_X)
#define isnan(_X) _isnan(_X)
#endif

#ifndef round(_x)
#define round(_x) floor(_x + 0.5)
#endif

inline int fabs(int _x) { return (int)fabs((float)(_x)); }

#ifndef strtof(_s, _t)
#define strtof(_s, _t) (float) strtod(_s, _t)
#endif

#endif



#include <errno.h>
#include <vector>
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

// NiftySeg Defines

#define SegPrecisionTYPE float
#define non_PV_numclass 5
#define PV_numbclass 7
#define max_numbclass 10
#define maxallowedpowerorder 6
#define redux_factor_for_bias 2
#define MaxMultispectalSize 6
#define MaxMultiLableClass 250

#define LabFusion_datatype float
#define classifier_datatype unsigned char

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
    SegPrecisionTYPE relax_factor;
    SegPrecisionTYPE relax_gauss_kernel;
    int  bias_order;
    float Bias_threshold;
    SegPrecisionTYPE  MRF_strength;
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
    bool flag_Outlierness;
    char * filename_out_outlier;
    bool flag_out_outlier;
    float OutliernessThreshold;
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
    SegPrecisionTYPE loglik;
    SegPrecisionTYPE oldloglik;
    int improv_phase;
    bool prior_relax;
    bool do_pv_modeling;
    bool pv_modeling_on;
    bool sg_delineation;
    bool out;
}FLAGS;


