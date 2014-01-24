# pragma once

//Global includes
#include <sys/types.h>

#if (defined(_WIN32) || defined(_WINDOWS)) && !defined(__CYGWIN__)
#include <dirent_win.h>
#include <float.h>
#include <time.h>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif



#if (defined(_WIN32) || defined(_WINDOWS)) && !defined(__CYGWIN__)
#define SEP "\\0"
#else
#define SEP "/\0"
#endif

template<typename T> inline bool isinf(T value)
{
    return std::numeric_limits<T>::has_infinity && value == std::numeric_limits<T>::infinity();
}
#ifndef isnan(_X)
#define isnan(_X) _isnan(_X)
#endif

#ifndef round(_x)
#define round(_x) floor(_x + 0.5)
#endif

inline int fabs(int _x)
{
    return (int)fabs((float)(_x));
}

#ifndef strtof(_s, _t)
#define strtof(_s, _t) (float) strtod(_s, _t)
#endif

#else  //IF NOT ON WINDOWS
#include <dirent.h>
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
#include <new>
#include <exception>

// NiftySeg Defines

#define SegPrecisionTYPE float
#define non_PV_numclass 5
#define PV_numbclass 20
#define max_numbclass 12
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



typedef struct
{
    long xsize;
    long ysize;
    long zsize;
    long usize;
    long tsize;
    long numclass;
    long numel;
    long numelmasked;
    long numelbias;
    float rescale_max[MaxMultispectalSize];
    float rescale_min[MaxMultispectalSize];
} ImageSize;


typedef struct SegParams
{

    SegPrecisionTYPE relax_factor;
    SegPrecisionTYPE relax_gauss_kernel;
    int  bias_order;
    float Bias_threshold;
    SegPrecisionTYPE  MRF_strength;
    int  maxIteration;
    int  minIteration;
    int  verbose_level;
    long numb_classes;
    bool flag_T1;
    bool flag_mask;
    bool flag_MRF;
    bool flag_out;
    bool flag_Bias;
    bool flag_PV_model;
    bool flag_SG_deli;
    bool flag_bc_out;
    bool flag_manual_priors;
    bool flag_priors4D;
    bool flag_Outlierness;
    char * filename_out_outlier;
    bool flag_out_outlier;
    float OutliernessThreshold;
    float OutliernessRatio;
    float * MAP_M;
    float * MAP_V;
    bool flag_MAP;
    bool aprox;
    char * filename_T1;
    char * filename_out;
    char * filename_bias;
    char * filename_mask;
    char ** filename_priors;

    SegParams( void )
    {
        relax_factor = 0;
        relax_gauss_kernel = 0;
        bias_order = 0;
        Bias_threshold = 0;
        MRF_strength = 0;
        maxIteration = 0;
        minIteration = 0;
        verbose_level = 0;
        numb_classes = 0;
        flag_T1 = 0;
        flag_mask = 0;
        flag_MRF = 0;
        flag_out = 0;
        flag_Bias = 0;
        flag_PV_model = 0;
        flag_SG_deli = 0;
        flag_bc_out = 0;
        flag_manual_priors = 0;
        flag_priors4D = 0;
        flag_Outlierness = 0;
        filename_out_outlier = 0;
        flag_out_outlier = 0;
        OutliernessThreshold = 0;
        OutliernessRatio = 0;
        MAP_M = 0;
        MAP_V = 0;
        flag_MAP = 0;
        aprox = 0;
        filename_T1 = 0;
        filename_out = 0;
        filename_bias = 0;
        filename_mask = 0;
        filename_priors = 0;
    }

    void Print( ostream &sout )
    {
        sout << "PARAMS: verbose_level: " << verbose_level << endl
             << "PARAMS: flag_T1: " << flag_T1 << endl
             << "PARAMS: filename_T1: " << (char *) ( ( filename_T1 ) ? filename_T1 : "0" ) << endl
             << "PARAMS: flag_out: " << flag_out << endl
             << "PARAMS: filename_out: " << (char *) ( ( filename_out ) ? filename_out : "0" ) << endl
             << "PARAMS: flag_mask: " << flag_mask << endl
             << "PARAMS: filename_mask: " << (char *) ( ( filename_mask ) ? filename_mask : "0" ) << endl
             << "PARAMS: numb_classes: " << numb_classes << endl
             << "PARAMS: flag_manual_priors: " << flag_manual_priors << endl;

        if ( filename_priors )
        {
            if ( flag_priors4D )
                sout << "PARAMS: filename_priors[4D]: " << (char *) ( ( filename_priors[0] ) ? filename_priors[0] : "0" ) << endl;
            else
                for(long classnum=0; classnum<numb_classes; classnum++)
                    sout << "PARAMS: filename_priors["
                         << classnum << "]: " << (char *) ( ( filename_priors[classnum] ) ? filename_priors[classnum] : "0" ) << endl;
        }
        else
            sout << "PARAMS: filename_priors: 0" << endl;

        sout << "PARAMS: flag_Outlierness: " << flag_Outlierness << endl
             << "PARAMS: OutliernessThreshold: " << OutliernessThreshold << endl
             << "PARAMS: OutliernessRatio: " << OutliernessRatio << endl
             << "PARAMS: flag_out_outlier: " << flag_out_outlier << endl
             << "PARAMS: filename_out_outlier: " << (char *) ( ( filename_out_outlier ) ? filename_out_outlier : "0" ) << endl
             << "PARAMS: flag_MAP: " << flag_MAP << endl;

        if ( MAP_M )
            for(long classnum=0; classnum<numb_classes; classnum++)
                sout << "PARAMS: MAP_M[" << classnum << "]: " << MAP_M[classnum] << endl;
        else
            sout << "PARAMS: MAP_M: " << MAP_M << endl;

        if ( MAP_V )
            for(long classnum=0; classnum<numb_classes; classnum++)
                sout << "PARAMS: MAP_V[" << classnum << "]: " << MAP_V[classnum] << endl;
        else
            sout << "PARAMS: MAP_V: " << MAP_V << endl;

        sout << "PARAMS: flag_PV_model: " << flag_PV_model << endl
             << "PARAMS: flag_SG_deli: " << flag_SG_deli << endl
             << "PARAMS: flag_Bias: " << flag_Bias << endl
             << "PARAMS: bias_order: " << bias_order << endl
             << "PARAMS: Bias_threshold: " << Bias_threshold << endl
             << "PARAMS: relax_factor: " << relax_factor << endl
             << "PARAMS: relax_gauss_kernel: " << relax_gauss_kernel << endl
             << "PARAMS: flag_bc_out: " << flag_bc_out << endl
             << "PARAMS: filename_bias: " << (char *) ( ( filename_bias ) ? filename_bias : "0" ) << endl
             << "PARAMS: flag_MRF: " << flag_MRF << endl
             << "PARAMS: MRF_strength: " << MRF_strength << endl
             << "PARAMS: maxIteration: " << maxIteration << endl
             << "PARAMS: minIteration: " << minIteration << endl
             << "PARAMS: aprox: " << aprox << endl;

    };
} SEG_PARAM;

typedef struct
{
    double loglik;
    SegPrecisionTYPE oldloglik;
    int improv_phase;
    bool prior_relax;
    bool do_pv_modeling;
    bool pv_modeling_on;
    bool sg_delineation;
    bool out;
} FLAGS;


