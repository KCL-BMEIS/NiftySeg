# pragma once
/** @file
 * @defgroup DEFINES_SEG_COMMON seg_common.h defines
 * @{  */
//Global includes

/// @brief  will provide many of the math defines and constants
#define _USE_MATH_DEFINES
#define NOMINMAX

#if ((defined(_WIN32) || defined(_WINDOWS)) && !defined(__CYGWIN__))
    #include <dirent_win.h>
    #include <float.h>
    #include <time.h>
    #include <malloc.h>

    #ifndef M_PI
        /// @brief M_PI as to be defined in windows
        #define M_PI (3.14159265358979323846)
    #endif

#if !(defined _MSC_VER && _MSC_VER >= 1800)
    template<typename T> inline bool isinf(T value)
    {
        return std::numeric_limits<T>::has_infinity && value == std::numeric_limits<T>::infinity();
    }

    #ifndef isnan(_X)
        /// @brief Define necessary to make it multi-platform compatible
        #define isnan(_X) _isnan(_X)
    #endif

    #ifndef round(_x)
        /// @brief Define necessary to make it multi-platform compatible, as Windows does not necessarily have a round function
        #define round(_x) floor(_x + 0.5)
    #endif
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
    #include <alloca.h>
#endif

#if (defined(_WIN32) || defined(_WINDOWS)) && !defined(__CYGWIN__)
    /// @brief Defines the system wide path separator for filenames
    #define SEP "\\"
#else
    /// @brief Defines the system wide path separator for filenames
    #define SEP "/\0"
#endif

#include <sys/types.h>
#include <errno.h>
#include <vector>
#include <cstdio>
#include <new>
#include <fstream>
#include <limits>
#include <cstdlib>
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "nifti1_io.h"
#include <new>
#include <exception>
#include <algorithm>
#ifdef _OPENMP
#include "omp.h"
#endif

// NiftySeg data defines
/// @brief Defines the overall data type used internally. As float is hardcoded in some places, please do not change this for now.
#define segPrecisionTYPE float
/// @brief Defines the data type for categorical labels in seg_LabFusion.
#define categoricalLabelType unsigned char

/// @brief Defines the number of non-PV classes used in seg_LoAd. As seg_LoAd will be deprecated, this will also disapear.
#define nonPVNumClass 5
/// @brief Defines the maximum number of classes the algorithm can infer. This is only used for memory and speed reasons. The value can be increase if necessary.
#define maxNumbClass 12
/// @brief Defines the maximum order of the polynomial for the Bias correction. This is only used for memory and speed reasons. The value can be increase if necessary.
#define maxAllowedBCPowerOrder 6
/// @brief Defines the subsampling size to estimate the bias field. The intensities/function are sampled only every redux_factor_for_bias voxels
#define reduxFactorForBias 2
/// @brief Defines the maximum number of multimodal images (_nu)
#define maxMultispectalSize 6
/// @brief Defines the maximum number of labels for seg_LabFusion
#define maxMultiLableClass 256



// seg_LoAD define (deprecated)
// These define the expected order of the tissues if the seg_LoAd model is used. This is deprecated.
#define WMclass 0
#define GMclass 1
#define CSFclass 2
#define dGMclass 3
#define iCSFclass 4
#define WMGMpvclass 5
#define GMCSFpvclass 6


#if ( defined ( __FreeBSD__ ) || !defined(_LARGEFILE64_SOURCE) )
     /// @brief Defines the readdir64 as readdir because different systems define it in different ways
     #define readdir64 readdir
     /// @brief Defines the dirent64 as dirent because different systems define it in different ways
     #define dirent64 dirent
 #endif
/** @} */

using namespace std;

/// @brief Structure defining the size of the image.
///
/// This structure will store the size of the image with and without mask, and also the intesity rescalings used to normalise the image.
/// This structure is used to make it easier to pass these values around between functions.
/// We could have used the nifti image structure instead, however, the nifty image structure limits the size of the input in x y and z to short, which is not enough some times.
typedef struct
{
    long xsize;
    long ysize;
    long zsize;
    long usize;
    long tsize;
    long numclass;
    long nummod;
    long numel;
    long numelmasked;
    long numelbias;
    float rescale_max[maxMultispectalSize];
    float rescale_min[maxMultispectalSize];
} ImageSize;

/// @brief Structure defining the flags of the seg_EM segmentation.
///
/// This is used to know the stage of the segmentation process. It could and should be integrated into the functions themselves.
typedef struct
{
    double loglik;
    segPrecisionTYPE oldloglik;
    int improv_phase;
    bool prior_relax;
    bool do_pv_modeling;
    bool pv_modeling_on;
    bool sg_delineation;
    bool out;
} seg_EM_Flags;



/// @brief Structure defining the parameters of the seg_EM segmentation, used mostly to comunicate with NiftyView.
///
/// This is used to pass the current segmentation parameters around. This extends seg_EM_Flags and is used to interact with NiftyView
/// As NiftyView is going to use CLI modules in the future (instead of actually inking against niftyseg) this will probably disapear in the near future.
typedef struct seg_EM_Params
{
    segPrecisionTYPE relax_factor;
    segPrecisionTYPE relax_gauss_kernel;
    int  bias_order;
    float Bias_threshold;
    segPrecisionTYPE  MRF_strength;
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

    seg_EM_Params( void )
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
}seg_EM_Params;


inline segPrecisionTYPE pow_int(const segPrecisionTYPE base,
                                int exp)
{
    if(exp==0)
    {
        return 1;
    }
    segPrecisionTYPE result = base;
    while (--exp)
    {
        result *= base;
    }
    return result;
}
