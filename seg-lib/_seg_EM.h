#pragma once
#include "NiftySegWinExportHeader.h"

#include "_seg_common.h"
#include "_seg_matrix.h"
extern "C++" {
#include "_seg_tools.h"
}

/** @class seg_EM
 *  @brief General purpouse EM segmentation class
 */
class NIFTYSEG_WINEXPORT seg_EM
{
protected:

    /// @brief A pointer to an externally defined input data. Do not delete this pointer internally.
    nifti_image*    InputImage;
    /// @brief A status bool assigning if the image has been read.
    bool    inputImage_status;
    /// @brief The level of verbose of the algorithm. 0=off, 1="balanced amount of info", 2="dump all possible info"
    int     verbose_level;
    /// @brief The filename of the Output result of the segmentation.
    string  filenameOut;

    /// @brief Number of dimensions of the image (2d=2, 3d=3, 4d=4)
    int     dimensions;
    /// @brief Dimention of the input image in X
    int     nx;
    /// @brief Dimention of the input image in Y
    int     ny;
    /// @brief Dimention of the input image in Z
    int     nz;
    /// @brief Dimention of the input image in T
    int     nt;
    /// @brief Dimention of the input image in U
    int     nu;
    /// @brief Spacing between voxels in X
    float   dx;
    /// @brief Spacing between voxels in X
    float   dy;
    /// @brief Spacing between voxels in X
    float   dz;
    /// @brief Number of elements (voxels) in the image
    int     numElements;
    /// @brief Current number of iterations
    int     iter;
    /// @brief The amount of regularisation of the off diagonal components of the covariance matrix
    float reg_factor;
    /// @brief A data structure to pass the dimensions of the image around more easily.
    ImageSize * CurrSizes;


    /// @brief A vector of concatenated class means, with K*D elements (K classes and D data dimensions),
    float*  M;
    /// @brief A vector of concatenated covariance matrices with K*D^2 elements (K classes times a D^2 covariance matrix per class),
    float*  V;
    /// @brief The current model Expectation. This should be the same size as the input image in x,y and z, but with K (number of classes) time points
    float*  Expec;
    /// @brief A memory preserving version of the population priors, with data only available for voxels within the input mask.
    float*  ShortPrior;
    /// @brief A integer array storing how to go from an index in a Short* array (size of the mask) to its corresponding index in the input data.
    int*    S2L;
    /// @brief A integer array doing the oposit from Short_2_Long_Indices,
    /// i.e. storing how to go from an index in the input data to its corresponding Short* array (size of the mask) location.
    int*    L2S;
    /// @brief The actual number of classes used in the mixture model.
    int     numberOfClasses;
    /// @brief The current log likelihood estimates
    double  loglik;
    /// @brief The log likelihood estimates in the previous iteration (used to test for convergence)
    double  oldloglik;
    /// @brief Current estimate of the ratio of convergence used as a stopping criteria
    float   ratio;
    /// @brief The Maximum number of iterations (overrides the covergenge criteria "ratio")
    int     maxIteration;
    /// @brief The Minimum number of iterations (overrides the covergenge criteria "ratio")
    int     minIteration;
    /// @brief The convergence criteria ratio for the LogLik increase
    float convCrit;

    /// @brief A float vector of size maxMultispectalSize containing the maximum intensity value used for the intensity rescaling
    float rescale_max[maxMultispectalSize];
    /// @brief A float vector of size maxMultispectalSize containing the minimum intensity value used for the intensity rescaling
    float rescale_min[maxMultispectalSize];

    // Mask
    /// @brief A pointer to the mask (region of interest). Do not delete this pointer if maskImage_status==1, delete if maskImage_status==2.
    nifti_image*  Mask;
    /// @brief A boolean defining if an image mask is being used.
    unsigned char maskImageStatus;
    /// @brief Defines the number of elements (voxels) whithin the mask, i.e. same as numel, but only within the mask ROI
    int     numElementsMasked;
    /// @brief Defines the number of elements (voxels) whithin the bias field image, i.e. numel rescaled by the subsampling factor reduxFactorForBias
    int     numelBias;

    // Priors Specific
    /// @brief A pointer to an externally defined set of priors. Do not delete this pointer internally.
    nifti_image*  Priors;
    /// @brief A boolean defining if population priors were used.
    bool    priorsStatus;

    // MRF Specific
    /// @brief A boolean defining if the MRF model is used
    bool    mrfStatus;
    /// @brief A float point defining the off diagonal energy term of the MRF energy matrix
    float   mrfStrength;
    /// @brief A float image of size numelmasked*numb_classes storing the MRF probability per class
    float*  MRF;
    /// @brief A float image storing the per class MRF probability
    float*  MRFBeta;
    /// @brief A float matrix used to actually store the MRF energy matrix. This allows for more complicated energy matrices than the offdiagonal ones defined by MRF_strength
    float*  MRFTransitionMatrix;

    /// @brief A float image of size numelmasked used to store the outlierness of the model as deffined in Koen's paper
    float * Outlierness;
    /// @brief A float pointer set to NULL of the outlierness is not beein used, or set to Outlierness if the outlierness is beeing used
    float * OutliernessUSE;
    /// @brief A boolean defining if the outlierness model is used
    bool    outliernessStatus;
    /// @brief A float storing the mahalonobis distance threshold as defined in Koen's paper
    float   outliernessThreshold;
    /// @brief Only start optimising the Outlieness part of the model if the convergence ratio is below Outlierness_ratio
    float   outliernessRatio;

    // BiasField Specific
    /// @brief A boolean defining if the BiasField model is used
    bool    biasFieldStatus;
    /// @brief An integer defining the order of the polynomial used to model the bias field
    int     biasFieldOrder;
    /// @brief A float image of size nt*numel containing the current estimate of the bias field per time point
    float*  BiasField;
    /// @brief A float vector of size ((BiasField_order+1)*(BiasField_order+2)/2*(BiasField_order+3)/3)*nu*nt containing all the polynomial basis' coefficients
    float*  biasFieldCoeficients;
    /// @brief Only start optimising the BiasField part of the model if the convergence ratio is below BiasField_ratio
    float   biasFieldRatio;

    // This section is specific for the LoAd algorithm, and only used under those conditions
    /// @brief An integer defining the current stage of the LoAd pipeline (see the original paper)
    int     stageLoAd;
    /// @brief A boolean defining if we are modeling partial volume by adding a few extra classes. Note that the classes have to be in the correct order for this to work.
    bool    pvModelStatus;
    /// @brief A boolean defining if the Sulci and Giri delineation improvement is beeing used
    bool    sgDelineationStatus;

    // This section is specific for the AdaPT algorithm, but can be used in combination with other parameters (see AdaPT paper)
    /// @brief A boolean defining if the priors relaxation should be used
    bool    relaxStatus;
    /// @brief A float factor defining the alpha and (1-alpha) of the relaxation, as stated in the AdaPT paper.
    float   relaxFactor;
    /// @brief A float defining the size of the Gaussian Kernel (in voxels) for prior relaxation regularisation
    float   relaxGaussKernelSize;

    /// @brief A boolean defining if the MAP regularisation should be used
    bool    mapStatus;
    /// @brief A float vector of size numb_classes containing the mean of the MAP prior per class (only works for monomodal data)
    float*  MAP_M;
    /// @brief A float vector of size numb_classes containing the variance of the MAP prior per class (only works for monomodal data)
    float*  MAP_V;



    // Private funcs
    void     CreateDiagonalMRFTransitionMatrix();
    void     CreateCurrSizes();
    void     InitializeAndNormalizeImageAndPriors();
    void     InitializeAndAllocate();
    void     InitializeMeansUsingIntensity();
    void     InitializeAndNormalizeNaNPriors();
    void     InitializeAndNormalizeImage();
    void     CreateShort2LongMatrix();
    void     CreateLong2ShortMatrix();

    void     RunMaximization();
    void     RunExpectation();
    void     RunPriorRelaxation();
    void     RunMRF();
    void     RunMRF3D();
    void     RunMRF2D();
    void     RunBiasField();
    void     RunBiasField3D();
    void     RunBiasField2D();


public:
    seg_EM (int _numb_classes,int _nu,int _nt);
    ~seg_EM ();

    void     Run_EM();

    // Setters of the options
    void    SetInputImage(nifti_image *);
    void    SetMaskImage(nifti_image *);
    void    SetPriorImage(nifti_image *_r);
    void    SetFilenameOut(char *_f);
    void    SetMAP(float *_M, float* _V);
    void    SetRegValue(float reg);
    void    SetRelaxation(float relax_factor,float relax_gauss_kernel);
    void    SetMRF(float MRF_strenght);
    void    SetOutlierness(float _OutliernessThreshold, float _OutliernessRatio);
    void    SetBiasField(int _BiasFieldOrder, float _BiasFieldRatio);
    void    SetMaximalIterationNumber(unsigned int numberiter);
    void    SetMinIterationNumber(unsigned int numberiter);
    void    SetConvergenceCriteria(float inConvCrit);
    void    SetVerbose(unsigned int verblevel);

    // Getters of the results
    float*   GetMeans();
    float*   GetSTD();
    nifti_image* GetResult();
    nifti_image* GetBiasCorrected(char * filename);
    nifti_image* GetOutlierness(char * filename);
};

