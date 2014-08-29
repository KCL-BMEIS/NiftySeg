/** @file seg_EM_CLIxml.h
 * @date 23/10/2012
 * @author John Hipwell
 * @brief Header file that contains the text to be output
 * for the NifTK command line interface (CLI) module and slicer extensions of seg_EM.
 */


std::string xml_segEM =
    "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"

    // Executable description
    // ~~~~~~~~~~~~~~~~~~~~~~

    "<executable>\n"

    "  <category>Segmentation.NiftySeg</category>\n"
    "  <title>EM-MRF Segmentation (seg_EM)</title>\n"
    "  <description><![CDATA[Module/executable for general purpose intensity based, "
    "Expectation-Maximisation, statistical image segmentation.]]></description>\n"
    "  <version>0.0.1</version>\n"
    "  <documentation-url>http://cmic.cs.ucl.ac.uk/home/software/</documentation-url>\n"
    "  <license>BSD</license>\n"
    "  <contributor>M. Jorge Cardoso, Marc Modat, Matt Clarkson and Sebastien Ourselin (UCL)</contributor>\n"

    // The mandatory input parameters
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    "  <parameters advanced=\"false\">\n"

    "    <label>Mandatory Input Parameters</label>\n"
    "    <description>These parameters are the minimum necessary for successful execution of the module.</description>\n"

    // Filename of the input image

    "    <image fileExtensions=\".nii,.nii.gz\">\n"
    "      <name>inputImageName</name>\n"
    "      <longflag>in</longflag>\n"
    "      <description>Filename of the input image to be segmented. The input image should be 2D, 3D or 4D images. 2D images should be on the XY plane. 4D images are segmented as if they were multimodal.</description>\n"
    "      <label>Input image</label>\n"
    "      <default>required</default>\n"
    "      <channel>input</channel>\n"
    "    </image>\n"

    // Filename of the output segmented image

    "    <image fileExtensions=\".nii,.nii.gz\">\n"
    "      <name>outputImageName</name>\n"
    "      <longflag>out</longflag>\n"
    "      <description>Filename of the output segmented image.</description>\n"
    "      <label>Output segmented image</label>\n"
    "      <default>SegmentedImage.nii</default>\n"
    "      <channel>output</channel>\n"
    "    </image>\n"

    // The number of classes to segment

    "    <integer>\n"
    "      <name>numb_classes</name>\n"
    "      <longflag>nopriors</longflag>\n"
    "      <description>The number of classes to segment</description>\n"
    "      <label>Either...\n   number of classes [maximum is 5]</label>\n"
    "      <default></default>\n"
    "      <constraints>\n"
    "        <minimum>0</minimum>\n"
    "        <maximum>5</maximum>\n"
    "        <step>1</step>\n"
    "      </constraints>\n"
    "    </integer>\n"

    // Filename of an image of 4D priors

    "    <image fileExtensions=\".nii,.nii.gz\">\n"
    "      <name>filename_priors</name>\n"
    "      <longflag>priors4D</longflag>\n"
    "      <description>Filename of a 4D image of priors stacked in the 4th dimension. All priors should be registered to the input image.</description>\n"
    "      <label>or...\n   priors (4D image)</label>\n"
    "      <default></default>\n"
    "      <channel>input</channel>\n"
    "    </image>\n"


    "  </parameters>\n"


    // General Options
    // ~~~~~~~~~~~~~~~

    "  <parameters advanced=\"true\">\n"

    "    <label>General Options</label>\n"
    "    <description>Optional input parameters</description>\n"

    // Verbose output

    "    <integer-enumeration>\n"
    "      <name>verbose_level</name>\n"
    "      <flag>v</flag>\n"
    "      <description>The level of information printed during execution of the module [0=off (default), 1=on, 2=debug].</description>\n"
    "      <label>Verbose output</label>\n"
    "      <default>0</default>\n"
    "      <element>0</element>\n"
    "      <element>1</element>\n"
    "      <element>2</element>\n"
    "    </integer-enumeration>\n"

    // Filename of an optional mask image

    "    <image fileExtensions=\".nii,.nii.gz\">\n"
    "      <name>filename_mask</name>\n"
    "      <longflag>mask</longflag>\n"
    "      <description>Filename of an optional mask for the input image.</description>\n"
    "      <label>Mask</label>\n"
    "      <default></default>\n"
    "      <channel>input</channel>\n"
    "    </image>\n"

    // The maximum number of iterations

    "    <integer>\n"
    "      <name>maxIteration</name>\n"
    "      <longflag>max_iter</longflag>\n"
    "      <description>The maximum number of iterations to perform (default is 100).</description>\n"
    "      <label>Number of iterations</label>\n"
    "      <default>100</default>\n"
    "      <constraints>\n"
    "        <minimum>1</minimum>\n"
    "        <maximum>1000000</maximum>\n"
    "        <step>1</step>\n"
    "      </constraints>\n"
    "    </integer>\n"

    // MRF prior strength

    "    <float>\n"
    "      <name>MRF_strength</name>\n"
    "      <longflag>mrf_beta</longflag>\n"
    "      <description>MRF prior strength [off = 0, max = 1] (default is 0.4).</description>\n"
    "      <label>MRF prior strength</label>\n"
    "      <default>0.4</default>\n"
    "      <constraints>\n"
    "        <minimum>0</minimum>\n"
    "        <maximum>1</maximum>\n"
    "        <step>0.01</step>\n"
    "      </constraints>\n"
    "    </float>\n"

#if 0 // (Rarely used and only relevant for multi-modal images)

    // Regularisation

    "    <float>\n"
    "      <name>regularization_amount</name>\n"
    "      <longflag>reg</longflag>\n"
    "      <description>Amount of regularization over the diagonal of the covariance matrix (greater than 1).</description>\n"
    "      <label>Regularisation</label>\n"
    "      <default></default>\n"
    "      <constraints>\n"
    "        <minimum>1</minimum>\n"
    "        <step>0.01</step>\n"
    "      </constraints>\n"
    "    </float>\n"
#endif

    "  </parameters>\n"


    // Bias Field
    // ~~~~~~~~~~

    "  <parameters advanced=\"true\">\n"

    "    <label>Bias Field Correction</label>\n"
    "    <description>Parameters which control correction of a background (MRI) bias field.</description>\n"

    // The bias field correction polynomial order

    "    <integer>\n"
    "      <name>bias_order</name>\n"
    "      <longflag>bc_order</longflag>\n"
    "      <description>Polynomial order to use for bias field correction [min 0, max 5] (default is 3).</description>\n"
    "      <label>Bias field polynomial order</label>\n"
    "      <default>3</default>\n"
    "      <constraints>\n"
    "        <minimum>0</minimum>\n"
    "        <maximum>5</maximum>\n"
    "        <step>1</step>\n"
    "      </constraints>\n"
    "    </integer>\n"

    // Bias field correction threshold

    "    <float>\n"
    "      <name>Bias_threshold</name>\n"
    "      <longflag>bc_thresh</longflag>\n"
    "      <description>Bias field correction will only run if the ratio of improvement is below this threshold (default is 0 i.e. OFF).</description>\n"
    "      <label>Ratio threshold</label>\n"
    "      <default>0</default>\n"
    "    </float>\n"

    // Filename of the output bias field correction

    "    <image fileExtensions=\".nii,.nii.gz\">\n"
    "      <name>filename_bias</name>\n"
    "      <longflag>bc_out</longflag>\n"
    "      <description>File to write the bias field corrected image to.</description>\n"
    "      <label>Output corrected image</label>\n"
    "      <default>BiasFieldCorrected.nii</default>\n"
    "      <channel>output</channel>\n"
    "    </image>\n"

    "  </parameters>\n"

    // Outlier Detection
    // ~~~~~~~~~~~~~~~~~

    "  <parameters advanced=\"true\">\n"

    "    <label>Outlier Detection</label>\n"
    "    <description>Parameters which control detection of classification outliers.</description>\n"

    // Outlier threshold

    "    <float>\n"
    "      <name>OutliernessThreshold</name>\n"
    "      <longflag>outlier_thresh</longflag>\n"
    "      <description>The Mahalanobis threshold for outlier detection [recommended between 3 and 7] (as in Van Leemput TMI 2003).</description>\n"
    "      <label>Outlier threshold</label>\n"
    "      <default>0.</default>\n"
    "      <constraints>\n"
    "        <minimum>0</minimum>\n"
    "        <step>0.2</step>\n"
    "      </constraints>\n"
    "    </float>\n"

// Outlier ratio

    "    <float>\n"
    "      <name>OutliernessRatio</name>\n"
    "      <longflag>outlier_ratio</longflag>\n"
    "      <description>Convergence ratio below which the outlier detection is going to be done [recommended 0.01] (as in Van Leemput TMI 2003).</description>\n"
    "      <label>Outlier ratio</label>\n"
    "      <default>0.01</default>\n"
    "      <constraints>\n"
    "        <minimum>0</minimum>\n"
    "        <step>0.2</step>\n"
    "      </constraints>\n"
    "    </float>\n"

// Outlier output filename

    "    <image fileExtensions=\".nii,.nii.gz\">\n"
    "      <name>filename_out_outlier</name>\n"
    "      <longflag>out_outlier</longflag>\n"
    "      <description>File to write the outlierness image to.</description>\n"
    "      <label>Output outlier image</label>\n"
    "      <default>Outliers.nii</default>\n"
    "      <channel>output</channel>\n"
    "    </image>\n"


    "  </parameters>\n"

    // Relaxation
    // ~~~~~~~~~~

    "  <parameters advanced=\"true\">\n"

    "    <label>Relaxation</label>\n"
    "    <description>Parameters which control relaxation of the priors.</description>\n"

    // Relaxation factor

    "    <float>\n"
    "      <name>relax_factor</name>\n"
    "      <longflag>rf_factor</longflag>\n"
    "      <description>Priors relaxation factor, [between 0 and 1] (recommended: 0.5).</description>\n"
    "      <label>Relaxation factor</label>\n"
    "      <default>0.</default>\n"
    "      <constraints>\n"
    "        <minimum>0</minimum>\n"
    "        <maximum>1</maximum>\n"
    "        <step>0.01</step>\n"
    "      </constraints>\n"
    "    </float>\n"

    // Gaussian regularization

    "    <float>\n"
    "      <name>relax_gauss_kernel</name>\n"
    "      <longflag>rf_gauss_kernel</longflag>\n"
    "      <description>Gaussian regularization standard deviation (for 3D only), must be greater than zero [recommended=2.0].</description>\n"
    "      <label>Gaussian regularization</label>\n"
    "      <default>0.</default>\n"
    "      <constraints>\n"
    "        <minimum>0</minimum>\n"
    "        <step>0.1</step>\n"
    "      </constraints>\n"
    "    </float>\n"

    "  </parameters>\n"


    // Relaxation
    // ~~~~~~~~~~

    "  <parameters advanced=\"true\">\n"

    "    <label>Maximum a Posteriori Formulation</label>\n"
    "    <description>Parameters which specify the semiconjugate priors over the mean of each class.</description>\n"

    "    <integer-vector>\n"
    "      <name>neighborhood</name>\n"
    "      <longflag>MAP_MV1</longflag>\n"
    "      <description>The mean and variance of the semiconjugate prior over the mean for class 1</description>\n"
    "      <label>Class 1 (mean,var)</label>\n"
    "      <default></default>\n"
    "    </integer-vector>\n"

    "    <integer-vector>\n"
    "      <name>neighborhood</name>\n"
    "      <longflag>MAP_MV2</longflag>\n"
    "      <description>The mean and variance of the semiconjugate prior over the mean for class 1</description>\n"
    "      <label>Class 2 (mean,var)</label>\n"
    "      <default></default>\n"
    "    </integer-vector>\n"

    "    <integer-vector>\n"
    "      <name>neighborhood</name>\n"
    "      <longflag>MAP_MV3</longflag>\n"
    "      <description>The mean and variance of the semiconjugate prior over the mean for class 1</description>\n"
    "      <label>Class 3 (mean,var)</label>\n"
    "      <default></default>\n"
    "    </integer-vector>\n"

    "    <integer-vector>\n"
    "      <name>neighborhood</name>\n"
    "      <longflag>MAP_MV4</longflag>\n"
    "      <description>The mean and variance of the semiconjugate prior over the mean for class 1</description>\n"
    "      <label>Class 4 (mean,var)</label>\n"
    "      <default></default>\n"
    "    </integer-vector>\n"

    "    <integer-vector>\n"
    "      <name>neighborhood</name>\n"
    "      <longflag>MAP_MV5</longflag>\n"
    "      <description>The mean and variance of the semiconjugate prior over the mean for class 1</description>\n"
    "      <label>Class 5 (mean,var)</label>\n"
    "      <default></default>\n"
    "    </integer-vector>\n"

    "  </parameters>\n"

    "</executable>\n"
    ;
