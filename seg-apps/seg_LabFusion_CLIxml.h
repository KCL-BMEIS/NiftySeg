/** @file seg_EM_CLIxml.h
 * @date 23/10/2012
 * @author John Hipwell
 * @brief Header file that contains the text to be output
 * for the NifTK command line interface (CLI) module and slicer extensions of seg_EM.
 */


std::string xml_segLabFusion =
    "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"

    // Executable description
    // ~~~~~~~~~~~~~~~~~~~~~~

    "<executable>\n"

    "  <category>Segmentation.NiftySeg</category>\n"
    "  <title>Label Fusion (seg_LabFusion)</title>\n"
    "  <description><![CDATA[Module/executable for label fusion.]]></description>\n"
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


    "  </parameters>\n"

    "</executable>\n"
    ;
