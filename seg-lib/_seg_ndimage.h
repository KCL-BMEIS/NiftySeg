#ifndef SEG_NDIMAGE_H
#define SEG_NDIMAGE_H

#include "_seg_common.h"

class ndimage
{

    bool is_initialized;

    // data size and type
    short datatype;
    short nbyper;
    float pixdim[8];
    int dim[8];
    int size_xyz;
    int size_all;

    // data presentation
    float scl_slope;
    float scl_inter;
    float cal_max;
    float cal_min;

    // data orientation
    short qform_code;
    short sform_code;
    float quatern_b;
    float quatern_c;
    float quatern_d;
    float qoffset_x;
    float qoffset_y;
    float qoffset_z;
    float sform_matrix[4][4];

    // data itself
    void * data; // if the image is not masked, then data has size this->size_all
    // if the image is masked, then data has size this->size_masked*(dim[4]*dim[5]*dim[6]*dim[7])

    // data mask
    bool is_masked;   // assignes of the image is masked or not
    bool size_masked; // number of elements equal to 1 within the mask
    bool * mask; // the mask, if exists, is always the size of dim[1]*dim[2]*dim[3]
    // i.e., is the same for all time points

    // data name

    char * filename;

public:
    ndimage(nifti_image * input_image); // converts input_image into a ndimage
    ndimage(nifti_image * input_image,bool destroy_input); //LOW_MEMORY: it will steal data from input
    ndimage(char * input_image_name);  //read from disk
    ndimage();  //just create structure but do not initialise
    ~ndimage();

    void dealocate();
};





#endif // SEG_NDIMAGE_H
