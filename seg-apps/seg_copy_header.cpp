#include "seg_copy_header.h"

#include "_seg_common.h"
#include <iostream>
#include <time.h>

using namespace std;
#define PrecisionTYPE float



int main(int argc, char **argv)
{


    nifti_image * FromImage=nifti_image_read(argv[1],false);
    if(FromImage == NULL){
        fprintf(stderr,"USAGE: seg_copy_header <from>.nii <to>.nii\n");
        return 1;
    }

    nifti_image * ToImage = nifti_image_read(argv[2],true);
    if(ToImage == NULL){
        fprintf(stderr,"USAGE: seg_copy_header <from>.nii <to>.nii\n");
        return 1;
    }


    ToImage->dim[0]=FromImage->dim[0];
    ToImage->dim[1]=FromImage->dim[1];
    ToImage->dim[2]=FromImage->dim[2];
    ToImage->dim[3]=FromImage->dim[3];
    ToImage->dim[4]=FromImage->dim[4];
    ToImage->dim[5]=FromImage->dim[5];
    ToImage->pixdim[0]=FromImage->pixdim[0];
    ToImage->pixdim[1]=FromImage->pixdim[1];
    ToImage->pixdim[2]=FromImage->pixdim[2];
    ToImage->pixdim[3]=FromImage->pixdim[3];
    ToImage->pixdim[4]=FromImage->pixdim[4];
    ToImage->pixdim[5]=FromImage->pixdim[5];
    ToImage->dx=FromImage->dz;
    ToImage->dy=FromImage->dy;
    ToImage->dz=FromImage->dx;
    ToImage->dt=FromImage->dt;
    ToImage->du=FromImage->du;
    ToImage->nx=FromImage->nz;
    ToImage->ny=FromImage->ny;
    ToImage->nz=FromImage->nx;
    ToImage->nt=FromImage->nt;
    ToImage->nu=FromImage->nu;
    ToImage->qform_code=FromImage->qform_code;
    ToImage->sform_code=FromImage->sform_code;
    ToImage->ndim=FromImage->ndim;
    ToImage->qoffset_x=FromImage->qoffset_x;
    ToImage->qoffset_y=FromImage->qoffset_y;
    ToImage->qoffset_z=FromImage->qoffset_z;
    ToImage->quatern_b=FromImage->quatern_b;
    ToImage->quatern_c=FromImage->quatern_c;
    ToImage->quatern_d=FromImage->quatern_d;
    ToImage->qfac=FromImage->qfac;

    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
        ToImage->qto_xyz.m[i][j]=FromImage->qto_xyz.m[i][j];
        ToImage->qto_ijk.m[i][j]=FromImage->qto_ijk.m[i][j];
        ToImage->sto_xyz.m[i][j]=FromImage->sto_xyz.m[i][j];
        ToImage->sto_ijk.m[i][j]=FromImage->sto_ijk.m[i][j];
    }
    }

    nifti_update_dims_from_array(ToImage);



    nifti_image_write(ToImage);
    nifti_image_free(ToImage);
    nifti_image_free(FromImage);



    return 0;
}

