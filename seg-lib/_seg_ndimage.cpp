#include "_seg_ndimage.h"

ndimage::ndimage(nifti_image * input_image) // creates a copy of the data
{
  this->is_initialized=1;
  // copy size
  this->filename = (char *) malloc(strlen(input_image->fname) + 1);
  strcpy(this->filename, input_image->fname);
  this->datatype=input_image->datatype;
  this->nbyper=input_image->nbyper;
  this->size_xyz=1;
  this->size_all=1;
  for(unsigned int i=0; i<8; i++){
      this->dim[i]=input_image->dim[i];
      this->pixdim[i]=input_image->pixdim[i];
      this->size_xyz*=(i>0 && i<4)?(input_image->dim[i]>0?input_image->dim[i]:1):1;
      this->size_all*=i>0?(input_image->dim[i]>0?input_image->dim[i]:1):1;
    }
  this->scl_slope=input_image->scl_slope;
  this->scl_inter=input_image->scl_inter;

  // copy orientation
  this->qform_code=input_image->qform_code;
  this->sform_code=input_image->sform_code;
  this->quatern_b=quatern_b;
  this->quatern_c=quatern_c;
  this->quatern_d=quatern_d;
  this->qoffset_x=qoffset_x;
  this->qoffset_y=qoffset_y;
  this->qoffset_z=qoffset_z;
  for(unsigned int i=0; i<4; i++){
      for(unsigned int j=0; j<4; j++){
          this->sform_matrix[i][j]=input_image->sto_xyz.m[i][j] ;
        }
    }
  this->data=calloc(this->size_all,this->nbyper); //allocate space for the data
  if(this->data==NULL){
      // Allocation failed... die gracefully (dont forget the try/catch statement)
      throw std::bad_alloc();
    }
  else{
      memcpy(this->data, input_image->data,this->size_all*this->nbyper); //copy of the data
      this->is_masked=0;
      this->size_masked=this->size_xyz;
      this->mask=NULL;
    }
}

ndimage::ndimage(nifti_image * input_image,bool destroy_input) // This method will take posetion of the data. Use in case of low memory contrains;
{

  // copy size
  this->is_initialized=1;
  this->filename =  (char *) malloc(strlen(input_image->fname) + 1);
  strcpy(this->filename, input_image->fname);
  this->datatype=input_image->datatype;
  this->nbyper=input_image->nbyper;
  this->size_xyz=1;
  this->size_all=1;
  for(unsigned int i=0; i<8; i++){
      this->dim[i]=input_image->dim[i];
      this->pixdim[i]=input_image->pixdim[i];
      this->size_xyz*=(i>0 && i<4)?(input_image->dim[i]>0?input_image->dim[i]:1):1;
      this->size_all*=i>0?(input_image->dim[i]>0?input_image->dim[i]:1):1;
    }
  this->scl_slope=input_image->scl_slope;
  this->scl_inter=input_image->scl_inter;

  // copy orientation
  this->qform_code=input_image->qform_code;
  this->sform_code=input_image->sform_code;
  this->quatern_b=quatern_b;
  this->quatern_c=quatern_c;
  this->quatern_d=quatern_d;
  this->qoffset_x=qoffset_x;
  this->qoffset_y=qoffset_y;
  this->qoffset_z=qoffset_z;
  for(unsigned int i=0; i<4; i++){
      for(unsigned int j=0; j<4; j++){
          this->sform_matrix[i][j]=input_image->sto_xyz.m[i][j] ;
        }
    }


  if(destroy_input){
      this->data=input_image->data;
      this->is_masked=0;
      this->size_masked=this->size_xyz;
      this->mask=NULL;
      input_image->data=NULL;
      nifti_image_free(input_image);
      input_image=NULL;
    }
  else{
      this->data=calloc(this->size_all,this->nbyper); //allocate space for the data
      if(this->data==NULL){
          // Allocation failed... die gracefully (dont forget the try/catch statement)
          throw std::bad_alloc();
        }
      else{
          memcpy(this->data, input_image->data,this->size_all*this->nbyper); //copy of the data
          this->is_masked=0;
          this->size_masked=this->size_xyz;
          this->mask=NULL;
        }
    }
}

ndimage::ndimage(char * input_image_name) // This method will take posetion of the data. Use in case of low memory contrains;
{

  this->is_initialized=1;
  nifti_image * input_image=nifti_image_read(input_image_name,true);
  this->filename =  (char *) malloc(strlen(input_image->fname) + 1);
  strcpy(this->filename, input_image->fname);
  if(input_image==NULL){
      // Read of image failed... die gracefully (dont forget the try/catch statement)
      throw std::bad_alloc();
    }
  else{
      // copy size
      this->datatype=input_image->datatype;
      this->nbyper=input_image->nbyper;
      this->size_xyz=1;
      this->size_all=1;
      for(unsigned int i=0; i<8; i++){
          this->dim[i]=input_image->dim[i];
          this->pixdim[i]=input_image->pixdim[i];
          this->size_xyz*=(i>0 && i<4)?(input_image->dim[i]>0?input_image->dim[i]:1):1;
          this->size_all*=i>0?(input_image->dim[i]>0?input_image->dim[i]:1):1;
        }
      this->scl_slope=input_image->scl_slope;
      this->scl_inter=input_image->scl_inter;

      // copy orientation
      this->qform_code=input_image->qform_code;
      this->sform_code=input_image->sform_code;
      this->quatern_b=quatern_b;
      this->quatern_c=quatern_c;
      this->quatern_d=quatern_d;
      this->qoffset_x=qoffset_x;
      this->qoffset_y=qoffset_y;
      this->qoffset_z=qoffset_z;
      for(unsigned int i=0; i<4; i++){
          for(unsigned int j=0; j<4; j++){
              this->sform_matrix[i][j]=input_image->sto_xyz.m[i][j] ;
            }
        }
      this->data=input_image->data;
      this->is_masked=0;
      this->size_masked=this->size_xyz;
      this->mask=NULL;
      input_image->data=NULL;
      nifti_image_free(input_image);
      input_image=NULL;
    }

}

ndimage::ndimage() // This method will create the structure but not initialise it;
{
  this->is_initialized=0;
  this->filename = NULL;
      // copy size
      this->datatype=0;
      this->nbyper=0;
      this->size_xyz=0;
      this->size_all=0;
      for(unsigned int i=0; i<8; i++){
          this->dim[i]=0;
          this->pixdim[i]=0;
          this->size_xyz*=1;
          this->size_all*=1;
        }
      this->scl_slope=0;
      this->scl_inter=0;

      // copy orientation
      this->qform_code=0;
      this->sform_code=0;
      this->quatern_b=0;
      this->quatern_c=0;
      this->quatern_d=0;
      this->qoffset_x=0;
      this->qoffset_y=0;
      this->qoffset_z=0;
      for(unsigned int i=0; i<4; i++){
          for(unsigned int j=0; j<4; j++){
              this->sform_matrix[i][j]=0;
            }
        }
      this->data=NULL;
      this->is_masked=0;
      this->size_masked=0;
      this->mask=0;
}

ndimage::~ndimage()
{
  this->is_initialized=0;
  // copy size
  if(this->data!=NULL){
      free(this->data);
      this->data=NULL;
    }
  if(this->mask!=NULL){
      free(this->mask);
      this->mask=NULL;
    }
}

void ndimage::dealocate()
{
  this->is_initialized=0;
  // copy size
  if(this->data!=NULL){
      free(this->data);
      this->data=NULL;
    }
  if(this->mask!=NULL){
      free(this->mask);
      this->mask=NULL;
    }
  if(this->filename!=NULL){
      free(this->filename);
      this->filename=NULL;
    }
}
