#include <iostream>
#include <time.h>
#include "_seg_common.h"
#include "_seg_tools.h"
#include "_seg_fill_lesions.h"
#include "_seg_fill_lesions_other.h"

using namespace std;
#define SegPrecisionTYPE float

void Usage(char *exec)
{
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    printf("Usage:\t%s -i <input> -l <lesionmask> -o <output> [options].\n\n",exec);
    printf("\t* * Optional parameters * *\n");
    printf("\t-dil\t<int>\t\tDilate the mask <int> times (in voxels, by default 0).\n");
    printf("\t-match\t<float>\t\tPercentage of minimum number of voxels between patches <float> (by default 0.5).\n");
    printf("\t-search\t<float>\t\tMinimum percentage of valid voxels in target patch <float> (by default 0).\n");
    printf("\t-smo\t<float>\t\tSmoothing by <float> (in minimal 6-neighbourhood voxels (by default 0.1)).\n");
    printf("\t-size\t<int>\t\tSearch regions size respect biggest patch size (by default 4).\n");
    printf("\t-cwf\t<float>\t\tPatch cardinality weighting factor (by default 2).\n");
    printf("\t-mask\t\t\tGive a binary mask with the valid search areas.\n");
    printf("\t-2D\t\t\tUses 2D patches in the Z axis, by default 3D.\n");
    printf("\t-other\t\t\tGuizard et al. (FIN 2015) method, it doesn't include the multiresolution/hierarchical inpainting part, this part needs to be done with some \n");
    printf("\t\t\t\texternal software such as reg_tools and reg_resample from NiftyReg. By default it uses the method presented in Prados et al. (Neuroimage 2016).\n");
    printf("\t-debug\t\t\tSave all intermidium files (by default OFF).\n");
    printf("\t-odt <datatype> \tSet output <datatype> (char, short, int, uchar, ushort, uint, float, double).\n");
    printf("\t-v\t\t\tVerbose (by default OFF).\n");
    #if defined (_OPENMP)
    printf("\t-omp <int>\t\tNumber of openmp threads [%d]\n",omp_get_max_threads());
    #endif
    #ifdef _GIT_HASH
    printf("\t--version\t\tPrint current source code git hash key and exit\n\t\t\t\t(%s)\n",_GIT_HASH);
    #endif
    printf("\tProgram based on the following publication:\n");
    printf("\tF. Prados, M. J. Cardoso, B. Kanber, O. Ciccarelli, R. Kapoor, C. A. M. Gandini Wheeler-Kingshott, S. Ourselin.\n");
    printf("\tA multi-time-point modality-agnostic patch-based method for lesion filling in multiple sclerosis.\n");
    printf("\tNeuroimage 139, 376-384 (2016).\n");
    printf("\n\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    return;
}

void no_memory ()
{
    cout << "Failed to allocate memory!\n";
    exit (1);
}

int main(int argc, char **argv)
{
    try
    {
        bool verbose=0;
        int datatypeoutput=NIFTI_TYPE_FLOAT32;
        double smo=0.1;
        double match=0.5;
        double expanding=0;
        int dil=0;
        bool debug=0;
        int size=4;
	int method=0;
	double k=2.0;
	bool patch2D=false;
        time_t start;
        time(&start);

        set_new_handler(no_memory);
        if (argc < 7)
        {
	    #ifdef _GIT_HASH
	    if(argc>=2 && strcmp(argv[1], "--version")==0)
	    {
	    	printf("%s\n",_GIT_HASH);
	    	return 0;
	    }
	    #endif
            Usage(argv[0]);
            return 0;
        }
        if(strcmp(argv[1], "-help")==0 || strcmp(argv[1], "-Help")==0 ||
                strcmp(argv[1], "-HELP")==0 || strcmp(argv[1], "-h")==0 ||
                strcmp(argv[1], "--h")==0 || strcmp(argv[1], "--help")==0)
        {
            Usage(argv[0]);
            return 0;
        }

        char * filename_in=NULL;
        char * filename_lesion_mask=NULL;
        string filename_out="";
        char * filename_mask=NULL;
	bool ok_input_file=false;
	bool ok_mask_file=false;
	bool ok_output_file=false;
	
        // We read the parameters for the algorithm
        for(long i=1; i<argc; i++)
        {
            if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
                    strcmp(argv[i], "-HELP")==0 || strcmp(argv[i], "-h")==0 ||
                    strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0)
            {
                Usage(argv[0]);
                return 0;
            }
	    else  if(strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "-input") == 0 || strcmp(argv[i], "-in") == 0)
            {
		filename_in=argv[++i];
		ok_input_file=true;
	    }
	    else  if(strcmp(argv[i], "-l") == 0 || strcmp(argv[i], "-lm") == 0 || strcmp(argv[i], "-lesion") == 0 || strcmp(argv[i], "-lesions") == 0)
            {
		filename_lesion_mask=argv[++i];
		ok_mask_file=true;
	    }
	    else  if(strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "-output") == 0 || strcmp(argv[i], "-out") == 0)
            {
		filename_out=argv[++i];
		ok_output_file=true;
	    }
            else  if(strcmp(argv[i], "-size") == 0)
            {
                string parser=argv[++i];
                size=(int)round(strtod(parser.c_str(),NULL));
            }
            else  if(strcmp(argv[i], "-debug") == 0)
            {
                debug=1;
            }
            else  if(strcmp(argv[i], "-dil") == 0)
            {
                string parser=argv[++i];
                dil=(int)round(strtod(parser.c_str(),NULL));
            }
            else  if(strcmp(argv[i], "-smo") == 0)
            {
                string parser=argv[++i];
                smo=strtod(parser.c_str(),NULL);
            }
            else  if(strcmp(argv[i], "-match") == 0)
            {
                string parser=argv[++i];
                match=strtod(parser.c_str(),NULL);
            }
            else  if(strcmp(argv[i], "-cwf") == 0)
            {
                string parser=argv[++i];
                k=strtod(parser.c_str(),NULL);
            }
            else  if(strcmp(argv[i], "-search") == 0)
            {
                string parser=argv[++i];
                expanding=strtod(parser.c_str(),NULL);
            }
	    else  if(strcmp(argv[i], "-2D") == 0)
            {
                patch2D=true;
            }
            else  if(strcmp(argv[i], "-other") == 0)
            {
                method=1;
            }
            else  if(strcmp(argv[i], "-mask") == 0 || strcmp(argv[i], "-m") == 0)
            {
                filename_mask=argv[++i];
            }
	    #if defined (_OPENMP)
            else if(strcmp(argv[i], "-omp")==0 || strcmp(argv[i], "--omp")==0)
            {
                omp_set_num_threads(atoi(argv[++i]));
            }
	    #endif
            // *********************  output data type  *************************
            else if(strcmp(argv[i], "-v") == 0)
            {
                verbose=1;
            }
            else if(strcmp(argv[i], "-odt") == 0)
            {
                string parser=argv[++i];
                if(parser.find("uchar")!=string::npos)
                {
                    datatypeoutput=NIFTI_TYPE_UINT8;
                }
                else if(parser.find("ushort")!=string::npos)
                {
                    datatypeoutput=NIFTI_TYPE_UINT16;
                }
                else if(parser.find("uint")!=string::npos)
                {
                    datatypeoutput=NIFTI_TYPE_UINT32;
                }
                else if(parser.find("char")!=string::npos)
                {
                    datatypeoutput=NIFTI_TYPE_INT8;
                }
                else if(parser.find("short")!=string::npos)
                {
                    datatypeoutput=NIFTI_TYPE_INT16;
                }
                else if(parser.find("int")!=string::npos)
                {
                    datatypeoutput=NIFTI_TYPE_INT32;
                }
                else if(parser.find("float")!=string::npos)
                {
                    datatypeoutput=NIFTI_TYPE_FLOAT32;
                }
                else if(parser.find("double")!=string::npos)
                {
                    datatypeoutput=NIFTI_TYPE_FLOAT64;
                }
                else
                {
                    cerr << "* Error datatype "<< parser << " is unknown"<<endl;
                    i=argc;
                }
            }
            else
            {
                cout << "Option "<< string(argv[i]) << " unkown"<<endl;
                i=argc;
                return 0;
            }
        }
	if (!ok_input_file) {
	    cerr<<"* Error input file not provided -i <filename>"<<endl;
	    Usage(argv[0]);
            return -1;
	}
	if (!ok_mask_file) {
	    cerr<<"* Error lesion mask file not provided -l <filename>"<<endl;
	    Usage(argv[0]);
            return -1;
	}
	if (!ok_output_file) {
	    cerr<<"* Error output file not provided -o <filename>"<<endl;
	    Usage(argv[0]);
            return -1;
	}
	
	// We read the two mains files
        nifti_image * InputImage=nifti_image_read(filename_in,true);
        if(InputImage == NULL)
        {
            fprintf(stderr,"* Error when reading the input image\n");
            return 1;
        }
        if(InputImage->datatype!=NIFTI_TYPE_FLOAT32)
        {
            seg_changeDatatype<SegPrecisionTYPE>(InputImage);
        }
        SegPrecisionTYPE * InputImagePtr = static_cast<SegPrecisionTYPE *>(InputImage->data);
        ImageSize * CurrSize = new ImageSize [1]();
        CurrSize->numel=(long)(InputImage->nx*InputImage->ny*InputImage->nz);
        CurrSize->xsize=InputImage->nx;
        CurrSize->ysize=InputImage->ny;
        CurrSize->zsize=InputImage->nz;
        CurrSize->usize=(InputImage->nu>1)?InputImage->nu:1;
        CurrSize->tsize=(InputImage->nt>1)?InputImage->nt:1;
	
	nifti_image * InputLesionMask=nifti_image_read(filename_lesion_mask,true);
        if(InputLesionMask == NULL)
        {
            fprintf(stderr,"* Error when reading the input lesion mask\n");
            return 1;
        }
        if(InputLesionMask->datatype!=NIFTI_TYPE_FLOAT32)
        {
            seg_changeDatatype<SegPrecisionTYPE>(InputLesionMask);
        }
	SegPrecisionTYPE * MaskPtr =  static_cast<SegPrecisionTYPE *>(InputLesionMask->data);
	
	// We read the brain mask only if it is given
	nifti_image * InputMask=NULL;
	if(filename_mask!=NULL) {
	     InputMask=nifti_image_read(filename_mask,true);
             if(InputMask == NULL)
             {
                 fprintf(stderr,"* Error when reading the input mask\n");
                 return 1;
             }
             if(InputMask->datatype!=NIFTI_TYPE_FLOAT32)
             {
                 seg_changeDatatype<SegPrecisionTYPE>(InputMask);
             }
             if(InputMask->nx!=InputImage->nx || InputMask->ny!=InputImage->ny || InputMask->nz!=InputImage->nz || InputMask->nt!=InputImage->nt) {
                 cerr << "ERROR: Mask "<< filename_mask << " is the wrong size  -  original = ( "<<InputImage->nx<<","
                      <<InputImage->ny<<","<<InputImage->nz<<","<<InputImage->nt<<","<<InputImage->nu<<" )  New image = ( "<<InputMask->nx<<","
                      <<InputMask->ny<<","<<InputMask->nz<<","<<InputMask->nt<<","<<InputMask->nu<<" )"<<endl;
		 return 0;
	     }
	}

        // We show the parameters
        if(verbose)
        {
            cout<<"Parameters values:"<<endl;
	    if(method==0) {
		cout<<"PARAM[method]=Prados et al. (Neuroimage 2016)"<<endl;
		if(patch2D) cout<<"PARAM[dimensionality]=2D"<<endl;
		else cout<<"PARAM[dimensionality]=3D"<<endl;
	    }
	    else {
		cout<<"PARAM[method]=Guizard et al. (FIN 2015)"<<endl;
		if(patch2D) cout<<"PARAM[dimensionality]=2D option NOT available for this method"<<endl;
		else cout<<"PARAM[dimensionality]=3D"<<endl;
	    }
	    cout<<"PARAM[input]="<<filename_in<<endl;
            cout<<"PARAM[les_mask]="<<filename_lesion_mask<<endl;
            if(filename_mask!=NULL) cout<<"PARAM[mask]="<<filename_mask<<endl;
	    cout<<"PARAM[output]="<<filename_out<<endl;
            cout<<"PARAM[dil]="<<dil<<endl;
            cout<<"PARAM[smo]="<<smo<<endl;
            cout<<"PARAM[match]="<<match<<endl;
            cout<<"PARAM[search]="<<expanding<<endl;
            cout<<"PARAM[size]="<<size<<endl;
	    cout<<"PARAM[card_weight_factor]="<<k<<endl;
            cout<<"PARAM[debug]="<<debug<<endl;
            cout<<"PARAM[odt/idt]="<<datatypeoutput<<"/"<<InputImage->datatype<<endl;
	    #if defined (_OPENMP)
	    cout<<"PARAM[threads]="<<omp_get_max_threads()<<endl;
	    #endif
            cout<<"PARAM[v]="<<verbose<<endl;         
        }        

        // We include the small empty voxels in the middle of the lesions at the mask
        if(verbose) cout <<"Filling gaps in the lesion mask" << endl;
        Close_Forground_ConnectComp<float,float>(MaskPtr,MaskPtr,CurrSize);
        if(dil>0) {
            if(verbose) cout <<"Dillating "<<dil<<" times the lesion mask" << endl;
            Dillate(MaskPtr,dil,CurrSize);
        }        

        // We are ready for filling lesion
        if(InputLesionMask->nx==InputImage->nx && InputLesionMask->ny==InputImage->ny && InputLesionMask->nz==InputImage->nz && InputLesionMask->nt==InputImage->nt)
        {
            if(verbose) cout<<"Fill lesions process starts"<<endl;
            if(method==0) {
		seg_fill_lesions <SegPrecisionTYPE> *fillLesions=new seg_fill_lesions<SegPrecisionTYPE>();
		fillLesions->setVerbose(verbose);
		fillLesions->setDebug(debug);
		fillLesions->setInputImage(InputImage);
		fillLesions->setInputLesionMask(InputLesionMask);
		if(filename_mask!=NULL) fillLesions->setInputMask(InputMask);
		fillLesions->setPatchSearchAreaSize(size);
		fillLesions->setPatchPercentage(match);
		fillLesions->setExpandingPercentage(expanding);
		fillLesions->setSmoothing(smo);
		fillLesions->setDimensionality(patch2D);
		fillLesions->setK(k);
		fillLesions->runIt();
	    }
	    else {
		if(patch2D) {
		    cout<<"PARAM[dimensionality]=2D option NOT available for this method"<<endl;
		    return -1;
		}
		seg_fill_lesions_other <SegPrecisionTYPE> *fillLesions=new seg_fill_lesions_other<SegPrecisionTYPE>();
		fillLesions->setVerbose(verbose);
		fillLesions->setDebug(debug);
		fillLesions->setInputImage(InputImage);
		fillLesions->setInputLesionMask(InputLesionMask);
		fillLesions->setPatchSize(4); //Paper says 9, then 4*2+1=9
		fillLesions->setSearchArea(10); // Paper says radius 10.
		fillLesions->runIt();
	    }
        }
        else
        {
            cout << "ERROR: Lesion mask "<< filename_lesion_mask << " is the wrong size  -  original = ( "<<InputImage->nx<<","
                 <<InputImage->ny<<","<<InputImage->nz<<","<<InputImage->nt<<","<<InputImage->nu<<" )  New image = ( "<<InputLesionMask->nx<<","
                <<InputLesionMask->ny<<","<<InputLesionMask->nz<<","<<InputLesionMask->nt<<","<<InputLesionMask->nu<<" )"<<endl;
	    return -1;
        }

        if(filename_out.find(string(".nii"))>0 || filename_out.find(string(".img")) || filename_out.find(string(".hdr"))>0)
        {
            // saving output
            if(verbose) cout<< "Saving results at: "<<filename_out<<endl;
            nifti_image * OutputImage = nifti_copy_nim_info(InputImage);
            OutputImage->datatype=datatypeoutput;
            nifti_set_filenames(OutputImage,filename_out.c_str(),0,0);
            OutputImage->dim[1]=InputImage->nx;
            OutputImage->dim[2]=InputImage->ny;
            OutputImage->dim[3]=InputImage->nz;
            OutputImage->dim[4]=OutputImage->nt=InputImage->nt;
            OutputImage->dim[5]=OutputImage->nu=InputImage->nu;
            OutputImage->dim[6]=OutputImage->nv=1;
            OutputImage->dim[7]=OutputImage->nw=1;
            OutputImage->dim[0]=3;
            OutputImage->dim[0]=(OutputImage->dim[4]>1?4:OutputImage->dim[0]);
            OutputImage->dim[0]=(OutputImage->dim[5]>1?5:OutputImage->dim[0]);
            OutputImage->dim[0]=(OutputImage->dim[6]>1?6:OutputImage->dim[0]);
            OutputImage->dim[0]=(OutputImage->dim[7]>1?7:OutputImage->dim[0]);

            if(verbose)
            {
                cout << "Output Dim = [ ";
                for(long i=0; i<8; i++)
                {
                    cout<<(float)OutputImage->dim[i];
                    if(i<7)
                    {
                        cout<<" , ";
                    }
                }
                cout<<" ] "<<endl;
                flush(cout);
            }
            long numElem=1;
            for(long i=1; i<8; i++)
            {
                if(OutputImage->dim[i]>0) {
                        numElem*=OutputImage->dim[i];
                }
            }
            nifti_update_dims_from_array(OutputImage);
            nifti_datatype_sizes(OutputImage->datatype,&OutputImage->nbyper,&OutputImage->swapsize);
            if(datatypeoutput==NIFTI_TYPE_UINT8)
            {
                OutputImage->data = (void *) calloc(numElem, sizeof(unsigned char));
                unsigned char * OutputImagePtr = static_cast<unsigned char *>(OutputImage->data);
                for(long i=0; i<(long)(numElem); i++)
                {
                    OutputImagePtr[i]=(unsigned char)round(InputImagePtr[i]);
                }
            }
            else if(datatypeoutput==NIFTI_TYPE_UINT16)
            {
                OutputImage->data = (void *) calloc(OutputImage->nvox, sizeof(unsigned short));
                unsigned short * OutputImagePtr = static_cast<unsigned short *>(OutputImage->data);
                for(long i=0; i<(long)(numElem); i++)
                {
                    OutputImagePtr[i]=(unsigned short)round(InputImagePtr[i]);
                }
            }
            else if(datatypeoutput==NIFTI_TYPE_UINT32)
            {
                OutputImage->data = (void *) calloc(numElem, sizeof(unsigned int));
                unsigned int * OutputImagePtr = static_cast<unsigned int *>(OutputImage->data);
                for(long i=0; i<(long)(numElem); i++)
                {
                    OutputImagePtr[i]=(unsigned int)round(InputImagePtr[i]);
                }
            }
            else if(datatypeoutput==NIFTI_TYPE_INT8)
            {
                OutputImage->data = (void *) calloc(numElem, sizeof(char));
                char * OutputImagePtr = static_cast<char *>(OutputImage->data);
                for(long i=0; i<(long)(numElem); i++)
                {
                    OutputImagePtr[i]=(char)round(InputImagePtr[i]);
                }
            }
            else if(datatypeoutput==NIFTI_TYPE_INT16)
            {
                OutputImage->data = (void *) calloc(numElem, sizeof(short));
                short * OutputImagePtr = static_cast<short *>(OutputImage->data);
                for(long i=0; i<(long)(numElem); i++)
                {
                    OutputImagePtr[i]=(short)round(InputImagePtr[i]);
                }
            }
            else if(datatypeoutput==NIFTI_TYPE_INT32)
            {
                OutputImage->data = (void *) calloc(numElem, sizeof(int));
                int * OutputImagePtr = static_cast<int *>(OutputImage->data);
                for(long i=0; i<(long)(numElem); i++)
                {
                    OutputImagePtr[i]=(int)round(InputImagePtr[i]);
                }
            }
            else if(datatypeoutput==NIFTI_TYPE_FLOAT32)
            {
                OutputImage->data = (void *) calloc(numElem, sizeof(float));
                float * OutputImagePtr = static_cast<float *>(OutputImage->data);
                for(long i=0; i<(long)(numElem); i++)
                {
                    OutputImagePtr[i]=(float)InputImagePtr[i];
                }
            }
            else if(datatypeoutput==NIFTI_TYPE_FLOAT64)
            {
                OutputImage->data = (void *) calloc(numElem, sizeof(double));
                double * OutputImagePtr = static_cast<double *>(OutputImage->data);
                for(long i=0; i<(long)(numElem); i++)
                {
                    OutputImagePtr[i]=(double)round(InputImagePtr[i]);
                }
            }
            nifti_image_write(OutputImage);
            nifti_image_free(OutputImage);
        }
        nifti_image_free(InputLesionMask);
        nifti_image_free(InputImage);
        time_t end;
        time( &end );
        int minutes = (int)floorf(float(end-start)/60.0f);
        int seconds = (int)(end-start - 60*minutes);
        if(verbose) cout<<"Lesions Filled in "<<minutes<<" min "<<seconds<<" sec"<<endl;
        if(verbose) cout<<"Have a good day !"<<endl;
    }

    catch(std::exception & e)
    {
        std::cerr << "Standard exception: " << e.what() << std::endl;
    }

    catch(...)
    {
        std::cerr << "Unhandled Exception: Something went wrong! Please report the error to f.carrasco"<<(char) 64<<"ucl.ac.uk" << std::endl;
    }
    return 0;
}

