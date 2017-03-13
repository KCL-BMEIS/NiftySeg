/**
 * @file seg_PatchMatch.cpp
 * @author Ferran Prados
 * @date 16/06/2016
 *
 * Copyright (c) 2016, University College London. All rights reserved.
 * Centre for Medical Image Computing (CMIC)
 * See the LICENSE.txt file in the nifty_seg root folder
 *
 */

#include <iostream>
#include <string>
#include <time.h>
#include "_seg_common.h"
#include "_seg_tools.h"
#include "_seg_PatchMatch.h"

using namespace std;

void Usage(char *exec)
{
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    printf("Usage:\t%s -i <input> -m <input_mask> -db <database> -o <output_result> [options].\n\n",exec);
    printf("\t* * Considerations * *\n\n");
    printf("\t\tFor an extended help, please read the NiftySeg wiki page at: http://cmictig.cs.ucl.ac.uk/wiki/index.php/NiftySeg\n");
    printf("\t\tThe database file is a text file and in each line we have a template file, a mask with the search region to consider and a file with the label to propagate.\n");
    printf("\t\tInput image, input mask, template images from database and masks from database must have the same 4D resolution (same number of XxYxZ voxels, modalities and/or time-points).\n");
    printf("\t\tLabel files from database must have the same 3D resolution (XxYxZ voxels) than input image but can have different number of volumes than the input image allowing to propagate multiple labels in the same execution.\n");
    printf("\t\tThe file should be as:\n\n");
    printf("\t\t/absolute_path/image_filename_1.[nii/nii.gz/hdr][DELIMITER]/absolute_path/mask_filename_1.[nii/nii.gz/hdr]/absolute_path/label_filename_1.[nii/nii.gz/hdr]\n");
    printf("\t\t/absolute_path/image_filename_2.[nii/nii.gz/hdr][DELIMITER]/absolute_path/mask_filename_2.[nii/nii.gz/hdr]/absolute_path/label_filename_2.[nii/nii.gz/hdr]\n");
    printf("\t\t...\n");
    printf("\t\t/absolute_path/image_filename_N.[nii/nii.gz/hdr][DELIMITER]/absolute_path/mask_filename_N.[nii/nii.gz/hdr]/absolute_path/label_filename_N.[nii/nii.gz/hdr]\n\n");
    printf("\t\tDELIMITER can be whitespace, ',', ';' or tab.\n\n");
    printf("\t* * Optional parameters * *\n\n");
    printf("\t-size\t\t\tPatch size, #voxels (by default 5).\n");
    printf("\t-cs\t\t\tConstrained search area size, number of times bigger than the patchsize (by default 4).\n");
    printf("\t-match\t\t\tNumber of better matching (by default 10).\n");
    printf("\t-pm\t\t\tNumber of patchmatch executions (by default 10). It should be equal or bigger than the number of better matching.\n");
    printf("\t-it\t\t\tNumber of iterations for the patchmatch algorithm (by default 5).\n");
    printf("\t-fill\t\t\tIt applies label fusion at all the voxels independtly if they are inside the input mask or not.\n");
    printf("\t-dist\t\t\tUsed distance (by default 0, SSD=0, LNCC=1).\n");
    printf("\t-debug\t\t\tSave all intermidium files (by default OFF).\n");
    printf("\t-odt <datatype> \tSet output <datatype> (char, short, int, uchar, ushort, uint, float, double).\n");
    printf("\t-v\t\t\tVerbose (by default OFF).\n");
    #if defined (_OPENMP)
    printf("\t-omp <int>\t\tNumber of openmp threads [%d]\n",omp_get_max_threads());
    #endif
    #ifdef _GIT_HASH
    printf("\t--version\t\tPrint current source code git hash key and exit\n\t\t\t\t(%s)\n",_GIT_HASH);
    #endif
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
        bool debug=0;
        int size=5;
        int cs=4;
        int better_match=10;
        int pm_threads=10;
        int it=5;
        int distance=0;
	bool filling=false;
        time_t start;
        time(&start);
        vector<string> imageFiles;
        vector<string> maskFiles;
        vector<string> outputFiles;
        char * filename_in=NULL;
        char * filename_mask=NULL;
        char * filename_database=NULL;
        string filename_out="";
	bool ok_input_file=false;
	bool ok_mask_file=false;
	bool ok_database_file=false;
	bool ok_output_file=false;
	
        set_new_handler(no_memory);
        if (argc <= 8)
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

        // We read extra parameters for the algorithm
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
	    else  if(strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "-mask") == 0)
            {
		filename_mask=argv[++i];
		ok_mask_file=true;
	    }
	    else  if(strcmp(argv[i], "-db") == 0 || strcmp(argv[i], "-database") == 0)
            {
		filename_database=argv[++i];
		ok_database_file=true;
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
                if(size%2!=1 || size<3) {
                    cerr << "ERROR: Wrong patch size "<< size << ", it has to be an odd number equal or bigger than 3."<<endl;
                    return -1;
                }
            }
            else if(strcmp(argv[i], "-debug") == 0)
            {
                debug=1;
            }
            else if(strcmp(argv[i], "-fill") == 0)
            {
                filling=true;
            }
            else if(strcmp(argv[i], "-cs") == 0)
            {
                string parser=argv[++i];
                cs=strtod(parser.c_str(),NULL);
                if(cs<2) {
                    cerr << "ERROR: Wrong number of the contrained search area size "<< cs << " is too small. At least 2."<<endl;
                    return -1;
                }
            }
            else if(strcmp(argv[i], "-match") == 0)
            {
                string parser=argv[++i];
                better_match=strtod(parser.c_str(),NULL);
                if(better_match<1) {
                    cerr << "ERROR: Wrong number of better match "<< better_match << " is not possible"<<endl;
                    return -1;
                }
            }
            else if(strcmp(argv[i], "-pm") == 0)
            {
                string parser=argv[++i];
                pm_threads=strtod(parser.c_str(),NULL);
                if(pm_threads<1) {
                    cerr << "ERROR: Wrong number of Patchmatch executions "<< pm_threads << " is not possible"<<endl;
                    return -1;
                }
            }
            else if(strcmp(argv[i], "-it") == 0)
            {
                string parser=argv[++i];
                it=strtod(parser.c_str(),NULL);
                if(it<1) {
                    cerr << "ERROR: Wrong number of iterations "<< it << " is not possible"<<endl;
                    return -1;
                }
            }
            else if(strcmp(argv[i], "-dist") == 0)
            {
                string parser=argv[++i];
                distance=(int)round(strtod(parser.c_str(),NULL));
                if(distance>1 || distance<0) {
                    cerr << "ERROR: Wrong distance "<< distance << " is unknown"<<endl;
                    return -1;
                }
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
                    cerr << "ERROR: Datatype "<< parser << " is unknown"<<endl;
                    i=argc;
		    return -1;
                }
            }
	    #ifdef _GIT_HASH
            else if( strcmp(argv[i], "--version")==0)
            {
                printf("%s\n",_GIT_HASH);
                return 0;
            }
	    #endif
            else
            {
                cerr << "Option "<< string(argv[i]) << " unkown"<<endl;
                i=argc;
                return -1;
            }
        }
	if (!ok_input_file) {
	    cerr<<"* Error input file not provided -i <filename>"<<endl;
	    Usage(argv[0]);
            return -1;
	}
	if (!ok_mask_file) {
	    cerr<<"* Error input file not provided -m <filename>"<<endl;
	    Usage(argv[0]);
            return -1;
	}
	if (!ok_database_file) {
	    cerr<<"* Error database file not provided -db <filename>"<<endl;
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
        nifti_image * InputMask=nifti_image_read(filename_mask,true);
        if(InputMask == NULL)
        {
            fprintf(stderr,"* Error when reading the input mask\n");
            return 1;
        }
        if(InputImage->datatype!=NIFTI_TYPE_FLOAT32)
        {
            seg_changeDatatype<SegPrecisionTYPE>(InputMask);
        }

        ImageSize * CurrSize = new ImageSize [1]();
        CurrSize->numel=(long)(InputImage->nx*InputImage->ny*InputImage->nz);
        CurrSize->xsize=InputImage->nx;
        CurrSize->ysize=InputImage->ny;
        CurrSize->zsize=InputImage->nz;
        CurrSize->usize=(InputImage->nu>1)?InputImage->nu:1;
        CurrSize->tsize=(InputImage->nt>1)?InputImage->nt:1;
	
        if(pm_threads<better_match) {
            cerr<< "ERROR: the number of patchmatch executions can't be less than the expected matching better results."<<endl;
            exit(-1);
        }
        if(!(InputMask->nx==InputImage->nx && InputMask->ny==InputImage->ny && InputMask->nz==InputImage->nz))
        {
            cerr << "ERROR: input mask "<< filename_mask << " has a different size  respect to input  image = ( "<<InputImage->nx<<","
                 <<InputImage->ny<<","<<InputImage->nz<<","<<" )  input mask = ( "<<InputMask->nx<<","
                <<InputMask->ny<<","<<InputMask->nz<<" )"<<endl;
            exit(-1);
        }

        // We show the parameters
        if(verbose)
        {
            cout<<"Parameters values:"<<endl;
            cout<<"PARAM[input]="<<filename_in<<endl;
            cout<<"PARAM[mask]="<<filename_mask<<endl;
            cout<<"PARAM[database]="<<filename_database<<endl;
            cout<<"PARAM[output]="<<filename_out<<endl;
            cout<<"PARAM[distance]="<<distance<<endl;
            cout<<"PARAM[patchmatch_executions]="<<pm_threads<<endl;
            cout<<"PARAM[iteratons]="<<it<<endl;
            cout<<"PARAM[better_match]="<<better_match<<endl;
            cout<<"PARAM[cs_size]="<<cs<<endl;
            cout<<"PARAM[patch_size]="<<size<<endl;
            if (filling) cout<<"PARAM[filling]=True"<<endl;
	    else cout<<"PARAM[filling]=False"<<endl;
	    cout<<"PARAM[debug]="<<debug<<endl;
            cout<<"PARAM[odt/idt]="<<datatypeoutput<<"/"<<InputImage->datatype<<endl;
	    cout<<"PARAM[v]="<<verbose<<endl;         
        }
        float *OutputResult=NULL;
        // We are ready for lesion detection

        if(verbose) cout<<"Loading database "<<filename_database<<endl;
        std::ifstream file(filename_database);
        std::string str,files,delimiter;
        std::getline(file, str);

        if(verbose) cout<<"Detecting delimiter: ";

        delimiter = " ";
        if(str.find(delimiter)!= std::string::npos)  delimiter = " ";
        else {
            delimiter = ",";
            if(str.find(delimiter)!= std::string::npos)  delimiter = ",";
            else {
                delimiter = ";";
                if(str.find(delimiter)!= std::string::npos)  delimiter = ";";
                else {
                    delimiter = "\t";
                    if(str.find(delimiter)!= std::string::npos)  delimiter = "\t";
                }
            }
        }
        if(verbose) cout<<"Delimiter detected: "<<delimiter<<endl;

        // We start again to read the file
        std::ifstream file2(filename_database);
        if(verbose) cout<<"Reading file:"<<endl;
        while (std::getline(file2, str))
        {
            files = str.substr(0, str.find(delimiter));
            if(verbose) std::cout<<files<<"|";
            imageFiles.push_back(files);
            str.erase(0, str.find(delimiter) + delimiter.length());


            files = str.substr(0, str.find(delimiter));
            if(verbose) std::cout<<files<<"|"<<endl;
            maskFiles.push_back(files);
            str.erase(0, str.find(delimiter) + delimiter.length());

            files = str.substr(0, str.find(delimiter));
            if(verbose) std::cout<<files<<"|"<<endl;
            outputFiles.push_back(files);
            str.erase(0, str.find(delimiter) + delimiter.length());
        }
        if(maskFiles.size()<=0 || imageFiles.size()<=0 || outputFiles.size()<=0) {
            cout<<"The database is empty, 0 files detected"<<endl;
            return 0;
        }

        if(verbose) cout<<"PatchMatch process starts"<<endl;
        if(verbose) cout<<"PatchMatching..."<<endl;
        seg_PatchMatch <SegPrecisionTYPE> *patchmatch=new seg_PatchMatch<SegPrecisionTYPE>();
        patchmatch->setVerbose(verbose);
        patchmatch->setDebug(debug);
        patchmatch->setInputImage(InputImage);
        patchmatch->setInputMask(InputMask);
        patchmatch->setInputImageDatabase(imageFiles);
        patchmatch->setInputMaskDatabase(maskFiles);
        patchmatch->setOutputFilesDatabase(outputFiles);
        patchmatch->setPatchMatchExecutions(pm_threads);
        patchmatch->setPatchMatchIterations(it);
        patchmatch->setConstrainedSearchAreaSize(cs);
        patchmatch->setBetterMatch(better_match);
        patchmatch->setPatchSize(size);
        patchmatch->setDistance(distance);
	patchmatch->setFilling(filling);
        patchmatch->runIt();
        OutputResult= patchmatch->getOutputResult();

        if(filename_out.find(string(".nii"))>0 || filename_out.find(string(".img")) || filename_out.find(string(".hdr"))>0)
        {
            // saving output
            if(verbose) cout<< "Saving results at: "<<filename_out<<" "<<outputFiles[0].c_str()<<endl;
            nifti_image * OutputImage = nifti_image_read(outputFiles[0].c_str(),true);
            OutputImage->datatype=datatypeoutput;
            nifti_set_filenames(OutputImage,filename_out.c_str(),0,0);
            float max=std::numeric_limits<float>::min();
            float min=std::numeric_limits<float>::max();
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
            if(verbose)
            {
                cout <<"Converting data "<<numElem<<endl;
            }
            if(datatypeoutput==NIFTI_TYPE_UINT8)
            {
                OutputImage->data = (void *) calloc(numElem, sizeof(unsigned char));
                unsigned char * OutputImagePtr = static_cast<unsigned char *>(OutputImage->data);
                for(long i=0; i<numElem; i++)
                {
                    max=max>OutputResult[i]?max:OutputResult[i];
                    min=min<OutputResult[i]?min:OutputResult[i];
                    OutputImagePtr[i]=(unsigned char)round(OutputResult[i]);
                }
            }
            else if(datatypeoutput==NIFTI_TYPE_UINT16)
            {
                OutputImage->data = (void *) calloc(numElem, sizeof(unsigned short));
                unsigned short * OutputImagePtr = static_cast<unsigned short *>(OutputImage->data);
                for(long i=0; i<numElem; i++)
                {
                    max=max>OutputResult[i]?max:OutputResult[i];
                    min=min<OutputResult[i]?min:OutputResult[i];
                    OutputImagePtr[i]=(unsigned short)round(OutputResult[i]);
                }
            }
            else if(datatypeoutput==NIFTI_TYPE_UINT32)
            {
                OutputImage->data = (void *) calloc(numElem, sizeof(unsigned int));
                unsigned int * OutputImagePtr = static_cast<unsigned int *>(OutputImage->data);
                for(long i=0; i<numElem; i++)
                {
                    max=max>OutputResult[i]?max:OutputResult[i];
                    min=min<OutputResult[i]?min:OutputResult[i];
                    OutputImagePtr[i]=(unsigned int)round(OutputResult[i]);
                }
            }
            else if(datatypeoutput==NIFTI_TYPE_INT8)
            {
                OutputImage->data = (void *) calloc(numElem, sizeof(char));
                char * OutputImagePtr = static_cast<char *>(OutputImage->data);
                for(long i=0; i<numElem; i++)
                {
                    max=max>OutputResult[i]?max:OutputResult[i];
                    min=min<OutputResult[i]?min:OutputResult[i];
                    OutputImagePtr[i]=(char)round(OutputResult[i]);
                }
            }
            else if(datatypeoutput==NIFTI_TYPE_INT16)
            {
                OutputImage->data = (void *) calloc(numElem, sizeof(short));
                short * OutputImagePtr = static_cast<short *>(OutputImage->data);
                for(long i=0; i<numElem; i++)
                {
                    max=max>OutputResult[i]?max:OutputResult[i];
                    min=min<OutputResult[i]?min:OutputResult[i];
                    OutputImagePtr[i]=(short)round(OutputResult[i]);
                }
            }
            else if(datatypeoutput==NIFTI_TYPE_INT32)
            {
                OutputImage->data = (void *) calloc(numElem, sizeof(int));
                int * OutputImagePtr = static_cast<int *>(OutputImage->data);
                for(long i=0; i<numElem; i++)
                {
                    max=max>OutputResult[i]?max:OutputResult[i];
                    min=min<OutputResult[i]?min:OutputResult[i];
                    OutputImagePtr[i]=(int)round(OutputResult[i]);
                }
            }
            else if(datatypeoutput==NIFTI_TYPE_FLOAT32)
            {
                OutputImage->data = (void *) calloc(numElem, sizeof(float));
                float * OutputImagePtr = static_cast<float *>(OutputImage->data);
                for(long i=0; i<numElem; i++)
                {
                    max=max>OutputResult[i]?max:OutputResult[i];
                    min=min<OutputResult[i]?min:OutputResult[i];
                    OutputImagePtr[i]=(float)OutputResult[i];
                }
            }
            else if(datatypeoutput==NIFTI_TYPE_FLOAT64)
            {
                OutputImage->data = (void *) calloc(numElem, sizeof(double));
                double * OutputImagePtr = static_cast<double *>(OutputImage->data);
                for(long i=0; i<numElem; i++)
                {
                    max=max>OutputResult[i]?max:OutputResult[i];
                    min=min<OutputResult[i]?min:OutputResult[i];
                    OutputImagePtr[i]=(double)round(OutputResult[i]);
                }
            }
            OutputImage->cal_max=max;
            OutputImage->cal_min=min;
            OutputImage->scl_inter=0;
            OutputImage->scl_slope=1;
            if(verbose)
            {
                cout <<"Writing image"<<endl;
            }
            nifti_image_write(OutputImage);
            nifti_image_free(OutputImage);
        }
        if(verbose)
        {
            cout <<"Free"<<endl;
        }
        nifti_image_free(InputImage);
        time_t end;
        time( &end );
        int minutes = (int)floorf(float(end-start)/60.0f);
        int seconds = (int)(end-start - 60*minutes);
        if(verbose) cout<<"PatchMatch done in "<<minutes<<" min "<<seconds<<" sec"<<endl;
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

