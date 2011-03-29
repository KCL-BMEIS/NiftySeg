/*
 *  _reg_tools_gpu.cu
 *  
 *
 *  Created by Marc Modat and Pankaj Daga on 24/03/2009.
 *  Copyright (c) 2009, University College London. All rights reserved.
 *  Centre for Medical Image Computing (CMIC)
 *  See the LICENSE.txt file in the nifty_reg root folder
 *
 */

#ifndef _REG_TOOLS_GPU_CU
#define _REG_TOOLS_GPU_CU

#include "_reg_blocksize_gpu.h"
#include "_reg_tools_kernels.cu"


void reg_voxelCentric2NodeCentric_gpu(	nifti_image *targetImage,
					nifti_image *controlPointImage,
					float4 **voxelNMIGradientArray_d,
					float4 **nodeNMIGradientArray_d,
					float weight)
{
	const int nodeNumber = controlPointImage->nx * controlPointImage->ny * controlPointImage->nz;
	const int voxelNumber = targetImage->nx * targetImage->ny * targetImage->nz;
	const int3 targetImageDim = make_int3(targetImage->nx, targetImage->ny, targetImage->nz);
	const int3 gridSize = make_int3(controlPointImage->nx, controlPointImage->ny, controlPointImage->nz);
	const float3 voxelNodeRatio_h = make_float3(
		controlPointImage->dx / targetImage->dx,
		controlPointImage->dy / targetImage->dy,
		controlPointImage->dz / targetImage->dz);

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_NodeNumber,&nodeNumber,sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_TargetImageDim,&targetImageDim,sizeof(int3)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_ControlPointImageDim,&gridSize,sizeof(int3)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_VoxelNodeRatio,&voxelNodeRatio_h,sizeof(float3)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_Weight,&weight,sizeof(float)));

	CUDA_SAFE_CALL(cudaBindTexture(0, gradientImageTexture, *voxelNMIGradientArray_d, voxelNumber*sizeof(float4)));

	const unsigned int Grid_reg_voxelCentric2NodeCentric = (unsigned int)ceil((float)nodeNumber/(float)Block_reg_voxelCentric2NodeCentric);
	dim3 B1(Block_reg_voxelCentric2NodeCentric,1,1);
	dim3 G1(Grid_reg_voxelCentric2NodeCentric,1,1);

	reg_voxelCentric2NodeCentric_kernel <<< G1, B1 >>> (*nodeNMIGradientArray_d);
	CUDA_SAFE_CALL(cudaThreadSynchronize());
#ifndef NDEBUG
    printf("[NiftyReg CUDA DEBUG] reg_voxelCentric2NodeCentric_gpu kernel: %s - Grid size [%i %i %i] - Block size [%i %i %i]\n",
	       cudaGetErrorString(cudaGetLastError()),G1.x,G1.y,G1.z,B1.x,B1.y,B1.z);
#endif
}

void reg_convertNMIGradientFromVoxelToRealSpace_gpu(	mat44 *sourceMatrix_xyz,
							nifti_image *controlPointImage,
							float4 **nodeNMIGradientArray_d)
{
	const int nodeNumber = controlPointImage->nx * controlPointImage->ny * controlPointImage->nz;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_NodeNumber,&nodeNumber,sizeof(int)));

    float4 *matrix_h;CUDA_SAFE_CALL(cudaMallocHost(&matrix_h, 3*sizeof(float4)));
	matrix_h[0] = make_float4(sourceMatrix_xyz->m[0][0], sourceMatrix_xyz->m[0][1], sourceMatrix_xyz->m[0][2], sourceMatrix_xyz->m[0][3]);
	matrix_h[1] = make_float4(sourceMatrix_xyz->m[1][0], sourceMatrix_xyz->m[1][1], sourceMatrix_xyz->m[1][2], sourceMatrix_xyz->m[1][3]);
	matrix_h[2] = make_float4(sourceMatrix_xyz->m[2][0], sourceMatrix_xyz->m[2][1], sourceMatrix_xyz->m[2][2], sourceMatrix_xyz->m[2][3]);
	float4 *matrix_d;
    CUDA_SAFE_CALL(cudaMalloc(&matrix_d, 3*sizeof(float4)));
	CUDA_SAFE_CALL(cudaMemcpy(matrix_d, matrix_h, 3*sizeof(float4), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaFreeHost((void *)matrix_h));
	CUDA_SAFE_CALL(cudaBindTexture(0, matrixTexture, matrix_d, 3*sizeof(float4)));
	
	const unsigned int Grid_reg_convertNMIGradientFromVoxelToRealSpace =
		(unsigned int)ceil((float)nodeNumber/(float)Block_reg_convertNMIGradientFromVoxelToRealSpace);
	dim3 B1(Grid_reg_convertNMIGradientFromVoxelToRealSpace,1,1);
	dim3 G1(Block_reg_convertNMIGradientFromVoxelToRealSpace,1,1);

	_reg_convertNMIGradientFromVoxelToRealSpace_kernel <<< G1, B1 >>> (*nodeNMIGradientArray_d);
	CUDA_SAFE_CALL(cudaThreadSynchronize());
#ifndef NDEBUG
    printf("[NiftyReg CUDA DEBUG] reg_convertNMIGradientFromVoxelToRealSpace: %s - Grid size [%i %i %i] - Block size [%i %i %i]\n",
	       cudaGetErrorString(cudaGetLastError()),G1.x,G1.y,G1.z,B1.x,B1.y,B1.z);
#endif
	CUDA_SAFE_CALL(cudaFree(matrix_d));
}


void reg_initialiseConjugateGradient(	float4 **nodeNMIGradientArray_d,
					float4 **conjugateG_d,
					float4 **conjugateH_d,
					int nodeNumber)
{
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_NodeNumber,&nodeNumber,sizeof(int)));
	CUDA_SAFE_CALL(cudaBindTexture(0, gradientImageTexture, *nodeNMIGradientArray_d, nodeNumber*sizeof(float4)));

	const unsigned int Grid_reg_initialiseConjugateGradient =
		(unsigned int)ceil((float)nodeNumber/(float)Block_reg_initialiseConjugateGradient);
	dim3 B1(Grid_reg_initialiseConjugateGradient,1,1);
	dim3 G1(Block_reg_initialiseConjugateGradient,1,1);

	reg_initialiseConjugateGradient_kernel <<< G1, B1 >>> (*conjugateG_d);
	CUDA_SAFE_CALL(cudaThreadSynchronize());
#ifndef NDEBUG
    printf("[NiftyReg CUDA DEBUG] reg_initialiseConjugateGradient: %s - Grid size [%i %i %i] - Block size [%i %i %i]\n",
	       cudaGetErrorString(cudaGetLastError()),G1.x,G1.y,G1.z,B1.x,B1.y,B1.z);
#endif
	CUDA_SAFE_CALL(cudaMemcpy(*conjugateH_d, *conjugateG_d, nodeNumber*sizeof(float4), cudaMemcpyDeviceToDevice));
}

void reg_GetConjugateGradient(	float4 **nodeNMIGradientArray_d,
				float4 **conjugateG_d,
				float4 **conjugateH_d,
				int nodeNumber)
{
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_NodeNumber,&nodeNumber,sizeof(int)));
	CUDA_SAFE_CALL(cudaBindTexture(0, conjugateGTexture, *conjugateG_d, nodeNumber*sizeof(float4)));
	CUDA_SAFE_CALL(cudaBindTexture(0, conjugateHTexture, *conjugateH_d, nodeNumber*sizeof(float4)));
	CUDA_SAFE_CALL(cudaBindTexture(0, gradientImageTexture, *nodeNMIGradientArray_d, nodeNumber*sizeof(float4)));

	// gam = sum((grad+g)*grad)/sum(HxG);
	const unsigned int Grid_reg_GetConjugateGradient1 = (unsigned int)ceil((float)nodeNumber/(float)Block_reg_GetConjugateGradient1);
	dim3 B1(Block_reg_GetConjugateGradient1,1,1);
	dim3 G1(Grid_reg_GetConjugateGradient1,1,1);

	float2 *sum_d;
    CUDA_SAFE_CALL(cudaMalloc(&sum_d, nodeNumber*sizeof(float2)));
	reg_GetConjugateGradient1_kernel <<< G1, B1 >>> (sum_d);
	CUDA_SAFE_CALL(cudaThreadSynchronize());
#ifndef NDEBUG
    printf("[NiftyReg CUDA DEBUG] reg_GetConjugateGradient1 kernel: %s - Grid size [%i %i %i] - Block size [%i %i %i]\n",
	       cudaGetErrorString(cudaGetLastError()),G1.x,G1.y,G1.z,B1.x,B1.y,B1.z);
#endif
    float2 *sum_h;CUDA_SAFE_CALL(cudaMallocHost(&sum_h, nodeNumber*sizeof(float2)));
	CUDA_SAFE_CALL(cudaMemcpy(sum_h,sum_d, nodeNumber*sizeof(float2),cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaFree(sum_d));
	double dgg = 0.0;
	double gg = 0.0;
	for(int i=0; i<nodeNumber; i++){
		dgg += sum_h[i].x;
		gg += sum_h[i].y;
	}
	float gam = (float)(dgg / gg);
	CUDA_SAFE_CALL(cudaFreeHost((void *)sum_h));

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_ScalingFactor,&gam,sizeof(float)));
	const unsigned int Grid_reg_GetConjugateGradient2 = (unsigned int)ceil((float)nodeNumber/(float)Block_reg_GetConjugateGradient2);
	dim3 B2(Block_reg_GetConjugateGradient2,1,1);
	dim3 G2(Grid_reg_GetConjugateGradient2,1,1);
	reg_GetConjugateGradient2_kernel <<< G2, B2 >>> (*nodeNMIGradientArray_d, *conjugateG_d, *conjugateH_d);
	CUDA_SAFE_CALL(cudaThreadSynchronize());
#ifndef NDEBUG
    printf("[NiftyReg CUDA DEBUG] reg_GetConjugateGradient2 kernel: %s - Grid size [%i %i %i] - Block size [%i %i %i]\n",
	       cudaGetErrorString(cudaGetLastError()),G1.x,G1.y,G1.z,B1.x,B1.y,B1.z);
#endif


}

float reg_getMaximalLength_gpu(	float4 **nodeNMIGradientArray_d,
				int nodeNumber)
{

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_NodeNumber,&nodeNumber,sizeof(int)));
	CUDA_SAFE_CALL(cudaBindTexture(0, gradientImageTexture, *nodeNMIGradientArray_d, nodeNumber*sizeof(float4)));

	// each thread extract the maximal value out of 128
	const int threadNumber = (int)ceil((float)nodeNumber/128.0f);
	const unsigned int Grid_reg_getMaximalLength = (unsigned int)ceil((float)threadNumber/(float)Block_reg_getMaximalLength);
	dim3 B1(Block_reg_getMaximalLength,1,1);
	dim3 G1(Grid_reg_getMaximalLength,1,1);

	float *all_d;
    CUDA_SAFE_CALL(cudaMalloc(&all_d, threadNumber*sizeof(float)));
	reg_getMaximalLength_kernel <<< G1, B1 >>> (all_d);
	CUDA_SAFE_CALL(cudaThreadSynchronize());
#ifndef NDEBUG
    printf("[NiftyReg CUDA DEBUG] reg_getMaximalLength kernel: %s - Grid size [%i %i %i] - Block size [%i %i %i]\n",
	       cudaGetErrorString(cudaGetLastError()),G1.x,G1.y,G1.z,B1.x,B1.y,B1.z);
#endif
    float *all_h;CUDA_SAFE_CALL(cudaMallocHost(&all_h, nodeNumber*sizeof(float)));
	CUDA_SAFE_CALL(cudaMemcpy(all_h, all_d, threadNumber*sizeof(float),cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaFree(all_d));
	double maxDistance = 0.0f;
	for(int i=0; i<threadNumber; i++) maxDistance = all_h[i]>maxDistance?all_h[i]:maxDistance;
	CUDA_SAFE_CALL(cudaFreeHost((void *)all_h));

	return (float)maxDistance;
}

void reg_updateControlPointPosition_gpu(nifti_image *controlPointImage,
					float4 **controlPointImageArray_d,
					float4 **bestControlPointPosition_d,
					float4 **nodeNMIGradientArray_d,
					float currentLength)
{
	const int nodeNumber = controlPointImage->nx * controlPointImage->ny * controlPointImage->nz;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_NodeNumber,&nodeNumber,sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_ScalingFactor,&currentLength,sizeof(float)));

	CUDA_SAFE_CALL(cudaBindTexture(0, controlPointTexture, *bestControlPointPosition_d, nodeNumber*sizeof(float4)));
	CUDA_SAFE_CALL(cudaBindTexture(0, gradientImageTexture, *nodeNMIGradientArray_d, nodeNumber*sizeof(float4)));

	const unsigned int Grid_reg_updateControlPointPosition = (unsigned int)ceil((float)nodeNumber/(float)Block_reg_updateControlPointPosition);
	dim3 B1(Block_reg_updateControlPointPosition,1,1);
	dim3 G1(Grid_reg_updateControlPointPosition,1,1);

	reg_updateControlPointPosition_kernel <<< G1, B1 >>> (*controlPointImageArray_d);
	CUDA_SAFE_CALL(cudaThreadSynchronize());
#ifndef NDEBUG
    printf("[NiftyReg CUDA DEBUG] reg_updateControlPointPosition kernel: %s - Grid size [%i %i %i] - Block size [%i %i %i]\n",
	       cudaGetErrorString(cudaGetLastError()),G1.x,G1.y,G1.z,B1.x,B1.y,B1.z);
#endif
}

void reg_gaussianSmoothing_gpu( nifti_image *image,
                                float4 **imageArray_d,
                                float sigma,
                                bool smoothXYZ[8])

{
    const int voxelNumber = image->nx * image->ny * image->nz;
    const int3 imageDim = make_int3(image->nx, image->ny, image->nz);

    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_ImageDim, &imageDim,sizeof(int3)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_VoxelNumber, &voxelNumber,sizeof(int3)));

    bool axisToSmooth[8];
    if(smoothXYZ==NULL){
        for(int i=0; i<8; i++) axisToSmooth[i]=true;
    }
    else{
        for(int i=0; i<8; i++) axisToSmooth[i]=smoothXYZ[i];
    }

    for(int n=1; n<4; n++){
        if(axisToSmooth[n]==true){
            float currentSigma;
            if(sigma>0) currentSigma=sigma/image->pixdim[n];
            else currentSigma=fabs(sigma); // voxel based if negative value
            int radius=(int)ceil(currentSigma*3.0f);
            if(radius>0){
                int kernelSize = 1+radius*2;
                float *kernel_h;
                CUDA_SAFE_CALL(cudaMallocHost(&kernel_h, kernelSize*sizeof(float)));
                float kernelSum=0;
                for(int i=-radius; i<=radius; i++){
                    kernel_h[radius+i]=(float)(exp( -(i*i)/(2.0*currentSigma*currentSigma)) / (currentSigma*2.506628274631)); // 2.506... = sqrt(2*pi)
                    kernelSum += kernel_h[radius+i];
                }
                for(int i=0; i<kernelSize; i++)
                    kernel_h[i] /= kernelSum;
                float *kernel_d;
                CUDA_SAFE_CALL(cudaMalloc(&kernel_d, kernelSize*sizeof(float)));
                CUDA_SAFE_CALL(cudaMemcpy(kernel_d, kernel_h, kernelSize*sizeof(float), cudaMemcpyHostToDevice));
                CUDA_SAFE_CALL(cudaFreeHost(kernel_h));

                float4 *smoothedImage;
                CUDA_SAFE_CALL(cudaMalloc(&smoothedImage,voxelNumber*sizeof(float4)));

                CUDA_SAFE_CALL(cudaBindTexture(0, convolutionKernelTexture, kernel_d, kernelSize*sizeof(float)));
                CUDA_SAFE_CALL(cudaBindTexture(0, gradientImageTexture, *imageArray_d, voxelNumber*sizeof(float4)));
                unsigned int Grid_reg_ApplyConvolutionWindow;
                dim3 B,G;
                switch(n){
                    case 1:
                        Grid_reg_ApplyConvolutionWindow =
                            (unsigned int)ceil((float)voxelNumber/(float)Block_reg_ApplyConvolutionWindowAlongX);
                        B=dim3(Block_reg_ApplyConvolutionWindowAlongX,1,1);
                        G=dim3(Grid_reg_ApplyConvolutionWindow,1,1);
                        _reg_ApplyConvolutionWindowAlongX_kernel <<< G, B >>> (smoothedImage, kernelSize);
                        CUDA_SAFE_CALL(cudaThreadSynchronize());
#ifndef NDEBUG
                        printf("[NiftyReg CUDA DEBUG] reg_ApplyConvolutionWindowAlongX_kernel: %s - Grid size [%i %i %i] - Block size [%i %i %i]\n",
                            cudaGetErrorString(cudaGetLastError()),G.x,G.y,G.z,B.x,B.y,B.z);
#endif
                        break;
                    case 2:
                        Grid_reg_ApplyConvolutionWindow =
                            (unsigned int)ceil((float)voxelNumber/(float)Block_reg_ApplyConvolutionWindowAlongY);
                        B=dim3(Block_reg_ApplyConvolutionWindowAlongY,1,1);
                        G=dim3(Grid_reg_ApplyConvolutionWindow,1,1);
                        _reg_ApplyConvolutionWindowAlongY_kernel <<< G, B >>> (smoothedImage, kernelSize);
                        CUDA_SAFE_CALL(cudaThreadSynchronize());
#ifndef NDEBUG
                        printf("[NiftyReg CUDA DEBUG] reg_ApplyConvolutionWindowAlongY_kernel: %s - Grid size [%i %i %i] - Block size [%i %i %i]\n",
                            cudaGetErrorString(cudaGetLastError()),G.x,G.y,G.z,B.x,B.y,B.z);
#endif
                        break;
                    case 3:
                        Grid_reg_ApplyConvolutionWindow =
                            (unsigned int)ceil((float)voxelNumber/(float)Block_reg_ApplyConvolutionWindowAlongZ);
                        B=dim3(Block_reg_ApplyConvolutionWindowAlongZ,1,1);
                        G=dim3(Grid_reg_ApplyConvolutionWindow,1,1);
                        _reg_ApplyConvolutionWindowAlongZ_kernel <<< G, B >>> (smoothedImage, kernelSize);
                        CUDA_SAFE_CALL(cudaThreadSynchronize());
#ifndef NDEBUG
                        printf("[NiftyReg CUDA DEBUG] reg_ApplyConvolutionWindowAlongZ_kernel: %s - Grid size [%i %i %i] - Block size [%i %i %i]\n",
                            cudaGetErrorString(cudaGetLastError()),G.x,G.y,G.z,B.x,B.y,B.z);
#endif
                        break;
                }
                CUDA_SAFE_CALL(cudaFree(kernel_d));
                CUDA_SAFE_CALL(cudaMemcpy(*imageArray_d, smoothedImage, voxelNumber*sizeof(float4), cudaMemcpyDeviceToDevice));
                CUDA_SAFE_CALL(cudaFree(smoothedImage));
            }
        }
    }

}


void reg_smoothImageForCubicSpline_gpu( nifti_image *image,
                                        float4 **imageArray_d,
                                        int *smoothingRadius)
{
    const int voxelNumber = image->nx * image->ny * image->nz;
    const int3 imageDim = make_int3(image->nx, image->ny, image->nz);

    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_ImageDim, &imageDim,sizeof(int3)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_VoxelNumber, &voxelNumber,sizeof(int3)));

    for(int n=0; n<3; n++){
        if(smoothingRadius[n]>0){
            int kernelSize = 1+smoothingRadius[n]*2;
            float *kernel_h;
            CUDA_SAFE_CALL(cudaMallocHost(&kernel_h, kernelSize*sizeof(float)));
//             float kernelSum=0;
            for(int i=-smoothingRadius[n]; i<=smoothingRadius[n]; i++){
                float coeff = fabs(2.0f*(float)(i)/smoothingRadius[n]);
                if(coeff<1.0f)  kernel_h[smoothingRadius[n]+i] = 2.0f/3.0f - coeff*coeff + 0.5f*coeff*coeff*coeff;
                else        kernel_h[smoothingRadius[n]+i] = -(coeff-2.0f)*(coeff-2.0f)*(coeff-2.0f)/6.0f;
//                 kernelSum += kernel_h[smoothingRadius[n]+i];
            }
//             for(int i=0; i<kernelSize; i++) kernel_h[i] /= kernelSum;
            float *kernel_d;
            CUDA_SAFE_CALL(cudaMalloc(&kernel_d, kernelSize*sizeof(float)));
            CUDA_SAFE_CALL(cudaMemcpy(kernel_d, kernel_h, kernelSize*sizeof(float), cudaMemcpyHostToDevice));
            CUDA_SAFE_CALL(cudaFreeHost(kernel_h));
            CUDA_SAFE_CALL(cudaBindTexture(0, convolutionKernelTexture, kernel_d, kernelSize*sizeof(float)));

            float4 *smoothedImage_d;
            CUDA_SAFE_CALL(cudaMalloc(&smoothedImage_d,voxelNumber*sizeof(float4)));

            CUDA_SAFE_CALL(cudaBindTexture(0, gradientImageTexture, *imageArray_d, voxelNumber*sizeof(float4)));

            unsigned int Grid_reg_ApplyConvolutionWindow;
            dim3 B,G;
            switch(n){
                case 0:
                    Grid_reg_ApplyConvolutionWindow =
                        (unsigned int)ceil((float)voxelNumber/(float)Block_reg_ApplyConvolutionWindowAlongX);
                    B=dim3(Block_reg_ApplyConvolutionWindowAlongX,1,1);
                    G=dim3(Grid_reg_ApplyConvolutionWindow,1,1);
                    _reg_ApplyConvolutionWindowAlongX_kernel <<< G, B >>> (smoothedImage_d, kernelSize);
                    CUDA_SAFE_CALL(cudaThreadSynchronize());
#ifndef NDEBUG
                    printf("[NiftyReg CUDA DEBUG] reg_ApplyConvolutionWindowAlongX_kernel: %s - Grid size [%i %i %i] - Block size [%i %i %i]\n",
                        cudaGetErrorString(cudaGetLastError()),G.x,G.y,G.z,B.x,B.y,B.z);
#endif
                    break;
                case 1:
                    Grid_reg_ApplyConvolutionWindow =
                        (unsigned int)ceil((float)voxelNumber/(float)Block_reg_ApplyConvolutionWindowAlongY);
                    B=dim3(Block_reg_ApplyConvolutionWindowAlongY,1,1);
                    G=dim3(Grid_reg_ApplyConvolutionWindow,1,1);
                    _reg_ApplyConvolutionWindowAlongY_kernel <<< G, B >>> (smoothedImage_d, kernelSize);
                    CUDA_SAFE_CALL(cudaThreadSynchronize());
#ifndef NDEBUG
                    printf("[NiftyReg CUDA DEBUG] reg_ApplyConvolutionWindowAlongY_kernel: %s - Grid size [%i %i %i] - Block size [%i %i %i]\n",
                        cudaGetErrorString(cudaGetLastError()),G.x,G.y,G.z,B.x,B.y,B.z);
#endif
                    break;
                case 2:
                    Grid_reg_ApplyConvolutionWindow =
                        (unsigned int)ceil((float)voxelNumber/(float)Block_reg_ApplyConvolutionWindowAlongZ);
                    B=dim3(Block_reg_ApplyConvolutionWindowAlongZ,1,1);
                    G=dim3(Grid_reg_ApplyConvolutionWindow,1,1);
                    _reg_ApplyConvolutionWindowAlongZ_kernel <<< G, B >>> (smoothedImage_d, kernelSize);
                    CUDA_SAFE_CALL(cudaThreadSynchronize());
#ifndef NDEBUG
                    printf("[NiftyReg CUDA DEBUG] reg_ApplyConvolutionWindowAlongZ_kernel: %s - Grid size [%i %i %i] - Block size [%i %i %i]\n",
                        cudaGetErrorString(cudaGetLastError()),G.x,G.y,G.z,B.x,B.y,B.z);
#endif
                    break;
            }
            CUDA_SAFE_CALL(cudaFree(kernel_d));
            CUDA_SAFE_CALL(cudaMemcpy(*imageArray_d, smoothedImage_d, voxelNumber*sizeof(float4), cudaMemcpyDeviceToDevice));
            CUDA_SAFE_CALL(cudaFree(smoothedImage_d));
        }
    }
}

#endif

