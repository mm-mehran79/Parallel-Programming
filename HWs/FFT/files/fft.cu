//ONLY MODIFY THIS FILE!
//YOU CAN MODIFY EVERYTHING IN THIS FILE!

#include "fft.h"

#define tx threadIdx.x
#define ty threadIdx.y
#define tz threadIdx.z

#define bx blockIdx.x
#define by blockIdx.y
#define bz blockIdx.z

#define Radix (4)
#define gridDim (1<<10)
#define blockDim (1<<6)
#define BETA (blockDim*1)
#define ALPHA (gridDim*1)
#define uint unsigned int

#define THREAD_PER_BLOCK_SORT 1024 
#define BLOCK_X  1024
#define THREAD_PER_BLOCK_FLY  256  

// you may define other parameters here!
// you may define other macros here!
// you may define other functions here!
// __global__ void sharedMemoryFft(float* x_r_d, float* x_i_d, const unsigned int N, const unsigned int M)
// {
// 	__shared__ float x_shared_r[BETA];
// 	__shared__ float x_shared_i[BETA];//could more than one element
	

	

// }

__global__ void gpuSort(float* x, float* x_temp, const unsigned int M)
{
	uint index = tx + ( by * BLOCK_X + bx  ) * THREAD_PER_BLOCK_SORT;
	uint reversed2bit = index;
	reversed2bit = (reversed2bit & 0xFFFF0000) >> 16 | (reversed2bit & 0x0000FFFF) << 16;
	reversed2bit = (reversed2bit & 0xFF00FF00) >>  8 | (reversed2bit & 0x00FF00FF) <<  8;
	reversed2bit = (reversed2bit & 0xF0F0F0F0) >>  4 | (reversed2bit & 0x0F0F0F0F) <<  4;
	reversed2bit = (reversed2bit & 0xAAAAAAAA) >>  2 | (reversed2bit & 0x33333333) <<  2;
	reversed2bit = reversed2bit >> (32-M);
	x_temp[index] = x[reversed2bit];
}

__global__ void gpuCopy(float* x, float* x_temp)
{
	uint index = tx + ( by * BLOCK_X + bx  ) * THREAD_PER_BLOCK_SORT;
	x[index] = x_temp[index];
}

__global__ void gpuButterfly(float* x_r_d, float* x_i_d, uint v)
{
	uint index = tx + ( by * BLOCK_X + bx  ) * THREAD_PER_BLOCK_FLY / 2;
	uint p = index % (1<<v);
	uint fl = 1<<v;
	uint indexNew = (index/fl)*(4*fl) + p;
	float const theta = (PI * p)/(2*fl);
	float cosinus, sinus;
	float x_i_temp, x_r_temp;


	float w0_real, w0_imag, w1_real, w1_imag, w2_real, w2_imag, w3_real, w3_imag;
	w0_real = x_r_d[indexNew];
	w0_imag = x_i_d[indexNew];
	
	indexNew += fl;
	x_r_temp = x_r_d[indexNew];
	x_i_temp = x_i_d[indexNew];
	sincosf(theta, &sinus, &cosinus);
	w1_real = cosinus * x_r_temp + sinus * x_i_temp;
	w1_imag = cosinus * x_i_temp - sinus * x_r_temp;

	indexNew += fl;
	float theta_temp = theta + theta;
	x_r_temp = x_r_d[indexNew];
	x_i_temp = x_i_d[indexNew];
	sincosf(theta_temp, &sinus, &cosinus);
	w2_real = cosinus * x_r_temp + sinus * x_i_temp;
	w2_imag = cosinus * x_i_temp - sinus * x_r_temp;

	indexNew += fl;
	theta_temp += theta;
	x_r_temp = x_r_d[indexNew];
	x_i_temp = x_i_d[indexNew];
	sincosf(theta_temp, &sinus, &cosinus);
	w3_real = cosinus * x_r_temp + sinus * x_i_temp;
	w3_imag = cosinus * x_i_temp - sinus * x_r_temp;

	x_r_d[index] = w0_real - w1_imag - w2_real + w3_imag;
	x_i_d[index] = w0_imag + w1_real - w2_imag - w3_real;
	index += fl;
	x_r_d[index] = w0_real - w1_real + w2_real - w3_real;
	x_i_d[index] = w0_imag - w1_imag + w2_imag - w3_imag;
	index += fl;
	x_r_d[index] = w0_real + w1_imag - w2_real - w3_imag;
	x_i_d[index] = w0_imag - w1_real - w2_imag + w3_real;
	index += fl;
	x_r_d[index] = w0_real + w1_real + w2_real + w3_real;
	x_i_d[index] = w0_imag + w1_imag + w2_imag + w3_imag;
}

//-----------------------------------------------------------------------------
void gpuKernel(float* x_r_d, float* x_i_d, const unsigned int N, const unsigned int M)
{
	float *d_temp;
	HANDLE_ERROR(cudaMalloc((void**) & d_temp, N * sizeof(float)));
	dim3 dimGrid_Swap(BLOCK_X,N/BLOCK_X/THREAD_PER_BLOCK_SORT); 
	dim3 dimBlock_Swap(THREAD_PER_BLOCK_SORT,1);
	dim3 dimGrid_Butterfly(BLOCK_X,N/4/BLOCK_X/(THREAD_PER_BLOCK_FLY/2)); 
	dim3 dimBlock_Butterfly(THREAD_PER_BLOCK_FLY/2,1);
	gpuSort<<< dimGrid_Swap,dimBlock_Swap >>>(x_r_d, d_temp, M);
	gpuCopy<<< dimGrid_Swap,dimBlock_Swap >>>(x_r_d, d_temp);
	gpuSort<<< dimGrid_Swap,dimBlock_Swap >>>(x_i_d, d_temp, M);
	gpuCopy<<< dimGrid_Swap,dimBlock_Swap >>>(x_i_d, d_temp);
	HANDLE_ERROR(cudaFree(d_temp));
	for(uint i = 0; i < M; i += 2)
		gpuButterfly<<< dimGrid_Butterfly ,dimBlock_Butterfly >>>(x_r_d,x_i_d,i);


	// In this function, both inputs and outputs are on GPU.
	// No need for cudaMalloc, cudaMemcpy or cudaFree.
	// This function does not run on GPU. 
	// You need to define another function and call it here for GPU execution.
	
}
