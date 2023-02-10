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

//Replace all sincosf with __sincosf for faster speed but lower accuracy (and vice versa)

// you may define other parameters here!
// you may define other macros here!
// you may define other functions here!
// __global__ void sharedMemoryFft(float* x_r_d, float* x_i_d, const unsigned int N, const unsigned int M)
// {
// 	__shared__ float x_shared_r[BETA];
// 	__shared__ float x_shared_i[BETA];//could more than one element
	

	

// }

__global__ void gpuSort(float* x, float* x_temp, const unsigned int N, const unsigned int M)
{
	// int id = tx + ( by * BLOCK_X + bx  ) * THREAD_PER_BLOCK_SORT;
	// int new_id = id;
	// int half = N/4;
	// int b = 0;

	// for(int i = 0 ; i<M; i++ ){
	// 	if( new_id % 4 == 0 ){
	// 		new_id = (new_id-b)/4 + b;
    // 	}
    // 	else if ( new_id % 4 == 1 ) {
	// 		new_id = (new_id-b)/4 + half + b;
	// 		b += half ;
    // 	}
	// 	else if ( new_id % 4 == 2 ) {
	// 		new_id = (new_id-b)/4 + 2*half + b;
	// 		b += 2*half ;
    // 	}
	// 	else{
	// 		new_id = (new_id-b)/4 + 3*half + b;
	// 		b += 3*half ;
    // 	}
	// 	half /= 4;
	// }
	// x_temp[new_id]=x[id];
	///////////////////// 2nd way:
	// int index = tx + ( by * BLOCK_X + bx  ) * THREAD_PER_BLOCK_SORT;
	// int reversed2bit = 0, index_copy = index;
	// for(int i = 0; i < M; i++)
	// {
	// 	reversed2bit *= 4;
	// 	reversed2bit += index_copy%4;
	// 	index_copy /= 4;
	// }
	// uint reversed2bit = index;
	// uint reversed2bit = (index & 0xFFFF0000U) >> 16 | (index & 0x0000FFFFU) << 16;
	// uint reversed2bit1 = (reversed2bit & 0xFF00FF00U) >>  8 | (reversed2bit & 0x00FF00FFU) <<  8;
	// uint reversed2bit2 = (reversed2bit1 & 0xF0F0F0F0U) >>  4 | (reversed2bit1 & 0x0F0F0F0FU) <<  4;
	// uint reversed2bit3 = (reversed2bit2 & 0xCCCCCCCCU) >>  2 | (reversed2bit2 & 0x33333333U) <<  2;
	// uint reversed2bit4 = reversed2bit3 >> (32-M);
	// x_temp[index] = x[reversed2bit4];
	//////////////////////////////// 4th way:
	uint index = tx + ( by * BLOCK_X + bx  ) * THREAD_PER_BLOCK_SORT;
	uint indexBitReversed = __brev(index);
	uint index2BitReversed = (indexBitReversed & 0xAAAAAAAAU) >>  1 | (indexBitReversed & 0x55555555U) <<  1;
	uint index2BitReversedOut = index2BitReversed >> (32U-M);
	x_temp[index2BitReversedOut] = x[index];
}

__global__ void gpuSortRADIX2(float* x, float* x_temp, const unsigned int N, const unsigned int M)
{
	uint index = tx + ( by * BLOCK_X + bx  ) * THREAD_PER_BLOCK_SORT;
	uint indexBitReversed = __brev(index);
	uint indexBitReversedOut = indexBitReversed >> (32U-M);
	x_temp[indexBitReversedOut] = x[index];
}

__global__ void gpuCopy(float* x, float* x_temp)
{
	uint index = tx + ( by * BLOCK_X + bx  ) * THREAD_PER_BLOCK_SORT;
	x[index] = x_temp[index];
}

__global__ void gpuButterfly(float* x_r_d, float* x_i_d, uint v)
{
	uint index = tx + ( by * BLOCK_X + bx  ) * THREAD_PER_BLOCK_FLY / 2;
	uint fl = 1<<v;
	uint p = index % fl;
	uint indexNew = (index/fl)*(4*fl) + p;
	index = indexNew;
	float const theta = (PI * p)/(2*fl);
	float cosinus, sinus;
	float x_i_temp, x_r_temp;


	float w0_real, w0_imag, w1_real, w1_imag, w2_real, w2_imag, w3_real, w3_imag;
	w0_real = x_r_d[indexNew];
	w0_imag = x_i_d[indexNew];
	
	indexNew += fl;
	x_r_temp = x_r_d[indexNew];
	x_i_temp = x_i_d[indexNew];
	__sincosf(theta, &sinus, &cosinus);
	w1_real = cosinus * x_r_temp + sinus * x_i_temp;
	w1_imag = cosinus * x_i_temp - sinus * x_r_temp;

	indexNew += fl;
	float theta_temp = theta + theta;
	x_r_temp = x_r_d[indexNew];
	x_i_temp = x_i_d[indexNew];
	__sincosf(theta_temp, &sinus, &cosinus);
	w2_real = cosinus * x_r_temp + sinus * x_i_temp;
	w2_imag = cosinus * x_i_temp - sinus * x_r_temp;

	indexNew += fl;
	theta_temp += theta;
	x_r_temp = x_r_d[indexNew];
	x_i_temp = x_i_d[indexNew];
	__sincosf(theta_temp, &sinus, &cosinus);
	w3_real = cosinus * x_r_temp + sinus * x_i_temp;
	w3_imag = cosinus * x_i_temp - sinus * x_r_temp;

	x_r_d[index] = w0_real + w1_real + w2_real + w3_real;
	x_i_d[index] = w0_imag + w1_imag + w2_imag + w3_imag;
	index += fl;
	x_r_d[index] = w0_real + w1_imag - w2_real - w3_imag;
	x_i_d[index] = w0_imag - w1_real - w2_imag + w3_real;
	index += fl;
	x_r_d[index] = w0_real - w1_real + w2_real - w3_real;
	x_i_d[index] = w0_imag - w1_imag + w2_imag - w3_imag;
	index += fl;
	x_r_d[index] = w0_real - w1_imag - w2_real + w3_imag;
	x_i_d[index] = w0_imag + w1_real - w2_imag - w3_real;
}

__global__ void gpuButterflyRADIX2(float* x_r_d, float* x_i_d, uint v)
{
	uint index = tx + ( by * BLOCK_X + bx  ) * THREAD_PER_BLOCK_FLY;
	uint fl = 1<<v;
	uint p = index % fl;
	uint indexNew = (index/fl)*(2*fl) + p;
	index = indexNew;
	float const theta = (PI * p)/(fl);
	float cosinus, sinus;
	float x_i_temp, x_r_temp;


	float w0_real, w0_imag, w1_real, w1_imag;
	w0_real = x_r_d[indexNew];
	w0_imag = x_i_d[indexNew];
	
	indexNew += fl;
	x_r_temp = x_r_d[indexNew];
	x_i_temp = x_i_d[indexNew];
	__sincosf(theta, &sinus, &cosinus);
	w1_real = cosinus * x_r_temp + sinus * x_i_temp;
	w1_imag = cosinus * x_i_temp - sinus * x_r_temp;

	x_r_d[index] = w0_real + w1_real;
	x_i_d[index] = w0_imag + w1_imag;
	index += fl;
	x_r_d[index] = w0_real - w1_real;
	x_i_d[index] = w0_imag - w1_imag;
}

//-----------------------------------------------------------------------------
void gpuKernel(float* x_r_d, float* x_i_d, const unsigned int N, const unsigned int M)
{
	float *d_temp;
	HANDLE_ERROR(cudaMalloc((void**) & d_temp, N * sizeof(float)));
	if (M%2)
	{
		//radix 2
		dim3 dimGrid_Swap(BLOCK_X,N/BLOCK_X/THREAD_PER_BLOCK_SORT); 
		dim3 dimBlock_Swap(THREAD_PER_BLOCK_SORT,1);
		gpuSortRADIX2<<< dimGrid_Swap,dimBlock_Swap >>>(x_r_d, d_temp, N, M);
		gpuCopy<<< dimGrid_Swap,dimBlock_Swap >>>(x_r_d, d_temp);

		dim3 dimGrid_Butterfly(BLOCK_X,N/2/BLOCK_X/(THREAD_PER_BLOCK_FLY)); 
		dim3 dimBlock_Butterfly(THREAD_PER_BLOCK_FLY,1);
		gpuSortRADIX2<<< dimGrid_Swap,dimBlock_Swap >>>(x_i_d, d_temp, N, M);
		gpuCopy<<< dimGrid_Swap,dimBlock_Swap >>>(x_i_d, d_temp);

		for(uint i = 0; i < M; i += 1)
			gpuButterflyRADIX2<<< dimGrid_Butterfly ,dimBlock_Butterfly >>>(x_r_d,x_i_d,i);
		HANDLE_ERROR(cudaFree(d_temp));
	}
	else
	{
		//radix 4
		dim3 dimGrid_Swap(BLOCK_X,N/BLOCK_X/THREAD_PER_BLOCK_SORT); 
		dim3 dimBlock_Swap(THREAD_PER_BLOCK_SORT,1);
		gpuSort<<< dimGrid_Swap,dimBlock_Swap >>>(x_r_d, d_temp, N, M);
		gpuCopy<<< dimGrid_Swap,dimBlock_Swap >>>(x_r_d, d_temp);
		
		dim3 dimGrid_Butterfly(BLOCK_X,N/4/BLOCK_X/(THREAD_PER_BLOCK_FLY/2)); 
		dim3 dimBlock_Butterfly(THREAD_PER_BLOCK_FLY/2,1);
		gpuSort<<< dimGrid_Swap,dimBlock_Swap >>>(x_i_d, d_temp, N, M);
		gpuCopy<<< dimGrid_Swap,dimBlock_Swap >>>(x_i_d, d_temp);
		
		
		for(uint i = 0; i < M; i += 2)
			gpuButterfly<<< dimGrid_Butterfly ,dimBlock_Butterfly >>>(x_r_d,x_i_d,i);
		HANDLE_ERROR(cudaFree(d_temp));
	}
	

		
	// In this function, both inputs and outputs are on GPU.
	// No need for cudaMalloc, cudaMemcpy or cudaFree.
	// This function does not run on GPU. 
	// You need to define another function and call it here for GPU execution.
	
}
