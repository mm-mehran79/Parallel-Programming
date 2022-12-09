// ONLY MODIFY THIS FILE!
// YOU CAN MODIFY EVERYTHING IN THIS FILE!

#include "bmm.h"

#define tx threadIdx.x
#define ty threadIdx.y
#define tz threadIdx.z

#define bx blockIdx.x
#define by blockIdx.y
#define bz blockIdx.z
#define _UI16_MAX 0xffff

#define TILEXSHIFT 4
#define TILEYSHIFT 4
#define DEPTHSHIFT 3
// TILEX and TILEY are used to set the number of threads in a CUDA block
#define TILEX (1 << TILEXSHIFT)
#define TILEY (1 << TILEYSHIFT)

// DEPH is used to set size of array in shared memmory
#define DEPTH (1 << DEPTHSHIFT)

// you may define other parameters here!
// you may define other macros here!
// you may define other functions here!

dim3 getDimGrid(const int m, const int n)
{
	dim3 dimGrid(n / TILEX, n / TILEY);
	return dimGrid;
}
dim3 getDimBlock(const int m, const int n)
{
	dim3 dimBlock(TILEX, TILEY);
	return dimBlock;
}
__global__ void kernelFunc(float *ad, float *bd, float *cd, const int m, const int n)
{
	__shared__ float As[TILEY][DEPTH];
	__shared__ float Bs[DEPTH][TILEX];
	float tempC = 0;
	unsigned short k, l;
	unsigned short i = ty + (by << TILEYSHIFT);
	unsigned short j = tx + (bx << TILEXSHIFT);
	for (k = 0; k < (n >> DEPTHSHIFT); k++)
	{
		__syncthreads();
		if ( (tx >> DEPTHSHIFT) == (k & (~(_UI16_MAX<<(TILEXSHIFT - DEPTHSHIFT)))) )			// tx / DEPTH == (k % TILEX) / DEPTH
			As[ty][tx & (~(_UI16_MAX << DEPTHSHIFT)) ] = mem2d(ad, m, i, k << DEPTHSHIFT + (ty & (~(_UI16_MAX << DEPTHSHIFT))) );
		if ( (ty >> DEPTHSHIFT) == (k & (~(_UI16_MAX<<(TILEYSHIFT - DEPTHSHIFT)))) )
			Bs[ ty & (~(_UI16_MAX << DEPTHSHIFT)) ][tx] = mem2d(bd, m, k << DEPTHSHIFT + (ty & (~(_UI16_MAX << DEPTHSHIFT))), j);
		__syncthreads();

		for (l = 0; l < DEPTH; l++)
		{
			tempC += As[ty][l] * Bs[l][tx];
		}
	}
	mem2d(cd, m, i, j) = tempC;

	// write your GPU kernel function here
}
