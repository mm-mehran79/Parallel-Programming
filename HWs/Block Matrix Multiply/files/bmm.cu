// ONLY MODIFY THIS FILE!
// YOU CAN MODIFY EVERYTHING IN THIS FILE!

#include "bmm.h"

#define tx threadIdx.x
#define ty threadIdx.y
#define tz threadIdx.z

#define bx blockIdx.x
#define by blockIdx.y
#define bz blockIdx.z

// TILEX and TILEY are used to set the number of threads in a CUDA block
#define TILEX 32
#define TILEY 16

// DEPH is used to set size of array in shared memmory
#define DEPTH 128

// you may define other parameters here!
// you may define other macros here!
//#define MOD(x, y) ((x) & ((1 << (y)) - 1))
// #define MOD(x, y) ((x) % ((1 << (y))))
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
	int k, l;
	int tyk, txk;
	const int i = ty + (by * TILEY);
	const int j = tx + (bx * TILEX);
	for (k = 0; k < n; k += DEPTH)
	{
		txk = k + tx;
		for (l = 0; l < DEPTH; l += TILEX)
			As[ty][tx + l] = mem2d(ad, m, i, txk + l);
		tyk = k + ty;
		for (l = 0; l < DEPTH; l += TILEY)
			Bs[ty + l][tx] = mem2d(bd, m, tyk + l, j);
		__syncthreads();
		
		for (l = 0; l < DEPTH; ++l)
		{
			tempC += As[ty][l] * Bs[l][tx];
		}
		__syncthreads();
	}
	mem2d(cd, m, i, j) = tempC;
}
