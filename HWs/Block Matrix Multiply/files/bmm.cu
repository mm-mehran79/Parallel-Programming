// ONLY MODIFY THIS FILE!
// YOU CAN MODIFY EVERYTHING IN THIS FILE!

#include "bmm.h"

#define tx threadIdx.x
#define ty threadIdx.y
#define tz threadIdx.z

#define bx blockIdx.x
#define by blockIdx.y
#define bz blockIdx.z

#define TILEXSHIFT 5
#define TILEYSHIFT 5
#define DEPTHSHIFT 7
// TILEX and TILEY are used to set the number of threads in a CUDA block
#define TILEX (1 << TILEXSHIFT)
#define TILEY (1 << TILEYSHIFT)

// DEPH is used to set size of array in shared memmory
#define DEPTH (1 << DEPTHSHIFT)

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
	unsigned short k, l;
	const unsigned short i = ty + (by << TILEYSHIFT);
	const unsigned short j = tx + (bx << TILEXSHIFT);
	for (k = 0; k < (n >> DEPTHSHIFT); k++)
	{
		__syncthreads();
		#if DEPTH == TILEX
			As[ty][tx] = mem2d(ad, m, i, (k<<DEPTHSHIFT) + tx);
		#elif DEPTH > TILEX
			for (l = 0; l < (1 << (DEPTHSHIFT-TILEXSHIFT)); l++)
				As[ty][tx + (l<<TILEXSHIFT)] = mem2d(ad, m, i, (k << DEPTHSHIFT) + tx + (l<<TILEXSHIFT));
		#else
			if ( (tx >> DEPTHSHIFT) == 0 )
				As[ty][tx] = mem2d(ad, m, i, (k << DEPTHSHIFT) + tx);
		#endif

		#if DEPTH == TILEY
			Bs[ty][tx] = mem2d(bd, m, (k<<DEPTHSHIFT) + ty, j);
		#elif DEPTH > TILEY
			for (l = 0; l < (1 << (DEPTHSHIFT-TILEYSHIFT)); l++)
				Bs[ty + (l<<TILEYSHIFT)][tx] = mem2d(bd, m, (k << DEPTHSHIFT) + ty + (l<<TILEYSHIFT), j);
		#else
			if ( (ty >> DEPTHSHIFT) == 0 )
				Bs[ty][tx] = mem2d(bd, m, (k << DEPTHSHIFT) + ty, j);
		#endif
		__syncthreads();
		for (l = 0; l < DEPTH; l++)
		{
			tempC += As[ty][l] * Bs[l][tx];
		}
	}
	mem2d(cd, m, i, j) = tempC;
}
