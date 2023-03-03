// Include your C header files here
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include "pth_msort.h"

#define uint unsigned int
#define FIRST_STAGE_THREAD 8U
#define LAST_NODE_THREAD 4U
#define RADIX_OF_SORTING 4U
#define LOG2_RADIX_OF_SORTING 2U

struct firstStageSortArgs
{
	int* values;
	int* sorted;
	uint start;
	uint end;
};
struct middleStageSortArgs
{
	int* values;
	int* sorted;
	uint start;
	uint middle;
	uint end;
};

struct lastStageSortArgs
{
	int* values;
	int* sorted;
	uint i_start, j_start;
	uint i_end, j_end;
	uint start;
};

void swap(int* a, int* b)
{
	int t = *a;
	*a = *b;
	*b = t;
}

uint partition(int arr[], uint l, uint h)
{
	int x = arr[h];
	uint i = (l - 1), j;
  
	for (j = l; j <= h - 1; j++) {
		if (arr[j] <= x) {
			i++;
			swap(&arr[i], &arr[j]);
		}
	}
	swap(&arr[i + 1], &arr[h]);
	return (i + 1);
}


void *firstStageSort(void *arg)
{
	struct firstStageSortArgs *inArgs = (struct firstStageSortArgs *) arg;
	uint start = inArgs->start;
	uint end = inArgs->end;
	int* values = inArgs->values;
	int* sorted = inArgs->sorted;
	
	/* // RADIX Sort::
	uint i, j, count[RADIX_OF_SORTING];
	int *temp;
	for (i = 0U; i < 32U; i+=LOG2_RADIX_OF_SORTING)
	{
		// uint count[RADIX_OF_SORTING] = {0};
		for (j = 0; j < RADIX_OF_SORTING; j++)
			count[j] = 0;
		for (j = start; j < end; j++)
			count[(values[j]>>i)%RADIX_OF_SORTING]++;
		for (j = 1; j < RADIX_OF_SORTING; j++)
			count[j] += count[j - 1];
		for (j = end - 1; j >= start; j--)
			sorted[count[(values[j]>>i)%RADIX_OF_SORTING]--] = values[j]; // Error : this line cause core dump
		
		temp = values;
		values = sorted;
		sorted = temp;
	} */

	// Quick sort:
	// "sorted" used as stack                    
	uint low = start, high = end - 1U;              // high and low are index for values; and top and start are used for sorted;
	uint *stack = (uint*)&sorted[start];			// initialize top of stack
	uint *stackStart = stack;
  
	// push initial values of l and h to stack
	*(stack++) = low;
	*(stack++) = high;
  
	// Keep popping from stack while is not empty
	while (stack != stackStart) {
		// Pop h and l
		high = *(--stack);
		low = *(--stack);
  
		// Set pivot element at its correct position
		// in stack array
		uint p = partition(values, low, high);
  
		// If there are elements on left side of pivot,
		// then push left side to stack
		if (p > low + 1U) {
			*(stack++) = low;
			*(stack++) = p - 1U;
		}
		// If there are elements on right side of pivot,
		// then push right side to stack
		if (p + 1U < high) {
			*(stack++) = p + 1U;
			*(stack++) = high;
		}
	}
	// swap two array: (it could be implemented better if i use pointer to array instead of array in pth_msort)
	uint i;
	for (i = start; i < end; i++) sorted[i] = values[i];
}

void *middleStageSort(void *arg)
{
	struct middleStageSortArgs *inArgs = (struct middleStageSortArgs *) arg;
	uint start = inArgs->start;
	uint middle = inArgs->middle;
	uint end = inArgs->end;
	int *values = inArgs->values;
	int *sorted = inArgs->sorted;
	// merge two sorted array:
	uint i = start, j = middle, k = start;    
	while ((i < middle) && (j < end))
	{
		if(values[i] < values[j])
			sorted[k++] = values[i++];                   
		else
			sorted[k++] = values[j++];                   
	}
	while (i < middle) sorted[k++] = values[i++];
	while (j < end) sorted[k++] = values[j++];
}

void *lastStageSort(void *arg)
{
	struct lastStageSortArgs *inArgs = (struct lastStageSortArgs *) arg;
	uint start = inArgs->start;
	uint i_start = inArgs->i_start;
	uint i_end = inArgs->i_end;
	uint j_start = inArgs->j_start;
	uint j_end = inArgs->j_end;
	int *values = inArgs->values;
	int *sorted = inArgs->sorted;
	// merge two sorted array:
	while (i_start < i_end && j_start < j_end)
		if(values[i_start] < values[j_start])
			sorted[start++] = values[i_start++];
		else
			sorted[start++] = values[j_start++];
	while (i_start < i_end) sorted[start++] = values[i_start++];
	while (j_start < j_end) sorted[start++] = values[j_start++];
}

uint binarySearch(const int *arr, uint lowIndex, uint highIndex, int searchQuery)
{
	while (lowIndex <= highIndex)
	{
		uint midIndex = lowIndex + (highIndex - lowIndex) / 2;
		if (arr[midIndex] == searchQuery)
			return midIndex;
		if (arr[midIndex] < searchQuery)
			lowIndex = midIndex + 1;
		else
			highIndex = midIndex - 1;
	}
	return lowIndex;        //arr[lowIndex] > searchQuery
}


void mergeSortParallel (const int* values, unsigned int N, int* sorted)
{
	// first stage sort (values => sorted)
	uint i, temp = 0;
	pthread_t handle[FIRST_STAGE_THREAD];
	struct firstStageSortArgs args[FIRST_STAGE_THREAD];
	
	for (i = 0U; i < FIRST_STAGE_THREAD; i++)
	{
		args[i].start = temp;
		temp += N/FIRST_STAGE_THREAD;
		args[i].end = temp;
		args[i].values = (int *) values;
		args[i].sorted = sorted;
		pthread_create(&handle[i], NULL, firstStageSort, (void *) &args[i]);
	}

	struct middleStageSortArgs middleArgs[(FIRST_STAGE_THREAD/2U)];
	temp = 0;
	// middle stage serial merge sort (sorted => values)
	for (i = 0U; i < (FIRST_STAGE_THREAD/2U); i++)
	{
		middleArgs[i].start = temp;
		temp += N/FIRST_STAGE_THREAD;
		middleArgs[i].middle = temp;
		temp += N/FIRST_STAGE_THREAD;
		middleArgs[i].end = temp;
		middleArgs[i].values = sorted;
		middleArgs[i].sorted = (int *)values;
	}
	
	for ( i = 0U; i < FIRST_STAGE_THREAD; i++)
	{
		pthread_join(handle[i], NULL);
	}
	for (i = 0U; i < (FIRST_STAGE_THREAD/2U); i++)
	{
		pthread_create(&handle[i], NULL, middleStageSort, (void *) &middleArgs[i]);
	}
	
	struct lastStageSortArgs lastArgs[LAST_NODE_THREAD];
	temp = 0U; // N/(LAST_NODE_THREAD * 2U);
	uint lowIndex = N / 2U;
	for (i = 0U; i < (FIRST_STAGE_THREAD/2U); i++) pthread_join(handle[i], NULL);

	for (i = 0U; i < LAST_NODE_THREAD; i++)
	{
		lastArgs[i].values = (int *)values;
		lastArgs[i].sorted = sorted;
		lastArgs[i].start = temp + (lowIndex - (N/2U));
		lastArgs[i].i_start = temp;
		temp += N/(LAST_NODE_THREAD * 2U);
		lastArgs[i].i_end = temp;
		lastArgs[i].j_start = lowIndex;
		if(i == (LAST_NODE_THREAD-1)) lastArgs[i].j_end = N;
		else
		{
			lowIndex = binarySearch(values, lowIndex, (N-1U), values[temp]);
			lastArgs[i].j_end = lowIndex;
		}
		pthread_create(&handle[i], NULL, lastStageSort, (void*) &lastArgs[i]);        
	}
	for (i = 0u; i < LAST_NODE_THREAD; i++) pthread_join(handle[i], NULL);
}
