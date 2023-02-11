// Include your C header files here
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include "pth_msort.h"

#define uint unsigned int
#define FIRST_STAGE_THREAD 4U
#define LAST_NODE_THREAD 4U
#define RADIX_OF_SORTING 4U

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


void *firstStageSort(void *arg)
{
    struct firstStageSortArgs *inArgs = (struct firstStageSortArgs *) arg;
    uint start = inArgs->start;
    uint end = inArgs->end;
    int* values = inArgs->values;
    int* sorted = inArgs->sorted;
    //const uint LENGTH = end - start;

    uint count[RADIX_OF_SORTING];
    uint i, j;
    int *temp;
    for (i = 0U; i < 32U; i+=RADIX_OF_SORTING)
    {
        uint count[RADIX_OF_SORTING] = {0};
        for (j = start; j < end; j++)
            count[(values[j]>>i)&RADIX_OF_SORTING]++;
        for (j = 1; j < RADIX_OF_SORTING; j++)
            count[j] += count[j - 1];
        for (j = end - 1; j >= 0; j--)
        {
            sorted[count[(values[j]>>i)&RADIX_OF_SORTING] - 1] = values[j];
            count[(values[j]>>i)&RADIX_OF_SORTING]--;
        }
        temp = values;
        values = sorted;
        sorted = temp;
    }
    

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
    while (i < middle && j < end)
        if(values[i] < values[j])
            sorted[k++] = values[i++];
        else
            sorted[k++] = values[j++];
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
    printf("0");
    uint i, temp = 0;
    pthread_t handle[FIRST_STAGE_THREAD];
    struct firstStageSortArgs args[FIRST_STAGE_THREAD];
    
    printf("1");
    for (i = 0U; i < FIRST_STAGE_THREAD; i++)
    {
        args[i].start = temp;
        temp += N/FIRST_STAGE_THREAD;
        args[i].end = temp;
        args[i].values = (int *) values;
        args[i].sorted = sorted;
        pthread_create(&handle[i], NULL, firstStageSort, (void *) &args[i]);
    }
    printf("2");

    struct middleStageSortArgs middleArgs[(FIRST_STAGE_THREAD/2)];
    temp = 0;
    // middle stage serial merge sort (sorted => values)
    for (i = 0U; i < (FIRST_STAGE_THREAD/2U); i++)
    {
        middleArgs[i].start = i * N/FIRST_STAGE_THREAD;
        temp += N/FIRST_STAGE_THREAD;
        middleArgs[i].middle = temp;
        temp += N/FIRST_STAGE_THREAD;
        middleArgs[i].end = temp;
        middleArgs[i].values = sorted;
        middleArgs[i].sorted = (int *)values;
    }
    
    printf("3");
    pthread_t middleHandle[(FIRST_STAGE_THREAD/2U)];
    for (i = 0U; i < FIRST_STAGE_THREAD; i++)
    {
        pthread_join(handle[i++], NULL);
        pthread_join(handle[i], NULL);
        pthread_create(&middleHandle[(i/2U)], NULL, middleStageSort, (void*) &middleArgs[i]);
    }
    
    printf("4");

    pthread_t lastNodeHandle[LAST_NODE_THREAD];
    struct lastStageSortArgs lastArgs[LAST_NODE_THREAD];
    temp = 0U; // N/(LAST_NODE_THREAD * 2U);
    uint lowIndex = N / 2U;
    for (i = 0U; i < (FIRST_STAGE_THREAD/2U); i++) 
        pthread_join(handle[i], NULL);

    
    printf("5");
    // last node sort (values)
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
        pthread_create(&lastNodeHandle[i], NULL, lastStageSort, (void*) &lastArgs[i]);        
    }
    
    printf("6");
    for (i = 0u; i < LAST_NODE_THREAD; i++)
        pthread_join(lastNodeHandle[i], NULL);
    
    printf("7");
}
