/*
 * Kyle Waas
 * Programming IV
 *
 * Hey Professor I just wanted to say that it worked for me so let me know if you have any questions. Also wanted to say that
 * I had a great semester!
 *
 * Sources:
 * https://www.youtube.com/watch?v=tjb2QQ9oPCE
 * https://www.youtube.com/watch?v=mNo9QlT_rbw
 * https://www.programmingsimplified.com/c/program/c-program-check-if-array-is-sorted-or-not
 *
 * Commands to compile and run:
 * mpicc -o mergesort mergesort.c
 * mpirun -np 10
 * ./mergesort
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#define true 1

int numbertasks;
int MERGESORT_CUTOFF = 3;

	//Taken from Class and changed to fit Comparable
	void swap(int x[], int a, int b) {
		int temp = x[a];
		x[a] = x[b];
		x[b] = temp;
	}


	//Taken from class an changed to fit Comparable
	int medianIndex(int x[], int a, int b, int c) {
		return (x[a] < x[b]) ? (x[b] < x[c] ? b : x[a] < x[c] ? c: a) :(x[b] > x[c] ? b : x[a] > x[c] ? c: a);
	}

	//Taken from class an changed to fit Comparable
	void realQuickSort(int x[], int start, int length) {

		//long startTime = System.currentTimeMillis();

		int medIndex;

		int lowIndex = start;
		int highIndex = start + length - 1;
		int midIndex = (lowIndex + highIndex) / 2;

		medIndex = medianIndex (x, lowIndex, highIndex, midIndex);

		int pivot = x[medIndex];

		int left = lowIndex;
		int right = start + length - 1;

		while (true) {
			while (left <= right && x[left] < pivot) {
				//NumberOfCompares5 ++;
				left++;
			}
			while (right >= left && x[right] > pivot) {
				//NumberOfCompares5++;
				right--;
			}
			if (left > right) {
				break;
			}
			swap (x, left, right);
			//NumberOfSwaps5 ++;
			left++;
			right--;
		}
		if (right + 1 - start > 1) {
			realQuickSort(x, start, right + 1 - start);
		}
		if (start + length - left > 1) {
			realQuickSort(x, left, start + length - left);
		}

		//long endTime = System.currentTimeMillis();
		//TimeToComplete5 = (long) ((endTime - startTime) * 0.001);
	}

	//Taken from class an changed to fit Comparable
	void quickSort(int x[]) {
		realQuickSort(x,0, 100000/numbertasks);
	}

		//Taken from class an changed to fit Comparable
	void realMergeSort(int source[], int dest[], int low, int high) {

		//long startTime = System.currentTimeMillis();

		int length = high - low;

		if (length < MERGESORT_CUTOFF) {
			for (int i = low; i < high; i++) {
				for (int j = i; j > low && dest[j-1] > dest[j]; j--) {
					swap(dest, j, j-1);
					//NumberOfSwaps4++;
				}
			}
			return;
		}

		int mid = (low + high) >> 1;

		realMergeSort(dest, source, low, mid);
		realMergeSort(dest, source, mid, high);
		//NumberOfCompares4++;

		int i = low;
		int leftPointer = low;
		int rightPointer = mid;

		for( ; i < high; i++) {
			//NumberOfCompares4++;
			if (rightPointer >= high || (leftPointer< mid && source[leftPointer] < source[rightPointer])) {
				dest[i] = source[leftPointer];
				leftPointer++;
			}
			else {
				dest[i] = source[rightPointer];
				rightPointer++;
			}
		}

		//long endTime = System.currentTimeMillis();
		//TimeToComplete4 = (long) ((endTime - startTime) * .001);
	}

	//Taken from class an changed to fit Comparable
	void mergeSort(int x[]) {

		int dest[100000];

		for(int i = 0; i < 100000; i++) {
			dest[i] = x[i];
		}

		//System.arraycopy(x,  0,  dest.toArray(), 0, x.size());
		realMergeSort(dest, x, 0, 100000);
	}

int isArraySorted(int s[], int n) {
  int a = 1, d = 1, i = 0;

  while ((a == 1 || d == 1) && i < n - 1) {
    if (s[i] < s[i+1])
      d = 0;
    else if (s[i] > s[i+1])
      a = 0;
    i++;
  }

  if (a == 1)
    return 1;
  else if (d == 1)
    return 2;
  else
    return 0;
}


int main (int argc, char *argv[]){

	int ScatterData[100000];

	int	taskid,	        // task ID - also used as seed number
	numtasks;       // number of tasks


	//Set seed for random number generator equal to task ID
	srandom (taskid);
	int index;
	for (index = 0; index < 100000; index++)
	{
		ScatterData[index] = rand();
	}

	// Obtain number of tasks and task ID
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

	numbertasks = numtasks;

	int ReceiveData[100000/numtasks];
	MPI_Scatter(ScatterData, 100000/numtasks, MPI_INT, &ReceiveData, 100000/numtasks, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	quickSort(ReceiveData);

	int GatheredData[100000];
	MPI_Gather(&ReceiveData, 100000/numtasks, MPI_INT, GatheredData,100000/numtasks, MPI_INT, 0, MPI_COMM_WORLD);

	mergeSort(GatheredData);

	int r = isArraySorted(GatheredData, 100000);

	if (r == 1)
		printf("The array is sorted in ascending order.\n");
	else
		printf("The array isn't sorted.\n");


	MPI_Finalize();

	return 0;
}
