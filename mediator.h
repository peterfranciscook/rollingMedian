#pragma once
//Copyright (c) 2011 ashelly.myopenid.com under <http://www.opensource.org/licenses/mit-license>

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>

typedef struct Mediator_t
{
	unsigned char* data;  //circular queue of values
	int* pos;   //index into `heap` for each value
	int* heap;  //max/median/min heap holding indexes into `data`.
	int   N;     //allocated size.
	int   idx;   //position in circular queue
	int   ct;    //count of items in queue
	int size; //size of data in bytes
	int (*cmp)(const void* key, const void* elt); // comparison function for type
	void (*mean)(const void* a, const void* b, void* c); // mean function for type
} Mediator;

/*--- Helper Functions ---*/
#define ItemLess(a,b)  ((a)<(b))
#define ItemMean(a,b)  (((a)+(b))/2)
#define minCt(m) (((m)->ct-1)/2) //count of items in minheap
#define maxCt(m) (((m)->ct)/2)   //count of items in maxheap 

// typed comparison functions:
// return -1 if x < y
// return  0 if x == y
// return  1 if x > y
// #define COMPARE(x,y) ((x)<(y)?-1:((x)==(y)?0:1))
int cmp_int8(const void* a, const void* b) {
	return (int)ItemLess(*(int8_t*)a, *(int8_t*)b);
}
int cmp_uint8(const void* a, const void* b) {
	return (int)ItemLess(*(uint8_t*)a, *(uint8_t*)b);
}
int cmp_int16(const void* a, const void* b) {
	return (int)ItemLess(*(int16_t*)a, *(int16_t*)b);
}
int cmp_uint16(const void* a, const void* b) {
	return (int)ItemLess(*(uint16_t*)a, *(uint16_t*)b);
}
int cmp_int32(const void* a, const void* b) {
	return (int)ItemLess(*(int32_t*)a, *(int32_t*)b);
}
int cmp_uint32(const void* a, const void* b) {
	return (int)ItemLess(*(uint32_t*)a, *(uint32_t*)b);
}
int cmp_int64(const void* a, const void* b) {
	return (int)ItemLess(*(int64_t*)a, *(int64_t*)b);
}
int cmp_uint64(const void* a, const void* b) {
	return (int)ItemLess(*(uint64_t*)a, *(uint64_t*)b);
}
int cmp_float(const void* a, const void* b) {
	return (int)ItemLess(*(float*)a, *(float*)b);
}
int cmp_double(const void* a, const void* b) {
	return (int)ItemLess(*(double*)a, *(double*)b);
}
// typed averaging functions
void mean_int8(const void* a, const void* b, void* c) {
	*(int8_t*)c = ItemMean(*(int8_t*)a, *(int8_t*)b);
}
void mean_uint8(const void* a, const void* b, void* c) {
	*(uint8_t*)c = ItemMean(*(uint8_t*)a, *(uint8_t*)b);
}
void mean_int16(const void* a, const void* b, void* c) {
	*(int16_t*)c = ItemMean(*(int16_t*)a, *(int16_t*)b);
}
void mean_uint16(const void* a, const void* b, void* c) {
	*(uint16_t*)c = ItemMean(*(uint16_t*)a, *(uint16_t*)b);
}
void mean_int32(const void* a, const void* b, void* c) {
	*(int32_t*)c = ItemMean(*(int32_t*)a, *(int32_t*)b);
}
void mean_uint32(const void* a, const void* b, void* c) {
	*(uint32_t*)c = ItemMean(*(uint32_t*)a, *(uint32_t*)b);
}
void mean_int64(const void* a, const void* b, void* c) {
	*(int64_t*)c = ItemMean(*(int64_t*)a, *(int64_t*)b);
}
void mean_uint64(const void* a, const void* b, void* c) {
	*(uint64_t*)c = ItemMean(*(uint64_t*)a, *(uint64_t*)b);
}
void mean_float(const void* a, const void* b, void* c) {
	*(float*)c = ItemMean(*(float*)a, *(float*)b);
}
void mean_double(const void* a, const void* b, void* c) {
	*(double*)c = ItemMean(*(double*)a, *(double*)b);
}

//returns 1 if heap[i] < heap[j]
int mmless(Mediator* m, int i, int j)
{
	int u = m->heap[i];
	int v = m->heap[j];
	int sz = m->size;
	const void* x = (m->data + u * sz);
	const void* y = (m->data + v * sz);
	return m->cmp(x, y);
}

//swaps items i&j in heap, maintains indexes
int mmexchange(Mediator* m, int i, int j)
{
	int t = m->heap[i];
	m->heap[i] = m->heap[j];
	m->heap[j] = t;
	m->pos[m->heap[i]] = i;
	m->pos[m->heap[j]] = j;
	return 1;
}

//swaps items i&j if i<j;  returns true if swapped
int mmCmpExch(Mediator* m, int i, int j)
{
	return (mmless(m, i, j) && mmexchange(m, i, j));
}

//maintains minheap property for all items below i/2.
void minSortDown(Mediator* m, int i)
{
	for (; i <= minCt(m); i *= 2) {
		if (i > 1 && i < minCt(m) && mmless(m, i + 1, i)) {
			++i;
		}
		if (!mmCmpExch(m, i, i / 2)) {
			break;
		}
	}
}

//maintains maxheap property for all items below i/2. (negative indexes)
void maxSortDown(Mediator* m, int i)
{
	for (; i >= -maxCt(m); i *= 2) {
		if (i<-1 && i > -maxCt(m) && mmless(m, i, i - 1)) {
			--i;
		}
		if (!mmCmpExch(m, i / 2, i)) {
			break;
		}
	}
}

//maintains minheap property for all items above i, including median
//returns true if median changed
int minSortUp(Mediator* m, int i)
{
	while ((i > 0) && mmCmpExch(m, i, i / 2)) {
		i /= 2;
	}
	return (i == 0);
}

//maintains maxheap property for all items above i, including median
//returns true if median changed
int maxSortUp(Mediator* m, int i)
{
	while ((i < 0) && mmCmpExch(m, i / 2, i)) {
		i /= 2;
	}
	return (i == 0);
}

/*--- Public Interface ---*/


//creates new Mediator: to calculate `nItems` running median. 
//mallocs single block of memory, caller must free.
Mediator* MediatorNew(int nItems, int dataSize,
	int (*cmp)(const void* key, const void* elt),
	void (*mean)(const void* a, const void* b, void* c))
{
	int size = sizeof(Mediator) + nItems * (dataSize + sizeof(int) * 2);
	//Mediator* m = malloc(size);
	Mediator* m = calloc(size, 1);
	if (m == NULL) {
		return m;
	}
	m->data = (unsigned char*)(m + 1);
	m->pos = (int*)(m->data + nItems * dataSize);
	m->heap = m->pos + nItems + (nItems / 2); //points to middle of storage.
	m->N = nItems;
	m->ct = m->idx = 0;
	m->size = dataSize;
	m->cmp = cmp;
	m->mean = mean;
	while (nItems--) {
		//set up initial heap fill pattern: median,max,min,max,...
		m->pos[nItems] = ((nItems + 1) / 2) * ((nItems & 1) ? -1 : 1);
		m->heap[m->pos[nItems]] = nItems;
	}
	return m;
}


//Inserts item, maintains median in O(lg nItems)
void MediatorInsert(Mediator* m, const void* v)
{
	// isNew : 1 if filling the buffer, 0 if the buffer is full
	// k : current position in circular queue
	// p : heap index for position k
	// old : value currently at data position k
	// u : pointer to current value of data at position k
	// v : pointer to new value to insert into data position k
	// isBigger  : 1 if v > data[k], 0 otherwise
	// isSmaller : 1 if v < data[k], 0 otherwise
	int isNew = (m->ct < m->N);
	int k = m->idx;
	int p = m->pos[k];
	int sz = m->size;
	void* u = (void*)(m->data + k * sz);
	int isBigger, isSmaller = 0;
	isBigger = m->cmp(u, v);
	if (!isBigger) {
		isSmaller = m->cmp(v, u);
	}

	// insert new value
	memcpy(u, v, sz);
	// advance index of circular queue
	m->idx = (k + 1) % (m->N);
	// increment buffer population count
	m->ct += isNew;

	if (p > 0) {
		//new item is in minHeap
		if (!isNew && isBigger) {
			minSortDown(m, p * 2);
		}
		else if (minSortUp(m, p)) {
			maxSortDown(m, -1);
		}
	}
	else if (p < 0) {
		//new item is in maxheap
		if (!isNew && isSmaller) {
			maxSortDown(m, p * 2);
		}
		else if (maxSortUp(m, p)) {
			minSortDown(m, 1);
		}
	}
	else {
		//new item is at median
		if (maxCt(m)) {
			maxSortDown(m, -1);
		}
		if (minCt(m)) {
			minSortDown(m, 1);
		}
	}
}

//returns median item (or average of 2 when item count is even)
void* MediatorMedian(Mediator* m, void* u)
{
	int k = m->heap[0];
	int sz = m->size;
	void* v = m->data + k * sz;
	if ((m->ct & 1) == 0) {
		int k = m->heap[-1];
		m->mean(v, m->data + k * sz, u);
		return u;
	}
	return v;
}


///*--- Test Code ---*/
//void PrintMaxHeap(Mediator* m)
//{
//	int i;
//	int sz = m->size;
//	if (maxCt(m)) {
//		int k = m->heap[-1];
//		printf("Max: %3d", *(int*)(m->data + k * sz));
//		//printf("Max: %3d", m->data[m->heap[-1]]);
//	}
//	for (i = 2; i <= maxCt(m); ++i) {
//		int k = m->heap[-i];
//		printf("|%3d ", *(int*)(m->data + k * sz));
//		//printf("|%3d ", m->data[m->heap[-i]]);
//		if (++i <= maxCt(m)) {
//			int k = m->heap[-i];
//			printf("%3d", *(int*)(m->data + k * sz));
//			//printf("%3d", m->data[m->heap[-i]]);
//		}
//	}
//	printf("\n");
//}
//void PrintMinHeap(Mediator* m)
//{
//	int i;
//	int sz = m->size;
//	if (minCt(m)) {
//		int k = m->heap[1];
//		printf("Min: %3d", *(int*)(m->data + k * sz));
//		//printf("Min: %3d", m->data[m->heap[1]]);
//	}
//	for (i = 2; i <= minCt(m); ++i) {
//		int k = m->heap[i];
//		printf("|%3d ", *(int*)(m->data + k * sz));
//		//printf("|%3d ", m->data[m->heap[i]]);
//		if (++i <= minCt(m)) {
//			int k = m->heap[i];
//			printf("%3d", *(int*)(m->data + k * sz));
//			//printf("%3d", m->data[m->heap[i]]);
//		}
//	}
//	printf("\n");
//}
//void ShowTree(Mediator* m)
//{
//	PrintMaxHeap(m);
//	printf("Mid: %3d\n", *(int*)(m->data + (m->heap[0]) * (m->size)));
//	//printf("Mid: %3d\n", m->data[m->heap[0]]);
//	PrintMinHeap(m);
//	printf("\n");
//}
//
//void MediatorTest_int8(N) {
//	int8_t v;
//	int sz = sizeof(int8_t);
//	Mediator* m = MediatorNew(N, sz);
//	void* u = malloc(sz);
//	if (u == NULL) {
//		return;
//	}
//	for (int i = 0; i < 128; i++) {
//		v = (int8_t)((rand() - (1 << 14)) % 127);
//		MediatorInsert(m, &v, cmp_int8);
//		const void* w = MediatorMedian(m, u, mean_int8);
//		// use this print routine to copy into matlab and verify result
//		for (int k = 0; k < (m->N); k++) {
//			printf("%3i ", *(int8_t*)(m->data + k * sz));
//		}
//		printf("%i\n", *(int8_t*)w);
//	}
//	// free memory
//	free(m);
//	free(u);
//}
//
//void MediatorTest_uint8(N) {
//	uint8_t v;
//	int sz = sizeof(uint8_t);
//	Mediator* m = MediatorNew(N, sz);
//	void* u = malloc(sz);
//	if (u == NULL) {
//		return;
//	}
//	for (int i = 0; i < 128; i++) {
//		v = (uint8_t)(rand() % 127);
//		MediatorInsert(m, &v, cmp_uint8);
//		const void* w = MediatorMedian(m, u, mean_uint8);
//		// use this print routine to copy into matlab and verify result
//		for (int k = 0; k < (m->N); k++) {
//			printf("%3i ", *(uint8_t*)(m->data + k * sz));
//		}
//		printf("%i\n", *(uint8_t*)w);
//	}
//	// free memory
//	free(m);
//	free(u);
//}
//
//void MediatorTest_int16(N) {
//	int16_t v;
//	int sz = sizeof(int16_t);
//	Mediator* m = MediatorNew(N, sz);
//	void* u = malloc(sz);
//	if (u == NULL) {
//		return;
//	}
//	for (int i = 0; i < 128; i++) {
//		v = (int16_t)((rand() - (1 << 14)) % 127);
//		MediatorInsert(m, &v, cmp_int16);
//		const void* w = MediatorMedian(m, u, mean_int16);
//		// use this print routine to copy into matlab and verify result
//		for (int k = 0; k < (m->N); k++) {
//			printf("%3i ", *(int16_t*)(m->data + k * sz));
//		}
//		printf("%i\n", *(int16_t*)w);
//	}
//	// free memory
//	free(m);
//	free(u);
//}
//
//void MediatorTest_uint16(N) {
//	uint16_t v;
//	int sz = sizeof(uint16_t);
//	Mediator* m = MediatorNew(N, sz);
//	void* u = malloc(sz);
//	if (u == NULL) {
//		return;
//	}
//	for (int i = 0; i < 128; i++) {
//		v = (uint16_t)(rand() % 127);
//		MediatorInsert(m, &v, cmp_uint16);
//		const void* w = MediatorMedian(m, u, mean_uint16);
//		// use this print routine to copy into matlab and verify result
//		for (int k = 0; k < (m->N); k++) {
//			printf("%3i ", *(uint16_t*)(m->data + k * sz));
//		}
//		printf("%i\n", *(uint16_t*)w);
//	}
//	// free memory
//	free(m);
//	free(u);
//}
//
//void MediatorTest_int32(N) {
//	int v;
//	int sz = sizeof(int);
//	Mediator* m = MediatorNew(N, sz);
//	void* u = malloc(sz);
//	if (u == NULL) {
//		return;
//	}
//	for (int i = 0; i < 128; i++) {
//		v = (rand() - (1 << 14)) % 127;
//		MediatorInsert(m, &v, cmp_int32);
//		const void* w = MediatorMedian(m, u, mean_int32);
//		// use this print routine to copy into matlab and verify result
//		for (int k = 0; k < (m->N); k++) {
//			printf("%3i ", *(int*)(m->data + k * sz));
//		}
//		printf("%i\n", *(int*)w);
//	}
//	// free memory
//	free(m);
//	free(u);
//}
//
//void MediatorTest_uint32(N) {
//	uint32_t v;
//	int sz = sizeof(uint32_t);
//	Mediator* m = MediatorNew(N, sz);
//	void* u = malloc(sz);
//	if (u == NULL) {
//		return;
//	}
//	for (int i = 0; i < 128; i++) {
//		v = (uint32_t)(rand() % 127);
//		MediatorInsert(m, &v, cmp_uint32);
//		const void* w = MediatorMedian(m, u, mean_uint32);
//		// use this print routine to copy into matlab and verify result
//		for (int k = 0; k < (m->N); k++) {
//			printf("%3i ", *(uint32_t*)(m->data + k * sz));
//		}
//		printf("%i\n", *(uint32_t*)w);
//	}
//	// free memory
//	free(m);
//	free(u);
//}
//
//void MediatorTest_int64(N) {
//	int64_t v;
//	int sz = sizeof(int64_t);
//	Mediator* m = MediatorNew(N, sz);
//	void* u = malloc(sz);
//	if (u == NULL) {
//		return;
//	}
//	for (int i = 0; i < 128; i++) {
//		v = (int64_t)((rand() - (1 << 14)) % 127);
//		MediatorInsert(m, &v, cmp_int64);
//		const void* w = MediatorMedian(m, u, mean_int64);
//		// use this print routine to copy into matlab and verify result
//		for (int k = 0; k < (m->N); k++) {
//			printf("%3i ", *(int64_t*)(m->data + k * sz));
//		}
//		printf("%i\n", *(int64_t*)w);
//	}
//	// free memory
//	free(m);
//	free(u);
//}
//
//void MediatorTest_uint64(N) {
//	uint64_t v;
//	int sz = sizeof(uint64_t);
//	Mediator* m = MediatorNew(N, sz);
//	void* u = malloc(sz);
//	if (u == NULL) {
//		return;
//	}
//	for (int i = 0; i < 128; i++) {
//		v = (uint64_t)(rand() % 127);
//		MediatorInsert(m, &v, cmp_uint64);
//		const void* w = MediatorMedian(m, u, mean_uint64);
//		// use this print routine to copy into matlab and verify result
//		for (int k = 0; k < (m->N); k++) {
//			printf("%3i ", *(uint64_t*)(m->data + k * sz));
//		}
//		printf("%i\n", *(uint64_t*)w);
//	}
//	// free memory
//	free(m);
//	free(u);
//}
//
//void MediatorTest_float(N) {
//	float v;
//	int sz = sizeof(float);
//	Mediator* m = MediatorNew(N, sz);
//	void* u = malloc(sz);
//	if (u == NULL) {
//		return;
//	}
//	for (int i = 0; i < 128; i++) {
//		v = (float)((rand() - (1 << 14)) % 127);
//		MediatorInsert(m, &v, cmp_float);
//		const void* w = MediatorMedian(m, u, mean_float);
//		// use this print routine to copy into matlab and verify result
//		for (int k = 0; k < (m->N); k++) {
//			printf("%3.1f ", *(float*)(m->data + k * sz));
//		}
//		printf("%3.1f\n", *(float*)w);
//	}
//	// free memory
//	free(m);
//	free(u);
//}
//
//void MediatorTest_double(N) {
//	double v;
//	int sz = sizeof(double);
//	Mediator* m = MediatorNew(N, sz);
//	void* u = malloc(sz);
//	if (u == NULL) {
//		return;
//	}
//	for (int i = 0; i < 128; i++) {
//		v = (double)((rand() - (1 << 14)) % 127);
//		MediatorInsert(m, &v, cmp_double);
//		const void* w = MediatorMedian(m, u, mean_double);
//		// use this print routine to copy into matlab and verify result
//		for (int k = 0; k < (m->N); k++) {
//			printf("%3.1f ", *(double*)(m->data + k * sz));
//		}
//		printf("%3.1f\n", *(double*)w);
//	}
//	// free memory
//	free(m);
//	free(u);
//}
//
//
//int mediatorTest(N)
//{
//	MediatorTest_int8(N);
//	printf("\n");
//	MediatorTest_uint8(N);
//	printf("\n");
//	MediatorTest_int16(N);
//	printf("\n");
//	MediatorTest_uint16(N);
//	printf("\n");
//	MediatorTest_int32(N);
//	printf("\n");
//	MediatorTest_uint32(N);
//	printf("\n");
//	MediatorTest_int64(N);
//	printf("\n");
//	MediatorTest_uint64(N);
//	printf("\n");
//	MediatorTest_float(N);
//	printf("\n");
//	MediatorTest_double(N);
//	return 1;
//}

void medFilt(const void* A, unsigned char* B, size_t M, size_t N, size_t R, size_t C, size_t sz,
	int (*cmp)(const void* key, const void* elt),
	void (*mean)(const void* a, const void* b, void* c))
{
	// NWPtr : pointer to upper-left corner of A
	// SWPtr : pointer to lower-left corner of A
	// NEPtr : pointer to upper-right corner of A
	// SEPtr : pointer to lower-right corner of A
	// M : # rows in A
	// N : # cols in A
	// R : # rows in filter window
	// C : # cols in filter window
	// sz: sizeof(Type)
	// I : # rows in clipped/edge filter window
	// J : # cols in clipped/edge filter window
	// K : # elements in filter window
	// m : current row in A
	// n : current col in A
	// i : current row in B (output index)
	// j : current col in B (output index)
	// u : mirror row position of i (?)
	// v : mirror col position of j (?)

	//size_t sz = sizeof(double);
	void* uPtr = malloc(sz); // memory location to hold mean in the event of an even window
	void* vPtr = malloc(sz); // memory location to hold mean in the event of an even window
	void* wPtr = malloc(sz); // memory location to hold mean in the event of an even window
	const void* w;
	const void* x;
	const void* y;
	unsigned char* NWPtr = (unsigned char*)A;
	//unsigned char* SWPtr;
	//unsigned char* NEPtr;
	//unsigned char* SEPtr;

	// on each iteration we'll work on 4 points simultaneously
	// (i,j), (M-1-i,j), (i, N-1-j), (M-1-i,N-1-j)

	// part 1: filter the left and right edges of the array
	for (size_t j = 0; j < C / 2; j++) {
		// suppose C = 4 -> J = 2 + {0,1} + 0 = {2, 3} -> n = {0,1}, {0,1,2}
		// suppose C = 5 -> J = 2 + {0,1} + 1 = {3, 4} -> n = {0,1,2}, {0,1,2,3}
		size_t I = R;
		size_t J = C / 2 + j + (C % 2);
		size_t K = I * J;

		// init and fill 4 mediators
		Mediator* mediatorPtrNW = MediatorNew(K, sz, cmp, mean);
		Mediator* mediatorPtrSW = MediatorNew(K, sz, cmp, mean);
		Mediator* mediatorPtrNE = MediatorNew(K, sz, cmp, mean);
		Mediator* mediatorPtrSE = MediatorNew(K, sz, cmp, mean);
		// ...

		// 1.1: initial fill of median heap/buffer(s)
		for (size_t m = 0; m < R / 2 + (R % 2); m++) {
			// this is a FIFO queue so need to iterate rows then cols
			size_t u = M - 1 - m;
			// 1.1.1 insert J elements into buffer 
			for (size_t n = 0; n < J; n++) {
				// k : offset in bytes to element (m,n) of A
				size_t v = N - 1 - n;
				size_t k = (n * M + m) * sz;
				// TODO: make pointer arithmetic better/easier/faster
				MediatorInsert(mediatorPtrNW, NWPtr + k);
				k = (n * M + u) * sz;
				MediatorInsert(mediatorPtrSW, NWPtr + k);
				k = (v * M + m) * sz;
				MediatorInsert(mediatorPtrNE, NWPtr + k);
				k = (v * M + u) * sz;
				MediatorInsert(mediatorPtrSE, NWPtr + k);
			}
		}

		// 1.2: compute median of window around elements (0,j), (M-1,j), (0,N-1-j), (M-1,N-1-j)
		// k : offset in bytes to element (0,j) of B
		size_t u = M - 1;
		size_t v = N - 1 - j;
		// top left corner
		size_t k = sz * M * j;
		w = MediatorMedian(mediatorPtrNW, uPtr);
		memcpy(B + k, w, sz);
		// bottom left corner
		k = sz * (M * j + u);
		w = MediatorMedian(mediatorPtrSW, uPtr);
		memcpy(B + k, w, sz);
		// top right corner
		k = sz * M * v;
		w = MediatorMedian(mediatorPtrNE, uPtr);
		memcpy(B + k, w, sz);
		// bottom right corner
		k = sz * (M * v + u);
		w = MediatorMedian(mediatorPtrSE, uPtr);
		memcpy(B + k, w, sz);
		// ...

		// 1.3: compute window median for rows 1 to M/2-1
		for (size_t m = R / 2 + (R % 2); m < M / 2 + R / 2 + (R % 2) - 1; m++) {
			// 1.3.1: insert A[m,0], A[m,1] ... A[m,J-1] into buffer
			size_t u = M - 1 - m;
			for (size_t n = 0; n < J; n++) {
				// k : offset in bytes to element (m,n) of A
				size_t v = N - 1 - n;
				size_t k = (n * M + m) * sz;
				MediatorInsert(mediatorPtrNW, NWPtr + k);
				k = (n * M + u) * sz;
				MediatorInsert(mediatorPtrSW, NWPtr + k);
				k = (v * M + m) * sz;
				MediatorInsert(mediatorPtrNE, NWPtr + k);
				k = (v * M + u) * sz;
				MediatorInsert(mediatorPtrSE, NWPtr + k);
			}
			// 1.3.2: compute median of window around elements (i,j), (M-1-i,j), (i,N-1-j), (M-1-i,N-1-j)
			// i : row index of B corresponding to centroid of current window
			//     an additional index shift is applied for even-row-length windows
			// k : offset in bytes to element (i,j) of B
			size_t i = m - R / 2 + !(R % 2); // questionable but seems to work on paper
			u = M - 1 - i;
			size_t v = N - 1 - j;
			// top left corner
			size_t k = sz * (M * j + i);
			w = MediatorMedian(mediatorPtrNW, uPtr);
			memcpy(B + k, w, sz);
			// bottom left corner
			k = sz * (M * j + u);
			w = MediatorMedian(mediatorPtrSW, uPtr);
			memcpy(B + k, w, sz);
			// top right corner
			k = sz * (M * v + i);
			w = MediatorMedian(mediatorPtrNE, uPtr);
			memcpy(B + k, w, sz);
			// bottom right corner
			k = sz * (M * v + u);
			w = MediatorMedian(mediatorPtrSE, uPtr);
			memcpy(B + k, w, sz);
		}
		// 1.4: fill in middle row if needed
		if (M % 2) {
			size_t m = M / 2 + R / 2 + (R % 2) - 1;
			size_t u = M - 1 - m;
			// 1.4.1: insert A[m,0], A[m,1] ... A[m,J-1] into buffer
			for (size_t n = 0; n < J; n++) {
				// k : offset in bytes to element (m,n) of A
				size_t v = N - 1 - n;
				size_t k = (n * M + m) * sz;
				MediatorInsert(mediatorPtrNW, NWPtr + k);
				k = (n * M + u) * sz;
				MediatorInsert(mediatorPtrSW, NWPtr + k);
				k = (v * M + m) * sz;
				MediatorInsert(mediatorPtrNE, NWPtr + k);
				k = (v * M + u) * sz;
				MediatorInsert(mediatorPtrSE, NWPtr + k);
			}
			// 1.4.2: compute median of window around elements (i,j), (M-1-i,j), (i,N-1-j), (M-1-i,N-1-j)
			// i : row index of B corresponding to centroid of current window
			//     an additional index shift is applied for even-row-length windows
			// k : offset in bytes to element (i,j) of B
			size_t i = m - R / 2 + !(R % 2); // questionable but seems to work on paper
			size_t v = N - 1 - j;

			// average the values from the top-left and bottom-left
			size_t k = sz * (M * j + i);
			x = MediatorMedian(mediatorPtrNW, uPtr);
			y = MediatorMedian(mediatorPtrSW, vPtr);
			mean(x, y, wPtr);
			memcpy(B + k, wPtr, sz);

			// average the values from the top-right and bottom-right
			k = sz * (M * v + i);
			x = MediatorMedian(mediatorPtrNE, uPtr);
			y = MediatorMedian(mediatorPtrSE, vPtr);
			mean(x, y, wPtr);
			memcpy(B + k, wPtr, sz);
		}
		// free memory
		free(mediatorPtrNW);
		free(mediatorPtrSW);
		free(mediatorPtrNE);
		free(mediatorPtrSE);
	}

	// part 2: filter the top and bottom edges of the array
	for (size_t i = 0; i < R / 2; i++) {
		size_t I = R / 2 + i + (R % 2);
		size_t J = C;
		size_t K = I * J;

		// init and fill 4 mediators
		Mediator* mediatorPtrNW = MediatorNew(K, sz, cmp, mean);
		Mediator* mediatorPtrSW = MediatorNew(K, sz, cmp, mean);
		Mediator* mediatorPtrNE = MediatorNew(K, sz, cmp, mean);
		Mediator* mediatorPtrSE = MediatorNew(K, sz, cmp, mean);
		// ...

		// 2.1: initial fill of median heap/buffer(s)
		for (size_t n = 0; n < J; n++) {
			// 2.1.1 insert I elements into buffer
			size_t v = N - 1 - n;
			for (size_t m = 0; m < I; m++) {
				// k : offset in bytes to element (m,n) of A
				size_t u = M - 1 - m;
				size_t k = (n * M + m) * sz;
				MediatorInsert(mediatorPtrNW, NWPtr + k);
				k = (n * M + u) * sz;
				MediatorInsert(mediatorPtrSW, NWPtr + k);
				k = (v * M + m) * sz;
				MediatorInsert(mediatorPtrNE, NWPtr + k);
				k = (v * M + u) * sz;
				MediatorInsert(mediatorPtrSE, NWPtr + k);
			}
		}
		// 2.2: compute median of window around elements (i,C/2), (M-1,C/2), (0,N-1-C/2), (M-1,N-1-C/2)
		size_t j = C / 2;
		size_t u = M - 1 - i;
		size_t v = N - 1 - j;
		// top left corner
		size_t k = sz * (M * j + i);
		w = MediatorMedian(mediatorPtrNW, uPtr);
		memcpy(B + k, w, sz);
		// bottom left corner
		k = sz * (M * j + u);
		w = MediatorMedian(mediatorPtrSW, uPtr);
		memcpy(B + k, w, sz);
		// top right corner
		k = sz * (v * M + i);
		w = MediatorMedian(mediatorPtrNE, uPtr);
		memcpy(B + k, w, sz);
		// bottom right corner
		k = sz * (v * M + u);
		w = MediatorMedian(mediatorPtrSE, uPtr);
		memcpy(B + k, w, sz);

		// 2.3: compute window median for columns C/2 to N/2-1
		for (size_t n = J; n < N / 2 + C / 2 + (C % 2) - 1; n++) {
			size_t v = N - 1 - n;
			for (size_t m = 0; m < I; m++) {
				// k : offset in bytes to element (m,n) of A
				size_t u = M - 1 - m;
				size_t k = (n * M + m) * sz;
				MediatorInsert(mediatorPtrNW, NWPtr + k);
				k = (n * M + u) * sz;
				MediatorInsert(mediatorPtrSW, NWPtr + k);
				k = (v * M + m) * sz;
				MediatorInsert(mediatorPtrNE, NWPtr + k);
				k = (v * M + u) * sz;
				MediatorInsert(mediatorPtrSE, NWPtr + k);
			}
			// j : col index of B corresponding to centroid of current window
			//     an additional index shift is applied for even-row-length windows
			// k : offset in bytes to element (i,j) of B
			size_t j = n - C / 2 + !(C % 2);
			size_t u = M - 1 - i;
			v = N - 1 - j;
			// top left corner
			size_t k = sz * (M * j + i);
			w = MediatorMedian(mediatorPtrNW, uPtr);
			memcpy(B + k, w, sz);
			// bottom left corner
			k = sz * (M * j + u);
			w = MediatorMedian(mediatorPtrSW, uPtr);
			memcpy(B + k, w, sz);
			// top right corner
			k = (v * M + i) * sz;
			w = MediatorMedian(mediatorPtrNE, uPtr);
			memcpy(B + k, w, sz);
			// bottom right corner
			k = (v * M + u) * sz;
			w = MediatorMedian(mediatorPtrSE, uPtr);
			memcpy(B + k, w, sz);
		}
		// 2.4: fill in middle col if needed
		if (N % 2) {
			size_t n = N / 2 + C / 2 + (C % 2) - 1;
			size_t v = N - 1 - n;
			for (size_t m = 0; m < I; m++) {
				// k : offset in bytes to element (m,n) of A
				size_t u = M - 1 - m;
				size_t k = (n * M + m) * sz;
				MediatorInsert(mediatorPtrNW, NWPtr + k);
				k = (n * M + u) * sz;
				MediatorInsert(mediatorPtrSW, NWPtr + k);
				k = (v * M + m) * sz;
				MediatorInsert(mediatorPtrNE, NWPtr + k);
				k = (v * M + u) * sz;
				MediatorInsert(mediatorPtrSE, NWPtr + k);
			}
			// j : col index of B corresponding to centroid of current window
			//     an additional index shift is applied for even-row-length windows
			// k : offset in bytes to element (i,j) of B
			size_t j = n - C / 2 + !(C % 2);
			size_t u = M - 1 - i;

			// average values from the top-left and top-right
			size_t k = sz * (M * j + i);
			x = MediatorMedian(mediatorPtrNW, uPtr);
			y = MediatorMedian(mediatorPtrNE, vPtr);
			mean(x, y, wPtr);
			memcpy(B + k, wPtr, sz);

			// average values from the bottom-left and bottom-right
			k = sz * (M * j + u);
			x = MediatorMedian(mediatorPtrSW, uPtr);
			y = MediatorMedian(mediatorPtrSE, vPtr);
			mean(x, y, wPtr);
			memcpy(B + k, wPtr, sz);
		}
		// free memory
		free(mediatorPtrNW);
		free(mediatorPtrSW);
		free(mediatorPtrNE);
		free(mediatorPtrSE);
	}

	// part 3: filter the middle of the array
	// init 4 mediators
	size_t K = R * C;
	Mediator* mediatorPtrNW = MediatorNew(K, sz, cmp, mean);
	Mediator* mediatorPtrSW = MediatorNew(K, sz, cmp, mean);
	Mediator* mediatorPtrNE = MediatorNew(K, sz, cmp, mean);
	Mediator* mediatorPtrSE = MediatorNew(K, sz, cmp, mean);
	// ...
	for (size_t i = R / 2; i < M / 2; i++) {
		// 3.1: initial fill of median heap/buffer(s)
		for (size_t n = 0; n < C; n++) {
			size_t v = N - 1 - n;
			for (size_t m = i - R / 2; m < (i + R / 2 + (R % 2)); m++) {
				size_t u = M - 1 - m;
				size_t k = (n * M + m) * sz;
				MediatorInsert(mediatorPtrNW, NWPtr + k);
				k = (n * M + u) * sz;
				MediatorInsert(mediatorPtrSW, NWPtr + k);
				k = (v * M + m) * sz;
				MediatorInsert(mediatorPtrNE, NWPtr + k);
				k = (v * M + u) * sz;
				MediatorInsert(mediatorPtrSE, NWPtr + k);
			}
		}
		// 3.2: compute median of window around elements (i,j), (M-1-1,j), (i,N-1-j), (M-1-1,N-1-j)
		size_t j = C / 2;
		size_t u = M - 1 - i;
		size_t v = N - 1 - j;
		// top right corner
		size_t k = (M * j + i) * sz;
		w = MediatorMedian(mediatorPtrNW, uPtr);
		memcpy(B + k, w, sz);
		// bottom left corner
		k = (j * M + u) * sz;
		w = MediatorMedian(mediatorPtrSW, uPtr);
		memcpy(B + k, w, sz);
		// top right corner
		k = (v * M + i) * sz;
		w = MediatorMedian(mediatorPtrNE, uPtr);
		memcpy(B + k, w, sz);
		// bottom right corner
		k = (v * M + u) * sz;
		w = MediatorMedian(mediatorPtrSE, uPtr);
		memcpy(B + k, w, sz);

		// 3.3: compute window median for columns C/2 to N/2-1
		for (size_t n = C; n < N / 2 + C / 2 + (C % 2) - 1; n++) {
			size_t v = N - 1 - n;
			for (size_t m = i - R / 2; m < (i + R / 2 + (R % 2)); m++) {
				size_t u = M - 1 - m;
				size_t k = (n * M + m) * sz;
				MediatorInsert(mediatorPtrNW, NWPtr + k);
				k = (n * M + u) * sz;
				MediatorInsert(mediatorPtrSW, NWPtr + k);
				k = (v * M + m) * sz;
				MediatorInsert(mediatorPtrNE, NWPtr + k);
				k = (v * M + u) * sz;
				MediatorInsert(mediatorPtrSE, NWPtr + k);
			}
			size_t j = n - C / 2 + !(C % 2);
			size_t u = M - 1 - i;
			v = N - 1 - j;
			// top left corner
			size_t k = (M * j + i) * sz;
			w = MediatorMedian(mediatorPtrNW, uPtr);
			memcpy(B + k, w, sz);
			// bottom left corner
			k = (j * M + u) * sz;
			w = MediatorMedian(mediatorPtrSW, uPtr);
			memcpy(B + k, w, sz);
			// top right corner
			k = (v * M + i) * sz;
			w = MediatorMedian(mediatorPtrNE, uPtr);
			memcpy(B + k, w, sz);
			// bottom right corner
			k = (v * M + u) * sz;
			w = MediatorMedian(mediatorPtrSE, uPtr);
			memcpy(B + k, w, sz);
		}
		// 3.4: fill in middle col if needed
		if (N % 2) {
			size_t n = N / 2 + C / 2 + (C % 2) - 1;
			size_t v = N - 1 - n;
			for (size_t m = i - R / 2; m < (i + R / 2 + (R % 2)); m++) {
				size_t u = M - 1 - m;
				size_t k = (n * M + m) * sz;
				MediatorInsert(mediatorPtrNW, NWPtr + k);
				k = (n * M + u) * sz;
				MediatorInsert(mediatorPtrSW, NWPtr + k);
				k = (v * M + m) * sz;
				MediatorInsert(mediatorPtrNE, NWPtr + k);
				k = (v * M + u) * sz;
				MediatorInsert(mediatorPtrSE, NWPtr + k);
			}
			size_t j = n - C / 2 + !(C % 2);
			size_t u = M - 1 - i;
			v = N - 1 - j;
			// average values from the top-left and top-right
			size_t k = (M * j + i) * sz;
			x = MediatorMedian(mediatorPtrNW, uPtr);
			y = MediatorMedian(mediatorPtrNE, vPtr);
			mean(x, y, wPtr);
			memcpy(B + k, wPtr, sz);

			// average values from the bottom-left and bottom-right
			k = (j * M + u) * sz;
			x = MediatorMedian(mediatorPtrSW, uPtr);
			y = MediatorMedian(mediatorPtrSE, vPtr);
			mean(x, y, wPtr);
			memcpy(B + k, wPtr, sz);
		}
	}
	// 3.5: fill in middle row if needed
	if (M % 2) {
		size_t i = M / 2;
		// 3.5.1: initial fill of median heap/buffer(s)
		for (size_t n = 0; n < C; n++) {
			size_t v = N - 1 - n;
			for (size_t m = i - R / 2; m < (i + R / 2 + (R % 2)); m++) {
				size_t u = M - 1 - m;
				size_t k = (n * M + m) * sz;
				MediatorInsert(mediatorPtrNW, NWPtr + k);
				k = (n * M + u) * sz;
				MediatorInsert(mediatorPtrSW, NWPtr + k);
				k = (v * M + m) * sz;
				MediatorInsert(mediatorPtrNE, NWPtr + k);
				k = (v * M + u) * sz;
				MediatorInsert(mediatorPtrSE, NWPtr + k);
			}
		}
		// 3.5.2: compute median of window around elements (i,j), (M-1-1,j), (i,N-1-j), (M-1-1,N-1-j)
		size_t j = C / 2;
		size_t u = M - 1 - i;
		size_t v = N - 1 - j;

		// average values from top-left and bottom-left
		size_t k = (M * j + i) * sz;
		x = MediatorMedian(mediatorPtrNW, uPtr);
		y = MediatorMedian(mediatorPtrSW, vPtr);
		mean(x, y, wPtr);
		memcpy(B + k, wPtr, sz);

		// average values from top-right and bottom-right
		k = (v * M + i) * sz;
		x = MediatorMedian(mediatorPtrNE, uPtr);
		y = MediatorMedian(mediatorPtrSE, vPtr);
		mean(x, y, wPtr);
		memcpy(B + k, wPtr, sz);

		// 3.5.3: compute window median for columns C/2 to N/2-1
		for (size_t n = C; n < N / 2 + C / 2 + (C % 2) - 1; n++) {
			size_t v = N - 1 - n;
			for (size_t m = i - R / 2; m < (i + R / 2 + (R % 2)); m++) {
				size_t u = M - 1 - m;
				size_t k = (n * M + m) * sz;
				MediatorInsert(mediatorPtrNW, NWPtr + k);
				k = (n * M + u) * sz;
				MediatorInsert(mediatorPtrSW, NWPtr + k);
				k = (v * M + m) * sz;
				MediatorInsert(mediatorPtrNE, NWPtr + k);
				k = (v * M + u) * sz;
				MediatorInsert(mediatorPtrSE, NWPtr + k);
			}
			size_t j = n - C / 2 + !(C % 2);
			size_t u = M - 1 - i;
			v = N - 1 - j;

			// average values from top-left and bottom-left
			size_t k = (M * j + i) * sz;
			x = MediatorMedian(mediatorPtrNW, uPtr);
			y = MediatorMedian(mediatorPtrSW, vPtr);
			mean(x, y, wPtr);
			memcpy(B + k, wPtr, sz);

			// average values from top-right and bottom-right
			k = (v * M + i) * sz;
			x = MediatorMedian(mediatorPtrNE, uPtr);
			y = MediatorMedian(mediatorPtrSE, vPtr);
			mean(x, y, wPtr);
			memcpy(B + k, wPtr, sz);
		}
		// 3.5.4 fill in center pixel if both N and M are odd
		if (N % 2) {
			size_t n = N / 2 + C / 2 + (C % 2) - 1;
			size_t v = N - 1 - n;
			for (size_t m = i - R / 2; m < (i + R / 2 + (R % 2)); m++) {
				size_t u = M - 1 - m;
				size_t k = (n * M + m) * sz;
				MediatorInsert(mediatorPtrNW, NWPtr + k);
				k = (n * M + u) * sz;
				MediatorInsert(mediatorPtrSW, NWPtr + k);
				k = (v * M + m) * sz;
				MediatorInsert(mediatorPtrNE, NWPtr + k);
				k = (v * M + u) * sz;
				MediatorInsert(mediatorPtrSE, NWPtr + k);
			}
			size_t j = n - C / 2 + !(C % 2);
			size_t u = M - 1 - i;
			v = N - 1 - j;
			// average values from the top-left and top-right
			size_t k = (M * j + i) * sz;
			x = MediatorMedian(mediatorPtrNW, uPtr);
			y = MediatorMedian(mediatorPtrNE, vPtr);
			mean(x, y, wPtr);
			memcpy(B + k, wPtr, sz);

			// average values from the bottom-left and bottom-right
			k = (j * M + u) * sz;
			x = MediatorMedian(mediatorPtrSW, uPtr);
			y = MediatorMedian(mediatorPtrSE, vPtr);
			mean(x, y, wPtr);

			// average the two averages
			mean(B + k, wPtr, uPtr);
			memcpy(B + k, uPtr, sz);
		}
	}
	// free memory
	free(mediatorPtrNW);
	free(mediatorPtrSW);
	free(mediatorPtrNE);
	free(mediatorPtrSE);
	free(uPtr);
	free(vPtr);
	free(wPtr);
	// ...
}