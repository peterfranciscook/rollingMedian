/*
* MATLAB/MEX wrapper for fast rolling median by ashelly
* based on min/max heap
* see https://gist.github.com/ashelly/5665911
* also https://stackoverflow.com/questions/5527437/rolling-median-in-c-turlach-implementation
* peter cook 2018
* to compile in MATLAB type "mex -V -O -largeArrayDims rollingMedian.c"
*/

#include <stdlib.h>
#include <stdint.h>
#include "mex.h"
#include "matrix.h"
#include "mediator.h"


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	// A : input array (M x N) or (M x N x 3)
	// M : # rows in A
	// N : # cols in A
	// R : # rows in filter window
	// C : # cols in filter window

	int (*cmp)(const void* key, const void* elt);
	void (*mean)(const void* a, const void* b, void* c);
	size_t sz = 8;
	
	// step 0a: get dimensions of input array
	const size_t* dims = mxGetDimensions(prhs[0]);
	size_t nDim = mxGetNumberOfDimensions(prhs[0]);
	size_t M = dims[0];
	size_t N = dims[1];
	
	// step 0b: check dimensions of input array
	if ((M < 1) || (N < 1)) {
		// throw error
		mexErrMsgTxt("Leading Dimensions of A Must be Nonzero");
	}

	// step 1a: get dimensions of filter window
	size_t R = (size_t)mxGetScalar(prhs[1]);
	size_t C = (size_t)mxGetScalar(prhs[2]);
	// step 1b: check dimensions of filter window
	// arbitrarily impose the condition that R < M / 2
	// unstable behavior less likely with this constraint. may remove in future
	if (R > M / 2) {
		mexErrMsgTxt("R must be less than M/2");
	}
	// arbitrarily impose the condition that C < N / 2
	// unstable behavior less likely with this constraint. may remove in future
	if (C > N / 2) {
		mexErrMsgTxt("C must be less than N/2");
	}
	
	// step 2: check class of input array
	if (!mxIsNumeric(prhs[0])) {
		mexErrMsgTxt("Input Array (A) Must be Numeric");
	}
	else if (mxIsComplex(prhs[0])) {
		mexErrMsgTxt("Input Array (A) Cannot be Complex");
	}
	
	// step 2: check class of input array (cont.)
	if (mxIsDouble(prhs[0])) {
		cmp = cmp_double;
		mean = mean_double;
		sz = sizeof(double);
		plhs[0] = mxCreateNumericArray(nDim, dims, mxDOUBLE_CLASS, mxREAL);
	}
	else if (mxIsSingle(prhs[0])) {
		cmp = cmp_float;
		mean = mean_float;
		sz = sizeof(float);
		plhs[0] = mxCreateNumericArray(nDim, dims, mxSINGLE_CLASS, mxREAL);
	}
	else if (mxIsUint8(prhs[0])) {
		cmp = cmp_uint8;
		mean = mean_uint8;
		sz = sizeof(uint8_t);
		plhs[0] = mxCreateNumericArray(nDim, dims, mxUINT8_CLASS, mxREAL);
	}
	else if (mxIsInt8(prhs[0])) {
		cmp = cmp_int8;
		mean = mean_int8;
		sz = sizeof(int8_t);
		plhs[0] = mxCreateNumericArray(nDim, dims, mxINT8_CLASS, mxREAL);
	}
	else if (mxIsUint16(prhs[0])) {
		cmp = cmp_uint16;
		mean = mean_uint16;
		sz = sizeof(uint16_t);
		plhs[0] = mxCreateNumericArray(nDim, dims, mxUINT16_CLASS, mxREAL);
	}
	else if (mxIsInt16(prhs[0])) {
		cmp = cmp_int16;
		mean = mean_int16;
		sz = sizeof(int16_t);
		plhs[0] = mxCreateNumericArray(nDim, dims, mxINT16_CLASS, mxREAL);
	}
	else if (mxIsUint32(prhs[0])) {
		cmp = cmp_uint32;
		mean = mean_uint32;
		sz = sizeof(uint32_t);
		plhs[0] = mxCreateNumericArray(nDim, dims, mxUINT32_CLASS, mxREAL);
	}	
	else if (mxIsInt32(prhs[0])) {
		cmp = cmp_int32;
		mean = mean_int32;
		sz = sizeof(int32_t);
		plhs[0] = mxCreateNumericArray(nDim, dims, mxINT32_CLASS, mxREAL);
	}
	else if (mxIsUint64(prhs[0])) {
		cmp = cmp_uint64;
		mean = mean_uint64;
		sz = sizeof(uint64_t);
		plhs[0] = mxCreateNumericArray(nDim, dims, mxUINT64_CLASS, mxREAL);
	}
	else if (mxIsInt64(prhs[0])) {
		cmp = cmp_int64;
		mean = mean_int64;
		sz = sizeof(int64_t);
		plhs[0] = mxCreateNumericArray(nDim, dims, mxINT64_CLASS, mxREAL);
	}
	else {
		// throw error, class not supported
		const char classString[] = "uint8_t, int8_t, uint16_t, int16_t, uint32_t, int, uint64_t, int64_t, float, double";
		const char errHead[] = "Unsupported Numeric Class. Numeric Classes Supported";
		char errString[_MAX_PATH];
		sprintf(errString, "%s : %s", errHead, classString);
		mexErrMsgTxt(errString);
		return;
	}
	
	// determine the number of dimensions of the input array
	// some common configurations (non-exhaustive)
	// M x N : most probable
	// M x N x 3 : RGB Image
	// M x N x 4 : RGB Image with transparency mask
	// M x N x ? : Stack of 1-Channel Images
	// M x N x 3|4 x ? : Stack of RGB Images
	// M x N x ? x ? x ? : Getting a little ridiculous
	size_t nFrame = 1;
	if (nDim > 2) {
		for (size_t k = 2; k < nDim; k++) {
			nFrame *= dims[k];
		}
	}
	size_t k = 0;
	const void* A = (const void*)mxGetData(prhs[0]);
	do {
		medFilt(A, (unsigned char*)mxGetData(plhs[0]) + (k * M * N * sz), M, N, R, C, sz, cmp, mean);
		k++;
	} while (k < nFrame);
	
}


///* gateway function */
//void mexFunction(int nlhs, mxArray* plhs[],
//	int nrhs, const mxArray* prhs[])
//{
//	//int32_t M, N, R, C, i, j, k, rowIdx, colIdx;
//	//double* A, * B;
//	//Mediator* m;
//
//	///* get input array dimensions M x N */
//	////dims = (mwSize *)mxGetDimensions(prhs[0]);
//	//M = (int32_t)mxGetM(prhs[0]);
//	//N = (int32_t)mxGetN(prhs[0]);
//
//	///* get filter window dimensions R x C */
//	//R = (int32_t)mxGetScalar(prhs[1]);
//	//C = (int32_t)mxGetScalar(prhs[2]);
//
//	///* pointer to input array */
//	//A = (double*)mxGetData(prhs[0]);
//
//	///* pointer to output array */
//	//plhs[0] = mxCreateNumericMatrix(M, N, mxDOUBLE_CLASS, mxREAL);
//	//B = (double*)mxGetData(plhs[0]);
//
//	//// Filter Window Passes
//	//// 1. First C/2-1 Columns, Rows 0 to M-R/2
//	//// 2. Columns C/2 to N-C/2-1, Rows 0 to M-R/2
//	//// 3. Columns N-C/2 to N-1, Rows 0 to M-R/2
//	//// 4. Rows M-R/2 to M-1, Columns 0 to N-C/2-1
//	//// 5. Rows M-R/2 to M-1, Columns N-C/2 to N-1
//	//// 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3
//	//// 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3
//	//// 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3
//	//// 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3
//	//// 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3
//	//// 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3
//	//// 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3
//	//// 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3
//	//// 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3
//	//// 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3
//	//// 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3
//	//// 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3
//	//// 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3
//	//// 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3
//	//// 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3
//	//// 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3
//	//// 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3
//	//// 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5
//	//// 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5
//	//// 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5
//	//// 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5
//	//// 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5
//
//	////loop over the first C/2 columns - each one will have a different window size
//	//for (colIdx = 0; colIdx < C / 2; colIdx++)
//	//{
//	//	//init & partially fill heap then compute median for first column element
//	//	m = MediatorNew(R * (1 + C / 2 + colIdx));
//	//	for (rowIdx = 0; rowIdx < R / 2 + 1; rowIdx++)
//	//	{
//	//		for (j = 0; j < C / 2 + colIdx + 1; j++)
//	//		{
//	//			k = j * M + rowIdx;
//	//			MediatorInsert(m, A[k]);
//	//		}
//	//	}
//	//	B[M * colIdx] = MediatorMedian(m);
//	//	////debug
//	//	//mexPrintf("i = %i\n", M*colIdx);
//	//	//ShowTree(m);
//	//	////end debug
//
//	//	// for rows 1 to M-R/2-1 add C/2+1+colIdx elements to heap & compute median
//	//	for (rowIdx = 1; rowIdx < M - R / 2; rowIdx++)
//	//	{
//	//		for (j = 0; j < C / 2 + colIdx + 1; j++)
//	//		{
//	//			k = j * M + rowIdx + R / 2;
//	//			MediatorInsert(m, A[k]);
//	//		}
//	//		B[M * colIdx + rowIdx] = MediatorMedian(m);
//	//		////debug
//	//		//mexPrintf("i = %i\n", M*colIdx + rowIdx);
//	//		//ShowTree(m);
//	//		////end debug
//	//	}
//
//	//	//this (PopOldest) seems to cause unpredictable memory access violations
//	//	////for rows M-R/2 to M-1 pop C/2+1+colIdx elements from heap & compute median
//	//	//for (rowIdx = M - R / 2; rowIdx < M; rowIdx++)
//	//	//{
//	//	//	for (j = 0; j < C / 2 + colIdx + 1; j++)
//	//	//	{
//	//	//		PopOldest(m);
//	//	//		mexPrintf("Popping %ith Element from Mediator\n", j);
//	//	//	}
//	//	//	B[M*colIdx + rowIdx] = MediatorMedian(m);
//	//	//}
//	//	mxFree(m);
//	//}
//
//	////loop over the middle columns: C/2 to N-1-C/2
//	//for (colIdx = C / 2; colIdx < N - C / 2; colIdx++)
//	//{
//	//	//init & partially fill heap then compute median for first column element
//	//	m = MediatorNew(R * C);
//	//	for (rowIdx = 0; rowIdx < R / 2 + 1; rowIdx++)
//	//	{
//	//		for (j = colIdx - C / 2; j < colIdx + C / 2 + 1; j++)
//	//		{
//	//			k = j * M + rowIdx;
//	//			MediatorInsert(m, A[k]);
//	//		}
//	//	}
//	//	B[M * colIdx] = MediatorMedian(m);
//
//	//	for (rowIdx = 1; rowIdx < M - R / 2; rowIdx++)
//	//	{
//	//		for (j = colIdx - C / 2; j < colIdx + C / 2 + 1; j++)
//	//		{
//	//			k = j * M + rowIdx + R / 2;
//	//			MediatorInsert(m, A[k]);
//	//		}
//
//	//		B[M * colIdx + rowIdx] = MediatorMedian(m);
//	//	}
//	//	mxFree(m);
//	//}
//
//	////loop over the last C/2 columns - each one will have a different window size
//	//for (colIdx = N - C / 2; colIdx < N; colIdx++)
//	//{
//	//	//init & partially fill heap then compute median for first column element
//	//	m = MediatorNew(R * (1 + C / 2 + colIdx));
//	//	for (rowIdx = 0; rowIdx < R / 2 + 1; rowIdx++)
//	//	{
//	//		for (j = colIdx - C / 2; j < N; j++)
//	//		{
//	//			k = j * M + rowIdx;
//	//			MediatorInsert(m, A[k]);
//	//		}
//	//	}
//	//	B[M * colIdx] = MediatorMedian(m);
//
//	//	// for rows 1 to M-R/2-1 add C/2+1+colIdx elements to heap & compute median
//	//	for (rowIdx = 1; rowIdx < M - R / 2; rowIdx++)
//	//	{
//	//		for (j = colIdx - C / 2; j < N; j++)
//	//		{
//	//			k = j * M + rowIdx + R / 2;
//	//			MediatorInsert(m, A[k]);
//	//		}
//	//		B[M * colIdx + rowIdx] = MediatorMedian(m);
//	//	}
//	//	mxFree(m);
//	//}
//
//
//	//for (rowIdx = M - R / 2; rowIdx < M; rowIdx++)
//	//{
//	//	//init & partially fill heap then compute median for first column element
//	//	m = MediatorNew(C * (M + R / 2 - rowIdx));
//
//	//	for (colIdx = 0; colIdx < C / 2 + 1; colIdx++)
//	//	{
//	//		for (i = rowIdx - R / 2; i < M; i++)
//	//		{
//	//			k = i + colIdx * M;
//	//			MediatorInsert(m, A[k]);
//	//		}
//	//	}
//	//	B[rowIdx] = MediatorMedian(m);
//
//	//	//between columns C/2+1 to N-C/2 add (M+R/2-rowIdx) elements to heap & compute median
//	//	for (colIdx = 1; colIdx < N - C / 2; colIdx++)
//	//	{
//	//		for (i = rowIdx - R / 2; i < M; i++)
//	//		{
//	//			k = i + (colIdx + C / 2) * M;
//	//			MediatorInsert(m, A[k]);
//	//		}
//	//		B[rowIdx + M * colIdx] = MediatorMedian(m);
//	//	}
//	//	mxFree(m);
//
//	//	//init & partially fill heap then compute median for last column element
//	//	m = MediatorNew(C * (M + R / 2 - rowIdx));
//	//	for (colIdx = N - 1; colIdx > N - C / 2 - 2; colIdx--)
//	//	{
//	//		for (i = rowIdx - R / 2; i < M; i++)
//	//		{
//	//			k = i + colIdx * M;
//	//			MediatorInsert(m, A[k]);
//	//		}
//	//	}
//	//	B[rowIdx + M * (N - 1)] = MediatorMedian(m);
//
//	//	//slide backwards between columns N-2 and N-C/2
//	//	for (colIdx = N - 2; colIdx > N - C / 2 - 1; colIdx--)
//	//	{
//	//		for (i = rowIdx - R / 2; i < M; i++)
//	//		{
//	//			k = i + (colIdx - C / 2) * M;
//	//			MediatorInsert(m, A[k]);
//	//		}
//	//		B[rowIdx + M * colIdx] = MediatorMedian(m);
//	//	}
//	//	mxFree(m);
//	//}
//
//	return;
//}


	//// part 1: filter the left and right edges of the array
	//// outer loop : columns 0 to C/2-1
	//for (int n = 0; n < C / 2; n++) {
	//	//int J = C / 2 + n + (C % 2); // not sure about the mod2, could also just be + 1 ?
	//	int J = C / 2 + n + 1;
	//	int I = R;
	//	int K = I * (J - 1);
	//	// init and fill 4 mediators
	//	//Mediator* mediatorPtrNW = MediatorNew(K, sz);
	//	//Mediator* mediatorPtrSW = MediatorNew(K, sz);
	//	//Mediator* mediatorPtrNE = MediatorNew(K, sz);
	//	//Mediator* mediatorPtrSE = MediatorNew(K, sz);
	//	
	//	for (int c = 0; c < J; c++) {
	//		int v = N - 1 - c;
	//		for (int m = 0; m < R / 2 + 1; m++) {
	//			int u = M - 1 - m;
	//			//MediatorInsert(mediatorPtrNW, (unsigned char*)A + (c * M + m) * sz, cmp);
	//			//MediatorInsert(mediatorPtrSW, (unsigned char*)A + (c * M + u) * sz, cmp);
	//			//MediatorInsert(mediatorPtrNE, (unsigned char*)A + (v * M + m) * sz, cmp);
	//			//MediatorInsert(mediatorPtrSE, (unsigned char*)A + (v * M + u) * sz, cmp);
	//		}
	//	}
	//	// pop median from heap and copy to B
	//	//w = MediatorMedian(mediatorPtrNW, u, mean);
	//	//memcpy(B + M * n * sz, w, sz);
	//	//w = MediatorMedian(mediatorPtrSW, u, mean);
	//	//memcpy(B + (M * (n + 1) - 1) * sz, w, sz);
	//	//w = MediatorMedian(mediatorPtrNE, u, mean);
	//	//memcpy(B + M * (N - n - 1) * sz, w, sz);
	//	//w = MediatorMedian(mediatorPtrSE, u, mean);
	//	//memcpy(B + (M * (N - n) - 1) * sz, w, sz);
	//	
	//	// loop over rows i = 1 to M/2-1
	//	for (int m = R / 2 + 1; m < M / 2 + R / 2; m++) {
	//		int u = M - 1 - m;
	//		for (int c = 0; c < J; c++) {
	//			int v = N - 1 - c;
	//			//MediatorInsert(mediatorPtrNW, A[c * M + m], cmp);
	//			//MediatorInsert(mediatorPtrSW, A[c * M + u], cmp);
	//			//MediatorInsert(mediatorPtrNE, A[v * M + m], cmp);
	//			//MediatorInsert(mediatorPtrSE, A[v * M + u], cmp);
	//		}
	//		int i = m - R / 2;
	//		int j = n;
	//		// B[j * M + i] = MediatorMedian(mediatorPtrNW, u, mean);
	//		// B[j * M + u] = MediatorMedian(mediatorPtrSW, u, mean);
	//		// B[v * M + i] = MediatorMedian(mediatorPtrNE, u, mean);
	//		// B[v * M + u] = MediatorMedian(mediatorPtrSE, u, mean);
	//	}
	//	// free pointers
	//	//free(mediatorPtrNW);
	//	//free(mediatorPtrSW);
	//	//free(mediatorPtrNE);
	//	//free(mediatorPtrSE);
	//}

//	// part 2: filter the top and bottom edges of the array
//	// outer loop : rows 0 to R / 2 - 1 
//for (int m = 0; m < R / 2; m++) {
//	//int I = R / 2 + m + (R % 2); // not sure about the mod2, could also just be + 1 ?
//	int I = R / 2 + m + 1;
//	int J = C;
//	int K = (I - 1) * J;
//	// init and fill 4 mediators
//	// Mediator* mediatorPtrNW = MediatorNew(K, sz);
//	// Mediator* mediatorPtrSW = MediatorNew(K, sz);
//	// Mediator* mediatorPtrNE = MediatorNew(K, sz);
//	// Mediator* mediatorPtrSE = MediatorNew(K, sz);
//	//
//	// inner loop 1: cols
//	for (int n = 0; n < C; n++) {
//		// inner loop 2: rows
//		int v = N - 1 - n;
//		for (int r = 0; r < I; r++) {
//			int u = M - 1 - r;
//			//MediatorInsert(mediatorPtrNW, A[n * M + r], cmp);
//			//MediatorInsert(mediatorPtrSW, A[n * M + u], cmp);
//			//MediatorInsert(mediatorPtrNE, A[v * M + r], cmp);
//			//MediatorInsert(mediatorPtrSE, A[v * M + u], cmp);
//		}
//	}
//	int i = m;
//	int j = C / 2;
//	int u = M - 1 - i;
//	int v = N - 1 - j;
//	// B[j * M + i] = MediatorMedian(mediatorPtrNW, u, mean);
//	// B[j * M + u] = MediatorMedian(mediatorPtrSW, u, mean);
//	// B[v * M + i] = MediatorMedian(mediatorPtrNE, u, mean);
//	// B[v * M + u] = MediatorMedian(mediatorPtrSE, u, mean);
//
//	// inner loop 1: cols C to N/2-1
//	for (int n = C; n < N / 2 + C / 2; n++) {
//		// inner loop 2: rows 0 to R/2 + m
//		int v = N - 1 - n;
//		for (int r = 0; r < I; r++) {
//			int u = M - 1 - r;
//			//MediatorInsert(mediatorPtrNW, A[n * M + r], cmp);
//			//MediatorInsert(mediatorPtrSW, A[n * M + u], cmp);
//			//MediatorInsert(mediatorPtrNE, A[v * M + r], cmp);
//			//MediatorInsert(mediatorPtrSE, A[v * M + u], cmp);
//		}
//		int i = m;
//		int j = n - C / 2 + !(C % 2); // not 100% sure about the mod2
//		//int j = n - C / 2;
//		int u = M - 1 - i;
//		v = N - 1 - j;
//		// B[j * M + i] = MediatorMedian(mediatorPtrNW, u, mean);
//		// B[j * M + u] = MediatorMedian(mediatorPtrSW, u, mean);
//		// B[v * M + i] = MediatorMedian(mediatorPtrNE, u, mean);
//		// B[v * M + u] = MediatorMedian(mediatorPtrSE, u, mean);
//	}
//	// free pointers
//	//free(mediatorPtrNW);
//	//free(mediatorPtrSW);
//	//free(mediatorPtrNE);
//	//free(mediatorPtrSE);
//}