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
