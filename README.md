# rollingMedian
C-MEX for 2D Rolling Median

MATLAB/MEX wrapper for fast rolling median by ashelly

based on min/max heap

see https://gist.github.com/ashelly/5665911
also https://stackoverflow.com/questions/5527437/rolling-median-in-c-turlach-implementation

to compile in MATLAB type "mex -V -O -largeArrayDims rollingMedian.c"

B = rollingMedian(A, R, C) Performs median filtering of the
matrix A in two dimensions with minimal edge effects and phase shift.
Inputs
------
A : Input Array
    Dimensions Allowed: (M x N), (M x N x ?), (M x N x ? x ?), ...
    As long as the leading dimensions of A (M & N) are nonzero, the
    filter will operate on all trailing dimensions. 
R : Filter Window Rows (1 < R < M / 2)
C : Filter Window Cols (1 < C < N / 2)

Outputs
-------
B : Output Array with the same dimensions and class as A.

Remarks
-------
rollingMedian uses a median-heap to compute the rolling median rather
than a sorting approach (i.e. sort all elements for each window). 
The time complexity of a sorting approach (for e.g. quicksort, mergesort) is 
O(M*N*R*C*log(R*C)). 
The time complexity of the median heap approach is O(M*N*log(R*C)).

Edge Effects
------------
The left and right edges (1) are filtered first using
successively wider filter windows for all pixels whose col index is less
than C/2. The top and bottom edges (2) are filtered second using
successively taller filter windows for all pixels whose row index is less
than R/2.

Phase Distortion
----------------
The algorithm operates on 4 pointers simultaneously (one for each of the
top-left, bottom-left, top-right, and bottom-right of the array) and
moves from the edges of the array inward. This creates a south-east phase
shift in the top-left quadrant, a north-east phase shift in the
bottom-left quadrant, a south-west phase shift for the top-right quadrant,
and a north-west phase shift in the bottom right quadrant. This may
create distortion at N/2 if C is even, and M/2 if R is even. If M or N is odd,
the median windows from both sides are advanced one row or col and the average
of both sides is used. 

Filter Window Passes
--------------------
```
1a: cols 0 to C/2-1, rows 0 to M/2-1
1b: cols 0 to C/2-1, rows M-1 to M-M/2 (reverse)
1c: cols N-1 to N-C/2 (reverse), rows 0 to M/2-1
1d: cols N-1 to N-C/2 (reverse), rows M-1 to M-M/2 (reverse)
1B: if M%2 : (cols 0 to C/2-1, row M/2) & (col N-1 to N-C/2 (reverse), row M/2)
2a: cols C/2 to N/2-1, rows 0 to R/2-1
2b: cols C/2 to N/2-1, rows M-1 to M-R/2 (reverse)
2c: cols N-C/2-1 to N-N/2 (reverse), rows 0 to R/2-1
2d: cols N-C/2-1 to N-N/2 (reverse), rows M-1 to M-R/2 (reverse)
2B: if N%2 : (col N/2, rows 0 to R/2-1) & (cols N/2, rows M-1 to M-R/2 (reverse))
3a: cols C/2 to N/2-1, rows R/2 to M/2-1
3b: cols C/2 to N/2-1, rows M-R/2-1 to M-M/2 (reverse)
3c: cols N-C/2-1 to N-N/2 (reverse), rows M/2 to M/2-1
3d: cols N-C/2-1 to N-N/2 (reverse), rows M-R/2-1 to M-M/2 (reverse)
3B: if N%2 : (col N/2, row R/2 to M/2-1) & (col N/2, rows M-R/2-1 to M-M/2 (reverse))
3C: if M%2 : (cols C/2 to N/2-1, row M/2) & (cols N-C/2-1 to N-N/2 (reverse), row M/2)
3D: if M%2 & N%2: average of 3B & 3C at (col N/2, row M/2)
   -> -> -> -> -> -> -> -> -> -> -> -> v  <- <- <- <- <- <- <- <- <- <- <- <-  
v  1  1  1  2  2  2  2  2  2  2  2  2 2B  2  2  2  2  2  2  2  2  2  1  1  1 v
|  1  1  1  2  2  2  2  2  2  2  2  2 2B  2  2  2  2  2  2  2  2  2  1  1  1 |
v  1  1  1  3  3  3  3  3  3  3  3  3 3B  3  3  3  3  3  3  3  3  3  1  1  1 v
|  1  1  1  3  3  3  3  3  3  3  3  3 3B  3  3  3  3  3  3  3  3  3  1  1  1 |
v  1  1  1  3  3  3  3  3  3  3  3  3 3B  3  3  3  3  3  3  3  3  3  1  1  1 v
|  1  1  1  3  3  3  3  3  3  3  3  3 3B  3  3  3  3  3  3  3  3  3  1  1  1 |
v  1  1  1  3  3  3  3  3  3  3  3  3 3B  3  3  3  3  3  3  3  3  3  1  1  1 v
|  1  1  1  3  3  3  3  3  3  3  3  3 3B  3  3  3  3  3  3  3  3  3  1  1  1 |
v  1  1  1  3  3  3  3  3  3  3  3  3 3B  3  3  3  3  3  3  3  3  3  1  1  1 v
|  1  1  1  3  3  3  3  3  3  3  3  3 3B  3  3  3  3  3  3  3  3  3  1  1  1 |
v  1  1  1  3  3  3  3  3  3  3  3  3 3B  3  3  3  3  3  3  3  3  3  1  1  1 v
-> 1B 1B 1B 3C3C 3C 3C 3C 3C 3C 3C 3C 3D 3C 3C 3C 3C 3C 3C 3C 3C 3C 1B 1B 1B<-
^  1  1  1  3  3  3  3  3  3  3  3  3 3B  3  3  3  3  3  3  3  3  3  1  1  1 ^
|  1  1  1  3  3  3  3  3  3  3  3  3 3B  3  3  3  3  3  3  3  3  3  1  1  1 |
^  1  1  1  3  3  3  3  3  3  3  3  3 3B  3  3  3  3  3  3  3  3  3  1  1  1 ^
|  1  1  1  3  3  3  3  3  3  3  3  3 3B  3  3  3  3  3  3  3  3  3  1  1  1 |
^  1  1  1  3  3  3  3  3  3  3  3  3 3B  3  3  3  3  3  3  3  3  3  1  1  1 ^
|  1  1  1  3  3  3  3  3  3  3  3  3 3B  3  3  3  3  3  3  3  3  3  1  1  1 |
^  1  1  1  3  3  3  3  3  3  3  3  3 3B  3  3  3  3  3  3  3  3  3  1  1  1 ^
|  1  1  1  3  3  3  3  3  3  3  3  3 3B  3  3  3  3  3  3  3  3  3  1  1  1 |
^  1  1  1  3  3  3  3  3  3  3  3  3 3B  3  3  3  3  3  3  3  3  3  1  1  1 ^
|  1  1  1  2  2  2  2  2  2  2  2  2 2B  2  2  2  2  2  2  2  2  2  1  1  1 |
^  1  1  1  2  2  2  2  2  2  2  2  2 2B  2  2  2  2  2  2  2  2  2  1  1  1 ^
   -> -> -> -> -> -> -> -> -> -> -> -> ^  <- <- <- <- <- <- <- <- <- <- <- <-  
```
Class Support
-------------
uint8, int8, uint16, int16, uint32, int32, uint64, int64, float, double

Peter Cook 2019
