#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/time.h>
#include<omp.h>
#include <immintrin.h>

#define likely(x) __builtin_expect(!!(x), 1) 
#define unlikely(x) __builtin_expect(!!(x), 0)

// Calculate sum of distance while combining different pivots. Complexity : O( n^2 )
double SumDistance(const int k, const int n, const int dim, double* coord, int* pivots){
    //double* rebuiltCoord = (double*)malloc(sizeof(double) * n * k);
    double* rebuiltCoord = (double*)_mm_malloc(sizeof(double) * n * k,32);
     //double* rebuiltCoord = (double*)mmap(NULL,sizeof(double) * n * k  , PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS , -1 , 0);
    double temp[5];
    // Rebuild coordinates. New coordinate of one point is its distance to each pivot.
        for(int ki=0; ki<k; ++ki){
            int pivoti = pivots[ki];
            for(int i=0; i<n; ++i){
                double distance = 0;
                for(int j=0; j<dim; ++j){
                    double temp = (coord[pivoti*dim+j] - coord[i*dim+j]);
                    distance += temp*temp;
               // distance += 
               // (coord[pivoti*dim + j] - coord[i*dim + j]) * (coord[pivoti*dim + j] - coord[i*dim + j]);
                }   
                rebuiltCoord[i*k + ki] = sqrt(distance);
            }
        }

    register double chebyshevSum = 0;
    //Calculate the sum of Chebyshev distance with rebuilt coordinates between every points
    //malloc temp
    for(int i=0; i<n; ++i){
        double *t;
        t = &(rebuiltCoord[i*k]);
        for(int j=i+1; j<n; ++j){
            double *t2;
            t2 = &(rebuiltCoord[j*k]);
            register double chebyshev = 0;        
            for(int ki=0; ki<k; ++ki){
                temp[ki] = fabs( *(t + ki) - *(t2 + ki));
            }
            for(int ki=0; ki<k; ++ki){
                if(temp[ki] > chebyshev){
                    chebyshev = temp[ki];
                }
            }
            chebyshevSum += chebyshev;
        }
    }
    //free(rebuiltCoord);
   _mm_free(rebuiltCoord);
   //munmap(rebuiltCoord, sizeof(double) * n * k);
    return chebyshevSum;
}


// Recursive function Combination() : combine pivots and calculate the sum of distance while combining different pivots.
// ki  : current depth of the recursion
// k   : number of pivots
// n   : number of points
// dim : dimension of metric space
// M   : number of combinations to store
// coord  : coordinates of points
// pivots : indexes of pivots
// maxDistanceSum  : the largest M distance sum
// maxDisSumPivots : the top M pivots combinations
// minDistanceSum  : the smallest M distance sum
// minDisSumPivots : the bottom M pivots combinations
void Combination(int ki, const int k, const int n, const int dim, const int M, double* coord, int* pivots,
                 double* maxDistanceSum, int* maxDisSumPivots, double* minDistanceSum, int* minDisSumPivots){
    if(ki==k-1){
        int start = pivots[ki-1]+1;
#pragma omp parallel
{        
#pragma omp for nowait schedule(dynamic)
        for(int i=start; i<n; i++){
        // pivots -> pivotp    
            int* pivotp = (int*)malloc(sizeof(int) * k);
            for(int i=0;i<k-1;i++){
                  pivotp[i] = pivots[i];
            }
            pivotp[k-1] = i;
            // Calculate sum of distance while combining different pivots.
            double distanceSum = SumDistance(k, n, dim, coord, pivotp);
#pragma omp critical
{
            // put data at the end of array
            maxDistanceSum[M] = distanceSum;
            minDistanceSum[M] = distanceSum;
            for(int kj=0; kj<k; kj++){
                maxDisSumPivots[M*k + kj] = pivotp[kj];
                minDisSumPivots[M*k + kj] = pivotp[kj];
            }
            // sort
            int a;
            for(a=M; a>0; a--){
                if(likely(maxDistanceSum[a] > maxDistanceSum[a-1])){
                    double temp = maxDistanceSum[a];
                    maxDistanceSum[a] = maxDistanceSum[a-1];
                    maxDistanceSum[a-1] = temp;
                    int kj;
                    for(kj=0; kj<k; kj++){
                        int temp = maxDisSumPivots[a*k + kj];
                        maxDisSumPivots[a*k + kj] = maxDisSumPivots[(a-1)*k + kj];
                        maxDisSumPivots[(a-1)*k + kj] = temp;
                    }
                }
                else{
                    break;
                }
            }
            for(a=M; a>0; a--){
                if(likely(minDistanceSum[a] < minDistanceSum[a-1])){
                    double temp = minDistanceSum[a];
                    minDistanceSum[a] = minDistanceSum[a-1];
                    minDistanceSum[a-1] = temp;
                    int kj;
                    for(kj=0; kj<k; kj++){
                        int temp = minDisSumPivots[a*k + kj];
                        minDisSumPivots[a*k + kj] = minDisSumPivots[(a-1)*k + kj];
                        minDisSumPivots[(a-1)*k + kj] = temp;
                    }
                }
                else{
                    break;
                }
            }
}
            free(pivotp);
        }
}
        return;
    }

    // Recursively call Combination() to combine pivots
    for(int i=pivots[ki-1]+1; i<n; i++){
        pivots[ki] = i;
        Combination(ki+1, k, n, dim, M, coord, pivots, maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots);
        /** Iteration Log : pivots computed, best pivots, max distance sum, min distance sum pivots, min distance sum
        *** You can delete the logging code. **/
    }
}

int main(int argc, char* argv[]){
    // filename : input file namespace
    char* filename = (char*)"uniformvector-2dim-5h.txt";
    if( argc==2 ) {
        filename = argv[1];
    }  else if(argc != 1) {
        printf("Usage: ./pivot <filename>\n");
        return -1;
    }
    // M : number of combinations to store
    const int M = 1000;
    // dim : dimension of metric space
    int dim;
    // n : number of points
    int n;
    // k : number of pivots
    int k;

    // Read parameter
    FILE* file = fopen(filename, "r");
    if( file == NULL ) {
        printf("%s file not found.\n", filename);
        return -1;
    }
    fscanf(file, "%d", &dim);
    fscanf(file, "%d", &n);
    fscanf(file, "%d", &k);
    printf("dim = %d, n = %d, k = %d\n", dim, n, k);

    // Start timing
    struct timeval start;

    // Read Data
    double* coord = (double*)_mm_malloc(sizeof(double)*dim*n, 32);
    // double* coord = (double*)_mm_malloc(sizeof(double) * dim * n, 32);
    int i;
    for(i=0; i<n; i++){
        int j;
        for(j=0; j<dim; j++){
            fscanf(file, "%lf", &coord[i*dim + j]);
        }
    }
    fclose(file);
    gettimeofday(&start, NULL);

    // maxDistanceSum : the largest M distance sum
   double* maxDistanceSum = (double*)malloc(sizeof(double) * (M+1));
    //double*  maxDistanceSum = (double*)_mm_malloc(sizeof(double) *(M+1), 32);
    for(i=0; i<M; i++){
        maxDistanceSum[i] = 0;
    }
    // maxDisSumPivots : the top M pivots combinations
    int* maxDisSumPivots = (int*)malloc(sizeof(int) * k * (M+1));

    // minDistanceSum : the smallest M distance sum
    double* minDistanceSum = (double*)malloc(sizeof(double) * (M+1));
    for(i=0; i<M; i++){
        minDistanceSum[i] = __DBL_MAX__;
    }
    // minDisSumPivots : the bottom M pivots combinations
    int* minDisSumPivots = (int*)malloc(sizeof(int) * k * (M+1));

    // temp : indexes of pivots with dummy array head
    int* temp = (int*)malloc(sizeof(int) * (k+1));
    temp[0] = -1;
    
    // Main loop. Combine different pivots with recursive function and evaluate them. Complexity : O( n^(k+2) )
    Combination(0, k, n, dim, M, coord, &temp[1], maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots);

    // End timing
    struct timeval end;
    gettimeofday (&end, NULL);
    printf("Using time : %f ms\n", (end.tv_sec-start.tv_sec)*1000.0+(end.tv_usec-start.tv_usec)/1000.0);

    // Store the result
    FILE* out = fopen("result.txt", "w");
    for(i=0; i<M; i++){
        int ki;
        for(ki=0; ki<k-1; ki++){
            fprintf(out, "%d ", maxDisSumPivots[i*k + ki]);
        }
        fprintf(out, "%d\n", maxDisSumPivots[i*k + k-1]);
    }
    for(i=0; i<M; i++){
        int ki;
        for(ki=0; ki<k-1; ki++){
            fprintf(out, "%d ", minDisSumPivots[i*k + ki]);
        }
        fprintf(out, "%d\n", minDisSumPivots[i*k + k-1]);
    }
    fclose(out);

    // Log
    int ki;
    printf("max : ");
    for(ki=0; ki<k; ki++){
        printf("%d ", maxDisSumPivots[ki]);
    }
    printf("%lf\n", maxDistanceSum[0]);
    printf("min : ");
    for(ki=0; ki<k; ki++){
        printf("%d ", minDisSumPivots[ki]);
    }
    printf("%lf\n", minDistanceSum[0]);
    // for(i=0; i<M; i++){
        // int ki;
        // for(ki=0; ki<k; ki++){
            // printf("%d\t", maxDisSumPivots[i*k + ki]);
        // }
        // printf("%lf\n", maxDistanceSum[i]);
    // }
    // for(i=0; i<M; i++){
        // int ki;
        // for(ki=0; ki<k; ki++){
            // printf("%d\t", minDisSumPivots[i*k + ki]);
        // }
        // printf("%lf\n", minDistanceSum[i]);
    // }
    return 0;
}

