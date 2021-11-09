#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>

void changeValues (int *a, int *b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

void shuffle(int *arr, int n) {
    for (int i = n-1; i > 0; i--)
    {
        int j = rand() % (i+1);
        changeValues(&arr[i], &arr[j]);
    }
}

int generate_random_code(int n, int dv, int dc, int *variable_lookup, int *check_lookup, bool *parity_check, int iterations, bool first_run, int seed) {
    if(first_run == true) {
        srand(seed);
        first_run = false;
    }
    if (iterations > 1000000) {
        return 0;
    }
    int k = n*(dc-dv)/dc;

    // Create check_lookup as array of 1,2,...,n*dv
    for(int i=0; i<n*dv; i++) {
        check_lookup[i] = i;
    }
    // Shuffle order of list to create random permutation
    shuffle(check_lookup, n*dv);
    // Floor divide by dv to get 
    for(int i=0; i<n*dv; i++) {
        check_lookup[i] /= dv;
    }

    // Check that no one check node is connected to the same variable node more than once, use same loop to start building parity_check matrix
    for(long long i=0; i<n-k; i++) {
        for(long long j=0; j<dc; j++) {
            for(int l=0; l<dc; l++) {
                if (j!=l && check_lookup[i*dc+j] == check_lookup[i*dc+l]) {
                    return generate_random_code(n,dv,dc,variable_lookup,check_lookup,parity_check,iterations+1, false, seed);
                }
            }

        }
    }

    int* indexes = (int*) malloc(n*sizeof(int));

    memset(indexes, 0, n*sizeof(int));

    for(long long i=0; i<n-k; i++) {
        for(long long j=0; j<dc; j++) {
            parity_check[i*n+check_lookup[i*dc+j]] = 1; 
        }
        for(int j=0; j<n; j++) {
            if(parity_check[i*n+j] == 1) {
                variable_lookup[j*dv+indexes[j]] = i;

                indexes[j] += 1;
            }
        }
    }
    
    return 1;
}

// int main() {
//     printf("Running main");
//     srand(time(0));
// }

// int main() {
//     printf("Starting test run new\n");
//     long long n = 100;
//     long long k = 50;
//     int* check_lookup = (int*) malloc(n*3*sizeof(int));
//     int* variable_lookup = (int*) malloc(n*3*sizeof(int));
//     printf("Assigning parity_check\n");
//     bool* parity_check = (bool*) malloc(n*k*sizeof(bool));
//     if(parity_check) {
//         printf("Assigned\n");
//     }
//     else{
//         printf("Error in allocation\n");
//     }
//     for(long long j=0; j<n; j++) {
//         for(long long l=0; l<k; l++) {
//             parity_check[j*k+l] = 0;
//         }
//     }
//     generate_random_code(n,3,6,variable_lookup, check_lookup, parity_check, 0, true);
//     printf("Done");

//     free(check_lookup);
//     free(variable_lookup);
//     free(parity_check);
//     return 0;
// }