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

int generate_random_code(int n, int dv, int dc, int *variable_lookup, int *check_lookup, int *sequence, bool *parity_check, int iterations, bool first_run, int seed) {
    if(first_run == true) {
        srand(time(NULL));
        first_run = false;
    }
    if (iterations > 10000) {
        return 0;
    }
    int k = n*(dc-dv)/dc;

    // Shuffle order of list to create random permutation
    shuffle(sequence, n*dv);
    // Floor divide by dv to get 
    for(int i=0; i<n*dv; i++) {
        check_lookup[i] = sequence[i]/dv;
    }

    // Check that no one check node is connected to the same variable node more than once, use same loop to start building parity_check matrix
    for(long long i=0; i<n-k; i++) {
        for(long long j=0; j<dc; j++) {
            for(int l=0; l<dc; l++) {
                if (j!=l && check_lookup[i*dc+j] == check_lookup[i*dc+l]) {
                    return generate_random_code(n,dv,dc,variable_lookup,check_lookup,sequence,parity_check,iterations+1, false, seed);
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