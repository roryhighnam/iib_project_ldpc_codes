#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

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

int generate_random_code(int n, int dv, int dc, int *variable_lookup, int *check_lookup, int *parity_check, int iterations) {

    if (iterations > 1000000) {
        return 0;
    }

    int k = n*(dc-dv)/dc;

    // Create random socket connections, and create variable_lookup and check_lookup
    // int* variable_lookup = (int*) malloc(n*dv*sizeof(int));
    // int* check_lookup = (int*) malloc(n*dv*sizeof(int));

    for(int i=0; i<n*dv; i++) {
        check_lookup[i] = i;
    }
    shuffle(check_lookup, n*dv);
    for(int i=0; i<n*dv; i++) {
        variable_lookup[i] = check_lookup[i] / dc;
    }
    for(int i=0; i<n*dv; i++) {
        check_lookup[i] /= dv;
    }
    // printf("Variable lookup:\n");
    // for(int i=0; i<n*dv; i++) {
    //     printf("%d\n", variable_lookup[i]);
    // }
    // printf("Check lookup:\n");
    // for(int i=0; i<n*dv; i++) {
    //     printf("%d\n", check_lookup[i]);
    // }

    

    // Check that no one check node is connected to the same variable node more than once
    for(int i=0; i<n-k; i++) {
        for(int j=0; j<dc; j++) {
            for(int l=0; l<dc; l++) {
                if (check_lookup[i*dc+j] == check_lookup[i*dc+l] && j!=l) {
                    return generate_random_code(n,dv,dc,variable_lookup,check_lookup,parity_check,iterations+1);
                }
            }
        }
    }
    // Check that no one check node is connected to the same variable node more than once
    // for(int i=0; i<n; i++) {
    //     for(int j=0; j<dv; j++) {
    //         for(int l=0; l<dv; l++) {
    //             if (variable_lookup[i*dv+j] == variable_lookup[i*dv+l] && j!=l) {
    //                 // printf("Duplicate");
    //                 return generate_random_code(n,dv,dc,variable_lookup,check_lookup,parity_check,iterations+1);
    //             }
    //         }
    //     }
    // }

    for(int i=0; i<n-k; i++) {
        for(int j=0; j<dc; j++) {
            parity_check[i*n+check_lookup[i*dc+j]] = 1; 
        }
    }
    return 1;
}

// int main() {
//     printf("Starting test run\n");
//     // int variable_lookup[3*10000];
//     // int check_lookup[3*10000];
//     int* check_lookup = (int*) malloc(10000*3*sizeof(int));
//     int* variable_lookup = (int*) malloc(10000*3*sizeof(int));
//     int* parity_check = (int*) malloc(10000*5000*sizeof(int));
//     for (int i=0; i<1000; i++) {
//         for(int j=0; j<10000; j++) {
//             for(int l=0; l<5000; l++) {
//                 parity_check[j*5000+l] = 0;
//             }
//         }
//         generate_random_code(10000,3,6,variable_lookup, check_lookup, parity_check, 0);
//         printf("i:%d\n",i); 
//     }
//     printf("Ending test run\n");
//     free(check_lookup);
//     free(variable_lookup);
//     free(parity_check);
// }