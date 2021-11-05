#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>

int message_passing(int *Mvc, int iterations, int *variable_to_check_list, int *check_to_variable_list, int *errors, int n, int k, int dv, int dc) {
    int it;

    int* Mcv = (int*) malloc(n*dv*sizeof(int));
    int total = 0;
    bool erasure = false;
    int check_node;

    for (it=0; it<iterations; it++) {

        // Update check->variable messages
        // For each check node
        for(int i=0; i<n-k; i++) {
            // For each variable connected to each check node
            for(int j=0; j<dc; j++) {
                // Sum all messages incoming to check node (minus variable we are updating)
                total = 0;
                erasure = false;
                for(int l=0; l<dc; l++) {
                    if (Mvc[check_to_variable_list[i*dc+l]] == 2 && l!=j) {
                        erasure = true;
                        break;
                    }
                    total += Mvc[check_to_variable_list[i*dc+l]];
                }
                total -= Mvc[check_to_variable_list[i*dc+j]];
                if (erasure==1) {
                    Mcv[i*dc+j] = 2;
                }
                else {
                    Mcv[i*dc+j] = total % 2;
                }
            }
        }

        int i;
        int j;
        int l;

        // Update variable->check messages
        for(i=0; i<n; i++) {
            // For each check connected to each variable node
            if (Mvc[i]==2) {
                for(j=0; j<dv; j++) {
                    check_node = variable_to_check_list[i*dv+j];

                    for(int l=0; l<dc; l++) {
                        if (i==check_to_variable_list[check_node*dc+l] && Mcv[check_node*dc+l] != 2) {
                            Mvc[i] = Mcv[check_node*dc+l];
                        }
                    }
                }
            }
        }


        // Check if erasure in messages
        erasure = false;
        for (int i=0; i<n; i++) {
            if (Mvc[i]==2) {
                erasure = true;
                errors[it] += 1;
            }
        }
        if (!erasure) {
            break;
        }
    }
    free(Mcv);
    return it;
}

// int test_message_passing(int *Mvc, int iterations, int *variable_to_check_list, int *check_to_variable_list, int *errors, int n, int k, int dv, int dc) {
//     int it;
//     int total = 0;
//     bool erasure = false;
//     bool valueFound = false;
//     int16_t check_node;
//     int i,j,l;

//     int8_t Mcv[n-k][n];
    
//     // int8_t Mcv[n*dv];
//     // double *Mcv;
//     // Mcv = calloc(n*dv, sizeof(int8_t));
//     printf("Allocated memory\n");

//     for (it=0; it<iterations; it++) {

//         printf("%d\n", it);

//         // Update check->variable messages
//         // For each check node
//         for(i=0; i<n-k; i++) {
//             // For each variable connected to each check node
//             for(j=0; j<dc; j++) {
//                 // Sum all messages incoming to check node (minus variable we are updating)
//                 total = 0;
//                 erasure = false;
//                 for(l=0; l<dc; l++) {
//                     if (Mvc[check_to_variable_list[i*dc+l]] == 2 && l!=j) {
//                         erasure = true;
//                         break;
//                     }
//                     total += Mvc[check_to_variable_list[i*dc+l]];
//                 }
//                 total -= Mvc[check_to_variable_list[i*dc+j]];
//                 if (erasure==1) {
//                     // Mcv[i*dc+j] = 2;
//                     Mcv[i][check_to_variable_list[i*dc+j]] = 2;
//                 }
//                 else {
//                     // Mcv[i*dc+j] = total % 2;
//                     Mcv[i][check_to_variable_list[i*dc+j]] = total % 2;
//                 }
//             }
//         }
//         printf("Done check to variable");
//         // bool valueFound = false;
//         // int16_t check_node;
//         // Update variable->check messages
//         for(i=0; i<n; i++) {
//             // For each check connected to each variable node
//             valueFound = false;
//             if (Mvc[i]==2) {
//                 for(j=0; i<dv; j++) {
//                     // check_node = variable_to_check_list[i*dv+j];
//                     // printf("%d\n", variable_to_check_list[i*dv+j]);
//                     // for(l=0; l<dc; l++) {
//                     //     if (i==check_to_variable_list[check_node+l] && Mcv[check_node*dc+l] != 2) {
//                     //         Mvc[i] = Mcv[check_node*dc+l];
//                     //     }
//                     // }
//                 }

//                 // for(int j=0; j<dv; j++) {
//                 //     // if (Mcv[variable_to_check_list[i*dv+j]*dc+]
//                 //     if (Mcv[variable_to_check_list[i*dv+j]][i] != 2) {
//                 //         Mvc[i] = Mcv[variable_to_check_list[i*dv+j]][i];
//                 //         break;
//                 //     }
//                 // }
//             }
//         }


//     //     // Check if erasure in messages
//     //     bool erasure = false;
//     //     for (int i=0; i<n; i++) {
//     //         if (Mvc[i]==2) {
//     //             erasure = true;
//     //             errors[it] += 1;
//     //         }
//     //     }
//     //     if (!erasure) {
//     //         break;
//     //     }
//     }

//     // printf("Iterations: %d\n", it);
//     // for (int i=0; i<n; i++) {
//     //     printf("%d\n", Mvc[i]);
//     // }
//     // return it;
//     return it;
// }


// int main(int argc, char *argv[]) {
//     int variable_to_check[12] = {0,1,0,1,0,2,1,2,0,2,1,2};
//     int check_to_variable[12] = {0,1,2,4,0,1,3,5,2,3,4,5};
//     int channel_output[6] = {1,1,2,0,1,0};
//     int decoded_message[6];
//     message_passing(channel_output, 1, variable_to_check, check_to_variable, 6, 3, 2, 4);
//     for (int i=0; i<6; i++) {
//         printf("%d\n", channel_output[i]);
//     }
// }