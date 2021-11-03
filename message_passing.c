#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

int message_passing(int8_t *Mvc, int iterations, int *variable_to_check_list, int *check_to_variable_list, int *errors, int n, int k, int dv, int dc) {
    int it;

    int Mcv[n-k][n];

    for (it=0; it<iterations; it++) {

        // Update check->variable messages
        // For each check node
        for(int i=0; i<n-k; i++) {
            // For each variable connected to each check node
            for(int j=0; j<dc; j++) {
                // Sum all messages incoming to check node (minus variable we are updating)
                int total = 0;
                bool erasure = false;
                for(int k=0; k<dc; k++) {
                    if (Mvc[check_to_variable_list[i*dc+k]] == 2 && k!=j) {
                        erasure = true;
                        break;
                    }
                    total += Mvc[check_to_variable_list[i*dc+k]];
                }
                total -= Mvc[check_to_variable_list[i*dc+j]];
                if (erasure==1) {
                    Mcv[i][check_to_variable_list[i*dc+j]] = 2;
                }
                else {
                    Mcv[i][check_to_variable_list[i*dc+j]] = total % 2;
                }
            }
        }

        // Updated variable->check messages
        for(int i=0; i<n; i++) {
            // For each check connected to each variable node
            if (Mvc[i]==2) {
                for(int j=0; j<dv; j++) {
                    if (Mcv[variable_to_check_list[i*dv+j]][i] != 2) {
                        Mvc[i] = Mcv[variable_to_check_list[i*dv+j]][i];
                        break;
                    }
                }
            }
        }

        // Check if erasure in messages
        bool erasure = false;
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

    // printf("Iterations: %d\n", it);
    // for (int i=0; i<n; i++) {
    //     printf("%d\n", Mvc[i]);
    // }
    return it;
}


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