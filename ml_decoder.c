#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>

void ml_decode(int *channel_output, bool *target, bool *parity_check, bool *remaining_parity_checks, int n, int dv, int dc) {
    int k = n*(dc-dv)/dc;

    int unknown_bit_index = 0;
    for(int i=0; i<n; i++) {
        if (channel_output[i] == 2) {
            for(int j=0; j<n-k; j++) {
                remaining_parity_checks[unknown_bit_index] = parity_check[j*n+i];
                unknown_bit_index += 1;
            }
        }
        else {
            // If channel output at current bit is not erasure, calculate target vector
            for(int j=0; j<n-k; j++) {
                target[j] = target[j] ^ (channel_output[i]&&parity_check[j*n+i]);
            }
        }
    }
    // printf("No of erasures: %d\n", unknown_bit_index/(n-k));

    // printf("\n");
    // for(int i=0; i<3; i++) {
    //     printf("%d\n", target[i]);
    // }
    // printf("\n");
    // for(int i=0; i<6; i++) {
    //     printf("%d\n", remaining_parity_checks[i]);
    // }
    return;
}

// int main() {
//     int channel_output[6] = {0,2,1,2,0,1};
//     bool target[4] = {0,0,0,0};
//     bool parity_check[18] = {1,1,0,1,0,1,0,1,1,1,1,0,1,0,1,0,1,1};
//     ml_decode(channel_output, target, parity_check, 6, 2, 4);
// }

