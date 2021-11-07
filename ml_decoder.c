#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>

bool mod2add(bool var1, bool var2) {
    if (!var1 && !var2) {
        return 0;
    }
    else {
        return !(var1 && var2);
    }
}

void ml_decode(int *channel_output, bool *parity_check, int n, int dv, int dc) {
    int k = n*(dc-dv)/dc;

    bool* remaining_parity_checks = (bool*) malloc(n*(n-k)*sizeof(bool));
    bool *target = (bool*)malloc((n-k)*sizeof(bool));

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
                target[j] = mod2add(target[j], channel_output[i]&&parity_check[j*n+i]);
            }
        }
    }
    printf("\n");
    for(int i=0; i<3; i++) {
        printf("%d\n", target[i]);
    }
    printf("\n");
    for(int i=0; i<6; i++) {
        printf("%d\n", remaining_parity_checks[i]);
    }
    free(remaining_parity_checks);
    free(target);
}

int main() {
    int channel_output[6] = {0,2,1,2,0,1};
    bool parity_check[18] = {1,1,0,1,0,1,0,1,1,1,1,0,1,0,1,0,1,1};
    ml_decode(channel_output, parity_check, 6, 2, 4);
}

