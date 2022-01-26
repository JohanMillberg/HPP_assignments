#include <stdio.h>
#include <stdlib.h>

int element(int n, int k) {
    int sum = 1;
    if (n == 0) {
        return sum;
    }
    for (int i = 0; i < k; i ++) {
        sum *= (n-i)/(i+1);
    }
    return sum;
}

int main(int argc, char *argv[]) {
    int i, j;
    int el;
    int n = atoi(argv[1]);
    printf("Main goes here, %d \n", n);

    for (j = 0; j < n; j++) {
        for (i = 0; i < j; i++) {
            el = element(j,n);
            printf("%d", el);
        }
        printf("\n");
    }

}