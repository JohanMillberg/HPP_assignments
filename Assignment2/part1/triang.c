#include <stdio.h>
#include <stdlib.h>

/*int element(int n, int k) {
    int product = 1;
    if (k > (n - k)) {
        k = n - k;
    }
    for (int i = 1; i <= k; i ++) {
        product *= (n+1-i)/(i);
    }
    return product;
}

int main(int argc, char *argv[]) {
    int i, j;
    int el;
    int n = atoi(argv[1]);

    for (i = 0; i < n; i++) {
        for (j = 0; j <= i; j++) {
            el = element(i,j);
            printf("%d", el);
        }
        printf("\n");
    }

}*/

int main(int argc, char *argv[]) {
    int n = atoi(argv[1]);
    int temp = 1;
    int **triangle = (int **) malloc(n*sizeof(int*));

    for (int i = 0; i < n; i++) {
        triangle[i] = (int *) malloc(i*sizeof(int));

        for (int j = 0; j <= i; j++) {
            if (i == 0 || j == 0) {
                temp = 1;
                triangle[i][j] = temp;
            }
            else {
                temp = temp * (i - j + 1)/j;
                triangle[i][j] = temp;
            }
            printf("%d ", triangle[i][j]);
        }
        printf("\n");
    }
}