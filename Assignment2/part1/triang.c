#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
    int n = atoi(argv[1]);
    int temp = 1;
    int **triangle = (int **) malloc(n*sizeof(int*)); // Allocates memory for a 2D array

    for (int i = 0; i < n; i++) {
        triangle[i] = (int *) malloc((i+1)*sizeof(int)); // Allocates memory for each row (array)

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
    for (int i = 0; i < n; i++) {
        free(triangle[i]);
        triangle[i] = NULL;
    }
    free(triangle);
    triangle = NULL;
}