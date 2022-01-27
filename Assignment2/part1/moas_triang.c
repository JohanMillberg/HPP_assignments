/* HPP - Assignment 2: Programming in C
    Part 1 - Pascal's Triangle
 */

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    int n = atoi(argv[1]); // index of Pascal's triangle
    int c = 1;
    
    // Allocate memory for 2D matrix
    int **arr = (int **)malloc(n * sizeof(int*));
    for (int i = 0; i < n; i++) {
        arr[i] = (int *)malloc((i+1) * sizeof(int));
    }
    
    // Assign binomial coefficients to matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            if (j == 0 || i == 0) {
                c = 1;
                arr[i][j] = c;
            }
            else {
                c = c * (i - j + 1)/j;
                arr[i][j] = c;
            }
        }
    }
    
    // Display results
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i >= j) {
                printf("%d ", arr[i][j]);
            }
            
        }
        printf("\n");
    }
    
    // Free the allocated memory
    for (int j = 0; j < n; j++) {
        free(arr[j]);
        arr[j] = NULL;
    }
    free(arr);
    arr = NULL;
    
    return 0;
}
