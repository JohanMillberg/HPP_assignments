#include <stdio.h>
#include <stdlib.h>

int main() {
    int integernumber;
    double doublenumber;
    char character;
    float floatnumber;

    FILE *fptr;
    fptr = fopen("little_bin_file", "r"); // Opens the file

    //Reads each entry into the corresponding variable and prints
    fread(&integernumber, sizeof(int), 1, fptr);
    printf("The integer is %d\n", integernumber);

    fread(&doublenumber, sizeof(double), 1, fptr);
    printf("The double is %lf\n", doublenumber);

    fread(&character, sizeof(char), 1, fptr);
    printf("The char is %c\n", character);

    fread(&floatnumber, sizeof(float), 1, fptr);
    printf("The float is %f \n", floatnumber);

    fclose(fptr); // Closes the file
}