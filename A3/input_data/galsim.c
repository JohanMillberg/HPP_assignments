#include <stdio.h>
#include <stdlib.h>

typedef struct particle {
    double pos_x, pos_y;
    double mass;
    double vel_x, vel_y;
    double brightness;
} particle_t;

/*
Method of reading data from .gal files was inspired by the function read_doubles_from_file
in the given compare_gal_files.c
*/
void set_initial_data(int N, char* filename, particle_t** particles) {
    int amount_lines = 6*N;
    double b[amount_lines];
    FILE* fptr = fopen(filename, "rb");

    fseek(fptr, 0L, SEEK_END);
    size_t fileSize = ftell(fptr);
    fseek(fptr, 0L, SEEK_SET);

    fread(b, sizeof(char), fileSize, fptr);

    fclose(fptr);

    for (int i = 0; i < N; i++) {

        (*particles)[i].pos_x = b[i*6+0];
        (*particles)[i].pos_y = b[i*6+1];
        (*particles)[i].mass = b[i*6+2];
        (*particles)[i].vel_x = b[i*6+3];
        (*particles)[i].vel_y = b[i*6+4];
        (*particles)[i].brightness = b[i*6+5];
    }
}

int main(int argc, char *argv[]) {
    // Implement 1 loop for time, one loop for each particle
    const int amount_steps = atoi(argv[3]);
    const double dt = atof(argv[4]);
    const double T = amount_steps*dt;

    const int N = atoi(argv[1]);
    char* filename = argv[2];


    particle_t* particles = (particle_t*) malloc(N* sizeof(*particles));

    set_initial_data(N, filename, &particles);
    for (int i = 0; i < N; i++) {
        printf("%lf\n", particles[i].pos_x);
    }

    free(particles);
}