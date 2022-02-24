#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include "graphics.h"


/*
Method of reading data from .gal files was inspired by the function read_doubles_from_file
in the given compare_gal_files.c
*/
void set_initial_data(int N, double** particles, const char* filename) {
    FILE* file = fopen(filename, "rb");
    fseek(file, 0L, SEEK_END);
    size_t fileSize = ftell(file);
    fseek(file, 0L, SEEK_SET);

    fread(*particles, sizeof(char), fileSize, file);

    fclose(file);
}

double get_timings() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
    return seconds;
}

int main(int argc, char *argv[]) {
    double time = get_timings();
    // Set input parameters
    const int N = atoi(argv[1]);
    const char* filename = argv[2];
    const int nsteps = atoi(argv[3]);
    const double delta_t = atof(argv[4]);
    const int graphics = atoi(argv[5]);
    const int windowWidth=800;

    // Declare variables
    double F_x, F_y, F_const;
    double r_x, r_y, r;
    double x, y, mass;
    double denom;
    const double G = (double) -100.0 / N;
    const double eps_0 = 0.001;

    double *particles = (double*) malloc(N*sizeof(double)*6);

    double *forces = (double*) malloc(2*N*sizeof(double));

    set_initial_data(N, &particles, filename);

    if (graphics == 1) {
        InitializeGraphics(argv[0],windowWidth,windowWidth);
        SetCAxes(0,1);
    }
    int t;
    int l;
    int i;
    int j;
    for (t = 0; t < nsteps; t++) {
        if (graphics == 1) {
            ClearScreen();
            for (int l = 0; l < N; l++) {
                DrawCircle(particles[l*6+0], particles[l*6+1], 1, 1, particles[l*6+2]*0.001, 0);
            }
            Refresh();
            usleep(2000);
        }
        for (l = 0; l < 2*N; l++) {
            forces[l] = 0;
        }

        for (i = 0; i < N; i++) {
            x = particles[i*6];
            y = particles[i*6 + 1];
            mass = particles[i*6 + 2];

            for (j = i; j < N; j++) {
                r_x = x - particles[j*6 + 0];
                r_y = y - particles[j*6 + 1];
                r = sqrt(r_x*r_x + r_y*r_y);

                denom = (r + eps_0)*(r + eps_0)*(r + eps_0);
                F_const = G * mass * (particles[j*6 + 2]/denom);
                F_x = F_const * r_x;
                F_y = F_const * r_y;

                forces[i*2 + 0] += F_x;
                forces[i*2 + 1] += F_y;
                forces[j*2 + 0] += -F_x;
                forces[j*2 + 1] += -F_y;

            }
            particles[i*6 + 3] = particles[i*6 + 3] + delta_t*(forces[i*2+0]/particles[i*6 + 2]);
            particles[i*6 + 4] = particles[i*6 + 4] + delta_t*(forces[i*2+1]/particles[i*6 + 2]);

            particles[i*6 + 0] = particles[i*6 + 0] + delta_t*particles[i*6 + 3];
            particles[i*6 + 1] = particles[i*6 + 1] + delta_t*particles[i*6 + 4];

        }
    }

    if (graphics == 1) {
        FlushDisplay();
        CloseDisplay();
    }

    FILE *ptr;

    ptr = fopen("results.gal", "wb");
    fwrite(particles, N*sizeof(double)*6, 1, ptr);
    fclose(ptr);

    free(particles);
    free(forces);
    printf("Galsim program took %7.3f wall seconds.\n", get_timings() - time);
    return 0;

}