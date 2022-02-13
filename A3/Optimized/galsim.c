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
    double r_x, r_y, r, r2;
    double denom;
    const double G = (double) -100.0 / N;
    const double eps_0 = 0.001;

    double *particles = (double*) malloc(N*sizeof(double)*6);
    double *particles_new = (double*) malloc(N*sizeof(double)*6);
    double *temp;

    set_initial_data(N, &particles, filename);
    memcpy(particles_new, particles, N*sizeof(double)*6);

    if (graphics == 1) {
        InitializeGraphics(argv[0],windowWidth,windowWidth);
        SetCAxes(0,1);
    }

    for (int t = 0; t < nsteps; t++) {
        if (graphics == 1) {
            ClearScreen();
            for (int l = 0; l < N; l++) {
                DrawCircle(particles[l*6+0], particles[l*6+1], 1, 1, particles[l*6+2]*0.001, 0);
            }
            Refresh();
            usleep(2000);
        }
        for (int i = 0; i < N; i++) {
            F_x = 0;
            F_y = 0;
            for (int j = 0; j < N; j++) {
                if (i != j) {
                    r_x = particles[i*6 + 0] - particles[j*6 + 0];
                    r_y = particles[i*6 + 1] - particles[j*6 + 1];
                    r2 = r_x*r_x + r_y*r_y;
                    r = sqrt(r2);

                    denom = (r + eps_0)*(r + eps_0)*(r + eps_0);
                    F_const = G * particles[i*6 + 2] * (particles[j*6 + 2]/denom);
                    F_x += F_const * r_x;
                    F_y += F_const * r_y;
                }
            }
            particles_new[i*6 + 3] = particles[i*6 + 3] + delta_t*(F_x/particles[i*6 + 2]);
            particles_new[i*6 + 4] = particles[i*6 + 4] + delta_t*(F_y/particles[i*6 + 2]);

            particles_new[i*6 + 0] = particles[i*6 + 0] + delta_t*particles_new[i*6 + 3];
            particles_new[i*6 + 1] = particles[i*6 + 1] + delta_t*particles_new[i*6 + 4];
        }
        temp = particles;
        particles = particles_new;
        particles_new = temp;
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
    free(particles_new);
    printf("Galsim program took %7.3f wall seconds.\n", get_timings() - time);
    return 0;

}