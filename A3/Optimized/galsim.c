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
int set_initial_data(int N, double** particles, const char* filename) {
    FILE* file = fopen(filename, "rb");
    if (!file) return 0;
    fseek(file, 0L, SEEK_END);
    size_t file_size = ftell(file);
    if (file_size != 6*N*sizeof(double)) return 0;
    fseek(file, 0L, SEEK_SET);

    fread(*particles, sizeof(char), file_size, file);

    fclose(file);
    return 1;
}

/*
get_timings() was inspired by the function get_wall_seconds() from Task 4 in Lab 5
*/
double get_timings() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double sec = tv.tv_sec + (double)tv.tv_usec / 1000000;
    return sec;
}

int main(int argc, char *argv[]) {
    if (argc < 6) {
        printf("Use following syntax to run program: ./galsim N file_name amount_steps step_length graphics_on \n");
        return 0;
    }
    double time = get_timings();

    // Set input parameters
    const int N = atoi(argv[1]);
    const char* filename = argv[2];
    const int nsteps = atoi(argv[3]);
    const double delta_t = atof(argv[4]);
    const int graphics = atoi(argv[5]);
    const int windowWidth=800;
    int successful;

    // Declare variables
    double F_x, F_y, F_const;
    double r_x, r_y, r;
    double x, y, mass;
    double denom;

    // Declare constants
    const double G = (double) -100.0 / N;
    const double eps_0 = 0.001;

    // Store the properties of all particles in the array particles
    double *particles = (double*) malloc(N*sizeof(double)*6);

    // The sum of the forces in the x and y direction for each particle is stored in forces
    double forces[2*N];

    successful = set_initial_data(N, &particles, filename);

    if (!successful) {
        printf("Error reading initial data file. \n");
        return 0;
    }

    if (graphics != 0) {
        InitializeGraphics(argv[0],windowWidth,windowWidth);
        SetCAxes(0,1);
    }
    // Declaring iteration variables
    unsigned int t;
    unsigned int l;
    unsigned int i;
    unsigned int j;
    unsigned int bl;

    // Amount of blocks and blocksize
    const int blocksize = 50;
    const int nBlocks = (2*N)/blocksize;
    int l_start;
    for (t = 0; t < nsteps; t++) {
        /*
            Draws all particles if graphics are enabled
        */
        if (graphics != 0) {
            ClearScreen();
            for (l = 0; l < N; l++) {
                DrawCircle(particles[l*6+0], particles[l*6+1], 1, 1, particles[l*6+2]*0.002, 0);
            }
            Refresh();
            usleep(2000);
        }

        //Utilizes cache blocking if 2*N is divisible by the block size
        if ((2*N) % blocksize == 0) {
            for (bl = 0; bl < nBlocks; bl ++) {
                l_start = bl*blocksize;
                for (l = l_start; l < (l_start + blocksize); l++) {
                    forces[l] = 0;
                }
            }
        }
        for (l = 0; l < 2*N; l++) {
            forces[l] = 0;
        }

        for (i = 0; i < N; i++) {
            x = particles[i*6];
            y = particles[i*6 + 1];
            mass = particles[i*6 + 2];

            /*
                Calculates all forces acting on particle i
            */
            for (j = i; j < N; j++) {
                r_x = x - particles[j*6 + 0];
                r_y = y - particles[j*6 + 1];
                r = sqrt(r_x*r_x + r_y*r_y);

                denom = (r + eps_0)*(r + eps_0)*(r + eps_0);
                F_const = G * mass * (particles[j*6 + 2]/denom);
                F_x = F_const * r_x;
                F_y = F_const * r_y;

                /*
                    Utilize the fact that the forces between the two particles are equal
                    but acting in the opposite direction, minimizing the needed iterations
                */
                forces[i*2 + 0] += F_x;
                forces[i*2 + 1] += F_y;
                forces[j*2 + 0] += -F_x;
                forces[j*2 + 1] += -F_y;

            }
            // Update the properties of the particle i
            particles[i*6 + 3] = particles[i*6 + 3] + delta_t*(forces[i*2+0]/particles[i*6 + 2]);
            particles[i*6 + 4] = particles[i*6 + 4] + delta_t*(forces[i*2+1]/particles[i*6 + 2]);

            particles[i*6 + 0] = particles[i*6 + 0] + delta_t*particles[i*6 + 3];
            particles[i*6 + 1] = particles[i*6 + 1] + delta_t*particles[i*6 + 4];

        }
    }

    if (graphics != 0) {
        FlushDisplay();
        CloseDisplay();
    }

    FILE *ptr;

    ptr = fopen("results.gal", "wb");
    fwrite(particles, N*sizeof(double)*6, 1, ptr); // Write all data to binary file
    fclose(ptr);

    free(particles);
    printf("Galsim program took %7.3f wall seconds.\n", get_timings() - time);
    return 0;

}