#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <pthread.h>
#include "graphics.h"

#define NUM_THREADS 4

typedef struct particle {
    double x, y, mass, vel_x, vel_y, brightness;
} particle_t;

const double eps_0 = 1E-3;
double* particles;

typedef struct thread_args {
    int start;
    double* forces;
    particle_t* particles;
    int particles_per_thread;
    double G;
    int N;
    
} thread_args_t;

/*
Method of reading data from .gal files was inspired by the function read_doubles_from_file
in the given compare_gal_files.c
*/

int set_initial_data(int N, particle_t** particles, const char* filename) {
    int i;
    double *buffer = (double*) malloc(N*sizeof(double)*6);
    FILE* file = fopen(filename, "rb");
    if (!file) {
        free(buffer);
        return 0;
    }
    fseek(file, 0L, SEEK_END);
    size_t file_size = ftell(file);
    if (file_size != 6*N*sizeof(double)) {
        free(buffer);
        return 0;
    }
    fseek(file, 0L, SEEK_SET);

    fread(buffer, sizeof(char), file_size, file);

    fclose(file);

    for (i = 0; i < N; i++) {
        (*particles)[i].x = buffer[i*6+0];
        (*particles)[i].y = buffer[i*6+1];
        (*particles)[i].mass = buffer[i*6+2];
        (*particles)[i].vel_x = buffer[i*6+3];
        (*particles)[i].vel_y = buffer[i*6+4];
        (*particles)[i].brightness = buffer[i*6+5];
    }
    free(buffer);
    return 1;
}

void* thread_function(void *arg) {
    
    // Set input parameters
    thread_args_t* args = (thread_args_t*) arg;
    double* forces = args->forces;
    particle_t* particles = args->particles;
    const int particles_per_thread = args->particles_per_thread;
    const int start = args->start;
    const double G = args->G;
    const int N = args->N;

    // Declare variables
    double pos_x, pos_y, m;
    double r_x, r_y, r;
    double F_const, F_x, F_y;
    double denom;
    
    // Declare iteration variables and counter
    unsigned int i;
    unsigned int j;
    int c = 0; // counter
    
    // Calculate all forces acting on particle i
    for (i = start; i < (start + particles_per_thread); i++) {
        pos_x = particles[i].x;
        pos_y = particles[i].y;
        m = particles[i].mass;
        
        F_x = 0;
        F_y = 0;
        
        for (j = 0; j < N; j++) {
            if (i != j) {
            r_x = pos_x - particles[j].x;
            r_y = pos_y - particles[j].y;
            r = sqrt(r_x*r_x + r_y*r_y);
            
            denom = (r + eps_0)*(r + eps_0)*(r + eps_0);
            F_const = G * m * (particles[j].mass/denom);
            F_x += F_const * r_x;
            F_y += F_const * r_y;
            }
        }
        
        forces[c*2 + 0] = F_x;
        forces[c*2 + 1] = F_y;
        c++;
    }

    return NULL;
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
        
    // Declare constants
    double G = (double) -100.0 / N;
    
    // Store the properties of all particles in the array particles
    particle_t* particles = (particle_t*) malloc(N*sizeof(particle_t)*6);
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
    unsigned int j;
    unsigned int t;
    unsigned int i;
    unsigned int l;
    
    // Thread related declarations
    pthread_t threads[NUM_THREADS];

    if (N % NUM_THREADS != 0) {
        printf("Number of elements is not divisible by number of threads.\n");
        return 0;
    }
    
    int particles_per_thread = N / NUM_THREADS;
    thread_args_t** arguments = (thread_args_t**) malloc(NUM_THREADS*sizeof(thread_args_t*));


    for (j = 0; j < NUM_THREADS; j++) {
        arguments[j] = malloc(sizeof(thread_args_t));
        arguments[j]->start = j*particles_per_thread;
        arguments[j]->forces =  malloc(sizeof(double)*2*particles_per_thread);
        arguments[j]->particles = malloc(sizeof(particle_t));
        arguments[j]->particles_per_thread = particles_per_thread;
        arguments[j]->G = G;
        arguments[j]->N = N;
    }
    
    for (t = 0; t < nsteps; t++) {
        if (graphics!= 0) {
            ClearScreen();
            for (l = 0; l < N; l++) {
                DrawCircle(particles[l].x, particles[l].y, 1, 1, particles[l].mass*0.002, 0);
            }
            Refresh();
            usleep(2000);
        }
        
        for (j = 0; j < NUM_THREADS; j++) {
            arguments[j]->particles = particles;
            pthread_create(&threads[j], NULL, thread_function, arguments[j]);
        }
        
        for (j = 0; j < NUM_THREADS; j++) {
            pthread_join(threads[j], NULL);
                    
            for (i = 0; i < particles_per_thread; i++) {
                particles[j*particles_per_thread+i].vel_x = particles[j*particles_per_thread+i].vel_x +
                            delta_t*arguments[j]->forces[i*2+0]/particles[j*particles_per_thread+i].mass;

                particles[j*particles_per_thread+i].vel_y = particles[j*particles_per_thread+i].vel_y +
                            delta_t*arguments[j]->forces[i*2+1]/particles[j*particles_per_thread+i].mass;
                        
                particles[j*particles_per_thread+i].x = particles[j*particles_per_thread+i].x +
                            delta_t * particles[j*particles_per_thread+i].vel_x;
                
                particles[j*particles_per_thread+i].y = particles[j*particles_per_thread+i].y +
                            delta_t * particles[j*particles_per_thread+i].vel_y;
                    }
            arguments[j]->particles = NULL;
                }

    }
    
    FILE *f_ptr;

    f_ptr = fopen("results.gal", "wb");
    fwrite(particles, N*sizeof(particle_t), 1, f_ptr); // Write all data to binary file
    fclose(f_ptr);
    
    for (j = 0; j < NUM_THREADS; j++) {
        free(arguments[j]->forces);
        free(arguments[j]->particles);
    }
    
    free(arguments);
    free(particles);
    printf("Galsim program took %7.3f wall seconds.\n", get_timings() - time);
    
    return 0;
}

