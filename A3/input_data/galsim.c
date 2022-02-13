#include <math.h>
#include <stdlib.h>
#include <stdio.h>
<<<<<<< HEAD

typedef struct particle {
    double pos_x, pos_y, vel_x, vel_y, mass, brightness;
    
} particle_t;

void set_initial_data(int N, particle_t** particle, const char* filename) {
    FILE* file = fopen(filename, "rb");
    
    fseek(file, 0L, SEEK_END);
    size_t fileSize = ftell(file);
    fseek(file, 0L, SEEK_SET);
    
=======
#include "graphics.h"

typedef struct particle {
    double pos_x, pos_y, vel_x, vel_y, mass, brightness;
} particle_t;

/*
Method of reading data from .gal files was inspired by the function read_doubles_from_file
in the given compare_gal_files.c
*/
void set_initial_data(int N, particle_t** particle, const char* filename) {
    FILE* file = fopen(filename, "rb");
    fseek(file, 0L, SEEK_END);
    size_t fileSize = ftell(file);
    fseek(file, 0L, SEEK_SET);

>>>>>>> origin/master
    double buffer[6*N];

    fread(buffer, sizeof(char), fileSize, file);

    for (int i = 0; i < N; i++) {
        (*particle)[i].pos_x = buffer[i*6 + 0];
        (*particle)[i].pos_y = buffer[i*6 + 1];
        (*particle)[i].mass = buffer[i*6 + 2];
        (*particle)[i].vel_x = buffer[i*6 + 3];
        (*particle)[i].vel_y = buffer[i*6 + 4];
        (*particle)[i].brightness = buffer[i*6 + 5];
    }
<<<<<<< HEAD
    
=======
>>>>>>> origin/master
    fclose(file);
}

void Force(int N, int i, particle_t *particle, double arr[]) {
<<<<<<< HEAD
    
    double F_const, r, r2, denom;
    double Fx = 0; double Fy = 0;
    const double G = -100/N;
    const double eps_0 = 0.001;

    for (int j = 0; j < (N-1); j++) {
        r2 = (particle[i].pos_x - particle[j].pos_x)*(particle[i].pos_x - particle[j].pos_x) + (particle[i].pos_y - particle[j].pos_y)*(particle[i].pos_y - particle[j].pos_y);
        r = sqrt(r2);
        
        denom = (r + eps_0)*(r + eps_0)*(r + eps_0);
        F_const = G*particle[i].mass * (particle[j].mass/denom);
        Fx += F_const*(particle[i].pos_x - particle[j].pos_x);
        Fy += F_const*(particle[i].pos_y - particle[j].pos_y);
    }
    
=======

    double F_const, r, r2, denom;
    double Fx = 0; double Fy = 0;
    double r_x, r_y;
    const double G = (double) -100.0/N; // Extremely important to cast to double
    const double eps_0 = 0.001;

    for (int j = 0; j < N; j++) {
        if (i != j) {
            r_x = particle[i].pos_x - particle[j].pos_x;
            r_y = particle[i].pos_y - particle[j].pos_y;
            r2 = r_x*r_x + r_y*r_y;
            r = sqrt(r2);

            denom = (r + eps_0)*(r + eps_0)*(r + eps_0);
            F_const = G*particle[i].mass * (particle[j].mass/denom);
            Fx += F_const*(r_x);
            Fy += F_const*(r_y);
        }
    }
>>>>>>> origin/master
    arr[0] = Fx;
    arr[1] = Fy;
}


int main(int argc, char *argv[]) {
<<<<<<< HEAD
  
=======

>>>>>>> origin/master
    // Set input parameters
    const int N = atoi(argv[1]);
    const char* filename = argv[2];
    const int nsteps = atoi(argv[3]);
    const double delta_t = atof(argv[4]);
<<<<<<< HEAD
    //const int graphics = atoi(argv[5]);
    
    // Set constants
    const double T = nsteps*delta_t;
    double Fx, Fy;
    
    // Declare variables
    
    // Declare time variables
    double vel_y_new, vel_x_new;
    double pos_x_new, pos_y_new;
    
    particle_t *particle = (particle_t*)malloc(N*sizeof(*particle));
    particle_t *particle_new = (particle_t*)malloc(N*sizeof(*particle));

    set_initial_data(N, &particle, filename);
    set_initial_data(N, &particle_new, filename);

    for(double t = 0; t < T; t+=delta_t) {
        for (int i = 0; i < N; i++) {
            double arr[2];
            
            Force(N, i, particle, arr);
            Fx = arr[0];
            Fy = arr[1];
            
            particle_new[i].vel_x = particle[i].vel_x + delta_t*(Fx/particle[i].mass);
            particle_new[i].vel_y = particle[i].vel_y + delta_t*(Fy/particle[i].mass);
            
            particle_new[i].pos_x = particle[i].pos_x + delta_t*particle_new[i].vel_x;
            particle_new[i].pos_y = particle[i].pos_y + delta_t*particle_new[i].vel_y;
        }
        
        particle = particle_new;
    }
    
    free(particle);
    free(particle_new);
    
    double buffer[6*N];
    
=======
    const int graphics = atoi(argv[5]);
    const int windowWidth=800;

    // Set constants
    double Fx, Fy;

    particle_t *particle = (particle_t*)malloc(N*sizeof(particle_t));
    particle_t *particle_new = (particle_t*)malloc(N*sizeof(particle_t));
    particle_t *temp;

    set_initial_data(N, &particle, filename);
    set_initial_data(N, &particle_new, filename);
    if (graphics == 1) {
        InitializeGraphics(argv[0],windowWidth,windowWidth);
        SetCAxes(0,1);
    }

    for(int t = 0; t < nsteps; t++) {
        if (graphics == 1) {
            ClearScreen();
            for (int l = 0; l < N; l++) {
                DrawCircle(particle[l].pos_x, particle[l].pos_y, 1, 1, particle[l].mass*0.001, 0);
            }
            Refresh();
            usleep(2000);
        }
        for (int i = 0; i < N; i++) {
            double arr[2];

            Force(N, i, particle, arr);
            Fx = arr[0];
            Fy = arr[1];

            particle_new[i].vel_x = particle[i].vel_x + delta_t*(Fx/particle[i].mass);
            particle_new[i].vel_y = particle[i].vel_y + delta_t*(Fy/particle[i].mass);

            particle_new[i].pos_x = particle[i].pos_x + delta_t*particle_new[i].vel_x;
            particle_new[i].pos_y = particle[i].pos_y + delta_t*particle_new[i].vel_y;

        }
        temp = particle;
        particle = particle_new;
        particle_new = temp;
        /*
        for(int k = 0; k < N; k++) {
            particle[k].vel_x = particle_new[k].vel_x;
            particle[k].vel_y = particle_new[k].vel_y;
            particle[k].pos_x = particle_new[k].pos_x;
            particle[k].pos_y = particle_new[k].pos_y;
        }
        */
    }
    if (graphics == 1) {
        FlushDisplay();
        CloseDisplay();
    }

    double buffer[6*N];

>>>>>>> origin/master
    for (int k = 0; k < N; k++) {
        buffer[6*k + 0] = particle[k].pos_x;
        buffer[6*k + 1] = particle[k].pos_y;
        buffer[6*k + 2] = particle[k].mass;
        buffer[6*k + 3] = particle[k].vel_x;
        buffer[6*k + 4] = particle[k].vel_y;
        buffer[6*k + 5] = particle[k].brightness;
    }
<<<<<<< HEAD
    
    FILE *ptr;
    
    ptr = fopen("results_ellipse.gal", "wb");
    fwrite(buffer, sizeof(buffer), 1, ptr);
    
    
    return 0;
}
=======

    FILE *ptr;

    ptr = fopen("results.gal", "wb");
    fwrite(buffer, sizeof(buffer), 1, ptr);
    fclose(ptr);

    free(particle);
    free(particle_new);

    return 0;
}
>>>>>>> origin/master
