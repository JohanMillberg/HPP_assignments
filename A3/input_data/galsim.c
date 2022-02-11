#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
    const int N = atoi(argv[1]);
    char* filename = argv[2];

    const int amount_steps = atoi(argv[3]);
    const double dt = atof(argv[4]);
    const double T = amount_steps*dt;

    const double G = 6.67408E-11;

    particle_t* particles = (particle_t*) malloc(N* sizeof(*particles));

    set_initial_data(N, filename, &particles);

    double mass_inverse;
    double gm;
    double F_common;
    double F_x;
    double F_y;
    double a_x;
    double a_y;
    double r;
    double r2;
    double epsilon = 10E-3;
    for (int t = 0; t < T; t+=dt) {
        for (int i = 0; i < N; i++) {
            // Update the attributes of each particle here, seperate function would mean many calls
            // Update each attribute here according to equations given. Calculating for each particle
            // is a parallell process.
            gm = G*particles[i].mass;
            F_x = 0;
            F_y = 0;
            mass_inverse = 1/particles[i].mass;
            for (int j = 0; j < N-1; j++) {
                r2 = (particles[i].pos_x-particles[j].pos_x)*(particles[i].pos_x-particles[j].pos_x)
                 + (particles[i].pos_y-particles[j].pos_y)*(particles[i].pos_y-particles[j].pos_y);
                r = sqrt(r2);

                F_common = gm * (particles[j].mass/((r+epsilon)*(r+epsilon)*(r+epsilon)));
                F_x += F_common * (particles[i].pos_x - particles[j].pos_x);
                F_y += F_common * (particles[i].pos_y - particles[j].pos_y);

                a_x = F_x*mass_inverse;
                a_y = F_y*mass_inverse;

                particles[i].vel_x = particles[i].vel_x + dt*a_x;
                particles[i].vel_y = particles[i].vel_y + dt*a_y;

                particles[i].pos_x = particles[i].pos_x + dt*particles[i].vel_x;
                particles[i].pos_y = particles[i].pos_y + dt*particles[i].vel_y;
             }
        }
    }

    for (int i = 0; i < N; i++) {
        printf("%lf\n", particles[i].pos_x);
    }

    free(particles);
}