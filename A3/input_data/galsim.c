#include <math.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct particle {
    double pos_x, pos_y, vel_x, vel_y, mass, brightness;
    
} particle_t;

void set_initial_data(int N, particle_t** particle, const char* filename) {
    FILE* file = fopen(filename, "rb");
    
    fseek(file, 0L, SEEK_END);
    size_t fileSize = ftell(file);
    fseek(file, 0L, SEEK_SET);
    
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
    
    fclose(file);
}

void Force(int N, int i, particle_t *particle, double arr[]) {
    
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
    
    arr[0] = Fx;
    arr[1] = Fy;
}


int main(int argc, char *argv[]) {
  
    // Set input parameters
    const int N = atoi(argv[1]);
    const char* filename = argv[2];
    const int nsteps = atoi(argv[3]);
    const double delta_t = atof(argv[4]);
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
    
    for (int k = 0; k < N; k++) {
        buffer[6*k + 0] = particle[k].pos_x;
        buffer[6*k + 1] = particle[k].pos_y;
        buffer[6*k + 2] = particle[k].mass;
        buffer[6*k + 3] = particle[k].vel_x;
        buffer[6*k + 4] = particle[k].vel_y;
        buffer[6*k + 5] = particle[k].brightness;
    }
    
    FILE *ptr;
    
    ptr = fopen("results_ellipse.gal", "wb");
    fwrite(buffer, sizeof(buffer), 1, ptr);
    
    
    return 0;
}
