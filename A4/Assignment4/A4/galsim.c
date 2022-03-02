#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <pthread.h>
#include "graphics.h"

#define NUM_THREADS 4

const double eps_0 = 1E-3;

typedef struct particle {
    double x, y, mass, vel_x, vel_y, brightness;
} particle_t;

typedef struct tree_node {
    double q_x, q_y, width, height;
    struct tree_node* topLeft;
    struct tree_node* topRight;
    struct tree_node* botLeft;
    struct tree_node* botRight;
    int amount_particles;
    particle_t* particles;
    double m;
    double m_x;
    double m_y;
} tree_node_t;

typedef struct thread_args {
    int start;
    //int end;
    double* forces;
    tree_node_t* root;
    int amount_of_particles;
    double theta_max;
    double G;
} thread_args_t;

tree_node_t* create_node(double x, double y, double width, double height, int N);
tree_node_t* make_node( tree_node_t* node_parent, int counter, int position_index);
void make_tree(tree_node_t* node, int graphics);
void delete_tree(tree_node_t* node);

tree_node_t* create_node(double x, double y, double width, double height, int N) {
    tree_node_t* node = (tree_node_t*) malloc(sizeof(tree_node_t));
    node->q_x = x; node->q_y = y; node->width = width; node->height = height;
    node->amount_particles = N;
    node->particles =(particle_t *) malloc(N*sizeof(particle_t));
    node->topLeft = NULL;
    node->topRight = NULL;
    node->botLeft = NULL;
    node->botRight = NULL;
    node->m = 0;
    node->m_x = 0;
    node->m_y = 0;
    return node;
}

tree_node_t* make_node(tree_node_t* node_parent, int counter, int position_index) {
    int i;
    int j = 0;
    int x_bool, y_bool;
    double x_parent, y_parent;
    double new_width = (node_parent->width)/2.0;
    double new_height = (node_parent->height)/2.0;
    double new_x, new_y;
    x_parent = node_parent->q_x;
    y_parent = node_parent->q_y;
    tree_node_t* node = NULL;

    switch (position_index) {
        case 0: {
            new_x = x_parent;
            new_y = y_parent;

            node = create_node(new_x, new_y, new_width, new_height, counter);
            for (i = 0; i < node_parent->amount_particles; i++) {
                x_bool = (node_parent->particles[i].x >= x_parent &&
                 node_parent->particles[i].x < x_parent + new_width);

                y_bool = (node_parent->particles[i].y >= y_parent &&
                 node_parent->particles[i].y < y_parent + new_height);
                if (x_bool && y_bool) {

                    node->particles[j].x = node_parent->particles[i].x;
                    node->particles[j].y = node_parent->particles[i].y;
                    node->particles[j].mass = node_parent->particles[i].mass;
                    node->m += node->particles[j].mass;
                    node->m_x += node->particles[j].x*node->particles[j].mass;
                    node->m_y += node->particles[j].y*node->particles[j].mass;
                    j++;
                }
            }
            node->m_x = node->m_x/node->m;
            node->m_y = node->m_y/node->m;
            break;
        }
        case 1: {
            new_x = x_parent + new_width;
            new_y = y_parent;

            node = create_node(new_x, new_y, new_width, new_height, counter);
            for (i = 0; i < node_parent->amount_particles; i++) {
                x_bool = (node_parent->particles[i].x >= x_parent + new_width &&
                 node_parent->particles[i].x <= x_parent + node_parent->width);

                y_bool = (node_parent->particles[i].y >= y_parent &&
                 node_parent->particles[i].y < y_parent + new_height);

                if (x_bool && y_bool) {
                    node->particles[j].x = node_parent->particles[i].x;
                    node->particles[j].y = node_parent->particles[i].y;
                    node->particles[j].mass = node_parent->particles[i].mass;
                    node->m += node->particles[j].mass;
                    node->m_x += node->particles[j].x*node->particles[j].mass;
                    node->m_y += node->particles[j].y*node->particles[j].mass;
                    j++;
                }
            }
            node->m_x = node->m_x/node->m;
            node->m_y = node->m_y/node->m;
            break;
        }
        case 2: {
            new_x = x_parent;
            new_y = y_parent + new_height;

            node = create_node(new_x, new_y, new_width, new_height, counter);
            for (i = 0; i < node_parent->amount_particles; i++) {
                x_bool = (node_parent->particles[i].x >= x_parent &&
                 node_parent->particles[i].x < x_parent + new_width);

                y_bool = (node_parent->particles[i].y >= y_parent + new_height &&
                 node_parent->particles[i].y <= y_parent + node_parent->height);
                if (x_bool && y_bool) {
                    node->particles[j].x = node_parent->particles[i].x;
                    node->particles[j].y = node_parent->particles[i].y;
                    node->particles[j].mass = node_parent->particles[i].mass;
                    node->m += node->particles[j].mass;
                    node->m_x += node->particles[j].x*node->particles[j].mass;
                    node->m_y += node->particles[j].y*node->particles[j].mass;
                    j++;
                }
            }
            node->m_x = node->m_x/node->m;
            node->m_y = node->m_y/node->m;
            break;
        }
        case 3: {
            new_x = x_parent + new_width;
            new_y = y_parent + new_height;

            node = create_node(new_x, new_y, new_width, new_height, counter);
            for (i = 0; i < node_parent->amount_particles; i++) {
                x_bool = (node_parent->particles[i].x >= x_parent + new_width &&
                 node_parent->particles[i].x <= x_parent + node_parent->width);
                y_bool = (node_parent->particles[i].y >= y_parent + new_height &&
                 node_parent->particles[i].y <= y_parent + node_parent->height);
                if (x_bool && y_bool) {
                    node->particles[j].x = node_parent->particles[i].x;
                    node->particles[j].y = node_parent->particles[i].y;
                    node->particles[j].mass = node_parent->particles[i].mass;
                    node->m += node->particles[j].mass;
                    node->m_x += node->particles[j].x*node->particles[j].mass;
                    node->m_y += node->particles[j].y*node->particles[j].mass;
                    j++;
                }
            }
            node->m_x = node->m_x/node->m;
            node->m_y = node->m_y/node->m;
            break;
    }

    }
    return node;
}

void make_tree(tree_node_t* node, int graphics) {
    int i;
    int bool_NW[2], bool_NE[2], bool_SW[2], bool_SE[2];
    int counter[4];
    int n = node->amount_particles;

    for (i = 0; i < 4; i++) {
        counter[i] = 0;
    }

    if (n > 1) {
            // Check if the particles are in this quadrant:
        for (i = 0; i < n; i++) {
            bool_NW[0] = (node->particles[i].x >= node->q_x &&
                node->particles[i].x < node->q_x + node->width/2.0);

            bool_NW[1] = (node->particles[i].y >= node->q_y &&
                node->particles[i].y < node->q_y + node->height/2.0);

            if (bool_NW[0] && bool_NW[1]) {
                counter[0]++;
            }

            bool_NE[0] = (node->particles[i].x >= node->q_x + node->width/2.0 &&
                node->particles[i].x <= node->q_x + node->width);

            bool_NE[1] = (node->particles[i].y >= node->q_y &&
                node->particles[i].y < node->q_y + node->height/2.0);
            if (bool_NE[0] && bool_NE[1]) {
                counter[1]++;
            }

            bool_SW[0] = (node->particles[i].x >= node->q_x &&
                node->particles[i].x < node->q_x + node->width/2.0);

            bool_SW[1] = (node->particles[i].y >= node->q_y + node->height/2.0 &&
                node->particles[i].y <= node->q_y + node->height);
            if (bool_SW[0] && bool_SW[1]) {
                counter[2]++;
            }

            bool_SE[0] = (node->particles[i].x >= node->q_x + node->width/2.0 &&
                node->particles[i].x <= node->q_x + node->width);

            bool_SE[1] = (node->particles[i].y >= node->q_y + node->height/2.0 &&
                node->particles[i].y <= node->q_y + node->height);
            if (bool_SE[0] && bool_SE[1]) {
                counter[3]++;
            }
        }
        if (counter[0] > 0) {
            node->topLeft = make_node(node, counter[0], 0);
            make_tree(node->topLeft, graphics);
        }

        if (counter[1] > 0) {
            node->topRight = make_node(node, counter[1], 1);
            make_tree(node->topRight, graphics);
        }

        if (counter[2] > 0) {
            node->botLeft = make_node(node, counter[2], 2);
            make_tree(node->botLeft, graphics);
        }

        if (counter[3] > 0) {
            node->botRight = make_node(node, counter[3], 3);
            make_tree(node->botRight, graphics);
        }
    }
    return;
}

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

void traverse_tree(tree_node_t* node, particle_t particle,
        double theta_max, double* forces_x, double*forces_y, double G) {
    double theta, F_const;
    double distance_x, distance_y, r;
    double denom;
    int leaf_nodes_exist;

    if (node == NULL) {
        return;
    }

    distance_x = (particle.x-node->m_x);
    distance_y = (particle.y-node->m_y);
    r = sqrt(distance_x*distance_x + distance_y*distance_y);
    if (r == 0) {
        return;
    }
    theta = node->width/r;

    leaf_nodes_exist = (node->topLeft != NULL || node->topRight != NULL || node->botLeft != NULL ||
        node->botRight != NULL);

    if (theta > theta_max && leaf_nodes_exist) {
        traverse_tree(node->topLeft, particle, theta_max, forces_x, forces_y, G);
        traverse_tree(node->topRight, particle, theta_max, forces_x, forces_y, G);
        traverse_tree(node->botLeft, particle, theta_max, forces_x, forces_y, G);
        traverse_tree(node->botRight, particle, theta_max, forces_x, forces_y, G);
    }

    else {
        //printf("Calculating force...\n");
        //Calculate x and y forces here
        denom = (r + eps_0)*(r + eps_0)*(r + eps_0);
        F_const = G * particle.mass * (node->m/denom);
        (*forces_x) += F_const * distance_x;
        (*forces_y) += F_const * distance_y;
    }
}

void delete_tree(tree_node_t* node) {
    if (node == NULL) {
        return;
    }

    delete_tree(node->topLeft);
    delete_tree(node->topRight);
    delete_tree(node->botLeft);
    delete_tree(node->botRight);

    free(node->particles);
    free(node);
}

void* thread_function(void* arg) {
    thread_args_t* args = (thread_args_t*) arg;
    tree_node_t* root = args->root;
    double* force = args->forces;
    int start = args->start;
    //int end = args->end;
    int amount = args->amount_of_particles;
    int i;
    double G = args->G;
    double theta_max = args->theta_max;
    double force_x, force_y;

    for (i = 0; i < amount; i++) {
        force_x = 0;
        force_y = 0;
        traverse_tree(root, root->particles[i+start], theta_max, &force_x, &force_y, G);
        force[i*2+0] = force_x;
        force[i*2+1] = force_y;
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
    const double theta_max = atof(argv[5]);
    const int graphics = atoi(argv[6]);
    const int windowWidth=800;
    int successful;

    // Declare constants
    const double G = (double) -100.0 / N;

    // Store the properties of all particles in the array particles
    particle_t *particles = (particle_t*) malloc(N*sizeof(particle_t));

    successful = set_initial_data(N, &particles, filename);

    if (!successful) {
        printf("Error reading initial data file. \n");
        return 0;
    }

    tree_node_t* root;

    if (graphics != 0) {
        InitializeGraphics(argv[0],windowWidth,windowWidth);
        SetCAxes(0,1);
    }
    // Declaring iteration variables
    unsigned int t;
    unsigned int i;
    unsigned int l;
    unsigned int j;

    //double forces_x;
    //double forces_y;

    // Thread related declarations
    pthread_t threads[NUM_THREADS];
    if (N % NUM_THREADS != 0) {
        printf("Number of elements is not divisible with number of threads.\n");
        return 0;
    }
    int particles_per_thread = N / NUM_THREADS;
    thread_args_t** arguments = (thread_args_t**) malloc(NUM_THREADS*sizeof(thread_args_t*));

    for (j = 0; j < NUM_THREADS; j++) {
        arguments[j] = malloc(sizeof(thread_args_t));
        arguments[j]->start = j*particles_per_thread;
        arguments[j]->amount_of_particles = particles_per_thread;
        arguments[j]->forces = malloc(sizeof(double)*2*particles_per_thread);
        arguments[j]->root = NULL;
        arguments[j]->theta_max = theta_max;
        arguments[j]->G = G;

    }

    for (t = 0; t < nsteps; t++) {
        if (graphics != 0) {
            ClearScreen();
            for (l = 0; l < N; l++) {
                DrawCircle(particles[l].x, particles[l].y, 1, 1, particles[l].mass*0.002, 0);
            }
            Refresh();
            usleep(2000);
        }

        root = create_node(0, 0, 1, 1, N);
        memcpy(root->particles, particles, N*sizeof(particle_t));
        make_tree(root, graphics);

        for (j = 0; j < NUM_THREADS; j++) {
            arguments[j]->root = root;
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
            arguments[j]->root = NULL;
        }

        delete_tree(root);
        root = NULL;
    }

    // Prints positions for debugging purposes
    /*
    for (i = 0; i < N; i++) {
        printf("(%lf, %lf)\n", particles[i].x, particles[i].y);
    }
    */
    if (graphics != 0) {
        FlushDisplay();
        CloseDisplay();
    }

    FILE *ptr;

    ptr = fopen("result.gal", "wb");
    fwrite(particles, N*sizeof(particle_t), 1, ptr); // Write all data to binary file
    fclose(ptr);

    free(particles);

    for (j = 0; j < NUM_THREADS; j++) {
        free(arguments[j]->forces);
        free(arguments[j]);
    }
    free(arguments);
    printf("Galsim program took %7.3f wall seconds.\n", get_timings() - time);
    return 0;

}
