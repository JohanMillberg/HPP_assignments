#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <pthread.h>
#include "graphics.h"

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
} tree_node_t;

void create_node(tree_node_t* node, double x, double y, double width, double height, int N);
tree_node_t* make_node(tree_node_t* node, tree_node_t* node_parent, int position_index);
void make_tree(tree_node_t* node);

void create_node(tree_node_t* node, double x, double y, double width, double height, int N) {
    node->q_x = x; node->q_y = y; node->width = width; node->height = height;
    node->amount_particles = N;
    node->particles = malloc(N*sizeof(particle_t));
    node->topLeft = NULL;
    node->topRight = NULL;
    node->botLeft = NULL;
    node->botRight = NULL;
}

tree_node_t* make_node(tree_node_t* node, tree_node_t* node_parent, int position_index) {
    int i;
    int counter = 0;
    int x_bool, y_bool;
    double x_parent, y_parent;
    particle_t* eligible_particles;
    double new_width = node_parent->width/2;
    double new_height = node_parent->height/2;
    double new_x, new_y;
    x_parent = node_parent->q_x;
    y_parent = node_parent->q_y;
    switch (position_index) {
        case 0: {
            new_x = x_parent;
            new_y = y_parent + new_height;
            // Check if the particles are in this quadrant:
            for (i = 0; i < node_parent->amount_particles; i++) {
                x_bool = (node_parent->particles[i].x > x_parent &&
                 node_parent->particles[i].x <= x_parent + new_width);
                y_bool = (node_parent->particles[i].y > y_parent + new_height &&
                 node_parent->particles[i].y <= y_parent + node_parent->height);
                if (x_bool && y_bool) {
                    counter++;
                }
            }

            create_node(node, new_x, new_y, new_width, new_height, counter);
            for (i = 0; i < node_parent->amount_particles; i++) {
                x_bool = (node_parent->particles[i].x > x_parent &&
                 node_parent->particles[i].x <= x_parent + new_width);
                y_bool = (node_parent->particles[i].y > y_parent + new_height &&
                 node_parent->particles[i].y <= y_parent + node_parent->height);
                if (x_bool && y_bool) {
                    node->particles[i].x = node_parent->particles[i].x;
                    node->particles[i].y = node_parent->particles[i].y;
                    node->particles[i].mass = node_parent->particles[i].mass;
                    node->particles[i].vel_x = node_parent->particles[i].vel_x;
                    node->particles[i].vel_y = node_parent->particles[i].vel_y;
                    node->particles[i].brightness = node_parent->particles[i].brightness;
                }
            }
            break;
        }
        case 1: {
            new_x = x_parent + new_width;
            new_y = y_parent + new_height;
            for (i = 0; i < node_parent->amount_particles; i++) {
                x_bool = (node_parent->particles[i].x > new_width &&
                 node_parent->particles[i].x <= x_parent + node_parent->width);
                y_bool = (node_parent->particles[i].y > y_parent + new_height &&
                 node_parent->particles[i].y <= y_parent + node_parent->height);
                if (x_bool && y_bool) {
                    counter++;
                }
            }
            create_node(node, new_x, new_y, new_width, new_height, counter);
            for (i = 0; i < node_parent->amount_particles; i++) {
                x_bool = (node_parent->particles[i].x > x_parent &&
                 node_parent->particles[i].x <= x_parent + new_width);
                y_bool = (node_parent->particles[i].y > y_parent + new_height &&
                 node_parent->particles[i].y <= y_parent + node_parent->height);
                if (x_bool && y_bool) {
                    node->particles[i].x = node_parent->particles[i].x;
                    node->particles[i].y = node_parent->particles[i].y;
                    node->particles[i].mass = node_parent->particles[i].mass;
                    node->particles[i].vel_x = node_parent->particles[i].vel_x;
                    node->particles[i].vel_y = node_parent->particles[i].vel_y;
                    node->particles[i].brightness = node_parent->particles[i].brightness;
                }
            }
            break;
        }
        case 2: {
            new_x = x_parent;
            new_y = y_parent;
            for (i = 0; i < node_parent->amount_particles; i++) {
                x_bool = (node_parent->particles[i].x > x_parent &&
                 node_parent->particles[i].x <= x_parent + new_width);
                y_bool = (node_parent->particles[i].y > y_parent &&
                 node_parent->particles[i].y <= y_parent + new_height);
                if (x_bool && y_bool) {
                    counter++;
                }
            }
            create_node(node, new_x, new_y, new_width, new_height, counter);
            for (i = 0; i < node_parent->amount_particles; i++) {
                x_bool = (node_parent->particles[i].x > x_parent &&
                 node_parent->particles[i].x <= x_parent + new_width);
                y_bool = (node_parent->particles[i].y > y_parent + new_height &&
                 node_parent->particles[i].y <= y_parent + node_parent->height);
                if (x_bool && y_bool) {
                    node->particles[i].x = node_parent->particles[i].x;
                    node->particles[i].y = node_parent->particles[i].y;
                    node->particles[i].mass = node_parent->particles[i].mass;
                    node->particles[i].vel_x = node_parent->particles[i].vel_x;
                    node->particles[i].vel_y = node_parent->particles[i].vel_y;
                    node->particles[i].brightness = node_parent->particles[i].brightness;
                }
            }
            break;
        }
        case 3: {
            new_x = x_parent + new_width;
            new_y = y_parent;
            for (i = 0; i < node_parent->amount_particles; i++) {
                x_bool = (node_parent->particles[i].x > new_width &&
                 node_parent->particles[i].x <= x_parent + node_parent->width);
                y_bool = (node_parent->particles[i].y > y_parent &&
                 node_parent->particles[i].y <= y_parent + new_height);
                if (x_bool && y_bool) {
                    counter++;
                }
            }
            create_node(node, new_x, new_y, new_width, new_height, counter);
            for (i = 0; i < node_parent->amount_particles; i++) {
                x_bool = (node_parent->particles[i].x > x_parent &&
                 node_parent->particles[i].x <= x_parent + new_width);
                y_bool = (node_parent->particles[i].y > y_parent + new_height &&
                 node_parent->particles[i].y <= y_parent + node_parent->height);
                if (x_bool && y_bool) {
                    node->particles[i].x = node_parent->particles[i].x;
                    node->particles[i].y = node_parent->particles[i].y;
                    node->particles[i].mass = node_parent->particles[i].mass;
                    node->particles[i].vel_x = node_parent->particles[i].vel_x;
                    node->particles[i].vel_y = node_parent->particles[i].vel_y;
                    node->particles[i].brightness = node_parent->particles[i].brightness;
                }
            }
            break;
    }


    }
    return node;
}

void make_tree(tree_node_t* node) {
    tree_node_t* new_node = malloc(sizeof(tree_node_t));
    int n = node->amount_particles;

    if (n > 1) {
        node->topLeft = malloc(sizeof(tree_node_t));
        new_node = make_node(node->topLeft, node, 0);
        make_tree(new_node);

        node->topRight = malloc(sizeof(tree_node_t));
        new_node = make_node(node->topRight, node, 1);
        make_tree(new_node);

        node->botLeft = malloc(sizeof(tree_node_t));
        new_node = make_node(node->botLeft, node, 2);
        make_tree(new_node);

        node->botRight = malloc(sizeof(tree_node_t));
        new_node = make_node(node->botRight, node, 3);
        make_tree(new_node);
    }
}

void print_tree(tree_node_t *node, int level) {
    if (node == NULL) {
        return;
    }
    for (int i = 0; i < level; i++) {
        printf(i == level-1 ? "|-" : " ");
    }
    printf("%d\n", node->amount_particles);
    print_tree(node->topLeft, level+1);
    print_tree(node->topRight, level+1);
    print_tree(node->botLeft, level+1);
    print_tree(node->botRight, level+1);
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

    // Declare variables
    double F_x, F_y, F_const;
    double r_x, r_y, r;
    double x, y, mass;
    double denom;

    // Declare constants
    const double G = (double) -100.0 / N;
    const double eps_0 = 0.001;

    // Store the properties of all particles in the array particles
    particle_t *particles = (particle_t*) malloc(N*sizeof(particle_t));

    // The sum of the forces in the x and y direction for each particle is stored in forces
    double forces[2*N];

    successful = set_initial_data(N, &particles, filename);

    tree_node_t* root = malloc(sizeof(tree_node_t));
    create_node(root, 0, 0, 1, 1, N);

    memcpy(root->particles, particles, N*sizeof(particle_t));
    make_tree(root);
    print_tree(root,0);
/*
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

    }

    FILE *ptr;

    ptr = fopen("result.gal", "wb");
    fwrite(particles, N*sizeof(double)*6, 1, ptr); // Write all data to binary file
    fclose(ptr);

    free(particles);
    printf("Galsim program took %7.3f wall seconds.\n", get_timings() - time);
    return 0;
*/
}