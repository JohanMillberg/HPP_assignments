#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <pthread.h>
#include "graphics.h"

const double eps_0 = 1E-3;

typedef struct tree_node {
    double q_x, q_y, width, height;
    struct tree_node* topLeft;
    struct tree_node* topRight;
    struct tree_node* botLeft;
    struct tree_node* botRight;
    int amount_particles;
    double* particles;
    //particle_t* particles;
    double m;
    double m_x;
    double m_y;
} tree_node_t;

/* ------------------------- Initiating Methods --------------------------------------------------------------*/
void create_node(tree_node_t* node, double x, double y, double width, double height, int N);
tree_node_t* make_node(tree_node_t* node, tree_node_t* node_parent, int counter, int position_index);
void make_tree(tree_node_t* node);
void delete_tree(tree_node_t* node);

/* -------------------------- Methods for Building Tree ----------------------------------------------------- */
void create_node(tree_node_t* node, double x, double y, double width, double height, int N) {
    node->q_x = x; node->q_y = y; node->width = width; node->height = height;
    node->amount_particles = N;
    //node->particles = (particle_t* )malloc(N*sizeof(particle_t));
    node->particles = (double*)malloc(N*sizeof(double*)*6);
    if (!node->particles) {
        printf("Error! Not enough memory!\n");
    }
    node->topLeft = NULL;
    node->topRight = NULL;
    node->botLeft = NULL;
    node->botRight = NULL;
    node->m = 0;
    node->m_x = 0;
    node->m_y = 0;
}

tree_node_t* make_node(tree_node_t* node, tree_node_t* node_parent, int counter, int position_index) {
    int i;
    int j = 0;
    int x_bool, y_bool;
    double x_parent, y_parent;
    double new_width = (node_parent->width)/2;
    double new_height = (node_parent->height)/2;
    double new_x, new_y;
    int n = node_parent->amount_particles;
    x_parent = node_parent->q_x;
    y_parent = node_parent->q_y;
    
    switch (position_index) {
        /* North West Quadrant*/
        case 0: {
            new_x = x_parent;
            new_y = y_parent;

            create_node(node, new_x, new_y, new_width, new_height, counter);
            for (i = 0; i < node_parent->amount_particles; i++) {
                x_bool = (node_parent->particles[i*6 + 0] >= x_parent &&
                 node_parent->particles[i*6 + 0] < x_parent + new_width);

                y_bool = (node_parent->particles[i*6 + 1] >= y_parent &&
                 node_parent->particles[i*6 + 1] < y_parent + new_height);
                if (x_bool && y_bool) {

                    node->particles[j*6 + 0] = node_parent->particles[i*6 + 0];
                    node->particles[j*6 + 1] = node_parent->particles[i*6 + 1];
                    node->particles[j*6 + 2] = node_parent->particles[i*6 + 2];
                    node->particles[j*6 + 3] = node_parent->particles[i*6 + 3];
                    node->particles[j*6 + 4] = node_parent->particles[i*6 + 4];
                    node->particles[j*6 + 5] = node_parent->particles[i*6 + 0];
                    node->m += node->particles[j*6 + 2];
                    node->m_x += node->particles[j*6 + 0]*node->particles[j*6 + 2];
                    node->m_y += node->particles[j*6 + 1]*node->particles[j*6 + 2];
                    j++;
                }
            }
            node->m_x = node->m_x/node->m;
            node->m_y = node->m_y/node->m;
            break;
        }
        /* North East Quadrant*/
        case 1: {
            new_x = x_parent + new_width;
            new_y = y_parent;
            
            create_node(node, new_x, new_y, new_width, new_height, counter);
            for (i = 0; i < node_parent->amount_particles; i++) {
                x_bool = (node_parent->particles[i*6 + 0] >= x_parent + new_width &&
                 node_parent->particles[i*6 + 0] <= x_parent + node_parent->width);

                y_bool = (node_parent->particles[i*6 + 1] >= y_parent &&
                 node_parent->particles[i*6 + 1] < y_parent + new_height);

                if (x_bool && y_bool) {
                    node->particles[j*6 + 0] = node_parent->particles[i*6 + 0];
                    node->particles[j*6 + 1] = node_parent->particles[i*6 + 1];
                    node->particles[j*6 + 2] = node_parent->particles[i*6 + 2];
                    node->particles[j*6 + 3] = node_parent->particles[i*6 + 3];
                    node->particles[j*6 + 4] = node_parent->particles[i*6 + 4];
                    node->particles[j*6 + 5] = node_parent->particles[i*6 + 0];
                    node->m += node->particles[j*6 + 2];
                    node->m_x += node->particles[j*6 + 0]*node->particles[j*6 + 2];
                    node->m_y += node->particles[j*6 + 1]*node->particles[j*6 + 2];
                    j++;
                }
            }
            node->m_x = node->m_x/node->m;
            node->m_y = node->m_y/node->m;
            break;
        }
        
        /* South West Quadrant*/
        case 2: {
            new_x = x_parent;
            new_y = y_parent + new_height;

            create_node(node, new_x, new_y, new_width, new_height, counter);
            for (i = 0; i < node_parent->amount_particles; i++) {
                x_bool = (node_parent->particles[i*6 + 0] >= x_parent &&
                 node_parent->particles[i*6 + 0] < x_parent + new_width);

                y_bool = (node_parent->particles[i*6 + 1] >= y_parent + new_height &&
                 node_parent->particles[i*6 + 1] <= y_parent + node_parent->height);
                if (x_bool && y_bool) {
                    node->particles[j*6 + 0] = node_parent->particles[i*6 + 0];
                    node->particles[j*6 + 1] = node_parent->particles[i*6 + 1];
                    node->particles[j*6 + 2] = node_parent->particles[i*6 + 2];
                    node->particles[j*6 + 3] = node_parent->particles[i*6 + 3];
                    node->particles[j*6 + 4] = node_parent->particles[i*6 + 4];
                    node->particles[j*6 + 5] = node_parent->particles[i*6 + 0];
                    node->m += node->particles[j*6 + 2];
                    node->m_x += node->particles[j*6 + 0]*node->particles[j*6 + 2];
                    node->m_y += node->particles[j*6 + 1]*node->particles[j*6 + 2];
                    j++;
                }
            }
            node->m_x = node->m_x/node->m;
            node->m_y = node->m_y/node->m;
            break;
        }
        
        /* South East Quadrant*/
        case 3: {
            new_x = x_parent + new_width;
            new_y = y_parent + new_height;

            create_node(node, new_x, new_y, new_width, new_height, counter);
            for (i = 0; i < node_parent->amount_particles; i++) {
                x_bool = (node_parent->particles[i*6 + 0] >= x_parent + new_width &&
                 node_parent->particles[i*6 + 0] <= x_parent + node_parent->width);
                y_bool = (node_parent->particles[i*6 + 1] >= y_parent + new_height &&
                 node_parent->particles[i*6 + 1] <= y_parent + node_parent->height);
                if (x_bool && y_bool) {
                    node->particles[j*6 + 0] = node_parent->particles[i*6 + 0];
                    node->particles[j*6 + 1] = node_parent->particles[i*6 + 1];
                    node->particles[j*6 + 2] = node_parent->particles[i*6 + 2];
                    node->particles[j*6 + 3] = node_parent->particles[i*6 + 3];
                    node->particles[j*6 + 4] = node_parent->particles[i*6 + 4];
                    node->particles[j*6 + 5] = node_parent->particles[i*6 + 0];
                    node->m += node->particles[j*6 + 2];
                    node->m_x += node->particles[j*6 + 0]*node->particles[j*6 + 2];
                    node->m_y += node->particles[j*6 + 1]*node->particles[j*6 + 2];
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

void make_tree(tree_node_t* node) {
    int i;
    int bool_NW[2], bool_NE[2], bool_SW[2], bool_SE[2];
    int counter[4];
    tree_node_t* new_node;
    int n = node->amount_particles;

    for (i = 0; i < 4; i++) {
        counter[i] = 0;
    }

    if (n > 1) {
        // Check if the particles are in this quadrant:
        /* Nort West Quadrant*/
        for (i = 0; i < n; i++) {
            bool_NW[0] = (node->particles[i*6 + 0] >= node->q_x &&
                node->particles[i*6 + 0] < node->q_x + node->width/2);

            bool_NW[1] = (node->particles[i*6 + 1] >= node->q_y &&
                node->particles[i*6 + 1] < node->q_y + node->height/2);

            if (bool_NW[0] && bool_NW[1]) {
                counter[0]++;
            }

            /* Nort East Quadrant*/
            bool_NE[0] = (node->particles[i*6 + 0] >= node->q_x + node->width/2 &&
                node->particles[i*6 + 0] <= node->q_x + node->width);

            bool_NE[1] = (node->particles[i*6 + 1] >= node->q_y &&
                node->particles[i*6 + 1] < node->q_y + node->height/2);
            if (bool_NE[0] && bool_NE[1]) {
                counter[1]++;
            }

            /* South West Quadrant*/
            bool_SW[0] = (node->particles[i*6 + 0] >= node->q_x &&
                node->particles[i*6 + 0] < node->q_x + node->width/2);

            bool_SW[1] = (node->particles[i*6 + 1] >= node->q_y + node->height/2 &&
                node->particles[i*6 + 1] <= node->q_y + node->height);
            if (bool_SW[0] && bool_SW[1]) {
                counter[2]++;
            }

            /* South East Quadrant*/
            bool_SE[0] = (node->particles[i*6 + 0] >= node->q_x + node->width/2 &&
                node->particles[i*6 + 0] <= node->q_x + node->width);

            bool_SE[1] = (node->particles[i*6 + 1] >= node->q_y + node->height/2 &&
                node->particles[i*6 + 1] <= node->q_y + node->height);
            if (bool_SE[0] && bool_SE[1]) {
                counter[3]++;
            }
        }
        /* Freeing particles */
        /* North West Quadrant*/
        if (counter[0] > 1) {
            node->topLeft = (tree_node_t*)malloc(sizeof(tree_node_t));
            new_node = make_node(node->topLeft, node, counter[0], 0);
            make_tree(new_node);
        }

        /* North East Quadrant*/
        if (counter[1] > 1) {
            node->topRight = (tree_node_t*) malloc(sizeof(tree_node_t));
            new_node = make_node(node->topRight, node, counter[1], 1);
            make_tree(new_node);
        }

        /* South West Quadrant*/
        if (counter[2] > 1) {
            node->botLeft = (tree_node_t*) malloc(sizeof(tree_node_t));
            new_node = make_node(node->botLeft, node, counter[2], 2);
            make_tree(new_node);
        }

        /* South East Quadrant*/
        if (counter[3] > 1) {
            node->botRight = (tree_node_t*) malloc(sizeof(tree_node_t));
            new_node = make_node(node->botRight, node, counter[3], 3);
            make_tree(new_node);
        }

    }

}

void print_tree(tree_node_t *node, int level, int quadrant) {
    if (node == NULL) {
        return;
    }
    for (int i = 0; i < level; i++) {
        printf(i == level-1 ? "|-" : " ");
    }
    printf("Level %d, quadrant %d\n", level, quadrant);
    printf("%d\n", node->amount_particles);
    print_tree(node->topLeft, level+1, 1);
    print_tree(node->topRight, level+1, 2);
    print_tree(node->botLeft, level+1, 3);
    print_tree(node->botRight, level+1, 4);
}

/*
Method of reading data from .gal files was inspired by the function read_doubles_from_file
in the given compare_gal_files.c
*/

int set_initial_data(int N, double** particles, const char* filename) {
    int i;
    FILE* file = fopen(filename, "rb");
    if (!file) {
        free(*particles);
        return 0;
    }
    fseek(file, 0L, SEEK_END);
    size_t file_size = ftell(file);
    if (file_size != 6*N*sizeof(double)) {
        return 0;
    }
    fseek(file, 0L, SEEK_SET);

    fread(*particles, sizeof(char), file_size, file);

    fclose(file);
    return 1;
}

/* Method for determining mass distribution and calculating forces */
void traverse_tree(tree_node_t* node, double x, double y, double mass,
        double theta_max, double* forces_x, double*forces_y, double G) {
    double theta, F_const;
    double distance_x, distance_y, r;
    double denom;

    if (node == NULL) {
        return;
    }

    distance_x = (x-node->m_x);
    distance_y = (y-node->m_y);

    r = sqrt(distance_x*distance_x + distance_y*distance_y);

    if (r == 0) {
        return;
    }
    theta = node->width/r;
    //printf("%lf\n", theta);

    if (theta > theta_max) {
        traverse_tree(node->topLeft, x, y, mass, theta_max, forces_x, forces_y, G);
        traverse_tree(node->topRight, x, y, mass, theta_max, forces_x, forces_y, G);
        traverse_tree(node->botLeft, x, y, mass, theta_max, forces_x, forces_y, G);
        traverse_tree(node->botRight, x, y, mass, theta_max, forces_x, forces_y, G);
    }

    else {
        //Calculate x and y forces here
        denom = (r + eps_0)*(r + eps_0)*(r + eps_0);
        F_const = G * mass * (node->m/denom);
        (*forces_x) += F_const * distance_x;
        (*forces_y) += F_const * distance_y;
    }
}

void delete_tree(tree_node_t* node) {
    
    if (node == NULL) {
        return;
    }

    if (node->amount_particles == 1) {
        free(node->particles);
        free(node);
        return;
    }
    
    delete_tree(node->topLeft);
    delete_tree(node->topRight);
    delete_tree(node->botLeft);
    delete_tree(node->botRight);

    //printf("Deleting node: %d\n", node->amount_particles);
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
    double x_pos, y_pos;
    double x, y, mass;
    double denom;

    // Declare constants
    const double G = (double) -100.0 / N;

    // Store the properties of all particles in the array particles
    double *particles = (double*) malloc(N*6*sizeof(double));

    // The sum of the forces in the x and y direction for each particle is stored in forces
    double forces[2*N];

    successful = set_initial_data(N, &particles, filename);

    if (!successful) {
        printf("Error reading initial data file. \n");
        return 0;
    }

    tree_node_t* root = (tree_node_t*)malloc(sizeof(tree_node_t));
    create_node(root, 0, 0, 1, 1, N);
    root->particles = particles;
    
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

    double forces_x;
    double forces_y;
    
    for (t = 0; t < nsteps; t++) {
        make_tree(root);

        for (i = 0; i < N; i++) {
            forces_x = 0;
            forces_y = 0;
            x_pos = root->particles[i*6 + 0];
            y_pos = root->particles[i*6 + 1];
            mass = root->particles[i*6 + 2];
            traverse_tree(root, x_pos, y_pos, mass, theta_max, &forces_x, &forces_y, G);
            forces[i*2 + 0] = forces_x;
            forces[i*2 + 1] = forces_y;
        }

        for (i = 0; i < N; i++) {
            root->particles[i*6 + 3] = root->particles[i*6 + 3] + delta_t*forces[i*2+0]/root->particles[i*6 + 2];
            root->particles[i*6 + 4] = root->particles[i*6 + 4] + delta_t*forces[i*2+1]/root->particles[i*6 + 2];

            root->particles[i*6 + 0] = root->particles[i*6 + 0] + delta_t * root->particles[i*6 + 3];
            root->particles[i*6 + 1] = root->particles[i*6 + 1] + delta_t * root->particles[i*6 + 4];
        }
        delete_tree(root);
    }

    for (i = 0; i < N; i++) {
        printf("Particle %d: x = %.3lf, y = %.3lf, mass = %.3lf\n", i, root->particles[i*6+0], root->particles[i*6 + 1], root->particles[i*6 + 2]);
    }
    FILE *ptr;

    ptr = fopen("result.gal", "wb");
    fwrite(particles, N*sizeof(double)*6, 1, ptr); // Write all data to binary file
    fclose(ptr);

    free(particles);
    free(root);
    printf("Galsim program took %7.3f wall seconds.\n", get_timings() - time);
    return 0;

}
