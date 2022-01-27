#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct node {
    int index;
    double min, max;
    struct node *next;
} node_t;

/*
    Creates a node using the given data
*/
node_t* create_node(int index, double min, double max) {
    node_t *node = malloc(sizeof(node_t));
    node->index = index;
    node->min = min;
    node->max = max;
    node->next = NULL;
    return node;
}

/*
    Function used to insert nodes into the linked list. The function also sorts the linked list
    automatically, by comparing the integers indices.
*/
void insert(node_t **head_ptr, node_t *node) {
    if (*head_ptr == NULL || node->index < (*head_ptr)->index) {
        node->next = *head_ptr;
        *head_ptr = node;
    }
    else if ((*head_ptr)->index == node->index) {
        (*head_ptr)->min = node->min;
        (*head_ptr)->max = node->max;
    }
    else {
        node_t *curr = *head_ptr;
        while (curr->next != NULL && curr->next->index <= node->index) {
            curr = curr->next;
        }
        if (curr->index != node->index) {
            node->next = curr->next;
            curr->next = node;
        }
        else {
            curr->min = node->min;
            curr->max = node->max;
        }
    }
}

/*
    Function used to delete a specific node given by its index.
*/
void delete_node(node_t **head_ptr, int index) {
    node_t *curr = *head_ptr;
    node_t *temp;

    if (curr == NULL) {
        printf("Database is empty. \n");
        return;
    }

    if (curr != NULL && curr->index == index) {
        *head_ptr = curr->next;
        free(curr);
        return;
    }

    while (curr != NULL && curr->index != index) {
        temp = curr;
        curr = curr->next;
    }

    if (curr == NULL) {
        printf("Index not in database.");
        return;
    }
    else {
        temp->next = curr->next;
        free(curr);
    }
}

/*
    Prints the linked list in a formatted way.
*/
void print(node_t *head_ptr) {
    node_t *curr = head_ptr;
    if (curr == NULL) {
        printf("Database is empty.\n");
    }
    else {
        printf("day \t min \t \t max\n");
        while (curr != NULL) {
            printf("%d \t %lf \t %lf\n", curr->index, curr->min, curr->max);
            curr = curr->next;
        }
    }
}

/*
    Deletes all remaining nodes as the program is closed.
*/
void delete_all(node_t **head_ptr) {
    node_t *prev = *head_ptr;
    node_t *curr;

    while (prev != NULL) {
        curr = prev->next;
        free(prev);
        prev = curr;

    }
}

int main() {
    int run = 1;
    int index, read_input;
    double min, max;
    char command;
    node_t *head = NULL;

    do {
        printf("Enter command: ");
        read_input = scanf(" %c", &command);

        if (read_input != 1) {
            printf("Invalid type in command. First argument must be a char.\n");
            continue;
        }

        if (command == 'A' || command == 'D'){
            read_input = scanf("%d", &index);

            if (read_input != 1 || index < 1 || index > 31) {
                printf("Invalid input in index. Index must be an integer between 1 and 31.\n");
                continue;
            }

            if (command == 'A') {
                read_input = scanf("%lf %lf", &min, &max);

                if (read_input != 2 || isnan(min) || isnan(max)) {
                    printf("Invalid input in min and max. Min and max must be real doubles.\n");
                    continue;
                }
            }
        }
        switch (command) {
            case 'A': {
                node_t *new_node = create_node(index, min, max);
                insert(&(head), new_node);
                break;
            }
            case 'D': {
                delete_node(&(head), index);
                break;
            }
            case 'P': {
                print(head);
                break;
            }
            case 'Q': {
                run = 0;
                delete_all(&(head));
                break;
            }
            default: {
                printf("Invalid command.\n");
            }
        }
    } while (run);
}