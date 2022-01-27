#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct node {
    int index;
    double min, max;
    struct node *next;
} node_t;

node_t* create_node(int index, double min, double max) {
    node_t *node = malloc(sizeof(node_t));
    node->index = index;
    node->min = min;
    node->max = max;
    node->next = NULL;
    return node;
}

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
        while (curr->next != NULL && curr->next->index < node->index) {
            curr = curr->next;
        }
        node->next = curr->next;
        curr->next = node;
    }
}

void delete_node(node_t *head_ptr, int index) {
    node_t *curr = head_ptr;
}

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

int main() {
    int run = 1;
    int index;
    double min, max;
    char command;
    node_t *head = NULL;

    do {
        printf("Enter command: ");
        scanf(" %c", &command);
        if (command == 'A' || command == 'D'){
            scanf("%d", &index);
            if (command == 'A') {
                scanf("%lf %lf", &min, &max);
            }
        }
        switch (command) {
            case 'A': {
                node_t *new_node = create_node(index, min, max);
                insert(&(head), new_node);
                break;
            }
            case 'D': {
                printf("D\n");
                break;
            }
            case 'P': {
                print(head);
                break;
            }
            case 'Q': {
                run = 0;
                break;
            }
        }
    } while (run);
}