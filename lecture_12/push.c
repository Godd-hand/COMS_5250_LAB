#include <stdio.h>
#include <stdlib.h>
#include "node.h"

void Push(double a, double b, node** top)
{
    if (*top == NULL)
    {
        *top = (node*)malloc(sizeof(struct node));
        (*top)->next = NULL;
        (*top)->a = a;
        (*top)->b = b;
        (*top)->position = 1;
    }
    else
    {
        node* temp; // Creating a new node and pushing it onto the stack
        temp = (node*)malloc(sizeof(struct node)); // allocating memory for the new node
        temp->next = *top; 
        temp->a = a;
        temp->b = b;
        temp->position = 1;
        *top = temp;

        node* ptr = (*top)->next; // Updating the position of all nodes in the stack after pushing a new node
        while (ptr != NULL)
        {
            ptr->position = ptr->position + 1;
            ptr = ptr->next;
        }
    }
}