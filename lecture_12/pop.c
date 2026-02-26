#include <stdlib.h>
#include "node.h"

void Pop(node** top, double* a, double* b)
{
    node* temp = *top;

    if (temp == NULL)
    { 
        return; 
    }
    else
    { 
        temp = temp->next; // 
    }
    
    *a = (*top)->a; // Retrieving the values of a and b from the top node before popping it
    *b = (*top)->b;
    free(*top);
    *top = temp;

    node* ptr = *top; // Updating the position of all nodes in the stack after popping the top node
    while (ptr != NULL)
    {
        ptr->position = ptr->position - 1;
        ptr = ptr->next;
    }
}