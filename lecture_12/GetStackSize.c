#include <stdlib.h>
#include "node.h"

void GetStackSize(node* top, int* stack_size)
{
    if (top == NULL)
    { 
        *stack_size = 0; 
        return; 
    }
    if (top->next == NULL)
    { 
        *stack_size = top->position; 
        return; 
    }
    GetStackSize(top->next, stack_size);
}