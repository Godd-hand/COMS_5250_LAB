#include "node.h"
#include <stddef.h>

void Peek(node* top, double* a, double* b)
{
    if (top != NULL) // Checking if the stack is not empty before peeking the top node
    {
        *a = top->a;
        *b = top->b;
    }

    
}