#include <stdio.h>
#include "node.h"

void PrintNode(node* top)
{
    printf(" %3i | %8.5f | %8.5f | %15p | %15p |\n", 
           top->position, top->a, top->b, (void*)top, (void*)top->next);
    
    if (top->next == NULL)
    { 
        return; 
    }
    PrintNode(top->next);
}

void DisplayStack(node* top)
{
    if (top == NULL)
    { 
        printf(" Stack is empty.\n"); 
        return; 
    }
    printf(" ------------------------------------------------------------------\n");
    printf(" |Pos|    a     |    b     |     Address     |      Next       |\n");
    printf(" ------------------------------------------------------------------\n");
    PrintNode(top);
    printf(" ------------------------------------------------------------------\n");
}