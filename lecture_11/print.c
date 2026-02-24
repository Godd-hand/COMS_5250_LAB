#include <stdio.h>
#include "poly.h"

void print_poly(poly_node* head)
{
    if (head == NULL)
    {
        printf("0\n");
        return;
    }
    
    poly_node* ptr = head;
    while (ptr != NULL)
    {
        printf("%.1fx^%d ", ptr->coeff, ptr->power);
        ptr = ptr->next;
        if (ptr != NULL)
        {
            printf("+ ");
        }
    }
    printf("\n");
}