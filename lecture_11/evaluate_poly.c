#include <math.h>
#include <stdlib.h>
#include "poly.h"

double eval_poly(poly_node* head, double x)
{
    double sum = 0.0;
    poly_node* ptr = head;
    
    while (ptr != NULL)
    {
        sum += ptr->coeff * pow(x, ptr->power);
        ptr = ptr->next;
    }
    
    return sum;
}