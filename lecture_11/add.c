/*
C file to implement the add_poly function for adding two polynomials represented as linked lists.
*/

#include <stdlib.h>
#include "poly.h"

poly_node* add_poly(poly_node* p1, poly_node* p2)
{
    poly_node* result = NULL;
    poly_node* ptr = p1;
    
    while (ptr != NULL)
    {
        insert_term(&result, ptr->coeff, ptr->power);
        ptr = ptr->next;
    }
    
    ptr = p2;
    while (ptr != NULL)
    {
        insert_term(&result, ptr->coeff, ptr->power);
        ptr = ptr->next;
    }
    
    return result;
}