#include <stdlib.h>
#include "poly.h"

void delete_poly(poly_node** head)
{
    poly_node* temp;
    while (*head != NULL)
    {
        temp = *head;
        *head = (*head)->next;
        free(temp);
    }
}