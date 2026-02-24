/*
C file to implement the insert_term function for a polynomial linked list.
*/

#include <stdio.h>
#include <stdlib.h>
#include "poly.h"

void insert_term(poly_node** head, double coeff, int power)
{
    if (coeff == 0.0) 
    {
        return;
    }

    poly_node* current = *head;
    poly_node* prev = NULL;

    while (current != NULL && current->power > power)
    {
        prev = current;
        current = current->next;
    }

    if (current != NULL && current->power == power)
    {
        current->coeff += coeff;
        
        if (current->coeff == 0.0)
        {
            if (prev == NULL)
            {
                *head = current->next;
            }
            else
            {
                prev->next = current->next;
            }
            free(current);
        }
        return;
    }

    poly_node* temp = (poly_node*)malloc(sizeof(poly_node));
    temp->coeff = coeff;
    temp->power = power;
    temp->next = current;

    if (prev == NULL)
    {
        *head = temp;
    }
    else
    {
        prev->next = temp;
    }
}
