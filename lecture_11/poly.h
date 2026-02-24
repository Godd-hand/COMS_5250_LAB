/*
Header file for polynomial linked list implementation.
*/

#ifndef __POLY_H__
#define __POLY_H__

typedef struct poly_node poly_node;
struct poly_node
{
    double coeff;
    int power;
    poly_node* next;
};

void insert_term(poly_node** head, double coeff, int power);
poly_node* add_poly(poly_node* p1, poly_node* p2);
double eval_poly(poly_node* head, double x);
void print_poly(poly_node* head);
void delete_poly(poly_node** head);

#endif