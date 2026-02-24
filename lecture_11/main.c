#include <stdio.h>
#include "poly.h"

int main()
{
    poly_node* p1 = NULL;
    poly_node* p2 = NULL;

    printf("Inserting terms for Poly 1...\n");
    insert_term(&p1, 3.0, 2);
    insert_term(&p1, 2.0, 1);
    insert_term(&p1, 5.0, 0);

    printf("Inserting terms for Poly 2...\n");
    insert_term(&p2, -3.0, 2); 
    insert_term(&p2, 4.0, 1);
    insert_term(&p2, -1.0, 0);
    insert_term(&p2, 7.0, 3); 

    printf("\nPoly 1: ");
    print_poly(p1);

    printf("Poly 2: ");
    print_poly(p2);

    poly_node* p3 = add_poly(p1, p2);
    
    printf("Sum:    ");
    print_poly(p3);

    double x = 2.0;
    double result = eval_poly(p3, x);
    printf("\nEvaluating Sum at x = %.1f\n", x);
    printf("Result: %.2f\n", result);

    delete_poly(&p1);
    delete_poly(&p2);
    delete_poly(&p3);

    return 0;
}