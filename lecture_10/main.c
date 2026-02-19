#include <stdio.h>
#include "quad_struct.h"

void print_quad_info(quad q, const char* name)
{
    printf("\n=== %s ===\n", name);
    printf(" Nodes:\n");
    printf("  #1: ( %7.2f, %7.2f )\n", q.node1.x, q.node1.y);
    printf("  #2: ( %7.2f, %7.2f )\n", q.node2.x, q.node2.y);
    printf("  #3: ( %7.2f, %7.2f )\n", q.node3.x, q.node3.y);
    printf("  #4: ( %7.2f, %7.2f )\n", q.node4.x, q.node4.y);
    printf(" Perimeter: %.4f\n", q.perimeter);
    printf(" Area:      %.4f\n", q.area);
    printf(" Angles (degrees): %.2f, %.2f, %.2f, %.2f\n", 
           q.angles[0], q.angles[1], q.angles[2], q.angles[3]);
    printf(" Sum of angles:    %.2f\n", q.angles[0] + q.angles[1] + q.angles[2] + q.angles[3]);
}

int main()
{
    // Test 1: Simple Square
    quad square;
    square.node1.x = 0.0; square.node1.y = 0.0;
    square.node2.x = 2.0; square.node2.y = 0.0;
    square.node3.x = 2.0; square.node3.y = 2.0;
    square.node4.x = 0.0; square.node4.y = 2.0;

    quad_perimeter(&square);
    quad_area(&square);
    quad_angles(&square);
    
    print_quad_info(square, "Test 1: Square");

    // Test 2: Complex Quadrilateral
    quad complex_q;
    complex_q.node1.x = 0.0; complex_q.node1.y = 0.0;
    complex_q.node2.x = 5.0; complex_q.node2.y = 1.0;
    complex_q.node3.x = 4.0; complex_q.node3.y = 4.0;
    complex_q.node4.x = 1.0; complex_q.node4.y = 3.0;

    quad_perimeter(&complex_q);
    quad_area(&complex_q);
    quad_angles(&complex_q);

    print_quad_info(complex_q, "Test 2: Complex Quadrilateral");

    return 0;
}