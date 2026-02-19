#include <math.h>
#include "quad_struct.h"

// Helper function to calculate distance between two points
double distance(point p1, point p2)
{
    return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
}

void quad_perimeter(quad* q)
{
    double d1 = distance(q->node1, q->node2);
    double d2 = distance(q->node2, q->node3);
    double d3 = distance(q->node3, q->node4);
    double d4 = distance(q->node4, q->node1);
    
    q->perimeter = d1 + d2 + d3 + d4;
}