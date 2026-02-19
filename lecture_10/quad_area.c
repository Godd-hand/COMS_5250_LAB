#include <math.h>
#include "quad_struct.h"

// Calculate triangle area using cross product methodology from the slides
double triangle_area(point n1, point n2, point n3)
{
    double u[3], v[3], w[3];
    
    u[0] = n3.x - n1.x;
    u[1] = n3.y - n1.y;
    u[2] = 0.0;

    v[0] = n2.x - n1.x;
    v[1] = n2.y - n1.y;
    v[2] = 0.0;

    w[2] = u[0]*v[1] - u[1]*v[0];
    
    return 0.5 * fabs(w[2]);
}

void quad_area(quad* q)
{
    // Split quadrilateral into two triangles: (1, 2, 3) and (1, 3, 4)
    double area1 = triangle_area(q->node1, q->node2, q->node3);
    double area2 = triangle_area(q->node1, q->node3, q->node4);
    
    q->area = area1 + area2;
}
