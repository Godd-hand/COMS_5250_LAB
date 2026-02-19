#include <math.h>
#include "quad_struct.h"

// Helper function for distance
double dist(point p1, point p2)
{
    return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
}

// Computes the inner angle at p2 formed by p1-p2 and p2-p3 using Law of Cosines
double angle_law_of_cosines(point p1, point p2, point p3)
{
    double a = dist(p1, p2);
    double b = dist(p2, p3);
    double c = dist(p1, p3); // Opposite side
    
    // Law of cosines: c^2 = a^2 + b^2 - 2ab*cos(theta)
    double cos_val = (a*a + b*b - c*c) / (2.0 * a * b);
    
    return acos(cos_val) * (180.0 / M_PI); // Convert radians to degrees
}

void quad_angles(quad* q)
{
    q->angles[0] = angle_law_of_cosines(q->node4, q->node1, q->node2); // Angle at node 1
    q->angles[1] = angle_law_of_cosines(q->node1, q->node2, q->node3); // Angle at node 2
    q->angles[2] = angle_law_of_cosines(q->node2, q->node3, q->node4); // Angle at node 3
    q->angles[3] = angle_law_of_cosines(q->node3, q->node4, q->node1); // Angle at node 4
}