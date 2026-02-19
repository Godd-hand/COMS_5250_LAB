/*
Structs defined in a header file
*/

#ifndef __QUAD_STRUCT_H__
#define __QUAD_STRUCT_H__

// Define a point struct to hold x, y coordinates
typedef struct point point;
struct point
{
    double x;
    double y;
};


// define the quadrilateral struct
typedef struct quad quad;
struct quad
{
    point node1;
    point node2;
    point node3;
    point node4;
    double perimeter;
    double area;
    double angles[4];
};

// Function declarations to be implemented in separate files
void quad_perimeter(quad* q);
void quad_area(quad* q);
void quad_angles(quad* q);

#endif
