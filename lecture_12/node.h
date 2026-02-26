#ifndef NODE_H
#define NODE_H

typedef struct node node;
struct node {
    int position;
    double a;
    double b;
    node* next;
};

void Push(double a, double b, node** top);
void Pop(node** top, double* a, double* b);
void Peek(node* top, double* a, double* b);
void DisplayStack(node* top);
void GetStackSize(node* top, int* stack_size);
void DeleteStack(node** top);

#endif