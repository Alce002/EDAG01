#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bnb.h"

#define EPSILON 1e-6


double** make_matrix(int m, int n){
    double **a;
    int i;
    a = calloc(m, sizeof(*a));
    for (i = 0; i < m; i++)
        a[i] = calloc(n, sizeof(**a));
    return a;
}

node_t* initial_node (int m, int n, double** a, double* b, double* c) {
    node_t *p = calloc(1, sizeof(node_t));
    p -> a   = make_matrix(m + 1, n + 1);
    p -> b   = calloc(m + 1, sizeof(double));
    p -> c   = calloc(n + 1, sizeof(double));
    p -> x   = calloc(n + 1, sizeof(double));
    p -> min = calloc(n, sizeof(double));
    p -> max = calloc(n, sizeof(double));
    p -> a = m;
    p -> a = n;
    for(int i = 0; i < n; i++) {
        p -> min[i] = -INFINITY;
        p -> max[i] = INFINITY;
    }
    return p;
}

node_t* extend(node_t* p, int m, int n, double** a, double* b, double* c, int k, double ak, double bk) {
    node_t *q = calloc(1, sizeof(node_t));
    int i, j;
    q -> k  = k;
    q -> ak = ak;
    q -> bk = bk;
    if(ak > 0 && p->max)
        q -> m = p -> m;
    else if(ak < 0 && p -> min[k] > 0)
        q-> m = p -> m;
    else
        q -> m = p->m + 1;
    q -> n   = p -> n;
    q -> h   = -1;
    q -> a   = make_matrix(q->m + 1, q-> n + 1);
    q -> b   = calloc(q->m + 1, sizeof(double));
    q -> c   = calloc(q->n + 1, sizeof(double));
    q -> x   = calloc(q->n + 1, sizeof(double));
    q -> min = calloc(n, sizeof(double));
    q -> max = calloc(n, sizeof(double));



}
int is_integer(double* xp) {
    double x = *xp;
    double r = round(x);
    if(abs(r-x) < EPSILON) {
        *xp = r;
        return 1;
    }
    return 0;
    
}

int integer(node_t* p) {
    int i;
    for(i = 0; i < p->n; i++) {
        if(!is_integer(&p->x[i]))
            return 0;
    }
    return 1;
}

void bound(node_t* p, int h, double* zp, double* x) {
    if(p->z > *zp) {
        *zp = p->z;
        
    }
}
/*int is_finite(double x) {
    if(x == NAN || abs(x) == INFINITY)
        return 0;
    else 
        return 1;
}*/
int branch(node_t* q, double z) {
    double min, max;
    if(q->z < z)
        return 0;
    for(int h = 0; h < q->n; h++) {
        if(!is_integer(&q->x[h]))
            min = 0;
        else
            min = q->min[h];
        max = q->max[h];
        if(floor(q->x[h]) < min || ceil(q->x[h]) > max )
            continue;
        q->h = h;
        q->xh = q-> x[h];

        for(int i = 0; i < q->n; i++)
            free(q->a[i]);
        free(q->a);
        free(q->b);
        free(q->c);
        free(q->x);
        return 1;
    }
    return 0;
}
void succ (node_t* p, int h, int m, int n, double** a, double* b, double* c, int k, double ak, double bk, double* zp, double* x) {
    node_t* q = extend(p, m, n, a, b, c, k, ak, bk);
    if(q == NULL)
        return;
    q->z = simplex(q->m, q->n, q->a, q->b, q->c, q->x, 0);
    if(isfinite(q->z)) {
        if(integer(q)) {
            bound(q, h, zp, x);
        } else if(branch(q, *zp)) {
            
            return;
        }
    }
}
double intopt(int m, int n, double** a, double* b, double* c, double* x) {
    
}