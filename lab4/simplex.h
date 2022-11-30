#ifndef SIMPLEX_H
#define SIMPLEX_H

#define EPSILON 1e-6

typedef struct simplex_t {
    int m; /* Constraints. */
    int n; /* Decision variables. */
    int *var; //[n+m+1]; /* 0..n - 1 are nonbasic. */
    double **a; // [m][n+1]; /* A . */
    double *b; // [m]; /* b . */
    double *x; // [n+1]; /* x . */
    double *c; // [n]; /* c . */
    double y; /* y. */
} simplex_t;

double simplex(int m, int n, double **a, double *b, double *c, double *x, double y); 
double** make_matrix(int m, int n);

#endif // SIMPLEX_H
