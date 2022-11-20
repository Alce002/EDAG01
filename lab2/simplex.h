#ifndef SIMPLEX_H
#define SIMPLEX_H

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

#endif // SIMPLEX_H