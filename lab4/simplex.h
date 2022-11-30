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

int init(simplex_t *s, int m, int n, double **a, double *b, double *c, double *x, double y, int *var);
int select_nonbasic(simplex_t *s);
int prepare(simplex_t *s, int k);
int initial(simplex_t *s, int m, int n, double **a, double *b, double *c, double *x, double y, int *var);
void pivot(simplex_t *s, int row, int col);
double xsimplex(int m, int n, double **a, double *b, double *c, double *x, double y, int *var, int h);
double simplex(int m, int n, double **a, double *b, double *c, double *x, double y);

#endif // SIMPLEX_H