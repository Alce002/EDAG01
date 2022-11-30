#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "simplex.h"

#define EPSILON 1e-6

void print_all(int m, int n, double **a, double *b, double *c, double y);

int init(simplex_t *s, int m, int n, double **a, double *b, double *c, double *x, double y, int *var);
int select_nonbasic(simplex_t *s);
int prepare(simplex_t *s, int k);
int initial(simplex_t *s, int m, int n, double **a, double *b, double *c, double *x, double y, int *var);
void pivot(simplex_t *s, int row, int col);
double xsimplex(int m, int n, double **a, double *b, double *c, double *x, double y, int *var, int h);
double simplex(int m, int n, double **a, double *b, double *c, double *x, double y);

int init(simplex_t *s, int m, int n, double **a, double *b, double *c, double *x, double y, int *var) {
    int i, k;
    (*s) = (simplex_t){m, n, var, a, b, x, c, y};
    if(s->var == NULL) {
        s->var = calloc(m+n+1,sizeof(*(s->var)));
        for(i = 0; i < m + n; i++) 
            s->var[i] = i;
    }
    for(k = 0, i = 1; i < m; i++)
        if (b[i] < b[k])
            k = i;
    return k;
}

int select_nonbasic(simplex_t *s) {
    int i;
    for(i = 0; i < s->n; i++) {
        //printf("c[%d]: %lf\n",i, s->c[i]);
        //printf("%lf\n", EPSILON);
        if(s->c[i] > EPSILON)
            return i;
    }
    return -1;
}

int prepare (simplex_t *s, int k) {
    int m = s->m;
    int n = s->n;
    int i;

    // make room for x_m + n at s.var[n] by moving s.var[n..n+m-1] one
    // step to the right.

    for(i = m+n; i > n; i--) 
        s->var[i] = s->var[i-1];
    
    s->var[n] = m + n;
    // add x_m + n to each constraint
    n++;
    
    for(i = 0; i < m; i++) 
        s->a[i][n-1] = -1;

    s->x = calloc(m+n, sizeof(*(s->x)));
    s->c = calloc(n, sizeof(*(s->c)));

    s->c[n-1] = -1;
    s->n = n;
    
    pivot(s,k,n-1);
}


int initial(simplex_t *s, int m, int n, double **a, double *b, double *c, double *x, double y, int *var) {
    int i, j, k;
    double w;
    k = init(s,m,n,a,b,c,x,y,var);
    if(b[k] >= 0)
        return 1;
    prepare(s, k);
    n = s->n;
    s->y = xsimplex(m,n,s->a, s->b, s->c, s->x, 0, s->var, 1);
    for(i = 0; i < m+n; i++) {
        if(s->var[i] == m+n-1) {
            if(abs(s->x[i])>EPSILON) {
                free(s->x);
                free(s->c);
                return 0;
            } else {
                break;
            }
        }
    }
    if(i >= n) {
        for(j = k = 0; k < n; k++) {
            if(abs(s->a[i-n][k]) > abs(s->a[i-n][j])) 
                j = k;
        }
        pivot(s,i-n,j);
        i = j;
    }
    if(i < n - 1) {
        k = s->var[i];
        s->var[i] = s->var[n-1];
        s->var[n-1] = k;
        for(k = 0; k < m; k++) {
            w = s->a[k][n-1];
            s->a[k][n-1] = s->a[k][i];
            s->a[k][i] = w;
        }
    }
    free(s->c);
    s->c = c;
    s->y = y;
    for(k = n-1; k < n+m-1; k++) 
        s->var[k] = s->var[k+1];
    n = s->n = s->n - 1;
    double *t = calloc(n, sizeof(double));
    for(k = 0; k < n; k++) {
        for (j = 0; j < n; j++) {
            if (k == s->var[j]) {
                t[j] = t[j] + s->c[k];
                goto next_k;
            }
        }
        for (j = 0; j < m; j++) {
            if(s->var[n+j] == k)
                break;
        }
        s->y = s->y + s->c[k] * s->b[j];
        for (i = 0; i < n; i++) 
            t[i] = t[i] - s->c[k] * s->a[j][i];
    next_k:;    
    }
    for(i = 0; i < n; i++)
        s->c[i] = t[i];
    free(t);
    free(s->x);
    return 1;
}

void pivot(simplex_t *s, int row, int col) {
    double **a = s->a;
    double *b = s->b;
    double *c = s->c;
    int m = s->m;
    int n = s->n;
    int i, j, t;

    t = s->var[col];
    s->var[col] = s->var[n + row];
    s->var[n + row] = t;
    s->y = s->y + c[col] * b[row] / a[row][col];

    for(i = 0; i < n; i++) {
        if(i != col) {
            c[i] = c[i] - c[col] * a[row][i] / a[row][col];
        }
    }

    c[col] = -c[col] / a[row][col];

	for(i = 0; i < m; i++) {
		if(i != row) {
			b[i] = b[i] - a[i][col] * b[row] / a[row][col];
		}
	}

    for(i = 0; i < m; i++) {
        if(i != row) {
            for(j = 0; j < n; j++) {
                if(j != col) {
                    a[i][j] = a[i][j] - a[i][col] * a[row][j] / a[row][col];
                }
            }
        }
    }

    for(i = 0; i < m; i++) {
        if(i != row) {
            a[i][col] = -a[i][col] / a[row][col];
        }
    }

    for(i = 0; i < n; i++) {
        if(i != col) {
            a[row][i] = a[row][i] / a[row][col];
        }
    }

    b[row] = b[row] / a[row][col];
    a[row][col] = 1 / a[row][col];

	print_all(s->m, s->n, s->a, s->b, s->c, s->y);
}

double xsimplex(int m, int n, double **a, double *b, double *c, double *x, double y, int *var, int h) {
    simplex_t s;
    int i, row, col;

    if(!initial(&s, m, n, a, b, c, x, y, var)) {
        free(s.var);
        return NAN;
    }

    while((col = select_nonbasic(&s)) >= 0) {
        row = -1;
        for(i = 0; i < m; i++) {
            if(a[i][col] > EPSILON &&
               (row < 0 || b[i] / a[i][col] < b[row] / a[row][col])) {
                row = i;
            }
        }
            
        if(row < 0) {
            free(s.var);
            return INFINITY;
        }
            
        pivot(&s, row, col);
    }

    if(h == 0) {
        for(i = 0; i < n; i++) {
            if(s.var[i] < n) {
                x[s.var[i]] = 0;
            }
        }

        for(i = 0; i < m; i++) {
            if(s.var[n + i] < n) {
                x[s.var[n + i]] = s.b[i];
            }
        }

        free(s.var);
    } else {
        for(i = 0; i < n; i++) {
            x[i] = 0;
        }

        for(i = n; i < n + m; i++) {
            x[i] = s.b[i - n];
        }
    }

    return s.y;
}

double simplex(int m, int n, double **a, double *b, double *c, double *x, double y) {
    return xsimplex(m, n, a, b, c, x, y, NULL, 0);
}

double** make_matrix(int m, int n){
    double **a;
    int i;
    a = calloc(m, sizeof(*a));
    for (i = 0; i < m; i++)
        a[i] = calloc(n, sizeof(**a));
    return a;
}

void print_all(int m, int n, double **a, double *b, double *c, double y) {
	printf("------------------------------\n");
	/* Print statements*/
    // m & n
    printf("m = %d, n = %d\n", m, n);

    // c
    printf("max z = ");
    for(int i=0; i < n-1; i++) {
        printf("%10.3lf x%d+", c[i], i);
    }
    printf("%10.3lf x%d\n", c[n-1], n-1);

    // matrix
    for(int y=0; y < m; y++) {
        for(int x=0; x < n-1; x++) {
            printf("%10.3lf x%d+", a[y][x], x);
        }
        printf("%10.3lf x%d \u2264 %10.3lf\n", a[y][n - 1], n-1, b[y]);
    }

	printf("y = %lf\n", y);

}

int main(int argc, char** argv) {
    int n, m;
    double *c, *b;
    double **a;

    scanf("%d %d", &m, &n);

    c = calloc(n, sizeof(*c));
    for(int i=0; i < n; i++) {
        scanf("%lf", &c[i]);
    }

    a = make_matrix(m, n + 1);
    for(int i = 0; i <m; i++) {
        for(int j = 0; j < n; j++) {
            scanf("%lf",&a[i][j]);
        }
    }

    b = calloc(m, sizeof(*b));
    for(int i=0; i < m; i++) {
        scanf("%lf", &b[i]);
    }

	print_all(m, n, a, b, c, 0);

    double *x = calloc(n + 1, sizeof(*x));

    double y = simplex(m, n, a, b, c, x, 0);

    print_all(m, n, a, b, c, y);
    printf("------------------------------\n");
    
    free(x);
    free(c);
    free(b);
    for(int i = 0; i < n; i++) {
        free(a[i]);
    }
    free(a);
    return 0;
} 
