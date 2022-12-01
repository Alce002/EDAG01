#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPSILON 1e-6

typedef struct node_t {
	int m;
	int n;
	int k;
	int h;
	double xh;
	double ak;
	double bk;
	double *min;
	double *max;
	double **a;
	double *b;
	double *x;
	double *c;
	double z;
	struct node_t *next;
} node_t;

double intopt(int m, int n, double **a, double *b, double *c, double *x);

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

void print_all(int m, int n, double **a, double *b, double *c, double y);

void pivot(simplex_t *s, int row, int col);
double xsimplex(int m, int n, double **a, double *b, double *c, double *x, double y, int *var, int h);

int init(simplex_t *s, int m, int n, double **a, double *b, double *c, double *x, double y, int *var) {
    int i, k;
    *s = (simplex_t){m, n, var, a, b, x, c, y};
    if(s->var == NULL) {
        s->var = calloc(m + n + 1, sizeof(*(s->var)));
        for(i = 0; i < m + n; i++) {
            s->var[i] = i;
        }
    }

    for(k = 0, i = 1; i < m; i++) {
        if(b[i] < b[k]) {
            k = i;
        }
    }

    return k;
}

int select_nonbasic(simplex_t *s) {
    int i;
    for(i = 0; i < s->n; i++) {
        if(s->c[i] > EPSILON) {
            return i;
        }
    }
    return -1;
}

void prepare(simplex_t *s, int k) {
    int m = s->m;
    int n = s->n;
    int i;

    for(i = m + n; i > n; i--) {
        s->var[i] = s->var[i - 1];
    }
    s->var[n] = m + n;

    n++;
    for(i = 0; i < m; i++) {
        s->a[i][n - 1] = -1;
    }

    s->x = calloc(m + n, sizeof(*s->x));
    s->c = calloc(n, sizeof(*s->c));
    s->c[n - 1] = -1;
    s->n = n;
    pivot(s, k, n - 1);
}

int initial(simplex_t *s, int m, int n, double **a, double *b, double *c, double *x, double y, int *var) {
    int i, j, k;
    double w;
    k = init(s, m, n, a, b, c, x, y, var);
    
    if(b[k] >= 0) {
        return 1;
    }

    prepare(s, k);

    n = s->n;
    s->y = xsimplex(m, n, s->a, s->b, s->c, s->x, 0, s->var, 1);

    for(i = 0; i < m + n; i++) {
        if(s->var[i] == m + n - 1) {
            if(fabs(s->x[i]) > EPSILON) {
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
            if(fabs(s->a[i - n][k]) > fabs(s->a[i - n][j])) {
                j = k;
            }
        }
        pivot(s, i - n, j);
        i = j;
    }

    if(i < n - 1) {
        k = s->var[i]; s->var[i] = s->var[n - 1]; s->var[n - 1] = k;
        for(k = 0; k < m; k++) {
            w = s->a[k][n - 1]; s->a[k][n - 1] = s->a[k][i]; s->a[k][i] = w;
        }
    }

    free(s->c);
    s->c = c;
    s->y = y;

    for(k = n - 1; k < n + m - 1; k++) {
        s->var[k] = s->var[k + 1];
    }

    n = s->n = s->n - 1;
    double *t = calloc(n, sizeof(*t));

    for(k = 0; k < n; k++) {
        for(j = 0; j < n; j++) {
            if(k == s->var[j]) {
                t[j] = t[j] + s->c[k];
                goto next_k;
            }
        }

        for(j = 0; j < m; j++) {
            if(s->var[n + j] == k) {
                break;
            }
        }

        s->y = s->y + s->c[k] * s->b[j];
        for(i = 0; i < n; i++) {
            t[i] = t[i] - s->c[k] * s->a[j][i];
        }
        next_k:;
    }

    for(i = 0; i < n; i++) {
        s->c[i] = t[i];
    }

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

	//print_all(s->m, s->n, s->a, s->b, s->c, s->y);
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

	//print_all(s.m, s.n, s.a, s.b, s.c, s.y);

    return s.y;
}

double simplex(int m, int n, double **a, double *b, double *c, double *x, double y) {
    //print_all(m, n, a, b, c, y);
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
	printf("------------------------------\n");
}


void free_node(node_t *p) {
	for(int i = 0; i < p->m + 1; i++) {
		free(p->a[i]);
	}
	free(p->a);
	free(p->b);
	free(p->c);
	free(p->x);
	free(p->min);
	free(p->max);
	free(p);
}

node_t* inital_node(int m, int n, double **a, double *b, double *c) {
	node_t *p = calloc(1, sizeof(node_t));
	p->a = make_matrix(m + 1, n + 1);
	p->b = calloc(m + 1, sizeof(double));
	p->c = calloc(n + 1, sizeof(double));
	p->x = calloc(n + 1, sizeof(double));
	p->min = calloc(n, sizeof(double));
	p->max = calloc(n, sizeof(double));
	p->m = m;
	p->n = n;
	p->next = NULL;
	// copy a, b, c
	int i, j;
	
	for(j = 0; j < m; j++) {
		for(i = 0; i < n; i++) {
			p->a[j][i] = a[j][i];
		}
	}

	for(j = 0; j < m; j++) {
		p->b[j] = b[j];
	}

	for(j = 0; j < n; j++) {
		p->c[j] = c[j];
	}

	for(i = 0; i < n; i++) {
		p->min[i] = -INFINITY;
		p->max[i] = INFINITY;
	}

	return p;
}

node_t* extend(node_t *p, int m, int n, double **a, double *b, double *c, int k, double ak, double bk) {
	node_t *q = calloc(1, sizeof(node_t));
	int i, j;
	q->k = k;
	q->ak = ak;
	q->bk = bk;
	if((ak > 0) && (p->max[k] < INFINITY)) {
		q->m = p->m;
	} else if((ak < 0) && (p->min[k] > 0)) {
		q->m = p->m;
	} else {
		q->m = p->m + 1;
	}

	q->n = p->n;
	q->h = -1;
	q->a = make_matrix(q->m + 1, q->n + 1);
	q->b = calloc(q->m + 1, sizeof(double));
	q->c = calloc(q->n + 1, sizeof(double));
	q->x = calloc(q->n + 1, sizeof(double));
	q->min = calloc(n, sizeof(double));
	q->max = calloc(n, sizeof(double));
	q->next = NULL;

	for(i = 0; i < n; i++) {
		q->min[i] = p->min[i];
		q->max[i] = p->max[i];
	}

	for(j = 0; j < m; j++) {
		for(i = 0; i < n; i++) {
			q->a[j][i] = a[j][i];
		}
	}

	for(j = 0; j < m; j++) {
		q->b[j] = b[j];
	}

	for(j = 0; j < n; j++) {
		q->c[j] = c[j];
	}

	if(ak > 0) {
		if((q->max[k] == INFINITY) || (bk < q->max[k])) {
			q->max[k] = bk;
		}
	} else if((q->min[k] == -INFINITY) || (-bk > q->min[k])) {
		q->min[k] = -bk;
	}
	
	for(i = m, j = 0; j < n; j++) {
		if(q->min[j] > -INFINITY) {
			q->a[i][j] = -1;
			q->b[i] = -q->min[j];
			i++;
		}

		if(q->max[j] < INFINITY) {
			q->a[i][j] = 1;
			q->b[i] = q->max[j];
			i++;
		}
	}

	return q;
}

int is_integer(double *xp) {
	double x = *xp;
	double r = lround(x);
	if(fabs(r - x) < EPSILON) {
		*xp = r;
		return 1;
	} else {
		return 0;
	}
}

int integer(node_t *p) {
	int i;
	for(i = 0; i < p->n; i++) {
		if(!is_integer(&p->x[i])) {
			return 0;
		}
	}
	return 1;
}

void bound(node_t *p, node_t **h, double *zp, double *x) {
	if(p->z > *zp) {
		*zp = p->z;
		for(int i = 0; i < p->n; i++) {
			x[i] = p->x[i];
		}
		// Cull h
		while(*h != NULL && (*h)->z < p->z) {
			node_t *old = *h;
			*h = (*h)->next;
			free_node(old);
		}
		
		if(*h != NULL) {
			node_t *current = *h;
			while(current->next != NULL) {
				if(current->next->z < p->z) {
					node_t *old = current->next;
					current->next = current->next->next;
					free_node(old);
				} else {
					current = current->next;
				}
			}
		}
	}
	//printf("Best z: %lf\n", *zp);
}

int branch(node_t *q, double z) {
	double min, max;
	if(q->z < z) return 0;

	for(int h = 0; h < q->n; h++) {
		if(!is_integer(&q->x[h])) {
			if(q->min[h] == -INFINITY) min = 0;
			else min = q->min[h];

			max = q->max[h];

			if((floor(q->x[h]) < min) || (ceil(q->x[h]) > max)) continue;
			
			q->h = h;
			q->xh = q->x[h];
			return 1;
		}
	}
	return 0;
}

void succ(node_t *p, node_t **h, int m, int n, double **a, double *b, double *c, int k, double ak, double bk, double *zp, double *x) {
	node_t *q = extend(p, m, n, a, b, c, k, ak, bk);
	if(q == NULL) return;
	q->z = simplex(q->m, q->n, q->a, q->b, q->c, q->x, 0);
	if(isfinite(q->z)) {
			//printf("Checking %lf\n", q->z);
			//for(int i = 0; i < n; i++)
			//	printf("x%d = %lf\n", i, q->x[i]);
			if(integer(q)) {
				bound(q, h, zp, x);
			} else if(branch(q, *zp)) {
			//	printf("h: %d, xh: %lf\n", q->h, q->xh);
				q->next = *h;
				*h = q;
				return;
			}
	}
	free_node(q);
}

double intopt(int m, int n, double **a, double *b, double *c, double *x) {
	node_t *p = inital_node(m, n, a, b, c);
	node_t *h = p;
	double z = -INFINITY;
	//printf("Begin intopt\n");
	p->z = simplex(p->m, p->n, p->a, p->b, p->c, p->x, 0);
	if(integer(p) || !isfinite(p->z)) {
		z = p->z;
		if(integer(p)) {
			for(int i = 0; i < p->n; i++) {
				x[i] = p->x[i];
			}
		}
		free_node(p);
		return z;
	}
	branch(p, z);
	while(h != NULL) {
		p = h;
		h = p->next;
		succ(p, &h, m, n, a, b, c, p->h, 1, floor(p->xh), &z, x);
		succ(p, &h, m, n, a, b, c, p->h, -1, -ceil(p->xh), &z, x);
		free_node(p);
	}

	if(z == -INFINITY) return NAN;
	else return z;
}
