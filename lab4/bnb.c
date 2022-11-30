#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "simplex.h"
#include "bnb.h"

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
	printf("Best z: %lf\n", *zp);
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
			printf("Checking %lf\n", q->z);
			for(int i = 0; i < n; i++)
				printf("x%d = %lf\n", i, q->x[i]);
			if(integer(q)) {
				bound(q, h, zp, x);
			} else if(branch(q, *zp)) {
				printf("h: %d, xh: %lf\n", q->h, q->xh);
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
	printf("Begin intopt\n");
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
