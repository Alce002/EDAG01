#ifndef BRANCH_AND_BOUND_H
#define BRANCH_AND_BOUND_H

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

#endif // BRANCH_AND_BOUND_H
