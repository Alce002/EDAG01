#ifndef BRANCH_AND_BOUND_H
#define BRANCH_AND_BOUND_H

typedef struct {
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
} node_t;

#endif // BRANCH_AND_BOUND_H
