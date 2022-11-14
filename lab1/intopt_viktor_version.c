#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double** make_matrix(int m, int n);
void print_everything(int m, int n, double *c, double *b, double **matrix);

int main() {
    int n, m;
    double *c, *b;
    double **matrix;

    scanf("%d %d",&n, &m);
    
    c = calloc(n, sizeof(*c));
    for(int i = 0; i < n; i++) {
        scanf("%lf", &c[i]);
    }

    matrix = make_matrix(m,n);
    for(int i = 0; i <m; i++) {
        for(int j = 0; j < n; j++) {
            scanf("%lf",&matrix[j][i]);
        }
    }

    b = calloc(m, sizeof(*b));
    for(int i = 0; i < n; i++) {
        scanf("%lf", &b[i]);
    }


    print_everything(n, m, c, b, matrix);
    
    free(c);
    free(b);

    for(int i = 0; i < n; i++) {
        free(matrix[i]);
    }
    free(matrix);

}

double** make_matrix(int m, int n){
    double** a;
    int i;
    a = calloc(m, sizeof(double*));
    for (i = 0; i < m; i += 1)
        a[i] = calloc(n, sizeof(double));
    return a;
}

void print_everything(int m, int n, double *c, double *b, double **matrix) {
    printf("m: %d n: %d\n",m,n);
    printf("max z = ");
    for(int i = 0; i < n; i++) {
        printf("%10.3lf",c[i]);
    }
    printf("\n");
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            printf("%10.3lf",matrix[j][i]);
        }
        printf("\n");
    }
    for(int i = 0; i < n; i++) {
        printf("%10.3lf",b[i]);
    }
    printf("\n");
}