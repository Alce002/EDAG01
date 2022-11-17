#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    int n, m;
    double *c, *b;
    double *matrix;

    scanf("%d %d", &m, &n);

    c = calloc(n, sizeof(*c));
    matrix = calloc(n * m, sizeof(*matrix));
    b = calloc(m, sizeof(*b));

    for(int i=0; i < n; i++) {
        scanf("%lf", &c[i]);
    }

    for(int y=0; y < m; y++) {
        for(int x=0; x < n; x++) {
            scanf("%lf", &matrix[x + n * y]);
        }
    }

    for(int i=0; i < m; i++) {
        scanf("%lf", &b[i]);
    }

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
            printf("%10.3lf x%d+", matrix[x + y * n], x);
        }
        printf("%10.3lf x%d \u2264 %10.3lf\n", matrix[n - 1 + y * n], n-1, b[y]);
    }

    free(c);
    free(b);
    free(matrix);
    return 0;
}    
