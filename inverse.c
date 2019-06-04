#include <stdio.h>
#include <math.h>
#define N 2

//Calculate the cofactor of mat[x][y] in temp[][]
void coFactor(double mat[N][N], double temp[N][N], int p, int q, int n) {
    int i=0, j=0;
    for (int row=0; row<n; row++) {
        for (int col = 0; col<n; col++) {
           if (row!=p && col!=q) {
                temp[i][j++] = mat[row][col];
                if (j==n-1) {
                    j=0;
                    i++;
                }
           }
        }
    }
}

/* Recursive function for finding determinant of matrix.
   n is current dimension of mat[][]. */
double determinant(double mat[N][N], int n) {
    double D=0;
    if (n==1)
        return mat[0][0];
    double temp[N][N];
    int sign = 1;
    for (int f=0;f<n;f++) {
        coFactor(mat, temp, 0, f, n);
        D += sign * mat[0][f] * determinant(temp, n-1);
        sign = -sign;
    }
    return D;
}

void inverseMat(double mat[N][N], double temp[N][N]) {
    int i=0, j=0;
    double det = determinant(mat, N);
    double blanck[N][N];
    for (i=0;i<N;i++) {
        for (j=0;j<N;j++) {
            coFactor(mat, blanck, j, i, N);
            temp[i][j] = pow(-1, i+j) * determinant(blanck, N-1);
            temp[i][j] /= det;
        }
    }
}

/* function for displaying the matrix */
void display(double mat[N][N], int row, int col)
{
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
            printf("  %lf", mat[i][j]);
        printf("\n");
    }
}

int main()
{
    /*double mat[N][N] = {{2., 1., 3.},
                     {5., 6., 4.},
                     {9., 8., 7.}};
*/
    double mat2[N][N];
    mat2[0][0] = 1.000000007451644;
    mat2[0][1] = -0.000000007451644;
    mat2[1][1] = 1.000000007451644;
    mat2[1][0] = -0.000000007451644;
/*
    printf("Determinant of the matrix is : %lf \n",
            determinant(mat, N));
    display(mat, 3, 3);
    double inv[N][N];
    inverseMat(mat, inv);
    display(inv, N, N);
*/
    display(mat2, N, N);
    printf("Determinant of the matrix is : %lf \n",
            determinant(mat2, N));
    double inv2[N][N];
    inverseMat(mat2, inv2);
    display(inv2, N, N);

    return 0;
}
