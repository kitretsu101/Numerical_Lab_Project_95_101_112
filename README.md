# Table of Contents

- [Solution of Linear Equations](#solution-of-linear-equations)
  - [Matrix Inversion](#matrix-inversion)
    - [Theory](#matrix-inversion-theory)
    - [Code](#matrix-inversion-code)
    - [Input](#matrix-inversion-input)
    - [Output](#matrix-inversion-output)
   
## Solution of Linear Equations
### Matrix Inversion
#### Matrix Inversion Theory
Given an n x n square matrix A, if there exists another matrix B such that AB = BA = I (where I is the identity matrix), then B is called the inverse matrix of A and is denoted by A<sup>-1</sup>. The general formula for finding the inverse matrix of A is,
A<sup>-1</sup> = $\frac{1}{det(A)}$ Adj(A)

det(A) is referred as the determinant of A. A determinant of a square matrix is a scaler valued function of the entries of that matrix. 

Minor

Given a row i and column j, a minor element of matrix A is the determinant of the smaller matrix found by excluding the specified row and column from the given matrix A.

Cofactor

Given a row i and column j, a cofactor element is the minor element multiplied by (-1)<sup>(i+j)</sup>.

Adjoint

An adjoint matrix of A is the transpose of the cofactor matrix of A. To get the adjoint matrix, we create the cofactor matrix of A by evaluating cofactor for each element a<sub>i,j</sub> and get the transpose of it.

To find the inverse of a matrix A, the matrix must be :
 - Square matrix
 - Non-singular. i.e. det(A) != 0. If the matrix is singular, then it is non-invertible.

#### Matrix Inversion Code
```cpp
#include <iostream>
#include <bits/stdc++.h>
#include <fstream>

using namespace std;

vector<vector<double>> get_cofactor_matrix(vector<vector<double>>& mat, int r, int c, int n)
{
    vector<vector<double>> cofac(n-1, vector<double>(n-1));
    int curr_r = 0;
    int curr_c = 0;
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            if((i == r || j == c))
                continue;
            else
            {
                cofac[curr_r][curr_c++] = mat[i][j];
                if(curr_c == n-1)
                {
                    curr_c = 0;
                    curr_r++;
                }
            }
        }
    }

    return cofac;
}

double get_determinant(vector<vector<double>>& mat, int n)
{
    if(n == 1) return mat[0][0];
    if(n == 2) return mat[0][0]*mat[1][1] - mat[1][0]*mat[0][1];

    double determinant = 0;
    int sign = 0;
    vector<vector<double>> cofac(n-1, vector<double>(n-1));

    for(int col = 0; col < n; col++)
    {
        cofac = get_cofactor_matrix(mat, 0, col, n);
        sign = ((col % 2) == 0) ? 1 : -1;
        determinant += sign * get_determinant(cofac, n-1) * mat[0][col];
    }
    return determinant;
}

vector<vector<double>> get_adjoint_matrix(vector<vector<double>> mat, int n)
{
    vector<vector<double>> cofac(n-1, vector<double>(n-1));
    vector<vector<double>> adj(n, vector<double>(n));
    int sign = 0;
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            cofac = get_cofactor_matrix(mat, i, j, n);
            sign = ((i+j)%2==0) ? 1 : -1;
            double det = get_determinant(cofac, n-1);
            adj[j][i] = sign * det;
        }
    }
    return adj;
}

int main()
{
    ifstream infile("input.txt");
    ofstream outfile("output.txt");
    int n;
    while(infile >> n)
    {
        stringstream ss;
        vector<vector<double>> mat(n, vector<double>(n));
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                infile >> mat[i][j];
            }
        }
        vector<vector<double>> adj = get_adjoint_matrix(mat, n);
        double det = get_determinant(mat, n);

        if(det == 0)
        {
            ss << "The matrix is singular. Not invertible" << endl;
            outfile << ss.str();
            cout << ss.str();
            continue;
        }

        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                adj[i][j] = adj[i][j] / (double) det;
            }
        }

        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                ss << adj[i][j] << " ";
            }
            ss << endl;
        }
        ss << endl;
        outfile << ss.str();
        cout << ss.str();
    }
    outfile.close();
    infile.close();
    //cout << get_determinant(mat, 3);
    return 0;
}
```

#### Matrix Inversion Input
```
3
1 2 4
4 5 9
7 8 9

3
1 2 3
4 5 6
7 8 9
```

#### Matrix Inversion Output
```
-1.8 0.933333 -0.133333 
1.8 -1.26667 0.466667 
-0.2 0.4 -0.2 

The matrix is singular. Not invertible
```

