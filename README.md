# Table of Contents

- [Solution of Linear Equations](#solution-of-linear-equations)
  - [Matrix Inversion](#matrix-inversion)
    - [Theory](#matrix-inversion-theory)
    - [Code](#matrix-inversion-code)
    - [Input](#matrix-inversion-input)
    - [Output](#matrix-inversion-output)
  - [Gauss Elimination Method](#gauss-elimination-method)
    - [Theory](#gauss-elimination-theory)
    - [Code](#gauss-elimination-code)
    - [Input](#gauss-elimination-input)
    - [Output](#gauss-elimination-output)
  - [Gauss Jordan Elimination Method](#gauss-jordan-elimination-method)
    - [Theory](#gauss-jordan-theory)
    - [Code](#gauss-jordan-code)
    - [Input](#gauss-jordan-input)
    - [Output](#gauss-jordan-output)
  - [LU Decomposition Method](#lu-decomposition-method)
    - [Theory](#lu-decomposition-theory)
    - [Code](#lu-decomposition-code)
    - [Input](#lu-decomposition-input)
    - [Output](#lu-decomposition-output)
- [Solution of Non-Linear Equations](#solution-of-non-linear-equations)
  - [Bisection Method](#bisection-method)
    - [Theory](#bisection-theory)
    - [Code](#bisection-code)
    - [Input](#bisection-input)
    - [Output](#bisection-output)
  - [False Position Method](#false-position-method)
    - [Theory](#false-position-theory)
    - [Code](#false-position-code)
    - [Input](#false-position-input)
    - [Output](#false-position-output)
  - [Newton Raphson Method](#newton-raphson-method)
    - [Theory](#newton-raphson-theory)
    - [Code](#newton-raphson-code)
    - [Input](#newton-raphson-input)
    - [Output](#newton-raphson-output)
  - [Secant Method](#secant-method)
    - [Theory](#secant-theory)
    - [Code](#secant-code)
    - [Input](#secant-input)
    - [Output](#secant-output)
- [Newton's Interpolation](#newtons-interpolation)
  - [Newton's Forward Interpolation](#newtons-forward-interpolation)
    - [Theory](#newtons-interpolation-theory)
    - [Code](#newtons-interpolation-code)
    - [Input](#newtons-interpolation-input)
    - [Output](#newtons-interpolation-output)
  - [Newton's Backward Interpolation](#newtons-backward-interpolation)
    - [Theory](#newtons-backward-interpolation-theory)
    - [Code](#newtons-backward-interpolation-code)
    - [Input](#newtons-backward-interpolation-input)
    - [Output](#newtons-backward-interpolation-output)
  - [Newton's Divided Difference](#newtons-divided-difference)
    - [Theory](#newtons-divided-difference-theory)
    - [Code](#newtons-divided-difference-code)
    - [Input](#newtons-divided-difference-input)
    - [Output](#newtons-divided-difference-output)
- [Numerical Integration](#numerical-integration)
  - [Simpson Method](#simpson-method)
    - [Theory](#simpson-method-theory)
    - [Code](#simpson-method-code)
    - [Input](#simpson-method-input)
    - [Output](#simpson-method-output)
- [Numerical Differentiation](#numerical-differentiation)
  - [Differentiation with Forward Interpolation](#differentiation-with-forward-interpolation)
    - [Theory](#differentiation-with-forward-interpolation-theory)
    - [Code](#differentiation-with-forward-interpolation-code)
    - [Input](#differentiation-with-forward-interpolation-input)
    - [Output](#differentiation-with-forward-interpolation-output)
  - [Differentiation with Backward Interpolation](#differentiation-with-backward-interpolation)
    - [Theory](#differentiation-with-backward-interpolation-theory)
    - [Code](#differentiation-with-backward-interpolation-code)
    - [Input](#differentiation-with-backward-interpolation-input)
    - [Output](#differentiation-with-backward-interpolation-output)
- [Ordinary Differential Equation](#ordinary-differential-equation)
  - [RK Method](#rk-method)
    - [Theory](#rk-method-theory)
    - [Code](#rk-method-code)
    - [Input](#rk-method-input)
    - [Output](#rk-method-output)
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

### Gauss Elimination Method
#### Gauss Elimination Theory
The Gauss Elimination Method is a numerical technique used to solve a system of linear equations. The main idea is to eliminate variables step by step and convert the system into an upper triangular matrix, from which the solution can be easily obtained using backward substitution.

It works on the augmented matrix [A | b] formed from the system A x = b where,
A is an n × n matrix,
x is the unknown vector,
b is the constant vector.

To get the upper triangular matrix forward elimination is have to be done. To do this,
For each pivot row k = 1 to n-1:
factor = a[i][k] / a[k][k].

To update the rows:
a[i][j] = a[i][j] - factor * a[k][j],
b[i]    = b[i]    - factor * b[k]
where:
i = k+1 to n,
j = k to n.
This process transforms the matrix into an upper triangular form.

Once the matrix is upper triangular, the unknowns are solved using backward sunstitution.
x[n] = b[n] / a[n][n],
x[i] = ( b[i] - Σ(a[i][j] * x[j]) ) / a[i][i]
where:
j = i+1 to n,
i = n-1, n-2, ..., 1.
The x vector is the solution of the system.

#### Gauss Elimination Code
```cpp
#include <bits/stdc++.h>
using namespace std;

int type(vector<vector<double>>&a,int n)
{
    int r=0,c=0;
    for(int i=0;i<n;i++){
        bool nz=false;
        for(int j=0;j<n;j++)
            if(fabs(a[i][j])>1e-9) nz=true;
        if(nz) r++;
        else if(fabs(a[i][n])>1e-9) return 0;
    }
    for(int j=0;j<n;j++){
        for(int i=0;i<n;i++)
            if(fabs(a[i][j])>1e-9){
                c++;
                break;
            }
    }
    if(r==c) return 1;
    return 2;
}

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    if(!fin){
        cout<<"input.txt not found"<<endl;
        return 0;
    }

    while(true)
    {
        int n;
        fin>>n;
        vector<vector<double>>a(n,vector<double>(n+1));
        for(int i=0;i<n;i++)
            for(int j=0;j<=n;j++)
                fin>>a[i][j];

        for(int k=0;k<n;k++){
            int mx=k;
            for(int i=k+1;i<n;i++)
                if(fabs(a[i][k])>fabs(a[mx][k])) mx=i;
            swap(a[k],a[mx]);

            if(fabs(a[k][k])<1e-9) continue;

            for(int i=k+1;i<n;i++){
                double f=a[i][k]/a[k][k];
                for(int j=k;j<=n;j++)
                    a[i][j]-=f*a[k][j];
            }
        }

        int t=type(a,n);

        fout<<"System size: "<<n<<endl;

        if(t==0){
            fout<<"Solution type: No solution"<<endl;
        }
        else if(t==2){
            fout<<"Solution type: Infinite solutions"<<endl;
        }
        else{
            fout<<"Solution type: Unique solution"<<endl;
            vector<double>x(n);
            for(int i=n-1;i>=0;i--){
                x[i]=a[i][n];
                for(int j=i+1;j<n;j++)
                    x[i]-=a[i][j]*x[j];
                x[i]/=a[i][i];
            }
            for(int i=0;i<n;i++)
                fout<<"x"<<i+1<<" = "<<x[i]<<endl;
        }

        char ch;
        fin>>ch;
        if(ch!='y' && ch!='Y') break;
        fout<<endl;
    }

    fin.close();
    fout.close();
    return 0;
}
```

#### Gauss Elimination Input
```
3
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3
y
2
1 1 2
2 2 5
y
2
1 1 2
2 2 4
n
```

#### Gauss Elimination Output
```
System size: 3
Solution type: Unique solution
x1 = 2
x2 = 3
x3 = -1

System size: 2
Solution type: No solution

System size: 2
Solution type: Infinite solutions
```

### Gauss Jordan Elimination Method
#### Gauss Jordan Theory
The Gauss–Jordan Method is an extension of Gaussian elimination. Instead of stopping at an upper triangular matrix, it continues elimination to convert the matrix into Reduced Row Echelon Form (RREF). This allows the solution to be obtained directly, without backward substitution.
Given:
A x = b where,
A is an n × n matrix,
x is the unknown vector,
b is the constant vector.

The augmented matrix [A | b] is transformed into:
[I | x],
where I is the identity matrix.

For each pivot element a[k][k]:
1. Normalizing the pivot row
   Row[k] = Row[k] / a[k][k].

2. Eliminating all other elements in the pivot column
a[i][j] = a[i][j] - a[i][k] * a[k][j]
for:
i ≠ k

After reaching RREF, the vector x would the solution of the system.

#### Gauss Jordan Code
```cpp
#include <bits/stdc++.h>
using namespace std;

int type(vector<vector<double>>&a,int n)
{
    int r=0,c=0;
    for(int i=0;i<n;i++){
        bool nz=false;
        for(int j=0;j<n;j++)
            if(fabs(a[i][j])>1e-9) nz=true;
        if(nz) r++;
        else if(fabs(a[i][n])>1e-9) return 0;
    }
    for(int j=0;j<n;j++){
        for(int i=0;i<n;i++)
            if(fabs(a[i][j])>1e-9){
                c++;
                break;
            }
    }
    if(r==c) return 1;
    return 2;
}

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    if(!fin){
        cout<<"input.txt not found"<<endl;
        return 0;
    }

    while(true)
    {
        int n;
        fin>>n;

        vector<vector<double>>a(n,vector<double>(n+1));
        for(int i=0;i<n;i++)
            for(int j=0;j<=n;j++)
                fin>>a[i][j];

        for(int k=0;k<n;k++){
            int mx=k;
            for(int i=k+1;i<n;i++)
                if(fabs(a[i][k])>fabs(a[mx][k])) mx=i;
            swap(a[k],a[mx]);

            if(fabs(a[k][k])<1e-9) continue;

            double div=a[k][k];
            for(int j=0;j<=n;j++)
                a[k][j]/=div;

            for(int i=0;i<n;i++){
                if(i==k) continue;
                double f=a[i][k];
                for(int j=0;j<=n;j++)
                    a[i][j]-=f*a[k][j];
            }
        }

        int t=type(a,n);

        fout<<"System size: "<<n<<endl;

        if(t==0){
            fout<<"Solution type: No solution"<<endl;
        }
        else if(t==2){
            fout<<"Solution type: Infinite solutions"<<endl;
        }
        else{
            fout<<"Solution type: Unique solution"<<endl;
            for(int i=0;i<n;i++)
                fout<<"x"<<i+1<<" = "<<a[i][n]<<endl;
        }

        char ch;
        fin>>ch;
        if(ch!='y' && ch!='Y') break;
        fout<<endl;
    }

    fin.close();
    fout.close();
    return 0;
}
```

#### Gauss Jordan Input
```
3
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3
y
2
1 1 2
2 2 5
y
2
1 1 2
2 2 4
n
```

#### Gauss Jordan Output
```
System size: 3
Solution type: Unique solution
x1 = 2
x2 = 3
x3 = -1

System size: 2
Solution type: No solution

System size: 2
Solution type: Infinite solutions
```

### LU Decomposition Method
#### LU Decomposition Theory
The LU Factorization Method is used to solve systems of linear equations by decomposing the coefficient matrix into a lower triangular matrix (L) and an upper triangular matrix (U).
It is especially efficient when solving the same system with multiple right-hand sides.

Given: A augmented matrix A x = b where,
A is an n × n matrix,
x is the unknown vector,
b is the constant vector.

Now, we have to factorize it such as: 
A = LU
where,
L is a lower triangular matrix, 
U is an upper triangular matrix.

Steps:
Step 1: Forward Substitution:
To do this we have to solve,
L y = b using, 
y[i] = b[i] - Σ(L[i][j] * y[j]).
where, j = 1 to i-1.

Step 2: Backward Substitution
To do this we have to solve,
U x = y using, 
x[n] = y[n] / U[n][n], 
x[i] = ( y[i] - Σ(U[i][j] * x[j]) ) / U[i][i].
where, j = i+1 to n

#### LU Decomposition Code
```cpp
#include <bits/stdc++.h>
using namespace std;

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int n;
    fin >> n;

    vector<vector<double>> A(n, vector<double>(n));
    vector<double> B(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            fin >> A[i][j];
        fin >> B[i];
    }

    vector<vector<double>> L(n, vector<double>(n, 0));
    vector<vector<double>> U(n, vector<double>(n, 0));

    for (int i = 0; i < n; i++) {
        for (int k = i; k < n; k++) {
            double sum = 0;
            for (int j = 0; j < i; j++)
                sum += L[i][j] * U[j][k];
            U[i][k] = A[i][k] - sum;
        }

        L[i][i] = 1;

        if (fabs(U[i][i]) < 1e-9) {
            fout << "The system has no unique solution.\n";
            return 0;
        }

        for (int k = i + 1; k < n; k++) {
            double sum = 0;
            for (int j = 0; j < i; j++)
                sum += L[k][j] * U[j][i];
            L[k][i] = (A[k][i] - sum) / U[i][i];
        }
    }

    vector<double> Y(n);
    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < i; j++)
            sum += L[i][j] * Y[j];
        Y[i] = B[i] - sum;
    }

    vector<double> X(n);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < n; j++)
            sum += U[i][j] * X[j];
        X[i] = (Y[i] - sum) / U[i][i];
    }

    fout << fixed << setprecision(3);

    fout << "Lower Triangular Matrix (L):\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            fout << setw(8) << L[i][j];
        fout << endl;
    }

    fout << "\nUpper Triangular Matrix (U):\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            fout << setw(8) << U[i][j];
        fout << endl;
    }

    fout << "\nSolution:\n";
    for (int i = 0; i < n; i++)
        fout << "x" << i + 1 << " = " << X[i] << endl;

    fout << "\nThe system has unique solution.\n";

    fin.close();
    fout.close();
    return 0;
}
```

#### LU Decomposition Input
```
3
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3
```

#### LU Decomposition Output
```
Lower Triangular Matrix (L):
   1.000   0.000   0.000
  -1.500   1.000   0.000
  -1.000   4.000   1.000

Upper Triangular Matrix (U):
   2.000   1.000  -1.000
   0.000   0.500   0.500
   0.000   0.000  -1.000

Solution:
x1 = 2.000
x2 = 3.000
x3 = -1.000

The system has unique solution.
```

## Solution of Non-Linear Equations
### Bisection Method
#### Bisection Theory
The Bisection Method is a simple and reliable numerical technique used to find a root of a nonlinear equation. The method works by repeatedly dividing an interval into two halves and selecting the subinterval that contains the root.It is based on the Intermediate Value Theorem.

Given a nonlinear equation: f(x) = 0.

Choose two initial values a and b such that:
f(a) * f(b) < 0.
This ensures that a root lies between a and b.
Steps:

1. Computing the midpoint:c = (a + b) / 2

2. Evaluating f(c),
   If f(a) * f(c) < 0, set b = c
   If f(c) * f(b) < 0, set a = c

This process is repeated until convergence where convergence condition is,
|b - a| < ε or |f(c)| < ε, 
where ε is the tolerance.

We can get all the roots of the plynomial by doing this process repeteadly in an interval.

#### Bisection Code
```cpp
#include <bits/stdc++.h>
using namespace std;

vector<double> a;
int n;

double f(double x) {
    double sum = 0;
    for (int i = 0; i <= n; i++) {
        sum += a[i] * pow(x, n - i);
    }
    return sum;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    fin >> n;
    a.resize(n + 1);

    for (int i = 0; i <= n; i++)
        fin >> a[i];

    double step, eps;
    fin >> step >> eps;

    double an = a[0];
    double an1 = (n >= 1) ? a[1] : 0;
    double an2 = (n >= 2) ? a[2] : 0;

    double xmax = sqrt((an1 * an1 - 2 * an * an2) / an);

    fout << fixed << setprecision(4);
    fout << "Search Interval: [" << -xmax << ", " << xmax << "]\n\n";

    vector<pair<double, double>> brackets;

    for (double x = -xmax; x < xmax; x += step) {
        if (f(x) * f(x + step) <= 0) {
            brackets.push_back({x, x + step});
        }
    }

    int rootCount = 1;

    for (auto br : brackets) {
        double left = br.first;
        double right = br.second;
        double mid;
        int iter = 0;

        while (true) {
            mid = (left + right) / 2.0;
            iter++;

            if (fabs(f(mid)) < eps)
                break;

            if (f(left) * f(mid) < 0)
                right = mid;
            else
                left = mid;
        }

        fout << "Root " << rootCount++
             << " = " << mid
             << "   Iterations = " << iter << endl;

        if (rootCount > n)
            break;
    }

    fin.close();
    fout.close();
    return 0;
}
```

#### Bisection Input
```
3
1 -6 11 -6
0.5
0.0001
```

#### Bisection Output
```
Search Interval: [-3.7417, 3.7417]

Root 1 = 1.0000   Iterations = 10
Root 2 = 2.0000   Iterations = 10
Root 3 = 3.0000   Iterations = 10
```

### False Position Method
#### False Position Theory
The False Position Method is an improvement over the bisection method. Instead of choosing the midpoint like Bisection method, it uses a linear interpolation between the two endpoints to approximate the root.
Given:
f(x) = 0 where f(x) is a polynomial of any degree.

Let, the initial values are a and b such that:
f(a) * f(b) < 0.

Formula:
Formula used for this method is:
c = (a * f(b) - b * f(a)) / (f(b) - f(a)).

Steps:

1. Computing c using the above formula.

2. Evaluating f(c),
   If f(a) * f(c) < 0, setting b = c,
   If f(c) * f(b) < 0, setting a = c.
This process if repeated until convergence.

Convergence Condition:
|f(c)| < ε or |b - a| < ε where ε is tolerance.

Like Bisection method, we can get all the roots of the polynomial by repeating the process inside a specific interval.

#### False Position Code
```cpp
#include <bits/stdc++.h>
using namespace std;

double f(double x, vector<double>& coef) {
    double res = 0;
    int n = coef.size();
    for (int i = 0; i < n; i++)
        res = res * x + coef[i];
    return res;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int n;
    fin >> n;

    vector<double> coef(n + 1);
    for (int i = 0; i <= n; i++)
        fin >> coef[i];

    double step, eps;
    fin >> step >> eps;

    double bound = sqrt((coef[1] * coef[1] - 2 * coef[0] * coef[2]) / coef[0]);
    fout << fixed << setprecision(4);
    fout << "Search Interval: [-" << bound << ", " << bound << "]\n\n";

    vector<pair<double, double>> intervals;

    for (double x = -bound; x < bound; x += step) {
        if (f(x, coef) * f(x + step, coef) < 0)
            intervals.push_back({x, x + step});
    }

    int rootCount = 1;
    for (auto it : intervals) {
        double a = it.first, b = it.second;
        double fa = f(a, coef), fb = f(b, coef);
        double c, fc;
        int iter = 0;

        do {
            c = (a * fb - b * fa) / (fb - fa);
            fc = f(c, coef);

            if (fa * fc < 0) {
                b = c;
                fb = fc;
            } else {
                a = c;
                fa = fc;
            }

            iter++;
        } while (fabs(fc) > eps);

        fout << "Root " << rootCount++ << " = " << c
             << "   Iterations = " << iter << endl;
    }

    if (intervals.empty())
        fout << "No real root found.\n";

    fin.close();
    fout.close();
    return 0;
}
```

#### False Position Input
```
3
1 -6 11 -6
0.5
0.0001
```

#### False Position Output
```
Search Interval: [-3.7417, 3.7417]

Root 1 = 1.0000   Iterations = 8
Root 2 = 1.9999   Iterations = 2
Root 3 = 3.0000   Iterations = 8
```

### Newton Raphson Method
#### Newton Raphson Theory
The Newton–Raphson Method is an iterative numerical technique used to find roots of nonlinear equations. It starts from an initial guess and improves the approximation using the derivative of the function.
The method is based on the tangent line approximation.
Given:
f(x) = 0 where f(x) is a polynomial.

The iterative formula for this method is:
x[n+1] = x[n] - f(x[n]) / f'(x[n]).

Convergence Condition:
The iteration is repeated until:
|x[n+1] - x[n]| < ε
where ε is the permissible error / tolerance.

If the process is repeated in a specific interval for a specific step size, all the roots can be determined for the corresponding polynomial.

#### Newton Raphson Code
```cpp
#include <bits/stdc++.h>
using namespace std;

double poly(const vector<double>&c,double x)
{
    double s=0;
    int n=c.size();
    for(int i=0;i<n;i++)s+=c[i]*pow(x,n-1-i);
    return s;
}

vector<double> derivative(const vector<double>&c)
{
    int n=c.size();
    vector<double>d;
    for(int i=0;i<n-1;i++)
        d.push_back(c[i]*(n-1-i));
    return d;
}

void printPoly(const vector<double>& c,ostream& out)
{
    int n=c.size();
    for(int i=0;i<n;i++)
    {
        if(c[i]==0)continue;
        if(i>0 && c[i]>0)out<<" + ";
        if(c[i]<0)out<<" - ";
        if(abs(c[i])!=1 || i==n-1)out<<abs(c[i]);
        if(i<n-1)
        {
            out<<"x";
            if(n-1-i>1)out<<"^"<<n-1-i;
        }
    }
    out<<endl;
}

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    if(!fin)
    {
        cout<<"Input file not found"<<endl;
        return 0;
    }

    int degree;
    fin>>degree;

    vector<double>coef(degree+1);
    for(int i=0;i<=degree;i++)fin>>coef[i];

    double a,b,step,tol;
    fin>>a>>b>>step>>tol;

    vector<double>dcoef=derivative(coef);
    vector<double>roots;

    for(double x=a;x<=b;x+=step)
    {
        double x0=x;

        for(int i=0;i<100;i++)
        {
            double f=poly(coef,x0);
            double df=poly(dcoef,x0);

            if(fabs(df)<1e-8)break;

            double x1=x0-f/df;

            if(fabs(x1-x0)<tol)
            {
                if(x1>=a && x1<=b)
                {
                    bool exists=false;
                    for(double r:roots)
                        if(fabs(r-x1)<tol) exists=true;

                    if(!exists) roots.push_back(x1);
                }
                break;
            }
            x0=x1;
        }
    }

    cout<<"Degree of Polynomial: "<<degree<<endl;
    fout<<"Degree of Polynomial: "<<degree<<endl;

    cout<<"Polynomial: ";
    fout<<"Polynomial: ";
    printPoly(coef,cout);
    printPoly(coef,fout);

    cout<<"Derivative: ";
    fout<<"Derivative: ";
    printPoly(dcoef,cout);
    printPoly(dcoef,fout);

    if(roots.empty())
    {
        cout<<"No roots found in given range"<<endl;
        fout<<"No roots found in given range"<<endl;
    }
    else
    {
        cout<<"Roots found:"<<endl;
        fout<<"Roots found:"<<endl;
        for(int i=0;i<roots.size();i++)
        {
            cout<<"Root "<<i+1<<": "<<roots[i]<<endl;
            fout<<"Root "<<i+1<<": "<<roots[i]<<endl;
        }
    }

    fin.close();
    fout.close();
    return 0;
}
```

#### Newton Raphson Input
```
4
1 0 -5 0 4
-5 5
0.5
0.0001
```

#### Newton Raphson Output
```
Degree of Polynomial: 4
Polynomial: x^4 - 5x^2 + 4
Derivative: 4x^3 - 10x
Roots found:
Root 1: -2
Root 2: -1
Root 3: 1
Root 4: 2
```

### Secant Method
#### Secant Theory
The Secant Method is an iterative numerical technique used to find the root of a non-linear equation. It is similar to the Newton–Raphson method, but does not require the derivative of the function.
Instead of using the slope of a tangent, it approximates the slope using a secant line passing through two nearby points on the function.

Given a nonlinear equation:
f(x) = 0 where f(x) is a polynomial.

The goal is to determine a root x* such that:
f(x*) = 0.

Steps:
1. Choosing two initial approximations x₀ and x₁ such that:
f(x₀) ≠ f(x₁).
2. Computing the next approximation using the secant formula:
   x[n+1] = x[n] - f(x[n]) * (x[n] - x[n-1])
                 / (f(x[n]) - f(x[n-1])).
3. Repeating the iteration until the convergence condition is satisfied.

Convergence Condition:
The iteration is stopped when:
|x[n+1] - x[n]| < ε..
where ε is the allowable error / tolerance.

If the process is repeated in a specific interval for a specific step size, all the roots can be determined for the corresponding polynomial.

#### Secant Code
```cpp
#include <bits/stdc++.h>
using namespace std;

double poly(const vector<double>&c,double x){
    double r=0;
    int n=c.size();
    for(int i=0;i<n;i++) r+=c[i]*pow(x,n-1-i);
    return r;
}

void printPoly(const vector<double>&c,ostream &out){
    int n=c.size();
    for(int i=0;i<n;i++){
        if(c[i]==0) continue;
        if(i>0&&c[i]>0) out<<" + ";
        if(c[i]<0) out<<" - ";
        if(abs(c[i])!=1||i==n-1) out<<abs(c[i]);
        if(i<n-1){
            out<<"x";
            if(n-1-i>1) out<<"^"<<n-1-i;
        }
    }
    out<<endl;
}

int main(){
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    if(!fin) return 0;

    int degree;
    fin>>degree;
    vector<double>coef(degree+1);
    for(int i=0;i<=degree;i++) fin>>coef[i];

    double minX,maxX,step,tol;
    fin>>minX>>maxX>>step>>tol;

    fout<<"Degree of Polynomial: "<<degree<<endl;
    fout<<"Polynomial: ";
    printPoly(coef,fout);

    vector<double>roots;
    double x=minX;

    while(x<maxX){
        double x0=x,x1=x+step;
        if(x1>maxX) x1=maxX;

        double f0=poly(coef,x0),f1=poly(coef,x1);

        if(f0*f1<=0){
            double curr=x0,prev=x1;
            double fcurr=poly(coef,curr),fprev=poly(coef,prev);
            double next;

            do{
                if(fcurr==fprev) break;
                next=curr-fcurr*(curr-prev)/(fcurr-fprev);
                prev=curr;
                fprev=fcurr;
                curr=next;
                fcurr=poly(coef,curr);
            }while(fabs(fcurr)>tol);

            bool exists=false;
            for(double r:roots) if(fabs(r-curr)<tol) exists=true;
            if(!exists) roots.push_back(curr);
        }

        x+=step;
    }

    if(roots.empty()) fout<<"No roots found in given range"<<endl;
    else{
        fout<<"Roots found:"<<endl;
        for(int i=0;i<roots.size();i++) fout<<"Root "<<i+1<<": "<<roots[i]<<endl;
    }

    fin.close();
    fout.close();
    return 0;
}
```

#### Secant Input
```
4
1 0 -5 0 4
-3.5 3.5
0.5
0.0001
```

#### Secant Output
```
Degree of Polynomial: 4
Polynomial: x^4 - 5x^2 + 4
Roots found:
Root 1: -2
Root 2: -1
Root 3: 1
Root 4: 2
```

## Newton's Interpolation
Interpolation is referred to as the technique to find the value of a function at a point that lies between two unknown points. It allows us a to construct a function from the given known data points and find out the values of the function at the intermediate data points. There are several methods of interpolation :
- Newton's Forward Interpolation
- Newton's Backward Interpolation
- Newton's Divided difference method
### Newton's Forward Interpolation
#### Newton's Forward Interpolation Theory
Newton's Forward Interpolation method is used when the given data points of the independent variable x are of equidistance and the value to be interpolated is near at the beginning of the table. 

Let the given known data points of the independent variable x are x<sub>0</sub>, x<sub>1</sub>, x<sub>2</sub>, ......, x<sub>n</sub>

Now, x<sub>1</sub> - x<sub>0</sub> = x<sub>2</sub> - x<sub>1</sub> = x<sub>n</sub> - x<sub>n-1</sub> = h

Let, u = $\frac{x-x0}{h}$

We also need to create a forward difference table. The forward difference is calculated as follows,

y<sub>1</sub> - y<sub>0</sub> = Δy<sub>0</sub>

Δy<sub>1</sub> - Δy<sub>0</sub> = Δ<sup>2</sup>y<sub>0</sub>

Δ<sup>2</sup>y<sub>1</sub> - Δ<sup>2</sup>y<sub>0</sub> = Δ<sup>3</sup>y<sub>0</sub>

.....

Newton's forward interpolation formula is given as,

y = y<sub>0</sub> + uΔy<sub>0</sub> + $\frac{u(u-1)}{2!}$ Δ<sup>2</sup>y<sub>0</sub> + $\frac{u(u-1)(u-2)}{3!}$ Δ<sup>3</sup>y<sub>0</sub> + ..... + $\frac{u(u-1)(u-2)...(u-n+1)}{n!}$ Δ<sup>n</sup>y<sub>0</sub> 

#### Newton's Forward Interpolation Code
```cpp
#include <iostream>
#include <bits/stdc++.h>
#include <fstream>

using namespace std;

int main()
{
    ifstream infile("input.txt");
    ofstream outfile("output.txt");


    int n;
    //cin >> n;
    while(infile >> n) {
    vector<vector<double>> y(n, vector<double>(n, 0));
    vector<double> x(n, 0);

    for(int i = 0; i < n; i++) {
        infile >> x[i] >> y[i][0];
    }

    for(int j = 1; j < n; j++) {
        for(int i = 0; i < n-j; i++) {
            y[i][j] = y[i+1][j-1] - y[i][j-1];
        }
    }

    double interpolate_val;
    infile >> interpolate_val;

    double u = (interpolate_val - x[0]) / (x[1] - x[0]);
    double sum = y[0][0];
    double u_sum = 1;
    for(int i = 0; i < n-1; i++) {
        u_sum *= (u - i) / (i+1);
        sum += u_sum * y[0][i+1];
    }

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            cout << y[i][j] << " ";
        }
        cout << endl;
    }

    stringstream os;
    os << "value of y at " << interpolate_val << " is : " << sum << endl;
    cout << os.str();
    outfile << os.str();
    }
}
```
#### Newton's Forward Interpolation Input
```
4
3 180
5 150
7 120
9 90
4

4
45 0.7071
50 0.766
55 0.8192
60 0.866
52
```

#### Newton's Forward Interpolation Output
```
value of y at 4 is : 165
value of y at 52 is : 0.788003
```
### Newton's Backward Interpolation
#### Newton's Backward Interpolation Theory
Newton's Backward Interpolation method is used when the given data points of the independent variable x are of equidistance and the value to be interpolated is near at the end of the table. 

Let the given known data points of the independent variable x are x<sub>0</sub>, x<sub>1</sub>, x<sub>2</sub>, ......, x<sub>n</sub>

Now, x<sub>1</sub> - x<sub>0</sub> = x<sub>2</sub> - x<sub>1</sub> = x<sub>n</sub> - x<sub>n-1</sub> = h

Let, v = $\frac{x-xn}{h}$

We also need to create a backward difference table. The backward difference is calculated as follows,

y<sub>n</sub> - y<sub>n-1</sub> = Δy<sub>n</sub>

Δy<sub>n</sub> - Δy<sub>n-1</sub> = Δ<sup>2</sup>y<sub>n</sub>

Δ<sup>2</sup>y<sub>n</sub> - Δ<sup>2</sup>y<sub>n-1</sub> = Δ<sup>3</sup>y<sub>n</sub>

.....

Newton's backward interpolation formula is given as,

y = y<sub>n</sub> + vΔy<sub>n</sub> + $\frac{v(v+1)}{2!}$ Δ<sup>2</sup>y<sub>n</sub> + $\frac{v(v+1)(v+2)}{3!}$ Δ<sup>3</sup>y<sub>n</sub> + ..... + $\frac{v(v+1)(v+2)...(v+n-1)}{n!}$ Δ<sup>n</sup>y<sub>n</sub> 

#### Newton's Backward Interpolation Code
```cpp
#include <iostream>
#include <bits/stdc++.h>
#include <fstream>

using namespace std;

int main()
{
    ifstream infile("input.txt");
    ofstream outfile("output.txt");


    int n;
    //cin >> n;
    while(infile >> n)
    {
        vector<vector<double>> y(n, vector<double>(n, 0));
        vector<double> x(n, 0);

        for(int i = 0; i < n; i++)
        {
            infile >> x[i] >> y[i][0];
        }

        for(int j = 1; j < n; j++)
        {
            for(int i = 0; i < n-j; i++)
            {
                y[i][j] = y[i+1][j-1] - y[i][j-1];
            }
        }

        double interpolate_val;
        infile >> interpolate_val;

        double u = (interpolate_val - x[n-1]) / (x[1] - x[0]);
        double sum = y[n-1][0];
        double u_sum = 1;
        for(int i = 0; i < n-1; i++)
        {
            u_sum *= (u + i) / (i+1);
            sum += u_sum * y[n-i-2][i+1];
        }

        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                cout << y[i][j] << " ";
            }
            cout << endl;
        }

        stringstream os;
        os << "value of y at " << interpolate_val << " is : " << sum << endl;
        cout << os.str();
        outfile << os.str();
    }
}
```
#### Newton's Backward Interpolation Input
```
5
1891 46
1901 66
1911 81
1921 93
1931 101
1925

5
24 28.06
28 30.19
32 32.75
36 34.94
40 40
33
```

#### Newton's Backward Interpolation Output
```
value of y at 1925 is : 96.8368
value of y at 33 is : 33.2747
```


### Newton's Divided Difference
#### Newton's Divided Difference Theory
Newton's divided difference method is used when the given data points of the independent variable x are not of equidistance.

Let the given known data points of the independent variable x are x<sub>0</sub>, x<sub>1</sub>, x<sub>2</sub>, ......, x<sub>n</sub>

Now, We also need to create a forward difference table for unequal intervals. The forward difference table is calculated as follows,

f[x<sub>i</sub>, x<sub>j</sub>] = $\frac{f(xi)-f(xj)}{xi-xj}$ 

f[x<sub>i</sub>, x<sub>j</sub>, x<sub>k</sub>] = $\frac{f[xi, xj]-f[xj, xk]}{xi-xk}$

....

Newton's divided difference formula is given as,

f(x<sub>n</sub>) = f(x<sub>0</sub>) + (x-x<sub>0</sub>)f[x<sub>1</sub>, x<sub>0</sub>] + (x-x<sub>0</sub>)(x-x<sub>1</sub>)f[x<sub>2</sub>, x<sub>1</sub>, x<sub>0</sub>] + .... + (x-x<sub>0</sub>)(x-x<sub>1</sub>)...(x-x<sub>n-1</sub>)f[x<sub>n</sub>, x<sub>n-1</sub>, x<sub>n-2</sub>, ......, x<sub>1</sub>, x<sub>0</sub>]

#### Newton's Divided Difference Code
```cpp
#include <iostream>
#include <bits/stdc++.h>
#include <fstream>

using namespace std;

int main()
{
    ifstream infile("input.txt");
    ofstream outfile("output.txt");


    int n;
    //cin >> n;
    while(infile >> n)
    {
        vector<vector<double>> y(n, vector<double>(n, 0));
        vector<double> x(n, 0);

        for(int i = 0; i < n; i++)
        {
            infile >> x[i] >> y[i][0];
        }

        for(int j = 1; j < n; j++)
        {
            for(int i = 0; i < n-j; i++)
            {
                y[i][j] = (y[i+1][j-1] - y[i][j-1]) / (double)(x[j+i] - x[i]);
            }
        }

        double interpolate_val;
        infile >> interpolate_val;

        //double u = (interpolate_val - x[n-1]) / (x[1] - x[0]);
        double sum = y[0][0];
        double u_sum = 1;
        for(int i = 0; i < n-1; i++)
        {
            u_sum *= (interpolate_val - x[i]);
            sum += u_sum * y[0][i+1];
        }

        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                cout << y[i][j] << " ";
            }
            cout << endl;
        }

        stringstream os;
        os << "value of y at " << interpolate_val << " is : " << sum << endl;
        cout << os.str();
        outfile << os.str();
    }
}
```
#### Newton's Divided Difference Input
```
4
1 0
4 1.386294
6 1.79175
5 1.609438
2
```

#### Newton's Divided Difference Output
```
value of y at 2 is : 0.628762
```

## Numerical Integration
### Simpson Method
#### Simpson Method Theory
Simpson’s Methods are numerical integration techniques used to approximate the value of a definite integral. They work by replacing the function with polynomial curves over small subintervals, which gives higher accuracy than the trapezoidal rule for smooth functions.

Simpson’s 1/3 Rule:

Simpson’s 1/3 Rule approximates the integrand using a second-degree polynomial over two consecutive subintervals. The number of subintervals n must be even.

Formula: 
∫[a to b] f(x) dx ≈ (h / 3) [
    f(x0)
  + 4 (f(x1) + f(x3) + ... + f(xn-1))
  + 2 (f(x2) + f(x4) + ... + f(xn-2))
  + f(xn)
]

Simpson’s 3/8 Rule:

Simpson’s 3/8 Rule approximates the integrand using a third-degree (cubic) polynomial over three consecutive subintervals. The number of subintervals n must be a multiple of 3.

Formula: 
∫[a to b] f(x) dx ≈ (3h / 8) [
    f(x0)
  + 3 (f(x1) + f(x2) + f(x4) + f(x5) + ...)
  + 2 (f(x3) + f(x6) + ...)
  + f(xn)
]

#### Simpson Method Code
```cpp
#include <bits/stdc++.h>
using namespace std;

double f(double x)
{
    double sin_x = sin(x);
    return pow(sin_x, 5) + 4.0 * pow(sin_x, 4) + 1.0;
}

double simp_one_third(double a, double b, int N)
{
    if (N % 2 != 0)
    {
        return NAN;
    }

    double h = (b - a) / N;
    double integral = f(a) + f(b);

    for (int i = 1; i < N; ++i)
    {
        double x = a + i * h;
        if (i % 2 == 0)
        {
            integral += 2.0 * f(x);
        }
        else
        {
            integral += 4.0 * f(x);
        }
    }

    return (h / 3.0) * integral;
}

double simpsons_three_eighths(double a, double b, int N)
{
    if (N % 3 != 0)
    {
        return NAN;
    }

    double h = (b - a) / N;
    double integral = f(a) + f(b);

    for (int i = 1; i < N; ++i)
    {
        double x = a + i * h;
        if (i % 3 == 0)
        {
            integral += 2.0 * f(x);
        }
        else
        {
            integral += 3.0 * f(x);
        }
    }

    return (3.0 * h / 8.0) * integral;
}

void run_simpsons_rules()
{
    double a, b;
    int N;

    ifstream inputFile("input.txt");
    if (!inputFile.is_open())
    {
        cerr << "Error: Could not open input.txt" << endl;
        return;
    }

    if (!(inputFile >> a >> b >> N))
    {
        cerr << "Error: Failed to read a, b, and N from input.txt. Check the file format." << endl;
        inputFile.close();
        return;
    }
    inputFile.close();

    double integral_1_3 = simp_one_third(a, b, N);
    double integral_3_8 = simpsons_three_eighths(a, b, N);

    ofstream outputFile("output.txt");
    if (!outputFile.is_open())
    {
        cerr << "Error: Could not open output.txt for writing." << endl;
        return;
    }

    outputFile << fixed << setprecision(10);

    outputFile << "--- Numerical Integration Results (Simpson's Rules) ---" << endl;
    outputFile << "Function: f(x) = sin^5(x) + 4*sin^4(x) + 1" << endl;
    outputFile << "Limits: a = " << a << ", b = " << b << endl;
    outputFile << "Number of Intervals (N) = " << N << endl;
    outputFile << "Interval width (h) = " << (b - a) / N << endl;
    outputFile << "--------------------------------------------------------" << endl;

    if (N % 2 == 0)
    {
        outputFile << "Simpson's 1/3 Rule (N is even):" << endl;
        outputFile << "Integral Estimate: " << integral_1_3 << endl;
    }
    else
    {
        outputFile << "Simpson's 1/3 Rule: N=" << N << " is not even. Calculation skipped." << endl;
    }

    outputFile << "--------------------------------------------------------" << endl;

    if (N % 3 == 0)
    {
        outputFile << "Simpson's 3/8 Rule (N is multiple of 3):" << endl;
        outputFile << "Integral Estimate: " << integral_3_8 << endl;
    }
    else
    {
        outputFile << "Simpson's 3/8 Rule: N=" << N << " is not a multiple of 3. Calculation skipped." << endl;
    }

    outputFile.close();

    cout << "Calculation complete. Results written to output.txt." << endl;

    cout << "\n--- Generated output.txt Content ---" << endl;
    ifstream mockOutput("output.txt");
    string line;
    while (getline(mockOutput, line))
    {
        cout << line << endl;
    }
    mockOutput.close();
}

int main()
{
    cout << fixed << setprecision(10);
    run_simpsons_rules();
    return 0;
}
```

#### Simpson Method Input
```
0.0 3.14159 6
```

#### Simpson Method Output
```
--- Numerical Integration Results (Simpson's Rules) ---
Function: f(x) = sin^5(x) + 4*sin^4(x) + 1
Limits: a = 0.0000000000, b = 3.1415900000
Number of Intervals (N) = 6
Interval width (h) = 0.5235983333
--------------------------------------------------------
Simpson's 1/3 Rule (N is even):
Integral Estimate: 8.9358309110
--------------------------------------------------------
Simpson's 3/8 Rule (N is multiple of 3):
Integral Estimate: 8.6610423801
```

## Ordinary Differential Equation
### RK Method
#### RK Method Theory
The Runge–Kutta (RK4) Method is a powerful numerical technique used to solve first-order ordinary differential equations. It improves accuracy by using four slope evaluations at each step.
It provides high accuracy without requiring higher derivatives.

Given:
dy/dx = f(x, y) with initial condition:
y(x₀) = y₀.

RK4 Formula: 
Let step size be h.
k1 = h * f(x[n], y[n]), 
k2 = h * f(x[n] + h/2, y[n] + k1/2), 
k3 = h * f(x[n] + h/2, y[n] + k2/2), 
k4 = h * f(x[n] + h,   y[n] + k3). 

Updating y: 
y[n+1] = y[n] + (1/6) * (k1 + 2*k2 + 2*k3 + k4).
Updating x: 
x[n+1] = x[n] + h.

#### RK Method Code
```cpp
#include <bits/stdc++.h>
using namespace std;

float dydx(float x,float y)
{
    return (x-y)/2;
}

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    if(!fin)
    {
        cout<<"input.txt not found"<<endl;
        return 0;
    }

    float x0,y0,x,h;
    fin>>x0>>y0>>x>>h;
    float X0=x0;

    cout<<"Input Read From File:"<<endl;
    cout<<"x0 = "<<x0<<endl;
    cout<<"y0 = "<<y0<<endl;
    cout<<"x  = "<<x<<endl;
    cout<<"h  = "<<h<<endl;

    int n=(int)((x-x0)/h);
    float y=y0;
    for(int i=1;i<=n;i++)
    {
        float k1=h*dydx(x0,y);
        float k2=h*dydx(x0+0.5*h,y+0.5*k1);
        float k3=h*dydx(x0+0.5*h,y+0.5*k2);
        float k4=h*dydx(x0+h,y+k3);
        y=y+(k1+2*k2+2*k3+k4)/6.0;
        x0+=h;
    }

    cout<<"The value of y at x is: "<<y<<endl;
    fout<<"Initial x0: "<<X0<<endl;
    fout<<"Initial y0: "<<y0<<endl;
    fout<<"Final x: "<<x<<endl;
    fout<<"Step h: "<<h<<endl;
    fout<<"The value of y at x is: "<<y<<endl;

    fin.close();
    fout.close();
    return 0;
}
```

#### RK Method Input
```
0
1
2
0.1
```

#### RK Method Output
```
Initial x0: 0
Initial y0: 1
Final x: 2
Step h: 0.1
The value of y at x is: 1.10364
```

## Numerical Differentiation

Numerical Differentiation is the process of approximating derivatives of a function when the function is available only at a finite number of discrete points instead of as an explicit formula.
In such cases, direct analytical differentiation is not possible, so an interpolating polynomial is constructed from the given data points and then differentiated.

### Differentiation with Forward Interpolation
#### Differentiation with Forward Interpolation Theory

The Forward Interpolation Differentiation Method is used when the value of the derivative is required near the beginning of the data table. This method is based on Newton’s Forward Interpolation Formula, which expresses the function in terms of forward differences.

After constructing the interpolation polynomial, it is differentiated to obtain approximate formulas for the first and second derivatives.

Forward Interpolation Formula: 
f(x) = y0
     + u Δy0
     + u(u−1)/2! Δ²y0
     + u(u−1)(u−2)/3! Δ³y0
     + ... 
where:
u = (x − x0) / h.

First Derivative Formula after differentiating the forward interpolation polynomial and evaluating at x = x0 (u = 0):

f'(x0) ≈ (1/h) [
    Δy0
  − (1/2) Δ²y0
  + (1/3) Δ³y0
  − (1/4) Δ⁴y0
  + ...
].


#### Differentiation with Forward Interpolation Code
```cpp
#include <bits/stdc++.h>
using namespace std;

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int n;
    fin >> n;

    vector<double> x(n), y(n);
    for (int i = 0; i < n; i++) fin >> x[i];
    for (int i = 0; i < n; i++) fin >> y[i];

    double h = x[1] - x[0];

    vector<vector<double>> diff(n, vector<double>(n));
    for (int i = 0; i < n; i++)
        diff[i][0] = y[i];

    for (int j = 1; j < n; j++)
        for (int i = 0; i < n - j; i++)
            diff[i][j] = diff[i + 1][j - 1] - diff[i][j - 1];

    fout << fixed << setprecision(6);
    fout << "Forward Difference Table:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n - i; j++)
            fout << diff[i][j] << " ";
        fout << "\n";
    }

    double derivative = 0;
    for (int i = 1; i < n; i++) {
        double term = diff[0][i];
        for (int j = 1; j < i; j++)
            term *= (double)j / (j + 1);
        derivative += (i % 2 ? term : -term);
    }

    derivative /= h;

    fout << "\nFirst derivative at x = " << x[0]
         << " is " << derivative << endl;

    return 0;
}
```

#### Differentiation with Forward Interpolation Input
```
5
0 1 2 3 4
1 2 4 8 16
```

#### Differentiation with Forward Interpolation Output
```
Forward Difference Table:
1.000000 1.000000 1.000000 1.000000 1.000000 
2.000000 2.000000 2.000000 2.000000 
4.000000 4.000000 4.000000 
8.000000 8.000000 
16.000000 

First derivative at x = 0.000000 is 0.583333
```

### Differentiation with Backward Interpolation
#### Differentiation with Backward Interpolation Theory

The Backward Interpolation Differentiation Method is used when the value of the derivative is required near the end of the data table. This method is based on Newton’s Backward Interpolation Formula, which expresses the function in terms of backward differences.

It provides better accuracy near the upper end of the tabulated values.

Backward Interpolation Formula: 
f(x) = yn
     + u ∇yn
     + u(u+1)/2! ∇²yn
     + u(u+1)(u+2)/3! ∇³yn
     + ... 
where:
u = (x − xn) / h.

First Derivative Formula after differentiating and evaluating at x = xn (u = 0):
f'(xn) ≈ (1/h) [
    ∇yn
  + (1/2) ∇²yn
  + (1/3) ∇³yn
  + (1/4) ∇⁴yn
  + ...
].

#### Differentiation with Backward Interpolation Code
```cpp
#include <bits/stdc++.h>
using namespace std;

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int n;
    fin >> n;

    vector<double> x(n), y(n);
    for (int i = 0; i < n; i++) fin >> x[i];
    for (int i = 0; i < n; i++) fin >> y[i];

    double h = x[1] - x[0];

    vector<vector<double>> diff(n, vector<double>(n));
    for (int i = 0; i < n; i++)
        diff[i][0] = y[i];

    for (int j = 1; j < n; j++)
        for (int i = n - 1; i >= j; i--)
            diff[i][j] = diff[i][j - 1] - diff[i - 1][j - 1];

    fout << fixed << setprecision(6);
    fout << "Backward Difference Table:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++)
            fout << diff[i][j] << " ";
        fout << "\n";
    }

    double derivative = 0;
    for (int i = 1; i < n; i++)
        derivative += diff[n - 1][i] / i;

    derivative /= h;

    fout << "\nFirst derivative at x = " << x[n - 1]
         << " is " << derivative << endl;

    return 0;
}
```

#### Differentiation with Backward Interpolation Input
```
6
1 2 3 4 5 6
1 8 27 64 125 216
```

#### Differentiation with Backward Interpolation Output
```
Backward Difference Table:
1.000000 
8.000000 7.000000 
27.000000 19.000000 12.000000 
64.000000 37.000000 18.000000 6.000000 
125.000000 61.000000 24.000000 6.000000 0.000000 
216.000000 91.000000 30.000000 6.000000 0.000000 0.000000 

First derivative at x = 6.000000 is 108.000000
```




