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
- [Solution of Non-Linear Equations](#solution-of-non-linear-equations)
  - [Newton Raphson Method](#newton-raphson-method)
    - [Theory](#newton-raphson-theory)
    - [Code](#newton-raphson-code)
    - [Input](#newton-raphson-input)
    - [Output](#newton-raphson-output)
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
1. Normalize the pivot row
Row[k] = Row[k] / a[k][k].

2. Eliminate all other elements in the pivot column
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

## Solution of Non-Linear Equations
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
where ε is the permissible error / tolarance.

If the process is repeated in a specific interval for a specific step size, all the roots can be determined for tha corresponding polynomial.

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




