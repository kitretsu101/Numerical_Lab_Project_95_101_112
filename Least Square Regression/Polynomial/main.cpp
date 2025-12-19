#include <bits/stdc++.h>
using namespace std;

void gauss_elimination(vector<vector<double>>& a, int n)
{
    for(int k=0; k<n; k++)
    {
        int mx=k;
        for(int i=k+1; i<n; i++)
            if(fabs(a[i][k])>fabs(a[mx][k])) mx=i;
        swap(a[k],a[mx]);

        if(fabs(a[k][k])<1e-9) continue;

        for(int i=k+1; i<n; i++)
        {
            double f=a[i][k]/a[k][k];
            for(int j=k; j<=n; j++)
                a[i][j]-=f*a[k][j];
        }
    }
}

vector<double> backward_sub(vector<vector<double>>& a, int n)
{
    vector<double>x(n);
    for(int i=n-1; i>=0; i--)
    {
        x[i]=a[i][n];
        for(int j=i+1; j<n; j++)
            x[i]-=a[i][j]*x[j];
        x[i]/=a[i][i];
    }
    return x;
}

int main()
{
    ifstream infile("input.txt");
    ofstream outfile("output.txt");
    stringstream ss;
    int n;
    while(infile >> n) {
        int order;
        infile >> order;
        int mat_size = order + 1;
        vector<double> x(n);
        vector<double> y(n);
        vector<double> sol(n);
        vector<vector<double>> mat(mat_size, vector<double>(mat_size + 1));

        for(int i = 0; i < n; i++) {
            infile >> x[i] >> y[i];
        }

        for(int i = 0; i < mat_size; i++) {
            for(int j = 0; j < mat_size; j++) {
                double sum = 0;
                for(int k = 0; k < n; k++) {
                    sum += pow(x[k], i+j);
                }
                mat[i][j] = sum;
            }
        }

        for(int i = 0; i < mat_size; i++) {
            double sum = 0;
            for(int j = 0; j < n; j++) {
                sum += pow(x[j], i) * y[j];
            }
            mat[i][mat_size] = sum;
        }

        gauss_elimination(mat, mat_size);
        sol = backward_sub(mat, mat_size);


        ss << "y = ";
        cout << "y = ";
        for(int i = 0; i < mat_size; i++) {
            cout << sol[i] << "x^" << i;
            ss << sol[i] << "x^" << i;
            if(i != mat_size - 1) {
                cout << " + ";
                ss << " + ";
            }
            else {
                cout << endl;
                ss << endl;
            }
        }

        outfile << ss.str();
//        for(int i = 0; i < mat_size; i++) {
//            for(int j = 0; j < mat_size + 1; j++) {
//                cout << mat[i][j] << " ";
//            }
//            cout << endl;
//        }


    }
    return 0;
}
