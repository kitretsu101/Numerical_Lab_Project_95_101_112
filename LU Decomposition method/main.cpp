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
