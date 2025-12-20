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
