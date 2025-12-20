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
