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
