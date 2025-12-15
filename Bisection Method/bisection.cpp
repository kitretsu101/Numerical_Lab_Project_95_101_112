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
